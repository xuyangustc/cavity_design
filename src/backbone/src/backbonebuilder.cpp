/*
 * backbonebuilder.cpp
 *
 *  Created on: 2018年1月9日
 *      Author: hyliu
 */

#include "backbone/backbonebuilder.h"
#include "dstl/randomengine.h"
#include "pdbstatistics/phipsidistr.h"
#include "backbone/closealoop.h"
using namespace NSPproteinrep;
using namespace NSPgeometry;
std::vector<BackBoneSite> BackBoneBuilder::buildforwardbackbone(int length,
		const BackBoneSite & nflankingsite, const std::vector<std::pair<int,int>> & helixregions,
		const std::vector<std::pair<int,int>> &strandregions, const std::set<int> & cissites){
	auto rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	BackBoneSite bsn=nflankingsite;
	std::vector<BackBoneSite> chain(length);
	for(int i=0;i<chain.size();++i){
		bool helix=false;
		bool strand=false;
		for(auto &r:helixregions){
			if(i>=r.first &&i<r.first+r.second) {
				helix=true;break;
			}
		}
		for(auto &r:strandregions){
			if(i>=r.first &&i<r.first+r.second) {
				strand=true;break;
			}
		}
		double phi,psi;
#pragma omp critical
		{
			if(helix) NSPpdbstatistics::PhiPsiDistr::helixdistr().randomphipsi(rng,&phi,&psi);
			else if(strand)NSPpdbstatistics::PhiPsiDistr::stranddistr().randomphipsi(rng,&phi,&psi);
			else NSPpdbstatistics::PhiPsiDistr::mixcoildistr().randomphipsi(rng,&phi,&psi);
		}
		BackBoneSite *bsp=&bsn;
		if(i>0) bsp=&chain[i-1];
		bool cispep=(cissites.find(i) != cissites.end());
		genbackbonesite(bsp,cispep,phi,psi,&chain[i]);
		if(i==0) {
			chain[0].chainid='A';
			chain[0].resseq=0;
		}
		chain[i].resname="ALA";
		if(cispep) chain[i].resname="PRO";
		if(helix) chain[i].sscode='H';
		else chain[i].sscode='C';
	}
	return chain;
}
std::vector<BackBoneSite> BackBoneBuilder::buildbackwardbackbone(int length,
		const BackBoneSite & cflankingsite, const std::vector<std::pair<int,int>> & helixregions,
		const std::vector<std::pair<int,int>> & strandregions, const std::set<int> & cissites){
	auto rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	BackBoneSite bsc=cflankingsite;
	std::vector<BackBoneSite> chain(length);
	for(int i=chain.size()-1;i>=0;--i){
		bool helix=false;
		bool strand=false;
		for(auto &r:helixregions){
			if(i>=r.first &&i<r.first+r.second) {
				helix=true;break;
			}
		}
		for(auto &r:strandregions){
			if(i>=r.first &&i<r.first+r.second) {
				strand=true;break;
			}
		}
		double phi,psi;
#pragma omp critical
		{
			if(helix) NSPpdbstatistics::PhiPsiDistr::helixdistr().randomphipsi(rng,&phi,&psi);
			else if(strand)NSPpdbstatistics::PhiPsiDistr::stranddistr().randomphipsi(rng,&phi,&psi);
			else NSPpdbstatistics::PhiPsiDistr::mixcoildistr().randomphipsi(rng,&phi,&psi);
		}
		BackBoneSite *bsm=&bsc;
		if(i<chain.size()-1) bsm=&chain[i+1];
		bool cispep=(cissites.find(i) != cissites.end());
		double omiga=180.0;
		if(cispep) omiga=0.0;
		NSPproteinrep::genprevbackbonesite(bsm,omiga,psi,phi,
				&chain[i]);
		if(i==chain.size()-1) {
			chain[i].chainid='A';
			chain[i].resseq=i;
			chain[i].resname="ALA";
		}
		if(helix) chain[i].sscode='H';
		else chain[i].sscode='C';
	}
	return chain;
}

std::vector<BackBoneSite> BackBoneBuilder::buildstrandat(int length,NSPgeometry::XYZ r0,
		NSPgeometry::XYZ direction,bool forward){
	BackBoneSite nsite;
	genbackbonesite(nullptr, false, 0.0,0.0, &nsite);
	std::vector<std::pair<int,int>> strandregions;
	strandregions.push_back(std::make_pair(0,length));
	std::vector<BackBoneSite> chain=buildforwardbackbone(length,nsite, std::vector<std::pair<int,int>>(),strandregions,std::set<int> ());
	movechainto(r0,direction,forward,chain);
	return chain;
}
std::vector<BackBoneSite> BackBoneBuilder::buildhelixat(int length,NSPgeometry::XYZ r0,
		NSPgeometry::XYZ direction,bool forward){
	BackBoneSite nsite;
	genbackbonesite(nullptr, false, 0.0,0.0, &nsite);
	std::vector<std::pair<int,int>> helixregions;
	helixregions.push_back(std::make_pair(0,length));
	std::vector<BackBoneSite> chain=buildforwardbackbone(length,nsite,
			helixregions,std::vector<std::pair<int,int>>(),std::set<int> ());
	movechainto(r0,direction,forward,chain);
	return chain;
}

void BackBoneBuilder::movechainto(NSPgeometry::XYZ r0,NSPgeometry::XYZ direction,
		bool forward, std::vector<BackBoneSite> & chain){
	NSPgeometry::XYZ rhead;
	NSPgeometry::XYZ rheadtotail;
	if (forward){
		rhead=chain[0].cacrd();
		rheadtotail=chain.back().cacrd()-rhead;
	} else {
		rhead=chain.back().cacrd();
		rheadtotail=chain[0].cacrd()-rhead;
	}
	NSPgeometry::XYZ trans=r0-rhead;
	for (auto &bs:chain){
		bs.translate(trans);
	}
	NSPgeometry::XYZ axis=NSPgeometry::cross(rheadtotail,direction);
	if(axis.squarednorm()>0.0){
		double angle=acos(NSPgeometry::dot(rheadtotail,direction)/(sqrt(rheadtotail.squarednorm()*direction.squarednorm())));
		NSPgeometry::Rotation rot(NSPgeometry::QuaternionCrd(axis,angle,1.0),r0);
		for(auto &bs:chain){
			bs.rotate(rot);
		}
	}
}
std::vector<std::shared_ptr<std::vector<BackBoneSite>>> BackBoneBuilder::buildlinkers(int length,
		const BackBoneSite &nflankingsite,
		const BackBoneSite &cflankingsite,
		const std::vector<std::pair<int,int>> & helixregions,
		const std::vector<std::pair<int,int>> & strandregions,
		const std::set<int> & cissites, int n12result){
	if(length==1) {
		bool cis0=cissites.find(0)!=cissites.end();
		double dnn2=4.0*4.0;
		XYZ nn=nflankingsite.getnextN();
		XYZ nnn=cflankingsite.ncrd();
		if((nn-nnn).squarednorm()>dnn2) return std::vector<std::shared_ptr<std::vector<BackBoneSite>>>();
		return buildlinker1(nflankingsite,cflankingsite,n12result,cis0);
	}
	if(length==2) {
		bool cis0=cissites.find(0)!=cissites.end();
		bool cis1=cissites.find(1)!=cissites.end();
		return buildlinker2(nflankingsite,cflankingsite,n12result,cis0,cis1);
	}
	std::vector<BackBoneSite> temploop(length+2);
	temploop[0]=nflankingsite;
	temploop.back()=cflankingsite;
	std::vector<BackBoneSite> middle=buildforwardbackbone(length,nflankingsite, helixregions,strandregions,cissites);
	temploop[1]=middle[0];
	std::vector<BackBoneSite> revmiddle=buildbackwardbackbone(1,cflankingsite, helixregions,strandregions,cissites);
	temploop[length]=revmiddle[0];
	std::vector<double> torsions;
	for(auto & bs:middle) {
		torsions.push_back(bs.phi());
		torsions.push_back(bs.psi());
		torsions.push_back(bs.omiga());
	}
	std::set<int> fixedsites;
	for(int i=0;i<length;++i){
		bool fixedi=false;
		for(auto &h:helixregions) if(i>=h.first && i<h.first+h.second) fixedi=true;
		for(auto &s:strandregions) if(i>=s.first && i<s.first+s.second) fixedi=true;
		if(fixedi) fixedsites.insert(i+1);
	}
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>> solutions;
#pragma omp critical
	{
		CloseALoop closer(temploop,1,length,torsions);
		solutions=closer.getsolutions(fixedsites);
	}
	return solutions;
}



std::vector<std::shared_ptr<std::vector<BackBoneSite>>> BackBoneBuilder::buildlinker1(
		const BackBoneSite &nflankingsite, const BackBoneSite &cflankingsite, int nr, bool cis0) {
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>> rts;
	NSPgeometry::XYZ nnn = cflankingsite.ncrd();
	double dis2 = 0.2*0.2;
	for(int i=0;i<nr;i++) {
		double stphi=NSPdstl::RandomEngine<>::getinstance().realrng(-180.0,180.0)();
		double stpsi=NSPdstl::RandomEngine<>::getinstance().realrng(-180.0,180.0)();
		BackBoneSite bs = nflankingsite.NextSite(cis0, stphi, stpsi);
		NSPgeometry::XYZ n = bs.ncrd();
		NSPgeometry::XYZ ca = bs.cacrd();
		NSPgeometry::XYZ c = bs.ccrd();
		NSPgeometry::XYZ o = bs.ocrd();
		NSPgeometry::XYZ nn = bs.getnextN();
		NSPgeometry::XYZ axis = ca -n;
		double angle = 5.0;
		NSPgeometry::QuaternionCrd qc(axis,angle);
		NSPgeometry::Rotation rot(qc,n);
		for(int j=0;j<72;j++) {
			rot.apply(&c);
			rot.apply(&o);
			rot.apply(&nn);
			LocalFrame lf0=make_localframe(ca,c,nn);
			XYZ lnn = lf0.global2localcrd(nn);
			LocalFrame lf1=make_localframe(ca,c,nnn);
			XYZ gnn = lf1.local2globalcrd(lnn);
			double d2=(gnn-nnn).squarednorm();
			//std::cout <<"distance2:\t" <<d2 <<std::endl;
			if(d2>dis2) continue;
			nn = gnn;
			o = lf0.global2localcrd(o);
			o = lf1.local2globalcrd(o);
			bs.data_[BackBoneSite::CCRD] = c.x_;
			bs.data_[BackBoneSite::CCRD+1] = c.y_;
			bs.data_[BackBoneSite::CCRD+2] = c.z_;
			bs.data_[BackBoneSite::OCRD] = o.x_;
			bs.data_[BackBoneSite::OCRD+1] = o.y_;
			bs.data_[BackBoneSite::OCRD+2] = o.z_;
			rts.push_back(std::shared_ptr < std::vector<BackBoneSite> > (new std::vector<BackBoneSite>));
			rts.back()->push_back(bs);
			break;
		}
	}
	return rts;
}

std::vector<std::shared_ptr<std::vector<BackBoneSite>>> BackBoneBuilder::buildlinker2(
		const BackBoneSite &nflankingsite, const BackBoneSite &cflankingsite, int nr, bool cis0, bool cis1) {
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>> rts;
	XYZ nnn = cflankingsite.ncrd();
	double dnn2=4.0*4.0;
	for(int i=0;i<nr;i++) {
		double stphi0=NSPdstl::RandomEngine<>::getinstance().realrng(-180.0,180.0)();
		double stpsi0=NSPdstl::RandomEngine<>::getinstance().realrng(-180.0,180.0)();
		BackBoneSite bs0 = nflankingsite.NextSite(cis0, stphi0, stpsi0);
		NSPgeometry::XYZ n0 = bs0.ncrd();
		NSPgeometry::XYZ ca0 = bs0.cacrd();
		NSPgeometry::XYZ c0 = bs0.ccrd();
		NSPgeometry::XYZ o0 = bs0.ocrd();
		NSPgeometry::XYZ nn0 = bs0.getnextN();
		NSPgeometry::XYZ axis = ca0 - n0;
		double angle = 5.0;
		NSPgeometry::QuaternionCrd qcphi(axis,angle);
		NSPgeometry::Rotation rotphi(qcphi,ca0);
		axis = c0 - ca0;
		NSPgeometry::QuaternionCrd qcpsi(axis,angle);
		NSPgeometry::Rotation rotpsi(qcpsi,ca0);
		for(int h=0;h<72;h++) {
			rotphi.apply(&c0);
			rotphi.apply(&o0);
			rotphi.apply(&nn0);
			bool find = false;
			for(int s=0;s<72;s++) {
				rotpsi.apply(&o0);
				rotpsi.apply(&nn0);
				double d2 = (nnn-nn0).squarednorm();
				if(d2>dnn2) continue;
				bs0.data_[BackBoneSite::CCRD] = c0.x_;
				bs0.data_[BackBoneSite::CCRD+1] = c0.y_;
				bs0.data_[BackBoneSite::CCRD+2] = c0.z_;
				bs0.data_[BackBoneSite::OCRD] = o0.x_;
				bs0.data_[BackBoneSite::OCRD+1] = o0.y_;
				bs0.data_[BackBoneSite::OCRD+2] = o0.z_;
				std::vector<std::shared_ptr<std::vector<BackBoneSite>>> r1 = buildlinker1(bs0,cflankingsite,1,cis1);
				if(r1.empty()) continue;
				rts.push_back(std::shared_ptr < std::vector<BackBoneSite> > (new std::vector<BackBoneSite>));
				rts.back()->push_back(bs0);
				rts.back()->push_back(r1[0]->at(0));
				find = true;
				break;
			}
			if(find) break;
		}
	}
	return rts;
}













NSPgeometry::Rotation BackBoneBuilder::transform(NSPgeometry::XYZ fixed1, NSPgeometry::XYZ fixed2,
		NSPgeometry::XYZ moved1, NSPgeometry::XYZ moved2, NSPgeometry::XYZ &translated) {
	translated = fixed1 - moved1;
	moved1 = moved1 + translated;
	moved2 = moved2 + translated;
	double ang = NSPgeometry::angle(fixed2,fixed1,moved2) * 180.0/3.14159265;
	NSPgeometry::XYZ md = moved2 - moved1;
	NSPgeometry::XYZ fd = fixed2 - fixed1;
	NSPgeometry::XYZ axis = NSPgeometry::cross(md,fd);
	NSPgeometry::QuaternionCrd qc(axis,ang);
	return NSPgeometry::Rotation(qc,fixed1);
}

std::vector<BackBoneSite> BackBoneBuilder::chainrotation(const BackBoneSite &bs, const std::vector<BackBoneSite> &chain) {
	std::vector<std::pair<int,int>> helix;
	std::vector<std::pair<int,int>> strand;
	std::set<int> cis;
	std::vector<BackBoneSite> ch = buildforwardbackbone(1,bs, helix,strand,cis);
	std::vector<BackBoneSite> newchain;
	NSPgeometry::XYZ fixedn(ch[0].data_[BackBoneSite::NCRD],ch[0].data_[BackBoneSite::NCRD+1],ch[0].data_[BackBoneSite::NCRD+2]);
	NSPgeometry::XYZ fixedca(ch[0].data_[BackBoneSite::CACRD],ch[0].data_[BackBoneSite::CACRD+1],ch[0].data_[BackBoneSite::CACRD+2]);
	NSPgeometry::XYZ movedn(chain[0].data_[BackBoneSite::NCRD],
			chain[0].data_[BackBoneSite::NCRD+1],chain[0].data_[BackBoneSite::NCRD+2]);
	NSPgeometry::XYZ movedca(chain[0].data_[BackBoneSite::CACRD],
			chain[0].data_[BackBoneSite::CACRD+1],chain[0].data_[BackBoneSite::CACRD+2]);
	NSPgeometry::XYZ translated;
	NSPgeometry::Rotation rt=transform(fixedn,fixedca,movedn,movedca,translated);
	for(int i=0;i<chain.size();i++) {
		BackBoneSite bs=chain[i];
		for(int k=BackBoneSite::NCRD;k<=BackBoneSite::OCRD;k+=3) {
			NSPgeometry::XYZ m(bs.data_[k],bs.data_[k+1],bs.data_[k+2]);
			m = m + translated;
			rt.apply(&m);
			bs.data_[k] = m.x_;
			bs.data_[k+1] = m.y_;
			bs.data_[k+2] = m.z_;
		}
		newchain.push_back(bs);
	}
	return newchain;
}

std::vector<BackBoneSite> BackBoneBuilder::getbackbone(const std::vector<std::pair<std::pair<int,int>,char>> &label,
		const std::set<int> &cissite, const std::vector<std::vector<BackBoneSite>>&fixed, const BackBoneSite &bsn) {
	std::vector<BackBoneSite> chain;
	std::vector<std::pair<int,int>> helix;
	std::vector<std::pair<int,int>> strand;
	//std::vector<std::pair<int,int>> lp;
	//std::vector<std::pair<int,int>> fx;
	int istart=0;
	std::set<int> cis;
	int lth=0;
	int nfix=0;
	BackBoneSite bsstart=bsn;
	for(int i=0;i<label.size();i++) {
		if(i==label.size()-1) {
			if(label[i].second=='f') {
				std::vector<BackBoneSite> ch = buildforwardbackbone(lth,bsstart,helix,strand,cis);
				for(BackBoneSite &bs:ch) chain.push_back(bs);
				std::vector<BackBoneSite> newch=chainrotation(ch.back(),fixed[nfix]);
				for(BackBoneSite &bs:newch) chain.push_back(bs);
				helix.clear();
				strand.clear();
				cis.clear();
				lth=0;
				nfix++;
				bsstart=chain.back();
			} else {
				if(label[i].second=='a') helix.push_back(label[i].first);
				else if(label[i].second=='b') strand.push_back(label[i].first);
				for(const auto & c:cissite) {
					if(c>=label[i].first.first && c<label[i].first.first+label[i].first.second)
						cis.insert(c-istart);
				}
				lth += label[i].first.second;
				std::vector<BackBoneSite> ch = buildforwardbackbone(lth,bsstart,helix,strand,cis);
				for(BackBoneSite &bs:ch) chain.push_back(bs);
			}
		} else if(label[i].second=='f') {
			std::vector<BackBoneSite> ch = buildforwardbackbone(lth,bsstart,helix,strand,cis);
			for(BackBoneSite &bs:ch) chain.push_back(bs);
			std::vector<BackBoneSite> newch=chainrotation(ch.back(),fixed[nfix]);
			for(BackBoneSite &bs:newch) chain.push_back(bs);
			helix.clear();
			strand.clear();
			cis.clear();
			istart=chain.size();
			lth=0;
			nfix++;
			bsstart=chain.back();
		} else if(label[i].second=='a') {
			helix.push_back(label[i].first);
			for(const auto & c:cissite) {
				if(c>=label[i].first.first && c<label[i].first.first+label[i].first.second)
					cis.insert(c-istart);
			}
			lth += label[i].first.second;
		} else if(label[i].second=='b') {
			strand.push_back(label[i].first);
			for(const auto & c:cissite) {
				if(c>=label[i].first.first && c<label[i].first.first+label[i].first.second)
					cis.insert(c-istart);
			}
			lth += label[i].first.second;
		} else if(label[i].second=='l') {
			for(const auto & c:cissite) {
				if(c>=label[i].first.first && c<label[i].first.first+label[i].first.second)
					cis.insert(c-istart);
			}
			lth += label[i].first.second;
		}
	}
	return chain;
}

std::vector<std::shared_ptr<std::vector<BackBoneSite>>> BackBoneBuilder::buildlinkers(
		const BackBoneSite & nflankingsite, const BackBoneSite & cflankingsite,
		const std::vector<std::pair<std::pair<int,int>,char>> &label,
		const std::vector<std::vector<BackBoneSite>> & fixedconf, const std::set<int> & cissites) {
	std::vector<BackBoneSite> chain=getbackbone(label, cissites, fixedconf, nflankingsite);
	/*{
		std::shared_ptr<std::vector<BackBoneSite>> ch1=std::make_shared<std::vector<BackBoneSite>>(chain);
		std::vector<std::shared_ptr<std::vector<BackBoneSite>>> sol;
		sol.push_back(ch1);
		return sol;
	}*/
	std::vector<std::pair<int,int>> helix;
	std::vector<std::pair<int,int>> strand;
	std::set<int> cis;
	std::vector<BackBoneSite> revmiddle=buildbackwardbackbone(1,cflankingsite, helix,strand,cis);
	std::vector<BackBoneSite> temploop{nflankingsite};
	for(int i=0;i<chain.size()-1;i++) temploop.push_back(chain[i]);
	temploop.push_back(revmiddle[0]);
	temploop.push_back(cflankingsite);
	CloseALoop closer(temploop, 1, chain.size(), chain);
	std::set<int> fixedsites;
	for(const auto & l:label) {
		if(l.second=='l') continue;
		int st=l.first.first;
		int lth=l.first.second;
		for(int i=st;i<st+lth;i++) fixedsites.insert(i);
	}
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>> solutions=closer.getsolutions(fixedsites);
	return solutions;
}

NSPgeometry::Rotation BackBoneBuilder::rot(NSPgeometry::XYZ fixed1, NSPgeometry::XYZ fixed2,
		NSPgeometry::XYZ fixed3, NSPgeometry::XYZ moved3) {
	NSPgeometry::XYZ ax = fixed1-fixed2;
	NSPgeometry::XYZ axf = fixed3-fixed2;
	NSPgeometry::XYZ axm = moved3-fixed2;
	NSPgeometry::XYZ p1=cross(ax,axf);
	NSPgeometry::XYZ p2=cross(ax,axm);
	double cos=p1.x_*p2.x_+p1.y_*p2.y_+p1.z_*p2.z_;
	cos = cos / sqrt(p1.squarednorm()*p2.squarednorm());
	cos = acos(cos);
	cos = cos *180.0 / 3.141592657;
	NSPgeometry::XYZ axis = fixed2-fixed1;
	NSPgeometry::QuaternionCrd qc(axis,cos);
	NSPgeometry::Rotation rt(qc,fixed1);

	NSPgeometry::XYZ m3 = moved3;
	rt.apply(&m3);

	if((fixed3-m3).squarednorm()>0.01) {
		//std::cout <<(fixed3-m3).squarednorm() <<std::endl;
		axis = fixed1-fixed2;
		NSPgeometry::QuaternionCrd qc1(axis,cos);
		rt = NSPgeometry::Rotation(qc1,fixed1);
		m3 = moved3;
		rt.apply(&m3);
		//std::cout <<'\t' <<(fixed3-m3).squarednorm() <<std::endl;
		//assert((fixed3-m3).squarednorm()<0.01);
	}
	return rt;
}

std::vector<FullSite> BackBoneBuilder::addsidechain(
		const std::vector<BackBoneSite>&bsch, const std::vector<FullSite>&orifs) {
	std::vector<FullSite> newfs=orifs;
	//for(const BackBoneSite &bs:bsch) newfs.push_back(make_fullsite(bs));
	for(int i=0;i<orifs.size();i++) {
		FullSite fs=orifs[i];
		std::map<std::string,NSPgeometry::XYZ> crds=fs.getcrds();
		NSPgeometry::XYZ translated;
		NSPgeometry::Rotation rt=transform(bsch[i].ncrd(),
				bsch[i].cacrd(),fs.getcrd("N"),fs.getcrd("CA"),translated);
		for(auto &c:crds) {
			c.second = c.second+translated;
			rt.apply(&(c.second));
		}
		/*{
			NSPgeometry::XYZ nm1=crds.at("N");
			NSPgeometry::XYZ cam1=crds.at("CA");
			NSPgeometry::XYZ nf1=bsch[i].ncrd();
			NSPgeometry::XYZ caf1=bsch[i].cacrd();
			double dn1 = (nm1-nf1).squarednorm();
			double dca1 = (cam1-caf1).squarednorm();
			std::cout <<dn1 <<'\t' <<dca1 <<std::endl;
			//assert(dn1<0.00001);
			//assert(dca1<0.00001);
		}*/
		//fs.changecrd(crds);
		rt=rot(bsch[i].ncrd(),bsch[i].cacrd(),bsch[i].ccrd(),crds.at("C"));
		for(auto &c:crds) {
			rt.apply(&(c.second));
		}
		/*{
			NSPgeometry::XYZ nm2=crds.at("N");
			NSPgeometry::XYZ cam2=crds.at("CA");
			NSPgeometry::XYZ cm2=crds.at("C");
			NSPgeometry::XYZ om2=crds.at("O");
			NSPgeometry::XYZ nf2=bsch[i].ncrd();
			NSPgeometry::XYZ caf2=bsch[i].cacrd();
			NSPgeometry::XYZ cf2=bsch[i].ccrd();
			NSPgeometry::XYZ of2=bsch[i].ocrd();
			double dn2 = (nm2-nf2).squarednorm();
			double dca2 = (cam2-caf2).squarednorm();
			double dc2 = (cm2-cf2).squarednorm();
			double do2 = (om2-of2).squarednorm();
			std::cout <<dn2 <<'\t' <<dca2 <<'\t' <<dc2 <<'\t' <<do2 <<std::endl;
			//assert(dn2.squarednorm()<0.00001);
			//assert(dca2.squarednorm()<0.00001);
			//assert(dc2.squarednorm()<0.00001);
			//assert(do2.squarednorm()<0.00001);
		}*/
		fs.changecrd(crds);
		newfs[i] = fs;
	}
	return newfs;
}











BackBoneBuilder::Region::Region(const std::string &label, const std::set<int> &cis) {
	length=label.size();
	char cl=label[0];
	int startsite=0;
	int endsite=0;
	for(int i=1;i<label.size();i++) {
		if(cl==label[i]) continue;
		startsite=endsite;
		endsite=i;
		if(cl=='a') helixregion.push_back({startsite,endsite});
		else if(cl=='b') strandregion.push_back({startsite,endsite});
		cl=label[i];
	}
	if(cl=='a') helixregion.push_back({startsite,endsite});
	else if(cl=='b') strandregion.push_back({startsite,endsite});
	cissites = cis;
}

std::vector<std::shared_ptr<std::vector<BackBoneSite>>> BackBoneBuilder::buildlinkers(
		const BackBoneSite & nflankingsite, const BackBoneSite & cflankingsite, const std::string & label,
		const std::vector<std::vector<BackBoneSite>> & fixedconf, const std::set<int> & cissites) {
	//label l-loop, a-helix,  b-strand, f-fixedsite
	std::vector<BackBoneSite> linker;
	linker.push_back(nflankingsite);
	std::string lb;
	std::set<int> cis;
	int startsite=0, endsite=0;
	int nfixed=0;
	for(int i=0;i<label.size();i++) {
		if(label[i]!='f') {
			lb.push_back(label[i]);
			continue;
		}
		startsite=endsite;
		endsite=i;
		//if(startsite == endsite) continue;
		for(auto c:cissites) if(c>=startsite && c<endsite) cis.insert(c-startsite);
		//build loop
		if(startsite != endsite) {
			Region rg(lb,cis);
			//length is 1 longer than origenal length, to link fixed part
			std::vector<BackBoneSite> middle=buildforwardbackbone(endsite-startsite+1,
					linker.back(),rg.helixregion,rg.strandregion,rg.cissites);
			for(BackBoneSite &bs:middle) linker.push_back(bs);
		}
		//place fixedregion
		std::vector<BackBoneSite> bss=fixedconf[nfixed];
		NSPgeometry::XYZ fixedn(linker.back().data_[BackBoneSite::NCRD],
				linker.back().data_[BackBoneSite::NCRD+1],linker.back().data_[BackBoneSite::NCRD+2]);
		NSPgeometry::XYZ fixedca(linker.back().data_[BackBoneSite::CACRD],
				linker.back().data_[BackBoneSite::CACRD+1],linker.back().data_[BackBoneSite::CACRD+2]);
//		NSPgeometry::XYZ fixedc(linker.back().data_[BackBoneSite::CCRD],
//				linker.back().data_[BackBoneSite::CCRD+1],linker.back().data_[BackBoneSite::CCRD+2]);
		NSPgeometry::XYZ movedn(bss[0].data_[BackBoneSite::NCRD],
				bss[0].data_[BackBoneSite::NCRD+1],bss[0].data_[BackBoneSite::NCRD+2]);
		NSPgeometry::XYZ movedca(bss[0].data_[BackBoneSite::CACRD],
				bss[0].data_[BackBoneSite::CACRD+1],bss[0].data_[BackBoneSite::CACRD+2]);
//		NSPgeometry::XYZ movedc(bss[0].data_[BackBoneSite::CCRD],
//				bss[0].data_[BackBoneSite::CCRD+1],bss[0].data_[BackBoneSite::CCRD+2]);
		NSPgeometry::XYZ translated;
		NSPgeometry::Rotation rt = transform(fixedn,fixedca,movedn,movedca,translated);
		linker.resize(linker.size()-1);
		for(int j=0;j<bss.size();j++) {
			for(int k=BackBoneSite::NCRD;k<=BackBoneSite::OCRD;k+=3) {
				NSPgeometry::XYZ m(bss[j].data_[k],bss[j].data_[k+1],bss[j].data_[k+2]);
				m = m + translated;
				rt.apply(&m);
				bss[j].data_[k] = m.x_;
				bss[j].data_[k+1] = m.y_;
				bss[j].data_[k+2] = m.z_;
			}
			linker.push_back(bss[j]);
		}
		i += fixedconf[nfixed].size()-1;
		lb.clear();
		cis.clear();
		nfixed++;
		startsite=endsite;
		endsite=i+1;
	}
	startsite=endsite;
	endsite=label.size();
	for(auto c:cissites) if(c>=startsite && c<endsite) cis.insert(c-startsite);
	//build loop
	if(startsite != endsite) {
		Region rg(lb,cis);
		std::vector<BackBoneSite> middle=buildforwardbackbone(endsite-startsite,
				linker.back(),rg.helixregion,rg.strandregion,rg.cissites);
		for(BackBoneSite &bs:middle) linker.push_back(bs);
	}
	std::vector<double> torsions;
	for(int i=1;i<linker.size();i++) {
		if(i==linker.size()-1) {
			torsions.push_back(linker[i].phi(linker[i-1]));
			torsions.push_back(360.0);
			torsions.push_back(180.0);
			continue;
		}
		torsions.push_back(linker[i].phi(linker[i-1]));
		torsions.push_back(linker[i].psi(linker[i+1]));
		torsions.push_back(linker[i].omiga(linker[i+1]));
	}
	std::vector<BackBoneSite> revmiddle=buildbackwardbackbone(1,cflankingsite,
			std::vector<std::pair<int,int>>(),std::vector<std::pair<int,int>>(),std::set<int>());
	linker.back() = revmiddle[0];
	linker.push_back(cflankingsite);
	assert(linker.size()==label.size()+2);
	CloseALoop closer(linker,1,label.size(),torsions);
	std::set<int> fixedsites;
	for(int i=0;i<label.size();++i){
		if(label[i]!='l') fixedsites.insert(i+1);
	}
	std::vector<std::shared_ptr<std::vector<BackBoneSite>>> solutions=closer.getsolutions(fixedsites);
	return solutions;
}

/*
std::vector<FullSite> BackBoneBuilder::addsidechain(
		const std::vector<BackBoneSite>&bsch, const std::map<int,std::vector<FullSite>>&orifs) {
	std::vector<FullSite> newfs;
	for(const BackBoneSite &bs:bsch) newfs.push_back(make_fullsite(bs));
	for(const auto &of:orifs) {
		for(int i=0;i<of.second.size();i++) {
			FullSite fs=of.second[i];
			NSPgeometry::XYZ translated;
			NSPgeometry::Rotation rt=transform(newfs[of.first+i].getcrd("N"),
					newfs[of.first+i].getcrd("C"),fs.getcrd("N"),fs.getcrd("C"),translated);
			std::map<std::string,NSPgeometry::XYZ> crds=fs.getcrds();
			for(auto &c:crds) {
				c.second = c.second+translated;
				rt.apply(&(c.second));
			}
			fs.changecrd(crds);
			newfs[of.first+i] = fs;
		}
	}
	return newfs;
}*/












