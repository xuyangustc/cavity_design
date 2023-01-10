/*
 * shape.cpp
 *
 *  Created on: Aug 18, 2021
 *      Author: xuyang
 */
#include "allatom/shape.h"
#include "allatom/basic_info.h"

using namespace NSPproteinrep;
using namespace NSPgeometry;
using namespace NSPallatom;

/*
 * ch from lf0 to lf1
 */
void Par2Helix::ChainCrdCovert(std::vector<Residue>&ch, LocalFrame &lf0, LocalFrame &lf1) {
	for(Residue &r:ch) {
		for(auto &p:r.rds()) {
			p.second.crd = lf0.global2localcrd(p.second.crd);
			p.second.crd = lf1.local2globalcrd(p.second.crd);
		}
	}
}
void Par2Helix::ChainRotate(std::vector<Residue>&ch, Rotation &rt) {
	for(Residue &r:ch) {
		for(auto &p:r.rds()) {
			rt.apply(&(p.second.crd));
		}
	}
}
void Par2Helix::ChainTranslate(std::vector<Residue>&ch, XYZ &trans) {
	for(Residue &r:ch) {
		for(auto &p:r.rds()) {
			p.second.crd = p.second.crd + trans;
		}
	}
}

std::vector<XYZ> Par2Helix::HelixTerminal(const std::vector<Residue>&rs) {
	int nres=7;
	if(rs.size()<nres+1) {
		std::cout <<"length of initial helix must be more than " <<nres+1 <<std::endl;
		exit(1);
	}
	XYZ nt(0.0,0.0,0.0);
	XYZ ct(0.0,0.0,0.0);
	int sz=rs.size();
	for(int i=0;i<nres;i++) {
		nt = nt + rs[i].rds().at("CA").crd;
		ct = ct + rs[sz-i-1].rds().at("CA").crd;
	}
	nt = nt /(double)nres;
	ct = ct /(double)nres;
	return {nt,ct};
}
XYZ Par2Helix::HelixAxis(const std::vector<Residue>&rs) {
	std::vector<XYZ> ts=HelixTerminal(rs);
	return ts[1]-ts[0];
}
/*
 * rotating axis from N 2 C
 */
void Par2Helix::ChainRotation(std::vector<Residue>&rs, double angle) {
	std::vector<XYZ> ts=HelixTerminal(rs);
	XYZ ax = ts[1] - ts[0];
	QuaternionCrd qc(ax,angle);
	Rotation rt(qc,ts[0]);
	ChainRotate(rs,rt);
}

std::vector<Residue> Par2Helix::MakeChain(std::vector<std::vector<double>>&dihs,
		std::vector<std::string>&seq) {
	assert(dihs.size()==seq.size());
	std::vector<BackBoneSite> chain(dihs.size());
	genbackbonesite(nullptr, false, dihs[0][0], dihs[0][1], &chain[0]);
	for(int i=1;i<dihs.size();i++) {
		genbackbonesite(&chain[i-1], false, dihs[i][0], dihs[i][1], &chain[i]);
	}
	std::vector<Residue> rs;
	for(int i=0;i<chain.size();i++) rs.push_back(Residue(chain[i],seq[i]));
	return rs;
}
std::vector<Residue> Par2Helix::InitialHelix(int lth, bool toup) {
	std::vector<double> phipsi{-57.0,-47.0};
	std::vector<std::vector<double>> dihs(lth,phipsi);
	std::vector<std::string> seq(lth,"GLY");
	std::vector<Residue> rs = MakeChain(dihs,seq);

	std::vector<XYZ> ts=HelixTerminal(rs);
	XYZ xy=rs[7].rds().at("CA").crd;
	LocalFrame lf0 = make_localframe(ts[0], ts[1], xy);
	XYZ localo(0.0,0.0,0.0);
	XYZ localx(0.0,0.0,1.0);
	XYZ localxy(0.0,1.0,1.0);
	if(!toup) {
		localx.z_ = -localx.z_;
		localxy.z_ = -localxy.z_;
	}
	LocalFrame lf1 = make_localframe(localo, localx, localxy);
	ChainCrdCovert(rs,lf0,lf1);
	return rs;
}

std::vector<Residue> Par2Helix::HelixGen() {
	std::vector<Residue> ch=InitialHelix(lth,toup);
	std::vector<XYZ> ts=HelixTerminal(ch);
	XYZ trans=axisA-ts[0];
	trans.z_ = 0.0;
	ChainTranslate(ch,trans);
	double Nhight=hightA;
	if(toup) Nhight += transdis;
	else Nhight -= transdis;
	trans.z_ = Nhight-ch[0].rds().at("CA").crd.z_;
	trans.x_ = 0.0;
	trans.y_ = 0.0;
	ChainTranslate(ch,trans);
	XYZ nres_crd = ch[nres].rds().at("CA").crd;
	XYZ localx = axisA;
	localx.z_ = nres_crd.z_;
	XYZ localy = localx;
	localy.z_ += 1.0;
	LocalFrame lf0 = make_localframe(localx, localy, nres_crd);
	XYZ localxy=axisB;
	localxy.z_ = nres_crd.z_;
	LocalFrame lf1 = make_localframe(localx, localy, localxy);
	ChainCrdCovert(ch,lf0,lf1);
	//std::cout <<22 <<std::endl;
	ChainRotation(ch,rotangle);
	//std::cout <<33 <<std::endl;
	return ch;
}










/*
 * helix must be vertical to XOY
 * fixedaxis: axis, Z-coordinate can be any value
 * helix: query helix
 * nres: return value based on (nres)th residue
 * return rotating angle of helix(axis from N 2 C)
 */
double Shape::Angle2Axis(XYZ fixedaxis, const std::vector<Residue>&helix, int nres) {
	XYZ fixedn = helix[nres].rds().at("CA").crd;
	std::vector<XYZ> helixaxis = Par2Helix::HelixTerminal(helix);
	XYZ fixedaxis1 = helixaxis[0];
	fixedaxis.z_ = fixedn.z_;
	fixedaxis1.z_ = fixedn.z_;
	//std::cout <<' ' <<fixedaxis.x_ <<' ' <<fixedaxis.y_ <<' ' <<fixedaxis.z_ <<std::endl;
	//std::cout <<' ' <<fixedaxis1.x_ <<' ' <<fixedaxis1.y_ <<' ' <<fixedaxis1.z_ <<std::endl;
	//std::cout <<' ' <<fixedn.x_ <<' ' <<fixedn.y_ <<' ' <<fixedn.z_ <<std::endl;
	double a = angle(fixedaxis, fixedaxis1, fixedn)*180.0/3.14159265358979323846;
	//std::cout <<a <<std::endl;

	XYZ fixedxy = fixedaxis1;
	if(helixaxis[1].z_ > helixaxis[0].z_) fixedxy.z_ += 1.0;
	else fixedxy.z_ -= 1.0;
	LocalFrame lf0 = make_localframe(fixedaxis1, fixedaxis, fixedxy);
	LocalFrame lf1 = make_localframe(XYZ(0.0,0.0,0.0), XYZ(0.0,1.0,0.0), XYZ(0.0,0.0,1.0));

	fixedn = lf0.global2localcrd(fixedn);
	fixedn = lf1.local2globalcrd(fixedn);
	if(fixedn.x_>0.0) a=360.0-a;

	if(a>359.9) return a-360.0;
	return a;
}

/*
 * rotangle, hightA, transdis, lth will be omitted
 * first par's hightA will be given value
 */
void Shape::ReadInputShapePar(std::string parfile) {
	std::vector<std::string> ls;
	std::ifstream ifs(parfile);
	std::string readline;
	while(std::getline(ifs,readline)) ls.push_back(readline);
	ifs.close();
	//hights.push_back(std::stod(ls[0]));
	std::vector<XYZ> axs;
	for(int i=1;i<ls.size();i++) {
		std::stringstream ss(ls[i]);
		XYZ c;
		double d;
		ss >>c.x_ >>c.y_ >>d;
		c.z_ =0.0;
		axs.push_back(c);
		hights.push_back(d);
	}
	for(int i=0;i<axs.size();i++) {
		Par2Helix p;
		p.axisA =axs[i];
		//p.transdis =0.0;
		if(i!=0) p.axisB =pars.back().axisA;
		p.nres =0;
		p.toup =true;
		pars.push_back(p);
	}
	pars[0].axisB =pars[1].axisA;
	pars[0].hightA =std::stod(ls[0]);
	if(hights[0]>pars[0].hightA) {
		for(int i=1;i<pars.size();i+=2) {
			pars[i].toup =false;
		}
	} else {
		for(int i=0;i<pars.size();i+=2) {
			pars[i].toup =false;
		}
	}
}

/*
 * Other 7 pars will be determined previously
 */
void Shape::FindBestLength(Par2Helix &par, double hightC) {
	double min=10000.0;
	std::vector<Residue> ch;
	//std::cout <<par.hightA <<' ' <<hightC <<std::endl;
	for(int i=8;i<50;i++) {
		par.lth =i;
		//std::cout <<i <<std::endl;
		ch =par.HelixGen();
		double delta =fabs(hightC-ch.back().rds().at("CA").crd.z_);
		//std::cout <<i <<' ' <<delta <<' ' <<hightC <<' '
		//		<<ch[0].rds().at("CA").crd.z_ <<' ' <<ch.back().rds().at("CA").crd.z_ <<std::endl;
		if(delta<min) {
			min=delta;
		} else {
			par.lth =i-1;
			break;
		}
	}
	//exit(1);
}






std::vector<XYZ> Shape::extractloopcrd(const std::vector<Residue>&lp) {
	std::vector<XYZ> cs;
	for(int j=0;j<5;j++) cs.push_back(lp[j].rds().at("CA").crd);
	int sz=lp.size();
	for(int j=0;j<5;j++) cs.push_back(lp[sz-5+j].rds().at("CA").crd);
	return cs;
}
void Shape::extractloopcrd() {
	//loopcrd.resize(looplib.size());
	for(int i=0;i<looplib.size();i++) {
		loopcrd.push_back(extractloopcrd(looplib[i]));
	}
}

std::vector<XYZ> Shape::extracthelixcrd(const std::vector<Residue>&ch0, int st,
		const std::vector<Residue>&ch1, int en) {
	int sz=ch0.size();
	std::vector<XYZ> cs;
	for(int i=0;i<5;i++) cs.push_back(ch0[sz-5+i-st].rds().at("CA").crd);
	for(int i=0;i<5;i++) cs.push_back(ch1[i+en].rds().at("CA").crd);
	return cs;
}

std::vector<std::vector<int>> Shape::LoopLibMatch(
		const std::vector<Residue>&ch0, const std::vector<Residue>&ch1) {
	std::vector<std::vector<int>> keys{{0,0},{0,1},{1,0},{0,2},{2,0}};
	for(int i=0;i<keys.size();i++) {
		std::vector<XYZ> crdh=extracthelixcrd(ch0,keys[i][0],ch1,keys[i][1]);
		std::vector<std::vector<int>> rts;
		//std::cout <<i;
		for(int j=0;j<looplib.size();j++) {
			NSPgeometry::QuatFit qf;
			double rmsd=sqrt(qf.setup(crdh,loopcrd[j]));
			//std::cout <<i <<' ' <<j <<' ' <<rmsd <<std::endl;
			if(rmsd>0.5) continue;
			//std::cout <<i <<' ' <<j <<' ' <<rmsd <<std::endl;
			std::vector<Residue> lp=looplib[j];
			RigidTransform rt = qf.getRigidTransform();
			for(Residue &r:lp) {
				for(auto &p:r.rds()) {
					rt.apply(&(p.second.crd));
				}
			}
			int sz=ch0.size();
			int st=-1, en=-1, stl, enl;
			//for(int k=0;k<5;k++) {
			for(int k=1;k<3;k++) {
				double dd = sqrt((ch0[sz-6+k-keys[i][0]].rds().at("C").crd
						-lp[k].rds().at("N").crd).squarednorm());
				//std::cout <<' ' <<dd;
				if(dd>1.99) continue;
				st=sz-6+k-keys[i][0];
				stl=k;
				break;
			}
			//std::cout <<std::endl;
			if(st==-1) continue;
			sz =lp.size();
			//for(int k=0;k<5;k++) {
			for(int k=1;k<3;k++) {
				double dd = sqrt((lp[sz-k-1].rds().at("C").crd
						-ch1[keys[i][1]+5-k].rds().at("N").crd).squarednorm());
				//std::cout <<' ' <<dd;
				if(dd>1.99) continue;
				en=keys[i][1]+5-k;
				enl=sz-k-1;
				break;
			}
			//std::cout <<std::endl;
			if(en==-1) continue;
			rts.push_back({j,keys[i][0],keys[i][1],st,en,stl,enl});
		}
		//std::cout <<std::endl;
		if(rts.empty()) continue;
		return rts;
	}
	//exit(1);
	return {};
}












std::vector<Residue> Shape::Merge(std::vector<std::vector<int>> nlp,
		std::map<int,std::string>&aafixed, std::vector<int>&lths) {
	std::vector<std::vector<Residue>> lps(chs.size()-1);
	std::vector<std::map<int,std::string>> aas(chs.size()-1);
	for(int i=nlp.size()-1;i>=0;i--) {
		std::vector<XYZ> csh=extracthelixcrd(chs[i],nlp[i][1],chs[i+1],nlp[i][2]);
		std::vector<Residue> lp=looplib[nlp[i][0]];
		std::vector<XYZ> csl=extractloopcrd(lp);
		NSPgeometry::QuatFit qf;
		double rmsd2=qf.setup(csh,csl);
		RigidTransform rt = qf.getRigidTransform();
		for(Residue &r:lp) {
			for(auto &p:r.rds()) {
				rt.apply(&(p.second.crd));
			}
		}
		std::map<int,std::string> aa;
		for(int j=5;j<lp.size()-5;j++) {
			std::string aan =lp[j].resname();
			if(aan!="PRO"&&aan!="GLY") continue;
			aa.insert({j-nlp[i][5],aan});
		}
		std::vector<Residue> lp1;
		for(int j=nlp[i][5];j<=nlp[i][6];j++) lp1.push_back(lp[j]);
		lps[i] =lp1;
		aas[i] =aa;
		lp.clear();
		for(int j=nlp[i][4];j<chs[i+1].size();j++) lp.push_back(chs[i+1][j]);
		chs[i+1]=lp;
		lp.clear();
		for(int j=0;j<=nlp[i][3];j++) lp.push_back(chs[i][j]);
		chs[i]=lp;
	}
	std::vector<Residue> ch=chs[0];
	lths.push_back(ch.size());
	for(int i=0;i<nlp.size();i++) {
		for(int j=0;j<lps[i].size();j++) {
			if(aas[i].find(j)!=aas[i].end()) aafixed.insert({ch.size(),aas[i].at(j)});
			ch.push_back(lps[i][j]);
		}
		lths.push_back(ch.size());
		for(int j=0;j<chs[i+1].size();j++) ch.push_back(chs[i+1][j]);
		lths.push_back(ch.size());
	}
	return ch;
}
double Shape::TransDis(double t) {
	t =NSPdstl::RandomEngine<>::getinstance().realrng(0.0,0.2)()-0.1+t;
	if(t<-0.75) t+=1.5;
	if(t>0.75) t-=1.5;
	return t;
}
double Shape::RotAngle(double r) {
	r =NSPdstl::RandomEngine<>::getinstance().realrng(0.0,2.0)()-1.0+r;
	if(r<0.0) r+=360.0;
	if(r>360.0) r-=360.0;
	return r;
}
std::vector<std::pair<int,std::string>> Shape::FindRange(int idx) {
	double rot =Angle2Axis(pars[idx+1].axisA, chs[idx], chs[idx].size()-1);
	//std::cout <<"rotangle: " <<rot <<std::endl;
	std::set<std::pair<int,std::string>> rs;
	for(auto &p0:allowed_par3) {
		double delta=fabs(rot-(double)(p0.first));
		if(delta>1.0) continue;
		for(auto &p1:p0.second) {
			for(const std::string &s:p1.second) {
				rs.insert({p1.first,s});
			}
		}
		//std::cout <<p0.first <<' ' <<rs.size() <<std::endl;
	}
	//std::cout <<rs.size() <<std::endl;
	std::vector<std::pair<int,std::string>> range1;
	for(auto &p:rs) range1.push_back(p);
	std::vector<int> n=RandomOrder(range1.size());
	//std::cout <<n.size() <<std::endl;
	//exit(1);
	std::vector<std::pair<int,std::string>> range;
	for(int i:n) {
		range.push_back(range1[i]);
		//std::cout <<' ' <<i;
	}
	//std::cout <<std::endl;
	//exit(1);
	return range;
}
bool Shape::AddHelix(int idx, std::vector<std::vector<int>> &lpinfo) {
	std::vector<std::pair<int,std::string>> range =FindRange(idx-1);
	if(range.empty()) return false;
	pars[idx].hightA =chs[idx-1].back().rds().at("CA").crd.z_;
	for(int i=0;i<range.size();i++) {
		if(i>=100) return false;
		std::cout <<'/';
		for(int i=0;i<=idx;i++) std::cout <<' ';
		std::cout <<idx <<' ' <<range.size() <<' ' <<i <<std::endl;
		pars[idx].rotangle =RotAngle((double)(range[i].first));
		pars[idx].transdis =TransDis(std::stod(range[i].second));
		//std::cout <<range[i].first <<' ' <<pars[idx].rotangle <<std::endl;
		//std::cout <<range[i].second <<' ' <<pars[idx].transdis <<std::endl;
		FindBestLength(pars[idx],hights[idx]);
		chs[idx] = pars[idx].HelixGen();
		std::vector<std::vector<int>> rts =LoopLibMatch(chs[idx-1],chs[idx]);
		//std::cout <<rts.size() <<std::endl;
		//exit(1);
		if(rts.empty()) continue;
		int n=NSPdstl::RandomEngine<>::getinstance().intrng(0,rts.size()-1)();
		lpinfo[idx-1] =rts[n];
		if(idx==pars.size()-1) return true;
		if(AddHelix(idx+1,lpinfo)) return true;
	}
	return false;
}

bool Shape::Clash(Residue r1, Residue r2, bool a1, bool a2) {
	std::set<std::string> mc{"CA","C","N","O"};
	std::vector<XYZ> c1;
	if(a1) {
		for(auto &p:r1.rds()) c1.push_back(p.second.crd);
	} else {
		for(auto s:mc) c1.push_back(r1.rds().at(s).crd);
	}
	std::vector<XYZ> c2;
	if(a2) {
		for(auto &p:r2.rds()) c2.push_back(p.second.crd);
	} else {
		for(auto s:mc) c2.push_back(r2.rds().at(s).crd);
	}
	double d2=9.0;
	for(XYZ &s1:c1) {
		for(XYZ &s2:c2) {
			double dis2=(s1-s2).squarednorm();
			if(dis2<d2) return true;
		}
	}
	return false;
}
bool Shape::Clash(std::vector<Residue>&ch, const std::map<int,std::string> &aas) {
	bool a1, a2;
	for(int i=0;i<ch.size();i++) {
		if(aas.find(i)==aas.end()) a1=false;
		else a1=true;
		for(int j=i+8;j<ch.size();j++) {
			if(aas.find(j)==aas.end()) a2=false;
			else a2=true;
			if(Clash(ch[i],ch[j],a1,a2)) return true;
		}
	}
	return false;
}

void Shape::ShapeGenRandom(std::string outdir, int nout) {
	chs.resize(pars.size());
	//if(outdir.back()!='/') outdir +='/';
	std::map<std::string,char> a31=AA_Atom::aa31();
	for(int op=0;op<nout;op++) {
		std::cout <<"Result: " <<op <<std::endl;
		std::vector<std::vector<int>> lpinfo(pars.size()-1);
		while(true) {
			std::cout <<"/ 0" <<std::endl;
			pars[0].rotangle =NSPdstl::RandomEngine<>::getinstance().realrng(0.0,360.0)();
			pars[0].transdis =NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.5)()-0.75;
			FindBestLength(pars[0],hights[0]);
			//std::cout <<"lth0: " <<pars[0].lth <<std::endl;
			chs[0] = pars[0].HelixGen();
			//exit(1);
			//std::cout <<22 <<std::endl;
			if(AddHelix(1,lpinfo)) break;
		}
		std::map<int,std::string> aas;
		std::vector<int> lths;
		std::vector<std::vector<Residue>> ch1{Merge(lpinfo,aas,lths)};
		if(Clash(ch1[0],aas)) {
			std::cout <<"Clash! " <<op <<std::endl;
			op--;
			continue;
		}
		std::ofstream ofs(outdir+"par8_"+std::to_string(op));
		for(int i=0;i<pars.size();i++) {
			pars[i].print(ofs);
		}
		ofs.close();
		ofs.open(outdir+"loop_"+std::to_string(op));
		for(auto &v:lpinfo) {
			for(int i:v) ofs <<i <<' ';
			ofs <<std::endl;
		}
		ofs.close();
		residue2pdb(outdir+"pdb_"+std::to_string(op),ch1);
		ofs.open(outdir+"resfile_"+std::to_string(op));
		ofs <<"default allButCys" <<std::endl;
		for(auto &p:aas) {
			ofs <<"A " <<p.first <<' ' <<a31.at(p.second) <<std::endl;
		}
		ofs.close();
		ofs.open(outdir+"length_"+std::to_string(op));
		for(int i:lths) ofs <<i <<std::endl;
		ofs.close();
	}
}


bool Shape::LoopLibMatch(const std::vector<Residue>&ch0,
		const std::vector<Residue>&ch1, double rmsdcut) {
	std::vector<std::vector<int>> keys{{0,0},{0,1},{1,0},{0,2},{2,0}};
	for(int i=0;i<keys.size();i++) {
		std::vector<XYZ> crdh=extracthelixcrd(ch0,keys[i][0],ch1,keys[i][1]);
		std::vector<std::vector<int>> rts;
		for(int j=0;j<looplib.size();j++) {
			NSPgeometry::QuatFit qf;
			double rmsd=sqrt(qf.setup(crdh,loopcrd[j]));
			if(rmsd<rmsdcut) return true;
		}
	}
	return false;
}
void Shape::HelixPairParLib(double dis, double trans, std::string outfile) {
	chs.resize(2);
	pars.resize(2);
	pars[0].axisA.x_ =pars[0].axisA.y_ =pars[0].axisA.z_ =0.0;
	pars[1].axisA =pars[0].axisA;
	pars[1].axisA.y_ =dis;
	pars[0].axisB =pars[1].axisA;
	pars[1].axisB =pars[0].axisA;
	pars[0].lth =pars[1].lth =20;
	pars[0].nres =19;
	pars[1].nres =0;
	pars[0].toup =true;
	pars[1].toup =false;
	pars[0].hightA =0.0;
	pars[0].transdis =0.0;
	pars[1].transdis =trans;
	std::ofstream ofs(outfile);
	for(int i=0;i<360;i++) {
		std::cout <<trans <<' ' <<i <<std::endl;
		pars[0].rotangle =(double)i;
		chs[0] =pars[0].HelixGen();
		pars[1].hightA =chs[0].back().rds().at("CA").crd.z_;
		for(int j=0;j<360;j++) {
			pars[1].rotangle =(double)j;
			chs[1] =pars[1].HelixGen();
			if(LoopLibMatch(chs[0],chs[1],0.6)) ofs <<i <<' ' <<j <<' ' <<trans <<std::endl;
		}
	}
	ofs.close();
}

void Shape::ReadAllowedPar3(std::string parfile) {
	std::string readline;
	std::ifstream ifs(parfile);
	while(std::getline(ifs,readline)) {
		std::stringstream ss(readline);
		int i, j;
		std::string t;
		ss >>i >>j >>t;
		if(allowed_par3.find(i)==allowed_par3.end()) {
			allowed_par3.insert({i,std::map<int,std::set<std::string>>()});
		}
		if(allowed_par3.at(i).find(j)==allowed_par3.at(i).end()) {
			allowed_par3.at(i).insert({j,std::set<std::string>()});
		}
		allowed_par3.at(i).at(j).insert(t);
	}
	ifs.close();
}



