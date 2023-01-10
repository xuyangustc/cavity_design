/*
 * genchain.cpp
 *
 *  Created on: 2018骞�1鏈�25鏃�
 *      Author: hyliu
 */
#include "sd/genchain.h"
#include "proteinrep/pdbrecord.h"
#include "sd/backboneff.h"
#include "proteinrep/intatomkey.h"
#include "sd/sidechainff.h"
#include "fullsite/structplus.h"
#include "dataio/splitstring.h"
using namespace NSPsd;
using namespace NSPproteinrep;
void NSPsd::definegenchaincontrol(std::string name,const std::vector<std::string> &controllines){
	std::map<std::string,double> doublepars;
	std::map<std::string,std::vector<std::string>> stringvecpars;
	std::map<std::string,std::vector<double>> doublevecpars;//{{"AtomSigmas",{3.8,3.8}}};
	std::map<std::string,std::string> stringpars{{"SCTypeFrom","seqstring"},
		{"RefPDBFile",""},{"TypeInCoil","UMC"},
		{"TypeInHelix","ZMC"},{"TypeInStrand","ZMC"},{"SeqString",""},
				{"AllFixedSites",""},{"MCFixedSites",""},{"OutPutPDBFile",""}};
	//std::map<std::string,std::vector<int>> intvecpars{{"AllFixedLocation",{}},{"MCFixedLocation",{}}};
	std::map<std::string,std::vector<int>> intvecpars;
	std::map<std::string,int>intpars;
	GenChainControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,
			stringvecpars,intvecpars);
	int nsuccess=GenChainControls::adjustvalues(name,controllines);
	if(nsuccess!= controllines.size()) {
		exit(1);
	}
	auto & pset=GenChainControls::getparameterset(name);
}
GenChain::GenChain(const std::string &controlname):
	controlname_(controlname){
	auto & pset=GenChainControls::getparameterset(controlname);
	std::string sctypefrom,refpdbfile,seqstring;
	pset.getval("SCTypeFrom",&sctypefrom);
	pset.getval("RefPDBFile",&refpdbfile);
	pset.getval("SeqString", &seqstring);
	if(!refpdbfile.empty()){
		refpdb_=std::shared_ptr<std::vector<std::vector<FullSite>>>(new
				std::vector<std::vector<FullSite>>);
		*refpdb_=readfullsitesfrompdb(refpdbfile);
		//sprefpdb_=std::shared_ptr<StructPlus>(new
		//		StructPlus(*refpdb_,true));
	} else if(sctypefrom.find("refpdb") != std::string::npos){
		std::cout <<"No RefPDBFile defined. Cannot obtain sidechain types from refpdb."<<std::endl;
		exit(1);
	}
	if(sctypefrom=="refpdb"){
		for(auto &c:*refpdb_){
			sctypes_.push_back(std::vector<std::string>());
			softsc_.push_back(std::vector<bool>());
			for(auto &aa:c){
				std::string upper;
				for(char a:aa.resname()) upper.push_back(toupper(a));
				if(upper==aa.resname()) softsc_.back().push_back(false);
				else softsc_.back().push_back(true);
				try {
					VSCType::getVSCType(upper);
				} catch (std::exception &e){
					std::cout <<"Undefined residue type in PDB:" <<aa.resname()<<std::endl;
					exit(1);
				}
				sctypes_.back().push_back(upper);
			}
		}
		if(!seqstring.empty()) {
			std::stringstream ss(seqstring);
			std::vector<std::string> vs;
			while(ss) {
				std::string s;
				ss >> s;
				vs.push_back(s);
			}
			std::vector<std::pair<std::vector<int>,std::string>> seq2change;
			for(int i=0;;i+=4) {
				if(i+3>=vs.size()) break;
				int c=std::stoi(vs[i]);
				int st=std::stoi(vs[i+1]);
				int en=std::stoi(vs[i+2]);
				if(en-st+1!=vs[i+3].size()) {
					std::cout <<"Wrong in SeqString!" <<std::endl;
					exit(1);
				}
				seq2change.push_back({{c,st},vs[i+3]});
			}
			for(int i=0;i<seq2change.size();i++) {
				int c=seq2change[i].first[0];
				int st=seq2change[i].first[1];
				for(int j=0;j<seq2change[i].second.size();j++) {
					std::string aa=VSCType::resnameof(seq2change[i].second[j]);
					sctypes_[c][st+j] = aa;
					(*refpdb_)[c][st+j].changeaa(aa);
				}
			}
		}
		if(!refpdbfile.empty()) sprefpdb_=std::shared_ptr<StructPlus>(new StructPlus(*refpdb_,true));
	}
	if(sctypefrom=="sstypeinrefpdb"){
		if(!refpdbfile.empty()) sprefpdb_=std::shared_ptr<StructPlus>(new StructPlus(*refpdb_,true));
		std::string htype,stype,ctype;
		pset.getval("TypeInHelix",&htype);
		pset.getval("TypeInStrand",&stype);
		pset.getval("TypeInCoil",&ctype);
		StructPlus & sp=*sprefpdb_;
		for(int cid=0;cid<refpdb_->size();++cid){
			sctypes_.push_back(std::vector<std::string>());
			softsc_.push_back(std::vector<bool>());
			for(int s=0;s<(*refpdb_)[cid].size();++s){
				char ss=sp.sstype(StructPlus::Position(cid,s));
				std::string sc;
				if(ss=='H') sc=htype;
				else if(ss=='E') sc=stype;
				else if(ss=='C') sc=ctype;
				std::string upper;
				for(char a:sc) upper.push_back(toupper(a));
				if(upper==sc) softsc_.back().push_back(false);
				else softsc_.back().push_back(true);
				try{
					VSCType::getVSCType(upper);
				} catch (std::exception &e){
					std::cout<<"Undefined sidechain type " <<sc<<std::endl;
					exit(1);
				}
				sctypes_.back().push_back(upper);
			}
		}
	}
	if(sctypefrom=="seqstring"){
		if(seqstring.empty()){
			std::cout<<"No SeqString defined."<<std::endl;
			exit(1);
		}
		std::vector<std::string> subseqs=NSPdataio::wordsInString(seqstring,",;:");
		sctypes_.assign(subseqs.size(),std::vector<std::string>());
		for(int c=0;c<subseqs.size();++c){
			softsc_.push_back(std::vector<bool>());
			for(int s=0;s<subseqs[c].size();++s){
				char aa=subseqs[c][s];
				char upper=toupper(aa);
				if(upper==aa) softsc_.back().push_back(false);
				else softsc_.back().push_back(true);
				sctypes_[c].push_back(VSCType::resnameof(upper));
			}
		}
		if(!refpdbfile.empty()) sprefpdb_=std::shared_ptr<StructPlus>(new StructPlus(*refpdb_,true));
	}
	//For the moment, the terminal residues cannot have any side chain torsional angles
	//This should be changed in later versions
/*	for(auto &c:sctypes_){
		if(!(VSCType::getVSCType(c[0]).rotameratoms.empty())){
			std::cout <<"Replacing first residue with ALA"<<std::endl;
			c[0]="ALA";
		}
		if(!(VSCType::getVSCType(c.back()).rotameratoms.empty())){
			std::cout <<"Replacing last residue with ALA"<<std::endl;
			c.back()="ALA";
		}
	}*/
}
std::string GenChain::seqstring() const {
	std::string seq;
	for(int c=0;c<sctypes_.size();++c){
		if(c>0) seq.push_back(';');
		for(auto & sc:sctypes_[c]){
			auto &vsc=VSCType::getVSCType(sc);
			seq.push_back(vsc.oneletter);
		}
	}
	return seq;
}
std::vector<double> GenChain::getcrd(const std::vector<std::vector<BackBoneSite>> &bssites,
		bool userefpdb) const {
	std::vector<NSPgeometry::XYZ> res;
	double deg = 3.14159265 / 180.0;
	assert(bssites.size()==sctypes_.size());
	for (int c=0;c<bssites.size();++c){
		assert(bssites[c].size()==sctypes_[c].size());
		for(int s=0;s<bssites[c].size();++s){
			const VSCType &vsc=VSCType::getVSCType(sctypes_[c][s]);
			const BackBoneSite &bs=bssites[c][s];
			std::vector<NSPgeometry::XYZ> r(vsc.nscatoms+4);
			r[0]=bs.ncrd();
			r[1]=bs.cacrd();
			r[vsc.nscatoms+2]=bs.ccrd();
			r[vsc.nscatoms+3]=bs.ocrd();
			std::vector<std::vector<std::pair<int,double>>> ics;
			if(userefpdb && refpdb_){
				if(sctypes_[c][s]==(*refpdb_)[c][s].resname()){
					ics=(*refpdb_)[c][s].internalcrds();
				}
			}
			if(ics.empty()){
				ics=vsc.internalcrds;
			}
			assert(ics.size()==vsc.nscatoms);
			int aidx=2;
			for(auto &ic:ics){
				r[aidx++]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
									r[ic[2].first],
									ic[0].second, ic[1].second*deg,ic[2].second*deg);
			}
			for(auto &x:r) res.push_back(x);
		}
	}
	std::vector<double> crd;
	for (auto &r : res) {
		for (int m = 0; m < 3; ++m)
			crd.push_back(r[m]);
	}
	return crd;
}
std::vector<double> GenChain::getcrd(const std::string &pdbfile) const {
	std::vector<std::vector<FullSite>> sites=readfullsitesfrompdb(pdbfile);
	return getcrd(sites);
}
std::vector<double> GenChain::getcrd(const
	std::vector<std::vector<FullSite>> &sites) const{
	std::vector<NSPgeometry::XYZ> res;
	double deg=3.14159265/180.0;
	for (int c=0;c<sites.size();++c){
		assert(sites[c].size()==sctypes_[c].size());
		for(int s=0;s<sites[c].size();++s){
			const VSCType &vsc=VSCType::getVSCType(sctypes_[c][s]);
			std::vector<NSPgeometry::XYZ> r(vsc.nscatoms+4);
			const FullSite &fs=sites[c][s];
			r[0]=fs.getcrd("N");
			r[1]=fs.getcrd("CA");
			r[vsc.nscatoms+2]=fs.getcrd("C");
			r[vsc.nscatoms+3]=fs.getcrd("O");
			bool sccomplete{true};
			for(int a=0;a<vsc.nscatoms;++a){
				if(!fs.hasatomcrd(vsc.atomnames[a])) {
					sccomplete=false;
					break;
				}
			}
			if( sccomplete){
				for(int a=0;a<vsc.nscatoms;++a){
					r[a+2]=fs.getcrd(vsc.atomnames[a]);
				}
			} else {
				auto & ics=vsc.internalcrds;
				assert(ics.size()==vsc.nscatoms);
				int aidx=2;
				for(auto &ic:ics){
					r[aidx++]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
								r[ic[2].first],
								ic[0].second, ic[1].second*deg,ic[2].second*deg);
						}
			}
			for(auto &x:r) res.push_back(x);
		}
	}
	std::vector<double> crd;
	for (auto &r : res) {
		for (int m = 0; m < 3; ++m)
			crd.push_back(r[m]);
	}
	return crd;
}
void GenChain::writepdb(const std::vector<double> & crd, std::ostream &os,double crdtoangstrom) const{
	int offset = 0;
	std::vector<PdbRecord> records;
	for(int c=0;c<sctypes_.size();++c){
		for (int i = 0; i < sctypes_[c].size(); ++i) {
			const VSCType &vsc=VSCType::getVSCType(sctypes_[c][i]);
			std::string resname=vsc.pdbname;
			auto nkey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "N", c, resname,0);
			NSPgeometry::XYZ ncrd(crd[offset], crd[offset + 1], crd[offset + 2]);
			offset += 3;
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(nkey, ncrd*crdtoangstrom, offset / 3 + 1));
			auto cakey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1,"CA", c, resname,0);
			NSPgeometry::XYZ cacrd(crd[offset], crd[offset + 1], crd[offset + 2]);
			offset += 3;
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(cakey,
					cacrd*crdtoangstrom, offset / 3 + 1));
			if (vsc.nscatoms>0) {
				for(int m=0;m<vsc.atomnames.size();++m){
					auto akey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, vsc.atomnames[m], c,
					resname,0);
					NSPgeometry::XYZ r(crd[offset], crd[offset + 1],
							crd[offset + 2]);
					offset += 3;
					records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(akey,
							r*crdtoangstrom, offset / 3 + 1));
				}
			}
			auto ckey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "C", c, resname,0);
			NSPgeometry::XYZ ccrd(crd[offset], crd[offset + 1], crd[offset + 2]);
			offset += 3;
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(ckey,
					ccrd*crdtoangstrom, offset / 3 + 1));
			auto okey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "O", c, resname,0);
			NSPgeometry::XYZ ocrd(crd[offset], crd[offset + 1], crd[offset + 2]);
			offset += 3;
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(okey,
					ocrd*crdtoangstrom, offset / 3 + 1));
		}
	}
	for (auto &r : records) {
		os << r.toString()<<std::endl;
	}
}
static void make_ff_chain_nophipsi(const std::vector<int> & nsites,
		ForceField &ff, const std::vector<std::vector<std::string>> &sctypes,
		const std::vector<std::vector<int>> &cissites,
		std::map<std::vector<int>,std::map<std::string,int>>&atomidinseq) {
	int offset = 0;
	int chainid = 0;
	int nscatoms = 0;
	std::vector<int> & stericatomtypes=ff.stericatomtypes();
	stericatomtypes.clear();
	for (int ns : nsites) {
		for (int i = 0; i < ns; ++i) {
			int nscatoms_p=nscatoms;
			const VSCType & vsc=VSCType::getVSCType(sctypes[chainid][i]);
			nscatoms = vsc.nscatoms;
			if (i > 0) {
				ff.addbond(offset - 2, offset, backboneff.b0_cn,
						2*KBT*backboneff.kb_cn);
				ff.addangle(offset - 3 - nscatoms_p, offset - 2, offset,
						backboneff.t0_cacn,2*KBT*KANG_FAC*backboneff.kt_cacn);
				ff.addangle(offset - 1, offset - 2, offset, backboneff.t0_ocn,
						2*KBT*KANG_FAC*backboneff.kt_ocn);
				ff.addangle(offset - 2, offset, offset + 1, backboneff.t0_cnca,
						2*KBT*KANG_FAC*backboneff.kt_cnca);
				ff.addimpdih(offset - 2, offset, offset - 3-nscatoms_p, offset - 1,
						backboneff.p_cncao, 2*KBT*KANG_FAC*backboneff.kp_cncao);
				double p0 = backboneff.p_cacnca;
				double p1 = backboneff.p_ocnca;
				for (auto s : cissites[chainid]) {
					if (i == s) {
						p0 = 0.0;
						p1 = 180.0;
						break;
					}
				}
				ff.addimpdih(offset - 3 - nscatoms_p, offset - 2, offset,
						offset + 1, p0, 2*KBT*KANG_FAC*backboneff.kp_cacnca);
				ff.addimpdih(offset - 1, offset - 2, offset, offset + 1, p1,
						2*KBT*KANG_FAC*backboneff.kp_ocnca);
				ff.adddih(offset - 4 - nscatoms_p, offset - 3 - nscatoms_p,
						offset - 2, offset);
				ff.adddih(offset - 2, offset, offset + 1,
						offset + 2 + nscatoms);
			}
			ff.addbsinchain(chainid, offset, offset + 1, offset + 2 + nscatoms,
					offset + 3 + nscatoms);
			std::vector<std::vector<int>> kaiatoms;
			for(int rsc:vsc.rotameratoms){
				kaiatoms.push_back(std::vector<int>());
				kaiatoms.back().push_back(rsc+offset+2);
				for(int i=0;i<3;++i){
					kaiatoms.back().push_back(vsc.internalcrds[rsc][i].first+offset);
				}
			}
/*			if(!kaiatoms.empty() &&(i==0 || i==ns-1)){
				std::cout <<"Current code does not support residues with side chain torsions at chain terminus"<<std::endl;
				exit(1);
			}*/
			if(ff.scinchains().size()<chainid+1) ff.scinchains().resize(chainid+1);
			ff.scinchains()[chainid].push_back(SCInChain(vsc.resname,kaiatoms,offset+2,nscatoms));
			ff.addbond(offset, offset + 1, backboneff.b0_nca,
					2*KBT*backboneff.kb_nca);
			ff.addbond(offset + 1, offset + 2 + nscatoms, backboneff.b0_cac,
					2*KBT*backboneff.kb_cac);
			ff.addbond(offset + 2+nscatoms, offset + 3 + nscatoms, backboneff.b0_co,
					2*KBT*backboneff.kb_co);
			ff.addangle(offset, offset + 1, offset + 2 + nscatoms,
					backboneff.t0_ncac, 2*KBT*KANG_FAC*backboneff.kt_ncac);
			ff.addangle(offset + 1, offset + 2 + nscatoms,
					offset + 3 + nscatoms, backboneff.t0_caco,
					2*KBT*KANG_FAC*backboneff.kt_caco);
			ff.adddih(offset, offset + 1, offset + 2 + nscatoms,
					offset + 3 + nscatoms);
			ff.setstericatom(offset, backboneff.sigma_n, backboneff.eps, true,
					false, backboneff.sigma_nhb);
			ff.setstericatom(offset + 1, backboneff.sigma_ca, backboneff.eps);
			stericatomtypes.push_back(VSCType::getstericatomtype(vsc.resname,"N"));
			stericatomtypes.push_back(VSCType::getstericatomtype(vsc.resname,"CA"));
			atomidinseq.insert({{chainid,i},{{"N",offset}}});
			atomidinseq.at({chainid,i}).insert({"CA",offset+1});
			if (nscatoms>0) {
				int aidx=2;
				for(auto & a:vsc.atomnames){
					int atype=VSCType::getstericatomtype(vsc.resname,a);
					double sigma=VSCType::packingatomtypes[atype].radius;
					int hbtype =VSCType::packingatomtypes[atype].hbtype;
					bool hbdonor=(hbtype==1 ||hbtype==3);
					bool hbacceptor=(hbtype==2||hbtype==3);
					double sigmahb=2.8;
					ff.setstericatom(offset + aidx++, sigma, 0.2,hbdonor,hbacceptor,sigmahb);
					atomidinseq.at({chainid,i}).insert({a,offset+aidx-1});
					stericatomtypes.push_back(atype);
				}
				for(int b=0;b<vsc.newbonds.size();++b){
					ff.addbond(offset+vsc.newbonds[b].first,offset+vsc.newbonds[b].second,
							vsc.b0[b],vsc.kb0[b]);
				}
				for(int a=0;a<vsc.newangles.size();++a){
						ff.addangle(offset+vsc.newangles[a][0],
								offset+vsc.newangles[a][1], offset+vsc.newangles[a][2],
								vsc.a0[a],vsc.ka0[a]);
					}

				for(int a=0;a<vsc.newimpdihs.size();++a){
						ff.addimpdih(offset+vsc.newimpdihs[a][0],
								offset+vsc.newimpdihs[a][1], offset+vsc.newimpdihs[a][2],
								offset+vsc.newimpdihs[a][3],
								vsc.imp0[a],vsc.kimp0[a]);
					}
				for(int a=0;a<vsc.newtorsions.size();++a){
					bool keep=true;
					for(auto k:vsc.newtorsions[a]){
						if(k<0 && i==0) keep=false;
						if(k>=nscatoms+4 && i==ns-1) keep=false;
					}
					if(!keep) continue;
					ff.adddih(offset+vsc.newtorsions[a][0],
								offset+vsc.newtorsions[a][1], offset+vsc.newtorsions[a][2],
								offset+vsc.newtorsions[a][3]);
				}
			}
			ff.setstericatom(offset + 2 + nscatoms, backboneff.sigma_c,
					backboneff.eps);
			ff.setstericatom(offset + 3 + nscatoms, backboneff.sigma_o,
					backboneff.eps, false, true, backboneff.sigma_ohb);
			stericatomtypes.push_back(VSCType::getstericatomtype(vsc.resname,"C"));
			stericatomtypes.push_back(VSCType::getstericatomtype(vsc.resname,"O"));
			atomidinseq.at({chainid,i}).insert({"C",offset + 2 + nscatoms});
			atomidinseq.at({chainid,i}).insert({"O",offset + 3 + nscatoms});
			offset += backboneff.natomspersites + nscatoms;
		}
		++chainid;
	}
}
std::vector<NSPallatom::ChemGroup> readcgs(std::string refpdbfile,
		std::map<std::vector<int>,std::map<std::string,int>> atomidinseq) {
	NSPallatom::PdbReader_xy pr;
	pr.init(refpdbfile,true);
	auto ress = pr.chains_new();
	NSPallatom::Atom3 a3(ress,12.0,6);
	a3.findatom3();
	a3.cover2NM();
	auto cs=a3.cgs();
	std::vector<NSPallatom::ChemGroup> cgs;
	for(auto &c:cs) {
		NSPallatom::ChemGroup g=*c;
		g.atomseqs.resize(g.resseqs.size(),0);
		std::vector<int> k{g.chainseq,g.resseq};
		for(int i=0;i<g.resseqs.size();i++) {
			g.atomseqs[i] = atomidinseq.at({g.chainseq,g.resseqs[i]}).at(g.crds_v[i].first);
		}
		cgs.push_back(g);
	}
	return cgs;
}
ForceField GenChain::make_forcefield(const std::string &controlname, std::vector<std::set<int>> notcis) {
	ForceField ff;
	std::vector<std::vector<int>> cissites;
	if(refpdb_){
		for(auto &chain:*refpdb_){
			std::vector<BackBoneSite> bs=backbone(chain);
			cissites.push_back(findcissites(bs));
		}
	}
	if (cissites.empty())
		cissites.assign(sctypes_.size(), std::vector<int>());
	for(int i=0;i<sctypes_.size();i++) {
		std::set<int> iscis;
		for(int ix:cissites[i]) iscis.insert(ix);
		for(int ix:cissites[i]) {
			if(notcis[i].find(ix)==notcis[i].end()) continue;
			iscis.erase(ix);
		}
		cissites[i].clear();
		for(int ix:iscis) cissites[i].push_back(ix);
	}
	int natoms = 0;
	std::vector<int>nsites;
	for (int i = 0; i < sctypes_.size(); ++i) {
		nsites.push_back(sctypes_[i].size());
		for (int j = 0; j < sctypes_[i].size(); ++j) {
			natoms += 4+VSCType::getVSCType(sctypes_[i][j]).nscatoms;
		}
	}
	ff.init(natoms);
	//std::map<std::vector<int>,std::map<std::string,int>> aiis;
	//make_ff_chain_nophipsi(nsites, ff, sctypes_, cissites, aiis);
	//atomidinseq=aiis;
	make_ff_chain_nophipsi(nsites, ff, sctypes_, cissites, atomidinseq);
	//for(auto &ad:atomidinseq) {
	//	for(int i:ad.first) std::cout <<i <<' ';
	//	std::cout <<std::endl;
	//	for(auto p:ad.second) std::cout <<' ' <<p.first <<' ' <<p.second <<std::endl;
	//}
	//exit(1);
	std::string conname=controlname.substr(0,controlname.size()-3) + "_genchain";
	auto &pset1 = SDControls::getparameterset(conname);
	std::string refpdbfile;
	pset1.getval("RefPDBFile",&refpdbfile);
	std::vector<NSPallatom::ChemGroup> cgs=readcgs(refpdbfile, atomidinseq);
	for(auto &cg:cgs) ff.addchemgroup(cg);

	auto & scs=ff.scinchains();
	for(int c=0;c<scs.size();++c){
		for(int p=0;p< scs[c].size();++p){
			scs[c][p].softpacking=softsc_[c][p];
		}
	}
	std::string seq;
	int nseq=0;
	std::map<std::string,char> aa20 {
		{"GLY",'G'},
		{"ALA",'A'},
		{"VAL",'V'},
		{"LEU",'L'},
		{"ILE",'I'},
		{"SER",'S'},
		{"THR",'T'},
		{"CYS",'C'},
		{"MET",'M'},
		{"ASP",'D'},
		{"GLU",'E'},
		{"ASN",'N'},
		{"GLN",'Q'},
		{"LYS",'K'},
		{"ARG",'R'},
		{"HIS",'H'},
		{"PRO",'P'},
		{"PHE",'F'},
		{"TYR",'Y'},
		{"TRP",'W'}
	};
	for(auto sc:sctypes_) for(auto s:sc) seq.push_back(aa20.at(s));
	int chainid = 0;
	for (int ns : nsites) {
		for (int i = 0; i < ns; ++i) {
			if (i == 0) {
				ff.addphipsi(i, chainid, true, false);
			} else if (i == ns - 1) {
				ff.addphipsi(i, chainid, false, true);
			} else {
				ff.addphipsi(i, chainid, false, false);
			}
			ff.changephipsiaa(nseq++, seq[nseq-1]);
		}
		++chainid;
	}
	ff.setexcln14();
	ff.usecontrols(controlname,this);
	return ff;
}

StochasticDynamics GenChain::make_sd(const std::string &controlname,
		const std::vector<std::vector<int>> &atomgroups, const std::vector<double> & temperatures) const{
	std::vector<double> masses;
	std::vector<double> gammas;
	double gamma;
	NSPdataio::ParameterSet &pset=SDControls::getparameterset(controlname);
	pset.getval("FrictionCoeff",&gamma);
	for(int c=0;c<sctypes_.size();++c){
		for (int i = 0; i < sctypes_[c].size(); ++i) {
			int natoms=4+VSCType::getVSCType(sctypes_[c][i]).nscatoms;
			for(int a=0;a<natoms;++a)
				masses.push_back(14.0);
			for (int d = 0; d < natoms; ++d)
				gammas.push_back(gamma);
		}
	}
	StochasticDynamics sd;
	double timestep;
	pset.getval("TimeStep",&timestep);
	sd.init(masses, gammas, timestep, atomgroups, temperatures,3);
	return sd;
}
ActiveSelections * GenChain::setactiveselections(const ForceField *ff){
	std::string allfixed;
	std::string mcfixed;
	auto & pset=GenChainControls::getparameterset(controlname_);
	pset.getval("AllFixedSites",&allfixed); //std::cout <<allfixed <<std::endl;
	pset.getval("MCFixedSites",&mcfixed); //std::cout <<mcfixed <<std::endl;
	acts_=std::shared_ptr<ActiveSelections>(new ActiveSelections(ff,allfixed,mcfixed));
	return acts_.get();
}
ActiveSelections * GenChain::setactiveselections(const ForceField *ff, std::string allfixed, std::string mcfixed){
	acts_=std::shared_ptr<ActiveSelections>(new ActiveSelections(ff,allfixed,mcfixed));
	return acts_.get();
}
SDRun GenChain::make_sdrun(const SDRun::SDRunIn &in, unsigned int seed, std::vector<std::set<int>> notcis){
	auto ff = std::shared_ptr < ForceField > (new ForceField);
	auto &pset = SDControls::getparameterset(in.sdcontrolname + "_sd");
	*ff=make_forcefield(in.sdcontrolname + "_ff", notcis);

	/*std::vector<NSPallatom::ChemGroup> cgs=ff->cgps();
	std::cout <<cgs.size() <<std::endl;
	for(int i=0;i<cgs.size();i++) {
		std::cout <<i <<std::endl;
		std::cout <<cgs[i].chainseq <<'\t' <<cgs[i].resseq <<'\t' <<cgs[i].resname <<std::endl;
		for(auto &p:cgs[i].crds_v) std::cout <<p.first <<' ';
		std::cout <<std::endl <<std::endl;
	}
	exit(1);*/

	auto shakebds = std::shared_ptr < ShakeBonds > (new ShakeBonds);
	int doshake;
	pset.getval("DoShake", &doshake);
	*shakebds = make_shakebonds(*ff);
	shakebds->seton(doshake != 0);
	std::vector<std::string> tgroups;
	pset.getval("TemperatureGroups",&tgroups);
	auto agrps=std::shared_ptr<std::vector<std::vector<int>>>(
			new std::vector<std::vector<int>>(tgroups.size()));
	int na=0;
	for(int i=0;i<tgroups.size();++i){
		if(tgroups[i]=="all") {
			for(int a=0;a<ff->natoms();++a) (*agrps)[i].push_back(a);
		} else	if(tgroups[i]=="mainchain"){
			(*agrps)[i]=ff->mainchainatoms();
		} else if(tgroups[i] =="sidechain"){
			(*agrps)[i]=ff->sidechainatoms();
		}  else if(tgroups[i] =="ssregions"){
			assert(sprefpdb_ !=nullptr);
			for(int c=0;c<sctypes_.size();++c){
				for(int p=0;p<sctypes_[c].size();++p){
					if(sprefpdb_->sstype(StructPlus::Position(c,p)) =='C') continue;
					std::vector<int> ia=ff->atomsinresidue(c,p);
					for(int a:ia) (*agrps)[i].push_back(a);
				}
			}
		} else if(tgroups[i]=="coilregions"){
			for(int c=0;c<sctypes_.size();++c){
				for(int p=0;p<sctypes_[c].size();++p){
					if(sprefpdb_->sstype(StructPlus::Position(c,p)) !='C') continue;
					std::vector<int> ia=ff->atomsinresidue(c,p);
					for(int a:ia) (*agrps)[i].push_back(a);
				}
			}
		} else {
			std::cout<<"Undefined TemperatureGroup Type: " <<tgroups[i]<<std::endl;
			exit(1);
		}
		na +=(*agrps)[i].size();
	}
	std::vector<double> temperatures;
	pset.getval("Temperatures",&temperatures);
	for(auto &t:temperatures) t*=KBT;
	assert (na==ff->natoms());
	auto sd = std::shared_ptr < StochasticDynamics > (new StochasticDynamics);
	*sd = make_sd(in.sdcontrolname + "_sd", *agrps,temperatures);
	SDRun sdrun(sd, ff, shakebds);
	sdrun.temperaturegroups()=agrps;
	sdrun.bathtemperatures()=temperatures;
	sdrun.shakeon() = doshake != 0;
	int nblsteps;
	pset.getval("NeighborListSteps", &nblsteps);
	sdrun.nblsteps() = nblsteps;
	sdrun.initrandomengine(seed);
	std::vector<bool> fixatoms(ff->natoms(),false);
	int  fixmainchain;
	pset.getval("FixMainChain",&fixmainchain);
	if(fixmainchain!=0){
		std::vector<int> mcatoms=ff->mainchainatoms();
		for(int a:mcatoms) fixatoms[a]=true;
	}
	SDRun::SDRunIn newin(*(in.crd),fixatoms);
	if(!acts_) this->setactiveselections(ff.get());
	sdrun.setactiveselections(acts_.get());
	sdrun.initstate(newin);
	return sdrun;
}
SDRun GenChain::make_sdrun(const SDRun::SDRunIn &in, unsigned int seed,
		std::string allfixed, std::string mcfixed, std::vector<std::set<int>> notcis){
	auto ff = std::shared_ptr < ForceField > (new ForceField);
	auto &pset = SDControls::getparameterset(in.sdcontrolname + "_sd");
	*ff=make_forcefield(in.sdcontrolname + "_ff", notcis);
	auto shakebds = std::shared_ptr < ShakeBonds > (new ShakeBonds);
	int doshake;
	pset.getval("DoShake", &doshake);
	*shakebds = make_shakebonds(*ff);
	shakebds->seton(doshake != 0);
	std::vector<std::string> tgroups;
	pset.getval("TemperatureGroups",&tgroups);
	auto agrps=std::shared_ptr<std::vector<std::vector<int>>>(
			new std::vector<std::vector<int>>(tgroups.size()));
	int na=0;
	for(int i=0;i<tgroups.size();++i){
		if(tgroups[i]=="all") {
			for(int a=0;a<ff->natoms();++a) (*agrps)[i].push_back(a);
		} else	if(tgroups[i]=="mainchain"){
			(*agrps)[i]=ff->mainchainatoms();
		} else if(tgroups[i] =="sidechain"){
			(*agrps)[i]=ff->sidechainatoms();
		}  else if(tgroups[i] =="ssregions"){
			assert(sprefpdb_ !=nullptr);
			for(int c=0;c<sctypes_.size();++c){
				for(int p=0;p<sctypes_[c].size();++p){
					if(sprefpdb_->sstype(StructPlus::Position(c,p)) =='C') continue;
					std::vector<int> ia=ff->atomsinresidue(c,p);
					for(int a:ia) (*agrps)[i].push_back(a);
				}
			}
		} else if(tgroups[i]=="coilregions"){
			for(int c=0;c<sctypes_.size();++c){
				for(int p=0;p<sctypes_[c].size();++p){
					if(sprefpdb_->sstype(StructPlus::Position(c,p)) !='C') continue;
					std::vector<int> ia=ff->atomsinresidue(c,p);
					for(int a:ia) (*agrps)[i].push_back(a);
				}
			}
		} else {
			std::cout<<"Undefined TemperatureGroup Type: " <<tgroups[i]<<std::endl;
			exit(1);
		}
		na +=(*agrps)[i].size();
	}
	std::vector<double> temperatures;
	pset.getval("Temperatures",&temperatures);
	for(auto &t:temperatures) t*=KBT;
	assert (na==ff->natoms());
	auto sd = std::shared_ptr < StochasticDynamics > (new StochasticDynamics);
	*sd = make_sd(in.sdcontrolname + "_sd", *agrps,temperatures);
	SDRun sdrun(sd, ff, shakebds);
	sdrun.temperaturegroups()=agrps;
	sdrun.bathtemperatures()=temperatures;
	sdrun.shakeon() = doshake != 0;
	int nblsteps;
	pset.getval("NeighborListSteps", &nblsteps);
	sdrun.nblsteps() = nblsteps;
	sdrun.initrandomengine(seed);
	std::vector<bool> fixatoms(ff->natoms(),false);
	int  fixmainchain;
	pset.getval("FixMainChain",&fixmainchain);
	if(fixmainchain!=0){
		std::vector<int> mcatoms=ff->mainchainatoms();
		for(int a:mcatoms) fixatoms[a]=true;
	}
	SDRun::SDRunIn newin(*(in.crd),fixatoms);
	if(!acts_) this->setactiveselections(ff.get(),allfixed,mcfixed);
	sdrun.setactiveselections(acts_.get());
	sdrun.initstate(newin);
	return sdrun;
}
void NSPsd::genchainreadcontrols(const std::string &filename,std::string name){
	NSPdataio::ControlFile cf;
	cf.readfile(filename);
	std::vector<std::string> sdcontrolines=cf.getcontrolines("SD");
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	std::vector<std::string> genchaincontrollines=cf.getcontrolines("GenChain");
	definesdcontrol(name+"_sd",sdcontrolines);
	defineforcefieldcontrol(name+"_ff",ffcontrolines);
	definegenchaincontrol(name+"_genchain",genchaincontrollines);
}
void NSPsd::genchainprintcontrols(std::string name,std::ostream &ofs){
	ofs<<"START GenChain"<<std::endl;
	GenChainControls::getparameterset(name+"_genchain").printparameters(ofs);
	ofs <<"START SD"<<std::endl;
	SDControls::getparameterset(name+"_sd").printparameters(ofs);
	ofs<<"END SD"<<std::endl;
	ofs<<"START ForceField"<<std::endl;
	ForceFieldControls::getparameterset(name+"_ff").printparameters(ofs);
	ofs<<"END ForceField"<<std::endl;
}
















std::vector<std::vector<FullSite>> GenChain::crd2fullsite(const std::vector<double> &crds, double nm2a) const {
	std::vector<std::vector<FullSite>> sites = *refpdb_;
	std::vector<NSPgeometry::XYZ> res;
	double deg=3.14159265/180.0;
	int count=0;
	for (int c=0;c<sites.size();++c){
		assert(sites[c].size()==sctypes_[c].size());
		for(int s=0;s<sites[c].size();++s){
			const VSCType &vsc=VSCType::getVSCType(sctypes_[c][s]);
			std::vector<std::string> atomnames{"N","CA"};
			std::vector<std::string> atomnames1=vsc.atomnames;
			for(auto &a:atomnames1) atomnames.push_back(a);
			atomnames.push_back("C");
			atomnames.push_back("O");
			std::map<std::string,NSPgeometry::XYZ> atomcrds;
			for(int i=0;i<atomnames.size();i++,count++) {
				NSPgeometry::XYZ c(crds[count*3]*nm2a,crds[count*3+1]*nm2a,crds[count*3+2]*nm2a);
				atomcrds.insert({atomnames[i],c});
			}
			sites[c][s].changecrd(atomcrds);
		}
	}
	return sites;
}

std::vector<std::vector<NSPallatom::Residue>> GenChain::crd2res(const std::vector<double> &crds,
		const std::vector<std::vector<NSPallatom::Residue>> &refres, double nm2a) {
	std::vector<std::vector<NSPallatom::Residue>> resss=refres;
	int count=0;
	for (int c=0;c<refres.size();++c){
		for(int r=0;r<refres[c].size();++r){
			const VSCType &vsc=VSCType::getVSCType(refres[c][r].resname());
			std::vector<std::string> atomnames{"N","CA"};
			std::vector<std::string> atomnames1=vsc.atomnames;
			for(auto &a:atomnames1) atomnames.push_back(a);
			atomnames.push_back("C");
			atomnames.push_back("O");
			for(int i=0;i<atomnames.size();i++,count++) {
				NSPgeometry::XYZ crd(crds[count*3]*nm2a,crds[count*3+1]*nm2a,crds[count*3+2]*nm2a);
				resss[c][r].rds().at(atomnames[i]).crd = crd;
			}
		}
	}
	return resss;
}
/*
double GenChain::calrmsd(std::vector<double> &crd1, std::vector<double> &crd2) {
	int natom=crd1.size()/3;
	int actatom=0;
	double r2=0;
	for(int i=0;i<natom;i++) {
		if(!(acts_->atomactive(i))) continue;
		actatom++;
		double dx=crd1[i*3]-crd2[i*3];
		double dy=crd1[i*3+1]-crd2[i*3+1];
		double dz=crd1[i*3+2]-crd2[i*3+2];
		r2 += dx*dx + dy*dy +dz*dz;
	}
	return sqrt(r2/actatom);
}

ActiveSelections * GenChain::setactiveselections_xy(const ForceField *ff){
	std::vector<int> allfixed;
	std::vector<int> mcfixed;
	auto & pset=GenChainControls::getparameterset(controlname_);
	pset.getval("AllFixedLocation",&allfixed);
	pset.getval("MCFixedLocation",&mcfixed);
	std::map<int, std::vector<int>> afs,mfs;
	int naf=allfixed.size()/3;
	for(int i=0;i<naf;i++) {
		std::vector<int> vis;
		vis.push_back(allfixed[i*3]);
		vis.push_back(allfixed[i*3+1]);
		vis.push_back(allfixed[i*3+2]);
		if(afs.find(vis[0])==afs.end()) afs.insert({vis[0],std::vector<int>()});
		for(int j=vis[1];j<vis[2];j++) afs.at(vis[0]).push_back(j);
	}
	int nmf=mcfixed.size()/3;
	for(int i=0;i<nmf;i++) {
		std::vector<int> vis;
		vis.push_back(mcfixed[i*3]);
		vis.push_back(mcfixed[i*3+1]);
		vis.push_back(mcfixed[i*3+2]);
		if(mfs.find(vis[0])==mfs.end()) mfs.insert({vis[0],std::vector<int>()});
		for(int j=vis[1];j<vis[2];j++) mfs.at(vis[0]).push_back(j);
	}
	acts_=std::shared_ptr<ActiveSelections>(new ActiveSelections(ff,afs,mfs));
	return acts_.get();
}*/

GenChain::GenChain(const std::string &controlname, const std::vector<std::vector<NSPproteinrep::FullSite>> &fsss):
	controlname_(controlname) {
	refpdb_=std::shared_ptr<std::vector<std::vector<FullSite>>>(new std::vector<std::vector<FullSite>>(fsss));
	sprefpdb_=std::shared_ptr<StructPlus>(new StructPlus(*refpdb_,true));
	for(auto &c:*refpdb_){
		sctypes_.push_back(std::vector<std::string>());
		softsc_.push_back(std::vector<bool>());
		for(auto &aa:c) {
			std::string upper;
			for(char a:aa.resname()) upper.push_back(toupper(a));
			if(upper==aa.resname()) softsc_.back().push_back(false);
			else softsc_.back().push_back(true);
			try {
				VSCType::getVSCType(upper);
			} catch (std::exception &e){
				std::cout <<"Undefined residue type in PDB:" <<aa.resname()<<std::endl;
				exit(1);
			}
			sctypes_.back().push_back(upper);
		}
	}
}

std::vector<double> GenChain::getcrd(const std::vector<std::vector<FullSite>> &sites,
		std::map<std::pair<int,int>,std::pair<int,int>>&atomseq) const{//chainseq, resdueseq, startseq, length
	std::vector<NSPgeometry::XYZ> res;
	int nres=0;
	double deg=3.14159265/180.0;
	for (int c=0;c<sites.size();++c){
		assert(sites[c].size()==sctypes_[c].size());
		for(int s=0;s<sites[c].size();++s){
			const VSCType &vsc=VSCType::getVSCType(sctypes_[c][s]);
			std::vector<NSPgeometry::XYZ> r(vsc.nscatoms+4);
			const FullSite &fs=sites[c][s];
			r[0]=fs.getcrd("N");
			r[1]=fs.getcrd("CA");
			r[vsc.nscatoms+2]=fs.getcrd("C");
			r[vsc.nscatoms+3]=fs.getcrd("O");
			bool sccomplete{true};
			for(int a=0;a<vsc.nscatoms;++a){
				if(!fs.hasatomcrd(vsc.atomnames[a])) {
					sccomplete=false;
					break;
				}
			}
			if( sccomplete){
				for(int a=0;a<vsc.nscatoms;++a){
					r[a+2]=fs.getcrd(vsc.atomnames[a]);
				}
			} else {
				auto & ics=vsc.internalcrds;
				assert(ics.size()==vsc.nscatoms);
				int aidx=2;
				for(auto &ic:ics){
					r[aidx++]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
								r[ic[2].first],
								ic[0].second, ic[1].second*deg,ic[2].second*deg);
						}
			}
			for(auto &x:r) res.push_back(x);
			atomseq[{c,s}]={nres,r.size()};
			nres+=r.size();
		}
	}
	std::vector<double> crd;
	for (auto &r : res) {
		for (int m = 0; m < 3; ++m)
			crd.push_back(r[m]);
	}
	return crd;
}

std::vector<double> GenChain::getcrd(const std::vector<std::vector<FullSite>> &sites,
		std::vector<std::pair<std::vector<int>,std::string>> &atomseq) {
	std::vector<NSPgeometry::XYZ> res;
	double deg=3.14159265/180.0;
	int na=0;
	for (int c=0;c<sites.size();++c){
		//assert(sites[c].size()==sctypes_[c].size());
		for(int s=0;s<sites[c].size();++s){
			const VSCType &vsc=VSCType::getVSCType(sites[c][s].resname());
			std::vector<NSPgeometry::XYZ> r(vsc.nscatoms+4);
			const FullSite &fs=sites[c][s];
			r[0]=fs.getcrd("N");
			r[1]=fs.getcrd("CA");
			r[vsc.nscatoms+2]=fs.getcrd("C");
			r[vsc.nscatoms+3]=fs.getcrd("O");
			bool sccomplete{true};
			for(int a=0;a<vsc.nscatoms;++a){
				if(!fs.hasatomcrd(vsc.atomnames[a])) {
					sccomplete=false;
					break;
				}
			}
			if( sccomplete){
				for(int a=0;a<vsc.nscatoms;++a){
					r[a+2]=fs.getcrd(vsc.atomnames[a]);
				}
			} else {
				auto & ics=vsc.internalcrds;
				assert(ics.size()==vsc.nscatoms);
				int aidx=2;
				for(auto &ic:ics){
					r[aidx++]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
								r[ic[2].first],
								ic[0].second, ic[1].second*deg,ic[2].second*deg);
						}
			}
			for(auto &x:r) res.push_back(x);
			std::vector<std::string> ans{"N","CA"};
			for(auto &s1:vsc.atomnames) ans.push_back(s1);
			ans.push_back("C");
			ans.push_back("O");
			assert(r.size()==ans.size());
			for(int i=0;i<r.size();i++) {
				atomseq.push_back({{c,s},ans[i]});
			}
		}
	}
	std::vector<double> crd;
	for (auto &r : res) {
		for (int m = 0; m < 3; ++m)
			crd.push_back(r[m]);
	}
	return crd;
}


