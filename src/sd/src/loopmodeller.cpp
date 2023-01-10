/*
 * loopmodeller.cpp
 *
 *  Created on: Jul 12, 2018
 *      Author: xuyang
 */
#include "sd/loopmodeller.h"
#include "backbone/backbonebuilder.h"
#include <iomanip>
using namespace NSPsd;
using namespace NSPgeometry;
using namespace NSPproteinrep;

std::map<char,std::string> AAAbbreviation13 {
	{'G',"GLY"},
	{'A',"ALA"},
	{'V',"VAL"},
	{'L',"LEU"},
	{'I',"ILE"},
	{'S',"SER"},
	{'T',"THR"},
	{'C',"CYS"},
	{'M',"MET"},
	{'D',"ASP"},
	{'E',"GLU"},
	{'N',"ASN"},
	{'Q',"GLN"},
	{'K',"LYS"},
	{'R',"ARG"},
	{'H',"HIS"},
	{'P',"PRO"},
	{'F',"PHE"},
	{'Y',"TYR"},
	{'W',"TRP"}
};

std::vector<NSPproteinrep::FullSite> TargetLoop::extractbaseconfiguration() const {
	int c1=startcp_.cid;
	int c2=endcp_.cid;
	int p1=startcp_.posi;
	int p2=endcp_.posi+1;
	std::vector<NSPproteinrep::FullSite> fss;
	if(c1==c2) {
		for(int i=p1;i!=p2;i++) fss.push_back(basechains_->at(c1)[i]);
	} else {
		for(int i=p1;i!=basechains_->at(c1).size();i++) fss.push_back(basechains_->at(c1)[i]);
		for(int i=0;i!=p2;i++) fss.push_back(basechains_->at(c2)[i]);
	}
	return fss;
}

bool TargetLoop::buildnewconformers() {
	std::set<int> cissites;
	for(int i=0;i<cissites.size();i++) {
		if(newsctypes_[i]=="PRO") cissites.insert(i);
	}
	NSPproteinrep::BackBoneSite bsn;
	NSPproteinrep::BackBoneSite bsc;
	std::vector<std::shared_ptr<std::vector<NSPproteinrep::BackBoneSite>>> bbss;
	if(startcp_.posi<=0) {
		bsc=basechains_->at(endcp_.cid)[endcp_.posi+1].getbackbonesite();
		std::shared_ptr<std::vector<NSPproteinrep::BackBoneSite>> bss = std::make_shared<std::vector<
				NSPproteinrep::BackBoneSite>>(NSPproteinrep::BackBoneBuilder::buildbackwardbackbone(
						newlength_,bsc,std::vector<std::pair<int,int>>(),std::vector<std::pair<int,int>>(),cissites));
		bbss.push_back(bss);
	} else if(endcp_.posi<=0 || endcp_.posi>=basechains_->size()-1) {
		bsn=basechains_->at(startcp_.cid)[startcp_.posi-1].getbackbonesite();
		std::shared_ptr<std::vector<NSPproteinrep::BackBoneSite>> bss = std::make_shared<std::vector<
				NSPproteinrep::BackBoneSite>>(NSPproteinrep::BackBoneBuilder::buildforwardbackbone(
						newlength_,bsn,std::vector<std::pair<int,int>>(),std::vector<std::pair<int,int>>(),cissites));
		bbss.push_back(bss);
	} else {
		bsn=basechains_->at(startcp_.cid)[startcp_.posi-1].getbackbonesite();
		bsc=basechains_->at(endcp_.cid)[endcp_.posi+1].getbackbonesite();
		bbss=NSPproteinrep::BackBoneBuilder::buildlinkers(newlength_,bsn,bsc,
				std::vector<std::pair<int,int>>(),std::vector<std::pair<int,int>>(),cissites);
		int ix=0;
		while(bbss.empty()) {
			bbss=NSPproteinrep::BackBoneBuilder::buildlinkers(newlength_,bsn,bsc,
					std::vector<std::pair<int,int>>(),std::vector<std::pair<int,int>>(),cissites);
			ix++;
			if(ix==10) return false;
		}
	}
	for(auto &bss:bbss) {
		assert(newsctypes_.size()==bss->size());
		std::shared_ptr<std::vector<NSPproteinrep::FullSite>> fss=std::shared_ptr<
				std::vector<NSPproteinrep::FullSite>>(new std::vector<NSPproteinrep::FullSite>());
		char cid=basechains_->at(startcp_.cid)[startcp_.posi].chainid();
		int pid=startcp_.posi;
		for(int i=0;i<bss->size();i++) {
			bss->at(i).chainid=cid;
			bss->at(i).resid=pid++;
			bss->at(i).resname=newsctypes_[i];
			fss->push_back(NSPproteinrep::make_fullsite(bss->at(i)));
			fss->back().resseq()=pid-1;
		}
		newconformers_.push_back(fss);
	}
	return true;
}
/*
void NSPsd::defineloopcontrol(std::string name,const std::vector<std::string> &controllines){
	std::map<std::string,double> doublepars{ {"SampTemperature",1.0}, {"SampPHIPSI",1.0}, {"SampLocal",1.0},
		{"SampPack",1.0}, {"SampSCConf",1.0}, {"SampSCPack",1.0}, {"SampLocalHB",1.0}, {"OptTemperature",1.0},
		{"OptPHIPSI",1.0}, {"OptLocal",1.0}, {"OptPack",1.0}, {"OptSCConf",1.0}, {"OptSCPack",1.0},
		{"OptLocalHB",1.0}, {"RMSDCut",3.0}, {"RMSDForce",10.0}, {"RMSDDifference",2.0}, {"EnergyDifference",2.0},
		{"SampRMSDCut",0.0}, {"LowRMSDCut",1.0}, {"HighRMSDCut",1.0}, {"EnergyCutOff",1.0}, {"SampEnergyProp",1.0} };
	std::map<std::string,std::vector<std::string>> stringvecpars{ {"LoopSeq",std::vector<std::string>()} };
	std::map<std::string,std::vector<double>> doublevecpars;
	std::map<std::string,std::string> stringpars{ {"InitialPDB",""}, {"ControlGenChain","sdffcontrol"},
		{"OptMethod","EM"}, {"RMSDCalculation","AllAtom"}, {"Pattern","AllAtom"}, {"SampMethod","RMSDCut"} };
	std::map<std::string,std::vector<int>> intvecpars{{"Loop",{0,0,0,0}}};
	std::map<std::string,int>intpars{ {"SampNumber",0}, {"SampMaxCycle",0}, {"SampMaxStep",0}, {"SampInterval",0},
		{"OptMaxStep",1000}, {"OptAlgorithm",1000}, {"OptInterval",1000} };
	LoopControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	int nsuccess=LoopControls::adjustvalues(name,controllines);
	if(nsuccess!= controllines.size()) {
		exit(1);
	}
}*/

void NSPsd::defineloopcontrol(std::string name,const std::vector<std::string> &controllines){
	std::map<std::string,double> doublepars{ {"EnergyDifference",0.0}, {"RMSDForce",0.0}, {"SampRMSDCut",0.0},
		{"SampWeightProp",0.0}, {"OptTemperature",0.0}, {"SampTemperature",0.0} };
	std::map<std::string,std::vector<std::string>> stringvecpars{ {"Loop",{}}, {"LoopSeq",{}} };
	std::map<std::string,std::vector<double>> doublevecpars;
	std::map<std::string,std::string> stringpars{ {"InitialPDB",""}, {"Pattern","AllAtom"}, {"ControlGenChain",""},
		{"OptMethod",""}, {"RMSDCalculation",""}, {"Path",""} };
	std::map<std::string,std::vector<int>> intvecpars{};
	std::map<std::string,int>intpars{ {"OptAlgorithm",0}, {"OptMaxStep",0}, {"OptInterval",0}, {"SampMaxStep",0},
		{"SampMaxCycle",0}, {"SampNumber",0} };
	LoopControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	int nsuccess=LoopControls::adjustvalues(name,controllines);
	if(nsuccess!= controllines.size()) {
		exit(1);
	}
}

void NSPsd::genchainreadcontrols_loop(const std::string &filename,std::string name){
	NSPdataio::ControlFile cf;
	cf.readfile(filename);
	std::vector<std::string> sdcontrolines=cf.getcontrolines("SD");
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	std::vector<std::string> genchaincontrollines=cf.getcontrolines("GenChain");
	std::vector<std::string> loopcontrollines=cf.getcontrolines("LOOP");
	definesdcontrol(name+"_sd",sdcontrolines);
	defineforcefieldcontrol(name+"_ff",ffcontrolines);
	definegenchaincontrol(name+"_genchain",genchaincontrollines);
	defineloopcontrol(name+"_loop",loopcontrollines);
}

bool ArchivedLoops::addmodel(const std::vector<double> &model, const std::vector<double> &energy,
		const std::vector<double> &rmsdrecord, std::set<int> &removed, int maxnum) {
	//true means model is added, false means not
	bool added{true};
	removed.clear();
	for(int i=0;i<confcrds_.size();i++) {
		if(rmsdrecord[i]<rmsdcut_) {
			removed.insert(i);
			if(energy[0]>energies_[i][0]) added=false;
		}
	}
	if(!added) return false;
	confcrds_.push_back(model);
	energies_.push_back(energy);

	//remove point in removed
	std::vector<int> rd;
	for(auto r:removed) rd.push_back(r);
	int ix=0;
	int newsize=confcrds_.size()-removed.size();
	int newadded;
	for(int i=newsize;i<confcrds_.size();i++) {
		if(removed.find(i)==removed.end()) {
			assert(rd[ix]<newsize);
			confcrds_[rd[ix]] = confcrds_[i];
			energies_[rd[ix]] = energies_[i];
			ix++;
			if(i==confcrds_.size()-1) newadded=rd[ix];
		}
	}
	confcrds_.resize(newsize);
	energies_.resize(newsize);
	if(newsize<=maxnum) return true;

	//remove point more than maxnum
	double bige=-100000.0;
	for(int i=0;i<newsize;i++) {
		if(energies_[i][0]>bige) {
			ix=i;
			bige=energies_[i][0];
		}
	}
	confcrds_[ix] = confcrds_.back();
	energies_[ix] = energies_.back();
	confcrds_.resize(newsize-1);
	energies_.resize(newsize-1);
	if(ix==newadded) return false;
	return true;
}

LoopModeller::LoopModeller(NSPdataio::ParameterSet &pset):pset_(pset) {
	std::vector<std::pair<char,std::pair<int,char>>> removedregion;
	std::vector<std::string> read_loop;
	pset_.getval("Loop",&read_loop);
	for(int i=0;i<read_loop.size();i+=3) {
		char chid=read_loop[i][0];
		int lid=std::stoi(read_loop[i+1]);
		char iid=read_loop[i+2][0];
		if(read_loop[i+2][0]=='_') iid=' ';
		removedregion.push_back({chid,{lid,iid}});
	}
	std::string initialpdb;
	pset_.getval("InitialPDB",&initialpdb);
	std::vector<std::pair<int,int>> loopregion;
	basechains_=NSPproteinrep::readfullsitesfrompdb_xy(initialpdb,removedregion,loopregion);
	std::string pattern;
	pset_.getval("Pattern",&pattern);
	if(pattern=="MCOnly") {
		for(auto & ch:basechains_) {
			for(NSPproteinrep::FullSite &fs:ch) {
				fs.covert2backbone();
			}
		}
	}
	pset_.getval("Path",&path_);
	if(path_.back()!='/') path_+='/';


	for(int i=0;i<loopregion.size();i+=2) {
		addtargetloop(loopregion[i].first,loopregion[i].second,loopregion[i+1].first,loopregion[i+1].second);
	}
	std::vector<std::string> read_seqs;
	pset_.getval("LoopSeq",&read_seqs);
	for(std::string &sq:read_seqs) {
		std::vector<std::string> str3;
		for(int i=0;i<sq.size();i++) {
			str3.push_back(AAAbbreviation13.at(sq[i]));
		}
		loopseq_.push_back(str3);
	}
	if(pattern=="MCOnly") {
		for(auto & ls:loopseq_) {
			for(std::string & s:ls) s="GLY";
		}
	}
	if(loopseq_.size()!=loopregion.size()/2) {
		std::cout <<"LOOP NUMBER IS NOT CORRESPONDED TO LOOP SEQUENCE NUM!!!" <<std::endl;
		exit(1);
	}

	std::string controlname;
	pset_.getval("ControlGenChain",&controlname);
	std::vector<std::vector<NSPproteinrep::FullSite>> initfsss;
	/*std::vector<std::set<int>> cisinloop;
	while(true) {
		std::vector<std::vector<bool>> newfixedsite;
		initfsss=getainstance(newfixedsite);
		if(cisinloop.empty()) {
			for(int i=0;i<initfsss.size();i++) {
				cisinloop.push_back(std::set<int>());
				for(int j=0;j<initfsss[i].size()-1;j++) {
					if(newfixedsite[i][j] && newfixedsite[i][j+1]) continue;
					if(initfsss[i][j+1].resname()=="PRO") cisinloop.back().insert(j+1);
				}
			}
		}
		bool allisright{true};
		for(int i=0;i<initfsss.size();i++){
			std::vector<BackBoneSite> bss=backbone(initfsss[i]);
			std::vector<int> cisch=findcissites(bss);
			for(int c:cisch) {
				if(c==0) continue;
				if(newfixedsite[i][c-1] && newfixedsite[i][c]) continue;
				if(cisinloop[i].find(c)!=cisinloop[i].end()) continue;
				allisright=false;
				std::cout <<"Initial Loop Has Wrong Omiga!" <<std::endl;
				break;
			}
			if(!allisright) break;
		}
		if(allisright) break;
	}*/
	std::vector<std::vector<bool>> newfixedsite;
	initfsss=getainstance(newfixedsite);
	std::vector<std::set<int>> notcisinloop;
	for(int i=0;i<initfsss.size();i++) {
		notcisinloop.push_back(std::set<int>());
		for(int j=0;j<initfsss[i].size()-1;j++) {
			if(newfixedsite[i][j] && newfixedsite[i][j+1]) continue;
			if(initfsss[i][j+1].resname()!="PRO") notcisinloop.back().insert(j+1);
		}
	}
	changeline(controlname+"_genchain");
	genchain_=std::make_shared<GenChain>(GenChain(controlname+"_genchain", initfsss));
	std::vector<double> initcrd = genchain_->getcrd(initfsss);
	for(double &c:initcrd) c *= A2NM;
	SDRun::SDRunIn sdrunin(initcrd,controlname);
	int rng = NSPdstl::RandomEngine<>::getinstance().intrng(1,1000)();
	sdrun_ = std::shared_ptr<SDRun>(new SDRun(genchain_->make_sdrun(sdrunin,rng,notcisinloop)));
	std::string method;
	pset_.getval("OptMethod",&method);
	if(method=="EM") {
		ff_=std::shared_ptr<ForceField>(new ForceField(genchain_->make_forcefield(controlname+"_ff",notcisinloop)));
		emrun_ = std::shared_ptr<EMRun> (new EMRun(genchain_.get(),ff_.get()));
	}
}

void LoopModeller::native() {
	std::vector<std::vector<NSPproteinrep::FullSite>> ntfs=basechains_;
	std::vector<double> crds=genchain_->getcrd(ntfs);
	for(auto &c:crds) c*=A2NM;
	std::vector<double> nates = sdrun_->confene(crds);
	std::vector<double> nativecrd=crds;
	std::vector<double> es=opt(crds);
	std::ofstream ofs(path_+"native_opt.pdb");
	genchain_->writepdb(crds,ofs,10.0);
	ofs.close();
	ofs.open(path_+"native.ene");
	for(double e:nates) ofs <<std::fixed <<std::setprecision(1) <<e <<'\t';
	ofs <<std::endl;
	for(double e:es) ofs <<std::fixed <<std::setprecision(1) <<e <<'\t';
	ofs <<std::fixed <<std::setprecision(3) <<ff_->getmcrmsd(crds,nativecrd,sdrun_->forceoff(),10.0) <<std::endl;
	ofs.close();
}

void LoopModeller::changeline(std::string controlname) {
	std::string line="AllFixedSites\t=\t";
	std::vector<std::vector<bool>> fixed;
	getainstance(fixed);
	for(int i=0;i<fixed.size();i++) {
		line = line + "chain" + std::to_string(i) + " ";
		for(int j=0;j<fixed[i].size();j++) if(fixed[i][j]) line = line + std::to_string(j) + " ";
	}
	auto & pset=GenChainControls::getparameterset(controlname);
	pset.readline(line);
}

void LoopModeller::remove_prev(std::map<std::pair<int,int>,std::vector<NSPproteinrep::FullSite>> &orifs, int ch, int po) {
	int lc,lp;
	for(auto & of:orifs) {
		if(of.first.first==ch && po>=of.first.second && po<=of.first.second+of.second.size()) {
			std::vector<NSPproteinrep::FullSite> fss;
			for(int i=0;i<of.second.size();i++) {
				if(of.first.second+i>po) fss.push_back(of.second[i]);
			}
			of.second=fss;
			lc=of.first.first;
			lp=of.first.second;
			//of.first.second=po+1;
			break;
		}
	}
	auto r=orifs.at({lc,lp});
	orifs.erase({lc,lp});
	orifs.insert({{lc,po+1},r});
}

void LoopModeller::remove_back(std::map<std::pair<int,int>,std::vector<NSPproteinrep::FullSite>> &orifs, int ch, int po) {
	for(auto & of:orifs) {
		if(of.first.first==ch && po>=of.first.second && po<=of.first.second+of.second.size()) {
			std::vector<NSPproteinrep::FullSite> fss;
			for(int i=0;i<of.second.size();i++) {
				if(of.first.second+i<po) fss.push_back(of.second[i]);
			}
			of.second=fss;
			break;
		}
	}
}

void LoopModeller::remove_middle(std::map<std::pair<int,int>,std::vector<NSPproteinrep::FullSite>> &orifs, int ch, int p1, int p2) {
	std::vector<NSPproteinrep::FullSite> fss;
	for(auto & of:orifs) {
		if(of.first.first==ch && p1>=of.first.second && p2<=of.first.second+of.second.size()) {
			std::vector<NSPproteinrep::FullSite> fssret;
			for(int i=0;i<of.second.size();i++) {
				if(of.first.second+i+1<p1) fssret.push_back(of.second[i]);
			}
			//of.second=fssret;
			//fss.clear();
			for(int i=0;i<of.second.size();i++) {
				if(of.first.second+i>=p2) fss.push_back(of.second[i]);
			}
			of.second=fssret;
			break;
		}
	}
	orifs.insert({{ch,p2+1},fss});
}

std::vector<std::vector<NSPproteinrep::FullSite>> LoopModeller::getainstance(std::vector<std::vector<bool>> &newfixed) {
	newfixed.clear();
	std::map<std::pair<int,int>,std::vector<NSPproteinrep::FullSite>> orifs;
	for(int i=0;i<basechains_.size();i++) orifs.insert({{i,0},basechains_[i]});
	std::map<std::vector<int>,std::vector<NSPproteinrep::FullSite>> loopchains;
	for(int i=0;i<loops_.size();i++) {
		int c1=loops_[i].sten()[0], c2=loops_[i].sten()[2];
		int p1=loops_[i].sten()[1], p2=loops_[i].sten()[3];
		if(c1==c2) {
			if(p1<=0) {
				remove_prev(orifs,c1,p2);
			} else if(p2<=0 || p2>=basechains_[c1].size()-1) {
				remove_back(orifs,c1,p1);
			} else {
				remove_middle(orifs,c1,p1,p2);
			}
		} else {
			remove_back(orifs,c1,p1);
			remove_prev(orifs,c2,p2);
		}
		std::vector<int> chpos;
		chpos.push_back(c1);
		chpos.push_back(p1);
		chpos.push_back(c2);
		chpos.push_back(p2);
		loopchains.insert({chpos,*(loops_[i].popconformer(loopseq_[i]))});
	}

	std::map<std::pair<int,int>,std::vector<bool>> fixori;
	for(auto & of:orifs) fixori.insert({of.first,std::vector<bool>(of.second.size(),true)});

	for(auto & ch:loopchains) {
		if(ch.first[1]<=0) {
			int oc=ch.first[0];
			int op=10000;
			for(auto &of:orifs) {
				if(of.first.first==oc) {
					if(of.first.second<op) op=of.first.second;
				}
			}
			std::vector<NSPproteinrep::FullSite> fss=orifs.at({oc,op});
			orifs.at({oc,op}).clear();
			fixori.at({oc,op}).clear();
			for(auto & fs:ch.second) {
				orifs.at({oc,op}).push_back(fs);
				fixori.at({oc,op}).push_back(false);
			}
			for(auto & fs:fss) {
				orifs.at({oc,op}).push_back(fs);
				fixori.at({oc,op}).push_back(true);
			}
		} else {
			std::vector<NSPproteinrep::FullSite> fss=ch.second;
			std::vector<bool> lb(fss.size(),false);
			int oc,op;
			for(auto &of:orifs) {
				if(of.first.first==ch.first[2] && of.first.second-1==ch.first[3]) {
					for(int i=0;i<of.second.size();i++) {
						fss.push_back(of.second[i]);
						lb.push_back(true);
					}
					oc=of.first.first;
					op=of.first.second;
					break;
				}
			}
			orifs.erase({oc,op});
			fixori.erase({oc,op});
			for(auto & of:orifs) {
				if(of.first.first==ch.first[0] && of.first.second+of.second.size()+1==ch.first[1]) {
					for(int i=0;i<fss.size();i++) {
						of.second.push_back(fss[i]);
						fixori.at(of.first).push_back(lb[i]);
					}
					break;
				}
			}
		}
	}

	std::vector<std::vector<NSPproteinrep::FullSite>> fsss;
	for(auto &of:orifs) {
		fsss.push_back(of.second);
		newfixed.push_back(fixori.at(of.first));
	}
	return fsss;
}

std::vector<double> LoopModeller::opt(std::vector<double> &crd) {
	static bool readed{false};
	static std::string method;
	static int opt_algorithm;
	static int opt_maxstep;
	static int check_interval;
	static double ene_diff;
	if(!readed) {
		readed = true;
		pset_.getval("OptMethod",&method);
		pset_.getval("OptMaxStep",&opt_maxstep);
		pset_.getval("OptAlgorithm",&opt_algorithm);
		pset_.getval("OptInterval",&check_interval);
		pset_.getval("EnergyDifference",&ene_diff);
	}
	if(method=="EM") {
		double ene;
		emrun_->run(opt_algorithm,opt_maxstep,&crd,&ene,check_interval,ene_diff);
		crd = emrun_->crdbck();
		return emrun_->potenergies();
	} else if(method=="SD") {
		sdrun_->resetstate(crd);
		sdrun_->run_opt(opt_maxstep,check_interval,ene_diff);
		crd = sdrun_->state().crd;
		return sdrun_->potenergies();
	} else {
		std::cout <<"Optimization Method Do Not Exist!" <<std::endl;
		exit(1);
	}
}

bool LoopModeller::knot(const std::vector<NSPproteinrep::FullSite> & chain) {
	std::vector<NSPgeometry::XYZ> cas;
	for(const auto &fs:chain) cas.push_back(fs.getcrd("CA"));
	std::vector<int> updown; //1-up, 2-down
	for(int i=0;i<cas.size()-1;i++) {
		for(int j=i+1;j<cas.size()-1;j++) {
			int ix=NSPgeometry::projection(cas[i],cas[i+1],cas[j],cas[j+1]);
			if(ix==0) continue;
			updown.push_back(ix);
		}
	}
	for(int i=0;i<updown.size()-2;i++) {
		if(updown[i]==1 && updown[i+1]==2 && updown[i+2]==1) return false;
		if(updown[i]==2 && updown[i+1]==1 && updown[i+2]==2) return false;
	}
	return true;
}

bool LoopModeller::topoisright(const std::vector<std::vector<NSPproteinrep::FullSite>> & conf) {
	static std::vector<std::pair<int,int>> loopregion;
	if(loopregion.empty()) {
		std::vector<std::pair<char,std::pair<int,char>>> removedregion;
		std::vector<std::string> read_loop;
		pset_.getval("Loop",&read_loop);
		for(int i=0;i<read_loop.size();i+=3) {
			char chid=read_loop[i][0];
			int lid=std::stoi(read_loop[i+1]);
			char iid=read_loop[i+2][0];
			if(read_loop[i+2][0]=='_') iid=' ';
			removedregion.push_back({chid,{lid,iid}});
		}
		std::string initialpdb;
		pset_.getval("InitialPDB",&initialpdb);
		NSPproteinrep::readfullsitesfrompdb_xy(initialpdb,removedregion,loopregion);
	}
	for(int i=0;i<loopregion.size();i+=2) {
		assert(loopregion[i].first==loopregion[i+1].first);
		std::vector<FullSite> fss;
		int cid=loopregion[i].first;
		for(int j=loopregion[i].second;j<=loopregion[i+1].second;j++) {
			fss.push_back(conf[cid][j]);
		}
		if(!knot(fss)) {
			std::cout <<"\tInitial Structure Has Wrong Topo!" <<std::endl;
			return false;
		}
	}
	return true;
}

void LoopModeller::makeinitialconf(std::vector<double> &crd, std::vector<double> &ene) {
	do {
		std::cout <<"Get Initial Structure ..." <<std::endl;
		std::vector<std::vector<NSPproteinrep::FullSite>> fsss=getainstance();
		crd = genchain_->getcrd(fsss);
		for (auto &c : crd) c *= A2NM;
		ene = opt(crd);
		break;
	} while(!topoisright(crd));
}
void LoopModeller::sample(const std::vector<double> &initcrd, const std::vector<double> &initene) {
	static bool readed{false};
	static int sample_maxstep{-1};
	static double force_index{0.0};
	static int maxcycle;
	static int nsample;
	static bool mconly{false};
	static double rmsdcut;
	if(!readed) {
		readed=true;
		pset_.getval("SampMaxStep",&sample_maxstep);
		pset_.getval("RMSDForce",&force_index);
		pset_.getval("SampMaxCycle",&maxcycle);
		pset_.getval("SampNumber",&nsample);
		std::string pattern;
		pset_.getval("RMSDCalculation",&pattern);
		if(pattern=="MCOnly") mconly = true;
		pset_.getval("SampRMSDCut",&rmsdcut);
		int lthloop=0;
		for(const auto &ls:loopseq_) lthloop += ls.size();
		rmsdcut = (double)lthloop*rmsdcut*A2NM;
	}

	std::vector<std::vector<double>> & confs = archivedloops_.confcrds();
	confs.clear();
	confs.push_back(initcrd);
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	std::vector<double> rdrd;
	std::vector<std::vector<double>> enessample;
	for(int c=0,ix=0; c<maxcycle && confs.size()<nsample; c++) {
		ix = rng.intrng(0,confs.size()-1)();
		if(!sdrun_->resetstate(confs[ix])) {
			std::cout <<"Initial State Shake Fail!" <<std::endl;
			if(confs.size()==1) return;
			continue;
		}
		if(!sdrun_->run_sample_rmsdcut(sample_maxstep,confs,force_index,rmsdcut,mconly,rdrd,enessample)) {
			if(confs.size()==1) return;
			continue;
		}
		confs.push_back(sdrun_->state().crd);
		std::cout <<"Get Conformer : " <<confs.size() <<std::endl;
	}
}

void LoopModeller::run() {
	std::vector<double> initcrd;
	std::vector<double> initene;
	double opttemp;
	pset_.getval("OptTemperature",&opttemp);
	sdrun_->changetemperature(opttemp);
	makeinitialconf(initcrd,initene);

	double prop;
	pset_.getval("SampWeightProp",&prop);
	sdrun_->changeweight(prop);
	double samptemp;
	pset_.getval("SampTemperature",&samptemp);
	sdrun_->changetemperature(samptemp);

	print(initcrd,initene,0);
	sample(initcrd,initene);

	while(archivedloops_.confcrds().size()==1) {
		makeinitialconf(initcrd,initene);
		print(initcrd,initene,0);
		sample(initcrd,initene);
	}

	typedef std::vector<double> VD;
	std::vector<VD> & confs = archivedloops_.confcrds();
	sdrun_->changeweight(1.0/prop);
	sdrun_->changetemperature(opttemp);
	for(int i=1;i<confs.size();i++) {
		std::cout <<"Optimization: " <<i <<std::endl;
		std::vector<double> es = opt(confs[i]);
		print(confs[i],es,i);
	}
	return;
}

void LoopModeller::print(const std::vector<double> &crd, const std::vector<double> &es, int i) {
	static std::vector<double> nativecrd;
	if(nativecrd.empty()) {
		nativecrd=genchain_->getcrd(basechains_);
		for(double &c:nativecrd) c *= A2NM;
	}
	std::ofstream ofs(path_+std::to_string(i)+".pdb");
	genchain_->writepdb(crd,ofs,10.0);
	ofs.close();
	if(i==0) ofs.open(path_+"record.dat");
	else ofs.open(path_+"record.dat",std::ofstream::app);
	ofs <<i <<'\t';
	for(double e:es) ofs <<std::fixed <<std::setprecision(1) <<e <<'\t';
	ofs <<std::fixed <<std::setprecision(3) <<ff_->getmcrmsd(crd,nativecrd,sdrun_->forceoff(),10.0) <<std::endl;
	ofs.close();
}



















/*
std::vector<std::vector<NSPproteinrep::FullSite>> LoopModeller::getainstance(std::vector<std::vector<bool>> &newfixed) {
	std::vector<std::vector<bool>> poslb;
	std::vector<bool> chlb;
	for(int i=0;i<basechains_.size();i++) {
		std::vector<bool> lb(basechains_[i].size(),true);
		poslb.push_back(lb);
		chlb.push_back(true);
	}
	std::vector<int> rds;
	pset_.getval("Loop",&rds);
	std::map<std::vector<int>,std::vector<NSPproteinrep::FullSite>> chains;
	int nlp=rds.size()/4;
	for(int i=0;i<nlp;i++) {
		int c1=rds[i*4], c2=rds[i*4+2];
		int p1=rds[i*4+1], p2=rds[i*4+3];
		if(c1==c2) for(int j=p1;j<=p2;j++) poslb[c1][j]=false;
		else {
			for(int j=p1;j<poslb[c1].size();j++) poslb[c1][j]=false;
			for(int j=0;j<=p2;j++) poslb[c2][j]=false;
		}
		std::vector<int> chpos;
		chpos.push_back(c1);
		chpos.push_back(p1);
		chpos.push_back(c2);
		chpos.push_back(p2);
		chains.insert({chpos,*(loops_[i].popconformer(loopseq_[i]))});
	}
	std::vector<std::vector<NSPproteinrep::FullSite>> newchain;
	while(!chains.empty()) {
		for(int i=0;i<basechains_.size();i++) {
			if(!chlb[i]) continue;
			if(!poslb[i][0]) continue;
			std::vector<NSPproteinrep::FullSite> fss;
			std::vector<bool> nfx;
			int ix=i;
			do {
				chlb[ix]=false;
				bool isend{true};
				for(int j=0;j<basechains_[ix].size();j++) {
					if(poslb[ix][j]) {
						fss.push_back(basechains_[ix][j]);
						nfx.push_back(true);
					}
					//std::vector<NSPproteinrep::FullSite> lpnew;
					std::vector<int> todelete;
					for(auto &ch:chains) {
						if(ix!=ch.first[0] || j!=ch.first[1]) continue;
						todelete=ch.first;
						//lpnew=ch.second;
						for(NSPproteinrep::FullSite &fs:ch.second) {
							fss.push_back(fs);
							nfx.push_back(false);
						}
						break;
					}
					if(todelete.empty()) continue;
					chains.erase(todelete);
					if(todelete[0]==todelete[2]) continue;
					ix=todelete[2];
					isend=false;
					break;
				}
				if(isend) break;
			} while(true);
			newchain.push_back(fss);
			newfixed.push_back(nfx);
		}
	}
	return newchain;
}*/


/*

void LoopModeller::nativetest(std::string fnene, std::string fnpdb, double prop) {
	std::vector<double> initcrd1;
	std::vector<double> initene1;
	makeinitialconf(initcrd1,initene1);
	double temp;
	pset_.getval("SampTemperature",&temp);
	sdrun_->changetemperature(temp);
	sdrun_->changeweight(prop);
	sample(initcrd1,initene1,false);
	std::vector<std::vector<double>> & confs = archivedloops_.confcrds();
	std::ofstream ofs;
	std::vector<double> crd2=genchain_->getcrd(basechains_);
	for(double &c:crd2) c *= A2NM;
	pset_.getval("OptTemperature",&temp);
	sdrun_->changetemperature(temp);
	sdrun_->changeweight(1.0/prop);
	for(int i=0;i<confs.size();i++) {
		std::vector<double> es = opt(confs[i],false);
		ofs.open("optene",std::ofstream::app);
		for(double e:es) ofs <<std::fixed <<std::setprecision(1) <<e <<'\t';
		ofs <<std::endl;
		ofs.close();
		ofs.open(std::to_string(i)+".pdb.opt");
		genchain_->writepdb(confs[i],ofs,10.0);
		ofs.close();
		ofs.open("mcrmsd",std::ofstream::app);
		ofs <<ff_->getmcrmsd(confs[i],crd2,sdrun_->forceoff(),10.0) <<std::endl;
		ofs.close();
		//ofs.open("allrmsd",std::ofstream::app);
		//ofs <<ff_->getallrmsd(confs[i],crd2,sdrun_->forceoff(),10.0) <<std::endl;
		//ofs.close();
	}
	return;
	std::vector<double> initcrd = genchain_->getcrd(basechains_);
	for(double &c:initcrd) c *= A2NM;
	std::vector<bool> forceoff = sdrun_->forceoff();
	std::shared_ptr<NeighborList> nbl = std::shared_ptr < NeighborList > (new NeighborList(initcrd, *ff_, &forceoff));
	std::vector<double> pot;
	ff_->forces(initcrd,*nbl,&pot,*(genchain_->getactiveselections()));
	std::vector<double> crd=initcrd;
	emrun_->nativecrd() = crd;
	emrun_->atomfixed() = sdrun_->forceoff();
	std::vector<double> fe=opt(crd);
	std::vector<double> enecurve = emrun_ ->enecurve();
	std::vector<double> rmsdmccurve = emrun_ ->rmsdmc();
	std::vector<double> rmsdallcurve = emrun_ ->rmsdall();

	std::ofstream ofs1(fnene);
	ofs <<"Native:"; for(double d:pot) ofs <<'\t' <<std::fixed <<std::setprecision(3) <<d;
	ofs <<std::endl;
	ofs <<"Stable:"; for(double d:fe) ofs <<'\t' <<std::fixed <<std::setprecision(3) <<d;
	ofs <<std::endl <<std::endl;
	ofs <<"RMSD_MC:\t" <<ff_->getmcrmsd(initcrd,crd,sdrun_->forceoff(),10.0) <<std::endl;
	ofs <<"RMSD_ALL:\t" <<ff_->getallrmsd(initcrd,crd,sdrun_->forceoff(),10.0) <<std::endl <<std::endl;
	ofs <<pot[0] <<std::endl;
	for(int i=0;i<enecurve.size();i++) ofs <<enecurve[i] <<'\t' <<rmsdmccurve[i] <<'\t' <<rmsdallcurve[i] <<std::endl;
	ofs.close();
	ofs.open(fnpdb);
	genchain_->writepdb(crd,ofs,10.0);
	ofs.close();
}

void LoopModeller::printresults(std::string dir, std::string fnene) {
	if(dir.back()!='/') dir+='/';
	std::ofstream ofs;
	std::vector<std::vector<double>> &crds=archivedloops_.confcrds();
	int ix=0;
	for(std::vector<double> &cs:crds) {
		ix++;
		std::string fn=dir+std::to_string(ix)+".pdb";
		ofs.open(fn);
		genchain_->writepdb(cs, ofs, 10.0);
		ofs.close();
	}
	std::vector<std::vector<double>> &enes=archivedloops_.energys();
	ofs.open(dir+fnene);
	for(int i=0;i<enes.size();i++) {
		for(double &e:enes[i]) ofs <<std::fixed <<std::setprecision(3) <<e <<'\t';
		ofs <<std::endl;
	}
	ofs.close();
}



void LoopModeller::sample_rmsd_ene(const std::vector<double> &initcrd, const std::vector<double> &initene) {
	static bool readed{false};
	static int sample_maxstep{-1};
	static double force_index{0.0};
	static int sample_record{-1};
	static double rmsdchangecutoff{0.0};
	static int maxcycle;
	static int nsample;
	static bool mconly{false};
	static std::string samplingmethod;
	static double rmsdcut;
	static double samptemp;
	static double lowrmsdcut;
	static double highrmsdcut;
	static double enecut;
	static double prop;
	if(!readed) {
		pset_.getval("SampMaxStep",&sample_maxstep);
		pset_.getval("RMSDForce",&force_index);
		pset_.getval("SampInterval",&sample_record);
		pset_.getval("RMSDDifference",&rmsdchangecutoff);
		rmsdchangecutoff *= A2NM;
		pset_.getval("SampMaxCycle",&maxcycle);
		pset_.getval("SampNumber",&nsample);
		std::string pattern;
		pset_.getval("RMSDCalculation",&pattern);
		if(pattern=="MCOnly") mconly = true;
		pset_.getval("SampMethod",&samplingmethod);
		pset_.getval("SampRMSDCut",&rmsdcut);
		int lthloop=0;
		for(const auto &ls:loopseq_) lthloop += ls.size();
		rmsdcut = (double)lthloop*rmsdcut*A2NM;
		pset_.getval("SampTemperature",&samptemp);
		pset_.getval("LowRMSDCut",&lowrmsdcut);
		pset_.getval("HighRMSDCut",&highrmsdcut);
		pset_.getval("EnergyCutOff",&enecut);
		pset_.getval("SampEnergyProp",&prop);
		lowrmsdcut = (double)lthloop*lowrmsdcut*A2NM;
		highrmsdcut = (double)lthloop*highrmsdcut*A2NM;
	}
	sdrun_->changetemperature(samptemp);
	sdrun_->changeweight(prop);

	std::vector<std::vector<double>> & confs = archivedloops_.confcrds();
	std::vector<std::vector<double>> & enes = archivedloops_.energys();
	confs.push_back(initcrd);
	enes.push_back(initene);
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	std::ofstream ofs;
	ofs.open("0.pdb");
	genchain_->writepdb(initcrd,ofs,10.0);
	ofs.close();
	for(int c=0,ix=0; c<maxcycle && confs.size()<nsample; c++) {
		ix = rng.intrng(0,confs.size()-1)();
		sdrun_->resetstate(confs[ix]);
		//if(!sdrun_->run_sample(sample_maxstep,confs,force_index,sample_record,rmsdchangecutoff,mconly)) continue;
		{
			std::vector<std::pair<int,double>> minrmsdrecord;
			std::vector<std::vector<double>> enesrecord;
			if(!sdrun_->run_sample_rmsdcut_ene(confs,lowrmsdcut,highrmsdcut,enecut,force_index,
					mconly,sample_maxstep,minrmsdrecord,enesrecord)) continue;
			ofs.open("sampene",std::ofstream::app);
			for(auto &e:sdrun_->potenergies()) {
				ofs <<std::fixed <<std::setprecision(1) <<e <<'\t';
			}
			ofs <<std::endl;
			ofs.close();
			ofs.open(std::to_string(confs.size())+".pdb");
			genchain_->writepdb(sdrun_->state().crd,ofs,10.0);
			ofs.close();
		}
		confs.push_back(sdrun_->state().crd);
		enes.push_back(sdrun_->potenergies());
	}
}



*/























































































/*
std::vector<XYZ> LoopModeller::buildbbsite(const BackBoneSite &base,
		bool backward, double phi, double psi, bool iscis) { // iscis is according to base
	XYZ bn=base.ncrd();
	XYZ bca=base.cacrd();
	XYZ bc=base.ccrd();
	std::vector<XYZ> newcrds(4); //order is N, CA, C, O
	if(backward) {
		newcrds[0] = InternaltoXYZ(); // plane angle
		if(iscis) {
			newcrds[1] = InternaltoXYZ(); // plane angle
		} else {
			newcrds[1] = InternaltoXYZ(); // plane angle
		}
		newcrds[2] = InternaltoXYZ(); // dihedral angle
		newcrds[3] = InternaltoXYZ(); // dihedral angle
	} else {

	}
}

void LoopModeller::sample_geometry() {
	double tole;
	double caca=3.8;
	//double can=1.5;
	//double cac=1.5;
	//double co=1.2;
	//double cn=1.3;
	std::vector<std::pair<double,double>> phipsi; // record assembles of PHIPSI
	std::vector<std::pair<double,double>> phipsigly; // only GLY could use this one
	std::vector<XYZ> fixedcrds;
	std::vector<std::pair<NSPproteinrep::BackBoneSite,NSPproteinrep::BackBoneSite>> startend;
	std::vector<std::vector<std::string>> lpseq; //attention: lpseq is not equal to loopseq_



	std::string fn;
	lowphipsi(fn,phipsi,phipsigly); //TODO : keep angle
	pset_.getval("SampToleration",tole);
	std::vector<std::vector<bool>> fixed;
	std::vector<std::vector<NSPproteinrep::FullSite>> fsss=getainstance(fixed);
	for(int i=0;i<fixed.size();i++) {
		for(int j=0;j<fixed[j].size();j++) {
			if(!fixed[i][j]) continue;
			std::map<std::string,NSPgeometry::XYZ> crds=fs.getcrds();
			for(auto &cd:crds) fixedcrds.push_back(cd.second);
		}
	}
	for(int i=0;i<fixed.size();i++) {
		for(int j=0;j<fixed[i].size();j++) {
			if(fixed[i][j]) continue;
			BackBoneSite bsst,bsen;
			if(j==0) bsst.resname="UNE"; // unexsited
			else bbst=fsss[i][j-1].getbackbonesite();
			std::vector<std::string> vs;
			int k=j;
			for(;k<fixed[i].size() && !fixed[i][k];k++) {
				vs.push_back(fsss[i][k].resname());
			}
			lpseq.push_back(vs);
			if(k==fixed.size()) bsen.resname="UNE";
			else bsen=fsss[i][k].getbackbonesite();
			startend.push_back({bsst,bsen});
			j=k;
		}
	}

	for(int i=0;i<lpseq.size();i++) {
		for(int j=0;j<lpseq[i].size();i++) {

		}
	}
}

*/






































































/*
void LoopModeller::run(bool out, std::string dir) {//default: archivedloops_ only has 1
	int sample_maxstep{-1};
	double force_index{0.0};
	int sample_record{-1};
	double rmsdchangecutoff{0.0};
	int opt_maxstep{-1};
	int opt_algorithm{-1};
	int nresults_{-1};
	int step_interval;

	pset_.getval("SampMaxStep",&sample_maxstep);
	pset_.getval("ForceIndex",&force_index);
	pset_.getval("SampRecord",&sample_record);
	pset_.getval("RMSDChangeCutOff",&rmsdchangecutoff);
	rmsdchangecutoff *= A2NM;
	pset_.getval("OptMaxStep",&opt_maxstep);
	pset_.getval("OptAlgorithm",&opt_algorithm);
	pset_.getval("NumberSample",&nresults_);
	pset_.getval("CheckInterval",&step_interval);

	double samptemp,opttemp;
	double samppp,optpp;
	double sampscconf,optscconf;
	double sampscpack,optscpack;
	double sampls,optls;
	double samppack,optpack;
	pset_.getval("SampTemperature",&samptemp);
	pset_.getval("OptTemperature",&opttemp);
	pset_.getval("SampPHIPSI",&samppp);
	pset_.getval("OptPHIPSI",&optpp);
	pset_.getval("SampLocal",&sampls);
	pset_.getval("OptLocal",&optls);
	pset_.getval("SampPack",&samppack);
	pset_.getval("OptPack",&optpack);
	pset_.getval("SampSCConf",&sampscconf);
	pset_.getval("OptSCConf",&optscconf);
	pset_.getval("SampSCPack",&sampscpack);
	pset_.getval("OptSCPack",&optscpack);

	std::string controlname;
	pset_.getval("ControlGenChain",&controlname);

	std::shared_ptr<ForceField> ff_=std::shared_ptr<ForceField>(new ForceField(genchain_->make_forcefield(controlname+"_ff")));
	ff_->changeweight(optpp,optls,optpack,optscconf,optscpack);
	genchain_->setactiveselections(ff_.get());

	std::ofstream recordstm,energystm,pdbstm;
	if(dir.back()!='/') dir+='/';
	std::string fnrecord=dir+"record.dat";
	std::string fnenergy=dir+"energy.dat";
	int numpdb=1;
	if(out) {
		recordstm.open(fnrecord);
		energystm.open(fnenergy);
	}

	//double ene11;
	while(true) {//TODO : write a function to automatic judge if reasonable
		std::cout <<"Get Initial Structure ..." <<std::endl;
		archivedloops_.confcrds().clear();
		archivedloops_.energys().clear();
		std::vector<std::vector<NSPproteinrep::FullSite>> fsss=getainstance();
		std::vector<double> initcrd = genchain_->getcrd(fsss);
		for (auto &c : initcrd) c *= A2NM;
		std::shared_ptr<EMRun> emrun=std::make_shared<EMRun>(EMRun(genchain_.get(),ff_.get()));
		double ene;
		emrun->run(opt_algorithm,opt_maxstep,&initcrd,&ene,step_interval);
		std::vector<double> enes=emrun->potenergies();
		std::set<int> removed;
		std::vector<double> rds;
		archivedloops_.addmodel(initcrd,enes,rds,removed,1);
		if(out) {
			recordstm <<"1 " <<removed.size();
			for(auto &r:removed) recordstm <<' ' <<r;
			recordstm <<std::endl;
			for(double &e:enes) energystm <<std::fixed <<std::setprecision(3) <<e <<'\t';
			energystm <<std::endl;
			pdbstm.open(dir+std::to_string(numpdb)+".pdb");
			genchain_->writepdb(initcrd,pdbstm,10.0);
			pdbstm.close();
			numpdb++;
		}
		break;
	}

	SDRun::SDRunIn sdrunin(archivedloops_.confcrds()[0],controlname);
	int rng = NSPdstl::RandomEngine<>::getinstance().intrng(1,1000)();
	std::shared_ptr<SDRun> sdrun_=std::make_shared<SDRun>(SDRun(genchain_->make_sdrun(sdrunin,rng)));
	sdrun_->changetemperature(samptemp);
	sdrun_->changeweight(samppp,sampls,samppack,sampscconf,sampscpack);

	int maxcycle;
	pset_.getval("MaxCycle",&maxcycle);
	int ix=0;
	int num=0;
	std::shared_ptr<EMRun> emrun=std::make_shared<EMRun>(EMRun(genchain_.get(),ff_.get()));
	while(archivedloops_.confcrds().size()<nresults_) {
		std::cout <<"Archived:   " <<archivedloops_.confcrds().size() <<' ' <<archivedloops_.energys().size() <<std::endl;
		std::cout <<"get decoy: " <<num++ <<std::endl;
		if(num==maxcycle) break;
		bool added{true};
		sdrun_->resetstate(archivedloops_.confcrds().at(ix));
		bool samplesuccess=sdrun_->runsteps_rmsd(sample_maxstep, archivedloops_.confcrds(),
				force_index, sample_record, rmsdchangecutoff);
		if(samplesuccess) {
			std::vector<double> crd=sdrun_->state().crd;
			double ene;
			//std::shared_ptr<EMRun> emrun=std::make_shared<EMRun>(EMRun(genchain_.get(),ff_.get()));
			emrun->run(opt_algorithm,opt_maxstep,&crd,&ene,step_interval);
			std::vector<double> enes=emrun->potenergies();
			//std::vector<double> enes=sdrun_->potenergies();
			std::set<int> removed;
			std::vector<double> rds;
			ff_->getmcrmsd(archivedloops_.confcrds(),crd,sdrun_->forceoff(),rds);
			bool added=archivedloops_.addmodel(crd,enes,rds,removed,nresults_);
			if(out) {
				if(added) recordstm <<"1 ";
				else recordstm <<"0 ";
				recordstm <<removed.size();
				for(auto &r:removed) recordstm <<' ' <<r;
				recordstm <<std::endl;
				for(double &e:enes) energystm <<std::fixed <<std::setprecision(3) <<e <<'\t';
				energystm <<std::endl;
				pdbstm.open(dir+std::to_string(numpdb)+".pdb");
				genchain_->writepdb(crd,pdbstm,10.0);
				pdbstm.close();
				numpdb++;
			}
			if(added) continue;
			ix++;
			if(ix==archivedloops_.confcrds().size()) ix=0;
		} else {
			//std::cout <<"\tshake failure! " <<num <<std::endl;
			ix++;
			if(ix==archivedloops_.confcrds().size()) ix=0;
		}
	}

	if(out) {
		recordstm.close();
		energystm.close();
		//recordstm.open(dir+"rmsd.dat");
		//std::vector<std::vector<double>> &cs=archivedloops_.confcrds();
		//for(const auto &c:cs) {
		//	std::vector<double> rds;
		//	ff_->getmcrmsd(cs,c,sdrun_->forceoff(),rds);
		//	for(double &d:rds) recordstm <<d <<' '; recordstm <<std::endl;
		//}
		//recordstm.close();
	}
}
*/






/*
bool ArchivedLoops::addmodel(std::vector<double> &model, const std::vector<double> &energy,
		int maxnum, GenChain *genchain, std::set<int> &ints) {
	//true means model is added, false means not
	//bool addnew{true};
	ints.clear();
	std::set<int> eps;
	for(int i=0;i<confcrds_.size();i++) {
		if(confcrds_[i].empty()) {
			eps.insert(i);
			continue;
		}
		if(genchain->calrmsd(confcrds_[i],model)>rmsdcut_) continue;
		if(energy[0]>energies_[i][0]) return false;
		ints.insert(i);
	}
	for(auto &i:ints) {
		confcrds_[i].clear();
		energies_[i].clear();
		eps.insert(i);
	}
	if(!eps.empty()) {
		auto it=eps.begin();
		confcrds_[*it] = model;
		energies_[*it] = energy;
		eps.erase(*it);
		return true;
	}
	if(confcrds_.size()==maxnum) {
		int ix=-1;
		double eh=-100000.0;
		for(int i=0;i<energies_.size();i++) {
			if(energies_[i][0]>eh) {
				ix=i;
				eh=energies_[i][0];
			}
		}
		if(eh>energy[0]) {
			confcrds_[ix] = model;
			energies_[ix] = energy;
			return true;
		}
		return false;
	}
	confcrds_.push_back(model);
	energies_.push_back(energy);
	return true;
}*/


//void LoopModeller::initopt() {
	/*std::vector<std::vector<NSPproteinrep::FullSite>> fsss=getainstance();
	std::vector<double> initcrd = genchain_->getcrd(fsss);
	for (auto &c : initcrd) c *= A2NM;
	std::string controlname;
	pset_.getval("ControlGenChain",&controlname);
	SDRun::SDRunIn sdrunin(initcrd,controlname);
	int rng = NSPdstl::RandomEngine<>::getinstance().intrng(1,1000)();
	SDRun sdrun=genchain_->make_sdrun(sdrunin,rng);
	ForceField ff=genchain_->make_forcefield(controlname+"_ff");
	genchain_->setactiveselections(&ff);
	EMRun emrun(genchain_.get(),&ff);*/

	/*int opt_maxstep{-1};
	int opt_algorithm{-1};

	pset_.getval("OptMaxStep",&opt_maxstep);
	pset_.getval("OptAlgorithm",&opt_algorithm);

	//ff_=std::shared_ptr<ForceField>(new ForceField(genchain_->make_forcefield(controlname+"_ff")));
	//ff_->changeweight(optpp,optls,optpack,optscconf,optscpack);
	//genchain_->setactiveselections(ff_.get());
	std::shared_ptr<EMRun> emrun__=std::make_shared<EMRun>(EMRun(genchain_.get(),ff_.get()));

	while(true) {//TODO : write a function to automatic judge if reasonable
		archivedloops_.confcrds().clear();
		archivedloops_.energys().clear();
		std::vector<std::vector<NSPproteinrep::FullSite>> fsss=getainstance();
		std::vector<double> initcrd = genchain_->getcrd(fsss);
		for (auto &c : initcrd) c *= A2NM;
		std::vector<double> gradient;
		double ene;
		emrun__->run(opt_algorithm,opt_maxstep,&initcrd,&ene,&gradient);
		std::vector<double> enes=emrun__->potenergies();
		archivedloops_.addmodel(initcrd,enes,1,genchain_.get());
		break;
	}*/
//}
/*
bool LoopModeller::getadecoy(std::vector<double> &initcrds, SDRun &sdrun, EMRun &emrun) {
	sdrun.resetstate(initcrds);
	//if(!sdrun.runsteps_rmsd(sample_maxstep, archivedloops_.confcrds(), force_index, sample_record, rmsdchangecutoff)) {
	//	std::cout <<"shake failure!" <<std::endl;
	//}
	std::vector<double> crd=sdrun.state().crd;
	std::vector<double> gradient;
	double ene;
	//emrun.run(opt_algorithm,opt_maxstep,&crd,&ene,&gradient);
	emrun.run(opt_algorithm,opt_maxstep,&initcrds,&ene,&gradient);
	std::vector<double> enes=emrun.potenergies();
	if(archivedloops_.addmodel(crd,enes,nresults_,genchain_.get())) return true;
	return false;
}

void LoopModeller::testarchived(std::string fncrd, std::string fnene) {
	std::vector<std::vector<NSPproteinrep::FullSite>> refpdb=NSPproteinrep::readfullsitesfrompdb(fncrd);
	std::vector<double> initcrd = genchain_->getcrd(refpdb);
	for (auto &c : initcrd) c *= A2NM;
	std::ifstream ifs(fnene);
	std::string line;
	std::vector<double> enes;
	while(std::getline(ifs,line)) enes.push_back(std::stod(line));
	ifs.close();
	archivedloops_.addmodel(initcrd,enes,100,genchain_.get());
}*/






