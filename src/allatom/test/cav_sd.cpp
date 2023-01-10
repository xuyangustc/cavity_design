/*
 * testchainsd.cpp
 *
 *  Created on: 2018骞�1鏈�25鏃�
 *      Author: hyliu
 */

#include "sd/genchain.h"
#include "sd/sidechainff.h"
#include "fullsite/fullsite.h"
#include "fullsite/structplus.h"
#include "sd/sdrun.h"
#include "sd/trajio.h"
#include "geometry/quatfit.h"
#include "sd/trajio.h"
#include "dstl/randomengine.h"
#include "sd/nnterm.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;
using namespace NSPsd;


struct TemperatureAnnealing{
	double t_high{0.0};
	double t_low{0.0};
	int nsteps_cycle{0},nsteps_high{0}, nsteps_drop{0};
	TemperatureAnnealing(){;}
	TemperatureAnnealing(double th,double tl, double cycle,double sh,double sd):
	t_high(th),t_low(tl),nsteps_cycle(cycle),nsteps_high(sh),nsteps_drop(sd){
		assert(sh+sd <= cycle);
	}
	double temperature(int step) const {
		step=step%nsteps_cycle;
		if(step<=nsteps_high) return t_high;
		else if(step>=(nsteps_high+nsteps_drop)) return t_low;
		else return t_high - (t_high-t_low)*(double)(step-nsteps_high)/(double)(nsteps_drop);
	}
};

/**Called back in SDRun, this functor monitor quantities during SD.
 * It can also be used to effect conditional earlier ending of the
 * calling SDRun (by returning true)
 *
 */
struct SDStepCallBack {
	/**
	 * called in SDRun at every step, after the force calculation
	 */
	bool operator()(SDRun &sdrun) const {
		if (sdrun.nstepsrun() == 0) {
			//std::cout<<"Pot energies at starting point: "<<std::endl;
			for(auto e:sdrun.potenergies()) std::cout<< " "<<e;
			std::cout<<std::endl;
		}
		//monitor some quantities
		if((sdrun.nstepsrun()+1) %100==0){
			double rg=NSPgeometry::radiusgyr(sdrun.state().crd)/A2NM;
			std::cout<<"Radius of gyration: " <<rg;
			std::vector<double> w(refcrd_.size()/3,0.0);
			std::vector<int> mcatoms=sdrun.ff()->mainchainatoms();
			for(auto i:mcatoms) w[i]=1.0;
			NSPgeometry::QuatFit qfmc,qfall;
			double rmsd2mc=qfmc.setup(sdrun.state().crd,refcrd_,w);
			double rmsd2=qfall.setup(sdrun.state().crd,refcrd_);
			std::cout <<" RMSD from ref, main chain " << sqrt(rmsd2mc)/A2NM
					<<" all atom " <<sqrt(rmsd2)/A2NM;
			std::cout<<std::endl;
		}
		if(trajstep_>0 &&(sdrun.nstepsrun()+1)%trajstep_==0){
			trajio_.writeframe(sdrun.state().crd);
			trajio_.writeframe(sdrun.potenergies());
		}
		if(doannealing_){
			sdrun.changetemperature(ta_.temperature(sdrun.nstepsrun()+1),tagroup_);
		}
		//no earlier ending of the calling SDRun
		return false;
	}
	void setrefcrd(const std::vector<double> &ref){
		refcrd_=ref;
	}
	void setannealing(const std::string &controlname){
		auto &pset=SDControls::getparameterset(controlname);
		std::vector<double> annealingscheme;
		pset.getval("AnnealingScheme",&annealingscheme);
		if(annealingscheme.empty()) return;
		doannealing_=true;
		int widx=0;
		tagroup_=(int) (annealingscheme[widx++]);
		double th=annealingscheme[widx++];
		double tl=annealingscheme[widx++];
		int nstep_cycle= (int) annealingscheme[widx++];
		int step_h=(int) annealingscheme[widx++];
		int step_d=(int) annealingscheme[widx++];
		ta_=TemperatureAnnealing(th,tl,nstep_cycle,step_h,step_d);
	}
	void createtrajfile(const std::string &controlname){
		auto &pset=SDControls::getparameterset(controlname);
		std::string trajfile;
		pset.getval("TrajFile",&trajfile);
		pset.getval("NTrajSteps", &trajstep_);
		if(!trajfile.empty() && trajstep_>0) trajio_.createwriter(trajfile);
		else trajstep_=-1;
	}
	SDStepCallBack(const std::string & controlname){
		setannealing(controlname);
		createtrajfile(controlname);
	}
	std::vector<double> refcrd_;
	bool doannealing_{false};
	TemperatureAnnealing ta_;
	int tagroup_{0};
	TrajIO trajio_;
	int trajstep_{-1};
};
void writecrd(const std::vector<double> &crd,std::ostream & os){
	os<<crd.size();
	for(auto &c:crd) os<<" "<<c;
}
bool readcrd(std::istream &is,std::vector<double> &crd){
	int nrecord;
	is>>nrecord;
	crd.resize(nrecord);
	for(int i=0;i<nrecord;++i) is >>crd[i];
	return(!is.fail());
}
int main(int argc, char **argv) {

	if(argc==1) {
		std::cout <<"Usage:" <<std::endl;
		std::cout <<"\t${parfile}" <<std::endl;
		return 1;
	}

	std::string parfile(argv[1]);
	//std::string pdbfile(argv[2]);

	std::string controlname="sdffcontrol";
	NSPdataio::ControlFile cf;
	cf.readfile(parfile);
	std::vector<std::string> sdcontrolines=cf.getcontrolines("SD");
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	std::vector<std::string> genchaincontrollines=cf.getcontrolines("GenChain");
	/*for(int i=0;i<genchaincontrollines.size();i++) {
		if(genchaincontrollines[i].substr(0,10)=="RefPDBFile") {
			genchaincontrollines[i] = "RefPDBFile\t=\t" + pdbfile;
			break;
		}
	}
	std::set<char> labch{'\t',' ',',',';',':'};
	bool startcheck{false};
	bool normsd{true};
	for(int i=0;i<ffcontrolines.size();i++) {
		if(ffcontrolines[i].substr(0,17)=="TotalRMSDRestrain") {
			for(int j=0;j<ffcontrolines[i].size();j++) {
				if(ffcontrolines[i][j]=='=') {
					startcheck=true;
					continue;
				}
				if(!startcheck) continue;
				if(labch.find(ffcontrolines[i][j])==labch.end()) {
					normsd=false;
					break;
				}
			}
			break;
		}
	}
	if(!normsd) {
		for(int i=0;i<ffcontrolines.size();i++) {
			if(ffcontrolines[i].substr(0,17)=="TotalRMSDRestrain") {
				//ffcontrolines[i] = "RefPDBFile\t=\t" + pdbfile;
				std::string line=ffcontrolines[i];
				std::stringstream ss(line);
				std::vector<std::string> vs(5);
				for(int j=0;j<vs.size();j++) ss >> vs[j];
				ffcontrolines[i] = vs[0] + "\t=\t" +pdbfile +"\t" + vs[2] +"\t" + vs[3] +"\t" + vs[4];
				break;
			}
		}
	}*/
	definesdcontrol(controlname+"_sd",sdcontrolines);
	defineforcefieldcontrol(controlname+"_ff",ffcontrolines);
	definegenchaincontrol(controlname+"_genchain",genchaincontrollines);

	//std::vector<std::string> temp_str;
	//NSPsd::Atom3NNTerm::gettheestimator(temp_str,temp_str);

	//std::ostream & ocont=std::cout;
	//genchainprintcontrols(controlname,ocont);

	std::string outf;
	SDControls::getparameterset(controlname+"_genchain").getval("OutPutPDBFile",&outf);
	std::string trajf;
	SDControls::getparameterset(controlname+"_sd").getval("TrajFile",&trajf);
	if(outf.empty() && trajf.empty()) {
		std::cout <<"No Trajectory & No OutPut!!!" <<std::endl;
		exit(1);
	}

	int seed;
	SDControls::getparameterset(controlname+"_sd").getval("RandomSeed",&seed);

	NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	GenChain genchain(controlname+"_genchain");
	std::vector<std::vector<FullSite>> fschains = *(genchain.getrefpdb());
	std::shared_ptr<std::vector<double>> initcrd;

	initcrd = std::shared_ptr<
				std::vector<double>>(new std::vector<double>(genchain.getcrd(fschains)));
	for (auto &c : *initcrd)
		c *= A2NM;

	SDRun::SDRunIn sdrunin(*initcrd,controlname);
	std::vector<std::set<int>> notcis(fschains.size());
	SDRun sdrun = genchain.make_sdrun(sdrunin,seed,notcis);
	int nsteps;
	SDControls::getparameterset(controlname+"_sd").getval("NumberStep",&nsteps);
	SDStepCallBack callback(controlname+"_sd");
	callback.setrefcrd(*initcrd);

	/*bool traj{false};
	int nit;
	if(argc==6) {
		traj=true;
		nit=std::stoi(std::string(argv[5]));
	}*/
	for(int i=0;i<nsteps/100;++i) {
		if(!(sdrun.runsteps(100,callback))){
			std::cout<<"Shake failure occurred."<<std::endl;
			exit(1);
		}
		double temp=sdrun.temperature();
		std::cout <<"run " << i <<": ";
		std::cout <<temp<<"\t";
		for(auto e:sdrun.potenergies()) std::cout << " "<<e;
		std::cout <<std::endl;
		/*std::ofstream ofs;
		std::string filename=outf;
		if(traj) {
			int nst=(i+1)*100;
			if(nst%nit!=0) continue;
			int j=nst/nit;
			filename = outf + "_" + std::to_string(j);
			ofs.open(filename.c_str());
			genchain.writepdb(sdrun.state().crd,ofs,1.0/A2NM);
			ofs.close();
		}*/
		/*else {
			std::ofstream ofs;
			ofs.open(filename.c_str());
			genchain.writepdb(sdrun.state().crd,ofs,1.0/A2NM);
			ofs.close();
		}*/
	}

	if(!outf.empty()) {
		std::ofstream ofs;
		ofs.open(outf.c_str());
		genchain.writepdb(sdrun.state().crd,ofs,1.0/A2NM);
		ofs.close();
	}
}



