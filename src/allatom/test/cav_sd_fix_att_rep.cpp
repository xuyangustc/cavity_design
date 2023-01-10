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
/*
 * fixed_points
 * atoms_seq
 * dis, f_
 */
std::vector<ForceField::Fixedpoint_Attraction_Repulsion> readfixedpointinfo(
		std::string parfile, ForceField::Site_Restrict_RMSD &sirerm,
		const std::map<std::vector<int>,std::map<std::string,int>>&atomseq) {
	std::vector<std::string> ls;
	std::ifstream ifs(parfile);
	std::string readline;
	while(std::getline(ifs,readline)) {
		if(readline.empty()) continue;
		if(readline[0]=='#') continue;
		ls.push_back(readline);
	}
	ifs.close();
	std::vector<ForceField::Fixedpoint_Attraction_Repulsion> pointinfo;
	int nfar=std::stoi(ls[0]);
	int i=1;
	//std::cout <<"FAR size: " <<nfar <<std::endl;
	//for(int i=0;i<ls.size();) {
	for(int nf=0;nf<nfar;nf++) {
		int nline=std::stoi(ls[i++]);
		ForceField::Fixedpoint_Attraction_Repulsion far;
		far.fixed_points.resize(nline);
		for(int j=0;j<nline;j++) {
			std::stringstream ss(ls[i++]);
			ss >>far.fixed_points[j].x_ >>far.fixed_points[j].y_ >>far.fixed_points[j].z_;
			far.fixed_points[j] = far.fixed_points[j]*0.1;
		}
		nline=std::stoi(ls[i++]);
		for(int j=0;j<nline;j++) {
			std::stringstream ss(ls[i++]);
			int cid;
			std::string rid;
			std::string lb;
			ss >>cid >>rid >>lb;
			bool found{false};
			for(char &c:rid) {
				if(c=='-') {
					c=' ';
					found=true;
					break;
				}
			}
			int r1, r2;
			if(found) {
				std::stringstream ss1(rid);
				ss1 >>r1 >>r2;
			} else {
				r1=r2=std::stoi(rid);
			}
			for(int k=r1;k<=r2;k++) {
				std::vector<int> key{cid,k};
				if(lb=="CA") {
					far.atoms_seq.push_back(atomseq.at(key).at("CA"));
				} else if(lb=="MC") {
					far.atoms_seq.push_back(atomseq.at(key).at("N"));
					far.atoms_seq.push_back(atomseq.at(key).at("CA"));
					far.atoms_seq.push_back(atomseq.at(key).at("C"));
					far.atoms_seq.push_back(atomseq.at(key).at("O"));
				} else if(lb=="AA") {
					for(const auto &p:atomseq.at(key)) {
						far.atoms_seq.push_back(p.second);
					}
				}
			}
		}
		//nline=std::stoi(ls[i++]);
		far.dis.resize(6);
		std::stringstream ss(ls[i++]);
		for(int j=0;j<6;j++) {
			ss >>far.dis[j];
			far.dis[j] *=0.1;
		}
		ss >>far.f_rep >>far.f_att;
		pointinfo.push_back(far);
		//std::cout <<nf <<std::endl;
		//std::cout <<'\t' <<pointinfo.back().fixed_points.size() <<std::endl;
		//std::cout <<'\t' <<pointinfo.back().atoms_seq.size() <<std::endl;
		//std::cout <<'\t' <<pointinfo.back().f_rep <<'\t' <<pointinfo.back().f_att <<std::endl;
		//std::cout <<'\t' <<pointinfo.back().dis[0] <<'\t' <<pointinfo.back().dis[1] <<'\t'
		//		<<pointinfo.back().dis[2] <<'\t' <<pointinfo.back().dis[3] <<'\t'
		//		<<pointinfo.back().dis[4] <<'\t' <<pointinfo.back().dis[5] <<std::endl;
	}
	if(i==ls.size()) return pointinfo;
	std::stringstream ssd(ls[i++]);
	ssd >>sirerm.dis_near >>sirerm.dis_far;
	sirerm.dis_near *= 0.1;
	sirerm.dis_far *= 0.1;
	std::stringstream ssf(ls[i++]);
	ssf >>sirerm.f_near >>sirerm.f_far;
	int nline =std::stoi(ls[i++]);
	for(int nl=0;nl<nline;nl++) {
		std::vector<XYZ> cs(4);
		for(int j=0;j<4;j++) {
			std::stringstream ss(ls[i++]);
			ss >>cs[j].x_ >>cs[j].y_ >>cs[j].z_;
			cs[j] = cs[j]*0.1;
		}
		sirerm.fixed.push_back(cs);
	}
	for(const auto &p:atomseq) {
		std::vector<int> v;
		v.push_back(p.second.at("N"));
		v.push_back(p.second.at("CA"));
		v.push_back(p.second.at("C"));
		v.push_back(p.second.at("O"));
		sirerm.nbb.push_back(v);
	}
	//std::cout <<"RMSD Restrict:" <<std::endl;
	//std::cout <<'\t' <<sirerm.fixed.size() <<std::endl;
	//std::cout <<'\t' <<sirerm.fixed[1][1].x_ <<'\t' <<sirerm.fixed[1][1].y_ <<'\t'
	//		<<sirerm.fixed[1][1].z_ <<std::endl;
	//std::cout <<'\t' <<sirerm.nbb.size() <<std::endl;
	//std::cout <<'\t' <<sirerm.f_near <<'\t' <<sirerm.f_far <<'\t'
	//		<<sirerm.dis_near <<'\t' <<sirerm.dis_far <<std::endl;
	//exit(1);
	return pointinfo;
}
int main(int argc, char **argv) {

	if(argc==1) {
		std::cout <<"Usage:" <<std::endl;
		std::cout <<"\t${parfile1(SD)}, ${parfile2}" <<std::endl;
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

	std::map<std::vector<int>,std::map<std::string,int>> atomseq=genchain.atomseq();
	std::string parfile2(argv[2]);
	ForceField::Site_Restrict_RMSD sirerm;
	std::vector<ForceField::Fixedpoint_Attraction_Repulsion> pointinfo
			= readfixedpointinfo(parfile2,sirerm,atomseq);

	/*std::cout <<pointinfo.size() <<std::endl;
	for(int i=0;i<pointinfo.size();i++) {
		std::cout <<i;
		for(int j:pointinfo[i].first.first) std::cout <<' ' <<j;
		std::cout <<std::endl;
		std::cout <<pointinfo[i].first.second.x_ <<' ' <<pointinfo[i].first.second.y_ <<' '
				<<pointinfo[i].first.second.z_ <<' ' <<pointinfo[i].second[0] <<' '
				<<pointinfo[i].second[1] <<' ' <<pointinfo[i].second[2] <<std::endl;
	}*/

	/*bool traj{false};
	int nit;
	if(argc==6) {
		traj=true;
		nit=std::stoi(std::string(argv[5]));
	}*/
	for(int i=0;i<nsteps/100;++i) {
		if(!(sdrun.runsteps_fixedpoint_attraction_repulsion(100,pointinfo,sirerm,callback))){
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



