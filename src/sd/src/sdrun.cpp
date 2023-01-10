/*
 * sdrun.cpp
 *
 *  Created on: 2017骞�11鏈�9鏃�
 *      Author: hyliu
 */

#include "sd/sdrun.h"
using namespace NSPsd;
/*
SDRun NSPsd::make_backbone_sdrun(const SDRun::SDRunIn &in,unsigned int seed) {
	int nsite = in.crd->size()/12;
	auto ff = std::shared_ptr < ForceField > (new ForceField);
	auto &pset=SDControls::getparameterset(in.sdcontrolname+"_sd");
	std::vector<int> nsites;
	nsites.push_back(nsite);
	std::vector<std::vector<int>> cissites;
	cissites.push_back(in.cissites);
	*ff = make_forcefield_backbone(nsites,cissites);
	ff->usecontrols(in.sdcontrolname+"_ff");
	auto shakebds = std::shared_ptr < ShakeBonds > (new ShakeBonds);
	int doshake;
	pset.getval("DoShake",&doshake);
	*shakebds = make_shakebonds(*ff);
	shakebds->seton(doshake!=0);
	auto sd = std::shared_ptr < StochasticDynamics > (new StochasticDynamics);
	*sd = make_sd_backbone(nsite,in.sdcontrolname+"_sd");
	SDRun sdrun(sd, ff, shakebds);
	sdrun.shakeon()=doshake!=0;
	int nblsteps;
	pset.getval("NeighborListSteps",&nblsteps);
	sdrun.nblsteps()=nblsteps;
	sdrun.initrandomengine(seed);
	sdrun.initstate(in);
	return sdrun;
}*/
void NSPsd::sdreadcontrols(const std::string &filename,std::string name){
	NSPdataio::ControlFile cf;
	cf.readfile(filename);
	std::vector<std::string> sdcontrolines=cf.getcontrolines("SD");
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	definesdcontrol(name+"_sd",sdcontrolines);
	defineforcefieldcontrol(name+"_ff",ffcontrolines);
}
void NSPsd::sdprintcontrols(std::string name,std::ostream &ofs){
	ofs <<"START SD"<<std::endl;
	SDControls::getparameterset(name+"_sd").printparameters(ofs);
	ofs<<"END SD"<<std::endl;
	ofs<<"START ForceField"<<std::endl;
	ForceFieldControls::getparameterset(name+"_ff").printparameters(ofs);
	ofs<<"END ForceField"<<std::endl;
}



int SDRun::runsteps_minE(int maxstep, int checkstep) {
	std::shared_ptr<NeighborList> nbl;
	std::shared_ptr<GroupNeighbor> gn;
	//std::cout <<0 <<std::endl;
	if(ff_->stericwght()>0.000001)
		nbl = std::shared_ptr < NeighborList > (new NeighborList(state_->crd, *ff_, &forceoff_));
	//std::cout <<"sd: " <<acts_->ff()->bsinchains().size() <<std::endl;
	//std::cout <<" " <<acts_->ff()->bsinchains()[0].size() <<std::endl;
	std::vector<double> forces = ff_->forces(state_->crd, *nbl, *gn, &potenergies_, *acts_);
	//std::cout <<2 <<std::endl;
	auto state_min = state_;
	std::vector<double> es_min=potenergies_;
	int nstep_min = 0;
	//std::cout <<"start\t" <<es_min[0] <<std::endl;
	//std::cout <<3 <<std::endl;
	for (int i = 0; i < maxstep; ++i) {
		//std::cout <<4 <<std::endl;
		if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces, shakebds_->shakefunction(), &vmasses_)) break;
		//std::cout <<5 <<std::endl;
		auto temp = state_;
		state_ = buffstate_;
		buffstate_ = temp;
		++nstepsrun_;
		if ((i+1) % nblsteps_ == 0) {
			if(ff_->stericwght()>0.000001)
				nbl = std::shared_ptr < NeighborList > (new NeighborList(state_->crd, *ff_, &forceoff_));
		}
		forces = ff_->forces(state_->crd, *nbl, *gn, &potenergies_, *acts_);
		if(potenergies_[0]<es_min[0]) {
			state_min = state_;
			es_min=potenergies_;
			nstep_min=i+1;
		}
		//std::cout <<i <<"\t" <<potenergies_[0] <<"\t" <<es_min[0] <<std::endl;
		if(i+1-nstep_min>=checkstep) {
			break;
		}
	}
	//exit(1);
	state_ = state_min;
	potenergies_ = es_min;
	return nstep_min;
}








/*
bool SDRun::smallchanged(std::vector<double> &records, int step_interval, double ene_diff) {
	int prev = records.size() - step_interval;
	if(prev < 0) return false;
	if(fabs(records.back()-records[prev])>ene_diff) return false;
	return true;
}

bool SDRun::run_sample_maxrmsd(int nsteps, const std::vector<std::vector<double>> &confs,
		double index, int step_interval, double rmsd_diff, bool mconly, std::vector<double> &rdrd,
		std::vector<std::vector<double>> &enes) {
	std::shared_ptr<NeighborList> nbl;
	int max_step{0};
	std::vector<double> max_ene;
	StochasticDynamics::State max_state;
	std::vector<double> records;
	bool get_max{false};

	for (int i = 0; i < nsteps; ++i) {
		if (i % nblsteps_ == 0)
			nbl = std::shared_ptr < NeighborList > (new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_->forces(state_->crd, *nbl, &potenergies_, *acts_);
		enes.push_back(potenergies_);
		double minrmsd=ff_->forces_rmsd(state_->crd, forces, forceoff_, confs, index, mconly);
		rdrd.push_back(minrmsd);
		//if(i%20==0) std::cout <<"Step = " <<i <<"\tMinRMSD = " <<minrmsd <<"\tEnergy = " <<potenergies_[0] <<std::endl;
		if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces, shakebds_->shakefunction(), &vmasses_)) {
			std::cout <<"Sampling:  Shake Fail!" <<std::endl;
			std::cout <<"\tNowStep = " <<i <<"\tMaxStep = " <<max_step <<std::endl;
			if(!records.empty()) {
				potenergies_ = max_ene;
				*state_ = max_state;
			}
			return false;
		}
		if(!get_max) {
			auto temp = state_;
			state_ = buffstate_;
			buffstate_ = temp;
			++nstepsrun_;
			if(records.empty() || minrmsd > records.back()) {
				records.push_back(minrmsd);
				max_step = i;
				max_ene = potenergies_;
				max_state = *state_;
			} else {
				records.push_back(records.back());
			}
			if(smallchanged(records,step_interval,rmsd_diff)) {
				std::cout <<"Sampling:  Meet Max RMSD!" <<std::endl;
				get_max = true;
			}
		}
		if(i+1==nsteps) {
			get_max = true;
			std::cout <<"Sampling:  Max Step!" <<std::endl;
		}
		if(get_max) {
			if(!records.empty()) {
				potenergies_ = max_ene;
				*state_ = max_state;
			}
			std::cout <<"\tNowStep = " <<i <<"\tMaxStep = " <<max_step <<std::endl;
			break;
		}
	}
	return true;
}

std::vector<double> SDRun::confene(const std::vector<double> &crds) {
	std::shared_ptr<NeighborList> nbl = std::shared_ptr < NeighborList > (new NeighborList(crds, *ff_, &forceoff_));
	ff_->forces(crds, *nbl, &potenergies_, *acts_);
	return potenergies_;
}

bool SDRun::run_opt(int nsteps, int step_interval, double ene_diff) {
	std::shared_ptr<NeighborList> nbl;
	int min_step{0};
	std::vector<double> min_ene;
	StochasticDynamics::State min_state;
	std::vector<double> records;
	bool get_min{false};

	for (int i = 0; i < nsteps; ++i) {
		if (i % nblsteps_ == 0)
			nbl = std::shared_ptr < NeighborList > (new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_->forces(state_->crd, *nbl, &potenergies_, *acts_);
		if(!sd_->leapstep(*state_, *buffstate_, *rng_, forces, shakebds_->shakefunction(), &vmasses_)) {
			std::cout <<"Optimization:  Shake Fail!" <<std::endl;
			std::cout <<"\tNowStep = " <<i <<"\tMinStep = " <<min_step <<std::endl;
			if(!records.empty()) {
				potenergies_ = min_ene;
				*state_ = min_state;
			}
			return false;
		}
		if(!get_min) {
			auto temp = state_;
			state_ = buffstate_;
			buffstate_ = temp;
			++nstepsrun_;
			if(records.empty() || potenergies_[0] < records.back()) {
				records.push_back(potenergies_[0]);
				min_step = i;
				min_ene = potenergies_;
				min_state = *state_;
			} else {
				records.push_back(records.back());
			}
			if(smallchanged(records,step_interval,ene_diff)) {
				std::cout <<"Optimization:  Meet Energy Minimum!" <<std::endl;
				get_min = true;
			}
		}
		if(i+1==nsteps) {
			get_min = true;
			std::cout <<"Optimization:  Max Step!" <<std::endl;
		}
		if(get_min) {
			if(!records.empty()) {
				potenergies_ = min_ene;
				*state_ = min_state;
			}
			std::cout <<"\tNowStep = " <<i <<"\tMinStep = " <<min_step <<std::endl;
			break;
		}
	}
	return true;
}

bool SDRun::resetstate(const std::vector<double> & initcrds) {
	state_ = sd_->make_initstate(initcrds, *rng_, &vmasses_);
	*buffstate_ = *state_;
	nstepsrun_ = 0;
	return sd_->shakestate(*state_, shakebds_->shakefunction(),true, &vmasses_);
}

bool SDRun::run_sample_rmsdcut(int nsteps, const std::vector<std::vector<double>> &confs,
		double index, double rmsdcut, bool mconly, std::vector<double> &rdrd,
		std::vector<std::vector<double>> &enes) {
	std::shared_ptr<NeighborList> nbl;
	for (int i = 0; i < nsteps; ++i) {
		if (i % nblsteps_ == 0)
			nbl = std::shared_ptr < NeighborList > (new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_->forces(state_->crd, *nbl, &potenergies_, *acts_);
		//enes.push_back(potenergies_);
		double minrmsd=ff_->forces_rmsd(state_->crd, forces, forceoff_, confs, index, rmsdcut, mconly);
		//rdrd.push_back(minrmsd);
		if(minrmsd>rmsdcut) {
			std::cout <<"Sampling:  Success!" <<std::endl;
			std::cout <<"\tNowStep = " <<i <<"\tMinRMSD: " <<minrmsd*10.0 <<std::endl;
			return true;
		}
		//return true;
		//if(i%20==0) std::cout <<"Step = " <<i <<"\tMinRMSD = " <<minrmsd <<"\tEnergy = " <<potenergies_[0] <<std::endl;
		if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces, shakebds_->shakefunction(), &vmasses_)) {
			std::cout <<"Sampling:  Shake Fail!" <<std::endl;
			std::cout <<"\tNowStep = " <<i <<std::endl;
			return false;
		}
		auto temp = state_;
		state_ = buffstate_;
		buffstate_ = temp;
		++nstepsrun_;
		if(i+1==nsteps) {
			std::cout <<"Sampling:  Max Step!" <<std::endl;
			std::cout <<"\tNowStep = " <<i <<"\tMinRMSD: " <<minrmsd*10.0 <<std::endl;
			return false;
		}
	}
}

bool SDRun::run_sample_rmsdcut_ene(const std::vector<std::vector<double>> &confs,
		double lowrmsdcut, double highrmsdcut, double enecut, double index, bool mconly, int nsteps,
		std::vector<std::pair<int,double>> &minrmsdrecord, std::vector<std::vector<double>> &enes) {
	std::shared_ptr<NeighborList> nbl;
	for (int i = 0; i < nsteps; ++i) {
		if (i % nblsteps_ == 0)
			nbl = std::shared_ptr < NeighborList > (new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_->forces(state_->crd, *nbl, &potenergies_, *acts_);
		enes.push_back(potenergies_);
		double staene = potenergies_[ForceField::ETOT] - potenergies_[ForceField::EBOND] -
				potenergies_[ForceField::EANG] - potenergies_[ForceField::EIMPDIH] -
				potenergies_[ForceField::ESTERIC] - potenergies_[ForceField::ESCPACKING];
		double rmsdcut=highrmsdcut;
		if(staene<enecut) rmsdcut=lowrmsdcut;
		int minstep;
		double minrmsd=ff_->forces_rmsd_ene(state_->crd, forces, forceoff_, confs, index, mconly, rmsdcut, minstep);
		minrmsdrecord.push_back({minstep,minrmsd});
		if(minrmsd>rmsdcut) {
			std::cout <<"Sampling:  Success!" <<std::endl;
			std::cout <<"\tNowStep = " <<i <<"\tMinRMSD: " <<minrmsd <<std::endl;
			return true;
		}
		//if(i%20==0) std::cout <<"Step = " <<i <<"\tMinRMSD = " <<minrmsd <<"\tEnergy = " <<potenergies_[0] <<std::endl;
		if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces, shakebds_->shakefunction(), &vmasses_)) {
			std::cout <<"Sampling:  Shake Fail!" <<std::endl;
			std::cout <<"\tNowStep = " <<i <<std::endl;
			return false;
		}
		auto temp = state_;
		state_ = buffstate_;
		buffstate_ = temp;
		++nstepsrun_;
		if(i+1==nsteps) {
			std::cout <<"Sampling:  Max Step!" <<std::endl;
			std::cout <<"\tNowStep = " <<i <<"\tMinRMSD: " <<minrmsd <<std::endl;
			return false;
		}
	}
}

bool SDRun::runsteps_beta_sheet(int nsteps, const std::vector<std::pair<int,int>> &hbpairs,
		double f_index_hb, const std::vector<std::pair<int,int>> &capairs, double f_index_ca) {
	std::shared_ptr<NeighborList> nbl;
	for (int i = 0; i < nsteps; ++i) {
		if (i % nblsteps_ == 0)
			nbl = std::shared_ptr < NeighborList
					> (new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_->forces_beta_sheet(state_->crd, *nbl,
				&potenergies_, *acts_, hbpairs, f_index_hb, capairs, f_index_ca);
		if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces,
				shakebds_->shakefunction(), &vmasses_))
			return false;
		auto temp = state_;
		state_ = buffstate_;
		buffstate_ = temp;
		++nstepsrun_;
	}
	return true;
}

bool SDRun::runsteps_strand_unbend(int nsteps, double f_index, const std::vector<std::pair<int,int>> &capairs) {
	std::shared_ptr<NeighborList> nbl;
	for (int i = 0; i < nsteps; ++i) {
		if (i % nblsteps_ == 0)
			nbl = std::shared_ptr < NeighborList
					> (new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_->forces_strand_unbend(state_->crd, *nbl,
				&potenergies_, *acts_, f_index, capairs);
		if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces,
				shakebds_->shakefunction(), &vmasses_))
			return false;
		auto temp = state_;
		state_ = buffstate_;
		buffstate_ = temp;
		++nstepsrun_;
	}
	return true;
}

bool SDRun::runsteps_strand_hb(int nsteps, const std::vector<std::pair<int,int>> &hbpairs,
		double discutoff, double f_index) {
	std::shared_ptr<NeighborList> nbl;
	for (int i = 0; i < nsteps; ++i) {
		if (i % nblsteps_ == 0)
			nbl = std::shared_ptr < NeighborList
					> (new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_->forces_strand_hb(state_->crd, *nbl,
				&potenergies_, *acts_, hbpairs, discutoff, f_index);
		if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces,
				shakebds_->shakefunction(), &vmasses_))
			return false;
		auto temp = state_;
		state_ = buffstate_;
		buffstate_ = temp;
		++nstepsrun_;
	}
	return true;
}

bool SDRun::runsteps_subentry(int nsteps,
		const std::vector<std::pair<std::vector<int>,std::string>> &atomseq,
		EnergyTerm &eneterm, const std::vector<int> &sizes) {
	std::shared_ptr<NeighborList> nbl;
	for (int i = 0; i < nsteps; ++i) {
		if (i % nblsteps_ == 0)
			nbl = std::shared_ptr < NeighborList
					> (new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_ -> forces_subentry(state_->crd, *nbl, &potenergies_, *acts_,
					atomseq, eneterm, sizes);
		if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces,
				shakebds_->shakefunction(), &vmasses_))
			return false;
		auto temp = state_;
		state_ = buffstate_;
		buffstate_ = temp;
		++nstepsrun_;
	}
	return true;
}
*/
void SDRun::energyfordecoy(int nsteps,
		const std::vector<std::pair<std::vector<int>,std::string>> &atomseq,
		EnergyTerm &eneterm, const std::vector<int> &sizes) {
	std::shared_ptr<NeighborList> nbl = std::shared_ptr < NeighborList
						> (new NeighborList(state_->crd, *ff_, &forceoff_));
	ff_ -> forces_subentry(state_->crd, *nbl, &potenergies_, *acts_, atomseq, eneterm, sizes);
}
void SDRun::printenergy() {
	std::shared_ptr<NeighborList> nbl;
	std::shared_ptr<GroupNeighbor> gn;
	nbl = std::shared_ptr < NeighborList > (new NeighborList(state_->crd, *ff_, &forceoff_));
	ff_->forces(state_->crd, *nbl, *gn, &potenergies_, *acts_);
}

bool SDRun::runsteps_subentry(int nsteps,
		const std::vector<std::pair<std::vector<int>,std::string>> &atomseq,
		EnergyTerm &eneterm, const std::vector<int> &sizes) {
	return true;
}

SDRun SDRunAssemble::FindSDRun(SDRunAssembleControl &sdrac) {
	if(sdruns.find(sdrac)!=sdruns.end()) {
		//std::cout <<"in" <<std::endl;
		//return sdruns.at(sdrac);
	}
	//std::cout <<"out" <<std::endl;

	std::vector<std::vector<std::map<std::string,int>>> serial;
	//ff
	std::string ffcontrolname = sdrac.controlname+"_ff";
	auto ff = std::shared_ptr < ForceField > (new ForceField);
	auto cissite_temp = sdrac.cissites;
	*ff = make_forcefield_allatom(sdrac.sequence, ffcontrolname, cissite_temp, serial);
	ff->setanalysismode(true);
	//sdrun.ff() = std::shared_ptr < ForceField > (new ForceField);
	//*(sdrun.ff()) = ff_temp;
	auto &pset = SDControls::getparameterset(sdrac.controlname + "_sd");
	//shakebd
	auto shakebds = std::shared_ptr < ShakeBonds > (new ShakeBonds);
	*shakebds = make_shakebonds(*ff);
	//sd
	std::vector<std::string> tgroups;
	pset.getval("TemperatureGroups",&tgroups);
	auto agrps=std::shared_ptr<std::vector<std::vector<int>>>(
			new std::vector<std::vector<int>>(tgroups.size()));
	int na=0;
	//lack: ssregions, coilregions
	for(int i=0;i<tgroups.size();++i){
		if(tgroups[i]=="all") {
			for(int a=0;a<ff->natoms();++a) (*agrps)[i].push_back(a);
		} else	if(tgroups[i]=="mainchain"){
			(*agrps)[i]=ff->mainchainatoms();
		} else if(tgroups[i] =="sidechain"){
			(*agrps)[i]=ff->sidechainatoms();
		} else {
			std::cout<<"Undefined TemperatureGroup Type: " <<tgroups[i]<<std::endl;
			exit(1);
		}
		na +=(*agrps)[i].size();
	}
	assert (na==ff->natoms());
	std::vector<double> temperatures;
	pset.getval("Temperatures",&temperatures);
	for(auto &t:temperatures) t*=KBT;
	std::vector<double> masses;
	std::vector<double> gammas;
	double gamma;
	pset.getval("FrictionCoeff",&gamma);
	for(int c=0;c<sdrac.sequence.size();++c){
		for (int i = 0; i < sdrac.sequence[c].size(); ++i) {
			int natoms=4+VSCType::getVSCType(sdrac.sequence[c][i]).nscatoms;
			for(int a=0;a<natoms;++a)
				masses.push_back(14.0);
			for (int d = 0; d < natoms; ++d)
				gammas.push_back(gamma);
		}
	}
	double timestep;
	pset.getval("TimeStep",&timestep);
	auto sd = std::shared_ptr < StochasticDynamics > (new StochasticDynamics);
	sd->init(masses, gammas, timestep, *agrps, temperatures,3);

	//sdruns.insert({sdrac,SDRun(sd, ff, shakebds)});
	SDRun sdrun(sd, ff, shakebds);
	//SDRun &sdrun = sdruns.at(sdrac);
	sdrun.serial() = serial;
	//shake
	int doshake;
	pset.getval("DoShake", &doshake);
	sdrun.shakebds()->seton(doshake != 0);
	//acts
	//auto acts = std::shared_ptr<ActiveSelections>(new ActiveSelections(ff.get(), sdrac.allfixed, sdrac.mcfixed));
	//sdrun.setactiveselections(acts.get());
	sdrac.acts.init(ff.get(), sdrac.allfixed, sdrac.mcfixed);
	sdrun.setactiveselections(&(sdrac.acts));
	//auto f1=acts->atomfixed();
	//std::cout <<f1.size() <<std::endl;
	//auto f2=sdrun.acts()->atomfixed();
	//std::cout <<f2.size() <<std::endl;
	//sdrun.acts()->init(ff, sdrac.allfixed, sdrac.mcfixed);
	//shakeon
	sdrun.shakeon() = doshake != 0;
	//others
	sdrun.temperaturegroups()=agrps;
	sdrun.bathtemperatures()=temperatures;
	int nblsteps;
	pset.getval("NeighborListSteps", &nblsteps);
	sdrun.nblsteps() = nblsteps;
	int seed=NSPdstl::RandomEngine<>::getinstance().intrng(0,999999)();
	sdrun.initrandomengine(seed);

	//sdruns.insert({sdrac,sdrun});
	//auto f21=sdrun.acts()->atomfixed();
	//std::cout <<f21.size() <<std::endl;
	return sdrun;
}

SDRun SDRunAssemble::SDRunInit(SDRunAssembleControl &sdrac,
		const std::vector<std::vector<NSPallatom::Residue>>&ress, double a2nm) {
	SDRun sdrun = FindSDRun(sdrac);
	//auto f2=sdrun.acts()->atomfixed();
	//std::cout <<f2.size() <<' ' <<sdrun.acts()->ff()->natoms() <<std::endl;
	auto &pset = SDControls::getparameterset(sdrac.controlname + "_sd");
	std::vector<bool> fixatoms(sdrun.ff()->natoms(),false);
	int fixmainchain;
	pset.getval("FixMainChain",&fixmainchain);
	if(fixmainchain!=0){
		std::vector<int> mcatoms=sdrun.ff()->mainchainatoms();
		for(int a:mcatoms) fixatoms[a]=true;
	}
	std::vector<double> crds(fixatoms.size()*3);
	std::vector<std::vector<std::map<std::string,int>>> serial=sdrun.serial();
	for(int i=0;i<serial.size();i++) {
		for(int j=0;j<serial[i].size();j++) {
			for(auto &p:serial[i][j]) {
				XYZ c=ress[i][j].rds().at(p.first).crd;
				c = c * a2nm;
				int n = p.second*3;
				crds[n] = c.x_;
				crds[n+1] = c.y_;
				crds[n+2] = c.z_;
			}
		}
	}
	SDRun::SDRunIn in(crds,fixatoms);
	//std::cout <<"st" <<std::endl;
	//auto f2=sdrun.acts()->atomfixed();
	//std::cout <<f2.size() <<std::endl;
	sdrun.initstate(in);
	return sdrun;
}


