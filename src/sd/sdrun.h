/*
 * SDRun.h
 *
 *  Created on: 2017骞�11鏈�9鏃�
 *      Author: hyliu
 */

#ifndef SD_SDRUN_H_
#define SD_SDRUN_H_
#include "sd/stochasticdynamics.h"
#include "sd/forcefield.h"
#include "sd/shakebonds.h"
#include "dstl/randomengine.h"
namespace NSPsd {
class SDRun;
class EmptySDCallBack{
public:
	bool operator() (SDRun &run) const { return false;}
};
class SDRun {
public:
	struct SDRunIn {
		const std::vector<double> *crd { nullptr };
		const std::vector<bool> *fixatom { nullptr };
		bool shakeinitcrd { true };
		SDRunIn(const std::vector<double> & initcrd,std::string controlname="") :
				crd(&initcrd),sdcontrolname(controlname) {
			;
		}
		SDRunIn(const std::vector<double> & initcrd,
				const std::vector<bool> &fixatoms) :
				crd(&initcrd), fixatom(&fixatoms) {
			;
		}
		std::vector<int> cissites;
		std::string sdcontrolname { "" };
		std::vector<std::string> sctypes;
	};
	SDRun() {
		;
	}
	SDRun(std::shared_ptr<StochasticDynamics> sd,
			std::shared_ptr<ForceField> ff,
			std::shared_ptr<ShakeBonds> shakebds) :
			sd_(sd), ff_(ff), shakebds_(shakebds) {
		;
	}
	/*	SDRun(const SDRun &run)  {
	 sd_=run.sd_;
	 ff_=run.ff_;
	 shakebds_=run.shakebds_;
	 }*/
	template<typename SEED>
	void initrandomengine(SEED seed) {
		rng_ = std::shared_ptr<NSPdstl::RandomEngine<>>(
				new NSPdstl::RandomEngine<>);
		rng_->init(seed);
		//std::cout << "First random number:" << rng_->realrng()() << std::endl;
	}
	std::shared_ptr<StochasticDynamics> sd() {
		return sd_;
	}
	std::shared_ptr<ForceField> ff() {
		return ff_;
	}
	std::shared_ptr<ShakeBonds> shakebds() {
		return shakebds_;
	}
	bool initstate(const SDRunIn & in) {
		auto &pset=SDControls::getparameterset(in.sdcontrolname+"_sd");
			vmasses_ = sd_->masses();
		int natoms = in.crd->size() / 3;
//		forceoff_.assign(natoms, false);
		forceoff_=acts_->atomfixed();
		nfixedcrd_=0;
		if(in.fixatom){
			for(int i=0;i<natoms;++i) forceoff_[i]= forceoff_[i] || in.fixatom->at(i);
		}
		for (int i = 0; i < natoms; ++i) {
				if (forceoff_[i]) {
					for (int idx = 3 * i; idx < 3 * i + 3; ++idx) {
						vmasses_[idx] = MASS_MAX;
						++nfixedcrd_;
					}
				}
			if(shakeon_) nfixedcrd_-=(nfixedcrd_/3-1);
		}
		state_ = sd_->make_initstate(*(in.crd), *rng_, &vmasses_);
		buffstate_ = std::shared_ptr < StochasticDynamics::State
				> (new StochasticDynamics::State);
		*buffstate_ = *state_;
		nstepsrun_ = 0;
		return sd_->shakestate(*state_, shakebds_->shakefunction(),
				in.shakeinitcrd, &vmasses_);
		//return true;
	}


	template<typename CALLBACK=EmptySDCallBack>
	bool runsteps(int nsteps, const CALLBACK &stepcallback=CALLBACK()) {
		std::shared_ptr<NeighborList> nbl;
		std::shared_ptr<GroupNeighbor> gn;
		for (int i = 0; i < nsteps; ++i) {
			if (i % nblsteps_ == 0) {
				if(ff_->stericwght()>0.000001)
				nbl = std::shared_ptr < NeighborList
						> (new NeighborList(state_->crd, *ff_, &forceoff_));
				//gn = std::shared_ptr < GroupNeighbor
				//		> (new GroupNeighbor(state_->crd, *ff_));
			}
			std::vector<double> forces = ff_->forces(state_->crd, *nbl, *gn,
					&potenergies_, *acts_);
			bool done = stepcallback(*this);
			//std::cout <<i <<' ' <<temperature() <<std::endl;
			if (done)
				return true;
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
	template<typename CALLBACK=EmptySDCallBack>
	bool runsteps_attract(int nsteps,const std::vector<std::pair<std::vector<int>,double>>&attrs,const CALLBACK &stepcallback=CALLBACK()) {
		std::shared_ptr<NeighborList> nbl;
		std::shared_ptr<GroupNeighbor> gn;
		for (int i = 0; i < nsteps; ++i) {
			if (i % nblsteps_ == 0) {
				if(ff_->stericwght()>0.000001)
				nbl = std::shared_ptr < NeighborList
						> (new NeighborList(state_->crd, *ff_, &forceoff_));
				//gn = std::shared_ptr < GroupNeighbor
				//		> (new GroupNeighbor(state_->crd, *ff_));
			}
			std::vector<double> forces = ff_->forces_attract(state_->crd, *nbl, *gn,&potenergies_, *acts_,attrs);
			bool done = stepcallback(*this);
			//std::cout <<i <<' ' <<temperature() <<std::endl;
			if (done)
				return true;
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
	template<typename CALLBACK=EmptySDCallBack>
	bool runsteps_fixedpoint_restrict(int nsteps, const std::vector<std::pair<std::pair<std::vector<int>,
				XYZ>,std::vector<double>>>&pointinfo, const CALLBACK &stepcallback=CALLBACK()) {
		std::shared_ptr<NeighborList> nbl;
		std::shared_ptr<GroupNeighbor> gn;
		for (int i = 0; i < nsteps; ++i) {
			if (i % nblsteps_ == 0) {
				if(ff_->stericwght()>0.000001)
				nbl = std::shared_ptr < NeighborList
						> (new NeighborList(state_->crd, *ff_, &forceoff_));
				//gn = std::shared_ptr < GroupNeighbor
				//		> (new GroupNeighbor(state_->crd, *ff_));
			}
			std::vector<double> forces = ff_->forces_fixedpoint_restrict(state_->crd,
					*nbl, *gn,&potenergies_, *acts_, pointinfo);
			bool done = stepcallback(*this);
			//std::cout <<i <<' ' <<temperature() <<std::endl;
			if (done)
				return true;
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
	template<typename CALLBACK=EmptySDCallBack>
	bool runsteps_fixedpoint_attraction_repulsion(int nsteps,
			const std::vector<ForceField::Fixedpoint_Attraction_Repulsion>&pointinfo,
			const ForceField::Site_Restrict_RMSD &sirerm, const CALLBACK &stepcallback=CALLBACK()) {
		std::shared_ptr<NeighborList> nbl;
		std::shared_ptr<GroupNeighbor> gn;
		for (int i = 0; i < nsteps; ++i) {
			if (i % nblsteps_ == 0) {
				if(ff_->stericwght()>0.000001)
				nbl = std::shared_ptr < NeighborList
						> (new NeighborList(state_->crd, *ff_, &forceoff_));
				//gn = std::shared_ptr < GroupNeighbor
				//		> (new GroupNeighbor(state_->crd, *ff_));
			}
			std::vector<double> forces = ff_->forces_fixedpoint_attraction_repulsion(state_->crd,
					*nbl, *gn,&potenergies_, *acts_, pointinfo, sirerm);
			bool done = stepcallback(*this);
			//std::cout <<i <<' ' <<temperature() <<std::endl;
			if (done)
				return true;
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
	StochasticDynamics::State & state() {
		return *state_;
	}
	const StochasticDynamics::State &state() const {
		return *state_;
	}
	int nstepsrun() const {
		return nstepsrun_;
	}
	double temperature() const {
		double nfreedof;
		if (shakeon_)
			nfreedof = state_->crd.size() -nfixedcrd_- shakebds_->nbonds();
		else
			nfreedof = state_->crd.size()-nfixedcrd_;
		return 2.0 * (sd_->ekin(*state_)) / (nfreedof * KBT);
	}
	const std::vector<double> &potenergies() const {
		return potenergies_;
	}
	bool & shakeon() {
		return shakeon_;
	}
	const bool & shakeon() const {
		return shakeon_;
	}
	int &nblsteps() {return nblsteps_;}
	const int & nblsteps() const {return nblsteps_;}
	void changetemperature(double kbt_new){
		double scale=kbt_new*KBT/sd_->kbT()[0];
		state_->scaletemp(scale);
		sd_->scaletemperatures(scale);
	}
	void changetemperature(double kbt_new,const std::vector<int> &atoms){
		double scale=kbt_new*KBT/sd_->kbT()[3*atoms[0]];
		state_->scaletemp(scale,atoms);
		sd_->scaletemperatures(scale,atoms,3);
	}
	void changetemperature(double kbt_new,int tgrp){
		changetemperature(kbt_new,temperaturegroups_->at(tgrp));
	}
	std::shared_ptr<std::vector<std::vector<int>>> & temperaturegroups(){
		return temperaturegroups_;}
	const std::shared_ptr<std::vector<std::vector<int>>> & temperaturegroups() const {
		return temperaturegroups_;}
	std::vector<double> &bathtemperatures() {return bathtemperatures_;}
	void setactiveselections(ActiveSelections *acts){
		acts_=acts;
	}



	void SymmetryOperation(const Rotation &rot, const std::vector<std::pair<int,int>>&corresp) {
		//if(!symmetryon_) return;
		std::vector<double> &cs=buffstate_->crd;
#pragma omp parallel for schedule(dynamic, 1)
		for(int i=0;i<corresp.size();i++) {
			int j=corresp[i].first*3;
			XYZ c(cs[j],cs[j+1],cs[j+2]);
			rot.apply(&c);
			j=corresp[i].second*3;
			cs[j] = c.x_;
			cs[j+1] = c.y_;
			cs[j+2] = c.z_;
		}
	}

	template<typename CALLBACK=EmptySDCallBack>
	bool runsteps_symmetry(int nsteps, const Rotation &rot,
			const std::vector<std::pair<int,int>>&corresp,
			const CALLBACK &stepcallback=CALLBACK()) {
		std::shared_ptr<NeighborList> nbl;
		std::shared_ptr<GroupNeighbor> gn;
		for (int i = 0; i < nsteps; ++i) {
			if (i % nblsteps_ == 0) {
				nbl = std::shared_ptr < NeighborList
						> (new NeighborList(state_->crd, *ff_, &forceoff_));
				//gn = std::shared_ptr < GroupNeighbor
				//		> (new GroupNeighbor(state_->crd, *ff_));
			}
			std::vector<double> forces = ff_->forces(state_->crd, *nbl, *gn,
					&potenergies_, *acts_);
			bool done = stepcallback(*this);
			//std::cout <<i <<' ' <<temperature() <<std::endl;
			if (done)
				return true;
			if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces,
					shakebds_->shakefunction(), &vmasses_))
				return false;
			SymmetryOperation(rot,corresp);
			auto temp = state_;
			state_ = buffstate_;
			buffstate_ = temp;
			++nstepsrun_;
		}
		return true;
	}



	bool smallchanged(std::vector<double> &records, int step_interval, double ene_diff);
	bool run_sample_maxrmsd(int nsteps, const std::vector<std::vector<double>> &confs,double index,
			int step_interval, double rmsd_diff, bool mconly, std::vector<double> &rdrd,
			std::vector<std::vector<double>> &enes);
	bool run_sample_rmsdcut(int nsteps, const std::vector<std::vector<double>> &confs, double index,
			double rmsdcut, bool mconly, std::vector<double> &rdrd, std::vector<std::vector<double>> &enes);
	bool run_sample_rmsdcut_ene(const std::vector<std::vector<double>> &confs,
			double lowrmsdcut, double highrmsdcut, double enecut, double index, bool mconly, int nsteps,
			std::vector<std::pair<int,double>> &minrmsdrecord, std::vector<std::vector<double>> &enes);
	bool run_opt(int nsteps, int step_interval, double ene_diff);
	bool resetstate(const std::vector<double> & initcrds);
	void changeweight(double &ppw, double &lsw, double &pkw, double &sccfw, double &scpkw, double &lhbw) {
		ff_->changeweight(ppw, lsw, pkw, sccfw, scpkw, lhbw);
	}
	void changeweight(double prop) {
		ff_->changeweight(prop);
	}
	std::vector<bool> forceoff() {return forceoff_;}
	std::vector<double> confene(const std::vector<double> &crds);
	bool runsteps_beta_sheet(int nsteps, const std::vector<std::pair<int,int>> &hbpairs,
			double f_index_hb, const std::vector<std::pair<int,int>> &capairs, double f_index_ca);
	bool runsteps_strand_unbend(int nsteps, double f_index, const std::vector<std::pair<int,int>> &capairs);
	bool runsteps_strand_hb(int nsteps, const std::vector<std::pair<int,int>> &hbpairs,
			double discutoff, double f_index);
	bool runsteps_subentry(int nsteps,
			const std::vector<std::pair<std::vector<int>,std::string>> &atomseq,
			EnergyTerm &eneterm, const std::vector<int> &sizes);
	void energyfordecoy(int nsteps,
			const std::vector<std::pair<std::vector<int>,std::string>> &atomseq,
			EnergyTerm &eneterm, const std::vector<int> &sizes);
	ActiveSelections* acts() {
		return acts_;
	}
	void printenergy();

	int runsteps_minE(int maxstep, int checkstep);

	std::vector<std::vector<std::map<std::string,int>>>& serial() {return serial_;}
private:
	std::shared_ptr<StochasticDynamics> sd_ { nullptr };
	std::shared_ptr<ForceField> ff_ { nullptr };
	std::shared_ptr<ShakeBonds> shakebds_ { nullptr };
	std::shared_ptr<StochasticDynamics::State> state_ { nullptr };
	std::shared_ptr<StochasticDynamics::State> buffstate_ { nullptr };
	std::shared_ptr<NSPdstl::RandomEngine<>> rng_ { nullptr };
	std::shared_ptr<std::vector<std::vector<int>>> temperaturegroups_{nullptr};
	ActiveSelections* acts_;
	std::vector<double> bathtemperatures_;
	std::vector<double> potenergies_;
	std::vector<double> vmasses_;
	std::vector<bool> forceoff_;
	int nfixedcrd_{0};
	int nstepsrun_ { 0 };
	int nblsteps_{50};
	bool shakeon_ { true };
	std::vector<std::vector<std::map<std::string,int>>> serial_;//chainseq, resseq, atomname, atomnumber
};

SDRun make_backbone_sdrun(const SDRun::SDRunIn &in, unsigned int randomseed);
void sdreadcontrols(const std::string &filename="sdffcontrols.par",std::string name="");
void sdprintcontrols(std::string name,std::ostream &ofs=std::cout);

SDRun make_sdrun_allatom(const std::vector<std::vector<std::string>> &scnames,
		const std::string &controlname, std::vector<std::vector<int>> &cissites);

struct SDRunAssembleControl {
public:
	std::vector<std::vector<std::string>> sequence;
	std::vector<std::vector<int>> cissites;
	std::map<int, std::vector<int>> allfixed;
	std::map<int, std::vector<int>> mcfixed;
	std::string controlname;
	ActiveSelections acts;
};
class Comp_SDRunAssembleControl {
public:
	bool operator() (SDRunAssembleControl s1, SDRunAssembleControl s2) {
		if(s1.sequence<s2.sequence) return true;
		else if(s1.sequence>s2.sequence) return false;
		else if(s1.cissites<s2.cissites) return true;
		else if(s1.cissites>s2.cissites) return false;
		else if(s1.allfixed<s2.allfixed) return true;
		else if(s1.allfixed>s2.allfixed) return false;
		else if(s1.mcfixed<s2.mcfixed) return true;
		else if(s1.mcfixed>s2.mcfixed) return false;
		else if(s1.controlname<s2.controlname) return true;
		return false;
	}
};
class SDRunAssemble {
public:
	SDRun FindSDRun(SDRunAssembleControl &sdrac);
	SDRun SDRunInit(SDRunAssembleControl &sdrac,
			const std::vector<std::vector<NSPallatom::Residue>>&ress, double a2nm);
private:
	std::map<SDRunAssembleControl,SDRun,Comp_SDRunAssembleControl> sdruns;
};
}

#endif /* SD_SDRUN_H_ */
