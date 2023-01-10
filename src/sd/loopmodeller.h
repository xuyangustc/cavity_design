/*
 * loopmodeller.h
 *
 *  Created on: 2018年7月10日
 *      Author: hyliu
 */

#ifndef LOOPMODELLER_H_
#define LOOPMODELLER_H_
#include "fullsite/fullsite.h"
#include "dataio/parameters.h"
#include "sd/genchain.h"
#include "sd/emrun.h"
#include <memory>

namespace NSPsd{
// a position in possibly multi-chain protein
struct ChainPosi{
	int cid{0}; //chain id
	int posi{-1}; //position id
	ChainPosi(){;}
	ChainPosi(int p):posi(p){;}
	ChainPosi(int c,int p):cid(c),posi(p){;}
	//make it comparable
	friend bool operator<(const ChainPosi &cp1, const ChainPosi &cp2){
		if(cp1.cid<cp2.cid) return true;
		else if(cp1.cid==cp2.cid) return cp1.posi<cp2.posi;
		else return false;
	}
};

// Specifies a part of the basechains_, which is to be substituted by the rebuilt loop
class TargetLoop{
public:
	TargetLoop(){;}
	TargetLoop(const std::vector<std::vector<NSPproteinrep::FullSite>>&basec,
			int pos1,int pos2): basechains_(&basec),
		startcp_(pos1),endcp_(pos2){;}
	TargetLoop(const std::vector<std::vector<NSPproteinrep::FullSite>>&basec,
			int cid1,int pos1,int cid2,int pos2):basechains_(&basec),
		startcp_(cid1,pos1),endcp_(cid2,pos2){;}
	friend bool operator<(const TargetLoop &tl1, const TargetLoop &tl2){
			return tl1.startcp_<tl2.startcp_;
		}
	bool isstartposi(int c,int p){
		return (c==startcp_.cid) &&(p==startcp_.posi);
	}
	bool isendposi(int c,int p){
		return (c==endcp_.cid) &&(p==endcp_.posi);
	}
	bool contains(int c, int p){
		if(isstartposi(c,p)) return true;
		ChainPosi cp(c,p);
		return (startcp_<cp && cp<endcp_);
	}
	//extract the structures in the basechain
	std::vector<NSPproteinrep::FullSite> extractbaseconfiguration() const;//TODO
	//retrieve a generated loop, remove it from the reservior of new conformers_
	std::shared_ptr<std::vector<NSPproteinrep::FullSite>> popconformer(const
	std::vector<std::string> &newsequence=std::vector<std::string>()){
		std::shared_ptr<std::vector<NSPproteinrep::FullSite>> res=getconformer(newsequence);
		if(res) newconformers_.resize(newconformers_.size()-1);
		return res;
	}
	//retrieves a generated loop
	std::shared_ptr<std::vector<NSPproteinrep::FullSite>> getconformer(
			const std::vector<std::string> &newsequence=std::vector<std::string>()){
			bool sequencechanged{false};
			if(!newsequence.empty()){
				if(newsctypes_!=newsequence) sequencechanged=true;
				newlength_=newsequence.size();
				newsctypes_=newsequence;
			}
			assert(newlength_>0);
			if(sequencechanged || newconformers_.empty()){
				if(!buildnewconformers()) return nullptr;
			}
			return newconformers_.back();
	}
/*	std::shared_ptr<std::vector<NSPproteinrep::FullSite>> getconformer(
			const std::vector<std::string> &newsequence=std::vector<std::string>()){
			if(!newsequence.empty()){
				newlength_=newsequence.size();
				newsctypes_=newsequence;
			}
			assert(newlength_>0);
			if(newconformers_.empty()){
				if(!buildnewconformers()) return nullptr;
			}
			return newconformers_.back();
	}*/
	std::vector<ChainPosi> getchainposi() {
		std::vector<ChainPosi> cps;
		cps.push_back(startcp_);
		cps.push_back(endcp_);
		return cps;
	}
	std::vector<int> sten() {
		std::vector<int> cp;
		cp.push_back(startcp_.cid);
		cp.push_back(startcp_.posi);
		cp.push_back(endcp_.cid);
		cp.push_back(endcp_.posi);
		return cp;
	}
private:
	//* starting and ending positions of the target loops in the basechains
	ChainPosi startcp_,endcp_;
	//* base chains
	const std::vector<std::vector<NSPproteinrep::FullSite>> *basechains_{nullptr};
	//* intended length of the new loop
	int newlength_{-1};
	//* intended amino acid sequences of the new loop
	std::vector<std::string> newsctypes_;
	//* stores loop conformations generated form backbone builder
	std::vector<std::shared_ptr<std::vector<NSPproteinrep::FullSite>>> newconformers_;
	//calls backbonebuilder to generate random but properly closed loop structures
	//the generated structures are stored in newconformers_
	bool buildnewconformers(); //TODO
};

//stores the previously optimized conformations of all the target loops.
class ArchivedLoops{
public:
	ArchivedLoops(){;}
	ArchivedLoops(double rcut):rmsdcut_(rcut){;}
	//try to add another optimized conformer to the archive
	//with restrict : RMSD between any two model must be larger than rmsdcut_ , total number can not go over maxnum
	bool addmodel(const std::vector<double> &model, const std::vector<double> &energy,
			const std::vector<double> &rmsdrecord, std::set<int> &removed, int maxnum);
	//no restrict
	bool addmodel(const std::vector<double> &model, const std::vector<double> &energy) {
		confcrds_.push_back(model);
		energies_.push_back(energy);
	}
	template<class GENCHAIN>
	std::vector<std::vector<std::vector<NSPproteinrep::FullSite>>> conformers(GENCHAIN &gench) {
		std::vector<std::vector<std::vector<NSPproteinrep::FullSite>>> confs;
		for(const auto &cc:confcrds_) {
			std::vector<std::vector<NSPproteinrep::FullSite>> sites;
			confs.push_back(gench.crd2fullsite(cc,10.0));
		}
		return confs;
	}
	std::vector<std::vector<double>> & confcrds() {
		return confcrds_;
	}
	std::vector<std::vector<double>> & energys() {
		return energies_;
	}
private:
	double rmsdcut_{0.0};
	//std::vector<std::shared_ptr<std::vector<NSPproteinrep::FullSite>>> conformers_;//A
	std::vector<std::vector<double>> energies_;
	std::vector<std::vector<double>> confcrds_;//NM
};

class LoopModeller{
public:
	LoopModeller(){;}
	LoopModeller(NSPdataio::ParameterSet &pset);
	std::vector<std::vector<NSPproteinrep::FullSite>> getainstance(std::vector<std::vector<bool>> &newfixed);
	std::vector<std::vector<NSPproteinrep::FullSite>> getainstance() {
		std::vector<std::vector<bool>> newfixed;
		return getainstance(newfixed);
	}
	std::vector<double> opt(std::vector<double> &crd);
	void makeinitialconf(std::vector<double> &crd, std::vector<double> &ene);
	void sample(const std::vector<double> &initcrd, const std::vector<double> &initene);
	void run();
	void nativetest(std::string fnene, std::string fnpdb, double prop);
	//void run(bool out=false, std::string dir="");
	void sample_rmsd_ene(const std::vector<double> &initcrd, const std::vector<double> &initene);
	void print(const std::vector<double> &crd, const std::vector<double> &es, int i);
	void native();


	std::vector<std::vector<NSPproteinrep::FullSite>>  & basechains() {
		return basechains_;
	}
	const std::vector<std::vector<NSPproteinrep::FullSite>>  & basechains() const {
		return basechains_;
	}
	const ArchivedLoops & archivedloops() const {
		return archivedloops_;
	}
	void printresults(std::string dir, std::string fnene);
private:
	std::vector<std::vector<NSPproteinrep::FullSite>> basechains_;
	std::shared_ptr<GenChain> genchain_;
	ArchivedLoops archivedloops_;
	std::vector<TargetLoop> loops_;
	std::vector<std::vector<std::string>> loopseq_;
	NSPdataio::ParameterSet pset_;
	//call from constructor
	void addtargetloop(int p1,int p2){
		loops_.push_back(TargetLoop(basechains_,p1,p2));
	}
	void addtargetloop(int c1, int p1,int c2,int p2){
		loops_.push_back(TargetLoop(basechains_,c1,p1,c2,p2));
	}
	std::shared_ptr<ForceField> ff_;
	std::shared_ptr<EMRun> emrun_;
	std::shared_ptr<SDRun> sdrun_;
	//EMRun emrun2;
	void remove_prev(std::map<std::pair<int,int>,std::vector<NSPproteinrep::FullSite>> &orifs, int ch, int po);
	void remove_back(std::map<std::pair<int,int>,std::vector<NSPproteinrep::FullSite>> &orifs, int ch, int po);
	void remove_middle(std::map<std::pair<int,int>,std::vector<NSPproteinrep::FullSite>> &orifs, int ch, int p1, int p2);
	void changeline(std::string controlname);

	bool knot(const std::vector<NSPproteinrep::FullSite> & chain);
	bool topoisright(const std::vector<std::vector<NSPproteinrep::FullSite>> & conf);
	bool topoisright(const std::vector<double> & crds) {
		std::vector<std::vector<NSPproteinrep::FullSite>> conf=genchain_->crd2fullsite(crds,10.0);
		return topoisright(conf);
	}
	std::string path_;
};
typedef NSPdataio::TypedMultiInstanceControls<LoopModeller> LoopControls;
void genchainreadcontrols_loop(const std::string &filename,std::string name);
void defineloopcontrol(std::string name,const std::vector<std::string> &controllines);
}



#endif /* LOOPMODELLER_H_ */
