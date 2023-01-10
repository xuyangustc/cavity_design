/*
 * genchain.h
 *
 *  Created on: 2018骞�1鏈�25鏃�
 *      Author: hyliu
 */

#ifndef SD_GENCHAIN_H_
#define SD_GENCHAIN_H_
#include "sd/forcefield.h"
#include "sd/sdrun.h"
#include "backbone/backbonesite.h"
#include "fullsite/fullsite.h"
#include "fullsite/structplus.h"
#include "allatom/pdbreader_xy.h"
#include <memory>
namespace NSPsd{
class GenChain{
public:
	GenChain(const std::string &controlname);
	std::vector<std::vector<std::string>> &sctypes(){return sctypes_;}
	const std::vector<std::vector<std::string>> &sctypes() const {return sctypes_;}
	std::string seqstring() const;
	void writepdb(const std::vector<double> & crd, std::ostream &os,double crdtoangstrom) const;
	std::vector<double> getcrd(const std::vector<std::vector<NSPproteinrep::BackBoneSite>> &sites,
			bool userefpdb=true) const;
	std::vector<double> getcrd(const std::string &pdbfile) const;
	std::vector<double> getcrd(const std::vector<std::vector<NSPproteinrep::FullSite>> & ) const;
	ForceField make_forcefield(const std::string &controlname, std::vector<std::set<int>> notcis=std::vector<std::set<int>>());
	StochasticDynamics make_sd(const std::string &controlname,
			const std::vector<std::vector<int>> &atomgroups,
			const std::vector<double> &temperatures) const;
	SDRun make_sdrun(const SDRun::SDRunIn &in, unsigned int seed, std::vector<std::set<int>> notcis=std::vector<std::set<int>>());
	SDRun make_sdrun(const SDRun::SDRunIn &in, unsigned int seed,
			std::string allfixed, std::string mcfixed,std::vector<std::set<int>> notcis=std::vector<std::set<int>>());
	ActiveSelections *setactiveselections(const ForceField *ff);
	ActiveSelections *setactiveselections(const ForceField *ff, std::string allfixed, std::string mcfixed);
	ActiveSelections *getactiveselections() const {
		return acts_.get();
	}



	std::vector<std::vector<NSPproteinrep::FullSite>> crd2fullsite(const std::vector<double> &crds, double nm2a) const;
	//double calrmsd(std::vector<double> &crd1, std::vector<double> &crd2);
	//ActiveSelections *setactiveselections_xy(const ForceField *ff);
	GenChain(const std::string &controlname, const std::vector<std::vector<NSPproteinrep::FullSite>> &fsss);
	std::vector<double> getcrd(double crd2nm=0.1) {
		std::vector<double> crds=getcrd(*refpdb_);
		for(double &c:crds) c *= crd2nm;
		return crds;
	}
	std::vector<double> getcrd(const std::vector<std::vector<NSPproteinrep::FullSite>> &sites,
			std::map<std::pair<int,int>,std::pair<int,int>>&atomseq) const;
	static std::vector<double> getcrd(const std::vector<std::vector<NSPproteinrep::FullSite>> &sites,
			std::vector<std::pair<std::vector<int>,std::string>> &atomseq);
	std::shared_ptr<std::vector<std::vector<NSPproteinrep::FullSite>>> getrefpdb() {return refpdb_;}
	static std::vector<std::vector<NSPallatom::Residue>> crd2res(const std::vector<double> &crds,
			const std::vector<std::vector<NSPallatom::Residue>> &refres, double nm2a);
	std::map<std::vector<int>,std::map<std::string,int>> atomseq() {
		return atomidinseq;
	}
private:
	std::vector<std::vector<std::string>> sctypes_;
	std::vector<std::vector<bool>> softsc_;
	std::string controlname_;
	std::shared_ptr<std::vector<std::vector<NSPproteinrep::FullSite>>> refpdb_{nullptr};
	std::shared_ptr<NSPproteinrep::StructPlus> sprefpdb_{nullptr};
	std::shared_ptr<ActiveSelections> acts_;
	std::map<std::vector<int>,std::map<std::string,int>> atomidinseq;//chainid, resid, atomname, atomseqinallchain,,,,,given value in make_forcefield
};

typedef NSPdataio::TypedMultiInstanceControls<GenChain> GenChainControls;
void definegenchaincontrol(std::string name,const std::vector<std::string> &controllines);
void genchainreadcontrols(const std::string &filename,std::string name);
void genchainprintcontrols(std::string name,std::ostream &ofs);
}




#endif /* SD_GENCHAIN_H_ */
