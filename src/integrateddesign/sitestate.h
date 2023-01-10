/*
 * siteproperties.h
 *
 *  Created on: 2018年4月2日
 *      Author: hyliu
 */
#include "designseq/StructureInfo.h"
#include "sd/genchain.h"
#ifndef SITESTATE_H_
#define SITESTATE_H_
namespace NSPproteinrep{
struct BackBoneState{
	double phi{360.0};
	double psi{360.0};
	double omiga{180.0};
	double sai{-1.0};
	int phipsistate{-1};
	char SSState{' '};
	int SASAState{-1};
};
struct SideChainState{
	std::string resiudetype{"UNK"};
	int rotamerstate{-1};
	std::vector<double> torsions;
};
struct SiteState{
	BackBoneState backbonestate;
	SideChainState sidechainstate;
};

std::vector<std::vector<SiteState>> sitestates(NSPdesignseq::PDB &pdb);
std::vector<std::vector<SiteState>> sitestates(const NSPsd::GenChain &genchain,
		const std::vector<double> &crd, double crdtoanstrom=1.0);

int degeneratedkai_loose(const std::string &resname);
int degeneratedkai_strict(const std::string &resname);
std::vector<double> kai_diffs(const SideChainState &s1,const SideChainState &s2,
		bool loose_degenerate=true);
}


#endif /* SITESTATE_H_ */
