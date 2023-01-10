/*
 * sitestate.cpp
 *
 *  Created on: 2018年4月2日
 *      Author: hyliu
 */

#include "integrateddesign/sitestate.h"
#include "integrateddesign/proteinreputil.h"
#include "designseq/SasaPSD.h"
using namespace NSPproteinrep;
using namespace NSPdesignseq;
using namespace NSPsd;

std::vector<std::vector<SiteState>> NSPproteinrep::sitestates(
		NSPdesignseq::PDB &pdb){
	NSPdesignseq::StructureInfo strinfo(&pdb);
	strinfo.updateTorsion();
	strinfo.updateSecondaryStructure();
	SasaPSD sasapsd;
	strinfo.updateSAI(&sasapsd);
	std::vector<ProteinChain *> &pchains=pdb.getChains();
	int nchains=pchains.size();
	std::vector<std::vector<SiteState>> result;
	for(int c=0;c<nchains;++c){
		int nres=pchains[c]->getResList().size();
		result.push_back(std::vector<SiteState>(nres));
	}
	int seqid=0;

	for(int c=0;c<nchains;++c){
		for(auto &state:result[c]){
			BackBoneState &bs=state.backbonestate;
			SideChainState &scs=state.sidechainstate;
			bs.phi=strinfo.getPhi(seqid);
			bs.psi=strinfo.getPsi(seqid);
			bs.omiga=strinfo.getOmg(seqid);
			bs.SSState=strinfo.getSS(seqid);
			bs.sai=strinfo.getSai(seqid);
			Residue *residue=strinfo.getResidue(seqid);
			scs.resiudetype=residue->getType();
			scs.torsions=sidechaintorsions(residue);
			seqid++;
		}
	}
	return result;
}
std::vector<std::vector<SiteState>> NSPproteinrep::sitestates(const GenChain &genchain,
		const std::vector<double> &crd, double crdtoanstrom){
		PDB pdb=makedesignseqPDB(genchain,crd,crdtoanstrom);
		return sitestates(pdb);
}
std::vector<double> NSPproteinrep::kai_diffs(const SideChainState &s1, const SideChainState & s2,
		bool loose_degenerate){
	std::vector<double> diffs;
	int nkai=s1.torsions.size()<s2.torsions.size()? s1.torsions.size():s2.torsions.size();
	for(int i=0;i<nkai;++i){
		double diff=s2.torsions.at(i)-s1.torsions.at(i);
		if(diff>180.0) diff-=360.0;
		else if(diff<=-180.0) diff+=360.0;
		if(diff<-90.0 || diff>90.0){
			int (*degen)(const std::string &);
			if (loose_degenerate) degen= &(NSPproteinrep::degeneratedkai_loose);
			else degen=&(NSPproteinrep::degeneratedkai_strict);
			if(i==degen(s1.resiudetype) || i==degen(s2.resiudetype)){
				if(diff<-90.0) diff+=180.0;
				else diff -=180.0;
			}
		}
		diffs.push_back(diff);
	}
	return diffs;
}
int NSPproteinrep::degeneratedkai_loose(const std::string &resname){
	if(resname=="ASP" || resname=="ASN" ||
			resname=="HIS"||
			resname=="PHE" ||
			resname=="TYR") return 1;
	else if(resname=="GLU" ||resname=="GLN") return 2;
	else return -1;
}
int NSPproteinrep::degeneratedkai_strict(const std::string &resname){
	if(resname=="ASP" ||
			resname=="PHE" ||
			resname=="TYR") return 1;
	else if(resname=="GLU") return 2;
	else return -1;
}

