/*
 * proteinreputil.cpp
 *
 *  Created on: 2018年4月2日
 *      Author: hyliu
 */

#include "integrateddesign/proteinreputil.h"
using namespace NSPproteinrep;
using namespace NSPdesignseq;
using namespace NSPsd;
PDB NSPproteinrep::makedesignseqPDB(const GenChain &genchain,
		const std::vector<double> & crd,
		double crdtoangstrom){
	std::ostringstream os_pdb;
	genchain.writepdb(crd,os_pdb,crdtoangstrom);
	std::istringstream is_pdb(os_pdb.str());
	return PDB(is_pdb,"NNNN");
}

std::vector<double> NSPproteinrep::sidechaintorsions(NSPdesignseq::Residue *residue){
	const VSCType & vsc=VSCType::getVSCType(residue->getType());
	const std::vector<int> & rotameratoms=vsc.rotameratoms;
	std::vector<double> result;
	const double rad=180.0/3.14159265;
	for(int l:rotameratoms){
		std::vector<NSPgeometry::XYZ> xyz;
		xyz.push_back(residue->getAtom(vsc.atomnames.at(l))->getCoord());
		for (int m=0;m<3;++m){
			int a=vsc.internalcrds[l][m].first;
			std::string aname;
			if(a==0) aname="N";
			else if(a==1) aname="CA";
			else aname=vsc.atomnames[a-2];
			xyz.push_back(residue->getAtom(aname)->getCoord());
		}
		double t=NSPgeometry::torsion(xyz[0],xyz[1],xyz[2],xyz[3])*rad;
		result.push_back(t);
	}
	return result;
}
