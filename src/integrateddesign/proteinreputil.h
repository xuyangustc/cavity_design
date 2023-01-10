/*
 * proteinreputil.h
 *
 *  Created on: 2018年4月2日
 *      Author: hyliu
 */

#ifndef PROTEINREPUTIL_H_
#define PROTEINREPUTIL_H_
#include "sd/genchain.h"
#include "designseq/StructureInfo.h"
namespace NSPproteinrep{
NSPdesignseq::PDB makedesignseqPDB(const NSPsd::GenChain &genchain,
		const std::vector<double> & crd,
		double crdtoangstrom=1.0);
std::vector<double> sidechaintorsions(NSPdesignseq::Residue *residue);

}




#endif /* PROTEINREPUTIL_H_ */
