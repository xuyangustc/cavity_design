/*
 * readpdb.h
 *
 *  Created on: Jun 11, 2018
 *      Author: xuyang
 */

#ifndef INTEGRATEDDESIGN_READPDB_H_
#define INTEGRATEDDESIGN_READPDB_H_

#include "fullsite/fullsite.h"
#include "proteinrep/pdbrecord.h"
#include "backbone/backbonesite.h"

namespace NSPproteinrep {
class PdbReader_xy {
public:
	typedef std::pair<int,char> ResKeyType; //residue number + insertion code
	typedef std::map<ResKeyType,std::vector<PdbRecord>> ResMapType;
	typedef std::map<char,ResMapType> RecordsType;
	PdbReader_xy() {;}
	void readpdb(std::string fn);
	PdbReader_xy(std::string fn, std::string pdbname=""):pdbid(pdbname) {
		readpdb(fn);
	}
	/*
	 * if MODEL, only read first MODEL
	 * read all atoms include "ATOM  " or "HETATM"
	 * read all atoms only from amino acid
	 * read all atoms which contain N,CA,C,O
	 * read atoms only from residue which are 20 standard amino acid
	 * read N,CA,C,O from backbone and discard incomplete backbone site
	 */
	bool model() {return model_;}
	RecordsType & getrecords() {return records_;}
	std::vector<std::vector<PdbRecord>> & noaa() {return aa_excluded;}
	std::vector<std::vector<BackBoneSite>> getbackbonesite();
private:
	std::string pdbid;
	bool model_{false};
	RecordsType records_;
	std::vector<std::vector<PdbRecord>> aa_excluded;
	//std::vector<PdbRecord> extractrecord(std::string fn);
};
}


#endif /* INTEGRATEDDESIGN_READPDB_H_ */
