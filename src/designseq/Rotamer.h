/*
 * Rotamer.h
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ROTAMER_H_
#define DESIGNSEQ_ROTAMER_H_

#include <string>
#include <vector>
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "designseq/AtomLib.h"


namespace NSPdesignseq {
using namespace std;

class Rotamer {

public:
	string rotName;
	string triName;
	vector<string> atomNameList;
	vector<NSPgeometry::XYZ> coordList;
	bool hasAromaticRing;
	NSPgeometry::XYZ normalVectorOfRing;


	Rotamer();
	Rotamer(const string& rotName, const string& triName);

	void addAtom(string atomName, NSPgeometry::XYZ localCoord);
	NSPgeometry::XYZ& getAtomCoord(const string& atomName);
	void updateLawOfRing();
	void buildSidechain(NSPgeometry::LocalFrame& cs, vector<NSPgeometry::XYZ>& xyzList);
	virtual ~Rotamer();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ROTAMER_H_ */
