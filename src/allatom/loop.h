/*
 * loop.h
 *
 *  Created on: Jun 13, 2020
 *      Author: xuyang
 */

#ifndef ALLATOM_LOOP_H_
#define ALLATOM_LOOP_H_

#include "allatom/secondary_structure.h"
#include "geometry/calculators.h"
#include "geometry/rotation.h"
#include "backbone/backbonebuilder.h"
#include "dstl/randomengine.h"
#include "dataio/parameters.h"
#include "sd/genchain.h"
#include "allatom/filter.h"
#include "dstl/searchtree.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
using namespace NSPproteinrep;
using namespace NSPgeometry;
using namespace NSPsd;

namespace NSPallatom {
class LOOP {
public:
	static std::vector<std::vector<Residue>> NewChains(const std::vector<std::vector<Residue>>&chains,
			const std::vector<std::vector<Residue>>&lps, const std::vector<std::vector<int>>&chsten,
			std::vector<std::vector<int>>&loop_added);
	static std::vector<std::vector<Residue>> MakeNLoop(const std::vector<std::vector<Residue>>&chains,
			int nloop, int c0, int st, int c1, int en, int lth, int maxtry=1000);
	static std::vector<std::vector<Residue>> MakeNLinker(const std::vector<std::vector<Residue>>&chains,
			const std::vector<std::pair<int,int>> & helixregions, //first is start position, second is length
			const std::vector<std::pair<int,int>> & strandregions,//first is start position, second is length
			int nloop, int c0, int st, int c1, int en, int lth, int maxtry=1000);
	static NSPproteinrep::BackBoneSite getbbs(const std::vector<std::vector<Residue>>&chains, int c, int r);
private:
	//
};
class LoopLengthAndLocation {
public:
	//SS, length, weight, tree
	//static std::map<std::string,std::map<int,std::pair<double,SearchTree::TREE>>> trees;
	static std::map<int,std::pair<double,SearchTree::TREE>> Build1Tree(std::string fn);
	static std::map<std::string,std::map<int,std::pair<double,SearchTree::TREE>>> BuildTree();
	//static std::map<std::string,std::map<int,double>> criterion;
	static std::map<std::string,std::map<int,double>> Criterion(double maxval=20.0);
	//static std::map<std::string,std::map<int,std::vector<double>>> discut;
	static std::map<std::string,std::map<int,std::vector<double>>> DisCut();
	//static std::vector<double> ranges;
	static std::vector<double> Ranges();
	//static void InitAll(double maxval);

	LoopLengthAndLocation(const std::vector<std::vector<Residue>>&chs,
			const std::vector<std::string>&ss):chains(chs), secstr(ss) {
		;
	}
	std::map<int,double> score(int c0, int r0, int c1, int r1, int minlth, int maxlth);
	std::map<int,double> criterionfilter(int c0, int r0, int c1, int r1, std::map<int,double>&base);
	std::vector<std::vector<Residue>> looplengthtest(int c0, int r0, int c1, int r1, int lth);
	std::vector<std::pair<std::vector<int>,std::pair<double,std::vector<std::vector<Residue>>>>> linkfilter(
				int c0, int r0, int c1, int r1, bool test1more, std::map<int,double>&base);
	//return: cst, rst, cen, rst, lth, {score,vector<result>}
	//sts: cst, resst, resen
	//ens: cen, resst, resen
	//first step is do distance filter
	std::map<std::vector<int>,std::pair<double,std::vector<std::vector<Residue>>>> lengthandlocation(
			std::vector<int>&sts, std::vector<int>&ens, int minlth, int maxlth, bool test1more);
private:
	std::vector<std::vector<Residue>> chains;
	std::vector<std::string> secstr;
	LocalFrame nHelix(int c, int r, XYZ &c0);
	LocalFrame cHelix(int c, int r, XYZ &c0);
	LocalFrame nStrand(int c, int r, XYZ &c0);
	LocalFrame cStrand(int c, int r, XYZ &c0);
	std::vector<double> AddElement(XYZ&c0, LocalFrame&lf0, XYZ&c1, LocalFrame&lf1);
};
}


#endif /* ALLATOM_LOOP_H_ */
