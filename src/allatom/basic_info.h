/*
 * basic_info.h
 *
 *  Created on: Nov 18, 2018
 *      Author: xuyang
 */

#ifndef ALLATOM_BASIC_INFO_H_
#define ALLATOM_BASIC_INFO_H_

#include "geometry/xyz.h"
#include <string>
#include <map>
#include <set>

namespace NSPallatom {
struct AA_Atom {
	static std::map<std::string,char> aa31();
	static std::map<char,std::string> aa13();
	static std::set<std::string> mainchain();
	static std::map<char,std::set<std::string>> sidechain();
	static std::map<std::string,std::set<std::string>> sidechain3();
	static std::map<char,std::map<std::string,double>> radius();
	static std::map<std::string,std::map<std::string,double>> radius3();
	static std::set<std::pair<char,std::string>> donors();
	static std::set<std::pair<char,std::string>> accepters();
	static std::map<std::pair<char,char>,int> BLOSUM62();
	static std::vector<NSPgeometry::XYZ> sphere256();
	static std::set<std::string> polarresidues();
	static std::map<std::string,std::map<std::string,std::set<std::vector<std::string>>>> sidechainHBunits();
	static std::map<std::string,std::vector<std::vector<std::string>>> groups(int index);
	static std::map<std::string,std::vector<std::set<std::string>>> groups();
	static std::map<std::string,std::vector<std::vector<std::string>>> groups_v();
	static std::map<std::string,std::vector<std::set<std::string>>> groups_center();
	static std::map<std::string,std::vector<std::vector<std::string>>> groups_crd();
	static std::map<std::string,std::vector<std::vector<std::string>>> groups_force();
	static std::map<std::vector<std::vector<std::string>>,double> atompairclashdis();
	static std::map<std::vector<std::string>,double> Weight_SAI();
	static std::map<std::string,std::vector<std::vector<std::string>>> Atom3();
	static std::map<std::string,std::map<std::string,std::vector<int>>> CADis();
	static std::vector<NSPgeometry::XYZ> PolyHedron20();
	static std::vector<NSPgeometry::XYZ> PolyHedron12();
};

//Given a vector {0,1,2,...,n-1} whose size is n, this function will return all kind of permutations, size is (n-1)!
std::vector<std::vector<int>> GetPermutation_norepeat(int n);
//Given a vector {0,1,2,...,n-1} whose size is n, this function will return all kind of permutations, size is n^n
std::vector<std::vector<int>> GetPermutation_repeat(int n);
//Given a vector<vector> {v0,v1,v2,v(n-1)}, v(n)={0,1,2,...,m}, extract 1 element from each vector
//the function will return all combination
//int represent size of subvector, so it can not be 0 or minus !!!!!
std::vector<std::vector<int>> ElementCombination(const std::vector<int> &vi);
//output is random order from 0 to n-1
std::vector<int> RandomOrder(int n);
}


#endif /* ALLATOM_BASIC_INFO_H_ */
