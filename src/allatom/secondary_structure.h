/*
 * secondary_structure.h
 *
 *  Created on: Nov 20, 2018
 *      Author: xuyang
 */

#ifndef ALLATOM_SECONDARY_STRUCTURE_H_
#define ALLATOM_SECONDARY_STRUCTURE_H_

#include "allatom/pdbreader_xy.h"
#include "allatom/basic_info.h"

namespace NSPallatom {
typedef std::pair<int,std::pair<int,int>> polypeptide;
struct Beta_Sheet {
	std::vector<polypeptide> sites;
	std::set<std::pair<std::pair<int,int>,std::pair<int,int>>> hbondpairs;//donor->accepter
};
class Secondary_Structure {
public:
	Secondary_Structure(const Peptide &p, int nh=4, int nb=3):pep_(p) {
		ss_divide();
		axis(nh,nb);
	}
	const std::vector<polypeptide> & helix() const {return alpha_helixes;}
	const std::vector<polypeptide> & strand() const {return beta_sheets;}
	const std::vector<std::vector<NSPgeometry::XYZ>> & h_axis() const {return helix_axis;}
	const std::vector<std::vector<NSPgeometry::XYZ>> & s_axis() const {return strand_axis;}
private:
	Peptide pep_;
	std::vector<std::vector<NSPgeometry::XYZ>> helix_axis;//center of 4 Ca
	std::vector<std::vector<NSPgeometry::XYZ>> strand_axis;//center of 2 Ca
	std::vector<polypeptide> alpha_helixes;
	std::vector<polypeptide> beta_sheets;
	void ss_divide();
	void axis(int nh, int nb);
	int hydrogenbond(const Residue &r1, const Residue &r2);
	std::set<std::pair<int,int>> helix_recognition(const std::vector<int> &ch_length,
			const std::set<std::pair<std::pair<int,int>,std::pair<int,int>>> &hbps);
	void strand_recognition(const std::vector<int> &ch_length,
			const std::map<std::pair<int,int>,std::set<std::pair<int,int>>> &n2c,
			const std::map<std::pair<int,int>,std::set<std::pair<int,int>>> &c2n);
};
}


#endif /* ALLATOM_SECONDARY_STRUCTURE_H_ */
