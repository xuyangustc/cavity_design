/*
 * confbuilder.h
 *
 *  Created on: Dec 11, 2018
 *      Author: xuyang
 */

#ifndef ALLATOM_CONFBUILDER_H_
#define ALLATOM_CONFBUILDER_H_
#include "allatom/secondary_structure.h"
#include "geometry/calculators.h"
#include "backbone/backbonebuilder.h"
#include "dstl/randomengine.h"
#include "dataio/parameters.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
using namespace NSPproteinrep;

namespace NSPallatom {
class InitialConfBuilder {
public:
	InitialConfBuilder(NSPdataio::ParameterSet &pset);
	enum SSType {Fixed, Strand, Helix, Loop};
	static void printpdb(std::string outfile, const std::vector<std::vector<Residue>>& outpep);

	void unbendstrand_MC(std::vector<NSPproteinrep::BackBoneSite> &chain,
			std::vector<std::pair<double,double>> &dihs, bool firstisfixed);
	std::vector<Residue> buildfirststrand(NSPproteinrep::BackBoneSite &st);
	std::vector<NSPproteinrep::BackBoneSite> pairedstrandN2C(NSPproteinrep::BackBoneSite &st,
			const std::vector<NSPproteinrep::BackBoneSite>&chain1, int lth);
	std::vector<NSPproteinrep::BackBoneSite> pairedstrandC2N(NSPproteinrep::BackBoneSite &st,
			const std::vector<NSPproteinrep::BackBoneSite>&chain1, int lth);
	std::vector<Residue> buildsecondstrand(const std::vector<Residue>&firstchain,
			const NSPproteinrep::BackBoneSite &firstst);
	void hbpaired_MC(BackBoneSite &st1, BackBoneSite &st2, BackBoneSite &bs1, BackBoneSite &bs2,
			bool firstisn2c, bool firstisinhb, bool parrallel, std::vector<double> &stdih);
	std::vector<double> hbpaired_ergodic(BackBoneSite &st1, BackBoneSite &st2, BackBoneSite &bs1, BackBoneSite &bs2,
			bool firstisn2c, bool firstisinhb, bool parrallel);
	std::vector<int> hbpaired_SD(std::string stfile, std::string controlfile,
			std::vector<std::pair<int,int>> &hb01no_ori, std::vector<std::pair<int,int>> &hb01on_ori,
			std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01no,
			std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01on);
	bool hbpaired_SD(std::string controlfile, std::string pdbfile, std::string outfile, std::vector<int> &acts,
			std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01no_ori,
			std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01on_ori);
	bool extendstrand(std::string stfile, std::string controlfile, std::string outfile);
	void sheet(std::string outfile);
	bool addloop();

	bool buildconf(std::string stfile, std::string controlfile, std::string outfile);
	void printpdb(std::string outfile);

	void builddomain_only4nanobody(std::string outfile);
private:
	std::vector<std::vector<Residue>> fixed_;
	std::vector<std::pair<int,std::vector<Residue>>> design_;
	std::vector<int> length_; // correspond with design_
	std::vector<int> extendpair_; //int correspond to serial in design_, size==2
	std::map<int,double> c2nphi_; // correspond to design_, size<=2
	//std::map<int,std::map<int,char>> betacontanctpair_; //trans - 't' , parrallel - 'p'
	std::vector<int> seqinsheet_; //direction is extend1->extend0
	std::vector<bool> isparrallel_; //size==seqinsheet.size()-1, do not contain the last one
	std::pair<double,double> strandphirange_;
	std::pair<double,double> strandpsirange_;
	std::string strandaa_;
	std::string helixaa_;
	std::string loopaa_;
	double prec_;
	int unbendstep_;
	double steplength_;
	double erglth_;
	double hblth_;
	int mcstep_;
	double mcsteplth_;
	double mcrange_;
	int sdstep_;
	double hbforce_;
	std::pair<int,XYZ> unbenddir_;

	std::vector<NSPproteinrep::BackBoneSite> buildchainC2N(NSPproteinrep::BackBoneSite &st, double phi,
			const std::vector<std::pair<double,double>>& dihs);
	std::vector<NSPproteinrep::BackBoneSite> buildchainN2C(NSPproteinrep::BackBoneSite &st,
			const std::vector<std::pair<double,double>>& dihs);
	void rotphiN(std::vector<NSPproteinrep::BackBoneSite> &chain, int n, double ang_add);
	void rotpsiN(std::vector<NSPproteinrep::BackBoneSite> &chain, int n, double ang_add);
	void rotphiC(std::vector<NSPproteinrep::BackBoneSite> &chain, int n, double ang_add);
	void rotpsiC(std::vector<NSPproteinrep::BackBoneSite> &chain, int n, double ang_add);
	void changephipsi(std::vector<NSPproteinrep::BackBoneSite> &chain, int n, double ang_add, bool isphi, bool movedisn);
};
typedef NSPdataio::TypedMultiInstanceControls<InitialConfBuilder> InitialConfBuilderControls;
void readextendpar(std::string parfile);
}


#endif /* ALLATOM_CONFBUILDER_H_ */
