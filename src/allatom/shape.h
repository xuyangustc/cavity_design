/*
 * shape.h
 *
 *  Created on: Aug 18, 2021
 *      Author: xuyang
 */

#ifndef ALLATOM_SHAPE_H_
#define ALLATOM_SHAPE_H_

#include "dstl/randomengine.h"
#include "backbone/backbonebuilder.h"
#include "allatom/pdbreader_xy.h"
#include "geometry/quatfit.h"

using namespace NSPgeometry;

namespace NSPallatom {
class Par2Helix {
public:
	//to generate HeilxA
	XYZ axisA; //axis of HelixA
	double hightA; //starting Z-coordinate of N-ter-Residue
	double transdis; //translating distance from N to C
	bool toup; //true(Z-coordinate of N2C are from low 2 high)
	int lth; //length of HelixA

	XYZ axisB; //reference axis
	double rotangle; //rotating angle based on axisB
	int nres; //rotangle based on (nres)th residue

	std::vector<Residue> HelixGen();
	void print(std::ostream &os) {
		os <<axisA.x_ <<' ' <<axisA.y_ <<' ' <<axisA.z_
				<<' ' <<hightA <<' ' <<transdis <<' ' <<toup <<' ' <<lth <<' '
				<<axisB.x_ <<' ' <<axisB.y_ <<' ' <<axisB.z_ <<' '
				<<rotangle <<' ' <<nres <<std::endl;
	}
	Par2Helix() {
		;
	}
	Par2Helix(std::string line) {
		init(line);
	}
	void init(std::string line) {
		std::stringstream ss(line);
		ss >>axisA.x_ >>axisA.y_ >>axisA.z_ >>hightA >>transdis >>toup >>lth
			>>axisB.x_ >>axisB.y_ >>axisB.z_ >>rotangle >>nres;
	}
	static std::vector<XYZ> HelixTerminal(const std::vector<Residue>&rs);
private:
	void ChainCrdCovert(std::vector<Residue>&ch, LocalFrame &lf0, LocalFrame &lf1);
	void ChainRotate(std::vector<Residue>&ch, Rotation &rt);
	void ChainTranslate(std::vector<Residue>&ch, XYZ &trans);
	XYZ HelixAxis(const std::vector<Residue>&rs);
	void ChainRotation(std::vector<Residue>&rs, double angle);
	std::vector<Residue> MakeChain(std::vector<std::vector<double>>&dihs, std::vector<std::string>&seq);
	std::vector<Residue> InitialHelix(int lth, bool toup);
};

class Shape {
public:
	std::vector<Par2Helix> pars;
	std::vector<double> hights;
	std::vector<std::vector<Residue>> chs;
	std::vector<std::vector<Residue>> looplib;
	std::vector<std::vector<XYZ>> loopcrd;
	std::map<int,std::map<int,std::set<std::string>>> allowed_par3;




	double Angle2Axis(XYZ fixedaxis, const std::vector<Residue>&helix, int nres);

	void ReadInputShapePar(std::string parfile);

	void FindBestLength(Par2Helix &par, double hightC);
	std::vector<XYZ> extractloopcrd(const std::vector<Residue>&lp);
	void extractloopcrd();
	std::vector<XYZ> extracthelixcrd(const std::vector<Residue>&ch0, int st,
			const std::vector<Residue>&ch1, int en);
	std::vector<std::vector<int>> LoopLibMatch(
			const std::vector<Residue>&ch0, const std::vector<Residue>&ch1);
	std::vector<Residue> Merge(std::vector<std::vector<int>> nlp,
			std::map<int,std::string>&aafixed, std::vector<int>&lths);
	double TransDis(double t);
	double RotAngle(double r);
	bool AddHelix(int idx, std::vector<std::vector<int>> &lpinfo);
	void ShapeGenRandom(std::string outdir, int nout);

	bool LoopLibMatch(const std::vector<Residue>&ch0,
			const std::vector<Residue>&ch1, double rmsdcut);
	void HelixPairParLib(double dis, double trans, std::string outfile);

	std::vector<std::pair<int,std::string>> FindRange(int idx);
	void ReadAllowedPar3(std::string parfile);
	bool Clash(Residue r1, Residue r2, bool a1, bool a2);
	bool Clash(std::vector<Residue>&ch, const std::map<int,std::string> &aas);
};
}


#endif /* ALLATOM_SHAPE_H_ */
