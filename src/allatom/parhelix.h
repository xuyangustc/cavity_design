/*
 * parhelix.h
 *
 *  Created on: 2020��5��23��
 *      Author: xuyang
 */

#ifndef ALLATOM_PARHELIX_H_
#define ALLATOM_PARHELIX_H_

#include "allatom/pdbreader_xy.h"
using namespace NSPgeometry;

namespace NSPallatom {

class ParHelix {
public:
	std::vector<XYZ> cacrd(int lth, int st=0);
	std::vector<Residue> helix(int lth, int st=0);
	const double pi{3.1415925};
	double rad0;
	double omg0;
	double phi0;
	double rad1;
	double omg1;
	double phi1;
	double alpha;
	double dr;
	double dz;
	void ideal() {
		double wl = -2.85*pi/180.0;
		double wh = 1.8*pi/180.0;
		if(omg0>wh) omg0 = wh;
		else if(omg0<wl) omg0 = wl;
		omg1 = 100.0*pi/180.0 - omg0;
		rad1 = 2.26;
		dr = 1.51;
		alpha = asin(rad0*omg0/dr);
	}
private:
	void transform(XYZ o0, XYZ x0, XYZ y0, XYZ o1, XYZ x1, XYZ y1, std::vector<XYZ>&cs);
};

class ParHelixBundle {
public:
private:
};

}


#endif /* ALLATOM_PARHELIX_H_ */
