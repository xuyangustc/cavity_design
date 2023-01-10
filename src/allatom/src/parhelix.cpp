/*
 * parhelix.cpp
 *
 *  Created on: 2020��5��23��
 *      Author: xuyang
 */
#include "allatom/parhelix.h"
#include "sd/sidechainff.h"
#include "allatom/basic_info.h"
#include "fullsite/fullsite.h"
#include "backbone/backbonesite.h"
using namespace NSPgeometry;
using namespace NSPallatom;
using namespace NSPproteinrep;

std::vector<XYZ> ParHelix::cacrd(int lth, int st) {
	double tana = tan(alpha);
	double phi00 = phi0 + dz*tana/rad0;
	double rca = rad1*cos(alpha);
	double rsa = rad1*sin(alpha);
	double az = omg0*rad0/tana;
	std::vector<XYZ> cs;
	for(int i=st;i<lth+st;i++) {
		double di = (double)i;
		double a0 = omg0*di+phi00;
		double c0 = cos(a0);
		double s0 = sin(a0);
		double a1 = omg1*di+phi1;
		double c1 = cos(a1);
		double s1 = sin(a1);
		double x = rad0*c0 + rad1*c0*c1 - rca*s0*s1;
		//std::cout <<x <<' ' <<rad0*c0 <<' ' <<rad1*c0*c1 <<' ' <<rca*s0*s1 <<std::endl;
		//exit(1);
		double y = rad0*s0 + rad1*s0*c1 + rca*c0*s1;
		double z = az*di - rsa*s1 + dz;
		cs.push_back(XYZ(x,y,z));
	}
	return cs;
}
/*
//abadon
std::vector<Residue> ParHelix::helix(int lth, int st) {
	std::vector<XYZ> cas = cacrd(lth, st);
	std::vector<std::vector<XYZ>> css;
	std::vector<XYZ> pu = BackBoneSite::PeptideUnit(); //CA, C, O, N, CA

	/*Residue r1;
	r1.chid() = 'A';
	r1.chid_new() = 'A';
	r1.rid() = 0;
	r1.rid_new() = 0;
	r1.resname() = "GLY";
	Residue r2 = r2;
	r2.rid() = 0;
	r2.rid_new() = 0;
	std::map<std::string,AtomDataInPDB>& rds1 = r1.rds();
	rds1.insert({"CA",AtomDataInPDB()});
	rds1.at("CA").crd = pu[0];
	rds1.insert({"C",AtomDataInPDB()});
	rds1.at("C").crd = pu[1];
	rds1.insert({"O",AtomDataInPDB()});
	rds1.at("O").crd = pu[2];
	std::map<std::string,AtomDataInPDB>& rds2 = r2.rds();
	rds2.insert({"N",AtomDataInPDB()});
	rds2.at("N").crd = pu[3];
	rds2.insert({"CA",AtomDataInPDB()});
	rds2.at("CA").crd = pu[4];
	return {r1,r2};*/

/*	LocalFrame lf0=make_localframe(pu[0], pu[4], pu[3]);
	for(XYZ &c:pu) c = lf0.global2localcrd(c);
	std::vector<XYZ> qs{pu[1],pu[2],pu[3]};
	//std::cout <<pu[0].x_ <<' ' <<pu[0].y_ <<' ' <<pu[0].z_ <<std::endl;
	//std::cout <<pu[1].x_ <<' ' <<pu[1].y_ <<' ' <<pu[1].z_ <<std::endl;
	//std::cout <<pu[2].x_ <<' ' <<pu[2].y_ <<' ' <<pu[2].z_ <<std::endl;
	//std::cout <<pu[3].x_ <<' ' <<pu[3].y_ <<' ' <<pu[3].z_ <<std::endl;
	//std::cout <<pu[4].x_ <<' ' <<pu[4].y_ <<' ' <<pu[4].z_ <<std::endl;
	for(int i=0;i<cas.size()-1;i++) {
		XYZ y1=cas[i+1];
		y1.z_ -= 1.0;
		//if(i>3 && i<cas.size()-4) {
		//	XYZ c0 = cas[i-1] + cas[i-2] + cas[i-3] + cas[i-4];
		//	XYZ c1 = cas[i+1] + cas[i+2] + cas[i+3] + cas[i+4];
		//	c0 = (c0-c1)/4.0;
		//	c0 = c0/sqrt(c0.squarednorm());
		//	y1 = c0;
		//}
		std::vector<XYZ> cs=qs;
		//transform(pu[0],pu[4],pu[3],cas[i],cas[i+1],y1,cs);
		//std::cout <<cs[0].x_ <<' ' <<cs[0].y_ <<' ' <<cs[0].z_ <<std::endl;
		LocalFrame lf1=make_localframe(cas[i], cas[i+1], y1);
		for(XYZ &c:cs) c = lf1.local2globalcrd(c);
		css.push_back(cs);
	}
	XYZ nc, cc, oc;
	IdealGeometries & igdat = IdealGeometries::getGlobalInstance();
	nc = NSPgeometry::InternaltoXYZ(cas[0], css[0][0], css[0][1], igdat.idealLength("N", "CA"),
			igdat.idealAngle("N", "CA", "C"), igdat.idealTorsion("N", "CA", "C","O"));
	double degree = 3.14159265 / 180.0;
	double phi=-60.0*degree;
	cc = NSPgeometry::InternaltoXYZ(cas.back(), css.back()[2], css.back()[0],
			igdat.idealLength("C", "CA"),igdat.idealAngle("C", "CA", "N"), phi);
	oc = NSPgeometry::InternaltoXYZ(cc,cas.back(), css[0][2], igdat.idealLength("O", "C"),
			igdat.idealAngle("O", "C", "CA"), igdat.idealTorsion("N", "CA", "C","O"));
	Residue res;
	res.chid() = 'A';
	res.chid_new() = 'A';
	res.rid() = 0;
	res.rid_new() = 0;
	res.resname() = "GLY";
	std::map<std::string,AtomDataInPDB>& rds = res.rds();
	rds.insert({"C",AtomDataInPDB()});
	rds.at("C").element = " C";
	rds.insert({"CA",AtomDataInPDB()});
	rds.at("CA").element = " C";
	rds.insert({"N",AtomDataInPDB()});
	rds.at("N").element = " N";
	rds.insert({"O",AtomDataInPDB()});
	rds.at("O").element = " O";
	std::vector<Residue> ch;
	for(int i=0;i<lth;i++) {
		Residue r = res;
		r.rid() = i;
		r.rid_new() = i;
		r.rds().at("CA").crd = cas[i];
		r.rds().at("N").crd = nc;
		r.rds().at("C").crd = cc;
		r.rds().at("O").crd = oc;
		if(i!=lth-1) {
			r.rds().at("C").crd = css[i][0];
			r.rds().at("O").crd = css[i][1];
		}
		if(i!=0) {
			r.rds().at("N").crd = css[i-1][2];
		}
		ch.push_back(r);
	}
	return ch;
}
*/

std::vector<Residue> ParHelix::helix(int lth, int st) {
	std::vector<XYZ> pu = BackBoneSite::PeptideUnit(); //CA, C, O, N, CA
	LocalFrame lf0=make_localframe(pu[0], pu[4], pu[3]);
	for(XYZ &c:pu) c = lf0.global2localcrd(c);

	std::vector<std::vector<XYZ>> css;
	std::vector<XYZ> cas = cacrd(lth+2, st-1);
	std::vector<XYZ> qs{pu[1],pu[2],pu[3]};
	for(int i=0;i<cas.size()-1;i++) {
		XYZ y1=cas[i+1];
		y1.z_ -= 1.0;
		std::vector<XYZ> cs=qs;
		LocalFrame lf1=make_localframe(cas[i], cas[i+1], y1);
		for(XYZ &c:cs) c = lf1.local2globalcrd(c);
		css.push_back(cs);
	}

	Residue res;
	res.chid() = 'A';
	res.chid_new() = 'A';
	res.rid() = 0;
	res.rid_new() = 0;
	res.resname() = "GLY";
	std::map<std::string,AtomDataInPDB>& rds = res.rds();
	rds.insert({"C",AtomDataInPDB()});
	rds.at("C").element = " C";
	rds.insert({"CA",AtomDataInPDB()});
	rds.at("CA").element = " C";
	rds.insert({"N",AtomDataInPDB()});
	rds.at("N").element = " N";
	rds.insert({"O",AtomDataInPDB()});
	rds.at("O").element = " O";

	std::vector<Residue> ch;
	for(int i=0;i<css.size()-1;i++) {
		Residue r = res;
		r.rid() = i;
		r.rid_new() = i;
		r.rds().at("CA").crd = cas[i+1];
		r.rds().at("N").crd = css[i][2];
		r.rds().at("C").crd = css[i+1][0];
		r.rds().at("O").crd = css[i+1][1];
		ch.push_back(r);
	}
	assert(ch.size()==lth);
	return ch;
}

//0->1
void ParHelix::transform(XYZ o0, XYZ x0, XYZ y0, XYZ o1, XYZ x1, XYZ y1, std::vector<XYZ>&cs) {
	LocalFrame lf0=make_localframe(o0, x0, y0);
	LocalFrame lf1=make_localframe(o1, x1, y1);
	for(auto &c:cs) {
		//std::cout <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<std::endl;
		c = lf0.global2localcrd(c);
		//std::cout <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<std::endl;
		c = lf1.local2globalcrd(c);
		//std::cout <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<std::endl;
		//exit(1);
	}
}

