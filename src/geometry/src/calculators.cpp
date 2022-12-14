/*
 * geometry.cpp
 *
 *  Created on: 2015年10月30日
 *      Author: hyliu
 */

#include "geometry/calculators.h"
#include "geometry/rotation.h"
#include <string>

using namespace NSPgeometry;
double NSPgeometry::distance(const XYZ & p1, const XYZ & p2,std::vector<XYZ> *deriv){
	XYZ r12=p2-p1;
	double r=sqrt(r12.squarednorm());
	r12=(1.0/r)*r12;
	deriv->clear();
	deriv->push_back(-r12);
	deriv->push_back(r12);
	return r;
}
double NSPgeometry::angle(const XYZ &p1, const XYZ &p2, const XYZ & p3) {
	XYZ p21 = p1 - p2;
	XYZ p23 = p3 - p2;
	double theta;
	try {
		double ctheta = dot(p21, p23)
				/ std::sqrt(p21.squarednorm() * p23.squarednorm());
		if(ctheta<-1.0) ctheta=-1.0;
		if(ctheta>1.0) ctheta=1.0;
		theta = std::acos(ctheta);
	} catch (std::exception &e) {
		std::string msg { "Error calculating angle" };
//		NSPgeometry::attendError(msg);
		theta = NULLANGLE; //error indicating value
	}
	return theta;
}
double NSPgeometry::angle(const XYZ &p1, const XYZ &p2, const XYZ &p3,
                          std::vector<XYZ> *deriv) {
    XYZ p21 = p1 - p2;
    XYZ p23 = p3 - p2;
    double r21 = std::sqrt(p21.squarednorm());
    double r23 = std::sqrt(p23.squarednorm());
    double costheta = dot(p21,p23) / (r21*r23);
    double sintheta = std::sqrt(1.0-costheta*costheta);
    double dthetadx = 0.0;
    if (sintheta > 1.0e-5) dthetadx = -1.0 / sintheta;
    deriv->clear();
    deriv->push_back(dthetadx*(1/r21)*(p23/r23-p21*(costheta/r21)));
    deriv->push_back(XYZ());
    deriv->push_back(dthetadx*(1/r23)*(p21/r21-p23*(costheta/r23)));
    (*deriv)[1]=-((*deriv)[0]+(*deriv)[2]);
    return std::acos(costheta);
}
double NSPgeometry::cos_angle(const XYZ &p1, const XYZ &p2, const XYZ & p3,std::vector<XYZ> *deriv) {
	XYZ p21 = p1 - p2;
	XYZ p23 = p3 - p2;
	double r21=std::sqrt(p21.squarednorm());
	double r23=std::sqrt(p23.squarednorm());
	double costheta=dot(p21,p23)/(r21*r23);
	deriv->clear();
	deriv->push_back((1/r21)*(p23/r23-p21*(costheta/r21)));
	deriv->push_back(XYZ());
	deriv->push_back((1/r23)*(p21/r21-p23*(costheta/r23)));
	(*deriv)[1]=-((*deriv)[0]+(*deriv)[2]);
	return costheta;
}
double NSPgeometry::torsion(const XYZ &p1, const XYZ &p2, const XYZ &p3,
		const XYZ &p4) {
	XYZ p21 = p1 - p2, p23 = p3 - p2, p34 = p4 - p3;
	XYZ n213 = cross(p21, p23), n324 = cross(p34, p23);
	try {
		double cphi = dot(n213, n324)
				/ std::sqrt(n213.squarednorm() * n324.squarednorm());
		if(cphi > 1.0) cphi=1.0;
		if(cphi <-1.0) cphi=-1.0;
		double d = dot(n324, p21);
		if (d >= 0)
			return std::acos(cphi);
		else
			return -std::acos(cphi);
	} catch (std::exception &e) {
		std::string msg { "Error calculation torsion" };
		//	molNSPgeometry::attendError(msg);
		return NULLANGLE; //error indicating value
	}
}
double NSPgeometry::torsion(const XYZ &p1, const XYZ &p2, const XYZ &p3, const XYZ &p4,std::vector<XYZ> *deriv){
	XYZ p21 = p1 - p2, p23 = p3 - p2, p43 = p3 - p4;
	XYZ n213 = cross(p21, p23), n324 = cross(p23, p43);
	double r23_2=p23.squarednorm();
	double rn213_2=n213.squarednorm();
	double rn324_2=n324.squarednorm();
	double r23=std::sqrt(r23_2);
	double rn213=std::sqrt(rn213_2);
	double rn324=std::sqrt(rn324_2);
	double cphi = dot(n213, n324)
				/ (rn213*rn324);
	if(cphi > 1.0) cphi=1.0;
	if(cphi <-1.0) cphi=-1.0;
	double d = dot(n324, p21);
	double phi;
	if (d >= 0)
			 phi=std::acos(cphi);
	else
			phi=-std::acos(cphi);
    double ki = -r23 / rn213_2;
    double kl =   r23 / rn324_2;
    double kj1 = dot(p21, p23) / r23_2 - 1.0;
    double kj2 = dot(p43, p23) / r23_2;
    deriv->resize(4);
    deriv->at(0)=-ki * n213;
    deriv->at(3)=-kl*n324;
    deriv->at(1)=kj1*deriv->at(0) - kj2*deriv->at(3);
    deriv->at(2)=-(deriv->at(0)+deriv->at(1)+deriv->at(3));
    return phi;
}
double NSPgeometry::rmsd(const std::vector<XYZ> & crd1, const std::vector<XYZ> & crd2) {
		double d=0.0;
		for(int i=0; i<crd1.size(); ++i){
			d += (crd1[i]-crd2[i]).squarednorm();
		}
		return sqrt(d/(double) crd1.size());
}
XYZ NSPgeometry::center(const std::vector<XYZ> & points,std::vector<double> weights) {
	if(weights.empty()) {
		for(int i=0; i<points.size(); ++i) weights.push_back(1.0);
	}
	XYZ sum;
	double wtot=0.0;
	for(int i=0; i<points.size();++i) {
		sum =sum+weights[i]*points[i];
		wtot +=weights[i];
	}
	return sum/wtot;
}
XYZ NSPgeometry::InternaltoXYZ(const XYZ & rk, const XYZ &rj, const XYZ & ri,
		double b, double theta, double phi) {
	XYZ rjk = rk - rj, rji = ri - rj;
	try {
		XYZ ejk = rjk / std::sqrt(rjk.squarednorm());
		rji= rji- dot(rji,ejk)*ejk;
		XYZ eji = rji/ std::sqrt(rji.squarednorm());
		XYZ rl = rk - (b * cos(theta)) * ejk + (b * sin(theta)) * eji;
		Rotation r;
		r.init(QuaternionCrd(ejk,phi*180.0/3.14159265),rk);
		r.apply(&rl);
//		std::cout <<"interaltoxyz result: " <<b <<"\t"<< theta<<"\t" <<phi<<"\t"
//				<< std::sqrt((rl-rk).squarednorm()) <<"\t"
//				<<angle(rl,rk,rj) <<"\t"<< torsion(rl,rk,rj,ri)<<std::endl;
		return rl;
	} catch (std::exception &e) {
		std::string msg { "Error generating XYZ from internal" };
//				molNSPgeometry::attendError(msg);
		return XYZ(); //returns origin on error
	}
}
XYZ NSPgeometry::InternaltoXYZ(const XYZ & rk, double b) {
	return rk + XYZ(b, 0.0, 0.0);
}
XYZ NSPgeometry::InternaltoXYZ(const XYZ & rk, const XYZ & rj, double b,
		double theta) {
	XYZ rjk = rk - rj;
	try {
		XYZ ejk = rjk / std::sqrt(rjk.squarednorm());
		XYZ ei;
		double nxy = std::sqrt(ejk.x_ * ejk.x_ + ejk.y_ * ejk.y_);
		if (nxy > 0) {
			ei.x_ = ejk.y_ / nxy;
			ei.y_ = -ejk.x_ / nxy;
			ei.z_ = 0.0;
		} else {
			double nyz = std::sqrt(ejk.y_ * ejk.y_ + ejk.z_ * ejk.z_);
			ei.x_ = 0.0;
			ei.y_ = ejk.z_ / nyz;
			ei.z_ = -ejk.y_ / nyz;
		}
		XYZ rl = rk - (b * cos(theta)) * ejk + (b * sin(theta)) * ei;
		return rl;
	} catch (std::exception &e) {
		std::string msg { "Error generating XYZ from internal" };
//		molNSPgeometry::attendError(msg);
		return XYZ(); //returns origin on error
	}
}

double NSPgeometry::radiusgyr(const std::vector<NSPgeometry::XYZ> & crd) {
	NSPgeometry::XYZ center = NSPgeometry::center(crd);
	double rg=0.0;
	for (auto &c : crd) {
		rg += (c - center).squarednorm();
	}
	rg = sqrt(rg / (double) crd.size());
	return rg;
}

double NSPgeometry::radiusgyr(const std::vector<double> & x) {
	std::vector<XYZ> crd;
	for(int i=0;i<x.size()/3;++i) crd.push_back(XYZ(x[3*i],x[3*i+1],x[3*i+2]));
	return radiusgyr(crd);
}





std::vector<double> NSPgeometry::eq_line(XYZ &c1, XYZ &c2) {
	std::vector<double> pars;
	pars.push_back(c1.x_);
	pars.push_back(c1.x_-c2.x_);
	pars.push_back(c1.y_);
	pars.push_back(c1.y_-c2.y_);
	pars.push_back(c1.z_);
	pars.push_back(c1.z_-c2.z_);
	return pars;
}

int NSPgeometry::projection(XYZ &c1, XYZ &c2, XYZ &c3, XYZ &c4) {
	XYZ c1p=c1;
	c1p.z_=0.0;
	XYZ c2p=c1;
	c2p.z_=0.0;
	XYZ c3p=c1;
	c3p.z_=0.0;
	XYZ c4p=c1;
	c4p.z_=0.0;
	std::vector<double> line1p=eq_line(c1p,c2p);
	std::vector<double> line2p=eq_line(c3p,c4p);
	double x0,y0;
	if(line1p[1]==0.0) {
		if(line2p[1]==0.0) return 0;//two lines are parallel
		x0 = line1p[0];
		y0 = line2p[2] + line2p[3] * (line1p[0]-line2p[0]) / line2p[1];
	} else if(line1p[3]==0.0) {
		if(line2p[3]==0.0) return 0;//two lines are parallel
		x0 = line2p[0] + line2p[1] * (line1p[2]-line2p[2]) / line2p[3];
		y0 = line1p[2];
	} else {
		double k=line2p[1]/line1p[1]-line2p[3]/line1p[3];
		if(k==0.0) return 0;//two lines are parallel
		double t2 = ((line2p[0]-line1p[0]) / line1p[1] - (line2p[2]-line1p[2]) / line1p[3]) / k;
		x0 = line2p[0] + line2p[1] * t2;
		y0 = line2p[2] + line2p[3] * t2;
	}
	if(x0<c1p.x_ && x0<c2p.x_) return 0;//no intersection
	if(y0<c1p.y_ && y0<c2p.y_) return 0;//no intersection
	if(x0<c3p.x_ && x0<c4p.x_) return 0;//no intersection
	if(y0<c3p.y_ && y0<c4p.y_) return 0;//no intersection
	double z1,z2;
	std::vector<double> line1=eq_line(c1,c2);
	if(line1[1]==0.0) {
		if(line1[3]==0.0) return 0;//line1 is vertical to plane XOY
		z1 = line1[4] + line1[5] * (y0 - line1[2]) / line1[3];
	} else {
		z1 = line1[4] + line1[5] * (x0 - line1[0]) / line1[1];
	}
	std::vector<double> line2=eq_line(c3,c4);
	if(line2[1]==0.0) {
		if(line2[3]==0.0) return 0;//line1 is vertical to plane XOY
		z2 = line2[4] + line2[5] * (y0 - line2[2]) / line2[3];
	} else {
		z2 = line2[4] + line2[5] * (x0 - line2[0]) / line2[1];
	}
	return 1?2:z1>z2;
}

















