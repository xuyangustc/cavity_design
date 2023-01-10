/*
 * polarenergy.h
 *
 *  Created on: Jun 28, 2018
 *      Author: xuyang
 */

#ifndef INTEGRATEDDESIGN_POLARENERGY_H_
#define INTEGRATEDDESIGN_POLARENERGY_H_

#include "integrateddesign/readpdb.h"
#include "dstl/randomengine.h"
#include "geometry/rotation.h"
#include "dstl/nbstree.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <map>
#include <cassert>
#include <cmath>
#include <set>
#include <memory>
#include <algorithm>
#include <iomanip>
using namespace NSPdstl;
using namespace NSPgeometry;

namespace NSPproteinrep {
void extractcrd_cc(std::string fnin, std::string fnout, double ccmax);
void extractcrd_on(std::string fnin, std::string fnout, double ccmax);
void transform_random_new(std::string fnin, std::string fnout,
		std::string fnoutnoclash, long numbg, double onmax, double ccmax, int resolution);
void randomchangebig(std::vector<XYZ> &crds, double rad, int center);
void randomchangesml(std::vector<XYZ> &crds, double rrad, double rang);
bool findsmldis(std::vector<XYZ> &c1, std::vector<XYZ> &c2, double &dis);
bool findsmldis(std::vector<XYZ> &crds, double &dis);
void transform_random_cc(std::string fnin, std::string fnout, std::string fnoutnoclash, int numrand, double ccmax);
void transform_random_on(std::string fnin, std::string fnout, std::string fnoutnoclash, int numrand, double ccmax);



void transform_random_sidechain(std::string fnin, std::string fnout, std::string fnoutnoclash, int ntimes, double radius);
void scgetparameter(std::string fnin, std::string fnout);

void random_check(std::string fnin, std::string fnlth, std::string fnang);
std::vector<long> random_seq(long lx);
void random_seq(std::string fnin, std::string fnout);
void getparameter(std::string fnin, std::string fnout);

class NeighborCountingScore {
public:
	void init(std::string fn);
	NeighborCountingScore(std::string fn) {
		init(fn);
	}
	void init(double lo, double li, double ri, double ro) {
		linner=li;rinner=ri;louter=lo;router=ro;
	}
	void settrainon(std::string fn) {
		probetrain=1;
		trainoutfile=fn;
	}
	void setrefon(std::string fn) {
		proberef=1;
		refoutfile=fn;
	}
	void buildtree(NBSTree<long> &tree, std::vector<std::vector<double>> &data,
			std::vector<long> &serial, std::vector<double> &ondis, std::string fn);
	void buildtree() {
		buildtree(traintree,traindata,trainserial,trainondis,traindatafile);
		buildtree(reftree,refdata,refserial,refondis,refdatafile);
	}
	double refscore(NBSQuery &query, NBSTree<long> &tree,
			std::vector<std::vector<double>> &data, long &totnum, long &num, bool ref, long seq);
	double trainscore(NBSQuery &query, NBSTree<long> &tree,
			std::vector<std::vector<double>> &data, long &totnum, long &num, bool ref, long seq);
	void getscore(NBSTree<long> &tree, std::vector<std::vector<double>> &data,
			std::vector<long> &serial, std::vector<double> &ondis, std::string fn, bool ref,long totnum);
	void getscore() {
		if(probetrain) getscore(traintree,traindata,trainserial,trainondis,trainoutfile,false,enetrainnum);
		if(proberef) getscore(reftree,refdata,refserial,refondis,refoutfile,true,enerefnum);
	}
	void testalpha(double alpha, std::string fn);
	void changereadfile(std::string trainin, std::string refin) {
		traindatafile=trainin;
		refdatafile=refin;
	}
	void changeoutfile(std::string trainout, std::string refout) {
		trainoutfile=trainout;
		refoutfile=refout;
	}
private:
	int ndim=6;
	std::vector<double> startval{-11.0,-11.0,-11.0,-11.0,-11.0,-11.0};
	std::vector<double> binwidth{0.001,0.001,0.001,0.001,0.001,0.001};
	double onprob1=2.0;
	double alphamin=0.001;
	double alphagrow=0.0005;
	double alphamax=0.01;
	double scoremin=1.0;
	double stadis0=4.0;
	double stadis1=4.25;
	double stadis2=4.5;
	long enetrainnum=50000;
	long enerefnum=10000;
	double linner;
	double rinner;
	double louter;
	double router;
	int probetrain=0;
	int proberef=0;
	std::string traindatafile;
	std::string refdatafile;
	std::string trainoutfile;
	std::string refoutfile;

	NBSTree<long> traintree;
	std::vector<std::vector<double>> traindata;
	std::vector<long> trainserial;
	std::vector<double> trainondis;
	NBSTree<long> reftree;
	std::vector<std::vector<double>> refdata;
	std::vector<long> refserial;
	std::vector<double> refondis;
};
void out_crdene(std::string enefile, std::string crdfile, std::string outfile, long trainnum, long refnum, double minsprop);
void getpar(std::vector<XYZ> &crds, double &con, double &cno, double &ncon, double &ocno, double &tor, double &lth);
void checksampletrain(std::string fnin, std::string fnout, int nsample);
void checksampleref(std::string fnin, std::string fnout, int nstart, int nperconf);
}


#endif /* INTEGRATEDDESIGN_POLARENERGY_H_ */
