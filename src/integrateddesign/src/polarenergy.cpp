/*
 * polarenergy.cpp
 *
 *  Created on: Jun 28, 2018
 *      Author: xuyang
 */
#include "integrateddesign/polarenergy.h"
#include "sd/nnterm.h"
using namespace NSPproteinrep;
using namespace NSPgeometry;
using namespace NSPdstl;
using namespace NSPsd;

void NSPproteinrep::extractcrd_cc(std::string fnin, std::string fnout, double ccmax) {
	double comax=1.5;
	double cnmax=1.5;
	PdbReader_xy pr;
	pr.readpdb(fnin);
	std::vector<std::vector<BackBoneSite>> bbsss=pr.getbackbonesite();
	std::ofstream ofs;
	ofs.open(fnout,std::ofstream::app);
	int ix=0;
	for(auto &bbss:bbsss) {
		for(int i=0;i<bbss.size()-1;i++) {
			double c1x=bbss[i].data_[BackBoneSite::CCRD];
			double c1y=bbss[i].data_[BackBoneSite::CCRD+1];
			double c1z=bbss[i].data_[BackBoneSite::CCRD+2];
			double o1x=bbss[i].data_[BackBoneSite::OCRD];
			double o1y=bbss[i].data_[BackBoneSite::OCRD+1];
			double o1z=bbss[i].data_[BackBoneSite::OCRD+2];
			double n1x=bbss[i+1].data_[BackBoneSite::NCRD];
			double n1y=bbss[i+1].data_[BackBoneSite::NCRD+1];
			double n1z=bbss[i+1].data_[BackBoneSite::NCRD+2];
			double dx=c1x-o1x;
			double dy=c1y-o1y;
			double dz=c1z-o1z;
			double dd=dx*dx+dy*dy+dz*dz;
			if(dd>comax*comax) {
				std::cout <<fnin <<' ' <<i <<" CO " <<sqrt(dd) <<std::endl;
				continue;
			}
			dx=c1x-n1x;
			dy=c1y-n1y;
			dz=c1z-n1z;
			dd=dx*dx+dy*dy+dz*dz;
			if(dd>cnmax*cnmax) {
				std::cout <<fnin <<' ' <<i <<" CN " <<sqrt(dd) <<std::endl;
				continue;
			}
			for(int j=3;j<6&&i+j<bbss.size()-1;j++) {
				double c2x=bbss[i+j].data_[BackBoneSite::CCRD];
				double c2y=bbss[i+j].data_[BackBoneSite::CCRD+1];
				double c2z=bbss[i+j].data_[BackBoneSite::CCRD+2];
				double o2x=bbss[i+j].data_[BackBoneSite::OCRD];
				double o2y=bbss[i+j].data_[BackBoneSite::OCRD+1];
				double o2z=bbss[i+j].data_[BackBoneSite::OCRD+2];
				double n2x=bbss[i+j+1].data_[BackBoneSite::NCRD];
				double n2y=bbss[i+j+1].data_[BackBoneSite::NCRD+1];
				double n2z=bbss[i+j+1].data_[BackBoneSite::NCRD+2];
				dx=c2x-o2x;
				dy=c2y-o2y;
				dz=c2z-o2z;
				dd=dx*dx+dy*dy+dz*dz;
				if(dd>comax*comax) {
					std::cout <<fnin <<' ' <<i+j <<" CO " <<sqrt(dd) <<std::endl;
					continue;
				}
				dx=c2x-n2x;
				dy=c2y-n2y;
				dz=c2z-n2z;
				dd=dx*dx+dy*dy+dz*dz;
				if(dd>cnmax*cnmax) {
					std::cout <<fnin <<' ' <<i+j <<" CN " <<sqrt(dd) <<std::endl;
					continue;
				}
				dx=c1x-c2x;
				dy=c1y-c2y;
				dz=c1z-c2z;
				dd=dx*dx+dy*dy+dz*dz;
				if(dd>ccmax*ccmax) continue;
				ofs <<std::fixed <<std::setprecision(3) <<c1x <<'\t'
						<<std::fixed <<std::setprecision(3) <<c1y <<'\t'
						<<std::fixed <<std::setprecision(3) <<c1z <<'\t'
						<<std::fixed <<std::setprecision(3) <<o1x <<'\t'
						<<std::fixed <<std::setprecision(3) <<o1y <<'\t'
						<<std::fixed <<std::setprecision(3) <<o1z <<'\t'
						<<std::fixed <<std::setprecision(3) <<n1x <<'\t'
						<<std::fixed <<std::setprecision(3) <<n1y <<'\t'
						<<std::fixed <<std::setprecision(3) <<n1z <<'\t'
						<<std::fixed <<std::setprecision(3) <<c2x <<'\t'
						<<std::fixed <<std::setprecision(3) <<c2y <<'\t'
						<<std::fixed <<std::setprecision(3) <<c2z <<'\t'
						<<std::fixed <<std::setprecision(3) <<o2x <<'\t'
						<<std::fixed <<std::setprecision(3) <<o2y <<'\t'
						<<std::fixed <<std::setprecision(3) <<o2z <<'\t'
						<<std::fixed <<std::setprecision(3) <<n2x <<'\t'
						<<std::fixed <<std::setprecision(3) <<n2y <<'\t'
						<<std::fixed <<std::setprecision(3) <<n2z <<std::endl;
				ix++;
			}
		}
	}
	ofs.close();
	std::cout <<"get pairs: " <<ix <<"   ";
}

void NSPproteinrep::extractcrd_on(std::string fnin, std::string fnout, double ccmax) {
	double comax=1.5;
	double cnmax=1.5;
	PdbReader_xy pr;
	pr.readpdb(fnin);
	std::vector<std::vector<BackBoneSite>> bbsss=pr.getbackbonesite();
	std::ofstream ofs;
	ofs.open(fnout,std::ofstream::app);
	int ix=0;
	for(auto &bbss:bbsss) {
		for(int i=0;i<bbss.size()-1;i++) {
			double c1x=bbss[i].data_[BackBoneSite::CCRD];
			double c1y=bbss[i].data_[BackBoneSite::CCRD+1];
			double c1z=bbss[i].data_[BackBoneSite::CCRD+2];
			double o1x=bbss[i].data_[BackBoneSite::OCRD];
			double o1y=bbss[i].data_[BackBoneSite::OCRD+1];
			double o1z=bbss[i].data_[BackBoneSite::OCRD+2];
			double n1x=bbss[i+1].data_[BackBoneSite::NCRD];
			double n1y=bbss[i+1].data_[BackBoneSite::NCRD+1];
			double n1z=bbss[i+1].data_[BackBoneSite::NCRD+2];
			double dx=c1x-o1x;
			double dy=c1y-o1y;
			double dz=c1z-o1z;
			double dd=dx*dx+dy*dy+dz*dz;
			if(dd>comax*comax) {
				std::cout <<fnin <<' ' <<i <<" CO " <<sqrt(dd) <<std::endl;
				continue;
			}
			dx=c1x-n1x;
			dy=c1y-n1y;
			dz=c1z-n1z;
			dd=dx*dx+dy*dy+dz*dz;
			if(dd>cnmax*cnmax) {
				std::cout <<fnin <<' ' <<i <<" CN " <<sqrt(dd) <<std::endl;
				continue;
			}
			for(int j=3;j<6&&i+j<bbss.size()-1;j++) {
				double c2x=bbss[i+j].data_[BackBoneSite::CCRD];
				double c2y=bbss[i+j].data_[BackBoneSite::CCRD+1];
				double c2z=bbss[i+j].data_[BackBoneSite::CCRD+2];
				double o2x=bbss[i+j].data_[BackBoneSite::OCRD];
				double o2y=bbss[i+j].data_[BackBoneSite::OCRD+1];
				double o2z=bbss[i+j].data_[BackBoneSite::OCRD+2];
				double n2x=bbss[i+j+1].data_[BackBoneSite::NCRD];
				double n2y=bbss[i+j+1].data_[BackBoneSite::NCRD+1];
				double n2z=bbss[i+j+1].data_[BackBoneSite::NCRD+2];
				dx=c2x-o2x;
				dy=c2y-o2y;
				dz=c2z-o2z;
				dd=dx*dx+dy*dy+dz*dz;
				if(dd>comax*comax) {
					std::cout <<fnin <<' ' <<i+j <<" CO " <<sqrt(dd) <<std::endl;
					continue;
				}
				dx=c2x-n2x;
				dy=c2y-n2y;
				dz=c2z-n2z;
				dd=dx*dx+dy*dy+dz*dz;
				if(dd>cnmax*cnmax) {
					std::cout <<fnin <<' ' <<i+j <<" CN " <<sqrt(dd) <<std::endl;
					continue;
				}
				dx=o1x-n2x;
				dy=o1y-n2y;
				dz=o1z-n2z;
				double don=dx*dx+dy*dy+dz*dz;
				dx=n1x-o2x;
				dy=n1y-o2y;
				dz=n1z-o2z;
				double dno=dx*dx+dy*dy+dz*dz;
				dd=don;
				if(dno<don) dd=dno;
				if(dd>ccmax*ccmax) continue;
				ofs <<std::fixed <<std::setprecision(3) <<c1x <<'\t'
						<<std::fixed <<std::setprecision(3) <<c1y <<'\t'
						<<std::fixed <<std::setprecision(3) <<c1z <<'\t'
						<<std::fixed <<std::setprecision(3) <<o1x <<'\t'
						<<std::fixed <<std::setprecision(3) <<o1y <<'\t'
						<<std::fixed <<std::setprecision(3) <<o1z <<'\t'
						<<std::fixed <<std::setprecision(3) <<n1x <<'\t'
						<<std::fixed <<std::setprecision(3) <<n1y <<'\t'
						<<std::fixed <<std::setprecision(3) <<n1z <<'\t'
						<<std::fixed <<std::setprecision(3) <<c2x <<'\t'
						<<std::fixed <<std::setprecision(3) <<c2y <<'\t'
						<<std::fixed <<std::setprecision(3) <<c2z <<'\t'
						<<std::fixed <<std::setprecision(3) <<o2x <<'\t'
						<<std::fixed <<std::setprecision(3) <<o2y <<'\t'
						<<std::fixed <<std::setprecision(3) <<o2z <<'\t'
						<<std::fixed <<std::setprecision(3) <<n2x <<'\t'
						<<std::fixed <<std::setprecision(3) <<n2y <<'\t'
						<<std::fixed <<std::setprecision(3) <<n2z <<std::endl;
				ix++;
			}
		}
	}
	ofs.close();
	std::cout <<"get pairs: " <<ix <<"   ";
}

void NSPproteinrep::randomchangebig(std::vector<XYZ> &crds, double rad, int center) {
	XYZ trans0=crds[center];
	for(int i=0;i<crds.size();i++) crds[i] = crds[i]-trans0;
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	XYZ ori(0.0,0.0,0.0);
	QuaternionCrd qc(rng.realrng(0.0,1.0),0);
	Rotation rt(qc,ori);
	for(int i=0;i<crds.size();i++) rt.apply(&(crds[i]));
	XYZ rand(rng.realrng(0.0,1.0),rad);
	for(int i=0;i<crds.size();i++) crds[i]=crds[i]+rand;
}

void NSPproteinrep::randomchangesml(std::vector<XYZ> &crds, double rrad, double rang) {
	XYZ fixed(0.0,0.0,0.0);
	for(int i=0;i<crds.size();i++) fixed=fixed+crds[i];
	double size=(double)(crds.size());
	fixed.x_ /= size;
	fixed.y_ /= size;
	fixed.z_ /= size;
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	XYZ axis(rng.realrng(0,1.0),1.0);
	double ang=rng.realrng(0,rang)();
	QuaternionCrd qc(axis,ang);
	Rotation rt(qc,fixed);
	for(int i=0;i<crds.size();i++) rt.apply(&(crds[i]));
	XYZ rand(rng.realrng(0.0,1.0),rrad);
	for(int i=0;i<crds.size();i++) crds[i]=crds[i]+rand;
}

bool NSPproteinrep::findsmldis(std::vector<XYZ> &c1, std::vector<XYZ> &c2, double &dis) {
	XYZ don=c1[1]-c2[2];
	XYZ dno=c1[2]-c2[1];
	double don2=don.squarednorm();
	double dno2=dno.squarednorm();
	if(don2>=dno2) {
		dis=sqrt(dno2);
		return false;
	} else {
		dis=sqrt(don2);
		return true;
	}
}

bool NSPproteinrep::findsmldis(std::vector<XYZ> &crds, double &dis) {
	XYZ don=crds[1]-crds[5];
	XYZ dno=crds[2]-crds[4];
	double don2=don.squarednorm();
	double dno2=dno.squarednorm();
	if(don2>=dno2) {
		dis=sqrt(dno2);
		return false;
	} else {
		dis=sqrt(don2);
		return true;
	}
}

double ccprob1=2.0;
double daa=2.5;
double dhb=2.0;

void NSPproteinrep::transform_random_cc(std::string fnin, std::string fnout, std::string fnoutnoclash, int numrand, double ccmax) {
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	std::ifstream ifs(fnin);
	std::ofstream ofs(fnout);
	std::ofstream nocl(fnoutnoclash);
	std::string line;
	XYZ ori(0.0,0.0,0.0);
	long linenum=0;
	long randomnum=0;
	long noclnum=0;
	//std::cout <<"Probability=1.0 until distance2: " <<ccprob1 <<std::endl;
	//std::cout <<"Clash between atoms less than distance2: " <<daa <<std::endl;
	//std::cout <<"Clash between HBatoms less than distance2: " <<dhb <<std::endl;
	while(std::getline(ifs,line)) {
		//if(linenum%10000==0) std::cout <<"generating random pairs: " <<linenum <<' ' <<randomnum <<' ' <<noclnum <<std::endl;
		linenum++;
		std::vector<XYZ> fixedcrd(3), randomcrd(3);
		std::istringstream iss(line);
		for(int i=0;i<3;i++) {
			iss >>fixedcrd[i].x_;
			iss >>fixedcrd[i].y_;
			iss >>fixedcrd[i].z_;
		}
		assert((fixedcrd[0]-fixedcrd[1]).squarednorm()<1.51*1.51);
		assert((fixedcrd[0]-fixedcrd[2]).squarednorm()<1.51*1.51);
		for(int i=0;i<3;i++) {
			iss >>randomcrd[i].x_;
			iss >>randomcrd[i].y_;
			iss >>randomcrd[i].z_;
		}
		assert((randomcrd[0]-randomcrd[1]).squarednorm()<1.51*1.51);
		assert((randomcrd[0]-randomcrd[2]).squarednorm()<1.51*1.51);
		XYZ fixedtrans=fixedcrd[0];
		XYZ randomtrans=randomcrd[0];
		for(int i=0;i<3;i++) {
			fixedcrd[i] = fixedcrd[i]-fixedtrans;
			randomcrd[i] = randomcrd[i]-randomtrans;
		}
		for(int i=0;i<numrand;i++) {
			randomchangebig(randomcrd,ccmax,0);
			double d0;
			findsmldis(fixedcrd,randomcrd,d0);
			if(rng.realrng(0.0,1.0)()>ccprob1/d0/d0) continue;
			randomnum++;
			for(int j=0;j<3;j++) {
				ofs <<std::fixed <<std::setprecision(3) <<fixedcrd[j].x_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<fixedcrd[j].y_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<fixedcrd[j].z_ <<'\t';
			}
			for(int j=0;j<3;j++) {
				ofs <<std::fixed <<std::setprecision(3) <<randomcrd[j].x_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<randomcrd[j].y_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<randomcrd[j].z_ <<'\t';
			}
			ofs <<std::endl;
			bool clash{false};
			for(int j=0;j<3;j++) {
				for(int k=0;k<3;k++) {
					double dist2=(fixedcrd[j]-randomcrd[k]).squarednorm();
					if(j==1&&k==2 || j==2&&k==1) {
						if(dist2<dhb) {
							clash=true;
							break;
						}
					} else if(dist2<daa) {
						clash=true;
						break;
					}
				}
				if(clash) break;
			}
			if(clash) continue;
			noclnum++;
			for(int j=0;j<3;j++) {
				nocl <<std::fixed <<std::setprecision(3) <<fixedcrd[j].x_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<fixedcrd[j].y_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<fixedcrd[j].z_ <<'\t';
			}
			for(int j=0;j<3;j++) {
				nocl <<std::fixed <<std::setprecision(3) <<randomcrd[j].x_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<randomcrd[j].y_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<randomcrd[j].z_ <<'\t';
			}
			nocl <<std::endl;
		}
	}
	ifs.close();
	ofs.close();
	nocl.close();
	std::cout <<"total random pairs: " <<linenum <<' ' <<randomnum <<' ' <<noclnum <<std::endl;
}

void NSPproteinrep::transform_random_on(std::string fnin, std::string fnout, std::string fnoutnoclash, int numrand, double ccmax) {
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	std::ifstream ifs(fnin);
	std::ofstream ofs(fnout);
	std::ofstream nocl(fnoutnoclash);
	std::string line;
	XYZ ori(0.0,0.0,0.0);
	long linenum=0;
	long randomnum=0;
	long noclnum=0;
	//std::cout <<"Probability=1.0 until distance2: " <<ccprob1 <<std::endl;
	//std::cout <<"Clash between atoms less than distance2: " <<daa <<std::endl;
	//std::cout <<"Clash between HBatoms less than distance2: " <<dhb <<std::endl;
	while(std::getline(ifs,line)) {
		//if(linenum%10000==0) std::cout <<"generating random pairs: " <<linenum <<' ' <<randomnum <<' ' <<noclnum <<std::endl;
		linenum++;
		std::istringstream iss(line);
		std::vector<XYZ> crds(6);
		for(int i=0;i<crds.size();i++) {
			iss >>crds[i].x_;
			iss >>crds[i].y_;
			iss >>crds[i].z_;
		}
		double ron;
		std::vector<XYZ> fixedcrd(3), randomcrd(3);
		if(findsmldis(crds,ron)) {
			fixedcrd[0]=crds[0];
			fixedcrd[1]=crds[1];
			fixedcrd[2]=crds[2];
			randomcrd[0]=crds[3];
			randomcrd[1]=crds[4];
			randomcrd[2]=crds[5];
		} else {
			fixedcrd[0]=crds[3];
			fixedcrd[1]=crds[4];
			fixedcrd[2]=crds[5];
			randomcrd[0]=crds[0];
			randomcrd[1]=crds[1];
			randomcrd[2]=crds[2];
		}
		assert((fixedcrd[0]-fixedcrd[1]).squarednorm()<1.51*1.51);
		assert((fixedcrd[0]-fixedcrd[2]).squarednorm()<1.51*1.51);
		assert((randomcrd[0]-randomcrd[1]).squarednorm()<1.51*1.51);
		assert((randomcrd[0]-randomcrd[2]).squarednorm()<1.51*1.51);
		XYZ fixedtrans=fixedcrd[1];
		XYZ randomtrans=randomcrd[2];
		for(int i=0;i<3;i++) {
			fixedcrd[i] = fixedcrd[i]-fixedtrans;
			randomcrd[i] = randomcrd[i]-randomtrans;
		}
		for(int i=0;i<numrand;i++) {
			randomchangebig(randomcrd,ccmax,2);
			if(!findsmldis(fixedcrd,randomcrd,ron)) {
				i--;
				continue;
			}
			if(rng.realrng(0.0,1.0)()>ccprob1/ron/ron) continue;
			randomnum++;
			for(int j=0;j<3;j++) {
				ofs <<std::fixed <<std::setprecision(3) <<fixedcrd[j].x_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<fixedcrd[j].y_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<fixedcrd[j].z_ <<'\t';
			}
			for(int j=0;j<3;j++) {
				ofs <<std::fixed <<std::setprecision(3) <<randomcrd[j].x_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<randomcrd[j].y_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<randomcrd[j].z_ <<'\t';
			}
			ofs <<std::endl;
			bool clash{false};
			for(int j=0;j<3;j++) {
				for(int k=0;k<3;k++) {
					double dist2=(fixedcrd[j]-randomcrd[k]).squarednorm();
					if(j==1&&k==2 || j==2&&k==1) {
						if(dist2<dhb) {
							clash=true;
							break;
						}
					} else if(dist2<daa) {
						clash=true;
						break;
					}
				}
				if(clash) break;
			}
			if(clash) continue;
			noclnum++;
			for(int j=0;j<3;j++) {
				nocl <<std::fixed <<std::setprecision(3) <<fixedcrd[j].x_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<fixedcrd[j].y_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<fixedcrd[j].z_ <<'\t';
			}
			for(int j=0;j<3;j++) {
				nocl <<std::fixed <<std::setprecision(3) <<randomcrd[j].x_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<randomcrd[j].y_ <<'\t'
						<<std::fixed <<std::setprecision(3) <<randomcrd[j].z_ <<'\t';
			}
			nocl <<std::endl;
		}
	}
	ifs.close();
	ofs.close();
	nocl.close();
	std::cout <<"total random pairs: " <<linenum <<' ' <<randomnum <<' ' <<noclnum <<std::endl;
}

void NSPproteinrep::transform_random_new(std::string fnin, std::string fnout,
		std::string fnoutnoclash, long numbg, double onmax, double ccmax, int resolution) {
	std::vector<std::vector<XYZ>> crds;
	std::ifstream ifs(fnin);
	std::string line;
	while(std::getline(ifs,line)) {
		std::istringstream iss(line);
		std::vector<XYZ> co(6);
		for(int i=0;i<co.size();i++) {
			iss >>co[i].x_;
			iss >>co[i].y_;
			iss >>co[i].z_;
		}
		XYZ trans=co[0];
		co[0]=co[0]-trans;
		co[1]=co[1]-trans;
		co[2]=co[2]-trans;
		crds.push_back(co);
	}
	ifs.close();
	std::cout <<"background__readlines: " <<crds.size() <<std::endl;
	long nperbin=(long)(((double)numbg)/(onmax-ccprob1)/((double)resolution))+1;
	std::map<int,long> binrecord;
	for(int i=(int)(ccprob1*(double)resolution);i<(int)(onmax*(double)resolution);i++) binrecord.insert({i,nperbin});

	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	std::ofstream ofs(fnout);
	std::ofstream nocl(fnoutnoclash);

	XYZ ori(0.0,0.0,0.0);
	long randomnum=-1;
	long acceptnum=0;
	long noclnum=0;
	while(!binrecord.empty()) {
		randomnum++;
		if(randomnum==crds.size()) {
			std::cout <<randomnum <<' ' <<acceptnum <<' ' <<noclnum <<std::endl;
			randomnum=0;
		}
		std::vector<XYZ> fixedcrd(3), randomcrd(3);
		fixedcrd[0]=crds[randomnum][0];
		fixedcrd[1]=crds[randomnum][1];
		fixedcrd[2]=crds[randomnum][2];
		randomcrd[0]=crds[randomnum][3];
		randomcrd[1]=crds[randomnum][4];
		randomcrd[2]=crds[randomnum][5];
		randomchangebig(randomcrd,ccmax,0);
		double ron;
		findsmldis(fixedcrd,randomcrd,ron);
		if(ron<ccprob1) continue;
		int hblth=(int)((double)resolution*ron);
		if(binrecord.find(hblth)==binrecord.end()) continue;
		acceptnum++;
		binrecord.at(hblth)--;
		if(binrecord.at(hblth)==0) binrecord.erase(hblth);
		for(int j=0;j<3;j++) {
			ofs <<std::fixed <<std::setprecision(3) <<fixedcrd[j].x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<fixedcrd[j].y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<fixedcrd[j].z_ <<'\t';
		}
		for(int j=0;j<3;j++) {
			ofs <<std::fixed <<std::setprecision(3) <<randomcrd[j].x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<randomcrd[j].y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<randomcrd[j].z_ <<'\t';
		}
		ofs <<std::endl;
		bool clash{false};
		for(int j=0;j<3;j++) {
			for(int k=0;k<3;k++) {
				double dist2=(fixedcrd[j]-randomcrd[k]).squarednorm();
				if(j==1&&k==2 || j==2&&k==1) {
					if(dist2<dhb) {
						clash=true;
						break;
					}
				} else if(dist2<daa) {
					clash=true;
					break;
				}
			}
			if(clash) break;
		}
		if(clash) continue;
		noclnum++;
		for(int j=0;j<3;j++) {
			nocl <<std::fixed <<std::setprecision(3) <<fixedcrd[j].x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<fixedcrd[j].y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<fixedcrd[j].z_ <<'\t';
		}
		for(int j=0;j<3;j++) {
			nocl <<std::fixed <<std::setprecision(3) <<randomcrd[j].x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<randomcrd[j].y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<randomcrd[j].z_ <<'\t';
		}
		nocl <<std::endl;
	}
	ofs.close();
	nocl.close();
	std::cout <<"total random pairs: " <<randomnum <<' ' <<acceptnum <<' ' <<noclnum <<std::endl;
}

void NSPproteinrep::random_check(std::string fnin, std::string fnlth, std::string fnang) {
	std::vector<std::vector<double>> vvd;
	std::ifstream ifs(fnin);
	std::string readline;
	int ix=0;
	while(std::getline(ifs,readline)) {
		std::istringstream iss(readline);
		double d;
		std::vector<double> vd;
		for(int i=0;i<18;i++) {
			iss >>d;
			vd.push_back(d);
		}
		vvd.push_back(vd);
		ix++;
		if(ix%1000000==0) std::cout <<"readline: " <<ix <<std::endl;
	}
	ifs.close();

	std::map<int,long> disdist;
	std::map<int,long> angdist;
	ix=0;
	for(auto & vd:vvd) {
		double dx=vd[15]-vd[3];
		double dy=vd[16]-vd[4];
		double dz=vd[17]-vd[5];
		double dis=dx*dx+dy*dy+dz*dz;
		dx=vd[12]-vd[6];
		dy=vd[13]-vd[7];
		dz=vd[14]-vd[8];
		double dis2=dx*dx+dy*dy+dz*dz;
		if(dis2<dis) dis=dis2;
		dis=sqrt(dis);
		int idis=(int)(dis*10.0);
		if(disdist.find(idis)==disdist.end()) disdist.insert({idis,0});
		disdist.at(idis)++;
		double dx1=vd[3]-vd[0];
		double dy1=vd[4]-vd[1];
		double dz1=vd[5]-vd[2];
		double dx2=vd[12]-vd[9];
		double dy2=vd[13]-vd[10];
		double dz2=vd[14]-vd[11];
		double cos1=dx1*dx2+dy1*dy2+dz1*dz2;
		double l1=dx1*dx1+dy1*dy1+dz1*dz1;
		double l2=dx2*dx2+dy2*dy2+dz2*dz2;
		double cos2=cos1*cos1/l1/l2;
		double costheta=sqrt(cos2);
		double ang=acos(costheta);
		int iang=(int)(ang*100.0);
		if(angdist.find(iang)==angdist.end()) angdist.insert({iang,0});
		angdist.at(iang)++;
		ix++;
		if(ix%1000000==0) std::cout <<"getdist: " <<ix <<std::endl;
	}
	std::ofstream ofs;
	ofs.open(fnlth);
	for(auto &dd:disdist) ofs <<dd.first <<' ' <<dd.second <<std::endl;
	ofs.close();
	ofs.open(fnang);
	for(auto &ad:angdist) ofs <<ad.first <<' ' <<sin((double)(ad.first)/100.0) <<' ' <<ad.second <<std::endl;
	ofs.close();
}

std::vector<long> NSPproteinrep::random_seq(long lx) {
	long lx1=(long)(sqrt((double)lx)+1.0);
	std::vector<std::vector<long>> vvl;
	std::vector<std::vector<int>> vvi;
	std::vector<long> vl;
	std::vector<int> vi;
	for(long l=0;l<lx;l++) {
		if((l+1)%lx1==0) {
			vvl.push_back(vl);
			vvi.push_back(vi);
			vl.clear();
			vi.clear();
		}
		vl.push_back(l);
		vi.push_back(0);
	}
	if(!vl.empty()) {
		vvl.push_back(vl);
		vvi.push_back(vi);
	}
	std::vector<long> seq;
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	long size=vvl.size();
	std::vector<long> ss;
	for(int i=0;i<vvl.size();i++) ss.push_back(vvl[i].size());
	long lll=0;
	while(true) {
		if(size==0) break;
		long lx1=rng.intrng(0,size-1)();
		long lx2=-1;
		long lx3;
		for(long l=0;l<vvl.size();l++) {
			if(ss[l]!=0) {
				lx2++;
				if(lx2==lx1) {
					ss[l]--;
					lx3=l;
					break;
				}
			}
		}
		long lx4=rng.intrng(0,ss[lx3])();
		if(ss[lx3]==0) size--;
		long lx5;
		long lx6=-1;
		for(long l=0;l<=vvl[lx3].size();l++) {
			if(vvi[lx3][l]==0) {
				lx6++;
				if(lx6==lx4) {
					vvi[lx3][l]=1;
					lx5=vvl[lx3][l];
					break;
				}
			}
		}
		seq.push_back(lx5);
		lll++;
		//if(lll%100000==0) std::cout <<lll <<std::endl;
	}
	return seq;
}

void NSPproteinrep::random_seq(std::string fnin, std::string fnout) {
	std::vector<std::string> vs;
	std::ifstream ifs(fnin);
	std::string line;
	while(std::getline(ifs,line)) vs.push_back(line);
	ifs.close();
	std::vector<long> seq=random_seq(vs.size());
	std::ofstream ofs(fnout);
	for(long &l:seq) ofs <<vs[l] <<std::endl;
	ofs.close();
}

void NSPproteinrep::getparameter(std::string fnin, std::string fnout) {
	std::ifstream ifs(fnin);
	std::ofstream ofs(fnout);
	std::string line;
	while(std::getline(ifs,line)) {
		std::istringstream iss(line);
		std::vector<XYZ> crds(6);
		for(int i=0;i<6;i++) {
			iss >>crds[i].x_ >>crds[i].y_ >>crds[i].z_;
		}
		double ron;
		if(!findsmldis(crds,ron)) {
			std::vector<XYZ> crds1(6);
			crds1[0]=crds[3];
			crds1[1]=crds[4];
			crds1[2]=crds[5];
			crds1[3]=crds[0];
			crds1[4]=crds[1];
			crds1[5]=crds[2];
			crds=crds1;
		}
		double con1=angle(crds[0],crds[1],crds[5]);
		double cno2=angle(crds[3],crds[5],crds[1]);
		double ncon1=torsion(crds[2],crds[0],crds[1],crds[5]);
		double ocno2=torsion(crds[4],crds[3],crds[5],crds[1]);
		double cost1=cos(con1);
		double sint1=sin(con1);
		double cosp1=cos(ncon1);
		double sinp1=sin(ncon1);
		double cost2=cos(cno2);
		double sint2=sin(cno2);
		double cosp2=cos(ocno2);
		double sinp2=sin(ocno2);
		ofs <<std::fixed <<std::setprecision(4) <<ron*cost1 <<'\t'
				<<std::fixed <<std::setprecision(4) <<ron*sint1*cosp1 <<'\t'
				<<std::fixed <<std::setprecision(4) <<ron*sint1*sinp1 <<'\t'
				<<std::fixed <<std::setprecision(4) <<ron*cost2 <<'\t'
				<<std::fixed <<std::setprecision(4) <<ron*sint2*cosp2 <<'\t'
				<<std::fixed <<std::setprecision(4) <<ron*sint2*sinp2 <<'\t'
				<<std::fixed <<std::setprecision(4) <<ron <<std::endl;
	}
	ifs.close();
	ofs.close();
}

void NeighborCountingScore::init(std::string fn) {
	std::ifstream ifs(fn);
	std::string line;
	while(std::getline(ifs,line)) {
		std::istringstream iss(line);
		std::string label,str;
		iss >>label >>str;
		if(label=="ndim") iss>>ndim;
		else if(label=="onprob1") iss >>onprob1;
		else if(label=="alphamin") iss >>alphamin;
		else if(label=="alphamax") iss >>alphamax;
		else if(label=="alphagrow") iss >>alphagrow;
		else if(label=="scoremin") iss >>scoremin;
		else if(label=="stadis0") iss >>stadis0;
		else if(label=="stadis1") iss >>stadis1;
		else if(label=="stadis2") iss >>stadis2;
		else if(label=="enetrainnum") iss >>enetrainnum;
		else if(label=="enerefnum") iss >>enerefnum;
		else if(label=="linner") iss >>linner;
		else if(label=="rinner") iss >>rinner;
		else if(label=="louter") iss >>louter;
		else if(label=="router") iss >>router;
		else if(label=="probetrain") iss >>probetrain;
		else if(label=="proberef") iss >>proberef;
		else if(label=="traindatafile") iss >>traindatafile;
		else if(label=="refdatafile") iss >>refdatafile;
		else if(label=="trainoutfile") iss >>trainoutfile;
		else if(label=="refoutfile") iss >>refoutfile;
		else if(label=="startval") {startval.resize(ndim);for(int i=0;i<ndim;i++) iss>>startval[i];}
		else if(label=="binwidth") {binwidth.resize(ndim);for(int i=0;i<ndim;i++) iss>>binwidth[i];}
		/*switch (label) {
		case "ndim" :{iss >>ndim;} break;
		case "startval" :{startval.resize(ndim);for(int i=0;i<ndim;i++) iss>>startval[i];} break;
		case "binwidth" :{binwidth.resize(ndim);for(int i=0;i<ndim;i++) iss>>binwidth[i];} break;
		case "onprob1" :{iss >>onprob1;} break;
		case "alphamin" :{iss >>alphamin;} break;
		case "alphamax" :{iss >>alphamax;} break;
		case "alphagrow" :{iss >>alphagrow;} break;
		case "scoremin" :{iss >>scoremin;} break;
		case "stadis0" :{iss >>stadis0;} break;
		case "stadis1" :{iss >>stadis1;} break;
		case "stadis2" :{iss >>stadis2;} break;
		case "totnum" :{iss >>totnum;} break;
		case "linner" :{iss >>linner;} break;
		case "rinner" :{iss >>rinner;} break;
		case "louter" :{iss >>louter;} break;
		case "router" :{iss >>router;} break;
		case "probetrain" :{iss >>probetrain;} break;
		case "proberef" :{iss >>proberef;} break;
		case "traindatafile" :{iss >>traindatafile;} break;
		case "refdatafile" :{iss >>refdatafile;} break;
		case "trainoutfile" :{iss >>trainoutfile;} break;
		case "refoutfile" :{iss >>refoutfile;} break;
		}*/
	}
	ifs.close();
	traintree.init(startval,binwidth);
	reftree.init(startval,binwidth);
}

void NeighborCountingScore::buildtree(NBSTree<long> &tree, std::vector<std::vector<double>> &data,
		std::vector<long> &serial, std::vector<double> &ondis, std::string fn) {
	std::ifstream ifs(fn);
	std::string line;
	long lx=-1;
	while(std::getline(ifs,line)) {
		lx++;
		if(lx%1000000==0) std::cout <<"readline: " <<lx <<std::endl;
		std::istringstream iss(line);
		std::vector<double> vd(ndim);
		for(int i=0;i<ndim;i++) iss >>vd[i];
		double on;
		iss >>on;
		if(vd[0]>=router || vd[0]<louter) continue;
		data.push_back(vd);
		serial.push_back(lx);
		ondis.push_back(on);
	}
	ifs.close();
	for(long l=0;l<data.size();l++) {
		if(l%1000000==0) std::cout <<"treebuilding: " <<l <<std::endl;
		tree.addpoint(l,data[l]);
	}
}

double NeighborCountingScore::trainscore(NBSQuery &query, NBSTree<long> &tree,
		std::vector<std::vector<double>> &data, long &totnum, long &num, bool ref, long seq) {
	std::vector<long> nbs;
	tree.findneighbors(query,&nbs);
	num=totnum=nbs.size();
	double sc=0.0;
	for(long &l:nbs) {
		if(l==seq && !ref) {
			num--;
			continue;
		}
		double s=1.0;
		bool inrange{true};
		for(int i=0;i<ndim;i++) {
			double diff = fabs(traindata[l][i]-query.vals[i]);
			if(diff<query.lcuts[i]/2.0) continue;
			if(diff>query.lcuts[i]) {
				inrange=false;
				s=0.0;
				break;
			}
			double dd = (query.lcuts[i]-diff) / query.lcuts[i]*2.0;
			s *= dd*dd;
		}
		if(!inrange) {
			num--;
			continue;
		}
		sc += s;
	}
	return sc;
}

double NeighborCountingScore::refscore(NBSQuery &query, NBSTree<long> &tree,
		std::vector<std::vector<double>> &data, long &totnum, long &num, bool ref, long seq) {
	std::vector<long> nbs;
	tree.findneighbors(query,&nbs);
	num=totnum=nbs.size();
	double sc=0.0;
	for(long &l:nbs) {
		if(l==seq && ref) {
			num--;
			continue;
		}
		double s=1.0;
		bool inrange{true};
		for(int i=0;i<ndim;i++) {
			double diff = fabs(refdata[l][i]-query.vals[i]);
			if(diff<query.lcuts[i]/2.0) continue;
			if(diff>query.lcuts[i]) {
				inrange=false;
				s=0.0;
				break;
			}
			double dd = (query.lcuts[i]-diff) / query.lcuts[i]*2.0;
			s *= dd*dd;
		}
		if(!inrange) {
			num--;
			continue;
		}
		if(refondis[l]<onprob1) sc += s;
		else sc += s*refondis[l]*refondis[l]/onprob1/onprob1;
	}
	return sc;
}

void NeighborCountingScore::getscore(NBSTree<long> &tree, std::vector<std::vector<double>> &data,
		std::vector<long> &serial, std::vector<double> &ondis, std::string fn, bool ref, long totnum) {
	std::ofstream ofs(fn);
	long getnumber=0;
	double alpha=alphamin;
	double refsc=-1.0;
	long refnum=-1;
	long totrefnum=-1;
	for(long l=0;l<data.size()&&getnumber<totnum;l++) {
		//if(ondis[l]>stadis2+0.5) continue;
		if(data[l][0]<linner || data[l][0]>=rinner) continue;
		getnumber++;
		std::vector<double> &dt=data[l];

		double lthhalf=ondis[l]*ondis[l]*alpha;
		NBSQuery query;
		for(int i=0;i<ndim;i++) {
			query.vals.push_back(dt[i]);
			//query.lcuts.push_back(dt[i]-lthhalf*2.0);
			//query.rcuts.push_back(dt[i]+lthhalf*2.0);
			//double lthhalf=dt[i]*dt[i]*alpha;
			query.lcuts.push_back(lthhalf*2.0);
			query.rcuts.push_back(lthhalf*2.0);
		}
		refsc=refscore(query,reftree,refdata,totrefnum,refnum,ref,l);
		if(refsc<scoremin) {
			while(refsc<scoremin) {
				if(alpha+alphagrow>alphamax) break;
				alpha+=alphagrow;
				lthhalf=ondis[l]*ondis[l]*alpha;
				query.vals.clear();
				query.lcuts.clear();
				query.rcuts.clear();
				for(int i=0;i<ndim;i++) {
					query.vals.push_back(dt[i]);
					//query.lcuts.push_back(dt[i]-lthhalf*2.0);
					//query.rcuts.push_back(dt[i]+lthhalf*2.0);
					//double lthhalf=dt[i]*dt[i]*alpha;
					query.lcuts.push_back(lthhalf*2.0);
					query.rcuts.push_back(lthhalf*2.0);
				}
				refsc = refscore(query,reftree,refdata,totrefnum,refnum,ref,l);
			}
		} else if(refsc>scoremin) {
			long refnum1=refnum;
			long totnum1=totrefnum;
			double refsc1=refsc;
			double alpha1=alpha;
			double lthhalf1=lthhalf;
			NBSQuery query1=query;
			while(refsc1>scoremin) {
				alpha=alpha1;
				query.vals=query1.vals;
				query.lcuts=query1.lcuts;
				query.rcuts=query1.rcuts;
				alpha1-=alphagrow;
				query1.vals.clear();
				query1.lcuts.clear();
				query1.rcuts.clear();
				lthhalf1=lthhalf;
				for(int i=0;i<ndim;i++) {
					query1.vals.push_back(dt[i]);
					//query1.lcuts.push_back(dt[i]-lthhalf1*2.0);
					//query1.rcuts.push_back(dt[i]+lthhalf1*2.0);
					//double lthhalf=dt[i]*dt[i]*alpha;
					query1.lcuts.push_back(lthhalf1*2.0);
					query1.rcuts.push_back(lthhalf1*2.0);
				}
				if(alpha1<alphamin) break;
				refsc=refsc1;
				refnum=refnum1;
				totrefnum=totnum1;
				refsc1 = refscore(query1,reftree,refdata,totnum1,refnum1,ref,l);
			}
		}
		long tottrainnum;
		long trainnum;
		double trainsc = trainscore(query,traintree,traindata,tottrainnum,trainnum,ref,l);
		ofs <<serial[l] <<'\t' <<totrefnum <<'\t' <<refnum <<'\t' <<tottrainnum <<'\t' <<trainnum <<'\t'
				<<std::fixed <<std::setprecision(6) <<refsc <<'\t'
				<<std::fixed <<std::setprecision(6) <<trainsc <<'\t'
				<<std::fixed <<std::setprecision(4) <<ondis[l] <<'\t'
				<<std::fixed <<std::setprecision(4) <<alpha <<std::endl;
	}
	ofs.close();
	std::cout <<"  score total line: " <<getnumber <<std::endl;
}

void NeighborCountingScore::testalpha(double alpha, std::string fn) {
	std::ofstream ofs(fn,std::ofstream::app);
	double sc;
	long totnum;
	long num;

	/*ofs <<"train" <<std::endl;
	for(long l=0;l<100;l++) {
		double lthhalf=trainondis[l]*trainondis[l]*alpha;
		NBSQuery query;
		for(int i=0;i<ndim;i++) {
			query.vals.push_back(traindata[l][i]);
			//query.lcuts.push_back(traindata[l][i]-lthhalf*2.0);
			//query.rcuts.push_back(traindata[l][i]+lthhalf*2.0);
			//double lthhalf=traindata[l][i]*traindata[l][i]*alpha;
			query.lcuts.push_back(lthhalf*2.0);
			query.rcuts.push_back(lthhalf*2.0);
		}
		sc=trainscore(query,traintree,traindata,totnum,num,false,l);
		ofs <<totnum <<'\t' <<num <<'\t' <<std::fixed <<std::setprecision(6) <<sc <<std::endl;
	}
	ofs <<"ref" <<std::endl;*/
	ofs <<"alpha = " <<alpha <<std::endl;
	std::set<long> sls;
	NSPdstl::RandomEngine<>::getinstance().reseed(876);
	auto & rng = NSPdstl::RandomEngine<>::getinstance().intrng(0,refdata.size()-1);
	for(long li=0;li<100;li++) {
		long l=rng();
		while(sls.find(l)!=sls.end()) l=rng();
		sls.insert(l);
		double lthhalf=refondis[l]*refondis[l]*alpha;
		NBSQuery query;
		for(int i=0;i<ndim;i++) {
			query.vals.push_back(refdata[l][i]);
			//query.lcuts.push_back(refdata[l][i]-lthhalf*2.0);
			//query.rcuts.push_back(refdata[l][i]+lthhalf*2.0);
			//double lthhalf=refdata[l][i]*refdata[l][i]*alpha;
			query.lcuts.push_back(lthhalf*2.0);
			query.rcuts.push_back(lthhalf*2.0);
		}
		sc=refscore(query,reftree,refdata,totnum,num,false,l);
		ofs <<totnum <<'\t' <<num <<'\t' <<std::fixed <<std::setprecision(6) <<sc <<std::endl;
	}
	ofs.close();
}

void NSPproteinrep::out_crdene(std::string enefile, std::string crdfile,
		std::string outfile, long trainnum, long refnum, double minsprop) {
	std::map<long,double> enes;
	std::ifstream ifs(enefile);
	std::string line;
	double dtn=(double)trainnum;
	double drn=(double)refnum;
	while(std::getline(ifs,line)) {
		std::istringstream iss(line);
		long sx;
		iss >>sx;
		double rsc,tsc,ron;
		for(int i=0;i<5;i++) iss >>rsc;
		if(rsc<0.00001) continue;
		iss >>tsc;
		iss >>ron;
		double wgh=1.0;
		if(ron>4.0 && ron<=4.25) wgh = -8.0*(ron-4.0)*(ron-4.0)+1.0;
		else if(ron<4.5 && ron>4.25) wgh = (ron-4.5)*(ron-4.5)*8.0;
		else if(ron>=4.5) wgh = 0.0;
		double prob=((tsc*drn/rsc/dtn+minsprop)/(1+minsprop))*wgh+1.0-wgh;
		//if(prob<=0.0) std::cout <<tsc <<' ' <<rsc <<' ' <<wgh <<' ' <<prob <<std::endl;
		double e=-log(prob);
		enes.insert({sx,e});
	}
	ifs.close();
	std::ofstream ofs(outfile);
	ifs.open(crdfile);
	long lx=-1;
	while(std::getline(ifs,line)) {
		lx++;
		if(enes.find(lx)==enes.end()) continue;
		ofs <<line <<' ' <<std::fixed <<std::setprecision(6) <<enes.at(lx) <<std::endl;
	}
	ifs.close();
	ofs.close();
}

void NSPproteinrep::getpar(std::vector<XYZ> &crds, double &con, double &cno, double &ncon, double &ocno, double &tor, double &lth) {
	XYZ n1,c1,o1,n2,c2,o2;
	double dhb1=(crds[1]-crds[5]).squarednorm();
	double dhb2=(crds[2]-crds[4]).squarednorm();
	if(dhb1<dhb2) {
		n1=crds[2];
		c1=crds[0];
		o1=crds[1];
		n2=crds[5];
		c2=crds[3];
		o2=crds[4];
		lth=sqrt(dhb1);
	} else {
		n1=crds[5];
		c1=crds[3];
		o1=crds[4];
		n2=crds[2];
		c2=crds[0];
		o2=crds[1];
		lth=sqrt(dhb2);
	}
	double deg=180.0/3.14159265;
	con=angle(c1,o1,n2)*deg;
	cno=angle(o1,n2,c2)*deg;
	ncon=torsion(n1,c1,o1,n2)*deg;
	ocno=torsion(o2,c2,n2,o1)*deg;
	tor=torsion(c1,o1,n2,c2)*deg;
}

void NSPproteinrep::checksampletrain(std::string fnin, std::string fnout, int nsample) {
	std::string line;
	std::ifstream ifs(fnin);
	std::ofstream ofs(fnout);
	//ofs <<std::fixed <<std::setprecision(2);
	int ix=0;
	LocalBbHBGeoNNTerm nnmodel;
	std::vector<int> atomids = { 0, 1, 2, 3, 4, 5 };
	std::vector<DvDxi> dvdxi;
	while(std::getline(ifs,line)) {
		std::istringstream iss(line);
		std::vector<XYZ> vc(6);
		for(int i=0;i<vc.size();i++) {
			iss >>vc[i].x_;
			iss >>vc[i].y_;
			iss >>vc[i].z_;
			vc[i].x_ *= 0.1;
			vc[i].y_ *= 0.1;
			vc[i].z_ *= 0.1;
		}
		nnmodel.setup(vc, atomids);
		double ene = nnmodel.outvalue(&dvdxi);
		double con,cno,ncon,ocno,tor,lth;
		getpar(vc,con,cno,ncon,ocno,tor,lth);
		if(lth>0.45) continue;
		ofs <<lth*10.0 <<'\t' <<con <<'\t' <<cno <<'\t' <<ncon <<'\t' <<ocno <<'\t' <<ene <<std::endl;
		ix++;
		if(ix==nsample) break;
	}
	ifs.close();
	ofs.close();
}

void NSPproteinrep::checksampleref(std::string fnin, std::string fnout, int nstart, int nstep) {
	std::vector<std::vector<XYZ>> crds;
	std::string line;
	std::ifstream ifs(fnin);
	while(std::getline(ifs,line)) {
		std::istringstream iss(line);
		std::vector<XYZ> vc(6);
		for(int i=0;i<vc.size();i++) {
			iss >>vc[i].x_;
			iss >>vc[i].y_;
			iss >>vc[i].z_;
			vc[i].x_ *= 0.1;
			vc[i].y_ *= 0.1;
			vc[i].z_ *= 0.1;
		}
		crds.push_back(vc);
		if(crds.size()==nstart*20) break;
	}
	ifs.close();
	std::ofstream ofs(fnout);
	//ofs <<std::fixed <<std::setprecision(2);
	LocalBbHBGeoNNTerm nnmodel;
	std::vector<int> atomids = { 0, 1, 2, 3, 4, 5 };
	std::vector<DvDxi> dvdxi;
	double ene;
	//NSPdstl::RandomEngine<>::getinstance().reseed(seed);
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	//int nsamp=0;
	for(int nsamp=0;nsamp<nstart;nsamp++) {
		std::vector<XYZ> crs=crds[nsamp];
		XYZ trans0=crs[0];
		for(int i=0;i<crs.size();i++) crs[i]=crs[i]-trans0;
		std::vector<XYZ> fixedcrd(3),randomcrd(3);
		for(int i=0;i<3;i++) {
			fixedcrd[i]=crs[i];
			randomcrd[i]=crs[i+3];
		}
		bool clash=true;
		while(clash){
			randomchangebig(randomcrd,0.45,0);
			clash=false;
			for(int j=0;j<3;j++) {
						for(int k=0;k<3;k++) {
							double dist2=(fixedcrd[j]-randomcrd[k]).squarednorm();
							if(j==1&&k==2 || j==2&&k==1) {
								if(dist2<0.041) {
									clash=true;
									break;
								}
							} else if(dist2<0.06251) {
								clash=true;
								break;
							}
						}
			}
		}
		for(int i=0;i<3;i++) crs[i+3]=randomcrd[i];
		nnmodel.setup(crs, atomids);
		double eneprev = nnmodel.outvalue(&dvdxi);
		ene=eneprev;
		double con,cno,ncon,ocno,tor,lth;
		for(int i=0;i<nstep;i++) {
			std::vector<XYZ> rand1=crs;
		    randomcrd[0]=crs[3];
		    randomcrd[1]=crs[4];
		    randomcrd[2]=crs[5];
			randomchangesml(randomcrd,0.05,5.0);
			for(int j=3;j<6;j++) crs[j]=randomcrd[j-3];
			//double con,cno,ncon,ocno,tor,lth;
			getpar(crs,con,cno,ncon,ocno,tor,lth);
			if(lth>0.45) {
				//i--;
				crs=rand1;
				continue;
			}
			bool clash{false};
			for(int j=0;j<3;j++) {
				for(int k=0;k<3;k++) {
					double dist2=(fixedcrd[j]-randomcrd[k]).squarednorm();
					if(j==1&&k==2 || j==2&&k==1) {
						if(dist2<0.04) {
							clash=true;
							break;
						}
					} else if(dist2<0.0625) {
						clash=true;
						break;
					}
				}
				if(clash) break;
			}
			ene=0.0;
			if(clash) {
				//i--;
				crs=rand1;
				continue;
				//ene=10.0;
			}
			nnmodel.setup(crs, atomids);
			double enenn = nnmodel.outvalue(&dvdxi);
			ene +=enenn;
			double delta=ene-eneprev;
			bool accept{false};
			if(delta<0.0) accept=true;
			else {
				if(rng.realrng(0.0,1.0)()<exp(-delta)) accept=true;
			}
			if(!accept) {
				//i--;
				crs=rand1;
				continue;
			}
			eneprev=ene;
			//double con,cno,ncon,ocno,tor,lth;
			//getpar(crs,con,cno,ncon,ocno,tor,lth);
			//ofs <<lth*10.0 <<'\t' <<con <<'\t' <<cno <<'\t' <<ncon <<'\t' <<ocno <<'\t' <<ene <<std::endl;
		}
		if(lth>0.45) {
			nsamp--;
			continue;
		}
		ofs <<lth*10.0 <<'\t' <<con <<'\t' <<cno <<'\t' <<ncon <<'\t' <<ocno <<'\t' <<ene <<std::endl;
		if(nsamp%100==0) std::cout <<nsamp <<std::endl;
	}
	ofs.close();
}




void NSPproteinrep::transform_random_sidechain(std::string fnin, std::string fnout,
		std::string fnoutnoclash, int ntimes, double radius) {
	std::vector<std::vector<XYZ>> nativecrds;
	std::ifstream ifs(fnin);
	std::string readline;
	while(std::getline(ifs,readline)) {
		std::vector<XYZ> cs(6);
		std::istringstream iss(readline);
		for(int i=0;i<cs.size();i++) {
			iss >>cs[i].x_;
			iss >>cs[i].y_;
			iss >>cs[i].z_;
		}
		XYZ trans=cs[0];
		cs[0]=cs[0]-trans;
		cs[1]=cs[1]-trans;
		cs[2]=cs[2]-trans;
		nativecrds.push_back(cs);
	}
	ifs.close();

	int num=ntimes*nativecrds.size();
	int nbins=(int)((radius-dhb)/0.1);
	long nperbin=num/nbins+1;
	std::map<int,long> binrecord;
	int binstart=(int)(dhb*10.0);
	assert(binstart>19);
	for(int i=0;i<nbins;i++) binrecord.insert({binstart+i,nperbin});
	assert(nbins*nperbin>=num);

	std::ofstream ofs(fnout);
	std::ofstream nocl(fnoutnoclash);
	XYZ ori(0.0,0.0,0.0);
	long randomnum=-1;
	long acceptnum=0;
	long noclnum=0;
	while(!binrecord.empty()) {
		randomnum++;
		if(randomnum==nativecrds.size()) {
			std::cout <<randomnum <<' ' <<acceptnum <<' ' <<noclnum <<std::endl;
			randomnum=0;
		}
		std::vector<XYZ> fixedcrd(3), randomcrd(3);
		fixedcrd[0]=nativecrds[randomnum][0];
		fixedcrd[1]=nativecrds[randomnum][1];
		fixedcrd[2]=nativecrds[randomnum][2];
		randomcrd[0]=nativecrds[randomnum][3];
		randomcrd[1]=nativecrds[randomnum][4];
		randomcrd[2]=nativecrds[randomnum][5];
		randomchangebig(randomcrd,radius,0);
		int hblth=(int)(sqrt((randomcrd[0]-fixedcrd[0]).squarednorm())*10);
		if(binrecord.find(hblth)==binrecord.end()) continue;
		acceptnum++;
		binrecord.at(hblth)--;
		if(binrecord.at(hblth)==0) binrecord.erase(hblth);
		for(int j=0;j<3;j++) {
			ofs <<std::fixed <<std::setprecision(3) <<fixedcrd[j].x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<fixedcrd[j].y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<fixedcrd[j].z_ <<'\t';
		}
		for(int j=0;j<3;j++) {
			ofs <<std::fixed <<std::setprecision(3) <<randomcrd[j].x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<randomcrd[j].y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<randomcrd[j].z_ <<'\t';
		}
		ofs <<std::endl;
		bool clash{false};
		for(int j=0;j<3;j++) {
			for(int k=0;k<3;k++) {
				double dist2=(fixedcrd[j]-randomcrd[k]).squarednorm();
				if(j==0&&k==0) {
					if(dist2<dhb) {
						clash=true;
						break;
					}
				} else if(dist2<daa) {
					clash=true;
					break;
				}
			}
			if(clash) break;
		}
		if(clash) continue;
		noclnum++;
		for(int j=0;j<3;j++) {
			nocl <<std::fixed <<std::setprecision(3) <<fixedcrd[j].x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<fixedcrd[j].y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<fixedcrd[j].z_ <<'\t';
		}
		for(int j=0;j<3;j++) {
			nocl <<std::fixed <<std::setprecision(3) <<randomcrd[j].x_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<randomcrd[j].y_ <<'\t'
					<<std::fixed <<std::setprecision(3) <<randomcrd[j].z_ <<'\t';
		}
		nocl <<std::endl;
	}
	ofs.close();
	nocl.close();
	std::cout <<"total random pairs: " <<randomnum <<' ' <<acceptnum <<' ' <<noclnum <<std::endl;
}

void NSPproteinrep::scgetparameter(std::string fnin, std::string fnout) {
	std::ifstream ifs(fnin);
	std::ofstream ofs(fnout);
	std::string line;
	while(std::getline(ifs,line)) {
		std::istringstream iss(line);
		std::vector<XYZ> crds(6);
		for(int i=0;i<6;i++) {
			iss >>crds[i].x_ >>crds[i].y_ >>crds[i].z_;
		}
		double rhb=sqrt((crds[0]-crds[3]).squarednorm());
		double ad1=angle(crds[1],crds[0],crds[3]);
		double ad2=angle(crds[4],crds[3],crds[0]);
		double aaad1=torsion(crds[2],crds[1],crds[0],crds[3]);
		double aaad2=torsion(crds[5],crds[4],crds[3],crds[0]);
		double cost1=cos(ad1);
		double sint1=sin(ad1);
		double cosp1=cos(aaad1);
		double sinp1=sin(aaad1);
		double cost2=cos(ad2);
		double sint2=sin(ad2);
		double cosp2=cos(aaad2);
		double sinp2=sin(aaad2);
		ofs <<std::fixed <<std::setprecision(4) <<rhb*cost1 <<'\t'
				<<std::fixed <<std::setprecision(4) <<rhb*sint1*cosp1 <<'\t'
				<<std::fixed <<std::setprecision(4) <<rhb*sint1*sinp1 <<'\t'
				<<std::fixed <<std::setprecision(4) <<rhb*cost2 <<'\t'
				<<std::fixed <<std::setprecision(4) <<rhb*sint2*cosp2 <<'\t'
				<<std::fixed <<std::setprecision(4) <<rhb*sint2*sinp2 <<'\t'
				<<std::fixed <<std::setprecision(4) <<rhb <<std::endl;
	}
	ifs.close();
	ofs.close();
}










