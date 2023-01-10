/*
 * confbuilder.cpp
 *
 *  Created on: Dec 11, 2018
 *      Author: xuyang
 */
#include "allatom/confbuilder.h"
#include "dataio/controlfile.h"
#include "dataio/parameters.h"
#include "backbone/backbonesite.h"
#include "geometry/localframe.h"
#include "sd/sdrun.h"
#include "sd/genchain.h"
#include <algorithm>
using namespace NSPallatom;
using namespace NSPproteinrep;
using namespace NSPgeometry;

std::vector<BackBoneSite> InitialConfBuilder::buildchainC2N(BackBoneSite &st, double phi,
		const std::vector<std::pair<double,double>>& dihs) {
	std::vector<BackBoneSite> chain(dihs.size());
	st.data_[BackBoneSite::PHI]=phi;
	genprevbackbonesite(&st, 180.0, dihs[0].second, dihs[0].first, &chain[0]);
	for(int i=1;i<chain.size();i++) {
		genprevbackbonesite(&chain[i-1], 180.0, dihs[i].second, dihs[i].first, &chain[i]);
	}
	std::vector<BackBoneSite> bbss;
	for(int i=chain.size()-1;i>=0;i--) {
		bbss.push_back(chain[i]);
	}
	return bbss;
}

std::vector<BackBoneSite> InitialConfBuilder::buildchainN2C(BackBoneSite &st,
		const std::vector<std::pair<double,double>>& dihs) {
	std::vector<BackBoneSite> bbss(dihs.size());
	st.resetpsi();
	st.data_[BackBoneSite::OMIGA]=180.0;
	genbackbonesite(&st, false, dihs[0].first, dihs[0].second, &bbss[0]);
	for(int i=1;i<bbss.size();i++) {
		genbackbonesite(&bbss[i-1], false, dihs[i].first, dihs[i].second, &bbss[i]);
	}
	return bbss;
}

void InitialConfBuilder::rotphiN(std::vector<BackBoneSite> &chain, int n, double ang_add) {
	XYZ axis=chain[n].ncrd()-chain[n].cacrd();
	QuaternionCrd qc(axis,ang_add);
	Rotation rt(qc,chain[n].cacrd());
	for(int i=0;i<n;i++) {
		chain[i].rotate(rt);
	}
}

void InitialConfBuilder::rotpsiN(std::vector<BackBoneSite> &chain, int n, double ang_add) {
	XYZ axis=chain[n].cacrd()-chain[n].ccrd();
	QuaternionCrd qc(axis,ang_add);
	Rotation rt(qc,chain[n].cacrd());
	XYZ ncrd=chain[n].ncrd();
	rt.apply(&ncrd);
	chain[n].data_[BackBoneSite::NCRD] = ncrd.x_;
	chain[n].data_[BackBoneSite::NCRD+1] = ncrd.y_;
	chain[n].data_[BackBoneSite::NCRD+2] = ncrd.z_;
	for(int i=0;i<n;i++) {
		chain[i].rotate(rt);
	}
}

void InitialConfBuilder::rotphiC(std::vector<BackBoneSite> &chain, int n, double ang_add) {
	XYZ axis=chain[n].cacrd()-chain[n].ncrd();
	QuaternionCrd qc(axis,ang_add);
	Rotation rt(qc,chain[n].cacrd());
	XYZ ccrd=chain[n].ccrd();
	rt.apply(&ccrd);
	chain[n].data_[BackBoneSite::CCRD] = ccrd.x_;
	chain[n].data_[BackBoneSite::CCRD+1] = ccrd.y_;
	chain[n].data_[BackBoneSite::CCRD+2] = ccrd.z_;
	XYZ ocrd=chain[n].ocrd();
	rt.apply(&ocrd);
	chain[n].data_[BackBoneSite::OCRD] = ocrd.x_;
	chain[n].data_[BackBoneSite::OCRD+1] = ocrd.y_;
	chain[n].data_[BackBoneSite::OCRD+2] = ocrd.z_;
	for(int i=n+1;i<chain.size();i++) {
		chain[i].rotate(rt);
	}
}

void InitialConfBuilder::rotpsiC(std::vector<BackBoneSite> &chain, int n, double ang_add) {
	XYZ axis=chain[n].ccrd()-chain[n].cacrd();
	QuaternionCrd qc(axis,ang_add);
	Rotation rt(qc,chain[n].cacrd());
	XYZ ocrd=chain[n].ocrd();
	rt.apply(&ocrd);
	chain[n].data_[BackBoneSite::OCRD] = ocrd.x_;
	chain[n].data_[BackBoneSite::OCRD+1] = ocrd.y_;
	chain[n].data_[BackBoneSite::OCRD+2] = ocrd.z_;
	for(int i=n+1;i<chain.size();i++) {
		chain[i].rotate(rt);
	}
}

void InitialConfBuilder::changephipsi(std::vector<BackBoneSite> &chain, int n, double ang_add, bool isphi, bool movedisn) {
	if(isphi) {
		if(movedisn) rotphiN(chain,n,ang_add);
		else rotphiC(chain,n,ang_add);
	} else {
		if(movedisn) rotpsiN(chain,n,ang_add);
		else rotpsiC(chain,n,ang_add);
	}
}

void InitialConfBuilder::unbendstrand_MC(std::vector<BackBoneSite> &chain, std::vector<std::pair<double,double>> &dihs, bool firstisfixed) {
	std::vector<std::pair<double,std::pair<double,double>>> das;
	for(auto &d:dihs) {
		das.push_back({d.first,strandphirange_});
		das.push_back({d.second,strandpsirange_});
	}
	std::vector<int> nang;
	if(firstisfixed) {
		for(int i=2;i<das.size()-1;i++) nang.push_back(i);
	} else {
		for(int i=1;i<das.size()-1;i++) nang.push_back(i);
	}
	auto rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	auto rng_i=NSPdstl::RandomEngine<>::getinstance().intrng(0,nang.size()-1);
	auto rng_i2=NSPdstl::RandomEngine<>::getinstance().intrng(0,nang.size()/2);
	auto rng_i01=NSPdstl::RandomEngine<>::getinstance().intrng(0,1);
	double dis=sqrt((unbenddir_.second-(firstisfixed?chain.back():chain[0]).cacrd()).squarednorm());
	for(int i=0;i<unbendstep_;i++) {
		//if(i%100==0) std::cout <<"UnBending Step: " <<i <<std::endl;
		std::vector<std::pair<int,double>> changed(rng_i01()==0?1:rng_i2());
		//assign changed
		std::set<int> lb;
		while(lb.size()<changed.size()) lb.insert(nang[rng_i()]);
		int nassign=0;
		for(int j:lb) {
			double a;
			do {
				a=rng()*steplength_*2.0-steplength_;
			} while(das[j].first+a<das[j].second.first || das[j].first+a>das[j].second.second);
			changed[nassign]={j,a};
			nassign++;
		}
		//assign changed end
		std::vector<BackBoneSite> chain1=chain;
		//change chain1
		for(int j=0;j<changed.size();j++) {
			changephipsi(chain1,changed[j].first/2,changed[j].second,changed[j].first%2==0,!firstisfixed);
		}
		//change chain1 end
		double dis1=sqrt((unbenddir_.second-(firstisfixed?chain1.back():chain1[0]).cacrd()).squarednorm());
		//assign dis and compare with the former
		bool accept=(dis1>dis || rng()<exp(dis1-dis));
		if(accept) {
			chain=chain1;
			dis=dis1;
			for(auto &ce:changed) das[ce.first].first+=ce.second;
		}
		//std::cout <<i <<'\t' <<changed.size() <<'\t' <<dis1 <<'\t' <<dis <<'\t' <<accept <<std::endl;
	}
	std::cout <<"UnBending Result: " <<std::endl;
	for(int i=0;i<das.size();i+=2) {
		std::cout <<'\t' <<das[i].first;
	}std::cout <<std::endl;
	for(int i=1;i<das.size();i+=2) {
		std::cout <<'\t' <<das[i].first;
	}std::cout <<std::endl;
}

std::vector<Residue> InitialConfBuilder::buildfirststrand(BackBoneSite &st) {
	int ix=extendpair_[unbenddir_.first];
	bool c2n;
	if(ix==0) {
		c2n=true;
		st=design_[ix+1].second[0].getbackbonesite();
	} else {
		if(design_[ix-1].first==Fixed) {
			c2n=false;
			st=design_[ix-1].second.back().getbackbonesite();
		} else {
			c2n=true;
			st=design_[ix+1].second[0].getbackbonesite();
		}
	}
	std::vector<BackBoneSite> chain;
	double defaultphi=-130.0;
	double defaultpsi=120.0;
	std::vector<std::pair<double,double>> dihs(length_[ix],{defaultphi,defaultpsi});
	std::vector<Residue> chs;
	if(c2n) {
		chain=buildchainC2N(st,defaultphi,dihs);
		chain.push_back(st);
		dihs.push_back({defaultphi,st.data_[BackBoneSite::PSI]});
		unbendstrand_MC(chain,dihs,false);
		for(int i=0;i<chain.size()-1;i++) {
			chain[i].resname="GLY";
			chs.push_back(Residue(chain[i]));
		}
	} else {
		std::vector<BackBoneSite> chain1=buildchainN2C(st,dihs);
		chain.push_back(st);
		for(auto &c:chain1) chain.push_back(c);
		dihs.resize(dihs.size()+1);
		for(int i=dihs.size()-1;i>0;i--) {
			dihs[i]=dihs[i-1];
		}
		st.resetpsi();
		dihs[0]={st.data_[BackBoneSite::PHI],st.data_[BackBoneSite::PSI]};
		unbendstrand_MC(chain,dihs,true);
		for(int i=1;i<chain.size();i++) {
			chain[i].resname="GLY";
			chs.push_back(Residue(chain[i]));
		}
	}
	return chs;
}

//build strand chain paired with chain1
//st paired with chain1[0] or chain1.back()
std::vector<BackBoneSite> InitialConfBuilder::pairedstrandN2C(BackBoneSite &st, const std::vector<BackBoneSite>&chain1, int lth) {
	double hbdis=2.9;
	st.resetpsi();
	st.data_[BackBoneSite::OMIGA]=180.0;
	std::vector<BackBoneSite> chain;
	XYZ ca0=st.cacrd(), ca1=chain1[0].cacrd(), ca2=chain1.back().cacrd();
	bool isparallel=(ca0-ca1).squarednorm()<(ca0-ca2).squarednorm();
	std::set<int> inhb;
	if(isparallel) {
		XYZ o1=st.ocrd(), n1=st.getnextN(), o2=chain1[0].ocrd(), n2=chain1[1].ncrd();
		if((n1-o2).squarednorm()<(n2-o1).squarednorm()) {
			for(int i=0;i<lth;i++) {
				if(i%2==0) inhb.insert(i);
			}
		} else {
			for(int i=0;i<lth;i++) {
				if(i%2==1) inhb.insert(i);
			}
		}
	} else {
		XYZ o1=st.ocrd(), n1=st.getnextN(), o2=chain1[chain1.size()-2].ocrd(), n2=chain1.back().ncrd();
		if((o1-n2).squarednorm()<(o2-n1).squarednorm()) {
			for(int i=0;i<lth;i++) {
				if(i%2==1) inhb.insert(i);
			}
		} else {
			for(int i=0;i<lth;i++) {
				if(i%2==0) inhb.insert(i);
			}
		}
	}
	for(int i=0;i<lth;i++) {
		std::cout <<"HB Paired Chain Building: " <<i <<std::endl;
		chain.push_back(BackBoneSite());
		double mh,ms;
		double mdis{10000.0};
		for(double h=strandphirange_.first;h<strandphirange_.second;h+=prec_) {
			for(double s=strandpsirange_.first;s<strandpsirange_.second;s+=prec_) {
				if(i==0) genbackbonesite(&st, false, mh, ms, &chain[i]);
				else genbackbonesite(&chain[i-1], false, mh, ms, &chain[i]);
				double delta;
				if(isparallel) {
					if(inhb.find(i)==inhb.end()) delta=fabs(sqrt((chain[i].getnextN()-chain1[i-1].ocrd()).squarednorm())-hbdis);
					else delta=fabs(sqrt((chain[i].ocrd()-chain1[i].ncrd()).squarednorm())-hbdis);
				} else {
					if(inhb.find(i)==inhb.end()) delta=fabs(sqrt((chain[i].getnextN()-
							chain1[chain1.size()-2-i].ocrd()).squarednorm())-hbdis);
					else delta=fabs(sqrt((chain[i].ocrd()-chain1[chain1.size()-i-2].ncrd()).squarednorm())-hbdis);
				}
				if(delta<mdis) {
					mh=h;
					ms=s;
					mdis=delta;
				}
			}
		}
		if(i==0) genbackbonesite(&st, false, mh, ms, &chain[i]);
		else genbackbonesite(&chain[i-1], false, mh, ms, &chain[i]);
		chain[i].resetpsi();
		chain[i].data_[BackBoneSite::OMIGA]=180.0;
	}
	return chain;
}

std::vector<BackBoneSite> InitialConfBuilder::pairedstrandC2N(BackBoneSite &st, const std::vector<BackBoneSite>&chain1, int lth) {
	double hbdis=2.9;
	std::vector<BackBoneSite> chain;
	XYZ ca0=st.cacrd(), ca1=chain1[0].cacrd(), ca2=chain1.back().cacrd();
	bool isparallel=(ca0-ca1).squarednorm()>(ca0-ca2).squarednorm();
	std::set<int> inhb;
	if(isparallel) {
		XYZ n1=st.ncrd(), o2=chain1[chain1.size()-2].ocrd();
		if((n1-o2).squarednorm()<3.5*3.5) {
			for(int i=lth-1;i>=0;i--) {
				if((lth-i)%2==0) inhb.insert(i);
			}
		} else {
			for(int i=lth-1;i>=0;i--) {
				if((lth-i)%2==1) inhb.insert(i);
			}
		}
	} else {
		XYZ n1=st.ncrd(), o2=chain1[chain1.size()-1].ocrd();
		if((n1-o2).squarednorm()<3.5*3.5) {
			for(int i=lth-1;i>=0;i--) {
				if((lth-i)%2==0) inhb.insert(i);
			}
		} else {
			for(int i=lth-1;i>=0;i--) {
				if((lth-i)%2==1) inhb.insert(i);
			}
		}
	}
	for(int i=0;i<lth;i++) {
		std::cout <<"HB Paired Chain Building: " <<i <<std::endl;
		chain.push_back(BackBoneSite());
		for(int i=chain.size()-1;i>0;i--) chain[i]=chain[i-1];
		double mh,ms;
		double mdis{10000.0};
		for(double h=strandphirange_.first;h<strandphirange_.second;h+=prec_) {
			for(double s=strandpsirange_.first;s<strandpsirange_.second;s+=prec_) {
				if(i==0) genprevbackbonesite_phiislast(&st, 180.0, ms, mh, &chain[0]);
				else genprevbackbonesite_phiislast(&chain[1], 180.0, ms, mh, &chain[0]);
				double delta;
				if(isparallel) {
					if(inhb.find(i)==inhb.end()) delta=fabs(sqrt((chain[0].getprevO()-chain1[chain1.size()-1-i].ncrd()).squarednorm())-hbdis);
					else delta=fabs(sqrt((chain[0].ncrd()-chain1[chain1.size()-3-i].ocrd()).squarednorm())-hbdis);
				} else {
					if(inhb.find(i)==inhb.end()) delta=fabs(sqrt((chain[0].getprevO()-chain1[chain1.size()-2-i].ncrd()).squarednorm())-hbdis);
					else delta=fabs(sqrt((chain[0].ncrd()-chain1[chain1.size()-i-2].ocrd()).squarednorm())-hbdis);
				}
				if(delta<mdis) {
					mh=h;
					ms=s;
					mdis=delta;
				}
			}
		}
		if(i==0) genprevbackbonesite_phiislast(&st, 180.0, ms, mh, &chain[0]);
		else genprevbackbonesite_phiislast(&chain[1], 180.0, ms, mh, &chain[0]);
	}
	return chain;
}

std::vector<Residue> InitialConfBuilder::buildsecondstrand(const std::vector<Residue>&firstchain, const BackBoneSite &firstst) {
	std::vector<BackBoneSite> chf;
	XYZ ca0=firstst.cacrd();
	XYZ ca1=firstchain[0].rds().at("CA").crd;
	XYZ ca2=firstchain.back().rds().at("CA").crd;
	if((ca0-ca1).squarednorm()<(ca0-ca2).squarednorm()) {
		chf.push_back(firstst);
		for(int i=0;i<firstchain.size();i++) chf.push_back(firstchain[i].getbackbonesite());
	} else {
		for(int i=0;i<firstchain.size();i++) chf.push_back(firstchain[i].getbackbonesite());
		chf.push_back(firstst);
	}
	int ix=extendpair_[unbenddir_.first==0?1:0];
	bool c2n;
	BackBoneSite st;
	if(ix==0) {
		c2n=true;
		st=design_[ix+1].second[0].getbackbonesite();
	} else {
		if(design_[ix-1].first==Fixed) {
			c2n=false;
			st=design_[ix-1].second.back().getbackbonesite();
		} else {
			c2n=true;
			st=design_[ix+1].second[0].getbackbonesite();
		}
	}

	std::vector<BackBoneSite> chain;
	std::vector<Residue> chs;
	if(c2n) {
		chain=pairedstrandC2N(st,chf,length_[ix]);
		for(int i=0;i<chain.size();i++) {
			chain[i].resname="GLY";
			chs.push_back(Residue(chain[i]));
		}
	} else {
		chain=pairedstrandC2N(st,chf,length_[ix]);
		for(int i=0;i<chain.size();i++) {
			chain[i].resname="GLY";
			chs.push_back(Residue(chain[i]));
		}
	}
	return chs;
}

void InitialConfBuilder::hbpaired_MC(BackBoneSite &st1, BackBoneSite &st2, BackBoneSite &bs1, BackBoneSite &bs2,
		bool firstisn2c, bool firstisinhb, bool parrallel, std::vector<double> &stdih) {
	int nstep{mcstep_};
	double dihrange{mcrange_};
	double steplength{mcsteplth_};
	double hbdis{2.9};
	std::vector<std::pair<double,std::pair<double,double>>> dih{{stdih[0],{stdih[0]-dihrange,stdih[0]+dihrange}},
		{stdih[1],{stdih[1]-dihrange,stdih[1]+dihrange}},{stdih[2],{stdih[2]-dihrange,stdih[2]+dihrange}},
		{stdih[3],{stdih[3]-dihrange,stdih[3]+dihrange}}};
	auto rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	auto rng_i=NSPdstl::RandomEngine<>::getinstance().intrng(0,dih.size()-1);
	auto rng_i01=NSPdstl::RandomEngine<>::getinstance().intrng(0,1);
	//XYZ po1=bs1.getprevO(), n1=bs1.ncrd(), o1=bs1.ocrd(), nn1=bs1.getnextN();
	//XYZ po2=bs2.getprevO(), n2=bs2.ncrd(), o2=bs2.ocrd(), nn2=bs2.getnextN();
	double dis=1000.0;
	for(int i=0;i<nstep;i++) {
		int nchanged=rng_i01()==0?1:rng_i();
		std::vector<double> changed(4,0.0);
		std::set<int> lb;
		while(lb.size()<nchanged) lb.insert(rng_i());
		int nassign=0;
		for(int j:lb) {
			double a;
			do {
				a=rng()*steplength*2.0-steplength;
			} while(dih[j].first+a<dih[j].second.first || dih[j].first+a>dih[j].second.second);
			changed[j]=a;
			nassign++;
		}
		BackBoneSite b1, b2;
		double dis1;
		if(firstisn2c) {
			genbackbonesite(&st1, false, dih[0].first+changed[0], dih[1].first+changed[1], &b1);
			if(parrallel) {
				genbackbonesite(&st2, false, dih[2].first+changed[2], dih[3].first+changed[3], &b2);
				if(firstisinhb) {
					dis1=fabs(sqrt((b1.getnextN()-b2.ocrd()).squarednorm())-hbdis);
				} else {
					dis1=fabs(sqrt((b1.ocrd()-b2.getnextN()).squarednorm())-hbdis);
				}
			} else {
				//genprevbackbonesite_phiislast(&st2, 180.0, dih[3].first+changed[3], dih[2].first+changed[2], &b2);
				genprevbackbonesite(&st2, 180.0, dih[3].first+changed[3], dih[2].first+changed[2], &b2);
				if(firstisinhb) {
					dis1=fabs(sqrt((b1.getnextN()-b2.getprevO()).squarednorm())-hbdis);
				} else {
					dis1=fabs(sqrt((b1.ocrd()-b2.ncrd()).squarednorm())-hbdis);
				}
			}
		} else {
			//genprevbackbonesite_phiislast(&st1, 180.0, dih[1].first+changed[1], dih[0].first+changed[0], &b1);
			genprevbackbonesite(&st1, 180.0, dih[1].first+changed[1], dih[0].first+changed[0], &b1);
			if(parrallel) {
				//genprevbackbonesite_phiislast(&st2, 180.0, dih[3].first+changed[3], dih[2].first+changed[2], &b2);
				genprevbackbonesite(&st2, 180.0, dih[3].first+changed[3], dih[2].first+changed[2], &b2);
				if(firstisinhb) {
					dis1=fabs(sqrt((b1.getprevO()-b2.ncrd()).squarednorm())-hbdis);
				} else {
					dis1=fabs(sqrt((b1.ncrd()-b2.getprevO()).squarednorm())-hbdis);
				}
			} else {
				genbackbonesite(&st2, false, dih[2].first+changed[2], dih[3].first+changed[3], &b2);
				if(firstisinhb) {
					dis1=fabs(sqrt((b1.getprevO()-b2.getnextN()).squarednorm())-hbdis);
				} else {
					dis1=fabs(sqrt((b1.ncrd()-b2.ocrd()).squarednorm())-hbdis);
				}
			}
		}
		//std::cout <<dis1 <<'\t' <<dih[0].first+changed[0] <<'\t' <<dih[1].first+changed[1]
		//		<<'\t' <<dih[2].first+changed[2] <<'\t' <<dih[3].first+changed[3] <<std::endl;
		//bool accept=(dis1<dis || rng()<exp(dis-dis1)/10.0);
		bool accept=(dis1<dis);
		if(accept) {
			bs1=b1;
			bs2=b2;
			dis=dis1;
			for(int i=0;i<changed.size();i++) dih[i].first += changed[i];
			//std::cout <<dis1 <<'\t' <<dih[0].first <<'\t' <<dih[1].first <<'\t' <<dih[2].first <<'\t' <<dih[3].first <<std::endl;
		}
	}
	//std::cout <<std::endl;
	//std::cout <<dih[0].first <<'\t' <<dih[1].first <<'\t' <<dih[2].first <<'\t' <<dih[3].first <<std::endl;
}



/*
void InitialConfBuilder::hbinextendpair(std::map<std::pair<int,int>,std::pair<int,int>> &hbn2o, bool &ex0stisinhb, bool &ex1stisinhb) {
	bool extend0isn2c{!(extendpair_[0]!=0 && design_[extendpair_[0]-1].first==Fixed)};
	bool extend1isn2c{!(extendpair_[1]!=0 && design_[extendpair_[1]-1].first==Fixed)};
	assert(extend1isn2c==(betacontanctpair_.at(extendpair_[0]).at(extendpair_[1])=='p'?extend0isn2c:(!extend0isn2c)));
	BackBoneSite st0{(!extend0isn2c)?design_[extendpair_[0]-1].second.back().getbackbonesite():
			design_[extendpair_[0]+1].second[0].getbackbonesite()};
	BackBoneSite st1{(!extend1isn2c)?design_[extendpair_[1]-1].second.back().getbackbonesite():
			design_[extendpair_[1]+1].second[0].getbackbonesite()};
	int l0=length_[extendpair_[0]], l1=length_[extendpair_[1]];
	assert(l0==l1);
	if(extend0isn2c) {
		if(extend1isn2c) {
			XYZ n0=st0.ncrd(), ca0=st0.cacrd(), n1=st1.ncrd(), ca1=st1.cacrd();
			bool st0isinhb{(ca0-n1).squarednorm()>(ca1-n0).squarednorm()};
			if(st0isinhb) {
				hbn2o.insert({{extendpair_[0]+1,0},{extendpair_[1],l1-1}});
				for(int i=0;i<l0-1;i++) {
					if(i%2==0) {
						hbn2o.insert({{extendpair_[1],l1-1-i},{extendpair_[0],l0-2-i}});
					} else {
						hbn2o.insert({{extendpair_[0],l0-1-i},{extendpair_[1],l1-2-i}});
					}
				}
			} else {
				hbn2o.insert({{extendpair_[1]+1,0},{extendpair_[0],l0-1}});
				for(int i=0;i<l0-1;i++) {
					if(i%2==0) {
						hbn2o.insert({{extendpair_[0],l0-1-i},{extendpair_[1],l1-2-i}});
					} else {
						hbn2o.insert({{extendpair_[1],l1-1-i},{extendpair_[0],l1-2-i}});
					}
				}
			}
			if(hbn2o.find({extendpair_[0],1})==hbn2o.end()) {
				ex0stisinhb=true;
				ex1stisinhb=false;
			} else {
				ex0stisinhb=false;
				ex1stisinhb=true;
			}
		} else {
			XYZ n0=st0.ncrd(), o1=st1.ocrd();
			bool st0isinhb{(n0-o1).squarednorm()<3.3*3.3};
			if(st0isinhb) {
				for(int i=1;i<l1;i+=2) {
					hbn2o.insert({{extendpair_[1],i},{extendpair_[0],l0-1-i}});
					hbn2o.insert({{extendpair_[0],l0-1-i},{extendpair_[1],i}});
				}
			} else {
				for(int i=0;i<l1;i+=2) {
					hbn2o.insert({{extendpair_[1],i},{extendpair_[0],l0-1-i}});
					hbn2o.insert({{extendpair_[0],l0-1-i},{extendpair_[1],i}});
				}
			}
			if(hbn2o.find({extendpair_[1],0})==hbn2o.end()) {
				ex0stisinhb=false;
				ex1stisinhb=false;
			} else {
				ex0stisinhb=true;
				ex1stisinhb=true;
			}
		}
	} else {
		if(extend1isn2c) {
			XYZ n0=st0.ncrd(), o1=st1.ocrd();
			bool st0isinhb{(n0-o1).squarednorm()<3.3*3.3};
			if(st0isinhb) {
				for(int i=1;i<l0;i+=2) {
					hbn2o.insert({{extendpair_[0],i},{extendpair_[1],l0-1-i}});
					hbn2o.insert({{extendpair_[1],l0-1-i},{extendpair_[0],i}});
				}
			} else {
				for(int i=0;i<l0;i+=2) {
					hbn2o.insert({{extendpair_[0],i},{extendpair_[1],l0-1-i}});
					hbn2o.insert({{extendpair_[1],l0-1-i},{extendpair_[0],i}});
				}
			}
			if(hbn2o.find({extendpair_[0],0})==hbn2o.end()) {
				ex0stisinhb=false;
				ex1stisinhb=false;
			} else {
				ex0stisinhb=true;
				ex1stisinhb=true;
			}
		} else {
			XYZ n0=st0.ncrd(), ca0=st0.cacrd(), n1=st1.ncrd(), ca1=st1.cacrd();
			bool st0isinhb{(ca0-n1).squarednorm()>(ca1-n0).squarednorm()};
			if(st0isinhb) {
				hbn2o.insert({{extendpair_[1],0},{extendpair_[0]-1,length_[extendpair_[0]-1]-1}});
				for(int i=0;i<l0-1;i++) {
					if(i%2==0) {
						hbn2o.insert({{extendpair_[0],i+1},{extendpair_[1],i}});
					} else {
						hbn2o.insert({{extendpair_[1],i+1},{extendpair_[0],i}});
					}
				}
			} else {
				hbn2o.insert({{extendpair_[0],0},{extendpair_[1]-1,length_[extendpair_[1]-1]-1}});
				for(int i=0;i<l0-1;i++) {
					if(i%2==0) {
						hbn2o.insert({{extendpair_[1],i+1},{extendpair_[0],i}});
					} else {
						hbn2o.insert({{extendpair_[0],i+1},{extendpair_[1],i}});
					}
				}
			}
			if(hbn2o.find({extendpair_[0],l0-1})==hbn2o.end()) {
				ex0stisinhb=false;
				ex1stisinhb=true;
			} else {
				ex0stisinhb=true;
				ex1stisinhb=false;
			}
		}
	}
}

void InitialConfBuilder::hbinstrandpair(std::map<std::pair<int,int>,std::pair<int,int>> &hbn2o,
		int l0, int l1, bool l0isn2c, bool l0stisinhb, bool l1isn2c, int i0, int i1) {
	int l2=l0>l1?l0:l1;
	if(l0isn2c) {
		if(l1isn2c) {
			if(l0stisinhb) {
				for(int i=0;i<l2;i++) {
					if(i%2==0) {
						if(i+1>=l1 || i>=l0) break;
						hbn2o.insert({{i1,i+1},{i0,i}});
					} else {
						if(i+1>=l0 || i>=l1) break;
						hbn2o.insert({{i0,i+1},{i1,i}});
					}
				}
			} else {
				for(int i=0;i<l2;i++) {
					if(i%2==0) {
						if(i+1>=l0 || i>=l1) break;
						hbn2o.insert({{i0,i+1},{i1,i}});
					} else {
						if(i+1>=l1 || i>=l0) break;
						hbn2o.insert({{i1,i+1},{i0,i}});
					}
				}
			}
		} else {
			int ix{l0stisinhb?0:1};
			for(int i=ix;i<l2;i+=2) {
				if(i>=l0 || i+1<l1) break;
				hbn2o.insert({{i0,i},{i1,l1-1-i}});
				hbn2o.insert({{i1,l1-1-i},{i0,i}});
			}
		}
	} else {
		if(l1isn2c) {
			int ix{l0stisinhb?0:1};
			for(int i=ix;i<l2;i+=2) {
				if(i>=l0 || i+1<l1) break;
				hbn2o.insert({{i0,l0-1-i},{i1,i}});
				hbn2o.insert({{i1,i},{i0,l0-1-i}});
			}
		} else {
			if(l0stisinhb) {
				for(int i=0;i<l2;i++) {
					if(i%2==0) {
						if(i+1>l0 || i+2>l1) break;
						hbn2o.insert({{i0,l0-1-i},{i1,l1-2-i}});
					} else {
						if(i+1>l1 || i+2>l0) break;
						hbn2o.insert({{i1,l1-1-i},{i0,l0-2-i}});
					}
				}
			} else {
				for(int i=0;i<l2;i++) {
					if(i%2==1) {
						if(i+1>l0 || i+2>l1) break;
						hbn2o.insert({{i0,l0-1-i},{i1,l1-2-i}});
					} else {
						if(i+1>l1 || i+2>l0) break;
						hbn2o.insert({{i1,l1-1-i},{i0,l0-2-i}});
					}
				}
			}
		}
	}
}

void InitialConfBuilder::hbnetinsheet() {
	std::map<std::pair<int,int>,std::pair<int,int>> hbn2o;
	bool ex0stisinhb, ex1stisinhb;
	hbinextendpair(hbn2o,ex0stisinhb,ex1stisinhb);
	int extend1;
	for(int i=0;i<seqinsheet_.size();i++) {
		if(seqinsheet_[i]!=extendpair_[1]) continue;
		extend1=i;
		break;
	}
	std::vector<bool> isn2c;
	bool l0stisinhb;
	for(int i=extend1-1;i>=0;i--) {
		hbinstrandpair(hbn2o, length_[seqinsheet_[i]], length_[seqinsheet_[i+1]], isn2c[i],
				l0stisinhb, isn2c[i+1], seqinsheet_[i], seqinsheet_[i+1]);
	}
}
*/
std::vector<double> InitialConfBuilder::hbpaired_ergodic(BackBoneSite &st1, BackBoneSite &st2, BackBoneSite &bs1, BackBoneSite &bs2,
		bool firstisn2c, bool firstisinhb, bool parrallel) {
	double prec{erglth_};
	double hbdis{hblth_};
	double mh1{0.0},ms1{0.0},mh2{0.0},ms2{0.0};
	std::vector<double> vdtemp{mh1,ms1,mh2,ms2};
	int nstore=1;
	std::vector<std::pair<double,std::pair<std::vector<double>,std::vector<BackBoneSite>>>> stores;
	for(int i=0;i<nstore;i++) {
		stores.push_back({1000.0,{vdtemp,std::vector<BackBoneSite>()}});
	}
	double dis{1000.0};
	BackBoneSite b1, b2;
	double dis1;
	for(double h1=strandphirange_.first;h1<strandphirange_.second;h1+=prec) {
		for(double s1=strandpsirange_.first;s1<strandpsirange_.second;s1+=prec) {
			if(firstisn2c) genbackbonesite(&st1, false, h1, s1, &b1);
			//else genprevbackbonesite_phiislast(&st1, 180.0, s1, h1, &b1);
			else genprevbackbonesite(&st1, 180.0, s1, h1, &b1);
			for(double h2=strandphirange_.first;h2<strandphirange_.second;h2+=prec) {
				for(double s2=strandpsirange_.first;s2<strandpsirange_.second;s2+=prec) {
					if(firstisn2c) {
						//genbackbonesite(&st1, false, h1, s1, &b1);
						if(parrallel) {
							genbackbonesite(&st2, false, h2, s2, &b2);
							if(firstisinhb) {
								dis1=fabs(sqrt((b1.getnextN()-b2.ocrd()).squarednorm())-hbdis);
							} else {
								dis1=fabs(sqrt((b1.ocrd()-b2.getnextN()).squarednorm())-hbdis);
							}
						} else {
							//genprevbackbonesite_phiislast(&st2, 180.0, s2, h2, &b2);
							genprevbackbonesite(&st2, 180.0, s2, h2, &b2);
							if(firstisinhb) {
								dis1=fabs(sqrt((b1.getnextN()-b2.getprevO()).squarednorm())-hbdis);
							} else {
								dis1=fabs(sqrt((b1.ocrd()-b2.ncrd()).squarednorm())-hbdis);
							}
						}
					} else {
						//genprevbackbonesite_phiislast(&st1, 180.0, s1, h1, &b1);
						if(parrallel) {
							//genprevbackbonesite_phiislast(&st2, 180.0, s2, h2, &b2);
							genprevbackbonesite(&st2, 180.0, s2, h2, &b2);
							if(firstisinhb) {
								dis1=fabs(sqrt((b1.getprevO()-b2.ncrd()).squarednorm())-hbdis);
							} else {
								dis1=fabs(sqrt((b1.ncrd()-b2.getprevO()).squarednorm())-hbdis);
							}
						} else {
							genbackbonesite(&st2, false, h2, s2, &b2);
							if(firstisinhb) {
								dis1=fabs(sqrt((b1.getprevO()-b2.getnextN()).squarednorm())-hbdis);
							} else {
								dis1=fabs(sqrt((b1.ncrd()-b2.ocrd()).squarednorm())-hbdis);
							}
						}
					}
					if(dis1<stores.back().first) {
						std::vector<double> vd{h1,s1,h2,s2};
						std::vector<BackBoneSite> bstemp{b1,b2};
						stores.back().first=dis1;
						stores.back().second.first=vd;
						stores.back().second.second=bstemp;
						std::sort(stores.begin(),stores.end(),[](
								std::pair<double,std::pair<std::vector<double>,std::vector<BackBoneSite>>> p1,
								std::pair<double,std::pair<std::vector<double>,std::vector<BackBoneSite>>> p2)->bool{
							return p1.first<p2.first;});
					}
				}
			}
		}
	}
	//std::cout <<"Result: " <<std::endl;
	//std::cout <<dis <<'\t' <<mh1 <<'\t' <<ms1 <<'\t' <<mh2 <<'\t' <<ms2 <<std::endl;
	int nchose=NSPdstl::RandomEngine<>::getinstance().intrng(0,nstore-1)();
	bs1=stores[nchose].second.second[0];
	bs2=stores[nchose].second.second[1];
	return stores[nchose].second.first;
}

std::vector<int> InitialConfBuilder::hbpaired_SD(std::string stfile, std::string controlfile,
		std::vector<std::pair<int,int>> &hb01no_ori, std::vector<std::pair<int,int>> &hb01on_ori,
		std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01no,
		std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01on) {
	std::vector<std::vector<Residue>> stchs;
	std::vector<Residue> chtemp;
	std::vector<std::vector<int>> activeparts(extendpair_.size());
	std::vector<int> acts(4);
	for(int i=0;i<design_.size();i++) {
		if(design_[i].second.empty()) {
			if(chtemp.empty()) continue;
			stchs.push_back(chtemp);
			chtemp.clear();
			continue;
		}
		int stnum=chtemp.size();
		for(auto &r:design_[i].second) chtemp.push_back(r);
		int ennum=chtemp.size();
		std::vector<int> vd{stchs.size(),stnum,ennum};
		if(i==extendpair_[0]) {
			activeparts[0]=vd;
			acts[0]=stchs.size();
			acts[1]=stnum;
		}
		if(i==extendpair_[1]) {
			activeparts[1]=vd;
			acts[2]=stchs.size();
			acts[3]=stnum;
		}
	}
	if(!chtemp.empty()) stchs.push_back(chtemp);
	for(auto &f:fixed_) stchs.push_back(f);
	for(int i=0;i<hb01no_ori.size();i++) {
		int c1=activeparts[0][0];
		int n1=activeparts[0][1]+hb01no_ori[i].first;
		int c2=activeparts[1][0];
		int n2=activeparts[1][1]+hb01no_ori[i].second;
		if(n1<0 || n1>=stchs[c1].size()) continue;
		if(n2<0 || n2>=stchs[c2].size()) continue;
		hb01no.push_back({{c1,n1},{c2,n2}});
	}
	for(int i=0;i<hb01on_ori.size();i++) {
		int c1=activeparts[0][0];
		int n1=activeparts[0][1]+hb01on_ori[i].first;
		int c2=activeparts[1][0];
		int n2=activeparts[1][1]+hb01on_ori[i].second;
		if(n1<0 || n1>=stchs[c1].size()) continue;
		if(n2<0 || n2>=stchs[c2].size()) continue;
		hb01on.push_back({{c1,n1},{c2,n2}});
	}
	/*std::ofstream ofs(hbfile);
	for(auto &h:hb01no)
		ofs <<h.first.first <<'\t' <<h.first.second <<"\tN\t" <<h.second.first <<'\t' <<h.second.second <<"\tO" <<std::endl;
	for(auto &h:hb01on)
		ofs <<h.first.first <<'\t' <<h.first.second <<"\tO\t" <<h.second.first <<'\t' <<h.second.second <<"\tN" <<std::endl;
	ofs.close();*/
	std::map<int,std::pair<int,int>> fxch;
	for(int i=0;i<stchs.size();i++) {
		if(i==activeparts[0][0] || i==activeparts[1][0]) {
			int j=i==activeparts[0][0]?0:1;
			if(activeparts[j][1]==0) {
				fxch.insert({i,{activeparts[j][2],stchs[i].size()-1}});
			} else {
				fxch.insert({i,{0,activeparts[j][1]-1}});
			}
		} else {
			fxch.insert({i,{0,stchs[i].size()}});
		}
	}
	std::vector<std::string> controllines;
	std::ifstream ifs(controlfile);
	std::string readline;
	while(std::getline(ifs,readline)) {
		if(readline.substr(0,10)=="RefPDBFile") readline="RefPDBFile\t=\t"+stfile;
		else if(readline.substr(0,13)=="AllFixedSites") {
			std::string l1{"AllFixedSites\t=\t"};
			for(auto &ch:fxch) {
				std::string l2="chain"+std::to_string(ch.first)+" "+std::to_string(ch.second.first)
						+"-"+std::to_string(ch.second.second)+" ";
				l1=l1+l2;
			}
			readline=l1;
		} else if(readline.substr(0,15)=="ExtendRestrains") {
			std::string l0=std::to_string(activeparts[0][0])+" "+std::to_string(activeparts[0][1])+" "
					+std::to_string(activeparts[0][1]==0?activeparts[0][2]:activeparts[0][2]-activeparts[0][1]-1);
			std::string l1=std::to_string(activeparts[1][0])+" "+std::to_string(activeparts[1][1])+" "
					+std::to_string(activeparts[1][1]==0?activeparts[1][2]:activeparts[1][2]-activeparts[1][1]-1);
			readline="ExtendRestrains\t=\t"+l0+" "+l1;
		}
		controllines.push_back(readline);
	}
	ifs.close();
	std::ofstream ofs(controlfile);
	for(auto &l:controllines) ofs <<l <<std::endl;
	ofs.close();
	printpdb(stfile,stchs);

	return acts;
}

bool InitialConfBuilder::hbpaired_SD(std::string controlname,
		std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01no_ori,
		std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01on_ori) {
	//build outfile successfully -> return true;
	int seed=NSPdstl::RandomEngine<>::getinstance().intrng(0,1000)();
	int step{sdstep_};
	double f_i_hb{hbforce_};
	double f_i_ca{0.0};
	//std::vector<std::vector<FullSite>> fschains=NSPproteinrep::readfullsitesfrompdb(pdbfile);
	std::vector<std::vector<FullSite>> fschains;
	std::vector<std::vector<Residue>> stchs;
	std::map<int,std::vector<int>> correspondindesign;
	std::vector<Residue> chtemp;
	for(int i=0;i<design_.size();i++) {
		if(design_[i].second.empty()) {
			if(chtemp.empty()) continue;
			stchs.push_back(chtemp);
			chtemp.clear();
			continue;
		}
		int stnum=chtemp.size();
		for(auto &r:design_[i].second) chtemp.push_back(r);
		int ennum=chtemp.size();
		correspondindesign.insert({i,{stchs.size(),stnum,ennum}});
	}
	if(!chtemp.empty()) stchs.push_back(chtemp);
	for(auto &f:fixed_) stchs.push_back(f);
	for(auto &st:stchs) {
		std::vector<FullSite> fss;
		for(auto &r:st) {
			FullSite fs;
			std::map<std::string, XYZ> cs;
			for(auto &c:r.rds()) {
				cs.insert({c.first,c.second.crd});
			}
			fs.changecrd(cs);
			fss.push_back(fs);
		}
		fschains.push_back(fss);
	}

	//std::string controlname="sdffcontrol";
	//NSPsd::genchainreadcontrols(controlfile,controlname);
	NSPsd::GenChain genchain(controlname+"_genchain");
	std::map<std::pair<int,int>,std::pair<int,int>> atomseq;
	std::vector<double> initcrd = genchain.getcrd(fschains,atomseq);
	for (auto &c : initcrd) c *= A2NM;
	std::vector<std::set<int>> notcis(fschains.size());
	NSPsd::ForceField ff=genchain.make_forcefield(controlname+"_ff",notcis);
	NSPsd::ActiveSelections *acts=genchain.setactiveselections(&ff);
	std::vector<bool> forceoff=acts->atomfixed();
	NSPsd::NeighborList nbl(initcrd,ff,&forceoff);
	NSPsd::SDRun::SDRunIn sdrunin(initcrd,controlname);
	NSPsd::SDRun sdrun = genchain.make_sdrun(sdrunin,seed,notcis);

	std::vector<std::pair<int,int>> hbps;
	for(auto &h:hb01no_ori) {
		int i1=atomseq.at(h.first).first;
		int i2=atomseq.at(h.second).first+atomseq.at(h.second).second-1;
		hbps.push_back({i1,i2});
	}
	for(auto &h:hb01on_ori) {
		int i1=atomseq.at(h.first).first+atomseq.at(h.first).second-1;
		int i2=atomseq.at(h.second).first;
		hbps.push_back({i1,i2});
	}
	std::vector<std::pair<int,int>> caps;
	for(int i=0;i<step/100;i++) {
		std::cout <<"SD: " <<i*100 <<std::endl;
		if(!(sdrun.runsteps_beta_sheet(100,hbps,f_i_hb,caps,f_i_ca))) {
			std::cout<<"Shake failure occurred."<<std::endl;
			return false;
		}
	}
	std::vector<std::vector<FullSite>> fsss=genchain.crd2fullsite(sdrun.state().crd, 1.0/A2NM);
	std::vector<std::vector<Residue>> chs1;
	for(auto &fss:fsss) {
		std::vector<Residue> rs;
		for(auto &fs:fss) rs.push_back(Residue(fs));
		chs1.push_back(rs);
	}

	for(int i=0;i<design_.size();i++) {
		if(design_[i].first==Fixed) continue;
		if(correspondindesign.find(i)==correspondindesign.end()) continue;
		std::vector<int> vi=correspondindesign.at(i);
		for(int j=vi[1];j<vi[2];j++) {
			design_[i].second[j-vi[1]]=chs1[vi[0]][j];
		}
	}
	return true;
}

bool InitialConfBuilder::hbpaired_SD(std::string controlfile, std::string pdbfile,
		std::string outfile, std::vector<int> &acts1,
		std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01no_ori,
		std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> &hb01on_ori) {
	//build outfile successfully -> return true;
	int seed=NSPdstl::RandomEngine<>::getinstance().intrng(0,1000)();
	int step{sdstep_};
	double f_i_hb{hbforce_};
	double f_i_ca{0.0};
	std::vector<std::vector<FullSite>> fschains=NSPproteinrep::readfullsitesfrompdb(pdbfile);
	std::string controlname="sdffcontrol";
	NSPsd::genchainreadcontrols(controlfile,controlname);
	NSPsd::GenChain genchain(controlname+"_genchain");
	std::map<std::pair<int,int>,std::pair<int,int>> atomseq;
	std::vector<double> initcrd = genchain.getcrd(fschains,atomseq);
	for (auto &c : initcrd) c *= A2NM;
	std::vector<std::set<int>> notcis(fschains.size());
	NSPsd::ForceField ff=genchain.make_forcefield(controlname+"_ff",notcis);
	NSPsd::ActiveSelections *acts=genchain.setactiveselections(&ff);
	std::vector<bool> forceoff=acts->atomfixed();
	NSPsd::NeighborList nbl(initcrd,ff,&forceoff);
	NSPsd::SDRun::SDRunIn sdrunin(initcrd,controlname);
	NSPsd::SDRun sdrun = genchain.make_sdrun(sdrunin,seed,notcis);

	std::vector<std::pair<int,int>> hbps;
	for(auto &h:hb01no_ori) {
		int i1=atomseq.at(h.first).first;
		int i2=atomseq.at(h.second).first+atomseq.at(h.second).second-1;
		hbps.push_back({i1,i2});
	}
	for(auto &h:hb01on_ori) {
		int i1=atomseq.at(h.first).first+atomseq.at(h.first).second-1;
		int i2=atomseq.at(h.second).first;
		hbps.push_back({i1,i2});
	}
	/*std::ifstream ifs(hbfile);
	std::string line;
	std::stringstream ss;
	while(std::getline(ifs,line)) {
		if(line[0]=='#') continue;
		int c1,c2,c3,c4;
		char noro1, noro2;
		ss << line;
		ss >>c1 >>c2 >>noro1 >>c3 >>c4 >>noro2;
		ss.clear();
		int i1=atomseq.at({c1,c2}).first;
		int i2=atomseq.at({c3,c4}).first;
		if(noro1=='O') i1=atomseq.at({c1,c2}).first+atomseq.at({c1,c2}).second-1;
		if(noro2=='O') i2=atomseq.at({c3,c4}).first+atomseq.at({c3,c4}).second-1;
		hbps.push_back({i1,i2});
	}
	ifs.close();*/

	std::vector<std::pair<int,int>> caps;
	/*ifs.open(cafile);
	while(std::getline(ifs,line)) {
		if(line[0]=='#') continue;
		int c1,c2,c3,c4;
		ss << line;
		ss >>c1 >>c2 >>c3 >>c4;
		ss.clear();
		int i1=atomseq.at({c1,c2}).first+1;
		int i2=atomseq.at({c3,c4}).first+1;
		caps.push_back({i1,i2});
	}
	ifs.close();*/

	for(int i=0;i<step/100;i++) {
		std::cout <<"SD: " <<i*100 <<std::endl;
		if(!(sdrun.runsteps_beta_sheet(100,hbps,f_i_hb,caps,f_i_ca))) {
			std::cout<<"Shake failure occurred."<<std::endl;
			return false;
		}
	}
	std::ofstream ofs(outfile);
	genchain.writepdb(sdrun.state().crd,ofs,1.0/A2NM);
	ofs.close();

	PdbReader_xy pr;
	pr.init(outfile);
	std::vector<std::vector<Residue>> chs1=pr.pep().chs();
	printpdb(outfile,chs1);

	int size0=design_[extendpair_[0]].second.size();
	int size1=design_[extendpair_[1]].second.size();
	for(int i=0;i<2;i++) {
		int n=extendpair_[i];
		int size=design_[n].second.size();
		for(int j=0;j<size;j++) {
			design_[n].second[j] = chs1[acts1[i*2]][acts1[i*2+1]+j];
			design_[n].second[j].changeaa(strandaa_);
		}
	}
	return true;
}

void InitialConfBuilder::extendstrand(std::vector<std::pair<int,int>> &hb01novi,
		std::vector<std::pair<int,int>> &hb01onvi, bool &firstendisinhb, bool &secondendisinhb) {
	bool n2c1;
	BackBoneSite st1;
	int ex1=extendpair_[0];
	if(ex1==0) {
		n2c1=false;
		st1=design_[ex1+1].second[0].getbackbonesite();
	} else {
		if(design_[ex1-1].first==Fixed) {
			n2c1=true;
			st1=design_[ex1-1].second.back().getbackbonesite();
		} else {
			n2c1=false;
			st1=design_[ex1+1].second[0].getbackbonesite();
		}
	}
	bool n2c2;
	BackBoneSite st2;
	int ex2=extendpair_[1];
	if(ex2==0) {
		n2c2=false;
		st2=design_[ex2+1].second[0].getbackbonesite();
	} else {
		if(design_[ex2-1].first==Fixed) {
			n2c2=true;
			st2=design_[1].second.back().getbackbonesite();
		} else {
			n2c2=false;
			st2=design_[ex2+1].second[0].getbackbonesite();
		}
	}
	assert(length_[ex1]==length_[ex2]);
	int lth=length_[ex1];
	bool parallel=n2c1==n2c2;
	std::set<int> firstinhb;
	if(parallel) {
		XYZ ca1=st1.cacrd(), o1=st1.ocrd(), ca2=st2.cacrd(), o2=st2.cacrd();
		if((o1-ca2).squarednorm()<(o2-ca1).squarednorm()) {
			for(int i=1;i<lth;i+=2) firstinhb.insert(i);
		} else {
			for(int i=0;i<lth;i+=2) firstinhb.insert(i);
		}
	} else {
		XYZ n1=st1.ncrd(), o2=st2.ocrd();
		if((n1-o2).squarednorm()<3.3*3.2) {
			for(int i=0;i<lth;i+=2) firstinhb.insert(i);
		} else {
			for(int i=1;i<lth;i+=2) firstinhb.insert(i);
		}
	}
	st1.resetpsi();
	st1.data_[BackBoneSite::OMIGA] = 180.0;
	if(c2nphi_.find(ex1)!=c2nphi_.end()) st1.data_[BackBoneSite::PHI]=c2nphi_.at(ex1);
	st2.resetpsi();
	st2.data_[BackBoneSite::OMIGA] = 180.0;
	if(c2nphi_.find(ex2)!=c2nphi_.end()) st2.data_[BackBoneSite::PHI]=c2nphi_.at(ex2);
	std::vector<BackBoneSite> bs1;
	std::vector<BackBoneSite> bs2;
	std::set<std::pair<int,int>> hb01no;
	std::set<std::pair<int,int>> hb01on;
	std::cout <<"building residue..." <<std::endl;
	for(int i=0;i<lth;i++) {
		//std::cout <<"building residue: " <<i <<std::endl;
		bs1.resize(i+1);
		bs2.resize(i+1);
		std::vector<double> vd;
		bool fihb=firstinhb.find(i)!=firstinhb.end();
		//std::cout <<"egodic:" <<std::endl;
		if(i==0) vd=hbpaired_ergodic(st1,st2,bs1[i],bs2[i],n2c1,fihb,parallel);
		else vd=hbpaired_ergodic(bs1[i-1],bs2[i-1],bs1[i],bs2[i],n2c1,fihb,parallel);
		//std::cout <<"MC" <<std::endl;
		if(i==0) hbpaired_MC(st1,st2,bs1[i],bs2[i],n2c1,fihb,parallel,vd);
		else hbpaired_MC(bs1[i-1],bs2[i-1],bs1[i],bs2[i],n2c1,fihb,parallel,vd);
		if(n2c1) {
			if(parallel) {
				if(fihb) {
					hb01no.insert({i+1,i});
					if(i==lth-1) {
						firstendisinhb=false;
						secondendisinhb=true;
					}
				} else {
					hb01on.insert({i,i+1});
				}
			} else {
				if(fihb) {
					hb01no.insert({i+1,i+1});
				} else {
					hb01on.insert({i,i});
				}
			}
		} else {
			if(parallel) {
				if(fihb) {
					hb01on.insert({i+1,i});
				} else {
					hb01no.insert({i,i+1});
				}
			} else {
				if(fihb) {
					hb01on.insert({i+1,i+1});
				} else {
					hb01no.insert({i,i});
				}
			}
		}
	}
	std::cout <<"build residue done!" <<std::endl;
	hb01novi.clear();
	for(auto h:hb01no) hb01novi.push_back(h);
	hb01onvi.clear();
	for(auto h:hb01on) hb01onvi.push_back(h);
	if(!n2c1) {
		std::vector<BackBoneSite> temp;
		for(int i=bs1.size()-1;i>=0;i--) temp.push_back(bs1[i]);
		bs1=temp;
		for(auto &h:hb01novi) h.first=lth-h.first-1;
		for(auto &h:hb01onvi) h.first=lth-h.first-1;
	}
	if(!n2c2) {
		std::vector<BackBoneSite> temp;
		for(int i=bs2.size()-1;i>=0;i--) temp.push_back(bs2[i]);
		bs2=temp;
		for(auto &h:hb01novi) h.second=lth-h.second-1;
		for(auto &h:hb01onvi) h.second=lth-h.second-1;
	}
	for(int i=0;i<bs1.size();i++) {
		design_[ex1].second.push_back(Residue(bs1[i],strandaa_));
	}
	for(int i=0;i<bs2.size();i++) {
		design_[ex2].second.push_back(Residue(bs2[i],strandaa_));
	}
	//std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> hb01no_new;
	//std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> hb01on_new;
	//std::vector<int> vis=hbpaired_SD(stfile, controlfile, hb01novi, hb01onvi, hb01no_new, hb01on_new);
	//return hbpaired_SD(controlfile, stfile, outfile,vis, hb01no_new, hb01on_new);
}

void InitialConfBuilder::sheet(std::string outfile) {
	XYZ ori, xaxis, xyplane;
	bool firstisfromN{true};
	bool secondisfromN{true};
	if(extendpair_[0]==0) {
		ori=design_[extendpair_[0]].second[0].rds().at("CA").crd;
		xaxis=design_[extendpair_[0]].second[2].rds().at("CA").crd;
	} else {
		if(design_[extendpair_[0]-1].first==Fixed) {
			ori=design_[extendpair_[0]].second.back().rds().at("CA").crd;
			xaxis=design_[extendpair_[0]].second[design_[extendpair_[0]].second.size()-3].rds().at("CA").crd;
			firstisfromN=false;
		} else {
			ori=design_[extendpair_[0]].second[0].rds().at("CA").crd;
			xaxis=design_[extendpair_[0]].second[2].rds().at("CA").crd;
		}
	}
	if(extendpair_[1]==0) {
		xyplane=design_[extendpair_[1]].second[0].rds().at("CA").crd;
	} else {
		if(design_[extendpair_[1]-1].first==Fixed) {
			xyplane=design_[extendpair_[1]].second.back().rds().at("CA").crd;
			secondisfromN=false;
		} else {
			xyplane=design_[extendpair_[1]].second[0].rds().at("CA").crd;
		}
	}
	std::vector<int> strandseq{extendpair_[0],extendpair_[1]};
	std::vector<bool> n2c{firstisfromN};
	int ixtemp1=extendpair_[0];
	int ixtemp2=extendpair_[1];
	bool n2ctemp=firstisfromN;
	while(strandseq.size()<betacontanctpair_.size()) {
		std::map<int,char> paired=betacontanctpair_.at(ixtemp1);
		int ix;
		bool finded{false};
		for(auto &p:paired) {
			if(ixtemp2!=p.first) {
				ix=p.first;
				finded=true;
				break;
			}
		}
		if(!finded) break;
		strandseq.resize(strandseq.size()+1);
		for(int i=strandseq.size()-1;i>0;i--) strandseq[i]=strandseq[i-1];
		strandseq[0]=ix;
		ixtemp2=ixtemp1;
		ixtemp1=ix;
		bool isparrallel=betacontanctpair_.at(ixtemp1).at(ixtemp2)=='p'?true:false;
		bool isfromN=isparrallel?(n2ctemp?true:false):(n2ctemp?false:true);
		n2ctemp=isfromN;
		n2c.resize(n2c.size()+1);
		for(int i=n2c.size()-1;i>0;i--) n2c[i]=n2c[i-1];
		n2c[0]=isfromN;
	}
	ixtemp2=extendpair_[0];
	ixtemp1=extendpair_[1];
	n2c.push_back(secondisfromN);
	n2ctemp=secondisfromN;
	while(strandseq.size()<betacontanctpair_.size()) {
		std::map<int,char> paired=betacontanctpair_.at(ixtemp1);
		int ix;
		for(auto &p:paired) {
			if(ixtemp2!=p.first) {
				ix=p.first;
				break;
			}
		}
		strandseq.push_back(ix);
		ixtemp2=ixtemp1;
		ixtemp1=ix;
		bool isparrallel=betacontanctpair_.at(ixtemp1).at(ixtemp2)=='p'?true:false;
		bool isfromN=isparrallel?(n2ctemp?true:false):(n2ctemp?false:true);
		n2c.push_back(isfromN);
		n2ctemp=isfromN;
	}
	int nori;
	for(int i=0;i<strandseq.size();i++) {
		if(strandseq[i]==extendpair_[0]) {
			nori=i;
			break;
		}
	}
	std::vector<XYZ> startpoints;
	double disstrand{5.0};
	for(int i=0;i<strandseq.size();i++) {
		double d=(double)(i-nori)*disstrand;
		startpoints.push_back(XYZ(0.0,d,0.0));
	}
	XYZ stranddir(1.0,0.0,0.0);
	LocalFrame lf=make_localframe(ori,xaxis,xyplane);
	stranddir=lf.local2globalcrd(stranddir)-ori;
	for(auto &sp:startpoints) {
		sp=lf.local2globalcrd(sp);
	}
	for(int i=0;i<strandseq.size();i++) {
		if(i==nori || i==nori+1) continue;
		int lth=length_[strandseq[i]];
		std::vector<BackBoneSite> ch=BackBoneBuilder::buildstrandat(lth, startpoints[i], stranddir, n2c[i]);
		for(auto &bs:ch) {
			bs.resname="GLY";
			design_[strandseq[i]].second.push_back(Residue(bs));
		}
	}
	//printpdb(outfile);

	double disstrandhelix{10.0};
	for(int i=2;i<design_.size()-2;i++) {
		if(design_[i].first!=Helix) continue;
		if(design_[i-2].first!=Strand || design_[i-1].first!=Loop || design_[i+1].first!=Loop || design_[i+2].first!=Strand) {
			std::cout <<"Helix Is Not In Two Strands!" <<std::endl;
			exit(1);
		}
		XYZ c1=design_[i-2].second.back().rds().at("CA").crd;
		XYZ c2=design_[i+2].second[0].rds().at("CA").crd;
		int seqinstrand, seqindesign;
		for(int j=0;j<strandseq.size();j++) {
			if(strandseq[j]==i-2 || strandseq[j]==i+2) {
				seqinstrand=j;
				seqindesign=strandseq[j];
				break;
			}
		}
		double dis1=disstrandhelix;
		if(seqindesign==i-2) {
			if(n2c[seqinstrand]) dis1=-disstrandhelix;
		} else {
			if(!n2c[seqinstrand]) dis1=-disstrandhelix;
		}
		c1=lf.global2localcrd(c1);
		c2=lf.global2localcrd(c2);
		c1.z_+=dis1;
		c2.z_+=dis1;
		c1=lf.local2globalcrd(c1);
		c2=lf.local2globalcrd(c2);
		c2=c2-c1;
		std::vector<BackBoneSite> bbss=BackBoneBuilder::buildhelixat(length_[i],c1,c2,true);
		for(auto &b:bbss) {
			b.resname=helixaa_;
			design_[i].second.push_back(Residue(b));
		}
	}
}

bool InitialConfBuilder::addloop() {
	int maxlooplth{10};
	int maxtrytime{10};
	std::vector<std::pair<int,int>> hrs;
	std::vector<std::pair<int,int>> srs;
	std::set<int> cis;
	for(int i=1;i<design_.size()-1;i++) {
		if(design_[i].first!=Loop) continue;
		BackBoneSite bsn=design_[i-1].second.back().getbackbonesite();
		BackBoneSite bsc=design_[i+1].second[0].getbackbonesite();
		std::vector<BackBoneSite> ch;
		for(int l=length_[i];l<maxlooplth;l++) {
			int ntime=-1;
			do {
				ntime++;
				auto bb=BackBoneBuilder::buildlinkers(l, bsn, bsc, hrs, srs, cis);
				if(bb.empty()) continue;
				int nchose=NSPdstl::RandomEngine<>::getinstance().intrng(0,bb.size()-1)();
				ch=*(bb[nchose]);
				break;
			} while(ntime<maxtrytime);
			if(!ch.empty()) break;
		}
		if(ch.empty()) {
			std::cout <<"can not find loop connection: " <<i <<std::endl;
			return false;
		}
		for(auto &h:ch) {
			h.resname = loopaa_;
			design_[i].second.push_back(Residue(h));
		}
	}
	return true;
}

void InitialConfBuilder::builddomain_only4nanobody(std::string outfile) {
	std::vector<std::vector<Residue>> chains(2);
	for(auto &r:design_[0].second) chains[0].push_back(r);
	for(auto &r:design_[1].second) {
		auto r1=r;
		r1.changeaa(strandaa_);
		chains[0].push_back(r1);
	}
	for(auto &r:design_[15].second) {
		auto r1=r;
		r1.changeaa(strandaa_);
		chains[1].push_back(r1);
	}
	for(auto &r:design_[16].second) chains[1].push_back(r);
	std::vector<Residue> rs0{chains[0][41],chains[0][42],chains[0][43],chains[0][44],chains[0][45],chains[0][46]};
	std::vector<Residue> rs1{chains[1][0],chains[1][1],chains[1][2],chains[1][3],chains[1][4],chains[1][5]};
	LocalFrame lf=make_localframe(rs0[5].rds().at("CA").crd,rs0[3].rds().at("CA").crd,rs1[0].rds().at("CA").crd);
	std::vector<std::pair<std::vector<Residue>,XYZ>> chs{{rs0,XYZ(0.0,-5.0,0.0)},{rs0,XYZ(0.0,-15.0,0.0)},
		{rs1,XYZ(0.0,-15.0,0.0)},{rs1,XYZ(0.0,5.0,0.0)}};
	XYZ ori=lf.local2globalcrd(XYZ(0.0,0.0,0.0));
	for(int i=0;i<chs.size();i++) {
		chs[i].second=lf.local2globalcrd(chs[i].second);
		for(int j=0;j<chs[i].first.size();j++) {
			auto & rds=chs[i].first[j].rds();
			for(auto &r:rds) r.second.crd=r.second.crd+chs[i].second-ori;
		}
	}
	XYZ st0(0.0,0.0,-10.0);
	st0=lf.local2globalcrd(st0);
	XYZ dir0(1.0,-1.0,0.0);
	dir0=lf.local2globalcrd(dir0)-ori;
	std::vector<BackBoneSite> bss0=BackBoneBuilder::buildhelixat(12,st0,dir0,true);
	XYZ st1(14.0,9.0,-10.0);
	st1=lf.local2globalcrd(st1);
	XYZ dir1(-3.0,-1.0,0.0);
	dir1=lf.local2globalcrd(dir1)-ori;
	std::vector<BackBoneSite> bss1=BackBoneBuilder::buildhelixat(12,st1,dir1,true);
	for(auto &ch:chs) chains.push_back(ch.first);
	std::vector<Residue> b0;
	for(auto &b:bss0) {
		b.resname=helixaa_;
		b0.push_back(Residue(b));
	}
	chains.push_back(b0);
	std::vector<Residue> b1;
	for(auto &b:bss1) {
		b.resname=helixaa_;
		b1.push_back(Residue(b));
	}
	chains.push_back(b1);
	//printpdb("temp",chains);

	std::vector<std::pair<BackBoneSite,BackBoneSite>> stenbs{
		{chains[0].back().getbackbonesite(),b0[0].getbackbonesite()},
		{b0.back().getbackbonesite(),chs[1].first[0].getbackbonesite()},
		{chs[1].first.back().getbackbonesite(),chs[2].first[0].getbackbonesite()},
		{chs[2].first.back().getbackbonesite(),chs[0].first[0].getbackbonesite()},
		{chs[0].first.back().getbackbonesite(),chs[3].first[0].getbackbonesite()},
		{chs[3].first.back().getbackbonesite(),b1[0].getbackbonesite()},
		{b1.back().getbackbonesite(),chains[1][0].getbackbonesite()}
	};
	std::vector<std::vector<Residue>> lps;
	std::vector<std::pair<int,int>> hrs;
	std::vector<std::pair<int,int>> srs;
	std::set<int> cis;
	//std::cout <<"loop start!" <<std::endl;
	for(int i=0;i<stenbs.size();i++) {
		//std::cout <<stenbs[i].first.cacrd().x_ <<'\t' <<stenbs[i].first.cacrd().y_ <<'\t' <<stenbs[i].first.cacrd().z_ <<std::endl;
		//std::cout <<stenbs[i].second.cacrd().x_ <<'\t' <<stenbs[i].second.cacrd().y_ <<'\t' <<stenbs[i].second.cacrd().z_ <<std::endl;
		std::vector<BackBoneSite> ch;
		for(int l=3;l<10;l++) {
			int ntime=-1;
			do {
				ntime++;
				auto bb=BackBoneBuilder::buildlinkers(l, stenbs[i].first, stenbs[i].second, hrs, srs, cis);
				if(bb.empty()) continue;
				ch=*(bb[0]);
				break;
			} while(ntime<10);
			if(!ch.empty()) break;
		}
		if(ch.empty()) {
			std::cout <<"can not find loop connection: " <<i <<std::endl;
			return;
		}
		std::vector<Residue> rs;
		for(auto &h:ch) {
			h.resname = loopaa_;
			rs.push_back(Residue(h));
		}
		lps.push_back(rs);
	}
	//std::cout <<"loop end!" <<std::endl;
	std::vector<Residue> finalchain;
	for(auto &r:chains[0]) finalchain.push_back(r);
	for(auto &r:lps[0]) finalchain.push_back(r);
	for(auto &r:b0) finalchain.push_back(r);
	for(auto &r:lps[1]) finalchain.push_back(r);
	for(auto &r:chs[1].first) finalchain.push_back(r);
	for(auto &r:lps[2]) finalchain.push_back(r);
	for(auto &r:chs[2].first) finalchain.push_back(r);
	for(auto &r:lps[3]) finalchain.push_back(r);
	for(auto &r:chs[0].first) finalchain.push_back(r);
	for(auto &r:lps[4]) finalchain.push_back(r);
	for(auto &r:chs[3].first) finalchain.push_back(r);
	for(auto &r:lps[5]) finalchain.push_back(r);
	for(auto &r:b1) finalchain.push_back(r);
	for(auto &r:lps[6]) finalchain.push_back(r);
	for(auto &r:chains[1]) finalchain.push_back(r);

	std::vector<std::vector<Residue>> newchain{finalchain};
	printpdb(outfile,newchain);
}

bool InitialConfBuilder::buildconf(std::string stfile, std::string controlfile, std::string outfile) {
	for(auto &d:design_) {
		if(d.first==Fixed) continue;
		d.second.clear();
	}
	/*BackBoneSite st;
	std::vector<Residue> es1=buildfirststrand(st);
	std::vector<Residue> es2=buildsecondstrand(es1,st);
	std::vector<std::vector<Residue>> rs{design_[0].second,es1,es2,design_.back().second};
	return rs;*/
	std::vector<std::pair<int,int>> hb01novi;
	std::vector<std::pair<int,int>> hb01onvi;
	std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> hb01no_new;
	std::vector<std::pair<std::pair<int,int>,std::pair<int,int>>> hb01on_new;
	std::string controlname="sdffcontrol";
	genchainreadcontrols(controlfile,controlname);
	if(!extendstrand(stfile, controlfile, outfile)) return false;

	//sheet(outfile);
	//if(!addloop()) return false;
	printpdb(outfile);
	return true;

	//build local frame
	//get strand
	//get helix
	//get loop
}

void InitialConfBuilder::printpdb(std::string outfile, const std::vector<std::vector<Residue>>& outpep) {
	std::ofstream ofs(outfile);
	int aid=0;
	for(int i=0;i<outpep.size();i++) {
		char cid='A'+i;
		for(int j=0;j<outpep[i].size();j++) {
			std::vector<PdbRecord> prs=outpep[i][j].getpdbrecords();
			for(PdbRecord &p:prs) {
				p.residueid = j;
				p.atomid = aid++;
				p.chainid = cid;
				ofs <<p.toString() <<std::endl;
			}
		}
	}
	ofs.close();
}

void InitialConfBuilder::printpdb(std::string outfile) {
	std::vector<std::vector<Residue>> stchs;
	std::vector<Residue> chtemp;
	for(int i=0;i<design_.size();i++) {
		if(design_[i].second.empty()) {
			if(chtemp.empty()) continue;
			stchs.push_back(chtemp);
			chtemp.clear();
			continue;
		}
		for(auto &r:design_[i].second) chtemp.push_back(r);
	}
	if(!chtemp.empty()) stchs.push_back(chtemp);
	for(auto &f:fixed_) stchs.push_back(f);
	printpdb(outfile,stchs);
}

InitialConfBuilder::InitialConfBuilder(NSPdataio::ParameterSet &pset) {
	std::string pdbfile;
	pset.getval("PDBFile",&pdbfile);
	PdbReader_xy pr;
	pr.init(pdbfile);
	std::vector<std::vector<Residue>> chs=pr.pep().chs();
	std::vector<std::vector<Residue>> chs_native=chs;

	std::vector<std::string> breakps;
	pset.getval("BreakPoint",&breakps);
	for(int i=0;i<breakps.size();i+=3) {
		char ch=breakps[i][0];
		std::string lab1=breakps[i+1];
		int n1;
		if(lab1!="_") n1=std::stoi(breakps[i+1]);
		std::string lab2=breakps[i+2];
		int n2;
		if(lab2!="_") n2=std::stoi(breakps[i+2]);

		int nch, r1{-1}, r2{10000};
		for(int a=0;a<chs.size();a++) {
			if(chs[a][0].chid()!=ch) continue;
			nch=a;
			for(int b=0;b<chs[a].size();b++) {
				if(lab1!="_" && chs[a][b].rid()==n1) r1=b;
				if(lab2!="_" && chs[a][b].rid()==n2) r2=b;
			}
			break;
		}
		std::vector<Residue> rss;
		for(int j=0;j<chs[nch].size();j++) {
			if(j>r1 && j<r2) continue;
			rss.push_back(chs[nch][j]);
		}
		chs[nch]=rss;
	}
	std::vector<std::pair<std::vector<Residue>,bool>> chs_new;
	for(auto &ch:chs) {
		std::vector<std::vector<Residue>> temp=Peptide::continuechain(ch);
		for(auto &t:temp) chs_new.push_back({t,true});
	}

	std::vector<std::string> ssseqs;
	pset.getval("SSSequence",&ssseqs);
	std::vector<std::pair<char,int>> phirds;
	for(int i=0;i<ssseqs.size();i++) {
		char c=ssseqs[i][0];
		int ii=std::stoi(ssseqs[i].substr(1));
		if(c=='l') {
			length_.push_back(ii);
			design_.push_back({InitialConfBuilder::Loop,std::vector<Residue>()});
		} else if(c=='h') {
			length_.push_back(ii);
			design_.push_back({InitialConfBuilder::Helix,std::vector<Residue>()});
		} else if(c=='e') {
			length_.push_back(ii);
			design_.push_back({InitialConfBuilder::Strand,std::vector<Residue>()});
		} else {
			for(int j=0;j<chs_new.size();j++) {
				if(chs_new[j].first[0].chid()==c && chs_new[j].first[0].rid()==ii ||
						chs_new[j].first.back().chid()==c && chs_new[j].first.back().rid()==ii) {
					chs_new[j].second=false;
					length_.push_back(chs_new[j].first.size());
					design_.push_back({InitialConfBuilder::Fixed,chs_new[j].first});
					break;
				}
			}
			if(i==0) continue;
			int lb=design_[design_.size()-2].first;
			if(lb==Loop || lb==Helix || lb==Strand) {
				for(int j=0;j<chs_native.size();j++) {
					for(int k=0;k<chs_native[j].size();k++) {
						if(chs_native[j][k].chid()==c && chs_native[j][k].rid()==ii) {
							double resphi=torsion(chs_native[j][k-1].rds().at("C").crd,chs_native[j][k].rds().at("N").crd,
									chs_native[j][k].rds().at("CA").crd,chs_native[j][k].rds().at("C").crd)*180.0/3.14159265;
							c2nphi_.insert({design_.size()-2,resphi});
							break;
						}
					}
				}
			}
		}
	}
	for(auto &ch:chs_new) {
		if(ch.second) fixed_.push_back(ch.first);
	}

	pset.getval("ExtendStrandPair",&extendpair_);
	std::vector<int> seqofsheet;
	std::string strandpairisparral;
	pset.getval("StrandHBPair",&seqofsheet);
	pset.getval("SheetDirection",&strandpairisparral);
	int ex0insheet, ex1insheet;
	for(int i=0;i<seqofsheet;i++) {
		if(seqofsheet[i]==extendpair_[0]) ex0insheet=i;
		if(seqofsheet[i]==extendpair_[1]) ex1insheet=i;
	}
	if(ex1insheet<ex0insheet) {

	} else {

	}


	std::vector<std::string> contanct;
	pset.getval("StrandHBPair",&contanct);
	for(int i=0;i<contanct.size();i+=3) {
		int i1=std::stoi(contanct[i]);
		int i2=std::stoi(contanct[i+1]);
		char c=contanct[i+2][0];
		if(design_[i1].first!=Strand || design_[i2].first!=Strand) {
			std::cout <<"Hydrogen-Bond Pair is not Strand :  " <<i/3 <<std::endl;
			exit(1);
		}
		if(betacontanctpair_.find(i1)==betacontanctpair_.end()) betacontanctpair_.insert({i1,std::map<int,char>()});
		if(betacontanctpair_.find(i2)==betacontanctpair_.end()) betacontanctpair_.insert({i2,std::map<int,char>()});
		betacontanctpair_.at(i1).insert({i2,c});
		betacontanctpair_.at(i2).insert({i1,c});
	}

	/*std::map<int,int> betainchain;
	int nbeta=0;
	for(int i=0;i<design_.size();i++) {
		if(design_[i].first!=InitialConfBuilder::Strand) continue;
		betainchain.insert({nbeta,i});
		nbeta++;
	}
	std::vector<std::string> contanct;
	pset.getval("StrandHBPair",&contanct);
	for(int i=0;i<contanct.size();i+=3) {
		int i1=betainchain.at(std::stoi(contanct[i]));
		int i2=betainchain.at(std::stoi(contanct[i+1]));
		char c=contanct[i+2][0];
		if(betacontanctpair_.find(i1)==betacontanctpair_.end()) betacontanctpair_.insert({i1,std::map<int,int>()});
		if(betacontanctpair_.find(i2)==betacontanctpair_.end()) betacontanctpair_.insert({i2,std::map<int,int>()});
		if(c=='t') {
			betacontanctpair_.at(i1).insert({i2,0});
			betacontanctpair_.at(i2).insert({i1,0});
		} else {
			betacontanctpair_.at(i1).insert({i2,1});
			betacontanctpair_.at(i2).insert({i1,1});
		}
	}*/

	/*std::vector<std::string> eps;
	pset.getval("ExtendStrandPair",&eps);
	for(int i=0;i<6;i+=3) {
		char c=eps[i][0];
		int n=std::stoi(eps[i+1]);
		int lth=std::stoi(eps[i+2]);
		for(int j=0;j<design_.size();j++) {
			if(design_[j].first!=Fixed) continue;
			const Residue & r1=design_[j].second[0];
			if(r1.chid()==c && r1.rid()==n) {
				extendpair_.push_back(j);
				break;
			}
			const Residue & r2=design_[j].second.back();
			if(r2.chid()==c && r2.rid()==n) {
				extendpair_.push_back(j);
				break;
			}
		}
	}*/


	std::vector<std::string> unbdir;
	pset.getval("UnBendDirection",&unbdir);
	unbenddir_.first = std::stoi(unbdir[0]);
	char chid=unbdir[1][0];
	int rid=std::stoi(unbdir[2]);
	std::string atomname=unbdir[3];
	for(int i=0;i<design_.size();i++) {
		if(design_[i].first!=Fixed) continue;
		for(int j=0;j<design_[i].second.size();j++) {
			const Residue & r=design_[i].second[j];
			if(r.chid()==chid && r.rid()==rid) {
				unbenddir_.second = r.rds().at(atomname).crd;
				break;
			}
		}
	}
	for(int i=0;i<fixed_.size();i++) {
		for(int j=0;j<fixed_[i].size();j++) {
			const Residue & r=fixed_[i][j];
			if(r.chid()==chid && r.rid()==rid) {
				unbenddir_.second = r.rds().at(atomname).crd;
				break;
			}
		}
	}

	std::vector<double> betaphi,betapsi;
	pset.getval("BetaPHI",&betaphi);
	pset.getval("BetaPSI",&betapsi);
	strandphirange_.first=betaphi[0];
	strandphirange_.second=betaphi[1];
	strandpsirange_.first=betapsi[0];
	strandpsirange_.second=betapsi[1];
	pset.getval("StrandAA",&strandaa_);
	pset.getval("HelixAA",&helixaa_);
	pset.getval("LoopAA",&loopaa_);
	pset.getval("Precision",&prec_);
	pset.getval("UnBendStep",&unbendstep_);
	pset.getval("StepLength",&steplength_);
	pset.getval("MCStepLength",&mcsteplth_);
	pset.getval("MCRange",&mcrange_);
	pset.getval("MCStep",&mcstep_);
	pset.getval("HydrogenBondLength",&hblth_);
	pset.getval("ErgodicLength",&erglth_);
	pset.getval("SDStep",&sdstep_);
	pset.getval("HBForce",&hbforce_);
}

void NSPallatom::readextendpar(std::string parfile) {
	NSPdataio::ControlFile cf;
	cf.readfile(parfile);
	std::vector<std::string> controllines=cf.getcontrolines("SSEXTEND");

	std::map<std::string,double> doublepars{{"Precision",0.0},{"StepLength",0.0},{"MCStepLength",0.0},
		{"MCRange",0.0},{"HydrogenBondLength",0.0},{"ErgodicLength",0.0}, {"HBForce",0.0}};
	std::map<std::string,std::vector<std::string>> stringvecpars{{"BreakPoint",{}},{"SSSequence",{}},
		{"StrandHBPair",{}},{"UnBendDirection",{}}};
	std::map<std::string,std::vector<double>> doublevecpars{{"BetaPHI",{}},{"BetaPSI",{}}};
	std::map<std::string,std::string> stringpars{{"PDBFile","pdb"},{"HelixAA","GLY"},{"StrandAA","GLY"},{"LoopAA","GLY"}};
	std::map<std::string,std::vector<int>> intvecpars{{"ExtendStrandPair",{}}};
	std::map<std::string,int>intpars{{"UnBendStep",1}, {"MCStep",0}, {"SDStep",0}};

	InitialConfBuilderControls::initdefaultkeyvals("SSEXTEND_1",doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	int nsuccess=InitialConfBuilderControls::adjustvalues("SSEXTEND_1",controllines);
	if(nsuccess!= controllines.size()) {
		exit(1);
	}
}

