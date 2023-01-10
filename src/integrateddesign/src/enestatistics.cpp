/*
 * enestatistics.cpp
 *
 *  Created on: 2018年4月6日
 *      Author: hyliu
 */

#include "integrateddesign/enestatistics.h"
#include "sd/nnterm.h"

using namespace NSPproteinrep;
using namespace NSPsd;
using namespace NSPgeometry;

void EneHistogram::accumulate(const std::vector<double> &esamples){
	for(auto e:esamples){
		int index=valtoindex(e);
		noccur_[index] +=1.0;
	}
}
std::vector<std::pair<double,double>> EneHistogram::makehistogram(){
	std::vector<std::pair<double,double>> result;
	double ntot=0.0;
	for(auto n:noccur_) ntot+=n;
	for(int i=0;i<nbins_;++i){
		double e=indextoval(i);
		double n=noccur_[i]/ntot;
		result.push_back(std::make_pair(e,n));
	}
	return result;
}

EneCompositions::EneCompositions(ForceField *ff, const std::vector<double> &crd){
	ff_=ff;
	std::vector<std::vector<BSInChain>> & bsinchains=ff->bsinchains();
	int nchains = bsinchains.size();
	phipsicodes_.clear();
	for (int i = 0; i < nchains; ++i) {
		phipsicodes_.push_back(makephipsicodes(crd, bsinchains[i]));
	}
	ephipsi_=0.0;
	phipsienergies_.clear();
	double rad=180.0/3.14159265;
	const NSPpdbstatistics::PhiPsiDistr &distr
		=NSPpdbstatistics::PhiPsiDistr::mixcoildistr();
	for(int i=0;i<nchains;++i){
		for(int s=0;s<bsinchains[i].size();++s){
			auto & codes = phipsicodes_[i][s];
			double ene,dedphi,dedpsi;
			if (s == 0) {
				ene = distr.intplene_psi(codes.psi * rad, &dedpsi);
			} else if (s == bsinchains[i].size()-1) {
				ene = distr.intplene_phi(codes.phi * rad, &dedphi);
			} else {
				PhiPsiNNTerm<PhiPsiNNTerm<>::MIXCOIL> phipsinnterm;
				std::vector<DvDxi> dedxi;
				ene = phipsinnterm.phipsiene(phipsicodes_[i], s, &dedxi);
			}
			phipsienergies_.push_back(ene);
			ephipsi_+=ene;
		}
	}
	els_=0.0;
	lsenergies_.clear();
	for(int c=0;c<nchains;++c){
		for (int i = 2; i < phipsicodes_[c].size() - 2; ++i) {
			LSNNTerm lsnnterm;
			std::vector<DvDxi> dedx;
			lsnnterm.setup(phipsicodes_[c], i);
			double e = lsnnterm.outvalue(&dedx);
			els_+=e;
			lsenergies_.push_back(e);
		}
	}
	esp_=0.0;
	sitepairenergies_.clear();
	std::vector<std::vector<SSCode>> sscodes;
	for (int i = 0; i < nchains; ++i) {
		sscodes.push_back(estimatess(phipsicodes_[i]));
	}
	for (int n = 0; n < nchains; ++n) {
			const std::vector<BSInChain> &bsn = bsinchains [n];
			for (int m = n; m < nchains; ++m) {
				const std::vector<BSInChain> &bsm = bsinchains[m];
				for (int i = 1; i < bsn.size() - 1; ++i) {
					int jstart = 1;
					int jend = bsm.size() - 1;
					if (m == n) {
						jstart = i + 1;
					}
					for (int j = jstart; j < jend; ++j) {
						if (m == n) {
							if (ff->sitepairoff(sscodes[n], i, j))
								continue;
						}
						SitePairNNTerm spterm;
						spterm.setup(crd, bsn, phipsicodes_[n], i, bsm,
								phipsicodes_[m], j);
						double rca = spterm.rca();
						const std::vector<DvDxi> &drcadx = spterm.drcadx();
						std::vector<DvDxi> dedx;
						double e = spterm.outvalue(&dedx);
						if(e==0) continue;
						esp_+=e;
						sitepairenergies_.push_back(e);
					}
				}
			}
		}
}

