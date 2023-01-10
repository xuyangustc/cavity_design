/*
 * enestatistics.h
 *
 *  Created on: 2018年4月6日
 *      Author: hyliu
 */

#ifndef ENESTATISTICS_H_
#define ENESTATISTICS_H_
#include "sd/forcefield.h"
namespace NSPproteinrep{
class EneHistogram{
public:
	EneHistogram(double emin, double emax, double step):
		emin_(emin),emax_(emax),step_(step){
		nbins_=(emax_-emin_+0.99*step_)/step_;
		emax_=emin_+step_*(double) nbins_;
		noccur_.assign(nbins_,0.0);}
	void accumulate(const std::vector<double> &esamples);
	std::vector<std::pair<double,double>> makehistogram();
private:
	double emin_;
	double emax_;
	double step_;
	int nbins_;
	std::vector<double> noccur_;
	std::vector<double> distr_;
	int valtoindex(double val){
		if(val<emin_) return 0;
		int index=(val-emin_)/step_;
		if(index>nbins_-1) index=nbins_-1;
		return index;
	}
	double indextoval(int i){
		return emin_+((double) i+0.5)*step_;
	}
};
class EneCompositions{
public:
	EneCompositions(NSPsd::ForceField *ff, const std::vector<double> &crd);
	const std::vector<double> &lsenergies() const {return lsenergies_;}
	const std::vector<double> &sitepairenergies() const {return sitepairenergies_;}
	const std::vector<double> &phipsienergies() const {return phipsienergies_;}
	double els_total() const{return els_;}
	double esp_total() const{return esp_;}
	double ephipsi_total() const {return ephipsi_;}
private:
		NSPsd::ForceField *ff_;
		std::vector<std::vector<NSPsd::PhiPsiCodes>> phipsicodes_;
		std::vector<double> lsenergies_;
		std::vector<double> sitepairenergies_;
		std::vector<double> phipsienergies_;
		double els_;
		double esp_;
		double ephipsi_;
};
}



#endif /* ENESTATISTICS_H_ */
