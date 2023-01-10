/*
 * nnvterm.h
 *
 *  Created on: 2017骞�12鏈�12鏃�
 *      Author: hyliu
 */

#ifndef SD_NNTERM_H_
#define SD_NNTERM_H_
#include "sd/bsinchain.h"
#include "dstl/nnregressionmodel.h"
#include "dataio/datapaths.h"
#include "allatom/pdbreader_xy.h"
//#define USE_CB
namespace NSPsd{

class NNTerm{
public:
	NNTerm(){;}
	virtual NSPdstl::CombinedEstimator &getestimator()=0;
	virtual void setup(const std::vector<PhiPsiCodes> & phipsicodes, int posi){
		abort(); //this function should not be called;
	}
	double outvalue(std::vector<DvDxi> *dvdxi);
	//get called when derived NNterm yield a single value
	virtual void setup(const std::vector<double> &crd,
			const std::vector<BSInChain> &bsinchain1,
			const std::vector<PhiPsiCodes> & phipsicodes1, int posi1,
			const std::vector<BSInChain> &bsinchain2,
						const std::vector<PhiPsiCodes> & phipsicodes2, int posi2){
		setup(phipsicodes1,posi1);
	}
	virtual std::vector<double> outvalues (std::vector<std::vector<DvDxi>> *dvdxi){
			std::vector<double> res;
			dvdxi->push_back(std::vector<DvDxi>());
			res.push_back(outvalue(&(dvdxi->back())));
			return res;
	}
	virtual ~NNTerm(){;}
protected:
	std::vector<double>  codevec_;
	std::vector<std::vector<DvDxi>> dcdx_;
};
template <int DISTRTYPE=-1>
class PhiPsiNNTerm :public NNTerm{
public:
	enum DistrType {COIL,MIXCOIL,GLY,PRO,PREPRO,HELIX,STRAND,TRANSPRO};
	static NSPdstl::CombinedEstimator & gettheestimator(){
		static NSPdstl::CombinedEstimator lsestimator;
		static bool estimatorread { false };
		if (!estimatorread) {
#ifdef _OPENMP
#pragma omp critical(phpsinnterm_global)
			{
				if(!estimatorread){
#endif
			std::string filename;
			if(DISTRTYPE==MIXCOIL){
//				filename = NSPdataio::datapath() + "mixcoilphipsi_nnmodel.dat";
				filename = NSPdataio::datafilename("mixcoilphipsi_nnmodel.dat");
			} else if(DISTRTYPE==GLY){
				filename = NSPdataio::datafilename("glyphipsi_nnmodel.dat");
			} else if(DISTRTYPE==PRO){
				filename = NSPdataio::datafilename("prophipsi_nnmodel.dat");
			}
			std::ifstream ifs;
			ifs.open(filename.c_str());
			lsestimator.setup(ifs);
			estimatorread = true;
#ifdef _OPENMP
#pragma omp flush
				}
			}
#endif
		}
		return lsestimator;
	}
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator();
	}
	double phipsiene(const std::vector<PhiPsiCodes> & phipsicodes, int posi,
			std::vector<DvDxi> *dedx){
		codevec_.clear();
		dcdx_.clear();
		for (int i=0;i<6;++i){
					codevec_.push_back(phipsicodes[posi].phicodes[i]);
					dcdx_.push_back(phipsicodes[posi].dphicodesdx[i]);
			}
		for (int i=0;i<6;++i){
					codevec_.push_back(phipsicodes[posi].psicodes[i]);
					dcdx_.push_back(phipsicodes[posi].dpsicodesdx[i]);
		}
		return outvalue(dedx);
	}
private:
};
class LSNNTerm :public NNTerm{
public:
	static NSPdstl::CombinedEstimator & gettheestimator();
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator();
	}
	void setup(const std::vector<PhiPsiCodes> & phipsicodes, int posi);
private:
};
class CBData{
public:
	CBData(){;}
	CBData(std::vector<NSPgeometry::XYZ> &crdncac);
	const NSPgeometry::XYZ &cbcrd() const{ return cbcrd_;}
	const std::vector<std::vector<NSPgeometry::XYZ>> &drcb() const {return drcb_;}
	std::vector<NSPgeometry::XYZ> distributederiv(const NSPgeometry::XYZ derivrcb) const;
private:
	NSPgeometry::XYZ cbcrd_;
	std::vector<std::vector<NSPgeometry::XYZ>> drcb_;
};

class SitePairNNTerm :public NNTerm{
public:
	static NSPdstl::CombinedEstimator & gettheestimator();
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator();
	}
	virtual void setup(const std::vector<double> &crd,
			const std::vector<BSInChain> &bsinchain1,
			const std::vector<PhiPsiCodes> & phipsicodes1, int posi1,
			const std::vector<BSInChain> &bsinchain2,
			const std::vector<PhiPsiCodes> & phipsicodes2, int posi2);
	void setupwithcbdata(const std::vector<double> &crd,
				const std::vector<BSInChain> &bsinchain1, const CBData & cbdata1,
				const std::vector<PhiPsiCodes> & phipsicodes1, int posi1,
				const std::vector<BSInChain> &bsinchain2, const CBData &cbdata2,
				const std::vector<PhiPsiCodes> & phipsicodes2, int posi2);
	double rca() const {return rca_;}
	const std::vector<DvDxi> & drcadx() const {return drcadx_;}
private:
	void add_dmcodevec(const std::vector<double> &crd,
			const BSInChain &bs1,const BSInChain &bs2);
	void add_dmcodevec(const std::vector<double> &crd,
			const BSInChain &bs1,const CBData & cbdata1,
			const BSInChain &bs2,const CBData &cbdata2);
	std::vector<double> encode_r(double r,std::vector<double> *dcdr);
	double rca_;
	std::vector<DvDxi> drcadx_;
};
class NN_SSTerm {
public:
	enum SSIDX{HELIX,STRAND,COIL,TERMINUS};
	static int sstype(std::vector<double> probabilities){
		if(probabilities.empty()) return TERMINUS;
		int sidx=0;
		for(int i=1;i<3;++i){
			if(probabilities[i]>probabilities[sidx])sidx=i;
		}
		return sidx;
	}
	static NSPdstl::L3SoftMaxClassifier & getclassifier();
	std::vector<double> probabilities(const std::vector<PhiPsiCodes> & phipsicodes,int posi,
			std::vector<std::vector<DvDxi>> *dp3dx);
	double probability(int ssidx,const std::vector<PhiPsiCodes> &phipsicodes,int posi,
			std::vector<DvDxi> *dpdx){
		std::vector<std::vector<DvDxi>> dp3dx;
		std::vector<double> p3=probabilities(phipsicodes,posi,&dp3dx);
		*dpdx=dp3dx[ssidx];
		return p3[ssidx];
	}
private:
	NSPdstl::L3SoftMaxClassifier classifier_;
	void setup(const std::string & filename);
};
class ConformerCode;
class NN_KaiTerm:public NNTerm{
public:
	static NSPdstl::CombinedEstimator &gettheestimator(const std::string &restype);
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator(resname_);
	}
	NN_KaiTerm(){;}
	NN_KaiTerm(const ConformerCode &conformercode);
private:
	std::string resname_;
};
//side chain conformer terms for terminus residues
class NN_KaiTerm_T:public NNTerm{
public:
	static NSPdstl::CombinedEstimator &gettheestimator(const std::string &restype);
	NSPdstl::CombinedEstimator &getestimator(){
		return gettheestimator(resname_);
	}
	NN_KaiTerm_T(){;}
	NN_KaiTerm_T(const ConformerCode &conformercode);
private:std::string resname_;
};

class LocalBbHBNNTerm : public NNTerm
{
public:
    static NSPdstl::CombinedEstimator & gettheestimator();
    NSPdstl::CombinedEstimator &getestimator(){
        return gettheestimator();
    }
    virtual void setup(const std::vector<std::vector<double>> &dmatrix,
                       const std::vector<int> &atomids);
    virtual void setup(const std::vector<NSPgeometry::XYZ> &crds,
                        const std::vector<int> &atomids);
    virtual void setup(const std::vector<double> &crd,
            const std::vector<BSInChain> &bsinchain1, int posi1,
            const std::vector<BSInChain> &bsinchain2, int posi2);
    double rnomin() const {return rnomin_;}
    const std::vector<DvDxi> & drcdx() const {return drcdx_;}
private:
    void add_dmcodevec(const std::vector<std::vector<double>> &dmatrix,
                       const std::vector<int> &atomids);
    void add_dmcodevec(const std::vector<NSPgeometry::XYZ> &crd,
                       const std::vector<int> &atomids);
    void add_dmcodevec(const std::vector<double> &crd,
            const BSInChain &bs1,const BSInChain &bs2);
    std::vector<double> encode_r(double r,std::vector<double> *dcdr);
    double rnomin_;
    std::vector<DvDxi> drcdx_;
};

class LocalBbHBGeoNNTerm : public NNTerm
{
public:
    static NSPdstl::CombinedEstimator & gettheestimator();
    NSPdstl::CombinedEstimator &getestimator(){
        return gettheestimator();
    }
    double outvalue(std::vector<DvDxi> *dvdxi) {
        double result = NNTerm::outvalue(dvdxi);
        for (auto & d : *dvdxi) {
            d.first = atomids_[d.first];
            d.second = srno_ * d.second;
        }
        for (auto & d : dsrdx_) {
            dvdxi->push_back(std::make_pair(atomids_[d.first], result * d.second));
        }
        return result * srno_;
    }
    double outvalue() {
    	std::vector<DvDxi> dvdxi;
        double result = NNTerm::outvalue(&dvdxi);
        return result * srno_;
    }
	virtual void setup(const std::vector<NSPgeometry::XYZ> &crds,
                        const std::vector<int> &atomids) {
         std::vector<NSPgeometry::XYZ> internalcrds;
         for (int id : atomids) {
             internalcrds.push_back(crds[id]);
         }
         setup(internalcrds);
         atomids_ = atomids;
	}
    virtual void setup(const std::vector<double> &crd,
            const std::vector<BSInChain> &bsinchain1, int posi1,
            const std::vector<BSInChain> &bsinchain2, int posi2);
    double rnomin() const {return rnomin_;}
    const std::vector<DvDxi> & drcdx() const {return drcdx_;}
    virtual void setup(const std::vector<NSPgeometry::XYZ> &crds);
private:
    std::vector<double> encode_r(double r, std::vector<double> *dcdr);
    //virtual void setup(const std::vector<NSPgeometry::XYZ> &crds);
    int calc_geometry(const std::vector<NSPgeometry::XYZ> &coords,
                      std::vector<double> &geometry,
                      std::vector<std::vector<DvDxi>> &dgdxi);
    void make_codevec(const std::vector<double> &geometry,
                      std::vector<std::vector<DvDxi>> &dgdxi);
    double radian2degree(const double rad);
    std::vector<int> atomids_;
    double srno_;
    std::vector<DvDxi> dsrdx_;
    double rnomin_;
    std::vector<DvDxi> drcdx_;
    inline double switch1d(double x0, double x1, double x, double *deriv) {
        *deriv = 0.0;
        if (x <= x0)
            return 1.0;
        else if (x >= x1)
            return 0.0;
        if (x <= x0 + 0.5 * (x1 - x0)) {
            double xs = (x - x0) / (x1 - x0);
            double fx = 1.0 - 2.0 * xs * xs;
            *deriv = -4.0 * xs/(x1-x0);
            return fx;
        } else {
            double xs = (x1 - x) / (x1 - x0);
            double fx = 2.0 * xs * xs;
            *deriv = -4.0 * xs/(x1-x0);
            return fx;
        }
    }
};

class Atom3NNTerm {
public:
	static NSPdstl::CombinedEstimator & gettheestimator(
			const std::vector<std::string>&ns1,const std::vector<std::string>&ns2);
    NSPdstl::CombinedEstimator &getestimator(
    		const std::vector<std::string>&ns1,const std::vector<std::string>&ns2){
        return gettheestimator(ns1,ns2);
    }
    static std::vector<int> GetDisIndex(const std::vector<std::string>&ns);
    static std::vector<int> GetCrdIndex(const std::vector<std::string>&ns);

    void setup(const std::vector<std::string>&ns1,const std::vector<NSPgeometry::XYZ> &crds1,
    		const std::vector<std::string>&ns2,const std::vector<NSPgeometry::XYZ> &crds2);
    double outvalue(std::vector<DvDxi> *dvdxi);
    double outvalue(std::vector<DvDxi> *dvdxi, std::vector<int>&idx) {
        double result = outvalue(dvdxi);
        std::map<int,int> mps;
        if(chord_) {
        	mps.insert({0,idx[3]});
        	mps.insert({1,idx[4]});
        	mps.insert({2,idx[5]});
        	mps.insert({3,idx[0]});
        	mps.insert({4,idx[1]});
        	mps.insert({5,idx[2]});
        } else {
        	mps.insert({0,idx[0]});
        	mps.insert({1,idx[1]});
        	mps.insert({2,idx[2]});
        	mps.insert({3,idx[3]});
        	mps.insert({4,idx[4]});
        	mps.insert({5,idx[5]});
        }
        for (auto & d : *dvdxi) {
            d.first = mps.at(d.first);
            d.second = d.second*w5565_;
        }
        dw5565_ *= result;
        for(DvDxi &dd :dwdr_) dvdxi->push_back(dw5565_*dd);
        /*int n=0;
        for(int i=0;i<dcdx_.size();i++) {
        	for(int j=0;j<dcdx_[i].size();j++) {
        		DvDxi dd = dvdxi->at(n);
        		dvdxi->at(n).first  = mps.at(dd.first);
        		dvdxi->at(n).second = dd.second*w5565_ + dw5565_*result*dcdx_[i][j].second;
        		n++;
        	}
        }*/
        return result*w5565_;
    }
private:
	std::vector<double>  codevec_;
	std::vector<std::vector<DvDxi>> dcdx_;
    bool chord_{false};
    std::vector<std::string> names1;
    std::vector<std::string> names2;
    void changeorder(const std::vector<NSPgeometry::XYZ> &crds1,
    		const std::vector<NSPgeometry::XYZ> &crds2);
    void setup(const std::vector<NSPgeometry::XYZ> &crds1, const std::vector<NSPgeometry::XYZ> &crds2);
    void encode_dis(double r, std::vector<DvDxi>&dgdxi);
    void encode_ang(double t, std::vector<DvDxi>&dgdxi);
    void encode_dih(double t1, std::vector<DvDxi>&dgdxi1, double t, std::vector<DvDxi>&dgdxi);
    void encode_dih2(double t1, std::vector<DvDxi>&dgdxi1,
    		double t2, std::vector<DvDxi>&dgdxi2, double t, std::vector<DvDxi>&dgdxi);
    double weight5565(double d, double &dd);
    double w5565_{1.0};
    double dw5565_{0.0};
    std::vector<DvDxi> dwdr_;
};

}


#endif /* SD_NNTERM_H_ */
