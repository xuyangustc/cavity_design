/*
 * nnterm.cpp
 *
 *  Created on: 2017骞�12鏈�12鏃�
 *      Author: hyliu
 */
#include "dataio/datapaths.h"
#include "geometry/calculators.h"
#include "sd/nnterm.h"
#include "sd/sidechainff.h"
#include "allatom/basic_info.h"

using namespace NSPsd;

NSPdstl::CombinedEstimator &LSNNTerm::gettheestimator() {
	static NSPdstl::CombinedEstimator lsestimator;
	static bool estimatorread { false };
	if (!estimatorread) {
#ifdef _OPENMP
#pragma omp critical(lsnnterm_global)
		{
			if(!estimatorread){
#endif
//		std::string filename = NSPdataio::datapath() + "ls_nnmodel.dat";
		std::string filename = NSPdataio::datafilename("ls_nnmodel.dat");
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

NSPdstl::CombinedEstimator &SitePairNNTerm::gettheestimator() {
	static NSPdstl::CombinedEstimator estimator;
	static bool estimatorread { false };
	if (!estimatorread) {
#ifdef _OPENMP
#pragma omp critical(spnnterm_global)
		{
			if(!estimatorread){
#endif
//		std::string filename = NSPdataio::datapath() + "sitepair_nnmodel.dat";
		std::string filename = NSPdataio::datafilename("sitepair_nnmodel.dat");
		std::ifstream ifs;
		ifs.open(filename.c_str());
		estimator.setup(ifs);
		estimatorread = true;
#ifdef _OPENMP
#pragma omp flush
			}
		}
#endif
	}
	return estimator;
}
NSPdstl::L3SoftMaxClassifier & NN_SSTerm::getclassifier() {
	static NSPdstl::L3SoftMaxClassifier classifier;
	static std::string filename {"ssclassifier.dat" };
	static bool classifierread{false};
	if (!classifierread){
#ifdef _OPENMP
#pragma omp critical(nnsterm_global)
		{
			if(!classifierread){
#endif
		std::ifstream ifs;
//		std::string fn=NSPdataio::datapath()+filename;
		std::string fn=NSPdataio::datafilename(filename);
		ifs.open(fn.c_str());
		classifier.setup(ifs);
		classifierread=true;
		ifs.close();
#ifdef _OPENMP
#pragma omp flush
			}
		}
#endif
	}
	return classifier;
}

double NNTerm::outvalue(std::vector<DvDxi> *dvdxi) {
	if (codevec_.empty())
		return 0.0;
	auto & estimator = getestimator();
	std::vector<double> dvdc;
	double res = estimator(codevec_, dvdc);
	dvdxi->clear();
	for (int i = 0; i < dvdc.size(); ++i) {
		for (int j = 0; j < dcdx_[i].size(); ++j) {
			dvdxi->push_back(dvdc[i] * dcdx_[i][j]);
		}
	}
	return res;
}

void LSNNTerm::setup(const std::vector<PhiPsiCodes> & phipsicodes, int posi) {
	codevec_.clear();
	dcdx_.clear();
	if (posi < 2 || posi >= phipsicodes.size() - 2)
		return;
	for (int i = posi - 2; i < posi + 3; ++i) {
		auto it = phipsicodes.begin() + i;
		if (i != posi - 2) {
			for (auto & c : it->phicodes)
				codevec_.push_back(c);
			for (auto &dc : it->dphicodesdx)
				dcdx_.push_back(dc);
		}
		if (i != posi + 2) {
			for (auto & c : it->psicodes)
				codevec_.push_back(c);
			for (auto &dc : it->dpsicodesdx)
				dcdx_.push_back(dc);
		}
	}
}

void SitePairNNTerm::setup(	const std::vector<double> &crd,
		const std::vector<BSInChain> &bsinchain1,
		const std::vector<PhiPsiCodes> & phipsicodes1, int posi1,
		const std::vector<BSInChain> &bsinchain2,
					const std::vector<PhiPsiCodes> & phipsicodes2, int posi2) {
	codevec_.clear();
	dcdx_.clear();
	if (posi1 < 1 || posi1 >= phipsicodes1.size() - 1)
		return;
	if (posi2 < 1 || posi2 >= phipsicodes2.size() - 1)
		return;
	std::vector<NSPgeometry::XYZ> deriv;
	rca_=NSPgeometry::distance(
			getxyz(crd,bsinchain1[posi1].caid),getxyz(crd,bsinchain2[posi2].caid),&deriv);
	drcadx_.clear();
	drcadx_.push_back(std::make_pair(bsinchain1[posi1].caid,deriv[0]));
	drcadx_.push_back(std::make_pair(bsinchain2[posi2].caid,deriv[1]));
	if(rca_>0.95)return;
	int nangcode = 6;
	for (int i = posi1 - 1; i < posi1 + 2; ++i) {
		auto it = phipsicodes1.begin() + i;
		if (i != posi1 - 1) {
			for (int i = 0; i < nangcode; ++i)
				codevec_.push_back(it->phicodes[i]);
			for (int i = 0; i < nangcode; ++i)
				dcdx_.push_back(it->dphicodesdx[i]);
		}
		if (i != posi1 + 1) {
			for (int i = 0; i < nangcode; ++i)
				codevec_.push_back(it->psicodes[i]);
			for (int i = 0; i < nangcode; ++i)
				dcdx_.push_back(it->dpsicodesdx[i]);
		}
	}
	for (int i = posi2 - 1; i < posi2 + 2; ++i) {
		auto it = phipsicodes2.begin() + i;
		if (i != posi2 - 1) {
			for (int i = 0; i < nangcode; ++i)
				codevec_.push_back(it->phicodes[i]);
			for (int i = 0; i < nangcode; ++i)
				dcdx_.push_back(it->dphicodesdx[i]);
		}
		if (i != posi2 + 1) {
			for (int i = 0; i < nangcode; ++i)
				codevec_.push_back(it->psicodes[i]);
			for (int i = 0; i < nangcode; ++i)
				dcdx_.push_back(it->dpsicodesdx[i]);
		}
	}
	add_dmcodevec(crd, bsinchain1[posi1], bsinchain2[posi2]);
}
void SitePairNNTerm::setupwithcbdata(	const std::vector<double> &crd,
		const std::vector<BSInChain> &bsinchain1, const CBData &cbdata1,
		const std::vector<PhiPsiCodes> & phipsicodes1, int posi1,
		const std::vector<BSInChain> &bsinchain2,  const CBData &cbdata2,
					const std::vector<PhiPsiCodes> & phipsicodes2, int posi2) {
	codevec_.clear();
	dcdx_.clear();
	if (posi1 < 1 || posi1 >= phipsicodes1.size() - 1)
		return;
	if (posi2 < 1 || posi2 >= phipsicodes2.size() - 1)
		return;
	std::vector<NSPgeometry::XYZ> deriv;
	rca_=NSPgeometry::distance(
			getxyz(crd,bsinchain1[posi1].caid),getxyz(crd,bsinchain2[posi2].caid),&deriv);
	drcadx_.clear();
	drcadx_.push_back(std::make_pair(bsinchain1[posi1].caid,deriv[0]));
	drcadx_.push_back(std::make_pair(bsinchain2[posi2].caid,deriv[1]));
	if(rca_>0.95)return;
	int nangcode = 6;
	for (int i = posi1 - 1; i < posi1 + 2; ++i) {
		auto it = phipsicodes1.begin() + i;
		if (i != posi1 - 1) {
			for (int i = 0; i < nangcode; ++i)
				codevec_.push_back(it->phicodes[i]);
			for (int i = 0; i < nangcode; ++i)
				dcdx_.push_back(it->dphicodesdx[i]);
		}
		if (i != posi1 + 1) {
			for (int i = 0; i < nangcode; ++i)
				codevec_.push_back(it->psicodes[i]);
			for (int i = 0; i < nangcode; ++i)
				dcdx_.push_back(it->dpsicodesdx[i]);
		}
	}
	for (int i = posi2 - 1; i < posi2 + 2; ++i) {
		auto it = phipsicodes2.begin() + i;
		if (i != posi2 - 1) {
			for (int i = 0; i < nangcode; ++i)
				codevec_.push_back(it->phicodes[i]);
			for (int i = 0; i < nangcode; ++i)
				dcdx_.push_back(it->dphicodesdx[i]);
		}
		if (i != posi2 + 1) {
			for (int i = 0; i < nangcode; ++i)
				codevec_.push_back(it->psicodes[i]);
			for (int i = 0; i < nangcode; ++i)
				dcdx_.push_back(it->dpsicodesdx[i]);
		}
	}
	add_dmcodevec(crd, bsinchain1[posi1],cbdata1,bsinchain2[posi2],cbdata2);
}


CBData::CBData(std::vector<NSPgeometry::XYZ> &crdncac) {
	static const double theta=109.5*3.14159265/180.0;
	static const double t=120.0*3.14159265/180.0;
	drcb_.assign(3,std::vector<NSPgeometry::XYZ>());
	for(int i=0;i<3;++i){
		for(int m=0;m<3;++m){
			crdncac[i][m] += 0.0001;
			NSPgeometry::XYZ cbp=NSPgeometry::InternaltoXYZ(crdncac[1],crdncac[2],
					crdncac[0],0.15,theta,t);
			crdncac[i][m]-= 0.0002;
			NSPgeometry::XYZ cbm=NSPgeometry::InternaltoXYZ(crdncac[1],crdncac[2],
						crdncac[0],0.15,theta,t);
			drcb_[i].push_back((cbp-cbm)/0.0002);
			crdncac[i][m]+=0.0001;
		}
	}
	cbcrd_=NSPgeometry::InternaltoXYZ(crdncac[1],crdncac[2],
			crdncac[0],0.15,theta,t);
}
std::vector<NSPgeometry::XYZ> CBData::distributederiv(const NSPgeometry::XYZ derivrcb)const {
	std::vector<NSPgeometry::XYZ> result;
	for(int i=0;i<3;++i){
		NSPgeometry::XYZ derivri(0,0,0);
		for(int m=0;m<3;++m){
			for(int k=0;k<3;++k){
				derivri[m] += derivrcb[k]*drcb_[i][m][k];
			}
		}
		result.push_back(derivri);
	}
	return result;
}
void SitePairNNTerm::add_dmcodevec(const std::vector<double> & crd,
		const BSInChain &bs1, const CBData & cbdata1,
		const BSInChain &bs2,const CBData &cbdata2) {
	std::vector<NSPgeometry::XYZ> crd1;
	std::vector<NSPgeometry::XYZ> crd2;
	std::vector<int> index1 = bs1.atomids();
	std::vector<int> index2 = bs2.atomids();
	for (auto id : index1)
		crd1.push_back(getxyz(crd, id));
	for (auto id : index2)
		crd2.push_back(getxyz(crd, id));
	const std::vector<std::vector<NSPgeometry::XYZ>> & drcb1=cbdata1.drcb();
	const std::vector<std::vector<NSPgeometry::XYZ>> & drcb2=cbdata2.drcb();
	crd1.push_back(cbdata1.cbcrd());
	index1.push_back(-1);
	crd2.push_back(cbdata2.cbcrd());
	index2.push_back(-2);
	std::vector<std::vector<double>> dm;
	std::vector<std::vector<std::vector<DvDxi>>>dmdx0;
	for (int i = 0; i < crd1.size(); ++i) {
		dm.push_back(std::vector<double>());
		dmdx0.push_back(std::vector<std::vector<DvDxi>>());
		auto &di = dm.back();
		auto & dmdxi = dmdx0.back();
		for (int j = 0; j < crd2.size(); ++j) {
			std::vector<NSPgeometry::XYZ> drdx;
			double r = distance(crd1[i], crd2[j], &drdx);
			dmdxi.push_back(std::vector<DvDxi>(2));
			dmdxi[j][0] = std::make_pair(index1[i], drdx[0]);
			dmdxi[j][1] = std::make_pair(index2[j], drdx[1]);
			di.push_back(r);
		}
	}
	std::vector<std::vector<DvDxi>> dcdx0;
	for (int i = 0; i < crd1.size(); ++i) {
		auto &di = dm[i];
		auto &dmdxi = dmdx0[i];
		for (int j = 0; j < crd2.size(); ++j) {
			if (i != 1 || j != 1) {
				codevec_.push_back(10.0 * (di[j] - dm[1][1])); //in angstrom
				dcdx0.push_back(std::vector<DvDxi>());
				auto & dcdx = dcdx0.back();
				dcdx.push_back(10.0 * dmdxi[j][0]);
				dcdx.push_back(10.0 * dmdxi[j][1]);
				dcdx.push_back(-10.0 * dmdx0[1][1][0]);
				dcdx.push_back(-10.0 * dmdx0[1][1][1]);
			}
			std::vector<double> dcdr;
			std::vector<double> rcode = encode_r(di[j], &dcdr);
			for (int l = 0; l < rcode.size(); ++l) {
				codevec_.push_back(rcode[l]);
//				codevec_.push_back(0.0);
				dcdx0.push_back(std::vector<DvDxi>());
				dcdx0.back().push_back(dcdr[l] * dmdxi[j][0]);
				dcdx0.back().push_back(dcdr[l] * dmdxi[j][1]);
//				dcdx_.back().push_back(0.0 * dmdxi[j][0]);
//				dcdx_.back().push_back(0.0 * dmdxi[j][1]);
			}
		}
	}
	for(auto &d:dcdx0){
		dcdx_.push_back(std::vector<DvDxi>());
		auto & db=dcdx_.back();
		for(auto &dx:d){
			if(dx.first>=0) {
				db.push_back(dx);
			}
			else if(dx.first==-1){
				for(int i=0;i<3;++i){
					NSPgeometry::XYZ dcdi(0,0,0);
					for(int m=0;m<3;++m)
						for(int k=0;k<3;++k) dcdi[m] += dx.second[k]*drcb1[i][m][k];
					db.push_back(std::make_pair(index1[i],dcdi));
				}
			} else if(dx.first==-2){
				for(int i=0;i<3;++i){
					NSPgeometry::XYZ dcdi(0,0,0);
					for(int m=0;m<3;++m)
						for(int k=0;k<3;++k) dcdi[m] += dx.second[k]*drcb2[i][m][k];
					db.push_back(std::make_pair(index2[i],dcdi));
				}
			}
		}
	}
}
void SitePairNNTerm::add_dmcodevec(const std::vector<double> & crd,
		const BSInChain &bs1, const BSInChain &bs2) {
	std::vector<NSPgeometry::XYZ> crd1;
	std::vector<NSPgeometry::XYZ> crd2;
	std::vector<int> index1 = bs1.atomids();
	std::vector<int> index2 = bs2.atomids();
	for (auto id : index1)
		crd1.push_back(getxyz(crd, id));
	for (auto id : index2)
		crd2.push_back(getxyz(crd, id));
	std::vector<std::vector<double>> dm;
	std::vector<std::vector<std::vector<DvDxi>>>dmdx0;
	for (int i = 0; i < crd1.size(); ++i) {
		dm.push_back(std::vector<double>());
		dmdx0.push_back(std::vector<std::vector<DvDxi>>());
		auto &di = dm.back();
		auto & dmdxi = dmdx0.back();
		for (int j = 0; j < crd2.size(); ++j) {
			std::vector<NSPgeometry::XYZ> drdx;
			double r = distance(crd1[i], crd2[j], &drdx);
			dmdxi.push_back(std::vector<DvDxi>(2));
			dmdxi[j][0] = std::make_pair(index1[i], drdx[0]);
			dmdxi[j][1] = std::make_pair(index2[j], drdx[1]);
			di.push_back(r);
		}
	}
	for (int i = 0; i < crd1.size(); ++i) {
		auto &di = dm[i];
		auto &dmdxi = dmdx0[i];
		for (int j = 0; j < crd2.size(); ++j) {
			if (i != 1 || j != 1) {
				codevec_.push_back(10.0 * (di[j] - dm[1][1])); //in angstrom
				dcdx_.push_back(std::vector<DvDxi>());
				auto & dcdx = dcdx_.back();
				dcdx.push_back(10.0 * dmdxi[j][0]);
				dcdx.push_back(10.0 * dmdxi[j][1]);
				dcdx.push_back(-10.0 * dmdx0[1][1][0]);
				dcdx.push_back(-10.0 * dmdx0[1][1][1]);
			}
			std::vector<double> dcdr;
			std::vector<double> rcode = encode_r(di[j], &dcdr);
			for (int l = 0; l < rcode.size(); ++l) {
				codevec_.push_back(rcode[l]);
//				codevec_.push_back(0.0);
				dcdx_.push_back(std::vector<DvDxi>());
				dcdx_.back().push_back(dcdr[l] * dmdxi[j][0]);
				dcdx_.back().push_back(dcdr[l] * dmdxi[j][1]);
//				dcdx_.back().push_back(0.0 * dmdxi[j][0]);
//				dcdx_.back().push_back(0.0 * dmdxi[j][1]);
			}
		}
	}
}

std::vector<double> SitePairNNTerm::encode_r(double r,
		std::vector<double> *dcdr) {
//parameters old
	//	static const std::vector<double> center { 0.25, 0.35, 0.45, 0.55,
//			0.65, 0.75,0.85,0.95};
//	static const std::vector<double> sigma2 { 0.005, 0.005, 0.005, 0.005,
//			0.005, 0.005, 0.005,0.01};
static const std::vector<double> center { 0.25, 0.35, 0.45, 0.55,
					0.65, 0.75,0.9,1.05};
static const std::vector<double> sigma2 { 0.005, 0.005, 0.005, 0.005,
					0.005, 0.005, 0.005,0.005};
	int ng = center.size();
	std::vector<double> w(ng, 0.);
	std::vector<double> dw(ng,0);
	dcdr->assign(ng, 0.0);
	double wtot = 0.0;
	double dwtot=0.0;
	for (int i = 0; i < ng; ++i) {
		double dr = r - center[i];
		double drs = dr / sigma2[i];
		double wt = exp(-dr * drs);
		dw[i] = -2.0 * drs * wt;
		w[i] = wt;
		wtot += wt;
		dwtot+=dw[i];
	}
	//fci=exp[-(r-ri)^2/sigmai_2]
	//wtot=sum(fci)
	//wi=fci/wtot
	for (int i = 0; i < ng; ++i) {
		w[i] /= wtot;
		(*dcdr)[i] =(dw[i]-w[i]*dwtot)/wtot;
	}
	return w;
}
std::vector<double> NN_SSTerm::probabilities(
		const std::vector<PhiPsiCodes> & phipsicodes, int posi,
		std::vector < std::vector<DvDxi>> *dp3dx) {
	if (posi < 3 || posi >= phipsicodes.size() - 3)
		return std::vector<double>();
	std::vector<double> codes;
	std::vector<std::vector<DvDxi>> dcdx;
	for (int p = posi - 3; p < posi + 4; ++p) {
		auto it = phipsicodes.begin() + p;
		if (p != posi - 3) {
			for (auto & c : it->phicodes)
				codes.push_back(c);
			for (auto &dc : it->dphicodesdx)
				dcdx.push_back(dc);
		}
		if (p != posi + 3) {
			for (auto & c : it->psicodes)
				codes.push_back(c);
			for (auto &dc : it->dpsicodesdx)
				dcdx.push_back(dc);
		}
	}
	std::vector<std::vector<double>> dpdc;
	std::vector<double> p3 = getclassifier()(codes, dpdc);
	dp3dx->assign(3, std::vector<DvDxi>());
	for (int i = 0; i < 3; ++i) {
		auto & dpdx = dp3dx->at(i);
		auto &dpidc = dpdc[i];
		for (int k = 0; k < dcdx.size(); ++k) {
			auto & dckdx = dcdx[k];
			double dp = dpidc[k];
			for (auto dx : dckdx) {
				dpdx.push_back(dp * dx);
			}
		}
	}
	return p3;
}
std::vector<std::string> kairesidues{"CYS","ASP","GLU","PHE","HIS","ILE",
	"LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"
};
NSPdstl::CombinedEstimator &NN_KaiTerm::gettheestimator(const std::string &restype){
	static std::map<std::string,NSPdstl::CombinedEstimator> estimators;
	static bool modelsread{false};
//	auto it=estimators.find(restype);
//	if(it !=estimators.end()) return it->second;
	if(modelsread) return estimators.at(restype);
#ifdef _OPENMP
#pragma omp critical(nnkaiterm_global)
	{
		if(!modelsread) {
#endif
	for(auto &res:kairesidues){
		std::string filename = NSPdataio::datafilename(res+"conf_nnmodel.dat");
		std::ifstream ifs;
		estimators.insert(std::make_pair(res,NSPdstl::CombinedEstimator()));
		ifs.open(filename.c_str());
		estimators.at(res).setup(ifs);
	}
	modelsread=true;
#ifdef _OPENMP
#pragma omp flush
		}
	}
#endif
	return estimators.at(restype);
}
NSPdstl::CombinedEstimator &NN_KaiTerm_T::gettheestimator(const std::string &restype){
	static std::map<std::string,NSPdstl::CombinedEstimator> estimators;
	static bool modelsread{false};
//	auto it=estimators.find(restype);
//	if(it !=estimators.end()) return it->second;
	if(modelsread) return estimators.at(restype);
#ifdef _OPENMP
#pragma omp critical(nnkai_t_term_global)
	{
		if(!modelsread) {
#endif
		for(auto & res:kairesidues){
			std::string filename = NSPdataio::datafilename(res+"_T_conf_nnmodel.dat");
			std::ifstream ifs;
			estimators.insert(std::make_pair(res,NSPdstl::CombinedEstimator()));
			ifs.open(filename.c_str());
			estimators.at(res).setup(ifs);
		}
		modelsread=true;
#ifdef _OPENMP
#pragma omp flush
		}
	}
#endif
	return estimators.at(restype);
}
NN_KaiTerm::NN_KaiTerm(const ConformerCode & conformercode){
	resname_=conformercode.restype;
	codevec_.clear();
	dcdx_.clear();
	PhiPsiCodes &phipsicodes=*(conformercode.phipsicodes);
	for(int i=0;i<6;++i){
		codevec_.push_back(phipsicodes.phicodes[i]);
		dcdx_.push_back(phipsicodes.dphicodesdx[i]);
	}
	for(int i=0;i<6;++i){
		codevec_.push_back(phipsicodes.psicodes[i]);
		dcdx_.push_back(phipsicodes.dpsicodesdx[i]);
	}
	for(int i=0;i<conformercode.torsioncodes.size();++i){
		const std::vector<double> &torsioncodes=conformercode.torsioncodes[i];
		for(int j=0;j<torsioncodes.size();++j){
			codevec_.push_back(torsioncodes[j]);
			dcdx_.push_back(conformercode.dtorsioncodesdx[i][j]);
		}
	}
}
NN_KaiTerm_T::NN_KaiTerm_T(const ConformerCode & conformercode){
	resname_=conformercode.restype;
	codevec_.clear();
	dcdx_.clear();
	for(int i=0;i<conformercode.torsioncodes.size();++i){
		const std::vector<double> &torsioncodes=conformercode.torsioncodes[i];
		for(int j=0;j<torsioncodes.size();++j){
			codevec_.push_back(torsioncodes[j]);
			dcdx_.push_back(conformercode.dtorsioncodesdx[i][j]);
		}
	}
}

NSPdstl::CombinedEstimator &LocalBbHBNNTerm::gettheestimator() {
    static NSPdstl::CombinedEstimator estimator;
    static bool estimatorread { false };
    if (!estimatorread) {
#ifdef _OPENMP
#pragma omp critical(lbbhbnnterm_global)
        {
            if(!estimatorread){
#endif
        std::string filename = NSPdataio::datafilename("localbbhb_nnmodel.dat");
        std::ifstream ifs;
        ifs.open(filename.c_str());
        estimator.setup(ifs);
        estimatorread = true;
#ifdef _OPENMP
#pragma omp flush
            }
        }
#endif
    }
    return estimator;
}

void LocalBbHBNNTerm::setup(const std::vector<std::vector<double>> &dmatrix,
                            const std::vector<int> &atomids) {
    codevec_.clear();
    dcdx_.clear();
    double rno = dmatrix[2][1];
    double ron = dmatrix[1][2];
    rnomin_ = std::min(rno, ron);
    if (rnomin_ > 0.45) return;
    add_dmcodevec(dmatrix, atomids);
}

void LocalBbHBNNTerm::setup(const std::vector<NSPgeometry::XYZ> &crds,
                            const std::vector<int> &atomids) {
    codevec_.clear();
    dcdx_.clear();
    double rno = NSPgeometry::distance(crds[2], crds[4]);
    double ron = NSPgeometry::distance(crds[1], crds[5]);
    rnomin_ = std::min(rno, ron);
    if (rnomin_ > 0.45) return;
    add_dmcodevec(crds, atomids);
}

void LocalBbHBNNTerm::setup(const std::vector<double> &crd,
        const std::vector<BSInChain> &bsinchain1, int posi1,
        const std::vector<BSInChain> &bsinchain2, int posi2) {
    codevec_.clear();
    dcdx_.clear();
    double rno = NSPgeometry::distance(
            getxyz(crd,bsinchain1[posi1].nid),getxyz(crd,bsinchain2[posi2].oid));
    double ron = NSPgeometry::distance(
            getxyz(crd,bsinchain1[posi1].oid),getxyz(crd,bsinchain2[posi2].nid));
    rnomin_ = std::min(rno, ron);
    if (rnomin_ > 0.45) return;
    add_dmcodevec(crd, bsinchain1[posi1], bsinchain2[posi2]);
}

std::vector<double> LocalBbHBNNTerm::encode_r(double r,
        std::vector<double> *dcdr) {
    static const std::vector<double> center { 0.25, 0.275, 0.3, 0.325,
                    0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65};
    static const std::vector<double> sigma2 { 0.00015, 0.00015, 0.00015, 0.00015,
                    0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};
    int ng = center.size();
    std::vector<double> w(ng, 0.);
    std::vector<double> dw(ng,0);
    dcdr->assign(ng, 0.0);
    double wtot = 0.0;
    double dwtot=0.0;
    for (int i = 0; i < ng; ++i) {
        double dr = r - center[i];
        double drs = dr / sigma2[i];
        double wt = exp(-dr * drs);
        dw[i] = -2.0 * drs * wt;
        w[i] = wt;
        wtot += wt;
        dwtot+=dw[i];
    }
    for (int i = 0; i < ng; ++i) {
        w[i] /= wtot;
        (*dcdr)[i] =(dw[i]-w[i]*dwtot)/wtot;
    }
    return w;
}

void LocalBbHBNNTerm::add_dmcodevec(const std::vector<std::vector<double>> &dmatrix,
                                    const std::vector<int> &atomids) {
    std::vector<int> id1;
    std::vector<int> id2;
    int natoms = dmatrix.size();
    for (int i = 0; i < natoms; i++) {
        id1.push_back(atomids[i]);
        id2.push_back(atomids[natoms+i]);
    }
    std::vector<std::vector<double>> dm;
    std::vector<std::vector<std::vector<DvDxi>>> dmdx0;
    for (int i = 0; i < dmatrix.size(); ++i) {
        dm.push_back(std::vector<double>());
        dmdx0.push_back(std::vector<std::vector<DvDxi>>());
        auto &di = dm.back();
        auto &dmdxi = dmdx0.back();
        for (int j = 0; j < dmatrix[0].size(); ++j) {
            std::vector<NSPgeometry::XYZ> drdx(2);
            dmdxi.push_back(std::vector<DvDxi>(2));
            dmdxi[j][0] = std::make_pair(id1[i], drdx[0]);
            dmdxi[j][1] = std::make_pair(id2[j], drdx[1]);
            di.push_back(dmatrix[i][j]);
        }
    }

    for (int i = 0; i < dm.size(); ++i) {
        auto &di = dm[i];
        auto &dmdxi = dmdx0[i];
        for (int j = 0; j < dm[0].size(); ++j) {
            if (i != 0 || j != 0) {
                codevec_.push_back(10.0 * (di[j] - dm[0][0])); //in angstrom
                dcdx_.push_back(std::vector<DvDxi>());
                auto & dcdx = dcdx_.back();
                dcdx.push_back(10.0 * dmdxi[j][0]);
                dcdx.push_back(10.0 * dmdxi[j][1]);
                dcdx.push_back(-10.0 * dmdx0[1][1][0]);
                dcdx.push_back(-10.0 * dmdx0[1][1][1]);
            }
            std::vector<double> dcdr;
            std::vector<double> rcode = encode_r(di[j], &dcdr);
            for (int l = 0; l < rcode.size(); ++l) {
                codevec_.push_back(rcode[l]);
                dcdx_.push_back(std::vector<DvDxi>());
                dcdx_.back().push_back(dcdr[l] * dmdxi[j][0]);
                dcdx_.back().push_back(dcdr[l] * dmdxi[j][1]);
            }
        }
    }
}

void LocalBbHBNNTerm::add_dmcodevec(const std::vector<NSPgeometry::XYZ> &crds,
                                    const std::vector<int> &atomids) {
    std::vector<NSPgeometry::XYZ> crd1;
    std::vector<int> id1;
    crd1.push_back(crds[0]);
    crd1.push_back(crds[1]);
    crd1.push_back(crds[2]);
    id1.push_back(atomids[0]);
    id1.push_back(atomids[1]);
    id1.push_back(atomids[2]);
    std::vector<NSPgeometry::XYZ> crd2;
    std::vector<int> id2;
    crd2.push_back(crds[3]);
    crd2.push_back(crds[4]);
    crd2.push_back(crds[5]);
    id2.push_back(atomids[3]);
    id2.push_back(atomids[4]);
    id2.push_back(atomids[5]);
    std::vector<std::vector<double>> dm;
    std::vector<std::vector<std::vector<DvDxi>>>dmdx0;
    for (int i = 0; i < crd1.size(); ++i) {
        dm.push_back(std::vector<double>());
        dmdx0.push_back(std::vector<std::vector<DvDxi>>());
        auto &di = dm.back();
        auto & dmdxi = dmdx0.back();
        for (int j = 0; j < crd2.size(); ++j) {
            std::vector<NSPgeometry::XYZ> drdx;
            double r = distance(crd1[i], crd2[j], &drdx);
            dmdxi.push_back(std::vector<DvDxi>(2));
            dmdxi[j][0] = std::make_pair(id1[i], drdx[0]);
            dmdxi[j][1] = std::make_pair(id2[j], drdx[1]);
            di.push_back(r);
        }
    }
    for (int i = 0; i < crd1.size(); ++i) {
        auto &di = dm[i];
        auto &dmdxi = dmdx0[i];
        for (int j = 0; j < crd2.size(); ++j) {
            if (i != 0 || j != 0) {
                codevec_.push_back(10.0 * (di[j] - dm[0][0])); //in angstrom
                dcdx_.push_back(std::vector<DvDxi>());
                auto & dcdx = dcdx_.back();
                dcdx.push_back(10.0 * dmdxi[j][0]);
                dcdx.push_back(10.0 * dmdxi[j][1]);
                dcdx.push_back(-10.0 * dmdx0[1][1][0]);
                dcdx.push_back(-10.0 * dmdx0[1][1][1]);
            }
            std::vector<double> dcdr;
            std::vector<double> rcode = encode_r(di[j], &dcdr);
            for (int l = 0; l < rcode.size(); ++l) {
                codevec_.push_back(rcode[l]);
                dcdx_.push_back(std::vector<DvDxi>());
                dcdx_.back().push_back(dcdr[l] * dmdxi[j][0]);
                dcdx_.back().push_back(dcdr[l] * dmdxi[j][1]);
            }
        }
    }
}

void LocalBbHBNNTerm::add_dmcodevec(const std::vector<double> &crd,
        const BSInChain &bs1, const BSInChain &bs2) {
    std::vector<NSPgeometry::XYZ> crds;
    std::vector<int> index1 = bs1.atomids();
    std::vector<int> index2 = bs2.atomids();
    crds.push_back(getxyz(crd, index1[2]));
    crds.push_back(getxyz(crd, index1[3]));
    crds.push_back(getxyz(crd, index1[0]));
    crds.push_back(getxyz(crd, index2[2]));
    crds.push_back(getxyz(crd, index2[3]));
    crds.push_back(getxyz(crd, index2[0]));
    std::vector<int> atomids;
    atomids.push_back(index1[0]);
    atomids.push_back(index1[1]);
    atomids.push_back(index1[2]);
    atomids.push_back(index2[0]);
    atomids.push_back(index2[1]);
    atomids.push_back(index2[2]);
    add_dmcodevec(crds, atomids);
}


NSPdstl::CombinedEstimator &LocalBbHBGeoNNTerm::gettheestimator() {
    static NSPdstl::CombinedEstimator estimator;
    static bool estimatorread { false };
    if (!estimatorread) {
#ifdef _OPENMP
#pragma omp critical(lbbhbgeonnterm_global)
        {
            if(!estimatorread){
#endif
        std::string filename = NSPdataio::datafilename("localbbhb_nnmodel.dat");
        std::ifstream ifs;
        ifs.open(filename.c_str());
        estimator.setup(ifs);
        estimatorread = true;
#ifdef _OPENMP
#pragma omp flush
            }
        }
#endif
    }
    return estimator;
}

void LocalBbHBGeoNNTerm::setup(const std::vector<NSPgeometry::XYZ> &crds) {
    codevec_.clear();
    dcdx_.clear();
    dsrdx_.clear();
    std::vector<double> geometry;
    std::vector<std::vector<DvDxi>> dgdxi;
    int result = calc_geometry(crds, geometry, dgdxi);
    if (result == 0) return;
    rnomin_ = geometry[0];
    double dsdr;
    srno_ = switch1d(0.4, 0.45, rnomin_, &dsdr);
    for (auto & dg : dgdxi[0]) {
        dsrdx_.push_back(dsdr*dg);
    }
    make_codevec(geometry, dgdxi);
}

void LocalBbHBGeoNNTerm::setup(const std::vector<double> &crd,
        const std::vector<BSInChain> &bsinchain1, int posi1,
        const std::vector<BSInChain> &bsinchain2, int posi2) {
    std::vector<NSPgeometry::XYZ> crds;
    std::vector<int> atomids;
    // posi 1
    int c1id = bsinchain1[posi1].cid;
    atomids.push_back(c1id);
    crds.push_back(getxyz(crd, c1id));
    int o1id = bsinchain1[posi1].oid;
    atomids.push_back(o1id);
    crds.push_back(getxyz(crd, o1id));
    int n1id = bsinchain1[posi1].nid;
    atomids.push_back(n1id);
    crds.push_back(getxyz(crd, n1id));
    // posi 2
    int c2id = bsinchain2[posi2].cid;
    atomids.push_back(c2id);
    crds.push_back(getxyz(crd, c2id));
    int o2id = bsinchain2[posi2].oid;
    atomids.push_back(o2id);
    crds.push_back(getxyz(crd, o2id));
    int n2id = bsinchain2[posi2].nid;
    atomids.push_back(n2id);
    crds.push_back(getxyz(crd, n2id));
    setup(crds);
}

std::vector<double> LocalBbHBGeoNNTerm::encode_r(double r,
        std::vector<double> *dcdr) {
    static const std::vector<double> center { 0.25, 0.275, 0.3, 0.325,
                    0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65};
    static const std::vector<double> sigma2 { 0.00015, 0.00015, 0.00015, 0.00015,
                    0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};
    int ng = center.size();
    std::vector<double> w(ng, 0.);
    std::vector<double> dw(ng,0);
    dcdr->assign(ng, 0.0);
    double wtot = 0.0;
    double dwtot=0.0;
    for (int i = 0; i < ng; ++i) {
        double dr = r - center[i];
        double drs = dr / sigma2[i];
        double wt = exp(-dr * drs);
        dw[i] = -2.0 * drs * wt;
        w[i] = wt;
        wtot += wt;
        dwtot+=dw[i];
    }
    for (int i = 0; i < ng; ++i) {
        w[i] /= wtot;
        (*dcdr)[i] =(dw[i]-w[i]*dwtot)/wtot;
    }
    return w;
}

double LocalBbHBGeoNNTerm::radian2degree(const double rad) {
    static const double PI = 3.14159265359;
    return rad / PI * 180.0;
}

int LocalBbHBGeoNNTerm::calc_geometry(const std::vector<NSPgeometry::XYZ> &coords,
                                      std::vector<double> &geometry,
                                      std::vector<std::vector<DvDxi>> &dgdxi) {
    geometry.clear();
    int ret = 0; // initialized as failure
    std::vector<int> abidx{0,1,2,3,4,5};
    std::vector<int> baidx{3,4,5,0,1,2};
    std::vector<NSPgeometry::XYZ> drdx1;
    std::vector<NSPgeometry::XYZ> drdx2;
    double n1o2dist = NSPgeometry::distance(coords[abidx[2]], coords[abidx[4]],&drdx1);
    double n2o1dist = NSPgeometry::distance(coords[baidx[2]], coords[baidx[4]],&drdx2);
    double dadistcutoff = 0.45; // in nm
    if (n1o2dist > dadistcutoff && n2o1dist > dadistcutoff) return ret;
    std::vector<NSPgeometry::XYZ> groupA;
    std::vector<NSPgeometry::XYZ> groupB;
    double dadist;
    double dangle; // angle around donor N
    double ddihed;
    double aangle; // angle around acceptor O
    double adihed;
    dgdxi.resize(5, std::vector<DvDxi>());
    std::vector<int> *idx;
    std::vector<NSPgeometry::XYZ> *drdx;
    if (n1o2dist < n2o1dist) {
        idx = &abidx;
        drdx = &drdx1;
        dadist = n1o2dist;
        ret = 1; // donor - acceptor
    } else {
        dadist = n2o1dist;
        idx = &baidx;
        drdx = &drdx2;
        ret = -1; // acceptor - donor
    }
    for (int i = 0; i < 3; i++) {
        groupA.push_back(coords[idx->at(i)]);
        groupB.push_back(coords[idx->at(i+3)]);
    }
    dgdxi[0].push_back(DvDxi(idx->at(2), drdx->at(0)));
    dgdxi[0].push_back(DvDxi(idx->at(4), drdx->at(1)));

    std::vector<NSPgeometry::XYZ> ddangdx;
    dangle = radian2degree(angle(groupA[0], groupA[2], groupB[1], &ddangdx)); // C1-N1-O2
    std::vector<NSPgeometry::XYZ> dddihdx;
    if (abs(dangle - 90.0) > 85.0) {
        ddihed = 0.0;
        dddihdx.assign(4,NSPgeometry::XYZ(0.0,0.0,0.0));
    } else {
        ddihed = radian2degree(torsion(groupA[1], groupA[0], groupA[2], groupB[1], &dddihdx)); // O1-C1=N1-O2
    }
    std::vector<NSPgeometry::XYZ> daangdx;
    aangle = radian2degree(angle(groupB[0], groupB[1], groupA[2], &daangdx)); // C2-O2-N1
    std::vector<NSPgeometry::XYZ> dadihdx;
    if (abs(aangle - 90.0) > 85.0) {
        adihed = 0.0;
        dadihdx.assign(4, NSPgeometry::XYZ(0.0, 0.0, 0.0));
    } else {
        adihed = radian2degree(torsion(groupB[2], groupB[0], groupB[1], groupA[2], &dadihdx)); // N2-C2=O2-N1
    }

    geometry.push_back(dadist);
    geometry.push_back(dangle);
    geometry.push_back(ddihed);
    geometry.push_back(aangle);
    geometry.push_back(adihed);

    dgdxi[1].push_back(DvDxi(idx->at(0),ddangdx[0]));
    dgdxi[1].push_back(DvDxi(idx->at(2),ddangdx[1]));
    dgdxi[1].push_back(DvDxi(idx->at(4),ddangdx[2]));
    /*std::vector<NSPgeometry::XYZ> temp;
    double ap=angle(groupA[0]+NSPgeometry::XYZ(0,0.001,0), groupA[2], groupB[1], &temp);
    double am=angle(groupA[0]+NSPgeometry::XYZ(0,-0.001,0), groupA[2], groupB[1], &temp);
    std::cout << ddangdx[0].y_  << ' ' << (ap - am) / 0.002<< std::endl;*/
    dgdxi[2].push_back(DvDxi(idx->at(1),dddihdx[0]));
    dgdxi[2].push_back(DvDxi(idx->at(0),dddihdx[1]));
    dgdxi[2].push_back(DvDxi(idx->at(2),dddihdx[2]));
    dgdxi[2].push_back(DvDxi(idx->at(4),dddihdx[3]));

    dgdxi[3].push_back(DvDxi(idx->at(3),daangdx[0]));
    dgdxi[3].push_back(DvDxi(idx->at(4),daangdx[1]));
    dgdxi[3].push_back(DvDxi(idx->at(2),daangdx[2]));

    dgdxi[4].push_back(DvDxi(idx->at(5),dadihdx[0]));
    dgdxi[4].push_back(DvDxi(idx->at(3),dadihdx[1]));
    dgdxi[4].push_back(DvDxi(idx->at(4),dadihdx[2]));
    dgdxi[4].push_back(DvDxi(idx->at(2),dadihdx[3]));

    return ret;
}

void LocalBbHBGeoNNTerm::make_codevec(const std::vector<double> &geometry,
                                      std::vector<std::vector<DvDxi>> &dgdxi) {
    codevec_.clear();
    dcdx_.clear();
    double dadist = geometry[0];
    static const double PI = 3.14159265953;
    std::vector<double> dcdr;
    std::vector<double> distcode = encode_r(dadist, &dcdr);
    int cidx = 0;
    for (double d : distcode) {
    	codevec_.push_back(d);
        dcdx_.push_back(std::vector<DvDxi>());
        for(auto &drdx : dgdxi[0]) {
            dcdx_.back().push_back(dcdr[cidx]*drdx);
        }
        cidx++;
    }
    int nang[3] = { 1, 2, 4 };
    for (int i = 1; i < 5; i++) {
        double t = geometry[i] * PI / 180.0;
        for (double na : nang) {
        	double c = cos(na*t);
            double s = sin(na*t);
            if (i % 2 == 1) {
                codevec_.push_back(c);
                dcdx_.push_back(std::vector<DvDxi>());
                double dcda = -na * s;
                for(auto &dadx:dgdxi[i]) dcdx_.back().push_back(dcda*dadx);
                codevec_.push_back(s);
                dcdx_.push_back(std::vector<DvDxi>());
                double dsda = na * c;
                for(auto & dadx:dgdxi[i]) dcdx_.back().push_back(dsda*dadx);
            } else {
                double t1 = geometry[i-1] * PI / 180.0;
                double ct1 = cos(t1);
                double st1 = sin(t1);
                codevec_.push_back(st1*c);
                dcdx_.push_back(std::vector<DvDxi>());
                double dcdt1 = ct1 * c;
                for(auto &dt1dx : dgdxi[i-1]) dcdx_.back().push_back(dcdt1*dt1dx);
                double dcdt = -na * s * st1;
                for(auto &dtdx : dgdxi[i]) dcdx_.back().push_back(dcdt*dtdx);
                codevec_.push_back(st1*s);
                dcdx_.push_back(std::vector<DvDxi>());
                dcdt1 = ct1 * s;
                for(auto &dt1dx:dgdxi[i-1]) dcdx_.back().push_back(dcdt1*dt1dx);
                dcdt = na * c * st1;
                for(auto &dtdx:dgdxi[i]) dcdx_.back().push_back(dcdt*dtdx);
            }
        }
    }
}



NSPdstl::CombinedEstimator &Atom3NNTerm::gettheestimator(
		const std::vector<std::string>&names1,const std::vector<std::string>&names2) {
	static std::map<std::vector<std::vector<std::string>>,NSPdstl::CombinedEstimator> a3estimator;
#pragma omp critical
	if (a3estimator.find({names1,names2})==a3estimator.end()) {
		NSPdstl::CombinedEstimator estimator;
		std::string filename="atom3/Atom3NNModel/";
		if(names2<names1) {
			filename = filename +names2[0] +"_" +names2[1] +"_" +names2[2] +"_" +names2[3]
						+"_" +names1[0] +"_" +names1[1] +"_" +names1[2] +"_" +names1[3];
		} else {
			filename = filename +names1[0] +"_" +names1[1] +"_" +names1[2] +"_" +names1[3]
						+"_" +names2[0] +"_" +names2[1] +"_" +names2[2] +"_" +names2[3];
		}
		filename = NSPdataio::datafilename(filename);
		std::ifstream ifs;
		ifs.open(filename.c_str());
		//std::cout <<filename <<std::endl;
		estimator.setup(ifs);
		ifs.close();
		//std::cout <<"done" <<std::endl;
		a3estimator.insert({{names1,names2},estimator});
		a3estimator.insert({{names2,names1},estimator});
	}
	return a3estimator.at({names1,names2});
}

/*
NSPdstl::CombinedEstimator &Atom3NNTerm::gettheestimator(const std::vector<std::string>&names1,const std::vector<std::string>&names2) {
	static std::map<std::vector<std::vector<std::string>>,NSPdstl::CombinedEstimator> a3estimator;
	if(names1.empty()) {
		std::vector<std::vector<std::string>> gs;
		std::map<std::string,std::vector<std::vector<std::string>>> a3=NSPallatom::AA_Atom::Atom3();
		for(auto &a:a3) {
			for(auto &v:a.second) {
				std::vector<std::string> vs{a.first};
				for(auto &s:v) {
					vs.push_back(s);
				}
				gs.push_back(vs);
			}
		}
		for(int i=0;i<gs.size();i++) {
			std::vector<std::string> n1=gs[i];
			for(int j=i;j<gs.size();j++) {
				std::vector<std::string> n2=gs[j];
				NSPdstl::CombinedEstimator estimator;
				std::string filename="atom3/Atom3NNModel/";
				if(n2<n1) {
					filename = filename +n2[0] +"_" +n2[1] +"_" +n2[2] +"_" +n2[3]
								+"_" +n1[0] +"_" +n1[1] +"_" +n1[2] +"_" +n1[3];
				} else {
					filename = filename +n1[0] +"_" +n1[1] +"_" +n1[2] +"_" +n1[3]
								+"_" +n2[0] +"_" +n2[1] +"_" +n2[2] +"_" +n2[3];
				}
				filename = NSPdataio::datafilename(filename);
				std::ifstream ifs;
				ifs.open(filename.c_str());
				//std::cout <<filename <<std::endl;
				estimator.setup(ifs);
				ifs.close();
				//std::cout <<"done" <<std::endl;
				a3estimator.insert({{n1,n2},estimator});
				a3estimator.insert({{n2,n1},estimator});
			}
		}
		std::vector<std::string> temp{"MCTRANS","N","C","O"};
		a3estimator.insert({{names1,names2},a3estimator.at({temp,temp})});
	}
	return a3estimator.at({names1,names2});
}*/

void Atom3NNTerm::changeorder(const std::vector<NSPgeometry::XYZ> &crds1,
		const std::vector<NSPgeometry::XYZ> &crds2) {
	if(names2<names1) {
		chord_=true;
		return;
	}
	/*else if(names1==names2) {
		double d1=(crds1[0]-crds2[2]).squarednorm();
		double d2=(crds1[2]-crds2[0]).squarednorm();
		if(d2<d1) {
			chord_=true;
			return;
		}
	}*/
}

void Atom3NNTerm::encode_dis(double r, std::vector<DvDxi>&dgdxi) { //r is NM
    //static const std::vector<double> center { 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    //	0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8,
	//	0.9, 1.0, 1.2, 1.4};
    //static const std::vector<double> sigma2 { 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015,
    //	0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015,
    //    0.0005, 0.0005, 0.00075, 0.00075};

//	static const std::vector<double> center { 0.2, 0.25, 0.3,  0.35,
//		0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
//		0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
//		1.0, 1.05, 1.1, 1.15, 1.2, 1.25,
//		1.3, 1.35, 1.4};
//	static const std::vector<double> sigma2 { 0.0005, 0.0005, 0.0005, 0.0005,
//		0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
//		0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
//		0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
//		0.0005, 0.0005, 0.0005};

//	static const std::vector<double> center { 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375,
//		0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675,
//		0.7, 0.725, 0.75, 0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975,
//		1.0, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275,
//		1.3, 1.325, 1.35, 1.375, 1.4};
//	static const std::vector<double> sigma2 { 0.00015, 0.00015, 0.00015, 0.00015,
//	 	0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015,
//		0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015,
//		0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015,
//		0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015,
//		0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015, 0.00015};


	//static const std::vector<double> center { 0.2,0.225,0.25,0.275, 0.3,0.325,0.35, 0.4,0.45, 0.5};
    //static const std::vector<double> sigma2 { 0.0005,0.0005,0.0005,0.0005, 0.0005,0.0005,0.0005,
    //	0.0005, 0.0005, 0.0005};

	//static const std::vector<double> center { 0.1,0.2, 0.25,0.3, 0.35,0.4, 0.45,0.5, 0.55,0.6,
	//		0.65,0.7, 0.75,0.8, 0.85,0.9, 0.95 };
    //static const std::vector<double> sigma2 { 0.0005,0.0005,0.0005,0.0005, 0.0005,0.0005,0.0005,0.0005,
    //	0.0005,0.0005,0.0005,0.0005, 0.0005,0.0005,0.0005,0.0005, 0.0005};

	//static const std::vector<double> center { 0.1,0.2, 0.25,0.275, 0.3,0.325, 0.35,0.4,
	//	0.45,0.5, 0.55,0.65, 0.75,0.85, 0.95 };
   // static const std::vector<double> sigma2 { 0.0005,0.0005,0.0005,0.0005, 0.0005,0.0005,0.0005,0.0005,
    //	0.0005,0.0005,0.0005,0.0005, 0.0005,0.0005,0.0005};

	//static const std::vector<double> center { 0.1,0.15,0.2,0.25,0.275, 0.3,0.325,0.35,0.4,0.45, 0.5,0.55,0.6,0.65,
	//	0.7,0.75,0.8,0.85, 0.9,0.95,1.0,1.1, 1.2,1.3,1.4,1.5 };
    //static const std::vector<double> sigma2 { 0.0005,0.0005,0.0005,0.0005, 0.0005,0.0005,0.0005,0.0005,
    // 	0.0005,0.0005,0.0005,0.0005, 0.0005,0.0005,0.0005,0.0005, 0.0005,0.0005,0.0005,0.0005,
	//	0.0005,0.0005,0.0005,0.0005 ,0.0005,0.0005 };

	static const std::vector<double> center { 0.25,0.35, 0.45,0.55, 0.7,0.95, 1.2 };
    static const std::vector<double> sigma2 { 0.005,0.005,0.005,0.005, 0.01125, 0.015,0.03 };
	//static const std::vector<double> center { 0.2,0.25,0.3,0.35, 0.4,0.45,0.5,0.55, 0.6,0.65,0.7,0.75,
	//	0.8,0.85,0.9,0.95,1.0,1.2 };
    //static const std::vector<double> sigma2 { 0.00125,0.00125,0.00125, 0.00125,0.00125,0.00125,
    //	0.00125,0.00125,0.00125,0.00125, 0.00125,0.00125,0.00125,0.00125, 0.00125,0.00125,0.00125,0.00125 };
    //static const std::vector<double> center { 0.25, 0.275, 0.3, 0.325,
    //                0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65};
    //static const std::vector<double> sigma2 { 0.00015, 0.00015, 0.00015, 0.00015,
    //                0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};

	//static const std::vector<double> center { 0.2, 0.25, 0.3, 0.35,
	//	0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
	//	0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05};
	//static const std::vector<double> sigma2 { 0.0005, 0.0005, 0.0005, 0.0005,
	//	0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005,
	//	0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};

	//static const std::vector<double> center { 0.0,0.05, 0.1,0.15, 0.2,0.25,
	//	0.3,0.35, 0.4,0.45, 0.5,0.55, 0.6,0.65, 0.7,0.75, 0.8,0.85, 0.9,0.95, 1.0 };
	//static const std::vector<double> sigma2 { 0.0005,0.0005, 0.0005,0.0005, 0.0005,0.0005,
	//	0.0005,0.0005, 0.0005,0.0005, 0.0005,0.0005, 0.0005,0.0005, 0.0005,0.0005,
	//	00.0005,0.0005, 0.0005,0.0005, 0.0005 };

	//unit is A
//	static const std::vector<double> center { 2.0, 2.5, 3.0, 3.5,
//		4.0, 4.5, 5.0, 5.5, 6.0, 6.5,
//		7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5};
//	static const std::vector<double> sigma2 { 0.05, 0.05, 0.05, 0.05,
//		0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
//		0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

    int ng = center.size();
    std::vector<double> w(ng, 0.0);
    std::vector<double> dw(ng,0.0);
    std::vector<double> dcdr(ng,0.0);
    double wtot = 0.0;
    double dwtot=0.0;
    for (int i = 0; i < ng; ++i) {
        double dr = r - center[i];
        double drs = dr / sigma2[i];
        double wt = exp(-dr * drs);
        dw[i] = -2.0 * drs * wt;
        w[i] = wt;
        wtot += wt;
        dwtot+=dw[i];
    }
    for (int i = 0; i < ng; ++i) {
        w[i] /= wtot;
        dcdr[i] =(dw[i]-w[i]*dwtot)/wtot;
    }
    int cidx = 0;
    for (double d : w) {
    	codevec_.push_back(d);
        dcdx_.push_back(std::vector<DvDxi>());
        for(auto &drdx : dgdxi) {
            dcdx_.back().push_back(dcdr[cidx]*drdx);
        }
        cidx++;
    }
}

void Atom3NNTerm::encode_ang(double t, std::vector<DvDxi>&dgdxi) {
	std::vector<double> nang{ 1.0, 2.0, 4.0 };
    for (double na : nang) {
    	double c = cos(na*t);
        double s = sin(na*t);
        codevec_.push_back(c);
        dcdx_.push_back(std::vector<DvDxi>());
        double dcda = -na * s;
        for(auto &dadx:dgdxi) dcdx_.back().push_back(dcda*dadx);
        codevec_.push_back(s);
        dcdx_.push_back(std::vector<DvDxi>());
        double dsda = na * c;
        for(auto & dadx:dgdxi) dcdx_.back().push_back(dsda*dadx);
    }
}

void Atom3NNTerm::encode_dih(double t1, std::vector<DvDxi>&dgdxi1, double t, std::vector<DvDxi>&dgdxi) {
	std::vector<double> nang{ 1.0, 2.0, 4.0 };
    double ct1 = cos(t1);
    double st1 = sin(t1);
    for (double na : nang) {
    	double c = cos(na*t);
        double s = sin(na*t);
        codevec_.push_back(st1*c);
        dcdx_.push_back(std::vector<DvDxi>());
        double dcdt1 = ct1 * c;
        for(auto &dt1dx : dgdxi1) dcdx_.back().push_back(dcdt1*dt1dx);
        double dcdt = -na * s * st1;
        for(auto &dtdx : dgdxi) dcdx_.back().push_back(dcdt*dtdx);
        codevec_.push_back(st1*s);
        dcdx_.push_back(std::vector<DvDxi>());
        dcdt1 = ct1 * s;
        for(auto &dt1dx:dgdxi1) dcdx_.back().push_back(dcdt1*dt1dx);
        dcdt = na * c * st1;
        for(auto &dtdx:dgdxi) dcdx_.back().push_back(dcdt*dtdx);
    }
}

void Atom3NNTerm::encode_dih2(double t1, std::vector<DvDxi>&dgdxi1,
		double t2, std::vector<DvDxi>&dgdxi2, double t, std::vector<DvDxi>&dgdxi) {
	std::vector<double> nang{ 1.0, 2.0, 4.0 };
    double ct1 = cos(t1);
    double st1 = sin(t1);
    double ct2 = cos(t2);
    double st2 = sin(t2);
    for (double na : nang) {
		double c = cos(na*t);
		double s = sin(na*t);
		codevec_.push_back(st1*st2*c);
		dcdx_.push_back(std::vector<DvDxi>());
		double dcdt1 = ct1 * c * st2;
		for(auto &dt1dx : dgdxi1) dcdx_.back().push_back(dcdt1*dt1dx);
		double dcdt2 = ct2 * c * st1;
		for(auto &dt2dx : dgdxi2) dcdx_.back().push_back(dcdt2*dt2dx);
		double dcdt = -na * s * st1 * st2;
		for(auto &dtdx : dgdxi) dcdx_.back().push_back(dcdt*dtdx);
		codevec_.push_back(st1*st2*s);
		dcdx_.push_back(std::vector<DvDxi>());
		dcdt1 = ct1 * s * st2;
		for(auto &dt1dx:dgdxi1) dcdx_.back().push_back(dcdt1*dt1dx);
		dcdt2 = ct2 * s * st1;
		for(auto &dt2dx:dgdxi2) dcdx_.back().push_back(dcdt2*dt2dx);
		dcdt = na * c * st1 * st2;
		for(auto &dtdx:dgdxi) dcdx_.back().push_back(dcdt*dtdx);
    }
}
/*
 *
void Atom3NNTerm::setup(const std::vector<NSPgeometry::XYZ> &crds1,
		const std::vector<NSPgeometry::XYZ> &crds2) {
	std::vector<NSPgeometry::XYZ> dd;
	std::vector<NSPgeometry::XYZ> da1;
	std::vector<NSPgeometry::XYZ> da2;
	std::vector<NSPgeometry::XYZ> dt1;
	std::vector<NSPgeometry::XYZ> dt2;
	std::vector<NSPgeometry::XYZ> dt;

	double dis=distance(crds1[0],crds2[0],&dd);
	double ang1=angle(crds1[1],crds1[0],crds2[0],&da1);
	double ang2=angle(crds2[1],crds2[0],crds1[0],&da2);
	double dih1=torsion(crds1[2],crds1[1],crds1[0],crds2[0],&dt1);
	double dih2=torsion(crds2[2],crds2[1],crds2[0],crds1[0],&dt2);
	double dih=torsion(crds1[1],crds1[0],crds2[0],crds2[1],&dt);

    std::vector<std::vector<DvDxi>> dgdxi(6);
	dgdxi[0].push_back(DvDxi(0, dd[0]));//0
	dgdxi[0].push_back(DvDxi(3, dd[1]));
	dgdxi[1].push_back(DvDxi(1, da1[0]));//1
	dgdxi[1].push_back(DvDxi(0, da1[1]));
	dgdxi[1].push_back(DvDxi(3, da1[2]));
	dgdxi[2].push_back(DvDxi(4, da2[0]));//2
	dgdxi[2].push_back(DvDxi(3, da2[1]));
	dgdxi[2].push_back(DvDxi(0, da2[2]));
	dgdxi[3].push_back(DvDxi(2, dt1[0]));//3
	dgdxi[3].push_back(DvDxi(1, dt1[1]));
	dgdxi[3].push_back(DvDxi(0, dt1[2]));
	dgdxi[3].push_back(DvDxi(3, dt1[3]));
	dgdxi[4].push_back(DvDxi(5, dt2[0]));//4
	dgdxi[4].push_back(DvDxi(4, dt2[1]));
	dgdxi[4].push_back(DvDxi(3, dt2[2]));
	dgdxi[4].push_back(DvDxi(0, dt2[3]));
	dgdxi[5].push_back(DvDxi(1, dt[0]));//5
	dgdxi[5].push_back(DvDxi(0, dt[1]));
	dgdxi[5].push_back(DvDxi(3, dt[2]));
	dgdxi[5].push_back(DvDxi(4, dt[3]));

    codevec_.clear();
    dcdx_.clear();
	encode_dis(dis,dgdxi[0]);
    encode_ang(ang1,dgdxi[1]);
    encode_ang(ang2,dgdxi[2]);
    encode_dih(ang1,dgdxi[1],dih1,dgdxi[3]);
    encode_dih(ang2,dgdxi[2],dih2,dgdxi[4]);
    encode_dih2(ang1,dgdxi[1],ang2,dgdxi[2],dih,dgdxi[5]);
}
*/
void Atom3NNTerm::setup(const std::vector<NSPgeometry::XYZ> &crds1,
		const std::vector<NSPgeometry::XYZ> &crds2) {
	codevec_.clear();
	dcdx_.clear();
	std::vector<int> v1 = GetCrdIndex(names1);
	std::vector<int> v2 = GetCrdIndex(names2);
	int n1=0;
	for(int i:v1) {
		int n2=0;
		for(int j:v2) {
			std::vector<NSPgeometry::XYZ> dd;
			double dis=distance(crds1[i],crds2[j],&dd);
			std::vector<DvDxi> dgdxi;
			dgdxi.push_back(DvDxi(n1, dd[0]));
			dwdr_.push_back(dgdxi.back());
			dgdxi.push_back(DvDxi(n2, dd[1]));
			dwdr_.push_back(dgdxi.back());
			encode_dis(dis,dgdxi);
			n2++;
		}
		n1++;
	}
}

//ns1,ns2: first is resname, then list of atomname
void Atom3NNTerm::setup(const std::vector<std::string>&ns1,const std::vector<NSPgeometry::XYZ> &crds1,
		const std::vector<std::string>&ns2,const std::vector<NSPgeometry::XYZ> &crds2) {
	names1=ns1;
	names2=ns2;
	changeorder(crds1,crds2);
	if(chord_) {
		names1=ns2;
		names2=ns1;
		setup(crds2,crds1);
	} else {
		setup(crds1,crds2);
	}
	std::vector<int> v1 = GetDisIndex(ns1);
	std::vector<int> v2 = GetDisIndex(ns2);
	double mind=10000.0;
	for(int i:v1) {
		for(int j:v2) {
			double d2=(crds1[i]-crds2[j]).squarednorm();
			if(d2<mind) mind=d2;
		}
	}
	mind = sqrt(mind);
	w5565_ = weight5565(mind,dw5565_);
}

double Atom3NNTerm::outvalue(std::vector<DvDxi> *dvdxi) {
	if (codevec_.empty())
		return 0.0;
	//std::vector<std::string> n1=names1;
	//std::vector<std::string> n2=names2;
	//if(chord_) {
	//	n1=names2;
	//	n2=names1;
	//}
	auto & estimator = getestimator(names1,names2);
	std::vector<double> dvdc;
	double res = estimator(codevec_, dvdc);
	dvdxi->clear();
	for (int i = 0; i < dvdc.size(); ++i) {
		for (int j = 0; j < dcdx_[i].size(); ++j) {
			dvdxi->push_back(dvdc[i] * dcdx_[i][j]);
		}
	}
	return res;
}

std::vector<int> Atom3NNTerm::GetDisIndex(const std::vector<std::string>&ns) {
	static std::map<std::vector<std::string>,std::vector<int>> idx;
	static std::map<std::string,std::vector<std::vector<std::string>>> aall=NSPallatom::AA_Atom::groups(0);
	static std::map<std::string,std::vector<std::vector<std::string>>> adis_temp=NSPallatom::AA_Atom::groups(2);
#pragma omp critical
	{
		if(idx.empty()) {
			std::map<std::string,std::vector<std::set<std::string>>> adis;
			for(auto &p:adis_temp) {
				std::vector<std::set<std::string>> vs;
				for(auto &v:p.second) {
					std::set<std::string> s;
					for(auto &ss:v) s.insert(ss);
					vs.push_back(s);
				}
				adis.insert({p.first,vs});
			}
			for(auto &p:aall) {
				for(int i=0;i<p.second.size();i++) {
					std::vector<std::string> v{p.first};
					std::vector<int> vi;
					auto &ad=adis.at(p.first)[i];
					for(int j=0;j<p.second[i].size();j++) {
						v.push_back(p.second[i][j]);
						if(ad.find(p.second[i][j])!=ad.end()) vi.push_back(j);
					}
					idx.insert({v,vi});
				}
			}
		}
	}
	/*std::cout <<idx.size() <<std::endl;
	for(auto &p:idx) {
		for(auto &s:p.first) std::cout <<s <<'\t';
		for(auto &s:p.second) std::cout <<s <<'\t';
		std::cout <<std::endl;
	}
	exit(1);*/
	return idx.at(ns);
}
std::vector<int> Atom3NNTerm::GetCrdIndex(const std::vector<std::string>&ns) {
	static std::map<std::vector<std::string>,std::vector<int>> idx;
	static std::map<std::string,std::vector<std::vector<std::string>>> aall=NSPallatom::AA_Atom::groups(0);
	static std::map<std::string,std::vector<std::vector<std::string>>> adis_temp=NSPallatom::AA_Atom::groups(1);
#pragma omp critical
	{
		if(idx.empty()) {
			std::map<std::string,std::vector<std::set<std::string>>> adis;
			for(auto &p:adis_temp) {
				std::vector<std::set<std::string>> vs;
				for(auto &v:p.second) {
					std::set<std::string> s;
					for(auto &ss:v) s.insert(ss);
					vs.push_back(s);
				}
				adis.insert({p.first,vs});
			}
			for(auto &p:aall) {
				for(int i=0;i<p.second.size();i++) {
					std::vector<std::string> v{p.first};
					std::vector<int> vi;
					auto &ad=adis.at(p.first)[i];
					for(int j=0;j<p.second[i].size();j++) {
						v.push_back(p.second[i][j]);
						if(ad.find(p.second[i][j])!=ad.end()) vi.push_back(j);
					}
					assert(vi.size()==3);
					idx.insert({v,vi});
				}
			}
		}
	}
	/*std::cout <<idx.size() <<std::endl;
	for(auto &p:idx) {
		for(auto &s:p.first) std::cout <<s <<'\t';
		for(auto &s:p.second) std::cout <<s <<'\t';
		std::cout <<std::endl;
	}
	exit(1);*/
	return idx.at(ns);
}

//NM
double Atom3NNTerm::weight5565(double d, double &dd) {
	d*=10.0;
	if(d>6.5) {
		dd=0.0;
		return 0.0;
	}
	if(d<5.5) {
		dd=0.0;
		return 1.0;
	}
	double delta;
	if(d<6.0) delta=5.5-d;
	else delta=d-6.5;
	dd = 4.0*delta;
	delta=delta*delta*2.0;
	if(d>6.0) return delta;
	else return 1.0-delta;
}

































