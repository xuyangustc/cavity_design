/*
 * forcefield_omp.cpp
 *
 *  Created on: 2018骞�4鏈�12鏃�
 *      Author: hyliu
 */
#ifdef _OPENMP
#include "sd/forcefield.h"
#include "sd/genchain.h"
#include "sd/activeselections.h"
#include "sd/nnterm.h"
#include "sd/localhbterm.h"
#include "sd/backboneff.h"
#include "fullsite/structplus.h"
#include "sd/shadowterm.h"
#include <cassert>
#include <iomanip>
#include <sstream>
#include <omp.h>
#include <algorithm>
using namespace NSPsd;
using namespace NSPgeometry;
GroupNeighbor::GroupNeighbor(const std::vector<double> &crd, const ForceField &ff) {
	double rcut=0.7; //NM
	double rcut2=rcut*rcut;
	int nmainchain=5;
	int nbond=6;
	std::vector<XYZ> xyz;
	for (int i = 0; i < crd.size(); i += 3) {
		xyz.push_back(XYZ(crd[i], crd[i + 1], crd[i + 2]));
	}
	const auto cgs=ff.cgps();
#pragma omp parallel for schedule(dynamic, 1)
	for(int i=0;i<cgs.size();i++) {
		for(int j=i+1;j<cgs.size();j++) {
			if(cgs[i].chainseq==cgs[j].chainseq&&cgs[i].resseq==cgs[j].resseq) continue;
			if(cgs[i].chainseq==cgs[j].chainseq) {
				int nm=fabs(cgs[i].resseq==cgs[j].resseq);
				if(nm<nmainchain) {
					if(!NSPallatom::GroupPair::farinteraction(cgs[i],cgs[j],nbond)) continue;
				}
			}
			bool isngh{false};
			for(int h=0;h<cgs[i].atomseqs.size();h++) {
				for(int k=0;k<cgs[j].atomseqs.size();k++) {
					XYZ c1=xyz[cgs[i].atomseqs[h]];
					XYZ c2=xyz[cgs[j].atomseqs[k]];
					if((c1-c2).squarednorm()>rcut2) continue;
					isngh=true;
					break;
				}
				if(isngh) break;
			}
			if(!isngh) continue;
#pragma omp critical
			{
				neighs.push_back({i,j});
			}
		}
	}
	/*for(int i=0;i<neighs.size();i++) {
		int m=neighs[i][0];
		int n=neighs[i][1];
		std::cout <<i <<'\t' <<m <<' ' <<n <<std::endl;
		std::cout <<cgs[m].chainseq <<' ' <<cgs[m].resseq <<std::endl;
		for(auto &p:cgs[m].crds_v) {
			std::cout <<' ' <<p.first <<' ' <<p.second.x_ <<' ' <<p.second.y_ <<' ' <<p.second.z_ <<std::endl;
		}
		std::cout <<cgs[n].chainseq <<' ' <<cgs[n].resseq <<std::endl;
		for(auto &p:cgs[n].crds_v) {
			std::cout <<' ' <<p.first <<' ' <<p.second.x_ <<' ' <<p.second.y_ <<' ' <<p.second.z_ <<std::endl;
		}
		std::cout <<std::endl;
	}*/
	//exit(1);
}
NeighborList::NeighborList(const std::vector<double> &crd, const ForceField &ff,
		const std::vector<bool> *forceoff) {
	std::vector<XYZ> xyz;
	for (int i = 0; i < crd.size(); i += 3) {
		xyz.push_back(XYZ(crd[i], crd[i + 1], crd[i + 2]));
	}
	int natoms = crd.size() / 3;
	int nthreads=omp_get_max_threads();
	std::vector<std::vector<std::vector<int>>> neighbors_thread(nthreads,
			std::vector<std::vector<int>>(natoms));

	double rcut2 = 1.00; //nm^2
	if(ff.scmmsteric()){
#pragma omp parallel
		{
			int myid=omp_get_thread_num();
			std::vector<std::vector<int>> & nbrs_t=neighbors_thread[myid];
#pragma omp for schedule(dynamic,1)
		for (int i = 0; i < natoms; ++i) {

			std::vector<int> & nbrsi = nbrs_t[i];
			if (forceoff) {
				if (forceoff->at(i))
					continue;
			}
			for (int j = 0; j < natoms; ++j) {
				if (i == j)
					continue;
				if (forceoff) {
					if (!forceoff->at(j) && (j < i))
						continue;
				} else if (j < i)
					continue;
				if (fabs(j - i) < MAXEXCL) {
					if (ff.excluded(i, j))
						continue;
					if (ff.is14(i, j))
						continue;
				}
				double diff2 = (xyz[i] - xyz[j]).squarednorm();
				if (diff2 < rcut2) {
					nbrsi.push_back(j);
					nbrs_t[j].push_back(i);
				}//rcut2
			} // j
		} //i
		}//parallel
	} else {
#pragma omp parallel
		{
			int myid=omp_get_thread_num();
			auto & nbrs_t=neighbors_thread[myid];
#pragma omp for schedule(dynamic,1)
		for(int m=0;m<ff.bsinchains().size();++m){
			const std::vector<BSInChain> & bsm=ff.bsinchains()[m];
			for(int i=0;i<bsm.size();++i){
				std::vector<int> ia=bsm[i].atomids();
				for(int n=m;n<ff.bsinchains().size();++n){
					const std::vector<BSInChain> &bsn=ff.bsinchains()[n];
					int jstart=0;
					int jend=bsn.size();
					if(n==m) {
						jstart=i+1;
						if(jstart >=jend) continue;
					}
					for(int j=jstart;j<jend;++j){
						std::vector<int> ja=bsn[j].atomids();
						for(auto ii:ia){
							for(auto jj:ja){
								if(forceoff)
									if(forceoff->at(ii) && forceoff->at(jj))continue;
								if(fabs(jj-ii) <MAXEXCL){
									if(ff.excluded(ii,jj)) continue;
									if(ff.is14(ii,jj)) continue;
								}
								double diff2=(xyz[ii]-xyz[jj]).squarednorm();
								if(diff2<rcut2){
									nbrs_t[ii].push_back(jj);
									nbrs_t[jj].push_back(ii);
								} //rcut2
							} //ja
						}//ia
					} //j
				} //n
			} //i
		} //m
		}//parallel
	} //scmmsteric
	neighbors.assign(natoms, std::vector<int>());
	for(auto &nbrs_t:neighbors_thread){
		for(int i=0;i<natoms;++i){
			int osize=neighbors[i].size();
			int nsize=osize+nbrs_t[i].size();
			neighbors[i].resize(nsize);
			std::copy(nbrs_t[i].begin(),nbrs_t[i].end(),neighbors[i].begin()+osize);
		}
	}
}
std::vector<double> ForceField::sideconfene(const std::vector<double> &crd,
		const std::vector<std::vector<PhiPsiCodes>> &phipsicodes,
		std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
		const ActiveSelections &acts) const {
	const std::vector<std::vector<ConformerCode>> & sccodes=
			acts.conformercodes();
	//std::vector<std::vector<ConformerCode>> sccodes=acts.conformercodes();
	int nchains = phipsicodes.size();
/*	for (int c = 0; c < nchains; ++c) {
		sccodes.push_back(
				makeconformercodes(crd, scinchains_[c], phipsicodes[c]));
	}*/
	int nthreads=omp_get_max_threads();
	std::vector<double> ene_thread(nthreads,0);
	int poffset=0;
	for (int c = 0; c < nchains; ++c) {
#pragma omp parallel
		{
			int myid=omp_get_thread_num();
			double &ene=ene_thread[myid];
			std::vector<NSPgeometry::XYZ> *xyzf=&(xyzf_thread[myid]);
#pragma omp for schedule(dynamic,1)
		for (int p = 0; p < sccodes[c].size(); ++p) {
			if(!acts.rotameractive(c,p)) continue;
			if (scinchains_[c][p].kaiatoms.empty())
				continue;
//			assert(p != 0 && p != sccodes[c].size() - 1); //current code does not yet handle terminals properly
			std::vector<DvDxi> dedx;
			double e;
			if(p==0 || p==sccodes[c].size()-1){
				NN_KaiTerm_T kaiterm(sccodes[c][p]);
				e=kaiterm.outvalue(&dedx);
			} else {
				NN_KaiTerm kaiterm(sccodes[c][p]);
				e = kaiterm.outvalue(&dedx);
			}
			ene += e;
			for (auto &d : dedx) {
				xyzf->at(d.first) = xyzf->at(d.first) - scconfwght_ * d.second;
			}
			if (eneanalysismode_) {
#pragma omp critical(eneanalysis)
				{
					enecomp_.escconf[poffset+p] = e * scconfwght_;
				}
			}
		}
		}//parallel
		poffset += sccodes[c].size();
	}
	for(auto &ene:ene_thread) ene *= scconfwght_;
	return ene_thread;
}
std::vector<double> ForceField::scpackingene(const std::vector<NSPgeometry::XYZ> &xyz,
		std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
		const ActiveSelections & acts) const {
	int nchains = bsinchains_.size();
	int minsep=1;
	int nthreads=omp_get_max_threads();
	std::vector<double> esum_thread (nthreads,0.0);
	double rcut2=1.21;
	for (int n = 0; n < nchains; ++n) {
		const std::vector<BSInChain> &bsn = bsinchains_[n];
		for (int m = n; m < nchains; ++m) {
			const std::vector<BSInChain> &bsm = bsinchains_[m];
#pragma omp parallel
		{
			int myid=omp_get_thread_num();
			double &esum=esum_thread[myid];
			std::vector<NSPgeometry::XYZ> *xyzf=&(xyzf_thread[myid]);
#pragma omp for schedule(dynamic,1)
			for (int i = 0; i < bsn.size(); ++i) {
				int jstart = 0;
				int jend = bsm.size();
				if (m == n) {
					jstart = i + minsep;
				}
				bool isc = (scinchains_[n][i].nscatoms >0);
				if (jstart > jend)
					continue;

				for (int j = jstart; j < jend; ++j) {
					if(!acts.sidechainpairactive(n,i,m,j)) continue;
					bool jsc = (scinchains_[m][j].nscatoms > 0);
					double rca2=(xyz[bsn[i].caid]-xyz[bsm[j].caid]).squarednorm();
					if(rca2>rcut2) continue;
					if (!(isc || jsc))
						continue;
					double sw =
							siteweights_[n][i] < siteweights_[m][j] ?
									siteweights_[n][i] : siteweights_[m][j];
					if (sw == 0.0)
						continue;
					sw*=scpackingweight_;
					int sep=1000;
					if(m==n){
						sep=j-i;
					}
					double esc=0.0;
					if(jsc){
						std::vector<DvDxi> dedx;
						bool terminal=(j==jend);
						esc +=mcscpackingenergy(xyz, stericatomtypes_,
								bsn[i],scinchains_[m][j],
								sep,&dedx,terminal);
						for(auto &d:dedx)
							(*xyzf)[d.first]=(*xyzf)[d.first]-sw*d.second;
					}
					if(isc){
						std::vector<DvDxi> dedx;
						bool terminal=(i==0);
						esc +=mcscpackingenergy(xyz, stericatomtypes_,
								bsm[j],scinchains_[n][i],
								-sep,&dedx,terminal);
						for(auto &d:dedx)
							(*xyzf)[d.first]=(*xyzf)[d.first]-sw*d.second;
					}
					if (isc && jsc) {
						std::vector<DvDxi> dedx;
						esc +=scscpackingenergy(xyz,stericatomtypes_,
								scinchains_[n][i],scinchains_[m][j],sep,&dedx);
						for (auto &d : dedx)
							(*xyzf)[d.first] = (*xyzf)[d.first] - sw * d.second;
					}
					esum += esc * sw;
				}
			}
		}//parallel
		}
	}
	return esum_thread;
}
std::vector<double> ForceField::sitepairenergy(const std::vector<double> & crd,
		const std::vector<std::vector<PhiPsiCodes>> & phipsicodes,
		const std::vector<std::vector<SSCode>> &sscodes,
		std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
		std::vector<double> &eattr_thread,
		const ActiveSelections & acts) const {
	int nthreads=omp_get_max_threads();
	std::vector<double> ene_thread(nthreads,0.0);
	eattr_thread.assign(nthreads, 0.0);
	int nchains = bsinchains_.size();
#ifdef DOCB
	std::vector<std::vector<CBData>> cbdata;
	for(int n=0;n<nchains;++n) {
		const std::vector<BSInChain> &bsn=bsinchains_[n];
		cbdata.push_back(std::vector<CBData>());
		std::vector<CBData> &ccb=cbdata.back();
		ccb.resize(bsn.size());
#pragma omp parallel for schedule(dynamic,1)
		for (int m=0;m<bsn.size();++m) {
			std::vector<NSPgeometry::XYZ> crdncac;
			std::vector<int> index = bsn[m].atomids();
			for (auto id : index)
			crdncac.push_back(getxyz(crd, id));
//test
			/*			CBData cbt(crdncac);
			 for(int ii=0;ii<3;++ii){
			 std::vector<NSPgeometry::XYZ> dcbdi(3);
			 for(int jj=0;jj<3;++jj){
			 crdncac[ii][jj] +=0.0001;
			 NSPgeometry::XYZ cb1=CBData(crdncac).cbcrd();
			 crdncac[ii][jj] -=0.0002;
			 NSPgeometry::XYZ cb2=CBData(crdncac).cbcrd();
			 dcbdi[jj]=(cb1-cb2)/0.0002;
			 crdncac[ii][jj] +=0.0001;
			 std::cout <<dcbdi[jj].toString() <<" "<< cbt.drcb()[ii][jj].toString() <<std::endl;
			 }
			 }
			 exit(1);*/
			ccb[m]=CBData(crdncac);
		}
	}
#endif
	for (int n = 0; n < nchains; ++n) {
		const std::vector<BSInChain> &bsn = bsinchains_[n];
		for (int m = n; m < nchains; ++m) {
			const std::vector<BSInChain> &bsm = bsinchains_[m];
#pragma omp parallel
			{
			int myid=omp_get_thread_num();
			double &ene=ene_thread[myid];
			double *eattr=&(eattr_thread[myid]);
			std::vector<NSPgeometry::XYZ> *xyzf=&(xyzf_thread[myid]);
#pragma omp for schedule(dynamic,1)
			for (int i = 1; i < bsn.size() - 1; ++i) {
				int jstart = 1;
				int jend = bsm.size() - 1;
				if (m == n) {
					jstart = i + 1;
				}

				for (int j = jstart; j < jend; ++j) {
					if(!acts.sitepairactive(n,i,m,j)) continue;
					if (m == n) {
						if (sitepairoff(sscodes[n], i, j))
							continue;
					}
					SitePairNNTerm spterm;
#ifdef USE_CB
					spterm.setupwithcbdata(crd,bsn,cbdata[n][i],phipsicodes[n],i,bsm,
							cbdata[m][j],phipsicodes[m],j);
#else
					spterm.setup(crd, bsn, phipsicodes[n], i, bsm,
							phipsicodes[m], j);
#endif
					double rca = spterm.rca();
					const std::vector<DvDxi> &drcadx = spterm.drcadx();
					double detrapdr;
					double wij = 1.0;
					if (!coreweights_.empty()) {
						wij = coreweights_[bsn[i].caid]
								* coreweights_[bsm[j].caid];
					}
					*eattr += wij * eattract_ * etrap(rca, &detrapdr);
					detrapdr *= wij * eattract_;
					for (auto &d : drcadx) {
						xyzf->at(d.first) = xyzf->at(d.first)
								- detrapdr * d.second;
					}
#ifdef CBPAIRENE
					std::vector<NSPgeometry::XYZ> drdxcb;
					double rcb=NSPgeometry::distance(cbdata[n][i].cbcrd(),
							cbdata[m][j].cbcrd(),&drdxcb);
					double dedrcb=0.0;
//					double ecbp=CBPairEne::cbforce(rcb+0.0001,&dedrcb);
//					double ecbm=CBPairEne::cbforce(rcb-0.0001,&dedrcb);
//					std::cout <<rcb<<" "<<(ecbp-ecbm)/0.0002 << " "<< dedrcb<<std::endl;
					double ecb=KBT*CBPairEne::cbforce(rcb,&dedrcb);
					/*					double ecb;
					 if (rcb>0.5){
					 ecb=KBT*0.5*100*(rcb-0.5)*(rcb-0.5);
					 dedrcb=100*(rcb-0.5);
					 }*/
//					if((ecbp-ecbm)/0.0002 !=dedrcb){ //for debug
//						ecb=0.0;
//						dedrcb=0.0;
//					}
					dedrcb *=KBT;
					if(dedrcb !=0.0 ) {
						std::vector<NSPgeometry::XYZ> f1=cbdata[n][i].
						distributederiv(dedrcb*drdxcb[0]);
						std::vector<NSPgeometry::XYZ> f2=cbdata[m][j].
						distributederiv(dedrcb*drdxcb[1]);
						std::vector<int> index1=bsn[i].atomids();
						std::vector<int> index2=bsm[j].atomids();
						for(int ii=0;ii<3;++ii) {
							xyzf->at(index1[ii]) =xyzf->at(index1[ii])-f1[ii];
							xyzf->at(index2[ii]) =xyzf->at(index2[ii])-f2[ii];
						}
						if(eneanalysismode_) {
#pragma omp critical(eneanalysis)
							{
								std::cout<<n <<":"<<i<<"-"<<m<<":"<<j<<" rcb "<<rcb<<" ecb "<<ecb<<std::endl;
							}
						}
					}
					ene +=ecb/sitepairweight_;
#endif
					double sw =
							siteweights_[n][i] < siteweights_[m][j] ?
									siteweights_[n][i] : siteweights_[m][j];
					std::vector<DvDxi> dedx;
					double e = spterm.outvalue(&dedx);

					ene += e * sw;
					if (eneanalysismode_) {
						e *= 0.5 * sitepairweight_;
						int ni = siteid(n, i);
						int mj = siteid(m, j);
#pragma omp critical(eneanalysis)
						{
							enecomp_.esitepair[ni][mj] = e;
							enecomp_.esitepair[mj][ni] = e;
							enecomp_.epacking[ni] += e;
							enecomp_.epacking[mj] += e;
						}
					}
					sw *= sitepairweight_;
					for (auto &d : dedx) {
						xyzf->at(d.first) = xyzf->at(d.first) - sw * d.second;
					}
				}
			}
			}//parallel
		}
	}
	for(auto & ene:ene_thread) ene *=sitepairweight_;
	return ene_thread;
}
std::vector<double> ForceField::stericenergy(const NeighborList &nbl,
		const std::vector<XYZ> &xyz, std::vector<std::vector<NSPgeometry::XYZ>>
		&xyzf_thread) const {

	int nthread=omp_get_max_threads();
	std::vector<double> ene_thread(nthread, 0.0);
#pragma omp parallel
	{
	std::vector<XYZ> pairf(2);
	int myid=omp_get_thread_num();
	double &ene=ene_thread[myid];
	std::vector<NSPgeometry::XYZ> *xyzf=&(xyzf_thread[myid]);
//	std::cout <<"NThreads: "<< omp_get_num_threads()<<std::endl;
#pragma omp for schedule(dynamic,1)
	for (int i = 0; i < natoms_; ++i) {
		for (auto j : nbl.neighbors[i]) {
			if (j <= i)
				continue;
			ene += StericAtom::pairenergy(xyz[i], xyz[j], stericatoms_[i],
					stericatoms_[j], &pairf);
			(*xyzf)[i] = (*xyzf)[i] + stericwght_ * pairf[0];
			(*xyzf)[j] = (*xyzf)[j] + stericwght_ * pairf[1];
		}
	}
	} //end parallel
	for(auto & ene:ene_thread) ene *=stericwght_;
	return ene_thread;
}
std::vector<double> ForceField::lsenergy(const std::vector<PhiPsiCodes> & phipsicodes,
		const std::vector<double> &siteweight,
		std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
		const std::vector<bool> &lsactive) const {
	int nthread=omp_get_max_threads();
	std::vector<double> els_thread(nthread,0.0);
#pragma omp parallel
	{
		int myid=omp_get_thread_num();
		double &els=els_thread[myid];
		std::vector<NSPgeometry::XYZ> *xyzf=&(xyzf_thread[myid]);
#pragma omp for schedule(dynamic,1)
	for (int i = 2; i < phipsicodes.size() - 2; ++i) {
		if(!lsactive[i]) continue;
		LSNNTerm lsnnterm;
		std::vector<DvDxi> dedx;
		lsnnterm.setup(phipsicodes, i);
		double w = siteweight[i] * lsweight_;
		double e = lsnnterm.outvalue(&dedx);
//		if (e > 0.0)
//			w = w * 4.0;
		els += e * w;
		if (eneanalysismode_) {
			enecomp_.els[enecomp_.siteoffset + i] = e * w;
		}
		for (auto &d : dedx) {
			(*xyzf)[d.first] = (*xyzf)[d.first] - w * d.second;
		}
	}
	}//end parallel
	return els_thread;
}
std::vector<double> ForceField::scpackingeneNN(const std::vector<NSPgeometry::XYZ> &xyz,
		std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
		const GroupNeighbor &gn) const {
	int nthreads=omp_get_max_threads();
	std::vector<double> esum_thread (nthreads,0.0);
	std::map<std::vector<int>,double> es_record;
#pragma omp parallel
	{
		int myid=omp_get_thread_num();
		double &esum=esum_thread[myid];
		std::vector<NSPgeometry::XYZ> *xyzf=&(xyzf_thread[myid]);
#pragma omp for schedule(dynamic,1)
		for(int i=0;i<gn.neighs.size();i++) {
			NSPsd::Atom3NNTerm a3n;
			std::vector<std::string> n1=chemgroups_[gn.neighs[i][0]].makelabel();
			std::vector<std::string> n2=chemgroups_[gn.neighs[i][1]].makelabel();
			std::vector<XYZ> c1,c2;
			for(int j:chemgroups_[gn.neighs[i][0]].atomseqs) c1.push_back(xyz[j]);
			for(int j:chemgroups_[gn.neighs[i][1]].atomseqs) c2.push_back(xyz[j]);
			a3n.setup(n1,c1,n2,c2);
			std::vector<DvDxi> dvdx;
			std::vector<int> idx;
			for(int j:chemgroups_[gn.neighs[i][0]].atomseqs) idx.push_back(j);
			for(int j:chemgroups_[gn.neighs[i][1]].atomseqs) idx.push_back(j);
			double e = a3n.outvalue(&dvdx,idx);
			//esum += e*scpackingweight_;
			esum += e*lsweight_;
			for (auto &d : dvdx) {
				//xyzf->at(d.first) = xyzf->at(d.first) - scpackingweight_ * d.second;
				xyzf->at(d.first) = xyzf->at(d.first) - lsweight_ * d.second;
			}
/*#pragma omp critical
			{
				//std::cout <<"NN" <<std::endl;
				int m=gn.neighs[i][0];
				int n=gn.neighs[i][1];
				if(m==0 && n==109) {
					for(auto n:n1) std::cout <<' ' <<n; std::cout <<std::endl;
					for(auto c:c1) std::cout <<' ' <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_; std::cout <<std::endl;
					for(auto n:n2) std::cout <<' ' <<n; std::cout <<std::endl;
					for(auto c:c2) std::cout <<' ' <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_; std::cout <<std::endl;
				}
				es_record.insert({{m,n},e});
				std::cout <<i <<'\t' <<m <<' ' <<n <<std::endl;
				std::cout <<e <<std::endl;
				//exit(1);
				std::cout <<chemgroups_[m].chainseq <<' ' <<chemgroups_[m].resseq <<' ' <<chemgroups_[m].resname <<std::endl;
				for(auto &p:chemgroups_[m].crds_v) {
					std::cout <<' ' <<p.first <<' ' <<p.second.x_ <<' ' <<p.second.y_ <<' ' <<p.second.z_ <<std::endl;
				}
				for(int j:chemgroups_[m].atomseqs) std::cout <<' ' <<j;
				std::cout <<std::endl;
				for(auto &n:n1) std::cout <<' ' <<n <<std::endl;
				for(XYZ &c:c1) {
					std::cout <<' ' <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<std::endl;
				}
				std::cout <<std::endl;

				std::cout <<chemgroups_[n].chainseq <<' ' <<chemgroups_[n].resseq <<' ' <<chemgroups_[n].resname <<std::endl;
				for(auto &p:chemgroups_[n].crds_v) {
					std::cout <<' ' <<p.first <<' ' <<p.second.x_ <<' ' <<p.second.y_ <<' ' <<p.second.z_ <<std::endl;
				}
				for(int j:chemgroups_[n].atomseqs) std::cout <<' ' <<j;
				std::cout <<std::endl;
				for(auto &n:n2) std::cout <<' ' <<n <<std::endl;
				for(XYZ &c:c2) {
					std::cout <<' ' <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<std::endl;
				}
				exit(1);
			}*/
		}
	}
	//for(auto &p:es_record) {
	//	std::cout <<p.first[0] <<' ' <<p.first[1] <<'\t' <<p.second <<std::endl;
	//}
	//exit(1);
	//for(auto &p:xyzf_thread) std::cout <<p.size() <<std::endl;
	/*std::vector<NSPgeometry::XYZ> fs=xyzf_thread[0];
	for(int i=1;i<xyzf_thread.size();i++) {
		for(int j=0;j<fs.size();j++) {
			fs[j] = fs[j] + xyzf_thread[i][j];
		}
	}
	for(auto &c:fs) {
		std::cout <<c.x_ <<'\t' <<c.y_ <<'\t' <<c.z_ <<std::endl;
	}
	for(double d:esum_thread) std::cout <<d <<std::endl;
	double t=esum_thread[0];
	for(int i=1;i<esum_thread.size();i++) t += esum_thread[i];
	std::cout <<t <<std::endl;
	exit(1);*/
	return esum_thread;
}
std::vector<double> ForceField::forces(const std::vector<double> &crd,
		const NeighborList &nbl, const GroupNeighbor &gn, std::vector<double> *potenergies,
//		const std::vector<bool> *forceoff) const {
		ActiveSelections &acts) const{
	if (eneanalysismode_) {
		int nsites = 0;
		for (auto &bs : bsinchains_) {
			nsites += bs.size();
		}
		enecomp_.init(nsites);
	}
	std::vector<XYZ> xyz;
	assert(natoms_ == crd.size() / 3);
	for (int i = 0; i < crd.size(); i += 3) {
		xyz.push_back(XYZ(crd[i], crd[i + 1], crd[i + 2]));
	}
	std::vector<XYZ> xyzf(natoms_, XYZ(0.0, 0.0, 0.0));
	int nthread=omp_get_max_threads();
	std::vector<std::vector<XYZ>> xyzf_thread(nthread,xyzf);
	potenergies->assign(ETERMS, 0.0);
	std::vector<bool> forceoff=acts.atomfixed();
	std::vector<std::vector<double>> pot_thread(ETERMS, std::vector<double>(nthread,0.0));
	pot_thread.at(EBOND) = covalentenergy(bondterms_, xyz, xyzf_thread, &forceoff);
	pot_thread.at(EANG) = covalentenergy(angleterms_, xyz, xyzf_thread, &forceoff);
	pot_thread.at(EIMPDIH) = covalentenergy(impdihterms_, xyz, xyzf_thread,
			&forceoff);
	int nchains = bsinchains_.size();
	/*std::vector<std::vector<PhiPsiCodes>> phipsicodes;
	for (int i = 0; i < nchains; ++i) {
		phipsicodes.push_back(makephipsicodes(crd, bsinchains_[i]));
	}
	std::vector<std::vector<SSCode>> sscodes;
	for (int i = 0; i < nchains; ++i) {
		sscodes.push_back(estimatess(phipsicodes[i]));
	}*/
	acts.calccodes(crd);
	const std::vector<std::vector<PhiPsiCodes>> & phipsicodes=acts.phipsicodes();
	const std::vector<std::vector<SSCode>> & sscodes=acts.sscodes();
	int nps=phipsis_.size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<nps;++i){
//	for (auto it = phipsis_.begin(); it != phipsis_.end(); ++it) {
		auto it=phipsis_.begin()+i;
		if(!acts.phipsiactive(it->chainid,it->siteid)) {
			continue;
		}
		double w = phipsiwght_ * siteweights_[it->chainid][it->siteid];
		int myid=omp_get_thread_num();
		double e = it->energy(phipsicodes[it->chainid], w, &(xyzf_thread[myid]));
		pot_thread.at(EPHIPSI)[myid] += e;
		if (eneanalysismode_) {
#pragma omp critical(eneanalysis)
			{
				enecomp_.ephipsi[i] = e;
			}
		}
	}
	if (lstermon_) {
		//pot_thread.at(ELOCALSTRUCTURE) = scpackingeneNN(xyz,xyzf_thread,gn);
		enecomp_.siteoffset = 0;
		for (int i = 0; i < nchains; ++i) {
			pot_thread.at(ELOCALSTRUCTURE) = lsenergy(phipsicodes[i],
					siteweights_[i], xyzf_thread,acts.lsactive(i));
			enecomp_.siteoffset += bsinchains_.size();
		}
	}
	if (sitepairtermon_) {
		std::vector<double> etrap_thread;
		pot_thread.at(ESITEPAIRS) = sitepairenergy(crd, phipsicodes, sscodes,
				xyzf_thread, etrap_thread,acts);
		for(int i=0;i<nthread;++i){
			pot_thread.at(ESTRUCTREST)[i] += etrap_thread[i];
		}
	}
	if (scconfwght_ != 0.0) {
		pot_thread.at(ESCCONF) = sideconfene(crd, phipsicodes, xyzf_thread,acts);
	}
	if (shadowweight_ != 0.0) {
		potenergies->at(ESHADOW) = shadowenergy(xyz, &xyzf);
	}
	for (auto it = structrestraints_.begin(); it != structrestraints_.end();
			++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	for (auto it = disrestraints_.begin(); it != disrestraints_.end(); ++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	potenergies->at(ERGRESTRAINT) += rgrestraint_.energy(xyz, &xyzf,
			coreweights_);
	for (auto it = ssrestraints_.begin(); it != ssrestraints_.end(); ++it) {
		potenergies->at(ESSRESTRAINT) += it->energy(sscodes, &xyzf);
	}
	if (stericwght_ > 0.0)
		pot_thread.at(ESTERIC) = stericenergy(nbl, xyz, xyzf_thread);
	//clock_t t1=clock();
	if (scpackingweight_>0.0)
		pot_thread.at(ESCPACKING) = scpackingene(xyz,xyzf_thread,acts);
	//pot_thread.at(ELOCALSTRUCTURE) = scpackingeneNN(xyz,xyzf_thread,gn);
	//clock_t t2=clock();
	//std::cout <<t2-t1 <<std::endl;
	if (localhbweight_ > 0.0)
		potenergies->at(ELOCALHB) = localbbhbenergy(xyz, &xyzf);
#pragma omp parallel for schedule(dynamic,1)
	for(int t=0;t<ETERMS;++t){
		for(int i=0;i<nthread;++i){
			potenergies->at(t) += pot_thread[t][i];
		}
	}
#pragma omp parallel for schedule(dynamic,1)
	for(int a=0;a<natoms_;++a){
		NSPgeometry::XYZ &fa=xyzf.at(a);
		for(int i=0;i<nthread;++i){
			fa =fa+xyzf_thread[i][a];
		}
	}
	potenergies->at(ETOT) = 0.0;
	for (int i = 1; i < ETERMS; ++i) {
		potenergies->at(ETOT) += potenergies->at(i);
	}
	std::vector<double> forces(3*natoms_,0.0);
#pragma omp parallel for schedule(dynamic,5)
//	for (auto it = xyzf.begin(); it != xyzf.end(); ++it) {
	for(int i=0;i<natoms_;++i){
		int i3=3*i;
		forces[i3]=xyzf[i].x_;
		forces[i3+1]=xyzf[i].y_;
		forces[i3+2]=xyzf[i].z_;
	}
	/*std::vector<double> fs1;
	double t1=0;
	for(int i=0;i<natoms_;++i){
		double f=sqrt(xyzf[i].squarednorm());
		fs1.push_back(f);
		t1 += f;
	}
	t1 /= (double)fs1.size();
	std::cout <<t1 <<std::endl;
	std::sort(fs1.begin(),fs1.end());
	for(double d:fs1) std::cout <<d <<std::endl;
	exit(1);*/
	//for(double d:*potenergies) std::cout <<d <<std::endl;
	//std::cout <<std::endl;
	//for(auto &c:xyzf) {
	//	std::cout <<c.x_ <<'\t' <<c.y_ <<'\t' <<c.z_ <<std::endl;
	//}
	//exit(1);
	return forces;
}


double getf(const std::vector<double>&vd, double dis) {
	if(dis<vd[0]) return 0.0;
	if(dis>vd[1]) return vd[2];
	double v01 = (vd[0]+vd[1])/2.0;
	double a = vd[1]-vd[0];
	a *= a;
	a = 2.0 * vd[2] / a;
	if(dis>v01) {
		double delta = dis-vd[0];
		return delta*delta*a;
	} else {
		double delta = vd[1] - dis;
		return vd[2] - a * delta * delta;
	}
}
/*
 * pointinfo:
 * assemble of atoms, fixedpoint, pars(start,end,forces);
 */
std::vector<double> ForceField::forces_fixedpoint_restrict(const std::vector<double> &crd,
		const NeighborList &nbl, const GroupNeighbor &gn, std::vector<double> *potenergies,
//		const std::vector<bool> *forceoff) const {
		ActiveSelections &acts, const std::vector<std::pair<std::pair<std::vector<int>,
		XYZ>,std::vector<double>>>&pointinfo) const{
	if (eneanalysismode_) {
		int nsites = 0;
		for (auto &bs : bsinchains_) {
			nsites += bs.size();
		}
		enecomp_.init(nsites);
	}
	std::vector<XYZ> xyz;
	assert(natoms_ == crd.size() / 3);
	for (int i = 0; i < crd.size(); i += 3) {
		xyz.push_back(XYZ(crd[i], crd[i + 1], crd[i + 2]));
	}
	std::vector<XYZ> xyzf(natoms_, XYZ(0.0, 0.0, 0.0));
	int nthread=omp_get_max_threads();
	std::vector<std::vector<XYZ>> xyzf_thread(nthread,xyzf);
	potenergies->assign(ETERMS, 0.0);
	std::vector<bool> forceoff=acts.atomfixed();
	std::vector<std::vector<double>> pot_thread(ETERMS, std::vector<double>(nthread,0.0));
	pot_thread.at(EBOND) = covalentenergy(bondterms_, xyz, xyzf_thread, &forceoff);
	pot_thread.at(EANG) = covalentenergy(angleterms_, xyz, xyzf_thread, &forceoff);
	pot_thread.at(EIMPDIH) = covalentenergy(impdihterms_, xyz, xyzf_thread,
			&forceoff);
	int nchains = bsinchains_.size();
	/*std::vector<std::vector<PhiPsiCodes>> phipsicodes;
	for (int i = 0; i < nchains; ++i) {
		phipsicodes.push_back(makephipsicodes(crd, bsinchains_[i]));
	}
	std::vector<std::vector<SSCode>> sscodes;
	for (int i = 0; i < nchains; ++i) {
		sscodes.push_back(estimatess(phipsicodes[i]));
	}*/
	acts.calccodes(crd);
	const std::vector<std::vector<PhiPsiCodes>> & phipsicodes=acts.phipsicodes();
	const std::vector<std::vector<SSCode>> & sscodes=acts.sscodes();
	int nps=phipsis_.size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<nps;++i){
//	for (auto it = phipsis_.begin(); it != phipsis_.end(); ++it) {
		auto it=phipsis_.begin()+i;
		if(!acts.phipsiactive(it->chainid,it->siteid)) {
			continue;
		}
		double w = phipsiwght_ * siteweights_[it->chainid][it->siteid];
		int myid=omp_get_thread_num();
		double e = it->energy(phipsicodes[it->chainid], w, &(xyzf_thread[myid]));
		pot_thread.at(EPHIPSI)[myid] += e;
		if (eneanalysismode_) {
#pragma omp critical(eneanalysis)
			{
				enecomp_.ephipsi[i] = e;
			}
		}
	}
	if (lstermon_) {
		//pot_thread.at(ELOCALSTRUCTURE) = scpackingeneNN(xyz,xyzf_thread,gn);
		enecomp_.siteoffset = 0;
		for (int i = 0; i < nchains; ++i) {
			pot_thread.at(ELOCALSTRUCTURE) = lsenergy(phipsicodes[i],
					siteweights_[i], xyzf_thread,acts.lsactive(i));
			enecomp_.siteoffset += bsinchains_.size();
		}
	}
	if (sitepairtermon_) {
		std::vector<double> etrap_thread;
		pot_thread.at(ESITEPAIRS) = sitepairenergy(crd, phipsicodes, sscodes,
				xyzf_thread, etrap_thread,acts);
		for(int i=0;i<nthread;++i){
			pot_thread.at(ESTRUCTREST)[i] += etrap_thread[i];
		}
	}
	if (scconfwght_ != 0.0) {
		pot_thread.at(ESCCONF) = sideconfene(crd, phipsicodes, xyzf_thread,acts);
	}
	if (shadowweight_ != 0.0) {
		potenergies->at(ESHADOW) = shadowenergy(xyz, &xyzf);
	}
	for (auto it = structrestraints_.begin(); it != structrestraints_.end();
			++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	for (auto it = disrestraints_.begin(); it != disrestraints_.end(); ++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	potenergies->at(ERGRESTRAINT) += rgrestraint_.energy(xyz, &xyzf,
			coreweights_);
	for (auto it = ssrestraints_.begin(); it != ssrestraints_.end(); ++it) {
		potenergies->at(ESSRESTRAINT) += it->energy(sscodes, &xyzf);
	}
	if (stericwght_ > 0.0)
		pot_thread.at(ESTERIC) = stericenergy(nbl, xyz, xyzf_thread);
	//clock_t t1=clock();
	if (scpackingweight_>0.0)
		pot_thread.at(ESCPACKING) = scpackingene(xyz,xyzf_thread,acts);
	//pot_thread.at(ELOCALSTRUCTURE) = scpackingeneNN(xyz,xyzf_thread,gn);
	//clock_t t2=clock();
	//std::cout <<t2-t1 <<std::endl;
	if (localhbweight_ > 0.0)
		potenergies->at(ELOCALHB) = localbbhbenergy(xyz, &xyzf);
#pragma omp parallel for schedule(dynamic,1)
	for(int t=0;t<ETERMS;++t){
		for(int i=0;i<nthread;++i){
			potenergies->at(t) += pot_thread[t][i];
		}
	}
#pragma omp parallel for schedule(dynamic,1)
	for(int a=0;a<natoms_;++a){
		NSPgeometry::XYZ &fa=xyzf.at(a);
		for(int i=0;i<nthread;++i){
			fa =fa+xyzf_thread[i][a];
		}
	}

#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<pointinfo.size();i++) {
		std::vector<int> vi = pointinfo[i].first.first;
		XYZ fxpt = pointinfo[i].first.second;
		//double fc = pointinfo[i].second;
		XYZ cen(0.0,0.0,0.0);
		for(int j:vi) cen = cen + XYZ(crd[j*3],crd[j*3+1],crd[j*3+2]);
		cen = cen / (double)(vi.size());
		XYZ fab = fxpt-cen;
		fab.z_ =0.0;
		double dfab = sqrt(fab.squarednorm());
		double fc = getf(pointinfo[i].second,dfab);
		fab = fab / dfab * fc;
#pragma omp critical
		{
			for(int j:vi) xyzf[j] = xyzf[j] + fab;
		}
	}

	potenergies->at(ETOT) = 0.0;
	for (int i = 1; i < ETERMS; ++i) {
		potenergies->at(ETOT) += potenergies->at(i);
	}
	std::vector<double> forces(3*natoms_,0.0);
#pragma omp parallel for schedule(dynamic,5)
//	for (auto it = xyzf.begin(); it != xyzf.end(); ++it) {
	for(int i=0;i<natoms_;++i){
		int i3=3*i;
		forces[i3]=xyzf[i].x_;
		forces[i3+1]=xyzf[i].y_;
		forces[i3+2]=xyzf[i].z_;
	}
	/*std::vector<double> fs1;
	double t1=0;
	for(int i=0;i<natoms_;++i){
		double f=sqrt(xyzf[i].squarednorm());
		fs1.push_back(f);
		t1 += f;
	}
	t1 /= (double)fs1.size();
	std::cout <<t1 <<std::endl;
	std::sort(fs1.begin(),fs1.end());
	for(double d:fs1) std::cout <<d <<std::endl;
	exit(1);*/
	//for(double d:*potenergies) std::cout <<d <<std::endl;
	//std::cout <<std::endl;
	//for(auto &c:xyzf) {
	//	std::cout <<c.x_ <<'\t' <<c.y_ <<'\t' <<c.z_ <<std::endl;
	//}
	//exit(1);
	return forces;
}
std::vector<double> ForceField::forces_attract(const std::vector<double> &crd,
		const NeighborList &nbl, const GroupNeighbor &gn, std::vector<double> *potenergies,
//		const std::vector<bool> *forceoff) const {
		ActiveSelections &acts, const std::vector<std::pair<std::vector<int>,double>>&attrs) const{
	if (eneanalysismode_) {
		int nsites = 0;
		for (auto &bs : bsinchains_) {
			nsites += bs.size();
		}
		enecomp_.init(nsites);
	}
	std::vector<XYZ> xyz;
	assert(natoms_ == crd.size() / 3);
	for (int i = 0; i < crd.size(); i += 3) {
		xyz.push_back(XYZ(crd[i], crd[i + 1], crd[i + 2]));
	}
	std::vector<XYZ> xyzf(natoms_, XYZ(0.0, 0.0, 0.0));
	int nthread=omp_get_max_threads();
	std::vector<std::vector<XYZ>> xyzf_thread(nthread,xyzf);
	potenergies->assign(ETERMS, 0.0);
	std::vector<bool> forceoff=acts.atomfixed();
	std::vector<std::vector<double>> pot_thread(ETERMS, std::vector<double>(nthread,0.0));
	pot_thread.at(EBOND) = covalentenergy(bondterms_, xyz, xyzf_thread, &forceoff);
	pot_thread.at(EANG) = covalentenergy(angleterms_, xyz, xyzf_thread, &forceoff);
	pot_thread.at(EIMPDIH) = covalentenergy(impdihterms_, xyz, xyzf_thread,
			&forceoff);
	int nchains = bsinchains_.size();
	/*std::vector<std::vector<PhiPsiCodes>> phipsicodes;
	for (int i = 0; i < nchains; ++i) {
		phipsicodes.push_back(makephipsicodes(crd, bsinchains_[i]));
	}
	std::vector<std::vector<SSCode>> sscodes;
	for (int i = 0; i < nchains; ++i) {
		sscodes.push_back(estimatess(phipsicodes[i]));
	}*/
	acts.calccodes(crd);
	const std::vector<std::vector<PhiPsiCodes>> & phipsicodes=acts.phipsicodes();
	const std::vector<std::vector<SSCode>> & sscodes=acts.sscodes();
	int nps=phipsis_.size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<nps;++i){
//	for (auto it = phipsis_.begin(); it != phipsis_.end(); ++it) {
		auto it=phipsis_.begin()+i;
		if(!acts.phipsiactive(it->chainid,it->siteid)) {
			continue;
		}
		double w = phipsiwght_ * siteweights_[it->chainid][it->siteid];
		int myid=omp_get_thread_num();
		double e = it->energy(phipsicodes[it->chainid], w, &(xyzf_thread[myid]));
		pot_thread.at(EPHIPSI)[myid] += e;
		if (eneanalysismode_) {
#pragma omp critical(eneanalysis)
			{
				enecomp_.ephipsi[i] = e;
			}
		}
	}
	if (lstermon_) {
		//pot_thread.at(ELOCALSTRUCTURE) = scpackingeneNN(xyz,xyzf_thread,gn);
		enecomp_.siteoffset = 0;
		for (int i = 0; i < nchains; ++i) {
			pot_thread.at(ELOCALSTRUCTURE) = lsenergy(phipsicodes[i],
					siteweights_[i], xyzf_thread,acts.lsactive(i));
			enecomp_.siteoffset += bsinchains_.size();
		}
	}
	if (sitepairtermon_) {
		std::vector<double> etrap_thread;
		pot_thread.at(ESITEPAIRS) = sitepairenergy(crd, phipsicodes, sscodes,
				xyzf_thread, etrap_thread,acts);
		for(int i=0;i<nthread;++i){
			pot_thread.at(ESTRUCTREST)[i] += etrap_thread[i];
		}
	}
	if (scconfwght_ != 0.0) {
		pot_thread.at(ESCCONF) = sideconfene(crd, phipsicodes, xyzf_thread,acts);
	}
	if (shadowweight_ != 0.0) {
		potenergies->at(ESHADOW) = shadowenergy(xyz, &xyzf);
	}
	for (auto it = structrestraints_.begin(); it != structrestraints_.end();
			++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	for (auto it = disrestraints_.begin(); it != disrestraints_.end(); ++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	potenergies->at(ERGRESTRAINT) += rgrestraint_.energy(xyz, &xyzf,
			coreweights_);
	for (auto it = ssrestraints_.begin(); it != ssrestraints_.end(); ++it) {
		potenergies->at(ESSRESTRAINT) += it->energy(sscodes, &xyzf);
	}
	if (stericwght_ > 0.0)
		pot_thread.at(ESTERIC) = stericenergy(nbl, xyz, xyzf_thread);
	//clock_t t1=clock();
	if (scpackingweight_>0.0)
		pot_thread.at(ESCPACKING) = scpackingene(xyz,xyzf_thread,acts);
	//pot_thread.at(ELOCALSTRUCTURE) = scpackingeneNN(xyz,xyzf_thread,gn);
	//clock_t t2=clock();
	//std::cout <<t2-t1 <<std::endl;
	if (localhbweight_ > 0.0)
		potenergies->at(ELOCALHB) = localbbhbenergy(xyz, &xyzf);
#pragma omp parallel for schedule(dynamic,1)
	for(int t=0;t<ETERMS;++t){
		for(int i=0;i<nthread;++i){
			potenergies->at(t) += pot_thread[t][i];
		}
	}
#pragma omp parallel for schedule(dynamic,1)
	for(int a=0;a<natoms_;++a){
		NSPgeometry::XYZ &fa=xyzf.at(a);
		for(int i=0;i<nthread;++i){
			fa =fa+xyzf_thread[i][a];
		}
	}
	//std::cout <<"size:\t" <<attrs.size() <<std::endl;
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<attrs.size();i++) {
		int a=attrs[i].first[0];
		int b=attrs[i].first[1];
		XYZ crd_a=XYZ(crd[a*3],crd[a*3+1],crd[a*3+2]);
		XYZ crd_b=XYZ(crd[b*3],crd[b*3+1],crd[b*3+2]);
		XYZ fab = crd_b-crd_a;
		fab = fab / sqrt(fab.squarednorm());
		fab = fab * attrs[i].second;
		XYZ fba = -fab;
#pragma omp critical
		{
			//std::cout <<attrs[i].second <<std::endl;
			//exit(1);
			xyzf[a] = xyzf[a] + fab;
			xyzf[b] = xyzf[b] + fba;
		}
	}
	potenergies->at(ETOT) = 0.0;
	for (int i = 1; i < ETERMS; ++i) {
		potenergies->at(ETOT) += potenergies->at(i);
	}
	std::vector<double> forces(3*natoms_,0.0);
#pragma omp parallel for schedule(dynamic,5)
//	for (auto it = xyzf.begin(); it != xyzf.end(); ++it) {
	for(int i=0;i<natoms_;++i){
		int i3=3*i;
		forces[i3]=xyzf[i].x_;
		forces[i3+1]=xyzf[i].y_;
		forces[i3+2]=xyzf[i].z_;
	}
	/*std::vector<double> fs1;
	double t1=0;
	for(int i=0;i<natoms_;++i){
		double f=sqrt(xyzf[i].squarednorm());
		fs1.push_back(f);
		t1 += f;
	}
	t1 /= (double)fs1.size();
	std::cout <<t1 <<std::endl;
	std::sort(fs1.begin(),fs1.end());
	for(double d:fs1) std::cout <<d <<std::endl;
	exit(1);*/
	//for(double d:*potenergies) std::cout <<d <<std::endl;
	//std::cout <<std::endl;
	//for(auto &c:xyzf) {
	//	std::cout <<c.x_ <<'\t' <<c.y_ <<'\t' <<c.z_ <<std::endl;
	//}
	//exit(1);
	return forces;
}



void ForceField::Site_Restrict_RMSD::getf(
		std::vector<XYZ>&xyzf, const std::vector<XYZ>&crds) const {
	if(fixed.empty()) return ;
	//std::cout <<"start restrict: " <<std::endl;
	//first, check if found site near fixed site
	//if yes, then record
	//std::vector<std::pair<int,std::vector<std::pair<int,double>>>> dissss;//fixed,site,rmsd
	std::vector<std::pair<std::vector<int>,double>> dissss; //fixed,site,rmsd
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<fixed.size();i++) {
		for(int j=0;j<nbb.size();j++) {
			double t=0.0;
			for(int k=0;k<4;k++) {
				t+=sqrt((crds[nbb[j][k]]-fixed[i][k]).squarednorm());
				if(t>dis_far) break;
			}
			if(t<dis_far) {
#pragma omp critical
				{
					dissss.push_back({{i,j},t});
				}
			}
		}
	}
	std::sort(dissss.begin(),dissss.end(),[](std::pair<std::vector<int>,double>p1,
			std::pair<std::vector<int>,double>p2)->bool{return p1.second<p2.second;});
	std::map<int,int> site_chosed; //first in nbb, second in fixed
	std::set<int> fixed_found;
	for(int i=0;i<dissss.size();i++) {
		if(site_chosed.find(dissss[i].first[1])!=site_chosed.end()) continue;
		if(fixed_found.find(dissss[i].first[0])!=fixed_found.end()) continue;
		site_chosed.insert({dissss[i].first[1],dissss[i].first[0]});
		fixed_found.insert(dissss[i].first[0]);
	}
	std::set<int> fixed_not_found;
	for(int i=0;i<fixed.size();i++) {
		if(fixed_found.find(i)==fixed_found.end()) fixed_not_found.insert(i);
	}
	//std::cout <<"Check Point 1" <<std::endl;
	//add force for all Calpha
	if(!fixed_not_found.empty()) {
#pragma omp parallel for schedule(dynamic,1)
		for(int i=0;i<nbb.size();i++) {
			if(site_chosed.find(i)!=site_chosed.end()) continue;
			int kc=nbb[i][1];
			for(int j:fixed_not_found) {
				xyzf[kc] = xyzf[kc] + attractforce(crds[kc],fixed[j][1]);
			}
		}
	}
	std::vector<std::vector<int>> site_chosed_v;
	for(auto i:site_chosed) site_chosed_v.push_back({i.first,i.second});
	//std::cout <<"Check Point 2" <<std::endl;
	//add force for chosed site
	if(!site_chosed_v.empty()) {
		std::cout <<site_chosed_v.size() <<std::endl;
#pragma omp parallel for schedule(dynamic,1)
		for(int i=0;i<site_chosed_v.size();i++) {
			int kf=site_chosed_v[i][1];
			for(int j=0;j<4;j++) {
				int kc=nbb[site_chosed_v[i][0]][j];
				xyzf[kc] = xyzf[kc] + attractforce(crds[kc],fixed[kf][j]);
			}
//#pragma omp critical
//			{
//				std::cout <<'\t' <<kf <<' ' <<site_chosed_v[i][0] <<std::endl;
//			}
		}
	}
	//std::cout <<"END restrict: " <<std::endl;
}
std::vector<double> ForceField::forces_fixedpoint_attraction_repulsion(
		const std::vector<double> &crd,const NeighborList &nbl, const GroupNeighbor &gn,
		std::vector<double> *potenergies,ActiveSelections &acts,
		const std::vector<Fixedpoint_Attraction_Repulsion> &fars,
		const Site_Restrict_RMSD &sirerm) const{
	if (eneanalysismode_) {
		int nsites = 0;
		for (auto &bs : bsinchains_) {
			nsites += bs.size();
		}
		enecomp_.init(nsites);
	}
	std::vector<XYZ> xyz;
	assert(natoms_ == crd.size() / 3);
	for (int i = 0; i < crd.size(); i += 3) {
		xyz.push_back(XYZ(crd[i], crd[i + 1], crd[i + 2]));
	}
	std::vector<XYZ> xyzf(natoms_, XYZ(0.0, 0.0, 0.0));
	int nthread=omp_get_max_threads();
	std::vector<std::vector<XYZ>> xyzf_thread(nthread,xyzf);
	potenergies->assign(ETERMS, 0.0);
	std::vector<bool> forceoff=acts.atomfixed();
	std::vector<std::vector<double>> pot_thread(ETERMS, std::vector<double>(nthread,0.0));
	pot_thread.at(EBOND) = covalentenergy(bondterms_, xyz, xyzf_thread, &forceoff);
	pot_thread.at(EANG) = covalentenergy(angleterms_, xyz, xyzf_thread, &forceoff);
	pot_thread.at(EIMPDIH) = covalentenergy(impdihterms_, xyz, xyzf_thread,
			&forceoff);
	int nchains = bsinchains_.size();
	/*std::vector<std::vector<PhiPsiCodes>> phipsicodes;
	for (int i = 0; i < nchains; ++i) {
		phipsicodes.push_back(makephipsicodes(crd, bsinchains_[i]));
	}
	std::vector<std::vector<SSCode>> sscodes;
	for (int i = 0; i < nchains; ++i) {
		sscodes.push_back(estimatess(phipsicodes[i]));
	}*/
	acts.calccodes(crd);
	const std::vector<std::vector<PhiPsiCodes>> & phipsicodes=acts.phipsicodes();
	const std::vector<std::vector<SSCode>> & sscodes=acts.sscodes();
	int nps=phipsis_.size();
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<nps;++i){
//	for (auto it = phipsis_.begin(); it != phipsis_.end(); ++it) {
		auto it=phipsis_.begin()+i;
		if(!acts.phipsiactive(it->chainid,it->siteid)) {
			continue;
		}
		double w = phipsiwght_ * siteweights_[it->chainid][it->siteid];
		int myid=omp_get_thread_num();
		double e = it->energy(phipsicodes[it->chainid], w, &(xyzf_thread[myid]));
		pot_thread.at(EPHIPSI)[myid] += e;
		if (eneanalysismode_) {
#pragma omp critical(eneanalysis)
			{
				enecomp_.ephipsi[i] = e;
			}
		}
	}
	if (lstermon_) {
		//pot_thread.at(ELOCALSTRUCTURE) = scpackingeneNN(xyz,xyzf_thread,gn);
		enecomp_.siteoffset = 0;
		for (int i = 0; i < nchains; ++i) {
			pot_thread.at(ELOCALSTRUCTURE) = lsenergy(phipsicodes[i],
					siteweights_[i], xyzf_thread,acts.lsactive(i));
			enecomp_.siteoffset += bsinchains_.size();
		}
	}
	if (sitepairtermon_) {
		std::vector<double> etrap_thread;
		pot_thread.at(ESITEPAIRS) = sitepairenergy(crd, phipsicodes, sscodes,
				xyzf_thread, etrap_thread,acts);
		for(int i=0;i<nthread;++i){
			pot_thread.at(ESTRUCTREST)[i] += etrap_thread[i];
		}
	}
	if (scconfwght_ != 0.0) {
		pot_thread.at(ESCCONF) = sideconfene(crd, phipsicodes, xyzf_thread,acts);
	}
	if (shadowweight_ != 0.0) {
		potenergies->at(ESHADOW) = shadowenergy(xyz, &xyzf);
	}
	for (auto it = structrestraints_.begin(); it != structrestraints_.end();
			++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	for (auto it = disrestraints_.begin(); it != disrestraints_.end(); ++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	potenergies->at(ERGRESTRAINT) += rgrestraint_.energy(xyz, &xyzf,
			coreweights_);
	for (auto it = ssrestraints_.begin(); it != ssrestraints_.end(); ++it) {
		potenergies->at(ESSRESTRAINT) += it->energy(sscodes, &xyzf);
	}
	if (stericwght_ > 0.0)
		pot_thread.at(ESTERIC) = stericenergy(nbl, xyz, xyzf_thread);
	//clock_t t1=clock();
	if (scpackingweight_>0.0)
		pot_thread.at(ESCPACKING) = scpackingene(xyz,xyzf_thread,acts);
	//pot_thread.at(ELOCALSTRUCTURE) = scpackingeneNN(xyz,xyzf_thread,gn);
	//clock_t t2=clock();
	//std::cout <<t2-t1 <<std::endl;
	if (localhbweight_ > 0.0)
		potenergies->at(ELOCALHB) = localbbhbenergy(xyz, &xyzf);
#pragma omp parallel for schedule(dynamic,1)
	for(int t=0;t<ETERMS;++t){
		for(int i=0;i<nthread;++i){
			potenergies->at(t) += pot_thread[t][i];
		}
	}
#pragma omp parallel for schedule(dynamic,1)
	for(int a=0;a<natoms_;++a){
		NSPgeometry::XYZ &fa=xyzf.at(a);
		for(int i=0;i<nthread;++i){
			fa =fa+xyzf_thread[i][a];
		}
	}

	for(int i=0;i<fars.size();i++) {
#pragma omp parallel for schedule(dynamic,1)
		for(int k=0;k<fars[i].atoms_seq.size();k++) {
			int j=fars[i].atoms_seq[k];
			XYZ c(crd[j*3],crd[j*3+1],crd[j*3+2]);
			for(XYZ fp:fars[i].fixed_points) {
				fp = fp-c;
				//if(fars[i].dis.size()==2) fp =-fp;
				double dfp=sqrt(fp.squarednorm());
				double f;
				if(!fars[i].getf(dfp,f)) fp =-fp;
				fp =fp/dfp*f;
#pragma omp critical
				{
					xyzf[j] = xyzf[j] + fp;
				}
			}
		}
	}
	sirerm.getf(xyzf, xyz);
	potenergies->at(ETOT) = 0.0;
	for (int i = 1; i < ETERMS; ++i) {
		potenergies->at(ETOT) += potenergies->at(i);
	}
	std::vector<double> forces(3*natoms_,0.0);
#pragma omp parallel for schedule(dynamic,5)
//	for (auto it = xyzf.begin(); it != xyzf.end(); ++it) {
	for(int i=0;i<natoms_;++i){
		int i3=3*i;
		forces[i3]=xyzf[i].x_;
		forces[i3+1]=xyzf[i].y_;
		forces[i3+2]=xyzf[i].z_;
	}
	/*std::vector<double> fs1;
	double t1=0;
	for(int i=0;i<natoms_;++i){
		double f=sqrt(xyzf[i].squarednorm());
		fs1.push_back(f);
		t1 += f;
	}
	t1 /= (double)fs1.size();
	std::cout <<t1 <<std::endl;
	std::sort(fs1.begin(),fs1.end());
	for(double d:fs1) std::cout <<d <<std::endl;
	exit(1);*/
	//for(double d:*potenergies) std::cout <<d <<std::endl;
	//std::cout <<std::endl;
	//for(auto &c:xyzf) {
	//	std::cout <<c.x_ <<'\t' <<c.y_ <<'\t' <<c.z_ <<std::endl;
	//}
	//exit(1);
	return forces;
}

#endif
