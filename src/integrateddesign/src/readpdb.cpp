/*
 * readpdb.cpp
 *
 *  Created on: Jun 11, 2018
 *      Author: xuyang
 */
#include "integrateddesign/readpdb.h"
using namespace NSPproteinrep;

void PdbReader_xy::readpdb(std::string fn) {
	std::vector<PdbRecord> prsinpdb;
	std::ifstream ifs(fn);
	std::string lineinpdb;
	while(std::getline(ifs,lineinpdb)) {
		std::string label=lineinpdb.substr(0,6);
		if(label=="ENDMDL") {
			model_=true;
			break;
		}
		if(label!="ATOM  " && label!="HETATM") continue;
		PdbRecord pr(lineinpdb);
		if(prsinpdb.empty()) {
			prsinpdb.push_back(pr);
			continue;
		}
		if(pr.atomname==prsinpdb.back().atomname &&
				pr.residueid==prsinpdb.back().residueid &&
				pr.insertionid==prsinpdb.back().insertionid)
			continue;
		prsinpdb.push_back(pr);
	}
	ifs.close();

	std::vector<std::vector<PdbRecord>> ress;
	std::vector<PdbRecord> prstemp;
	int residtemp=prsinpdb[0].residueid;
	char insertionidtemp=prsinpdb[0].insertionid;
	for(auto &pr:prsinpdb) {
		if(pr.residueid!=residtemp || pr.insertionid!=insertionidtemp) {
			ress.push_back(prstemp);
			prstemp.clear();
			residtemp=pr.residueid;
			insertionidtemp=pr.insertionid;
		}
		prstemp.push_back(pr);
	}
	ress.push_back(prstemp);

	std::vector<std::vector<std::vector<PdbRecord>>> chains;
	PdbRecord cprev;
	std::vector<std::vector<PdbRecord>> chaintemp;
	for(auto &res:ress) {
		int nix=-1;
		int cix=-1;
		bool N,CA,C,O;
		N=CA=C=O=false;
		for(int i=0;i<res.size();i++) {
			if(res[i].atomname=="N") {nix=i;N=true;continue;}
			if(res[i].atomname=="CA") {CA=true;continue;}
			if(res[i].atomname=="C") {cix=i;C=true;continue;}
			if(res[i].atomname=="O") {O=true;continue;}
		}
		if(!N || !CA || !C || !O) {
			aa_excluded.push_back(res);
		} else {
			if(!chaintemp.empty()) {
				double dx=cprev.x-res[nix].x;
				double dy=cprev.y-res[nix].y;
				double dz=cprev.z-res[nix].z;
				if(dx*dx+dy*dy+dz*dz >4.0) {
					chains.push_back(chaintemp);
					chaintemp.clear();
				}
			}
			chaintemp.push_back(res);
			cprev=res[cix];
		}
	}
	if(!chaintemp.empty()) chains.push_back(chaintemp);

	char chainidtemp='A';
	for(auto & ch:chains) {
		records_.insert({chainidtemp,ResMapType()});
		for(int i=0;i<ch.size();i++) {
			residtemp=ch[i][0].residueid;
			insertionidtemp=ch[i][0].insertionid;
			records_.at(chainidtemp).insert({{residtemp,insertionidtemp},ch[i]});
		}
		chainidtemp++;
	}
}

std::vector<std::vector<BackBoneSite>> PdbReader_xy::getbackbonesite() {
	std::vector<std::vector<BackBoneSite>> bbsss;
	std::vector<BackBoneSite> bbss;
	BackBoneSite bbs;
	bbs.pdbid=pdbid;
	for(auto & rds:records_) {
		bbs.chainid=rds.first;
		int seqtemp=0;
		for(auto & rmt:rds.second) {
			bbs.resname=rmt.second[0].residuename;
			bbs.resid=rmt.second[0].residueid;
			bbs.resseq=seqtemp;
			seqtemp++;
			for(PdbRecord &pr:rmt.second) {
				if(pr.atomname=="N") {
					bbs.data_[BackBoneSite::NCRD]=pr.x;
					bbs.data_[BackBoneSite::NCRD+1]=pr.y;
					bbs.data_[BackBoneSite::NCRD+2]=pr.z;
				} else if(pr.atomname=="CA") {
					bbs.data_[BackBoneSite::CACRD]=pr.x;
					bbs.data_[BackBoneSite::CACRD+1]=pr.y;
					bbs.data_[BackBoneSite::CACRD+2]=pr.z;
				} else if(pr.atomname=="C") {
					bbs.data_[BackBoneSite::CCRD]=pr.x;
					bbs.data_[BackBoneSite::CCRD+1]=pr.y;
					bbs.data_[BackBoneSite::CCRD+2]=pr.z;
				} else if(pr.atomname=="O") {
					bbs.data_[BackBoneSite::OCRD]=pr.x;
					bbs.data_[BackBoneSite::OCRD+1]=pr.y;
					bbs.data_[BackBoneSite::OCRD+2]=pr.z;
				}
			}
			bbss.push_back(bbs);
		}
		bbsss.push_back(bbss);
		bbss.clear();
	}
	return bbsss;
}














