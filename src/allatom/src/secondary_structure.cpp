/*
 * secondary_structure.cpp
 *
 *  Created on: Nov 20, 2018
 *      Author: xuyang
 */
#include "allatom/secondary_structure.h"
#include <algorithm>
using namespace NSPallatom;

int Secondary_Structure::hydrogenbond(const Residue &r1, const Residue &r2) {
	//0-no hb, 1-r1->r2, 2-r2->r1, 3-both
	static double nrad=1.7;
	static double orad=1.6;
	static std::map<std::string,std::map<std::string,double>> rads;
	if(rads.empty()) {
		std::map<char,std::map<std::string,double>> r1=AA_Atom::radius();
		std::map<char,std::string> a13=AA_Atom::aa13();
		for(auto &r:r1) rads.insert({a13.at(r.first),r.second});
	}
	//std::cout <<rads.size() <<std::endl;
	//for(auto &r:rads) {
	//	std::cout <<r.first <<' ' <<r.second.size() <<std::endl;
	//	for(auto &p:r.second) std::cout <<'\t' <<p.first <<'\t' <<p.second <<std::endl;
	//}
	bool hb1{false};
	bool hb2{false};
	XYZ n1=r1.rds().at("N").crd;
	XYZ o1=r1.rds().at("O").crd;
	XYZ n2=r2.rds().at("N").crd;
	XYZ o2=r2.rds().at("O").crd;
	double rn1=nrad;
	double ro1=orad;
	double rn2=nrad;
	double ro2=orad;//std::cout <<r1.resname() <<std::endl;
	if(rads.find(r1.resname())!=rads.end()) {
		rn1=rads.at(r1.resname()).at("N");//std::cout <<'\t' <<rn1 <<std::endl;
		ro1=rads.at(r1.resname()).at("O");//std::cout <<'\t' <<ro1 <<std::endl;
	}
	if(rads.find(r2.resname())!=rads.end()) {
		rn2=rads.at(r2.resname()).at("N");
		ro2=rads.at(r2.resname()).at("O");
	}
	if((n1-o2).squarednorm()<(rn1+ro2)*(rn1+ro2)) hb1=true;
	if((n2-o1).squarednorm()<(rn2+ro1)*(rn2+ro1)) hb2=true;
	if(hb1&&hb2) return 3;
	if(!hb1&&hb2) return 2;
	if(hb1&&!hb2) return 1;
	return 0;
}

void Secondary_Structure::ss_divide() {
	//first find helix, then find strand
	//minium continuous number: helix-5, strand-3
	const std::vector<std::vector<Residue>> chains=pep_.chs_new();
	std::vector<int> ch_length;
	std::vector<std::pair<int,int>> lbs;
	for(int i=0;i<chains.size();i++) {
		for(int j=0;j<chains[i].size();j++) lbs.push_back({i,j});
		ch_length.push_back(chains[i].size());
	}
	std::set<std::pair<std::pair<int,int>,std::pair<int,int>>> hbps;//donor->accepter
	for(int i=0;i<lbs.size();i++) {//std::cout <<i <<' ' <<lbs.size() <<std::endl;
		int c1=lbs[i].first;
		int r1=lbs[i].second;
		for(int j=i+2;j<lbs.size();j++) {
			int c2=lbs[j].first;
			int r2=lbs[j].second;
			int l=hydrogenbond(chains[c1][r1],chains[c2][r2]);
			if(l==1) hbps.insert({{c1,r1},{c2,r2}});
			if(l==2) hbps.insert({{c2,r2},{c1,r1}});
			if(l==3) {
				hbps.insert({{c1,r1},{c2,r2}});
				hbps.insert({{c2,r2},{c1,r1}});
			}
		}
	}
	std::set<std::pair<int,int>> hlab=helix_recognition(ch_length,hbps);
	std::map<std::pair<int,int>,std::set<std::pair<int,int>>> n2c;
	std::map<std::pair<int,int>,std::set<std::pair<int,int>>> c2n;
	for(const auto &p:hbps) {
		if(hlab.find(p.first)!=hlab.end() && hlab.find(p.second)!=hlab.end()) continue;
		if(n2c.find(p.first)==n2c.end()) n2c.insert({p.first,std::set<std::pair<int,int>>()});
		if(c2n.find(p.second)==c2n.end()) c2n.insert({p.second,std::set<std::pair<int,int>>()});
		n2c.at(p.first).insert(p.second);
		c2n.at(p.second).insert(p.first);
	}
	strand_recognition(ch_length,n2c,c2n);
}

void Secondary_Structure::axis(int nh, int nb) {
	const auto & chs=pep_.chs_new();
	std::vector<std::vector<XYZ>> cas(chs.size());
	for(int i=0;i<cas.size();i++) {
		for(int j=0;j<chs[i].size();j++) {
			cas[i].push_back(chs[i][j].rds().at("CA").crd);
		}
	}
	for(auto & p:alpha_helixes) {
		//std::cout <<alpha_helixes.size() <<'\t' <<alpha_helixes[0].first <<'\t'
		//<<alpha_helixes[0].second.first <<'\t' <<alpha_helixes[0].second.second <<std::endl;
		//std::cout <<cas.size() <<'\t' <<cas[0].size() <<std::endl;
		helix_axis.push_back(std::vector<XYZ>());
		int c=p.first;
		for(int i=p.second.first;i<p.second.second-nh+1;i++) {
			double x=0.0, y=0.0, z=0.0;
			for(int j=0;j<nh;j++) {
				x+=cas[c][i+j].x_;
				y+=cas[c][i+j].y_;
				z+=cas[c][i+j].z_;
			}
			//std::cout <<cas[c][i+6].x_ <<std::endl;
			//double x=cas[c][i].x_+cas[c][i+1].x_+cas[c][i+2].x_+cas[c][i+3].x_;//+cas[c][i+4].x_+cas[c][i+5].x_+cas[c][i+6].x_+cas[c][i+7].x_;
			//double y=cas[c][i].y_+cas[c][i+1].y_+cas[c][i+2].y_+cas[c][i+3].y_;//+cas[c][i+4].y_+cas[c][i+5].y_+cas[c][i+6].y_+cas[c][i+7].y_;
			//double z=cas[c][i].z_+cas[c][i+1].z_+cas[c][i+2].z_+cas[c][i+3].z_;//+cas[c][i+4].z_+cas[c][i+5].z_+cas[c][i+6].z_+cas[c][i+7].z_;
			helix_axis.back().push_back(XYZ(x,y,z)/(double)nh);
		}
	}
	for(auto & p:beta_sheets) {
		strand_axis.push_back(std::vector<XYZ>());
		int c=p.first;
		for(int i=p.second.first;i<p.second.second-nb+1;i++) {
			double x=0.0, y=0.0, z=0.0;
			for(int j=0;j<nb;j++) {
				x+=cas[c][i+j].x_;
				y+=cas[c][i+j].y_;
				z+=cas[c][i+j].z_;
			}
			//double x=cas[c][i].x_+cas[c][i+1].x_+cas[c][i+2].x_;
			//double y=cas[c][i].y_+cas[c][i+1].y_+cas[c][i+2].y_;
			//double z=cas[c][i].z_+cas[c][i+1].z_+cas[c][i+2].z_;
			strand_axis.back().push_back(XYZ(x,y,z)/(double)nb);
		}
	}
}

std::set<std::pair<int,int>> Secondary_Structure::helix_recognition(
		const std::vector<int> &ch_length,
		const std::set<std::pair<std::pair<int,int>,std::pair<int,int>>> &hbps) {
	std::vector<polypeptide> h1;
	for(int i=0;i<ch_length.size();i++) {
		int j0;
		for(int j=0;j+4<ch_length[i];j++) {
			if(hbps.find({{i,j+4},{i,j}})==hbps.end()) continue;
			j0=j;
			for(;j+4<ch_length[i];j++) {
				if(hbps.find({{i,j+4},{i,j}})==hbps.end()) break;
				//if(j+4==chains[i].size()) break;
			}
			//if(j-j0<2) continue;
			h1.push_back({i,{j0,j+4}});
			//for(int j1=j0;j1<j+5;j1++) hlab.insert({i,j1});
		}
	}
	std::set<std::pair<int,int>> hlab;
	for(int i=0;i<h1.size();i++) {
		polypeptide p1=h1[i];
		for(i++;i<h1.size();i++) {
			if(p1.first!=h1[i].first) break;
			if(h1[i].second.first+2>p1.second.second) break;
			p1.second.second=h1[i].second.second;
		}
		i--;
		alpha_helixes.push_back(p1);
		for(int j=p1.second.first;j<p1.second.second;j++) hlab.insert({p1.first,j});
	}
	return hlab;
}

void Secondary_Structure::strand_recognition(const std::vector<int> &ch_length,
		const std::map<std::pair<int,int>,std::set<std::pair<int,int>>> &n2c,
		const std::map<std::pair<int,int>,std::set<std::pair<int,int>>> &c2n) {
	std::vector<polypeptide> strands;
	std::vector<polypeptide> nter;
	std::vector<polypeptide> cter;
	for(int i=0;i<ch_length.size();i++) {
		for(int j=0;j<ch_length[i];j++) {
			bool n1=n2c.find({i,j})!=n2c.end();
			bool c1=c2n.find({i,j})!=c2n.end();
			bool n2=n2c.find({i,j+2})!=n2c.end();
			bool c2=c2n.find({i,j+2})!=c2n.end();
			if(!c1 || !n2) continue;
			if(!n1 && !c2) continue;
			polypeptide ppep1{i,{j,j+2}};
			std::set<std::pair<int,int>> nc1;
			std::set<std::pair<int,int>> cn1=c2n.at({i,j});
			std::set<std::pair<int,int>> nc2=n2c.at({i,j+2});
			std::set<std::pair<int,int>> cn2;
			if(!c2) {
				nc1=n2c.at({i,j});
			} else if(!n1) {
				cn2=c2n.at({i,j+2});
			} else {
				nc1=n2c.at({i,j});
				cn2=c2n.at({i,j+2});
			}
			for(auto &p2:cn1) {
				for(auto &p3:nc2) {
					if(p2.first!=p3.first) continue;
					int dp=p3.second-p2.second;
					if(dp==-2) {//trans-parrallel
						polypeptide ppep2={p2.first,{p3.second,p2.second}};
						//std::cout <<ppep1.first <<'\t' <<ppep1.second.first <<'\t' <<ppep1.second.second <<'\t'
						//		<<ppep2.first <<'\t' <<ppep2.second.first <<'\t' <<ppep2.second.second <<std::endl;
						if(nc1.find(p2)==nc1.end()) {
							if(cn2.find(p3)==cn2.end()) continue;
							nter.push_back(ppep1);
							cter.push_back(ppep2);
						} else if(cn2.find(p3)==cn2.end()) {
							cter.push_back(ppep1);
							nter.push_back(ppep2);
						} else {
							strands.push_back(ppep1);
							strands.push_back(ppep2);
							//std::cout <<ppep1.first <<'\t' <<ppep1.second.first <<'\t' <<ppep1.second.second <<'\t'
							//		<<ppep2.first <<'\t' <<ppep2.second.first <<'\t' <<ppep2.second.second <<std::endl;
						}
					} else if(dp==0) {//parrallel
						polypeptide ppep2={p2.first,{p2.second,p2.second+2}};
						polypeptide ppep_1={p2.first,{p2.second-2,p2.second}};
						std::cout <<ppep1.first <<'\t' <<ppep1.second.first <<'\t' <<ppep1.second.second <<'\t'
								<<ppep_1.first <<'\t' <<ppep_1.second.first <<'\t' <<ppep_1.second.second <<'\t'
								<<ppep2.first <<'\t' <<ppep2.second.first <<'\t' <<ppep2.second.second <<std::endl;
						bool hb2{false};
						if(n2c.find({i,j+4})!=n2c.end()) {
							if(n2c.at({i,j+4}).find({p2.first,p2.second+2})==n2c.at({i,j+4}).end()) {
								hb2=true;
							}
						}
						bool hb1{false};
						if(c2n.find({i,j-2})!=c2n.end()) {
							if(c2n.at({i,j-2}).find({p2.first,p2.first-2})==c2n.at({i,j-2}).end()) {
								hb1=true;
							}
						}
						bool hb11=nc1.find({p2.first,p2.second-2})==nc1.end();
						bool hb22=cn2.find({p2.first,p2.second+2})==cn2.end();
						if(hb11 && hb22) strands.push_back(ppep1);
						else if(!hb11 && hb22) nter.push_back(ppep1);
						else if(hb11 && !hb22) cter.push_back(ppep1);
						if(hb11) {
							if(hb1) strands.push_back(ppep_1);
							else nter.push_back(ppep_1);
						}
						if(hb22) {
							if(hb2) strands.push_back(ppep2);
							else cter.push_back(ppep2);
						}
					}
				}
			}
		}
	}
	std::sort(strands.begin(),strands.end(),[](polypeptide p1,polypeptide p2)->bool{
		if(p1.first<p2.first) return true;
		else if((p1.first>p2.first)) return false;
		else if(p1.second.first<p2.second.first) return true;
		else if(p1.second.first>p2.second.first) return false;
		else if(p1.second.second<p2.second.second) return true;
		else return false;
	});
	for(int i=0;i<strands.size();i++) {
		polypeptide p=strands[i];//std::cout <<p.first <<'\t' <<p.second.first <<'\t' <<p.second.second <<std::endl;
		int j=i+1;
		for(;j<strands.size();j++) {
			if(strands[j].first!=strands[i].first) break;
			if(strands[j].second.first>strands[i].second.second+1) break;
			p.second.second=strands[j].second.second;
		}
		beta_sheets.push_back(p);
		i=j-1;
	}
	for(auto &p:beta_sheets) {
		for(auto &t:nter) {
			if(t.first!=p.first) continue;
			if(t.second.second!=p.second.first) continue;
			p.second.first=t.second.first;
			break;
		}
		for(auto &t:cter) {
			if(t.first!=p.first) continue;
			if(t.second.first!=p.second.second) continue;
			p.second.second=t.second.second;
			break;
		}
	}
	for(auto &p:beta_sheets) {
		p.second.second += 1;
	}
}







