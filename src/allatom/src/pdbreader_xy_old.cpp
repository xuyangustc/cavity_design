/*
 * pdbreader_xy.cpp
 *
 *  Created on: Aug 13, 2018
 *      Author: xuyang
 */
#include "../pdbreader_xy_old.h"

#include "dstl/randomengine.h"
using namespace NSPallatom;
using namespace NSPproteinrep;
using namespace NSPgeometry;


std::map<std::string,char> AA20 {
	{"GLY",'G'},
	{"ALA",'A'},
	{"VAL",'V'},
	{"LEU",'L'},
	{"ILE",'I'},
	{"SER",'S'},
	{"THR",'T'},
	{"CYS",'C'},
	{"MET",'M'},
	{"ASP",'D'},
	{"GLU",'E'},
	{"ASN",'N'},
	{"GLN",'Q'},
	{"LYS",'K'},
	{"ARG",'R'},
	{"HIS",'H'},
	{"PRO",'P'},
	{"PHE",'F'},
	{"TYR",'Y'},
	{"TRP",'W'}
};

std::vector<std::vector<std::string>> AtomPair::crdmainchain{
	{"O","C","CA"},{"C","O","CA"},{"CA","C","N"},{"N","CA","C"}
};

std::map<std::string,std::vector<std::vector<std::string>>> AtomPair::crdsidechain{
	{"GLY",{  }},
	{"ALA",{ {"C","CA","CB"} }},
	{"SER",{ {"C","CA","CB","OG"} }},
	{"CYS",{ {"C","CA","CB","SG"} }},
	{"MET",{ {"C","CA","CB","CG","SD","CE"} }},
	{"LYS",{ {"C","CA","CB","CG","CD","CE","NZ"} }},

	{"VAL",{ {"C","CA","CB","CG1"},{"CA","CB","CG2"} }},
	{"THR",{ {"C","CA","CB","OG1"},{"CA","CB","CG2"} }},
	{"LEU",{ {"C","CA","CB","CG","CD1"},{"CB","CG","CD2"} }},
	{"ASP",{ {"C","CA","CB","CG","OD1"},{"CB","CG","OD2"} }},
	{"ASN",{ {"C","CA","CB","CG","OD1"},{"CB","CG","ND2"} }},
	{"GLU",{ {"C","CA","CB","CG","CD","OE1"},{"CG","CD","OE2"} }},
	{"GLN",{ {"C","CA","CB","CG","CD","OE1"},{"CG","CD","NE2"} }},
	{"ARG",{ {"C","CA","CB","CG","CD","NE","CZ","NH1"},{"NE","CZ","NH2"} }},
	{"ILE",{ {"C","CA","CB","CG1","CD1"},{"CA","CB","CG2"} }},

	{"PRO",{ {"C","CA","CB","CG","CD"} }},
	{"HIS",{ {"C","CA","CB","CG","ND1","CE1","NE2","CD2"} }},
	{"PHE",{ {"C","CA","CB","CG","CD1","CE1","CZ"},{"CB","CG","CD2","CE2"} }},
	{"TYR",{ {"C","CA","CB","CG","CD1","CE1","CZ","OH"},{"CB","CG","CD2","CE2"} }},
	{"TRP",{ {"C","CA","CB","CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"} }},
};

std::set<std::string> AtomPair::getresidues() {
	std::set<std::string> ress;
	for(auto sc:crdsidechain) {
		ress.insert(sc.first);
	}
	return ress;
}
std::set<std::string> AtomPair::residues=getresidues();

std::set<std::string> AtomPair::getmainchains() {
	std::set<std::string> mc;
	for(auto ss:crdmainchain) {
		for(auto s:ss) mc.insert(s);
	}
	return mc;
}
std::set<std::string> AtomPair::mainchains=getmainchains();

std::map<std::string,std::set<std::string>> AtomPair::getsidechains() {
	std::map<std::string,std::set<std::string>> scs;
	for(auto sc:crdsidechain) {
		std::set<std::string> as;
		for(auto ss:sc.second) {
			for(auto s:ss) as.insert(s);
		}
		scs.insert({sc.first,as});
	}
	return scs;
}
std::map<std::string,std::set<std::string>> AtomPair::sidechains=getsidechains();

std::set<AtomPair::AtomName> AtomPair::getatomnames() {
	std::set<AtomPair::AtomName> ans;
	for(auto cs:AtomPair::crdsidechain) {
		for(auto cm:AtomPair::crdmainchain) {
			for(std::string s:cm) ans.insert({cs.first,s});
		}
		for(auto sd:cs.second) {
			for(std::string s:sd) ans.insert({cs.first,s});
		}
	}
	//std::cout <<"get neighbors: " <<ngs.size() <<std::endl;
	return ans;
}
std::set<AtomPair::AtomName> AtomPair::atomnames=getatomnames();

std::map<AtomPair::AtomName,std::pair<AtomPair::AtomName,AtomPair::AtomName>>
		AtomPair::getneighbors() {
	std::map<AtomPair::AtomName,std::pair<AtomPair::AtomName,AtomPair::AtomName>> ngs;
	AtomPair::AtomName cen;
	AtomPair::AtomName fi;
	AtomPair::AtomName se;
	for(auto cs:AtomPair::crdsidechain) {
		for(auto cm:AtomPair::crdmainchain) {
			cen={cs.first,cm[0]};
			fi={cs.first,cm[1]};
			se={cs.first,cm[2]};
			ngs.insert({cen,{fi,se}});
		}
		for(auto sd:cs.second) {
			for(int i=2;i<sd.size();i++) {
				cen={cs.first,sd[i]};
				fi={cs.first,sd[i-1]};
				se={cs.first,sd[i-2]};
				ngs.insert({cen,{fi,se}});
			}
		}
	}
	//std::cout <<"get neighbors: " <<ngs.size() <<std::endl;
	return ngs;
}
std::map<AtomPair::AtomName,std::pair<AtomPair::AtomName,
		AtomPair::AtomName>> AtomPair::neighbors=getneighbors();

std::map<AtomPair::AtomName,double> AtomPair::getradius() {
	std::map<AtomPair::AtomName,double> rads;
	std::string fn=NSPdataio::datafilename("resLib.dat");
	std::ifstream ifs(fn);
	std::string line;
	while(std::getline(ifs,line)) {
		std::istringstream iss(line);
		std::string s1,s2;
		int nline;
		iss >>s1 >>s2 >>nline;
		for(int i=0;i<nline;i++) {
			std::getline(ifs,line);
			std::istringstream ss(line);
			std::string s3,s4;
			double d;
			ss >>s3 >>s4 >>d;
			rads.insert({{s2,s4},d});
		}
	}
	ifs.close();
	return rads;
}
std::map<AtomPair::AtomName,double> AtomPair::radius=getradius();

std::set<AtomPair::AtomName> AtomPair::getdonors() {
	std::set<AtomPair::AtomName> ats;
	auto scs=crdsidechain;
	for(auto &sc:scs) {
		if(sc.first!="PRO") ats.insert({sc.first,"N"});
		for(auto &s:sc.second) for(auto &c:s) {if(c[0]=='N'||c[0]=='O') ats.insert({sc.first,c});}
	}
	return ats;
}
std::set<AtomPair::AtomName> AtomPair::donors=getdonors();

std::set<AtomPair::AtomName> AtomPair::getaccepters() {
	std::set<AtomPair::AtomName> ats;
	auto scs=crdsidechain;
	for(auto &sc:scs) {
		ats.insert({sc.first,"O"});
		for(auto &s:sc.second) for(auto &c:s) {if(c[0]=='N'||c[0]=='O') ats.insert({sc.first,c});}
	}
	return ats;
}
std::set<AtomPair::AtomName> AtomPair::accepters=getaccepters();

std::map<std::pair<char,char>,int> BLOSUM62::getmatrix() {
	std::map<std::pair<char,char>,int> mt;
	std::map<char,std::vector<int>> dt;
	std::string fn=NSPdataio::datafilename("blosum.mat");
	std::ifstream ifs(fn);
	std::string line;
	std::stringstream ss;
	for(int i=0;i<20;i++) {
		std::getline(ifs,line);
		ss << line;
		char c;
		ss >> c;
		std::vector<int> vi(i+1);
		for(int j=0;j<vi.size();j++) ss >> vi[j];
		ss.clear();
		dt.insert({c,vi});
	}
	std::getline(ifs,line);
	ss << line;
	std::vector<char> ress(20);
	for(int i=0;i<20;i++) ss >> ress[i];
	for(auto &d:dt) {
		std::vector<int> & vi=d.second;
		for(int i=0;i<vi.size();i++) {
			mt[{d.first,ress[i]}] = vi[i];
		}
	}
	ifs.close();
	return mt;
}
std::map<std::pair<char,char>,int> BLOSUM62::matrix=getmatrix();

void Residue::init(const std::vector<NSPproteinrep::PdbRecord> &prs) {
	chainid_ = prs[0].chainid;
	resid_ = prs[0].residueid;
	resname_ = prs[0].residuename;
	insertionid_ = prs[0].insertionid;
	if(AA20.find(resname_)==AA20.end()) {
		for(const PdbRecord &p:prs) crds_.insert({p.atomname,NSPgeometry::XYZ(p.x,p.y,p.z)});
		return ;
	}
	auto &vsc=NSPsd::VSCType::getVSCType(resname_);
	std::vector<std::string> scatoms=vsc.atomnames;
	std::set<std::string> atoms;
	for(std::string & s:scatoms) atoms.insert(s);
	atoms.insert("N");
	atoms.insert("CA");
	atoms.insert("C");
	atoms.insert("O");
	for(const PdbRecord &p:prs) {
		if(atoms.find(p.atomname)==atoms.end()) continue;
		crds_.insert({p.atomname,NSPgeometry::XYZ(p.x,p.y,p.z)});
	}
}

NSPproteinrep::FullSite Residue::Tofullsite() {
	FullSite fs;
	fs.changecrd(crds_);
	fs.resname()=resname_;
	fs.chainid()=chainid_;
	fs.resid()=resid_;
	fs.resseq()=resseq_;
	fs.insertionid()=insertionid_;
	return fs;
}

std::map<std::pair<std::string,std::string>,std::vector<XYZ>> Residue::getneighbor() {
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> nghs;
	for(auto c:crds_) {
		AtomPair::AtomName an={resname_,c.first};
		if(AtomPair::neighbors.find(an)==AtomPair::neighbors.end()) {
			std::cout <<"Atoms beyond 167 when get neighbor!" <<std::endl;
			std::cout <<'\t' <<resname_ <<'\t' <<c.first <<std::endl;
			exit(1);
		}
		auto fise=AtomPair::neighbors.at({resname_,c.first});
		std::string finame=fise.first.second;
		if(crds_.find(finame)==crds_.end()) continue;
		std::string sename=fise.second.second;
		if(crds_.find(sename)==crds_.end()) continue;
		std::vector<XYZ> nghcrd;
		nghcrd.push_back(crds_.at(c.first));
		nghcrd.push_back(crds_.at(finame));
		nghcrd.push_back(crds_.at(sename));
		nghs.insert({{resname_,c.first},nghcrd});
	}
	return nghs;
}

std::vector<XYZ> Residue::covert20(std::string resn) {
	std::vector<XYZ> cs;
	if(AA20.find(resname_)!=AA20.end()) return cs;
	if(resname_=="MSE") {
		resname_="MET";
		if(crds_.find("SE")!=crds_.end()) {
			XYZ c=crds_.at("SE");
			crds_.erase("SE");
			crds_.insert({"SD",c});
		}
		std::map<std::string,XYZ> crds2;
		std::set<std::string> metsc=AtomPair::sidechains.at("MET");
		metsc.insert("N");
		metsc.insert("CA");
		metsc.insert("C");
		metsc.insert("O");
		for(auto c:crds_) {
			if(metsc.find(c.first)==metsc.end()) cs.push_back(c.second);
			else crds2.insert(c);
		}
		crds_ = crds2;
		return cs;
	}
	if(resn.empty()) return cs;
	resname_ = resn;
	std::map<std::string,XYZ> newc;
	for(auto &c:crds_) {
		if(c.first=="N" || c.first=="CA" || c.first=="C" || c.first=="O") newc.insert(c);
		else cs.push_back(c.second);
	}
	crds_=newc;
	return cs;
}

void Residue::printrecords(int &atominitnum, std::ostream &os) const {
	PdbRecord p;
	p.label = "ATOM";
	p.chainid = 'A' + chainseq_;
	p.residueid = resseq_;
	std::vector<std::string> rds;
	for(const auto &c:crds_) {
		PdbRecord p1=p;
		p1.namesymbol = c.first.substr(0,1);
		p1.namemodifier = c.first.substr(1);
		p1.atomname = c.first;
		p1.residuename = resname_;
		p1.atomid = atominitnum++;
		p1.x = c.second.x_;
		p1.y = c.second.y_;
		p1.z = c.second.z_;
		p1.elementname[1] = c.first[0];
		os <<p1.toString() <<std::endl;
	}
}

void PdbReader_xy::init(std::string fn, std::string out20) {
	std::vector<std::vector<PdbRecord>> prs;
	std::ifstream ifs(fn);
	std::string line;
	std::vector<PdbRecord> ps;
	while(std::getline(ifs,line)) {
		std::string subl = line.substr(0,6);
		if(subl=="ENDMDL") break;
		if(subl=="ATOM  " || subl=="HETATM") {
			if(ps.empty()) ps.push_back(PdbRecord(line));
			else {
				PdbRecord p(line);
				PdbRecord p2=ps.back();
				if(p2.chainid!=p.chainid || p2.residueid!=p.residueid || p2.insertionid!=p.insertionid) {
					prs.push_back(ps);
					ps.clear();
				}
				ps.push_back(p);
			}
		}
	}
	if(!ps.empty()) prs.push_back(ps);
	ifs.close();

	for(int i=0;i<prs.size();i++) {
		if(prs[i][0].residuename=="HOH") {
			waters_.push_back(prs[i][0]);
			continue;
		}
		bool n,ca,c,o;
		n=ca=c=o=false;
		for(PdbRecord &p:prs[i]) {
			if(p.atomname=="N") n=true;
			else if(p.atomname=="CA") ca=true;
			else if(p.atomname=="C") c=true;
			else if(p.atomname=="O") o=true;
		}
		if(!n||!ca||!c||!o) {
			ligands_.push_back(prs[i]);
			continue;
		}
		Residue res;
		res.init(prs[i]);
		if(protein_.empty()) {
			protein_.push_back(Chain());
		} else {
			Residue & resprev=protein_.back().back();
			NSPgeometry::XYZ crdc=resprev.crds().at("C");
			NSPgeometry::XYZ crdn=res.crds().at("N");
			double bondlength = NSPgeometry::distance(crdc,crdn);
			if(bondlength>2.0) protein_.push_back(Chain());
		}
		res.chainseq() = protein_.size()-1;
		res.resseq() = protein_.back().size();
		protein_.back().push_back(res);
		if(AA20.find(res.resname())==AA20.end()) {
			beyond20s_.insert({res.chainseq(),res.resseq()});
		}
	}
	for(Chain &ch:protein_) {
		for(Residue &r:ch) {
			std::string resn=r.resname();
			std::vector<XYZ> cs=r.covert20(out20);
			if(cs.empty()) continue;
			std::vector<PdbRecord> prs;
			for(XYZ &c:cs) prs.push_back(getarecord(resn,"C",c));
		}
	}
}

XYZ PdbReader_xy::center_pep() {
	XYZ center(0.0,0.0,0.0);
	int natom=0;
	center = NSPgeometry::XYZ(0.0,0.0,0.0);
	for(const Chain &ch:protein_) {
		for(const Residue &res:ch) {
			std::string resname = res.resname();
			std::map<std::string,XYZ> cs=res.crds();
			for(auto & c:cs) {
				center.x_ += c.second.x_;
				center.y_ += c.second.y_;
				center.z_ += c.second.z_;
				natom++;
			}
		}
	}
	center.x_ /= (double)(natom);
	center.y_ /= (double)(natom);
	center.z_ /= (double)(natom);
	return center;
}

double PdbReader_xy::maxrad_pep(XYZ &center) {
	double maxrad=0.0;
	for(const Chain &ch:protein_) {
		for(const Residue &res:ch) {
			std::map<std::string,XYZ> cs=res.crds();
			for(auto & c:cs) {
				double d2=(center-c.second).squarednorm();
				if(d2>maxrad) maxrad=d2;
			}
		}
	}
	return sqrt(maxrad);
}

std::map<std::pair<std::string,std::string>,std::vector<XYZ>> PdbReader_xy::neighbors() {
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> nghs;
	for(Chain &ch:protein_) {
		for(Residue &r:ch) {
			if(AA20.find(r.resname())==AA20.end()) continue;
			std::map<std::pair<std::string,std::string>,std::vector<XYZ>> ns=r.getneighbor();
			//std::cout <<"Residue Neighbor is : " <<ns.size() <<std::endl;
			for(auto &n:ns) {
				if(nghs.find(n.first)==nghs.end()) nghs.insert({n.first,std::vector<XYZ>()});
				for(XYZ &c:n.second) nghs.at(n.first).push_back(c);
			}
			//std::cout <<"Total Neighbor is : " <<nghs.size() <<std::endl;
		}
	}
	return nghs;
}

std::map<std::pair<std::string,std::string>,std::vector<XYZ>> PdbReader_xy::peps() {
	std::map<std::pair<std::string,std::string>,std::vector<XYZ>> cs;
	for(Chain &ch:protein_) {
		for(Residue &r:ch) {
			std::string resname=r.resname();
			if(AA20.find(resname)==AA20.end()) continue;
			std::map<std::string,XYZ> rescrds=r.crds();
			for(auto &a:rescrds) {
				std::pair<std::string,std::string> p{resname,a.first};
				if(cs.find(p)==cs.end()) cs.insert({p,std::vector<XYZ>()});
				cs.at(p).push_back(a.second);
			}
		}
	}
	return cs;
}

int PdbReader_xy::print(std::ostream &os) {
	std::vector<std::vector<NSPproteinrep::FullSite>> conf;
	for(int i=0;i<protein_.size();i++) {
		conf.push_back(std::vector<NSPproteinrep::FullSite>());
		auto &ch =conf.back();
		for(int j=0;j<protein_[i].size();j++) {
			ch.push_back(protein_[i][j].Tofullsite());
		}
	}
	writetopdb(conf,os);
	for(auto & lg:ligands_) {
		for(auto & l:lg) os <<l.toString() <<std::endl;
	}
	int i=0;
	for(auto & w:waters_) {
		os <<w.toString() <<std::endl;
		if(w.atomid>i) i=w.atomid;
	}
	return i;
}

NSPproteinrep::PdbRecord PdbReader_xy::getarecord(std::string resname, std::string atomname, const XYZ &c) {
	PdbRecord p;
	p.label = "ATOM";
	p.residuename = resname;
	p.namesymbol = atomname.substr(0,1);
	p.namemodifier = atomname.substr(1);
	p.atomid = 0;
	p.residueid = 0;
	p.chainid = 'A';
	p.elementname[1] = atomname[0];
	p.x = c.x_;
	p.y = c.y_;
	p.z = c.z_;
	return p;
}

void PdbReader_xy::addss(std::string pdbfile) {
	std::map<char,std::map<int,std::map<std::string,std::vector<std::pair<char,bool>>>>> sstype;
	//<chainid,<seq,<restye,<sstype,readed>>>>, vector exist because of insertion number
	std::string cmd = "stride "+pdbfile;
	FILE *pp = popen(cmd.c_str(),"r");
	if(!pp) {
		std::cout <<"making strfile wrong!" <<std::endl;
		exit(1);
	}
	char tmp[100];
	while(!feof(pp)) {
		fgets(tmp,100,pp);
		std::string line(tmp);
		if(line.substr(0,3)!="ASG") continue;
		char cid=line[9];
		int seq=std::stoi(line.substr(11,4));
		std::string res=line.substr(5,3);
		char ss=line[24];
		sstype[cid][seq][res].push_back({ss,false});
	}
	pclose(pp);

	for(Chain &ch:protein_) {
		for(Residue &res:ch) {
			auto & pcb = sstype[res.chid()][res.resid()][res.resname()];
			bool finded{false};
			for(auto & cb:pcb) {
				if(cb.second) continue;
				res.ssid() = cb.first;
				finded = true;
			}
			if(!finded) {
				std::cout <<pdbfile <<" Can Not Find Residue: " <<res.chid() <<' ' <<res.resid() <<' ' <<res.resname() <<std::endl;
				continue;
			}
		}
	}
}

std::vector<NSPgeometry::XYZ> Grid::getwater(const NSPgeometry::XYZ &cen, double rad, double wtrad) {
	std::vector<XYZ> waters;
	int ninline=(int)(rad/wtrad)+1;

	std::vector<XYZ> vol;
	for(int x=-1;x<ninline;x++) {
		for(int y=-1;y<ninline;y++) {
			for(int z=-1;z<ninline;z++) {
				XYZ c(cen.x_-rad+wtrad*2.0*(double)x,
						cen.y_-rad+wtrad*2.0*(double)y,
						cen.z_-rad+wtrad*2.0*(double)z);
				vol.push_back(c);
			}
		}
	}

	for(int i=0;i<vol.size();i++) {
		//double d2=(vol[i]-cen).squarednorm();
		//if(d2>rad*rad) continue;
		waters.push_back(vol[i]);
	}
	return waters;
}

std::vector<Square> Grid::getrange(const NSPgeometry::XYZ &cen, double rad) {
	rad += Square_Length;
	std::vector<XYZ> waters;//=getwater(cen,rad,Square_Length);
	int ninline=(int)(rad/Square_Length)+1;

	std::vector<XYZ> vol;
	for(int x=-1;x<ninline;x++) {
		for(int y=-1;y<ninline;y++) {
			for(int z=-1;z<ninline;z++) {
				XYZ c(cen.x_-rad+Square_Length*2.0*(double)x,
						cen.y_-rad+Square_Length*2.0*(double)y,
						cen.z_-rad+Square_Length*2.0*(double)z);
				vol.push_back(c);
			}
		}
	}

	for(int i=0;i<vol.size();i++) {
		double d2=(vol[i]-cen).squarednorm();
		if(d2>rad*rad) continue;
		waters.push_back(vol[i]);
	}

	std::vector<Square> rgs;
	for(XYZ &c:waters) rgs.push_back(Square(c.x_, c.y_, c.z_));
	return rgs;
}

void MixModel::init(PdbReader_xy &p) {
	center_=p.center_pep();
	double maxrad=p.maxrad_pep(center_);
	waters_=Grid::getwater(center_,maxrad+rad_add,rad_water);
	//std::vector<Square> rs=Grid::getrange(center_,maxrad+rad_add*2.0);//make sure rs contain all waters
	//for(Square &r:rs) ranges_.insert({r,Content()});
	//partion(true,0.0);
	nghs_=p.neighbors();//std::cout <<"Neighbor size is : " <<nghs_.size() <<std::endl;
	crds_=p.peps();//std::cout <<"Crd size is : " <<crds_.size() <<std::endl;
	pdbid_ = p.pdbid();
}

void MixModel::microenvironment(std::string outpath, bool calatom, bool calwater) {
	if(outpath.back()!='/') outpath+='/';
	std::ofstream ofs;
	for(auto & ns:nghs_) {
		std::string nghres=ns.first.first;
		std::string nghcen=ns.first.second;
		if(AtomPair::neighbors.find(ns.first)==AtomPair::neighbors.end()) {
			std::cout <<"Neighbor is not in the definition!" <<std::endl;
			std::cout <<'\t' <<ns.first.first <<'\t' <<ns.first.second <<std::endl;
			exit(1);
		}
		auto nghfise=AtomPair::neighbors.at(ns.first);
		std::string nghfi=nghfise.first.second;
		std::string nghse=nghfise.second.second;
		std::vector<XYZ> & nghcrds = ns.second;
		if(calatom) {
			for(auto &cs:crds_) {
				std::vector<std::pair<std::vector<XYZ>,std::vector<XYZ>>> envs;
				std::vector<XYZ> & envcrds = cs.second;
				for(int i=0;i<nghcrds.size();i+=3) {
					std::vector<XYZ> p1;
					p1.push_back(nghcrds[i]);
					p1.push_back(nghcrds[i+1]);
					p1.push_back(nghcrds[i+2]);
					std::vector<XYZ> p2;
					for(int j=0;j<envcrds.size();j++) {
						double x=envcrds[j].x_-nghcrds[i].x_;
						if(x<-rad_add || x>rad_add) continue;
						double y=envcrds[j].y_-nghcrds[i].y_;
						if(y<-rad_add || y>rad_add) continue;
						double z=envcrds[j].z_-nghcrds[i].z_;
						if(z<-rad_add || z>rad_add) continue;
						double d2=x*x+y*y+z*z;
						if(d2<0.001 || d2>rad_add*rad_add) continue;
						p2.push_back(envcrds[j]);
					}
					if(p2.empty()) continue;
					envs.push_back({p1,p2});
				}
				if(envs.empty()) continue;
				ofs.open(outpath+nghres+"_"+nghcen+"_"+cs.first.first+"_"+cs.first.second,std::ofstream::app);
				for(auto &ev:envs) {
					ofs <<"START\t" <<(ev.second.size()+3)*3 <<'\t' <<pdbid_ <<std::endl;
					ofs <<ev.first[0].x_ <<std::endl;
					ofs <<ev.first[0].y_ <<std::endl;
					ofs <<ev.first[0].z_ <<std::endl;
					ofs <<ev.first[1].x_ <<std::endl;
					ofs <<ev.first[1].y_ <<std::endl;
					ofs <<ev.first[1].z_ <<std::endl;
					ofs <<ev.first[2].x_ <<std::endl;
					ofs <<ev.first[2].y_ <<std::endl;
					ofs <<ev.first[2].z_ <<std::endl;
					/*PdbRecord p;
					p=PdbReader_xy::getarecord(nghres,nghcen,ev.first[0]);
					ofs <<p.toString() <<std::endl;
					p=PdbReader_xy::getarecord(nghres,nghfi,ev.first[1]);
					ofs <<p.toString() <<std::endl;
					p=PdbReader_xy::getarecord(nghres,nghse,ev.first[2]);
					ofs <<p.toString() <<std::endl;*/
					for(int i=0;i<ev.second.size();i++) {
						//p=PdbReader_xy::getarecord(cs.first.first,cs.first.second,ev.second[i]);
						//ofs <<p.toString() <<std::endl;
						ofs <<ev.second[i].x_ <<std::endl;
						ofs <<ev.second[i].y_ <<std::endl;
						ofs <<ev.second[i].z_ <<std::endl;
					}
				}
				ofs.close();
			}
		}
		if(calwater) {
			std::vector<std::pair<std::vector<XYZ>,std::vector<XYZ>>> envs;
			for(int i=0;i<nghcrds.size();i+=3) {
				std::vector<XYZ> p1;
				p1.push_back(nghcrds[i]);
				p1.push_back(nghcrds[i+1]);
				p1.push_back(nghcrds[i+2]);
				std::vector<XYZ> p2;
				for(int j=0;j<waters_.size();j++) {
					double x=waters_[j].x_-nghcrds[i].x_;
					if(x<-rad_add || x>rad_add) continue;
					double y=waters_[j].y_-nghcrds[i].y_;
					if(y<-rad_add || y>rad_add) continue;
					double z=waters_[j].z_-nghcrds[i].z_;
					if(z<-rad_add || z>rad_add) continue;
					double d2=x*x+y*y+z*z;
					double dmin=Water_Radius+AtomPair::radius.at(ns.first);
					if(d2<0.001 || d2>rad_add*rad_add) continue;
					p2.push_back(waters_[j]);
				}
				if(p2.empty()) continue;
				envs.push_back({p1,p2});
			}
			if(envs.empty()) continue;
			ofs.open(outpath+nghres+"_"+nghcen+"_HOH_O",std::ofstream::app);

			for(auto &ev:envs) {
				ofs <<"START\t" <<(ev.second.size()+3)*3 <<'\t' <<pdbid_ <<std::endl;
				ofs <<ev.first[0].x_ <<std::endl;
				ofs <<ev.first[0].y_ <<std::endl;
				ofs <<ev.first[0].z_ <<std::endl;
				ofs <<ev.first[1].x_ <<std::endl;
				ofs <<ev.first[1].y_ <<std::endl;
				ofs <<ev.first[1].z_ <<std::endl;
				ofs <<ev.first[2].x_ <<std::endl;
				ofs <<ev.first[2].y_ <<std::endl;
				ofs <<ev.first[2].z_ <<std::endl;
				/*PdbRecord p;
				p=PdbReader_xy::getarecord(nghres,nghcen,ev.first[0]);
				ofs <<p.toString() <<std::endl;
				p=PdbReader_xy::getarecord(nghres,nghfi,ev.first[1]);
				ofs <<p.toString() <<std::endl;
				p=PdbReader_xy::getarecord(nghres,nghse,ev.first[2]);
				ofs <<p.toString() <<std::endl;*/
				for(int i=0;i<ev.second.size();i++) {
					//p=PdbReader_xy::getarecord("HOH","O",ev.second[i]);
					//ofs <<p.toString() <<std::endl;
					ofs <<ev.second[i].x_ <<std::endl;
					ofs <<ev.second[i].y_ <<std::endl;
					ofs <<ev.second[i].z_ <<std::endl;
				}
			}
			ofs.close();
		}
	}
}

void MixModel::rot_random() {
	auto & rng=NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0);
	QuaternionCrd qc(rng,0);
	Rotation rt(qc,center_);
	for(auto & cs:nghs_) {
		for(XYZ & c:cs.second) {
			//double d1=(c-center_).squarednorm();
			//std::cout <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<' ';
			rt.apply(&c);
			//std::cout <<c.x_ <<' ' <<c.y_ <<' ' <<c.z_ <<"\t\t";
			//double d2=(c-center_).squarednorm();
			//assert(fabs(d1-d2)<0.001);
			//std::cout <<d1 <<' ' <<d2 <<'\t';
		}//std::cout <<std::endl;
	}
	for(auto & cs:crds_) {
		for(XYZ & c:cs.second) rt.apply(&c);
	}
}

void MixModel::statistic(int ntime, std::string outpath, bool calatom, bool calwater) {
	microenvironment(outpath,calatom,false);
	if(!calwater) return;
	for(int i=0;i<ntime;i++) {
		//std::cout <<"State : " <<i <<std::endl;
		rot_random();
		microenvironment(outpath,false,calwater);
	}
}

std::map<MixModel::APS,std::vector<std::vector<XYZ>>> MixModel::statistic(int ntime, bool calatom, bool calwater) {
	std::map<APS,std::vector<std::vector<XYZ>>> atompairs;
	microenvironment(atompairs,calatom,false);
	if(!calwater) return atompairs;
	for(int i=0;i<ntime;i++) {
		//std::cout <<"State : " <<i <<std::endl;
		rot_random();
		microenvironment(atompairs,false,calwater);
	}
	return atompairs;
}

void MixModel::add_square(const XYZ &crd, XYZ &ori) {
	//XYZ ori(0.0,0.0,0.0);
	double half=Square_Length/2.0;
	while(ori<crd || crd<ori) {
		if(crd.x_<ori.x_-half) ori.x_ -= Square_Length;
		else if(crd.x_>=ori.x_+half) ori.x_ += Square_Length;
		if(crd.y_<ori.y_-half) ori.y_ -= Square_Length;
		else if(crd.y_>=ori.y_+half) ori.y_ += Square_Length;
		if(crd.z_<ori.z_-half) ori.z_ -= Square_Length;
		else if(crd.z_>=ori.z_+half) ori.z_ += Square_Length;
		//std::cout <<'\t' <<ori.x_ <<'\t' <<ori.y_ <<'\t' <<ori.z_ <<std::endl;
		//std::cout <<'\t' <<crd.x_ <<'\t' <<crd.y_ <<'\t' <<crd.z_ <<std::endl;
	}//std::cout <<"out " <<std::endl;
	ranges_.insert({ori,Content()});
}

void MixModel::partion(bool calwater, double buff_rad) {
	if(calwater) {
		XYZ ori(0.0,0.0,0.0);
		for(auto &r:ranges_) r.second.wats.clear();
		for(const XYZ & w:waters_) {
			if(ranges_.find(w)==ranges_.end()) {
				add_square(w,ori);
			}
			ranges_.at(w).wats.push_back(w);
		}
	}

	for(auto &r:ranges_) r.second.nghs.clear();
	//XYZ ori(0.0,0.0,0.0);
	for(const auto &nh:nghs_) {
		for(int i=0;i<nh.second.size();i+=3) {
			//if(ranges_.find(nh.second[i])==ranges_.end()) add_square(nh.second[i],ori);
			auto & rangenghs = ranges_.at(nh.second[i]).nghs;
			if(rangenghs.find(nh.first)==rangenghs.end()) rangenghs.insert({nh.first,std::vector<XYZ>()});
			rangenghs.at(nh.first).push_back(nh.second[i]);
			rangenghs.at(nh.first).push_back(nh.second[i+1]);
			rangenghs.at(nh.first).push_back(nh.second[i+2]);
		}
	}

	for(auto &r:ranges_) r.second.w_label.assign(r.second.wats.size(),0); //no clash
	double atom_max_rad=2.0;
	double max_rad=atom_max_rad+Water_Radius;
	int n_add=2;
	double l=Square_Length/2.0;
	while(l<max_rad) {
		l += Square_Length;
		n_add++;
	}
	//int num=0;
	for(auto &rg:ranges_) {//std::cout <<ranges_.size() <<'\t' <<num++ <<std::endl;
		auto & ns=rg.second.nghs;
		for(auto & ng:ns) {
			double dis=AtomPair::radius.at(ng.first)+Water_Radius-buff_rad;
			for(int i=0;i<ng.second.size();i+=3) {
				waterisclash(ng.second[i],n_add,dis);
			}
		}
	}
}

void MixModel::waterisclash(const XYZ & cen, int n_add, double dis) {
	Square st(cen.x_-(double)n_add*Square_Length,cen.y_-(double)n_add*Square_Length,cen.z_-(double)n_add*Square_Length);
	std::vector<Square> ss;
	int nlayer=n_add*2-1;
	for(int x=0;x<nlayer;x++) {
		for(int y=0;y<nlayer;y++) {
			for(int z=0;z<nlayer;z++) {
				ss.push_back(Square(st.x_+(double)x*Square_Length,st.y_+(double)y*Square_Length,st.z_+(double)z*Square_Length));
			}
		}
	}
	bool clash{false};
	for(auto &s:ss) {
		if(ranges_.find(s)==ranges_.end()) continue;
		auto & ws=ranges_.at(s).wats;
		auto & wl=ranges_.at(s).w_label;
		for(int i=0;i<ws.size();i++) {
			if(wl[i]==1) continue;
			double dx=ws[i].x_-cen.x_;
			if(dx>dis || dx<-dis) continue;
			double dy=ws[i].y_-cen.y_;
			if(dy>dis || dy<-dis) continue;
			double dz=ws[i].z_-cen.z_;
			if(dz>dis || dz<-dis) continue;
			double d2=(ws[i]-cen).squarednorm();
			if(d2<dis*dis) wl[i]=1;
		}
	}
}

void MixModel::nearregion(std::vector<XYZ> &atoms, const XYZ &cen, double max_rad, int n_add, double atomrad, double buffrad) {
	Square st(cen.x_-(double)n_add*Square_Length,cen.y_-(double)n_add*Square_Length,cen.z_-(double)n_add*Square_Length);
	std::vector<Square> ss;
	int nlayer=n_add*2-1;
	for(int x=0;x<nlayer;x++) {
		for(int y=0;y<nlayer;y++) {
			for(int z=0;z<nlayer;z++) {
				ss.push_back(Square(st.x_+(double)x*Square_Length,st.y_+(double)y*Square_Length,st.z_+(double)z*Square_Length));
			}
		}
	}
	double dm=max_rad+buffrad;
	double dm2=dm*dm;
	//double da=atomrad+Water_Radius-buffrad;
	//double da2=da*da;
	///*
	for(auto &s:ss) {
		if(ranges_.find(s)==ranges_.end()) continue;
		auto & ws=ranges_.at(s).wats;
		auto & wl=ranges_.at(s).w_label;
		for(int i=0;i<ws.size();i++) {
			if(wl[i]==1) continue;
			auto &w=ws[i];
			double dx=w.x_-cen.x_;
			if(dx>dm || dx<-dm) continue;
			double dy=w.y_-cen.y_;
			if(dy>dm || dy<-dm) continue;
			double dz=w.z_-cen.z_;
			if(dz>dm || dz<-dm) continue;
			double d22=(w-cen).squarednorm();
			//if(d22<da2) continue;
			if(d22>dm2) continue;
			atoms.push_back(w);
		}
	}//*/

	/*
	for(XYZ &w:waters_) {
		double dx=w.x_-cen.x_;
		if(dx>dm || dx<-dm) continue;
		double dy=w.y_-cen.y_;
		if(dy>dm || dy<-dm) continue;
		double dz=w.z_-cen.z_;
		if(dz>dm || dz<-dm) continue;
		double d22=(w-cen).squarednorm();
		if(d22<da2) continue;
		if(d22>dm2) continue;
		atoms.push_back(w);
	}*/
}

void MixModel::watermenv(std::map<APS,std::vector<std::vector<XYZ>>> &atompairs, double max_rad, double buffrad) {
	int n_add=3;
	double l=Square_Length/2.0;
	while(l<max_rad) {
		l += Square_Length;
		n_add++;
	}
	for(auto & ns:nghs_) {
		AtomPair::AtomName cen=ns.first;
		if(AtomPair::neighbors.find(cen)==AtomPair::neighbors.end()) {
			std::cout <<"Neighbor is not in the definition!" <<std::endl;
			std::cout <<'\t' <<cen.first <<'\t' <<cen.second <<std::endl;
			exit(1);
		}
		auto nghfise=AtomPair::neighbors.at(cen);
		std::string nghfi=nghfise.first.second;
		std::string nghse=nghfise.second.second;
		std::vector<XYZ> & nghcrds = ns.second;
		std::vector<std::vector<XYZ>> envs;//first 3 are cen, fi, se, others are waters
		for(int i=0;i<nghcrds.size();i+=3) {
			std::vector<XYZ> p1;
			p1.push_back(nghcrds[i]);
			p1.push_back(nghcrds[i+1]);
			p1.push_back(nghcrds[i+2]);
			nearregion(p1,nghcrds[i],max_rad,n_add,AtomPair::radius.at(cen),buffrad);
			envs.push_back(p1);
		}
		APS aps={cen,{"HOH","O"}};
		if(atompairs.find(aps)==atompairs.end()) atompairs.insert({aps,std::vector<std::vector<XYZ>>()});
		for(auto &es:envs) {
			atompairs.at(aps).push_back(es);
		}
	}
}


bool MixModel::getrt(const XYZ &crd, const XYZ &axis, Rotation &rt) {
	XYZ ori(0.0,0.0,0.0);
	double ang = angle(crd,ori,axis) * 180.0/3.14159265;//std::cout <<"\t" <<ang <<std::endl;
	if(ang>-0.000001 && ang<0.000001) return false;
	XYZ ax=cross(crd,axis);
	QuaternionCrd qc(ax,ang);
	rt=Rotation(qc,ori);

	XYZ c1=crd;
	rt.apply(&c1);
	double nang=angle(c1,ori,axis) * 180.0/3.14159265;//std::cout <<nang <<std::endl;
	if(nang<0.0001 && nang>-0.0001) return true;
	ax=-ax;
	qc=QuaternionCrd(ax,ang);
	rt=Rotation(qc,ori);
	c1=crd;
	rt.apply(&c1);
	nang=angle(c1,ori,axis) * 180.0/3.14159265;//std::cout <<nang <<std::endl;
	if(!(nang<0.0001 && nang>-0.0001)) {
		std::cout <<"Bad Data" <<std::endl;
		exit(1);
	}
	return true;
}

bool MixModel::getrt(const XYZ &crd, const XYZ &axis, const XYZ &plane, Rotation &rt) {
	XYZ ori(0.0,0.0,0.0);
	XYZ f1=cross(axis,plane);
	XYZ f2=cross(axis,crd);
	double ang = angle(f1,ori,f2) * 180.0/3.14159265;//std::cout <<"\t" <<ang <<std::endl;
	if(ang>-0.000001 && ang<0.000001) return false;
	XYZ ax=axis;
	QuaternionCrd qc(ax,ang);
	rt=Rotation(qc,ori);

	XYZ c1=crd;
	rt.apply(&c1);
	f2=cross(axis,c1);
	double nang=angle(f1,ori,f2) * 180.0/3.14159265;//std::cout <<nang <<std::endl;
	if(nang<0.0001 && nang>-0.0001) return true;
	ax=-ax;
	qc=QuaternionCrd(ax,ang);
	rt=Rotation(qc,ori);
	c1=crd;
	rt.apply(&c1);
	f2=cross(axis,c1);
	nang=angle(f1,ori,f2) * 180.0/3.14159265;//std::cout <<nang <<std::endl;
	if(!(nang<0.0001 && nang>-0.0001)) {
		std::cout <<"Bad Data" <<std::endl;
		exit(1);
	}
	return true;
}

void MixModel::rot(std::map<MixModel::APS,std::vector<std::vector<NSPgeometry::XYZ>>> &atompairs) {
	XYZ xmins(-1.0,0.0,0.0);
	XYZ ymins(0.0,-1.0,0.0);
	for(auto &ap:atompairs) {
		for(auto &cs:ap.second) {
			XYZ ori=cs[0];
			for(XYZ &c:cs) c=c-ori;
			Rotation rt;
			if(getrt(cs[1],xmins,rt)) {
				for(XYZ &c:cs) rt.apply(&c);
			}
			if(getrt(cs[2],xmins,ymins,rt)) {
				for(XYZ &c:cs) rt.apply(&c);
			}
		}
	}
}

void MixModel::microenvironment(std::map<APS,std::vector<std::vector<XYZ>>> &atompairs, bool calatom, bool calwater) {
	for(auto & ns:nghs_) {
		AtomPair::AtomName cen=ns.first;
		if(AtomPair::neighbors.find(cen)==AtomPair::neighbors.end()) {
			std::cout <<"Neighbor is not in the definition!" <<std::endl;
			std::cout <<'\t' <<cen.first <<'\t' <<cen.second <<std::endl;
			exit(1);
		}
		auto nghfise=AtomPair::neighbors.at(cen);
		std::string nghfi=nghfise.first.second;
		std::string nghse=nghfise.second.second;
		std::vector<XYZ> & nghcrds = ns.second;
		if(calatom) {
			for(auto &cs:crds_) {
				std::vector<std::pair<std::vector<XYZ>,std::vector<XYZ>>> envs;
				std::vector<XYZ> & envcrds = cs.second;
				for(int i=0;i<nghcrds.size();i+=3) {
					std::vector<XYZ> p1;
					p1.push_back(nghcrds[i]);
					p1.push_back(nghcrds[i+1]);
					p1.push_back(nghcrds[i+2]);
					std::vector<XYZ> p2;
					for(int j=0;j<envcrds.size();j++) {
						double x=envcrds[j].x_-nghcrds[i].x_;
						if(x<-rad_add || x>rad_add) continue;
						double y=envcrds[j].y_-nghcrds[i].y_;
						if(y<-rad_add || y>rad_add) continue;
						double z=envcrds[j].z_-nghcrds[i].z_;
						if(z<-rad_add || z>rad_add) continue;
						double d2=x*x+y*y+z*z;
						if(d2<0.001 || d2>rad_add*rad_add) continue;
						p2.push_back(envcrds[j]);
					}
					if(p2.empty()) continue;
					envs.push_back({p1,p2});
				}
				if(envs.empty()) continue;
				APS aps={cen,cs.first};
				if(atompairs.find(aps)==atompairs.end()) atompairs.insert({aps,std::vector<std::vector<XYZ>>()});
				for(auto &es:envs) {
					std::vector<XYZ> totcrds=es.first;
					for(XYZ & c:es.second) totcrds.push_back(c);
					atompairs.at(aps).push_back(totcrds);
				}
			}
		}
		if(calwater) {
			std::vector<std::pair<std::vector<XYZ>,std::vector<XYZ>>> envs;
			for(int i=0;i<nghcrds.size();i+=3) {
				std::vector<XYZ> p1;
				p1.push_back(nghcrds[i]);
				p1.push_back(nghcrds[i+1]);
				p1.push_back(nghcrds[i+2]);
				std::vector<XYZ> p2;
				for(int j=0;j<waters_.size();j++) {
					double x=waters_[j].x_-nghcrds[i].x_;
					if(x<-rad_add || x>rad_add) continue;
					double y=waters_[j].y_-nghcrds[i].y_;
					if(y<-rad_add || y>rad_add) continue;
					double z=waters_[j].z_-nghcrds[i].z_;
					if(z<-rad_add || z>rad_add) continue;
					double d2=x*x+y*y+z*z;
					double dmin=Water_Radius+AtomPair::radius.at(ns.first);
					if(d2<0.001 || d2>rad_add*rad_add) continue;
					p2.push_back(waters_[j]);
				}
				if(p2.empty()) continue;
				envs.push_back({p1,p2});
			}
			if(envs.empty()) continue;
			APS aps={cen,{"HOH","O"}};
			if(atompairs.find(aps)==atompairs.end()) atompairs.insert({aps,std::vector<std::vector<XYZ>>()});
			for(auto &es:envs) {
				std::vector<XYZ> totcrds=es.first;
				for(XYZ & c:es.second) totcrds.push_back(c);
				atompairs.at(aps).push_back(totcrds);
			}
		}
	}
}

void MixModel::print(std::string fn, std::map<APS,std::vector<std::vector<XYZ>>> &atompairs, bool fixedcrd) {
	std::ofstream ofs(fn);
	for(auto &ap:atompairs) {
		if(fixedcrd) {
			for(auto &cs:ap.second) {
				ofs <<ap.first.first.first <<' ' <<ap.first.first.second <<'\t' <<ap.first.second.first
						<<' ' <<ap.first.second.second <<std::endl;
				ofs <<cs.size()*3 <<std::endl;
				for(XYZ &c:cs) {
					ofs <<c.x_ <<std::endl;
					ofs <<c.y_ <<std::endl;
					ofs <<c.z_ <<std::endl;
				}
			}
		} else {
			ofs <<ap.first.first.first <<' ' <<ap.first.first.second <<'\t' <<ap.first.second.first
					<<' ' <<ap.first.second.second <<std::endl;
			std::vector<double> ds;
			for(auto &cs:ap.second) {
				for(int i=3;i<cs.size();i++) {
					ds.push_back(cs[i].x_);
					ds.push_back(cs[i].y_);
					ds.push_back(cs[i].z_);
				}
			}
			ofs <<ds.size() <<std::endl;
			for(double &d:ds) ofs <<d <<std::endl;
		}
	}
	ofs.close();
}

std::vector<std::pair<PdbReader_xy::OneAtom,PdbReader_xy::OneAtom>> PdbReader_xy::hbpairs() const {
	std::vector<std::pair<OneAtom,OneAtom>> hbps;
	std::map<AtomPair::AtomName,double> rads=AtomPair::radius;
	std::set<AtomPair::AtomName> drs1=AtomPair::donors;
	std::set<AtomPair::AtomName> acs1=AtomPair::accepters;
	std::map<std::string,std::set<std::string>> drs;
	std::map<std::string,std::set<std::string>> acs;
	for(auto &a:drs1) {
		if(drs.find(a.first)==drs.end()) drs.insert({a.first,std::set<std::string>()});
		drs.at(a.first).insert(a.second);
	}
	for(auto &a:acs1) {
		if(acs.find(a.first)==acs.end()) acs.insert({a.first,std::set<std::string>()});
		acs.at(a.first).insert(a.second);
	}
	struct OneRes {
		int chainid;
		int resid;
		Residue res;
		OneRes(int c,int r,Residue r2):chainid(c),resid(r),res(r2) {;}
	};
	std::vector<OneRes> ress;
	for(int i=0;i<protein_.size();i++) {
		for(int j=0;j<protein_[i].size();j++) {
			ress.push_back(OneRes(i,j,protein_[i][j]));
		}
	}
	std::set<std::string> polarscs;
	for(auto &ac:acs) if(ac.second.size()!=1 && ac.first!="SER" && ac.first!="THR") polarscs.insert(ac.first);
	for(int i=0;i<ress.size();i++) {
		std::string rn1=ress[i].res.resname();
		std::map<std::string,XYZ>& cs1=ress[i].res.crds();
		if(polarscs.find(rn1)!=polarscs.end()) {
			for(auto &c1:cs1) {
				if(drs1.find({rn1,c1.first})==drs1.end() && acs1.find({rn1,c1.first})==acs1.end()) continue;
				if(c1.first=="N") continue;
				if(c1.first=="O") continue;
				double d2=(c1.second-cs1.at("N")).squarednorm();
				double dc=rads.at({rn1,c1.first});
				double dn=rads.at({rn1,"N"});
				double d2n=(dc+dn)*(dc+dn);
				if(d2<d2n) hbps.push_back({OneAtom(ress[i].chainid,ress[i].resid,"N"),
					OneAtom(ress[i].chainid,ress[i].resid,c1.first)});
				d2=(c1.second-cs1.at("O")).squarednorm();
				double doa=rads.at({rn1,"O"});
				double d2o=(dc+doa)*(dc+doa);
				if(d2<d2o) hbps.push_back({OneAtom(ress[i].chainid,ress[i].resid,"O"),
					OneAtom(ress[i].chainid,ress[i].resid,c1.first)});
			}
		}
		for(int j=i+1;j<ress.size();j++) {
			std::string rn2=ress[j].res.resname();
			std::map<std::string,XYZ>& cs2=ress[j].res.crds();
			for(auto &c1:cs1) {
				AtomPair::AtomName an1={rn1,c1.first};
				bool bd1=drs1.find(an1)==drs1.end();
				bool ba1=acs1.find(an1)==acs1.end();
				if(bd1 && ba1) continue;
				double r1=rads.at({rn1,c1.first});
				for(auto &c2:cs2) {
					if(j-i==1) {
						if(c1.first=="O" && c2.first=="N" || c1.first=="N" && c2.first=="O") {
							double dd=(cs1.at("C")-cs2.at("N")).squarednorm();
							if(dd<4.0) continue;
						}
					}
					AtomPair::AtomName an2={rn2,c2.first};
					bool bd2=drs1.find(an2)==drs1.end();
					bool ba2=acs1.find(an2)==acs1.end();
					if((bd1||ba2)&&(bd2||ba1)) continue;
					double r2=rads.at({rn2,c2.first});
					double dr2=(r2+r1)*(r2+r1);
					double d2=(c1.second-c2.second).squarednorm();
					if(d2>dr2) continue;
					hbps.push_back({OneAtom(ress[i].chainid,ress[i].resid,c1.first),
							OneAtom(ress[j].chainid,ress[j].resid,c2.first)});
				}
			}
		}
	}
	return hbps;
}

std::vector<std::string> PdbReader_xy::seqs() {
	std::vector<std::string> sqs;
	for(Chain &ch:protein_) {
		std::string s;
		for(Residue &r:ch) {
			char c='X';
			if(AA20.find(r.resname())!=AA20.end()) c=AA20.at(r.resname());
			s.push_back(c);
		}
		sqs.push_back(s);
	}
	return sqs;
}














/*



void PdbReader_xy::getpar(std::vector<XYZ>&fixeds,std::vector<double>&rads,XYZ &center,double &maxdis,double wr,double ar) const {
	static std::map<std::pair<std::string,std::string>,double> rds;
	if(rds.empty()) {
		std::string fn=NSPdataio::datafilename("resLib.dat");
		std::ifstream ifs(fn);
		std::string line;
		while(std::getline(ifs,line)) {
			std::istringstream iss(line);
			std::string s1,s2;
			int nline;
			iss >>s1 >>s2 >>nline;
			for(int i=0;i<nline;i++) {
				std::getline(ifs,line);
				std::istringstream ss(line);
				std::string s3,s4;
				double d;
				ss >>s3 >>s4 >>d;
				rds.insert({{s2,s4},d});
			}
		}
		ifs.close();
	}
	for(const Chain &ch:protein_) {
		for(const Residue &res:ch) {
			std::string resname = res.resname();
			std::map<std::string,XYZ> cs=res.crds();
			for(auto & c:cs) {
				if(rds.find({resname,c.first})==rds.end()) continue;
				rads.push_back(rds.at({resname,c.first}));
				fixeds.push_back(c.second);
			}
		}
	}
	for(const auto & prs:ligands_) {
		for(const PdbRecord &p:prs) {
			fixeds.push_back(XYZ(p.x,p.y,p.z));
			rads.push_back(ar);
		}
	}
	for(const PdbRecord &p:waters_) {
		fixeds.push_back(XYZ(p.x,p.y,p.z));
		rads.push_back(wr);
	}
	maxdis = 0.0;
	center = NSPgeometry::XYZ(0.0,0.0,0.0);
	std::vector<XYZ> &prds = fixeds;
	for(int i=0;i<prds.size();i++) {
		center.x_ += prds[i].x_;
		center.y_ += prds[i].y_;
		center.z_ += prds[i].z_;
		for(int j=i+1;j<prds.size();j++) {
			double x = prds[i].x_-prds[j].x_;
			double y = prds[i].y_-prds[j].y_;
			double z = prds[i].z_-prds[j].z_;
			double d2 = x*x + y*y + z*z;
			if(d2>maxdis) maxdis = d2;
		}
	}
	center.x_ /= (double)(prds.size());
	center.y_ /= (double)(prds.size());
	center.z_ /= (double)(prds.size());
	maxdis = sqrt(maxdis)/2.0;
}

void PdbReader_xy::getpar_grid(std::vector<XYZ>&fixeds,std::vector<double>&rads,
		XYZ &center,double &maxdis) const {
	static std::map<std::pair<std::string,std::string>,double> rds;
	if(rds.empty()) {
		std::string fn=NSPdataio::datafilename("resLib.dat");
		std::ifstream ifs(fn);
		std::string line;
		while(std::getline(ifs,line)) {
			std::istringstream iss(line);
			std::string s1,s2;
			int nline;
			iss >>s1 >>s2 >>nline;
			for(int i=0;i<nline;i++) {
				std::getline(ifs,line);
				std::istringstream ss(line);
				std::string s3,s4;
				double d;
				ss >>s3 >>s4 >>d;
				rds.insert({{s2,s4},d});
			}
		}
		ifs.close();
	}
	int natom=0;
	maxdis = 0.0;
	center = NSPgeometry::XYZ(0.0,0.0,0.0);
	for(const Chain &ch:protein_) {
		for(const Residue &res:ch) {
			std::string resname = res.resname();
			std::map<std::string,XYZ> cs=res.crds();
			for(auto & c:cs) {
				if(rds.find({resname,c.first})==rds.end()) continue;
				rads.push_back(rds.at({resname,c.first}));
				fixeds.push_back(c.second);
				center.x_ += c.second.x_;
				center.y_ += c.second.y_;
				center.z_ += c.second.z_;
				natom++;
			}
		}
	}
	center.x_ /= (double)(natom);
	center.y_ /= (double)(natom);
	center.z_ /= (double)(natom);
	for(const auto & prs:ligands_) {
		for(const PdbRecord &p:prs) {
			fixeds.push_back(XYZ(p.x,p.y,p.z));
			rads.push_back(Atom_Radius);
		}
	}

	std::vector<XYZ> &prds = fixeds;
	for(int i=0;i<prds.size();i++) {
		for(int j=i+1;j<prds.size();j++) {
			double x = prds[i].x_-prds[j].x_;
			double y = prds[i].y_-prds[j].y_;
			double z = prds[i].z_-prds[j].z_;
			double d2 = x*x + y*y + z*z;
			if(d2>maxdis) maxdis = d2;
		}
	}
	maxdis = sqrt(maxdis)/2.0;
}

void MixModel::rewater() {
	iswt.resize(waters_.size(),true);
	for(XYZ &w:waters_) {
		for(int i=0;i<fixeds_.size();i++) {
			double d=sqrt((w-fixeds_[i]).squarednorm());
			if(d<rads_[i]) {
				iswt[i]=false;
				break;
			}
		}
	}
}




















void NSPallatom::definesdwatercontrol(std::string name,const std::vector<std::string> &controllines){
	std::map<std::string,double> doublepars{ {"Water_Radius",0.0}, {"Atom_Radius",0.0}, {"TimeStep",0.0},
			{"Force_Center",0.0}, {"Gamma",0.0}, {"Radius_Add",0.0} };
	std::map<std::string,std::vector<std::string>> stringvecpars{  };
	std::map<std::string,std::vector<double>> doublevecpars{ {"Temperature",{}} };
	std::map<std::string,std::string> stringpars{ {"PDBFile",""}, {"OutPath",""} };
	std::map<std::string,std::vector<int>> intvecpars{{"Opt_Step",{}}};
	std::map<std::string,int>intpars{ {"Run_Step",0}, {"Neighbor_Step",0}, {"Trajectory_Step",0}, {"NResults",0} };
	SDWaterControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,stringvecpars,intvecpars);
	int nsuccess=SDWaterControls::adjustvalues(name,controllines);
	if(nsuccess!= controllines.size()) {
		exit(1);
	}
}

void NSPallatom::readcontrols_sdwater(const std::string &filename,std::string name){
	NSPdataio::ControlFile cf;
	cf.readfile(filename);
	std::vector<std::string> sdcontrollines=cf.getcontrolines("SDWater");
	definesdwatercontrol(name+"_sdwater",sdcontrollines);
}




SDWater::SDWater(NSPdataio::ParameterSet &pset):pset_(pset) {
	std::string pdbfile;
	pset_.getval("PDBFile",&pdbfile);
	pr_.init(pdbfile);
	double wr,ar;
	pset_.getval("Water_Radius",&wr);
	pset_.getval("Atom_Radius",&ar);
	pr_.getpar(fixeds_,rads_,center_,maxrad_,wr,ar);
	add_water();
	sd_ = std::shared_ptr<NSPsd::StochasticDynamics>(new NSPsd::StochasticDynamics);
	double gamma;
	pset_.getval("Gamma",&gamma);
	std::vector<double> gammas;
	std::vector<double> mas;
	for(int i=0;i<fixeds_.size()+waters_.size();i++) {
		masses_.push_back(18.0);
		masses_.push_back(18.0);
		masses_.push_back(18.0);
		mas.push_back(18.0);
		gammas.push_back(gamma);
		forceoff_.push_back(false);
		if(i<fixeds_.size()) forceoff_.back()=true;
	}
	double timestep;
	pset_.getval("TimeStep",&timestep);
	std::vector<double> temp;
	pset_.getval("Temperature",&temp);
	temp[0] *= KBT;
	sd_->init(mas,gammas,timestep,temp[0],3);
	std::vector<double> initcrds;
	for(const XYZ & c:fixeds_) {
		initcrds.push_back(c.x_);
		initcrds.push_back(c.y_);
		initcrds.push_back(c.z_);
	}
	for(const XYZ & c:waters_) {
		initcrds.push_back(c.x_);
		initcrds.push_back(c.y_);
		initcrds.push_back(c.z_);
	}
	for(double & d:initcrds) d *= 0.1;
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	state_ = sd_->make_initstate(initcrds, rng, &masses_);
	buffstate_ = std::shared_ptr < NSPsd::StochasticDynamics::State > (new NSPsd::StochasticDynamics::State);
	*buffstate_ = *state_;
	potenes_.resize(2);
}

void SDWater::add_water() {
	double rad_add;
	double atom_rad;
	double water_rad;
	pset_.getval("Radius_Add",&rad_add);
	pset_.getval("Atom_Radius",&atom_rad);
	pset_.getval("Water_Radius",&water_rad);
	double bigrad=maxrad_+rad_add;
	double totvol = bigrad*bigrad*bigrad;
	double atomvol = water_rad*water_rad*water_rad;
	int nwaters = (int)(totvol/atomvol)-fixeds_.size();
	int numinline = (int)(maxrad_/(water_rad+0.1));
	int nlayer = nwaters / (6*numinline*numinline) + 1;

	std::cout <<"Fixed Number Is : " <<fixeds_.size() <<std::endl;
	std::cout <<"Water Number Is : " <<numinline*numinline*6*nlayer <<std::endl;
	std::cout <<"Total Number Is : " <<numinline*numinline*6*nlayer+fixeds_.size() <<std::endl;
	std::cout <<"Min Total Number Is : " <<(int)(totvol/atomvol) <<std::endl;
	std::cout <<"Max Radius is : " <<maxrad_ <<std::endl;
	std::cout <<"Water Number In One Line Is : " <<numinline <<std::endl;
	std::cout <<"Number Of Water Layer Is : " <<nlayer <<std::endl;

	for(int i=0;i<numinline;i++) {
		double li = (water_rad+0.1)*2.0*(double)i-maxrad_;
		for(int j=0;j<numinline;j++) {
			double lj = (water_rad+0.1)*2.0*(double)j-maxrad_;
			for(int r=0;r<nlayer;r++) {
				double rad = maxrad_ + (double)r*2.0*(water_rad+0.1) + 1.0;
				waters_.push_back(XYZ(center_.x_-rad,center_.y_+li,center_.z_+lj));
				waters_.push_back(XYZ(center_.x_+rad,center_.y_+li,center_.z_+lj));
				waters_.push_back(XYZ(center_.x_+li,center_.y_-rad,center_.z_+lj));
				waters_.push_back(XYZ(center_.x_+li,center_.y_+rad,center_.z_+lj));
				waters_.push_back(XYZ(center_.x_+li,center_.y_+lj,center_.z_-rad));
				waters_.push_back(XYZ(center_.x_+li,center_.y_+lj,center_.z_+rad));
			}
		}
	}
	assert(waters_.size()==numinline*numinline*6*nlayer);
}

double SDWater::stericene(const XYZ &c1,const XYZ &c2,std::vector<XYZ>*fs) {
	fs->assign(2, (0.0, 0.0, 0.0));
	double r = distance(c1, c2, fs);
	double r2 = r * r;
	double r3 = r2 * r;
	double r5 = r2 * r3;
	double r6 = r3 * r3;
	double r7 = r * r6;
	double r12 = r6 * r6;
	double r13 = r * r12;
	double sigma = 3.0*A2NM;
	double eps = 0.4;
	double sigma3 = sigma * sigma * sigma;
	double sigma6 = sigma3 * sigma3;
	double sigma12 = sigma6 * sigma6;
	double ene = 4.0 * eps * (sigma12 / r12 - sigma6 / r6) + eps;
	double dedr = -4.0 * eps * (-12.0 * sigma12 / r13 + 6.0 * sigma6 / r7);
	if (dedr > 500.0) dedr = 500.0;
	(*fs)[0] = dedr * (*fs)[0];
	(*fs)[1] = dedr * (*fs)[1];
	return ene;
}

std::vector<double> SDWater::forces(const std::vector<double> &cs, const std::set<std::pair<int,int>>&nbl, double f_cen, double discut) {
	double ene=0;
	std::vector<XYZ> xyzf(cs.size(),XYZ(0.0,0.0,0.0));
	for(const auto &nl:nbl) {
		if(forceoff_[nl.first] && forceoff_[nl.second]) continue;
		std::vector<XYZ> cf;
		XYZ c1=XYZ(cs[nl.first*3],cs[nl.first*3+1],cs[nl.first*3+2]);
		XYZ c2=XYZ(cs[nl.second*3],cs[nl.second*3+1],cs[nl.second*3+2]);
		ene+=stericene(c1,c2,&cf);
		if(!forceoff_[nl.first]) xyzf[nl.first] = xyzf[nl.first] + cf[0];
		if(!forceoff_[nl.second]) xyzf[nl.second] = xyzf[nl.second] + cf[1];
	}
	potenes_[0]=ene;
	ene=0.0;
	NSPgeometry::XYZ cen=center_*0.1;
	for(int i=0;i<forceoff_.size();i++) {
		if(forceoff_[i]) continue;
		XYZ c(cs[i*3],cs[i*3+1],cs[i*3+2]);
		XYZ f=cen-c;
		double d2=f.squarednorm();
		if(d2<discut*discut) continue;
		double d=sqrt(d2);
		f = f / d * f_cen;
		xyzf[i] = xyzf[i]+f;
		ene += f_cen * (d-discut);
	}
	potenes_[1] = ene;
	std::vector<double> fs;
	for(XYZ & c:xyzf) {
		fs.push_back(c.x_);
		fs.push_back(c.y_);
		fs.push_back(c.z_);
	}
	return fs;
}

std::set<std::pair<int,int>> SDWater::getneighbor(double rcut2) {
	std::set<std::pair<int,int>> nbl;
	std::vector<double> & crds = state_->crd;
	for(int m=0;m<crds.size()/3;m++) {
		for(int n=m+1;n<crds.size()/3;n++) {
			if(forceoff_[m] && forceoff_[n]) continue;
			double x=crds[m*3]-crds[n*3];
			double y=crds[m*3+1]-crds[n*3+1];
			double z=crds[m*3+2]-crds[n*3+2];
			double d2 = x*x + y*y + z*z;
			if(d2>rcut2) continue;
			nbl.insert({m,n});
		}
	}
	return nbl;
}

void SDWater::runstep(int nstep, bool isopt, int itvl, std::string path) {
	double rcut2=1.0;//nm^2
	std::set<std::pair<int,int>> nbl;
	auto & rng = NSPdstl::RandomEngine<>::getinstance();
	std::vector<bool> fcff;
	for(bool b:forceoff_) {
		fcff.push_back(b);
		fcff.push_back(b);
		fcff.push_back(b);
	}
	double f_cen;
	pset_.getval("Force_Center",&f_cen);
	double dcut;
	pset_.getval("Radius_Add",&dcut);
	double fd = (dcut+maxrad_)*0.1;
	for(int i=0;i<nstep;i++) {
		if(i%20==0) nbl=getneighbor(rcut2);
		std::vector<double> fs = forces(state_->crd,nbl,f_cen,fd);
		if(!isopt) {
			if((i+1)%itvl==0) {
				std::cout <<"Run Step : " <<i+1 <<std::endl;
				std::vector<double> &crds = state_->crd;
				for(int i=3*fixeds_.size(),j=0;i<crds.size();i+=3,j++) {
					waters_[j].x_ = crds[i] * 10.0;
					waters_[j].y_ = crds[i+1] * 10.0;
					waters_[j].z_ = crds[i+2] * 10.0;
				}
				//printtraj(path+std::to_string((i+1)/itvl)+".pdb");
				printdis(path+std::to_string((i+1)/itvl)+".dis");
				printene(path+"ene1.dat",path+"ene2.dat");
			}
		}
		sd_->leapstep(*state_,*buffstate_,rng,fs,&masses_,fcff);
		auto temp = state_;
		state_ = buffstate_;
		buffstate_ = temp;
	}
	std::vector<double> &crds = state_->crd;
	for(int i=3*fixeds_.size(),j=0;i<crds.size();i+=3,j++) {
		waters_[j].x_ = crds[i] * 10.0;
		waters_[j].y_ = crds[i+1] * 10.0;
		waters_[j].z_ = crds[i+2] * 10.0;
	}
}

void SDWater::move2center() {
	int nstep;
	pset_.getval("Run_Step",&nstep);
	int itvl;
	pset_.getval("Trajectory_Step",&itvl);
	std::string outpath;
	pset_.getval("OutPath",&outpath);
	if(outpath.back()!='/') outpath += '/';
	std::vector<double> temp;
	pset_.getval("Temperature",&temp);
	double scale=temp[0]*KBT/sd_->kbT()[0];
	state_->scaletemp(scale);
	sd_->scaletemperatures(scale);
	runstep(nstep, false, itvl, outpath);
}

void SDWater::opt() {
	std::vector<int> step;
	pset_.getval("Opt_Step",&step);
	std::string outpath;
	pset_.getval("OutPath",&outpath);
	if(outpath.back()!='/') outpath += '/';
	std::vector<double> temp;
	pset_.getval("Temperature",&temp);
	int nresult;
	pset_.getval("NResult",&nresult);
	double scale;
	for(int i=0;i<nresult;i++) {
		scale=temp[0]*KBT/sd_->kbT()[0];
		state_->scaletemp(scale);
		sd_->scaletemperatures(scale);
		runstep(step[0], true, 0, outpath);
		scale=temp[1]*KBT/sd_->kbT()[0];
		state_->scaletemp(scale);
		sd_->scaletemperatures(scale);
		runstep(step[1], true, 0, outpath);
		printtraj(outpath+std::to_string(i)+".pdb.opt");
	}
}

void SDWater::printtraj(std::string fn) {
	std::ofstream ofs(fn);
	int nw=pr_.print(ofs);
	int nr=0;
	for(NSPgeometry::XYZ &c:waters_) {
		if(nw==9998) nw=-1;
		if(nr==9999) nr=0;
		NSPproteinrep::PdbRecord p;
		p.label = "HETATM";
		p.namesymbol = "O";
		p.residuename = "HOH";
		p.chainid = 'A';
		p.atomid = ++nw;
		p.residueid = nr++;
		p.elementname[1] = 'O';
		p.x = c.x_;
		p.y = c.y_;
		p.z = c.z_;
		ofs <<p.toString() <<std::endl;
	}
	ofs.close();
}

void SDWater::printdis(std::string fn) {
	std::ofstream ofs(fn);
	std::vector<double> ds=diswater();
	for(double &d:ds) ofs <<d <<std::endl;
	ofs.close();
}

std::vector<double> SDWater::diswater() {
	std::vector<double> dis;
	for(XYZ & c:waters_) {
		dis.push_back(NSPgeometry::distance(c,center_));
	}
	return dis;
}

void SDWater::printene(std::string fn1, std::string fn2) {
	std::ofstream ofs;
	ofs.open(fn1,std::ofstream::app);
	ofs <<potenes_[0] <<std::endl;
	ofs.close();
	ofs.open(fn2,std::ofstream::app);
	ofs <<potenes_[1] <<std::endl;
	ofs.close();
}





void MixModel::getneighbor() {
	std::vector<Residue> pep=pr_.pep();
	auto ngh=AtomPair::neighbors;
	for(Residue &res:pep) {
		std::string rn=res.resname();
		std::map<std::string,XYZ> cs=res.crds();
		for(auto &c:cs) {
			auto nh=ngh.at({rn,c.first});
			std::string fi=nh.first.second;
			std::string se=nh.second.second;
			if(cs.find(fi)==cs.end()) continue;
			if(cs.find(se)==cs.end()) continue;
			nhbs_.at({rn,c.first}).push_back(AtomNeighbor());
			AtomPair &an=nhbs_.at({rn,c.first}).back();
			an.resname = rn;
			an.cen = c.first;
			an.fi = fi;
			an.se = se;
			std::vector<XYZ> &ancs=an.crds;
			std::vector<std::pair<std::string,std::string>> &anan=an.names;
			for(int i=0;i<pep.size();i++) {
				std::string nghresname=pep[i].resname();
				std::map<std::string,XYZ> nghcs=pep[i].crds();
				for(auto & nc:nghcs) {
					if(nc.second-c.second).squarednorm()>
				}
			}
			for(int i=0;i<waters_.size();i++) {

			}
		}
	}
}


*/





