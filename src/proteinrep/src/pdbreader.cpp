/*
 * pdbreader.cpp
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */

#include "proteinrep/pdbreader.h"
#include "dataio/inputlines.h"
#include "dataio/datapaths.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/filesystem.hpp>

using namespace NSPproteinrep;
namespace BFS = boost::filesystem;
void PdbReader::readpdb(const std::string & filename) {
	NSPdataio::TextLines lines;
	lines.init(filename);
	readpdb(lines.lines());
}

void PdbReader::readpdb(std::vector<std::string> & lines) {
	for (auto & line : lines) {
		if (line.substr(0, 6) != "ATOM  " && line.substr(0, 6) != "HETATM")
			continue;
		addRecord(PdbRecord(line));
	}
	mappdbkeyint_=std::shared_ptr<MapPdbKeyInt>(new MapPdbKeyInt(records_));
	for( auto &c:records_) {
		aminoacidsequences_.insert(std::make_pair(c.first,std::vector<std::string>()));
		std::vector<std::string>&seq=aminoacidsequences_.at(c.first);
		int posi=0;
		for(auto &r:c.second){
//			if(r.second[0].label != "ATOM") continue;
			assert(posi++ == mappdbkeyint_->posiNumber(r.first,c.first));
			seq.push_back(r.second[0].residuename);
		}
	}
}

void PdbReader::addRecord(const PdbRecord & record) {
	char chainid = record.chainid;
	if (records_.find(chainid) == records_.end())
		records_.insert(std::make_pair(chainid, ResMapType()));
	ResMapType & resmap = records_.at(chainid);
	ResKeyType reskey = std::make_pair(record.residueid, record.insertionid);
	if (resmap.find(reskey) == resmap.end())
		resmap.insert(std::make_pair(reskey, std::vector<PdbRecord>()));
	resmap.at(reskey).push_back(record);
}

MapPdbKeyInt::MapPdbKeyInt(const typename PdbReader::RecordsType & records) {
	mapchainidint_.init(records);
	mapreskeyint_.resize(records.size(),
			NSPdstl::MapKeyInt<typename PdbReader::ResKeyType>());
	for (auto & c : records)
		mapreskeyint_.at(mapchainidint_.keynumber(c.first)).init(c.second);
}
std::string pdbfilename(const PDBModelID &mid) {
	std::string res;
	res.resize(mid.pdbid.size());
	std::transform(mid.pdbid.begin(), mid.pdbid.end(), res.begin(), ::tolower);
	res += ".pdb" + std::to_string(mid.biounit);
	return res;
}
std::string modelfilename(const PDBModelID &mid) {
	std::string res;
	res = mid.pdbid + "_model" + std::to_string(mid.model) + ".pdb"+std::to_string(mid.biounit);
	return res;
}
std::string NSPproteinrep::extractpdbmodel(const PDBModelID & mid) {
	std::string pdbpath = NSPdataio::downloadedpdbpath();
	std::string filename = pdbfilename(mid);
	std::string fullname = pdbpath + filename + ".gz";
	BFS::path fullpath(fullname);
	if (!BFS::exists(fullpath) || BFS::is_empty(fullpath)) {
		std::string cmd;
		cmd =
				"wget -q ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/"
						+ filename + ".gz";
//#ifndef NDEBUG
		std::cout << "downloading " << filename << ".gz from internet...\n";
//#endif
		FILE *pp = popen(cmd.c_str(), "w");
		if (pp == NULL) {
			std::cout << "downloading unsuccessful\n";
			return std::string();
		}
		std::fclose(pp);
		if (fullname != (filename + ".gz")) {
			std::cout <<"moving downloaded file \n";
			cmd = "mv " + filename + ".gz " + fullname;
			pp=popen(cmd.c_str(), "w");
			std::fclose(pp);
		}
	}
	std::string cmd;
	cmd = "gunzip -c " + fullname + " > " + filename;
	std::cout <<cmd<<std::endl;
	FILE *pp=popen(cmd.c_str(), "w");
	std::fclose(pp);
	return extractpdbmodel(filename, mid);
}
std::string NSPproteinrep::extractpdbmodel(const std::string & filename,
		const PDBModelID &mid) {
	std::string mfile = modelfilename(mid);
	std::ofstream ofs;
	ofs.open(mfile.c_str());
	NSPdataio::TextLines lines;
	lines.init(filename);
	bool inmodel = false;
	if (mid.model == 0)
		inmodel = true;
	int natoms = 0;
	for (auto &line : lines) {
		if (line.substr(0, 6) != "ATOM  " && line.substr(0, 6) != "HETATM"
				&& line.substr(0, 6) != "ENDMDL"
				&& line.substr(0, 5) != "MODEL") {
			ofs << line << std::endl;
			continue;
		}
		if (line.substr(0, 5) == "MODEL") {
			inmodel = false;
			int num = std::stoi(line.substr(6));
			if (num == mid.model)
				inmodel = true;
			continue;
		}
		if (line.substr(0, 6) == "ENDMDL") {
			inmodel = false;
			continue;
		}
		if (inmodel) {
			ofs << line << std::endl;
			++natoms;
		}
	}
	if (natoms > 0)
		return mfile;
	else
		return std::string();
}
