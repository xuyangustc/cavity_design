/*
 *aaconformer.h
 *
 *  Created on: 2017年6月23日
 *      Author: hyliu
 */

#ifndef PROTEINREP_AACONFORMER_H_
#define PROTEINREP_AACONFORMER_H_
#include "geometry/xyz.h"
#include "geometry/localframe.h"
#include "proteinrep/pdbreader.h"
#include <memory>
namespace NSPproteinrep {

struct AAConformer {
	static std::map<std::string, std::vector<std::string>> sidechainatoms;
	static std::vector<std::string> backboneatoms;
	static std::vector<std::string> mainchainatoms;
	std::string residuename { "" };
	std::map<std::string, NSPgeometry::XYZ> globalcrd;
	std::map<std::string, NSPgeometry::XYZ> localcrd;
	NSPgeometry::LocalFrame localframe;
	static NSPgeometry::LocalFrame calclocalframe(
			const std::map<std::string, NSPgeometry::XYZ> & gcrd);
	std::map<std::string, NSPgeometry::XYZ> & calclocalcrd();
	const std::map<std::string, NSPgeometry::XYZ> &getlocalcrd() const {
		return localcrd;
	}
	std::map<std::string, NSPgeometry::XYZ> & calcglobalcrd(
			const std::map<std::string, NSPgeometry::XYZ> &refglobal) {
		localframe = calclocalframe(refglobal);
		return calcglobalcrd();
	}
	std::map<std::string, NSPgeometry::XYZ> & calcglobalcrd();
	const std::map<std::string, NSPgeometry::XYZ> & getglobalcrd() const {
		return globalcrd;
	}
	bool mainchaincrdcomplete()const {
		for (auto &m : mainchainatoms) {
			if (globalcrd.find(m) == globalcrd.end())
				return false;
		}
		return true;
	}
	bool sidechaincrdcomplete() const {
		for (auto & a : sidechainatoms.at(residuename)) {
			if (globalcrd.find(a) == globalcrd.end())
				return false;
		}
		return true;
	}
	bool crdcomplete() const {
		return mainchaincrdcomplete() && sidechaincrdcomplete();
	}
	bool isnaturalaa() const {
		return (sidechainatoms.find(residuename) != sidechainatoms.end());
	}
	bool isaa() const {
		return isnaturalaa() || mainchaincrdcomplete();
		}
	void removeoxtorother(const std::string & name="OXT"){
		globalcrd.erase(name);
		localcrd.erase(name);
	}
	int removeatomsnotlisted();
	double phi(const AAConformer &cp) const {
		return NSPgeometry::torsion(cp.globalcrd.at("C"), globalcrd.at("N"),
				globalcrd.at("CA"), globalcrd.at("C"));
	}
	double psi(const AAConformer &cn) const {
		return NSPgeometry::torsion(globalcrd.at("N"), globalcrd.at("CA"),
				globalcrd.at("C"), cn.globalcrd.at("N"));
	}
	double omega_n(const AAConformer &cp) const {
		return NSPgeometry::torsion(cp.globalcrd.at("CA"), cp.globalcrd.at("C"),
				globalcrd.at("N"), globalcrd.at("CA"));
	}
	double omega_c(const AAConformer &cn) const {
		return NSPgeometry::torsion(globalcrd.at("CA"), globalcrd.at("C"),
				cn.globalcrd.at("N"), cn.globalcrd.at("CA"));
	}
	double dihedral(const std::string &a1, const std::string &a2,
			const std::string &a3, const std::string &a4) {
		if (!globalcrd.empty()) {
			return NSPgeometry::torsion(globalcrd.at(a1), globalcrd.at(a2),
					globalcrd.at(a3), globalcrd.at(a4));
		} else {
			return NSPgeometry::torsion(localcrd.at(a1), localcrd.at(a2),
					localcrd.at(a3), localcrd.at(a4));
		}
	}
	std::vector<PdbRecord> make_pdbrecords(char chainid, int resid,
			char insertionid, int atomid0);
};

double sidechainRMSD2_local(const AAConformer &c1, const AAConformer &c2);
double RMSD2_global(const AAConformer &c1, const AAConformer &c2);
double backboneRMSD2(const AAConformer &c1, const AAConformer &c2);

double distance(const AAConformer &c1, const AAConformer &c2,
		const std::string & atom1 = "CA", const std::string &atom2 = "CA");
bool isbackbone(const std::string &atomname);
bool ismainchain(const std::string &atomname);
AAConformer make_aaconformer(const std::vector<PdbRecord> & residuerecords,
		std::vector<AAConformer> *altconformers = nullptr);
struct AAConformersInModel {
	PDBModelID pdbmodelid;
	std::shared_ptr<MapPdbKeyInt> mappdbkeyint { nullptr };
	std::vector<std::vector<AAConformer>> conformers;
	std::vector<std::vector<std::vector<AAConformer>>>altconformers;
	std::string getsequence(int chainnumber) const;
	void readpdbfile(const std::string &filename);
	const std::vector<AAConformer> & conformersinchain(char chainid=' ') const
	{
		if( chainid != ' ' && mappdbkeyint) {
			return conformers.at(mappdbkeyint->chainNumber(chainid));
		}
		else return conformers[0];
	};
};
}

#endif /* PROTEINREP_AACONFORMER_H_ */
