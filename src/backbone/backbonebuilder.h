/*
 * backbonebuilder.h
 *
 *  Created on: 2018年1月9日
 *      Author: hyliu
 */

#ifndef BACKBONE_BACKBONEBUILDER_H_
#define BACKBONE_BACKBONEBUILDER_H_
#include "backbone/backbonesite.h"
#include "fullsite/fullsite.h"
#include <set>
#include <memory>
namespace NSPproteinrep{

class BackBoneBuilder {
public:
	/** build a segment of strand with given starting(or ending) position and direction.
	 * @para length:  number of residues forming the strand
	 * @para r0: position of the CA atom of the head residue
	 * @para direction:  Direction from the head CA position to the tail CA position.
	 * @para forward: if true, the head residue is the first residue. Otherwise it is the last residue.
	 */
	static std::vector<BackBoneSite> buildstrandat(int length,NSPgeometry::XYZ r0, NSPgeometry::XYZ direction,bool forward);
	/** build a segment of helix with given starting(or ending) position and direction.
		 * @para length:  number of residues forming the helix
		 * @para r0: position of the CA atom of the head residue
		 * @para direction:  Direction from the head CA position to the tail CA position.
		 * @para forward: if true, the head residue is the first residue. Otherwise it is the last residue.
		 */
	static std::vector<BackBoneSite> buildhelixat(int length, NSPgeometry::XYZ r0,
			NSPgeometry::XYZ direction_,bool forward);
	/**build a segment of backbone, starting from the N-terminal.
	 * @para length: number of residues in the segment
	 * @para nflankingsite: a backbone site right before the segment to build.
	 * @para helixregions: helix regions in the built segment. Each region is specified a pair of position and length.
	 * @para strandregions: strand regions in the build segment. Each region is specified as a pair of position and length.
	 * @para cissites: positions with cis-peptide bonds.
	 */
	static std::vector<BackBoneSite> buildforwardbackbone(int length,
			const BackBoneSite & nflankingsite, const std::vector<std::pair<int,int>> & helixregions,
			const std::vector<std::pair<int,int>> & strandregions, const std::set<int> & cissites);
	/**build a segment of backbone, starting from the C-terminal. The positions are still indexed from the N to the C terminal
	 * @para length: number of residues in the segment
	 * @para nflankingsite: a backbone site right before the segment to build.
	 * @para helixregions: helix regions in the built segment. Each region is specified a pair of position and length.
	 * @para strandregions: strand regions in the build segment. Each region is specified as a pair of position and length.
	 * @para cissites: positions with cis-peptide bonds.
	 */
	static std::vector<BackBoneSite> buildbackwardbackbone(int length,
			const BackBoneSite & cflankingsite, const std::vector<std::pair<int,int>> & helixregions,
			const std::vector<std::pair<int,int>> & strandregions, const std::set<int> & cissites);
	/** build a set of segments that can close a gap between a N terminal and a C terminalbackbone site.
	 *@para length:length of each segments.
	 *@para nflankingsite: n-terminal site. Not included in the segment to return.
	 *@para cflankingsite: c-terminal site. Not included in the segment to return.
	 *@para helixregions: helix regions in the built segment.
	 *@para strand regions: strand regions in the buildsegment
	 *@para cissites: positions of cis peptide bonds.
	 *@para n12result: will be used when length==1 or length==2
	 */
	static  std::vector<std::shared_ptr<std::vector<BackBoneSite>>> buildlinkers(int length,const BackBoneSite &nflankingsite,
			const BackBoneSite &clankingsite,
			const std::vector<std::pair<int,int>> & helixregions,
			const std::vector<std::pair<int,int>> & strandregions,
			const std::set<int> & cissites, int n12result=1);
	/**move a segment to a given head CA position, and rotate it so its head-to-tail CA-CA vector
	 * points along direction.
	 * If foward, head is the first residue and tail is the last residue.
	 * Otherwise the first residue is head and the last residue is tail.
	 */
	static void movechainto(NSPgeometry::XYZ r0,NSPgeometry::XYZ direction,bool forward,
			std::vector<BackBoneSite> &chain);



	static  std::vector<std::shared_ptr<std::vector<BackBoneSite>>> buildlinker1(
			const BackBoneSite &nflankingsite, const BackBoneSite &clankingsite, int nr, bool cis0);
	static  std::vector<std::shared_ptr<std::vector<BackBoneSite>>> buildlinker2(
			const BackBoneSite &nflankingsite, const BackBoneSite &clankingsite, int nr, bool cis0, bool cis1);




	static NSPgeometry::Rotation rot(NSPgeometry::XYZ fixed1, NSPgeometry::XYZ fixed2,
			NSPgeometry::XYZ fixed3, NSPgeometry::XYZ moved3);
	static NSPgeometry::Rotation transform(NSPgeometry::XYZ fixed1, NSPgeometry::XYZ fixed2,
			NSPgeometry::XYZ moved1, NSPgeometry::XYZ moved2, NSPgeometry::XYZ &translated);
	static std::vector<BackBoneSite> chainrotation(const BackBoneSite &bs, const std::vector<BackBoneSite> &chain);
	static std::vector<BackBoneSite> getbackbone(const std::vector<std::pair<std::pair<int,int>,char>> &label,
			const std::set<int> &cissite, const std::vector<std::vector<BackBoneSite>>&fixed, const BackBoneSite &bsn);
	static std::vector<std::shared_ptr<std::vector<BackBoneSite>>> buildlinkers(
			const BackBoneSite & nflankingsite, const BackBoneSite & cflankingsite,
			const std::vector<std::pair<std::pair<int,int>,char>> &label,
			const std::vector<std::vector<BackBoneSite>> & fixedconf, const std::set<int> & cissites);







	struct Region {
		Region(const std::string &label, const std::set<int> &cis);
		int length;
		std::vector<std::pair<int,int>> helixregion;
		std::vector<std::pair<int,int>> strandregion;
		std::set<int> cissites;
	};
	static std::vector<std::shared_ptr<std::vector<BackBoneSite>>> buildlinkers(
			const BackBoneSite & nflankingsite, const BackBoneSite & cflankingsite, const std::string & label,
			const std::vector<std::vector<BackBoneSite>> & fixedconf, const std::set<int> & cissites);
	static std::vector<FullSite> addsidechain(
			const std::vector<BackBoneSite>&bsch, const std::vector<FullSite>&orifs);
	static std::vector<FullSite> addsidechain(
			const std::vector<BackBoneSite>&bsch, const std::map<int,std::vector<FullSite>>&orifs);
private:
};


}



#endif /* BACKBONE_BACKBONEBUILDER_H_ */
