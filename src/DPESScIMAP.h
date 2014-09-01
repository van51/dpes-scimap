/*
	(c) 2014 Evangelos Anagnostopoulos
	This code is licensed under the GNU General Public License, version 2 (GPL-2.0)
	(see LICENSE for details)
 */

#ifndef __DPESSCIMAP_H__
#define __DPESSCIMAP_H__

#include "Configuration.h"
#include "LRMSD.h"
#include "PDBLoader.h"
#include <vector>
#include <utility>
#include <string>

typedef std::pair<int, double> nnScore;

class DPESScIMAP {
public:	
	DPESScIMAP(std::string filename, int numOfConfs, int numOfAtoms, int numOfPivots,
			int numOfEuclidNN, int numOfRMSDNN, int numOfLandmarks);

	~DPESScIMAP() {
		for (size_t i=0; i<confs.size(); i++) {
			delete confs[i];
		}
	}

	void run();

	void selectPivots(Configurations& pivots);

	double** projectConfigurations(Configurations& pivots);

	void filterNearestNeighbors(int** euclideanNNs);
private:
	int numberOfConfigurations;
	int numberOfAtoms;
	int numberOfPivots;
	int numberOfEuclideanNN;
	int numberOfRMSDNN;
	int numberOfLandmarks;
	Configurations confs;
};
#endif
