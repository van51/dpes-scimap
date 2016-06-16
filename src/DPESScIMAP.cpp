/*
	(c) 2014 Evangelos Anagnostopoulos
	This code is licensed under the GNU General Public License, version 2 (GPL-2.0)
	(see LICENSE for details)
 */

#define TAPKEE_WITH_ARPACK
#include <tapkee/tapkee.hpp>
#include <ANN/ANN.h>
#include "DPESScIMAP.h" 
#include <set>
#include <climits>
#include <fstream>
#include <algorithm>

using namespace std;

bool scoreComparison(nnScore s1, nnScore s2) {
	return s1.second < s2.second;
}

struct MyDistanceCallback {
	double distance(Configuration* l, Configuration* r) { 
		return LRMSD::computeRMSD(l, r);
	} 
} callbackDistance; 


DPESScIMAP::DPESScIMAP(string filename, int numOfConfs, int numOfAtoms, int numOfPivots,
	int numOfEuclidNN, int numOfRMSDNN, int numOfLandmarks) : 
		numberOfConfigurations(numOfConfs), numberOfAtoms(numOfAtoms), numberOfPivots(numOfPivots),
		numberOfEuclideanNN(numOfEuclidNN), numberOfRMSDNN(numOfRMSDNN), numberOfLandmarks(numOfLandmarks) {

	PDBLoader pdbLoader(filename, numberOfConfigurations, numberOfAtoms);
	confs.reserve(numberOfConfigurations);
	pdbLoader.getConfigurationsFromFile(confs);
}

void DPESScIMAP::run() {
	Configurations pivots;
	pivots.reserve(numberOfPivots);
	selectPivots(pivots);
	double** euclideanProjection = projectConfigurations(pivots);
	pivots.clear();
	/** NN panw ston projectedConfigurations pinaka */
	int** nns = new int*[numberOfConfigurations];
	for (int i=0; i<numberOfConfigurations; i++)
		nns[i] = new int[numberOfEuclideanNN];

	ANNkd_tree* kdTree = new ANNkd_tree(euclideanProjection, numberOfConfigurations, numberOfPivots);
	ANNdistArray dists = new ANNdist[numberOfEuclideanNN];

	for(int i = 0; i < numberOfConfigurations; ++i) {
		kdTree->annkSearch(euclideanProjection[i], numberOfEuclideanNN, nns[i], dists, 0);
	}
	delete kdTree;
	annClose();
	delete euclideanProjection;

	filterNearestNeighbors(nns);
	for (int i=0; i<numberOfConfigurations; i++) {
		tapkee::linkNeighbors.push_back(vector<int>());
		for (int j=0; j<numberOfRMSDNN; j++) {
			tapkee::linkNeighbors[i].push_back(nns[i][j]);
		}
		delete[] nns[i];
	}
	delete[] nns;
	tapkee::TapkeeOutput to = tapkee::initialize().withParameters((tapkee::method=tapkee::Isomap, tapkee::target_dimension=2,
				tapkee::num_neighbors=numberOfLandmarks, tapkee::eigen_method=tapkee::Arpack)).withDistance(callbackDistance).embedRange(confs.begin(), confs.end());

	ofstream output("out.pdo", std::ofstream::out);
	for (int i=0; i<numberOfConfigurations; i++) {
		output << (i+1) << "\t" << to.embedding(i,0) << endl;
	}
	output.close();
}
void DPESScIMAP::filterNearestNeighbors(int** euclideanNNs) {
	vector<nnScore> scores;
	scores.resize(numberOfEuclideanNN);
	for (int i=0; i<numberOfConfigurations; i++) {
		for (int j=0; j<numberOfEuclideanNN; j++) {
			scores[j] = nnScore(euclideanNNs[i][j], LRMSD::computeRMSD(confs[i], confs[euclideanNNs[i][j]]));
		}
		sort(scores.begin(), scores.end(), scoreComparison);
		for (int j=0; j<numberOfRMSDNN; j++) {
			euclideanNNs[i][j] = scores[j].first;
		}
		//euclideanNNs[i].resize(numberOfRMSDNN);
	}	
}

void DPESScIMAP::selectPivots(Configurations& pivots) {
	set<int> pivotIndices;

	int initialPivot = rand() % numberOfConfigurations;
	pivotIndices.insert(initialPivot);

	while (static_cast<int>(pivots.size())<numberOfPivots) {
		double maxRMSD = -1.0;
		int maxRMSDpivot = -1;
		for (int i=0; i<numberOfConfigurations; i++) {
			if (pivotIndices.find(i)!=pivotIndices.end()) {
				double minRMSD = INFINITY;
				for (set<int>::iterator it=pivotIndices.begin(); it!=pivotIndices.end(); it++) {
					double rmsd = LRMSD::computeRMSD(confs[i], confs[(*it)]);
					if (rmsd<minRMSD) {
						minRMSD = rmsd;
					}
				}
				if (minRMSD>maxRMSD) {
					maxRMSD = minRMSD;
					maxRMSDpivot = i;
				}
			}
		}
		pivotIndices.insert(maxRMSDpivot);
		pivots.push_back(confs[maxRMSDpivot]);
	}

}

double** DPESScIMAP::projectConfigurations(Configurations& pivots) {
	//matrix* projectedConfs = new matrix(numberOfConfigurations, numberOfPivots);
	double** projectedConfs = new double*[numberOfConfigurations];
	for (int i=0; i<numberOfConfigurations; i++)
		projectedConfs[i] = new double[numberOfPivots];

	for (int i=0; i<numberOfConfigurations; i++) {
		for (int j=0; j<numberOfPivots; j++) {
			//((*projectedConfs)(i,j)) = LRMSD::computeRMSD(confs[i], pivots[j]);
			projectedConfs[i][j] = LRMSD::computeRMSD(confs[i], pivots[j]);
		}
	}
	return projectedConfs;
}

