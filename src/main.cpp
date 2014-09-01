/*
	(c) 2014 Evangelos Anagnostopoulos
	This code is licensed under the GNU General Public License, version 2 (GPL-2.0)
	(see LICENSE for details)
 */

#include <iostream>
#include <string>
#include <fstream>
#include "DPESScIMAP.h"

using namespace std;

int main(int argc, char* argv[]) {
	iostream::sync_with_stdio(false);

	int numberOfPivots = 80;
	int numberOfConfigurations = 164;
	int numberOfAtoms = 369;
	int k = 163;
	int l = 163;
	string filename = "../simulation_data/abeta_trajectory";
	string outFilename = "out.pdo";
	if (argc > 1) {
		ifstream input(argv[1]);
		input >> numberOfPivots;
		input >> numberOfConfigurations;
		input >> numberOfAtoms;
		input >> k;
		input >> l;
		input >> filename;
		input >> outFilename;
		input.close();
	}
	DPESScIMAP alg(filename, numberOfConfigurations, numberOfAtoms, numberOfPivots, k, l, l);
	alg.run();
}
