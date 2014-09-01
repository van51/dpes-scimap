/*
	(c) 2014 Evangelos Anagnostopoulos
	This code is licensed under the GNU General Public License, version 2 (GPL-2.0)
	(see LICENSE for details)
 */

#include <fstream>
#include <sstream>
#include "PDBLoader.h"
#include <iostream>

void PDBLoader::getConfigurationsFromFile(Configurations& confs) {
	std::ifstream file(filename.c_str());

	std::string line;
	std::getline(file, line);
	int numberOfConfsRead = 0;
	while (numberOfConfsRead<numOfConfs) {
		Configuration* configuration = new Configuration(numOfAtoms);
		int currentAtom = 0;
		while (std::getline(file, line)) {
			std::istringstream parsePoints(line.substr(27, std::string::npos));
			double x, y, z;
			parsePoints >> x >> y >> z;
			configuration->set(currentAtom, 0, x);
			configuration->set(currentAtom, 1, y);
			configuration->set(currentAtom, 2, z);
			if (++currentAtom==numOfAtoms)
				break;
		}
		confs.push_back(configuration);
		std::getline(file, line); // read END
		numberOfConfsRead++;
	}
}

