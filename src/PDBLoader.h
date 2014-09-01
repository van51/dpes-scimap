/*
	(c) 2014 Evangelos Anagnostopoulos
	This code is licensed under the GNU General Public License, version 2 (GPL-2.0)
	(see LICENSE for details)
 */

#ifndef __PDBLOADER_H__
#define __PDBLOADER_H__

#include <vector>
#include "Configuration.h"
#include <string>

class PDBLoader {
public:	

	PDBLoader(std::string file, int numberOfConfigurations, int numberOfAtoms) : filename(file),
		numOfConfs(numberOfConfigurations), numOfAtoms(numberOfAtoms) {} 

	void getConfigurationsFromFile(Configurations& confs);
private:
	std::string filename;
	int numOfConfs;
	int numOfAtoms;
};

#endif
