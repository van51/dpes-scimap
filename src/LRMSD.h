/*
	(c) 2014 Evangelos Anagnostopoulos
	This code is licensed under the GNU General Public License, version 2 (GPL-2.0)
	(see LICENSE for details)
 */

#ifndef __LRMSD_H__
#define __LRMSD_H__

#include <iostream>
#include <cmath>

#include "Configuration.h"

typedef Eigen::Matrix<double,1,Eigen::Dynamic> RowVector;
class LRMSD {
public:
	
	static double computeRMSD(Configuration* conf1, Configuration* conf2);
};

#endif
