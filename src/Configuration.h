/*
	(c) 2014 Evangelos Anagnostopoulos
	This code is licensed under the GNU General Public License, version 2 (GPL-2.0)
	(see LICENSE for details)
 */

#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

#include <vector>
#include <eigen3/Eigen/Dense>

typedef Eigen::MatrixXd matrix;

class Configuration {
public:	
	Configuration(int numberOfPoints) : points(numberOfPoints, 3) {
		id = counter++;
	}

	double& operator() (int i, int j) {
		return points(i, j);
	}

	void set(int i, int j, double value) {
		points(i, j) = value;
	}

	double get(int i, int j) {
		return points(i, j);
	}

	matrix& getMatrix() {
		return points;
	}

	int getNumberOfAtoms() {
		return points.rows();
	}

	int getID() {
		return id;
	}

	static int counter;

private:
	matrix points;
	int id;
};

typedef std::vector<Configuration* > Configurations;

#endif // __CONFIGURATION_H__
