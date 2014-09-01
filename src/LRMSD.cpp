/*
	(c) 2014 Evangelos Anagnostopoulos
	This code is licensed under the GNU General Public License, version 2 (GPL-2.0)
	(see LICENSE for details)
 */

#include "LRMSD.h"
#include "eigen3/Eigen/SVD"
#include "eigen3/Eigen/LU"

double LRMSD::computeRMSD(Configuration* conf1, Configuration* conf2) {
	RowVector conf1Mean = conf1->getMatrix().colwise().mean();
	RowVector conf2Mean = conf2->getMatrix().colwise().mean();

	matrix centeredX = conf1->getMatrix().rowwise() - conf1Mean;
	matrix centeredY = conf2->getMatrix().rowwise() - conf2Mean;
	Eigen::JacobiSVD<matrix> svd(centeredY.transpose() * centeredX, Eigen::ComputeThinU | Eigen::ComputeThinV);
	matrix Q = svd.matrixU() * svd.matrixV().transpose();
	if (Q.determinant() < 0) 
		Q(2,2) = -Q(2,2);

	centeredY = centeredY * Q;

	return 1 / sqrt(conf1->getNumberOfAtoms()) * (centeredX - centeredY).norm();
}
