#define EXPORTING
#include "Analysis/FixationEstabilityAnalysis.h"
#include "Analysis/Numerics.h"
#include <Eigen/Core>
#include <Eigen/SVD>

namespace Analysis {
	FixationEstabilityAnalysis::FixationEstabilityAnalysis() {
		doneProcessing = false;

	}
	void FixationEstabilityAnalysis::LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& yL, const std::vector<double>& xR, const std::vector<double>& yR, bool processNow) {
		this->t = t;
		this->xL = xL;
		this->yL = yL;
		this->xR = xR;
		this->yR = yR;
		doneProcessing = false;
		if (processNow) {
			Process();
		}
	}
	void FixationEstabilityAnalysis::LoadVelocities(const std::vector<double>& vxL, const std::vector<double>& vyL, const std::vector<double>& vxR, const std::vector<double>& vyR) {

		this->vxL = vxL;
		this->vyL = vyL;
		this->vxR = vxR;
		this->vyR = vyR;

	}
	void FixationEstabilityAnalysis::Process() {
		if (xL.size() != vxL.size()) {
			firstDerivative(xL.data(), vxL.data(), xL.size());
			firstDerivative(yL.data(), vyL.data(), yL.size());
			firstDerivative(xR.data(), vxR.data(), xR.size());
			firstDerivative(yR.data(), vyR.data(), yR.size());
		}
		std::vector<double> vl(vxL.size()), vr(vxR.size()), al(vxL.size()), ar(vxR.size());
		polarCoordinates(vxL.data(), vyL.data(), vl.data(), al.data(), vxL.size());
		polarCoordinates(vxR.data(), vyR.data(), vr.data(), ar.data(), vxR.size());
		harmonicVelocities[0] = nanharmmean(vl.data(), vxL.size());
		harmonicVelocities[1] = nanharmmean(vr.data(), vxR.size());

		GetDispersionFromTraces(xL, yL, variances[0], mainComponentDirections[0]);
		GetDispersionFromTraces(xR, yR, variances[1], mainComponentDirections[1]);
		entropies[0] = entropy2D(xL.data(), yL.data(), xL.size(), -3, -3, 3, 3, 30, 30);
		entropies[1] = entropy2D(xR.data(), yR.data(), xR.size(), -3, -3, 3, 3, 30, 30);

		doneProcessing = true;
	}

	void FixationEstabilityAnalysis::GetEntropies(double& left, double& right) const {
		left = entropies[0];
		right = entropies[1];
	}
	void FixationEstabilityAnalysis::GetVariances(double& lh, double& lv, double& rh, double& rv) const {
		lh = variances[0][0];
		lv = variances[0][1];
		rh = variances[1][0];
		rv = variances[1][1];
	}
	void FixationEstabilityAnalysis::GetHarmonicVelocities(double& left, double& right) const {
		left = harmonicVelocities[0];
		right = harmonicVelocities[1];
	}
	void FixationEstabilityAnalysis::GetMainAxes(FixationEstabilityAnalysis::MainComponentDirectionsType& left, FixationEstabilityAnalysis::MainComponentDirectionsType& right) const {
		left = this->mainComponentDirections[0];
		right = this->mainComponentDirections[1];
	}

	void FixationEstabilityAnalysis::Clear() {
		vxL.resize(0);
		vyL.resize(0);
		vxR.resize(0);
		vyR.resize(0);
		xL.resize(0);
		yL.resize(0);
		xR.resize(0);
		yR.resize(0);
		doneProcessing = false;


	}

	void FixationEstabilityAnalysis::GetDispersionFromTraces(const std::vector<double>& x, const std::vector<double>& y, MainComponentsType& lambdas, MainComponentDirectionsType& directions) {
		using namespace Eigen;

		//pack data in Eigen matrix, removing nans
		unsigned long countNan = 0;
		std::vector<double> xN, yN;
		for (std::size_t k = 0; k < x.size(); k++) {
			if (isnan(x[k]) || isnan(y[k])) {
				countNan++;
			}
			else {
				xN.push_back(x[k]);
				yN.push_back(y[k]);
			}
		}
		Map<VectorXd> vx(xN.data(), Index(xN.size()));
		Map<VectorXd> vy(yN.data(), Index(yN.size()));

		MatrixXd m(xN.size(), 2);
		m.col(0) = vx;
		m.col(1) = vy;

		//calculate svd decomposition
		JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
		
		lambdas[0] = svd.singularValues()[0] * svd.singularValues()[0] / double(xN.size() - 1);
		lambdas[1] = svd.singularValues()[1] * svd.singularValues()[1] / double(xN.size() - 1);
		directions[0][0] = svd.matrixV().col(0)[0];
		directions[0][1] = svd.matrixV().col(0)[1];
		directions[1][0] = svd.matrixV().col(1)[0];
		directions[1][1] = svd.matrixV().col(1)[1];



	}

	void FixationEstabilityAnalysis::GetBCEA(double& lBCEA, double& rBCEA) const {

		double lvx, lvy, rvx, rvy;
		GetVariances(lvx, lvy, rvx, rvy);
		lBCEA = 2.291*M_PI*std::sqrt(lvx)*std::sqrt(lvy); //since the variances are calculated over the main components (not X and Y), Pearson correlation between componenets is zero and the equation simplifies to this
		rBCEA = 2.291*M_PI*std::sqrt(rvx)*std::sqrt(rvy);


	}
}