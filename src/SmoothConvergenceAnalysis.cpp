#define EXPORTING
#include "Analysis/SmoothConvergenceAnalysis.h"
#include <memory>
#include <algorithm>
#include <numeric>
#include "Analysis\Algorithms.h"
#include "ceres/ceres.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

namespace Analysis {

	
	SmoothConvergenceAnalysis::SmoothConvergenceAnalysis() { doneProcessing = false; };
	void SmoothConvergenceAnalysis::LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& xR, bool processNow) {
		this->t = t;
		this->xL = xL;
		this->xR = xR;
		if (processNow)
			Process();
	}
	void SmoothConvergenceAnalysis::SetVelocity(double v) {
		this->v = v;

	}
	void SmoothConvergenceAnalysis::SetStartingPosition(double d0) {
		this->d0 = d0;
	}
	double SmoothConvergenceAnalysis::GetApproachingConvergenceDistance() const {
		return 0.0;
	}
	unsigned int SmoothConvergenceAnalysis::GetApproachingMisses() const {
		return 0;
	}
	bool SmoothConvergenceAnalysis::Process() {

		d.resize(t.size());
		std::size_t n = 1;
		d[0] = d0;
		unsigned int direction = v > 0 ? 1 : 0;
		std::vector<std::vector<double> > distances(2);
		std::vector<std::vector<double> >vergence(2);
		double velocity = v;
		std::vector<int> periodDirection;
		std::vector<std::size_t> periodStartSample;
		
		for (std::size_t k = 0; k < tDirectionChanges.size(); k++) {
			periodStartSample.push_back(n);
			periodDirection.push_back(direction);
			while (t[n] < tDirectionChanges[k] & n < d.size() - 1) {
				d[n] = d[n - 1] + velocity*0.004;
				n++;

			}
			velocity *= -1.0;
			direction = (direction + 1) % 2;
		}
		verg.resize(xL.size());
		std::transform(xL.begin(), xL.end(), xR.begin(), verg.begin(), std::minus<double>());

		for (std::size_t k = 0; k < periodStartSample.size()-1; k++) { //remove the first and last periods since they are incomplete
			auto startV = verg.begin() + periodStartSample[k];
			auto endV = verg.begin() + periodStartSample[k + 1] - 1;
			std::vector<double> trialVergence(startV, endV);
			auto startD = d.begin() + periodStartSample[k];
			auto endD = d.begin() + periodStartSample[k + 1] - 1;
			std::vector<double> dTrial(startD, endD);
			double baseline = 0.0;
			double a = 2.0;
			double b = 3.0;
			double c = 15.0;
			double G = 10.0;
			trialFitting(trialVergence, dTrial, a, b, c, baseline, G);
			if (periodDirection[k] == 0) {
				approachingDistances.push_back(dTrial);
				approachingVergences.push_back(trialVergence);
				approachingNPC.push_back(c);
				approachingSlopes.push_back(b / (2.0*a));
			}
			else {
				distancingDistances.push_back(dTrial);
				distancingVergences.push_back(trialVergence);
				distancingNPC.push_back(c);
				distancingSlopes.push_back(b / (2.0*a));
			}
		}

		doneProcessing = true;
		return true;
	}

	void SmoothConvergenceAnalysis::trialFitting(const std::vector<double>& v, const std::vector<double>& d, double& a, double& b, double& c, double& baseline, double& G) {
		using namespace ceres;

		Problem fitting;

		for (std::size_t k = 1; k < v.size(); k++) {
			if (!std::isnan(v[k])) {
				CostFunction* cost_function_approaching =
					new AutoDiffCostFunction<typename GeneralizedBellFunctionResidual, 1, 1, 1, 1, 1, 1>
					(new GeneralizedBellFunctionResidual(d[k], v[k]));
				fitting.AddResidualBlock(cost_function_approaching, new CauchyLoss(1.0), &a, &b, &c, &baseline, &G);
				//fittingApproaching.AddResidualBlock(cost_function_approaching, NULL, &a, &b, &c, &baseline, &G);
				fitting.SetParameterLowerBound(&a, 0, 0.0);
				fitting.SetParameterLowerBound(&b, 0, 0.0);
				fitting.SetParameterLowerBound(&c, 0, 0.0);
				fitting.SetParameterLowerBound(&G, 0, 0.0);
				fitting.SetParameterUpperBound(&G, 0, 50.0);


			}
		}
		Solver::Options solverOptions;
		solverOptions.max_num_iterations = 300;
		solverOptions.linear_solver_type = ceres::DENSE_QR;
		solverOptions.minimizer_progress_to_stdout = false;

		Solver::Summary summary;
		Solve(solverOptions, &fitting, &summary);
		
	}

	void SmoothConvergenceAnalysis::GetApproachingWaveforms(std::vector<std::vector<double> >& d, std::vector<std::vector<double> >& v) const {
		d = approachingDistances;
		v = approachingVergences;
	}

	bool SmoothConvergenceAnalysis::Clear() {
		doneProcessing = false;
		return true;
	}

	void SmoothConvergenceAnalysis::SetDirectionChangeTimestamps(std::vector<double> t) {
		this->tDirectionChanges = t;
	}

	double SmoothConvergenceAnalysis::GetMedianApproachingNPC() const {
		return nanmedian(approachingNPC.data(), approachingNPC.size());
	}

	double SmoothConvergenceAnalysis::GetMedianDistancingNPC() const {
		return nanmedian(distancingNPC.data(), distancingNPC.size());
	}

	Json::Value SmoothConvergenceAnalysis::ToJSONObject() {
		if (!doneProcessing) return "";
		Json::Value root;   // starts as "null"; will contain the root value after parsing


		root["Approaching"]["NPC"] = GetMedianApproachingNPC();
		root["Distancing"]["NPC"] = GetMedianDistancingNPC();


		return root;
	}
}