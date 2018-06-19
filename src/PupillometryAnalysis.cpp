#define EXPORTING
#include "Analysis/PupillometryAnalysis.h"
#include <memory>
#include <algorithm>
#include <numeric>
#include "Analysis\Algorithms.h"
#include "ceres/ceres.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif
#define COS225 0.924

namespace Analysis {

	
	const unsigned int PupillometryAnalysis::logitNComponents = 3;
		
		void PupillometryAnalysis::LoadTrials(const TaskDataType& data, const TaskEventsType& events, bool processNow) {
			this->data = data;
			this->events = events;

			if (processNow) {
				ProcessAll();
			}

			
		}

		

		bool PupillometryAnalysis::ProcessAll() {
			using namespace ceres;
			
			double tFlash = events[0].inTrialTS*1000;
			double tFlashEnd = events[1].inTrialTS*1000;
			Problem fittingLeft, fittingRight;

			std::vector<double> aLeft = { -2.0, 5.0, 0.0 };
			std::vector<double> aRight = { -2.0, 5.0, 0.0 };
			double bLeft = 5.0;
			double bRight = 5.0;
			std::vector<double> TLeft = { 0.5, 1.5, 3.0 };
			std::vector<double> TRight = { 0.5, 1.5,3.0 };
			std::vector<double> DLeft = { 0.1, 1.0, 0.1 };
			std::vector<double> DRight = { 0.1, 1.0, 0.1 };
			timeAroundFlash.resize(data.at("t").size());
			std::transform(data.at("t").begin(), data.at("t").end(), timeAroundFlash.begin(), [tFlash](double x) -> double {return x - tFlash; });

			for (std::size_t k = 1000; k < data.at("t").size(); k++) {
				if (!std::isnan(data.at("lp")[k])) {
					CostFunction* cost_function_left =
						new AutoDiffCostFunction<typename InverseLogitModelResidual, 1, logitNComponents, logitNComponents, logitNComponents, 1>
						(new InverseLogitModelResidual(logitNComponents, timeAroundFlash[k] / 1000.0, data.at("lp")[k]));
					//fittingLeft.AddResidualBlock(cost_function_left,new CauchyLoss(0.5), (double*)aLeft.data(), (double*)TLeft.data(), (double*)DLeft.data(), &bLeft);
					fittingLeft.AddResidualBlock(cost_function_left, NULL, (double*)aLeft.data(), (double*)TLeft.data(), (double*)DLeft.data(), &bLeft);
					for (std::size_t i = 0; i < logitNComponents; i++) {
						fittingLeft.SetParameterLowerBound(TLeft.data(), i, 0.0);
						fittingLeft.SetParameterLowerBound(DLeft.data(), i, 0.0);
					}

				}
				if (!std::isnan(data.at("rp")[k])) {
					CostFunction* cost_function_right =
						new AutoDiffCostFunction<typename InverseLogitModelResidual, 1, logitNComponents, logitNComponents, logitNComponents, 1>
						(new InverseLogitModelResidual(logitNComponents, timeAroundFlash[k] / 1000.0, data.at("rp")[k]));
					//fittingRight.AddResidualBlock(cost_function_right, new CauchyLoss(0.5), (double*)aRight.data(), (double*)TRight.data(), (double*)DRight.data(), &bRight);
					fittingRight.AddResidualBlock(cost_function_right, NULL, (double*)aRight.data(), (double*)TRight.data(), (double*)DRight.data(), &bRight);
					for (std::size_t i = 0; i < logitNComponents; i++) {
						fittingRight.SetParameterLowerBound(TRight.data(), i, 0.0);
						fittingRight.SetParameterLowerBound(DRight.data(), i, 0.0);
					}

				}
			}
			Solver::Options solverOptions;
			solverOptions.max_num_iterations = 300;
			solverOptions.linear_solver_type = ceres::DENSE_QR;
			solverOptions.minimizer_progress_to_stdout = false;

			Solver::Summary summaryLeft, summaryRight;
			Solve(solverOptions, &fittingLeft, &summaryLeft);
			//std::cout << summaryLeft.BriefReport() << std::endl;
			Solve(solverOptions, &fittingRight, &summaryRight);
			//std::cout << summaryRight.BriefReport() << std::endl;

			
			std::vector<double> yModelLeft, yModelRight;

			InverseLogitModelResidual::CalculateModelValues(aLeft, TLeft, DLeft, bLeft, timeAroundFlash, yModelLeft);
			InverseLogitModelResidual::CalculateModelValues(aRight, TRight, DRight, bRight, timeAroundFlash, yModelRight);

			Problem leftMinimizer, rightMinimizer;
			CostFunction* left_cost = new AutoDiffCostFunction<typename InverseLogitModelMinimizer, 1, 1>
				(new InverseLogitModelMinimizer(aLeft, TLeft, DLeft, bLeft));
			CostFunction* right_cost = new AutoDiffCostFunction<typename InverseLogitModelMinimizer, 1, 1>
				(new InverseLogitModelMinimizer(aRight, TRight, DRight, bRight));

			double minLeft, minRight;

			timeToMaxContractionLeft = 1.0;
			timeToMaxContractionRight = 1.0;
			double troughLeft;
			leftMinimizer.AddResidualBlock(left_cost, NULL, &timeToMaxContractionLeft);
			rightMinimizer.AddResidualBlock(right_cost, NULL, &timeToMaxContractionRight);
			Solve(solverOptions, &leftMinimizer, &summaryLeft);
			Solve(solverOptions, &rightMinimizer, &summaryRight);

			InverseLogitModelResidual::CalculateModelValues(aLeft, TLeft, DLeft, bLeft, timeToMaxContractionLeft, minLeft);
			InverseLogitModelResidual::CalculateModelValues(aRight, TRight, DRight, bRight, timeToMaxContractionRight, minRight);

			maxContractionLeft = std::abs(bLeft - minLeft);
			maxContractionRight = std::abs(bRight - minRight);

			auto minPeakItL = std::find_if(timeAroundFlash.begin(), timeAroundFlash.end(), [this](const double& t)-> bool {return t >= this->timeToMaxContractionLeft; });
			auto minPeakItR = std::find_if(timeAroundFlash.begin(), timeAroundFlash.end(), [this](const double& t)-> bool {return t >= this->timeToMaxContractionRight; });

			for (auto itL = minPeakItL; itL != timeAroundFlash.end(); ++itL) {
				std::size_t index = itL - timeAroundFlash.begin();
				if (data.at("glx")[index] >= maxContractionLeft + (bLeft - maxContractionLeft) / 2.0) {
					timeToRecoveryLeft = *itL - timeToMaxContractionLeft;
					break;
				}
			}

			for (auto itR = minPeakItR; itR != timeAroundFlash.end(); ++itR) {
				std::size_t index = itR - timeAroundFlash.begin();
				if (data.at("grx")[index] >= maxContractionRight + (bRight - maxContractionRight) / 2.0) {
					timeToRecoveryRight = *itR - timeToMaxContractionRight;
					break;
				}
			}

			doneProcessing = true;
			
			return true;
		}

		bool PupillometryAnalysis::Clear() {
			return true;
		}

		PupillometryAnalysis::PupillometryAnalysis() {
			
			doneProcessing = false;

		}


		void PupillometryAnalysis::GetMaxConstriction(double& left, double& right) const {
			left = maxContractionLeft;
			right = maxContractionRight;
		}
		void PupillometryAnalysis::GetTimeToMaxConstriction(double& left, double& right) const {
			left = timeToMaxContractionLeft;
			right = timeToMaxContractionRight;
		}

		void PupillometryAnalysis::GetTimeToRecovery(double& left, double& right) const {
			left = timeToRecoveryLeft;
			right = timeToRecoveryRight;
		}

		void PupillometryAnalysis::GetTimeAroundFlash(std::vector<double>& t) const {
			t = timeAroundFlash;
		}

		Json::Value PupillometryAnalysis::ToJSONObject() {
			Json::Value root;

			root["LeftTimeToMaxContraction"] = timeToMaxContractionLeft;
			root["RightTimeToMaxContraction"] = timeToMaxContractionRight;
			root["LeftContraction"] = maxContractionLeft;
			root["RightContraction"] = maxContractionRight;
			root["LeftTimeToHalfRecovery"] = timeToRecoveryLeft;
			root["RightTimeToHalfRecovery"] = timeToRecoveryRight;
			return root;
		}
}