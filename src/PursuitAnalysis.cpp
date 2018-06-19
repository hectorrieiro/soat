#define EXPORTING
#include "Analysis/PursuitAnalysis.h"
#include <memory>
#include <algorithm>
#include <numeric>
#include "Analysis\Algorithms.h"
#include "ceres/ceres.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

namespace Analysis {

	
		PursuitAnalysis::PursuitAnalysis() {
			doneProcessing = false;
			f = 0.2;
			A = 20.0;
		}
		void PursuitAnalysis::LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& yL, const std::vector<double>& xR, const std::vector<double>& yR, bool processNow) {
			doneProcessing = false;
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

		void PursuitAnalysis::SetPursuitAmplitude(double A) {
			this->A = A;
		}

		void PursuitAnalysis::SetPursuitFrequency(double f) {
			this->f = f;
		}

		bool PursuitAnalysis::Process() {

			using namespace ceres;
			Problem fittingLeft, fittingRight;
			//skip the first second or so
			auto it = std::upper_bound(t.begin(), t.end(), 1000.0);
			std::size_t initIdx = it - t.begin();
			double Aleft = 20.0;
			double Aright = 20.0;
			double fleft = this->f;
			double fright = this->f;
			double phaseleft = 0.0;
			double phaseright = 0.0;
			double cleft = 0.0;
			double cright = 0.0;
			unsigned int nLeft = 0;
			unsigned int nRight = 0;
			for (std::size_t k = initIdx; k < t.size(); k++) {
				if (!std::isnan(xL[k])) {
					CostFunction* cost_function_left = new AutoDiffCostFunction<typename SinusoidalResidual, 1, 1, 1, 1, 1>(new SinusoidalResidual(t[k] / 1000.0, xL[k]));
					fittingLeft.AddResidualBlock(cost_function_left, NULL, &Aleft, &fleft, &phaseleft, &cleft);
					cleft += xL[k];
					nLeft++;
				}
				if (!std::isnan(xR[k])) {
					CostFunction* cost_function_right = new AutoDiffCostFunction<typename SinusoidalResidual, 1, 1, 1, 1, 1>(new SinusoidalResidual(t[k] / 1000.0, xR[k]));
					fittingRight.AddResidualBlock(cost_function_right, NULL, &Aright, &fright, &phaseright, &cright);
					cright += xR[k];
					nRight++;
				}
			}
			cleft = cleft / double(nLeft);
			cright = cright / double(nRight);
			Solver::Options solverOptions;
			solverOptions.max_num_iterations = 75;
			solverOptions.linear_solver_type = ceres::DENSE_QR;
			solverOptions.minimizer_progress_to_stdout = true;
			Solver::Summary summaryLeft, summaryRight;
			Solve(solverOptions, &fittingLeft, &summaryLeft);
			Solve(solverOptions, &fittingRight, &summaryRight);
			std::vector<double> modelData(t.size());
			for (std::size_t k = 0; k < t.size(); k++) {
				modelData[k] = -A*::sin(2.0*M_PI*f * t[k] / 1000.0);
			}

			double vMaxStim = 2.0 * M_PI*f*A;
			velGain = std::abs((2.0 * M_PI*fleft*Aleft + 2.0 * M_PI*fright*Aright) / 2.0 / vMaxStim);

			//std::unique_ptr<EngbertKlieglThreshold> thres = std::make_unique<EngbertKlieglThreshold>();
		//	thres->SetLambda(6);
			std::unique_ptr<HardThresholdType> thres = std::make_unique<HardThresholdType>();
			thres->SetThreshold(vMaxStim*3.0);
			VelocityBasedSaccadeDetectionAlgorithm det;
			det.SetThreshold(&*thres);
			det.SetMinimumSaccadeDuration(10);
			det.SetMaximumSaccadeSize(30.0);
		
			det.SetTaskType(SaccadeDetectionAlgorithm::FIXATION);
			std::vector<std::size_t> dataIndices(xL.size());
			std::iota(std::begin(dataIndices), std::end(dataIndices), 0); //0 is the starting number

			vecLeft = det.GetSaccadesFromTrace(const_cast<double*>(t.data()), const_cast<double*>(xL.data()), const_cast<double*>(yL.data()), xL.size(), dataIndices);
			vecRight = det.GetSaccadesFromTrace(const_cast<double*>(t.data()), const_cast<double*>(xR.data()), const_cast<double*>(yR.data()), xR.size(), dataIndices);

			det.RemoveMonoculars(vecLeft, vecRight);
			//det.RemoveBinocularOvershoots(vecLeft, vecRight);
			SWJVector swjs = SWJDetection(vecLeft, vecRight);
			freqGain = (fleft + fright) / (this->f*2.0);
			phase = (phaseleft + phaseright) / 2.0*180.0 / M_PI;
			double duration = (t[t.size() - 1] - t[0]) / 1000.0;
			saccRate = double(vecLeft.size()) / duration;
			swjRate = double(swjs.size()) / duration;
			saccAmplitude = 0.0;
			swjAmplitude = 0.0;
			for (std::size_t k = 0; k < vecLeft.size(); k++) {
				saccAmplitude += (vecLeft[k].magnitude + vecRight[k].magnitude) / 2.0;
			}
			saccAmplitude = saccAmplitude / double(vecLeft.size());
			for (std::size_t k = 0; k < swjs.size(); k++) {
				swjAmplitude += (swjs[k].firstMagnitude + swjs[k].secondMagnitude) / 2.0;
			}
			swjAmplitude = swjAmplitude / double(swjs.size());

			//get asymmetry measure:
			// -calulate left and right areas
			// the equation is in (0,1] if the area towards the right is bigger, in [-1,0) if the left is bigger
			double areaLL = 0.0;
			double areaRL = 0.0;
			double areaLR = 0.0;
			double areaRR = 0.0;
			double meanLeft = nanmean(xL.data(), xL.size());
			double meanRight = nanmean(xR.data(), xR.size());
			for (std::size_t k = 1; k < t.size(); k++) {
				if ((xL[k]-meanLeft) > 0) {
					areaRL += std::abs(xL[k]-meanLeft) * (t[k] - t[k - 1]) / 1000.0;
					areaRR += std::abs(xR[k]-meanRight) * (t[k] - t[k - 1]) / 1000.0;
				}
				else if ((xL[k] - meanLeft) < 0) {
					areaLL += std::abs(xL[k] - meanLeft) * (t[k] - t[k - 1]) / 1000.0;
					areaLR += std::abs(xR[k] - meanRight) * (t[k] - t[k - 1]) / 1000.0;
				}
			}

			double asymLeft = (areaRL - areaLL)*M_PI*fleft / Aleft;
			double asymRight = (areaRR - areaLR)*M_PI*fright / Aright;
			pursuitAsymmetry = (asymLeft + asymRight) / 2.0;
			double totalSaccadeDuration = 0.0;
			for (std::size_t k = 0; k < vecLeft.size(); k++) {
				totalSaccadeDuration += (vecLeft[k].duration + vecRight[k].duration) / 2.0;
			}

			fractionTimeInPursuit = 1.0 - totalSaccadeDuration / 1000.0 / duration;


			doneProcessing = true;
			return true;

		}

		double PursuitAnalysis::GetPeakVelocityGain() const {
			if (!doneProcessing) return 0.0;
			return velGain;
		}
		double PursuitAnalysis::GetFrequencyGain() const {
			if (!doneProcessing) return 0.0;
			return freqGain;
		}
		double PursuitAnalysis::GetPhase() const {
			if (!doneProcessing) return 0.0;
			return phase;
		}
		SaccadeVector PursuitAnalysis::GetLeftSaccades() const {
			if (!doneProcessing) return SaccadeVector(0);
			return vecLeft;
		}
		SaccadeVector PursuitAnalysis::GetRightSaccades() const {
			if (!doneProcessing) return SaccadeVector(0);
			return vecRight;
		}

		double PursuitAnalysis::GetSaccadeRate() const {
			return saccRate;
		}
		double PursuitAnalysis::GetSWJRate() const {
			return swjRate;
		}
		double PursuitAnalysis::GetSaccadeAmplitude() const {
			return saccAmplitude;
		}
		double PursuitAnalysis::GetSWJAmplitude() const {
			return swjAmplitude;
		}
		double PursuitAnalysis::GetAsymmetry() const {
			return pursuitAsymmetry;
		}

		bool PursuitAnalysis::Clear() {
			return true;
		}

		Json::Value PursuitAnalysis::ToJSONObject() {
			Json::Value root;

			root["Frequency"] = f;
			root["Amplitude"] = A;
			root["VelocityGain"] = velGain;
			root["FrequencyGain"] = freqGain;
			root["SaccadeRate"] = saccRate;
			root["SaccadeAmplitude"] = saccAmplitude;
			root["FractionTimeInPursuit"] = fractionTimeInPursuit;


			return root;
		}

	
}