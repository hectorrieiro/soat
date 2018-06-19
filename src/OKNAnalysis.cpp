#define EXPORTING
#include "Analysis/OKNAnalysis.h"
#include <memory>
#include <algorithm>
#include <numeric>
#include "Analysis\Algorithms.h"
#include "ceres/ceres.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

namespace Analysis {

	
	OKNAnalysis::OKNAnalysis() {};
	void OKNAnalysis::LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& yL, const std::vector<double>& xR, const std::vector<double>& yR, bool processNow) {
		this->t = t;
		this->xL = xL;
		this->xR = xR;
		this->yL = yL;
		this->yR = yR;
		if (processNow)
			Process();
	}
	

	bool OKNAnalysis::Process() {

		std::unique_ptr<EngbertKlieglThreshold> thres = std::make_unique<EngbertKlieglThreshold>();
		thres->SetLambda(6.0);
		VelocityBasedSaccadeDetectionAlgorithm det;
		det.SetThreshold(&*thres);
		det.SetMinimumSaccadeDuration(10);
		std::vector<std::size_t> dataIndices(xL.size());
		std::iota(std::begin(dataIndices), std::end(dataIndices), 0); //0 is the starting number

		vecLeft = det.GetSaccadesFromTrace(const_cast<double*>(t.data()), const_cast<double*>(xL.data()), const_cast<double*>(yL.data()), xL.size(), dataIndices);
		vecRight = det.GetSaccadesFromTrace(const_cast<double*>(t.data()), const_cast<double*>(xR.data()), const_cast<double*>(yR.data()), xR.size(), dataIndices);

		det.RemoveMonoculars(vecLeft, vecRight);
		det.RemoveBinocularOvershoots(vecLeft, vecRight);
		std::vector<double> mags;
		std::vector<double> durs;
		for (std::size_t k = 1; k < vecLeft.size(); k++) {
			double mL = std::sqrt(pow(vecLeft[k].startX - vecLeft[k - 1].endX, 2.0) + pow(vecLeft[k].startY - vecLeft[k - 1].endY, 2.0));
			double mR = std::sqrt(pow(vecRight[k].startX - vecRight[k - 1].endX, 2.0) + pow(vecRight[k].startY - vecRight[k - 1].endY, 2.0));
			mags.push_back((mL + mR) / 2.0);
			double dL = vecLeft[k].startTimestamp - vecLeft[k - 1].endTimestamp;
			double dR = vecRight[k].startTimestamp - vecRight[k - 1].endTimestamp;
			durs.push_back((dL + dR) / 2.0);
		}
		
		magnitude = nanmedian(mags.data(), mags.size());
		duration = nanmedian(durs.data(), durs.size());
		
		return true;
	}

	

	bool OKNAnalysis::Clear() {
		return true;
	}

	double OKNAnalysis::GetSlowPhasesDuration() const {
		return duration;
	}

	double OKNAnalysis::GetSlowPhasesMagnitude() const {
		return magnitude;
	}

	Json::Value OKNAnalysis::ToJSONObject() {
		Json::Value root;

		root["SlowPhaseDuration"] = duration;
		root["SlowPhaseMagnitude"] = magnitude;


		return root;
	}
}