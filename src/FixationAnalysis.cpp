#define EXPORTING
#include "Analysis/FixationAnalysis.h"
#include "Analysis/Numerics.h"
#include <algorithm>
#include "Analysis/FixationEstabilityAnalysis.h"
#include "Analysis/MicrosaccadeAnalysis.h"

namespace Analysis {

	FixationAnalysis::FixationAnalysis() {
		m_FixationEstabilityAnalysis = std::make_unique<FixationEstabilityAnalysis>();
		m_microsaccadeAnalysis = std::make_unique<MicrosaccadeAnalysis>();
	
		doneProcessing = false;

	}

	void FixationAnalysis::LoadTrials(const ExperimentResultsType& data, const ExperimentEventsType& events, bool processNow) {
		this->data = data;
		this->events = events;

		if (processNow) {
			ProcessAll();
		}

	}

	
	bool FixationAnalysis::ProcessAll() {

		const unsigned int paddingSize = 250;
		
		numPaddingSamples = 0;
		trialEndTimestamps.resize(data.size());
		trialStartTimestamps.resize(data.size());
		for (std::size_t i = 0; i < data.size(); i++) {
			double meanLx = nanmean(data[i].at("glx").data(), data[i].at("glx").size());
			double meanLy = nanmean(data[i].at("gly").data(), data[i].at("gly").size());
			double meanRx = nanmean(data[i].at("grx").data(), data[i].at("grx").size());
			double meanRy = nanmean(data[i].at("gry").data(), data[i].at("gry").size());

			std::for_each(data[i].at("glx").begin(), data[i].at("glx").end(), [this, meanLx](double x) {this->lx.insert(this->lx.end(), x - meanLx); });
			std::for_each(data[i].at("gly").begin(), data[i].at("gly").end(), [this, meanLy](double x) {this->ly.insert(this->ly.end(), x - meanLy); });
			std::for_each(data[i].at("grx").begin(), data[i].at("grx").end(), [this, meanRx](double x) {this->rx.insert(this->rx.end(), x - meanRx); });
			std::for_each(data[i].at("gry").begin(), data[i].at("gry").end(), [this, meanRy](double x) {this->ry.insert(this->ry.end(), x - meanRy); });

			lx.insert(lx.end(), paddingSize, std::nan(""));
			rx.insert(rx.end(), paddingSize, std::nan(""));
			ly.insert(ly.end(), paddingSize, std::nan(""));
			ry.insert(ry.end(), paddingSize, std::nan(""));

			std::vector<double> adjustedT = data[i].at("t");
			double lastTimestamp = i == 0 ? 0 : *(t.end() - 1) + 4;
			std::for_each(adjustedT.begin(), adjustedT.end(), [&lastTimestamp](double& ts) {ts = ts + lastTimestamp; });
			trialStartTimestamps[i] = *adjustedT.begin();
			trialEndTimestamps[i] = *(adjustedT.end() - 1);
			t.insert(t.end(), adjustedT.begin(), adjustedT.end());
			lastTimestamp = *(t.end() - 1);
			auto it = t.insert(t.end(), paddingSize, 0.0);
			std::generate_n(it, paddingSize, [&lastTimestamp] {lastTimestamp = lastTimestamp + 4; return lastTimestamp; });
			numPaddingSamples += paddingSize;
		}

		std::vector<double> dlx(lx.size()), dly(lx.size()), drx(lx.size()), dry(lx.size());
		firstDerivative(lx.data(), dlx.data(), lx.size());
		firstDerivative(ly.data(), dly.data(), ly.size());
		firstDerivative(rx.data(), drx.data(), rx.size());
		firstDerivative(ry.data(), dry.data(), ry.size());
		m_FixationEstabilityAnalysis->LoadTraces(t, lx, ly, rx, ry, false);
		m_FixationEstabilityAnalysis->LoadVelocities(dlx, dly, drx, dry);
		m_FixationEstabilityAnalysis->Process();
		m_microsaccadeAnalysis->LoadTraces(t, lx, ly, rx, ry, true);

		doneProcessing = true;
		return true;

	}
	const MicrosaccadeAnalysis& FixationAnalysis::GetMicrosaccadeAnalysis() const {
		if (doneProcessing) 
			return *m_microsaccadeAnalysis;
		else {
			//throw exception?
		}
	}

	const FixationEstabilityAnalysis& FixationAnalysis::GetFixationEstabilityAnalysis() const {
		if (doneProcessing)
		return *m_FixationEstabilityAnalysis;
		else {
			//throw exception?
		}

	}
	void FixationAnalysis::Clear() {
		doneProcessing = false;

	}

	unsigned int FixationAnalysis::GetNumberOfTrials() const {
		if (!doneProcessing) return 0;
		return trialEndTimestamps.size();
	}

	double FixationAnalysis::GetRecordingDuration() const {
		return 1.0/250.0*double(lx.size() - numPaddingSamples);
	}

	void FixationAnalysis::GetTracesForTrial(unsigned int trialNumber, std::vector<double>& t, std::vector<double>& xL, std::vector<double>& yL, std::vector<double>& xR, std::vector<double>& yR) const {
		if (!doneProcessing || trialNumber > trialStartTimestamps.size()) return;
		auto itStart = std::find(this->t.begin(), this->t.end(), trialStartTimestamps[trialNumber]);
		auto itEnd = std::find(this->t.begin(), this->t.end(), trialEndTimestamps[trialNumber]);

		std::size_t startIndex = itStart - this->t.begin();
		std::size_t endIndex = itEnd - this->t.begin();
		t.resize(endIndex - startIndex-1);
		xL.resize(endIndex - startIndex-1);
		xR.resize(endIndex - startIndex-1);
		yL.resize(endIndex - startIndex-1);
		yR.resize(endIndex - startIndex-1);
		std::copy(itStart, itEnd-1, t.begin());
		std::copy(lx.begin() + startIndex, lx.begin() + endIndex-1, xL.begin());
		std::copy(ly.begin() + startIndex, ly.begin() + endIndex-1, yL.begin());
		std::copy(rx.begin() + startIndex, rx.begin() + endIndex-1, xR.begin());
		std::copy(ry.begin() + startIndex, ry.begin() + endIndex-1, yR.begin());

	}

	void FixationAnalysis::GetFullTraces(std::vector<double>& t, std::vector<double>& xL, std::vector<double>& yL, std::vector<double>& xR, std::vector<double>& yR) const {
		t = this->t;
		xL = lx;
		yL = ly;
		xR = rx;
		yR = ry;
	}

	Json::Value FixationAnalysis::ToJSONObject() {
		Json::Value root;

		root["Microsaccades"]["MicrosaccadeRate"] = m_microsaccadeAnalysis->GetNumberOfMicrosaccades()/GetRecordingDuration();
		root["Microsaccades"]["MicrosaccadeMagnitude"] = m_microsaccadeAnalysis->GetMeanMagnitude();
		std::vector<double> mags, pkvel;
		m_microsaccadeAnalysis->GetMagnitudePeakVelocityMainSequence(mags, pkvel);
		root["Microsaccades"]["MicrosaccadeMainSequenceSlope"] = linearSlope(mags, pkvel);
		root["Microsaccades"]["SWJRate"] = m_microsaccadeAnalysis->GetNumberOfSWJs()/GetRecordingDuration();
		MicrosaccadeAnalysis::SaccadeParameterDistribution swjMag, nonSwjMag;
		MicrosaccadeAnalysis::SaccadeParameterDistribution swjVert, nonSwjVert;
		m_microsaccadeAnalysis->GetSWJMagnitude(swjMag, nonSwjMag);
		m_microsaccadeAnalysis->GetSWJVerticality(swjVert, nonSwjVert);
		root["Microsaccades"]["SWJMagnitude"] = swjMag.mean;
		root["Microsaccades"]["NonSWJMicrosaccadeVerticalComponent"] = nonSwjVert.mean;
		root["Microsaccades"]["SWJMicrosaccadeVerticalComponent"] = swjVert.mean;
		root["Microsaccades"]["NonSWJMicrosaccadeMagnitude"] = nonSwjMag.mean;
		root["Microsaccades"]["SWJMicrosaccadeMagnitude"] = swjMag.mean;
		double harmVelLeft, harmVelRight;
		m_FixationEstabilityAnalysis->GetHarmonicVelocities(harmVelLeft, harmVelRight);
		root["HarmonicVelocity"] = (harmVelLeft + harmVelRight) / 2.0;
		double leftEnt, rightEnt;
		m_FixationEstabilityAnalysis->GetEntropies(leftEnt, rightEnt);
		root["GazeEntropy"] = (leftEnt + rightEnt)/2.0;
		double leftBCEA, rightBCEA;
		
		m_FixationEstabilityAnalysis->GetBCEA(leftBCEA, rightBCEA);

		root["BCEA"] = (leftBCEA+rightBCEA)/2.0;
		return root;
	}

}