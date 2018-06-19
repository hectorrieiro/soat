#define EXPORTING
#include "Analysis/MicrosaccadeAnalysis.h"
#include <memory>
#include <algorithm>
#include <numeric>
#include "Analysis\Algorithms.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

namespace Analysis {
	MicrosaccadeAnalysis::MicrosaccadeAnalysis() {
		doneProcessing = false;


	}
	void MicrosaccadeAnalysis::LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& yL, const std::vector<double>& xR, const std::vector<double>& yR, bool processNow) {
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
	
	void MicrosaccadeAnalysis::SetMaximumMagnitudeThreshold(double) {

	}
	double MicrosaccadeAnalysis::GetMaximumMagnitudeThreshold() const {
		return 5.0;
	}
	bool MicrosaccadeAnalysis::Process() {
		
		std::unique_ptr<EngbertKlieglThreshold> thres = std::make_unique<EngbertKlieglThreshold>();
		thres->SetLambda(5.0);
		VelocityBasedSaccadeDetectionAlgorithm det;
		det.SetThreshold(&*thres);
		det.SetMinimumSaccadeDuration(6);
		std::vector<std::size_t> dataIndices(xL.size());
		std::iota(std::begin(dataIndices), std::end(dataIndices), 0); //0 is the starting number

		vecLeft = det.GetSaccadesFromTrace(const_cast<double*>(t.data()), const_cast<double*>(xL.data()), const_cast<double*>(yL.data()), xL.size(), dataIndices);
		vecRight = det.GetSaccadesFromTrace(const_cast<double*>(t.data()), const_cast<double*>(xR.data()), const_cast<double*>(yR.data()), xR.size(), dataIndices);

		det.RemoveMonoculars(vecLeft, vecRight);
		det.RemoveBinocularOvershoots(vecLeft, vecRight);
		vels.resize(vecLeft.size());
		mags.resize(vecLeft.size());
		dirs.resize(vecLeft.size());
		for (std::size_t k = 0; k < vels.size(); k++) {
			vels[k] = (vecLeft[k].peakVelocity + vecRight[k].peakVelocity) / 2.0;
			mags[k] = (vecLeft[k].magnitude + vecRight[k].magnitude) / 2.0;
			dirs[k] = (vecLeft[k].direction + vecRight[k].direction) / 2.0;
			bool changeAngleSign = (vecLeft[k].direction > 90 && vecLeft[k].direction < 180 && vecRight[k].direction < -90 && vecRight[k].direction > -180) ||
				(vecLeft[k].direction < -90 && vecLeft[k].direction > -180 && vecRight[k].direction > 90 && vecRight[k].direction < 180);
			if (changeAngleSign) {
				dirs[k] = dirs[k] > 0.0 ? 180.0 - dirs[k] : -180.0 - dirs[k];
			}


		}

		swjs = SWJDetection(vecLeft, vecRight);
		doneProcessing = true;
		return true;
	}

	void MicrosaccadeAnalysis::GetPolarHistogram(std::vector<double>& histCenters, std::vector<double>& N, unsigned int nBins) const {
		histCenters.resize(nBins);
		std::vector<double> histEdges(nBins+1);
		double binWidth = 360.0 / double(nBins);
		double edge = -190.0;
		double center = -185.0;
		std::generate_n(histEdges.begin(), nBins+1, [&edge, binWidth] {edge = edge + binWidth; return edge; });
		std::generate_n(histCenters.begin(), nBins, [&center, binWidth] {center = center + binWidth; return center; });
		N = HistogramFromEdges(dirs, histEdges);
		N.push_back(N[0]);
		histCenters.push_back(histCenters[0]);
		std::for_each(N.begin(), N.end(), [this](double& x) { x = x / double(this->vecLeft.size()); });

	}

	void MicrosaccadeAnalysis::GetMagnitudeHistogram(std::vector<double>& histCenters, std::vector<double>& N, unsigned int nBins, double maxMagnitude) const {
		histCenters.resize(nBins);
		std::vector<double> histEdges(nBins+1);
		double binWidth = maxMagnitude / double(nBins);
		double minMag = -binWidth;
		double minCenter = -binWidth / 2.0;
		std::generate_n(histEdges.begin(), nBins+1, [&minMag, binWidth] {minMag += binWidth; return minMag; });
		std::generate_n(histCenters.begin(), nBins, [&minCenter, binWidth] {minCenter += binWidth; return minCenter; });
		N = HistogramFromEdges(mags, histEdges);
		std::for_each(N.begin(), N.end(), [this](double& x) {x = x / double(this->vecLeft.size()); });
	}

	unsigned int MicrosaccadeAnalysis::GetNumberOfMicrosaccades() const {
		return mags.size();
	}

	unsigned int MicrosaccadeAnalysis::GetNumberOfSWJs() const {
		return swjs.size();
	}


	double MicrosaccadeAnalysis::GetMeanMagnitude() const {
		return nanmean(mags.data(), mags.size());
	}

	double MicrosaccadeAnalysis::GetSTDMagnitude() const {
		return std::sqrt(nanvar(mags.data(), GetMeanMagnitude(), mags.size()));
	}

	bool MicrosaccadeAnalysis::Clear() {
		doneProcessing = false;
		return true;

	}

	BinocularSaccades MicrosaccadeAnalysis::GetMicrosaccadesInInterval(double tStart, double tEnd) const {
		BinocularSaccades out;
		SaccadeVector outLeft;
		SaccadeVector outRight;
		for (std::size_t k = 0; k < vecLeft.size(); k++) {
			if (vecLeft[k].startTimestamp >= tStart) {
				if (vecLeft[k].startTimestamp > tEnd) break;
				outLeft.push_back(vecLeft[k]);
				outRight.push_back(vecRight[k]);
			}
		}
		out[LEFT] = outLeft;
		out[RIGHT] = outRight;
		return out;
	}

	void MicrosaccadeAnalysis::GetMagnitudePeakVelocityMainSequence(std::vector<double>& mags, std::vector<double>& peakVel) const {

		mags = this->mags;
		peakVel = this->vels;
	}
	

	void MicrosaccadeAnalysis::GetSWJMagnitude(SaccadeParameterDistribution& swj, SaccadeParameterDistribution& nonSwj) const {
		std::vector<double> swjMags, nonSwjMags;
		for (std::size_t k = 0; k < vecLeft.size(); k++) {
			if (vecLeft[k].inSWJ) {
				swjMags.push_back(vecLeft[k].magnitude);
			}
			else {
				nonSwjMags.push_back(vecLeft[k].magnitude);
			}
		}

		swj = GetDistributionForParameterArray(swjMags);
		nonSwj = GetDistributionForParameterArray(nonSwjMags);

	}

	void MicrosaccadeAnalysis::GetSWJPeakVelocity(SaccadeParameterDistribution& swj, SaccadeParameterDistribution& nonSwj) const {
		std::vector<double> swjMags, nonSwjMags;
		for (std::size_t k = 0; k < vecLeft.size(); k++) {
			if (vecLeft[k].inSWJ) {
				swjMags.push_back(vecLeft[k].peakVelocity);
			}
			else {
				nonSwjMags.push_back(vecLeft[k].peakVelocity);
			}
		}

		swj = GetDistributionForParameterArray(swjMags);
		nonSwj = GetDistributionForParameterArray(nonSwjMags);

	}

	void MicrosaccadeAnalysis::GetSWJVerticality(SaccadeParameterDistribution& swj, SaccadeParameterDistribution& nonSwj) const {
		std::vector<double> swjMags, nonSwjMags;
		for (std::size_t k = 0; k < vecLeft.size(); k++) {
			if (vecLeft[k].inSWJ) {
				swjMags.push_back(std::abs(std::sin(vecLeft[k].direction/180.0*M_PI)));
			}
			else {
				nonSwjMags.push_back(std::abs(std::sin(vecLeft[k].direction / 180.0*M_PI)));
			}
		}

		swj = GetDistributionForParameterArray(swjMags);
		nonSwj = GetDistributionForParameterArray(nonSwjMags);

	}

	double MicrosaccadeAnalysis::GetMainSequenceSlope() const {
		double slope = 0.0;
		for (std::size_t k = 0; k < mags.size(); k++) {
			slope += vels[k] / mags[k];
		}

		slope /= mags.size();
		return slope;

	}
	double MicrosaccadeAnalysis::GetMedianVertComponent() const {
		std::vector<double> cmp(dirs.size());
		for (std::size_t k = 0; k < dirs.size(); k++) {
			cmp[k] = abs(tan(dirs[k] / 180.0*M_PI));
		}

		
		return nanmedian(cmp.data(), cmp.size());

	}

	MicrosaccadeAnalysis::SaccadeParameterDistribution MicrosaccadeAnalysis::GetDistributionForParameterArray(const std::vector<double> x) const {
		SaccadeParameterDistribution str;
		if (x.size() == 0) {
			return str;
		}
		
		str.mean = nanmean(x.data(), x.size());
		std::vector<double> xsorted(x.size());
		std::copy(x.begin(), x.end(), xsorted.begin());
		std::sort(xsorted.begin(), xsorted.end());
		str.median = xsorted.size() % 2 == 0 ? (xsorted[xsorted.size() / 2] + xsorted[xsorted.size() / 2 - 1]) / 2.0 : xsorted[xsorted.size() / 2];
		str.std = nanvar(x.data(), str.mean, x.size());
		str.std = std::sqrt(str.std*str.std);
		str.sem = str.std / std::sqrt(x.size());
		double pct2_5Idx = double(x.size()) / 40.0;
		double pct97_5Idx = 39.0 / 40.0*double(x.size());
		double pct25Idx = double(x.size()) / 4.0;
		double pct75Idx = 3.0 / 4.0*double(x.size());
		double wLow = std::floor(pct2_5Idx) / (std::floor(pct2_5Idx) + std::ceil(pct2_5Idx));
		double wHigh = std::ceil(pct2_5Idx) / (std::floor(pct2_5Idx) + std::ceil(pct2_5Idx));
		str.confidence95[0] = wLow*xsorted[std::floor(pct2_5Idx)] + wHigh*xsorted[std::ceil(pct2_5Idx)];
		wLow = std::floor(pct97_5Idx) / (std::floor(pct97_5Idx) + std::ceil(pct97_5Idx));
		wHigh = std::ceil(pct97_5Idx) / (std::floor(pct97_5Idx) + std::ceil(pct97_5Idx));
		str.confidence95[1] = wLow*std::floor(pct97_5Idx) + wHigh*std::ceil(pct97_5Idx);
		wLow = std::floor(pct25Idx) / (std::floor(pct25Idx) + std::ceil(pct25Idx));
		wHigh = std::ceil(pct25Idx) / (std::floor(pct25Idx) + std::ceil(pct25Idx));
		str.confidence50[0] = wLow*xsorted[std::floor(pct25Idx)] + wHigh*xsorted[std::ceil(pct25Idx)];
		wLow = std::floor(pct75Idx) / (std::floor(pct75Idx) + std::ceil(pct75Idx));
		wHigh = std::ceil(pct75Idx) / (std::floor(pct75Idx) + std::ceil(pct75Idx));
		str.confidence50[1] = wLow*xsorted[std::floor(pct75Idx)] + wHigh*xsorted[std::min(std::size_t(std::ceil(pct75Idx)), xsorted.size()-1)];

		return str;
	}
}