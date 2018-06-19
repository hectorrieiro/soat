#ifndef MICROSACCADEANALYSIS_H
#define MICROSACCADEANALYSIS_H
#include "Analysis\DataTypes.h"
#include "Analysis/Algorithms.h"
#include <array>

namespace Analysis {

	class DllExport MicrosaccadeAnalysis {
	public:
		typedef struct {
			std::array<double, 2> confidence95 = { std::nan(""),std::nan("") };
			std::array<double, 2> confidence50 = { std::nan(""),std::nan("") };
			double mean = std::nan("");
			double median = std::nan("");
			double std = std::nan("");
			double sem = std::nan("");
		} SaccadeParameterDistribution;

		MicrosaccadeAnalysis();
		void LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& yL, const std::vector<double>& xR, const std::vector<double>& yR, bool processNow = false);
		void SetMaximumMagnitudeThreshold(double);
		double GetMaximumMagnitudeThreshold() const;
		void GetPolarHistogram(std::vector<double>& histCenters, std::vector<double>& N, unsigned int nBins = 36) const;
		void GetMagnitudeHistogram(std::vector<double>& histCenters, std::vector<double>& N, unsigned int nBins = 25, double maxMagnitude = 5.0) const;
		unsigned int GetNumberOfMicrosaccades() const;
		double GetMainSequenceSlope() const;
		double GetMedianVertComponent() const;
		double GetMeanMagnitude() const;
		double GetSTDMagnitude() const;
		BinocularSaccades GetMicrosaccadesInInterval(double tStart, double tEnd) const;
		void GetMagnitudePeakVelocityMainSequence(std::vector<double>& mags, std::vector<double>& peakVel) const;
		void GetSWJMagnitude(SaccadeParameterDistribution&, SaccadeParameterDistribution&) const;
		void GetSWJPeakVelocity(SaccadeParameterDistribution&, SaccadeParameterDistribution&) const;
		void GetSWJVerticality(SaccadeParameterDistribution&, SaccadeParameterDistribution&) const;
		unsigned int GetNumberOfSWJs() const;
		bool Process();
		bool Clear();

	private:
		std::vector<double> t;
		std::vector<double> xL;
		std::vector<double> yL;
		std::vector<double> xR;
		std::vector<double> yR;
		SaccadeVector vecLeft, vecRight;
		SWJVector swjs;
		std::vector<double> vels;
		std::vector<double> mags;
		std::vector<double> dirs;
		bool doneProcessing;
		double maxMagnitudeThreshold;
		unsigned long int totalPaddingSamples;
		SaccadeParameterDistribution GetDistributionForParameterArray(const std::vector<double>) const;
	};

}

#endif
