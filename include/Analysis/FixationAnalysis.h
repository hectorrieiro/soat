#ifndef FIXATIONANALYSIS_H
#define FIXATIONANALYSIS_H

#include "Analysis\DataTypes.h"
#include <memory>

namespace Analysis {

	class MicrosaccadeAnalysis;

	class FixationEstabilityAnalysis;

	class DllExport FixationAnalysis : public JSONifiable {
	public:
		FixationAnalysis();
		void LoadTrials(const ExperimentResultsType& data, const ExperimentEventsType& events, bool processNow = false);
		bool ProcessAll();
		const MicrosaccadeAnalysis& GetMicrosaccadeAnalysis() const;
		const FixationEstabilityAnalysis& GetFixationEstabilityAnalysis() const;
		double GetRecordingDuration() const;
		unsigned int GetNumberOfTrials() const;
		void GetTracesForTrial(unsigned int trialNumber, std::vector<double>& t, std::vector<double>& xL, std::vector<double>& yL, std::vector<double>& xR, std::vector<double>& yR) const;
		void GetFullTraces(std::vector<double>& t, std::vector<double>& xL, std::vector<double>& yL, std::vector<double>& xR, std::vector<double>& yR) const;
		void Clear();
		virtual Json::Value ToJSONObject() override;

	private:
		bool doneProcessing;
		std::unique_ptr<MicrosaccadeAnalysis> m_microsaccadeAnalysis;

		std::unique_ptr<FixationEstabilityAnalysis> m_FixationEstabilityAnalysis;
		ExperimentResultsType data;
		ExperimentEventsType events;
		std::vector<double> t, lx, rx, ly, ry;
		unsigned long int numPaddingSamples;
		std::vector<double> trialStartTimestamps;
		std::vector<double> trialEndTimestamps;

	};
	

}

#endif