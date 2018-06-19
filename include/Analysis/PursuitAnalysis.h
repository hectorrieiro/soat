#ifndef PURSUITANALYSIS_H
#define PURSUITANALYSIS_H

#include "Analysis/Algorithms.h"
#include "Analysis/DataTypes.h"
#include <array>

namespace Analysis {

	class DllExport PursuitAnalysis : public JSONifiable {
	public:

		PursuitAnalysis();
		void LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& yL, const std::vector<double>& xR, const std::vector<double>& yR, bool processNow = false);
		void SetPursuitFrequency(double f);
		void SetPursuitAmplitude(double A);
		double GetPeakVelocityGain() const;
		double GetFrequencyGain() const;
		double GetPhase() const;
		double GetSaccadeRate() const;
		double GetSWJRate() const;
		double GetSaccadeAmplitude() const;
		double GetSWJAmplitude() const;
		double GetAsymmetry() const;
		SaccadeVector GetLeftSaccades() const;
		SaccadeVector GetRightSaccades() const;
		bool Process();
		bool Clear();

		Json::Value ToJSONObject() override;

	private:
		std::vector<double> t;
		std::vector<double> xL;
		std::vector<double> yL;
		std::vector<double> xR;
		std::vector<double> yR;
		bool doneProcessing;
		double f;
		double A;
		SaccadeVector vecLeft;
		SaccadeVector vecRight;
		double freqGain, velGain, phase, saccRate, swjRate, saccAmplitude, swjAmplitude, pursuitAsymmetry;
		double fractionTimeInPursuit;
	};

}

#endif
