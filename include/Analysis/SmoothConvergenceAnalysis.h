#ifndef SMOOTHCONVERGENCEANALYSIS_H
#define SMOOTHCONVERGENCEANALYSIS_H

#include "Analysis/Algorithms.h"
#include "Analysis/DataTypes.h"
#include <array>

namespace Analysis {

	class DllExport SmoothConvergenceAnalysis : public JSONifiable {
	public:

		SmoothConvergenceAnalysis();
		void LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& xR,  bool processNow = false);
		void SetVelocity(double v);
		void SetStartingPosition(double d0);
		void SetDirectionChangeTimestamps(std::vector<double> t);
		double GetApproachingConvergenceDistance() const;
		unsigned int GetApproachingMisses() const;
		void GetApproachingWaveforms(std::vector<std::vector<double> >& d, std::vector<std::vector<double> >& v) const;
		double GetMedianApproachingNPC() const;
		double GetMedianDistancingNPC() const;
		bool Process();
		bool Clear();
		virtual Json::Value ToJSONObject() override;

	private:
		std::vector<double> t;
		std::vector<double> xL;
		std::vector<double> xR;
		std::vector<double> tDirectionChanges;
		bool doneProcessing;
		double v;
		double d0;
		SaccadeVector vecLeft;
		SaccadeVector vecRight;
		double freqGain, velGain, phase;
		std::vector<double> verg;
		std::vector<double> d;
		static void trialFitting(const std::vector<double>& v, const std::vector<double>& d, double& a, double& b, double& c, double& baseline, double& G);
		std::vector<std::vector<double> > approachingVergences, approachingDistances, distancingVergences, distancingDistances;
		std::vector<double> approachingSlopes, approachingNPC, distancingSlopes, distancingNPC;
		
	};

}

#endif
