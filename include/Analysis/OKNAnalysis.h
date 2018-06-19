#ifndef OKNANALYSIS_H
#define OKNANALYSIS_H

#include "Analysis/Algorithms.h"
#include "Analysis\DataTypes.h"
#include <array>

namespace Analysis {

	class DllExport OKNAnalysis : public JSONifiable{
	public:

		OKNAnalysis();
		void LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& yL, const std::vector<double>& xR, const std::vector<double>& yR, bool processNow = false);
		double GetSlowPhasesMagnitude() const;
		double GetSlowPhasesDuration() const;
		bool Process();
		bool Clear();
		Json::Value ToJSONObject() override;

	private:
		std::vector<double> t;
		std::vector<double> xL;
		std::vector<double> xR;
		std::vector<double> yL;
		std::vector<double> yR;
		double duration;
		double magnitude;
		SaccadeVector vecLeft, vecRight;
		
	};

}

#endif
