#ifndef PUPILLOMETRYANALYSIS_H
#define PUPILLOMETRYANALYSIS_H

#include "Analysis/Algorithms.h"
#include "Analysis/DataTypes.h"
#include <array>

namespace Analysis {

	class DllExport PupillometryAnalysis : public JSONifiable{
	public:
		
		PupillometryAnalysis();

		void LoadTrials(const TaskDataType& data, const TaskEventsType& events, bool processNow = false);
		
		void GetMaxConstriction(double& left, double& right) const;
		void GetTimeToMaxConstriction(double& left, double& right) const;
		void GetTimeToRecovery(double& left, double& right) const;
		bool ProcessAll();
		bool Clear();
		void GetTimeAroundFlash(std::vector<double>& t) const;
		virtual Json::Value ToJSONObject() override;

	private:
		bool doneProcessing;
		TaskDataType data;
		TaskEventsType events;
		static const unsigned int logitNComponents;
		double maxContractionLeft;
		double maxContractionRight;
		double timeToMaxContractionLeft, timeToMaxContractionRight;
		std::vector<double> timeAroundFlash;
		double timeToRecoveryLeft, timeToRecoveryRight;
		
	};

}

#endif
