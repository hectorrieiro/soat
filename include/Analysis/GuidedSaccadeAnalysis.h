#ifndef GUIDEDSACCADEANALYSIS_H
#define GUIDEDSACCADEANALYSIS_H

#include "Analysis/Algorithms.h"
#include "Analysis/DataTypes.h"
#include <array>

namespace Analysis {

	class DllExport GuidedSaccadeAnalysis : public JSONifiable {
	public:
		typedef unsigned int SaccadeCueProperties;
		typedef enum SaccadeCuePropertiesEnum_: unsigned int { NONE = 0, TO_CENTER = 1, TO_PERIPHERY = 2, LEFT = 4, RIGHT = 8, UP = 16, DOWN = 32, ALL = 0xFFFFFF} SaccadeFilterType;
		
		GuidedSaccadeAnalysis();

		void LoadTrials(const TaskDataType& data, const TaskEventsType& events, bool processNow = false);
		
		void GetMagnitudeGain(double&, double&, unsigned int filter = ALL) const;
		void GetAverageLatency(double&, double&, unsigned int filter = ALL) const;
		double GetFractionMultistep(unsigned int filter = ALL) const;
		unsigned int GetNumberOfMisses(unsigned int filter = ALL) const;
		SaccadeVector GetLeftSaccades(unsigned int filter = ALL) const;
		SaccadeVector GetRightSaccades(unsigned int filter = ALL) const;
		unsigned int GetNumberUncued(unsigned int filter = ALL) const;
		double GetMainSequenceSlope(unsigned int filter = ALL) const;
		double GetMultistepMagnitudeFrac(unsigned int filter = ALL) const;
		std::vector<double> GetSaccadicWaveforms(std::vector<std::vector<double> >& lw, std::vector<std::vector<double> > & rw, unsigned int filter = ALL, bool amplitudeNormalized = false);
		bool ProcessAll();
		bool Clear();
		Json::Value ToJSONObject() override;

	private:
		void SetTargetPositions();
		TaskDataType data;
		TaskEventsType events;
		std::vector<double> t;
		std::vector<double> xL;
		std::vector<double> yL;
		std::vector<double> xR;
		std::vector<double> yR;
		bool doneProcessing;
		SaccadeVector vecLeft;
		SaccadeVector vecRight;
		std::vector<unsigned int> saccadeCues;
		std::vector<double> cuedMagnitudes;
		std::vector<double> targetT, targetX, targetY;
		std::vector<double> latenciesLeft, latenciesRight;
		std::vector<double> magnitudeGainLeft, magnitudeGainRight;
		std::vector<double> peakVelocitiesLeft, peakVelocitiesRight;
		std::vector<double> magnitudesLeft, magnitudesRight;
		std::vector<bool> multistep;
		std::vector<double> multistepFractionMagnitude;
		std::vector<unsigned int> numberOfSaccadesPerTarget;
		std::vector<bool> hitOrMiss;
		std::vector<unsigned int> numberUncued;
		std::size_t waveformDuration;
		std::vector<std::vector<double> > leftWaveforms, rightWaveforms;
		std::vector<double> slopes;

		std::vector<std::size_t> filteredIndices(unsigned int filter) const;
		
	};

}

#endif
