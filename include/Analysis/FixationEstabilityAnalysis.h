#ifndef FIXATIONESTABILITYANALYSIS_H
#define FIXATIONESTABILITYANALYSIS_H
#include <array>
#include "Analysis/DataTypes.h"



namespace Analysis {

	class DllExport FixationEstabilityAnalysis {
	public:
		typedef std::array<double, 2> BinocularDataPointType;
		typedef std::array<std::array<double, 2>, 2> MainComponentDirectionsType;
		typedef std::array<double, 2> MainComponentsType;
		FixationEstabilityAnalysis();
		void LoadTraces(const std::vector<double>& t, const std::vector<double>& xL, const std::vector<double>& yL, const std::vector<double>& xR, const std::vector<double>& yR, bool processNow = false);
		void LoadVelocities(const std::vector<double>& vxL, const std::vector<double>& vyL, const std::vector<double>& vxR, const std::vector<double>& vyR);
		void Process();
		void GetEntropies(double& left, double& right) const;
		void GetVariances(double& lh, double& lv, double& rh, double& rv) const;
		void GetHarmonicVelocities(double& left, double& right) const;
		void GetMainAxes(MainComponentDirectionsType&, MainComponentDirectionsType&) const;
		void GetBCEA(double& lBCEA, double& rBCEA) const;
		void Clear();
	private:
		std::vector<double> t;
		std::vector<double> xL;
		std::vector<double> yL;
		std::vector<double> xR;
		std::vector<double> yR;
		std::vector<double> vxL;
		std::vector<double> vyL;
		std::vector<double> vxR;
		std::vector<double> vyR;
		bool doneProcessing;
		typedef std::array<double, 2> BinocularDataPointType;
		typedef std::array<std::array<double, 2>, 2> MainComponentDirectionsType;
		typedef std::array<double, 2> MainComponentsType;
		BinocularDataPointType entropies;
		BinocularDataPointType harmonicVelocities;
		std::array<MainComponentDirectionsType, 2>  mainComponentDirections;
		std::array<MainComponentsType, 2> variances;

		void GetDispersionFromTraces(const std::vector<double>&, const std::vector<double>&, MainComponentsType&, MainComponentDirectionsType&);


	};

	
}

#endif
