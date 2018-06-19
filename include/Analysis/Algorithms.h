#ifndef REPORTSALGORITHMS_H
#define REPORTSALGORITHMS_H

#include "Analysis\DataTypes.h"
#include "Analysis\Numerics.h"
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <fstream>

#define M_PI 3.14159265

namespace Analysis {

	typedef struct {
		double startTimestamp;
		double endTimestamp;
		double startX;
		double startY;
		double endX;
		double endY;
		double duration;
		double magnitude;
		double direction;
		double peakVelocity;
		double meanVelocity;
		bool hasOvershoot;
		bool inSWJ = false;
		std::size_t swjIndex;
		double overshootMagnitude;
		double overshootDuration;
		double overshootMeanVelocity;
		double overshootPeakVelocity;
		double overshootStartX;
		double overshootStartY;
		double overshootStartTimestamp;
		double overshootEndTimestamp;
		double overshootEndX;
		double overshootEndY;
	} SaccadeMonocularProperties;

	typedef enum {SWJ_LEFT, SWJ_RIGHT, SWJ_BOTH} SWJEyeEnum;
	typedef struct {
		SWJEyeEnum eye;
		std::size_t firstSaccadeIndex;
		std::size_t secondSaccadeIndex;
		double dirDifference;
		double isi;
		double firstDirection;
		double secondDirection;
		double firstMagnitude;
		double secondMagnitude;
		double swjScore;

	} SWJProperties;

	typedef std::vector<SWJProperties> SWJVector;
	typedef std::vector<SaccadeMonocularProperties> SaccadeVector;
	typedef enum {LEFT, RIGHT} SaccadeEyeEnum;
	typedef std::map<SaccadeEyeEnum, SaccadeVector> BinocularSaccades;


	class BlinkRemover {
	public:
		BlinkRemover() : outputBlinkValue(std::nan("")), inputBlinkValue(std::nan("")), bufferSamples(50), interpolateMissingSamples(true), minimumBlinkDuration(15) {};
		void SetInputBlinkValue(double);
		void SetOutputBlinkValue(double);
		void SetBufferSamples(unsigned int);
		void SetInterpolateMissingSamples(bool);
		void SetMinimumBlinkDuration(unsigned int samples);
		unsigned int operator()(double* x, double*y, std::size_t size);
	private:
		double outputBlinkValue;
		double inputBlinkValue;
		unsigned int bufferSamples;
		bool interpolateMissingSamples;
		unsigned int minimumBlinkDuration;
	};

	class VelocityThresholdType {
	public:
		virtual void AdjustThreshold(const double* vx, const double* vy, std::size_t size, std::vector<std::size_t> fixationIndices = {}) = 0;
		virtual bool AboveThreshold(const double* vx, const double* vy, std::size_t) = 0;
		virtual bool AboveReducedThreshold(const double* vx, const double* vy, std::size_t) = 0;
	};

	class HardThresholdType : public VelocityThresholdType {
	public:
		void SetThreshold(double t);
		void AdjustThreshold(const double* vx, const double* vy, std::size_t size, std::vector<std::size_t> fixationIndices) override;
		bool AboveThreshold(const double* vx, const double* vy, std::size_t);
		bool AboveReducedThreshold(const double* vx, const double* vy, std::size_t);
	private:
		double _t;
	};

	class EngbertKlieglThreshold : public VelocityThresholdType {
	public:
		void SetLambda(double l);
		void AdjustThreshold(const double* vx, const double* vy, std::size_t size, std::vector<std::size_t> fixationIndices) override;
		bool AboveThreshold(const double* vx, const double* vy, std::size_t) override;
		bool AboveReducedThreshold(const double* vx, const double* vy, std::size_t);
	private:
		double l;
		double thresholdX, thresholdY;
	};

	class SaccadeDetectionAlgorithm {
	public:
		typedef enum { FIXATION, PURSUIT } TaskType;
		SaccadeDetectionAlgorithm() : _detectOvershoots(true), minISI(20), doFilterBinocular(true), minDuration(15), magnitudeMinThreshold(0.2), magnitudeMaxThreshold(5.0), taskType(FIXATION) {};
		void DetectOvershoots(bool);
		void SetMinISI(double isi);
		void SetThresholds(double min, double max);
		void FilterBinocular(bool);
		void SetMinimumSaccadeDuration(double ms);
		void SetMaximumSaccadeSize(double deg);
		void SetTaskType(TaskType t);
		
		virtual SaccadeVector GetSaccadesFromTrace(double* t, double* x, double* y, std::size_t size, std::vector<std::size_t>& fixationIndices) = 0;
		void RemoveMonoculars(SaccadeVector&, SaccadeVector&);
		void RemoveOvershoots(SaccadeVector&);
		void RemoveBinocularOvershoots(SaccadeVector&, SaccadeVector&);
	protected:
		bool _detectOvershoots;
		double minISI;
		bool doFilterBinocular;
		double minDuration;
		double magnitudeMaxThreshold;
		double magnitudeMinThreshold;
		TaskType taskType;

		static int TimingOverlap(const SaccadeMonocularProperties& s1, const SaccadeMonocularProperties& s2);
		
	};

	DllExport SWJVector SWJDetection(SaccadeVector& vecLeft, SaccadeVector& vecRight);

    class VelocityBasedSaccadeDetectionAlgorithm : public SaccadeDetectionAlgorithm {
	public:
		VelocityBasedSaccadeDetectionAlgorithm() = default;
		void SetThreshold(VelocityThresholdType* t) {
			this->thresholdFcn = t;
		};

		SaccadeVector GetSaccadesFromTrace(double* t, double* x, double* y, std::size_t size, std::vector<std::size_t>& fixationIndices) {

			//blink detection
			BlinkRemover br;
			const double samplingFrequency = 250;
			unsigned int numBlinks = br(x, y, size);
			
			//velocity calculation
			std::vector<double> rho(size);
			std::vector<double> theta(size);
			std::vector<double> dx(size);
			std::vector<double> dy(size);
			if (taskType == FIXATION) {
				firstDerivative(x, dx.data(), size, samplingFrequency);
				firstDerivative(y, dy.data(), size, samplingFrequency);
			//	polarCoordinates(dx.data(), dy.data(), rho.data(), theta.data(), size);

			}
			else if (taskType == PURSUIT) {
				dx[0] = 0.0;
				dy[0] = 0.0;
				for (std::size_t k = 1; k < size; k++) {
					dx[k] = (x[k] - x[k - 1])*samplingFrequency;
					dy[k] = (y[k] - y[k - 1])*samplingFrequency;
				}
				const std::vector<double> B = { 0.97098210991349054310717292537447065115,
					- 2.91248634174777709660020263982005417347,
					2.91248634174777709660020263982005417347,
					- 0.97098210991349054310717292537447065115 };
				const std::vector<double> A = { 1.0,
					- 2.940643021552723723743838490918278694153,
					2.883487624020730777374410536140203475952,
					- 0.942806257749079668073477478174027055502 };

				ZeroPhaseIIRFilter H;
				H.SetDenominator(A);
				H.SetNumerator(B);
				std::vector<double> dxF(size);
				std::vector<double> dyF(size);
				H(dx, dxF);
				H.ResetFilter();
				H(dy, dyF);
				dx = dxF;
				std::copy(dxF.begin(), dxF.end(), dy.begin());
				//dy.resize(dxF.size());
				//dy.assign(dy.size(), 0.0);
			
				
			}
			//threshold calculation
			
			this->thresholdFcn->AdjustThreshold(dx.data(), dy.data(), size, fixationIndices);
			//peak finding
			SaccadeVector sacc = findPeaks(t, dx.data(), dy.data(), x, y, size);
			return sacc;


		}

		
	protected:
		VelocityThresholdType* thresholdFcn;
		
		SaccadeVector findPeaks(const double* t, const double* x, const double* y, const double* posx, const double* posy, std::size_t size) {
			std::size_t k = 0;
			std::size_t startIdx = 0;
			std::size_t endIdx = 0;
			SaccadeVector vec;
			if (size == 0) return vec;
			while (thresholdFcn->AboveThreshold(x, y, k)) {//ignore the peak if eye was moving at high velocity at the beginning, just move the cursor to a subthreshold position
				k++;
			}
			while (k < size) {
				if (thresholdFcn->AboveThreshold(x, y, k)) {
					
					//look for starting point
					std::size_t nStart = k-1;
					
					while (thresholdFcn->AboveReducedThreshold(x, y, nStart)) {
						nStart--;
					}
					startIdx = nStart + 1;
					
					//look for next subthreshold value, and collect the max value in the interval
					std::size_t sampleMax;
					double maxValue = -1.0;
					std::size_t n = k + 1;
					while ((n < size) && thresholdFcn->AboveReducedThreshold(x,y,n)) {
						double v = std::sqrt(x[n] * x[n] + y[n] * y[n]);
						if (v > maxValue) {
							maxValue = v;
							sampleMax = n;
						}
						n++;
					}
					k = n;
					//if reached the ending, ignore the peak
					if (n == size) continue;
					endIdx = n;
					SaccadeMonocularProperties sacc;
					sacc.duration = t[endIdx] - t[startIdx];
					if (sacc.duration < minDuration) continue;
					sacc.endTimestamp = t[endIdx];
					sacc.startTimestamp = t[startIdx];
					sacc.peakVelocity = maxValue;
					sacc.startX = posx[startIdx];
					sacc.startY = posy[startIdx];
					sacc.endX = posx[endIdx];
					sacc.endY = posy[endIdx];
					sacc.magnitude = std::sqrt((sacc.endX - sacc.startX)*(sacc.endX - sacc.startX) + (sacc.endY - sacc.startY)*(sacc.endY - sacc.startY));
					if (sacc.magnitude > magnitudeMaxThreshold || sacc.magnitude < magnitudeMinThreshold) continue;
					sacc.direction = atan2(sacc.endY - sacc.startY, sacc.endX - sacc.startX) / M_PI*180.0;
					sacc.duration = sacc.endTimestamp - sacc.startTimestamp;
					vec.push_back(sacc);


				}
				k++;
			}

			return vec;

		}
	};
	

	class DllExport NeweyWestRegression {
	public:
		NeweyWestRegression() = default;
		void SetY(const std::vector<double>& y);
		void SetX(const std::vector<std::vector<double> > x);
		void SetMaxLag(unsigned int l);
		void Solve();
		std::vector<double> GetPredictors() const;
		std::vector<double> GetPredictorsSem() const;
		double GetRsquared() const;
	private:
		typedef Eigen::MatrixXd MatrixType;
		typedef Eigen::VectorXd VectorType;
		double _rsquared;
		std::vector<double> _b;
		std::vector<double> _s2;
		unsigned int maxLag;
		MatrixType _X;
		VectorType _y;
	};
}

#undef M_PI

#endif
