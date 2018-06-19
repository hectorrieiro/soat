#define EXPORTING

#include <math.h>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>

#include <iostream>
#include "Analysis/Algorithms.h"
#include "Analysis/Numerics.h"

#ifndef M_PI
#define M_PI 3.14159265
#endif

#define COMPARE(val, comp) std::isnan((comp)) ? std::isnan((val)) : (val) == (comp)


namespace Analysis {

	

	void BlinkRemover::SetInputBlinkValue(double val) {
		inputBlinkValue = val;
	}
	void BlinkRemover::SetOutputBlinkValue(double val) {
		outputBlinkValue = val;
	}
	void BlinkRemover::SetBufferSamples(unsigned int n) {
		bufferSamples = n;
	}
	void BlinkRemover::SetMinimumBlinkDuration(unsigned int samples) {
		minimumBlinkDuration = samples;
	}
	void BlinkRemover::SetInterpolateMissingSamples(bool yesNo) {
		interpolateMissingSamples = yesNo;
	}

	unsigned int BlinkRemover::operator()(double* x, double*y, std::size_t size) { 
		
		unsigned int numBlinksRemoved = 0;

		std::size_t k = 0;
		while (k < size) {
			if (COMPARE(x[k], inputBlinkValue) || COMPARE(y[k], inputBlinkValue)) {
				std::size_t n = k + 1;
				while (n < size && (COMPARE(x[n], inputBlinkValue) || COMPARE(y[n], inputBlinkValue))) {
					n++;
				}
				if (n == size || k == 0 || (n - k + 1) >= minimumBlinkDuration) {
					std::size_t minIdx = std::min(k - bufferSamples, std::numeric_limits<std::size_t>::max());
					std::size_t maxIdx = std::min(n + bufferSamples, size - 1);
					for (std::size_t i = minIdx; i <= maxIdx; i++) {
						x[i] = outputBlinkValue;
						y[i] = outputBlinkValue;
					}
					k = n + bufferSamples + 1;
					numBlinksRemoved++;
					continue;
				}
				else if (n != size && interpolateMissingSamples) {
					double xStartValue = x[k - 1];
					double yStartValue = y[k - 1];
					double xEndValue = x[n];
					double yEndValue = y[n];
					double xStep = (xEndValue - xStartValue) / double(n - k + 1);
					double yStep = (yEndValue - yStartValue) / double(n - k + 1);
					for (std::size_t i = k; i < n; i++) {
						x[i] = xStep * (i - k + 1);
						y[i] = yStep * (i - k + 1);
					}
					k = n + 1;
				}
				
			}
			else {
				k++;
			}
		}
		return numBlinksRemoved;
	}


	void HardThresholdType::SetThreshold(double t) {
		this->_t = t;
	};

	void HardThresholdType::AdjustThreshold(const double* vx, const double *vy, std::size_t size, std::vector<std::size_t> fixationIndices) {};
	
	bool HardThresholdType::AboveThreshold(const double *vx, const double* vy, std::size_t k) { 
		
		
		return _t < std::sqrt(vx[k]*vx[k]+vy[k]*vy[k]); 
	}

	bool HardThresholdType::AboveReducedThreshold(const double *vx, const double* vy, std::size_t k) {


		return (_t/3.0) < std::sqrt(vx[k] * vx[k] + vy[k] * vy[k]);
	}
	void EngbertKlieglThreshold::SetLambda(double _l) {
		this->l = _l;
	}

	void EngbertKlieglThreshold::AdjustThreshold(const double* vx, const double *vy, std::size_t size, std::vector<std::size_t> fixationIndices) {
		double mx = nanmedian(vx, size, fixationIndices);
		double my = nanmedian(vy, size, fixationIndices);
		double stdX = std::sqrt(nanmedvar(vx, mx, size, fixationIndices));
		double stdY = std::sqrt(nanmedvar(vy, my, size, fixationIndices));
		this->thresholdX = abs(mx) + l*stdX;
		this->thresholdY = abs(my) + l*stdY;
	}
	bool EngbertKlieglThreshold::AboveThreshold(const double* vx, const double* vy, std::size_t k) {
		
		

		return (pow(vx[k] / thresholdX, 2.0) + pow(vy[k] / thresholdY, 2.0)) > 1.0;
	
	};

	bool EngbertKlieglThreshold::AboveReducedThreshold(const double* vx, const double* vy, std::size_t k) {



		return (pow(vx[k] / thresholdX, 2.0) + pow(vy[k] / thresholdY, 2.0)) > 1.0;

	};

	void SaccadeDetectionAlgorithm::DetectOvershoots(bool b) {
		_detectOvershoots = b;
	};

	
	void SaccadeDetectionAlgorithm::SetMinISI(double isi) {
		this->minISI = isi;
	};
	
	void SaccadeDetectionAlgorithm::FilterBinocular(bool yesNo) {
		this->doFilterBinocular = yesNo;
	};

	void SaccadeDetectionAlgorithm::SetMinimumSaccadeDuration(double d) {
		this->minDuration = d;
	}

	void SaccadeDetectionAlgorithm::SetMaximumSaccadeSize(double d) {
		this->magnitudeMaxThreshold = d;
	}

	void SaccadeDetectionAlgorithm::RemoveMonoculars(SaccadeVector& left, SaccadeVector& right) {

		SaccadeVector leftBinocc;
		SaccadeVector rightBinocc;

		std::size_t startingRightIndex = 0;
		for (std::size_t k = 0; k < left.size(); k++) {
			for (std::size_t j = startingRightIndex; j < right.size(); j++) {
				if (TimingOverlap(left[k], right[j])) {
					leftBinocc.push_back(left[k]);
					startingRightIndex = j + 1;
					break;
				}
				if (right[j].startTimestamp > left[k].endTimestamp) {
					break;
				}
			}
		}

		std::size_t startingLeftIndex = 0;
		for (std::size_t k = 0; k < right.size(); k++) {
			for (std::size_t j = startingLeftIndex; j < left.size(); j++) {
				if (TimingOverlap(right[k], left[j])) {
					rightBinocc.push_back(right[k]);
					startingLeftIndex = j + 1;
					break;
				}
				if (left[j].startTimestamp > right[k].endTimestamp) {
					break;
				}
			}
		}

		left = leftBinocc;
		right = rightBinocc;

	}

	void SaccadeDetectionAlgorithm::RemoveOvershoots(SaccadeVector& v) {
		std::vector<SaccadeMonocularProperties> vo;
		if (v.empty()) return;
		vo.push_back(v[0]);
		for (std::size_t k = 1; k < v.size(); k++) {
			if (v[k].startTimestamp - v[k - 1].endTimestamp > 20) {
				vo.push_back(v[k]);
			}
		}

		v = vo;

	}

	void SaccadeDetectionAlgorithm::SetTaskType(TaskType t) {
		this->taskType = t;
	}

	void SaccadeDetectionAlgorithm::RemoveBinocularOvershoots(SaccadeVector& v1, SaccadeVector& v2) {
		std::vector<SaccadeMonocularProperties> v1n, v2n;
		if (v1.empty()) return;
		v1n.push_back(v1[0]);
		v2n.push_back(v2[0]);
		for (std::size_t k = 1; k < v1.size(); k++) { 
			bool isOvershoot = (v1[k].startTimestamp - v1[k - 1].endTimestamp < 20) || (v2[k].startTimestamp - v2[k - 1].endTimestamp < 20);
			if (!isOvershoot)  {
				v1n.push_back(v1[k]);
				v2n.push_back(v2[k]);
			}
		}

		v1 = v1n;
		v2 = v2n;
	}


	int SaccadeDetectionAlgorithm::TimingOverlap(const SaccadeMonocularProperties& s1, const SaccadeMonocularProperties& s2) {

		if (s1.startTimestamp <= s2.startTimestamp && s1.endTimestamp >= s2.startTimestamp) {
			return 1;
		}
		else if (s1.startTimestamp <= s2.endTimestamp && s1.endTimestamp >= s2.endTimestamp) {
			return 1;
		}
		else if (s1.startTimestamp >= s2.startTimestamp && s1.endTimestamp <= s2.endTimestamp) {
			return 1;
		}
		else {
			return 0;
		}

	}
	
	void SaccadeDetectionAlgorithm::SetThresholds(double min, double max) {

		this->magnitudeMaxThreshold = max;
		this->magnitudeMinThreshold = min;

	}

	DllExport SWJVector SWJDetection(SaccadeVector& vecLeft, SaccadeVector& vecRight) {
		SWJVector swjs(0);
		if (vecLeft.size() == 0) return swjs;
		std::vector<double> scores(vecLeft.size() - 1);
		const std::vector<double> Dmu = { 180.0, 180.0 };
		const std::vector<double> Dsigma = { 30.0, 7.0 };
		const std::vector<double> Dw = { 0.4, 0.6 };
		const std::vector<double> Mmu = { 0.0, 0.0 };
		const std::vector<double> Msigma = { 0.39, 0.16 };
		const std::vector<double> Mw = { 0.4, 0.6 };
		

		for (std::size_t k = 0; k < vecLeft.size() - 1; k++) {
			double isi = vecLeft[k + 1].startTimestamp - vecLeft[k].endTimestamp;
			double dMag = (vecLeft[k + 1].magnitude - vecLeft[k].magnitude) / (vecLeft[k + 1].magnitude + vecLeft[k].magnitude);
			double dDir = vecLeft[k + 1].direction - vecLeft[k].direction;

			if (std::cos(dDir / 180.0*M_PI) > 0) {
				scores[k] = 0;
				continue;
			}
			while (dDir < 0.0) {
				dDir += 360.0;
			}
			double FD = GaussianMixtureCDF(dDir, Dmu, Dsigma, Dw);
			double fd = (dDir > 180.0) ? 1.0 - FD : FD;
			double FM = GaussianMixtureCDF(std::abs(dMag), Mmu, Msigma, Mw);
			double fm = 1.0 - FM;
			double FI = ExGaussianCDF(isi, 120.0, 60.0, 180.0);
			double fi = (isi >= 200.0) ? 1 - FI : FI;
			scores[k] = fd*fi*fm;

			if (scores[k] >= 0.0014) {
				SWJProperties props;
				props.eye = SWJ_BOTH;
				props.swjScore = scores[k];
				props.dirDifference = dDir;
				props.firstDirection = vecLeft[k].direction;
				props.firstMagnitude = vecLeft[k].magnitude;
				props.firstSaccadeIndex = k;
				props.isi = isi;
				props.secondDirection = vecLeft[k + 1].direction;
				props.secondMagnitude = vecLeft[k + 1].magnitude;
				props.secondSaccadeIndex = k + 1;

				swjs.push_back(props);

				vecLeft[k].inSWJ = true;
				vecRight[k].inSWJ = true;
				vecLeft[k + 1].inSWJ = true;
				vecRight[k + 1].inSWJ = true;
				vecRight[k].swjIndex = swjs.size() - 1;
				vecRight[k].swjIndex = swjs.size() - 1;


			}


		}
		return swjs;

	}


	void NeweyWestRegression::SetY(const std::vector<double>& y) {
		_y = MatrixType::Zero(y.size(), 1);
		for (std::size_t i = 0; i < y.size(); i++) {
			_y[i] = y[i];
		}

	}
	void NeweyWestRegression::SetX(const std::vector<std::vector<double> > x) {
		_X = MatrixType::Zero(x[0].size(), x.size());
		for (std::size_t m = 0; m < x.size(); m++) {
			for (std::size_t n = 0; n < x[m].size(); n++) {
				_X(n,m) = x[m][n];
			}
		}
	}
	void NeweyWestRegression::SetMaxLag(unsigned int l) {
		this->maxLag = l;
	}

	void NeweyWestRegression::Solve() {
		//solve least-squares problem
		auto solution = _X.colPivHouseholderQr().solve(_y).eval();
		
		_b.resize(solution.rows());
		VectorType v = solution;
		for (std::size_t n = 0; n < solution.rows(); n++) {
			_b[n] = v(n);
		}
		//residual calculation
		auto e = _y - _X*solution;

		//determination coefficient
		_rsquared = 1.0 - e.squaredNorm() / _y.squaredNorm();

		//standard error calculation
		auto XtX = (_X.transpose()*_X).eval();
		MatrixType ohm = MatrixType::Zero(_X.cols(), _X.cols());
		
		for (std::size_t i = 0; i < _X.rows(); i++) {
			auto C = _X.row(i).transpose()*_X.row(i);
			double e2 = e[i] * e[i];
			auto A = e2*C;
			ohm += A;
		}

		ohm *= double(_X.rows()) / double(_X.rows() - _X.cols());
		MatrixType ohmK = MatrixType::Zero(_X.cols(), _X.cols());
		for (int l = 1; l <= maxLag; l++) {
			double lambda = 1.0 - double(l) / double(maxLag + 1);
			MatrixType A = MatrixType::Zero(_X.cols(), _X.cols());
			for (int t = l+1; t < _X.rows(); t++) {
				auto XX = _X.row(t).transpose()*_X.row(t - l);
				A = e[t] * e[t - l] * ( XX - XX.transpose());
			}
			ohmK += lambda*A;
		}

		ohm += double(_X.rows()) / double(_X.rows() - _X.cols())*ohmK;
		auto S = XtX.inverse()*ohm*XtX.inverse(); //this works because I am assuming small number of estimator (comlumns of _X)
		_s2.resize(_X.cols());
		for (std::size_t i = 0; i < _X.cols(); i++) {
			_s2[i] = std::sqrt(std::abs(S(i,i)));
		}

	}
	std::vector<double> NeweyWestRegression::GetPredictors() const {
		return _b;
	}
	std::vector<double> NeweyWestRegression::GetPredictorsSem() const {
		return _s2;
	}
	double NeweyWestRegression::GetRsquared() const {
		return _rsquared;
	}
	
}