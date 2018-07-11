#ifndef REPORTSNUMERICS_H
#define REPORTSNUMERICS_H
#include <cstddef>
#include "Analysis\DataTypes.h"
#include <vector>
#include <algorithm>
#include <numeric>
#ifndef M_PI
#define M_PI 3.14159265359
#endif
namespace Analysis {

	DllExport void firstDerivative(const double* __restrict x, double* __restrict xprime, const long size, double sampleRate = 250);
	DllExport void firstAndSecondDerivatives(const double* __restrict x, double* __restrict xprime, double* __restrict xsecond, const long size, double sampleRate = 250);
	DllExport void polarCoordinates(const double* __restrict x, const double* __restrict y, double* __restrict rho, double* __restrict theta, const long int size);
	DllExport double nanmean(const double* __restrict x, const long int size);
	DllExport double nanvar(const double* __restrict x, const double mean, const long int size);
	DllExport double nanmedian(const double* __restrict x, const long int size);
	DllExport double nanmedvar(const double* __restrict x, const double median, const long int size);
	DllExport double nanmean(const double* __restrict x, const long int size, std::vector<std::size_t>& idx);
	DllExport double nanharmmean(const double* __restrict x, const long int size);
	DllExport double nanvar(const double* __restrict x, const double mean, const long int size, std::vector<std::size_t>& idx);
	DllExport double nanmedian(const double* __restrict x, const long int size, std::vector<std::size_t>& idx);
	DllExport double nanmedvar(const double* __restrict x, const double median, const long int size, std::vector<std::size_t>& idx);
	DllExport double nanmax(const double* __restrict x, const long int size);
	DllExport std::vector<double> HistogramFromEdges(const std::vector<double>& data, const std::vector<double>& edges);
	DllExport std::pair<double, double> eigenCrossCorrelation(std::vector<double>&, std::vector<double>& xCorrInputVecSecond);
	DllExport std::size_t delayEstimation(const std::vector<double>& ref, const std::vector<double> data);
	DllExport double entropy2D(const double*  __restrict x, const double* __restrict y, const long int size, const double minX, const double minY, 
		const double maxX, const double maxY, const unsigned int numBinsX, const unsigned int numBinsY);
	DllExport double GaussianCDF(double x, double mu, double sigma);
	DllExport double GaussianMixtureCDF(double x, std::vector<double> mu, std::vector<double> sigma, std::vector<double> w);
	DllExport double ExGaussianCDF(double x, double mu, double sigma, double tao);
	DllExport double linearSlope(const std::vector<double>& x, const std::vector<double>& y);
	DllExport inline bool NanCompareMin(double a, double b) {
		if (isnan(a) && isnan(b)) return false;
		if (isnan(a)) return false;
		if (isnan(b)) return true;
		return a < b;
	}
	DllExport void SubtractMean(double* x, const long sz);
	DllExport inline bool NanCompareMax(double a, double b) {
		if (isnan(a) && isnan(b)) return false;
		if (isnan(a)) return true;
		if (isnan(b)) return false;
		return a < b;
	}

	struct SinusoidalResidual {
		SinusoidalResidual(double t, double y) : t_(t), y_(y) {};
		template <typename T>
		inline bool operator()(const T* const A, const T* const f, const T* const phase, const T* const C, T* residual)  const {
			residual[0] = T(y_) - (A[0] * sin(2.0 * M_PI * f[0] * t_ + phase[0]) + C[0]);
			return true;
		}

	private:
		const double t_;
		const double y_;
	};

	struct GeneralizedBellFunctionResidual {
		GeneralizedBellFunctionResidual(double x, double y) : x_(x), y_(y) {};
		template <typename T>
		inline bool operator()(const T* a, const T* b, const T* c, const T* base, const T* G, T* residual) const {
			
			T den = 1.0 + ceres::pow(ceres::abs((x_ - *c) / *a), 2.0 * (*b));
			residual[0] = T(y_) - (*base + *G / den);
			return true;
		}
	private:
		const double x_;
		const double y_;
	};

	struct InverseLogitModelResidual {
		InverseLogitModelResidual(unsigned int nComponents, double t, double y) : nComponents_(nComponents), t_(t), y_(y), k_(4.6), noOverlapConstraint(false) {};
		template <typename T>
		inline bool operator()(const T* const a, const T* const Tp, const T* const D, const T* const b, T* residual) const {
			if (noOverlapConstraint) {
				for (std::size_t k = 1; k < nComponents_; k++) {
					if ((Tp[k] - Tp[k - 1]) <= k_*(D[k] + D[k - 1])) return false; //orthogonality conditions, using the trick described in ceres documentation
				}
			}
			T value = b[0];
			for (std::size_t k = 0; k < nComponents_; k++) {
				value += a[k] * 1.0 / (1.0 + ceres::exp(-(t_ - Tp[k]) / D[k]));//need to use ceres::exp instead of std::exp or the compiler chokes on templates arguments
			}

			residual[0] = T(y_) - value;
			return true;
		}

		static inline void CalculateModelValues(const std::vector<double>&a, const std::vector<double>&T, const std::vector<double>& D,
			double b, const std::vector<double>& t, std::vector<double>& y) {
			y.resize(t.size());
			for (std::size_t i = 0; i < t.size(); i++) {
				y[i] = b;
				for (std::size_t k = 0; k < a.size(); k++) {
					y[i] += a[k] * 1.0 / (1.0 + std::exp(-(t[i] - T[k]) / D[k]));
				}
			}
		}

		static inline void CalculateModelValues(const std::vector<double>&a, const std::vector<double>&T, const std::vector<double>& D,
			double b, double t, double& y) {


			y = b;
			for (std::size_t k = 0; k < a.size(); k++) {
				y += a[k] * 1.0 / (1.0 + std::exp(-(t - T[k]) / D[k]));
			}

		}

		static inline void CalculateModelValues(double a, double T, double  D, double b, const std::vector<double>& t, std::vector<double>& y) {
			y.resize(t.size());
			for (std::size_t i = 0; i < t.size(); i++) {
				y[i] = b + a * 1.0 / (1.0 + std::exp(-(t[i] - T) / D));
			}
		}

		static inline void CalculateModelValues(double a, double T, double  D, double b, double t, double& y) {


			y = b + a * 1.0 / (1.0 + std::exp(-(t - T) / D));

		}
	private:
		const double t_;
		const double y_;
		const double k_;
		const bool noOverlapConstraint;
		unsigned int nComponents_;

	};

	struct InverseLogitModelMinimizer {
		InverseLogitModelMinimizer(std::vector<double>& a, std::vector<double>& T, std::vector<double>& D, double b) : a_(a), T_(T), D_(D), b_(b) {}

		template <typename T>
		inline bool operator()(const T* const t, T* residual) const {
			T value(b_);
			for (std::size_t k = 0; k < a_.size(); k++) {
				value += a_[k] / (1.0 + ceres::exp(-(t[0] - T_[k]) / D_[k]));
			}
			residual[0] = value;
			return true;
		}

	private:
		std::vector<double>& a_;
		std::vector<double>& T_;
		std::vector<double>& D_;
		double b_;
	};


	class LinearFilter {
	public:
		LinearFilter();
		virtual void SetNumerator(const std::vector<double>& B);
		virtual void SetDenominator(const std::vector<double>& A);
		virtual void SetInitialConditions(const std::vector<double>& Z);
		virtual void ResetFilter();
		virtual void operator()(std::vector<double>& X, std::vector<double>& Y) = 0;
	protected:
		std::vector<double> A, B, Z;

	};

	class ZeroPhaseIIRFilter : public LinearFilter {
	public:
		ZeroPhaseIIRFilter();
		virtual void operator()(std::vector<double>& X, std::vector<double>& Y);
	private:
		typedef std::vector<double> vectord;
		typedef std::vector<int> vectori;
		void filter(const vectord &X, vectord &Y);
		void add_index_range(vectori &indices, int beg, int end, int inc = 1);
		void add_index_const(vectori &indices, int value, size_t numel);
		void append_vector(vectord &vec, const vectord &tail);
		vectord subvector_reverse(const vectord &vec, int idx_end, int idx_start);
		inline int max_val(const vectori& vec) {
			return std::max_element(vec.begin(), vec.end())[0];
		}

	};

	


}

#endif
