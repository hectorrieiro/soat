#define EXPORTING
#include "Analysis/Numerics.h"
#include <math.h>
#include <algorithm>
#include "Eigen/Dense"
#include "unsupported\Eigen\FFT"
#include <fstream>

namespace Analysis {
	void firstDerivative(const double* __restrict x, double* __restrict xprime, const long size, double sampleRate) {
		const double coeff1 = 126.0/1188.0;
		const double coeff2 = 193.0 / 1188.0;
		const double coeff3 = 142.0 / 1188.0;
		const double coeff4 = -86.0 / 1188.0;
		switch (size) {
		case 0:
			break;
		case 1:
			xprime[0] = 0.0;
			break;
		case 2:
			xprime[0] = 0.0;
			xprime[1] = (x[1] - x[2])*sampleRate;
			break;
		case 3:
			xprime[0] = 0.5*x[1]*sampleRate;
			xprime[1] = (-0.5*x[0] + 0.5*x[2])*sampleRate;
			xprime[2] = -0.5*x[1] * sampleRate;
			break;
		default:
			xprime[0] = 0.0;
			xprime[1] = (-coeff4*x[3]-coeff3*x[2]-coeff2*x[1]-coeff1*x[0] + coeff1*x[2] + coeff2*x[3] + coeff3*x[4] + coeff4*x[5])*sampleRate;
			xprime[2] = (-coeff4*x[2] - coeff3*x[1] - coeff2*x[0] - coeff1*x[1] + coeff1*x[3] + coeff2*x[4] + coeff3*x[5] + coeff4*x[6])*sampleRate;
			xprime[3] = (-coeff4*x[1] - coeff3*x[0] - coeff2*x[1] - coeff1*x[2] + coeff1*x[4] + coeff2*x[5] + coeff3*x[6] + coeff4*x[7])*sampleRate;
#pragma loop(hint_parallel(0))
#pragma loop(ivdep) 
			for (long k = 4; k < size - 4; ++k) {
				xprime[k] = -coeff4*x[k - 4];
				xprime[k] += -coeff3*x[k - 3];
				xprime[k] += -coeff2*x[k - 2];
				xprime[k] += -coeff1*x[k - 1];
				xprime[k] += coeff1*x[k + 1];
				xprime[k] += coeff2*x[k + 2];
				xprime[k] += coeff3*x[k + 3];
				xprime[k] += coeff4*x[k + 4];
				xprime[k] *= sampleRate;
			}
			xprime[size - 4] = (-coeff4*x[size - 8] - coeff3*x[size - 7] - coeff2*x[size - 6] - coeff1*x[size - 5] + coeff1*x[size - 3] + coeff2*x[size - 2] + coeff3*x[size - 1] + coeff4*x[size - 1])*sampleRate;
			xprime[size - 3] = (-coeff4*x[size - 7] - coeff3*x[size - 6] - coeff2*x[size - 5] - coeff1*x[size - 4] + coeff1*x[size - 2] + coeff2*x[size - 1] + coeff3*x[size - 1] + coeff4*x[size - 2])*sampleRate;
			xprime[size - 2] = (-coeff4*x[size - 6] - coeff3*x[size - 5] - coeff2*x[size - 4] - coeff1*x[size - 3] + coeff1*x[size - 1] + coeff2*x[size - 2] + coeff3*x[size - 3] + coeff4*x[size - 4])*sampleRate;
			xprime[size - 1] = 0.0;

		}

	}

	DllExport void firstAndSecondDerivatives(const double* __restrict x, double* __restrict xprime, double* __restrict xsecond, const long size, double sampleRate) {
		const double coeff11 = 2.0 / 3.0;
		const double coeff12 = 1.0 / 12.0;
		const double coeff20 = -5.0 / 2.0;
		const double coeff21 = 4.0 / 3.0;
		const double coeff22 = -1.0 / 12.0;
		switch (size) {
		case 0:
			break;
		case 1:
			xprime[0] = 0.0;
			xsecond[0] = 0.0;
			break;
		case 2:
			xprime[0] = 0.0;
			xprime[1] = (x[1] - x[2])*sampleRate;
			xsecond[0] = 0.0;
			xsecond[1] = 0.0;
			break;
		case 3:
			xprime[0] = 0.5*x[1] * sampleRate;
			xprime[1] = (-0.5*x[0] + 0.5*x[2])*sampleRate;
			xprime[2] = (-0.5*x[1])*sampleRate;
			xsecond[0] = (-2.0*x[0] + x[1])*sampleRate;
			xsecond[1] = (x[0] - 2.0*x[1] + x[2])*sampleRate;
			xsecond[2] = (x[1] - 2.0*x[2])*sampleRate;
			break;
		default:
			xprime[0] = (coeff11*x[1] - coeff12*x[2])*sampleRate;
			xprime[1] = (-coeff11*x[0] + coeff11*x[2] - coeff12*x[3])*sampleRate;
			xsecond[0] = (coeff20*x[0] + coeff21*x[1] + coeff22*x[2])*sampleRate;
			xsecond[1] = (coeff21*x[0] + coeff20*x[1] + coeff21*x[2] + coeff22*x[3])*sampleRate;
#pragma loop(hint_parallel(0))
#pragma loop(ivdep) 
			for (long k = 2; k < size - 2; ++k) {
				xprime[k] = (coeff12*x[k - 2] - coeff11*x[k - 1] + coeff11*x[k + 1] - coeff12*x[k + 2])*sampleRate;
				xsecond[k] = (coeff22*x[k - 2] + coeff21*x[k - 1] + coeff20*x[k] + coeff21*x[k + 1] + coeff22*x[k + 2])*sampleRate;
			}
			xprime[size - 2] = (coeff12*x[size - 4] - coeff11*x[size - 3] + coeff11*x[size - 1])*sampleRate;
			xprime[size - 1] = (coeff12*x[size - 3] - coeff11*x[size - 2])*sampleRate;
			xsecond[size - 2] = (coeff22*x[size - 4] + coeff21*x[size - 3] + coeff20*x[size - 2] + coeff21*x[size - 1])*sampleRate;
			xsecond[size - 1] = (coeff22*x[size - 3] + coeff21*x[size - 2] + coeff20*x[size - 1])*sampleRate;

		}

	}
	DllExport void polarCoordinates(const double* __restrict x, const double* __restrict y, double* __restrict rho, double* __restrict theta, const long int size) {
#pragma loop(hint_parallel(0))

		for (long k = 0; k < size; ++k) {
			rho[k] = sqrt(x[k] * x[k] + y[k] * y[k]);
			theta[k] = atan2(y[k], x[k]);
		}
	}

	double nanmean(const double* __restrict x, const long int size, std::vector<std::size_t>& idx) {
		if (idx.empty()) return nanmean(x, size);

		std::size_t n = 0;
		double acc = 0.0;
#pragma loop(hint_parallel(0))
#pragma loop(ivdep) 
		for (std::size_t k = 0; k < idx.size(); k++) {
			if (!isnan(x[idx[k]])) {
				acc += x[idx[k]];
				n++;
			}
		}

		return acc / double(n);

	}
	double nanvar(const double* __restrict x, const double mean, const long int size, std::vector<std::size_t>& idx) {
		if (idx.empty()) return nanvar(x, mean, size);

		std::size_t n = 0;
		double acc = 0.0;
#pragma loop(hint_parallel(0))
#pragma loop(ivdep) 
		for (std::size_t k = 0; k < idx.size(); k++) {
			if (!isnan(x[idx[k]])) {
				acc += (x[idx[k]] - mean)*(x[idx[k]] - mean);
				n++;
			}
		}

		return acc / double(n);
	}
	double nanmedian(const double* __restrict x, const long int size, std::vector<std::size_t>& idx) {
		if (idx.empty()) return nanmedian(x, size);
		std::vector<double> xf(idx.size());
		for (std::size_t k = 0; k < idx.size(); k++) {
			xf[k] = x[idx[k]];
		}

		return nanmedian(xf.data(), idx.size());

	}
	double nanmedvar(const double* __restrict x, const double median, const long int size, std::vector<std::size_t>& idx) {
		if (idx.empty()) return nanmedvar(x, median, size);
		std::vector<double> xf(idx.size());
		for (std::size_t k = 0; k < idx.size(); k++) {
			xf[k] = x[idx[k]];
		}

		return nanmedvar(xf.data(), median, idx.size());
	}

	DllExport double nanmean(const double* __restrict x, const long int size) {
		std::size_t n = 0;
		double acc = 0.0;
#pragma loop(hint_parallel(0))
#pragma loop(ivdep) 
		for (std::size_t k = 0; k < size; k++) {
			if (!isnan(x[k])) {
				acc += x[k];
				n++;
			}
		}

		return acc / double(n);

	}

	DllExport double nanharmmean(const double* __restrict x, const long int size) {
		std::size_t n = 0;
		double acc = 0.0;
#pragma loop(hint_parallel(0))
#pragma loop(ivdep) 
		for (std::size_t k = 0; k < size; k++) {
			if (!isnan(x[k]) && std::abs(x[k]) > 0.000001) {
				acc += 1.0 / x[k];
				n++;
			}
		}

		return double(n) / acc;
	}

	DllExport double nanvar(const double* __restrict x, const double mean, const long int size) {
		std::size_t n = 0;
		double acc = 0.0;
#pragma loop(hint_parallel(0))
#pragma loop(ivdep) 
		for (std::size_t k = 0; k < size; k++) {
			if (!isnan(x[k])) {
				acc += (x[k]-mean)*(x[k]-mean);
				n++;
			}
		}

		return acc / double(n);

	}

	DllExport double nanmedian(const double* __restrict x, const long int size) {

		std::vector<double> v(size);
		std::copy(x, x + size, v.begin());

		auto it = std::partition(v.begin(), v.end(), isnan<double>);
		std::sort(it, v.end());

		int numEl = v.end()-it;
		double med;
		if (numEl == 0) return 0.0;
		if (numEl == 1) return *it;
		if (numEl == 2) return (*it + *(it + 1)) / 2.0;
		if (numEl % 2 == 0) {
			med = (*(it + numEl / 2) + *(it + numEl / 2 + 1)) / 2.0;
		}
		else {
			med = *(it + std::size_t(numEl / 2.0));
		}

		return med;
	}
	double nanmax(const double* __restrict x, const long int size) {
		std::vector<double> v(size);
		std::copy(x, x + size, v.begin());

		auto it = std::partition(v.begin(), v.end(), isnan<double>);
		

		int numEl = v.end() - it;
		return *(std::max_element(it, v.end()));
	}


	DllExport double nanmedvar(const double* __restrict x, const double median, const long int size) {

		std::vector<double> v(size);
		for (std::size_t k = 0; k < size; k++) {
			v[k] = x[k] - median;
			v[k] *= v[k];
		}

		auto it = std::partition(v.begin(), v.end(), isnan<double>);
		std::sort(it, v.end());

		int numEl = v.end() - it;
		double med;
		if (numEl == 0) return 0.0;
		if (numEl == 1) return *it;
		if (numEl == 2) return (*it + *(it + 1)) / 2.0;
		if (numEl % 2 == 0) {
			med = (*(it + numEl / 2) + *(it+numEl / 2 + 1)) / 2.0;
		}
		else {
			med = *(it + std::size_t(numEl / 2.0));
		}

		return med;

	}
	std::vector<double> HistogramFromEdges(const std::vector<double>& data, const std::vector<double>& edges) {
		std::vector<double> out(edges.size() - 1);
		for (std::size_t k = 0; k < data.size(); k++) {
			for (std::size_t j = 0; j < edges.size() - 1; j++) {
				if (data[k] >= edges[j] && data[k] < edges[j + 1]) {
					out[j]++;
					break;
				}
			}
		}

		return out;
	}

	std::pair<double, double> eigenCrossCorrelation(std::vector<double>& xCorrInputVecFirst
		, std::vector<double>& xCorrInputVecSecond)
	{
		Eigen::FFT<double> fft;
		int N = std::max(xCorrInputVecFirst.size(), xCorrInputVecSecond.size());

		//Compute the FFT size as the "next power of 2" of the input vector's length (max)
		int b = ceil(log2(2.0 * N - 1));
		int fftsize = pow(2, b);
		int end = fftsize - 1;
		int maxlag = N - 1;
		size_t firstSize = xCorrInputVecFirst.size();
		size_t secondSize = xCorrInputVecSecond.size();

		//Zero Padd
		for (int i = xCorrInputVecFirst.size(); i < fftsize; i++)
		{
			xCorrInputVecFirst.push_back(0);
		}

		for (int i = xCorrInputVecSecond.size(); i < fftsize; i++)
		{
			xCorrInputVecSecond.push_back(0);
		}

		std::vector<std::complex<double> > freqvec;
		std::vector<std::complex<double> > freqvec2;

		//FFT for freq domain to both vectors
		fft.fwd(freqvec, xCorrInputVecFirst);
		fft.fwd(freqvec2, xCorrInputVecSecond);

		//Main step of cross corr
		for (int i = 0; i < fftsize; i++)
		{
			freqvec[i] = freqvec[i] * std::conj(freqvec2[i]);
		}

		std::vector<double> result;
		fft.inv(result, freqvec);

		//Will get rid of extra zero padding and move minus lags to beginning without copy
		std::vector<double> result2(std::make_move_iterator(result.begin() + end - maxlag + 1),
			std::make_move_iterator(result.end()));

		result2.insert(result2.end(), make_move_iterator(result.begin())
			, make_move_iterator(result.begin() + maxlag));


		auto minMaxRange = std::minmax_element(result2.begin(),
			result2.end());

		//Will take back the changes which made in input vector
		if (xCorrInputVecFirst.size() != firstSize)
			xCorrInputVecFirst.resize(firstSize);
		if (xCorrInputVecSecond.size() != secondSize)
			xCorrInputVecSecond.resize(secondSize);


		//Return val
		auto resultIndex = ((minMaxRange.second - result2.begin()) - N + 1);
		//std::cout << "Max element at: " << resultIndex << std::endl;
		auto maxValue = result[minMaxRange.second - result.begin()];
		return std::make_pair(resultIndex, maxValue);
	}

	std::size_t delayEstimation(const std::vector<double>& ref, const std::vector<double> data) {
		
		std::vector<double> normRef(ref.size());
		double refMean = nanmean(ref.data(), ref.size());
		double refStd = std::sqrt(nanvar(ref.data(), refMean, ref.size()));
		std::transform(ref.begin(), ref.end(), normRef.begin(), [refMean, refStd](double x) -> double {return (x - refMean) / refStd; });
		std::vector<double> normData(data.size());
		double dataMean = nanmean(data.data(), data.size());
		double dataStd = std::sqrt(nanvar(data.data(), dataMean, data.size()));
		std::transform(data.begin(), data.end(), normData.begin(), [dataMean, dataStd](double x) -> double {return (x - dataMean) / dataStd; });
		//normData.erase(normData.begin(), normData.begin() + 1000);
		//normRef.erase(normRef.begin(), normRef.begin() + 1000);
		std::vector<double> xcorr(2*normRef.size()+1);
		for (long k = 0; k < normRef.size(); k++) {
			std::size_t numSamples = k;
			double acc = 0.0;
			for (long m = 0; m < normData.size(); m++) {
				if ((m + k) < normRef.size()) {
					if (!isnan(normRef[m + k]) && !isnan(normData[m])) {
						acc += normRef[m + k] * normData[m];
						numSamples++;
					}
				}
			}
			xcorr[k] = acc / double(numSamples);
		}
		{
			double acc = 0.0;
			std::size_t numSamples = 0;
			for (long m = 0; m < normData.size(); m++) {

				if (!isnan(normRef[m]) && !isnan(normData[m])) {
					acc += normRef[m] * normData[m];
					numSamples++;
				}

			}
			xcorr[normRef.size()] = acc / double(numSamples);
		}
		for (long k = normRef.size() + 1; k < xcorr.size(); k++) {
			std::size_t numSamples = k - normRef.size();
			double acc = 0.0;
			for (long m = 0; m < normData.size(); m++) {
				if ((m - k + long(normRef.size())) >= 0) {
					if (!isnan(normRef[m - k + long(normRef.size())]) && !isnan(normData[m])) {
						acc += normRef[m - k + long(normRef.size())] * normData[m];
						numSamples++;
					}
				}
			}
			xcorr[k] = acc / double(numSamples);

		}

		std::ofstream of("c://down//xcorr.csv", std::ios::trunc);
		for (std::size_t k = 0; k < normRef.size(); k++) {
			of << ref[k] << "," << data[k] << "," << xcorr[k] << std::endl;
		}
		of.close();

		return std::min_element(xcorr.begin(), xcorr.end()) - (xcorr.begin());
	}

	DllExport double entropy2D(const double*  __restrict x, const double* __restrict y, const long int size, double minX, double minY,
		double maxX, double maxY, const unsigned int numBinsX, const unsigned int numBinsY) {

		double binWidthX = (maxX - minX) / numBinsX;
		double binWidthY = (maxY - minY) / numBinsY;

		unsigned int count = 0;
		std::vector<std::vector<unsigned int>> nB(numBinsX);
		for (std::size_t k = 0; k < numBinsX; k++) {
			nB[k].resize(numBinsY, 0);
		}

		//create 2d histogram
		for (std::size_t k = 0; k < size; k++) {
			unsigned int i, j;
			if (x[k] < minX || x[k] > maxX || y[k] < minY || y[k] > maxY || isnan(x[k]) || isnan(y[k])) continue;
			count++;
			for (i = 0; i < numBinsX - 1; i++) {
				if (x[k] >= minX + (double)i*binWidthX && x[k] < minX + (double)(i + 1)*binWidthX) {
					break;
				}
			}
			for (j = 0; j < numBinsY - 1; j++) {
				if (y[k] >= minY + (double)j*binWidthY && y[k] < minY + (double)(j + 1)*binWidthY) {
					break;
				}
			}
			nB[i][j]++;
		}

		//estimate entropy
		double Hxy = 0.0;
		for (std::size_t i = 0; i < numBinsX; i++) {
			for (std::size_t j = 0; j < numBinsY; j++) {
				if (nB[i][j] == 0) continue;
				double p = double(nB[i][j]) / double(count);
				Hxy += -p*std::log2(p);
			}
		}

	
		return Hxy;
		

		
	}

	double linearSlope(const std::vector<double>& x, const std::vector<double>& y) {
		const auto n = x.size();
		const auto s_x = std::accumulate(x.begin(), x.end(), 0.0);
		const auto s_y = std::accumulate(y.begin(), y.end(), 0.0);
		const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
		const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
		const auto a = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
		return a;
	}

	double GaussianCDF(double x, double mu, double sigma) {

		return 0.5*(1 + std::erf((x - mu) / (sigma*std::sqrt(2.0))));

	}
	double GaussianMixtureCDF(double x, std::vector<double> mu, std::vector<double> sigma, std::vector<double> w) {
		double acc = 0.0;
 
		for (std::size_t k = 0; k < w.size(); k++) {
			acc += w[k] * GaussianCDF(x, mu[k], sigma[k]);
		}
		return acc;
	}

	double ExGaussianCDF(double x, double mu, double sigma, double tao) {
		double u = tao*(x - mu);
		double v = tao*sigma;
		double cdf = GaussianCDF(u, 0.0, v);

		cdf = cdf - std::exp(-u + v*v/2.0 + std::log(GaussianCDF(u, v*v, v)));
		return cdf;
	}


	LinearFilter::LinearFilter() {}
	void LinearFilter::SetNumerator(const std::vector<double>& B) {
		this->B = B;
		ResetFilter();
	}
	void LinearFilter::SetDenominator(const std::vector<double>& A) {
		this->A = A;
		ResetFilter();
	}
	void LinearFilter::SetInitialConditions(const std::vector<double>& Z) {
		if (Z.size() != this->Z.size()) return; //exception?
		this->Z.assign(Z.begin(), Z.end());
	}
	void LinearFilter::ResetFilter() {
		Z.resize(std::max(A.size(), B.size()));
		Z.assign(Z.size(), 0.0);
	}
	
	ZeroPhaseIIRFilter::ZeroPhaseIIRFilter() {}
	void ZeroPhaseIIRFilter::operator()(std::vector<double>& X, std::vector<double>& Y) {
		using namespace Eigen;

		int len = X.size();     // length of input
		int na = A.size();
		int nb = B.size();
		int nfilt = (nb > na) ? nb : na;
		int nfact = 3 * (nfilt - 1); // length of edge transients

		if (len <= nfact)
		{
			throw std::domain_error("Input data too short! Data must have length more than 3 times filter order.");
		}

		// set up filter's initial conditions to remove DC offset problems at the
		// beginning and end of the sequence
		B.resize(nfilt, 0);
		A.resize(nfilt, 0);

		vectori rows, cols;
		//rows = [1:nfilt-1           2:nfilt-1             1:nfilt-2];
		add_index_range(rows, 0, nfilt - 2);
		if (nfilt > 2)
		{
			add_index_range(rows, 1, nfilt - 2);
			add_index_range(rows, 0, nfilt - 3);
		}
		//cols = [ones(1,nfilt-1)         2:nfilt-1          2:nfilt-1];
		add_index_const(cols, 0, nfilt - 1);
		if (nfilt > 2)
		{
			add_index_range(cols, 1, nfilt - 2);
			add_index_range(cols, 1, nfilt - 2);
		}
		// data = [1+a(2)         a(3:nfilt)        ones(1,nfilt-2)    -ones(1,nfilt-2)];

		auto klen = rows.size();
		vectord data;
		data.resize(klen);
		data[0] = 1 + A[1];  int j = 1;
		if (nfilt > 2)
		{
			for (int i = 2; i < nfilt; i++)
				data[j++] = A[i];
			for (int i = 0; i < nfilt - 2; i++)
				data[j++] = 1.0;
			for (int i = 0; i < nfilt - 2; i++)
				data[j++] = -1.0;
		}

		vectord leftpad = subvector_reverse(X, nfact, 1);
		double _2x0 = 2 * X[0];
		std::transform(leftpad.begin(), leftpad.end(), leftpad.begin(), [_2x0](double val) {return _2x0 - val; });

		vectord rightpad = subvector_reverse(X, len - 2, len - nfact - 1);
		double _2xl = 2 * X[len - 1];
		std::transform(rightpad.begin(), rightpad.end(), rightpad.begin(), [_2xl](double val) {return _2xl - val; });

		double y0;
		vectord signal1, signal2, zi;

		signal1.reserve(leftpad.size() + X.size() + rightpad.size());
		append_vector(signal1, leftpad);
		append_vector(signal1, X);
		append_vector(signal1, rightpad);

		// Calculate initial conditions
		MatrixXd sp = MatrixXd::Zero(max_val(rows) + 1, max_val(cols) + 1);
		for (size_t k = 0; k < klen; ++k)
		{
			sp(rows[k], cols[k]) = data[k];
		}
		auto bb = VectorXd::Map(B.data(), B.size());
		auto aa = VectorXd::Map(A.data(), A.size());
		MatrixXd zzi = (sp.inverse() * (bb.segment(1, nfilt - 1) - (bb(0) * aa.segment(1, nfilt - 1))));
		zi.resize(zzi.size());

		// Do the forward and backward filtering
		y0 = signal1[0];
		std::transform(zzi.data(), zzi.data() + zzi.size(), zi.begin(), [y0](double val) { return val*y0; });
		filter(signal1, signal2);
		std::reverse(signal2.begin(), signal2.end());
		y0 = signal2[0];
		std::transform(zzi.data(), zzi.data() + zzi.size(), zi.begin(), [y0](double val) { return val*y0; });
		filter(signal2, signal1);
		Y = subvector_reverse(signal1, signal1.size() - nfact - 1, nfact);
	}

	void ZeroPhaseIIRFilter::add_index_range(vectori &indices, int beg, int end, int inc)
	{
		for (int i = beg; i <= end; i += inc)
		{
			indices.push_back(i);
		}
	}

	void ZeroPhaseIIRFilter::add_index_const(vectori &indices, int value, size_t numel)
	{
		while (numel--)
		{
			indices.push_back(value);
		}
	}

	void ZeroPhaseIIRFilter::append_vector(vectord &vec, const vectord &tail)
	{
		vec.insert(vec.end(), tail.begin(), tail.end());
	}

	ZeroPhaseIIRFilter::vectord ZeroPhaseIIRFilter::subvector_reverse(const vectord &vec, int idx_end, int idx_start)
	{
		vectord result(&vec[idx_start], &vec[idx_end + 1]);
		std::reverse(result.begin(), result.end());
		return result;
	}

	
	

	void ZeroPhaseIIRFilter::filter(const vectord &X, vectord &Y)
	{
		if (A.empty())
		{
			throw std::domain_error("The feedback filter coefficients are empty.");
		}
		if (std::all_of(A.begin(), A.end(), [](double coef) { return coef == 0; }))
		{
			throw std::domain_error("At least one of the feedback filter coefficients has to be non-zero.");
		}
		if (A[0] == 0)
		{
			throw std::domain_error("First feedback coefficient has to be non-zero.");
		}

		// Normalize feedback coefficients if a[0] != 1;
		auto a0 = A[0];
		if (a0 != 1.0)
		{
			std::transform(A.begin(), A.end(), A.begin(), [a0](double v) { return v / a0; });
			std::transform(B.begin(), B.end(), B.begin(), [a0](double v) { return v / a0; });
		}

		size_t input_size = X.size();
		size_t filter_order = std::max(A.size(), B.size());
		B.resize(filter_order, 0);
		A.resize(filter_order, 0);
		Y.resize(input_size);

		const double *x = &X[0];
		const double *b = &B[0];
		const double *a = &A[0];
		double *z = &Z[0];
		double *y = &Y[0];

		for (size_t i = 0; i < input_size; ++i)
		{
			size_t order = filter_order - 1;
			while (order)
			{
 				if (i >= order)
				{
					z[order - 1] = b[order] * x[i - order] - a[order] * y[i - order] + z[order];
				}
				--order;
			}
			y[i] = b[0] * x[i] + z[0];
		}
		Z.resize(filter_order - 1);
	}


}