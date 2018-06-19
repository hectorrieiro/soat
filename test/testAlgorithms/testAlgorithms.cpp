#include "testAlgorithms.h"
#include "Analysis\Algorithms.h"
#include <iostream>
#include <chrono>
#include <random>

testAlgorithms::testAlgorithms() {

}


testAlgorithms::~testAlgorithms() {

}

void testAlgorithms::SetUp() {

}

void testAlgorithms::TearDown() {

}

TEST_F(testAlgorithms, SinglePredictorRegressionNoError) {
	Analysis::NeweyWestRegression model;
	
	std::vector<double> y = { 3, 5, 7, 9, 11 };
	std::vector<double> x0 = { 1, 1, 1, 1, 1 };
	std::vector<double> x1 = { 0, 1, 2, 3, 4 };
	std::vector<std::vector<double>> X;
	X.push_back(x0);
	X.push_back(x1);
	model.SetX(X);
	model.SetY(y);
	model.SetMaxLag(0);
	model.Solve();
	std::vector<double> r = model.GetPredictors();
	ASSERT_EQ(r.size(), 2);
	EXPECT_NEAR(r[0], 3.0, 0.001);
	EXPECT_NEAR(r[1], 2.0, 0.001);
	EXPECT_NEAR(model.GetRsquared(), 1.0, 0.001);


}

TEST_F(testAlgorithms, SinglePredictorRegressionError) {
	Analysis::NeweyWestRegression model;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::normal_distribution<double> distribution(0.0, 2.0);
	unsigned int numSamples = 25000;
	std::vector<double> y(numSamples);
	std::vector<double> x0;
	x0.assign(numSamples, 1.0);
	std::vector<double> x1(numSamples);
	std::iota(x1.begin(), x1.end(), 0.0);
	for (std::size_t i = 0; i < numSamples; i++) {
		y[i] += 3.0*x0[i] + 2*x1[i]+distribution(generator);
	}
	std::vector<std::vector<double>> X;
	X.push_back(x0);
	X.push_back(x1);
	model.SetX(X);
	model.SetY(y);
	model.SetMaxLag(0);
	model.Solve();
	std::vector<double> r = model.GetPredictors();
	std::vector<double> e = model.GetPredictorsSem();
	ASSERT_EQ(r.size(), 2);
	ASSERT_EQ(e.size(), 2);
	std::cerr << "Estimated parameters: ["<< r[0]<<", "<<r[1]<<"]" << std::endl;
	std::cerr << "Estimation SEM: ["<<e[0] << ", " << e[1] << "]" << std::endl;
	EXPECT_NEAR(r[0], 3.0, 0.1);
	EXPECT_NEAR(r[1], 2.0, 0.1);
	EXPECT_NEAR(model.GetRsquared(), 1.0, 0.1);


}

TEST_F(testAlgorithms, SinglePredictorRegressionErrorCorrelated) {
	Analysis::NeweyWestRegression model;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::normal_distribution<double> distribution(0.0, 2.0);
	unsigned int numSamples = 25000;
	std::vector<double> y(numSamples);
	std::vector<double> x0;
	x0.assign(numSamples, 1.0);
	std::vector<double> x1(numSamples);
	std::iota(x1.begin(), x1.end(), 0.0);
	std::vector<double> errors(numSamples);
	
	for (std::size_t i = 0; i < numSamples; i++) {

		errors[i] = distribution(generator);
		if (i > 0) {
			errors[i] += 0.3*errors[i - 1];
		}
		y[i] += 3.0*x0[i] + 2 * x1[i] + errors[i];
	}
	std::vector<std::vector<double>> X;
	X.push_back(x0);
	X.push_back(x1);
	model.SetX(X);
	model.SetY(y);
	model.SetMaxLag(10);
	model.Solve();
	std::vector<double> r = model.GetPredictors();
	std::vector<double> e = model.GetPredictorsSem();
	ASSERT_EQ(r.size(), 2);
	ASSERT_EQ(e.size(), 2);
	std::cerr << "Estimated parameters: [" << r[0] << ", " << r[1] << "]" << std::endl;
	std::cerr << "Estimation SEM: [" << e[0] << ", " << e[1] << "]" << std::endl;
	EXPECT_NEAR(r[0], 3.0, 0.1);
	EXPECT_NEAR(r[1], 2.0, 0.1);
	EXPECT_NEAR(model.GetRsquared(), 1.0, 0.1);


}


	