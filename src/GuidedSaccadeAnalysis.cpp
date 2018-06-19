#define EXPORTING
#include "Analysis/GuidedSaccadeAnalysis.h"
#include <memory>
#include <algorithm>
#include <numeric>
#include "Analysis\Algorithms.h"
#include "ceres/ceres.h"

#ifndef M_PI
#define M_PI 3.14159265359
#endif
#define COS225 0.924

namespace Analysis {

	
		
		void GuidedSaccadeAnalysis::LoadTrials(const TaskDataType& data, const TaskEventsType& events, bool processNow) {
			this->data = data;
			this->events = events;
			SetTargetPositions();
			if (processNow) {
				ProcessAll();
			}

			
		}

		

		bool GuidedSaccadeAnalysis::ProcessAll() {

			double meanLx = nanmean(data.at("glx").data(), data.at("glx").size());
			double meanLy = nanmean(data.at("gly").data(), data.at("gly").size());
			double meanRx = nanmean(data.at("grx").data(), data.at("grx").size());
			double meanRy = nanmean(data.at("gry").data(), data.at("gry").size());
			double t0 = data.at("t")[0];

			xL.reserve(data.at("glx").size());
			xR.reserve(data.at("glx").size());
			yL.reserve(data.at("glx").size());
			yR.reserve(data.at("glx").size());
			t.reserve(data.at("t").size());
			std::for_each(data.at("glx").begin(), data.at("glx").end(), [this, meanLx](double x) {this->xL.insert(this->xL.end(), x - meanLx); });
			std::for_each(data.at("gly").begin(), data.at("gly").end(), [this, meanLy](double x) {this->yL.insert(this->yL.end(), x - meanLy); });
			std::for_each(data.at("grx").begin(), data.at("grx").end(), [this, meanRx](double x) {this->xR.insert(this->xR.end(), x - meanRx); });
			std::for_each(data.at("gry").begin(), data.at("gry").end(), [this, meanRy](double x) {this->yR.insert(this->yR.end(), x - meanRy); });
			std::for_each(data.at("t").begin(), data.at("t").end(), [this, t0](double x) {this->t.insert(this->t.end(), x - t0); });

			std::unique_ptr<EngbertKlieglThreshold> thres = std::make_unique<EngbertKlieglThreshold>();
			thres->SetLambda(4.0);
			VelocityBasedSaccadeDetectionAlgorithm det;
			det.SetThreshold(&*thres);
			det.SetMinimumSaccadeDuration(20);
			det.SetMaximumSaccadeSize(30);
			det.SetTaskType(VelocityBasedSaccadeDetectionAlgorithm::FIXATION);
			std::vector<std::size_t> dataIndices(xL.size());
			std::iota(std::begin(dataIndices), std::end(dataIndices), 0); //0 is the starting number

			vecLeft = det.GetSaccadesFromTrace(const_cast<double*>(t.data()), const_cast<double*>(xL.data()), const_cast<double*>(yL.data()), xL.size(), dataIndices);
			vecRight = det.GetSaccadesFromTrace(const_cast<double*>(t.data()), const_cast<double*>(xR.data()), const_cast<double*>(yR.data()), xR.size(), dataIndices);

			det.RemoveMonoculars(vecLeft, vecRight);
			//det.RemoveBinocularOvershoots(vecLeft, vecRight);
			

	
			for (std::size_t k = 1; k < targetT.size(); k++) {
				SaccadeVector vL;
				SaccadeVector vR;
				double tMin = 100;
				double tMax = 800;
				hitOrMiss[k - 1] = false;
				for (std::size_t j = 0; j < this->vecLeft.size(); j++) {
					bool leftIn = (vecLeft[j].startTimestamp > targetT[k] + tMin) && (vecLeft[j].startTimestamp < targetT[k] + tMax);
					bool rightIn = (vecRight[j].startTimestamp > targetT[k] + tMin) && (vecRight[j].startTimestamp < targetT[k] + tMax);

					bool leftUncued, rightUncued;
					if (k < targetT.size() - 1) {
						leftUncued = (vecLeft[j].startTimestamp > targetT[k] + tMax) && (vecLeft[j].startTimestamp < targetT[k + 1] +tMin);
						rightUncued = (vecRight[j].startTimestamp > targetT[k] + tMax) && (vecRight[j].startTimestamp < targetT[k + 1] + tMin);
					}
					else {
						leftUncued = (vecLeft[j].startTimestamp > targetT[k] + tMax);
						rightUncued = (vecRight[j].startTimestamp > targetT[k] + tMax);
					}
					if (leftIn || rightIn) {
						vL.push_back(vecLeft[j]);
						vR.push_back(vecRight[j]);
					}
					else if (leftUncued || rightUncued) {
						if ((vecLeft[j].magnitude + vecRight[j].magnitude) / 2.0 > 3.0) {
							numberUncued[k-1]++;
						}
					}
				}



				tMin = 0;
				tMax = waveformDuration;
				auto initIt = std::find_if(t.begin(), t.end(), [this, tMin, k](const double& x) {return x >= this->targetT[k] + tMin; });
				auto endIt = std::find_if(t.begin(), t.end(), [this, tMax, k](const double& x) {return x > this->targetT[k] + tMax; });
				auto initD = std::distance(t.begin(), initIt);
				auto endD = std::distance(t.begin(), endIt);
				//fix this

				leftWaveforms[2 * (k - 1)].resize(150);
				leftWaveforms[2 * (k - 1) + 1].resize(150);
				rightWaveforms[2 * (k - 1)].resize(150);
				rightWaveforms[2 * (k - 1) + 1].resize(150);

				std::copy(xL.begin() + initD, xL.begin() + initD + 150 - 1, leftWaveforms[2 * (k - 1)].begin());
				std::copy(yL.begin() + initD, yL.begin() + initD + 150 - 1, leftWaveforms[2 * (k - 1) + 1].begin());
				std::copy(xR.begin() + initD, xR.begin() + initD + 150 - 1, rightWaveforms[2 * (k - 1)].begin());
				std::copy(yR.begin() + initD, yR.begin() + initD + 150 - 1, rightWaveforms[2 * (k - 1) + 1].begin());

			
				double xL0 = leftWaveforms[2 * (k - 1)][0];
				double yL0 = leftWaveforms[2 * (k - 1) +1][0];
				double xR0 = rightWaveforms[2 * (k - 1)][0];
				double yR0 = rightWaveforms[2 * (k - 1) +1][0];
				std::for_each(leftWaveforms[2 * (k - 1)].begin(), leftWaveforms[2 * (k - 1)].end(), [xL0](double& x) {x = x - xL0; });
				std::for_each(leftWaveforms[2 * (k - 1) +1].begin(), leftWaveforms[2 * (k - 1) + 1].end(), [yL0](double& x) {x = x - yL0; });
				std::for_each(rightWaveforms[2 * (k - 1)].begin(), rightWaveforms[2 * (k - 1)].end(), [xR0](double& x) {x = x - xR0; });
				std::for_each(rightWaveforms[2 * (k - 1) +1].begin(), rightWaveforms[2 * (k - 1) + 1].end(), [yR0](double& x) {x = x - yR0; });

	
		
				double totalMagnitude = 0.0;
				double multistepMagnitude = 0.0;
				double pvl = -1.0;
				double pvr = -1.0;
				double ml = 0.0;
				double mr = 0.0;
				int lastSaccade = -1;
				for (std::size_t j = 0; j < vL.size(); j++) {
					double dotL = (vL[j].endX - vL[j].startX)*(targetX[k] - targetX[k - 1]) + (vL[j].endY - vL[j].startY)*(targetY[k] - targetY[k - 1]);
					double dotR = (vR[j].endX - vR[j].startX)*(targetX[k] - targetX[k - 1]) + (vR[j].endY - vR[j].startY)*(targetY[k] - targetY[k - 1]);
					double cosL = dotL / vL[0].magnitude / cuedMagnitudes[k - 1];
					double cosR = dotR / vR[0].magnitude / cuedMagnitudes[k - 1];
					if (std::max(cosL, cosR) > COS225) {
						
						if (hitOrMiss[k - 1] == true && lastSaccade >= 0) { //not the first saccade detected
							if (vL[j].startTimestamp - vL[j - lastSaccade].endTimestamp < 300) {
								multistep[k - 1] = true;
								multistepMagnitude += (vL[j].magnitude + vR[j].magnitude) / 2.0;
							}
							else {
								break;
							}
						}
						if (hitOrMiss[k-1] == false) { // first saccade detected in the right direction
							magnitudeGainLeft[k - 1] = vL[j].magnitude / cuedMagnitudes[k - 1];
							magnitudeGainRight[k - 1] = vR[j].magnitude / cuedMagnitudes[k - 1];
							latenciesLeft[k - 1] = vL[j].startTimestamp - targetT[k];
							latenciesRight[k - 1] = vR[j].startTimestamp - targetT[k];
							hitOrMiss[k - 1] = true;
							slopes[k - 1] = ((vL[j].peakVelocity / vL[j].magnitude) + (vR[j].peakVelocity / vR[j].magnitude)) / 2.0;
						}
						lastSaccade = j;
						totalMagnitude += (vL[j].magnitude + vR[j].magnitude) / 2.0;
						
						if (vL[j].peakVelocity > pvl) {
							pvl = vL[j].peakVelocity;
							ml = vL[j].magnitude;
						}
						if (vR[j].peakVelocity > pvr) {
							pvr = vR[j].peakVelocity;
							mr = vR[j].magnitude;
						}
					}
				}
				peakVelocitiesLeft[k - 1] = pvl;
				peakVelocitiesRight[k - 1] = pvr;
				magnitudesLeft[k - 1] = ml;
				magnitudesRight[k - 1] = mr;
				if (multistep[k - 1]) {
					multistepFractionMagnitude[k - 1] = multistepMagnitude / totalMagnitude;
				}
				
			}
			

			doneProcessing = true;
			return true;

		}

		bool GuidedSaccadeAnalysis::Clear() {
			return true;
		}

		GuidedSaccadeAnalysis::GuidedSaccadeAnalysis() {
			waveformDuration = 600;
			doneProcessing = false;

		}
		void GuidedSaccadeAnalysis::SetTargetPositions() {

			targetT.clear();
			targetX.clear();
			targetY.clear();

			for (std::size_t k = 0; k < events.size(); k++) {
				targetT.push_back(events[k].inTrialTS*1000);
				double xt, yt, zt;
				if (events[k].message.find("CENTER") != std::string::npos) {
					xt = 0.0;
					yt = 0.0;
					zt = 0.0;
				}
				else {
					std::size_t pos = events[k].message.find(':');
					std::size_t closePos = events[k].message.find('|');
					std::string coord = events[k].message.substr(pos + 1, closePos - pos - 2);
					std::istringstream ss(coord);
					ss >> xt;
					ss >> yt;
					ss >> zt;
				}
				targetX.push_back(atan2(xt, 0.57)*180.0 / M_PI);
				targetY.push_back(atan2(yt, 0.57)*180.0 / M_PI);

			}
			
			cuedMagnitudes.resize(targetT.size() - 1);
			saccadeCues.resize(targetT.size() - 1, GuidedSaccadeAnalysis::NONE);
			hitOrMiss.resize(targetT.size() - 1, false);
			magnitudeGainLeft.resize(targetT.size() - 1, std::nan(""));
			magnitudeGainRight.resize(targetT.size() - 1, std::nan(""));
			latenciesLeft.resize(targetT.size() - 1, std::nan(""));
			latenciesRight.resize(targetT.size() - 1, std::nan(""));
			multistep.resize(targetT.size() - 1, false);
			multistepFractionMagnitude.resize(targetT.size() - 1, std::nan(""));
			leftWaveforms.resize((targetT.size() - 1) * 2);
			rightWaveforms.resize((targetT.size() - 1) * 2);
			peakVelocitiesLeft.resize(targetT.size() - 1, std::nan(""));
			peakVelocitiesRight.resize(targetT.size() - 1, std::nan(""));
			magnitudesLeft.resize(targetT.size() - 1, std::nan(""));
			magnitudesRight.resize(targetT.size() - 1, std::nan(""));
			numberUncued.resize(targetT.size() - 1, 0);
			slopes.resize(targetT.size() - 1, std::nan(""));

			for (std::size_t k = 1; k < targetT.size(); k++) {

				cuedMagnitudes[k - 1] = std::sqrt((targetX[k] - targetX[k - 1])*(targetX[k] - targetX[k - 1]) + (targetY[k] - targetY[k - 1])*(targetY[k] - targetY[k - 1]));

				if (targetX[k] > targetX[k - 1]) {
					saccadeCues[k - 1] |= GuidedSaccadeAnalysis::RIGHT;
				}
				else if (targetX[k] < targetX[k - 1]) {
					saccadeCues[k - 1] |= GuidedSaccadeAnalysis::LEFT;
				}

				if (targetY[k] > targetY[k - 1]) {
					saccadeCues[k - 1] |= GuidedSaccadeAnalysis::UP;
				}
				else if (targetY[k] < targetY[k - 1]) {
					saccadeCues[k - 1] |= GuidedSaccadeAnalysis::DOWN;
				}

				if (targetX[k] == 0.0 && targetY[k] == 0.0) {
					saccadeCues[k - 1] |= GuidedSaccadeAnalysis::TO_CENTER;
				}
				else if ((targetX[k] != 0.0 || targetY[k] != 0.0)) {
					saccadeCues[k - 1] |= GuidedSaccadeAnalysis::TO_PERIPHERY;
				}

			}

		}
		void GuidedSaccadeAnalysis::GetMagnitudeGain(double& left, double& right, unsigned int filter) const {
			std::vector<std::size_t> indices = filteredIndices(filter);
			std::vector<double> leftFilteredMagnitudes, rightFilteredMagnitudes;
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				leftFilteredMagnitudes.push_back(magnitudeGainLeft[*it]);
				rightFilteredMagnitudes.push_back(magnitudeGainRight[*it]);
			}

			left = nanmean(leftFilteredMagnitudes.data(), leftFilteredMagnitudes.size());
			right = nanmean(rightFilteredMagnitudes.data(), rightFilteredMagnitudes.size());

		}
		void GuidedSaccadeAnalysis::GetAverageLatency(double& left, double& right, unsigned int filter) const {
			std::vector<std::size_t> indices = filteredIndices(filter);
			std::vector<double> leftFilteredLatencies, rightFilteredLatencies;
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				leftFilteredLatencies.push_back(latenciesLeft[*it]);
				rightFilteredLatencies.push_back(latenciesRight[*it]);
			}

			left = nanmean(leftFilteredLatencies.data(), leftFilteredLatencies.size());
			right = nanmean(rightFilteredLatencies.data(), rightFilteredLatencies.size());

		}
		double GuidedSaccadeAnalysis::GetFractionMultistep(unsigned int filter) const {
			std::vector<std::size_t> indices = filteredIndices(filter);
			unsigned int count = 0;
			unsigned int totalSaccades = 0;
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				if (hitOrMiss[*it]) {
					totalSaccades++;
					if(multistep[*it]) count++;
				}
			}

			return double(count)/double(totalSaccades);

		}

		double GuidedSaccadeAnalysis::GetMainSequenceSlope(unsigned int filter) const {
			std::vector<std::size_t> indices = filteredIndices(filter);
			unsigned int count = 0;
			double acc = 0;
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				if (hitOrMiss[*it]) {
					acc += slopes[*it];
					count++;
				}
			}

			return acc / double(count);

		}
		double GuidedSaccadeAnalysis::GetMultistepMagnitudeFrac(unsigned int filter) const {
			std::vector<std::size_t> indices = filteredIndices(filter);
			unsigned int count = 0;
			double acc = 0;
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				if (hitOrMiss[*it] && multistep[*it]) {
					acc += multistepFractionMagnitude[*it];
					count++;
				}
			}

			return acc / double(count);

		}

		SaccadeVector GuidedSaccadeAnalysis::GetLeftSaccades(unsigned int filter) const {
			return vecLeft;


		}
		SaccadeVector GuidedSaccadeAnalysis::GetRightSaccades(unsigned int filter) const {
			return vecRight;
		}
		std::vector<double> GuidedSaccadeAnalysis::GetSaccadicWaveforms(std::vector<std::vector<double> >& lw, std::vector<std::vector<double> > & rw, unsigned int filter, bool amplitudeNormalized) {
			std::vector<std::size_t> indices = filteredIndices(filter);
			lw.resize(0);
			rw.resize(0);
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				
				lw.push_back(leftWaveforms[*it * 2]);
				lw.push_back(leftWaveforms[*it * 2 + 1]);
				rw.push_back(rightWaveforms[*it * 2]);
				rw.push_back(rightWaveforms[*it * 2 + 1]);
				if (amplitudeNormalized) {
					double amplitudeX = targetX[*it+1] - targetX[*it];
					double amplitudeY = targetY[*it + 1] - targetY[*it];
					std::for_each((*(lw.end() - 2)).begin(), (*(lw.end() - 2)).end(), [amplitudeX](double& x) {x = x / amplitudeX; });
					std::for_each((*(lw.end() - 1)).begin(), (*(lw.end() - 1)).end(), [amplitudeY](double& x) {x = x / amplitudeY; });
					std::for_each((*(rw.end() - 2)).begin(), (*(rw.end() - 2)).end(), [amplitudeX](double& x) {x = x / amplitudeX; });
					std::for_each((*(rw.end() - 1)).begin(), (*(rw.end() - 1)).end(), [amplitudeY](double& x) {x = x / amplitudeY; });


				}
			}

			for (std::size_t k = 1; k < lw.size(); k++) {
				if (lw[k].size() > lw[0].size()) lw[k].pop_back();
				else if (lw[k].size() < lw.size()) lw[k].push_back(std::nan(""));
			}
			for (std::size_t k = 1; k < rw.size(); k++) {
				if (rw[k].size() > rw[0].size()) rw[k].pop_back();
				else if (rw[k].size() < rw.size()) rw[k].push_back(std::nan(""));
			}

			std::vector<double> timestamps(lw[0].size());
			std::generate(timestamps.begin(), timestamps.end(), [n = -4]() mutable {n += 4; return n; });

			return timestamps;

		}
		
		unsigned int GuidedSaccadeAnalysis::GetNumberOfMisses(unsigned int filter) const {
			std::vector<std::size_t> indices = filteredIndices(filter);
			unsigned int count = 0;
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				if (!hitOrMiss[*it]) count++;
			}

			return count;
		}

		std::vector<std::size_t>  GuidedSaccadeAnalysis::filteredIndices(unsigned int filter) const {
			std::vector<std::size_t> out;
			for (std::size_t k = 0; k < this->saccadeCues.size(); k++) {
				if (filter & saccadeCues[k]) {
					out.push_back(k);
				}
			}

			return out;

		}

		unsigned int GuidedSaccadeAnalysis::GetNumberUncued(unsigned int filter) const {
			std::vector<std::size_t> indices = filteredIndices(filter);
			unsigned int count = 0;
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				count += numberUncued[*it];
			}

			return count;

		}

		Json::Value GuidedSaccadeAnalysis::ToJSONObject() {
			Json::Value root;
			double leftLatency, rightLatency;
			GetAverageLatency(leftLatency, rightLatency, ALL);
			root["AverageLatency"] = (leftLatency + rightLatency) / 2.0;
			double leftGain, rightGain;
			GetMagnitudeGain(leftGain, rightGain, ALL);
			root["AverageGain"] = (leftGain + rightGain) / 2.0;
			root["NumberOfMisses"] = GetNumberOfMisses(ALL);
			root["NumberOfMultistep"] = GetFractionMultistep(ALL);
			root["MultistepFirstGain"] = GetMultistepMagnitudeFrac(ALL);
			root["MainSequence"] = GetMainSequenceSlope(ALL);
			root["NumberUncuedSaccades"] = GetNumberUncued(ALL);

			GetAverageLatency(leftLatency, rightLatency, TO_PERIPHERY);
			root["Unpredictable"]["AverageLatency"] = (leftLatency + rightLatency) / 2.0;
			GetMagnitudeGain(leftGain, rightGain, TO_PERIPHERY);
			root["Unpredictable"]["AverageGain"] = (leftGain + rightGain) / 2.0;
			root["Unpredictable"]["NumberOfMisses"] = GetNumberOfMisses(TO_PERIPHERY);
			root["Unpredictable"]["NumberOfMultistep"] = GetFractionMultistep(TO_PERIPHERY);
			root["Unpredictable"]["MultistepFirstGain"] = GetMultistepMagnitudeFrac(TO_PERIPHERY);
			root["Unpredictable"]["MainSequence"] = GetMainSequenceSlope(TO_PERIPHERY);
			
			GetAverageLatency(leftLatency, rightLatency, TO_CENTER);
			root["Predictable"]["AverageLatency"] = (leftLatency + rightLatency) / 2.0;
			GetMagnitudeGain(leftGain, rightGain, TO_CENTER);
			root["Predictable"]["AverageGain"] = (leftGain + rightGain) / 2.0;
			root["Predictable"]["NumberOfMisses"] = GetNumberOfMisses(TO_CENTER);
			root["Predictable"]["NumberOfMultistep"] = GetFractionMultistep(TO_CENTER);
			root["Predictable"]["MultistepFirstGain"] = GetMultistepMagnitudeFrac(TO_CENTER);
			root["Predictable"]["MainSequence"] = GetMainSequenceSlope(TO_CENTER);

			return root;
		}

	
}