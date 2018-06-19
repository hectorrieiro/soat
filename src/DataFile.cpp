#define EXPORTING
#include "json/value.h"
#include "json/reader.h"
#include "Analysis/DataFile.h"
#include "Analysis/DataTypes.h"
#include "Analysis/HDF5.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>


#define M_PI 3.1415926

namespace Analysis {

	class SampleToDva {
	public:
		virtual void operator()(double& x, double& y) {};
		virtual ~SampleToDva() = default;
	};

	class EyeTechDSSampleToDVA : public SampleToDva {
	public:
	EyeTechDSSampleToDVA (int _pxHoriz, int _pxVert, double _distance, double _dimHoriz, double _dimVert) {
		pxHoriz = _pxHoriz;
		pxVert = _pxVert;
		distance = _distance;
		dimHoriz = _dimHoriz;
		dimVert = _dimVert;
	}
	void operator()(double& x, double& y) override {
		x = atan2((x - 50.0) / 100.0*dimHoriz, distance) / M_PI*180.0;
		y = atan2((y - 50.0) / 100.0*dimVert, distance) / M_PI*180.0;
	}
	virtual ~EyeTechDSSampleToDVA() = default;
	private:
		int pxHoriz, pxVert;
		double distance, dimHoriz, dimVert;
	};

	class SMIETHMDSampleToDVA : public SampleToDva {
	public:
		SMIETHMDSampleToDVA() = default;

		void SMIETHMDSampleToDVA::operator()(double& x, double& y) override {
			x = atan2(x - 1080.0, virtualDistanceX) / M_PI*180.0;
			y = -atan2(y - 600.0, virtualDistanceY) / M_PI*180.0;
		}
		virtual ~SMIETHMDSampleToDVA() = default;
	private:
		const double SMIETHMDSampleToDVA::virtualDistanceX = 1080.0 / tan(111.9 / 2.0 / 180 * M_PI);
		const double SMIETHMDSampleToDVA::virtualDistanceY = 600.0 / tan(105.6 / 2.0 / 180 * M_PI);
	};
	
	DataFile::DataFile(const std::string& filename) {
		this->filename = filename;
	}

	bool DataFile::Process() {
		bool success = false;
		taskNames.clear();
		results.clear();
		events.clear();
		if (isHDFFile()) {
			success = ProcessHDFFile();
		}
		else {
			success = ProcessTextFile();
		}
		
		return success;
	}

	bool DataFile::isHDFFile() {

		return Util::HDF5DataFile::IsValidHDF5(filename);
	}

	bool DataFile::ProcessHDFFile() {
		
		//read necessary information from headers
		Util::HDF5DataFile file(filename);
		std::string configurationString = file.GetConfiguration();
		
		
		ProcessExperimentalSetup(configurationString);

		std::string recordingInfo = file.GetRecordingInfo();
		ProcessRecordingInfo(recordingInfo);
		
		//now segment data and events
		results.clear();
		events.clear();
		std::vector<std::string> eventDescriptions;
		std::vector<double> eventTimestamps;
		file.GetEvents(eventTimestamps, eventDescriptions);
		std::size_t currentEventProcessing = 0;
		std::vector<double> recordingTimestamps;
		file.GetTimestamps(recordingTimestamps);
		std::vector<double> recordingLX, recordingLY, recordingRX, recordingRY, recordingLP, recordingRP;
		file.GetLeftGaze(recordingLX, recordingLY);
		file.GetRightGaze(recordingRX, recordingRY);
		file.GetLeftPupilSize(recordingLP);
		file.GetRightPupilSize(recordingRP);

		while (currentEventProcessing < eventDescriptions.size()) {

			boost::algorithm::trim(eventDescriptions[currentEventProcessing]);
			std::string line = eventDescriptions[currentEventProcessing];
			if (line.find("REPEATING LAST TEST") != line.npos) {
				if (!results.empty()) {
					results.pop_back();
					events.pop_back();
					taskNames.pop_back();
				}
				continue;

			}
			std::size_t pos = line.find("START");

			if (pos != line.npos) { //start of a trial
				std::vector<double> timestamps;
				double t0 = eventTimestamps[currentEventProcessing];
				double softT0;
				std::string s = eventDescriptions[currentEventProcessing];
				softT0 = std::stod(s.substr(s.find_last_of(':') + 1)) / 1000.0;
				Analysis::TaskDataType data;
				Analysis::TaskEventsType taskEvents;
				std::string taskType;
				getTaskTypeAndParameters(line, taskType, data);
				taskNames.push_back(taskType);
				double tEnd;
				//look for task end event and collect events along the way
				bool reachedEnd = false;
				while (!reachedEnd) {
					currentEventProcessing++;
					line = eventDescriptions[currentEventProcessing];
					boost::algorithm::trim(line);
					std::size_t pos = line.find("END");
					if (pos != line.npos) { //end of a trial
						reachedEnd = true;
						tEnd = eventTimestamps[currentEventProcessing];
					}
					else { //process events inside the trial
						Analysis::EventType ev;
						std::string st = eventDescriptions[currentEventProcessing];
						double softT = std::stod(st.substr(st.find_last_of(':') + 1)) / 1000.0;

						ev.timestamp = (eventTimestamps[currentEventProcessing] - t0) / 1000.0;
						ev.inTrialTS = (softT - softT0) / 1000.0;
						ev.message = eventDescriptions[currentEventProcessing];
						taskEvents.push_back(ev);
					}
				}
				std::vector<double> t, glx, gly, grx, gry, lp, rp;
				//find indices for data in trial
				auto itBegin = std::find_if(recordingTimestamps.begin(), recordingTimestamps.end(), [t0](const double& x) -> bool {return x >= t0; });
				auto itEnd = std::find_if(recordingTimestamps.begin(), recordingTimestamps.end(), [tEnd](const double& x)->bool {return x >= tEnd; });

				std::size_t numSamples = itEnd - itBegin;
				std::size_t indexBegin = itBegin - recordingTimestamps.begin();
				std::size_t indexEnd = itEnd - recordingTimestamps.begin();
				t.resize(numSamples);

				glx.resize(numSamples);
				grx.resize(numSamples);
				gly.resize(numSamples);
				gry.resize(numSamples);
				lp.resize(numSamples);
				rp.resize(numSamples);
				std::size_t n = 0;
				for (std::size_t k = indexBegin; k < indexEnd; k++) {
					t[n] = (recordingTimestamps[k] - t0) / 1000.0;
					double x = recordingLX[k];
					double y = recordingLY[k];
					(*dvaConversion)(x, y);
					glx[n] = x;
					gly[n] = y;
					x = recordingRX[k];
					y = recordingRY[k];
					(*dvaConversion)(x, y);
					grx[n] = x;
					gry[n] = y;
					lp[n] = recordingLP[k];
					rp[n] = recordingRP[k];
					n++;
				}

				data["t"] = t;
				data["glx"] = glx;
				data["gly"] = gly;
				data["grx"] = grx;
				data["gry"] = gry;
				data["lp"] = lp;
				data["rp"] = rp;
				results.push_back(data);
				events.push_back(taskEvents);

			}

			currentEventProcessing++;
		}

		return true;
	}

	bool DataFile::ProcessTextFile() {
		std::ifstream f;
		f.open(filename);

		std::string line;
		std::getline(f, line);
		boost::trim(line);
		Json::Value config;


		while (!f.eof()) {


			if (line.compare("[CONFIGURATION]") == 0) {
				std::string configurationString("");
				std::getline(f, line);
				boost::trim(line);
				while (line[0] != '[') {
					configurationString.append(line);
					std::getline(f, line);
					boost::trim(line);
				}

				ProcessExperimentalSetup(configurationString);

				continue;


			}
			else if (line.compare("[TASKS]") == 0) {

			}
			else if (line.compare("[DATA HEADER]") == 0) {
				std::string stringval("");
				std::getline(f, line);
				boost::trim(line);
				while (line[0] != '[') {
					stringval.append(line);
					std::getline(f, line);
					boost::trim(line);
				}
				ProcessRecordingInfo(stringval);
				continue;
			}
			else if (line.compare("[START DATA]") == 0) {

				readTextData(f, taskNames, results, events);

			}
			std::getline(f, line);
			boost::trim(line);
		}

		f.close();


		return true;
	}

	void DataFile::getTaskTypeAndParameters(const std::string& line, std::string& taskType, Analysis::TaskDataType& parameters) {

		using boost::algorithm::trim;
		using namespace boost;
		std::size_t pos = line.find("START");
		std::size_t dotPos = line.find(".", pos);
		taskType = line.substr(pos + 5, dotPos - (pos + 5));
		trim(taskType);
		if (std::strcmp(taskType.c_str(), "OKN") == 0) {
			std::size_t barPos = line.find("|");
			std::string speedStr = line.substr(dotPos + 1, barPos - dotPos - 3);
			trim(speedStr);
			parameters["speed"] = std::vector<double>(1, std::stod(speedStr));
		}
		else if (std::strcmp(taskType.c_str(), "FIXATION") == 0) {

		}
		else if (std::strcmp(taskType.c_str(), "GUIDED SACCADE") == 0) {

		}
		else if (std::strcmp(taskType.c_str(), "PURSUIT") == 0) {
			std::size_t movementDotPos = line.find('.', dotPos + 1);
			std::string movementString = line.substr(dotPos + 1, movementDotPos - 1 - dotPos);
			trim(movementString);
			std::vector<double> paramVector(6, 0.0);
			paramVector[0] = std::strcmp(movementString.c_str(), "LINEAR") == 0 ? 1.0 : 0.0;
			std::size_t aux = line.find('.', movementDotPos + 1);
			std::size_t speedDotPos = line.find('.', aux + 1);
			std::string speedStr = line.substr(movementDotPos + 1, speedDotPos - 1 - movementDotPos);
			trim(speedStr);
			paramVector[1] = std::stod(speedStr);
			std::size_t directionDot = line.find('.', speedDotPos + 1);
			std::string directionStr = line.substr(speedDotPos + 1, directionDot - speedDotPos - 1);
			trim(directionStr);
			if (std::strcmp(directionStr.c_str(), "HORIZONTAL") == 0) {
				paramVector[2] = 0.0;
			}
			else if (std::strcmp(directionStr.c_str(), "VERTICAL") == 0) {
				paramVector[2] = 1.0;
			}
			else if (std::strcmp(directionStr.c_str(), "CONVERGENCE") == 0) {
				paramVector[2] = 2.0;
			}

			std::size_t colonPos = line.find(':', directionDot + 1);
			std::size_t barPos = line.find('|', colonPos + 1);
			std::string posString = line.substr(colonPos + 1, barPos - colonPos - 1);
			char_separator<char> sep(", ");
			tokenizer<char_separator<char>> tokens(posString, sep);
			std::size_t idx = 3;
			for (const auto& t : tokens) {
				double val;
				trim(const_cast<std::string&>(t));
				val = stod(t);
				paramVector[idx++] = val;
			}


			parameters["motionParameters"] = paramVector;

		}
		else if (std::strcmp(taskType.c_str(), "PUPIL") == 0) {

		}
	}

	bool DataFile::ProcessExperimentalSetup(const std::string& stringval) {

		Json::Value config;

		std::stringstream ss(stringval);
		ss >> config;
		Json::Value setup = config["ExperimentalSetup"];

		if (strcmp(setup["EyeTracker"].asString().c_str(), "SMIETHMDTracker") == 0) {
			dvaConversion = std::make_unique<SMIETHMDSampleToDVA>();
		}
		else if (strcmp(setup["EyeTracker"].asString().c_str(), "EyeTechDSTracker") == 0) {
			dvaConversion = std::make_unique<EyeTechDSSampleToDVA>(setup["HorizontalScreenSizePX"].asInt(),
				setup["VerticalScreenSizePX"].asInt(),
				setup["SubjectDistanceToScreenCM"].asDouble(),
				setup["HorizontalScreenSizeCM"].asDouble(),
				setup["VerticalScreenSizeCM"].asDouble());
		}
		else {
			dvaConversion = std::make_unique<SampleToDva>();
		}
		return true;
	}

	bool DataFile::ProcessRecordingInfo(const std::string& info) {
		std::stringstream infoStream(info);
		std::string line;
		std::getline(infoStream, line);
		boost::algorithm::trim(line);

		while (!infoStream.eof()) {
			using namespace boost;
			if (line.size() == 0) continue;
			char_separator<char> sep("=");
			tokenizer<char_separator<char>> tokens(line, sep);
			std::vector<std::string> components;
			for (const auto& t : tokens) {
				std::string str = t.c_str();
				trim(str);
				components.push_back(str);
			}
			if (components.size() != 2) continue;
			if (components[0].compare("Recording started") == 0) {
				this->testDate = components[1];
			}
			else if (components[0].compare("Annotation") == 0) {
				this->testNotes = components[1];
			}
			else if (components[0].compare("Tracker") == 0) {
				if (components[1].compare("SMIETHMD") == 0) {
					dvaConversion = std::make_unique<SMIETHMDSampleToDVA>();
				}
			}
			else if (components[0].compare("Subject ID") == 0) {
				this->subjectID = components[1];
			}
			else if (components[0].compare("Subject DOB") == 0) {
				this->subjectDOB = components[1];
			}
			std::string line;
			std::getline(infoStream, line);
			boost::algorithm::trim(line);

		}

		return true;

	}

	void DataFile::readTextData(std::istream& input, std::vector<std::string>& taskTypes, Analysis::ExperimentResultsType& results, Analysis::ExperimentEventsType& events) {
		using namespace std;
		using namespace boost;

		//first line should contain the column headers, use this to pick the necessary columns
		string line;
		getline(input, line);
		trim(line);

		char_separator<char> sep(", ");
		tokenizer<char_separator<char>> tokens(line, sep);
		int timestampColumn, gazeRXColumn, gazeRYColumn, gazeLXColumn, gazeLYColumn, leftPupilColumn, rightPupilColumn;
		int counter = 0;
		for (const auto& t : tokens) {
			if (t.compare("Timestamp") == 0) {
				timestampColumn = counter;
			}
			else if (t.compare("LEFT_GAZE_X") == 0) {
				gazeLXColumn = counter;
			}
			else if (t.compare("LEFT_GAZE_Y") == 0) {
				gazeLYColumn = counter;
			}
			else if (t.compare("RIGHT_GAZE_X") == 0) {
				gazeRXColumn = counter;
			}
			else if (t.compare("RIGHT_GAZE_Y") == 0) {
				gazeRYColumn = counter;
			}
			else if (t.compare("LEFT_PUPIL_DIAMETER") == 0) {
				leftPupilColumn = counter;
			}
			else if (t.compare("RIGHT_PUPIL_DIAMETER") == 0) {
				rightPupilColumn = counter;
			}
			counter++;
		}

		//loop looking for tasks
		//for each task, read the relevant data and generate the plot

		while (!input.eof()) {
			getline(input, line);
			trim(line);
			if (line.find("REPEATING LAST TEST") != line.npos) {
				if (!results.empty()) {
					results.pop_back();
					events.pop_back();
					taskTypes.pop_back();
				}
				continue;

			}
			std::size_t pos = line.find("START");

			if (pos != line.npos) { //start of a trial
				std::vector<double> timestamps;
				std::vector<double> glx, grx, gly, gry, lp, rp;
				glx.clear();
				grx.clear();
				gly.clear();
				gry.clear();
				lp.clear();
				rp.clear();
				timestamps.clear();
				std::size_t commaPos = line.find(",");
				double t0 = atof(line.substr(0, commaPos).c_str());
				std::string s = line;
				double softT0 = std::stod(s.substr(s.find_last_of(':') + 1)) / 1000.0;
				Analysis::TaskDataType data;
				Analysis::TaskEventsType taskEvents;
				std::string taskType;
				getTaskTypeAndParameters(line, taskType, data);

				taskTypes.push_back(taskType);
				bool taskEnd = false;
				while (!taskEnd) {
					double offsetTSInTrial = 0.0;
					double offsetTS = 0.0;
					getline(input, line);
					trim(line);
					std::size_t pos = line.find("END");
					if (pos != line.npos) { //end of a trial
						taskEnd = true;
						continue;
					}
					char_separator<char> sep(", ");
					tokenizer<char_separator<char>> tokens(line, sep);
					int counter = 0;
					bool goodLine = true;
					double timestamp;
					double a, b, c, d, e, f;
					for (const auto& t : tokens) {
						double val;
						char* strend;
						val = std::strtod(t.c_str(), &strend);
						if (counter == timestampColumn) {
							timestamp = (val - t0) / 1000.0;
						}
						else if (val == 0 & strend == t.c_str()) {//it's an  event
							goodLine = false;
							Analysis::EventType ev;

							ev.timestamp = timestamp;
							std::vector<std::string> tok;
							std::copy(tokens.begin(), tokens.end(), std::back_inserter<std::vector<std::string> >(tok));
							tok.erase(tok.begin());
							std::stringstream imploded;
							std::copy(tok.begin(), tok.end(), std::ostream_iterator<std::string>(imploded, " "));

							ev.message = imploded.str();
							std::string st = imploded.str();
							double softT = std::stod(st.substr(st.find_last_of(':') + 1)) / 1000.0;
							ev.inTrialTS = (softT - softT0) / 1000.0;
							if (taskEvents.size() > 0 && ev.inTrialTS < taskEvents.back().inTrialTS) {
								//need to correct the timestamp
								ev.inTrialTS = (double(std::numeric_limits<unsigned long>::max()) / 1000.0 + softT - softT0) / 1000.0;
							}
							taskEvents.push_back(ev);
							break;
						}
						else { //it's a data sample
							if (counter == gazeLXColumn) {
								a = (val);
							}
							else if (counter == gazeLYColumn) {
								b = (val);
							}
							else if (counter == gazeRXColumn) {
								c = (val);
							}
							else if (counter == gazeRYColumn) {
								d = (val);
							}
							else if (counter == leftPupilColumn) {
								e = (val);
							}
							else if (counter == rightPupilColumn) {
								f = (val);
							}
						}

						counter++;

					}
					if (goodLine) {
						timestamps.push_back(timestamp);
						(*dvaConversion)(a, b);
						(*dvaConversion)(c, d);
						glx.push_back(a);
						gly.push_back(b);
						grx.push_back(c);
						gry.push_back(d);
						lp.push_back(e);
						rp.push_back(f);

					}


				}
				data["t"] = timestamps;

				data["glx"] = glx;
				data["gly"] = gly;
				data["grx"] = grx;
				data["gry"] = gry;
				data["lp"] = lp;
				data["rp"] = rp;
				results.push_back(data);
				events.push_back(taskEvents);
			}
		}
	}

	DataFile::~DataFile() = default;


}