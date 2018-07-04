#define EXPORTING
#include "json/value.h"
#include "json/reader.h"
#include "Analysis/DataFile.h"
#include "Analysis/DataTypes.h"
#include "Analysis/HDF5.h"
#include <H5Cpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <memory>


#define M_PI 3.1415926

namespace Analysis {
	namespace IO {

		

		class EyeTechDSSampleToDVA : public SampleToDva {
		public:
			EyeTechDSSampleToDVA(int _pxHoriz, int _pxVert, double _distance, double _dimHoriz, double _dimVert) {
				pxHoriz = _pxHoriz;
				pxVert = _pxVert;
				distance = _distance;
				dimHoriz = _dimHoriz;
				dimVert = _dimVert;
			}
			void convert(double& x, double& y) override {
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

			void SMIETHMDSampleToDVA::convert(double& x, double& y) override {
				x = atan2(x - 1080.0, virtualDistanceX) / M_PI*180.0;
				y = -atan2(y - 600.0, virtualDistanceY) / M_PI*180.0;
			}
			virtual ~SMIETHMDSampleToDVA() = default;
		private:
			const double SMIETHMDSampleToDVA::virtualDistanceX = 1080.0 / tan(111.9 / 2.0 / 180 * M_PI);
			const double SMIETHMDSampleToDVA::virtualDistanceY = 600.0 / tan(105.6 / 2.0 / 180 * M_PI);
		};

		bool RawDataFile::isRawTextFile(const std::string& filename) {
			return true;// do something here
		}

		RawDataFile::RawDataFile() = default;

		RawDataFile* RawDataFile::Read(const std::string& filename) {
			
			RawDataFile* file = new RawDataFile();
			file->filename = filename;

			bool success = false;
			file->taskNames.clear();
			file->results.clear();
			file->events.clear();
			try {
				if (isRawHDFFile(filename)) {
					success = file->ProcessHDFFile();
				}
				else if (isRawTextFile(filename)) {
					success = file->ProcessTextFile();
				}
				else {
					delete file;
					return nullptr;
				}
			}
			catch (...) {//ugly
				delete file;
				return nullptr;
			}
			
			return file;
		}

		bool RawDataFile::isRawHDFFile(const std::string& filename) {
			if (!Analysis::IO::IsValidHDF5(filename)) return false;
			return true;
		}

		bool RawDataFile::ProcessHDFFile() {

			//read necessary information from headers
			HDF5RawDataFile file(filename);
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
							Analysis::H5EventType ev;
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
						dvaConversion->convert(x, y);
						glx[n] = x;
						gly[n] = y;
						x = recordingRX[k];
						y = recordingRY[k];
						dvaConversion->convert(x, y);
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

		bool RawDataFile::ProcessTextFile() {
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
						stringval.append("\n");
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

		void RawDataFile::getTaskTypeAndParameters(const std::string& line, std::string& taskType, Analysis::TaskDataType& parameters) {

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

		bool RawDataFile::ProcessExperimentalSetup(const std::string& stringval) {

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

		bool RawDataFile::ProcessRecordingInfo(const std::string& info) {
			if (info.empty()) return false;
			std::stringstream infoStream(info);
			std::string line;

			using namespace boost;
			do {
				
				std::getline(infoStream, line);
				boost::algorithm::trim(line);
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
			} while (!infoStream.eof());
			 

			return true;

		}

		void RawDataFile::readTextData(std::istream& input, std::vector<std::string>& taskTypes, Analysis::ExperimentResultsType& results, Analysis::ExperimentEventsType& events) {
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
								Analysis::H5EventType ev;

								ev.timestamp = timestamp;
								std::vector<std::string> tok;
								std::copy(tokens.begin(), tokens.end(), std::back_inserter<std::vector<std::string> >(tok));
								tok.erase(tok.begin());
								std::stringstream imploded;
								std::copy(tok.begin(), tok.end(), std::ostream_iterator<std::string>(imploded, " "));

								ev.message = imploded.str();
								std::string st = imploded.str();
								double softT = std::stod(st.substr(st.find_last_of(':') + 1)) / 1000.0;
								//TODO: work on this part, it is messy
								ev.inTrialTS = (softT - softT0)/1000.0;
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
							dvaConversion->convert(a, b);
							dvaConversion->convert(c, d);
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

		void RawDataFile::GetRawData(std::vector<std::string>& taskNames, ExperimentResultsType& data, ExperimentEventsType& events) {

			taskNames = this->taskNames;
			data = this->results;
			events = this->events;
		}
		void RawDataFile::GetRecordingInfo(std::string& testDate, std::string& testNotes, std::string& subjectID, std::string& subjectDOB) {
			testDate = this->testDate;
			testNotes = this->testNotes;
			subjectID = this->subjectID;
			subjectDOB = this->subjectDOB;
		}
		bool RawDataFile::ExportAsHDF5(const std::string& filename) {
			return false;

		}
		bool RawDataFile::ExportAsText(const std::string& filename) {
			return false;
		}
		
		RawDataFile::~RawDataFile() = default;

		const std::vector<std::string> SegmentedDataFile::EyePositionsTableLabels = { "Timestamp", "Left Eye X", "Left Eye Y", "Right Eye X", "Right Eye Y", "Left Pupil Diameter", "Right Pupil Diameter" };
		
		
		SegmentedDataFile* RawDataFile::Segment(const std::string& output) {

			SegmentedDataFile* outFile = nullptr;
			using namespace H5;

			H5File h5f(output, H5F_ACC_TRUNC);

			DataSpace fileTypeAttrSpace (H5S_SCALAR);
			DataType fileTypeAttrType(PredType::STD_B8LE);
			Attribute fileTypeAttr = h5f.createAttribute("File Type", fileTypeAttrType, fileTypeAttrSpace);
			FileTypeHDF fileTypeAttrVal = ANALYSIS_FILE;
			fileTypeAttr.write(fileTypeAttrType, &fileTypeAttrVal);
			StrType stringType(0, H5T_VARIABLE);

			hsize_t dims1[] = { 1 };
			DataSpace idAttrSpace(1, dims1);
			Attribute idAttr = h5f.createAttribute(SegmentedFileNames.RootAttributes.subjectID, stringType, idAttrSpace);
			idAttr.write(stringType, this->subjectID);
			DataSpace dobAttrSpace(1, dims1);
			Attribute dobAttr = h5f.createAttribute(SegmentedFileNames.RootAttributes.subjectDOB, stringType, dobAttrSpace);
			dobAttr.write(stringType, this->subjectDOB);
			DataSpace dateAttrSpace(1, dims1);
			Attribute dateAttr = h5f.createAttribute(SegmentedFileNames.RootAttributes.recordingDate, stringType, dateAttrSpace);
			dateAttr.write(stringType, this->testDate);
			DataSpace notesAttrSpace(1, dims1);
			Attribute notesAttr = h5f.createAttribute(SegmentedFileNames.RootAttributes.recordingNotes, stringType, notesAttrSpace);
			notesAttr.write(stringType, this->testNotes);

			Group tasksGroup = h5f.createGroup(SegmentedFileNames.Tasks);
			ExperimentResultsType fixationData;
			ExperimentEventsType fixationEvents;
			std::vector<std::string> dotPosColNames;
			dotPosColNames.push_back(SegmentedFileNames.DotPositionColumnNames.timestamp);
			dotPosColNames.push_back(SegmentedFileNames.DotPositionColumnNames.X);
			dotPosColNames.push_back(SegmentedFileNames.DotPositionColumnNames.Y);
			dotPosColNames.push_back(SegmentedFileNames.DotPositionColumnNames.Z);
			std::vector<std::string> gazePosColNames;
			gazePosColNames.push_back(SegmentedFileNames.GazePositionColumnNames.timestamp);
			gazePosColNames.push_back(SegmentedFileNames.GazePositionColumnNames.LeftX);
			gazePosColNames.push_back(SegmentedFileNames.GazePositionColumnNames.LeftY);
			gazePosColNames.push_back(SegmentedFileNames.GazePositionColumnNames.RightX);
			gazePosColNames.push_back(SegmentedFileNames.GazePositionColumnNames.RightY);
			gazePosColNames.push_back(SegmentedFileNames.GazePositionColumnNames.LeftPupil);
			gazePosColNames.push_back(SegmentedFileNames.GazePositionColumnNames.RightPupil);
			
			Group fixationGroup;
			Group horizontalSaccadesGroup;
			Group verticalSaccadesGroup;
			Group smoothPursuitGroup;
			Group smoothConvergenceGroup;
			Group stepConvergenceGroup;
			Group oknGroup;
			Group pupillometryGroup;

			for (std::size_t k = 0; k < taskNames.size(); k++) {
				std::string groupName, dsName;
				Group grp, trialGroup;
				std::vector<DotPositionType> dotPositions;
				if (std::strcmp(taskNames[k].c_str(), "FIXATION") == 0) {
					grp = OpenGroup(SegmentedFileNames.TaskNames.Fixation, tasksGroup);
					
					GetTrialGroupAndDatasetNames(grp, groupName, dsName);
					trialGroup = grp.createGroup(groupName);
					DotPositionType dp;
					dp.timestamp = 0;
					dp.xPosition = 0;
					dp.yPosition = 0;
					dp.zPosition = 0;
					dotPositions.push_back(dp);
				}
				else if (std::strcmp(taskNames[k].c_str(), "GUIDED SACCADE") == 0) {
					dotPositions = extractSaccadePositionsFromEvents(events[k]);
					std::string saccadeTask;
					switch (getSaccadeTaskType(events[k])) {
					case HORIZONTAL: 
						saccadeTask = SegmentedFileNames.TaskNames.HorizontalSaccades;
						break;
					case VERTICAL:
						saccadeTask = SegmentedFileNames.TaskNames.VerticalSaccades;
						break;
					case CONVERGENCE:
						continue;
						break;
					}
					grp = OpenGroup(saccadeTask, tasksGroup);
					GetTrialGroupAndDatasetNames(grp, groupName, dsName);
					dsName += groupName;
					trialGroup = grp.createGroup(groupName);
				}
				else if (std::strcmp(taskNames[k].c_str(), "PURSUIT") == 0) {
					if (std::abs(results[k]["motionParameters"][2] - 2.0) < 0.1) { //smooth convergence
						grp = OpenGroup(SegmentedFileNames.TaskNames.SmoothConvergence, tasksGroup);
					
						GetTrialGroupAndDatasetNames(grp, groupName, dsName);
						trialGroup = grp.createGroup(groupName);

						DotPositionType dp;
						dp.timestamp = 0;
						dp.xPosition = 0;
						dp.yPosition = 0;
						dp.zPosition = results[k]["motionParameters"].back();
						dotPositions.push_back(dp);
						for (std::size_t j = 0; j < events[k].size(); j++) {
							DotPositionType dcp;
							dcp.timestamp = events[k][j].inTrialTS*1000.0;
							std::size_t colonPos = events[k][j].message.find_first_of(':');
							std::size_t barPos = events[k][j].message.find_first_of("|||");
							double displ = std::stod(events[k][j].message.substr(colonPos + 1, barPos - colonPos - 1).c_str());
							dcp.xPosition = 0;
							dcp.yPosition = 0;
							dcp.zPosition = dp.zPosition + displ;
							dotPositions.push_back(dcp);
						}
					}
					else {

						grp = OpenGroup(SegmentedFileNames.TaskNames.SmoothPursuit, tasksGroup);
						GetTrialGroupAndDatasetNames(grp, groupName, dsName);
						trialGroup = grp.createGroup(groupName);

						DotPositionType dp;
						dp.timestamp = 0;
						dp.xPosition = 0;
						dp.yPosition = 0;
						dp.zPosition = 0;
						dotPositions.push_back(dp);
						
					}
					DataType velType(H5T_NATIVE_DOUBLE);
					DataSpace velSpace(H5S_SCALAR);
					Attribute velAttr = trialGroup.createAttribute("Velocity", velType, velSpace);
					velAttr.write(velType, &results[k]["motionParameters"][1]);

				}
				else if (std::strcmp(taskNames[k].c_str(), "PUPIL") == 0) {
					grp = OpenGroup(SegmentedFileNames.TaskNames.Pupillometry, tasksGroup);
					GetTrialGroupAndDatasetNames(grp, groupName, dsName);
					trialGroup = grp.createGroup(groupName);
					DotPositionType dp;
					dp.timestamp = events[k][0].inTrialTS*1000.0;
					dotPositions.push_back(dp);
					dp.timestamp = events[k][1].inTrialTS*1000.0;
					dotPositions.push_back(dp);

					
				}
				else if (std::strcmp(taskNames[k].c_str(), "OKN") == 0) {
					grp = OpenGroup(SegmentedFileNames.TaskNames.OKN, tasksGroup);
					GetTrialGroupAndDatasetNames(grp, groupName, dsName);
					trialGroup = grp.createGroup(groupName);
					DataType velType(H5T_NATIVE_DOUBLE);
					DataSpace velSpace(H5S_SCALAR);
					Attribute velAttr = trialGroup.createAttribute("Velocity", velType, velSpace);
					velAttr.write(velType, &results[k]["speed"][0]);

					

				}
				else {
					continue;
				}
				WriteDotPositions(dotPositions, dsName, trialGroup);
				CopyGazeDataToDataset(results[k], dsName, trialGroup);
			}



			h5f.close();
			outFile = SegmentedDataFile::Read(output);
			return outFile;


		}

		bool RawDataFile::WriteDotPositions(const std::vector<DotPositionType>& dotPositions, const std::string& dsName, H5::Group& trialGroup) {
			
			hsize_t numRows = dotPositions.size();

			if (numRows == 0) return false;

			H5::DataType doubleType(H5T_NATIVE_DOUBLE);
			
			H5::CompType h5Type(sizeof(DotPositionType));
			h5Type.insertMember(SegmentedFileNames.DotPositionColumnNames.timestamp, HOFFSET(DotPositionType, timestamp), doubleType);
			h5Type.insertMember(SegmentedFileNames.DotPositionColumnNames.X, HOFFSET(DotPositionType, xPosition), doubleType);
			h5Type.insertMember(SegmentedFileNames.DotPositionColumnNames.Y, HOFFSET(DotPositionType, yPosition), doubleType);
			h5Type.insertMember(SegmentedFileNames.DotPositionColumnNames.Z, HOFFSET(DotPositionType, zPosition), doubleType);

			hsize_t dims[] = { numRows };
			H5::DataSpace space(1, dims);

			H5::DataSet ds = trialGroup.createDataSet((dsName+"_DotPositions").c_str(), h5Type, space);
			ds.write(dotPositions.data(), h5Type);
			ds.close();

			return true;
		}

		H5::Group RawDataFile::OpenGroup(const std::string& name, H5::Group& baseGroup) {
			H5::Group grp;
			try {
				grp = baseGroup.openGroup(name);
			}
			catch (H5::Exception& e) {
				grp = baseGroup.createGroup(name);
			}

			return grp;
		}

		void RawDataFile::GetTrialGroupAndDatasetNames(const H5::Group& grp, std::string& groupName, std::string& dsName) {
			hsize_t numTrials = grp.getNumObjs();
			
			groupName = SegmentedFileNames.TrialPrefix;
			groupName += "_";
			std::stringstream ss;
			ss.fill('0');
			ss.width(3);
			ss << numTrials + 1;
			groupName += ss.str();
			dsName = "Data_";
			dsName += groupName;
		}

		bool RawDataFile::CopyGazeDataToDataset(const TaskDataType& data, const std::string& datasetName, H5::Group& location, bool compression) {
			
			hsize_t numRows = data.at("t").size();
			std::vector<GazeSampleType> samples(numRows);
			for (std::size_t k = 0; k < numRows; k++) {
				GazeSampleType s;
				s.timestamp = data.at("t")[k];
				s.lx = data.at("glx")[k];
				s.ly = data.at("gly")[k];
				s.rx = data.at("grx")[k];
				s.ry = data.at("gry")[k];
				s.lp = data.at("lp")[k];
				s.rp = data.at("rp")[k];
				samples[k] = s;
			}
			
			H5::CompType h5Type(sizeof(GazeSampleType));
			H5::DataType doubleType(H5T_NATIVE_DOUBLE);
			h5Type.insertMember(SegmentedFileNames.GazePositionColumnNames.timestamp, HOFFSET(GazeSampleType, timestamp), doubleType);
			h5Type.insertMember(SegmentedFileNames.GazePositionColumnNames.LeftX, HOFFSET(GazeSampleType, lx), doubleType);
			h5Type.insertMember(SegmentedFileNames.GazePositionColumnNames.LeftY, HOFFSET(GazeSampleType, ly), doubleType);
			h5Type.insertMember(SegmentedFileNames.GazePositionColumnNames.RightX, HOFFSET(GazeSampleType, rx), doubleType);
			h5Type.insertMember(SegmentedFileNames.GazePositionColumnNames.RightY, HOFFSET(GazeSampleType, ry), doubleType);
			h5Type.insertMember(SegmentedFileNames.GazePositionColumnNames.LeftPupil, HOFFSET(GazeSampleType, lp), doubleType);
			h5Type.insertMember(SegmentedFileNames.GazePositionColumnNames.RightPupil, HOFFSET(GazeSampleType, rp), doubleType);

			hsize_t dims[] = { numRows };
			hsize_t chunk_dims[1] = { 512 };

			H5::DataSpace space(1, dims);
			std::unique_ptr<H5::DSetCreatPropList> props = std::make_unique<H5::DSetCreatPropList>();
			if (compression) {
				props->setChunk(1, chunk_dims);
				props->setDeflate(9);
			}
			
			H5::DataSet ds = location.createDataSet(datasetName.c_str(), h5Type, space, *props);
			ds.write(samples.data(), h5Type);
			ds.close();

			return true;
		}

		std::vector<RawDataFile::DotPositionType> RawDataFile::extractSaccadePositionsFromEvents(const TaskEventsType& events) {

			std::vector<DotPositionType> out;

			for (std::size_t k = 0; k < events.size(); k++) {
				DotPositionType dp;
				dp.timestamp = events[k].inTrialTS*1000.0;

				std::size_t pos = events[k].message.find(':');
				std::size_t closePos = events[k].message.find('|');
				std::string coord = events[k].message.substr(pos + 1, closePos - pos - 2);
				if (coord.find("CENTER") == std::string::npos) {
					std::istringstream ss(coord);
					ss >> dp.xPosition;
					ss >> dp.yPosition;
					ss >> dp.zPosition;
					dp.xPosition = atan2(dp.xPosition, 0.57)*180.0 / M_PI;
					dp.yPosition = atan2(dp.yPosition, 0.57)*180.0 / M_PI;
					dp.zPosition = -0.57;
				}
				else {
					dp.xPosition = 0.0;
					dp.yPosition = 0.0;
					dp.zPosition = -0.57;
				}
				out.push_back(dp);
			}
			return out;

		}

		RawDataFile::SaccadeTaskType RawDataFile::getSaccadeTaskType(const Analysis::TaskEventsType& events) const {

			std::vector<double> dotX(events.size());
			std::vector<double> dotY(events.size());
			std::vector<double> dotZ(events.size());
			bool xChanged = false;
			bool yChanged = false;
			bool zChanged = false;

			for (std::size_t k = 0; k < events.size(); k++) {
				std::size_t pos = events[k].message.find('z=');
				if (pos != std::string::npos) {
					std::size_t posClosing = events[k].message.find(')');
					std::string z = events[k].message.substr(pos + 2, posClosing - pos - 2);
					dotX[k] = 0.0;
					dotY[k] = 0.0;
					dotZ[k] = std::stod(z);
				}
				else {
					double xt, yt, zt;
					std::size_t pos = events[k].message.find(':');
					std::size_t closePos = events[k].message.find('|');
					std::string coord = events[k].message.substr(pos + 1, closePos - pos - 2);
					std::istringstream ss(coord);
					ss >> xt;
					ss >> yt;
					ss >> zt;
					dotX[k] = atan2(xt, 0.57)*180.0 / M_PI;
					dotY[k] = -atan2(yt, 0.57)*180.0 / M_PI;
					dotZ[k] = zt;

				}
				if (k > 0 && dotZ[k] != dotZ[k - 1]) {
					return CONVERGENCE;
				}
				else if (k > 0 && dotX[k] != dotX[k - 1]) {
					return HORIZONTAL;
				}
				else if (k > 0 && dotY[k] != dotY[k - 1]) {
					return VERTICAL;
				}
			}

		}

		SegmentedDataFile::SegmentedDataFile() {

		}

		SegmentedDataFile* SegmentedDataFile::Read(const std::string& filename) {
			if (!IsSegmentedHDFFile(filename)) return nullptr;
			SegmentedDataFile* f = new SegmentedDataFile;
			f->hdfFile = std::make_unique<H5::H5File>(filename, H5F_ACC_RDWR);
			return f;
		}

		bool SegmentedDataFile::DoAnalyses() {
			return true;
		}

		bool SegmentedDataFile::IsSegmentedHDFFile(const std::string& filename) {
			if (!IsValidHDF5(filename)) return false;
			H5::H5File f(filename, H5F_ACC_RDONLY);
			if (!f.attrExists("File type")) return false;
			H5::Attribute fileTypeAttr = f.openAttribute("File type");
			H5::DataType fileTypeDataType = fileTypeAttr.getDataType();

			short attrVal = 0;
			fileTypeAttr.read(fileTypeDataType, &attrVal);

			if (attrVal != FileTypeHDF::ANALYSIS_FILE) return false;
			
			return true;
		}

		SegmentedDataFile::~SegmentedDataFile() {
			hdfFile->close();
		}
		
	}

	


}