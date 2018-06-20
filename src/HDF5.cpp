#define EXPORTING
#include "Analysis/HDF5.h"
#include "H5Cpp.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <cstdio>

using namespace H5;

namespace Analysis {

	bool HDF5DataFile::IsValidHDF5(const std::string& file) {
		bool isH5;
		try {
			isH5 = H5::H5File::isHdf5(file);
		}
		catch (H5::FileIException& e) {
			isH5 = false;
		}
		return isH5;
	}

	bool HDF5DataFile::Text2HDF5(const std::string& origin, const std::string& dest) {
		std::ifstream inFile;
		try {
			H5File destFile(dest, H5F_ACC_TRUNC);


			StrType stringType(0, H5T_VARIABLE);
			CompType headerType(sizeof(SimpleTableType));
			headerType.insertMember("Title", HOFFSET(SimpleTableType, title), stringType);
			headerType.insertMember("Content", HOFFSET(SimpleTableType, content), stringType);
	
			CompType eventType(sizeof(EventType));
			eventType.insertMember("Timestamp", HOFFSET(EventType, timestamp), H5T_NATIVE_DOUBLE);
			eventType.insertMember("Text", HOFFSET(EventType, text), stringType);
	
			
	
			
			Group dataGroup = destFile.createGroup("Data");
	
			
			inFile.open(origin);
	
			std::string line;
			std::getline(inFile, line);
			std::vector<SimpleTableType> headers;
			std::vector<EventType> events;
			std::vector<std::vector<double> > data;
			std::vector<std::string> fieldNames;
		
			do {

				if (line.compare("[START DATA]") == 0) {

					std::getline(inFile, line);
					boost::trim(line);
					boost::char_separator<char> sep(", ");
					boost::tokenizer<boost::char_separator<char>> tokens(line, sep);


					for (const auto& t : tokens) {
						fieldNames.push_back(t);
					}

					data.resize(fieldNames.size());
					std::getline(inFile, line);
					boost::trim(line);
					while (!inFile.eof()) {

						//decide if line is event or data
						double t;
						boost::char_separator<char> sepToks(",");
						boost::tokenizer<boost::char_separator<char>> datTokens(line, sepToks);
						std::vector<std::string> toks;
						for (const auto& tk : datTokens) {
							toks.push_back(tk);
						}

						std::istringstream iss(toks[1]);
						double x;
						iss >> x;
						bool isEvent = iss.fail();

						if (isEvent) {
							EventType e;
							e.timestamp = std::stod(toks[0]);
							std::string s("");
							for (std::size_t k = 1; k < toks.size(); k++) {
								s += toks[k];
							}

							e.text = new char[s.size() + 1];
							strncpy(e.text, s.c_str(), s.size() + 1);
							events.push_back(e);
						}
						else {
							for (std::size_t k = 0; k < data.size(); k++) {
								data[k].push_back(std::stod(toks[k]));
							}
						}
						std::getline(inFile, line);
						boost::trim(line);
					}
				}
				else {
					if (line[0] == '[') {
						SimpleTableType t;
						std::string title(line);
						boost::erase_all(title, "[");
						boost::erase_all(title, "]");
						boost::erase_all(title, "\r");
						boost::erase_all(title, "\n");

						t.title = new char[title.size() + 1];
						strcpy_s(t.title, title.size() + 1, title.data());

						std::stringstream content;
						std::getline(inFile, line);
						boost::algorithm::trim(line);
						while (line[0] != '[') {

							if (line.size() != 0) {
								content << line << std::endl; //do i need to append the endline??
							}
							std::getline(inFile, line);
							boost::algorithm::trim(line);
						}

						t.content = new char[content.str().size() + 1];
						strcpy_s(t.content, content.str().size() + 1, content.str().data());

						headers.push_back(t);



						continue;
					}
				}

				std::getline(inFile, line);
				boost::trim(line);
			} while (!inFile.eof());

			hsize_t headerLength = headers.size();
			hsize_t headerRank = 1;
			hsize_t headerDim[1] = { headerLength };
			hsize_t eventsLength = events.size();
			hsize_t eventsRank = 1;
			hsize_t eventsDim[1] = { eventsLength };
			hsize_t dataLength = data[0].size();
			hsize_t dataRank = 1;
			hsize_t dataDim[1] = { dataLength };

			DataSpace eventsSpace(eventsRank, eventsDim);
			DataSet eventsDS = destFile.createDataSet("Events", eventType, eventsSpace);
			eventsDS.write(events.data(), eventType);

			DataSpace dataSpace(dataRank, dataDim);
			hsize_t chunk_dims[1] = { 512 };
			DSetCreatPropList *props = new DSetCreatPropList;
			props->setChunk(1, chunk_dims);
			props->setDeflate(9);
			for (std::size_t k = 0; k < data.size(); k++) {

				DataSet ds = dataGroup.createDataSet(fieldNames[k], H5T_NATIVE_DOUBLE, dataSpace, *props);
				ds.write(data[k].data(), H5T_NATIVE_DOUBLE);

			}
			DataSpace headerSpace(headerRank, headerDim);
			DataSet headerDS = destFile.createDataSet("Header", headerType, headerSpace);
			headerDS.write(headers.data(), headerType);
			//std::for_each(headers.begin(), headers.end(), [] (SimpleTableType& v) { delete[] v.title; delete[] v.content; });
			headerSpace.close();
			headerDS.close();
			destFile.close();

		}
		catch (...) { //this is horrible, but it works
			inFile.close();
			std::remove(dest.c_str());
			return false;
		}
		


		return true;
	}

	bool HDF5DataFile::HDF52Text(const std::string& origin, const std::string& dest) {

		std::ofstream destFile(dest, std::ios::trunc);

		H5File inFile(origin.c_str(), H5F_ACC_RDONLY);
		DataSet headers = inFile.openDataSet("Header");

		CompType datasetType(headers);
		DataSpace headersSpace = headers.getSpace();

		std::size_t numHeaders = headersSpace.getSimpleExtentNdims()*headersSpace.getSelectNpoints();
		std::size_t headersBufferSize = numHeaders*datasetType.getNmembers();

		std::vector<char*> rdata(headersBufferSize);

		headers.read((void*)rdata.data(), datasetType);

		for (unsigned i = 0; i < numHeaders; i++) {
			destFile << "[" << rdata[i * 2] << "]" << std::endl;
			destFile << rdata[i * 2 + 1] << std::endl;
			delete[] rdata[i * 2];
			delete[] rdata[i * 2 + 1];

		} /* end for */

		headers.close();

		destFile << "[START DATA]" << std::endl;
		DataSet eventsSet = inFile.openDataSet("Events");
		DataSpace eventsSpace = eventsSet.getSpace();
		CompType eventType(eventsSet);
		std::size_t numEvents = eventsSpace.getSimpleExtentNdims()*eventsSpace.getSelectNpoints();

		std::vector<EventType> eventData(numEvents);
		eventsSet.read((void*)eventData.data(), eventType);




		Group dataGroup = inFile.openGroup("Data");
		std::size_t numObjects = dataGroup.getNumObjs();
		std::vector<std::vector<double>> data(1);
		std::vector<std::string> columnNames(1);
		columnNames[0] = "Timestamp";
		DataSet timestamps = dataGroup.openDataSet("Timestamp"); //there needs to be a timestamp dataset
		std::size_t numSamples = timestamps.getSpace().getSelectNpoints();
		data[0].resize(numSamples);
		timestamps.read((void*)data[0].data(), H5T_NATIVE_DOUBLE);
		std::size_t columnIndex = 1;
		for (hsize_t k = 0; k < numObjects; k++) {
			H5G_obj_t objType = dataGroup.getObjTypeByIdx(k);
			if (objType == H5G_DATASET) {
				std::string datasetName = dataGroup.getObjnameByIdx(k);
				DataSet ds = dataGroup.openDataSet(datasetName);
				if (datasetName.compare("Timestamp") == 0) continue;
				columnNames.push_back(datasetName);
				std::vector<double> samples(numSamples);
				ds.read(samples.data(), H5T_NATIVE_DOUBLE);
				data.push_back(samples);
				columnIndex++;
			}
		}

		for (std::size_t k = 0; k < columnNames.size(); k++) {
			if (k > 0) destFile << ", ";
			destFile << columnNames[k];
		}
		destFile << std::endl;
		std::size_t eventCount = 0;
		std::streamsize precisionOrig = destFile.precision();
		destFile.precision(24);
		for (std::size_t k = 0; k < data[0].size(); k++) {
			if (data[0][k] >= eventData[eventCount].timestamp) {
				destFile << eventData[eventCount].timestamp << "," << eventData[eventCount].text << std::endl;
				delete[] eventData[eventCount].text;
				eventCount++;
			}
			destFile << data[0][k];
			for (std::size_t n = 1; n < data.size(); n++) {
				destFile << ",";
				destFile << data[n][k];
			}
			destFile << std::endl;

		}



		destFile.close();
		inFile.close();
		return true;

	}

	HDF5DataFile::HDF5DataFile(const std::string& filename) {
		if (!IsValidHDF5(filename)) return;
		file = std::make_unique<H5File>(filename.c_str(), H5F_ACC_RDONLY);
		dataGroup = std::make_unique<Group>(file->openGroup("Data"));
	}

	HDF5DataFile::~HDF5DataFile() {
		if (file != nullptr) {
			dataGroup->close();
			file->close();
		}


	}

	std::string HDF5DataFile::GetConfiguration() {
		return GetHeaderContent("CONFIGURATION");

	}
	std::string HDF5DataFile::GetTasks() {
		return GetHeaderContent("TASKS");

	}
	std::string HDF5DataFile::GetRecordingInfo() {
		return GetHeaderContent("DATA HEADER");
	}

	std::string HDF5DataFile::GetHeaderContent(const std::string& headerName) {
		DataSet headers = file->openDataSet("Header");

		CompType datasetType(headers);
		DataSpace headersSpace = headers.getSpace();

		std::size_t numHeaders = headersSpace.getSimpleExtentNdims()*headersSpace.getSelectNpoints();
		std::size_t headersBufferSize = numHeaders*datasetType.getNmembers();

		std::vector<char*> rdata(headersBufferSize);

		headers.read((void*)rdata.data(), datasetType);
		std::string result;
		for (unsigned i = 0; i < numHeaders; i++) {
			if (strcmp(rdata[i * 2], headerName.c_str()) == 0) {
				result = rdata[i * 2 + 1];
			}
		} /* end for */

		headers.close();
		return result;

	}

	void HDF5DataFile::GetEvents(std::vector<double>& t , std::vector<std::string>& text) {
		DataSet eventsSet = file->openDataSet("Events");
		DataSpace eventsSpace = eventsSet.getSpace();
		CompType eventType(eventsSet);
		std::size_t numEvents = eventsSpace.getSimpleExtentNdims()*eventsSpace.getSelectNpoints();

		std::vector<EventType> eventData(numEvents);
		eventsSet.read((void*)eventData.data(), eventType);
		t.resize(eventData.size());
		text.resize(eventData.size());
		for (std::size_t k = 0; k < eventData.size(); k++) {
			t[k] = eventData[k].timestamp;
			text[k] = eventData[k].text;
			delete[] eventData[k].text;
		}
	}

	void HDF5DataFile::GetDataWithName(const std::string& name, std::vector<double>& out) {
	
		std::size_t numObjects = dataGroup->getNumObjs();
		DataSet timestamps = dataGroup->openDataSet("Timestamp"); //there needs to be a timestamp dataset
		std::size_t numSamples = timestamps.getSpace().getSelectNpoints();
		timestamps.close();
		out.resize(numSamples);
		for (hsize_t k = 0; k < numObjects; k++) {
			H5G_obj_t objType = dataGroup->getObjTypeByIdx(k);
			if (objType == H5G_DATASET) {
				std::string datasetName = dataGroup->getObjnameByIdx(k);
				if (datasetName.compare(name) != 0) continue;
				DataSet ds = dataGroup->openDataSet(datasetName);
				ds.read(out.data(), H5T_NATIVE_DOUBLE);
				ds.close();
				break;
			}
		}
	}

	bool HDF5DataFile::HasDataWithName(const std::string& name) {
		std::size_t numObjects = dataGroup->getNumObjs();
		
		for (hsize_t k = 0; k < numObjects; k++) {
			H5G_obj_t objType = dataGroup->getObjTypeByIdx(k);
			if (objType == H5G_DATASET) {
				std::string datasetName = dataGroup->getObjnameByIdx(k);
				if (datasetName.compare(name) == 0) return true;
			}
		}
		return false;
	}

	void HDF5DataFile::GetTimestamps(std::vector<double>& out) {
		
		DataSet timestamps = dataGroup->openDataSet("Timestamp"); //there needs to be a timestamp dataset
		std::size_t numSamples = timestamps.getSpace().getSelectNpoints();
		out.resize(numSamples);
		timestamps.read((void*)out.data(), H5T_NATIVE_DOUBLE);
	}
		
		
	bool HDF5DataFile::HasLeftGaze() {
		return HasDataWithName("LEFT_GAZE_X");

	}
	void HDF5DataFile::GetLeftGaze(std::vector<double>& lx, std::vector<double>& ly) {
		GetDataWithName("LEFT_GAZE_X", lx);
		GetDataWithName("LEFT_GAZE_Y", ly);
	}
	bool HDF5DataFile::HasRightGaze() {
		return HasDataWithName("RIGHT_GAZE_X");
	}
	void HDF5DataFile::GetRightGaze(std::vector<double>& rx, std::vector<double>& ry) {
		GetDataWithName("RIGHT_GAZE_X", rx);
		GetDataWithName("RIGHT_GAZE_Y", ry);
	}
	bool HDF5DataFile::HasLeftPupilSize() {
		return HasDataWithName("LEFT_PUPIL_DIAMETER");
	}
	void HDF5DataFile::GetLeftPupilSize(std::vector<double>& lp) {
		GetDataWithName("LEFT_PUPIL_DIAMETER", lp);
	}
	bool HDF5DataFile::HasRightPupilSize(){
		return HasDataWithName("RIGHT_PUPIL_DIAMETER");
	}

	void HDF5DataFile::GetRightPupilSize(std::vector<double>& rp) {
		GetDataWithName("RIGHT_PUPIL_DIAMETER", rp);
	}

}
