#ifndef DATAFILEANALYSIS_H
#define DATAFILEANALYSIS_H

#include "Analysis/DataTypes.h"
#include <memory>

namespace Analysis {
	class DllExport DataFile {
	public:
		DataFile(const std::string& filename);
		bool Process();
		bool isHDFFile();

	private:
		virtual ~DataFile();
		void getTaskTypeAndParameters(const std::string& line, std::string& taskType, TaskDataType& parameters);
		std::string filename;
		bool ProcessTextFile();
		bool ProcessHDFFile();
		bool ProcessExperimentalSetup(const std::string& setup);
		bool ProcessRecordingInfo(const std::string& info);
		void readTextData(std::istream& input, std::vector<std::string>& taskTypes, ExperimentResultsType& results, ExperimentEventsType& events);

		std::unique_ptr<SampleToDva> dvaConversion;
		ExperimentResultsType results;
		ExperimentEventsType events;
		std::vector<std::string> taskNames;
		std::string testDate;
		std::string testNotes;
		std::string subjectID;
		std::string subjectDOB;

	};
}

#endif
