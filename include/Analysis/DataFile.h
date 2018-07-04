#ifndef DATAFILEANALYSIS_H
#define DATAFILEANALYSIS_H

#include "Analysis/DataTypes.h"
#include <memory>

namespace H5 {
	class H5File;
	class Group;
}

namespace Analysis {
	namespace IO {
		class DllExport SegmentedDataFile;

		typedef enum:short{ RAW_FILE = 0x00, ANALYSIS_FILE = 0x01 } FileTypeHDF;
		class DllExport RawDataFile {
		public:
			static RawDataFile* Read(const std::string& filename);
			void GetRawData(std::vector<std::string>&, ExperimentResultsType&, ExperimentEventsType&);
			void GetRecordingInfo(std::string& testDate, std::string& testNotes, std::string& subjectID, std::string& subjectDOB);
			bool ExportAsHDF5(const std::string& filename);
			bool ExportAsText(const std::string& filename);
			SegmentedDataFile* Segment(const std::string&);
			virtual ~RawDataFile();
			static bool isRawHDFFile(const std::string& filename);
			static bool isRawTextFile(const std::string& filename);

		private:
			RawDataFile();
			typedef enum { HORIZONTAL, VERTICAL, CONVERGENCE } SaccadeTaskType;
			typedef struct {
				double timestamp;
				double xPosition;
				double yPosition;
				double zPosition;
			} DotPositionType;
			typedef struct {
				double timestamp;
				double lx;
				double ly;
				double rx;
				double ry;
				double lp;
				double rp;
			} GazeSampleType;
			static std::vector<DotPositionType> extractSaccadePositionsFromEvents(const TaskEventsType&);
			static void GetTrialGroupAndDatasetNames(const H5::Group& grp, std::string& groupName, std::string& dsName);
			static bool CopyGazeDataToDataset(const TaskDataType& data, const std::string& datasetName, H5::Group& location, bool compression = true);
			static bool WriteDotPositions(const std::vector<DotPositionType>& dotPositions, const std::string& dsName, H5::Group& trialGroup);
			H5::Group OpenGroup(const std::string& name, H5::Group& baseGroup);
			SaccadeTaskType RawDataFile::getSaccadeTaskType(const Analysis::TaskEventsType& events) const;
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

		static const struct {
			struct {
				const char* subjectID = "Subject ID";
				const char* subjectDOB = "Subject DOB";
				const char* recordingDate = "Recording Date";
				const char* recordingNotes = "Recording Notes";
			} RootAttributes;
			
			const char* Tasks = "Tasks";
			const char* TrialPrefix = "Trial";
			const char* Analysis = "Analysis";
			const char* Figures = "Figures";

			struct {
				const char* Fixation = "Fixation";
				const char* HorizontalSaccades = "Horizontal Saccades";
				const char* VerticalSaccades = "Vertical Saccades";
				const char* SmoothConvergence = "Smooth Convergence";
				const char* SmoothPursuit = "Smooth Pursuit";
				const char* OKN = "OKN";
				const char* Pupillometry = "Pupillometry";
			} TaskNames;
			struct {
				const char* StartingDotPosition = "Starting Dot Position";
				const char* Velocity = "Velocity";
				const char* Amplitude = "Amplitude";
				const char* Orientation = "Orientation";
				const char* Horizontal = "Horizontal";
				const char* Vertical = "Vertical";
				const char* Direction = "Direction";
			} TaskAttributes;
			struct {
				const char* timestamp = "Timestamp";
				const char* X = "X";
				const char* Y = "Y";
				const char* Z = "Z";

			} DotPositionColumnNames;

			struct {
				const char* timestamp = "Timestamp";
				const char* LeftX = "LeftX";
				const char* LeftY = "LeftYX";
				const char* RightX = "RightX";
				const char* RightY = "RightY";
				const char* LeftPupil = "Left pupil";
				const char* RightPupil = "Right pupil";
			} GazePositionColumnNames;
			
		} SegmentedFileNames;
		
		class DllExport SegmentedDataFile {
		public:
			
			static SegmentedDataFile* Read(const std::string& filename);
			bool DoAnalyses();
			
			virtual ~SegmentedDataFile();
			static bool IsSegmentedHDFFile(const std::string&);
		private:
			static const std::vector<std::string >  EyePositionsTableLabels;
			std::unique_ptr<H5::H5File> hdfFile;
			SegmentedDataFile();
			
		};
	}
}

#endif
