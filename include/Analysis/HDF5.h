#ifndef UTIL_HDF5_H
#define UTIL_HDF5_H

#if (defined (_MSC_VER))

#ifdef EXPORTING
#define DllExport __declspec(dllexport)
#else
#define DllExport __declspec(dllimport)
#endif

#elif (defined (__clang__))
#define DllExport
#endif


#include <string>
#include <vector>
#include <memory>



namespace H5 {
	class H5File;
	class Group;

}

namespace Analysis {

	namespace IO {

	typedef struct SimpleTableType_ {
		char* title;
		char* content;

	} H5SimpleTableType;

	typedef struct {
		double timestamp;
		char* text;
	} H5EventType;

	bool IsValidHDF5(const std::string& file);

	class DllExport HDF5RawDataFile {
	public:
		static bool Text2HDF5(const std::string& origin, const std::string& dest);
		static bool HDF52Text(const std::string& origin, const std::string& dest);
		
		HDF5RawDataFile(const std::string& filename);
		std::string GetConfiguration();
		std::string GetTasks();
		std::string GetRecordingInfo();
		void GetEvents(std::vector<double>&, std::vector<std::string>&);
		void GetTimestamps(std::vector<double>&);
		bool HasLeftGaze();
		void GetLeftGaze(std::vector<double>&, std::vector<double>&);
		bool HasRightGaze();
		void GetRightGaze(std::vector<double>&, std::vector<double>&);
		bool HasLeftPupilSize();
		void GetLeftPupilSize(std::vector<double>&);
		bool HasRightPupilSize();
		void GetRightPupilSize(std::vector<double>&);
		virtual ~HDF5RawDataFile();
		

	private:
		std::unique_ptr<H5::H5File> file;
		std::unique_ptr<H5::Group> dataGroup;
		std::string GetHeaderContent(const std::string& headerName);
		void GetDataWithName(const std::string& name, std::vector<double>& out);
		bool HasDataWithName(const std::string&);

	};

	class DllExport HDF5SegmentedDataFile {
		static bool SegmentRawFile(const HDF5RawDataFile&, const std::string&, bool doAnalyses = false);
	}; //HDF5SegmentedDataFile
	
	
	}//namespace IO

}
#endif
