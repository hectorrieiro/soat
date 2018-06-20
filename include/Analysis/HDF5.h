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

	class DllExport HDF5DataFile {
	public:
		static bool Text2HDF5(const std::string& origin, const std::string& dest);
		static bool HDF52Text(const std::string& origin, const std::string& dest);
		static bool IsValidHDF5(const std::string& file);

		HDF5DataFile(const std::string& filename);
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
		
		virtual ~HDF5DataFile();

	private:

		std::unique_ptr<H5::H5File> file;
		std::unique_ptr<H5::Group> dataGroup;

		std::string GetHeaderContent(const std::string& headerName);
		void GetDataWithName(const std::string& name, std::vector<double>& out);
		bool HasDataWithName(const std::string&);


		typedef struct SimpleTableType_{
			char* title;
			char* content;
			
		} SimpleTableType;

		typedef struct {
			double timestamp;
			char* text;
		} EventType;
	};

}
#endif
