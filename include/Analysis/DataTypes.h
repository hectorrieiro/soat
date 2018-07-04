#ifndef ANALYSISDATATYPES_H
#define ANALYSISDATATYPES_H

#if (defined (_MSC_VER))

#ifdef EXPORTING
#define DllExport __declspec(dllexport)
#else
#define DllExport __declspec(dllimport)
#endif

#else
#define DllExport
#endif


#include <string>
#include <vector>
#include <map>
#include <json/value.h>

namespace Analysis {
	class SampleToDva {
	public:
		virtual void convert(double& x, double& y) {};
		virtual ~SampleToDva() = default;
	};
	typedef struct {
		double timestamp;
		double inTrialTS;
		std::string message;
	} H5EventType;
	typedef std::vector<H5EventType> TaskEventsType;
	typedef std::vector<TaskEventsType> ExperimentEventsType;
	typedef std::map<const std::string, std::vector<double> > TaskDataType;
	typedef std::vector<TaskDataType> ExperimentResultsType;

	class DllExport JSONifiable {
	public:
		virtual Json::Value ToJSONObject() = 0;
	};
}

#endif