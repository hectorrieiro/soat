#ifndef ANALYSISEXPORT_H
#define ANALYSISEXPORT_H

#if (defined (_MSC_VER))

#ifdef EXPORTING
#define DllExport __declspec(dllexport)
#else
#define DllExport __declspec(dllimport)
#endif

#else
#define DllExport
#endif





extern "C" {

	DllExport bool ImportAndProcessDatafile(const char* inFile, const char* outFile);
	DllExport bool AnalyzeSegmentedFile(const char* inFile);

}

	
	


#endif
