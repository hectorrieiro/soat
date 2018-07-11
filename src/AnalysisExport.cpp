#define EXPORTING
#include "Analysis/AnalysisExport.h"
#include "Analysis/DataTypes.h"
#include "Analysis/DataFile.h"
#include "Analysis/HDF5.h"
using namespace Analysis;
using namespace Analysis::IO;

typedef Analysis::IO::RawDataFile RawType;
typedef Analysis::IO::SegmentedDataFile SegmentedType;
typedef std::unique_ptr<SegmentedType> SegmentedPtr;
typedef std::unique_ptr<RawType> RawPtr;

bool ImportAndProcessDatafile(const char* _inFile, const char* _outFile) {
	
	

	RawPtr f(RawType::Read(_inFile));
	if (f == nullptr) return false;
	SegmentedPtr fs(f->Segment(_outFile));
	if (fs == nullptr) return false;
	bool res = fs->DoAnalyses();

	return res;
	

}

bool AnalyzeSegmentedFile(const char* file) {
	SegmentedPtr f(SegmentedType::Read(file));
	if (f == nullptr) return false;
	bool res = f->DoAnalyses();
	return res;
}
