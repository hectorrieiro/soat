
#include "Analysis/DataFile.h"
#include <boost\filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

int main() {
	using namespace boost::filesystem;
	path p("c:\\down\\seg");
	std::string out("c:\\down\\seg\\");
	typedef std::vector<path> vec;
	vec v;
	copy(directory_iterator(p), directory_iterator(), back_inserter(v));
	
	  for (vec::const_iterator it(v.begin()); it != v.end(); it++) {
		if (boost::algorithm::ends_with(it->string(), "dat")) {
		
				Analysis::IO::RawDataFile* f = Analysis::IO::RawDataFile::Read(it->string());
				f->Segment(out + it->filename().generic_string() + ".h5");
			
			
		}

	}
	
	
}