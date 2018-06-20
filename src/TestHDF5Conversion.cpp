
#include "Analysis/HDF5.h"
#include <boost\filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

int main() {
	using namespace boost::filesystem;
	path p("c:\\down\\sparcc");
	std::string out("c:\\down\\sparcc\\h5\\");
	typedef std::vector<path> vec;
	vec v;
	copy(directory_iterator(p), directory_iterator(), back_inserter(v));
	
	  for (vec::const_iterator it(v.begin()); it != v.end(); it++) {
		if (boost::algorithm::ends_with(it->string(), "dat")) {
			try {
				Analysis::HDF5DataFile::Text2HDF5(it->string(), out + it->filename().generic_string());
			}
			catch (...) {

			}
			
		}

	}
	
	
}