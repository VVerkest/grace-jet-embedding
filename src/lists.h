#ifndef lists__h
#define lists__h

#include <vector>
using namespace std;
// need functions to read in the runId lists
vector<int> read_file_list(const char* list, int col=0, bool sort=false);

#endif
