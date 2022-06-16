#include "lists.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <iostream>

using namespace std;

vector<int> read_file_list(const char* file_name, int col, bool do_sort){
    ifstream f_in;
    f_in.open(file_name);
    vector<int> vec {};
    if (!f_in.is_open()) {
        cout << " fatal: couldn't open input file " << file_name;
        return vec;
    }
    string line;
    while (getline(f_in, line)){
        if (line.rfind("#",0) == 0) continue; // it begins with # and is a comment

        istringstream is{ line };
        int i_val;
        for (int i{0};i<=col;++i) {
            is >> i_val;
            /* cout << " i: " << i << " val " << i_val << endl; */

        }
        vec.push_back(i_val);
    }
    f_in.close();

    // sort the list
    if (do_sort) sort(vec.begin(), vec.end());
    return vec;
};
