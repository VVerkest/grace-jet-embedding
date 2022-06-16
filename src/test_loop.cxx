
#include "events.h"
#include <sstream>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

using namespace std;
void test_loop(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"test_loop\"" << endl;
    istringstream options ( _options );
    int n_options = 0;
    string arg;
    // options >> arg;
    while (options >> arg) {
        cout    << " Option " << n_options << ":  " << arg << endl;
        dat.log << " Option " << n_options << ":  " << arg << endl;
        ++n_options;
    }

    // Histogram declarations here:
    // TH1D hg {"hg", "a;b;c", 10, 0., 1.};


    // Run loop here:
    while (dat.next()) {
        // cout << " Doing nothing in event " << dat.jentry << endl;
        // dat.log  << " Doing nothing in event " << dat.jentry << endl;
        // hg.Fill( x );
    }

    // Write histograms here
    // hg.Write();


    // Wrap-up work here:
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
