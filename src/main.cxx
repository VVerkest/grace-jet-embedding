
#include "events.h"
#include <sstream>
#include <fstream>
#include "TFile.h"
// #include <algorithm>
using namespace std;

int main(int nargs, char** argv) {
    /*
     * arguments:
     *   1: number of events
     *   2: input root file list name
     *   3: which program to run
     *   4: output base-name
     *   5: optional input
     */

    int n_events          {  (nargs>1) ? atoi(argv[1]) : -1 };
    string inp_list       {  (nargs>2) ? argv[2] : "in-lists/list_test.list" };
    string which_loop     {  (nargs>3) ? argv[3] : "jetEmbedding_loop" };
    string o_name_tag     {  (nargs>4) ? argv[4] : "jetEmbedding_loop" };

    ostringstream collect;
    for (int i{5};i<nargs;++i) {
        string arg {argv[i]};
        collect << arg << " ";
    }

    ofstream log;
    log.open((o_name_tag + ".log").c_str());
    log << "Starting output."  << endl
        << "Command line input:" << endl;
    for (int i{0}; i<nargs; ++i) log << "arg("<<i<<")  " << argv[i] << endl;
    log << endl << endl;

    events my_events{log, n_events, inp_list};

    
    TFile fout  { (o_name_tag+".root").c_str(), "recreate" };

    // run the loop
    cout << " Looking for function: " << which_loop << endl;
    if (which_loop == "empty-loop") {
        // TAG: empty-loop
   } else if (which_loop == "sys_err") {
        sys_err(my_events, collect.str());
   } else if (which_loop == "eta_match") {
        eta_match(my_events, collect.str());
   } else if (which_loop == "jetEmbedding_loop") {
        jetEmbedding_loop(my_events, collect.str());
   } else if (which_loop == "five_rooResF") {
        five_rooResF(my_events, collect.str());
   } else if (which_loop == "jet_draw") {
        jet_draw(my_events, collect.str());
   } else if (which_loop == "bins70") {
        bins70(my_events, collect.str());
   } else if (which_loop == "sane_bins") {
        sane_bins(my_events, collect.str());
   } else if (which_loop == "crazy_limit") {
        crazy_limit(my_events, collect.str());
   } else if (which_loop == "large_ones") {
        large_ones(my_events, collect.str());
   } else if (which_loop == "list_ids") {
        list_ids(my_events, collect.str());
   } else if (which_loop == "crazy") {
        crazy(my_events, collect.str());
   } else if (which_loop == "rat_rooResF") {
        rat_rooResF(my_events, collect.str());
   } else if (which_loop == "cut_rooResF") {
        cut_rooResF(my_events, collect.str());
   } else if (which_loop == "big_rooResF") {
        big_rooResF(my_events, collect.str());
   } else if (which_loop == "rooResF_check") {
        rooResF_check(my_events, collect.str());
   } else if (which_loop == "rooResF") {
        rooResF(my_events, collect.str());
   } else if (which_loop == "test_loop") {
        test_loop(my_events, collect.str());
    } else {
        cout << " -  Fatal error: couldn't find loop \"" << which_loop << "\"" << endl
            << "     ->  Terminating program." << endl;
    }

    fout.Close();
    cout << " Finished TFile:  " << endl;
    cout << "scp grace:/home/djs232/w2021/jet-embedding/" 
         << (o_name_tag+".root").c_str() << " ." << endl;
    log.close();
};
