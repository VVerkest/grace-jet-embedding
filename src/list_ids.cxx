#include "events.h"
#include <sstream>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TreeObj.h"
#include "friend_fnc.h"
#include "ioClass.h"
/* #include "vzZDCx_correlator.h" */

#include "ioXsec_pAu2015.h"
#include "ioJetMatcher.h"
#include "oiJetMaker.h"
/* #include "ioBins.h" */

#include "RooUnfoldResponse.h"


using namespace std;
void list_ids(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"list_ids\"" << endl;
    istringstream options ( _options );
    map<int,map<int,short>> data;
    ioMsgTree msg{false};
    msg.slurp_file("src/list_ids.cxx");
    msg.write();

    while (dat.next()) {
        int runId = dat.mu_event->runId;
        int eventId = dat.mu_event->eventId;
        if (data.count(runId)==0) {
            data[runId] = {};
            data[runId][eventId] = 1;
        } else {
            if (data[runId].count(eventId) == 0)
                data[runId][eventId] = 1;
            else
                data[runId][eventId]++;
        }
    }
    // make the output tree
    TTree tree {"ids","tree of runId, eventId, nEvents"};
    int runId, eventId; 
    short nEvents;
    tree.Branch("runId",&runId);
    tree.Branch("eventId",&eventId);
    tree.Branch("nCopies",&nEvents);

    for (auto& runs : data) {
        runId = runs.first;
        for (auto &events : runs.second) {
            eventId = events.first;
            nEvents = events.second;
            tree.Fill();
        }
    }
    tree.Write();

    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
