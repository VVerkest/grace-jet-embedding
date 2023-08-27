#ifndef friend_fnc__h
#define friend_fnc__h

#include "events.h"
#include "ioClass.h"
#include "TreeObj.h"

std::pair<double,double> calc_transEA(events& dat, double phi, vector<int>& bad_towers); 
// calc rho_Nch, rhopTEt

bool cut30_event(events& dat); // return true if there is a track or tower 30 GeV or above

struct Iter_GoodCorrTowers { // will iterate through all good powers 
                             // above 200 MeV
    TClonesArray* tca;
    mupicoTower*  ptr;
    ioIntSet&     bad_towers;
    int           index;
    int           n_towers;

    Iter_GoodCorrTowers  (TClonesArray* _tca, ioIntSet& _bad_towers);
    Iter_GoodCorrTowers begin();
    Iter_GoodCorrTowers end();
    void operator++();
    void skip_bad_towers();
    mupicoTower& operator*();
    bool operator!=(const Iter_GoodCorrTowers& rhs);
};
struct Iter_GoodTowerHits { // Like GoodCorrTowers but
                            // doesn't skip Et<200 (to be used with Et, to EtCorr()
    TClonesArray* tca;
    mupicoTower*  ptr;
    ioIntSet&     bad_towers;
    int           index;
    int           n_towers;

    Iter_GoodTowerHits  (TClonesArray* _tca, ioIntSet& _bad_towers);
    Iter_GoodTowerHits begin();
    Iter_GoodTowerHits end();
    void operator++();
    void skip_bad_towers();
    mupicoTower& operator*();
    bool operator!=(const Iter_GoodTowerHits& rhs);
};


struct get_good_towers {
    get_good_towers(events& _dat, 
            const char* bad_tower_list=
            "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/bad_tower_iter.list");
    events& dat;
    ioIntSet bad_towers;
    mupicoTower* operator()();
    bool next();
    void reset();
    int index;
    int n_towers;
    mupicoTower* tower;
    double Et;
};

struct Iter_GoodTracks { // will iterate through all good powers 
                         // above 200 MeV
    TClonesArray* tca;
    mupicoTrack*  ptr;
    int           index;
    int           n_tracks;

    Iter_GoodTracks  (TClonesArray* _tca);
    Iter_GoodTracks begin();
    Iter_GoodTracks end();
    void operator++();
    void skip_bad_tracks();
    mupicoTrack& operator*();
    bool operator!=(const Iter_GoodTracks& rhs);
};

double get_track_eff(double trPt, TH1D *hEffic);

int get_zdcX_bin(double zdcx, const double* zdcx_bins);

vector<fastjet::PseudoJet> gather_charged_UE(events& dat, double phi, const ioIntSet& bad_towers);

int get_track_eta_bin( double eta );

#endif
