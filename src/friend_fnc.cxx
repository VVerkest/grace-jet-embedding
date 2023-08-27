#include "friend_fnc.h"
#include "io_fnc.h"

std::pair<double,double> calc_transEA(events& dat, double phi, vector<int>& bad_towers) {
    const double rhoArea { (IO_piless1-1)*2. };
    double Nch{0};
    double sumpTEt{0};
    for (auto tower : dat.iter_tower){
        if (!io_isAbsTransPhi(tower.phi,phi)) continue;
        if (binary_search(bad_towers.begin(),bad_towers.end(), tower.towerID)) continue;
        sumpTEt += tower.EtCorr();
    }
    for (auto track : dat.iter_track) {
        if (!io_isAbsTransPhi(track.phi,phi)) continue;
        Nch += 1;
        sumpTEt += track.pt;
    }
    return { Nch, sumpTEt / rhoArea };
};

bool cut30_event(events& dat) {
    for (auto tower : dat.iter_tower)  if (tower.Et >= 30.) return true;
    for (auto track : dat.iter_track)  if (track.pt >= 30.) return true;
    for (auto track : dat.iter_mcTr)   if (track.pt >= 30.) return true;
    for (auto track : dat.iter_mcNeut) if (track.pt >= 30.) return true;
    return false;
};

get_good_towers::get_good_towers(events& _dat, const char* bad_tower_list) :
    dat{_dat}, bad_towers{bad_tower_list}, index{0}, n_towers {0} 
{};
bool get_good_towers::next() {
    if (index == 0) n_towers = dat.size_tower(); 
    while (index<n_towers) {
        tower = dat.get_tower(index);
        if (bad_towers(tower->towerID)) { ++index; continue; }
        Et = (tower->Et - tower->Et_hadroncorr);
        if (Et<0.2) { 
            ++index; 
            continue; 
        }
        break; 
    }
    if (index == n_towers) { 
        index=0; 
        return false; 
    }
    ++index;
    return true;
};
void get_good_towers::reset() { index=0; };
mupicoTower* get_good_towers::operator()(){ return tower; };

Iter_GoodCorrTowers::Iter_GoodCorrTowers
(TClonesArray* _tca, ioIntSet& _bad_towers) : 
        tca{_tca}, 
        bad_towers{_bad_towers},
        n_towers{tca->GetEntriesFast()},
        index{0}
{};
void Iter_GoodCorrTowers::skip_bad_towers() {
    while ( 
          index != n_towers 
       && (
              bad_towers.has(ptr->towerID)
           || ptr->EtCorr() < 0.2
          )
    ){
        ptr=(mupicoTower*)tca->UncheckedAt(++index);
    }
};
Iter_GoodCorrTowers Iter_GoodCorrTowers::begin() {
    Iter_GoodCorrTowers iter {tca, bad_towers};
    iter.ptr = (mupicoTower*) tca->UncheckedAt(0);
    iter.skip_bad_towers();
    return iter;
};
Iter_GoodCorrTowers Iter_GoodCorrTowers::end() {
    Iter_GoodCorrTowers iter {tca, bad_towers};
    iter.index = iter.n_towers;
    return iter;
};
void Iter_GoodCorrTowers::operator++() {
    ptr=(mupicoTower*)tca->UncheckedAt(++index);
    skip_bad_towers();
};
mupicoTower& Iter_GoodCorrTowers::operator*() {return *ptr;};
bool Iter_GoodCorrTowers::operator!=(const Iter_GoodCorrTowers& rhs) 
{ return index!=n_towers;};



Iter_GoodTowerHits::Iter_GoodTowerHits
(TClonesArray* _tca, ioIntSet& _bad_towers) : 
        tca{_tca}, 
        bad_towers{_bad_towers},
        n_towers{tca->GetEntriesFast()},
        index{0}
{};
void Iter_GoodTowerHits::skip_bad_towers() {
    while ( 
          index != n_towers 
       && bad_towers.has(ptr->towerID)
    ){
        ptr=(mupicoTower*)tca->UncheckedAt(++index);
    }
};
Iter_GoodTowerHits Iter_GoodTowerHits::begin() {
    Iter_GoodTowerHits iter {tca, bad_towers};
    iter.ptr = (mupicoTower*) tca->UncheckedAt(0);
    iter.skip_bad_towers();
    return iter;
};
Iter_GoodTowerHits Iter_GoodTowerHits::end() {
    Iter_GoodTowerHits iter {tca, bad_towers};
    iter.index = iter.n_towers;
    return iter;
};
void Iter_GoodTowerHits::operator++() {
    ptr=(mupicoTower*)tca->UncheckedAt(++index);
    skip_bad_towers();
};
mupicoTower& Iter_GoodTowerHits::operator*() {return *ptr;};
bool Iter_GoodTowerHits::operator!=(const Iter_GoodTowerHits& rhs) 
{ return index!=n_towers;};

Iter_GoodTracks::Iter_GoodTracks
(TClonesArray* _tca) :
        tca{_tca}, 
        n_tracks{tca->GetEntriesFast()},
        index{0}
{};
void Iter_GoodTracks::skip_bad_tracks() {
    while ( 
          index != n_tracks 
       && (!ptr->pass_cuts || (0.52*ptr->nHitsPoss > ptr->nHitsFit) || ptr->dcaXYZ>=3.)
    ){
        ptr=(mupicoTrack*)tca->UncheckedAt(++index);
    }
};
Iter_GoodTracks Iter_GoodTracks::begin() {
    Iter_GoodTracks iter {tca};
    iter.ptr = (mupicoTrack*) tca->UncheckedAt(0);
    iter.skip_bad_tracks();
    return iter;
};
Iter_GoodTracks Iter_GoodTracks::end() {
    Iter_GoodTracks iter {tca};
    iter.index = iter.n_tracks;
    return iter;
};
void Iter_GoodTracks::operator++() {
    ptr=(mupicoTrack*)tca->UncheckedAt(++index);
    skip_bad_tracks();
};
mupicoTrack& Iter_GoodTracks::operator*() {return *ptr;};
bool Iter_GoodTracks::operator!=(const Iter_GoodTracks& rhs) 
{ return index!=n_tracks;};

double get_track_eff(double trPt, TH1D *hEffic) {
    int ptBin = hEffic->FindBin(trPt);
    return hEffic->GetBinContent(ptBin);
};

int get_zdcX_bin(double zdcx, const double* bin_ZDCx) {
    if (zdcx>20000. || zdcx<5000.) { cout<<"zdcx must be between 5k and 20k"<<endl; }
    for (int i=0; i<5; ++i) {
        if ( zdcx>=bin_ZDCx[i] && zdcx<=bin_ZDCx[i+1] ) { return i; }
    }
    return 9999;
};

vector<fastjet::PseudoJet> gather_charged_UE(events& dat, double phi, const ioIntSet& bad_towers) {
    vector<fastjet::PseudoJet> chargedUE;
    for (auto& track : Iter_GoodTracks(dat.tca_track)) {
//        if (!is_phi_trans(track.phi,phi) || fabs(track.eta)>0.9) continue;
        if (!is_phi_trans(track.phi,phi)) continue;
        if(track.pt<0.2) {continue;}//TEMPORARY!!!
      else {
            fastjet::PseudoJet current;
            current.reset_momentum_PtYPhiM( track.pt, track.eta, track.phi, 0. );
            chargedUE.push_back( current );
        }
    }
    return chargedUE;
};

int get_track_eta_bin( double eta ){
    double eta20lo[20] = { -1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    double eta20hi[20] = { -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. };
    if (fabs(eta)>0.9) {
        cout<<"eta must be a value between -0.9 and 0.9 inclusive"<<endl;
//    if (fabs(eta)>1.) {
//        cout<<"eta must be a value between -1 and 1 inclusive"<<endl;
        return 9999;
    }
    else {
        for (int i=0; i<20; ++i) {
            if (eta>=eta20lo[i] && eta<=eta20hi[i]) { return i; }
        }
    }
    return 9999;
};


