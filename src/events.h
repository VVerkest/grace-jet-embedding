
#ifndef events_h
#define events_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"
#include "TClonesArray.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include "TreeObj.h"
#include "MemTimeProgression.h"


// template to iterate through TClonesArray members
// (will be used with members that return values with the add_class_lib.py loop)
// will allow iteration such as ` for (auto track : dat.iter_track() ) { /* do stuff */ };
template <class T> struct iterTCA {
    TClonesArray* tca;
    T* ptr;
    iterTCA  (TClonesArray* _tca) : tca{_tca} {};

    int index{0};
    iterTCA begin() {
        iterTCA iter {tca};
        iter.ptr = (T*) tca->UncheckedAt(0);
        return iter;
    };
    iterTCA end() {
        iterTCA iter {tca};
        iter.ptr=(T*)tca->UncheckedAt(tca->GetEntriesFast());
        return iter;
    };
    void operator++() {ptr=(T*)tca->UncheckedAt(++index);};
    T& operator*() {return *ptr;};
    bool operator!=(const iterTCA& rhs) { return ptr!=rhs.ptr;};
};


class events {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   long long int nevents;
   ofstream &log;
   MemTimeProgression stats;

//:TAG START:  Array Sizes
/* No need to use fixed array sizes because TClonesArray* are used instead of Arrays */
//:TAG END: Array Sizes

//:TAG START: Leaf Types
  //
//  note:  TOjects have been substituted in from header files:
// ${HOME}/AN_common/include/TreeObj.h 
// The following TObjects are present in the TTree and should be used.
//
// TObjects used in TTree + members found in above headers:
//  - JetwArea : pt eta phi area index_track index_tower
//  - mupicoEventHeader : runId eventId ZDCx vz BBC_Ein BBC_Eout BBC_Win BBC_Wout vzVpd ranking ZdcSumAdcEast ZdcSumAdcWest refMult vx vy
//  - embTrack : geantId id pt eta phi
//  - embNeutPart : geantId pt eta phi
//  - mupicoTrack : pt eta phi dcaXY dcaXYZ TOF_match BEMC_match towerID towerEt nHitsFit nHitsPoss nHitsDedx pass_cuts
//  - mupicoTower : Et eta phi Et_hadroncorr towerID
//
// TClonesArrays branches found in TTree (with their associated members):
// Note: a 'pseudo-iterator is provided for each TClonesArray. See leafs below.//  - mc_Cjet : pt eta phi area index_track index_tower
// Note: a 'pseudo-iterator is provided for each TClonesArray. See leafs below.//  - mc_Fjet : pt eta phi area index_track index_tower
// Note: a 'pseudo-iterator is provided for each TClonesArray. See leafs below.//  - Cjet : pt eta phi area index_track index_tower
// Note: a 'pseudo-iterator is provided for each TClonesArray. See leafs below.//  - Fjet : pt eta phi area index_track index_tower
// Note: a 'pseudo-iterator is provided for each TClonesArray. See leafs below.//  - mcTr : geantId id pt eta phi
// Note: a 'pseudo-iterator is provided for each TClonesArray. See leafs below.//  - mcNeut : geantId pt eta phi
// Note: a 'pseudo-iterator is provided for each TClonesArray. See leafs below.//  - track : pt eta phi dcaXY dcaXYZ TOF_match BEMC_match towerID towerEt nHitsFit nHitsPoss nHitsDedx pass_cuts
// Note: a 'pseudo-iterator is provided for each TClonesArray. See leafs below.//  - tower : Et eta phi Et_hadroncorr towerID
//
// Branches of TObjects found in TTree (with their associated members):
//  - mu_event : runId eventId ZDCx vx vy vz BBC_Ein BBC_Eout BBC_Win BBC_Wout vzVpd ranking ZdcSumAdcEast ZdcSumAdcWest
//
// Best matches (according to members) of Tree branches to TObject types"
// Note: if any of these "first place guesses" are wrong, fix them by hand!
// "Branch name" -> (class-name):(#members in TTree not in header, #members in header not in TTree) ...
// (TCloneArray Branches)
//  - mc_Cjet -> JetwArea:(0, 0) embNeutPart:(3, 1) embTrack:(3, 2) mupicoTower:(4, 3) mupicoTrack:(3, 10) mupicoEventHeader:(6, 15)
//  - mc_Fjet -> JetwArea:(0, 0) embNeutPart:(3, 1) embTrack:(3, 2) mupicoTower:(4, 3) mupicoTrack:(3, 10) mupicoEventHeader:(6, 15)
//  - Cjet -> JetwArea:(0, 0) embNeutPart:(3, 1) embTrack:(3, 2) mupicoTower:(4, 3) mupicoTrack:(3, 10) mupicoEventHeader:(6, 15)
//  - Fjet -> JetwArea:(0, 0) embNeutPart:(3, 1) embTrack:(3, 2) mupicoTower:(4, 3) mupicoTrack:(3, 10) mupicoEventHeader:(6, 15)
//  - mcTr -> embTrack:(0, 0) embNeutPart:(1, 0) JetwArea:(2, 3) mupicoTower:(3, 3) mupicoTrack:(2, 10) mupicoEventHeader:(5, 15)
//  - mcNeut -> embNeutPart:(0, 0) embTrack:(0, 1) JetwArea:(1, 3) mupicoTower:(2, 3) mupicoTrack:(1, 10) mupicoEventHeader:(4, 15)
//  - track -> mupicoTrack:(0, 0) embNeutPart:(10, 1) embTrack:(10, 2) mupicoTower:(10, 2) JetwArea:(10, 3) mupicoEventHeader:(13, 15)
//  - tower -> mupicoTower:(0, 0) embNeutPart:(3, 2) embTrack:(3, 3) JetwArea:(3, 4) mupicoTrack:(2, 10) mupicoEventHeader:(5, 15)
// (TObject Branches)
//  - mu_event -> mupicoEventHeader:(0, 1) embNeutPart:(14, 4) embTrack:(14, 5) mupicoTower:(14, 5) JetwArea:(14, 6) mupicoTrack:(14, 13)
    UInt_t              mc_Cjet_njets;

    //--accessor and iterator for TClonesArray* (tca) mc_Cjet
    TClonesArray        *tca_mc_Cjet {new TClonesArray("JetwArea")};
    JetwArea            *get_mc_Cjet(int=-1);
    iterTCA<JetwArea>    iter_mc_Cjet { tca_mc_Cjet };
    int                  size_mc_Cjet() { return tca_mc_Cjet->GetEntriesFast(); };

    Float_t             mc_Cjet_rho;
    Float_t             mc_Cjet_rho_sigma;
    UInt_t              mc_Fjet_njets;

    //--accessor and iterator for TClonesArray* (tca) mc_Fjet
    TClonesArray        *tca_mc_Fjet {new TClonesArray("JetwArea")};
    JetwArea            *get_mc_Fjet(int=-1);
    iterTCA<JetwArea>    iter_mc_Fjet { tca_mc_Fjet };
    int                  size_mc_Fjet() { return tca_mc_Fjet->GetEntriesFast(); };

    Float_t             mc_Fjet_rho;
    Float_t             mc_Fjet_rho_sigma;
    UInt_t              Cjet_njets;

    //--accessor and iterator for TClonesArray* (tca) Cjet
    TClonesArray        *tca_Cjet {new TClonesArray("JetwArea")};
    JetwArea            *get_Cjet(int=-1);
    iterTCA<JetwArea>    iter_Cjet { tca_Cjet };
    int                  size_Cjet() { return tca_Cjet->GetEntriesFast(); };

    Float_t             Cjet_rho;
    Float_t             Cjet_rho_sigma;
    UInt_t              Fjet_njets;

    //--accessor and iterator for TClonesArray* (tca) Fjet
    TClonesArray        *tca_Fjet {new TClonesArray("JetwArea")};
    JetwArea            *get_Fjet(int=-1);
    iterTCA<JetwArea>    iter_Fjet { tca_Fjet };
    int                  size_Fjet() { return tca_Fjet->GetEntriesFast(); };

    Float_t             Fjet_rho;
    Float_t             Fjet_rho_sigma;
    mupicoEventHeader   *mu_event;
    Short_t             pthat_bin;

    //--accessor and iterator for TClonesArray* (tca) mcTr
    TClonesArray        *tca_mcTr {new TClonesArray("embTrack")};
    embTrack            *get_mcTr(int=-1);
    iterTCA<embTrack>    iter_mcTr { tca_mcTr };
    int                  size_mcTr() { return tca_mcTr->GetEntriesFast(); };


    //--accessor and iterator for TClonesArray* (tca) mcNeut
    TClonesArray        *tca_mcNeut {new TClonesArray("embNeutPart")};
    embNeutPart         *get_mcNeut(int=-1);
    iterTCA<embNeutPart>  iter_mcNeut { tca_mcNeut };
    int                  size_mcNeut() { return tca_mcNeut->GetEntriesFast(); };


    //--accessor and iterator for TClonesArray* (tca) track
    TClonesArray        *tca_track {new TClonesArray("mupicoTrack")};
    mupicoTrack         *get_track(int=-1);
    iterTCA<mupicoTrack>  iter_track { tca_track };
    int                  size_track() { return tca_track->GetEntriesFast(); };


    //--accessor and iterator for TClonesArray* (tca) tower
    TClonesArray        *tca_tower {new TClonesArray("mupicoTower")};
    mupicoTower         *get_tower(int=-1);
    iterTCA<mupicoTower>  iter_tower { tca_tower };
    int                  size_tower() { return tca_tower->GetEntriesFast(); };

    UShort_t            EastBBC[24];
    UShort_t            WestBBC[24];
    UShort_t            ZdcSmdEastHorizontal[8];
    UShort_t            ZdcSmdEastVertical[8];
    UShort_t            ZdcSmdWestHorizontal[8];
    UShort_t            ZdcSmdWestVertical[8];

//:TAG END: Leaf Types

//:TAG START: Declare Branches
       // List of branches
    TBranch        *b_mc_Cjet_njets;   //!
    TBranch        *b_mc_Cjet_rho;   //!
    TBranch        *b_mc_Cjet_rho_sigma;   //!
    TBranch        *b_mc_Fjet_njets;   //!
    TBranch        *b_mc_Fjet_rho;   //!
    TBranch        *b_mc_Fjet_rho_sigma;   //!
    TBranch        *b_Cjet_njets;   //!
    TBranch        *b_Cjet_rho;   //!
    TBranch        *b_Cjet_rho_sigma;   //!
    TBranch        *b_Fjet_njets;   //!
    TBranch        *b_Fjet_rho;   //!
    TBranch        *b_Fjet_rho_sigma;   //!
    TBranch        *b_pthat_bin;   //!
    TBranch        *b_s;   //!
    TBranch        *b_s_0;   //!
    TBranch        *b_s_1;   //!
    TBranch        *b_s_2;   //!
    TBranch        *b_s_3;   //!
    TBranch        *b_s_4;   //!

//:TAG END: Declare Branches

   map<int,bool*>      trigger_map;
   bool                has_trigger(int);
   bool                has_trigger_all(vector<int>);
   bool                has_trigger_any(vector<int>);

   events(ofstream& log, int n_events, TString inlist);
   virtual ~events();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // values for running the loop
   bool next();
   Long64_t nentries;
   Long64_t jentry;

   // Make friend functions to fill in the runs
   // i.e. friend void runLoop(events&, string);
   // TAG: start-friend-functions
   friend void sys_err(events&, string);
   friend void eta_match(events&, string);
   friend void thesis_emb(events&, string);
   friend void five_rooResF(events&, string);
   friend void jet_draw(events&, string);
   friend void bins70(events&, string);
   friend void sane_bins(events&, string);
   friend void crazy_limit(events&, string);
   friend void large_ones(events&, string);
   friend void list_ids(events&, string);
   friend void crazy(events&, string);
   friend void rat_rooResF(events&, string);
   friend void cut_rooResF(events&, string);
   friend void big_rooResF(events&, string);
   friend void rooResF_check(events&, string);
   friend void rooResF(events&, string);
   friend void test_loop(events&, string);

};
#endif
