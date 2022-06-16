
#define events_cxx
#include "events.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


events::events(ofstream& _log, int _n_events, TString inp_file) :
    nevents{_n_events},
    log {_log},
    stats {10000},
    jentry {-1}
{
    TChain* tree  = new TChain("events");
    if (inp_file.EndsWith(".root")) {
            cout << " Adding input file: " << inp_file << endl;
        tree->Add(inp_file.Data());
    } else if (inp_file.EndsWith(".list")) {
        string line;
        ifstream list;
        list.open(inp_file.Data());
        while (getline(list, line)){
            cout << " Adding input file: " << line << endl;
            tree->Add(line.c_str());
        }
        list.close();
    }
    /* cout << " has " << tree->GetEntries() << " in tree" << endl; */
    /* cout << " has " << tree->GetEntriesFast() << " in tree" << endl; */
    Init(tree);

    //Set the trigger map

    nentries = tree->GetEntries();
    if (nevents == -1) nevents = nentries;
}
events::~events()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t events::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t events::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void events::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //:TAG START: Set Branches
   fChain->SetMakeClass(0); // note if is SetBranchAddress(1) then TObject's
                         // couldn't be read from the tree
    mu_event        = nullptr;
    // List of branches
    fChain->SetBranchAddress("mc_Cjet_njets", &mc_Cjet_njets, &b_mc_Cjet_njets);
    fChain->SetBranchAddress("mc_Cjet", &tca_mc_Cjet);
    fChain->SetBranchAddress("mc_Cjet_rho", &mc_Cjet_rho, &b_mc_Cjet_rho);
    fChain->SetBranchAddress("mc_Cjet_rho_sigma", &mc_Cjet_rho_sigma, &b_mc_Cjet_rho_sigma);
    fChain->SetBranchAddress("mc_Fjet_njets", &mc_Fjet_njets, &b_mc_Fjet_njets);
    fChain->SetBranchAddress("mc_Fjet", &tca_mc_Fjet);
    fChain->SetBranchAddress("mc_Fjet_rho", &mc_Fjet_rho, &b_mc_Fjet_rho);
    fChain->SetBranchAddress("mc_Fjet_rho_sigma", &mc_Fjet_rho_sigma, &b_mc_Fjet_rho_sigma);
    fChain->SetBranchAddress("Cjet_njets", &Cjet_njets, &b_Cjet_njets);
    fChain->SetBranchAddress("Cjet", &tca_Cjet);
    fChain->SetBranchAddress("Cjet_rho", &Cjet_rho, &b_Cjet_rho);
    fChain->SetBranchAddress("Cjet_rho_sigma", &Cjet_rho_sigma, &b_Cjet_rho_sigma);
    fChain->SetBranchAddress("Fjet_njets", &Fjet_njets, &b_Fjet_njets);
    fChain->SetBranchAddress("Fjet", &tca_Fjet);
    fChain->SetBranchAddress("Fjet_rho", &Fjet_rho, &b_Fjet_rho);
    fChain->SetBranchAddress("Fjet_rho_sigma", &Fjet_rho_sigma, &b_Fjet_rho_sigma);
    fChain->SetBranchAddress("mu_event", &mu_event);
    fChain->SetBranchAddress("pthat_bin", &pthat_bin, &b_pthat_bin);
    fChain->SetBranchAddress("mcTr", &tca_mcTr);
    fChain->SetBranchAddress("mcNeut", &tca_mcNeut);
    fChain->SetBranchAddress("track", &tca_track);
    fChain->SetBranchAddress("tower", &tca_tower);
    fChain->SetBranchAddress("EastBBC", EastBBC, &b_s);
    fChain->SetBranchAddress("WestBBC", WestBBC, &b_s_0);
    fChain->SetBranchAddress("ZdcSmdEastHorizontal", ZdcSmdEastHorizontal, &b_s_1);
    fChain->SetBranchAddress("ZdcSmdEastVertical", ZdcSmdEastVertical, &b_s_2);
    fChain->SetBranchAddress("ZdcSmdWestHorizontal", ZdcSmdWestHorizontal, &b_s_3);
    fChain->SetBranchAddress("ZdcSmdWestVertical", ZdcSmdWestVertical, &b_s_4);

   //:TAG END: Set Branches
   Notify();
}

Bool_t events::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void events::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t events::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

bool events::next() {
    jentry++;
    if (jentry >= nevents) {
        stats.set_get_stats();
        log  << " Final stats: " << stats.stats << endl;
        return false;
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) {
        cout << " Breaking out of loop at jentry on failure to read ientry: " << jentry << endl;
        stats.set_get_stats();
        log  << " Final stats: " << stats.stats << endl;
        return false;
    }
    fChain->GetEntry(jentry);
    if (stats.call()) log << stats.stats << endl;
    return true;
}

//:TAG START: Coda Functions
// Two  psuedo-iterator functions for TClonesArray(JetwArea) *tca_mc_Cjet
JetwArea* events::get_mc_Cjet(int i) {
    return (JetwArea*) tca_mc_Cjet->UncheckedAt(i);
};

// Two  psuedo-iterator functions for TClonesArray(JetwArea) *tca_mc_Fjet
JetwArea* events::get_mc_Fjet(int i) {
    return (JetwArea*) tca_mc_Fjet->UncheckedAt(i);
};

// Two  psuedo-iterator functions for TClonesArray(JetwArea) *tca_Cjet
JetwArea* events::get_Cjet(int i) {
    return (JetwArea*) tca_Cjet->UncheckedAt(i);
};

// Two  psuedo-iterator functions for TClonesArray(JetwArea) *tca_Fjet
JetwArea* events::get_Fjet(int i) {
    return (JetwArea*) tca_Fjet->UncheckedAt(i);
};

// Two  psuedo-iterator functions for TClonesArray(embTrack) *tca_mcTr
embTrack* events::get_mcTr(int i) {
    return (embTrack*) tca_mcTr->UncheckedAt(i);
};

// Two  psuedo-iterator functions for TClonesArray(embNeutPart) *tca_mcNeut
embNeutPart* events::get_mcNeut(int i) {
    return (embNeutPart*) tca_mcNeut->UncheckedAt(i);
};

// Two  psuedo-iterator functions for TClonesArray(mupicoTrack) *tca_track
mupicoTrack* events::get_track(int i) {
    return (mupicoTrack*) tca_track->UncheckedAt(i);
};

// Two  psuedo-iterator functions for TClonesArray(mupicoTower) *tca_tower
mupicoTower* events::get_tower(int i) {
    return (mupicoTower*) tca_tower->UncheckedAt(i);
};

//:TAG END: Coda Functions

// warning, below function will seg-fault if trigger is not in map
bool  events::has_trigger(int i_trig) { return *(trigger_map[i_trig]); };
bool  events::has_trigger_all(vector<int> triggers) { 
    for (auto T : triggers)  if (! *(trigger_map[T])) return false;
    return true;
};
bool  events::has_trigger_any(vector<int> triggers) { 
    for (auto T : triggers)  if (*(trigger_map[T])) return true;
    return false;
};
