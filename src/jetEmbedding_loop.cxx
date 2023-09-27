
#include "events.h"
#include <sstream>
#include <fstream>

#include <iomanip>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "friend_fnc.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

const double bin_leadPt[51] = {4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54.};
const double bin_iBBCEsum[11] = { 0.,  2802.38, 5460.32, 8420.66, 11722, 15410.1, 19584.9, 24371, 30105.3, 37666.3, 64000};
const double bin_ZDCx[6] = {5000., 8000., 11000., 14000., 17000., 20000.};
const double ZDCxWeight[5] = {1, 1.7291, 2.0717, 2.64997, 3.8812};

double calc_absDphi(double a, double b) {
    double val = abs(a - b);
    while (val > M_PI) val = abs(val - 2*M_PI);
    return val;
};

bool trigger_matches_jet(double trig_eta, double trig_phi, double lead_eta, double lead_phi) {
    const auto dphi = calc_absDphi(lead_phi,trig_phi);
    const auto deta = lead_eta-trig_eta;
    return  (dphi>(M_PI-0.4)) || ((dphi*dphi+deta*deta)<0.16);
};

int get_zdcX_bin(double zdcx) {
    if (zdcx>20000. || zdcx<5000.) { cout<<"zdcx must be between 5k and 20k"<<endl; }
    for (int i=0; i<5; ++i) {
        if ( zdcx>=bin_ZDCx[i] && zdcx<=bin_ZDCx[i+1] ) { return i; }
    }
    return 9999;
};

bool is_phi_p4match(double phi_A, double phi_B) {
    double delta = TMath::Abs(phi_A - phi_B);
    while (delta > M_PI) delta -= M_PI;
    if (delta < 0.4) return true;
    if (delta > (M_PI-0.4)) return true;
};

bool matched_jets( double phi_A, double phi_B, double eta_A, double eta_B ) {
    double dPhi = TMath::Abs(phi_A - phi_B);
    while (dPhi > 2.*M_PI) dPhi -= 2.*M_PI;
    double dEta = TMath::Abs(eta_A - eta_B);
    double dR = sqrt( dPhi*dPhi + dEta*dEta );
    if (dR<0.4) { return true; }
};

double pyth6_10x_Xsec(int i) { // here "i" indicates the pt-hat bin, which comes from dat.pthat_bin
    
    const array<double, 9> PYTHIA6_XSEC {
            0.1604997783173907, 0.0279900193730690, 0.006924398431,
            0.0028726057079642, 0.0005197051748372, 0.0000140447879818,
            0.0000006505378525, 0.0000000345848665, 0.0000000016149182 };
    const array<double, 9> PYTHIA6_EMB_NEVENTS { 3920155, 2101168, 1187058, 1695141, 4967075, 1797387, 260676, 261926, 262366};
    
    return PYTHIA6_XSEC[i] / PYTHIA6_EMB_NEVENTS[i];
};


using namespace std;
using namespace fastjet;

void jetEmbedding_loop(events& dat, string _options) {
    // Start of
    cout << " Running fn \"jetEmbedding_loop\"" << endl;
    istringstream options ( _options );
    int n_options = 0;
    string arg;
    // options >> arg;
    while (options >> arg) {
        cout    << " Option " << n_options << ":  " << arg << endl;
        dat.log << " Option " << n_options << ":  " << arg << endl;
        ++n_options;
    }


    // list of good events and bad towers
    ioIntMap  runid     { "in-data/good_run.list" , 1, 0 ,false };
    ioIntSet  badtow    { "in-data/bad_tower.list", 0, false };

    int triggerid { 500206 };
    const double vzCut = 10.0;   // |Vz|<=10 cm
    const double VertexZDiffCut = 6.0;
    
    JetDefinition jet_def(antikt_algorithm, 0.4);     //  JET DEFINITION
    Selector jetEtaSelector = SelectorAbsEtaMax( 0.6 );
    
    // Histogram declarations here:
    TH3D *hMatchedJets = new TH3D("hMatchedJets",";part-level leading jet p_{T} [GeV/c];det-level leading jet p_{T} [GeV/c];EA_{BBC}",50,bin_leadPt,50,bin_leadPt,10,bin_iBBCEsum);
    TH2D *hMissed = new TH2D("hMissed",";missed part-level leading jet p_{T} [GeV/c];EA_{BBC}",50,bin_leadPt,10,bin_iBBCEsum);
    TH2D *hFake = new TH2D("hFake",";fake det-level leading jet p_{T} [GeV/c];EA_{BBC}",50,bin_leadPt,10,bin_iBBCEsum);
    TH3D *hEventInfo = new TH3D("hEventInfo",";trigger tower E_{T} [GeV];det-level leading jet p_{T} [GeV/c];part-level leading jet p_{T} [GeV/c]",50,0.,50.,50,0.,50.,50,0.,50.);
    // maybe I will turn these into arrays or sparses for the jet pt-hat bins...
    
    vector<PseudoJet> jet_inputs, rawJets;

    while (dat.next()) { // / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
        
        // dat.log  << " Doing nothing in event " << dat.jentry << endl;
        if (!runid.has(dat.mu_event->runId)) continue; // keep only good runs
        if (fabs(dat.mu_event->vz)>vzCut) { continue; } // cut on vz
        if (fabs(dat.mu_event->vz - dat.mu_event->vzVpd)>VertexZDiffCut) { continue; } // cut on vz diff

        if (dat.mu_event->ZDCx<5000. || dat.mu_event->ZDCx>20000.)  continue; // cut on ZDCx

//        for (auto track : Iter_GoodTracks(dat.tca_track)) {        }
//        for (auto track : dat.iter_mcTr){        }
//        for (auto part : dat.iter_mcNeut){            }
        double trigEt = 0., trigEta, trigPhi;
        for (auto tower : Iter_GoodTowerHits(dat.tca_tower,badtow)) {
            if (tower.Et > trigEt && tower.Et < 30. && tower.Et > 4.) {
                trigEt = tower.Et;
                trigEta = tower.eta;
                trigPhi = tower.phi;
            }
        }
        if ( trigEt<4. ){ continue; } // take the highest energy offline trigger tower >4GeV
        
        jet_inputs.clear(), rawJets.clear();
        
        for (auto& track : Iter_GoodTracks(dat.tca_track)) {  // LOOP OVER TRACKS TO CLUSTER
            if(track.pt>0.2 && track.pt<30.) {
                fastjet::PseudoJet current;
                current.reset_momentum_PtYPhiM( track.pt, track.eta, track.phi, 0. );
                jet_inputs.push_back( current );
            }
        }
        for (auto tower : Iter_GoodCorrTowers(dat.tca_tower,badtow)) {
            if(tower.EtCorr()>0.2 && tower.EtCorr()<30.) {
                fastjet::PseudoJet current;
                current.reset_momentum_PtYPhiM( tower.EtCorr(), tower.eta, tower.phi, 0. );
                jet_inputs.push_back( current );
            }
        }
        
        ClusterSequence jetCluster( jet_inputs, jet_def );
        rawJets = sorted_by_pt( jetEtaSelector( jetCluster.inclusive_jets() ) );
        
        if (rawJets.size()==0 && dat.tca_mc_Fjet->GetEntriesFast()==0) { continue; } // REQUIRE JET
        
        int zdcxbin = get_zdcX_bin(dat.mu_event->ZDCx);

        double xsec_wt = pyth6_10x_Xsec( dat.pthat_bin );
        
        double mc_leadPt = 0., mc_leadEta, mc_leadPhi;
        for (auto jet : dat.iter_mc_Fjet){
            if (jet.pt > mc_leadPt && jet.pt > 4.) { // take highest-pT particle-level jet
                mc_leadPt = jet.pt;
                mc_leadEta = jet.eta;
                mc_leadPhi = jet.phi;
            }
        }
        
        double leadPt = 0., leadEta, leadPhi;
        if (rawJets.size()>0) {
            if (rawJets[0].pt() > 4.) { // take highest-pT jet
                leadPt = rawJets[0].pt();
                leadEta = rawJets[0].eta();
                leadPhi = rawJets[0].phi();
            }
        }
        if ( leadPt<4. || rawJets.size()==0 ) {
            if ( mc_leadPt<4. ) { continue; }
            hMissed->Fill(mc_leadPt, dat.mu_event->BBC_Ein, xsec_wt*ZDCxWeight[zdcxbin]); // if part jet and no det jet: fill missed
            continue;
        }
        /*else {
            cout<<"Veronica should revisit this after a snack"<<endl;
            break;
        }*/
//        cout<<trigPhi<<" \t "<<leadPhi<<endl;
//        if( !is_phi_p4match( trigPhi, leadPhi) ) { continue; }  // require trigger in det-level leading jet (or recoil region)
        if (!trigger_matches_jet(trigEta, trigPhi, leadEta, leadPhi)) {continue;}
        
        if ( mc_leadPt<4. ) {
            hFake->Fill(leadPt, dat.mu_event->BBC_Ein, xsec_wt*ZDCxWeight[zdcxbin]);         // FAKE JET
            continue;
        }
        
        if ( matched_jets( leadPhi, mc_leadPhi, leadEta, mc_leadEta ) ) {
            hMatchedJets->Fill( mc_leadPt, leadPt, dat.mu_event->BBC_Ein, xsec_wt*ZDCxWeight[zdcxbin] );        // FILL RESPONSE
        }
        else {            // FILL MISS AND MATCH
            hMissed->Fill(mc_leadPt, dat.mu_event->BBC_Ein, xsec_wt*ZDCxWeight[zdcxbin]);
            hFake->Fill(leadPt, dat.mu_event->BBC_Ein, xsec_wt*ZDCxWeight[zdcxbin]);
        }

    } // / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
    

    // Write histograms here
    hMatchedJets->Write();
    hMissed->Write();
    hFake->Write();
    hEventInfo->Write();

    // Wrap-up work here:
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
