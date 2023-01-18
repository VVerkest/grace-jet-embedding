
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

#include "ioClass.h"


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

    const double bin_leadPt[51] = {4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54.};
    const double bin_iBBCEsum[11] = { 0., 2767.15, 5397.98, 8333.35, 11610.5, 15280.9, 19440.2, 24219.7, 29959, 37534.5, 64000. }; //{0., 3559.12, 6735.12, 10126.1, 13752.1, 17669.1, 21948.1, 26718.1, 32283.1, 39473.1, 64000.};

    // list of good events and bad towers
    ioIntMap  runid     { "in-data/good_run.list" , 1, 0 ,false };
    ioIntSet  badtow    { "in-data/bad_tower.list", 0, false };

    int triggerid { 500206 };
    const double vzCut = 10.0;   // |Vz|<=10 cm
    const double VertexZDiffCut = 6.0;

    // Histogram declarations here:
    TH3D *hMatchedJets = new TH3D("hMatchedJets",";part-level leading jet p_{T} [GeV/c];det-level leading jet p_{T} [GeV/c];EA_{BBC}",50,bin_leadPt,50,bin_leadPt,10,bin_iBBCEsum);
    TH2D *hMissed = new TH2D("hMissed",";missed part-level leading jet p_{T} [GeV/c];EA_{BBC}",50,bin_leadPt,10,bin_iBBCEsum);
    TH2D *hFake = new TH2D("hFake",";fake det-level leading jet p_{T} [GeV/c];EA_{BBC}",50,bin_leadPt,10,bin_iBBCEsum);
    // maybe I will turn these into arrays or sparses for the jet pt-hat bins...
    

    while (dat.next()) { // / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
        
        // dat.log  << " Doing nothing in event " << dat.jentry << endl;
        if (!runid.has(dat.mu_event->runId)) continue; // keep only good runs
        if (fabs(dat.mu_event->vz)>vzCut) { continue; } // cut on vz
        if (fabs(dat.mu_event->vz - dat.mu_event->vzVpd)>VertexZDiffCut) { continue; } // cut on vz diff

//        for (auto track : Iter_GoodTracks(dat.tca_track)) {        }
//        for (auto track : dat.iter_mcTr){        }
//        for (auto part : dat.iter_mcNeut){            }
        double trigEt = 0., trigEta, trigPhi;
        for (auto tower : Iter_GoodCorrTowers(dat.tca_tower,badtow)) {
            if (tower.Et > trigEt && tower.Et < 30. && tower.Et > 4.) {
                trigEt = tower.Et;
                trigEta = tower.eta;
                trigPhi = tower.phi;
            }
        }
        if ( trigEt<4. ){ continue; } // take the highest energy offline trigger tower >4GeV
        if (dat.tca_Fjet->GetEntriesFast()==0 && dat.tca_mc_Fjet->GetEntriesFast()==0) { continue; }
        
        double leadPt = 0., leadEta, leadPhi;
        for (auto jet : dat.iter_Fjet){
            if (jet.pt > leadPt && jet.pt > 4.) { // take highest-pT jet
                leadPt = jet.pt;
                leadEta = jet.eta;
                leadPhi = jet.phi;
            }
        }
        double xsec_wt = pyth6_10x_Xsec( dat.pthat_bin );
        if ( leadPt<4. ) {
            double mc_leadPt = 0., mc_leadEta, mc_leadPhi;
            for (auto jet : dat.iter_mc_Fjet){
                if (jet.pt > mc_leadPt && jet.pt > 4.) { // take highest-pT particle-level jet
                    mc_leadPt = jet.pt;
                    mc_leadEta = jet.eta;
                    mc_leadPhi = jet.phi;
                }
            }
            if ( mc_leadPt<4. ) { continue; }
            hMissed->Fill(mc_leadPt, dat.mu_event->BBC_Ein, xsec_wt); // if part jet and no det jet: fill missed
            continue;
        }
//        cout<<trigPhi<<" \t "<<leadPhi<<endl;
        if( !is_phi_p4match( trigPhi, leadPhi) ) { continue; }  // require trigger in det-level leading jet (or recoil region)
        
        double mc_leadPt = 0., mc_leadEta, mc_leadPhi;
        for (auto jet : dat.iter_mc_Fjet){
            if (jet.pt > mc_leadPt && jet.pt > 4.) { // take highest-pT particle-level jet
                mc_leadPt = jet.pt;
                mc_leadEta = jet.eta;
                mc_leadPhi = jet.phi;
            }
        }
        if ( mc_leadPt<4. ) {
            hFake->Fill(leadPt, dat.mu_event->BBC_Ein, xsec_wt);         // FAKE JET
            continue;
        }
        
        if ( matched_jets( leadPhi, mc_leadPhi, leadEta, mc_leadEta ) ) {
            hMatchedJets->Fill( mc_leadPt, leadPt, dat.mu_event->BBC_Ein, xsec_wt );        // FILL RESPONSE
        }
        else {            // FILL MISS AND MATCH
            hMissed->Fill(mc_leadPt, dat.mu_event->BBC_Ein, xsec_wt);
            hFake->Fill(leadPt, dat.mu_event->BBC_Ein, xsec_wt);
        }

    } // / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
    

    // Write histograms here
    hMatchedJets->Write();
    hMissed->Write();
    hFake->Write();

    // Wrap-up work here:
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
