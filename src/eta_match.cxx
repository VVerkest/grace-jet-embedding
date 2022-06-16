#include "events.h"
#include <sstream>
#include <fstream>


#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TreeObj.h"
#include "friend_fnc.h"
#include "ioClass.h"
/* #include "vzZDCx_correlator.h" */

#include "ioXsec_pAu2015.h"
/* #include "ioJetMatcher.h" */
#include "ioJetMatcherArray.h"
#include "oiJetMaker.h"
/* #include "ioBins.h" */

#include "ioTHnSparse.h"
#include "RooUnfoldResponse.h"


using namespace std;
void eta_match(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"eta_match\"" << endl;
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

    ioMsgTree msg{false};
    msg.slurp_file("src/eta_match.cxx");
    msg.write();

    get_good_towers tower{dat};

    oiJetMaker T_F_jets   {{{"calc_areas",0,"jetrap",10.0}}};
    const char* edge_file = "/gpfs/loomis/home.grace/djs232/root_macros/io_lib"
                            "/pAu2015_common/bin_edges.txt";
    ioXsec Xsec   { "/gpfs/loomis/home.grace/djs232"
        "/root_macros/io_lib/pAu2015_common/pAuXsection.txt" };

    ioBinVec bin_eta {{-2.,-2, 40, 2.}};
    ioBinVec pt_bins  {{ 0., 8.,10.,10.,18, 100. }};
    cout << " pt_bins " << pt_bins << endl;
    cout << " bin_eta " << bin_eta << endl;
    TH2D AJ_eta { "AJ_eta",";pT_{lead};pT_{rec}",pt_bins, pt_bins, pt_bins, pt_bins,  };
    ioIntSet  badtow  { 
        "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/bad_tower_iter.list",0,false};
    
    const double one_eighth_pi { TMath::Pi()/8. };
    const double seven_eighth_pi { 7*TMath::Pi()/8. };
    while (dat.next()) {
        if (cut30_event(dat)) continue;
        int pthatbin { dat.pthat_bin };

        // make the T-full jets
        T_F_jets.reset();
        for (auto track : dat.iter_mcTr){
            T_F_jets.add_particle(track.pt, track.eta, track.phi);
        }
        for (auto part : dat.iter_mcNeut){
            T_F_jets.add_particle(part.pt, part.eta, part.phi);
        }
        T_F_jets.cluster_jets();

        mupicoTower* trigger=nullptr;
        double trigger_pt{0};
        for (auto& tower : Iter_GoodTowerHits{dat.tca_tower, badtow}) {
            if (tower.Et > trigger_pt) {
                trigger_pt = tower.Et;
                trigger = &tower;
            }
        }
        if (trigger_pt < 8) continue;
        double phi_trig { trigger->phi };

        double pt_tr_side = 0;
        double pt_rec  = 0;

        if (trigger == nullptr) continue;
        auto& jets = T_F_jets.pseudojets;
        for (auto jet : jets) {
            if (pt_tr_side != 0 && pt_rec != 0) break;
            double pt { jet.perp() };
            if (TMath::Abs(jet.eta()) >0.6) continue;
            double absDphi = (io_absDphi(jet.phi(), phi_trig));
            if (pt_tr_side==0 && (absDphi < one_eighth_pi)) {
                pt_tr_side = pt;
            }
            if (pt_rec==0  && (absDphi > seven_eighth_pi)) {
                pt_rec = pt;
            }
        };
        AJ_eta.Fill( pt_tr_side, pt_rec, Xsec.Xsec(pthatbin) );
    }
    AJ_eta.Write();
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
