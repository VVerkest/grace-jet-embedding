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
/* #include "ioJetMatcherGoodBins.h" */
#include "oiJetMaker.h"
/* #include "ioBins.h" */

#include "RooUnfoldResponse.h"


using namespace std;
void bins70(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"bins70\"" << endl;
    istringstream options ( _options );
    int n_options = 0;
    string arg;
    // options >> arg;
    while (options >> arg) {
        cout    << " Option " << n_options << ":  " << arg << endl;
        dat.log << " Option " << n_options << ":  " << arg << endl;
        ++n_options;
    }

    /* cout << " a0 " << endl; */
    // Histogram declarations here:
    // TH1D hg {"hg", "a;b;c", 10, 0., 1.};
    const char* edge_file = "/gpfs/loomis/home.grace/djs232/root_macros/io_lib"
                            "/pAu2015_common/bin_edges.txt";
    const char* goodbinfile = "/gpfs/loomis/home.grace/djs232/root_macros/io_lib"
                            "/pAu2015_common/good_bins.txt";
    const char* goodbinfinal = "/gpfs/loomis/home.grace/djs232/root_macros/io_lib"
                            "/pAu2015_common/good_bins_iter.txt";

    ioMsgTree msg{false};
    msg.slurp_file("src/bins70.cxx");
    msg.write();
    // normalized per event
    // cut the corner cases

    ioIntSet  badtow  { 
        "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/bad_tower_iter.list",0,false};

    /* cout << " a1 " << endl; */
    Iter_GoodCorrTowers iter_goodtowers { dat.tca_tower, badtow };
    Iter_GoodTowerHits  iter_trighits   { dat.tca_tower, badtow };
    /* cout << " a2 " << endl; */

    oiJetMaker meas_jets   {{{"calc_areas",0}}};
    oiJetMaker true_jets   {{{"calc_areas",0}}};
    /* cout << " a3 " << endl; */

    vector<long int> cnt_0   (9,0); // cnt initially
    vector<long int> cnt_1   (9,0); // cnt after cutting on 30
    vector<long int> cnt_2   (9,0); // cnt after cutting too high Measured
    vector<long int> cnt_3   (9,0); // cnt after matching

    /* ioJetMatcherGoodBins matcher { */ 
    /*     "raw", edge_file, "jetpt_70byones", "jetpt_70byones", "", {{"debug",true}} }; */
    /* ioJetMatcherGoodBins matcher_trim { */ 
    /*     "trim", edge_file, "jetpt_70byones", "jetpt_70byones", goodbinfile, {{"debug",true}} }; */
    /* ioJetMatcherGoodBins matcher_final { */ 
    /*     "final", edge_file, "jetpt_70byones", "jetpt_70byones", goodbinfinal, {{"debug",true, */
    /*     "ratio_AtoB",0.02}} }; */
    /* /1* cout << " a4 " << endl; *1/ */
    
    /* ioBinVec jbin {{ -1., -1., 124, 30. }}; */
    /* ioXsec Xsec   { "/gpfs/loomis/home.grace/djs232" */
    /*     "/root_macros/io_lib/pAu2015_common/pAuXsection.txt" }; */
    /* ioBinVec pthbins { Xsec.pthatbins }; */

    /* // keep track of the third hardest jets in embedding */
    /* TH2D* third_jets = new TH2D( "third_jets", "third_jets;pT of third True jet; #hat{#it{p}}_{T,bin}", */
    /*     jbin, jbin,  pthbins, pthbins ); */

    /* /1* cout << " a6 " << endl; *1/ */
    /* int z {0}; */
    /* TH1D ZDCx  { "zdcX",";ZDCx;N_{events}",    100,0.,35000.}; */
    /* while (dat.next()) { */
    /*     ++cnt_0[dat.pthat_bin]; */
    /*     if (cut30_event(dat)) continue; */
    /*     int pthatbin { dat.pthat_bin }; */
    /*     ++cnt_1[pthatbin]; */
    /*     ++cnt_2[pthatbin]; */

    /*     // make full measured jets */
    /*     meas_jets.reset(); */
    /*     for (auto track : dat.iter_track) { */
    /*         meas_jets.add_particle(track.pt, track.eta, track.phi); */
    /*     } */
    /*     for (auto tower : iter_goodtowers) { */
    /*         meas_jets.add_particle(tower.EtCorr(), tower.eta, tower.phi); */
    /*     } */
    /*     meas_jets.cluster_jets(); */

    /*     // make the T-full jets */
    /*     true_jets.reset(); */
    /*     for (auto track : dat.iter_mcTr){ */
    /*         true_jets.add_particle(track.pt, track.eta, track.phi); */
    /*     } */
    /*     for (auto part : dat.iter_mcNeut){ */
    /*         true_jets.add_particle(part.pt, part.eta, part.phi); */
    /*     } */
    /*     true_jets.cluster_jets(); */

    /*     // see how close to the weighting things stay... */
    /*     ZDCx.Fill(dat.mu_event->ZDCx,Xsec.Xsec(pthatbin)); */

    /*     auto third_pt { meas_jets.pseudojets.size() > 2 ? meas_jets.pseudojets[2].perp() : -1. }; */
    /*     double pthat_val { (0.5)*(Xsec.pthatbins[pthatbin]+ */
    /*                    Xsec.pthatbins[pthatbin+1]) }; */
    /*     third_jets->Fill(third_pt,pthat_val); */

    /*     // fill in the jetmatchers */
    /*     for (auto& jet : meas_jets.pseudojets) { */
    /*         /1* cout << " add: RECO: " << jet.perp() << endl; *1/ */
    /*         matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp()); */
    /*         matcher_trim.addjet_reco(jet.eta(), jet.phi(), jet.perp()); */
    /*         matcher_final.addjet_reco(jet.eta(), jet.phi(), jet.perp()); */
    /*     } */
    /*     for (auto& jet : true_jets.pseudojets) { */
    /*         /1* cout << " add: MC: " << jet.perp() << endl; *1/ */
    /*         matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp()); */
    /*         matcher_trim.addjet_MC(jet.eta(), jet.phi(), jet.perp()); */
    /*         matcher_final.addjet_MC(jet.eta(), jet.phi(), jet.perp()); */
    /*     } */

    /*     matcher.do_matching(pthatbin); */
    /*     matcher_trim.do_matching(pthatbin); */
    /*     if (!matcher_final.do_matching(pthatbin)) ++cnt_3[pthatbin]; */
    /* } */
    /* matcher_final.write(); */
    /* matcher_trim.write(); */
    /* matcher.write(); */
    /* ostringstream os; */
    /* os << " Runs cut for (A) >30GeV {track,tower}" << endl */
    /*    << " Runs cut for (B) Truth jet outlier" << endl */
    /*    << "              (C) measured or fake pT out of bounds" << endl */
    /*    << "              total: A||B||C " << endl */
    /*    << " Data is \"cut\"(no.lost:rat.lost) " << endl */
    /*    << Form(" %8s  %-21s %-21s %-21s %-21s","start",">30GeV","T-outlier","pT-bounds","total") << endl; */
    /* for (int i{0}; i<cnt_0.size(); ++i) { */
    /*     int i0 {(int) cnt_0[i]}; */
    /*     int i1 {(int) cnt_1[i]}; */
    /*     int i2 {(int) cnt_2[i]}; */
    /*     int i3 {(int) cnt_3[i]}; */
    /*     os << Form("%8i %s %s %s %s",i0, io_cutdiff(i0,i1), */
    /*             io_cutdiff(i1,i2), io_cutdiff(i2,i3), io_cutdiff(i0, i3)) << endl; */
    /* } */
    /* cout << os.str(); */
    /* dat.log << os.str(); */

    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
