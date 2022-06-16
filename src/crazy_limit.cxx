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
void crazy_limit(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"crazy_limit\"" << endl;
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
    const char* edge_file = "/gpfs/loomis/home.grace/djs232/root_macros/io_lib"
                            "/pAu2015_common/bin_edges.txt";

    ioMsgTree msg{false};
    msg.slurp_file("src/crazy_limit.cxx");
    msg.write();
    // normalized per event
    // cut the corner cases

    ioIntSet  badtow  { 
        "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/bad_tower_iter.list",0,false};

    /* int* runid   = &dat.mu_event->runId; */
    /* int* eventid = &dat.mu_event->eventId; */
    get_good_towers tower{dat};
    oiJetMaker meas_jets   {{{"calc_areas",0}}};
    oiJetMaker true_jets   {{{"calc_areas",0}}};
    vector<double> M_jetpt_limits = ioReadValVec("/gpfs/loomis/home.grace/djs232/w2021"
        "/pAu2015_common/M_jet_limit.txt");
    vector<long int> cnt_0   (9,0); // cnt initially
    vector<long int> cnt_1   (9,0); // cnt after cutting on 30
    vector<long int> cnt_2   (9,0); // cnt after cutting too high Measured
    vector<long int> cnt_3   (9,0); // cnt after matching

    ioXsec Xsec   { "/gpfs/loomis/home.grace/djs232"
        "/root_macros/io_lib/pAu2015_common/pAuXsection.txt" };
    const char* JES_JER_file { 
        "/gpfs/loomis/home.grace/djs232/root_macros/io_lib/pAu2015_common/JES_JER.txt" };
    const char* JES_JER_pthb { 
        "/gpfs/loomis/home.grace/djs232/root_macros/io_lib/pAu2015_common/JER_JET_v_pthatbin.txt" };

    ioJetMatcher matcher { 
        "F_ones", Xsec, edge_file, 
        "jetpt_bigbyones", "jetpt_bigbyones",
    {{
        "pthb_Mlimit_file", JES_JER_pthb,
        "JER_limit", 3.,
        "match_fake_limit",8.,
        "fake_limit", 8.,
        //
        /* "JES_JER_file", JES_JER_file, */
        /* "JER_limit", 3., */
        /* "pt_true","jet_pt", */
        /* "match_fake_limit",8., */
        /* "fake_limit", 8., */
        /* "ratio_AtoB", 0.10, */
        //
        "debug",true
    }}};
    
    // find the pt bounds of outliers
    const char* limit_file = "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/"
        "embedding_pt_limits.txt";
    vector<double> limitpt_FullTrue { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_truth"}}) };

    ioBinVec jbin {{ -1., -1., 124, 30. }};
    ioBinVec pthbins { Xsec.pthatbins };

    // keep track of the third hardest jets in embedding
    TH2D* third_jets = new TH2D( "third_jets", "third_jets;pT of third True jet; #hat{#it{p}}_{T,bin}",
        jbin, jbin,  pthbins, pthbins );

    int z {0};
    TH1D ZDCx  { "zdcX",";ZDCx;N_{events}",    100,0.,35000.};
    while (dat.next()) {
        ++cnt_0[dat.pthat_bin];
        if (cut30_event(dat)) continue;
        ++cnt_1[dat.pthat_bin];

        int pthatbin { dat.pthat_bin };

        // make full measured jets
        meas_jets.reset();
        for (auto track : dat.iter_track) {
            meas_jets.add_particle(track.pt, track.eta, track.phi);
        }
        while (tower.next()) {
            meas_jets.add_particle(tower.Et, tower()->eta, tower()->phi);
        }
        meas_jets.cluster_jets();

        // make the T-full jets
        true_jets.reset();
        for (auto track : dat.iter_mcTr){
            true_jets.add_particle(track.pt, track.eta, track.phi);
        }
        for (auto part : dat.iter_mcNeut){
            true_jets.add_particle(part.pt, part.eta, part.phi);
        }

        // cut on full jet bounds
        true_jets.cluster_jets();
        if (true_jets.pseudojets.size()>0 
                && true_jets.pseudojets[0].perp()>limitpt_FullTrue[pthatbin]) continue;
        ++cnt_2[dat.pthat_bin];
       
        // see how close to the weighting things stay...
        ZDCx.Fill(dat.mu_event->ZDCx,Xsec.Xsec(dat.pthat_bin));

        auto third_pt { meas_jets.pseudojets.size() > 2 ? meas_jets.pseudojets[2].perp() : -1. };
        double pthat_val { (0.5)*(Xsec.pthatbins[dat.pthat_bin]+
                       Xsec.pthatbins[dat.pthat_bin+1]) };
        third_jets->Fill(third_pt,pthat_val);

        // fill in the jetmatchers
        for (auto& jet : meas_jets.pseudojets) {
            matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
        }
        for (auto& jet : true_jets.pseudojets) {
            matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
        }
        if (!matcher.do_matching(dat.pthat_bin)) ++cnt_3[dat.pthat_bin];
    }
    matcher.write();
    ostringstream os;
    os << " Runs cut for (A) >30GeV {track,tower}" << endl
       << " Runs cut for (B) Truth jet outlier" << endl
       << "              (C) measured or fake pT out of bounds" << endl
       << "              total: A||B||C " << endl
       << " Data is \"cut\"(no.lost:rat.lost) " << endl
       << Form(" %8s  %-21s %-21s %-21s %-21s","start",">30GeV","T-outlier","pT-bounds","total") << endl;
    for (int i{0}; i<cnt_0.size(); ++i) {
        int i0 {(int) cnt_0[i]};
        int i1 {(int) cnt_1[i]};
        int i2 {(int) cnt_2[i]};
        int i3 {(int) cnt_3[i]};
        os << Form("%8i %s %s %s %s",i0, io_cutdiff(i0,i1),
                io_cutdiff(i1,i2), io_cutdiff(i2,i3), io_cutdiff(i0, i3)) << endl;
    }
    cout << os.str();
    dat.log << os.str();

    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
