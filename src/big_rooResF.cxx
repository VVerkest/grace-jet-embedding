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
void big_rooResF(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"big_rooResF\"" << endl;
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
    msg.slurp_file("src/big_rooResF.cxx");
    msg.write();
    // normalized per event

    ioIntSet  badtow  { 
        "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/bad_tower_iter.list",0,false};

    /* int* runid   = &dat.mu_event->runId; */
    /* int* eventid = &dat.mu_event->eventId; */
    get_good_towers tower{dat};
    oiJetMaker M_F_jets   {{{"calc_areas",0}}};
    oiJetMaker T_F_jets   {{{"calc_areas",0}}};
    vector<double> M_jetpt_limits = ioReadValVec("/gpfs/loomis/home.grace/djs232/w2021"
        "/pAu2015_common/M_jet_limit.txt");
    vector<long int> cnt_pthatbins (9,0);
    vector<long int> precut_cnt    (9,0);
    vector<long int> midcut_cnt    (9,0);
    auto bbc_limits = ioReadValVec ( edge_file, "BBCES_4bin" );

    // Jet Matchers for F/C  : True/Measured x Full/Charged
    ioXsec Xsec   { "/gpfs/loomis/home.grace/djs232"
        "/root_macros/io_lib/pAu2015_common/pAuXsection.txt" };
    /* ioJetMatcher F_matcher { "F", Xsec, edge_file, "pt_by1Raghav_M", "pt_by1Raghav_T" }; */
    /* ioJetMatcher C_matcher { "C", Xsec, edge_file, "pt_by1Raghav_M", "pt_by1Raghav_T" }; */
    ioJetMatcher F_matcher { "F", Xsec, edge_file, "jetpt_bigbyones", "jetpt_bigbyones" };
    ioJetMatcher C_matcher { "C", Xsec, edge_file, "jetpt_bigbyones", "jetpt_bigbyones" };
    ioJetMatcher matcher_hiEA { 
        "match_hiEA", Xsec, edge_file, 
        "jetpt_bigbyones", "jetpt_bigbyones" };
    ioJetMatcher matcher_loEA { 
        "match_loEA", Xsec, edge_file, 
        "jetpt_bigbyones", "jetpt_bigbyones" };

    // find the pt bounds of outliers
    const char* limit_file = "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/"
        "embedding_pt_limits.txt";
    vector<double> limitpt_FullTrue { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_truth"}}) };
    vector<double> limitpt_FullMeas { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_measured"}}) };
    vector<double> limitpt_ChTrue { ioReadValVec(limit_file,   {{"tag","emb_C_jet_bound_truth"}}) };
    vector<double> limitpt_ChMeas { ioReadValVec(limit_file,   {{"tag","emb_C_jet_bound_measured"}}) };

    TH1D ZDCx  { "zdcX",";ZDCx;N_{events}",    100,0.,35000.};
    while (dat.next()) {

        /* cout << "size : " << dat.size_mc_Fjet() << "  " << dat.Fjet_njets << endl; */
        /* continue; */

        precut_cnt[dat.pthat_bin] += 1;
        /* cout << " dat " << z++ << endl; */
        if (cut30_event(dat)) continue;
        midcut_cnt[dat.pthat_bin] += 1;

        int pthatbin { dat.pthat_bin };
        // cut on full jets
        if (dat.size_mc_Cjet()>0 && dat.get_mc_Cjet(0)->pt > limitpt_ChTrue[pthatbin]) continue;
        /* if (dat.size_Cjet()>0 && dat.get_Cjet(0)->pt > limitpt_ChMeas[pthatbin]) continue; */


        // make full measured jets
        M_F_jets.reset();
        for (auto track : Iter_GoodTracks(dat.tca_track)) {
            M_F_jets.add_particle(track.pt, track.eta, track.phi);
        }
        while (tower.next()) {
            M_F_jets.add_particle(tower.Et, tower()->eta, tower()->phi);
        }
        M_F_jets.cluster_jets();

        // make the T-full jets
        T_F_jets.reset();
        for (auto track : dat.iter_mcTr){
            T_F_jets.add_particle(track.pt, track.eta, track.phi);
        }
        for (auto part : dat.iter_mcNeut){
            T_F_jets.add_particle(part.pt, part.eta, part.phi);
        }

        // cut on full jet bounds
        T_F_jets.cluster_jets();
        if (T_F_jets.pseudojets.size()>0 
                && T_F_jets.pseudojets[0].perp()>limitpt_FullTrue[pthatbin]) continue;
       
        // see how close to the weighting things stay...
        cnt_pthatbins[dat.pthat_bin] += 1;
        ZDCx.Fill(dat.mu_event->ZDCx,Xsec.Xsec(dat.pthat_bin));

        // fill in the jetmatchers
        for (auto& jet : M_F_jets.pseudojets) {
            F_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
        }
        for (auto& jet : T_F_jets.pseudojets) {
            F_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
        }
        F_matcher.do_matching(dat.pthat_bin); 

        if (dat.mu_event->BBC_Ein > bbc_limits[3]) {
            for (auto& jet : M_F_jets.pseudojets) {
                matcher_hiEA.addjet_reco(jet.eta(), jet.phi(), jet.perp());
            }
            for (auto& jet : T_F_jets.pseudojets) {
                matcher_hiEA.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            }
            matcher_hiEA.do_matching(dat.pthat_bin);
        } else if (dat.mu_event->BBC_Ein < bbc_limits[2]) {
            for (auto& jet : M_F_jets.pseudojets) {
                matcher_loEA.addjet_reco(jet.eta(), jet.phi(), jet.perp());
            }
            for (auto& jet : T_F_jets.pseudojets) {
                matcher_loEA.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            }
            matcher_loEA.do_matching(dat.pthat_bin);
        }

        // fill charged jets (already are pre-made in embedding,
        // and don't change with bad tower list)
        for (auto jet : dat.iter_mc_Cjet) {
            C_matcher.addjet_MC(jet.eta, jet.phi, jet.pt);
        }
        for (auto jet : dat.iter_Cjet) {
            C_matcher.addjet_reco(jet.eta, jet.phi, jet.pt);
        }
        C_matcher.do_matching(dat.pthat_bin);
    }
    F_matcher.write();
    C_matcher.write();
    matcher_loEA.write();
    matcher_hiEA.write();
    cout << "   no-cut   cut-30   cut-OL Ncut-30(percent) Ncut-OL(percent)" << endl;

    /* dat.log << Form(" %-8s          ->  %-7s %-7s  %-8s", */
            /* "pre-cut","post-cut","delta","ratio-cut") << endl; */
    /* cout    << Form(" %-8s          ->  %-7s %-7s  %-8s", */
            /* "pre-cut","post-cut","delta","ratio-cut") << endl; */
    for (int i{0}; i<cnt_pthatbins.size(); ++i) {
        int pre  {(int)precut_cnt[i]};
        int mid  {(int)midcut_cnt[i]};
        int post {(int)cnt_pthatbins[i]};
        cout << Form(" %8i %8i %8i   %8i(%4.2f)   %8i(%4.2f)",
            pre, mid, post, (pre-mid), float(pre-mid)/(pre),
                            (mid-post),float(mid-post)/mid) << endl;
        /* dat.log << Form(" %-8i ->%7i->  %-7i  %-7i  %8.2f", */
                /* pre, mid, post, (post-pre),(1.-float(post)/pre)) << endl; */
        /* cout << Form(" %-8i ->%7i->  %-7i  %-7i  %8.2f", */
                /* pre, mid, post, (post-pre),(1.-float(post)/pre)) << endl; */
    }
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
