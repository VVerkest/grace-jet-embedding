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
void rat_rooResF(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"rat_rooResF\"" << endl;
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
    msg.slurp_file("src/rat_rooResF.cxx");
    msg.write();
    // normalized per event
    // cut the corner cases

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
    vector<long int> bounds_cut_cnt    (9,0);
    vector<long int> bounds_keep_cnt    (9,0);

    // Jet Matchers for F/C  : True/Measured x Full/Charged
    ioXsec Xsec   { "/gpfs/loomis/home.grace/djs232"
        "/root_macros/io_lib/pAu2015_common/pAuXsection.txt" };
    const char* file_bounds { 
        "/gpfs/loomis/home.grace/djs232/root_macros/io_lib/pAu2015_common/JES_JER.txt" };
    /* ioJetMatcherX F_matcher { "F", Xsec, edge_file, "pt_by1Raghav_M", "pt_by1Raghav_T" }; */
    /* ioJetMatcherX C_matcher { "C", Xsec, edge_file, "pt_by1Raghav_M", "pt_by1Raghav_T" }; */

    ioJetMatcher F_ones_match { 
        "F_ones", Xsec, edge_file, 
        "jetpt_byones_M", "jetpt_byones_T",
    {{
        "match_bounds_file", file_bounds, 
        "pt_true", "pt_jet_bin_upbound",
        "pt_measured", "pt_jet+4*JER+9.0",
        "ratio_AtoB", 0.984 
    }}};

    ioJetMatcher F_prelim_match { 
        "F_prelim", Xsec, edge_file, 
        "jetpt_prelim_M", "jetpt_prelim_T",
    {{
        "match_bounds_file", file_bounds, 
        "pt_true", "pt_jet_bin_upbound",
        "pt_measured", "pt_jet+4*JER+9.0",
        "ratio_AtoB", 0.984 
    }}};
    /* cout << " LION size: " << F_matcher.out_of_match_bounds.size << endl; */
    /* for (auto i : F_matcher.out_of_match_bounds.X) cout << " X: " << i << endl; */
    /* for (auto i : F_matcher.out_of_match_bounds.Y) cout << " Y: " << i << endl; */

    ioJetMatcher C_matcher { "C", Xsec, edge_file, "jetpt_bigbyones", "jetpt_bigbyones" };

    // find the pt bounds of outliers
    const char* limit_file = "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/"
        "embedding_pt_limits.txt";
    vector<double> limitpt_FullTrue { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_truth"}}) };
    vector<double> limitpt_FullMeas { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_measured"}}) };
    vector<double> limitpt_ChTrue { ioReadValVec(limit_file,   {{"tag","emb_C_jet_bound_truth"}}) };
    vector<double> limitpt_ChMeas { ioReadValVec(limit_file,   {{"tag","emb_C_jet_bound_measured"}}) };

    int z {0};
    TH1D ZDCx  { "zdcX",";ZDCx;N_{events}",    100,0.,35000.};
    while (dat.next()) {

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
        for (auto track : dat.iter_track) {
            M_F_jets.add_particle(track.pt, track.eta, track.phi);
        }
        while (tower.next()) {
            M_F_jets.add_particle(tower.Et, tower()->eta, tower()->phi);
        }
        M_F_jets.cluster_jets();

        // cut on full jets outside of limits
        /* if (M_F_jets.pseudojets.size()>0 */ 
                /* && M_F_jets.pseudojets[0].perp()>limitpt_FullMeas[pthatbin]) continue; */

        // cut events with cut in full measured jet limits
        /* double limit { M_jetpt_limits[dat.pthat_bin] }; */
        /* bool skip_event{false}; */
        /* for (auto jet : M_F_jets.pseudojets) { */
        /*     if (jet.perp() >= limit) { */
        /*         skip_event=true; */
        /*         break; */
        /*     } */
        /* } */
        /* if (skip_event) continue; */


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
            F_ones_match.addjet_reco(jet.eta(), jet.phi(), jet.perp());
            F_prelim_match.addjet_reco(jet.eta(), jet.phi(), jet.perp());
        }
        for (auto& jet : T_F_jets.pseudojets) {
            F_ones_match.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            F_prelim_match.addjet_MC(jet.eta(), jet.phi(), jet.perp());
        }
        F_prelim_match.do_matching(dat.pthat_bin);
        if (F_ones_match.do_matching(dat.pthat_bin)) {
            bounds_cut_cnt[dat.pthat_bin] += 1;
        } else {
            bounds_keep_cnt[dat.pthat_bin] += 1;
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
    F_ones_match.write();
    F_prelim_match.write();
    C_matcher.write();
    cout << "   no-cut   cut-30   cut-OL Ncut-30(percent) Ncut-OL(percent)" << endl;
    dat.log << "   no-cut   cut-30   cut-OL Ncut-30(percent) Ncut-OL(percent)" << endl;

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
        dat.log << Form(" %8i %8i %8i   %8i(%4.2f)   %8i(%4.2f)",
            pre, mid, post, (pre-mid), float(pre-mid)/(pre),
                            (mid-post),float(mid-post)/mid) << endl;
    }

    cout << " Cut for out of bounds: " << endl;
    cout << "        N-cut N-final    ratio cut" << endl;
    for (int i{0}; i<cnt_pthatbins.size(); ++i) {
        int post = cnt_pthatbins[i];
        int out  = bounds_cut_cnt[i];
        cout << Form("%10i %10i %10.2f", 
            post, out, double(out)/(out+post)) << endl;
    }

    dat.log << " Cut for out of bounds: " << endl;
    dat.log << "        N-cut N-final    ratio cut" << endl;
    for (int i{0}; i<cnt_pthatbins.size(); ++i) {
        int post = cnt_pthatbins[i];
        int out  = bounds_cut_cnt[i];
        dat.log << Form("%10i %10i %10.2f", 
            post, out, double(out)/(out+post)) << endl;
    }

    dat.log << " Cut for out of bounds: " << endl;
    dat.log << "    N-keep      N_cut  ratio cut" << endl;
    cout << " Cut for out of bounds: " << endl;
    cout << "        N-keep N-cup    ratio cut" << endl;
    for (int i{0}; i<cnt_pthatbins.size(); ++i) {
        int keep = bounds_keep_cnt[i];
        int cut  = bounds_cut_cnt[i];
        dat.log << Form("%10i %10i %10.4f", 
            keep, cut, double(cut)/(keep+cut)) << endl;
        cout    << Form("%10i %10i %10.4f", 
            keep, cut, double(cut)/(keep+cut)) << endl;
    }

    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
