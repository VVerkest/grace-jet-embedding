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
/* #include "ioJetMatcher.h" */
#include "ioJetMatcher100.h"
#include "oiJetMaker.h"
/* #include "ioBins.h" */

#include "ioTHnSparse.h"
#include "RooUnfoldResponse.h"

// use the same settings as thesis_emb, but only for the single set of
// binning (the one used for the data) with : 1. 5% track eff. cut, 2. 30 M/S cut (A to B)
// 
// uncertainties to add: (from Isaac Mooney's jet analysis note)
// 1. IP2(6)Iteration: iterations in Bayesian unfolding instead of 4 (fine, no problem! -- not done here)
// 2. TS - tower scale uncertainty of 3.8%
// 3. TU - tracking efficiency uncertainty 4%
// 4. HC50 -- hadronic subtraction uncertainty, use 50% instead of 100%
// 5. GS -- generator smearing of response matrix (10%?)

using namespace std;
void sys_err(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"sys_err\"" << endl;
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
    ioMsgTree msg{false};
    msg.slurp_file("src/sys_err.cxx");
    msg.write();
    // normalized per event

    ioIntSet  badtow  { 
        "/vast/palmer/home.grace/djs232/w2021/pAu2015_common/bad_tower_iter.list",0,false};

    /* int* runid   = &dat.mu_event->runId; */
    /* int* eventid = &dat.mu_event->eventId; */
    get_good_towers tower{dat};

    // regular jets
    oiJetMaker T_jets   {{{"calc_areas",0}}}; // truth jets
    oiJetMaker M_jets   {{{"calc_areas",0}}}; // measured jets
    // TS offset jet
    oiJetMaker TS_M_jets   {{{"calc_areas",0}}};
    // TU offset jet
    oiJetMaker TU_M_jets   {{{"calc_areas",0}}};
    // HC50 offset jet
    oiJetMaker HC50_M_jets   {{{"calc_areas",0}}};

    const char* Mlimit_file = "/vast/palmer/home.grace/djs232/root_macros/io_lib/pAu2015_common/JER_JET_v_pthatbin.txt";
    /* ioPtrDbl lim_ptbins   {"thesis_data/JER_JER.list","jet_pt"     }; */
    /* ioPtrDbl lim_M_limits {"thesis_data/JER_JER.list","pt_M_limit" }; */
    /* ioXYbounder {lim_ptbins.vec, lim_M_limits.vec}; */

    const char* edge_file = "/vast/palmer/home.grace/djs232/root_macros/io_lib"
                            "/pAu2015_common/bin_edges.txt";
    auto bbc_limits = ioReadValVec ( edge_file, "BBCES_4bin" );

    /* ioJetSpectraSparse sparse(edge_file); */

    // Jet Matchers for F/C  : True/Measured x Full/Charged
    ioXsec Xsec   { "/vast/palmer/home.grace/djs232"
        "/root_macros/io_lib/pAu2015_common/pAuXsection.txt" };


    const char* bin_file = "./thesis_data/thesis_bins.txt";
    /* vector<string> bin_names{ "pt"}; */
    /* vector<string> bin_tags { "one_pt" }; */
    /* vector<string> bin_names{ "EAbbc_ratio","Etbinning","BbcTpc","EAtpc_short","one"}; */
    /* vector<string> bin_tags { "dphi_EA_bbc_ratio", "EAbbc_Etbinning", */ 
                              /* "BbcTpc_comp", "draw_dphi_EAtpc_short","one_pt"}; */
    vector<string> bin_names{ "one"};
    vector<string> bin_tags_M { "one_pt"};
    vector<string> bin_tags_T { "one_pt"};

    ioJetMatcher100 F_matcher    { &Xsec };
    ioJetMatcher100 TU_matcher   { &Xsec,"TU_" };
    ioJetMatcher100 TS_matcher   { &Xsec,"TS_" };
    ioJetMatcher100 HC50_matcher { &Xsec,"HC50_" };

    // find the pt bounds of outliers
    const char* limit_file = "/vast/palmer/home.grace/djs232/w2021/pAu2015_common/"
        "embedding_pt_limits.txt";
    vector<double> limitpt_FullTrue { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_truth"}}) };
    vector<double> limitpt_FullMeas { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_measured"}}) };
    vector<double> limitpt_ChTrue   { ioReadValVec(limit_file, {{"tag","emb_C_jet_bound_truth"}}) };
    vector<double> limitpt_ChMeas   { ioReadValVec(limit_file, {{"tag","emb_C_jet_bound_measured"}}) };

    /* cout << " z2 " << endl; */
    TH1D ZDCx  { "zdcX",";ZDCx;N_{events}", 100,0.,35000.};
    TRandom3 _rand{};
    
    /* cout << " z3 " << endl; */
    while (dat.next()) {
        /* precut_cnt[dat.pthat_bin] += 1; */
        if (cut30_event(dat)) continue;
        /* midcut_cnt[dat.pthat_bin] += 1; */

        int pthatbin { dat.pthat_bin };

        M_jets.reset();
        TU_M_jets.reset();
        TS_M_jets.reset();
        HC50_M_jets.reset();

        //add tracks
        /* cout << " Uniform : " << _rand.Uniform() << endl; */
        /* int n = 0; */
        /* for (auto track : Iter_GoodTracks(dat.tca_track)) { */
        /*     cout << " track("<<n<<")  DCA: " << track.dcaXY << endl; */
        /*     n += 1; */
        /*     M_jets.add_particle(track.pt, track.eta, track.phi); */
        /*     if (_rand.Uniform() > 0.04) TU_M_jets.add_particle(track.pt, track.eta, track.phi); */
        /*     TS_M_jets.add_particle(track.pt, track.eta, track.phi); */
        /*     HC50_M_jets.add_particle(track.pt, track.eta, track.phi); */
        /* } */
        /* continue; */
        
        //add towers
        for (auto tower : Iter_GoodCorrTowers(dat.tca_tower,badtow)) {
            double EtCorr { tower.EtCorr() };
            M_jets.add_particle   (EtCorr, tower.eta, tower.phi);
            TU_M_jets.add_particle(EtCorr, tower.eta, tower.phi);
        }
        // need separate loop for towers that would otherwise be skipped
        for (auto tower : Iter_GoodTowerHits(dat.tca_tower,badtow)) {
            double Et = (tower.Et - 0.5*tower.Et_hadroncorr);
            if (Et>=0.2) HC50_M_jets.add_particle(Et, tower.eta, tower.phi);
        };
        for (auto tower : Iter_GoodTowerHits(dat.tca_tower,badtow)) {
            double Et { _rand.Gaus(1,0.038) * tower.Et - tower.Et_hadroncorr };
            if (Et>=0.2) TS_M_jets.add_particle(Et, tower.eta, tower.phi);
        };

        M_jets.cluster_jets();
        TU_M_jets.cluster_jets();
        TS_M_jets.cluster_jets();
        HC50_M_jets.cluster_jets();

        // make the T-full jets
        T_jets.reset();
        for (auto track : dat.iter_mcTr){
            T_jets.add_particle(track.pt, track.eta, track.phi);
        }
        for (auto part : dat.iter_mcNeut){
            T_jets.add_particle(part.pt, part.eta, part.phi);
        }
        /* cout << " z1 " << endl; */

        // cut on full jet bounds
        T_jets.cluster_jets();
        /* for (auto& jet : M_jets.pseudojets) { */
        /*     no_cull.addjet_reco(jet.eta(), jet.phi(), jet.perp()); */
        /* } */
        /* for (auto& jet : T_jets.pseudojets) { */
        /*     no_cull.addjet_MC(jet.eta(), jet.phi(), jet.perp()); */
        /* } */
        /* no_cull.do_matching(dat.pthat_bin); */ 

        if (dat.size_mc_Cjet()>0 && dat.get_mc_Cjet(0)->pt > limitpt_ChTrue[pthatbin]) continue;
        if (T_jets.pseudojets.size()>0 
                && T_jets.pseudojets[0].perp()>limitpt_FullTrue[pthatbin]) continue;
       


        // see how close to the weighting things stay...
        /* cnt_pthatbins[dat.pthat_bin] += 1; */
        ZDCx.Fill(dat.mu_event->ZDCx,Xsec.Xsec(dat.pthat_bin));

        // fill in the jetmatchers -- measured jets
        for (auto& jet : M_jets.pseudojets) {
            F_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
        }
        for (auto& jet : TU_M_jets.pseudojets) {
            TU_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
        }
        for (auto& jet : TS_M_jets.pseudojets) {
            TS_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
        }
        for (auto& jet : HC50_M_jets.pseudojets) {
            HC50_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
        }

        // 
        for (auto& jet : T_jets.pseudojets) {
            F_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            TU_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            TS_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            HC50_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
        }
        F_matcher.do_matching(dat.pthat_bin); 
        TU_matcher.do_matching(dat.pthat_bin); 
        TS_matcher.do_matching(dat.pthat_bin); 
        HC50_matcher.do_matching(dat.pthat_bin); 
    }
    /* sparse.write(); */
    F_matcher.write();
    TU_matcher.write();
    TS_matcher.write();
    HC50_matcher.write();

    cout << "   no-cut   cut-30   cut-OL Ncut-30(percent) Ncut-OL(percent)" << endl;
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
