#include "events.h"
#include <sstream>
#include <fstream>


#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TreeObj.h"
#include "friend_fnc.h"
#include "TRandom3.h"

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
void thesis_emb(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"thesis_emb\"" << endl;
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
    msg.slurp_file("src/thesis_emb.cxx");
    msg.write();
    // normalized per event

    ioIntSet  badtow  { 
        "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/bad_tower_iter.list",0,false};

    /* int* runid   = &dat.mu_event->runId; */
    /* int* eventid = &dat.mu_event->eventId; */
    get_good_towers tower{dat};

    oiJetMaker M_F_jets   {{{"calc_areas",0}}};
    oiJetMaker T_F_jets   {{{"calc_areas",0}}};

    oiJetMaker M_F_jets_5p   {{{"calc_areas",0}}}; // jets with 5% track efficiency penalty
    

    const char* Mlimit_file = "/gpfs/loomis/home.grace/djs232/root_macros/io_lib/pAu2015_common/JER_JET_v_pthatbin.txt";
    /* ioPtrDbl lim_ptbins   {"thesis_data/JER_JER.list","jet_pt"     }; */
    /* ioPtrDbl lim_M_limits {"thesis_data/JER_JER.list","pt_M_limit" }; */
    /* ioXYbounder {lim_ptbins.vec, lim_M_limits.vec}; */

    const char* edge_file = "/gpfs/loomis/home.grace/djs232/root_macros/io_lib"
                            "/pAu2015_common/bin_edges.txt";
    auto bbc_limits = ioReadValVec ( edge_file, "BBCES_4bin" );

    ioJetSpectraSparse sparse   (edge_file);
    ioJetSpectraSparse sparse_5p(edge_file,"5p");

    // Jet Matchers for F/C  : True/Measured x Full/Charged
    ioXsec Xsec   { "/gpfs/loomis/home.grace/djs232"
        "/root_macros/io_lib/pAu2015_common/pAuXsection.txt" };


    const char* bin_file = "./thesis_data/thesis_bins.txt";
    vector<string> bin_names{ "EAbbc_ratio","Etbinning","BbcTpc","EAtpc_short","one"};
    vector<string> bin_tags { "dphi_EA_bbc_ratio", "EAbbc_Etbinning", 
                              "BbcTpc_comp", "draw_dphi_EAtpc_short","one_pt"};
    /* ioJetMatcherArray C_matcher { "C", Xsec, edge_file, bin_names, bin_tags, Mlimit_file }; */
    /* ioJetMatcherArray F_one { "one", Xsec, bin_file, {{"cull"}}, {{"one_pt"}}, Mlimit_file }; */
    /* F_one.write_9 = true; */
    /* F_one.cut_high_sigma {5.}; */
    /* F_one.cut_high_sigma_offset {8.}; */
    /* F_one.cull_n = 20; */
    /* ioJetMatcherArray no_cull { "nocull", Xsec, bin_file, {{"one_pt"}}, {{"one_pt"}}, Mlimit_file }; */
    /* no_cull.cull_n = 10; */
    /* no_cull.write_9 = true; */
    /* no_cull.apply_Mlimit = true; */
    ioJetMatcherArray trig_matcher   { "trig",  Xsec, bin_file, bin_names, bin_tags, bin_tags, Mlimit_file };
    ioJetMatcherArray rec_matcher    { "recoil",  Xsec, bin_file, bin_names, bin_tags, bin_tags, Mlimit_file };
    ioJetMatcherArray cull_matcher   { "cull",  Xsec, bin_file, bin_names, bin_tags, bin_tags, Mlimit_file };
    ioJetMatcherArray F_matcher      { "all",   Xsec, bin_file, bin_names, bin_tags, bin_tags, Mlimit_file };
    ioJetMatcherArray F_matcher_hiEA { "hiEA",  Xsec, bin_file, bin_names, bin_tags, bin_tags, Mlimit_file };
    ioJetMatcherArray F_matcher_loEA { "loEA",  Xsec, bin_file, bin_names, bin_tags, bin_tags, Mlimit_file };

    TH2D trig_dphi { "trig_dphi", ";|#phi_{trigger}-#phi_{jet}|;#it{p}_{T,jet}",
        50, 0., TMath::Pi(), 50,8.,58.};
    
    trig_matcher.write_9 = true;
    /* trig_matcher.cut_high_sigma = 5.; */
    /* trig_matcher.cut_high_sigma_offset = 8.; */
    trig_matcher.cull_n = 10;

    rec_matcher.write_9 = true;
    /* rec_matcher.cut_high_sigma = 5.; */
    /* rec_matcher.cut_high_sigma_offset = 8.; */
    rec_matcher.cull_n = 10;

    cull_matcher.write_9 = true;
    cull_matcher.cut_high_sigma = 5.;
    cull_matcher.cut_high_sigma_offset = 8.;
    cull_matcher.cull_n = 20;

    F_matcher.write_9 = true;
    F_matcher.cut_high_sigma = 5.;
    F_matcher.cut_high_sigma_offset = 8.;
    F_matcher.cull_n = 20;

    F_matcher_loEA.write_9 = true;
    F_matcher_loEA.cut_high_sigma = 5.;
    F_matcher_loEA.cut_high_sigma_offset= 8.;
    F_matcher_loEA.cull_n = 20;

    F_matcher_hiEA.write_9 = true;
    F_matcher_hiEA.cut_high_sigma = 5.;
    F_matcher_hiEA.cut_high_sigma_offset = 8.;
    F_matcher_hiEA.cull_n = 20;


    // find the pt bounds of outliers
    const char* limit_file = "/gpfs/loomis/home.grace/djs232/w2021/pAu2015_common/"
        "embedding_pt_limits.txt";
    vector<double> limitpt_FullTrue { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_truth"}}) };
    vector<double> limitpt_FullMeas { ioReadValVec(limit_file, {{"tag","emb_F_jet_bound_measured"}}) };
    vector<double> limitpt_ChTrue   { ioReadValVec(limit_file, {{"tag","emb_C_jet_bound_truth"}}) };
    vector<double> limitpt_ChMeas   { ioReadValVec(limit_file, {{"tag","emb_C_jet_bound_measured"}}) };

    TH1D ZDCx  { "zdcX",";ZDCx;N_{events}",    100,0.,35000.};

    TRandom3 r3{};

    while (dat.next()) {
        /* precut_cnt[dat.pthat_bin] += 1; */
        if (cut30_event(dat)) continue;
        /* midcut_cnt[dat.pthat_bin] += 1; */

        int pthatbin { dat.pthat_bin };

        // make full measured jets
        M_F_jets.reset();
        M_F_jets_5p.reset();
        for (auto track : Iter_GoodTracks(dat.tca_track)) {
            M_F_jets.add_particle(track.pt, track.eta, track.phi);
            if (r3.Uniform()<=0.60) M_F_jets_5p.add_particle(track.pt, track.eta, track.phi);
        }
        for (auto tower : Iter_GoodCorrTowers(dat.tca_tower,badtow)) {
            M_F_jets.add_particle(tower.EtCorr(), tower.eta, tower.phi);
            M_F_jets_5p.add_particle(tower.EtCorr(), tower.eta, tower.phi);
        }
        M_F_jets.cluster_jets();
        M_F_jets_5p.cluster_jets();

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
        /* for (auto& jet : M_F_jets.pseudojets) { */
        /*     no_cull.addjet_reco(jet.eta(), jet.phi(), jet.perp()); */
        /* } */
        /* for (auto& jet : T_F_jets.pseudojets) { */
        /*     no_cull.addjet_MC(jet.eta(), jet.phi(), jet.perp()); */
        /* } */
        /* no_cull.do_matching(dat.pthat_bin); */ 

        if (dat.size_mc_Cjet()>0 && dat.get_mc_Cjet(0)->pt > limitpt_ChTrue[pthatbin]) continue;
        if (T_F_jets.pseudojets.size()>0 
                && T_F_jets.pseudojets[0].perp()>limitpt_FullTrue[pthatbin]) continue;
       


        // see how close to the weighting things stay...
        /* cnt_pthatbins[dat.pthat_bin] += 1; */
        ZDCx.Fill(dat.mu_event->ZDCx,Xsec.Xsec(dat.pthat_bin));

        // fill in the jetmatchers
        for (auto& jet : M_F_jets.pseudojets) {
            F_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
            /* F_one.addjet_reco(jet.eta(), jet.phi(), jet.perp()); */
        }
        for (auto& jet : T_F_jets.pseudojets) {
            F_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            trig_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            rec_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            cull_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            /* F_one.addjet_MC(jet.eta(), jet.phi(), jet.perp()); */
        }
        F_matcher.do_matching(dat.pthat_bin); 
        /* F_one.do_matching(dat.pthat_bin); */ 
        if (dat.mu_event->BBC_Ein > bbc_limits[3]) {
            for (auto& jet : M_F_jets.pseudojets) {
                F_matcher_hiEA.addjet_reco(jet.eta(), jet.phi(), jet.perp());
            }
            for (auto& jet : T_F_jets.pseudojets) {
                F_matcher_hiEA.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            }
            F_matcher_hiEA.do_matching(dat.pthat_bin);
        } else if (dat.mu_event->BBC_Ein < bbc_limits[2]) {
            for (auto& jet : M_F_jets.pseudojets) {
                F_matcher_loEA.addjet_reco(jet.eta(), jet.phi(), jet.perp());
            }
            for (auto& jet : T_F_jets.pseudojets) {
                F_matcher_loEA.addjet_MC(jet.eta(), jet.phi(), jet.perp());
            }
            F_matcher_loEA.do_matching(dat.pthat_bin);
        }
        mupicoTower*  trigger=nullptr;
        double max_pt{0};
        ioMinMaxPtr filter{};
        for (auto& tower : Iter_GoodTowerHits{dat.tca_tower, badtow}) {
            if (tower.Et > max_pt) {
                max_pt = tower.Et;
                trigger = &tower;
            }
            /* filter(tower.Et, &tower); */
        }
        trigger = (max_pt > 8) ? trigger : nullptr;
        const double pieighth = TMath::Pi()/8.;
        const double receight = 7*TMath::Pi()/8.;
        /* trigger = (mupicoTower*) filter.max_ptr; */
        if (trigger != nullptr && trigger->Et > 8.) {
            double trig_phi = trigger->phi;
            sparse.weight = Xsec.Xsec(pthatbin);
            sparse.fill_trig(dat.mu_event->BBC_Ein, 1., trigger->Et,
                    dat.mu_event->ZDCx, 0);

            sparse_5p.weight = Xsec.Xsec(pthatbin);
            sparse_5p.fill_trig(dat.mu_event->BBC_Ein, 1., trigger->Et,
                    dat.mu_event->ZDCx, 0);
            for (auto& jet : M_F_jets.pseudojets) {
                double ADphi = TMath::Abs(io_dphi(trig_phi,jet.phi()));
                sparse.fill_jetpt_absDphi(jet.perp(), ADphi);
                if (ADphi<=pieighth)        trig_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
                else if (ADphi >= receight) rec_matcher .addjet_reco(jet.eta(), jet.phi(), jet.perp());
                if (ADphi<=pieighth) cull_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
                trig_dphi.Fill(ADphi, jet.perp());
                /* cout << "ADphi : " << ADphi << " " << jet.perp() << " eta: " << jet.eta() << endl; */
            }
            for (auto& jet : M_F_jets_5p.pseudojets) {
                double ADphi = TMath::Abs(io_dphi(trig_phi,jet.phi()));
                sparse_5p.fill_jetpt_absDphi(jet.perp(), ADphi);
            }
        }
        trig_matcher.do_matching(dat.pthat_bin);
        rec_matcher.do_matching(dat.pthat_bin);
        cull_matcher.do_matching(dat.pthat_bin);


        // fill charged jets (already are pre-made in embedding,
        // and don't change with bad tower list)
        /* for (auto jet : dat.iter_mc_Cjet) { */
        /*     C_matcher.addjet_MC(jet.eta, jet.phi, jet.pt); */
        /* } */
        /* for (auto jet : dat.iter_Cjet) { */
        /*     C_matcher.addjet_reco(jet.eta, jet.phi, jet.pt); */
        /* } */
        /* C_matcher.do_matching(dat.pthat_bin); */
    }
    F_matcher.write();
    trig_matcher.write();
    rec_matcher.write();
    if (false) { // FIXME
        cull_matcher.write();
        F_matcher_loEA.write();
        F_matcher_hiEA.write();
    }
    sparse.write();
    sparse_5p.write();
    trig_dphi.Write();
    /* cout << "   no-cut   cut-30   cut-OL Ncut-30(percent) Ncut-OL(percent)" << endl; */

    /* dat.log << Form(" %-8s          ->  %-7s %-7s  %-8s", */
            /* "pre-cut","post-cut","delta","ratio-cut") << endl; */
    /* cout    << Form(" %-8s          ->  %-7s %-7s  %-8s", */
            /* "pre-cut","post-cut","delta","ratio-cut") << endl; */
    /* for (int i{0}; i<cnt_pthatbins.size(); ++i) { */
        /* int pre  {(int)precut_cnt[i]}; */
        /* int mid  {(int)midcut_cnt[i]}; */
        /* int post {(int)cnt_pthatbins[i]}; */
        /* cout << Form(" %8i %8i %8i   %8i(%4.2f)   %8i(%4.2f)", */
        /*     pre, mid, post, (pre-mid), float(pre-mid)/(pre), */
        /*                     (mid-post),float(mid-post)/mid) << endl; */
        /* dat.log << Form(" %-8i ->%7i->  %-7i  %-7i  %8.2f", */
                /* pre, mid, post, (post-pre),(1.-float(post)/pre)) << endl; */
        /* cout << Form(" %-8i ->%7i->  %-7i  %-7i  %8.2f", */
                /* pre, mid, post, (post-pre),(1.-float(post)/pre)) << endl; */
    /* } */
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
