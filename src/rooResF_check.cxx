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
void rooResF_check(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"rooResF_check\"" << endl;
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
    msg.slurp_file("src/rooResF_check.cxx");
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

    // Jet Matchers for F/C  : True/Measured x Full/Charged
    ioXsec Xsec   { "/gpfs/loomis/home.grace/djs232"
        "/root_macros/io_lib/pAu2015_common/pAuXsection.txt" };
    ioJetMatcher F_matcher { "F", Xsec, edge_file, "prelim_Raghav_M", "prelim_Raghav_T" };
    ioJetMatcher C_matcher { "C", Xsec, edge_file, "prelim_Raghav_M", "prelim_Raghav_T" };

    ioBinVec Tbins {edge_file,"prelim_Raghav_T"};
    ioBinVec Mbins {edge_file,"prelim_Raghav_M"};
    TH1D hg_T {"T","test truth;pt;N",Tbins,Tbins};
    TH1D hg_M {"M","test measure;pt;N",Mbins,Mbins};
    TH2D hg_MT {"MT",";pt M;pt T",Mbins,Mbins,Tbins,Tbins};
    TH1D A_T {"A_T","test truth;pt;N",Tbins,Tbins};
    TH1D A_M {"A_M","test measure;pt;N",Mbins,Mbins};
    TH2D A_MT {"A_MT",";pt M;pt T",Mbins,Mbins,Tbins,Tbins};
    TH1D B_T {"B_T","test truth;pt;N",Tbins,Tbins};
    TH1D B_M {"B_M","test measure;pt;N",Mbins,Mbins};
    TH2D B_MT {"B_MT",";pt M;pt T",Mbins,Mbins,Tbins,Tbins};

    /* ioBinVeC test01 {edge_file,"01_M"}; */
    /* cout << " test_01_M " << test01 << endl; */

    /* ioBinVeC test01 {edge_file,"01_T"}; */
    /* cout << " test_01_T " << test01 << endl; */
    /* return; */

    int z {0};
    TH1D ZDCx  { "zdcX",";ZDCx;N_{events}",    100,0.,35000.};
    bool is_A = true;
    while (dat.next()) {
        precut_cnt[dat.pthat_bin] += 1;
        /* cout << " dat " << z++ << endl; */
        if (cut30_event(dat)) continue;
        ZDCx.Fill(dat.mu_event->ZDCx,Xsec.Xsec(dat.pthat_bin));

        // make full measured jets
        M_F_jets.reset();
        for (auto track : dat.iter_track) {
            M_F_jets.add_particle(track.pt, track.eta, track.phi);
        }
        while (tower.next()) {
            M_F_jets.add_particle(tower.Et, tower()->eta, tower()->phi);
        }
        M_F_jets.cluster_jets();

        // cut events with cut in full measured jet limits
        double limit { M_jetpt_limits[dat.pthat_bin] };
        bool skip_event{false};
        for (auto jet : M_F_jets.pseudojets) {
            if (jet.perp() >= limit) {
                skip_event=true;
                break;
            }
        }
        if (skip_event) continue;
        is_A = !is_A;

        // do custom check build
        double W { Xsec.Xsec(dat.pthat_bin) };
        auto matched = vector<bool>(M_F_jets.pseudojets.size(),false);
        for (auto MC : dat.iter_mc_Fjet) {
            int i{-1};
            for (auto rec : M_F_jets.pseudojets) {
                ++i;

                hg_T.Fill(MC.pt,W);
                if (is_A) A_T.Fill(MC.pt,W);
                else      B_T.Fill(MC.pt,W);

                if (matched[i]) continue;
                if (io_R2(MC.eta, MC.phi, rec.eta(), rec.phi())<0.16) {
                    hg_MT.Fill(rec.perp(), MC.pt, W);
                    if (is_A) A_MT.Fill(rec.perp(), MC.pt, W);
                    else      B_MT.Fill(rec.perp(), MC.pt, W);
                    matched[i] = true;
                    continue;
                }
            }
        }
        for (auto rec : M_F_jets.pseudojets) {
            hg_M.Fill(rec.perp(),W);
            if (is_A) A_M.Fill(rec.perp(),W);
            else      B_M.Fill(rec.perp(),W);
        }
        
        // see how close to the weighting things stay...
        cnt_pthatbins[dat.pthat_bin] += 1;

        // make the T-full jets
        T_F_jets.reset();
        for (auto track : dat.iter_mcTr){
            T_F_jets.add_particle(track.pt, track.eta, track.phi);
        }
        for (auto part : dat.iter_mcNeut){
            T_F_jets.add_particle(part.pt, part.eta, part.phi);
        }
        T_F_jets.cluster_jets();

        // fill in the jetmatchers
        for (auto& jet : M_F_jets.pseudojets) {
            F_matcher.addjet_reco(jet.eta(), jet.phi(), jet.perp());
        }
        for (auto& jet : T_F_jets.pseudojets) {
            F_matcher.addjet_MC(jet.eta(), jet.phi(), jet.perp());
        }
        F_matcher.do_matching(dat.pthat_bin); 

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

    hg_T.Write();
    hg_M.Write();
    hg_MT.Write();

    A_T.Write();
    A_M.Write();
    A_MT.Write();

    B_T.Write();
    B_M.Write();
    B_MT.Write();


    dat.log << Form(" %-8s ->  %-7s  %-7s  %-8s",
            "pre-cut","post-cut","delta","ratio-cut") << endl;
    for (int i{0}; i<cnt_pthatbins.size(); ++i) {
        int pre  {(int)precut_cnt[i]};
        int post {(int)cnt_pthatbins[i]};
        dat.log << Form(" %-8i ->  %-7i  %-7i  %8.2f",
                pre, post, (post-pre),(1.-float(post)/pre)) << endl;
    }
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
