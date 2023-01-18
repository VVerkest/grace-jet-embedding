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


/*
 * Use to draw jet events.
 * 1. Make TH2D histogram for 
 *      Pythia dijet
 *      Non-Pythia dijet
 *      Together
 * 2. Add ghosts manually for all hisogram center of bins
 * 3. Cluster the three sets of jets
 * 4. Draw out to out.root
 */

using namespace std;
void jet_draw(events& dat, string _options) {
    // Start of 
    cout << " Running fn \"jet_draw\"" << endl;
    istringstream options ( _options );
    int n_options = 0;
    string arg;
    // options >> arg;
    vector<string> v_algo {"antikt","kt","cambridge"};
    vector<string> v_system {"emb","pAu","both"};
    vector<string> v_JTC {"jet","track","cells"};

    ioBinVec phibins {{ 0., 0.,  180, 6.283185307 }};
    ioBinVec etabins {{ -1, -1,  30, 1. }};

    map<string, oiJetMaker> jet_makers;
    map<string, TH2D>       hgrams;
    for (auto algo : v_algo)
    for (auto system : v_system) {
        /* cout << "algo: " << algo << endl; */
        jet_makers[Form("%s_%s", algo.c_str(), system.c_str())] 
            = oiJetMaker({{"jet_def",algo.c_str(),"calc_areas",1,"ghost_R",0.0005,
                    "jet_R",0.4,"min_jet_pt",0.,"jetrap",2.}});
        for (auto JTC: v_JTC) {
            string name = Form("%s_%s_%s", algo.c_str(), system.c_str(), JTC.c_str());
            /* cout << " name: " << name << endl; */
            hgrams[name] = {name.c_str(), name.c_str(), etabins, etabins, phibins, phibins };
            /* cout << "again: " << hgrams[name].GetName() << endl; */
        }
    }

    TH2D& hg = hgrams["CA_emb_jet"];

    while (dat.next()) {
        if (dat.size_mcTr() < 7) continue;

        map<int,int> emb_to_track;
        set<int> is_emb_track; 

        int i_emb{0};
        for (auto& track : dat.iter_mcTr) {
            is_emb_track.insert(track.id);

            auto tr = dat.get_track(track.id);
            if (track.id > -1) {
                emb_to_track[i_emb] = track.id;
            /*     cout << Form( */
            /*             "* eta(%6.2f,%6.2f,%6.2f) phi(%6.2f,%6.2f,%6.2f) pt(%6.2f,%6.2f,%6.2f)", */
            /*             track.eta, tr->eta, track.eta-tr->eta, */
            /*             track.phi, tr->phi, track.phi-tr->phi, */
            /*             track.pt,  tr->pt,  track.pt-tr->pt) << endl; */
            }
            ++i_emb;
        }

        /* for (auto& key : emb_to_track) { */
        /*     auto emb = dat.get_mcTr(key.first); */
        /*     auto tr  = dat.get_track(key.second); */
        /*     cout << Form( */
        /*             "eta(%6.2f,%6.2f,%6.2f) phi(%6.2f,%6.2f,%6.2f) pt(%6.2f,%6.2f,%6.2f)", */
        /*             emb->eta, tr->eta, emb->eta-tr->eta, */
        /*             emb->phi, tr->phi, emb->phi-tr->phi, */
        /*             emb->pt,  tr->pt, emb->pt-tr->pt) << endl; */
        /* } */
        // count pAu tracks:
        int n_cnt{0};
        int k{0};
        for (auto& tr : dat.iter_track) if (is_emb_track.count(k++) == 0) n_cnt++;

        /* if (dat.size_track() < 40) continue; */

        cout << Form(" emb: %2i  both: %2i  pAu: %2i   pT-lead: %4.2f",
            dat.size_mcTr(), dat.size_track(), n_cnt, dat.get_track(0)->pt) << endl;
        /* continue; */
        /* cout << " MC tracks: " << dat.size_mcTr() << " tracks: " << dat.size_track() << endl; */
        /* continue; */
        /* if (dat.size_mcTr() != 39 || dat.size_track() != 57) continue; */
        /* if (dat.size_mcTr() != 23 || dat.size_track() != 39) continue; */
        /* if (dat.size_mcTr() != 32 || dat.size_track() != 50 || (int)(dat.get_track(0)->pt)!=9) continue; */
        /* if (dat.size_mcTr() != 21 || dat.size_track() != 32 || (int)(dat.get_track(0)->pt)!=9) continue; */
        if (dat.size_mcTr() != 25 || dat.size_track() != 34 || (int)(dat.get_track(0)->pt)!=8) continue;
        /* if (dat.size_mcTr() != 13 || dat.size_track() != 27 || (int)(dat.get_track(0)->pt)!=11) continue; */


        

        // fill jet makers and cluster the jets
        /* set<int> is_emb_track; */ 
        for (auto track : dat.iter_mcTr) {
            if (track.id < 0) continue;
            /* is_emb_track.insert(track.id); */
            /* cout << " is_emb_track: " << track.id << " " << track.eta << " " << track.phi << " " << track.pt <<  endl; */
            for (auto algo : v_algo) {
                jet_makers[Form("%s_emb",algo.c_str())]
                    .add_particle(track.pt, track.eta, track.phi);
                /* cout << " tracks: " << jet_makers[Form("%s_emb",algo.c_str())].particles.size() << endl; */
            }
        }

        int i_track{0};
        int _cnt {0};
        for (auto track : dat.iter_track) {
            /* cout << " adding particle: " << track.pt << " " << track.eta << " " << track.phi << endl; */
            for (auto algo : v_algo) {
                jet_makers[Form("%s_both",algo.c_str())] 
                .add_particle(track.pt, track.eta, track.phi);
                /* cout << " 00000 " << jet_makers[Form("%s_both",algo.c_str())].particles.size() << endl; */
            }
            if (is_emb_track.count(i_track)==0) {
                for (auto algo : v_algo) {
                    jet_makers[Form("%s_pAu",algo.c_str())]
                    .add_particle(track.pt, track.eta, track.phi);
                    /* cout << " 11111 " << jet_makers[Form("%s_pAu",algo.c_str())].particles.size() << "  - > " << i_track << endl; */
                }
            }
            ++i_track;
            /* cout << " ITRACK " << i_track << endl; */
        } 

        // print how many
    for (auto algo : v_algo)
    for (auto system : v_system) 
        /* cout << " constituents: " << algo << " " << system << " " << jet_makers[ */
            /* Form("%s_%s",algo.c_str(),system.c_str())].particles.size() << endl; */

        /* cout << " YES " << endl; */
        // print out the number of particles of each type -- scoping work
        /* return; */


        bool first {true};
        for (auto algo : v_algo)
        for (auto system : v_system) {
            auto& maker = jet_makers[Form("%s_%s", algo.c_str(), system.c_str())];
            cout << "maker with " << maker.particles.size() << endl;
            maker.cluster_jets();

/* if (system == "emb") cout << " SIZE0: " << maker.particles.size() << endl; */
            /* cout << " njets-> " << maker.cluster_jets() << endl; */
/* if (system == "emb") cout << " SIZE1: " << maker.particles.size() << endl; */
/* if (system == "emb") cout << " jets: " << maker.pseudojets.size() << endl; */

            // fill the histograms
            auto& hgram_jet   = hgrams[Form("%s_%s_jet",   algo.c_str(), system.c_str())];
            auto& hgram_track = hgrams[Form("%s_%s_track", algo.c_str(), system.c_str())];
            auto& hgram_cells = hgrams[Form("%s_%s_cells", algo.c_str(), system.c_str())];
            for (auto& jet : maker.pseudojets) {
                double jet_pt = jet.is_pure_ghost() ? 0. : jet.perp();
                /* cout << " */ 
                if (!jet.is_pure_ghost()) hgram_jet.Fill(jet.eta(), jet.phi(), jet_pt);
                for (auto& track : jet.constituents()) {
                    int index = hgram_cells.FindBin(track.eta(),track.phi());
                    hgram_cells.SetBinContent(index, jet_pt);
                    if (!track.is_pure_ghost()) {
                         hgram_track.Fill(track.eta(), track.phi(), track.perp());
                    }
                }
            }
        }
        /* for (auto track : dat.iter_mcTr) { */
        /*     cout << " id: " << track.id << endl; */
        /*     if (track.id > 0) */ 
        /*     cout << " pt: " << track.pt << " phi: " << track.phi << " id: " << track.id */ 
        /*          << " pt-tracks: " << dat.get_track(track.id)->pt << " " << dat.get_track(track.id)->phi << endl; */
        /* } */
        /* break; */
        //break;
    }

    for (auto JTC :  v_JTC)
    for (auto algo : v_algo) 
    for (auto system : v_system) 
        hgrams[Form("%s_%s_%s", algo.c_str(), system.c_str(),JTC.c_str())].Write();


    ioPads pads {{{0.05,0.15, 0.42},{0.42,0.69},{0.69,0.69,.96,1.}},
                 {{0.,0.10, 0.36},{0.36,0.62},{0.62,0.62,.88,1.}}, 400, 800};

    gStyle->SetPalette(kCMYK);
    ioOptMap padopt {{"Title","","yAxisRangeLo",0.1,"yAxisRangeHi",2*TMath::Pi(),
                                  "zAxisRangeLo",0.2,"zAxisRangeHi",17.,
                                  "xAxisRangeLo",-0.80,"xAxisRangeHi",0.80,
                                    "xAxisNdivisions",5,
    "yAxisTitle","#phi","xAxisTitle","#eta","yAxisTitleOffset",3.75,
    "xAxisTitleOffset",2.5}};
    pads(2)->SetLogz();
    io_fmt(&hgrams["antikt_emb_cells"],padopt+ioOptMap{{"yAxisTitle","#phi(anti-k_{T})"}})->Draw("colz");
    pads(5)->SetLogz();
    io_fmt(&hgrams["antikt_pAu_cells"],padopt)->Draw("colz");
    pads(8)->SetLogz();
    io_fmt(&hgrams["antikt_both_cells"],padopt)->Draw("colz");

    pads(1)->SetLogz();
    io_fmt(&hgrams["kt_emb_cells"],padopt+ioOptMap{{"yAxisTitle","#phi(k_{T})"}})->Draw("colz");
    pads(4)->SetLogz();
    io_fmt(&hgrams["kt_pAu_cells"],padopt)->Draw("colz");
    pads(7)->SetLogz();
    io_fmt(&hgrams["kt_both_cells"],padopt)->Draw("colz");

    pads(0)->SetLogz();
    io_fmt(&hgrams["cambridge_emb_cells"],padopt+ioOptMap{{"yAxisTitle","#phi(Cambridge/Achen)",
            "xAxisTitle","#eta(#it{pp}_{embed})"}})->Draw("colz");
    pads(3)->SetLogz();
    io_fmt(&hgrams["cambridge_pAu_cells"],padopt+ioOptMap{{"xAxisTitle","#eta(#it{p}+Au)"}})->Draw("colz");
    pads(6)->SetLogz();
    io_fmt(&hgrams["cambridge_both_cells"],padopt+ioOptMap{{"xAxisTitle","#eta(Both)"}})->Draw("colz");

    pads.stamp("`date` `pwd`/$0");
    pads.canvas->SaveAs("pads.pdf");


    /* dat.log  << " Done running events" << endl; */
    /* cout     << " Done running events" << endl; */
}
