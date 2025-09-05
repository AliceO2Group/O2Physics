/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* 
Study for the antineutron CEX interaction background
author: Fabiola Lugo 
*/

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Logger.h"
#include "Framework/AnalysisDataModel.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "DetectorsBase/Propagator.h"
#include "TMCProcess.h"
#include <optional>
//ROOT
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include <cmath>

using namespace o2;
using namespace o2::framework;

struct AntineutronTask {
    HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
    OutputObj<TTree> cexTree{"cex_tree"};   // CEX
    OutputObj<TTree> bgTree{"bg_tree"};     // background
        // === Branch buffers ===
        struct PairBuffers {
          // MC (pair)
          double mc_pair_p, mc_pair_pt, mc_pair_pz;
          double mc_dplane;
          double mc_angle_deg;
          double mc_vtx_x, mc_vtx_y, mc_vtx_z;
          // TRK (pair, fitter)
          double trk_pair_p, trk_pair_pt, trk_pair_pz;
          double trk_angle_deg;
          double trk_vtxfit_dcaPair;
          double trk_vtxfit_R;
          double trk_vtxfit_distToPV;
          double trk_vtxfit_secVtx_x, trk_vtxfit_secVtx_y, trk_vtxfit_secVtx_z;
          // fit quality
          double vtxfit_chi2;
          int    vtxfit_status;
          int    nCand;
          double vtxfit_dX, vtxfit_dY, vtxfit_dZ; // (fit − MC)
          double vtxfit_d3D;
          // proton trk
          double p_trk_p, p_trk_px, p_trk_py, p_trk_pz, p_trk_eta, p_trk_tpcSignal;
          int    p_trk_nClsITS;
          // antiproton trk
          double antip_trk_p, antip_trk_px, antip_trk_py, antip_trk_pz, antip_trk_eta, antip_trk_tpcSignal;
          int    antip_trk_nClsITS;
          // Meta
          int    mother_pdg, colId, p_id, antip_id;
        } buf;
    
    void init(InitContext const&)
    {
        // Primary vertex
        histos.add("hVx", "Primary vertex X;X (cm);Entries", kTH1F, {{100, -5., 5.}});
        histos.add("hVy", "Primary vertex Y;Y (cm);Entries", kTH1F, {{100, -5., 5.}});
        histos.add("hVz", "Primary vertex Z;Z (cm);Entries", kTH1F, {{200, -20., 20.}});

        // Primary antineutrons
        histos.add("antin_p", "Total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("antin_px", "p_{x};p_{x} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("antin_py", "p_{y};p_{y} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("antin_pz", "p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("antin_eta","Pseudorapidity;#eta;Entries", kTH1F, {{100, -10., 10.}});
        histos.add("antin_p_ITScuts", "Momentum with ITS cuts;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

        // Primary neutrons
        histos.add("n_p",  "Total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("n_px", "p_{x};p_{x} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("n_py", "p_{y};p_{y} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("n_pz", "p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("n_eta","Pseudorapidity;#eta;Entries", kTH1F, {{100, -10., 10.}});
        histos.add("n_p_ITScuts", "Momentum with ITS cuts;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

        // Primary antiprotons
        histos.add("antip_p",  "Total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("antip_px", "p_{x};p_{x} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("antip_py", "p_{y};p_{y} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("antip_pz", "p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("antip_eta","Pseudorapidity;#eta;Entries", kTH1F, {{100, -10., 10.}});
        histos.add("antip_p_ITScuts", "Momentum with ITS cuts;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

        // Primary protons
        histos.add("p_p", "Total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("p_px", "p_{x};p_{x} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("p_py", "p_{y};p_{y} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("p_pz", "p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("p_eta","Pseudorapidity;#eta;Entries", kTH1F, {{100, -10., 10.}});
        histos.add("p_p_ITScuts", "Momentum with ITS cuts;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

        // test (MC)
        histos.add("antip_test", "Secondary antiprotons;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("primmom_test", "Secondary particles;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

        // CEX pair from antineutron (MC)
        histos.add("cex_pairmc_p", "CEX pair total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("cex_pairmc_pt", "CEX pair p_{T};p_{T} (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("cex_pairmc_pz", "CEX pair p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("cex_pairmc_dplane","CEX pair d_{plane};d_{plane} (cm);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("cex_pairmc_angle", "Pair opening angle;Angle (°);Entries", kTH1F, {{180, 0., 180.}});
        histos.add("cex_pairmc_vtx", "MC CEX pair vertex;X (cm);Y (cm)", kTH2F, {{100, -50., 50.}, {100, -50., 50.}});
        histos.add("cex_pairmc_vtxz", "MC secondary vertex Z;Z (cm);Entries", kTH1F, {{200, -60., 60.}});
        histos.add("cex_pairmc_pITScuts","CEX pair momentum (ITS cuts);|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

        // CEX pair normalized to antineutron (MC)
        histos.add("cexn_pairmc_p", "Pair p / antineutron p;p/p_{#bar{n}};Entries", kTH1F, {{100, 0., 2.}});
        histos.add("cexn_pairmc_pt", "Pair p_{T} / antineutron p_{T};p_{T}/p_{T,#bar{n}};Entries", kTH1F, {{100, 0., 2.}});
        histos.add("cexn_pairmc_pz", "Pair p_{z} / antineutron p_{z};p_{z}/p_{z,#bar{n}};Entries", kTH1F, {{100, -2., 2.}});

        // BG pair (not from antineutron) (MC)
        histos.add("cexbg_pairmc_p", "Background pair momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("cexbg_pairmc_pt", "Background pair p_{T};p_{T} (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("cexbg_pairmc_pz", "Background pair p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
        histos.add("cexbg_pairmc_dplane","Background d_{plane};d_{plane} (cm);Entries", kTH1F, {{100, 0., 10.}});
        histos.add("cexbg_pairmc_angle", "Background opening angle;Angle (°);Entries", kTH1F, {{180, 0., 180.}});
        histos.add("cexbg_pairmc_vtx", "Background pair vertex;X (cm);Y (cm)", kTH2F, {{100, -50., 50.}, {100, -50., 50.}});
        histos.add("cexbg_pairmc_vtxz", "Background secondary vertex Z;Z (cm);Entries", kTH1F, {{200, -60., 60.}});
        histos.add("cexbg_pairmc_pITScuts","Background momentum (ITS cuts);|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

        // CEX pair from antineutron (TRK)
        histos.add("cex_pairtrk_angle", "Pair opening angle (tracks);Angle (°);Entries", kTH1F, {{180, 0., 180.}});
        histos.add("cex_pairtrk_p", "Pair momentum (tracks);|p| (GeV/c);Entries", kTH1F, {{120, 0., 12.}});
        histos.add("cex_pairtrk_pt", "Pair p_{T} (tracks);p_{T} (GeV/c);Entries", kTH1F, {{120, 0., 12.}});
        histos.add("cex_pairtrk_pz", "Pair p_{z} (tracks);p_{z} (GeV/c);Entries", kTH1F, {{120, -12., 12.}});
        histos.add("cex_pairtrk_vtxfit_dcaPair", "DCA between tracks at PCA;DCA (cm);Entries", kTH1F, {{200, 0., 10.}});
        histos.add("cex_pairtrk_vtxfit_R", "Secondary-vertex radius (PCA);R (cm);Entries", kTH1F, {{200, 0., 60.}});
        histos.add("cex_pairtrk_vtxfit_distToPV", "Distance from secondary vertex to PV;dist (cm);Entries", kTH1F, {{240, 0., 120.}});
        histos.add("cex_pairtrk_vtxfit_secVtxXY", "Secondary vertex (PCA);X (cm);Y (cm)", kTH2F, {{200, -60., 60.}, {200, -60., 60.}});
        histos.add("cex_pairtrk_vtxfit_secVtxZ", "Secondary vertex Z (PCA);Z (cm);Entries", kTH1F, {{240, -60., 60.}});

        // BG pair (not from antineutron) (TRK)
        histos.add("cexbg_pairtrk_angle", "Background opening angle (tracks);Angle (°);Entries", kTH1F, {{180, 0., 180.}});
        histos.add("cexbg_pairtrk_p", "Pair momentum (tracks);|p| (GeV/c);Entries", kTH1F, {{120, 0., 12.}});
        histos.add("cexbg_pairtrk_pt", "Pair p_{T} (tracks);p_{T} (GeV/c);Entries", kTH1F, {{120, 0., 12.}});
        histos.add("cexbg_pairtrk_pz", "Pair p_{z} (tracks);p_{z} (GeV/c);Entries", kTH1F, {{120, -12., 12.}});
        histos.add("cexbg_pairtrk_vtxfit_dcaPair", "DCA between tracks at PCA;DCA (cm);Entries", kTH1F, {{200, 0., 10.}});
        histos.add("cexbg_pairtrk_vtxfit_R", "Secondary-vertex radius (PCA);R (cm);Entries", kTH1F, {{200, 0., 60.}});
        histos.add("cexbg_pairtrk_vtxfit_distToPV", "Distance from secondary vertex to PV;dist (cm);Entries", kTH1F, {{240, 0., 120.}});
        histos.add("cexbg_pairtrk_vtxfit_secVtxXY", "Secondary vertex (PCA);X (cm);Y (cm)", kTH2F, {{200, -60., 60.}, {200, -60., 60.}});
        histos.add("cexbg_pairtrk_vtxfit_secVtxZ", "Secondary vertex Z (PCA);Z (cm);Entries", kTH1F, {{240, -60., 60.}});

        // Vertex fit (DCAFitter2 / PCA)
        histos.add("vtxfit_chi2", "DCAFitter2 #chi^{2};#chi^{2};Entries", kTH1F, {{200, 0., 100.}});
        histos.add("vtxfit_status", "Fit status (0=OK);code;Entries", kTH1I, {{10, 0., 10.}});
        histos.add("vtxfit_mc_dX", "SV residual X (fit - MC);#Delta X (cm);Entries", kTH1F, {{400, -20., 20.}});
        histos.add("vtxfit_mc_dY", "SV residual Y (fit - MC);#Delta Y (cm);Entries", kTH1F, {{400, -20., 20.}});
        histos.add("vtxfit_mc_dZ", "SV residual Z (fit - MC);#Delta Z (cm);Entries", kTH1F, {{400, -20., 20.}});
        histos.add("vtxfit_mc_d3D", "SV distance |fit - MC|;#Delta r (cm);Entries",   kTH1F, {{300, 0., 30.}});
        
        //Trees
        cexTree.setObject(new TTree("cex_pairs", "CEX (antineutron daughters)"));
        bgTree.setObject(new TTree("bg_pairs",  "Background pairs"));
        // Helper
        auto defineBranches = [&](o2::framework::OutputObj<TTree>& t) {
          // MC
          t->Branch("mc_pair_p",    &buf.mc_pair_p);
          t->Branch("mc_pair_pt",   &buf.mc_pair_pt);
          t->Branch("mc_pair_pz",   &buf.mc_pair_pz);
          t->Branch("mc_dplane",    &buf.mc_dplane);
          t->Branch("mc_angle_deg", &buf.mc_angle_deg);
          t->Branch("mc_vtx_x",     &buf.mc_vtx_x);
          t->Branch("mc_vtx_y",     &buf.mc_vtx_y);
          t->Branch("mc_vtx_z",     &buf.mc_vtx_z);
          // TRK (par)
          t->Branch("trk_pair_p",   &buf.trk_pair_p);
          t->Branch("trk_pair_pt",  &buf.trk_pair_pt);
          t->Branch("trk_pair_pz",  &buf.trk_pair_pz);
          t->Branch("trk_angle_deg",&buf.trk_angle_deg);
          t->Branch("trk_vtxfit_dcaPair",      &buf.trk_vtxfit_dcaPair);
          t->Branch("trk_vtxfit_R",            &buf.trk_vtxfit_R);
          t->Branch("trk_vtxfit_distToPV",     &buf.trk_vtxfit_distToPV);
          t->Branch("trk_vtxfit_secVtx_x",     &buf.trk_vtxfit_secVtx_x);
          t->Branch("trk_vtxfit_secVtx_y",     &buf.trk_vtxfit_secVtx_y);
          t->Branch("trk_vtxfit_secVtx_z",     &buf.trk_vtxfit_secVtx_z);
          // Fit/residuales
          t->Branch("vtxfit_chi2",   &buf.vtxfit_chi2);
          t->Branch("vtxfit_status", &buf.vtxfit_status);
          t->Branch("nCand",         &buf.nCand);
          t->Branch("vtxfit_dX",     &buf.vtxfit_dX);
          t->Branch("vtxfit_dY",     &buf.vtxfit_dY);
          t->Branch("vtxfit_dZ",     &buf.vtxfit_dZ);
          t->Branch("vtxfit_d3D",    &buf.vtxfit_d3D);
          // Tracks individuales
          t->Branch("p_trk_p",         &buf.p_trk_p);
          t->Branch("p_trk_px",        &buf.p_trk_px);
          t->Branch("p_trk_py",        &buf.p_trk_py);
          t->Branch("p_trk_pz",        &buf.p_trk_pz);
          t->Branch("p_trk_eta",       &buf.p_trk_eta);
          t->Branch("p_trk_tpcSignal", &buf.p_trk_tpcSignal);
          t->Branch("p_trk_nClsITS",   &buf.p_trk_nClsITS);
          t->Branch("antip_trk_p",         &buf.antip_trk_p);
          t->Branch("antip_trk_px",        &buf.antip_trk_px);
          t->Branch("antip_trk_py",        &buf.antip_trk_py);
          t->Branch("antip_trk_pz",        &buf.antip_trk_pz);
          t->Branch("antip_trk_eta",       &buf.antip_trk_eta);
          t->Branch("antip_trk_tpcSignal", &buf.antip_trk_tpcSignal);
          t->Branch("antip_trk_nClsITS",   &buf.antip_trk_nClsITS);
          // Meta
          t->Branch("mother_pdg", &buf.mother_pdg);
          t->Branch("colId",      &buf.colId);
          t->Branch("p_id",       &buf.p_id);
          t->Branch("antip_id",   &buf.antip_id);
        };

        defineBranches(cexTree);
        defineBranches(bgTree);
    }
    //Check available tables in the AOD, specifically TracksIU, TracksCovIU
    using TracksWCovMc = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels, o2::aod::TracksCovIU>;
    
    static o2::track::TrackParCov makeTPCovFromAOD(const TracksWCovMc::iterator& tr)
    {
      using o2::track::TrackParCov;
      TrackParCov tpcov;
      // Local state: x, alpha, y, z, snp, tgl, q/pt
      float par[5] = {tr.y(), tr.z(), tr.snp(), tr.tgl(), tr.signed1Pt()};
      tpcov.set(tr.x(), tr.alpha(), par);

      // Covariance matrix (15 terms) in O2 order
      std::array<float, 15> cov = {
          tr.cYY(),       tr.cZY(),       tr.cZZ(),
          tr.cSnpY(),     tr.cSnpZ(),     tr.cSnpSnp(),
          tr.cTglY(),     tr.cTglZ(),     tr.cTglSnp(),
          tr.cTglTgl(),   tr.c1PtY(),     tr.c1PtZ(),
          tr.c1PtSnp(),   tr.c1PtTgl(),   tr.c1Pt21Pt2()
      };
      tpcov.setCov(cov);
      return tpcov;
    }
    
    void process(aod::McCollisions const& cols, aod::McParticles const& particles, TracksWCovMc const& tracks) {
        double pvtx_x = 0;
        double pvtx_y = 0;
        double pvtx_z = 0;
        for (auto const& col : cols) {
            const auto colId= col.globalIndex();
            if (std::isfinite(col.posX()) && std::isfinite(col.posY()) && std::isfinite(col.posZ())) {
                pvtx_x = col.posX();
                pvtx_y = col.posY();
                pvtx_z = col.posZ();
                histos.fill(HIST("hVx"), pvtx_x);
                histos.fill(HIST("hVy"), pvtx_y);
                histos.fill(HIST("hVz"), pvtx_z);
            }

            for (auto& particle : particles) {
                const auto procEnum = particle.getProcess();
                const bool isSecondaryFromMaterial = (!particle.producedByGenerator()) && (procEnum == kPHadronic || procEnum == kPHInhelastic);
                if (particle.mcCollisionId() != colId) continue;

                // Primary antineutrons
                if (particle.pdgCode() == -2112 && particle.isPhysicalPrimary()) {
                    histos.fill(HIST("antin_p"), particle.p());
                    histos.fill(HIST("antin_px"), particle.px());
                    histos.fill(HIST("antin_py"), particle.py());
                    histos.fill(HIST("antin_pz"), particle.pz());
                    histos.fill(HIST("antin_eta"), particle.eta());
                    if (TMath::Abs(particle.eta()) < 1.5 && TMath::Abs(particle.vz())< 5.3) histos.fill(HIST("antin_p_ITScuts"),particle.p());
                }
                // Primary neutrons
                if (particle.pdgCode() == 2112 && particle.isPhysicalPrimary()) {
                    histos.fill(HIST("n_p"), particle.p());
                    histos.fill(HIST("n_px"), particle.px());
                    histos.fill(HIST("n_py"), particle.py());
                    histos.fill(HIST("n_pz"), particle.pz());
                    histos.fill(HIST("n_eta"), particle.eta());
                    if (TMath::Abs(particle.eta()) < 1.5 && TMath::Abs(particle.vz())< 5.3) histos.fill(HIST("n_p_ITScuts"),particle.p());
                }
                // Primary antiprotons
                if (particle.pdgCode() == -2212 && particle.isPhysicalPrimary()) {
                    histos.fill(HIST("antip_p"), particle.p());
                    histos.fill(HIST("antip_px"), particle.px());
                    histos.fill(HIST("antip_py"), particle.py());
                    histos.fill(HIST("antip_pz"), particle.pz());
                    histos.fill(HIST("antip_eta"), particle.eta());
                    if (TMath::Abs(particle.eta()) < 1.5 && TMath::Abs(particle.vz())< 5.3) histos.fill(HIST("antip_p_ITScuts"),particle.p());
                }
                // Primary protons
                if (particle.pdgCode() == 2212 && particle.isPhysicalPrimary()) {
                    histos.fill(HIST("p_p"), particle.p());
                    histos.fill(HIST("p_px"), particle.px());
                    histos.fill(HIST("p_py"), particle.py());
                    histos.fill(HIST("p_pz"), particle.pz());
                    histos.fill(HIST("p_eta"), particle.eta());
                    if (TMath::Abs(particle.eta()) < 1.5 && TMath::Abs(particle.vz())< 5.3) histos.fill(HIST("p_p_ITScuts"),particle.p());
                }
                
                if (particle.pdgCode() == -2212 && isSecondaryFromMaterial && !particle.mothersIds().empty()) {
                    for (auto mother : particle.mothers_as<aod::McParticles>()) {
                        if ((mother.isPhysicalPrimary())) histos.fill(HIST("primmom_test"), particle.p());
                    }
                }
                
                //Seconday antiprotons from material
                if (particle.pdgCode() == -2212 && isSecondaryFromMaterial && !particle.mothersIds().empty()) {
                    histos.fill(HIST("antip_test"), particle.p());
                    //Primary mother
                    for (auto mother : particle.mothers_as<aod::McParticles>()) { 
                        if (!(mother.isPhysicalPrimary())) continue;
                        double mother_pt = mother.pt();
                        double mother_pz = mother.pz();
                        double mother_vz = mother.vz();
                        double mother_p = mother.p();
                        double mother_eta = mother.eta();
                        int mother_pdg = mother.pdgCode();
                        
                        double antip_vx = particle.vx();
                        double antip_vy = particle.vy();
                        double antip_vz = particle.vz();
                        double antip_p = particle.p();
                        double antip_px = particle.px();
                        double antip_py = particle.py();
                        double antip_pz = particle.pz();
                        double antip_e = particle.e();
                        double antip_eta = particle.eta();
                        int antip_id = particle.globalIndex();
                        
                        //Selection conditions
                        const double R = TMath::Sqrt(antip_vx*antip_vx + antip_vy*antip_vy);
                        //Config for ITS
                        //if(3.9<=R && R<=43.0 && TMath::Abs(antip_vz)<=48.9){
                        //Config for ITS2
                        if(2.2<=R && R<=39.0 && TMath::Abs(antip_vz)<=39.0){
                            if (TMath::Abs(mother_eta) < 1.5 && TMath::Abs(mother_vz) < 5.3)
                            {
                                //Pion minus veto
                                bool pion_m = false;
                                for (auto& particle1 : particles) {
                                    const auto proc1Enum = particle1.getProcess();
                                    const bool isSecondaryFromMaterial1 = (!particle1.producedByGenerator()) && (proc1Enum == kPHadronic || proc1Enum == kPHInhelastic);
                                    if (particle1.mcCollisionId() != colId) continue;
                                    if (particle1.pdgCode() != -211 || !isSecondaryFromMaterial1 || particle1.mothersIds().empty()) continue;
                                    bool hasPrimaryMother_pim = false;
                                    for (auto mother : particle1.mothers_as<aod::McParticles>()) {
                                        if (mother.isPhysicalPrimary()) {
                                            hasPrimaryMother_pim = true;
                                            break;
                                        }
                                    }
                                    if (!hasPrimaryMother_pim) continue;
                                    double pim_vx = particle1.vx();
                                    double pim_vy = particle1.vy();
                                    double pim_vz = particle1.vz();
                                    if(pim_vx == antip_vx && pim_vy == antip_vy && pim_vz == antip_vz){
                                        pion_m = true;
                                        break;
                                    }
                                }
                                //Pion plus veto
                                bool pion_p = false;
                                for (auto& particle2 : particles) {
                                    if (particle2.mcCollisionId() != colId) continue;
                                    const auto proc2Enum = particle2.getProcess();
                                    const bool isSecondaryFromMaterial2 = (!particle2.producedByGenerator()) && (proc2Enum == kPHadronic || proc2Enum == kPHInhelastic);
                                    if (particle2.pdgCode() != 211 || !isSecondaryFromMaterial2 || particle2.mothersIds().empty()) continue;
                                    bool hasPrimaryMother_pip = false;
                                    for (auto mother : particle2.mothers_as<aod::McParticles>()) {
                                        if (mother.isPhysicalPrimary()) {
                                            hasPrimaryMother_pip = true;
                                            break;
                                        }
                                    }
                                    if (!hasPrimaryMother_pip) continue;
                                    double pip_vx = particle2.vx();
                                    double pip_vy = particle2.vy();
                                    double pip_vz = particle2.vz();
                                    if(pip_vx == antip_vx && pip_vy == antip_vy && pip_vz == antip_vz){
                                        pion_p = true;
                                        break;
                                    }
                                }
                                if (pion_p == false && pion_m == false){
                                    // CEX selection
                                    double dplane = 10;
                                    double dplane_tmp = 0;
                                    double p = 0;
                                    double p_tmp = 0;
                                    double pcex_px = 0;
                                    double pcex_py = 0;
                                    double pcex_pz = 0;
                                    double pcex_p = 0;
                                    double e = 0;
                                    double e_tmp = 0;
                                    double pcex_e = 0;
                                    int k_plane = -1;
                                    int k_e = -1;
                                    int k_p = -1;
                                    //Secondary proton from material
                                    for (auto& particle3 : particles) {
                                        if (particle3.mcCollisionId() != colId) continue;
                                        const auto proc3Enum = particle3.getProcess();
                                        const bool isSecondaryFromMaterial3 = (!particle3.producedByGenerator()) && (proc3Enum == kPHadronic || proc3Enum == kPHInhelastic);
                                        if (particle3.pdgCode() != 2212 || !isSecondaryFromMaterial3 || particle3.mothersIds().empty()) continue;
                                        bool hasPrimaryMother_p = false;
                                        for (auto mother : particle3.mothers_as<aod::McParticles>()) {
                                            if (mother.isPhysicalPrimary()) {
                                                hasPrimaryMother_p = true;
                                                break;
                                            }
                                        }
                                        if (!hasPrimaryMother_p) continue;
                                        double p_vx = particle3.vx();
                                        double p_vy = particle3.vy();
                                        double p_vz = particle3.vz();
                                        double p_p = particle3.p();
                                        double p_px = particle3.px();
                                        double p_py = particle3.py();
                                        double p_pz = particle3.pz();
                                        double p_e = particle3.e();
                                        double p_eta = particle3.eta();
                                        
                                        if(p_vx == antip_vx && p_vy == antip_vy && p_vz == antip_vz){
                                            bool shareMother = false;
                                            {
                                              const auto& momsAp = particle.mothersIds();   //antiproton
                                              const auto& momsP  = particle3.mothersIds();  //proton
                                              for (auto ida : momsAp) {
                                                for (auto idp : momsP) {
                                                    if (ida == idp) {
                                                        shareMother = true;
                                                        break;
                                                    }
                                                }
                                                if (shareMother) break;
                                              }
                                            }
                                            if (!shareMother) continue;
                                            //dplane_tmp = (p_py*antip_pz - p_pz*antip_py)*(pvtx_x-antip_vx) + (p_pz*antip_px - p_px*antip_pz)*(pvtx_y-antip_vy) + (p_px*antip_py - p_py*antip_px)*(pvtx_z-antip_vz);
                                            double nx = (p_py*antip_pz - p_pz*antip_py);
                                            double ny = (p_pz*antip_px - p_px*antip_pz);
                                            double nz = (p_px*antip_py - p_py*antip_px);
                                            double rx = (pvtx_x - antip_vx);
                                            double ry = (pvtx_y - antip_vy);
                                            double rz = (pvtx_z - antip_vz);
                                            double denom = nx*nx + ny*ny + nz*nz;
                                            if (denom > 0.) {
                                              dplane_tmp = TMath::Abs(nx*rx + ny*ry + nz*rz) / TMath::Sqrt(denom);
                                            } else {
                                              dplane_tmp = 1e9;
                                            }
                                            if(TMath::Abs(dplane_tmp) < TMath::Abs(dplane))
                                            {
                                                k_plane = particle3.globalIndex();
                                                dplane = dplane_tmp;
                                            }
                                            
                                            e_tmp = antip_e + p_e;
                                            if(TMath::Abs(e_tmp) > TMath::Abs(e))
                                            {
                                                k_e = particle3.globalIndex();
                                                e = e_tmp;
                                                pcex_e = p_e;
                                            }
                                            
                                            p_tmp = TMath::Sqrt(TMath::Power((p_px+antip_px),2)+TMath::Power((p_py+antip_py),2)+TMath::Power((p_pz+antip_pz),2));
                                            if(TMath::Abs(p_tmp) > TMath::Abs(p))
                                            {
                                                k_p = particle3.globalIndex();
                                                p = p_tmp;
                                                pcex_p =  p_p;
                                                pcex_px =  p_px;
                                                pcex_py =  p_py;
                                                pcex_pz =  p_pz;
                                            }
                                            if(k_plane == k_e && k_plane == k_p && k_plane >= 0 ){
                                                int p_id = k_plane;
                                                TVector3 pVecProton = TVector3(pcex_px, pcex_py, pcex_pz);
                                                TVector3 pVecAntiproton = TVector3(antip_px, antip_py, antip_pz);
                                                TVector3 total_mc_pVec = pVecProton + pVecAntiproton;
                                                double cex_pairmc_p = total_mc_pVec.Mag();
                                                double cex_pairmc_pt = total_mc_pVec.Pt();
                                                double cex_pairmc_pz =  pcex_pz+antip_pz;
                                                double mcangleRad = pVecProton.Angle(pVecAntiproton);
                                                double mcangleDeg = mcangleRad * TMath::RadToDeg();
                                                
                                                //Antineutron mother
                                                if (mother_pdg == -2112) {
                                                    //CEX pair
                                                    histos.fill(HIST("cex_pairmc_p"), cex_pairmc_p);
                                                    histos.fill(HIST("cex_pairmc_pt"), cex_pairmc_pt);
                                                    histos.fill(HIST("cex_pairmc_pz"), cex_pairmc_pz);
                                                    histos.fill(HIST("cex_pairmc_dplane"), dplane);
                                                    histos.fill(HIST("cex_pairmc_angle"), mcangleDeg);
                                                    histos.fill(HIST("cex_pairmc_vtx"), antip_vx, antip_vy);
                                                    histos.fill(HIST("cex_pairmc_vtxz"), antip_vz);
                                                    if (TMath::Abs(mother_eta) < 0.9 && TMath::Abs(mother_vz)< 5.3)  histos.fill(HIST("cex_pairmc_pITScuts"),cex_pairmc_p);
                                                    //CEX pair normalized
                                                    if (mother_p  != 0) histos.fill(HIST("cexn_pairmc_p"),  cex_pairmc_p /mother_p);
                                                    if (mother_pt != 0) histos.fill(HIST("cexn_pairmc_pt"), cex_pairmc_pt/mother_pt);
                                                    if (mother_pz != 0) histos.fill(HIST("cexn_pairmc_pz"), cex_pairmc_pz/mother_pz);
                                                }
                                                //BG mother
                                                if (mother_pdg != -2112) {
                                                    //CEX pair
                                                    histos.fill(HIST("cexbg_pairmc_p"), cex_pairmc_p);
                                                    histos.fill(HIST("cexbg_pairmc_pt"), cex_pairmc_pt);
                                                    histos.fill(HIST("cexbg_pairmc_pz"), cex_pairmc_pz);
                                                    histos.fill(HIST("cexbg_pairmc_dplane"), dplane);
                                                    histos.fill(HIST("cexbg_pairmc_angle"), mcangleDeg);
                                                    histos.fill(HIST("cexbg_pairmc_vtx"), antip_vx, antip_vy);
                                                    histos.fill(HIST("cexbg_pairmc_vtxz"), antip_vz);
                                                    if (TMath::Abs(mother_eta) < 0.9 && TMath::Abs(mother_vz)< 5.3)  histos.fill(HIST("cexbg_pairmc_pITScuts"),cex_pairmc_p);
                                                }
                                                
                                                //Detector signal
                                                bool antip_layers = false;
                                                bool antip_hasTrack = false;
                                                double antip_trk_px = 0.;
                                                double antip_trk_py = 0.;
                                                double antip_trk_pz = 0.;
                                                double antip_trk_p = 0.;
                                                double antip_trk_eta = 0.;
                                                double antip_trk_tpcSignal = 0;
                                                //int antip_trk_nClsTPC = 0;
                                                int antip_trk_nClsITS = 0;
                                                
                                                bool p_layers = false;
                                                bool p_hasTrack = false;
                                                double p_trk_px = 0.;
                                                double p_trk_py = 0.;
                                                double p_trk_pz = 0.;
                                                double p_trk_p = 0.;
                                                double p_trk_eta = 0.;
                                                double p_trk_tpcSignal = 0;
                                                //int p_trk_nClsTPC = 0;
                                                int p_trk_nClsITS = 0;
                                                
                                                for (auto& track : tracks) {
                                                    if (!track.has_mcParticle()) continue;
                                                    const auto& mc = track.mcParticle();
                                                    uint8_t itsMap = track.itsClusterMap();
                                                    //Config for ITS1
                                                    /*bool hitSPD = (itsMap & 0x3) != 0; // bits 0 (SPD L1) & 1 (SPD L2)
                                                        bool hitSDD  = (itsMap & 0xC) != 0; // bits 2–3
                                                        bool hitSSD  = (itsMap & 0x30) != 0; // bits 4–5
                                                        bool layer_condition = (hitSDD || hitSSD) && !hitSPD;*/
                                                    //Config for ITS2
                                                    bool hitL0 = (itsMap & (1u<<0)) != 0;
                                                    bool hitL1 = (itsMap & (1u<<1)) != 0;
                                                    bool hitL2 = (itsMap & (1u<<2)) != 0;
                                                    bool hitL3 = (itsMap & (1u<<3)) != 0;
                                                    bool hitL4 = (itsMap & (1u<<4)) != 0;
                                                    bool hitL5 = (itsMap & (1u<<5)) != 0;
                                                    bool hitL6 = (itsMap & (1u<<6)) != 0;
                                                    bool hitIB = (hitL0 || hitL1 || hitL2);
                                                    bool hitOuter = (hitL3 || hitL4 || hitL5 || hitL6);
                                                    int nITS = track.itsNCls();
                                                    bool layer_condition = (!hitIB) && hitOuter && (nITS >= 2);
                                                    
                                                    if (mc.globalIndex() == antip_id) {
                                                        antip_trk_p = track.p();
                                                        antip_trk_px = track.px();
                                                        antip_trk_py = track.py();
                                                        antip_trk_pz = track.pz();
                                                        antip_trk_eta = track.eta();
                                                        antip_trk_tpcSignal = track.tpcSignal();
                                                        //antip_trk_nClsTPC = track.tpcNCls();
                                                        antip_trk_nClsITS = track.itsNCls();
                                                        antip_hasTrack = true;
                                                        if (layer_condition) antip_layers = true;
                                                    }
                                                    else if (mc.globalIndex() == p_id) {
                                                        p_trk_p = track.p();
                                                        p_trk_px = track.px();
                                                        p_trk_py = track.py();
                                                        p_trk_pz = track.pz();
                                                        p_trk_eta = track.eta();
                                                        p_trk_tpcSignal = track.tpcSignal();
                                                        //p_trk_nClsTPC = track.tpcNCls();
                                                        p_trk_nClsITS = track.itsNCls();
                                                        p_hasTrack = true;
                                                        if (layer_condition) p_layers = true;
                                                    }
                                                }
                                                
                                                if (p_hasTrack == true && antip_hasTrack == true){
                                                    TVector3 pVecProton_trk(p_trk_px, p_trk_py, p_trk_pz);
                                                    TVector3 AntipVecProton_trk(antip_trk_px, antip_trk_py, antip_trk_pz);
                                                    TVector3 total_trk_pVec = pVecProton_trk + AntipVecProton_trk;
                                                    double trkangleRad = AntipVecProton_trk.Angle(pVecProton_trk);
                                                    double trkangleDeg = trkangleRad * TMath::RadToDeg();
                                                    if(mother_pdg == -2112) histos.fill(HIST("cex_pairtrk_angle"), trkangleDeg);
                                                    if(mother_pdg != -2112) histos.fill(HIST("cexbg_pairtrk_angle"), trkangleDeg);
                                                    
                                                    // ==== Secondary vertex via central DCA vertexer (DCAFitter2) ====
                                                    using o2::vertexing::DCAFitter2;
                                                    constexpr float BZ_TESLA = 0.5f;

                                                    DCAFitter2 fitter(/*bz=*/BZ_TESLA, /*useAbsDCA=*/true, /*propagateToPCA=*/true);
                                                    fitter.setBz(BZ_TESLA);
                                                    
                                                    //float bz = o2::base::Propagator::Instance()->getNominalBz(); // en kGauss
                                                    //DCAFitter2 fitter(bz, /*useAbsDCA=*/true, /*propagateToPCA=*/true);
                                                    fitter.setMaxR(45.f);        // cm
                                                    fitter.setMaxDZIni(4.f);     // cm
                                                    fitter.setMaxDXYIni(4.f);    // cm
                                                    fitter.setMaxChi2(50.f);
                                                    fitter.setPropagateToPCA(true);
                                                    
                                                    std::optional<TracksWCovMc::iterator> pRow, apRow;
                                                    for (auto& tr : tracks) {
                                                        if (!tr.has_mcParticle()) continue;
                                                        const auto& mc = tr.mcParticle();
                                                        if (mc.globalIndex() == antip_id) apRow = tr;
                                                        if (mc.globalIndex() == p_id) pRow  = tr;
                                                        if (pRow && apRow) break;
                                                    }
                                                    if (!(pRow && apRow)) {
                                                    }
                                                    else {
                                                        //TrackParCov
                                                        auto trP  = makeTPCovFromAOD(*pRow);
                                                        auto trAP = makeTPCovFromAOD(*apRow);
                                                        int nCand = fitter.process(trP, trAP);
                                                        auto status = fitter.getFitStatus();
                                                        histos.fill(HIST("vtxfit_status"), (int)status);
                                                        
                                                        if (nCand > 0 && (status == DCAFitter2::FitStatus::Converged || status == DCAFitter2::FitStatus::MaxIter)) {
                                                            //Secondary vertex (commom PCA) [x,y,z] cm
                                                            auto vtx = fitter.getPCACandidatePos();
                                                            const double secX = vtx[0];
                                                            const double secY = vtx[1];
                                                            const double secZ = vtx[2];
                                                            
                                                            //DCA of the pair in the  PCA (equivalent to minDCA)
                                                            fitter.propagateTracksToVertex();
                                                            auto tp0 = fitter.getTrackParamAtPCA(0);
                                                            auto tp1 = fitter.getTrackParamAtPCA(1);
                                                            const auto p0 = tp0.getXYZGlo();
                                                            const auto p1 = tp1.getXYZGlo();
                                                            const double x0 = p0.X(), y0 = p0.Y(), z0 = p0.Z();
                                                            const double x1 = p1.X(), y1 = p1.Y(), z1 = p1.Z();
                                                            const double dcaPair = std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1));
                                                            
                                                            if (mother_pdg == -2112) histos.fill(HIST("cex_pairtrk_vtxfit_dcaPair"), dcaPair);
                                                            if (mother_pdg != -2112) histos.fill(HIST("cexbg_pairtrk_vtxfit_dcaPair"), dcaPair);
                                                            
                                                            //Cuts
                                                            //if (dcaPair > 10.0) continue;
                                                            //if (trkangleDeg > 60) continue;
                                                            //if (TMath::Abs(dplane) > 2) continue;
                                                            if (antip_layers == true && p_layers == true){
                                                                double cex_pairtrk_p = total_trk_pVec.Mag();
                                                                double cex_pairtrk_pt = total_trk_pVec.Pt();
                                                                double cex_pairtrk_pz = p_trk_pz + antip_trk_pz;
                                                                const double radius = std::hypot(secX, secY);
                                                                const double dxPV = secX - pvtx_x;
                                                                const double dyPV = secY - pvtx_y;
                                                                const double dzPV = secZ - pvtx_z;
                                                                const double distToPrimary = std::sqrt(dxPV*dxPV + dyPV*dyPV + dzPV*dzPV);
                                                                
                                                                if(mother_pdg == -2112){
                                                                    histos.fill(HIST("cex_pairtrk_p"), cex_pairtrk_p);
                                                                    histos.fill(HIST("cex_pairtrk_pt"), cex_pairtrk_pt);
                                                                    histos.fill(HIST("cex_pairtrk_pz"), cex_pairtrk_pz);
                                                                    histos.fill(HIST("cex_pairtrk_vtxfit_R"), radius);
                                                                    histos.fill(HIST("cex_pairtrk_vtxfit_distToPV"), distToPrimary);
                                                                    histos.fill(HIST("cex_pairtrk_vtxfit_secVtxXY"), secX, secY);
                                                                    histos.fill(HIST("cex_pairtrk_vtxfit_secVtxZ"), secZ);
                                                                }
                                                                else{
                                                                    histos.fill(HIST("cexbg_pairtrk_p"), cex_pairtrk_p);
                                                                    histos.fill(HIST("cexbg_pairtrk_pt"), cex_pairtrk_pt);
                                                                    histos.fill(HIST("cexbg_pairtrk_pz"), cex_pairtrk_pz);
                                                                    histos.fill(HIST("cexbg_pairtrk_vtxfit_R"), radius);
                                                                    histos.fill(HIST("cexbg_pairtrk_vtxfit_distToPV"), distToPrimary);
                                                                    histos.fill(HIST("cexbg_pairtrk_vtxfit_secVtxXY"), secX, secY);
                                                                    histos.fill(HIST("cexbg_pairtrk_vtxfit_secVtxZ"), secZ);
                                                                }
                                                                
                                                                const float chi2 = fitter.getChi2AtPCACandidate();
                                                                histos.fill(HIST("vtxfit_chi2"), chi2);
                                                                const double dX = secX - antip_vx;
                                                                const double dY = secY - antip_vy;
                                                                const double dZ = secZ - antip_vz;
                                                                const double d3D = std::sqrt(dX*dX + dY*dY + dZ*dZ);
                                                                histos.fill(HIST("vtxfit_mc_dX"), dX);
                                                                histos.fill(HIST("vtxfit_mc_dY"), dY);
                                                                histos.fill(HIST("vtxfit_mc_dZ"), dZ);
                                                                histos.fill(HIST("vtxfit_mc_d3D"), d3D);
                                                                
                                                                // Fill Trees
                                                                buf.mc_pair_p  = cex_pairmc_p;
                                                                buf.mc_pair_pt = cex_pairmc_pt;
                                                                buf.mc_pair_pz = cex_pairmc_pz;
                                                                buf.mc_dplane = dplane;
                                                                buf.mc_angle_deg = mcangleDeg;
                                                                buf.mc_vtx_x = antip_vx;
                                                                buf.mc_vtx_y = antip_vy;
                                                                buf.mc_vtx_z = antip_vz;
                                                                buf.trk_pair_p = cex_pairtrk_p;
                                                                buf.trk_pair_pt = cex_pairtrk_pt;
                                                                buf.trk_pair_pz = cex_pairtrk_pz;
                                                                buf.trk_angle_deg = trkangleDeg;
                                                                buf.trk_vtxfit_dcaPair = dcaPair;
                                                                buf.trk_vtxfit_R = radius;
                                                                buf.trk_vtxfit_distToPV = distToPrimary;
                                                                buf.trk_vtxfit_secVtx_x = secX;
                                                                buf.trk_vtxfit_secVtx_y = secY;
                                                                buf.trk_vtxfit_secVtx_z = secZ;
                                                                buf.vtxfit_chi2 = chi2;
                                                                buf.vtxfit_status = (int)status;
                                                                buf.nCand = nCand;
                                                                buf.vtxfit_dX = secX - antip_vx;
                                                                buf.vtxfit_dY = secY - antip_vy;
                                                                buf.vtxfit_dZ = secZ - antip_vz;
                                                                buf.vtxfit_d3D = std::sqrt(buf.vtxfit_dX*buf.vtxfit_dX + buf.vtxfit_dY*buf.vtxfit_dY + buf.vtxfit_dZ*buf.vtxfit_dZ);
                                                                buf.p_trk_p = p_trk_p;
                                                                buf.p_trk_px = p_trk_px;
                                                                buf.p_trk_py = p_trk_py;
                                                                buf.p_trk_pz = p_trk_pz;
                                                                buf.p_trk_eta = p_trk_eta;
                                                                buf.p_trk_tpcSignal = p_trk_tpcSignal;
                                                                buf.p_trk_nClsITS = p_trk_nClsITS;
                                                                buf.antip_trk_p = antip_trk_p;
                                                                buf.antip_trk_px = antip_trk_px;
                                                                buf.antip_trk_py = antip_trk_py;
                                                                buf.antip_trk_pz = antip_trk_pz;
                                                                buf.antip_trk_eta = antip_trk_eta;
                                                                buf.antip_trk_tpcSignal = antip_trk_tpcSignal;
                                                                buf.antip_trk_nClsITS = antip_trk_nClsITS;
                                                                buf.mother_pdg = mother_pdg;
                                                                buf.colId = colId;
                                                                buf.p_id = p_id;
                                                                buf.antip_id = antip_id;

                                                                if (mother_pdg == -2112) cexTree->Fill();
                                                                else bgTree->Fill();
                                                            }
                                                        }
                                                    }
                                                    // ==== end DCAFitter2 ====
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
    return WorkflowSpec{adaptAnalysisTask<AntineutronTask>(ctx)};
}

