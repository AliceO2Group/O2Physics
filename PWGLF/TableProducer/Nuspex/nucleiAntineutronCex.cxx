// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file nucleiAntineutronCex.cxx
/// \brief Analysis task for antineutron detection through cex interactions
/// \author Fabiola Lugo
///

#include <PWGLF/DataModel/LFAntinCexTables.h>

#include <Common/DataModel/PIDResponseITS.h>

#include <CommonConstants/MathConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <TMCProcess.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>

using namespace o2;
using namespace o2::framework;
using o2::constants::math::Rad2Deg;

struct NucleiAntineutronCex {
  // Slicing per colision
  Preslice<aod::McParticles> perMcByColl = aod::mcparticle::mcCollisionId;
  // Check available tables in the AOD, specifically TracksIU, TracksCovIU
  using TracksWCovMc = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels, o2::aod::TracksCovIU>;

  // === Cut values ===
  static constexpr double kIts2MinR = 2.2;    // ITS2 min radius [cm]
  static constexpr double kIts2MaxR = 39.0;   // ITS2 max radius [cm]
  static constexpr double kIts2MaxVz = 39.0;  // ITS2 max |vz| [cm]
  static constexpr double kAccMaxEta = 1.2;   // acceptance |eta|
  static constexpr double kAccMaxVz = 5.3;    // acceptance |vz| [cm]
  static constexpr double kStrictEta = 0.9;   // tighter eta cut
  static constexpr double kInitDplane = 10.0; // init dplane
  static constexpr double kHuge = 1e9;        // fallback for bad denom
  static constexpr int kMinItsHits = 2;
  static constexpr double kVtxTol = 1e-4;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Produces<aod::AntinCexPairs> outPairs;

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
    histos.add("antin_eta", "Pseudorapidity;#eta;Entries", kTH1F, {{100, -10., 10.}});
    histos.add("antin_p_ITScuts", "Momentum with ITS cuts;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

    // Primary neutrons
    histos.add("n_p", "Total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("n_px", "p_{x};p_{x} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("n_py", "p_{y};p_{y} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("n_pz", "p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("n_eta", "Pseudorapidity;#eta;Entries", kTH1F, {{100, -10., 10.}});
    histos.add("n_p_ITScuts", "Momentum with ITS cuts;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

    // Primary antiprotons
    histos.add("antipP", "Total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("antipPx", "p_{x};p_{x} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("antipPy", "p_{y};p_{y} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("antipPz", "p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("antipEta", "Pseudorapidity;#eta;Entries", kTH1F, {{100, -10., 10.}});
    histos.add("antipP_ITScuts", "Momentum with ITS cuts;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

    // Primary protons
    histos.add("pP", "Total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("pPx", "p_{x};p_{x} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("pPy", "p_{y};p_{y} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("pPz", "p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("pEta", "Pseudorapidity;#eta;Entries", kTH1F, {{100, -10., 10.}});
    histos.add("pP_ITScuts", "Momentum with ITS cuts;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

    // test (MC)
    histos.add("antip_test", "Secondary antiprotons;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

    // CEX pair from antineutron (MC)
    histos.add("cexPairMcP", "CEX pair total momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("cexPairMcPt", "CEX pair p_{T};p_{T} (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("cexPairMcPz", "CEX pair p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("cex_pairmcDplane", "CEX pair d_{plane};d_{plane} (cm);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("cex_pairmc_angle", "Pair opening angle;Angle (°);Entries", kTH1F, {{180, 0., 180.}});
    histos.add("cex_pairmc_vtx", "MC CEX pair vertex;X (cm);Y (cm)", kTH2F, {{100, -50., 50.}, {100, -50., 50.}});
    histos.add("cex_pairmc_vtxz", "MC secondary vertex Z;Z (cm);Entries", kTH1F, {{200, -60., 60.}});
    histos.add("cexPairMcPITScuts", "CEX pair momentum (ITS cuts);|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

    // CEX pair normalized to antineutron (MC)
    histos.add("cexn_pairmc_p", "Pair p / antineutron p;p/p_{#bar{n}};Entries", kTH1F, {{100, 0., 2.}});
    histos.add("cexn_pairmc_pt", "Pair p_{T} / antineutron p_{T};p_{T}/p_{T,#bar{n}};Entries", kTH1F, {{100, 0., 2.}});
    histos.add("cexn_pairmc_pz", "Pair p_{z} / antineutron p_{z};p_{z}/p_{z,#bar{n}};Entries", kTH1F, {{100, -2., 2.}});

    // BG pair (not from antineutron) (MC)
    histos.add("cexbg_pairmc_p", "Background pair momentum;|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("cexbg_pairmc_pt", "Background pair p_{T};p_{T} (GeV/c);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("cexbg_pairmc_pz", "Background pair p_{z};p_{z} (GeV/c);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("cexbg_pairmcDplane", "Background d_{plane};d_{plane} (cm);Entries", kTH1F, {{100, 0., 10.}});
    histos.add("cexbg_pairmc_angle", "Background opening angle;Angle (°);Entries", kTH1F, {{180, 0., 180.}});
    histos.add("cexbg_pairmc_vtx", "Background pair vertex;X (cm);Y (cm)", kTH2F, {{100, -50., 50.}, {100, -50., 50.}});
    histos.add("cexbg_pairmc_vtxz", "Background secondary vertex Z;Z (cm);Entries", kTH1F, {{200, -60., 60.}});
    histos.add("cexbg_pairmc_pITScuts", "Background momentum (ITS cuts);|p| (GeV/c);Entries", kTH1F, {{100, 0., 10.}});

    // CEX pair from antineutron (TRK)
    histos.add("cex_pairtrk_angle", "Pair opening angle (tracks);Angle (°);Entries", kTH1F, {{180, 0., 180.}});
    histos.add("cexPairTrkP", "Pair momentum (tracks);|p| (GeV/c);Entries", kTH1F, {{120, 0., 12.}});
    histos.add("cexPairTrkPt", "Pair p_{T} (tracks);p_{T} (GeV/c);Entries", kTH1F, {{120, 0., 12.}});
    histos.add("cexPairTrkPz", "Pair p_{z} (tracks);p_{z} (GeV/c);Entries", kTH1F, {{120, -12., 12.}});
    histos.add("cex_pairtrkVtxfitDcaPair", "DCA between tracks at PCA;DCA (cm);Entries", kTH1F, {{200, 0., 10.}});
    histos.add("cex_pairtrkVtxfitR", "Secondary-vertex radius (PCA);R (cm);Entries", kTH1F, {{200, 0., 60.}});
    histos.add("cex_pairtrkVtxfitDistToPv", "Distance from secondary vertex to PV;dist (cm);Entries", kTH1F, {{240, 0., 120.}});
    histos.add("cex_pairtrk_vtxfit_secVtxXY", "Secondary vertex (PCA);X (cm);Y (cm)", kTH2F, {{200, -60., 60.}, {200, -60., 60.}});
    histos.add("cex_pairtrk_vtxfit_secVtxZ", "Secondary vertex Z (PCA);Z (cm);Entries", kTH1F, {{240, -60., 60.}});

    // BG pair (not from antineutron) (TRK)
    histos.add("cexbg_pairtrk_angle", "Background opening angle (tracks);Angle (°);Entries", kTH1F, {{180, 0., 180.}});
    histos.add("cexbg_pairtrk_p", "Pair momentum (tracks);|p| (GeV/c);Entries", kTH1F, {{120, 0., 12.}});
    histos.add("cexbg_pairtrk_pt", "Pair p_{T} (tracks);p_{T} (GeV/c);Entries", kTH1F, {{120, 0., 12.}});
    histos.add("cexbg_pairtrk_pz", "Pair p_{z} (tracks);p_{z} (GeV/c);Entries", kTH1F, {{120, -12., 12.}});
    histos.add("cexbg_pairtrkVtxfitDcaPair", "DCA between tracks at PCA;DCA (cm);Entries", kTH1F, {{200, 0., 10.}});
    histos.add("cexbg_pairtrkVtxfitR", "Secondary-vertex radius (PCA);R (cm);Entries", kTH1F, {{200, 0., 60.}});
    histos.add("cexbg_pairtrkVtxfitDistToPv", "Distance from secondary vertex to PV;dist (cm);Entries", kTH1F, {{240, 0., 120.}});
    histos.add("cexbg_pairtrk_vtxfit_secVtxXY", "Secondary vertex (PCA);X (cm);Y (cm)", kTH2F, {{200, -60., 60.}, {200, -60., 60.}});
    histos.add("cexbg_pairtrk_vtxfit_secVtxZ", "Secondary vertex Z (PCA);Z (cm);Entries", kTH1F, {{240, -60., 60.}});

    // Vertex fit (DCAFitter2 / PCA)
    histos.add("vtxfitChi2", "DCAFitter2 #chi^{2};#chi^{2};Entries", kTH1F, {{200, 0., 100.}});
    histos.add("vtxfitStatus", "Fit status (0=OK);code;Entries", kTH1I, {{10, 0., 10.}});
    histos.add("vtxfit_mc_dX", "SV residual X (fit - MC);#Delta X (cm);Entries", kTH1F, {{400, -20., 20.}});
    histos.add("vtxfit_mc_dY", "SV residual Y (fit - MC);#Delta Y (cm);Entries", kTH1F, {{400, -20., 20.}});
    histos.add("vtxfit_mc_dZ", "SV residual Z (fit - MC);#Delta Z (cm);Entries", kTH1F, {{400, -20., 20.}});
    histos.add("vtxfit_mc_d3D", "SV distance |fit - MC|;#Delta r (cm);Entries", kTH1F, {{300, 0., 30.}});

    // ITS PID (protons / antiprotons, reconstructed tracks)
    histos.add("pItsNsigmaPr", "ITS n#sigma (p hyp., proton);n#sigma_{ITS}(p);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("apItsNsigmaPr", "ITS n#sigma (p hyp., antiproton);n#sigma_{ITS}(p);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("pItsPidValid", "ITS PID valid flag (proton);PidValid;Entries", kTH1F, {{2, 0., 2.}});
    histos.add("apItsPidValid", "ITS PID valid flag (antiproton);PidValid;Entries", kTH1F, {{2, 0., 2.}});
    histos.add("pTgl", "tgl (proton track);tgl;Entries", kTH1F, {{100, -2., 2.}});
    histos.add("apTgl", "tgl (antiproton track);tgl;Entries", kTH1F, {{100, -2., 2.}});
    histos.add("pItsNsigmaPr_bg", "ITS n#sigma (p hyp., proton);n#sigma_{ITS}(p);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("apItsNsigmaPr_bg", "ITS n#sigma (p hyp., antiproton);n#sigma_{ITS}(p);Entries", kTH1F, {{100, -10., 10.}});
    histos.add("pItsPidValid_bg", "ITS PID valid flag (proton);PidValid;Entries", kTH1F, {{2, 0., 2.}});
    histos.add("apItsPidValid_bg", "ITS PID valid flag (antiproton);PidValid;Entries", kTH1F, {{2, 0., 2.}});
    histos.add("pTgl_bg", "tgl (proton track);tgl;Entries", kTH1F, {{100, -2., 2.}});
    histos.add("apTgl_bg", "tgl (antiproton track);tgl;Entries", kTH1F, {{100, -2., 2.}});
  }

  static o2::track::TrackParCov makeTPCovFromAOD(const TracksWCovMc::iterator& tr)
  {
    using o2::track::TrackParCov;
    TrackParCov tpcov;
    // Local state: x, alpha, y, z, snp, tgl, q/pt
    float par[5] = {tr.y(), tr.z(), tr.snp(), tr.tgl(), tr.signed1Pt()};
    tpcov.set(tr.x(), tr.alpha(), par);

    // Covariance matrix (15 terms) in O2 order
    std::array<float, 15> cov = {
      tr.cYY(), tr.cZY(), tr.cZZ(),
      tr.cSnpY(), tr.cSnpZ(), tr.cSnpSnp(),
      tr.cTglY(), tr.cTglZ(), tr.cTglSnp(),
      tr.cTglTgl(), tr.c1PtY(), tr.c1PtZ(),
      tr.c1PtSnp(), tr.c1PtTgl(), tr.c1Pt21Pt2()};
    tpcov.setCov(cov);
    return tpcov;
  }

  void process(aod::McCollisions const& cols, aod::McParticles const& particles, TracksWCovMc const& tracks)
  {
    double pvtxX = 0;
    double pvtxY = 0;
    double pvtxZ = 0;
    for (auto const& col : cols) {
      const auto colId = col.globalIndex();
      auto mcPartsThis = particles.sliceBy(perMcByColl, colId);

      if (std::isfinite(col.posX()) && std::isfinite(col.posY()) && std::isfinite(col.posZ())) {
        pvtxX = col.posX();
        pvtxY = col.posY();
        pvtxZ = col.posZ();
        histos.fill(HIST("hVx"), pvtxX);
        histos.fill(HIST("hVy"), pvtxY);
        histos.fill(HIST("hVz"), pvtxZ);
      }

      for (const auto& particle : mcPartsThis) {

        // Primary antineutrons
        if (particle.pdgCode() == -kNeutron && particle.isPhysicalPrimary()) {
          histos.fill(HIST("antin_p"), particle.p());
          histos.fill(HIST("antin_px"), particle.px());
          histos.fill(HIST("antin_py"), particle.py());
          histos.fill(HIST("antin_pz"), particle.pz());
          histos.fill(HIST("antin_eta"), particle.eta());
          if (std::abs(particle.eta()) < kAccMaxEta && std::abs(particle.vz()) < kAccMaxVz)
            histos.fill(HIST("antin_p_ITScuts"), particle.p());
        }
        // Primary neutrons
        if (particle.pdgCode() == kNeutron && particle.isPhysicalPrimary()) {
          histos.fill(HIST("n_p"), particle.p());
          histos.fill(HIST("n_px"), particle.px());
          histos.fill(HIST("n_py"), particle.py());
          histos.fill(HIST("n_pz"), particle.pz());
          histos.fill(HIST("n_eta"), particle.eta());
          if (std::abs(particle.eta()) < kAccMaxEta && std::abs(particle.vz()) < kAccMaxVz)
            histos.fill(HIST("n_p_ITScuts"), particle.p());
        }
        // Primary antiprotons
        if (particle.pdgCode() == -kProton && particle.isPhysicalPrimary()) {
          histos.fill(HIST("antipP"), particle.p());
          histos.fill(HIST("antipPx"), particle.px());
          histos.fill(HIST("antipPy"), particle.py());
          histos.fill(HIST("antipPz"), particle.pz());
          histos.fill(HIST("antipEta"), particle.eta());
          if (std::abs(particle.eta()) < kAccMaxEta && std::abs(particle.vz()) < kAccMaxVz)
            histos.fill(HIST("antipP_ITScuts"), particle.p());
        }
        // Primary protons
        if (particle.pdgCode() == kProton && particle.isPhysicalPrimary()) {
          histos.fill(HIST("pP"), particle.p());
          histos.fill(HIST("pPx"), particle.px());
          histos.fill(HIST("pPy"), particle.py());
          histos.fill(HIST("pPz"), particle.pz());
          histos.fill(HIST("pEta"), particle.eta());
          if (std::abs(particle.eta()) < kAccMaxEta && std::abs(particle.vz()) < kAccMaxVz)
            histos.fill(HIST("pP_ITScuts"), particle.p());
        }

        // Seconday antiprotons from material
        const auto procEnum = particle.getProcess();
        const bool isSecondaryFromMaterial = (!particle.producedByGenerator()) && (procEnum == kPHadronic || procEnum == kPHInhelastic);
        if (particle.pdgCode() != -kProton || !isSecondaryFromMaterial || particle.mothersIds().empty())
          continue;
        histos.fill(HIST("antip_test"), particle.p());

        // Primary mother
        bool hasPrimaryMotherAntip = false;
        double motherPt = 0.0;
        double motherPz = 0.0;
        double motherVz = 0.0;
        double motherP = 0.0;
        double motherEta = 0.0;
        int motherPdg = 0;

        for (const auto& mother : particle.mothers_as<aod::McParticles>()) {
          if (mother.isPhysicalPrimary()) {
            hasPrimaryMotherAntip = true;
            motherPt = mother.pt();
            motherPz = mother.pz();
            motherVz = mother.vz();
            motherP = mother.p();
            motherEta = mother.eta();
            motherPdg = mother.pdgCode();
            break;
          }
        }
        if (!hasPrimaryMotherAntip)
          continue;

        double antipVx = particle.vx();
        double antipVy = particle.vy();
        double antipVz = particle.vz();
        double antipPx = particle.px();
        double antipPy = particle.py();
        double antipPz = particle.pz();
        double antipE = particle.e();
        int antipId = particle.globalIndex();

        // Selection conditions: Produced in the ITS
        const double r = std::sqrt(antipVx * antipVx + antipVy * antipVy);
        // Config for ITS
        // if(3.9<=r && r<=43.0 && std::abs(antipVz)<=48.9){
        // Config for ITS2
        if (r < kIts2MinR || r > kIts2MaxR || std::abs(antipVz) > kIts2MaxVz)
          continue;
        if (std::abs(motherEta) >= kAccMaxEta || std::abs(motherVz) >= kAccMaxVz)
          continue;

        // Pion minus veto
        bool pionMinus = false;
        for (const auto& particle1 : mcPartsThis) {
          const auto proc1Enum = particle1.getProcess();
          const bool isSecondaryFromMaterial1 = (!particle1.producedByGenerator()) && (proc1Enum == kPHadronic || proc1Enum == kPHInhelastic);
          if (particle1.mcCollisionId() != colId)
            continue;
          if (particle1.pdgCode() != kPiMinus || !isSecondaryFromMaterial1 || particle1.mothersIds().empty())
            continue;
          bool hasPrimaryMotherPim = false;
          for (const auto& mother : particle1.mothers_as<aod::McParticles>()) {
            if (mother.isPhysicalPrimary()) {
              hasPrimaryMotherPim = true;
              break;
            }
          }
          if (!hasPrimaryMotherPim)
            continue;
          double pimVx = particle1.vx();
          double pimVy = particle1.vy();
          double pimVz = particle1.vz();
          if (std::abs(pimVx - antipVx) < kVtxTol && std::abs(pimVy - antipVy) < kVtxTol && std::abs(pimVz - antipVz) < kVtxTol) {
            pionMinus = true;
            break;
          }
        }

        // Pion plus veto
        bool pionPlus = false;
        for (const auto& particle2 : mcPartsThis) {
          if (particle2.mcCollisionId() != colId)
            continue;
          const auto proc2Enum = particle2.getProcess();
          const bool isSecondaryFromMaterial2 = (!particle2.producedByGenerator()) && (proc2Enum == kPHadronic || proc2Enum == kPHInhelastic);
          if (particle2.pdgCode() != kPiPlus || !isSecondaryFromMaterial2 || particle2.mothersIds().empty())
            continue;
          bool hasPrimaryMotherPip = false;
          for (const auto& mother : particle2.mothers_as<aod::McParticles>()) {
            if (mother.isPhysicalPrimary()) {
              hasPrimaryMotherPip = true;
              break;
            }
          }
          if (!hasPrimaryMotherPip)
            continue;
          double pipVx = particle2.vx();
          double pipVy = particle2.vy();
          double pipVz = particle2.vz();
          if (std::abs(pipVx - antipVx) < kVtxTol && std::abs(pipVy - antipVy) < kVtxTol && std::abs(pipVz - antipVz) < kVtxTol) {
            pionPlus = true;
            break;
          }
        }

        if (pionPlus || pionMinus)
          continue;

        // CEX selection
        double dplane = kInitDplane;
        double dplaneTmp = 0;
        double p = 0;
        double pTmp = 0;
        double pcexPx = 0;
        double pcexPy = 0;
        double pcexPz = 0;
        double e = 0;
        double eTmp = 0;
        int k_plane = -1;
        int k_e = -1;
        int k_p = -1;

        // Secondary proton from material
        for (const auto& particle3 : mcPartsThis) {
          if (particle3.mcCollisionId() != colId)
            continue;
          const auto proc3Enum = particle3.getProcess();
          const bool isSecondaryFromMaterial3 = (!particle3.producedByGenerator()) && (proc3Enum == kPHadronic || proc3Enum == kPHInhelastic);
          if (particle3.pdgCode() != kProton || !isSecondaryFromMaterial3 || particle3.mothersIds().empty())
            continue;
          bool hasPrimaryMotherP = false;
          for (const auto& mother : particle3.mothers_as<aod::McParticles>()) {
            if (mother.isPhysicalPrimary()) {
              hasPrimaryMotherP = true;
              break;
            }
          }
          if (!hasPrimaryMotherP)
            continue;
          double protonVx = particle3.vx();
          double protonVy = particle3.vy();
          double protonVz = particle3.vz();
          double pPx = particle3.px();
          double pPy = particle3.py();
          double pPz = particle3.pz();
          double pE = particle3.e();
          if (std::abs(protonVx - antipVx) < kVtxTol && std::abs(protonVy - antipVy) < kVtxTol && std::abs(protonVz - antipVz) < kVtxTol) {
            // Same mother
            bool shareMother = false;
            const auto& momsAp = particle.mothersIds(); // antiproton
            const auto& momsP = particle3.mothersIds(); // proton
            for (const auto& ida : momsAp) {
              for (const auto& idp : momsP) {
                if (ida == idp) {
                  shareMother = true;
                  break;
                }
              }
              if (shareMother)
                break;
            }
            if (!shareMother)
              continue;

            // CEX proton selection
            // dplaneTmp = (pPy*antipPz - pPz*antipPy)*(pvtxX-antipVx) + (pPz*antipPx - pPx*antipPz)*(pvtxY-antipVy) + (pPx*antipPy - pPy*antipPx)*(pvtxZ-antipVz);
            double nx = (pPy * antipPz - pPz * antipPy);
            double ny = (pPz * antipPx - pPx * antipPz);
            double nz = (pPx * antipPy - pPy * antipPx);
            double rx = (pvtxX - antipVx);
            double ry = (pvtxY - antipVy);
            double rz = (pvtxZ - antipVz);
            double denom = nx * nx + ny * ny + nz * nz;
            if (denom > 0.) {
              dplaneTmp = std::abs(nx * rx + ny * ry + nz * rz) / std::sqrt(denom);
            } else {
              dplaneTmp = kHuge;
            }
            if (std::abs(dplaneTmp) < std::abs(dplane)) {
              k_plane = particle3.globalIndex();
              dplane = dplaneTmp;
            }

            eTmp = antipE + pE;
            if (std::abs(eTmp) > std::abs(e)) {
              k_e = particle3.globalIndex();
              e = eTmp;
            }

            pTmp = std::sqrt(std::pow((pPx + antipPx), 2) + std::pow((pPy + antipPy), 2) + std::pow((pPz + antipPz), 2));
            if (std::abs(pTmp) > std::abs(p)) {
              k_p = particle3.globalIndex();
              p = pTmp;
              pcexPx = pPx;
              pcexPy = pPy;
              pcexPz = pPz;
            }
          }
        }

        if (k_plane == k_e && k_plane == k_p && k_plane >= 0) {
          int pId = k_plane;
          TVector3 pVecProton = TVector3(pcexPx, pcexPy, pcexPz);
          TVector3 pVecAntiproton = TVector3(antipPx, antipPy, antipPz);
          TVector3 total_mc_pVec = pVecProton + pVecAntiproton;
          double cexPairMcP = total_mc_pVec.Mag();
          double cexPairMcPt = total_mc_pVec.Pt();
          double cexPairMcPz = pcexPz + antipPz;
          double mcangleRad = pVecProton.Angle(pVecAntiproton);
          double mcangleDeg = mcangleRad * Rad2Deg;

          // Antineutron mother
          if (motherPdg == -kNeutron) {
            // CEX pair
            histos.fill(HIST("cexPairMcP"), cexPairMcP);
            histos.fill(HIST("cexPairMcPt"), cexPairMcPt);
            histos.fill(HIST("cexPairMcPz"), cexPairMcPz);
            histos.fill(HIST("cex_pairmcDplane"), dplane);
            histos.fill(HIST("cex_pairmc_angle"), mcangleDeg);
            histos.fill(HIST("cex_pairmc_vtx"), antipVx, antipVy);
            histos.fill(HIST("cex_pairmc_vtxz"), antipVz);
            if (std::abs(motherEta) < kStrictEta && std::abs(motherVz) < kAccMaxVz)
              histos.fill(HIST("cexPairMcPITScuts"), cexPairMcP);
            // CEX pair normalized
            if (motherP != 0)
              histos.fill(HIST("cexn_pairmc_p"), cexPairMcP / motherP);
            if (motherPt != 0)
              histos.fill(HIST("cexn_pairmc_pt"), cexPairMcPt / motherPt);
            if (motherPz != 0)
              histos.fill(HIST("cexn_pairmc_pz"), cexPairMcPz / motherPz);
          }
          // BG mother
          if (motherPdg != -kNeutron) {
            // CEX pair
            histos.fill(HIST("cexbg_pairmc_p"), cexPairMcP);
            histos.fill(HIST("cexbg_pairmc_pt"), cexPairMcPt);
            histos.fill(HIST("cexbg_pairmc_pz"), cexPairMcPz);
            histos.fill(HIST("cexbg_pairmcDplane"), dplane);
            histos.fill(HIST("cexbg_pairmc_angle"), mcangleDeg);
            histos.fill(HIST("cexbg_pairmc_vtx"), antipVx, antipVy);
            histos.fill(HIST("cexbg_pairmc_vtxz"), antipVz);
            if (std::abs(motherEta) < kStrictEta && std::abs(motherVz) < kAccMaxVz)
              histos.fill(HIST("cexbg_pairmc_pITScuts"), cexPairMcP);
          }

          // Detector signal
          bool antipLayers = false;
          bool antipHasTrack = false;
          double antipTrkPx = 0.;
          double antipTrkPy = 0.;
          double antipTrkPz = 0.;
          double antipTrkP = 0.;
          double antipTrkEta = 0.;
          double antipTrkTpcSignal = 0;
          // int antip_trk_nClsTPC = 0;
          int antipTrkNClsIts = 0;
          uint16_t apItsMap = 0;
          float pTrkItsNSigmaPr = -999.f;
          int8_t pTrkItsPidValid = 0;
          float pTrkTgl = 0.f;

          bool pLayers = false;
          bool pHasTrack = false;
          double pTrkPx = 0.;
          double pTrkPy = 0.;
          double pTrkPz = 0.;
          double pTrkP = 0.;
          double pTrkEta = 0.;
          double pTrkTpcSignal = 0;
          // int p_trk_nClsTPC = 0;
          int pTrkNClsIts = 0;
          uint16_t pItsMap = 0;
          float antipTrkItsNSigmaPr = -999.f;
          int8_t antipTrkItsPidValid = 0;
          float antipTrkTgl = 0.f;

          o2::aod::ITSResponse itsResponse;

          for (const auto& track : tracks) {
            if (!track.has_mcParticle())
              continue;
            const auto& mc = track.mcParticle();
            if (mc.mcCollisionId() != colId)
              continue;
            uint8_t itsMap = track.itsClusterMap();
            // Config for ITS1
            /*bool hitSPD = (itsMap & 0x3) != 0; // bits 0 (SPD L1) & 1 (SPD L2)
             bool hitSDD  = (itsMap & 0xC) != 0; // bits 2–3
             bool hitSSD  = (itsMap & 0x30) != 0; // bits 4–5
             bool layerCondition = (hitSDD || hitSSD) && !hitSPD;*/
            // Config for ITS2
            bool hitL0 = (itsMap & (1u << 0)) != 0;
            bool hitL1 = (itsMap & (1u << 1)) != 0;
            bool hitL2 = (itsMap & (1u << 2)) != 0;
            bool hitL3 = (itsMap & (1u << 3)) != 0;
            bool hitL4 = (itsMap & (1u << 4)) != 0;
            bool hitL5 = (itsMap & (1u << 5)) != 0;
            bool hitL6 = (itsMap & (1u << 6)) != 0;
            bool hitIB = (hitL0 || hitL1 || hitL2);
            bool hitOuter = (hitL3 || hitL4 || hitL5 || hitL6);
            int nITS = track.itsNCls();
            bool layerCondition = (!hitIB) && hitOuter && (nITS >= kMinItsHits);

            if (mc.globalIndex() == antipId) {
              antipTrkP = track.p();
              antipTrkPx = track.px();
              antipTrkPy = track.py();
              antipTrkPz = track.pz();
              antipTrkEta = track.eta();
              antipTrkTpcSignal = track.tpcSignal();
              // antip_trk_nClsTPC = track.tpcNCls();
              antipTrkNClsIts = track.itsNCls();
              antipTrkTgl = track.tgl();
              const auto nsigmaITSantip = itsResponse.nSigmaITS<o2::track::PID::Proton>(track);
              antipTrkItsNSigmaPr = static_cast<float>(nsigmaITSantip);
              antipTrkItsPidValid = std::isfinite(nsigmaITSantip) ? 1 : 0;
              antipHasTrack = true;
              apItsMap = static_cast<uint16_t>(track.itsClusterMap());
              antipLayers = (apItsMap != 0);
              if (layerCondition)
                antipLayers = true;
              if (motherPdg == -kNeutron) {
                histos.fill(HIST("apItsNsigmaPr"), antipTrkItsNSigmaPr);
                histos.fill(HIST("apItsPidValid"), antipTrkItsPidValid);
                histos.fill(HIST("apTgl"), antipTrkTgl);
              }
              if (motherPdg != -kNeutron) {
                histos.fill(HIST("apItsNsigmaPr_bg"), antipTrkItsNSigmaPr);
                histos.fill(HIST("apItsPidValid_bg"), antipTrkItsPidValid);
                histos.fill(HIST("apTgl_bg"), antipTrkTgl);
              }
            } else if (mc.globalIndex() == pId) {
              pTrkP = track.p();
              pTrkPx = track.px();
              pTrkPy = track.py();
              pTrkPz = track.pz();
              pTrkEta = track.eta();
              pTrkTpcSignal = track.tpcSignal();
              // p_trk_nClsTPC = track.tpcNCls();
              pTrkNClsIts = track.itsNCls();
              pTrkTgl = track.tgl();
              const auto nsigmaITSp =
                itsResponse.nSigmaITS<o2::track::PID::Proton>(track);
              pTrkItsNSigmaPr = static_cast<float>(nsigmaITSp);
              pTrkItsPidValid = std::isfinite(nsigmaITSp) ? 1 : 0;
              pHasTrack = true;
              pItsMap = static_cast<uint16_t>(track.itsClusterMap());
              pLayers = (pItsMap != 0);
              if (layerCondition)
                pLayers = true;
              if (motherPdg == -kNeutron) {
                histos.fill(HIST("pItsNsigmaPr"), pTrkItsNSigmaPr);
                histos.fill(HIST("pItsPidValid"), pTrkItsPidValid);
                histos.fill(HIST("pTgl"), pTrkTgl);
              }
              if (motherPdg != -kNeutron) {
                histos.fill(HIST("pItsNsigmaPr_bg"), pTrkItsNSigmaPr);
                histos.fill(HIST("pItsPidValid_bg"), pTrkItsPidValid);
                histos.fill(HIST("pTgl_bg"), pTrkTgl);
              }
            }
          }
          if (!(pHasTrack && antipHasTrack))
            continue;

          TVector3 pVecProton_trk(pTrkPx, pTrkPy, pTrkPz);
          TVector3 AntipVecProton_trk(antipTrkPx, antipTrkPy, antipTrkPz);
          TVector3 total_trk_pVec = pVecProton_trk + AntipVecProton_trk;
          double trkangleRad = AntipVecProton_trk.Angle(pVecProton_trk);
          double trkangleDeg = trkangleRad * Rad2Deg;
          if (motherPdg == -kNeutron)
            histos.fill(HIST("cex_pairtrk_angle"), trkangleDeg);
          if (motherPdg != -kNeutron)
            histos.fill(HIST("cexbg_pairtrk_angle"), trkangleDeg);

          // ==== Secondary vertex via central DCA vertexer (DCAFitter2) ====
          using o2::vertexing::DCAFitter2;
          constexpr float kBzTesla = 0.5f;
          DCAFitter2 fitter(/*bz=*/kBzTesla, /*useAbsDCA=*/true, /*propagateToPCA=*/true);
          fitter.setBz(kBzTesla);
          // float bz = o2::base::Propagator::Instance()->getNominalBz(); // en kGauss
          // DCAFitter2 fitter(bz, /*useAbsDCA=*/true, /*propagateToPCA=*/true);
          fitter.setMaxR(45.f);     // cm
          fitter.setMaxDZIni(4.f);  // cm
          fitter.setMaxDXYIni(4.f); // cm
          fitter.setMaxChi2(50.f);
          fitter.setPropagateToPCA(true);
          std::optional<TracksWCovMc::iterator> pRow, apRow;
          for (const auto& tr : tracks) {
            if (!tr.has_mcParticle())
              continue;
            const auto& mc = tr.mcParticle();
            if (mc.globalIndex() == antipId)
              apRow = tr;
            if (mc.globalIndex() == pId)
              pRow = tr;
            if (pRow && apRow)
              break;
          }
          if (pRow && apRow) {
            // TrackParCov
            auto trP = makeTPCovFromAOD(*pRow);
            auto trAP = makeTPCovFromAOD(*apRow);
            int nCand = fitter.process(trP, trAP);
            auto status = fitter.getFitStatus();
            histos.fill(HIST("vtxfitStatus"), static_cast<int>(status));
            if (nCand > 0 && (status == DCAFitter2::FitStatus::Converged || status == DCAFitter2::FitStatus::MaxIter)) {
              // Secondary vertex (commom PCA) [x,y,z] cm
              auto vtx = fitter.getPCACandidatePos();
              const double secX = vtx[0];
              const double secY = vtx[1];
              const double secZ = vtx[2];
              // DCA of the pair in the  PCA (equivalent to minDCA)
              fitter.propagateTracksToVertex();
              auto tp0 = fitter.getTrackParamAtPCA(0);
              auto tp1 = fitter.getTrackParamAtPCA(1);
              const auto p0 = tp0.getXYZGlo();
              const auto p1 = tp1.getXYZGlo();
              const double x0 = p0.X(), y0 = p0.Y(), z0 = p0.Z();
              const double x1 = p1.X(), y1 = p1.Y(), z1 = p1.Z();
              const double dcaPair = std::sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1));

              if (motherPdg == -kNeutron)
                histos.fill(HIST("cex_pairtrkVtxfitDcaPair"), dcaPair);
              if (motherPdg != -kNeutron)
                histos.fill(HIST("cexbg_pairtrkVtxfitDcaPair"), dcaPair);

              if (!(antipLayers && pLayers))
                continue;
              double cexPairTrkP = total_trk_pVec.Mag();
              double cexPairTrkPt = total_trk_pVec.Pt();
              double cexPairTrkPz = pTrkPz + antipTrkPz;
              const double radius = std::hypot(secX, secY);
              const double dxPv = secX - pvtxX;
              const double dyPv = secY - pvtxY;
              const double dzPv = secZ - pvtxZ;
              const double distToPrimary = std::sqrt(dxPv * dxPv + dyPv * dyPv + dzPv * dzPv);

              const TVector3 pv2sv(secX - pvtxX, secY - pvtxY, secZ - pvtxZ);
              const double pairPointingAngleDeg = pv2sv.Angle(total_trk_pVec) * Rad2Deg;

              const double pP = pVecProton_trk.Mag();
              const double pAP = AntipVecProton_trk.Mag();
              const double ptP = pVecProton_trk.Pt();
              const double ptAP = AntipVecProton_trk.Pt();

              const double denomP = std::max(1e-9, pP + pAP);
              const double denomPt = std::max(1e-9, ptP + ptAP);

              const float pairPBalance = std::abs(pP - pAP) / denomP;
              const float pairPtBalance = std::abs(ptP - ptAP) / denomPt;

              const float pairQ = (pVecProton_trk - AntipVecProton_trk).Mag();

              // Trk - MC
              const float dPairP = cexPairTrkP - cexPairMcP;
              const float dPairPt = cexPairTrkPt - cexPairMcPt;
              const float dPairPz = cexPairTrkPz - cexPairMcPz;
              const float dOpenAngle = trkangleDeg - mcangleDeg;

              // Closest ITS layer: Radius need to be checked
              static const std::array<double, 7> rLayers = {2.2, 2.8, 3.6, 19.6, 24.0, 29.0, 35.0};
              int16_t svNearestLayerId = -1;
              float svDeltaRToLayer = 1e9f;
              for (int i = 0; i < static_cast<int>(rLayers.size()); ++i) {
                const float dR = static_cast<float>(std::abs(radius - rLayers[i]));
                if (dR < svDeltaRToLayer) {
                  svDeltaRToLayer = dR;
                  svNearestLayerId = static_cast<int16_t>(i);
                }
              }

              if (motherPdg == -kNeutron) {
                histos.fill(HIST("cexPairTrkP"), cexPairTrkP);
                histos.fill(HIST("cexPairTrkPt"), cexPairTrkPt);
                histos.fill(HIST("cexPairTrkPz"), cexPairTrkPz);
                histos.fill(HIST("cex_pairtrkVtxfitR"), radius);
                histos.fill(HIST("cex_pairtrkVtxfitDistToPv"), distToPrimary);
                histos.fill(HIST("cex_pairtrk_vtxfit_secVtxXY"), secX, secY);
                histos.fill(HIST("cex_pairtrk_vtxfit_secVtxZ"), secZ);
              } else {
                histos.fill(HIST("cexbg_pairtrk_p"), cexPairTrkP);
                histos.fill(HIST("cexbg_pairtrk_pt"), cexPairTrkPt);
                histos.fill(HIST("cexbg_pairtrk_pz"), cexPairTrkPz);
                histos.fill(HIST("cexbg_pairtrkVtxfitR"), radius);
                histos.fill(HIST("cexbg_pairtrkVtxfitDistToPv"), distToPrimary);
                histos.fill(HIST("cexbg_pairtrk_vtxfit_secVtxXY"), secX, secY);
                histos.fill(HIST("cexbg_pairtrk_vtxfit_secVtxZ"), secZ);
              }

              const float chi2 = fitter.getChi2AtPCACandidate();
              histos.fill(HIST("vtxfitChi2"), chi2);
              const double dx = secX - antipVx;
              const double dy = secY - antipVy;
              const double dz = secZ - antipVz;
              const double d3d = std::sqrt(dx * dx + dy * dy + dz * dz);
              histos.fill(HIST("vtxfit_mc_dX"), dx);
              histos.fill(HIST("vtxfit_mc_dY"), dy);
              histos.fill(HIST("vtxfit_mc_dZ"), dz);
              histos.fill(HIST("vtxfit_mc_d3D"), d3d);

              const bool isCex = (motherPdg == -kNeutron);

              const float vtxfitDX = secX - antipVx;
              const float vtxfitDY = secY - antipVy;
              const float vtxfitDZ = secZ - antipVz;
              const float vtxfitD3D = std::sqrt(vtxfitDX * vtxfitDX + vtxfitDY * vtxfitDY + vtxfitDZ * vtxfitDZ);

              const uint32_t selMask = 0u;

              outPairs(
                isCex,
                motherPdg,
                colId,
                pId,
                antipId,

                cexPairMcP,
                cexPairMcPt,
                cexPairMcPz,
                dplane,
                mcangleDeg,
                antipVx,
                antipVy,
                antipVz,

                cexPairTrkP,
                cexPairTrkPt,
                cexPairTrkPz,
                trkangleDeg,
                dcaPair,
                radius,
                distToPrimary,
                secX,
                secY,
                secZ,

                chi2,
                static_cast<int>(status),
                nCand,
                vtxfitDX,
                vtxfitDY,
                vtxfitDZ,
                vtxfitD3D,

                pTrkP,
                pTrkPx,
                pTrkPy,
                pTrkPz,
                pTrkEta,
                pTrkTpcSignal,
                pTrkNClsIts,

                antipTrkP,
                antipTrkPx,
                antipTrkPy,
                antipTrkPz,
                antipTrkEta,
                antipTrkTpcSignal,
                antipTrkNClsIts,

                selMask,

                pairPointingAngleDeg,
                pairPBalance,
                pairPtBalance,
                pairQ,

                dPairP,
                dPairPt,
                dPairPz,
                dOpenAngle,

                svNearestLayerId,
                svDeltaRToLayer,

                pItsMap,
                apItsMap,
                static_cast<int8_t>(pLayers ? 1 : 0),
                static_cast<int8_t>(antipLayers ? 1 : 0),

                pvtxZ,

                pTrkItsNSigmaPr,
                pTrkItsPidValid,
                pTrkTgl,

                antipTrkItsNSigmaPr,
                antipTrkItsPidValid,
                antipTrkTgl);
            }
          }
          // ==== end DCAFitter2 ====
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<NucleiAntineutronCex>(ctx)};
}
