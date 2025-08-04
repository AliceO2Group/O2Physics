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

/// \file   sigmaminustask.cxx
/// \brief Example of a simple task for the analysis of the Sigma-minus
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                             aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa,
                             aod::pidTOFFullPi, aod::pidTOFFullPr, aod::pidTOFFullKa>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSel>;

struct sigmaminustask {

  // Output Tables
  Produces<aod::SlimKinkCands> outputDataTable;
  Produces<aod::SlimKinkCandsMC> outputDataTableMC;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSigmaMinus{"sigmaminus", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaPi{"cutNSigmaPi", 4, "NSigmaTPCPion"};
  Configurable<float> cutRapMotherMC{"cutRapMotherMC", 1.0f, "Rapidity cut for mother Sigma in MC"};

  Configurable<bool> fillOutputTree{"fillOutputTree", true, "If true, fill the output tree with Kink candidates"};

  Preslice<aod::KinkCands> mPerCol = aod::track::collisionId;

  float radToDeg = o2::constants::math::Rad2Deg;

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec ptAxis{100, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nSigmaPiAxis{60, -30, 30, "n#sigma_{#pi}"};
    const AxisSpec nSigmaPrAxis{60, -30, 30, "n#sigma_{p}"};
    const AxisSpec sigmaMassAxis{100, 1.1, 1.4, "m (GeV/#it{c}^{2})"};
    const AxisSpec xiMassAxis{100, 1.2, 1.6, "m_{#Xi} (GeV/#it{c}^{2})"};
    const AxisSpec pdgAxis{10001, -5000, 5000, "PDG code"};
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};
    const AxisSpec dcaMothAxis{100, 0, 1, "DCA [cm]"};
    const AxisSpec dcaDaugAxis{200, 0, 20, "DCA [cm]"};

    const AxisSpec ptResolutionAxis{100, -0.5, 0.5, "#it{p}_{T}^{rec} - #it{p}_{T}^{gen} (GeV/#it{c})"};
    const AxisSpec massResolutionAxis{100, -0.1, 0.1, "m_{rec} - m_{gen} (GeV/#it{c}^{2})"};

    const AxisSpec boolAxis{2, -0.5, 1.5, "Boolean value (0=false, 1=true)"};
    const AxisSpec filtersAxis{10, -0.5, 9.5, "Filter index"};

    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    // Sigma-minus reconstruction
    rSigmaMinus.add("h2MassSigmaMinusPt", "h2MassSigmaMinusPt", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2SigmaMassVsXiMass", "h2SigmaMassVsXiMass", {HistType::kTH2F, {xiMassAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2NSigmaTPCPiPt", "h2NSigmaTPCPiPt", {HistType::kTH2F, {ptAxis, nSigmaPiAxis}});
    rSigmaMinus.add("h2DCAMothPt", "h2DCAMothPt", {HistType::kTH2F, {ptAxis, dcaMothAxis}});
    rSigmaMinus.add("h2DCADaugPt", "h2DCADaugPt", {HistType::kTH2F, {ptAxis, dcaDaugAxis}});

    if (doprocessMC) {
      // Add MC histograms if needed
      rSigmaMinus.add("h2MassPtMCRec", "h2MassPtMCRec", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
      rSigmaMinus.add("h2MassPtMCGen", "h2MassPtMCGen", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});

      rSigmaMinus.add("h2MassResolution", "h2MassResolution", {HistType::kTH2F, {ptAxis, massResolutionAxis}});
      rSigmaMinus.add("h2PtResolution", "h2PtResolution", {HistType::kTH2F, {ptAxis, ptResolutionAxis}});
      
      rSigmaMinus.add("h2NSigmaTOFPiPt", "h2NSigmaTOFPiPt", {HistType::kTH2F, {ptAxis, nSigmaPiAxis}});
      rSigmaMinus.add("h2NSigmaTOFPrPt", "h2NSigmaTOFPrPt", {HistType::kTH2F, {ptAxis, nSigmaPrAxis}});

      // BC ID comparison histograms
      rSigmaMinus.add("hMcCollIdCoherence", "McCollId (coll == daug)", {HistType::kTH1F, {boolAxis}});
      rSigmaMinus.add("h2CollId_BCId", "(McCollId coherence) vs (EvSelBC == McBC)", {HistType::kTH2F, {boolAxis, boolAxis}});
      rSigmaMinus.add("h2BCId_comp1", "(BC == McBC) vs (BC == EvSelBC)", {HistType::kTH2F, {boolAxis, boolAxis}});
      rSigmaMinus.add("h2BCId_comp2", "(McBC == EvSelBC) vs (BC == EvSelBC)", {HistType::kTH2F, {boolAxis, boolAxis}});
    }

    if (doprocessFindable) {
      // Add findable Sigma histograms
      rSigmaMinus.add("h2MassPtFindable", "h2MassPtFindable", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
      rSigmaMinus.add("hFilterIndex", "hFilterIndex", {HistType::kTH1F, {filtersAxis}});
    }
  }

  void processData(CollisionsFull::iterator const& collision, aod::KinkCands const& KinkCands, TracksFull const&)
  {
    if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
      return;
    }
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    for (const auto& kinkCand : KinkCands) {
      auto dauTrack = kinkCand.trackDaug_as<TracksFull>();

      if (std::abs(dauTrack.tpcNSigmaPi()) > cutNSigmaPi) {
        continue;
      }

      rSigmaMinus.fill(HIST("h2MassSigmaMinusPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2SigmaMassVsXiMass"), kinkCand.mXiMinus(), kinkCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2NSigmaPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tpcNSigmaPi());
      rSigmaMinus.fill(HIST("h2DCAMothPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.dcaMothPv());
      rSigmaMinus.fill(HIST("h2DCADaugPt"), kinkCand.mothSign() * kinkCand.ptDaug(), kinkCand.dcaDaugPv());

      if (fillOutputTree) {
        outputDataTable(kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx(),
                        kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth(),
                        kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug(),
                        kinkCand.dcaMothPv(), kinkCand.dcaDaugPv(), kinkCand.dcaKinkTopo(),
                        kinkCand.mothSign(),
                        dauTrack.tpcNSigmaPi(), dauTrack.tpcNSigmaPr(), dauTrack.tpcNSigmaKa(),
                        dauTrack.tofNSigmaPi(), dauTrack.tofNSigmaPr(), dauTrack.tofNSigmaKa());
      }
    }
  }
  PROCESS_SWITCH(sigmaminustask, processData, "Data processing", true);

  void processMC(CollisionsFullMC const& collisions, aod::KinkCands const& KinkCands, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC,  TracksFull const&, aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
        continue;
      }

      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      auto kinkCandPerColl = KinkCands.sliceBy(mPerCol, collision.globalIndex());

      for (const auto& kinkCand : kinkCandPerColl) {

        auto dauTrack = kinkCand.trackDaug_as<TracksFull>();
        auto mothTrack = kinkCand.trackMoth_as<TracksFull>();
        if (dauTrack.sign() != mothTrack.sign()) {
          LOG(info) << "Skipping kink candidate with opposite sign daughter and mother: " << kinkCand.globalIndex();
          continue; // Skip if the daughter has the opposite sign as the mother
        }
        if (std::abs(dauTrack.tpcNSigmaPi()) > cutNSigmaPi) {
          continue;
        }

        // histograms filled with all kink candidates
        rSigmaMinus.fill(HIST("h2MassSigmaMinusPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
        rSigmaMinus.fill(HIST("h2SigmaMassVsXiMass"), kinkCand.mXiMinus(), kinkCand.mSigmaMinus());
        rSigmaMinus.fill(HIST("h2NSigmaTPCPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tpcNSigmaPi());

        // do MC association
        auto mcLabSigma = trackLabelsMC.rawIteratorAt(mothTrack.globalIndex());
        auto mcLabPiDau = trackLabelsMC.rawIteratorAt(dauTrack.globalIndex());
        if (mcLabSigma.has_mcParticle() && mcLabPiDau.has_mcParticle()) {
          auto mcTrackSigma = mcLabSigma.mcParticle_as<aod::McParticles>();
          auto mcTrackPiDau = mcLabPiDau.mcParticle_as<aod::McParticles>();
          if (!mcTrackPiDau.has_mothers()) {
            continue;
          }
          for (auto& piMother : mcTrackPiDau.mothers_as<aod::McParticles>()) {
            if (piMother.globalIndex() != mcTrackSigma.globalIndex()) {
              continue;
            }
            if (std::abs(mcTrackSigma.pdgCode()) != 3112 && std::abs(mcTrackSigma.pdgCode()) != 3222) {
              continue;
            }
            if (std::abs(mcTrackPiDau.pdgCode()) != 211 && std::abs(mcTrackPiDau.pdgCode()) != 2212) {
              continue;
            }

            float MotherMassMC = std::sqrt(piMother.e() * piMother.e() - piMother.p() * piMother.p());
            float MotherpTMC = piMother.pt();
            float deltaXMother = mcTrackPiDau.vx() - piMother.vx();
            float deltaYMother = mcTrackPiDau.vy() - piMother.vy();
            float decayRadiusMC = std::sqrt(deltaXMother * deltaXMother + deltaYMother * deltaYMother);

            // Check coherence of MCcollision Id for daughter MCparticle and reconstructed collision
            bool mcCollisionIdCheck = false;
            if (collision.has_mcCollision()) {
              mcCollisionIdCheck = collision.mcCollision().globalIndex() == mcTrackPiDau.mcCollisionId();
            }
            // Check bunch crossing ID coherence
            auto mcCollision = mcTrackPiDau.template mcCollision_as<aod::McCollisions>();
            bool BCId_vs_MCBCId = collision.bcId() == mcCollision.bcId();
            bool BCId_vs_EvSel = collision.bcId() == collision.foundBCId();
            bool EvSel_vs_MCBCId = collision.foundBCId() == mcCollision.bcId();

            rSigmaMinus.fill(HIST("hMcCollIdCoherence"), static_cast<int>(mcCollisionIdCheck));
            rSigmaMinus.fill(HIST("h2CollId_BCId"), static_cast<int>(mcCollisionIdCheck), static_cast<int>(EvSel_vs_MCBCId));
            rSigmaMinus.fill(HIST("h2BCId_comp1"), static_cast<int>(BCId_vs_MCBCId), static_cast<int>(BCId_vs_EvSel));
            rSigmaMinus.fill(HIST("h2BCId_comp2"), static_cast<int>(EvSel_vs_MCBCId), static_cast<int>(BCId_vs_EvSel));

            rSigmaMinus.fill(HIST("h2MassPtMCRec"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
            rSigmaMinus.fill(HIST("h2MassResolution"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus() - MotherMassMC);
            rSigmaMinus.fill(HIST("h2PtResolution"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.ptMoth() - MotherpTMC);
            rSigmaMinus.fill(HIST("h2DCAMothPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.dcaMothPv());
            rSigmaMinus.fill(HIST("h2DCADaugPt"), kinkCand.mothSign() * kinkCand.ptDaug(), kinkCand.dcaDaugPv());

            if (std::abs(mcTrackPiDau.pdgCode()) == 211) {
              rSigmaMinus.fill(HIST("h2NSigmaTOFPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tofNSigmaPi());
            } else if (std::abs(mcTrackPiDau.pdgCode()) == 2212) {
              rSigmaMinus.fill(HIST("h2NSigmaTOFPrPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tofNSigmaPr());
            }

            // fill the output table with Mc information
            if (fillOutputTree) {
              outputDataTableMC(kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx(),
                                kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth(),
                                kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug(),
                                kinkCand.dcaMothPv(), kinkCand.dcaDaugPv(), kinkCand.dcaKinkTopo(),
                                kinkCand.mothSign(),
                                dauTrack.tpcNSigmaPi(), dauTrack.tpcNSigmaPr(), dauTrack.tpcNSigmaKa(),
                                dauTrack.tofNSigmaPi(), dauTrack.tofNSigmaPr(), dauTrack.tofNSigmaKa(),
                                mcTrackSigma.pdgCode(), mcTrackPiDau.pdgCode(),
                                MotherpTMC, MotherMassMC, decayRadiusMC, mcCollisionIdCheck);
            }
          }
        } // MC association and selection
      } // kink cand loop
    } // collision loop

    // Loop over all generated particles to fill MC histograms
    for (const auto& mcPart : particlesMC) {
      if ((std::abs(mcPart.pdgCode()) != 3112 && std::abs(mcPart.pdgCode()) != 3222) || std::abs(mcPart.y()) > cutRapMotherMC) { // only sigma mothers and rapidity cut
        continue;
      }
      if (!mcPart.has_daughters()) {
        continue; // Skip if no daughters
      }
      bool hasSigmaDaughter = false;
      int daug_pdg = 0;
      std::array<float, 3> secVtx;
      std::array<float, 3> momDaug;
      for (const auto& daughter : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(daughter.pdgCode()) == 211 || std::abs(daughter.pdgCode()) == 2212) { // Pi or proton daughter
          hasSigmaDaughter = true;
          secVtx = {daughter.vx(), daughter.vy(), daughter.vz()};
          momDaug = {daughter.px(), daughter.py(), daughter.pz()};
          daug_pdg = daughter.pdgCode();
          break; // Found a daughter, exit loop
        }
      }
      if (!hasSigmaDaughter) {
        continue; // Skip if no pi/proton daughter found
      }
      float mcMass = std::sqrt(mcPart.e() * mcPart.e() - mcPart.p() * mcPart.p());
      float mcDecayRadius = std::sqrt((secVtx[0] - mcPart.vx()) * (secVtx[0] - mcPart.vx()) + (secVtx[1] - mcPart.vy()) * (secVtx[1] - mcPart.vy()));
      int sigmaSign = mcPart.pdgCode() > 0 ? 1 : -1; // Determine the sign of the Sigma
      rSigmaMinus.fill(HIST("h2MassPtMCGen"), sigmaSign * mcPart.pt(), mcMass);

      // Fill output table with non reconstructed MC candidates
      if (fillOutputTree) {
        outputDataTableMC(-999, -999, -999,
                          -999, -999, -999,
                          -999, -999, -999,
                          -999, -999, -999,
                          sigmaSign,
                          -999, -999, -999,
                          -999, -999, -999,
                          mcPart.pdgCode(), daug_pdg,
                          mcPart.pt(), mcMass, mcDecayRadius, false);
      }
    }
  }

  PROCESS_SWITCH(sigmaminustask, processMC, "MC processing", false);

  void processFindable(TracksFull const& tracks, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const&, CollisionsFullMC const&)
  {
    for (const auto& motherTrack : tracks) {
      // Check if mother is Sigma in MC
      auto mcLabelMother = trackLabelsMC.rawIteratorAt(motherTrack.globalIndex());
      if (!mcLabelMother.has_mcParticle()) {
        continue;
      }
      auto mcMother = mcLabelMother.mcParticle_as<aod::McParticles>();
      if (std::abs(mcMother.pdgCode()) != 3112 && std::abs(mcMother.pdgCode()) != 3222) {
        continue;
      }

      for (const auto& daughterTrack : tracks) {
        // Check if daughter is pi/proton
        auto mcLabelDaughter = trackLabelsMC.rawIteratorAt(daughterTrack.globalIndex());
        if (!mcLabelDaughter.has_mcParticle()) {
          continue;
        }
        auto mcDaughter = mcLabelDaughter.mcParticle_as<aod::McParticles>();
        if (std::abs(mcDaughter.pdgCode()) != 211 && std::abs(mcDaughter.pdgCode()) != 2212) {
          continue;
        }
        
        // Verify the MC mother-daughter relationship
        bool isValidPair = false;
        if (mcDaughter.has_mothers()) {
          for (const auto& mother : mcDaughter.mothers_as<aod::McParticles>()) {
            if (mother.globalIndex() == mcMother.globalIndex()) {
              isValidPair = true;
              break;
            }
          }
        }
        if (!isValidPair) {
          continue; 
        }

        float mcMass = std::sqrt(mcMother.e() * mcMother.e() - mcMother.p() * mcMother.p());
        int sigmaSign = mcMother.pdgCode() > 0 ? 1 : -1;
        rSigmaMinus.fill(HIST("h2MassPtFindable"), sigmaSign * motherTrack.pt(), mcMass);

        // Define filter index and progressively apply kinkbuilder cuts to track pairs
        int filterIndex = 0;
        rSigmaMinus.fill(HIST("hFilterIndex"), filterIndex);

        // 1 - tracks with right ITS, TPC, TOF signals
        if (motherTrack.has_collision() && motherTrack.hasITS() && !motherTrack.hasTPC() && !motherTrack.hasTOF() &&
            daughterTrack.hasITS() && daughterTrack.hasTPC()) {
          filterIndex += 1;
          rSigmaMinus.fill(HIST("hFilterIndex"), filterIndex);
        } else {
          continue; 
        }
        
        // 2, 3 - mother track ITS properties
        if (motherTrack.itsNCls() < 6 &&
            motherTrack.itsNClsInnerBarrel() == 3 && motherTrack.itsChi2NCl() < 36) {
          filterIndex += 1;
          rSigmaMinus.fill(HIST("hFilterIndex"), filterIndex);
        }
        if (filterIndex > 1 && motherTrack.pt() > 0.5) {
          filterIndex += 1;
          rSigmaMinus.fill(HIST("hFilterIndex"), filterIndex);
        } else {
          continue; 
        }

        // 4 - daughter track ITS+TPC properties
        if (daughterTrack.itsNClsInnerBarrel() == 0 && daughterTrack.itsNCls() < 4 &&
            daughterTrack.tpcNClsCrossedRows() > 0.8 * daughterTrack.tpcNClsFindable() && daughterTrack.tpcNClsFound() > 80) {
          filterIndex += 1;
          rSigmaMinus.fill(HIST("hFilterIndex"), filterIndex);
        } else {
          continue; 
        }
        
        // 5 - geometric cuts: eta
        if (std::abs(motherTrack.eta()) < 1.0 && std::abs(daughterTrack.eta()) < 1.0) {
          filterIndex += 1;
          rSigmaMinus.fill(HIST("hFilterIndex"), filterIndex);
        } else {
          continue; 
        }
        // 6 - geometric cuts: phi
        if (std::abs(motherTrack.phi() - daughterTrack.phi()) * radToDeg < 100.0) {
          filterIndex += 1;
          rSigmaMinus.fill(HIST("hFilterIndex"), filterIndex);
        } else {
          continue; 
        }

        // 7 - collision selection
        auto collision = motherTrack.template collision_as<CollisionsFullMC>();
        if (!(std::abs(collision.posZ()) > cutzvertex || !collision.sel8())) {
          filterIndex += 1;
          rSigmaMinus.fill(HIST("hFilterIndex"), filterIndex);
        }
        
      }
    }
  }

  PROCESS_SWITCH(sigmaminustask, processFindable, "Findable Sigma processing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<sigmaminustask>(cfgc)};
}
