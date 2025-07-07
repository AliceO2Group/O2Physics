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

/// \file   spectraKinkPiKa.cxx
/// \brief Example of a simple task for the analysis of the muon from Kaon pion using kink topology
/// \author sandeep dudi sandeep.dudi@cern.ch

#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TPDGCode.h"
#include "TVector3.h"
#include <TMath.h>
#include <TString.h>

#include <cstdlib>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCMu>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSel>;
struct spectraKinkPiKa {
  Service<o2::framework::O2DatabasePDG> pdg;
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rpiKkink{"rpiKkink", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaPi{"cutNSigmaPi", 4, "NSigmaTPCPion"};
  Configurable<float> cutNSigmaKa{"cutNSigmaKa", 4, "NSigmaTPCKaon"};
  Configurable<float> rapCut{"rapCut", 0.8, "rapCut"};

  Configurable<int> pid{"pidMother", 321, ""};
  Configurable<bool> d0pid{"dopid", 0, ""};

  Preslice<aod::KinkCands> mPerCol = aod::track::collisionId;

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec ptAxis{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec qtAxis{2000, 0, 2, "#it{q}_{T} (GeV/#it{c})"};
    const AxisSpec kinkAxis{200, 0, 4, "#theta"};
    const AxisSpec etaAxis{200, -5.0, 5.0, "#eta"};
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};

    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    rpiKkink.add("h2_dau_pt_vs_eta_rec", "pt_vs_eta_dau", {HistType::kTH2F, {ptAxis, etaAxis}});
    rpiKkink.add("h2_moth_pt_vs_eta_rec", "pt_vs_eta_moth", {HistType::kTH2F, {ptAxis, etaAxis}});
    rpiKkink.add("h2_pt_moth_vs_dau_rec", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

    rpiKkink.add("h2_qt", "qt", {HistType::kTH1F, {qtAxis}});
    rpiKkink.add("h2_qt_vs_pt", "qt_pt", {HistType::kTH2F, {qtAxis, ptAxis}});

    rpiKkink.add("h2_kink_angle", "kink angle", {HistType::kTH1F, {kinkAxis}});

    // pion
    rpiKkink.add("h2_dau_pt_vs_eta_rec_pion", "pt_vs_eta_dau", {HistType::kTH2F, {ptAxis, etaAxis}});
    rpiKkink.add("h2_moth_pt_vs_eta_rec_pion", "pt_vs_eta_moth", {HistType::kTH2F, {ptAxis, etaAxis}});
    rpiKkink.add("h2_pt_moth_vs_dau_rec_pion", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

    rpiKkink.add("h2_qt_pion", "qt", {HistType::kTH1F, {qtAxis}});
    rpiKkink.add("h2_qt_vs_ptpion", "qt_pt", {HistType::kTH2F, {qtAxis, ptAxis}});
    rpiKkink.add("h2_kink_angle_pion", "kink angle", {HistType::kTH1F, {kinkAxis}});

    if (doprocessMC) {
      rpiKkink.add("h2_dau_pt_vs_eta_gen", "pt_vs_eta_dau", {HistType::kTH2F, {ptAxis, etaAxis}});
      rpiKkink.add("h2_moth_pt_vs_eta_gen", "pt_vs_eta_moth", {HistType::kTH2F, {ptAxis, etaAxis}});
      rpiKkink.add("h2_pt_moth_vs_dau_gen", "pt_moth_vs_dau", {HistType::kTH2F, {ptAxis, ptAxis}});

      rpiKkink.add("h2_qt_gen", "qt", {HistType::kTH1F, {qtAxis}});
      rpiKkink.add("h2_qt_rec", "qt", {HistType::kTH1F, {qtAxis}});
      rpiKkink.add("h2_kink_angle_gen", "kink angle", {HistType::kTH1F, {kinkAxis}});
    }
  }

  void processData(CollisionsFull::iterator const& collision, aod::KinkCands const& KinkCands, TracksFull const&)
  {
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector v1;

    if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8() || !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }

    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    for (const auto& kinkCand : KinkCands) {
      auto dauTrack = kinkCand.trackDaug_as<TracksFull>();
      auto mothTrack = kinkCand.trackMoth_as<TracksFull>();
      bool kaon = false;
      bool pion = false;
      if (std::abs(mothTrack.tpcNSigmaKa()) < cutNSigmaKa) {
        kaon = true;
      }
      if (std::abs(mothTrack.tpcNSigmaPi()) < cutNSigmaPi) {
        pion = true;
      }
      if (!kaon && !pion)
        continue;
      v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassPionCharged);
      v1.SetCoordinates(dauTrack.px(), dauTrack.py(), dauTrack.pz(), o2::constants::physics::MassMuon);

      float pMoth = v0.P();
      float pDaug = v1.P();
      float spKink = mothTrack.px() * dauTrack.px() + mothTrack.py() * dauTrack.py() + mothTrack.pz() * dauTrack.pz();
      float kinkangle = std::acos(spKink / (pMoth * pDaug));
      if (kaon) {
        rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec"), v0.Pt(), v0.Eta());
        rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec"), v1.Pt(), v1.Eta());
        rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec"), v0.Pt(), v1.Pt());
        rpiKkink.fill(HIST("h2_kink_angle"), kinkangle);
      }
      if (pion) {
        rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec_pion"), v0.Pt(), v0.Eta());
        rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec_pion"), v1.Pt(), v1.Eta());
        rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec_pion"), v0.Pt(), v1.Pt());
        rpiKkink.fill(HIST("h2_kink_angle_pion"), kinkangle);
      }
      TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
      // Compute transverse component
      TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
      double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)

      if (kaon) {
        rpiKkink.fill(HIST("h2_qt"), ptd);
        rpiKkink.fill(HIST("h2_qt_vs_pt"), ptd, v1.Pt());
      }
      if (pion) {
        rpiKkink.fill(HIST("h2_qt_pion"), ptd);
        rpiKkink.fill(HIST("h2_qt_vs_ptpion"), ptd, v1.Pt());
      }
    }
  }
  PROCESS_SWITCH(spectraKinkPiKa, processData, "Data processing", true);

  void processMC(CollisionsFullMC const& collisions, aod::KinkCands const& KinkCands, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, TracksFull const&)
  {
    for (const auto& collision : collisions) {
      ROOT::Math::PxPyPzMVector v0;
      ROOT::Math::PxPyPzMVector v1;
      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8() || !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
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
        bool kaon = false;
        bool pion = false;
        if (std::abs(mothTrack.tpcNSigmaKa()) < cutNSigmaKa) {
          kaon = true;
        }
        if (std::abs(mothTrack.tpcNSigmaPi()) < cutNSigmaPi) {
          pion = true;
        }
        if (!kaon && !pion)
          continue;

        v0.SetCoordinates(mothTrack.px(), mothTrack.py(), mothTrack.pz(), o2::constants::physics::MassPionCharged);
        v1.SetCoordinates(dauTrack.px(), dauTrack.py(), dauTrack.pz(), o2::constants::physics::MassMuon);

        float pMoth = v0.P();
        float pDaug = v1.P();
        float spKink = mothTrack.px() * dauTrack.px() + mothTrack.py() * dauTrack.py() + mothTrack.pz() * dauTrack.pz();
        float kinkangle = std::acos(spKink / (pMoth * pDaug));

        rpiKkink.fill(HIST("h2_moth_pt_vs_eta_rec"), v0.Pt(), v0.Eta());
        rpiKkink.fill(HIST("h2_dau_pt_vs_eta_rec"), v1.Pt(), v1.Eta());
        rpiKkink.fill(HIST("h2_pt_moth_vs_dau_rec"), v0.Pt(), v1.Pt());
        rpiKkink.fill(HIST("h2_kink_angle"), kinkangle);

        TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
        // Compute transverse component
        TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
        double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)

        rpiKkink.fill(HIST("h2_qt"), ptd);

        // do MC association
        auto mcLabMoth = trackLabelsMC.rawIteratorAt(mothTrack.globalIndex());
        auto mcLabDau = trackLabelsMC.rawIteratorAt(dauTrack.globalIndex());
        if (mcLabMoth.has_mcParticle() && mcLabDau.has_mcParticle()) {
          auto mcTrackMoth = mcLabMoth.mcParticle_as<aod::McParticles>();
          auto mcTrackDau = mcLabDau.mcParticle_as<aod::McParticles>();
          if (!mcTrackDau.has_mothers()) {
            continue;
          }
          for (const auto& piMother : mcTrackDau.mothers_as<aod::McParticles>()) {
            if (piMother.globalIndex() != mcTrackMoth.globalIndex()) {
              continue;
            }
            if (std::abs(mcTrackMoth.pdgCode()) != pid || std::abs(mcTrackDau.pdgCode()) != kMuonPlus) {
              continue;
            }
            // rpiKkink.fill(HIST("h2MassPtMCRec"), kinkCand.ptMoth(),   v1.Pt());
            rpiKkink.fill(HIST("h2_qt_rec"), ptd);
          }
        }
      }
    }
    for (const auto& mcPart : particlesMC) {
      ROOT::Math::PxPyPzMVector v0;
      ROOT::Math::PxPyPzMVector v1;

      if (!d0pid && (std::abs(mcPart.pdgCode()) != pid || std::abs(mcPart.y()) > rapCut)) {
        continue;
      }
      if (d0pid && (std::abs(mcPart.pdgCode()) != kD0 || std::abs(mcPart.pdgCode()) != kDPlus || std::abs(mcPart.pdgCode()) != kDStar || std::abs(mcPart.y()) > rapCut)) {
        continue;
      }

      if (!mcPart.has_daughters()) {
        continue; // Skip if no daughters
      }
      bool hasKaonpionDaughter = false;
      for (const auto& daughter : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(daughter.pdgCode()) == kMuonPlus) { // muon PDG code
          hasKaonpionDaughter = true;
          v1.SetCoordinates(daughter.px(), daughter.py(), daughter.pz(), o2::constants::physics::MassMuon);
          break; // Found a muon daughter, exit loop
        }
      }
      if (!hasKaonpionDaughter) {
        continue; // Skip if no muon daughter found
      }
      if (pid == kKPlus) {
        v0.SetCoordinates(mcPart.px(), mcPart.py(), mcPart.pz(), o2::constants::physics::MassKaonCharged);
      }

      if (pid == kPiPlus) {
        v0.SetCoordinates(mcPart.px(), mcPart.py(), mcPart.pz(), o2::constants::physics::MassPionCharged);
      }
      if (d0pid) {
        v0.SetCoordinates(mcPart.px(), mcPart.py(), mcPart.pz(), o2::constants::physics::MassD0);
      }

      float pMoth = v0.P();
      float pDaug = v1.P();
      float spKink = v0.Px() * v1.Px() + v0.Py() * v1.Py() + v0.Pz() * v1.Pz();
      float kinkangle = std::acos(spKink / (pMoth * pDaug));

      //  std::cout<< kinkCand.ptMoth()<<" check   "<<v0.Pt()<<std::endl;
      rpiKkink.fill(HIST("h2_moth_pt_vs_eta_gen"), v0.Pt(), v0.Eta());
      rpiKkink.fill(HIST("h2_dau_pt_vs_eta_gen"), v1.Pt(), v1.Eta());
      rpiKkink.fill(HIST("h2_pt_moth_vs_dau_gen"), v0.Pt(), v1.Pt());
      rpiKkink.fill(HIST("h2_kink_angle_gen"), kinkangle);

      TVector3 pdlab(v1.Px(), v1.Py(), v1.Pz());
      // Compute transverse component
      TVector3 motherDir(v0.Px(), v0.Py(), v0.Pz());
      double ptd = pdlab.Perp(motherDir); // or p_d_lab.Mag() * sin(theta)
      rpiKkink.fill(HIST("h2_qt_gen"), ptd);
    }
  }
  PROCESS_SWITCH(spectraKinkPiKa, processMC, "MC processing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<spectraKinkPiKa>(cfgc)};
}
