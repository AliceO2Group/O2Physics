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
/// \file resonanceTreeCreator.cxx
/// \brief Produces a TTree with machine learning variables for resonances in the LF group
/// \author Stefano Cannito (stefano.cannito@cern.ch)

#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace resomlcandidates
{
// Multiplicity class
DECLARE_SOA_COLUMN(MultClass, multClass, float); //! Event multiplicity class
// Daughter 1
DECLARE_SOA_COLUMN(PtDaughter1, ptdaughter1, float);     //! Transverse momentum of daughter1 (GeV/c)
DECLARE_SOA_COLUMN(PDaughter1, pdaughter1, float);       //! Momentum of daughter1 (GeV/c)
DECLARE_SOA_COLUMN(PhiDaughter1, phiDaughter1, float);   //! Azimuthal angle of daughter1 (rad)
DECLARE_SOA_COLUMN(EtaDaughter1, etaDaughter1, float);   //! Pseudorapidity of daughter1
DECLARE_SOA_COLUMN(YDaughter1, yDaughter1, float);       //! Rapidity of daughter1
DECLARE_SOA_COLUMN(DCAxyDaughter1, dcaDaughter1, float); //! DCA of daughter1 to primary vertex (cm)
DECLARE_SOA_COLUMN(DCAzDaughter1, dcaZDaughter1, float); //! DCA of daughter1 to primary vertex in z (cm)
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);       //! TPC Nsigma separation for daughter1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);       //! TPC Nsigma separation for daughter1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);       //! TOF Nsigma separation for daughter1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);       //! TOF Nsigma separation for daughter1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPi1, nSigTpcTofPi1, float); //! TPC and TOF combined Nsigma separation for daughter1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKa1, nSigTpcTofKa1, float); //! TPC and TOF combined Nsigma separation for daughter1 with kaon mass hypothesis
// Daughter 2
DECLARE_SOA_COLUMN(PtDaughter2, ptdaughter2, float);     //! Transverse momentum of daughter2 (GeV/c)
DECLARE_SOA_COLUMN(PDaughter2, pdaughter2, float);       //! Momentum of daughter2 (in GeV/c)
DECLARE_SOA_COLUMN(PhiDaughter2, phiDaughter2, float);   //! Azimuthal angle of daughter2 (rad)
DECLARE_SOA_COLUMN(EtaDaughter2, etaDaughter2, float);   //! Pseudorapidity of daughter2
DECLARE_SOA_COLUMN(YDaughter2, yDaughter2, float);       //! Rapidity of daughter2
DECLARE_SOA_COLUMN(DCAxyDaughter2, dcaDaughter2, float); //! DCA of daughter2 to primary vertex (cm)
DECLARE_SOA_COLUMN(DCAzDaughter2, dcaZDaughter2, float); //! DCA of daughter2 to primary vertex in z (cm)
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);       //! TPC Nsigma separation for daughter2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKa2, nSigTpcKa2, float);       //! TPC Nsigma separation for daughter2 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);       //! TOF Nsigma separation for daughter2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKa2, nSigTofKa2, float);       //! TOF Nsigma separation for daughter2 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPi2, nSigTpcTofPi2, float); //! TPC and TOF combined Nsigma separation for daughter2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKa2, nSigTpcTofKa2, float); //! TPC and TOF combined Nsigma separation for daughter2 with kaon mass hypothesis
// Candidate
DECLARE_SOA_COLUMN(M, m, float);                //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);              //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Phi, phi, float);            //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);            //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Y, y, float);                //! Rapidity of candidate
DECLARE_SOA_COLUMN(Sign, sign, int8_t);         //! Sign of the candidate
DECLARE_SOA_COLUMN(IsTruePhi, isTruePhi, bool); //! Flag to indicate if the candidate is a phi meson
} // namespace resomlcandidates

DECLARE_SOA_TABLE(ResoCandidates, "AOD", "RESOCANDIDATES",
                  resomlcandidates::M,
                  resomlcandidates::Pt,
                  resomlcandidates::P,
                  resomlcandidates::Phi,
                  resomlcandidates::Eta,
                  resomlcandidates::Y,
                  resomlcandidates::Sign,
                  resomlcandidates::IsTruePhi);

DECLARE_SOA_TABLE(ResoMLCandidates, "AOD", "RESOMLCANDIDATES",
                  resomlcandidates::MultClass,
                  resomlcandidates::PtDaughter1,
                  resomlcandidates::PDaughter1,
                  resomlcandidates::PhiDaughter1,
                  resomlcandidates::EtaDaughter1,
                  resomlcandidates::YDaughter1,
                  resomlcandidates::DCAxyDaughter1,
                  resomlcandidates::DCAzDaughter1,
                  resomlcandidates::NSigTpcPi1,
                  resomlcandidates::NSigTpcKa1,
                  resomlcandidates::NSigTofPi1,
                  resomlcandidates::NSigTofKa1,
                  resomlcandidates::NSigTpcTofPi1,
                  resomlcandidates::NSigTpcTofKa1,
                  resomlcandidates::PtDaughter2,
                  resomlcandidates::PDaughter2,
                  resomlcandidates::PhiDaughter2,
                  resomlcandidates::EtaDaughter2,
                  resomlcandidates::YDaughter2,
                  resomlcandidates::DCAxyDaughter2,
                  resomlcandidates::DCAzDaughter2,
                  resomlcandidates::NSigTpcPi2,
                  resomlcandidates::NSigTpcKa2,
                  resomlcandidates::NSigTofPi2,
                  resomlcandidates::NSigTofKa2,
                  resomlcandidates::NSigTpcTofPi2,
                  resomlcandidates::NSigTpcTofKa2,
                  resomlcandidates::M,
                  resomlcandidates::Pt,
                  resomlcandidates::P,
                  resomlcandidates::Phi,
                  resomlcandidates::Eta,
                  resomlcandidates::Y,
                  resomlcandidates::Sign);

namespace resomlselection
{
DECLARE_SOA_COLUMN(PhiBDTScore, gammaBDTScore, float);
} // namespace resomlselection

DECLARE_SOA_TABLE(ResoPhiMLSelection, "AOD", "RESOPHIMLSELECTION",
                  resomlselection::PhiBDTScore);
} // namespace o2::aod

struct resonanceTreeCreator {
  // Production of the TTree
  Produces<aod::ResoCandidates> resoCandidates;
  Produces<aod::ResoMLCandidates> resoMLCandidates;

  // Configurables for track selection
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0f, "Cut on the signed transverse momentum to select positive or negative tracks"};

  // Configurables for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 0.5f, "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10.0f, "Maximum pt for the application of the downsampling factor"};

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;

  // Defining the type of the tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi,
                               aod::pidTPCFullKa, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  Partition<FullTracks> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<FullTracks> negTracks = aod::track::signed1Pt < cfgCutCharge;

  Partition<FullMCTracks> posMCTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<FullMCTracks> negMCTracks = aod::track::signed1Pt < cfgCutCharge;

  Preslice<aod::Tracks> perColl = aod::track::collisionId;
  Preslice<aod::McParticles> perMCColl = aod::mcparticle::mcCollisionId;

  // Necessary to flag INEL>0 events in GenMC
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // Cache for manual slicing
  SliceCache cache;

  enum ParticleType {
    Pi,
    Ka,
    Pr
  };

  // Constants
  double massPi = o2::constants::physics::MassPiPlus;
  double massKa = o2::constants::physics::MassKPlus;

  void init(InitContext&)
  {
  }

  // Combine Nsigma values from TPC and TOF
  template <int massHypo, typename T>
  float combineNSigma(const T& track)
  {
    float nSigmaTPC, nSigmaTOF;
    switch (massHypo) {
      case Pi:
        nSigmaTPC = track.tpcNSigmaPi();
        nSigmaTOF = track.tofNSigmaPi();
        break;
      case Ka:
        nSigmaTPC = track.tpcNSigmaKa();
        nSigmaTOF = track.tofNSigmaKa();
        break;
      case Pr:
        nSigmaTPC = track.tpcNSigmaPr();
        nSigmaTOF = track.tofNSigmaPr();
        break;
      default:
        break;
    }

    static constexpr float defaultNSigmaTolerance = .1f;
    static constexpr float defaultNSigma = -999.f + defaultNSigmaTolerance;

    if (nSigmaTPC > defaultNSigma && nSigmaTOF > defaultNSigma)
      return std::sqrt(0.5f * std::pow(nSigmaTPC, 2) + std::pow(nSigmaTOF, 2));
    if (nSigmaTPC > defaultNSigma)
      return std::abs(nSigmaTPC);
    if (nSigmaTOF > defaultNSigma)
      return std::abs(nSigmaTOF);
    return nSigmaTPC;
  }

  // Reconstruct the candidate 4-momentum from two daughter tracks
  template <typename T>
  ROOT::Math::PxPyPzMVector recMother(const T& track1, const T& track2, float masstrack1, float masstrack2)
  {
    ROOT::Math::PxPyPzMVector daughter1(track1.px(), track1.py(), track1.pz(), masstrack1); // set the daughter1 4-momentum
    ROOT::Math::PxPyPzMVector daughter2(track2.px(), track2.py(), track2.pz(), masstrack2); // set the daughter2 4-momentum
    ROOT::Math::PxPyPzMVector mother = daughter1 + daughter2;                               // calculate the mother 4-momentum

    return mother;
  }

  template <typename T1, typename T2>
  void fillCandidateTree4ML(const T1& collision, const T2& track1, const T2& track2, float masstrack1, float masstrack2
                            /*std::optional<std::reference_wrapper<const aod::McParticles>> mcParticles = std::nullopt*/)
  {
    auto tpctofPi1 = combineNSigma<Pi>(track1);
    auto tpctofKa1 = combineNSigma<Ka>(track1);
    auto tpctofPi2 = combineNSigma<Pi>(track2);
    auto tpctofKa2 = combineNSigma<Ka>(track2);

    ROOT::Math::PxPyPzMVector recCandidate = recMother(track1, track2, masstrack1, masstrack2);

    if (downSampleBkgFactor < 1.) {
      float pseudoRndm = track1.pt() * 1000. - static_cast<int64_t>(track1.pt() * 1000);
      if (recCandidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor)
        return;
    }

    resoMLCandidates(collision.centFT0M(),
                     track1.pt(), track1.p(), track1.phi(), track1.eta(), track1.rapidity(masstrack1), track1.dcaXY(), track1.dcaZ(),
                     track1.tpcNSigmaPi(), track1.tpcNSigmaKa(), track1.tofNSigmaPi(), track1.tofNSigmaKa(), tpctofPi1, tpctofKa1,
                     track2.pt(), track2.p(), track2.phi(), track2.eta(), track2.rapidity(masstrack2), track2.dcaXY(), track2.dcaZ(),
                     track2.tpcNSigmaPi(), track2.tpcNSigmaKa(), track2.tofNSigmaPi(), track2.tofNSigmaKa(), tpctofPi2, tpctofKa2,
                     recCandidate.M(), recCandidate.Pt(), recCandidate.P(), recCandidate.Phi(),
                     recCandidate.Eta(), recCandidate.Rapidity(), track1.sign() + track2.sign());
  }

  template <typename T>
  bool isMCPhi(const T& track1, const T& track2, const aod::McParticles& mcParticles)
  {
    if (!track1.has_mcParticle() || !track2.has_mcParticle())
      return false; // Skip filling if no MC truth is available for both tracks

    auto mcTrack1 = mcParticles.rawIteratorAt(track1.mcParticleId());
    auto mcTrack2 = mcParticles.rawIteratorAt(track2.mcParticleId());

    if (mcTrack1.pdgCode() != PDG_t::kKPlus || !mcTrack1.isPhysicalPrimary())
      return false; // Skip filling if the first track is not a primary K+
    if (mcTrack2.pdgCode() != PDG_t::kKMinus || !mcTrack2.isPhysicalPrimary())
      return false; // Skip filling if the second track is not a primary K-

    const auto mcTrack1MotherIndexes = mcTrack1.mothersIds();
    const auto mcTrack2MotherIndexes = mcTrack2.mothersIds();

    for (const auto& mcTrack1MotherIndex : mcTrack1MotherIndexes) {
      for (const auto& mcTrack2MotherIndex : mcTrack2MotherIndexes) {
        if (mcTrack1MotherIndex != mcTrack2MotherIndex)
          continue;

        const auto mother = mcParticles.rawIteratorAt(mcTrack1MotherIndex);
        if (mother.pdgCode() == o2::constants::physics::Pdg::kPhi)
          return true;
      }
    }
    return false;
  }

  void processData4ML(SelCollisions::iterator const& collision, FullTracks const&)
  {
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (const auto& track1 : posThisColl) {
      for (const auto& track2 : negThisColl) {
        // Fill the ResoMLCandidates table with candidates in Data
        fillCandidateTree4ML(collision, track1, track2, massKa, massKa);
      }
    }
  }

  PROCESS_SWITCH(resonanceTreeCreator, processData4ML, "Fill ResoMLCandidates in Data", true);

  void processMC4ML(SimCollisions::iterator const& collision, FullMCTracks const&, aod::McParticles const& mcParticles)
  {
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (const auto& track1 : posThisColl) {
      for (const auto& track2 : negThisColl) {
        if (fillOnlySignal && !isMCPhi(track1, track2, mcParticles))
          return; // Skip filling if only signal is requested and not a phi in MC truth
        if (fillOnlyBackground && isMCPhi(track1, track2, mcParticles))
          return; // Skip filling if only background is requested and a phi in MC truth

        // Fill the ResoMLCandidates table with candidates in MC
        fillCandidateTree4ML(collision, track1, track2, massKa, massKa);
      }
    }
  }

  PROCESS_SWITCH(resonanceTreeCreator, processMC4ML, "Fill ResoMLCandidates in MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<resonanceTreeCreator>(cfgc)};
}
