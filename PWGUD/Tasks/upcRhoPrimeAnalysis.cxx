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
/// \file   upcRhoPrimeAnalysis.cxx
/// \brief  Task for analysis of rho prime in UPCs using UD tables (from SG producer).
/// \author Cesar Omar Ramirez Alvarez (cesar.ramirez@cern.ch), Autonomous University of Puebla

#include "PWGUD/DataModel/UDTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace ROOT::Math;

// Define UD tables
using UDtracks = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
using UDCollisions = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>;
using UDMcCollisions = aod::UDMcCollisions;
using UDMcParticles = aod::UDMcParticles;
using UDtracksMC = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels>;
using UDCollisionsMC = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDMcCollsLabels>;

namespace o2::aod
{
namespace fourpi
{
// Columns RD / MC reco
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);                        // Run number
DECLARE_SOA_COLUMN(M, m, double);                                         // System invariant mass
DECLARE_SOA_COLUMN(Pt, pt, double);                                       // System pT
DECLARE_SOA_COLUMN(Eta, eta, double);                                     // System pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, double);                                     // System azimuthal angle
DECLARE_SOA_COLUMN(PosX, posX, double);                                   // Vertex X position
DECLARE_SOA_COLUMN(PosY, posY, double);                                   // Vertex Y position
DECLARE_SOA_COLUMN(PosZ, posZ, double);                                   // Vertex Z position
DECLARE_SOA_COLUMN(TotalCharge, totalCharge, int);                        // Real total charge of the 4 selected tracks
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);        // FT0A amplitude
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);        // FT0C amplitude
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);        // FV0A amplitude
DECLARE_SOA_COLUMN(NumContrib, numContrib, int32_t);                      // Number of PV contributors
DECLARE_SOA_COLUMN(Sign, sign, std::vector<int>);                         // Track charges
DECLARE_SOA_COLUMN(TrackPt, trackPt, std::vector<float>);                 // Track pT
DECLARE_SOA_COLUMN(TrackEta, trackEta, std::vector<float>);               // Track eta
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, std::vector<float>);               // Track phi
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, std::vector<float>);         // TPC nSigma electron
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, std::vector<float>);         // TPC nSigma pion
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, std::vector<float>);         // TPC nSigma kaon
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, std::vector<float>);         // TPC nSigma proton
DECLARE_SOA_COLUMN(TrackID, trackID, std::vector<int>);                   // Track index within the system
DECLARE_SOA_COLUMN(IsReconstructedWithUPC, isReconstructedWithUPC, bool); // UPC reconstruction mode flag
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);                              // ZNA time
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);                              // ZNC time
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);              // ZNA energy
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);              // ZNC energy
DECLARE_SOA_COLUMN(IsChargeZero, isChargeZero, bool);                     // TotalCharge == 0
DECLARE_SOA_COLUMN(OccupancyInTime, occupancyInTime, int);                // Occupancy
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, double);                   // Hadronic interaction rate
} // namespace fourpi

DECLARE_SOA_TABLE(SYSTEMTREE, "AOD", "SystemTree",
                  fourpi::RunNumber, fourpi::M, fourpi::Pt, fourpi::Eta, fourpi::Phi,
                  fourpi::PosX, fourpi::PosY, fourpi::PosZ, fourpi::TotalCharge,
                  fourpi::TotalFT0AmplitudeA, fourpi::TotalFT0AmplitudeC, fourpi::TotalFV0AmplitudeA,
                  fourpi::NumContrib,
                  fourpi::Sign, fourpi::TrackPt, fourpi::TrackEta, fourpi::TrackPhi,
                  fourpi::TPCNSigmaEl, fourpi::TPCNSigmaPi, fourpi::TPCNSigmaKa, fourpi::TPCNSigmaPr,
                  fourpi::TrackID, fourpi::IsReconstructedWithUPC,
                  fourpi::TimeZNA, fourpi::TimeZNC, fourpi::EnergyCommonZNA, fourpi::EnergyCommonZNC,
                  fourpi::IsChargeZero, fourpi::OccupancyInTime, fourpi::HadronicRate);

namespace mcgen4pi
{
// Columns MC gen
DECLARE_SOA_COLUMN(McMotherPdg, mcMotherPdg, int);             // Common mother PDG
DECLARE_SOA_COLUMN(McMotherPt, mcMotherPt, float);             // Mother pT
DECLARE_SOA_COLUMN(McMotherPhi, mcMotherPhi, float);           // Mother phi
DECLARE_SOA_COLUMN(McMotherMass, mcMotherMass, float);         // Mother invariant mass
DECLARE_SOA_COLUMN(McMotherRapidity, mcMotherRapidity, float); // Mother rapidity
DECLARE_SOA_COLUMN(McTotalCharge, mcTotalCharge, int);         // Real total charge of the 4 generated particles

// Per-particle info
DECLARE_SOA_COLUMN(McTrackPdg, mcTrackPdg, int[4]);
DECLARE_SOA_COLUMN(McTrackPt, mcTrackPt, float[4]);
DECLARE_SOA_COLUMN(McTrackEta, mcTrackEta, float[4]);
DECLARE_SOA_COLUMN(McTrackPhi, mcTrackPhi, float[4]);
DECLARE_SOA_COLUMN(McTrackSign, mcTrackSign, int[4]);
DECLARE_SOA_COLUMN(McTrackIsPrimary, mcTrackIsPrimary, int[4]);

// Generated vertex
DECLARE_SOA_COLUMN(McPosX, mcPosX, float);
DECLARE_SOA_COLUMN(McPosY, mcPosY, float);
DECLARE_SOA_COLUMN(McPosZ, mcPosZ, float);

// Generated collision index
DECLARE_SOA_COLUMN(McCollisionIndex, mcCollisionIndex, int);

// Real run number for MC gen
DECLARE_SOA_COLUMN(McRunNumber, mcRunNumber, int);

DECLARE_SOA_COLUMN(RecoIndex, recoIndex, int);
} // namespace mcgen4pi

// MC reco<->gen match table
DECLARE_SOA_TABLE(FourPiMcMatchTree, "AOD", "FOURPIMCMATCH",
                  fourpi::IsReconstructedWithUPC,
                  mcgen4pi::McMotherPdg, mcgen4pi::McMotherPt, mcgen4pi::McMotherPhi,
                  mcgen4pi::McMotherMass, mcgen4pi::McMotherRapidity, mcgen4pi::McTotalCharge,
                  mcgen4pi::McTrackPdg, mcgen4pi::McTrackPt, mcgen4pi::McTrackEta, mcgen4pi::McTrackPhi,
                  mcgen4pi::McTrackSign, mcgen4pi::McTrackIsPrimary,
                  mcgen4pi::McPosX, mcgen4pi::McPosY, mcgen4pi::McPosZ, mcgen4pi::McCollisionIndex,
                  mcgen4pi::McRunNumber, mcgen4pi::RecoIndex);

// MC All generated collisions
DECLARE_SOA_TABLE(FourPiMcGenAllTree, "AOD", "FOURPIMCGALL",
                  mcgen4pi::McMotherPdg, mcgen4pi::McMotherPt, mcgen4pi::McMotherPhi,
                  mcgen4pi::McMotherMass, mcgen4pi::McMotherRapidity, mcgen4pi::McTotalCharge,
                  mcgen4pi::McTrackPdg, mcgen4pi::McTrackPt, mcgen4pi::McTrackEta, mcgen4pi::McTrackPhi,
                  mcgen4pi::McTrackSign, mcgen4pi::McTrackIsPrimary,
                  mcgen4pi::McPosX, mcgen4pi::McPosY, mcgen4pi::McPosZ, mcgen4pi::McCollisionIndex,
                  mcgen4pi::McRunNumber);
} // namespace o2::aod

struct upcRhoPrimeAnalysis {
  Produces<aod::SYSTEMTREE> systemTree;
  Produces<aod::FourPiMcMatchTree> fourPiMcMatchTree;
  Produces<aod::FourPiMcGenAllTree> fourPiMcGenAllTree;

  // System selection configuration
  Configurable<double> systemYCut{"systemYCut", 0.5, "Max Rapidity of rho prime"};
  Configurable<double> systemPtCut{"systemPtCut", 0.1, "Min Pt of rho prime"};
  Configurable<double> systemMassMinCut{"systemMassMinCut", 0.8, "Min Mass of rho prime"};
  Configurable<double> systemMassMaxCut{"systemMassMaxCut", 2.2, "Max Mass of rho prime"};
  Configurable<double> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};

  // Event selection configuration
  Configurable<float> vZCut{"vZCut", 10.0, "Cut on vertex Z position"};
  Configurable<int> numPVContrib{"numPVContrib", 4, "Number of PV contributors"};
  Configurable<float> fv0Cut{"fv0Cut", 50.0, "FV0 amplitude cut"};
  Configurable<float> ft0aCut{"ft0aCut", 50.0, "FT0A amplitude cut"};
  Configurable<float> ft0cCut{"ft0cCut", 50.0, "FT0C amplitude cut"};
  Configurable<float> zdcCut{"zdcCut", 0.0, "ZDC energy cut"};
  Configurable<bool> sbpCut{"sbpCut", true, "SBP cut"};
  Configurable<bool> itsROFbCut{"itsROFbCut", true, "ITS ROFb cut"};
  Configurable<bool> vtxITSTPCcut{"vtxITSTPCcut", true, "Vertex ITS-TPC cut"};
  Configurable<bool> tfbCut{"tfbCut", true, "TFB cut"};
  Configurable<bool> specifyGapSide{"specifyGapSide", true, "specify gap side for SG/DG produced data"};
  Configurable<int> gapSide{"gapSide", 2, "gap side for SG produced data"};

  // Track selection configuration
  Configurable<bool> useOnlyPVtracks{"useOnlyPVtracks", true, "Use only PV tracks"};
  Configurable<float> tpcChi2NClsCut{"tpcChi2NClsCut", 5.0, "TPC chi2/N clusters cut"};
  Configurable<float> itsChi2NClsCut{"itsChi2NClsCut", 36.0, "ITS chi2/N clusters cut"};
  Configurable<float> nSigmaTPCcut{"nSigmaTPCcut", 5.0, "TPC nSigma cut"};
  Configurable<float> dcaXYcut{"dcaXYcut", 0, "dcaXY cut"};
  Configurable<float> dcaZcut{"dcaZcut", 2, "dcaZ cut"};
  Configurable<int> minTPCFindableClusters{"minTPCFindableClusters", 70, "Minimum number of findable TPC clusters"};

  // Optional generatorId filter
  Configurable<int> genId{"genId", -1, "generator ID; -1 = no filter"};

  // Define histogram registry for RD / MC reco
  HistogramRegistry registry{
    "registry",
    {// Event flow histograms
     {"Events/Flow", "Event flow;Cut;Counts", {HistType::kTH1D, {{11, 0, 11}}}},
     {"Events/FlowDetailed", "Detailed event flow;Cut;Counts", {HistType::kTH1D, {{13, 0, 13}}}},
     {"Events/hRecoMode", "Reconstruction mode;;Counts", {HistType::kTH1D, {{2, 0, 2}}}},
     {"Events/VertexZ", "Vertex Z;z (cm);Counts", {HistType::kTH1F, {{200, -20, 20}}}},
     {"Events/NumContrib", "Number of contributors;N_{contrib};Counts", {HistType::kTH1F, {{100, 0, 100}}}},
     {"Events/FV0Amplitude", "FV0 amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
     {"Events/FT0AmplitudeA", "FT0A amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
     {"Events/FT0AmplitudeC", "FT0C amplitude;Amplitude;Counts", {HistType::kTH1F, {{200, 0, 200}}}},
     {"Events/ZDCEnergy", "ZDC energy;Energy (TeV);Counts", {HistType::kTH1F, {{200, 0, 2}}}},

     // Track quality histograms
     {"Tracks/Pt", "Track p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{200, 0, 2}}}},
     {"Tracks/Eta", "Track #eta;#eta;Counts", {HistType::kTH1F, {{200, -2, 2}}}},
     {"Tracks/TPCNSigmaPi", "TPC n#sigma for #pi;n#sigma;Counts", {HistType::kTH1F, {{200, -10, 10}}}},
     {"Tracks/TPCChi2NCl", "TPC #chi^{2}/N_{cls};#chi^{2}/N_{cls};Counts", {HistType::kTH1F, {{200, 0, 20}}}},
     {"Tracks/ITSChi2NCl", "ITS #chi^{2}/N_{cls};#chi^{2}/N_{cls};Counts", {HistType::kTH1F, {{200, 0, 50}}}},
     {"Tracks/RejectionReasons", "Track rejection reasons;Reason;Counts", {HistType::kTH1F, {{15, 0, 15}}}},
     {"Tracks/DCASpectrum", "Track DCA spectrum;DCA (cm);Counts", {HistType::kTH1F, {{100, 0, 5}}}},
     {"Tracks/ChargeDistribution", "Track charge distribution;Charge;Counts", {HistType::kTH1F, {{3, -1.5, 1.5}}}},
     {"Tracks/TPCClusters", "TPC clusters findable;N_{clusters};Counts", {HistType::kTH1F, {{100, 0, 200}}}},
     {"Tracks/NGoodTracksPerEvent", "Good tracks per event (before the ==4 cut);N;Counts", {HistType::kTH1F, {{11, -0.5, 10.5}}}},

     // System kinematics histograms
     {"System/hM", ";m (GeV/#it{c}^{2});counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"System/hPt", ";p_{T} (GeV/#it{c});counts", {HistType::kTH1F, {{1000, 0.0, 1.1}}}},
     {"System/hEta", ";#eta;counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
     {"System/hPhi", ";#phi;counts", {HistType::kTH1F, {{180, 0.0, 6.28}}}},
     {"System/hY", ";y;counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
     {"System/hTotalChargeBefore", "Total charge before M/Pt/Y cuts;Q;counts", {HistType::kTH1F, {{9, -4.5, 4.5}}}},
     // -4=----, -2=---+ (3-,1+), 0=+-+- (2+,2-), +2=+++- (3+,1-), +4=++++
     {"System/hTotalCharge", "System total charge (4 tracks, after M/Pt/Y cuts);Q;counts", {HistType::kTH1F, {{9, -4.5, 4.5}}}},
     {"System/hMVsTotalChargeBefore", "Invariant mass vs charge combination (before M/Pt/Y cuts);m (GeV/#it{c}^{2});Q", {HistType::kTH2F, {{1000, 0.0, 10.0}, {9, -4.5, 4.5}}}},
     {"System/hMVsTotalCharge", "Invariant mass vs charge combination (after M/Pt/Y cuts);m (GeV/#it{c}^{2});Q", {HistType::kTH2F, {{1000, 0.0, 10.0}, {9, -4.5, 4.5}}}},

     // Comparison histograms
     {"Cuts/MBefore", "Mass before cuts;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0, 10}}}},
     {"Cuts/MAfter", "Mass after cuts;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0, 10}}}},
     {"Cuts/PtBefore", "p_{T} before cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0, 1.1}}}},
     {"Cuts/PtAfter", "p_{T} after cuts;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0, 1.1}}}}}};

  // Define histogram registry for MC analysis
  HistogramRegistry mcRegistry{
    "mcRegistry",
    {// Event-level MC histograms
     {"MC/Events/hAllEvents", "All MC Events", {HistType::kTH1F, {{1, 0, 1}}}},
     {"MC/Events/hAccepted", "Accepted MC Events", {HistType::kTH1F, {{1, 0, 1}}}},
     {"MC/Events/hVertexZ", "MC Vertex Z;z (cm);Counts", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
     {"MC/Events/hNPrimaries", "Number of primary particles;N;Counts", {HistType::kTH1F, {{10, -0.5, 9.5}}}},

     // Track-level MC histograms
     {"MC/Tracks/hPt", "MC Track p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{200, 0.0, 2.0}}}},
     {"MC/Tracks/hEta", "MC Track #eta;#eta;Counts", {HistType::kTH1F, {{200, -2.0, 2.0}}}},
     {"MC/Tracks/hPhi", "MC Track #phi;#phi;Counts", {HistType::kTH1F, {{200, 0.0, 6.28}}}},
     {"MC/Tracks/hPdgCode", "PDG codes;PDG code;Counts", {HistType::kTH1F, {{2000, -1000, 1000}}}},

     // System-level (mother) MC histograms
     {"MC/System/hM", "MC Invariant Mass;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/System/hPt", "MC p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 1.1}}}},
     {"MC/System/hY", "MC Rapidity;y;Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
     {"MC/System/hMvsPt", "MC Mass vs p_{T};m (GeV/c^{2});p_{T} (GeV/c)", {HistType::kTH2F, {{1000, 0.0, 10.0}, {1000, 0.0, 1.1}}}},
     {"MC/System/hMvsY", "MC Mass vs Rapidity;m (GeV/c^{2});y", {HistType::kTH2F, {{1000, 0.0, 10.0}, {200, -2.0, 2.0}}}},
     {"MC/System/hTotalCharge", "Generated total charge;Q;Counts", {HistType::kTH1F, {{9, -4.5, 4.5}}}},

     // Basic summary
     {"MC/Summary/hEventCounter", "Generated vs reconstructed summary;;Counts", {HistType::kTH1D, {{2, 0, 2}}}},
     {"MC/Summary/hMatchStatus", "Reco<->Gen match status;;Counts", {HistType::kTH1D, {{3, 0, 3}}}},
     {"MC/Summary/hRecoMode", "Reconstruction mode;;Counts", {HistType::kTH1D, {{2, 0, 2}}}},
     {"MC/Summary/hRecoTotalCharge", "Reco total charge (MC events, after M/Pt/Y cuts);Q;Counts", {HistType::kTH1D, {{9, -4.5, 4.5}}}},

     // Rho prime specific histograms
     {"MC/RhoPrime/hFound", "Rho Prime Found;Found;Counts", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
     {"MC/RhoPrime/hMass", "Rho Prime Mass;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/RhoPrime/hMassUPC", "Rho Prime Mass UPC;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/RhoPrime/hMassSTD", "Rho Prime Mass STD;m (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/RhoPrime/hPt", "Rho Prime p_{T};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 1.1}}}},

     // Control histograms
     {"MC/Control/hNDaughtersOfMother", "Number of daughters of the common mother;N;Counts", {HistType::kTH1F, {{10, 0, 10}}}},
     {"MC/Control/hMotherPdg", "PDG of the common mother found;PDG;Counts", {HistType::kTH1F, {{200, 0, 40000}}}},
     {"MC/Control/hTracksWithMcParticle", "Reco tracks with an associated mcParticle (of 4);N;Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}},

     // MC reco<->gen match histograms
     {"MC/Match/hRecoEvents", "Reco events entering the match check;;Counts", {HistType::kTH1F, {{1, 0, 1}}}},
     {"MC/Match/hMatchedGenM", "Matched candidates - generated M;m_{gen} (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
     {"MC/Match/hMatchedGenPt", "Matched candidates - generated p_{T};p_{T,gen} (GeV/c);Counts", {HistType::kTH1F, {{1000, 0.0, 1.1}}}},
     {"MC/Match/hMatchedGenY", "Matched candidates - generated y;y_{gen};Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
     {"MC/Match/hRecoVsGenM", "Reco vs Gen M;m_{gen} (GeV/c^{2});m_{reco} (GeV/c^{2})", {HistType::kTH2F, {{500, 0.0, 5.0}, {500, 0.0, 5.0}}}}}};

  // Helper functions for kinematic calculations
  static float pt(float px, float py) { return std::sqrt(px * px + py * py); }

  static float eta(float px, float py, float pz)
  {
    if (std::abs(pz) > 1e-10) {
      float p = std::sqrt(px * px + py * py + pz * pz);
      return 0.5f * std::log((p + pz) / (p - pz));
    }
    return 0.0f;
  }

  static float phi(float px, float py)
  {
    if (std::abs(px) > 1e-10 || std::abs(py) > 1e-10) {
      return std::atan2(py, px);
    }
    return 0.0f;
  }

  // Generic "common mother"
  struct MotherInfo {
    bool found = false;
    int pdg = 0;
    float px = 0, py = 0;
  };

  // 4 particles share the same direct mother
  template <typename T>
  MotherInfo findCommonMotherGeneric(const std::vector<T>& daughters)
  {
    MotherInfo info;
    if (daughters.empty() || !daughters[0].has_mothers()) {
      return info;
    }
    auto firstMothers = daughters[0].template mothers_as<UDMcParticles>();
    if (firstMothers.begin() == firstMothers.end()) {
      return info;
    }
    auto motherIt = firstMothers.begin();
    int64_t motherGlobalIndex = motherIt->globalIndex();
    for (size_t i = 1; i < daughters.size(); i++) {
      if (!daughters[i].has_mothers()) {
        return info;
      }
      auto iMothers = daughters[i].template mothers_as<UDMcParticles>();
      if (iMothers.begin() == iMothers.end() || iMothers.begin()->globalIndex() != motherGlobalIndex) {
        return info;
      }
    }
    const auto& mother = *motherIt;
    info.found = true;
    info.pdg = mother.pdgCode();
    info.px = mother.px();
    info.py = mother.py();
    return info;
  }

  // Charge sign from PDG code
  static int signFromPdg(int pdgCode)
  {
    if (pdgCode == 211) {
      return 1;
    }
    if (pdgCode == -211) {
      return -1;
    }
    return (pdgCode > 0) ? 1 : (pdgCode < 0) ? -1
                                             : 0;
  }

  void init(InitContext&)
  {
    // Configure event flow histogram labels
    auto hFlow = registry.get<TH1>(HIST("Events/Flow"));
    hFlow->GetXaxis()->SetBinLabel(1, "All events");
    hFlow->GetXaxis()->SetBinLabel(2, "ITS-TPC cut");
    hFlow->GetXaxis()->SetBinLabel(3, "SBP cut");
    hFlow->GetXaxis()->SetBinLabel(4, "ITS ROFb cut");
    hFlow->GetXaxis()->SetBinLabel(5, "TFB cut");
    hFlow->GetXaxis()->SetBinLabel(6, "Gap Side cut");
    hFlow->GetXaxis()->SetBinLabel(7, "ZDC energy cut");
    hFlow->GetXaxis()->SetBinLabel(8, "PV contrib cut");
    hFlow->GetXaxis()->SetBinLabel(9, "Z vtx cut");
    hFlow->GetXaxis()->SetBinLabel(10, "4 tracks cut");
    hFlow->GetXaxis()->SetBinLabel(11, "System cuts (M,Pt,Y)");

    auto hFlowDetailed = registry.get<TH1>(HIST("Events/FlowDetailed"));
    hFlowDetailed->GetXaxis()->SetBinLabel(1, "All events");
    hFlowDetailed->GetXaxis()->SetBinLabel(2, "vtxITSTPC");
    hFlowDetailed->GetXaxis()->SetBinLabel(3, "sbp");
    hFlowDetailed->GetXaxis()->SetBinLabel(4, "itsROFb");
    hFlowDetailed->GetXaxis()->SetBinLabel(5, "tfb");
    hFlowDetailed->GetXaxis()->SetBinLabel(6, "gapSide");
    hFlowDetailed->GetXaxis()->SetBinLabel(7, "FV0 < cut");
    hFlowDetailed->GetXaxis()->SetBinLabel(8, "FT0A < cut");
    hFlowDetailed->GetXaxis()->SetBinLabel(9, "FT0C < cut");
    hFlowDetailed->GetXaxis()->SetBinLabel(10, "ZDC energy < cut");
    hFlowDetailed->GetXaxis()->SetBinLabel(11, "numContrib == 4");
    hFlowDetailed->GetXaxis()->SetBinLabel(12, "posZ < cut");
    hFlowDetailed->GetXaxis()->SetBinLabel(13, "4 tracks");

    // Readable labels for Events/hRecoMode
    auto hRecoMode = registry.get<TH1>(HIST("Events/hRecoMode"));
    hRecoMode->GetXaxis()->SetBinLabel(1, "STD");
    hRecoMode->GetXaxis()->SetBinLabel(2, "UPC");

    auto hChargeBefore = registry.get<TH1>(HIST("System/hTotalChargeBefore"));
    hChargeBefore->GetXaxis()->SetBinLabel(1, "----");
    hChargeBefore->GetXaxis()->SetBinLabel(3, "---+");
    hChargeBefore->GetXaxis()->SetBinLabel(5, "+-+-");
    hChargeBefore->GetXaxis()->SetBinLabel(7, "+++-");
    hChargeBefore->GetXaxis()->SetBinLabel(9, "++++");

    auto hCharge = registry.get<TH1>(HIST("System/hTotalCharge"));
    hCharge->GetXaxis()->SetBinLabel(1, "----");
    hCharge->GetXaxis()->SetBinLabel(3, "---+");
    hCharge->GetXaxis()->SetBinLabel(5, "+-+-");
    hCharge->GetXaxis()->SetBinLabel(7, "+++-");
    hCharge->GetXaxis()->SetBinLabel(9, "++++");

    auto hMChargeBefore = registry.get<TH2>(HIST("System/hMVsTotalChargeBefore"));
    hMChargeBefore->GetYaxis()->SetBinLabel(1, "----");
    hMChargeBefore->GetYaxis()->SetBinLabel(3, "---+");
    hMChargeBefore->GetYaxis()->SetBinLabel(5, "+-+-");
    hMChargeBefore->GetYaxis()->SetBinLabel(7, "+++-");
    hMChargeBefore->GetYaxis()->SetBinLabel(9, "++++");

    auto hMCharge = registry.get<TH2>(HIST("System/hMVsTotalCharge"));
    hMCharge->GetYaxis()->SetBinLabel(1, "----");
    hMCharge->GetYaxis()->SetBinLabel(3, "---+");
    hMCharge->GetYaxis()->SetBinLabel(5, "+-+-");
    hMCharge->GetYaxis()->SetBinLabel(7, "+++-");
    hMCharge->GetYaxis()->SetBinLabel(9, "++++");

    auto hReject = registry.get<TH1>(HIST("Tracks/RejectionReasons"));
    hReject->GetXaxis()->SetBinLabel(1, "All tracks");
    hReject->GetXaxis()->SetBinLabel(2, "isPVContributor");
    hReject->GetXaxis()->SetBinLabel(3, "hasITS && hasTPC");
    hReject->GetXaxis()->SetBinLabel(4, "pT > 0.1");
    hReject->GetXaxis()->SetBinLabel(5, "tpcChi2NCl < cut");
    hReject->GetXaxis()->SetBinLabel(6, "itsChi2NCl < cut");
    hReject->GetXaxis()->SetBinLabel(7, "tpcNClsFindable > cut");
    hReject->GetXaxis()->SetBinLabel(8, "|tpcNSigmaPi| < cut");
    hReject->GetXaxis()->SetBinLabel(9, "|eta| < cut");
    hReject->GetXaxis()->SetBinLabel(10, "|dcaZ| < cut");
    hReject->GetXaxis()->SetBinLabel(11, "|dcaXY| < cut");
    hReject->GetXaxis()->SetBinLabel(12, "Accepted tracks");

    auto hPdg = mcRegistry.get<TH1>(HIST("MC/Tracks/hPdgCode"));
    if (hPdg) {
      int bin211 = hPdg->GetXaxis()->FindBin(211);
      int binM211 = hPdg->GetXaxis()->FindBin(-211);
      int bin30113 = hPdg->GetXaxis()->FindBin(30113);
      if (bin211 > 0 && bin211 <= hPdg->GetNbinsX())
        hPdg->GetXaxis()->SetBinLabel(bin211, "#pi^{+}");
      if (binM211 > 0 && binM211 <= hPdg->GetNbinsX())
        hPdg->GetXaxis()->SetBinLabel(binM211, "#pi^{-}");
      if (bin30113 > 0 && bin30113 <= hPdg->GetNbinsX())
        hPdg->GetXaxis()->SetBinLabel(bin30113, "rho'");
    }

    auto hRhoFound = mcRegistry.get<TH1>(HIST("MC/RhoPrime/hFound"));
    hRhoFound->GetXaxis()->SetBinLabel(1, "Not Found");
    hRhoFound->GetXaxis()->SetBinLabel(2, "Found");

    auto hMcCharge = mcRegistry.get<TH1>(HIST("MC/System/hTotalCharge"));
    hMcCharge->GetXaxis()->SetBinLabel(1, "----");
    hMcCharge->GetXaxis()->SetBinLabel(3, "---+");
    hMcCharge->GetXaxis()->SetBinLabel(5, "+-+-");
    hMcCharge->GetXaxis()->SetBinLabel(7, "+++-");
    hMcCharge->GetXaxis()->SetBinLabel(9, "++++");

    auto hMcSummary = mcRegistry.get<TH1>(HIST("MC/Summary/hEventCounter"));
    hMcSummary->GetXaxis()->SetBinLabel(1, "Generated (denominator, gen all)");
    hMcSummary->GetXaxis()->SetBinLabel(2, "Reconstructed + matched (numerator)");

    auto hMatchStatus = mcRegistry.get<TH1>(HIST("MC/Summary/hMatchStatus"));
    hMatchStatus->GetXaxis()->SetBinLabel(1, "No mcParticle on some track");
    hMatchStatus->GetXaxis()->SetBinLabel(2, "Has mcParticle but no common mother");
    hMatchStatus->GetXaxis()->SetBinLabel(3, "Match found");

    auto hMcRecoMode = mcRegistry.get<TH1>(HIST("MC/Summary/hRecoMode"));
    hMcRecoMode->GetXaxis()->SetBinLabel(1, "STD");
    hMcRecoMode->GetXaxis()->SetBinLabel(2, "UPC");

    auto hMcRecoCharge = mcRegistry.get<TH1>(HIST("MC/Summary/hRecoTotalCharge"));
    hMcRecoCharge->GetXaxis()->SetBinLabel(1, "----");
    hMcRecoCharge->GetXaxis()->SetBinLabel(3, "---+");
    hMcRecoCharge->GetXaxis()->SetBinLabel(5, "+-+-");
    hMcRecoCharge->GetXaxis()->SetBinLabel(7, "+++-");
    hMcRecoCharge->GetXaxis()->SetBinLabel(9, "++++");

    if (doprocessDataCols && (doprocessMcCols || doprocessMcGenAll)) {
      LOGP(fatal,
           "Invalid config: processDataCols must not run together with processMcCols/processMcGenAll "
           "(it would duplicate SystemTree rows). Disable one of the two groups in your configuration.json.");
    }
  }

  // Collision + track selection
  template <typename ColType, typename TracksType>
  std::vector<typename TracksType::iterator> selectAndFillSystemTree(ColType const& collision, TracksType const& tracks, bool& outFilled)
  {
    outFilled = false;
    std::vector<typename TracksType::iterator> goodTracks;

    registry.fill(HIST("Events/Flow"), 0);
    registry.fill(HIST("Events/FlowDetailed"), 0);

    registry.fill(HIST("Events/VertexZ"), collision.posZ());
    registry.fill(HIST("Events/NumContrib"), collision.numContrib());
    registry.fill(HIST("Events/FV0Amplitude"), collision.totalFV0AmplitudeA());
    registry.fill(HIST("Events/FT0AmplitudeA"), collision.totalFT0AmplitudeA());
    registry.fill(HIST("Events/FT0AmplitudeC"), collision.totalFT0AmplitudeC());
    registry.fill(HIST("Events/ZDCEnergy"), collision.energyCommonZNA());
    registry.fill(HIST("Events/ZDCEnergy"), collision.energyCommonZNC());

    if (collision.vtxITSTPC() != vtxITSTPCcut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/Flow"), 1);
    registry.fill(HIST("Events/FlowDetailed"), 1);

    if (collision.sbp() != sbpCut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/Flow"), 2);
    registry.fill(HIST("Events/FlowDetailed"), 2);

    if (collision.itsROFb() != itsROFbCut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/Flow"), 3);
    registry.fill(HIST("Events/FlowDetailed"), 3);

    if (collision.tfb() != tfbCut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/Flow"), 4);
    registry.fill(HIST("Events/FlowDetailed"), 4);

    if (specifyGapSide && collision.gapSide() != gapSide) {
      return goodTracks;
    }

    registry.fill(HIST("Events/Flow"), 5);
    registry.fill(HIST("Events/FlowDetailed"), 5);

    if (collision.totalFV0AmplitudeA() > fv0Cut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/FlowDetailed"), 6);

    if (collision.totalFT0AmplitudeA() > ft0aCut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/FlowDetailed"), 7);

    if (collision.totalFT0AmplitudeC() > ft0cCut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/FlowDetailed"), 8);

    if (collision.energyCommonZNA() > zdcCut || collision.energyCommonZNC() > zdcCut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/Flow"), 6);
    registry.fill(HIST("Events/FlowDetailed"), 9);

    if (collision.numContrib() != numPVContrib) {
      return goodTracks;
    }
    registry.fill(HIST("Events/Flow"), 7);
    registry.fill(HIST("Events/FlowDetailed"), 10);

    if (std::abs(collision.posZ()) > vZCut) {
      return goodTracks;
    }
    registry.fill(HIST("Events/Flow"), 8);
    registry.fill(HIST("Events/FlowDetailed"), 11);

    // --- Track selection: up to 4 good tracks, charge combination ---
    goodTracks.reserve(4);
    for (const auto& track : tracks) {
      registry.fill(HIST("Tracks/RejectionReasons"), 0);

      if (useOnlyPVtracks && !track.isPVContributor()) {
        registry.fill(HIST("Tracks/RejectionReasons"), 1);
        continue;
      }
      if (!track.hasITS() || !track.hasTPC()) {
        registry.fill(HIST("Tracks/RejectionReasons"), 2);
        continue;
      }

      registry.fill(HIST("Tracks/Pt"), track.pt());
      registry.fill(HIST("Tracks/Eta"), eta(track.px(), track.py(), track.pz()));
      registry.fill(HIST("Tracks/TPCNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("Tracks/TPCChi2NCl"), track.tpcChi2NCl());
      registry.fill(HIST("Tracks/ITSChi2NCl"), track.itsChi2NCl());
      registry.fill(HIST("Tracks/DCASpectrum"), std::hypot(track.dcaXY(), track.dcaZ()));
      registry.fill(HIST("Tracks/ChargeDistribution"), track.sign());
      registry.fill(HIST("Tracks/TPCClusters"), track.tpcNClsFindable());

      if (track.pt() <= 0.1f) {
        registry.fill(HIST("Tracks/RejectionReasons"), 3);
        continue;
      }
      if (track.tpcChi2NCl() > tpcChi2NClsCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 4);
        continue;
      }
      if (track.itsChi2NCl() > itsChi2NClsCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 5);
        continue;
      }
      if (track.tpcNClsFindable() < minTPCFindableClusters) {
        registry.fill(HIST("Tracks/RejectionReasons"), 6);
        continue;
      }
      if (std::abs(track.tpcNSigmaPi()) > nSigmaTPCcut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 7);
        continue;
      }
      if (std::abs(eta(track.px(), track.py(), track.pz())) > etaCut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 8);
        continue;
      }
      if (std::abs(track.dcaZ()) > dcaZcut) {
        registry.fill(HIST("Tracks/RejectionReasons"), 9);
        continue;
      }
      float maxDCAxy = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
      if (dcaXYcut == 0 && std::fabs(track.dcaXY()) > maxDCAxy) {
        registry.fill(HIST("Tracks/RejectionReasons"), 10);
        continue;
      }

      registry.fill(HIST("Tracks/RejectionReasons"), 11);
      goodTracks.push_back(track);
      if (goodTracks.size() == 4) {
        break;
      }
    }

    // Basic control
    registry.fill(HIST("Tracks/NGoodTracksPerEvent"), goodTracks.size());

    if (goodTracks.size() != 4) {
      return goodTracks;
    }
    registry.fill(HIST("Events/Flow"), 9);
    registry.fill(HIST("Events/FlowDetailed"), 12);

    // Real total charge
    int totalCharge = 0;
    for (const auto& track : goodTracks) {
      totalCharge += track.sign();
    }
    bool isChargeZero = (totalCharge == 0);
    registry.fill(HIST("System/hTotalChargeBefore"), totalCharge);

    PxPyPzMVector fourPionSystem;
    for (const auto& track : goodTracks) {
      fourPionSystem += PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged);
    }

    registry.fill(HIST("Cuts/MBefore"), fourPionSystem.M());
    registry.fill(HIST("Cuts/PtBefore"), fourPionSystem.Pt());
    registry.fill(HIST("System/hMVsTotalChargeBefore"), fourPionSystem.M(), totalCharge);

    if (fourPionSystem.M() < systemMassMinCut || fourPionSystem.M() > systemMassMaxCut) {
      return goodTracks;
    }
    if (fourPionSystem.Pt() > systemPtCut) {
      return goodTracks;
    }
    if (std::abs(fourPionSystem.Rapidity()) > systemYCut) {
      return goodTracks;
    }

    registry.fill(HIST("Cuts/MAfter"), fourPionSystem.M());
    registry.fill(HIST("Cuts/PtAfter"), fourPionSystem.Pt());
    registry.fill(HIST("System/hM"), fourPionSystem.M());
    registry.fill(HIST("System/hPt"), fourPionSystem.Pt());
    registry.fill(HIST("System/hEta"), fourPionSystem.Eta());
    registry.fill(HIST("System/hPhi"), fourPionSystem.Phi() + o2::constants::math::PI);
    registry.fill(HIST("System/hY"), fourPionSystem.Rapidity());
    registry.fill(HIST("System/hTotalCharge"), totalCharge);
    registry.fill(HIST("System/hMVsTotalCharge"), fourPionSystem.M(), totalCharge);

    std::vector<float> trackPts, trackEtas, trackPhis;
    std::vector<int> trackSigns, trackIDs;
    std::vector<float> tpcNSigmasEl, tpcNSigmasPi, tpcNSigmasKa, tpcNSigmasPr;
    for (size_t i = 0; i < goodTracks.size(); i++) {
      const auto& track = goodTracks[i];
      trackPts.push_back(track.pt());
      trackEtas.push_back(eta(track.px(), track.py(), track.pz()));
      trackPhis.push_back(phi(track.px(), track.py()));
      trackSigns.push_back(track.sign());
      tpcNSigmasEl.push_back(track.tpcNSigmaEl());
      tpcNSigmasPi.push_back(track.tpcNSigmaPi());
      tpcNSigmasKa.push_back(track.tpcNSigmaKa());
      tpcNSigmasPr.push_back(track.tpcNSigmaPr());
      trackIDs.push_back(static_cast<int>(i));
    }

    bool isReconstructedWithUPC = (collision.flags() == 1);

    registry.fill(HIST("Events/Flow"), 10);
    registry.fill(HIST("Events/hRecoMode"), isReconstructedWithUPC ? 1 : 0);

    systemTree(
      collision.runNumber(),
      fourPionSystem.M(), fourPionSystem.Pt(), fourPionSystem.Rapidity(), fourPionSystem.Phi(),
      collision.posX(), collision.posY(), collision.posZ(),
      totalCharge,
      collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(),
      collision.numContrib(),
      trackSigns, trackPts, trackEtas, trackPhis,
      tpcNSigmasEl, tpcNSigmasPi, tpcNSigmasKa, tpcNSigmasPr,
      trackIDs, isReconstructedWithUPC,
      collision.timeZNA(), collision.timeZNC(), collision.energyCommonZNA(), collision.energyCommonZNC(),
      isChargeZero, collision.occupancyInTime(), collision.hadronicRate());

    outFilled = true;
    return goodTracks;
  } // end selectAndFillSystemTree

  int getMcRunNumber(aod::BCs const& bcs)
  {
    if (bcs.size() == 0) {
      return -1;
    }
    auto bc = bcs.begin();
    return bc.runNumber();
  }

  // processDataCols: RD or MC reco
  void processDataCols(UDCollisions::iterator const& collision, UDtracks const& tracks)
  {
    bool filled = false;
    selectAndFillSystemTree(collision, tracks, filled);
  }
  PROCESS_SWITCH(upcRhoPrimeAnalysis, processDataCols, "Process real data or MC reco", true);

  // processMcCols: MC only
  void processMcCols(UDCollisionsMC::iterator const& collision, UDtracksMC const& tracks, UDMcParticles const&, UDMcCollisions const&, aod::BCs const& bcs)
  {
    if (genId != -1) {
      if (!collision.has_udMcCollision() || collision.template udMcCollision_as<UDMcCollisions>().generatorsID() != genId) {
        return;
      }
    }

    bool filled = false;
    auto selTrks = selectAndFillSystemTree(collision, tracks, filled);
    if (!filled) {
      return;
    }

    bool isReconstructedWithUPC = (collision.flags() == 1);

    mcRegistry.fill(HIST("MC/Match/hRecoEvents"), 0);
    mcRegistry.fill(HIST("MC/Summary/hRecoMode"), isReconstructedWithUPC ? 1 : 0);

    int recoTotalCharge = 0;
    for (const auto& track : selTrks) {
      recoTotalCharge += track.sign();
    }
    mcRegistry.fill(HIST("MC/Summary/hRecoTotalCharge"), recoTotalCharge);

    // Defaults (pdg=0 / kinematics=-999 => "no match")
    int motherPdg = 0;
    float motherPt = -999, motherPhi = -999, motherMass = -999, motherRap = -999;
    int mcTotalCharge = 0;
    int trackPdgs[4] = {0, 0, 0, 0};
    float trackPts[4] = {-999, -999, -999, -999};
    float trackEtas[4] = {-999, -999, -999, -999};
    float trackPhis[4] = {-999, -999, -999, -999};
    int trackSigns[4] = {0, 0, 0, 0};
    int isPrimary[4] = {0, 0, 0, 0};
    float mcPosX = -999, mcPosY = -999, mcPosZ = -999;
    int mcCollisionIndex = -1;
    int mcRunNumber = getMcRunNumber(bcs);

    // Basic control
    int nTracksWithMcParticle = 0;
    for (const auto& track : selTrks) {
      if (track.has_udMcParticle()) {
        nTracksWithMcParticle++;
      }
    }
    mcRegistry.fill(HIST("MC/Control/hTracksWithMcParticle"), nTracksWithMcParticle);

    // All 4 reco tracks need an associated mcParticle to look for the mother
    std::vector<typename UDMcParticles::iterator> mcParts;
    bool allHaveMcParticle = true;
    for (const auto& track : selTrks) {
      if (!track.has_udMcParticle()) {
        allHaveMcParticle = false;
        break;
      }
      mcParts.push_back(track.template udMcParticle_as<UDMcParticles>());
    }

    MotherInfo mi;
    if (allHaveMcParticle) {
      mi = findCommonMotherGeneric(mcParts);
    }

    if (!allHaveMcParticle) {
      mcRegistry.fill(HIST("MC/Summary/hMatchStatus"), 0);
    } else if (!mi.found) {
      mcRegistry.fill(HIST("MC/Summary/hMatchStatus"), 1);
    } else {
      mcRegistry.fill(HIST("MC/Summary/hMatchStatus"), 2);
      mcRegistry.fill(HIST("MC/Summary/hEventCounter"), 1);
    }

    if (allHaveMcParticle) {
      for (size_t i = 0; i < mcParts.size() && i < 4; i++) {
        trackPdgs[i] = mcParts[i].pdgCode();
        trackPts[i] = pt(mcParts[i].px(), mcParts[i].py());
        trackEtas[i] = eta(mcParts[i].px(), mcParts[i].py(), mcParts[i].pz());
        trackPhis[i] = phi(mcParts[i].px(), mcParts[i].py());
        trackSigns[i] = signFromPdg(trackPdgs[i]);
        isPrimary[i] = mcParts[i].isPhysicalPrimary() ? 1 : 0;
        mcTotalCharge += trackSigns[i];
      }

      auto mcCollision = mcParts[0].template udMcCollision_as<UDMcCollisions>();
      mcPosX = mcCollision.posX();
      mcPosY = mcCollision.posY();
      mcPosZ = mcCollision.posZ();
      mcCollisionIndex = static_cast<int>(mcCollision.globalIndex());
    }

    if (mi.found) {
      motherPdg = mi.pdg;
      motherPt = pt(mi.px, mi.py);
      motherPhi = phi(mi.px, mi.py);
      PxPyPzMVector genSystem;
      for (const auto& p : mcParts) {
        genSystem += PxPyPzMVector(p.px(), p.py(), p.pz(), o2::constants::physics::MassPionCharged);
      }
      motherMass = genSystem.M();
      motherRap = genSystem.Rapidity();

      mcRegistry.fill(HIST("MC/Control/hMotherPdg"), motherPdg);
      mcRegistry.fill(HIST("MC/Match/hMatchedGenM"), motherMass);
      mcRegistry.fill(HIST("MC/Match/hMatchedGenPt"), motherPt);
      mcRegistry.fill(HIST("MC/Match/hMatchedGenY"), motherRap);

      PxPyPzMVector recoSystem;
      for (const auto& track : selTrks) {
        recoSystem += PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged);
      }
      mcRegistry.fill(HIST("MC/Match/hRecoVsGenM"), motherMass, recoSystem.M());
      if (motherPdg == 30113) { // only when the mother found is actually the rho prime
        if (isReconstructedWithUPC) {
          mcRegistry.fill(HIST("MC/RhoPrime/hMassUPC"), recoSystem.M());
        } else {
          mcRegistry.fill(HIST("MC/RhoPrime/hMassSTD"), recoSystem.M());
        }
      }
    }

    fourPiMcMatchTree(
      isReconstructedWithUPC,
      motherPdg, motherPt, motherPhi, motherMass, motherRap, mcTotalCharge,
      trackPdgs, trackPts, trackEtas, trackPhis, trackSigns, isPrimary,
      mcPosX, mcPosY, mcPosZ, mcCollisionIndex, mcRunNumber,
      systemTree.lastIndex());
  }
  PROCESS_SWITCH(upcRhoPrimeAnalysis, processMcCols, "Match MC reco<->gen (MC only)", false);

  // processMcGenAll: MC All generated collisions
  void processMcGenAll(UDMcCollisions::iterator const& mcCollision, UDMcParticles const& mcParticles, aod::BCs const& bcs)
  {
    if (genId != -1 && mcCollision.generatorsID() != genId) {
      return;
    }

    int mcRunNumber = getMcRunNumber(bcs);

    mcRegistry.fill(HIST("MC/Events/hAllEvents"), 0);
    mcRegistry.fill(HIST("MC/Events/hVertexZ"), mcCollision.posZ());

    std::vector<decltype(mcParticles.begin())> primaries;
    for (const auto& part : mcParticles) {
      if (part.isPhysicalPrimary()) {
        primaries.push_back(part);
      }
    }
    mcRegistry.fill(HIST("MC/Events/hNPrimaries"), primaries.size());

    if (primaries.size() != 4) {
      return;
    }
    mcRegistry.fill(HIST("MC/Summary/hEventCounter"), 0);

    int trackPdgs[4] = {0, 0, 0, 0};
    float trackPts[4] = {0, 0, 0, 0};
    float trackEtas[4] = {0, 0, 0, 0};
    float trackPhis[4] = {0, 0, 0, 0};
    int trackSigns[4] = {0, 0, 0, 0};
    int isPrimary[4] = {1, 1, 1, 1};
    int mcTotalCharge = 0;

    PxPyPzMVector genSystem;
    for (size_t i = 0; i < 4; i++) {
      const auto& p = primaries[i];
      trackPdgs[i] = p.pdgCode();
      trackPts[i] = pt(p.px(), p.py());
      trackEtas[i] = eta(p.px(), p.py(), p.pz());
      trackPhis[i] = phi(p.px(), p.py());
      trackSigns[i] = signFromPdg(trackPdgs[i]);
      mcTotalCharge += trackSigns[i];

      mcRegistry.fill(HIST("MC/Tracks/hPt"), trackPts[i]);
      mcRegistry.fill(HIST("MC/Tracks/hEta"), trackEtas[i]);
      mcRegistry.fill(HIST("MC/Tracks/hPhi"), trackPhis[i] + o2::constants::math::PI); // same 0-2pi convention as System/hPhi
      mcRegistry.fill(HIST("MC/Tracks/hPdgCode"), trackPdgs[i]);

      genSystem += PxPyPzMVector(p.px(), p.py(), p.pz(), o2::constants::physics::MassPionCharged);
    }

    // Generic common-mother search
    int motherPdg = 0;
    float motherPt = -999, motherPhi = -999, motherMass = -999, motherRap = -999;
    MotherInfo mi = findCommonMotherGeneric(primaries);
    mcRegistry.fill(HIST("MC/Control/hNDaughtersOfMother"), mi.found ? 4 : 0);
    if (mi.found) {
      motherPdg = mi.pdg;
      motherPt = pt(mi.px, mi.py);
      motherPhi = phi(mi.px, mi.py);
      motherMass = genSystem.M();
      motherRap = genSystem.Rapidity();
      mcRegistry.fill(HIST("MC/Control/hMotherPdg"), motherPdg);
    }

    mcRegistry.fill(HIST("MC/System/hM"), genSystem.M());
    mcRegistry.fill(HIST("MC/System/hPt"), genSystem.Pt());
    mcRegistry.fill(HIST("MC/System/hY"), genSystem.Rapidity());
    mcRegistry.fill(HIST("MC/System/hMvsPt"), genSystem.M(), genSystem.Pt());
    mcRegistry.fill(HIST("MC/System/hMvsY"), genSystem.M(), genSystem.Rapidity());
    mcRegistry.fill(HIST("MC/System/hTotalCharge"), mcTotalCharge);

    if (motherPdg == 30113) {
      mcRegistry.fill(HIST("MC/RhoPrime/hFound"), 1);
      mcRegistry.fill(HIST("MC/RhoPrime/hMass"), motherMass);
      mcRegistry.fill(HIST("MC/RhoPrime/hPt"), motherPt);
    } else {
      mcRegistry.fill(HIST("MC/RhoPrime/hFound"), 0);
    }

    fourPiMcGenAllTree(
      motherPdg, motherPt, motherPhi, motherMass, motherRap, mcTotalCharge,
      trackPdgs, trackPts, trackEtas, trackPhis, trackSigns, isPrimary,
      mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
      static_cast<int>(mcCollision.globalIndex()), mcRunNumber);

    mcRegistry.fill(HIST("MC/Events/hAccepted"), 0);
  }
  PROCESS_SWITCH(upcRhoPrimeAnalysis, processMcGenAll, "All generated collisions (MC only)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<upcRhoPrimeAnalysis>(cfgc, TaskName{"upc-rho-prime-analysis"})};
}
