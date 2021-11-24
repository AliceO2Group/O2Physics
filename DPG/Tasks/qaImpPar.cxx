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
/// \author Mattia Faggin <mattia.faggin@cern.ch>, Padova University and INFN
///
/// Event selection: o2-analysis-timestamp --aod-file AO2D.root -b | o2-analysis-event-selection -b | ---> not working with Run2 converted data/MC
/// Track selection: o2-analysis-trackextension | o2-analysis-trackselection | ---> add --isRun3 1 with Run 3 data/MC (then global track selection works)
/// PID: o2-analysis-pid-tpc-full | o2-analysis-pid-tof-full
/// Working configuration (2021 Oct 20th): o2-analysis-trackextension -b --aod-file ./AO2D.root | o2-analysis-trackselection -b --isRun3 1 | o2-analysis-pid-tpc-full -b | o2-analysis-pid-tof-full -b | o2-analysis-pp-qa-impact-parameter -b

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h" // for propagation to primary vertex
#include "Common/Core/MC.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDResponse.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/TrackSelection.h"

#include "iostream"

using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// QA task for impact parameter distribution monitoring
struct QaImpactPar {

  /// Input parameters
  ConfigurableAxis binningImpPar{"binningImpPar", {200, -500.f, 500.f}, "Impact parameter binning"};
  //Configurable<int> numberContributorsMin{"numberContributorsMin", 0, "Minimum number of contributors for the primary vertex"};
  Configurable<float> zVtxMax{"zVtxMax", 10.f, "Maximum value for |z_vtx|"};
  //Configurable<int> keepOnlyGlobalTracks{"keepOnlyGlobalTracks", 1, "Keep only global tracks or not"};
  Configurable<float> ptMin{"ptMin", 0.1f, "Minimum track pt [GeV/c]"};
  Configurable<float> nSigmaTPCPionMin{"nSigmaTPCPionMin", -99999.f, "Minimum nSigma value in TPC, pion hypothesis"};
  Configurable<float> nSigmaTPCPionMax{"nSigmaTPCPionMax", 99999.f, "Maximum nSigma value in TPC, pion hypothesis"};
  Configurable<float> nSigmaTPCKaonMin{"nSigmaTPCKaonMin", -99999.f, "Minimum nSigma value in TPC, kaon hypothesis"};
  Configurable<float> nSigmaTPCKaonMax{"nSigmaTPCKaonMax", 99999.f, "Maximum nSigma value in TPC, kaon hypothesis"};
  Configurable<float> nSigmaTPCProtonMin{"nSigmaTPCProtonMin", -99999.f, "Minimum nSigma value in TPC, proton hypothesis"};
  Configurable<float> nSigmaTPCProtonMax{"nSigmaTPCProtonMax", 99999.f, "Maximum nSigma value in TPC, proton hypothesis"};
  Configurable<float> nSigmaTOFPionMin{"nSigmaTOFPionMin", -99999.f, "Minimum nSigma value in TOF, pion hypothesis"};
  Configurable<float> nSigmaTOFPionMax{"nSigmaTOFPionMax", 99999.f, "Maximum nSigma value in TOF, pion hypothesis"};
  Configurable<float> nSigmaTOFKaonMin{"nSigmaTOFKaonMin", -99999.f, "Minimum nSigma value in TOF, kaon hypothesis"};
  Configurable<float> nSigmaTOFKaonMax{"nSigmaTOFKaonMax", 99999.f, "Maximum nSigma value in TOF, kaon hypothesis"};
  Configurable<float> nSigmaTOFProtonMin{"nSigmaTOFProtonMin", -99999.f, "Minimum nSigma value in TOF, proton hypothesis"};
  Configurable<float> nSigmaTOFProtonMax{"nSigmaTOFProtonMax", 99999.f, "Maximum nSigma value in TOF, proton hypothesis"};

  /// Selections with Filter (from o2::framework::expressions)
  // Primary vertex |z_vtx|<XXX cm
  Filter collisionZVtxFilter = nabs(o2::aod::collision::posZ) < zVtxMax;
  // Global tracks
  // with Run 3 data/MC enable '--isRun3 1' option
  Filter globalTrackFilter = (o2::aod::track::isGlobalTrack == (uint8_t) true); /// filterbit 4 track selections + tight DCA cuts
  // Pt selection
  Filter ptMinFilter = o2::aod::track::pt > ptMin;

  /// Histogram registry (from o2::framework)
  HistogramRegistry histograms{"HistogramsImpParQA"};

  /// init function - declare and define histograms
  void init(InitContext&)
  {
    // Primary vertex
    const AxisSpec collisionZAxis{100, -20.f, 20.f, "Z (cm)"};
    const AxisSpec collisionNumberContributorAxis{1000, 0, 1000, "Number of contributors"};

    histograms.add("vertexZ", "", kTH1D, {collisionZAxis});
    histograms.add("numberContributors", "", kTH1D, {collisionNumberContributorAxis});

    // tracks
    const AxisSpec trackPtAxis{100, 0.f, 10.f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec trackEtaAxis{40, -2.f, 2.f, "#it{#eta}"};
    const AxisSpec trackPhiAxis{24, 0.f, TMath::TwoPi(), "#varphi"};
    const AxisSpec trackImpParRPhiAxis{binningImpPar, "#it{d}_{r#it{#varphi}} (#mum)"};
    const AxisSpec trackImpParZAxis{binningImpPar, "#it{d}_{z} (#mum)"};
    const AxisSpec trackNSigmaTPCPionAxis{20, -10.f, 10.f, "Number of #sigma TPC #pi^{#pm}"};
    const AxisSpec trackNSigmaTPCKaonAxis{20, -10.f, 10.f, "Number of #sigma TPC K^{#pm}"};
    const AxisSpec trackNSigmaTPCProtonAxis{20, -10.f, 10.f, "Number of #sigma TPC proton"};
    const AxisSpec trackNSigmaTOFPionAxis{20, -10.f, 10.f, "Number of #sigma TOF #pi^{#pm}"};
    const AxisSpec trackNSigmaTOFKaonAxis{20, -10.f, 10.f, "Number of #sigma TOF K^{#pm}"};
    const AxisSpec trackNSigmaTOFProtonAxis{20, -10.f, 10.f, "Number of #sigma TOF proton"};

    histograms.add("pt", "", kTH1D, {trackPtAxis});
    histograms.add("h4ImpPar", "", kTHnD, {trackPtAxis, trackImpParRPhiAxis, trackEtaAxis, trackPhiAxis});
    histograms.add("h4ImpParZ", "", kTHnD, {trackPtAxis, trackImpParZAxis, trackEtaAxis, trackPhiAxis});
    histograms.add("h4ImpPar_Pion", "", kTHnD, {trackPtAxis, trackImpParRPhiAxis, trackEtaAxis, trackPhiAxis});
    histograms.add("h4ImpParZ_Pion", "", kTHnD, {trackPtAxis, trackImpParZAxis, trackEtaAxis, trackPhiAxis});
    histograms.add("h4ImpPar_Kaon", "", kTHnD, {trackPtAxis, trackImpParRPhiAxis, trackEtaAxis, trackPhiAxis});
    histograms.add("h4ImpParZ_Kaon", "", kTHnD, {trackPtAxis, trackImpParZAxis, trackEtaAxis, trackPhiAxis});
    histograms.add("h4ImpPar_Proton", "", kTHnD, {trackPtAxis, trackImpParRPhiAxis, trackEtaAxis, trackPhiAxis});
    histograms.add("h4ImpParZ_Proton", "", kTHnD, {trackPtAxis, trackImpParZAxis, trackEtaAxis, trackPhiAxis});
    histograms.add("hNSigmaTPCPion", "", kTH2D, {trackPtAxis, trackNSigmaTPCPionAxis});
    histograms.add("hNSigmaTPCKaon", "", kTH2D, {trackPtAxis, trackNSigmaTPCKaonAxis});
    histograms.add("hNSigmaTPCProton", "", kTH2D, {trackPtAxis, trackNSigmaTPCProtonAxis});
    histograms.add("hNSigmaTOFPion", "", kTH2D, {trackPtAxis, trackNSigmaTOFPionAxis});
    histograms.add("hNSigmaTOFKaon", "", kTH2D, {trackPtAxis, trackNSigmaTOFKaonAxis});
    histograms.add("hNSigmaTOFProton", "", kTH2D, {trackPtAxis, trackNSigmaTOFProtonAxis});
    histograms.add("hNSigmaTPCPion_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTPCPionAxis});
    histograms.add("hNSigmaTPCKaon_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTPCKaonAxis});
    histograms.add("hNSigmaTPCProton_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTPCProtonAxis});
    histograms.add("hNSigmaTOFPion_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTOFPionAxis});
    histograms.add("hNSigmaTOFKaon_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTOFKaonAxis});
    histograms.add("hNSigmaTOFProton_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTOFProtonAxis});
  }

  /// o2::aod::EvSels makes the execution crash, with the following error message:
  /// [240108:bc-selection-task]: [17:22:28][WARN] CCDB: Did not find an alien token; Cannot serve objects located on alien://
  /// [240108:bc-selection-task]: [17:22:28][ERROR] Requested resource does not exist: http://alice-ccdb.cern.ch/EventSelection/TriggerAliases/1511123421601/
  /// [240108:bc-selection-task]: [17:22:28][FATAL] Trigger aliases are not available in CCDB for run=282341 at timestamp=1511123421601
  //void process(o2::soa::Filtered<o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>>::iterator& collision,
  void process(o2::soa::Filtered<o2::aod::Collisions>::iterator& collision,
               o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksExtended, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>> const& tracks)
  {
    //o2::dataformats::DCA dca;
    // FIXME: get this from CCDB
    //constexpr float magneticField{5.0};      // in kG
    constexpr float toMicrometers = 10000.f; // Conversion from [cm] to [mum]
    const bool isPIDPionApplied = ((nSigmaTPCPionMin > -10.001 && nSigmaTPCPionMax < 10.001) || (nSigmaTOFPionMin > -10.001 && nSigmaTOFPionMax < 10.001));
    const bool isPIDKaonApplied = ((nSigmaTPCKaonMin > -10.001 && nSigmaTPCKaonMax < 10.001) || (nSigmaTOFKaonMin > -10.001 && nSigmaTOFKaonMax < 10.001));
    const bool isPIDProtonApplied = ((nSigmaTPCProtonMin > -10.001 && nSigmaTPCProtonMax < 10.001) || (nSigmaTOFProtonMin > -10.001 && nSigmaTOFProtonMax < 10.001));

    /// trigger selection (remove for the moment, need to join with o2::aod::EvSels)
    //if (useINT7Trigger) {
    //    // from Tutorial/src/multiplicityEventTrackSelection.cxx
    //    if (!collision.alias()[kINT7]) {
    //        return;
    //    }
    //    if (!collision.sel7()) {
    //        return;
    //    }
    //}

    histograms.fill(HIST("vertexZ"), collision.posZ());
    histograms.fill(HIST("numberContributors"), collision.numContrib());

    /// loop over tracks
    float pt = -999.f;
    float impParRPhi = -999.f;
    float impParZ = -999.f;
    float tpcNSigmaPion = -999.f;
    float tpcNSigmaKaon = -999.f;
    float tpcNSigmaProton = -999.f;
    float tofNSigmaPion = -999.f;
    float tofNSigmaKaon = -999.f;
    float tofNSigmaProton = -999.f;
    for (const auto& track : tracks) {

      /// Using the Filter instead
      ///if ((keepOnlyGlobalTracks) && (!track.isGlobalTrack())) {
      ///  /// not a global track (FB 4 with tight DCA cuts)
      ///  continue;
      ///}

      pt = track.pt();
      tpcNSigmaPion = track.tpcNSigmaPi();
      tpcNSigmaKaon = track.tpcNSigmaKa();
      tpcNSigmaProton = track.tpcNSigmaPr();
      tofNSigmaPion = track.tofNSigmaPi();
      tofNSigmaKaon = track.tofNSigmaKa();
      tofNSigmaProton = track.tofNSigmaPr();

      histograms.fill(HIST("pt"), pt);
      histograms.fill(HIST("hNSigmaTPCPion"), pt, tpcNSigmaPion);
      histograms.fill(HIST("hNSigmaTPCKaon"), pt, tpcNSigmaKaon);
      histograms.fill(HIST("hNSigmaTPCProton"), pt, tpcNSigmaProton);
      histograms.fill(HIST("hNSigmaTOFPion"), pt, tofNSigmaPion);
      histograms.fill(HIST("hNSigmaTOFKaon"), pt, tofNSigmaKaon);
      histograms.fill(HIST("hNSigmaTOFProton"), pt, tofNSigmaProton);

      // propagation to primary vertex for DCA
      // "crude" method not working with Run 3 MC productions
      //if (getTrackParCov(track).propagateToDCA(getPrimaryVertex(track.collision()), magneticField, &dca, 100.)) {

      /// propagation ok! Retrieve impact parameter
      // PR "Run 3 DCA extraction #187" required - correct calculation of DCAxy of tracks propagated to the PV
      impParRPhi = toMicrometers * track.dcaXY(); //dca.getY();
      impParZ = toMicrometers * track.dcaZ();     //dca.getY();

      /// all tracks
      histograms.fill(HIST("h4ImpPar"), pt, impParRPhi, track.eta(), track.phi());
      histograms.fill(HIST("h4ImpParZ"), pt, impParZ, track.eta(), track.phi());

      if (isPIDPionApplied && nSigmaTPCPionMin < tpcNSigmaPion && tpcNSigmaPion < nSigmaTPCPionMax && nSigmaTOFPionMin < tofNSigmaPion && tofNSigmaPion < nSigmaTOFPionMax) {
        /// PID selected pions
        histograms.fill(HIST("h4ImpPar_Pion"), pt, impParRPhi, track.eta(), track.phi());
        histograms.fill(HIST("h4ImpParZ_Pion"), pt, impParZ, track.eta(), track.phi());
        histograms.fill(HIST("hNSigmaTPCPion_afterPID"), pt, tpcNSigmaPion);
        histograms.fill(HIST("hNSigmaTOFPion_afterPID"), pt, tofNSigmaPion);
      }
      if (isPIDKaonApplied && nSigmaTPCKaonMin < tpcNSigmaKaon && tpcNSigmaKaon < nSigmaTPCKaonMax && nSigmaTOFKaonMin < tofNSigmaKaon && tofNSigmaKaon < nSigmaTOFKaonMax) {
        /// PID selected kaons
        histograms.fill(HIST("h4ImpPar_Kaon"), pt, impParRPhi, track.eta(), track.phi());
        histograms.fill(HIST("h4ImpParZ_Kaon"), pt, impParZ, track.eta(), track.phi());
        histograms.fill(HIST("hNSigmaTPCKaon_afterPID"), pt, tpcNSigmaKaon);
        histograms.fill(HIST("hNSigmaTOFKaon_afterPID"), pt, tofNSigmaKaon);
      }
      if (isPIDProtonApplied && nSigmaTPCProtonMin < tpcNSigmaProton && tpcNSigmaProton < nSigmaTPCProtonMax && nSigmaTOFProtonMin < tofNSigmaProton && tofNSigmaProton < nSigmaTOFProtonMax) {
        /// PID selected Protons
        histograms.fill(HIST("h4ImpPar_Proton"), pt, impParRPhi, track.eta(), track.phi());
        histograms.fill(HIST("h4ImpParZ_Proton"), pt, impParZ, track.eta(), track.phi());
        histograms.fill(HIST("hNSigmaTPCProton_afterPID"), pt, tpcNSigmaProton);
        histograms.fill(HIST("hNSigmaTOFProton_afterPID"), pt, tofNSigmaProton);
      }
      //}
    }
  }
};

/// QA task for impact parameter distribution monitoring
struct QaImpactParMC {
  /// Input parameters
  //Configurable<int> numberContributorsMin{"numberContributorsMin", 0, "Minimum number of contributors for the primary vertex"};
  Configurable<float> zVtxMaxMC{"zVtxMaxMC", 10.f, "Maximum value for |z_vtx| (MC)"};
  //Configurable<int> keepOnlyGlobalTracksMC{"keepOnlyGlobalTracksMC", 1, "Keep only global tracks or not (MC)"};
  Configurable<float> ptMinMC{"ptMinMC", 0.1f, "Minimum track pt [GeV/c] (MC)"};

  /// Selections with Filter (from o2::framework::expressions)
  // Primary vertex |z_vtx|<XXX cm
  Filter collisionZVtxFilter = nabs(o2::aod::collision::posZ) < zVtxMaxMC;
  // Global tracks
  // with Run 3 data/MC enable the '--isRun3 1' option
  Filter globalTrackFilter = (o2::aod::track::isGlobalTrack == (uint8_t) true); /// filterbit 4 track selections + tight DCA cuts
  // Pt selection
  Filter ptMinFilter = o2::aod::track::pt > ptMinMC;

  /// Histogram registry (from o2::framework)
  HistogramRegistry histograms{"HistogramsImpParQAMC"};

  /// init function - declare and define histograms
  void init(InitContext&)
  {
    // Primary vertex
    const AxisSpec collisionZAxis{100, -20.f, 20.f, "Z [cm]"};
    const AxisSpec collisionNumberContributorAxis{1000, 0, 1000, "Number of contributors"};

    histograms.add("vertexZ", "", kTH1D, {collisionZAxis});
    histograms.add("numberContributors", "", kTH1D, {collisionNumberContributorAxis});
    histograms.add("vertexZ_MCColl", "", kTH1D, {collisionZAxis});

    // tracks
    const AxisSpec trackPtAxis{100, 0.f, 10.f, "#it{p}_{T} [GeV/#it{c}]"};
    const AxisSpec trackEtaAxis{40, -2.f, 2.f, "#it{#eta}"};
    const AxisSpec trackPhiAxis{24, 0.f, TMath::TwoPi(), "#varphi"};
    const AxisSpec trackImpParRPhiAxis{200, -500, 500, "#it{d}_{r#it{#varphi}} [#mum]"};
    const AxisSpec trackImpParZAxis{200, -500, 500, "#it{d}_{z} [#mum]"};
    const AxisSpec trackPDGAxis{3, 0.5f, 3.5f, "species (1: pi, 2: K, 3: p)"};

    histograms.add("pt", "", kTH1D, {trackPtAxis});
    histograms.add("h3ImpPar_PhysPrimary", "", kTHnD, {trackPtAxis, trackImpParRPhiAxis, trackPDGAxis});
    histograms.add("h3ImpParZ_PhysPrimary", "", kTHnD, {trackPtAxis, trackImpParZAxis, trackPDGAxis});
    histograms.add("h3ImpPar_MCvertex_PhysPrimary", "", kTHnD, {trackPtAxis, trackImpParRPhiAxis, trackPDGAxis});
    histograms.add("h3ImpParZ_MCvertex_PhysPrimary", "", kTHnD, {trackPtAxis, trackImpParZAxis, trackPDGAxis});
  }

  //void process(o2::soa::Filtered<o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::McCollisionLabels>>::iterator& collision,
  void process(const o2::soa::Filtered<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>>::iterator& collision,
               const o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksExtended, o2::aod::McTrackLabels>>& tracks,
               const o2::aod::McCollisions&,
               const o2::aod::McParticles& mcParticles) // this Join should ensure to run over all the MC matched tracks
  {
    //o2::dataformats::DCA dca;
    // FIXME: get this from CCDB
    //constexpr float magneticField{5.0};      // in kG
    constexpr float toMicrometers = 10000.f; // Conversion from [cm] to [mum]

    /// trigger selection
    //if (useINT7TriggerMC) {
    //    // from Tutorial/src/multiplicityEventTrackSelection.cxx
    //    if (!collision.alias()[kINT7]) {
    //        return;
    //    }
    //    if (!collision.sel7()) {
    //        return;
    //    }
    //}

    const auto mccollision = collision.mcCollision();

    histograms.fill(HIST("vertexZ"), collision.posZ());
    histograms.fill(HIST("numberContributors"), collision.numContrib());
    histograms.fill(HIST("vertexZ_MCColl"), mccollision.posZ());

    auto PDGtoIndex = [](const int pdg) {
      switch (pdg) {
        case 211: // pion
          return 1;
        case 321: // kaon
          return 2;
        case 2212: // proton
          return 3;
        default: // not identified
          return 0;
      }
    };

    /// loop over tracks
    float impParRPhi = -999.f;
    float impParZ = -999.f;
    //o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      o2::base::GeometryManager::isGeometryLoaded();
      o2::base::GeometryManager::loadGeometry();
      o2::base::Propagator::initFieldFromGRP();
      auto matLUTFile = o2::base::NameConf::getMatLUTFileName();
      if (o2::utils::Str::pathExists(matLUTFile)) {
        auto* lut = o2::base::MatLayerCylSet::loadFromFile(matLUTFile);
        o2::base::Propagator::Instance()->setMatLUT(lut);
      }
    }
    for (const auto& track : tracks) {

      /// Using the Filter instead
      ///if ((keepOnlyGlobalTracksMC) && (!track.isGlobalTrack())) {
      ///  /// not a global track (FB 4 + tight DCA cuts)
      ///  continue;
      ///}

      histograms.fill(HIST("pt"), track.pt());
      const auto mcparticle = track.mcParticle();
      if (MC::isPhysicalPrimary(mcparticle)) {
        impParRPhi = toMicrometers * track.dcaXY(); // from TracksExtended
        impParZ = toMicrometers * track.dcaZ();     // from TracksExtended
        histograms.fill(HIST("h3ImpPar_PhysPrimary"), track.pt(), impParRPhi, PDGtoIndex(std::abs(mcparticle.pdgCode())));
        histograms.fill(HIST("h3ImpParZ_PhysPrimary"), track.pt(), impParZ, PDGtoIndex(std::abs(mcparticle.pdgCode())));
      }

      // propagation to primary vertex for DCA
      // NB: do not use 'track.collisions()' if the o2::aod::Collisions are joined with McCollisionLabels
      // "crude" method not working with Run 3 MC productions
      //if (getTrackParCov(track).propagateToDCA(getPrimaryVertex(/*track.collision()*/ collision), magneticField, &dca, 100.)) {
      std::array<float, 2> dca{1e10f, 1e10f};

      //}
      // Reconstructed vertex
      //if (getTrackParCov(track).propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, o2::base::Propagator::Instance()->getNominalBz(), &dca, 100.)) {
      // MC vertex
      if (getTrackParCov(track).propagateParamToDCA({mccollision.posX(), mccollision.posY(), mccollision.posZ()}, o2::base::Propagator::Instance()->getNominalBz(), &dca, 100.)) {

        /// propagation ok! Retrieve impact parameter
        // PR "Run 3 DCA extraction #187" required - correct calculation of DCAxy of tracks propagated to the PV
        impParRPhi = toMicrometers * dca[0]; //track.dcaXY(); //dca.getY()
        impParZ = toMicrometers * dca[1];    //track.dcaZ(); //dca.getZ()

        /// MC matching - physical primaries
        //const auto mcparticle = track.mcParticle();
        if (MC::isPhysicalPrimary(mcparticle)) {
          histograms.fill(HIST("h3ImpPar_MCvertex_PhysPrimary"), track.pt(), impParRPhi, PDGtoIndex(std::abs(mcparticle.pdgCode())));
          histograms.fill(HIST("h3ImpParZ_MCvertex_PhysPrimary"), track.pt(), impParZ, PDGtoIndex(std::abs(mcparticle.pdgCode())));

          // retrieve information about the mother mcparticle
          //const auto indexMother0 = mcparticle.mother0Id();
          //const auto indexMother1 = mcparticle.mother1Id();
          //std::cout << std::endl << "=========================================" << std::endl;
          //std::cout << "   indexMother0 " << indexMother0 << ", indexMother1 " << indexMother1 << std::endl;
          //if(mcparticle.has_mother0()) {
          //  std::cout << "===> mcparticle PDG " << mcparticle.pdgCode() << std::endl;
          //  std::cout << "===> mothers:\n\t";
          //  std::cout << "pdg mother 0: " << mcParticles.iteratorAt(indexMother0).pdgCode();
          //  if(mcparticle.has_mother1()){
          //    for(auto index=(indexMother0+1); index<=indexMother1; index++) {
          //      std::cout << " " << mcParticles.iteratorAt(index).pdgCode();
          //    }
          //  }
          //  std::cout << std::endl;
          //}
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w{
    adaptAnalysisTask<QaImpactPar>(cfgc, TaskName{"qa-impact-par"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    w.push_back(adaptAnalysisTask<QaImpactParMC>(cfgc, TaskName{"qa-impact-par-mc"}));
  }
  return w;
}