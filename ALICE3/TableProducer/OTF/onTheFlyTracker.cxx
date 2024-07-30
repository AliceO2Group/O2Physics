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

/// \file onTheFlyTracker.cxx
///
/// \brief LUT-based on-the-fly analysis task-level tracking
///
/// This task allows for the calculation of aod::collisions and aod::Tracks in a synthetic manner,
/// smearing MC particles with very configurable settings. This will allow for the usage of
/// custom LUTs (obtained through separate studies) and the subsequent estimate of the performance
/// of a future detector even in very statistics-hungry analyses.
///
/// \author David Dobrigkeit Chinellato, UNICAMP
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>, UniBo
/// \author Roberto Preghenella preghenella@bo.infn.it
///

#include <utility>

#include <TGeoGlobalMagField.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include <TPDGCode.h>
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "SimulationDataFormat/InteractionSampler.h"
#include "Field/MagneticField.h"

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"

using namespace o2;
using namespace o2::framework;

struct OnTheFlyTracker {
  Produces<aod::Collisions> collisions;
  Produces<aod::McCollisionLabels> collLabels;
  Produces<aod::StoredTracks> tracksPar;
  Produces<aod::TracksExtension> tracksParExtension;
  Produces<aod::StoredTracksCov> tracksParCov;
  Produces<aod::TracksCovExtension> tracksParCovExtension;
  Produces<aod::McTrackLabels> tracksLabels;
  Produces<aod::TracksDCA> tracksDCA;
  Produces<aod::TracksDCACov> tracksDCACov;
  Produces<aod::CollisionsAlice3> collisionsAlice3;
  Produces<aod::TracksAlice3> TracksAlice3;

  // optionally produced, empty (to be tuned later)
  Produces<aod::StoredTracksExtra> tracksExtra; // base table, extend later
  Produces<aod::TrackSelection> trackSelection;
  Produces<aod::TrackSelectionExtension> trackSelectionExtension;

  Configurable<int> seed{"seed", 0, "TGenPhaseSpace seed"};
  Configurable<float> magneticField{"magneticField", 20.0f, "magnetic field in kG"};
  Configurable<float> maxEta{"maxEta", 1.5, "maximum eta to consider viable"};
  Configurable<float> multEtaRange{"multEtaRange", 0.8, "eta range to compute the multiplicity"};
  Configurable<float> minPt{"minPt", 0.1, "minimum pt to consider viable"};
  Configurable<bool> enableLUT{"enableLUT", false, "Enable track smearing"};
  Configurable<bool> enableNucleiSmearing{"enableNucleiSmearing", false, "Enable smearing of nuclei"};
  Configurable<bool> enablePrimaryVertexing{"enablePrimaryVertexing", true, "Enable primary vertexing"};
  Configurable<bool> interpolateLutEfficiencyVsNch{"interpolateLutEfficiencyVsNch", true, "interpolate LUT efficiency as f(Nch)"};
  Configurable<bool> treatXi{"treatXi", false, "Manually decay Xi^{-} and fill tables with daughters"};

  Configurable<bool> populateTracksDCA{"populateTracksDCA", true, "populate TracksDCA table"};
  Configurable<bool> populateTracksDCACov{"populateTracksDCACov", false, "populate TracksDCACov table"};
  Configurable<bool> populateTracksExtra{"populateTracksExtra", false, "populate TracksExtra table (legacy)"};
  Configurable<bool> populateTrackSelection{"populateTrackSelection", false, "populate TrackSelection table (legacy)"};

  Configurable<bool> processUnreconstructedTracks{"processUnreconstructedTracks", false, "process (smear) unreco-ed tracks"};
  Configurable<bool> doExtraQA{"doExtraQA", false, "do extra 2D QA plots"};
  Configurable<bool> extraQAwithoutDecayDaughters{"extraQAwithoutDecayDaughters", false, "remove decay daughters from qa plots (yes/no)"};
  Configurable<bool> doXiQA{"doXiQA", false, "QA plots for when treating Xi"};

  Configurable<std::string> lutEl{"lutEl", "lutCovm.el.dat", "LUT for electrons"};
  Configurable<std::string> lutMu{"lutMu", "lutCovm.mu.dat", "LUT for muons"};
  Configurable<std::string> lutPi{"lutPi", "lutCovm.pi.dat", "LUT for pions"};
  Configurable<std::string> lutKa{"lutKa", "lutCovm.ka.dat", "LUT for kaons"};
  Configurable<std::string> lutPr{"lutPr", "lutCovm.pr.dat", "LUT for protons"};
  Configurable<std::string> lutDe{"lutDe", "lutCovm.de.dat", "LUT for deuterons"};
  Configurable<std::string> lutTr{"lutTr", "lutCovm.tr.dat", "LUT for tritons"};
  Configurable<std::string> lutHe3{"lutHe3", "lutCovm.he3.dat", "LUT for Helium-3"};

  Configurable<std::string> lutPi0{"lutPi0", "lutCovm.pi.20kG.rmin20.geometry_v0.dat", "LUT for pions without layer 0"};
  Configurable<std::string> lutPi1{"lutPi1", "lutCovm.pi.20kG.rmin20.geometry_v1.dat", "LUT for pions without layer 1"};
  Configurable<std::string> lutPi2{"lutPi2", "lutCovm.pi.20kG.rmin20.geometry_v2.dat", "LUT for pions without layer 2"};
  Configurable<std::string> lutPi3{"lutPi3", "lutCovm.pi.20kG.rmin20.geometry_v3.dat", "LUT for pions without layer 3"};
  Configurable<std::string> lutPi4{"lutPi4", "lutCovm.pi.20kG.rmin20.geometry_v4.dat", "LUT for pions without layer 4"};
  Configurable<std::string> lutPi5{"lutPi5", "lutCovm.pi.20kG.rmin20.geometry_v5.dat", "LUT for pions without layer 5"};
  Configurable<std::string> lutPr0{"lutPr0", "lutCovm.pr.20kG.rmin20.geometry_v0.dat", "LUT for protons without layer 0"};
  Configurable<std::string> lutPr1{"lutPr1", "lutCovm.pr.20kG.rmin20.geometry_v1.dat", "LUT for protons without layer 1"};
  Configurable<std::string> lutPr2{"lutPr2", "lutCovm.pr.20kG.rmin20.geometry_v2.dat", "LUT for protons without layer 2"};
  Configurable<std::string> lutPr3{"lutPr3", "lutCovm.pr.20kG.rmin20.geometry_v3.dat", "LUT for protons without layer 3"};
  Configurable<std::string> lutPr4{"lutPr4", "lutCovm.pr.20kG.rmin20.geometry_v4.dat", "LUT for protons without layer 4"};
  Configurable<std::string> lutPr5{"lutPr5", "lutCovm.pr.20kG.rmin20.geometry_v5.dat", "LUT for protons without layer 5"};

  ConfigurableAxis axisMomentum{"axisMomentum", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "#it{p} (GeV/#it{c})"};
  ConfigurableAxis axisNVertices{"axisNVertices", {20, -0.5, 19.5}, "N_{vertices}"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {100, -0.5, 99.5}, "N_{contributors}"};
  ConfigurableAxis axisVertexZ{"axisVertexZ", {40, -20, 20}, "vertex Z (cm)"};
  ConfigurableAxis axisDCA{"axisDCA", {400, -200, 200}, "DCA (#mum)"};
  ConfigurableAxis axisX{"axisX", {250, -50, 200}, "track X (cm)"};
  ConfigurableAxis axisRadius{"axisRadius", {55, 0.01, 100}, "decay radius"};

  using PVertex = o2::dataformats::PrimaryVertex;

  // Class to hold the track information for the O2 vertexing
  class TrackAlice3 : public o2::track::TrackParCov
  {
    using TimeEst = o2::dataformats::TimeStampWithError<float, float>;

   public:
    TrackAlice3() = default;
    ~TrackAlice3() = default;
    TrackAlice3(const TrackAlice3& src) = default;
    TrackAlice3(const o2::track::TrackParCov& src, const int64_t label, const float t = 0, const float te = 1, bool decayDauInput = false) : o2::track::TrackParCov(src), mcLabel{label}, timeEst{t, te}, isDecayDau(decayDauInput) {}
    const TimeEst& getTimeMUS() const { return timeEst; }
    int64_t mcLabel;
    TimeEst timeEst; ///< time estimate in ns
    bool isDecayDau;
  };

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer
  o2::delphes::DelphesO2TrackSmearer mSmearer;

  // Xi daughters uses different smearers
  o2::delphes::DelphesO2TrackSmearer mSmearer0;
  o2::delphes::DelphesO2TrackSmearer mSmearer1;
  o2::delphes::DelphesO2TrackSmearer mSmearer2;
  o2::delphes::DelphesO2TrackSmearer mSmearer3;
  o2::delphes::DelphesO2TrackSmearer mSmearer4;
  o2::delphes::DelphesO2TrackSmearer mSmearer5;

  // For processing and vertexing
  std::vector<TrackAlice3> tracksAlice3;
  std::vector<TrackAlice3> ghostTracksAlice3;
  std::vector<o2::InteractionRecord> bcData;
  o2::steer::InteractionSampler irSampler;
  o2::vertexing::PVertexer vertexer;

  void init(o2::framework::InitContext&)
  {
    if (enableLUT) {
      std::map<int, const char*> mapPdgLut;
      const char* lutElChar = lutEl->c_str();
      const char* lutMuChar = lutMu->c_str();
      const char* lutPiChar = lutPi->c_str();
      const char* lutKaChar = lutKa->c_str();
      const char* lutPrChar = lutPr->c_str();

      LOGF(info, "Will load electron lut file ..: %s", lutElChar);
      LOGF(info, "Will load muon lut file ......: %s", lutMuChar);
      LOGF(info, "Will load pion lut file ......: %s", lutPiChar);
      LOGF(info, "Will load kaon lut file ......: %s", lutKaChar);
      LOGF(info, "Will load proton lut file ....: %s", lutPrChar);

      mapPdgLut.insert(std::make_pair(11, lutElChar));
      mapPdgLut.insert(std::make_pair(13, lutMuChar));
      mapPdgLut.insert(std::make_pair(211, lutPiChar));
      mapPdgLut.insert(std::make_pair(321, lutKaChar));
      mapPdgLut.insert(std::make_pair(2212, lutPrChar));

      if (enableNucleiSmearing) {
        const char* lutDeChar = lutDe->c_str();
        const char* lutTrChar = lutTr->c_str();
        const char* lutHe3Char = lutHe3->c_str();
        mapPdgLut.insert(std::make_pair(1000010020, lutDeChar));
        mapPdgLut.insert(std::make_pair(1000010030, lutTrChar));
        mapPdgLut.insert(std::make_pair(1000020030, lutHe3Char));
      }
      for (auto e : mapPdgLut) {
        if (!mSmearer.loadTable(e.first, e.second)) {
          LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
        }
      }
      // interpolate efficiencies if requested to do so
      mSmearer.interpolateEfficiency(static_cast<bool>(interpolateLutEfficiencyVsNch));

      // smear un-reco'ed tracks if asked to do so
      mSmearer.skipUnreconstructed(static_cast<bool>(!processUnreconstructedTracks));

      if (treatXi) {
        std::map<int, const char*> mapPdgLut0;
        const char* lutPiChar0 = lutPi0->c_str();
        const char* lutPrChar0 = lutPr0->c_str();
        LOGF(info, "Load more pion lut files .....: %s", lutPiChar0);
        LOGF(info, "Load more proton lut files ...: %s", lutPrChar0);
        mapPdgLut0.insert(std::make_pair(211, lutPiChar0));
        mapPdgLut0.insert(std::make_pair(2212, lutPrChar0));

        std::map<int, const char*> mapPdgLut1;
        const char* lutPiChar1 = lutPi1->c_str();
        const char* lutPrChar1 = lutPr1->c_str();
        LOGF(info, "Load more pion lut files .....: %s", lutPiChar1);
        LOGF(info, "Load more proton lut files ...: %s", lutPrChar1);
        mapPdgLut1.insert(std::make_pair(211, lutPiChar1));
        mapPdgLut1.insert(std::make_pair(2212, lutPrChar1));

        std::map<int, const char*> mapPdgLut2;
        const char* lutPiChar2 = lutPi2->c_str();
        const char* lutPrChar2 = lutPr2->c_str();
        LOGF(info, "Load more pion lut files .....: %s", lutPiChar2);
        LOGF(info, "Load more proton lut files ...: %s", lutPrChar2);
        mapPdgLut2.insert(std::make_pair(211, lutPiChar2));
        mapPdgLut2.insert(std::make_pair(2212, lutPrChar2));

        std::map<int, const char*> mapPdgLut3;
        const char* lutPiChar3 = lutPi3->c_str();
        const char* lutPrChar3 = lutPr3->c_str();
        LOGF(info, "Load more pion lut files .....: %s", lutPiChar3);
        LOGF(info, "Load more proton lut files ...: %s", lutPrChar3);
        mapPdgLut3.insert(std::make_pair(211, lutPiChar3));
        mapPdgLut3.insert(std::make_pair(2212, lutPrChar3));

        std::map<int, const char*> mapPdgLut4;
        const char* lutPiChar4 = lutPi4->c_str();
        const char* lutPrChar4 = lutPr4->c_str();
        LOGF(info, "Load more pion lut files .....: %s", lutPiChar4);
        LOGF(info, "Load more proton lut files ...: %s", lutPrChar4);
        mapPdgLut4.insert(std::make_pair(211, lutPiChar4));
        mapPdgLut4.insert(std::make_pair(2212, lutPrChar4));

        std::map<int, const char*> mapPdgLut5;
        const char* lutPiChar5 = lutPi5->c_str();
        const char* lutPrChar5 = lutPr5->c_str();
        LOGF(info, "Load more pion lut files .....: %s", lutPiChar5);
        LOGF(info, "Load more proton lut files ...: %s", lutPrChar5);
        mapPdgLut5.insert(std::make_pair(211, lutPiChar5));
        mapPdgLut5.insert(std::make_pair(2212, lutPrChar5));

        for (auto e : mapPdgLut0) {
          if (!mSmearer0.loadTable(e.first, e.second)) {
            LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
          }
        }
        for (auto e : mapPdgLut1) {
          if (!mSmearer1.loadTable(e.first, e.second)) {
            LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
          }
        }
        for (auto e : mapPdgLut2) {
          if (!mSmearer2.loadTable(e.first, e.second)) {
            LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
          }
        }
        for (auto e : mapPdgLut3) {
          if (!mSmearer3.loadTable(e.first, e.second)) {
            LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
          }
        }
        for (auto e : mapPdgLut4) {
          if (!mSmearer4.loadTable(e.first, e.second)) {
            LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
          }
        }
        for (auto e : mapPdgLut5) {
          if (!mSmearer5.loadTable(e.first, e.second)) {
            LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
          }
        }

        // interpolate efficiencies if requested to do so
        mSmearer0.interpolateEfficiency(static_cast<bool>(interpolateLutEfficiencyVsNch));
        mSmearer1.interpolateEfficiency(static_cast<bool>(interpolateLutEfficiencyVsNch));
        mSmearer2.interpolateEfficiency(static_cast<bool>(interpolateLutEfficiencyVsNch));
        mSmearer3.interpolateEfficiency(static_cast<bool>(interpolateLutEfficiencyVsNch));
        mSmearer4.interpolateEfficiency(static_cast<bool>(interpolateLutEfficiencyVsNch));
        mSmearer5.interpolateEfficiency(static_cast<bool>(interpolateLutEfficiencyVsNch));

        // smear un-reco'ed tracks if asked to do so
        mSmearer0.skipUnreconstructed(static_cast<bool>(!processUnreconstructedTracks));
        mSmearer1.skipUnreconstructed(static_cast<bool>(!processUnreconstructedTracks));
        mSmearer2.skipUnreconstructed(static_cast<bool>(!processUnreconstructedTracks));
        mSmearer3.skipUnreconstructed(static_cast<bool>(!processUnreconstructedTracks));
        mSmearer4.skipUnreconstructed(static_cast<bool>(!processUnreconstructedTracks));
        mSmearer5.skipUnreconstructed(static_cast<bool>(!processUnreconstructedTracks));
      }
    }

    // Basic QA
    histos.add("hPtGenerated", "hPtGenerated", kTH1F, {axisMomentum});
    histos.add("hPtGeneratedEl", "hPtGeneratedEl", kTH1F, {axisMomentum});
    histos.add("hPtGeneratedPi", "hPtGeneratedPi", kTH1F, {axisMomentum});
    histos.add("hPtGeneratedKa", "hPtGeneratedKa", kTH1F, {axisMomentum});
    histos.add("hPtGeneratedPr", "hPtGeneratedPr", kTH1F, {axisMomentum});
    histos.add("hPtReconstructed", "hPtReconstructed", kTH1F, {axisMomentum});
    histos.add("hPtReconstructedEl", "hPtReconstructedEl", kTH1F, {axisMomentum});
    histos.add("hPtReconstructedPi", "hPtReconstructedPi", kTH1F, {axisMomentum});
    histos.add("hPtReconstructedKa", "hPtReconstructedKa", kTH1F, {axisMomentum});
    histos.add("hPtReconstructedPr", "hPtReconstructedPr", kTH1F, {axisMomentum});

    // Collision QA
    histos.add("hPVz", "hPVz", kTH1F, {axisVertexZ});
    histos.add("hLUTMultiplicity", "hLUTMultiplicity", kTH1F, {axisMultiplicity});
    histos.add("hSimMultiplicity", "hSimMultiplicity", kTH1F, {axisMultiplicity});
    histos.add("hRecoMultiplicity", "hRecoMultiplicity", kTH1F, {axisMultiplicity});

    if (doExtraQA) {
      histos.add("h2dVerticesVsContributors", "h2dVerticesVsContributors", kTH2F, {axisMultiplicity, axisNVertices});
      histos.add("hRecoVsSimMultiplicity", "hRecoVsSimMultiplicity", kTH2F, {axisMultiplicity, axisMultiplicity});
      histos.add("h2dDCAxy", "h2dDCAxy", kTH2F, {axisMomentum, axisDCA});

      histos.add("hSimTrackX", "hSimTrackX", kTH1F, {axisX});
      histos.add("hRecoTrackX", "hRecoTrackX", kTH1F, {axisX});
      histos.add("hTrackXatDCA", "hTrackXatDCA", kTH1F, {axisX});
    }

    if (doXiQA) {
      histos.add("hGenXi", "hGenXi", kTH2F, {axisRadius, axisMomentum});
      histos.add("hRecoXi", "hRecoXi", kTH2F, {axisRadius, axisMomentum});

      histos.add("hGenPiFromXi", "hGenPiFromXi", kTH2F, {axisRadius, axisMomentum});
      histos.add("hGenPiFromL0", "hGenPiFromL0", kTH2F, {axisRadius, axisMomentum});
      histos.add("hGenPrFromL0", "hGenPrFromL0", kTH2F, {axisRadius, axisMomentum});
      histos.add("hRecoPiFromXi", "hRecoPiFromXi", kTH2F, {axisRadius, axisMomentum});
      histos.add("hRecoPiFromL0", "hRecoPiFromL0", kTH2F, {axisRadius, axisMomentum});
      histos.add("hRecoPrFromL0", "hRecoPrFromL0", kTH2F, {axisRadius, axisMomentum});
    }

    LOGF(info, "Initializing magnetic field to value: %.3f kG", static_cast<float>(magneticField));
    o2::parameters::GRPMagField grpmag;
    grpmag.setFieldUniformity(true);
    grpmag.setL3Current(30000.f * (magneticField / 5.0f));
    auto field = grpmag.getNominalL3Field();
    o2::base::Propagator::initFieldFromGRP(&grpmag);

    auto fieldInstance = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (!fieldInstance) {
      LOGF(fatal, "Failed to set up magnetic field! Stopping now!");
    }

    // Cross-check
    LOGF(info, "Check field at (0, 0, 0): %.1f kG, nominal: %.1f", static_cast<float>(fieldInstance->GetBz(0, 0, 0)), static_cast<float>(field));

    LOGF(info, "Initializing empty material cylinder LUT - could be better in the future");
    o2::base::MatLayerCylSet* lut = new o2::base::MatLayerCylSet();
    lut->addLayer(200, 200.1, 2, 1.0f, 100.0f);
    LOGF(info, "MatLayerCylSet::optimizePhiSlices()");
    lut->optimizePhiSlices();
    LOGF(info, "Setting lut now...");
    o2::base::Propagator::Instance()->setMatLUT(lut);

    irSampler.setInteractionRate(100);
    irSampler.setFirstIR(o2::InteractionRecord(0, 0));
    irSampler.init();

    vertexer.setValidateWithIR(kFALSE);
    vertexer.setBunchFilling(irSampler.getBunchFilling());
    vertexer.init();
  }

  /// Function to decay the xi
  /// \param particle the particle to decay
  /// \param decayDaughters the address of resulting daughters
  /// \param xiDecayVertex the address of the xi decay vertex
  /// \param l0DecayVertex the address of the l0 decay vertex
  template <typename McParticleType>
  void decayParticle(McParticleType particle, std::vector<TLorentzVector>& decayDaughters, std::vector<double>& xiDecayVertex, std::vector<double>& l0DecayVertex)
  {
    std::vector<double> xiVelocity(3);
    std::vector<double> l0Velocity(3);
    std::vector<double> xiMomentum = {particle.px(), particle.py(), particle.pz()};
    std::vector<double> xiProductionVertex = {particle.vx(), particle.vy(), particle.vz()};
    double bz = magneticField * (-0.1); // To tesla (sign?!)
    TRandom3 rand;
    rand.SetSeed(seed);
    double u = rand.Uniform(0, 1);
    for (int i = 0; i < 3; i++) {
      xiVelocity[i] = xiMomentum[i] / particle.e();
    }
    double xi_v_tot = particle.p() / particle.e();
    double xi_ctau = 4.91 / 100; // xi
    double xi_charge = -1.6022e-19;
    double speedOfLight = 3e+8;
    double xi_v_xy = speedOfLight * sqrt(xiVelocity[0] * xiVelocity[0] + xiVelocity[1] * xiVelocity[1]);
    double xi_lifetime = (-xi_ctau * log(1 - u)) / xi_v_tot;
    double xi_kgMass = 1.32171 * 1.78e-27;
    double xi_r = (xi_kgMass * xi_v_xy) / (xi_charge * bz);
    double xi_theta = ((xi_lifetime / speedOfLight) * xi_charge * bz) / xi_kgMass;
    double xi_phi = TMath::ATan2(xiVelocity[1], xiVelocity[0]);
    xiDecayVertex.push_back(xiProductionVertex[0] + xi_r * (std::sin(xi_theta) * std::cos(xi_phi) - (1 - std::cos(xi_theta)) * std::sin(xi_phi)));
    xiDecayVertex.push_back(xiProductionVertex[1] + xi_r * ((1 - std::cos(xi_theta)) * std::cos(xi_phi) + std::sin(xi_theta) * std::sin(xi_phi)));
    xiDecayVertex.push_back(xiProductionVertex[2] + xi_lifetime * xiVelocity[2]);
    std::vector<double> xiDaughters = {1.115, 0.139};
    double px = (xi_v_xy * std::cos(xi_theta + xi_phi) * particle.e()) / speedOfLight;
    double py = (xi_v_xy * std::sin(xi_theta + xi_phi) * particle.e()) / speedOfLight;

    TLorentzVector xi(px, py, particle.pz(), particle.e());
    TGenPhaseSpace xiDecay;
    xiDecay.SetDecay(xi, 2, xiDaughters.data());
    xiDecay.Generate();
    decayDaughters.push_back(*xiDecay.GetDecay(1));

    double l0_ctau = 7.89 / 100; // lambda
    double l0_v_tot = xiDecay.GetDecay(0)->P() / xiDecay.GetDecay(0)->E();
    std::vector<double> l0Daughters = {0.139, 0.938};
    l0Velocity[0] = xiDecay.GetDecay(0)->Px() / xiDecay.GetDecay(0)->E();
    l0Velocity[1] = xiDecay.GetDecay(0)->Py() / xiDecay.GetDecay(0)->E();
    l0Velocity[2] = xiDecay.GetDecay(0)->Pz() / xiDecay.GetDecay(0)->E();
    double l0_lifetime = (-l0_ctau * log(1 - u)) / l0_v_tot;
    for (int i = 0; i < 3; i++) {
      l0DecayVertex.push_back(xiDecayVertex[i] + l0_lifetime * l0Velocity[i]);
    }

    TLorentzVector l0 = *xiDecay.GetDecay(0);
    TGenPhaseSpace l0Decay;
    l0Decay.SetDecay(l0, 2, l0Daughters.data());
    l0Decay.Generate();
    decayDaughters.push_back(*l0Decay.GetDecay(0));
    decayDaughters.push_back(*l0Decay.GetDecay(1));
  }

  /// Function to convert a TLorentzVector into a perfect Track
  /// \param pdgCode particle pdg
  /// \param particle the particle to convert (TLorentzVector)
  /// \param productionVertex where the particle was produced
  /// \param o2track the address of the resulting TrackParCov
  void convertTLorentzVectorToO2Track(int pdgCode, TLorentzVector particle, std::vector<double> productionVertex, o2::track::TrackParCov& o2track)
  {
    auto pdgInfo = pdgDB->GetParticle(pdgCode);
    int charge = 0;
    if (pdgInfo != nullptr) {
      charge = pdgInfo->Charge() / 3;
    }
    std::array<float, 5> params;
    std::array<float, 15> covm = {0.};
    float s, c, x;
    o2::math_utils::sincos(static_cast<float> particle.Phi(), s, c);
    o2::math_utils::rotateZInv(static_cast<float> productionVertex[0], static_cast<float> productionVertex[1], x, params[0], s, c);
    params[1] = static_cast<float> productionVertex[2];
    params[2] = 0;
    auto theta = 2. * std::atan(std::exp(-particle.PseudoRapidity()));
    params[3] = 1. / std::tan(theta);
    params[4] = charge / particle.Pt();

    // Initialize TrackParCov in-place
    new (&o2track)(o2::track::TrackParCov)(x, particle.Phi(), params, covm);
  }

  /// Function to convert a McParticle into a perfect Track
  /// \param particle the particle to convert (mcParticle)
  /// \param o2track the address of the resulting TrackParCov
  template <typename McParticleType>
  void convertMCParticleToO2Track(McParticleType& particle, o2::track::TrackParCov& o2track)
  {
    auto pdgInfo = pdgDB->GetParticle(particle.pdgCode());
    int charge = 0;
    if (pdgInfo != nullptr) {
      charge = pdgInfo->Charge() / 3;
    }
    std::array<float, 5> params;
    std::array<float, 15> covm = {0.};
    float s, c, x;
    o2::math_utils::sincos(particle.phi(), s, c);
    o2::math_utils::rotateZInv(particle.vx(), particle.vy(), x, params[0], s, c);
    params[1] = particle.vz();
    params[2] = 0.; // since alpha = phi
    auto theta = 2. * std::atan(std::exp(-particle.eta()));
    params[3] = 1. / std::tan(theta);
    params[4] = charge / particle.pt();

    // Initialize TrackParCov in-place
    new (&o2track)(o2::track::TrackParCov)(x, particle.phi(), params, covm);
  }

  float dNdEta = 0.f; // Charged particle multiplicity to use in the efficiency evaluation
  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    tracksAlice3.clear();
    ghostTracksAlice3.clear();
    bcData.clear();

    o2::dataformats::DCA dcaInfo;
    o2::dataformats::VertexBase vtx;

    // generate collision time
    auto ir = irSampler.generateCollisionTime();

    // First we compute the number of charged particles in the event
    dNdEta = 0.f;
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.eta()) > multEtaRange) {
        continue;
      }
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      const auto pdg = std::abs(mcParticle.pdgCode());
      if (pdg != kElectron && pdg != kMuonMinus && pdg != kPiPlus && pdg != kKPlus && pdg != kProton) {
        if (!treatXi) {
          continue;
        } else if (pdg != 3312) {
          continue;
        }
      }
      const auto& pdgInfo = pdgDB->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        LOG(warning) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        continue;
      }
      if (pdgInfo->Charge() == 0) {
        continue;
      }
      dNdEta += 1.f;
    }

    dNdEta /= (multEtaRange * 2.0f);
    uint32_t multiplicityCounter = 0;
    histos.fill(HIST("hLUTMultiplicity"), dNdEta);
    gRandom->SetSeed(seed);

    for (const auto& mcParticle : mcParticles) {
      double xiDecayRadius2D = 0;
      double l0DecayRadius2D = 0;
      std::vector<TLorentzVector> decayProducts;
      std::vector<double> xiDecayVertex, l0DecayVertex;
      std::vector<double> layers = {0.50, 1.20, 2.50, 3.75, 7.00, 12.0, 20.0};
      if (treatXi) {
        if (mcParticle.pdgCode() == 3312) {
          decayParticle(mcParticle, decayProducts, xiDecayVertex, l0DecayVertex);
          xiDecayRadius2D = sqrt(xiDecayVertex[0] * xiDecayVertex[0] + xiDecayVertex[1] * xiDecayVertex[1]) * 100;
          l0DecayRadius2D = sqrt(l0DecayVertex[0] * l0DecayVertex[0] + l0DecayVertex[1] * l0DecayVertex[1]) * 100;
        }
      }

      const auto pdg = std::abs(mcParticle.pdgCode());
      if (!mcParticle.isPhysicalPrimary()) {
        if (!treatXi) {
          continue;
        } else if (pdg != 3312) {
          continue;
        }
      }
      if (pdg != kElectron && pdg != kMuonMinus && pdg != kPiPlus && pdg != kKPlus && pdg != kProton) {
        if (!treatXi) {
          continue;
        } else if (pdg != 3312) {
          continue;
        }
      }
      if (std::fabs(mcParticle.eta()) > maxEta) {
        continue;
      }

      histos.fill(HIST("hPtGenerated"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 11)
        histos.fill(HIST("hPtGeneratedEl"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 211)
        histos.fill(HIST("hPtGeneratedPi"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 321)
        histos.fill(HIST("hPtGeneratedKa"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 2212)
        histos.fill(HIST("hPtGeneratedPr"), mcParticle.pt());

      if (doXiQA && mcParticle.pdgCode() == 3312) {
        histos.fill(HIST("hGenXi"), xiDecayRadius2D, mcParticle.pt());
        histos.fill(HIST("hGenPiFromXi"), xiDecayRadius2D, decayProducts[0].Pt());
        histos.fill(HIST("hGenPiFromL0"), l0DecayRadius2D, decayProducts[1].Pt());
        histos.fill(HIST("hGenPrFromL0"), l0DecayRadius2D, decayProducts[2].Pt());
      }

      if (mcParticle.pt() < minPt) {
        continue;
      }

      bool isDecayDaughter = false;
      if (mcParticle.getProcess() == 4)
        isDecayDaughter = true;

      multiplicityCounter++;
      const float t = (ir.timeInBCNS + gRandom->Gaus(0., 100.)) * 1e-3;
      std::vector<o2::track::TrackParCov> xiDaughterTrackParCovs(3);
      std::vector<bool> isReco(3);
      std::vector<o2::delphes::DelphesO2TrackSmearer> smearer = {mSmearer0, mSmearer1, mSmearer2, mSmearer3, mSmearer4, mSmearer5};
      if (treatXi && mcParticle.pdgCode() == 3312) {
        if (xiDecayRadius2D > 20) {
          continue;
        }

        convertTLorentzVectorToO2Track(-211, decayProducts[0], xiDecayVertex, xiDaughterTrackParCovs[0]);
        convertTLorentzVectorToO2Track(-211, decayProducts[1], l0DecayVertex, xiDaughterTrackParCovs[1]);
        convertTLorentzVectorToO2Track(2212, decayProducts[2], l0DecayVertex, xiDaughterTrackParCovs[2]);

        // Map daughter to smearer
        int firstSmearerIndex = -1;
        int secondSmearerIndex = -1;
        for (unsigned i = 0; i < layers.size(); i++) {
          if (xiDecayRadius2D > layers[i]) {
            firstSmearerIndex = i;
          }
          if (l0DecayRadius2D > layers[i]) {
            secondSmearerIndex = i;
          }
        }
        if (firstSmearerIndex > 5) {
          isReco[0] = false;
        } else if (firstSmearerIndex == -1) {
          isReco[0] = mSmearer.smearTrack(xiDaughterTrackParCovs[0], 211, dNdEta);
        } else {
          isReco[0] = smearer[firstSmearerIndex].smearTrack(xiDaughterTrackParCovs[0], 211, dNdEta);
        }
        if (secondSmearerIndex > 5) {
          isReco[1] = false;
          isReco[2] = false;
        } else if (secondSmearerIndex == -1) {
          isReco[1] = mSmearer.smearTrack(xiDaughterTrackParCovs[1], 211, dNdEta);
          isReco[2] = mSmearer.smearTrack(xiDaughterTrackParCovs[2], 2212, dNdEta);
        } else {
          isReco[1] = smearer[secondSmearerIndex].smearTrack(xiDaughterTrackParCovs[1], 211, dNdEta);
          isReco[2] = smearer[secondSmearerIndex].smearTrack(xiDaughterTrackParCovs[2], 2212, dNdEta);
        }
        for (int i = 0; i < 3; i++) {
          if (decayProducts[i].Pt() < minPt) {
            isReco[i] = false;
          }
          if (!isReco[i] && !processUnreconstructedTracks) {
            continue;
          }
          if (TMath::IsNaN(xiDaughterTrackParCovs[i].getZ())) {
            continue;
          }
          if (isReco[i]) {
            tracksAlice3.push_back(TrackAlice3{xiDaughterTrackParCovs[i], mcParticle.globalIndex(), t, 100.f * 1e-3, true});
          } else {
            ghostTracksAlice3.push_back(TrackAlice3{xiDaughterTrackParCovs[i], mcParticle.globalIndex(), t, 100.f * 1e-3, true});
          }
        }

        if (doXiQA && mcParticle.pdgCode() == 3312) {
          if (isReco[0] && isReco[1] && isReco[2])
            histos.fill(HIST("hRecoXi"), xiDecayRadius2D, mcParticle.pt());
          if (isReco[0])
            histos.fill(HIST("hRecoPiFromXi"), xiDecayRadius2D, decayProducts[0].Pt());
          if (isReco[1])
            histos.fill(HIST("hRecoPiFromL0"), l0DecayRadius2D, decayProducts[1].Pt());
          if (isReco[2])
            histos.fill(HIST("hRecoPrFromL0"), l0DecayRadius2D, decayProducts[2].Pt());
        }
        continue; // Not filling the tables with the xi itself
      }

      o2::track::TrackParCov trackParCov;
      convertMCParticleToO2Track(mcParticle, trackParCov);

      if (doExtraQA) {
        histos.fill(HIST("hSimTrackX"), trackParCov.getX());
      }

      bool reconstructed = mSmearer.smearTrack(trackParCov, mcParticle.pdgCode(), dNdEta);
      if (!reconstructed && !processUnreconstructedTracks) {
        continue;
      }
      if (TMath::IsNaN(trackParCov.getZ())) {
        // capture rare smearing mistakes / corrupted tracks
        continue;
      }

      // Base QA (note: reco pT here)
      histos.fill(HIST("hPtReconstructed"), trackParCov.getPt());
      if (TMath::Abs(mcParticle.pdgCode()) == 11)
        histos.fill(HIST("hPtReconstructedEl"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 211)
        histos.fill(HIST("hPtReconstructedPi"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 321)
        histos.fill(HIST("hPtReconstructedKa"), mcParticle.pt());
      if (TMath::Abs(mcParticle.pdgCode()) == 2212)
        histos.fill(HIST("hPtReconstructedPr"), mcParticle.pt());

      if (doExtraQA) {
        histos.fill(HIST("hRecoTrackX"), trackParCov.getX());
      }

      // populate vector with track if we reco-ed it
      if (reconstructed) {
        tracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), t, 100.f * 1e-3, isDecayDaughter});
      } else {
        ghostTracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), t, 100.f * 1e-3, isDecayDaughter});
      }
    }

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // Calculate primary vertex with tracks from this collision
    // data preparation
    o2::vertexing::PVertex primaryVertex;

    if (enablePrimaryVertexing) {
      std::vector<o2::MCCompLabel> lblTracks;
      std::vector<o2::vertexing::PVertex> vertices;
      std::vector<o2::vertexing::GIndex> vertexTrackIDs;
      std::vector<o2::vertexing::V2TRef> v2tRefs;
      std::vector<o2::MCEventLabel> lblVtx;
      lblVtx.emplace_back(mcCollision.globalIndex(), 1);
      std::vector<o2::dataformats::GlobalTrackID> idxVec; // store IDs

      idxVec.reserve(tracksAlice3.size());
      for (unsigned i = 0; i < tracksAlice3.size(); i++) {
        lblTracks.emplace_back(tracksAlice3[i].mcLabel, mcCollision.globalIndex(), 1, false);
        idxVec.emplace_back(i, o2::dataformats::GlobalTrackID::ITS); // let's say ITS
      }

      // Calculate vertices
      const int n_vertices = vertexer.process(tracksAlice3, // track array
                                              idxVec,
                                              gsl::span<o2::InteractionRecord>{bcData},
                                              vertices,
                                              vertexTrackIDs,
                                              v2tRefs,
                                              gsl::span<const o2::MCCompLabel>{lblTracks},
                                              lblVtx);

      if (n_vertices < 1) {
        return; // primary vertex not reconstructed
      }

      // Find largest vertex
      int largestVertex = 0;
      for (Int_t iv = 1; iv < n_vertices; iv++) {
        if (vertices[iv].getNContributors() > vertices[largestVertex].getNContributors()) {
          largestVertex = iv;
        }
      }
      primaryVertex = vertices[largestVertex];
      if (doExtraQA) {
        histos.fill(HIST("h2dVerticesVsContributors"), primaryVertex.getNContributors(), n_vertices);
      }
    } else {
      primaryVertex.setXYZ(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
    }
    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*

    // debug / informational
    histos.fill(HIST("hSimMultiplicity"), multiplicityCounter);
    histos.fill(HIST("hRecoMultiplicity"), tracksAlice3.size());
    histos.fill(HIST("hPVz"), primaryVertex.getZ());

    if (doExtraQA) {
      histos.fill(HIST("hRecoVsSimMultiplicity"), multiplicityCounter, tracksAlice3.size());
    }

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // populate collisions
    collisions(-1, // BC is irrelevant in synthetic MC tests for now, could be adjusted in future
               primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
               primaryVertex.getSigmaX2(), primaryVertex.getSigmaXY(), primaryVertex.getSigmaY2(),
               primaryVertex.getSigmaXZ(), primaryVertex.getSigmaYZ(), primaryVertex.getSigmaZ2(),
               0, primaryVertex.getChi2(), primaryVertex.getNContributors(),
               0, 0);
    collLabels(mcCollision.globalIndex(), 0);
    collisionsAlice3(dNdEta);
    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // populate tracks
    for (const auto& trackParCov : tracksAlice3) {
      // Fixme: collision index could be changeable
      aod::track::TrackTypeEnum trackType = aod::track::Track;

      if (populateTracksDCA) {
        float dcaXY = 1e+10, dcaZ = 1e+10;
        o2::track::TrackParCov trackParametrization(trackParCov);
        if (trackParametrization.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
          dcaXY = dcaInfo.getY();
          dcaZ = dcaInfo.getZ();
        }
        if (doExtraQA && (!extraQAwithoutDecayDaughters || (extraQAwithoutDecayDaughters && !trackParCov.isDecayDau))) {
          histos.fill(HIST("h2dDCAxy"), trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
          histos.fill(HIST("hTrackXatDCA"), trackParametrization.getX());
        }
        tracksDCA(dcaXY, dcaZ);
        if (populateTracksDCACov) {
          tracksDCACov(dcaInfo.getSigmaY2(), dcaInfo.getSigmaZ2());
        }
      }
      tracksPar(collisions.lastIndex(), trackType, trackParCov.getX(), trackParCov.getAlpha(), trackParCov.getY(), trackParCov.getZ(), trackParCov.getSnp(), trackParCov.getTgl(), trackParCov.getQ2Pt());
      tracksParExtension(trackParCov.getPt(), trackParCov.getP(), trackParCov.getEta(), trackParCov.getPhi());

      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tracksParCov(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                   std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tracksParCovExtension(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                            trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                            trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                            trackParCov.getSigma1Pt2());
      tracksLabels(trackParCov.mcLabel, 0);

      // populate extra tables if required to do so
      if (populateTracksExtra) {
        tracksExtra(0.0f, (uint32_t)0, (uint8_t)0, (uint8_t)0,
                    (int8_t)0, (int8_t)0, (uint8_t)0, (uint8_t)0,
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
      }
      if (populateTrackSelection) {
        trackSelection((uint8_t)0, false, false, false, false, false, false);
        trackSelectionExtension(false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);
      }
      TracksAlice3(true);
    }
    // populate ghost tracks
    for (const auto& trackParCov : ghostTracksAlice3) {
      // Fixme: collision index could be changeable
      aod::track::TrackTypeEnum trackType = aod::track::Track;

      if (populateTracksDCA) {
        float dcaXY = 1e+10, dcaZ = 1e+10;
        o2::track::TrackParCov trackParametrization(trackParCov);
        if (trackParametrization.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
          dcaXY = dcaInfo.getY();
          dcaZ = dcaInfo.getZ();
        }
        if (doExtraQA && (!extraQAwithoutDecayDaughters || (extraQAwithoutDecayDaughters && !trackParCov.isDecayDau))) {
          histos.fill(HIST("h2dDCAxy"), trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
          histos.fill(HIST("hTrackXatDCA"), trackParametrization.getX());
        }
        tracksDCA(dcaXY, dcaZ);
        if (populateTracksDCACov) {
          tracksDCACov(dcaInfo.getSigmaY2(), dcaInfo.getSigmaZ2());
        }
      }

      tracksPar(collisions.lastIndex(), trackType, trackParCov.getX(), trackParCov.getAlpha(), trackParCov.getY(), trackParCov.getZ(), trackParCov.getSnp(), trackParCov.getTgl(), trackParCov.getQ2Pt());
      tracksParExtension(trackParCov.getPt(), trackParCov.getP(), trackParCov.getEta(), trackParCov.getPhi());

      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tracksParCov(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                   std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tracksParCovExtension(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                            trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                            trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                            trackParCov.getSigma1Pt2());
      tracksLabels(trackParCov.mcLabel, 0);

      // populate extra tables if required to do so
      if (populateTracksExtra) {
        tracksExtra(0.0f, (uint32_t)0, (uint8_t)0, (uint8_t)0,
                    (int8_t)0, (int8_t)0, (uint8_t)0, (uint8_t)0,
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
      }
      if (populateTrackSelection) {
        trackSelection((uint8_t)0, false, false, false, false, false, false);
        trackSelectionExtension(false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);
      }
      TracksAlice3(false);
    }
  }
};

/// Extends TracksExtra if necessary
struct onTheFlyTrackerInitializer {
  Spawns<aod::TracksExtra> tracksExtra;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<OnTheFlyTracker>(cfgc),
    adaptAnalysisTask<onTheFlyTrackerInitializer>(cfgc)};
}
