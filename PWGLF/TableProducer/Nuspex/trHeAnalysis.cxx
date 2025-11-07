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
/// \file trHeAnalysis.cxx
///
/// \brief Triton and Helion Analysis on pp Data
///
/// \author Matthias Herzer <matthias.herzer@cern.ch>, Goethe University Frankfurt
///
#include "PWGLF/DataModel/LFNucleiTables.h"
#include "PWGLF/DataModel/LFPIDTOFGenericTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/Utils/pidTOFGeneric.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TF1.h>

#include <string>
#include <vector>

namespace o2::aod
{
namespace h3_data
{
DECLARE_SOA_COLUMN(TPt, tPt, float);
DECLARE_SOA_COLUMN(TEta, tEta, float);
DECLARE_SOA_COLUMN(TPhi, tPhi, float);
DECLARE_SOA_COLUMN(TCharge, tCharge, int8_t);
DECLARE_SOA_COLUMN(TP, tP, float);
DECLARE_SOA_COLUMN(TH3DeDx, tH3DeDx, float);
DECLARE_SOA_COLUMN(TnSigmaTpc, tnSigmaTpc, float);
DECLARE_SOA_COLUMN(TTofSignalH3, tTofSignalH3, float);
DECLARE_SOA_COLUMN(TDcaXY, tDcaXY, float);
DECLARE_SOA_COLUMN(TDcaZ, tDcaZ, float);
DECLARE_SOA_COLUMN(TSigmaYX, tSigmaYX, float);
DECLARE_SOA_COLUMN(TSigmaXYZ, tSigmaXYZ, float);
DECLARE_SOA_COLUMN(TSigmaZ, tSigmaZ, float);
DECLARE_SOA_COLUMN(TnTpcCluster, tnTpcCluster, int);
DECLARE_SOA_COLUMN(TnItsCluster, tnItsCluster, int);
DECLARE_SOA_COLUMN(TTpcChi2NCl, tTpcChi2NCl, float);
DECLARE_SOA_COLUMN(TItsChi2NCl, tItsChi2NCl, float);
DECLARE_SOA_COLUMN(TRigidity, tRigidity, float);
DECLARE_SOA_COLUMN(TItsClusterSize, tItsClusterSize, float);
DECLARE_SOA_COLUMN(THasTof, tHasTof, bool);
DECLARE_SOA_COLUMN(TDetectorMap, tDetectorMap, int8_t);
} // namespace h3_data
DECLARE_SOA_TABLE(H3Data, "AOD", "h3_data", h3_data::TPt, h3_data::TEta,
                  h3_data::TPhi, h3_data::TCharge, h3_data::TH3DeDx,
                  h3_data::TnSigmaTpc, h3_data::TTofSignalH3, h3_data::TDcaXY,
                  h3_data::TDcaZ, h3_data::TSigmaYX, h3_data::TSigmaXYZ,
                  h3_data::TSigmaZ, h3_data::TnTpcCluster,
                  h3_data::TnItsCluster, h3_data::TTpcChi2NCl,
                  h3_data::TItsChi2NCl, h3_data::TRigidity,
                  h3_data::TItsClusterSize, h3_data::THasTof, h3_data::TDetectorMap);
namespace he_data
{
DECLARE_SOA_COLUMN(TPt, tPt, float);
DECLARE_SOA_COLUMN(TEta, tEta, float);
DECLARE_SOA_COLUMN(TPhi, tPhi, float);
DECLARE_SOA_COLUMN(TCharge, tCharge, int8_t);
DECLARE_SOA_COLUMN(TP, tP, float);
DECLARE_SOA_COLUMN(THeDeDx, tHeDeDx, float);
DECLARE_SOA_COLUMN(TnSigmaTpc, tnSigmaTpc, float);
DECLARE_SOA_COLUMN(TTofSignalHe, tTofSignalHe, float);
DECLARE_SOA_COLUMN(TDcaXY, tDcaXY, float);
DECLARE_SOA_COLUMN(TDcaZ, tDcaZ, float);
DECLARE_SOA_COLUMN(TSigmaYX, tSigmaYX, float);
DECLARE_SOA_COLUMN(TSigmaXYZ, tSigmaXYZ, float);
DECLARE_SOA_COLUMN(TSigmaZ, tSigmaZ, float);
DECLARE_SOA_COLUMN(TnTpcCluster, tnTpcCluster, int);
DECLARE_SOA_COLUMN(TnItsCluster, tnItsCluster, int);
DECLARE_SOA_COLUMN(TTpcChi2NCl, tTpcChi2NCl, float);
DECLARE_SOA_COLUMN(TItsChi2NCl, tItsChi2NCl, float);
DECLARE_SOA_COLUMN(TRigidity, tRigidity, float);
DECLARE_SOA_COLUMN(TItsClusterSize, tItsClusterSize, float);
DECLARE_SOA_COLUMN(THasTof, tHasTof, bool);
DECLARE_SOA_COLUMN(TDetectorMap, tDetectorMap, int8_t);
} // namespace he_data
DECLARE_SOA_TABLE(HeData, "AOD", "he_data", he_data::TPt, he_data::TEta,
                  he_data::TPhi, he_data::TCharge, he_data::THeDeDx,
                  he_data::TnSigmaTpc, he_data::TTofSignalHe, he_data::TDcaXY,
                  he_data::TDcaZ, he_data::TSigmaYX, he_data::TSigmaXYZ,
                  he_data::TSigmaZ, he_data::TnTpcCluster,
                  he_data::TnItsCluster, he_data::TTpcChi2NCl,
                  he_data::TItsChi2NCl, he_data::TRigidity,
                  he_data::TItsClusterSize, he_data::THasTof, he_data::TDetectorMap);
} // namespace o2::aod
namespace
{
const int nBetheParams = 6;
const int nParticles = 2;
static const std::vector<std::string> particleNames{"triton", "helion"};
static const std::vector<int> particlePdgCodes{
  o2::constants::physics::kTriton, o2::constants::physics::kHelium3};
static const std::vector<float> particleMasses{
  o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3};
static const std::vector<int> particleCharge{1, 2};
static const std::vector<float> particleChargeFactor{2.3, 2.55};
static const std::vector<std::string> betheBlochParNames{
  "p0", "p1", "p2", "p3", "p4", "resolution"};
constexpr float BetheBlochDefault[nParticles][nBetheParams]{
  {0.248753, 3.58634, 0.0167065, 2.29194, 0.774344,
   0.07}, // triton
  {0.0274556, 18.3054, 3.99987e-05, 3.17219, 11.1775,
   0.07}}; // Helion
} // namespace
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using TracksFull =
  soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
            o2::aod::TracksDCA, aod::pidTOFmass, aod::pidTOFbeta,
            aod::pidTPCLfFullTr, aod::pidTPCLfFullHe,
            aod::TOFSignal, aod::TrackSelectionExtension,
            o2::aod::EvTimeTOFFT0ForTrack>;

class Particle
{
 public:
  TString name;
  int pdgCode;
  float mass;
  int charge;
  float resolution;
  float chargeFactor;
  std::vector<float> betheParams;
  static constexpr int NNumBetheParams = 5;

  Particle(const std::string name_, int pdgCode_, float mass_, int charge_,
           LabeledArray<float> bethe, float chargeFactor_)
  {
    name = TString(name_);
    pdgCode = pdgCode_;
    mass = mass_;
    charge = charge_;
    chargeFactor = chargeFactor_;

    resolution =
      bethe.get(name, "resolution"); // Access the "resolution" parameter

    betheParams.clear();
    for (int i = 0; i < NNumBetheParams; ++i) {
      betheParams.push_back(bethe.get(name, i));
    }
  }
};

struct TrHeAnalysis {
  Produces<o2::aod::H3Data> h3Data;
  Produces<o2::aod::HeData> heData;
  HistogramRegistry histos{
    "Histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  std::vector<Particle> particles;
  Configurable<bool> enableTr{"enableTr", true, "Flag to enable triton analysis."};
  Configurable<bool> enableHe{"enableHe", true, "Flag to enable helium-3 analysis."};
  Configurable<bool> cfgRigidityCorrection{"cfgRigidityCorrection", true, "Enable Rigidity correction"};
  ConfigurableAxis binsDeDx{"binsDeDx", {600, 0.f, 3000.f}, ""};
  ConfigurableAxis binsBeta{"binsBeta", {120, 0.0, 1.2}, ""};
  ConfigurableAxis binsDca{"binsDca", {400, -1.f, 1.f}, ""};
  ConfigurableAxis binsSigmaTpc{"binsSigmaTpc", {1000, -100, 100}, ""};
  ConfigurableAxis binsSigmaTof{"binsSigmaTof", {1000, -100, 100}, ""};
  ConfigurableAxis binsMassTr{"binsMassTr", {250, -2.5, 2.5f}, ""};
  ConfigurableAxis binsMassHe{"binsMassHe", {300, -3., 3.f}, ""};
  // Set the event selection cuts
  struct : ConfigurableGroup {
    Configurable<bool> useSel8{"useSel8", true, "Use Sel8 for run3 Event Selection"};
    Configurable<bool> tvxTrigger{"tvxTrigger", false, "Use TVX for Event Selection (default w/ Sel8)"};
    Configurable<bool> removeTfBorder{"removeTfBorder", false, "Remove TimeFrame border (default w/ Sel8)"};
    Configurable<bool> removeItsRofBorder{"removeItsRofBorder", false, "Remove ITS Read-Out Frame border (default w/ Sel8)"};
  } evselOptions;

  Configurable<bool> cfgTPCPidMethod{"cfgTPCPidMethod", false, "Using own or built in bethe parametrization"}; // false for built in
  Configurable<int> cfgMassMethod{"cfgMassMethod", 0, "0: Using built in 1: mass calculated with beta 2: mass calculated with the event time"};
  Configurable<bool> cfgEnableItsClusterSizeCut{"cfgEnableItsClusterSizeCut", false, "Enable ITS cluster size cut"};
  Configurable<bool> cfgEnableTofMassCut{"cfgEnableTofMassCut", false, "Enable TOF mass cut"};
  Configurable<float> cfgTofMassCutPt{"cfgTofMassCutPt", 1.6f, "Pt value for which the TOF-cut starts to be used"};
  // Set the multiplity event limits
  Configurable<float> cfgLowMultCut{"cfgLowMultCut", 0.0f, "Accepted multiplicity percentage lower limit"};
  Configurable<float> cfgHighMultCut{"cfgHighMultCut", 100.0f, "Accepted multiplicity percentage higher limit"};

  // Set the z-vertex event cut limits
  Configurable<float> cfgHighCutVertex{"cfgHighCutVertex", 10.0f, "Accepted z-vertex upper limit"};
  Configurable<float> cfgLowCutVertex{"cfgLowCutVertex", -10.0f, "Accepted z-vertex lower limit"};

  // Set the quality cuts for tracks
  Configurable<bool> rejectFakeTracks{"rejectFakeTracks", false, "Flag to reject ITS-TPC fake tracks (for MC)"};
  Configurable<float> cfgCutItsClusters{"cfgCutItsClusters", -1.f, "Minimum number of ITS clusters"};
  Configurable<float> cfgCutTpcXRows{"cfgCutTpcXRows", -1.f, "Minimum number of crossed TPC rows"};
  Configurable<float> cfgCutTpcClusters{"cfgCutTpcClusters", 40.f, "Minimum number of found TPC clusters"};
  Configurable<int> nItsLayer{"nItsLayer", 0, "ITS Layer (0-6)"};
  Configurable<float> cfgCutTpcCrRowToFindableCl{"cfgCutTpcCrRowToFindableCl", 0.8f, "Minimum ratio of crossed rows to findable cluster in TPC"};
  Configurable<float> cfgCutMaxChi2TpcH3{"cfgCutMaxChi2TpcH3", 4.f, "Maximum chi2 per cluster for TPC"};
  Configurable<float> cfgCutMaxChi2ItsH3{"cfgCutMaxChi2ItsH3", 36.f, "Maximum chi2 per cluster for ITS"};
  Configurable<float> cfgCutMaxChi2TpcHe{"cfgCutMaxChi2TpcHe", 4.f, "Maximum chi2 per cluster for TPC"};
  Configurable<float> cfgCutMaxChi2ItsHe{"cfgCutMaxChi2ItsHe", 36.f, "Maximum chi2 per cluster for ITS"};
  Configurable<bool> cfgCutTpcRefit{"cfgCutTpcRefit", 1, "TPC refit "};
  Configurable<bool> cfgCutItsRefit{"cfgCutItsRefit", 1, "ITS refit"};
  Configurable<float> cfgCutMaxItsClusterSizeHe{"cfgCutMaxItsClusterSizeHe", 4.f, "Maximum ITS Cluster Size for He "};
  Configurable<float> cfgCutMinItsClusterSizeHe{"cfgCutMinItsClusterSizeHe", 1.f, "Minimum ITS Cluster Size for He"};
  Configurable<float> cfgCutMaxItsClusterSizeH3{"cfgCutMaxItsClusterSizeH3", 4.f, "Maximum ITS Cluster Size for Tr"};
  Configurable<float> cfgCutMinItsClusterSizeH3{"cfgCutMinItsClusterSizeH3", 1.f, "Minimum ITS Cluster Size for Tr"};
  Configurable<float> cfgCutMinTofMassH3{"cfgCutMinTofMassH3", 5.f, "Minimum Tof mass H3"};
  Configurable<float> cfgCutMaxTofMassH3{"cfgCutMaxTofMassH3", 11.f, "Maximum TOF mass H3"};
  Configurable<float> cfgMaxRigidity{"cfgMaxRigidity", 10.f, "Maximum rigidity value"};
  Configurable<float> cfgMaxPt{"cfgMaxPt", 10.f, "Maximum pT value"};

  // Set the kinematic and PID cuts for tracks
  struct : ConfigurableGroup {
    Configurable<float> pCut{"pCut", 0.6f, "Value of the p selection for spectra (default 0.3)"};
    Configurable<float> etaCut{"etaCut", 0.8f, "Value of the eta selection for spectra (default 0.8)"};
    Configurable<float> yLowCut{"yLowCut", -1.0f, "Value of the low rapidity selection for spectra (default -1.0)"};
    Configurable<float> yHighCut{"yHighCut", 1.0f, "Value of the high rapidity selection for spectra (default 1.0)"};
  } kinemOptions;

  struct : ConfigurableGroup {
    Configurable<float> nsigmaTPCTr{"nsigmaTPCTr", 5.f, "Value of the Nsigma TPC cut for tritons"};
    Configurable<float> nsigmaTPCHe{"nsigmaTPCHe", 5.f, "Value of the Nsigma TPC cut for helium-3"};
  } nsigmaTPCvar;
  Configurable<LabeledArray<float>> cfgBetheBlochParams{"cfgBetheBlochParams", {BetheBlochDefault[0], nParticles, nBetheParams, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec dedxAxis{binsDeDx, "d#it{E}/d#it{x} A.U."};
    const AxisSpec betaAxis{binsBeta, "TOF #beta"};
    const AxisSpec dcaxyAxis{binsDca, "DCAxy (cm)"};
    const AxisSpec dcazAxis{binsDca, "DCAz (cm)"};
    const AxisSpec massTrAxis{binsMassTr, ""};
    const AxisSpec massHeAxis{binsMassHe, ""};
    const AxisSpec sigmaTPCAxis{binsSigmaTpc, ""};
    const AxisSpec sigmaTOFAxis{binsSigmaTof, ""};

    histos.add<TH1>("histogram/pT",
                    "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts",
                    HistType::kTH1F, {{500, 0., 10.}});
    histos.add<TH1>("histogram/p", "Track momentum; p (GeV/#it{c}); counts",
                    HistType::kTH1F, {{500, 0., 10.}});
    histos.add<TH2>("histogram/TPCsignVsTPCmomentum",
                    "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC "
                    "<-dE/dx> (a.u.)",
                    HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
    histos.add<TH2>(
      "histogram/TOFbetaVsP",
      "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta",
      HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    histos.add<TH2>("histogram/H3/H3-TPCsignVsTPCmomentum",
                    "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC "
                    "<-dE/dx> (a.u.)",
                    HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
    histos.add<TH2>(
      "histogram/H3/H3-TOFbetaVsP",
      "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta",
      HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    histos.add<TH2>("histogram/He/He-TPCsignVsTPCmomentum",
                    "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC "
                    "<-dE/dx> (a.u.)",
                    HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
    histos.add<TH2>(
      "histogram/He/He-TOFbetaVsP",
      "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta",
      HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    histos.add<TH1>("event/eventSelection", "eventSelection", HistType::kTH1D,
                    {{7, -0.5, 6.5}});
    auto h = histos.get<TH1>(HIST("event/eventSelection"));
    h->GetXaxis()->SetBinLabel(1, "Total");
    h->GetXaxis()->SetBinLabel(2, "TVX trigger cut");
    h->GetXaxis()->SetBinLabel(3, "TF border cut");
    h->GetXaxis()->SetBinLabel(4, "ITS ROF cut");
    h->GetXaxis()->SetBinLabel(5, "TVX + TF + ITS ROF");
    h->GetXaxis()->SetBinLabel(6, "Sel8 cut");
    h->GetXaxis()->SetBinLabel(7, "Z-vert Cut");
    histos.add<TH1>("histogram/cuts", "cuts", HistType::kTH1D,
                    {{13, -0.5, 12.5}});
    auto hCuts = histos.get<TH1>(HIST("histogram/cuts"));
    hCuts->GetXaxis()->SetBinLabel(1, "total");
    hCuts->GetXaxis()->SetBinLabel(2, "p cut");
    hCuts->GetXaxis()->SetBinLabel(3, "eta cut");
    hCuts->GetXaxis()->SetBinLabel(4, "TPC cluster");
    hCuts->GetXaxis()->SetBinLabel(5, "ITS clsuter");
    hCuts->GetXaxis()->SetBinLabel(6, "TPC crossed rows");
    hCuts->GetXaxis()->SetBinLabel(7, "max chi2 ITS");
    hCuts->GetXaxis()->SetBinLabel(8, "max chi2 TPC");
    hCuts->GetXaxis()->SetBinLabel(9, "crossed rows over findable cluster");
    hCuts->GetXaxis()->SetBinLabel(10, "TPC refit");
    hCuts->GetXaxis()->SetBinLabel(11, "ITS refit");
    hCuts->GetXaxis()->SetBinLabel(12, "ITS cluster size");
    hCuts->GetXaxis()->SetBinLabel(13, "TOF mass cut");
    for (int i = 0; i < nParticles; i++) {
      particles.push_back(Particle(particleNames.at(i), particlePdgCodes.at(i),
                                   particleMasses.at(i), particleCharge.at(i),
                                   cfgBetheBlochParams, particleChargeFactor.at(i)));
    }
  }
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& event,
               TracksFull const& tracks)
  {
    bool trRapCut = kFALSE;
    bool heRapCut = kFALSE;
    histos.fill(HIST("event/eventSelection"), 0);
    if ((event.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
        (event.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
        (event.selection_bit(aod::evsel::kIsTriggerTVX))) {
      histos.fill(HIST("event/eventSelection"), 4);
    }
    if (evselOptions.useSel8 && !event.sel8())
      return;
    histos.fill(HIST("event/eventSelection"), 5);
    if (event.posZ() < cfgLowCutVertex || event.posZ() > cfgHighCutVertex)
      return;
    histos.fill(HIST("event/eventSelection"), 6);
    if (cfgTPCPidMethod) {
      for (const auto& track : tracks) {
        trRapCut =
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Triton)) >
            kinemOptions.yLowCut &&
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Triton)) <
            kinemOptions.yHighCut;
        heRapCut =
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)) >
            kinemOptions.yLowCut &&
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)) <
            kinemOptions.yHighCut;
        histos.fill(HIST("histogram/cuts"), 0);
        if (std::abs(track.p()) < kinemOptions.pCut) {
          histos.fill(HIST("histogram/cuts"), 1);
          continue;
        }
        if (std::abs(track.eta()) > kinemOptions.etaCut) {
          histos.fill(HIST("histogram/cuts"), 2);
          continue;
        }
        if (track.pt() > cfgMaxPt || getRigidity(track) > cfgMaxRigidity) {
          continue;
        }
        if (track.tpcNClsFound() < cfgCutTpcClusters) {
          histos.fill(HIST("histogram/cuts"), 3);
          continue;
        }
        if (track.itsNCls() < cfgCutItsClusters) {
          histos.fill(HIST("histogram/cuts"), 4);
          continue;
        }
        if (track.tpcNClsCrossedRows() < cfgCutTpcXRows) {
          histos.fill(HIST("histogram/cuts"), 5);
          continue;
        }
        if (track.tpcCrossedRowsOverFindableCls() <= cfgCutTpcCrRowToFindableCl) {
          histos.fill(HIST("histogram/cuts"), 8);
          continue;
        }
        if (cfgCutTpcRefit) {
          if (!track.passedTPCRefit()) {
            histos.fill(HIST("histogram/cuts"), 9);
            continue;
          }
        }
        if (cfgCutItsRefit) {
          if (!track.passedITSRefit()) {
            histos.fill(HIST("histogram/cuts"), 10);
            continue;
          }
        }
        histos.fill(HIST("histogram/pT"), track.pt());
        histos.fill(HIST("histogram/p"), track.p());
        histos.fill(HIST("histogram/TPCsignVsTPCmomentum"),
                    getRigidity(track) * track.sign(),
                    track.tpcSignal());
        histos.fill(HIST("histogram/TOFbetaVsP"),
                    getRigidity(track) * track.sign(),
                    track.beta());
        if (enableTr && trRapCut) {
          if (std::abs(getTPCnSigma(track, particles.at(0))) <
              nsigmaTPCvar.nsigmaTPCTr) {
            if (track.itsChi2NCl() > cfgCutMaxChi2ItsH3) {
              histos.fill(HIST("histogram/cuts"), 6);
              continue;
            }
            if (track.tpcChi2NCl() > cfgCutMaxChi2TpcH3) {
              histos.fill(HIST("histogram/cuts"), 7);
              continue;
            }
            if (cfgEnableItsClusterSizeCut) {
              if (getMeanItsClsSize(track) / std::cosh(track.eta()) <= cfgCutMinItsClusterSizeHe ||
                  getMeanItsClsSize(track) / std::cosh(track.eta()) >= cfgCutMaxItsClusterSizeHe) {
                histos.fill(HIST("histogram/cuts"), 12);
                continue;
              }
            }
            if (cfgEnableTofMassCut && track.pt() > cfgTofMassCutPt) {
              if (getMass(track) < cfgCutMinTofMassH3 || getMass(track) > cfgCutMaxTofMassH3) {
                histos.fill(HIST("histogram/cuts"), 13);
                continue;
              }
            }
            histos.fill(HIST("histogram/H3/H3-TPCsignVsTPCmomentum"),
                        getRigidity(track) * track.sign(),
                        track.tpcSignal());
            histos.fill(HIST("histogram/H3/H3-TOFbetaVsP"),
                        getRigidity(track) * track.sign(),
                        track.beta());
            float tPt = track.pt();
            float tEta = track.eta();
            float tPhi = track.phi();
            int8_t tCharge = track.sign();
            float tH3DeDx = track.tpcSignal();
            float tnSigmaTpc = getTPCnSigma(track, particles.at(0));
            float tTofSignalH3 = getMass(track);
            float tDcaXY = track.dcaXY();
            float tDcaZ = track.dcaZ();
            float tSigmaYX = track.sigmaY();
            float tSigmaXYZ = track.sigmaSnp();
            float tSigmaZ = track.sigmaZ();
            int tnTpcCluster = track.tpcNClsFound();
            int tnItsCluster = track.itsNCls();
            float tTpcChi2NCl = track.tpcChi2NCl();
            float tItsChi2NCl = track.itsChi2NCl();
            float tRigidity = getRigidity(track);
            float tItsClusterSize =
              getMeanItsClsSize(track) / std::cosh(track.eta());
            bool tHasTof = track.hasTOF();
            int8_t tDetectorMap = track.detectorMap();
            h3Data(tPt, tEta, tPhi, tCharge, tH3DeDx, tnSigmaTpc, tTofSignalH3,
                   tDcaXY, tDcaZ, tSigmaYX, tSigmaXYZ, tSigmaZ, tnTpcCluster,
                   tnItsCluster, tTpcChi2NCl, tItsChi2NCl, tRigidity,
                   tItsClusterSize, tHasTof, tDetectorMap);
          }
        }
        if (enableHe && heRapCut) {
          if (std::abs(getTPCnSigma(track, particles.at(1))) <
              nsigmaTPCvar.nsigmaTPCHe) {
            if (track.itsChi2NCl() > cfgCutMaxChi2ItsHe) {
              histos.fill(HIST("histogram/cuts"), 6);
              continue;
            }
            if (track.tpcChi2NCl() > cfgCutMaxChi2TpcHe) {
              histos.fill(HIST("histogram/cuts"), 7);
              continue;
            }
            if (cfgEnableItsClusterSizeCut) {
              if (getMeanItsClsSize(track) / std::cosh(track.eta()) <= cfgCutMinItsClusterSizeHe ||
                  getMeanItsClsSize(track) / std::cosh(track.eta()) >= cfgCutMaxItsClusterSizeHe) {
                histos.fill(HIST("histogram/cuts"), 12);
                continue;
              }
            }
            histos.fill(HIST("histogram/He/He-TPCsignVsTPCmomentum"),
                        getRigidity(track) * track.sign(),
                        track.tpcSignal());
            histos.fill(HIST("histogram/He/He-TOFbetaVsP"),
                        getRigidity(track) * track.sign(),
                        track.beta());
            float tPt = track.pt();
            float tEta = track.eta();
            float tPhi = track.phi();
            int8_t tCharge = 2.f * track.sign();
            float tHeDeDx = track.tpcSignal();
            float tnSigmaTpc = getTPCnSigma(track, particles.at(1));
            float tTofSignalHe = getMass(track);
            float tDcaXY = track.dcaXY();
            float tDcaZ = track.dcaZ();
            float tSigmaYX = track.sigmaY();
            float tSigmaXYZ = track.sigmaSnp();
            float tSigmaZ = track.sigmaZ();
            int tnTpcCluster = track.tpcNClsFound();
            int tnItsCluster = track.itsNCls();
            float tTpcChi2NCl = track.tpcChi2NCl();
            float tItsChi2NCl = track.itsChi2NCl();
            float tRigidity = getRigidity(track);
            float tItsClusterSize =
              getMeanItsClsSize(track) / std::cosh(track.eta());
            bool tHasTof = track.hasTOF();
            int8_t tDetectorMap = track.detectorMap();
            heData(tPt, tEta, tPhi, tCharge, tHeDeDx, tnSigmaTpc, tTofSignalHe,
                   tDcaXY, tDcaZ, tSigmaYX, tSigmaXYZ, tSigmaZ, tnTpcCluster,
                   tnItsCluster, tTpcChi2NCl, tItsChi2NCl, tRigidity,
                   tItsClusterSize, tHasTof, tDetectorMap);
          }
        }
      }
    }
    if (!cfgTPCPidMethod) {
      for (const auto& track : tracks) {
        trRapCut =
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Triton)) >
            kinemOptions.yLowCut &&
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Triton)) <
            kinemOptions.yHighCut;
        heRapCut =
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)) >
            kinemOptions.yLowCut &&
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)) <
            kinemOptions.yHighCut;
        histos.fill(HIST("histogram/cuts"), 0);
        if (std::abs(track.p()) < kinemOptions.pCut) {
          histos.fill(HIST("histogram/cuts"), 1);
          continue;
        }
        if (std::abs(track.eta()) > kinemOptions.etaCut) {
          histos.fill(HIST("histogram/cuts"), 2);
          continue;
        }
        if (track.pt() > cfgMaxPt || getRigidity(track) > cfgMaxRigidity) {
          continue;
        }
        if (track.tpcNClsFound() < cfgCutTpcClusters) {
          histos.fill(HIST("histogram/cuts"), 3);
          continue;
        }
        if (track.itsNCls() < cfgCutItsClusters) {
          histos.fill(HIST("histogram/cuts"), 4);
          continue;
        }
        if (track.tpcNClsCrossedRows() < cfgCutTpcXRows) {
          histos.fill(HIST("histogram/cuts"), 5);
          continue;
        }
        if (track.tpcCrossedRowsOverFindableCls() <= cfgCutTpcCrRowToFindableCl) {
          histos.fill(HIST("histogram/cuts"), 8);
          continue;
        }
        if (cfgCutTpcRefit) {
          if (!track.passedTPCRefit()) {
            histos.fill(HIST("histogram/cuts"), 9);
            continue;
          }
        }
        if (cfgCutItsRefit) {
          if (!track.passedITSRefit()) {
            histos.fill(HIST("histogram/cuts"), 10);
            continue;
          }
        }
        histos.fill(HIST("histogram/pT"), track.pt());
        histos.fill(HIST("histogram/p"), track.p());
        histos.fill(HIST("histogram/TPCsignVsTPCmomentum"),
                    getRigidity(track) * (1.f * track.sign()),
                    track.tpcSignal());
        histos.fill(HIST("histogram/TOFbetaVsP"),
                    track.p() * (1.f * track.sign()), track.beta());
        if (enableTr && trRapCut) {
          if (std::abs(track.tpcNSigmaTr()) < nsigmaTPCvar.nsigmaTPCTr) {
            if (track.itsChi2NCl() > cfgCutMaxChi2ItsH3) {
              histos.fill(HIST("histogram/cuts"), 6);
              continue;
            }
            if (track.tpcChi2NCl() > cfgCutMaxChi2TpcH3) {
              histos.fill(HIST("histogram/cuts"), 7);
              continue;
            }
            if (cfgEnableItsClusterSizeCut) {
              if (getMeanItsClsSize(track) / std::cosh(track.eta()) <= cfgCutMinItsClusterSizeH3 ||
                  getMeanItsClsSize(track) / std::cosh(track.eta()) >= cfgCutMaxItsClusterSizeH3) {
                histos.fill(HIST("histogram/cuts"), 12);
                continue;
              }
            }
            if (cfgEnableTofMassCut && track.pt() > cfgTofMassCutPt) {
              if (getMass(track) < cfgCutMinTofMassH3 || getMass(track) > cfgCutMaxTofMassH3) {
                histos.fill(HIST("histogram/cuts"), 13);
                continue;
              }
            }
            histos.fill(HIST("histogram/H3/H3-TPCsignVsTPCmomentum"),
                        getRigidity(track) * (1.f * track.sign()),
                        track.tpcSignal());
            histos.fill(HIST("histogram/H3/H3-TOFbetaVsP"),
                        track.p() * (1.f * track.sign()),
                        track.beta());
            float tPt = track.pt();
            float tEta = track.eta();
            float tPhi = track.phi();
            int8_t tCharge = track.sign();
            float tH3DeDx = track.tpcSignal();
            float tnSigmaTpc = track.tpcNSigmaTr();
            float tTofSignalH3 = getMass(track);
            float tDcaXY = track.dcaXY();
            float tDcaZ = track.dcaZ();
            float tSigmaYX = track.sigmaY();
            float tSigmaXYZ = track.sigmaSnp();
            float tSigmaZ = track.sigmaZ();
            int tnTpcCluster = track.tpcNClsFound();
            int tnItsCluster = track.itsNCls();
            float tTpcChi2NCl = track.tpcChi2NCl();
            float tItsChi2NCl = track.itsChi2NCl();
            float tRigidity = getRigidity(track);
            float tItsClusterSize =
              getMeanItsClsSize(track) / std::cosh(track.eta());
            bool tHasTof = track.hasTOF();
            int8_t tDetectorMap = track.detectorMap();
            h3Data(tPt, tEta, tPhi, tCharge, tH3DeDx, tnSigmaTpc, tTofSignalH3,
                   tDcaXY, tDcaZ, tSigmaYX, tSigmaXYZ, tSigmaZ, tnTpcCluster,
                   tnItsCluster, tTpcChi2NCl, tItsChi2NCl, tRigidity,
                   tItsClusterSize, tHasTof, tDetectorMap);
          }
        }
        if (enableHe && heRapCut) {
          if (std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe) {
            if (track.itsChi2NCl() > cfgCutMaxChi2ItsHe) {
              histos.fill(HIST("histogram/cuts"), 6);
              continue;
            }
            if (track.tpcChi2NCl() > cfgCutMaxChi2TpcHe) {
              histos.fill(HIST("histogram/cuts"), 7);
              continue;
            }
            if (cfgEnableItsClusterSizeCut) {
              if (getMeanItsClsSize(track) / std::cosh(track.eta()) <= cfgCutMinItsClusterSizeHe ||
                  getMeanItsClsSize(track) / std::cosh(track.eta()) >= cfgCutMaxItsClusterSizeHe) {
                histos.fill(HIST("histogram/cuts"), 12);
                continue;
              }
            }
            histos.fill(HIST("histogram/He/He-TPCsignVsTPCmomentum"),
                        getRigidity(track) * track.sign(),
                        track.tpcSignal());
            histos.fill(HIST("histogram/He/He-TOFbetaVsP"),
                        getRigidity(track) * track.sign(),
                        track.beta());
            float tPt = track.pt();
            float tEta = track.eta();
            float tPhi = track.phi();
            int8_t tCharge = 2.f * track.sign();
            float tHeDeDx = track.tpcSignal();
            float tnSigmaTpc = track.tpcNSigmaHe();
            float tTofSignalHe = getMass(track);
            float tDcaXY = track.dcaXY();
            float tDcaZ = track.dcaZ();
            float tSigmaYX = track.sigmaY();
            float tSigmaXYZ = track.sigmaSnp();
            float tSigmaZ = track.sigmaZ();
            int tnTpcCluster = track.tpcNClsFound();
            int tnItsCluster = track.itsNCls();
            float tTpcChi2NCl = track.tpcChi2NCl();
            float tItsChi2NCl = track.itsChi2NCl();
            float tRigidity = getRigidity(track);
            float tItsClusterSize =
              getMeanItsClsSize(track) / std::cosh(track.eta());
            bool tHasTof = track.hasTOF();
            int8_t tDetectorMap = track.detectorMap();
            heData(tPt, tEta, tPhi, tCharge, tHeDeDx, tnSigmaTpc, tTofSignalHe,
                   tDcaXY, tDcaZ, tSigmaYX, tSigmaXYZ, tSigmaZ, tnTpcCluster,
                   tnItsCluster, tTpcChi2NCl, tItsChi2NCl, tRigidity,
                   tItsClusterSize, tHasTof, tDetectorMap);
          }
        }
      }
    }
  }

  template <class T>
  float getTPCnSigma(T const& track, Particle const& particle)
  {
    const float rigidity = getRigidity(track);
    if (!track.hasTPC())
      return -999;

    float expBethe{betheBlochAleph(particle, rigidity)};
    float expSigma{expBethe * particle.resolution};
    float sigmaTPC =
      static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
    return sigmaTPC;
  }

  template <class T>
  float betheBlochAleph(Particle const& particle, T const& rigidity)
  {
    double bg = particle.charge * rigidity / particle.mass;
    double beta = bg / std::sqrt(1. + bg * bg);
    double aa = std::pow(beta, particle.betheParams[3]);
    double bb = std::pow(1. / bg, particle.betheParams[4]);
    if ((particle.betheParams[2] + bb) <= 0)
      return 0;
    bb = std::log(particle.betheParams[2] + bb);
    return std::pow(particle.charge, particle.chargeFactor) * 50 * (particle.betheParams[1] - aa - bb) * particle.betheParams[0] / aa;
  }

  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    constexpr int NNumLayers = 8;
    constexpr int NBitsPerLayer = 4;
    constexpr int NBitMask = (1 << NBitsPerLayer) - 1;

    int sum = 0, n = 0;
    for (int i = 0; i < NNumLayers; i++) {
      int clsSize = (track.itsClusterSizes() >> (NBitsPerLayer * i)) & NBitMask;
      sum += clsSize;
      if (clsSize) {
        n++;
      }
    }
    return n > 0 ? static_cast<float>(sum) / n : 0.f;
  }
  template <class T>
  float getRigidity(T const& track)
  {
    if (!cfgRigidityCorrection)
      return track.tpcInnerParam();
    bool hePID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    return hePID ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
  }
  template <class T>
  float getMass(const T& track)
  {
    if (cfgMassMethod == 0) {
      float m = track.mass();
      return m * m;
    }
    if (cfgMassMethod == 1) {
      const float beta = track.beta();
      const float rigidity = getRigidity(track);
      float gamma = 1.f / std::sqrt(1.f - beta * beta);
      float mass = rigidity / std::sqrt(gamma * gamma - 1.f);
      return mass * mass;
    }
    if (cfgMassMethod == 2) {
      const float rigidity = getRigidity(track);
      float tofStartTime = track.evTimeForTrack();
      float tofTime = track.tofSignal();
      constexpr float CInCmPs = 2.99792458e-2f;
      float length = track.length();
      float time = tofTime - tofStartTime;
      if (time > 0.f && length > 0.f) {
        float beta = length / (CInCmPs * time);
        float gamma = 1.f / std::sqrt(1.f - beta * beta);
        float mass = rigidity / std::sqrt(gamma * gamma - 1.f);
        return mass * mass;
      }
      return -1.f;
    }
    return -1.f;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TrHeAnalysis>(cfgc),
  };
}
