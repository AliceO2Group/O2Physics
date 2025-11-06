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
//
// Nuclei spectra analysis task
// ========================
//
// Executable + dependencies:
//
// Data (run3):
// o2-analysis-lf-nuclei-spectra, o2-analysis-timestamp
// o2-analysis-pid-tof-base, o2-analysis-multiplicity-table, o2-analysis-event-selection
// (to add flow: o2-analysis-qvector-table, o2-analysis-centrality-table)

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFSlimNucleiTables.h"
#include "PWGLF/Utils/nucleiUtils.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Tools/TrackTuner.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector4D.h"
#include "TRandom3.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct nucleiFlowTree {
  enum {
    kProton = BIT(0),
    kDeuteron = BIT(1),
    kTriton = BIT(2),
    kHe3 = BIT(3),
    kHe4 = BIT(4),
    kHasTOF = BIT(5),
    kHasTRD = BIT(6),
    kIsAmbiguous = BIT(7), /// just a placeholder now
    kITSrof = BIT(8),
    kIsPhysicalPrimary = BIT(9), /// MC flags starting from the second half of the short
    kIsSecondaryFromMaterial = BIT(10),
    kIsSecondaryFromWeakDecay = BIT(11) /// the last 4 bits are reserved for the PID in tracking
  };

  Produces<o2::aod::NucleiTable> nucleiTable;
  Produces<o2::aod::NucleiTableMC> nucleiTableMC;
  Produces<o2::aod::NucleiTableFlow> nucleiTableFlow;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};

  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Centrality estimator (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3)"};
  Configurable<float> cfgCMrapidity{"cfgCMrapidity", 0.f, "Rapidity of the center of mass (only for p-Pb)"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutTpcMom{"cfgCutTpcMom", 0.2f, "Minimum TPC momentum for tracks"};
  Configurable<float> cfgCutRapidityMin{"cfgCutRapidityMin", -0.5, "Minimum rapidity for tracks"};
  Configurable<float> cfgCutRapidityMax{"cfgCutRapidityMax", 0.5, "Maximum rapidity for tracks"};
  Configurable<bool> cfgCutOnReconstructedRapidity{"cfgCutOnReconstructedRapidity", false, "Cut on reconstructed rapidity"};
  Configurable<float> cfgCutNclusITS{"cfgCutNclusITS", 5, "Minimum number of ITS clusters"};
  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};
  Configurable<float> cfgCutPtMinTree{"cfgCutPtMinTree", 0.2f, "Minimum track transverse momentum for tree saving"};
  Configurable<float> cfgCutPtMaxTree{"cfgCutPtMaxTree", 15.0f, "Maximum track transverse momentum for tree saving"};

  Configurable<LabeledArray<int>> cfgEventSelections{"cfgEventSelections", {nuclei::EvSelDefault[0], 8, 1, nuclei::eventSelectionLabels, nuclei::eventSelectionTitle}, "Event selections"};

  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {nuclei::bbMomScalingDefault[0], 5, 2, nuclei::names, nuclei::chargeLabelNames}, "TPC Bethe-Bloch momentum scaling for light nuclei"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {nuclei::betheBlochDefault[0], 5, 6, nuclei::names, nuclei::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgNsigmaTPC{"cfgNsigmaTPC", {nuclei::nSigmaTPCdefault[0], 5, 2, nuclei::names, nuclei::nSigmaConfigName}, "TPC nsigma selection for light nuclei"};
  Configurable<LabeledArray<double>> cfgDCAcut{"cfgDCAcut", {nuclei::DCAcutDefault[0], 5, 2, nuclei::names, nuclei::nDCAConfigName}, "Max DCAxy and DCAz for light nuclei"};
  Configurable<LabeledArray<double>> cfgDownscaling{"cfgDownscaling", {nuclei::DownscalingDefault[0], 5, 1, nuclei::names, nuclei::DownscalingConfigName}, "Fraction of kept candidates for light nuclei"};
  Configurable<LabeledArray<int>> cfgTreeConfig{"cfgTreeConfig", {nuclei::TreeConfigDefault[0], 5, 2, nuclei::names, nuclei::treeConfigNames}, "Filtered trees configuration"};

  ConfigurableAxis cfgNITSClusBins{"cfgNITSClusBins", {3, 4.5, 7.5}, "N ITS clusters binning"};
  ConfigurableAxis cfgNTPCClusBins{"cfgNTPCClusBins", {3, 89.5, 159.5}, "N TPC clusters binning"};

  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};

  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // CCDB options
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgCCDBurl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfgZorroCCDBpath{"cfgZorroCCDBpath", "/Users/m/mpuccio/EventFiltering/OTS/", "path to the zorro ccdb objects"};
  int mRunNumber = 0;
  float mBz = 0.f;

  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime>;

  // Configurable Harmonics index
  Configurable<int> cfgHarmonics{"cfgHarmonics", 2, "Harmonics index for flow analysis"};

  // Collisions with chentrality
  using CollWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentNTPVs, aod::CentFT0Cs>::iterator;

  // Flow analysis
  using CollWithEP = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::EPCalibrationTables>::iterator;

  using CollWithQvec = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorFV0AVecs, aod::QvectorTPCallVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs>::iterator;

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  float computeEventPlane(float y, float x)
  {
    return 0.5 * std::atan2(y, x);
  }

  template <class Tcoll>
  bool eventSelectionWithHisto(Tcoll& collision)
  {
    spectra.fill(HIST("hEventSelections"), 0);

    if (cfgEventSelections->get(nuclei::evSel::kTVX) && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kTVX + 1);

    if (cfgEventSelections->get(nuclei::evSel::kZvtx) && std::abs(collision.posZ()) > cfgCutVertex) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kZvtx + 1);

    if (cfgEventSelections->get(nuclei::evSel::kTFborder) && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kTFborder + 1);

    if (cfgEventSelections->get(nuclei::evSel::kITSROFborder) && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kITSROFborder + 1);

    if (cfgEventSelections->get(nuclei::evSel::kNoSameBunchPileup) && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kNoSameBunchPileup + 1);

    if (cfgEventSelections->get(nuclei::evSel::kIsGoodZvtxFT0vsPV) && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kIsGoodZvtxFT0vsPV + 1);

    if (cfgEventSelections->get(nuclei::evSel::kIsGoodITSLayersAll) && !collision.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kIsGoodITSLayersAll + 1);

    if constexpr (
      requires {
        collision.triggereventep();
      }) {
      if (cfgEventSelections->get(nuclei::evSel::kIsEPtriggered) && !collision.triggereventep()) {
        return false;
      }
      spectra.fill(HIST("hEventSelections"), nuclei::evSel::kIsEPtriggered + 1);
    }

    float centrality = getCentrality(collision);
    spectra.fill(HIST("hCentrality"), centrality);

    return true;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fHe");
      zorro.populateHistRegistry(spectra, bc.runNumber());
    }
    auto timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(nuclei::lut);
    mBz = static_cast<float>(grpmag->getNominalL3Field());
    LOGF(info, "Retrieved GRP for timestamp %ull (%i) with magnetic field of %1.2f kZG", timestamp, mRunNumber, mBz);
  }

  void init(o2::framework::InitContext&)
  {
    zorroSummary.setObject(zorro.getZorroSummary());
    zorro.setBaseCCDBPath(cfgZorroCCDBpath.value);
    ccdb->setURL(cfgCCDBurl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    spectra.add("hEventSelections", "hEventSelections", {HistType::kTH1D, {{nuclei::evSel::kNevSels + 1, -0.5f, static_cast<float>(nuclei::evSel::kNevSels) + 0.5f}}});
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(1, "all");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kTVX + 2, "TVX");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kZvtx + 2, "Zvtx");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kTFborder + 2, "TFborder");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kITSROFborder + 2, "ITSROFborder");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kNoSameBunchPileup + 2, "kNoSameBunchPileup");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kIsGoodZvtxFT0vsPV + 2, "isGoodZvtxFT0vsPV");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kIsGoodITSLayersAll + 2, "IsGoodITSLayersAll");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kIsEPtriggered + 2, "IsEPtriggered");

    spectra.add("hCentrality", "hCentrality", HistType::kTH1D, {{100, 0., 100., "Centrality (%)"}});

    spectra.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., 20., "z position (cm)"}});
    spectra.add("hTpcSignalData", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("hTpcSignalDataSelected", "Specific energy loss for selected particles", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("hTofSignalData", "TOF beta", HistType::kTH2F, {{500, 0., 5., "#it{p} (GeV/#it{c})"}, {750, 0, 1.5, "TOF #beta"}});

    for (int iS{0}; iS < nuclei::Species::kNspecies; ++iS) {
      for (int iMax{0}; iMax < 2; ++iMax) {
        nuclei::pidCutTPC[iS][iMax] = cfgNsigmaTPC->get(iS, iMax); // changed pidCut to pidCutTPC so that it compiles TODO: check if it is correct
      }
    }

    nuclei::lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
  }

  template <typename Tcoll>
  float getCentrality(Tcoll const& collision)
  {
    float centrality = 1.;
    if constexpr (o2::aod::HasCentrality<Tcoll>) {
      if (cfgCentralityEstimator == nuclei::centDetectors::kFV0A) {
        centrality = collision.centFV0A();
      } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0M) {
        centrality = collision.centFT0M();
      } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0A) {
        centrality = collision.centFT0A();
      } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0C) {
        centrality = collision.centFT0C();
      } else {
        LOG(warning) << "Centrality estimator not valid. Possible values: (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3). Centrality set to 1.";
      }
    }
    return centrality;
  }

  template <typename Tcoll, typename Ttrks>
  void fillDataInfo(Tcoll const& collision, Ttrks const& tracks)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    if (cfgSkimmedProcessing) {
      zorro.isSelected(bc.globalBC()); /// Just let Zorro do the accounting
    }
    gRandom->SetSeed(bc.timestamp());

    spectra.fill(HIST("hRecVtxZData"), collision.posZ());

    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

    const double bgScalings[5][2]{
      {nuclei::charges[0] * cfgMomentumScalingBetheBloch->get(0u, 0u) / nuclei::masses[0], nuclei::charges[0] * cfgMomentumScalingBetheBloch->get(0u, 1u) / nuclei::masses[0]},
      {nuclei::charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / nuclei::masses[1], nuclei::charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / nuclei::masses[1]},
      {nuclei::charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / nuclei::masses[2], nuclei::charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / nuclei::masses[2]},
      {nuclei::charges[3] * cfgMomentumScalingBetheBloch->get(3u, 0u) / nuclei::masses[3], nuclei::charges[3] * cfgMomentumScalingBetheBloch->get(3u, 1u) / nuclei::masses[3]},
      {nuclei::charges[4] * cfgMomentumScalingBetheBloch->get(3u, 0u) / nuclei::masses[4], nuclei::charges[4] * cfgMomentumScalingBetheBloch->get(3u, 1u) / nuclei::masses[4]}};

    for (auto& track : tracks) { // start loop over tracks
      if (std::abs(track.eta()) > cfgCutEta ||
          track.tpcInnerParam() < cfgCutTpcMom ||
          track.itsNCls() < cfgCutNclusITS ||
          track.tpcNClsFound() < cfgCutNclusTPC ||
          track.tpcNClsCrossedRows() < 70 ||
          track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
          track.tpcChi2NCl() > 4.f ||
          track.itsChi2NCl() > 36.f) {
        continue;
      }
      // temporary fix: tpcInnerParam() returns the momentum in all the software tags before
      bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
      float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();

      spectra.fill(HIST("hTpcSignalData"), correctedTpcInnerParam * track.sign(), track.tpcSignal());
      const int iC{track.sign() < 0};

      bool selectedTPC[5]{false}, goodToAnalyse{false};
      std::array<float, 5> nSigmaTPC;

      for (int iS{0}; iS < nuclei::Species::kNspecies; ++iS) {

        double expBethe{tpc::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScalings[iS][iC]), cfgBetheBlochParams->get(iS, 0u), cfgBetheBlochParams->get(iS, 1u), cfgBetheBlochParams->get(iS, 2u), cfgBetheBlochParams->get(iS, 3u), cfgBetheBlochParams->get(iS, 4u))};

        double expSigma{expBethe * cfgBetheBlochParams->get(iS, 5u)};

        nSigmaTPC[iS] = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);

        selectedTPC[iS] = (nSigmaTPC[iS] > nuclei::pidCutTPC[iS][0] && nSigmaTPC[iS] < nuclei::pidCutTPC[iS][1]);

        goodToAnalyse = goodToAnalyse || selectedTPC[iS];
      }
      if (!goodToAnalyse) {
        continue;
      }

      setTrackParCov(track, mTrackParCov);
      mTrackParCov.setPID(track.pidForTracking());

      std::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, mTrackParCov, mBz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

      float beta{o2::pid::tof::Beta::GetBeta(track)};
      spectra.fill(HIST("hTpcSignalDataSelected"), correctedTpcInnerParam * track.sign(), track.tpcSignal());
      spectra.fill(HIST("hTofSignalData"), correctedTpcInnerParam, beta);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      uint16_t flag = static_cast<uint16_t>((track.pidForTracking() & 0xF) << 12);
      std::array<float, 5> tofMasses{-3.f, -3.f, -3.f, -3.f, -3.f};
      bool fillTree{true}; // set to true and never used again
      bool fillDCAHist{false};
      bool correctPV{false};
      bool isSecondary{false};
      bool fromWeakDecay{false};

      if (track.hasTOF()) {
        flag |= kHasTOF;
      }
      if (track.hasTRD()) {
        flag |= kHasTRD;
      }
      if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        flag |= kITSrof;
      }
      for (int iS{0}; iS < nuclei::Species::kNspecies; ++iS) {
        bool selectedTOF{false};
        if (std::abs(dcaInfo[1]) > cfgDCAcut->get(iS, 1)) {
          continue;
        }
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> fvector{mTrackParCov.getPt() * nuclei::charges[iS], mTrackParCov.getEta(), mTrackParCov.getPhi(), nuclei::masses[iS]};
        if (selectedTPC[iS]) {
          if (track.hasTOF()) {
            selectedTOF = true; /// temporarly skipped
            float charge{1.f + static_cast<float>(iS == 3 || iS == 4)};
            tofMasses[iS] = correctedTpcInnerParam * charge * std::sqrt(1.f / (beta * beta) - 1.f) - nuclei::masses[iS];
          }
          if (cfgTreeConfig->get(iS, 1u) && !selectedTOF) {
            continue;
          }
          bool setPartFlag = cfgTreeConfig->get(iS, 0u);
          if (setPartFlag) {
            if (cfgDownscaling->get(iS) < 1. && gRandom->Rndm() > cfgDownscaling->get(iS)) {
              continue;
            }
            flag |= BIT(iS);
          }
        }
      }
      if (flag & (kProton | kDeuteron | kTriton | kHe3 | kHe4) /*|| doprocessMC*/) { /// ignore PID pre-selections for the MC
        if constexpr (requires {
                        collision.psiFT0A();
                      }) {
          nuclei::candidates_flow.emplace_back(NucleusCandidateFlow{
            collision.centFV0A(),
            collision.centFT0M(),
            collision.centFT0A(),
            collision.centFT0C(),
            collision.psiFT0A(),
            collision.psiFT0C(),
            collision.psiTPC(),
            collision.psiTPCL(),
            collision.psiTPCR(),
            collision.qFT0A(),
            collision.qFT0C(),
            collision.qTPC(),
            collision.qTPCL(),
            collision.qTPCR(),
          });
        } else if constexpr (requires {
                               collision.qvecFT0AImVec()[cfgHarmonics - 2];
                             }) {
          nuclei::candidates_flow.emplace_back(NucleusCandidateFlow{
            collision.centFV0A(),
            collision.centFT0M(),
            collision.centFT0A(),
            collision.centFT0C(),
            computeEventPlane(collision.qvecFT0AImVec()[cfgHarmonics - 2], collision.qvecFT0AReVec()[cfgHarmonics - 2]),
            computeEventPlane(collision.qvecFT0CImVec()[cfgHarmonics - 2], collision.qvecFT0CReVec()[cfgHarmonics - 2]),
            computeEventPlane(collision.qvecTPCallImVec()[cfgHarmonics - 2], collision.qvecTPCallReVec()[cfgHarmonics - 2]),
            computeEventPlane(collision.qvecTPCnegImVec()[cfgHarmonics - 2], collision.qvecTPCnegReVec()[cfgHarmonics - 2]),
            computeEventPlane(collision.qvecTPCposImVec()[cfgHarmonics - 2], collision.qvecTPCposReVec()[cfgHarmonics - 2]),
            std::hypot(collision.qvecFT0AImVec()[cfgHarmonics - 2], collision.qvecFT0AReVec()[cfgHarmonics - 2]),
            std::hypot(collision.qvecFT0CImVec()[cfgHarmonics - 2], collision.qvecFT0CReVec()[cfgHarmonics - 2]),
            std::hypot(collision.qvecTPCallImVec()[cfgHarmonics - 2], collision.qvecTPCallReVec()[cfgHarmonics - 2]),
            std::hypot(collision.qvecTPCnegImVec()[cfgHarmonics - 2], collision.qvecTPCnegReVec()[cfgHarmonics - 2]),
            std::hypot(collision.qvecTPCposImVec()[cfgHarmonics - 2], collision.qvecTPCposReVec()[cfgHarmonics - 2])});
        }
        if (flag & kTriton) {
          if (track.pt() < cfgCutPtMinTree || track.pt() > cfgCutPtMaxTree || track.sign() > 0)
            continue;
        }
        nuclei::candidates.emplace_back(NucleusCandidate{
          static_cast<int>(track.globalIndex()), static_cast<int>(track.collisionId()), (1 - 2 * iC) * mTrackParCov.getPt(), mTrackParCov.getEta(), mTrackParCov.getPhi(),
          correctedTpcInnerParam, beta, collision.posZ(), collision.numContrib(), dcaInfo[0], dcaInfo[1], track.tpcSignal(), track.itsChi2NCl(), track.tpcChi2NCl(), track.tofChi2(),
          nSigmaTPC, tofMasses, fillTree, fillDCAHist, correctPV, isSecondary, fromWeakDecay, flag, track.tpcNClsFindable(), static_cast<uint8_t>(track.tpcNClsCrossedRows()), track.itsClusterMap(),
          static_cast<uint8_t>(track.tpcNClsFound()), static_cast<uint8_t>(track.tpcNClsShared()), static_cast<uint8_t>(track.itsNCls()), static_cast<uint32_t>(track.itsClusterSizes())});
      }
    } // end loop over tracks
  }

  void processDataFlow(CollWithEP const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
    nuclei::candidates_flow.clear();
    if (!eventSelectionWithHisto(collision)) {
      return;
    }
    fillDataInfo(collision, tracks);
    for (auto& c : nuclei::candidates) {
      nucleiTable(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.nContrib, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.TOFchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.TPCnClsShared, c.clusterSizesITS);
    }
    for (auto& c : nuclei::candidates_flow) {
      nucleiTableFlow(c.centFV0A, c.centFT0M, c.centFT0A, c.centFT0C, c.psiFT0A, c.psiFT0C, c.psiTPC, c.psiTPCl, c.psiTPCr, c.qFT0A, c.qFT0C, c.qTPC, c.qTPCl, c.qTPCr);
    }
  }
  PROCESS_SWITCH(nucleiFlowTree, processDataFlow, "Data analysis with flow", true);

  void processDataFlowAlternative(CollWithQvec const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
    nuclei::candidates_flow.clear();
    if (!eventSelectionWithHisto(collision)) {
      return;
    }
    fillDataInfo(collision, tracks);
    for (auto& c : nuclei::candidates) {
      nucleiTable(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.nContrib, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.TOFchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.TPCnClsShared, c.clusterSizesITS);
    }
    for (auto& c : nuclei::candidates_flow) {
      nucleiTableFlow(c.centFV0A, c.centFT0M, c.centFT0A, c.centFT0C, c.psiFT0A, c.psiFT0C, c.psiTPC, c.psiTPCl, c.psiTPCr, c.qFT0A, c.qFT0C, c.qTPC, c.qTPCl, c.qTPCr);
    }
  }
  PROCESS_SWITCH(nucleiFlowTree, processDataFlowAlternative, "Data analysis with flow - alternative framework", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiFlowTree>(cfgc, TaskName{"nuclei-flow-trees"})};
}
