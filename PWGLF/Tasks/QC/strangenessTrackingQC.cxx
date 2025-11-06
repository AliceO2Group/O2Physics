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

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
// #include "PWGHF/Core/PDG.h"
#include "PWGLF/DataModel/LFNonPromptCascadeTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{

};

struct miniCasc {
  bool fillOmega;
  float pt;
  float eta;
  float phi;
  float radius;
  float massOmega;
  float massXi;
  float dcaXYCasc;
  float dcaXYTracked;
};

struct strangenessTrackingQC {

  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels>;

  Configurable<int> setting_materialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};

  Configurable<float> cascsetting_dcaCascDaughters{"casc_setting_dcaV0daughters", 0.1f, "DCA between the V0 daughters"};
  Configurable<float> cascsetting_cosPA{"casc_setting_cosPA", 0.995f, "Cosine of the pointing angle of the V0"};
  Configurable<float> cascsetting_massWindowXi{"casc_setting_massWindowXi", 0.01f, "Mass window for the Xi"};

  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> cfgGRPpath{"cfgGRPpath", "GLO/GRP/GRP", "Path of the grp file"};

  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};

  ConfigurableAxis ptBins{"ptBins", {200, -10.f, 10.f}, "Binning for #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis dcaBins{"dcaBins", {1e3, -0.1, 0.1}, "Binning for DCA (cm)"};
  ConfigurableAxis decayRadBins{"decayRadBins", {100, 0.f, 40.f}, "Binning for decay radius (cm)"};
  ConfigurableAxis omegaMassBins{"omegaMassBins", {125, 1.650, 1.700}, "Invariant mass (GeV/#it{c}^{2})"};
  ConfigurableAxis xiMassBins{"xiMassBins", {125, 1.296, 1.346}, "Invariant mass (GeV/#it{c}^{2})"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float bz = 0.f;
  o2::vertexing::DCAFitterN<2> m_fitter;

  HistogramRegistry registry{
    "registry",
    {
      {"omegaMass", "; Mass (GeV/#it{c}^{2}); Counts", {HistType::kTH1F, {{125, 1.650, 1.700}}}},
      {"xiMass", "; Mass (GeV/#it{c}^{2}); Counts", {HistType::kTH1F, {{125, 1.296, 1.346}}}},
      {"omegaMassTracked", "; Mass (GeV/#it{c}^{2}); Counts", {HistType::kTH1F, {{125, 1.650, 1.700}}}},
      {"xiMassTracked", "; Mass (GeV/#it{c}^{2}); Counts", {HistType::kTH1F, {{125, 1.296, 1.346}}}},

      {"omegaHist", "; #it{p}_{T} (GeV/#it{c}); Radius (cm); Mass", {HistType::kTH3F, {{ptBins, decayRadBins, omegaMassBins}}}},
      {"xiHist", "; #it{p}_{T} (GeV/#it{c}); Radius (cm); Mass", {HistType::kTH3F, {{ptBins, decayRadBins, xiMassBins}}}},
      {"xiHistTracked", "; #it{p}_{T} (GeV/#it{c}); Radius (cm); Mass", {HistType::kTH3F, {{ptBins, decayRadBins, xiMassBins}}}},
      {"omegaHistTracked", "; #it{p}_{T} (GeV/#it{c}); Radius (cm); Mass", {HistType::kTH3F, {{ptBins, decayRadBins, omegaMassBins}}}},

      {"xiDCAvsPt", "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", {HistType::kTH2F, {{ptBins, dcaBins}}}},
      {"omegaDCAvsPt", "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", {HistType::kTH2F, {{ptBins, dcaBins}}}},
      {"xiDCAvsPtTracked", "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", {HistType::kTH2F, {{ptBins, dcaBins}}}},
      {"omegaDCAvsPtTracked", "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", {HistType::kTH2F, {{ptBins, dcaBins}}}},
    }};

  template <typename T>
  float dcaToPV(o2::dataformats::VertexBase& PV, T& trackParCov)
  {
    auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(setting_materialCorrection.value);
    o2::dataformats::DCA impactParameterTrk;
    o2::base::Propagator::Instance()->propagateToDCA(PV, trackParCov, bz, 2.f, matCorr, &impactParameterTrk);
    return impactParameterTrk.getY();
  }

  template <typename T>
  bool qualityTrackSelection(const T& track)
  {
    if (std::abs(track.eta()) > 0.9) {
      return false;
    }
    if (track.tpcNClsFound() < cfgCutNclusTPC) {
      return false;
    }
    return true;
  }

  float computeMassMother(const float massA, const float massB, const std::array<float, 3>& momA, const std::array<float, 3>& momB) const
  {
    float eA = std::hypot(massA, std::hypot(momA[0], momA[1], momA[2]));
    float eB = std::hypot(massB, std::hypot(momB[0], momB[1], momB[2]));
    float momTot = std::hypot(momA[0] + momB[0], momA[1] + momB[1], momA[2] + momB[2]);
    float eMother = eA + eB;
    return std::sqrt(eMother * eMother - momTot * momTot);
  }

  template <class TCasc>
  bool buildCascade(TCasc const& casc, CollisionCandidates::iterator const& collision, aod::V0s const&, TrackCandidates const&, miniCasc& miniCasc)
  {
    const auto& v0 = casc.template v0_as<aod::V0s>();
    const auto& bachelor = casc.template bachelor_as<TrackCandidates>();
    const auto& ptrack = v0.template posTrack_as<TrackCandidates>();
    const auto& ntrack = v0.template negTrack_as<TrackCandidates>();
    if (!qualityTrackSelection(ptrack) || !qualityTrackSelection(ntrack) || !qualityTrackSelection(bachelor)) {
      return false;
    }
    const auto& protonTrack = bachelor.sign() > 0 ? ntrack : ptrack;
    const auto& pionTrack = bachelor.sign() > 0 ? ptrack : ntrack;
    if (std::abs(protonTrack.tpcNSigmaPr()) > 3 || std::abs(pionTrack.tpcNSigmaPi()) > 3) {
      return false;
    }
    auto primaryVertex = getPrimaryVertex(collision);
    std::array<float, 3> pvPos = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};

    float cascCpa = -1;
    float cascDauDCA = -1;

    std::array<float, 3> cascMom;
    std::array<float, 3> v0Mom;
    std::array<float, 3> bachelorMom;

    // track propagation
    o2::track::TrackParCov trackParCovV0;
    o2::track::TrackPar trackParV0;
    o2::track::TrackPar trackParBachelor;
    o2::track::TrackParCov trackParCovCasc;
    if (m_fitter.process(getTrackParCov(pionTrack), getTrackParCov(protonTrack))) {
      trackParCovV0 = m_fitter.createParentTrackParCov(0); // V0 track retrieved from p and pi daughters
      if (m_fitter.process(trackParCovV0, getTrackParCov(bachelor))) {
        trackParV0 = m_fitter.getTrackParamAtPCA(0);
        trackParBachelor = m_fitter.getTrackParamAtPCA(1);
        trackParV0.getPxPyPzGlo(v0Mom);
        trackParBachelor.getPxPyPzGlo(bachelorMom);
        trackParCovCasc = m_fitter.createParentTrackParCov();
        trackParCovCasc.getPxPyPzGlo(cascMom);
        cascCpa = RecoDecay::cpa(pvPos, m_fitter.getPCACandidate(), cascMom);
        cascDauDCA = std::sqrt(std::abs(m_fitter.getChi2AtPCACandidate()));
      } else {
        return false;
      }
    } else {
      return false;
    }

    if (cascCpa < cascsetting_cosPA) {
      return false;
    }

    if (cascDauDCA > cascsetting_dcaCascDaughters) {
      return false;
    }

    int chargeFactor = bachelor.sign() > 0 ? 1 : -1;
    miniCasc.pt = chargeFactor * std::hypot(cascMom[0], cascMom[1]);
    miniCasc.massOmega = computeMassMother(constants::physics::MassLambda0, constants::physics::MassKaonCharged, v0Mom, bachelorMom);
    miniCasc.massXi = computeMassMother(constants::physics::MassLambda0, constants::physics::MassPionCharged, v0Mom, bachelorMom);
    miniCasc.fillOmega = false;
    if (TMath::Abs(miniCasc.massXi - constants::physics::MassXiMinus) > cascsetting_massWindowXi && std::abs(bachelor.tpcNSigmaKa()) < 3) {
      miniCasc.fillOmega = true;
    }

    miniCasc.dcaXYCasc = dcaToPV(primaryVertex, trackParCovCasc);
    auto svPos = m_fitter.getPCACandidate();
    miniCasc.radius = std::hypot(svPos[0], svPos[1]);

    return true;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    auto timestamp = bc.timestamp();

    if (o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(cfgGRPpath, timestamp)) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = grpo->getNominalL3Field();
    } else if (o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgGRPmagPath, timestamp)) {
      o2::base::Propagator::initFieldFromGRP(grpmag);
      bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(debug) << "bz = " << bz;
    } else {
      LOG(fatal) << "Got nullptr from CCDB for path " << cfgGRPmagPath << " of object GRPMagField and " << cfgGRPpath << " of object GRPObject for timestamp " << timestamp;
    }
  }

  void init(InitContext const&)
  {
    mRunNumber = 0;
    bz = 0;

    if (static_cast<o2::base::Propagator::MatCorrType>(setting_materialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
      LOG(info) << "Setting material correction LUT";
      o2::base::Propagator::Instance(true)->setMatLUT(lut);
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    m_fitter.setPropagateToPCA(true);
    m_fitter.setMaxR(200.);
    m_fitter.setMinParamChange(1e-3);
    m_fitter.setMinRelChi2Change(0.9);
    m_fitter.setMaxDZIni(4);
    m_fitter.setMaxDXYIni(4);
    m_fitter.setMaxChi2(1e9);
    m_fitter.setUseAbsDCA(true);
    m_fitter.setWeightedFinalPCA(false);
    // int mat{static_cast<int>(setting_materialCorrection)};
    // m_fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));
  }

  void process(CollisionCandidates const&, aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades, aod::V0s const& v0s, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {

    for (const auto& trackedCascade : trackedCascades) {
      miniCasc miniCasc;
      const auto& casc = trackedCascade.cascade();
      auto collision = trackedCascade.collision_as<CollisionCandidates>();
      if (!collision.sel8() || std::abs(collision.posZ()) > 10) {
        continue;
      }
      initCCDB(collision.bc_as<aod::BCsWithTimestamps>());
      m_fitter.setBz(bz);
      if (buildCascade(casc, collision, v0s, tracks, miniCasc)) {

        // compute the dca of the tracked cascade
        const auto& track = trackedCascade.track_as<TrackCandidates>();
        auto trackCovTrk = getTrackParCov(track);

        auto primaryVertex = getPrimaryVertex(collision);
        miniCasc.dcaXYTracked = dcaToPV(primaryVertex, trackCovTrk);
        // fill the histograms
        if (miniCasc.fillOmega) {
          registry.fill(HIST("omegaMassTracked"), miniCasc.massOmega);
          registry.fill(HIST("omegaDCAvsPtTracked"), miniCasc.pt, miniCasc.dcaXYTracked);
          registry.fill(HIST("omegaHistTracked"), miniCasc.pt, miniCasc.radius, miniCasc.massOmega);
        }
        registry.fill(HIST("xiMassTracked"), miniCasc.massXi);
        registry.fill(HIST("xiDCAvsPtTracked"), miniCasc.pt, miniCasc.dcaXYTracked);
        registry.fill(HIST("xiHistTracked"), miniCasc.pt, miniCasc.radius, miniCasc.massXi);
      }
    }

    for (auto& cascade : cascades) {
      miniCasc miniCasc;
      auto collision = cascade.collision_as<CollisionCandidates>();
      if (!collision.sel8() || std::abs(collision.posZ()) > 10) {
        continue;
      }
      initCCDB(collision.bc_as<aod::BCsWithTimestamps>());
      m_fitter.setBz(bz);
      if (buildCascade(cascade, collision, v0s, tracks, miniCasc)) {
        if (miniCasc.fillOmega) {
          registry.fill(HIST("omegaMass"), miniCasc.massOmega);
          registry.fill(HIST("omegaDCAvsPt"), miniCasc.pt, miniCasc.dcaXYCasc);
          registry.fill(HIST("omegaHist"), miniCasc.pt, miniCasc.radius, miniCasc.massOmega);
        }
        registry.fill(HIST("xiMass"), miniCasc.massXi);
        registry.fill(HIST("xiDCAvsPt"), miniCasc.pt, miniCasc.dcaXYCasc);
        registry.fill(HIST("xiHist"), miniCasc.pt, miniCasc.radius, miniCasc.massXi);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessTrackingQC>(cfgc)};
}
