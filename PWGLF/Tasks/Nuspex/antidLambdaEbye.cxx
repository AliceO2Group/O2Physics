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

#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCDe, aod::pidTPCPr, aod::pidTPCPi, aod::TOFSignal, aod::TOFEvTime>;

namespace
{
// constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
constexpr double betheBlochDefault[1][6]{{-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNamesBB{"d"};
static const std::vector<std::string> pidHypotheses{"Electron", "Muon", "Pion", "Kaon", "Proton", "Deuteron", "Triton", "He3", "Alpha", "Pion0", "Photon", "K0", "Lambda", "HyperTriton", "Hyperhydrog4", "XiMinus", "OmegaMinus"};
} // namespace

struct antidLambdaEbye {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::pid::tof::Beta<TracksFull::iterator> responseBeta;

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.03f, o2::constants::physics::MassLambda0 + 0.03f}, "binning for the lambda invariant-mass"};
  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};

  ConfigurableAxis cosPaAxis{"cosPaAxis", {1e3, 0.95f, 1.00f}, "binning for the cosPa axis"};
  ConfigurableAxis radiusAxis{"radiusAxis", {1e3, 0.f, 100.f}, "binning for the radius axis"};
  ConfigurableAxis dcaV0daughAxis{"dcaV0daughAxis", {2e2, 0.f, 2.f}, "binning for the dca of V0 daughters"};
  ConfigurableAxis dcaDaughPvAxis{"dcaDaughPvAxis", {1e3, -10.f, 10.f}, "binning for the dca of positive daughter to PV"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNamesBB, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  ConfigurableAxis tpcNsigmaAxis{"tpcNsigmaAxis", {100, -5.f, 5.f}, "tpc nsigma axis"};
  // ConfigurableAxis tofNsigmaAxis{"tofNsigmaAxis", {100, -5.f, 5.f}, "tof nsigma axis"};
  ConfigurableAxis tofMassAxis{"tofMassAxis", {1000, 0., 3.f}, "tof mass axis"};
  ConfigurableAxis momAxis{"momAxis", {60., 0.f, 3.f}, "momentum axis binning"};
  ConfigurableAxis momAxisFine{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};
  ConfigurableAxis momResAxis{"momResAxis", {1.e2, -1.f, 1.f}, "momentum resolution binning"};
  ConfigurableAxis tpcAxis{"tpcAxis", {4.e2, 0.f, 4.e3f}, "tpc signal axis binning"};
  ConfigurableAxis tofAxis{"tofAxis", {1.e3, 0.f, 1.f}, "tof signal axis binning"};
  ConfigurableAxis trackingPidAxis{"trackingPidAxis", {static_cast<double>(pidHypotheses.size()), 0, static_cast<double>(pidHypotheses.size())}, "tracking pid hypothesis binning"};

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};

  Configurable<float> antidPtMin{"antidPtMin", 0.8f, "minimum antideuteron pT (GeV/c)"};
  Configurable<float> antidPtTof{"antidPtTof", 1.0f, "antideuteron pT to switch to TOF pid (GeV/c) "};
  Configurable<float> antidPtMax{"antidPtMax", 1.8f, "maximum antideuteron pT (GeV/c)"};

  Configurable<float> lambdaPtMin{"lambdaPtMin", 0.5f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 3.0f, "maximum (anti)lambda pT (GeV/c)"};

  Configurable<float> antidNcrossedRows{"antidNcrossedRows", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> antidNclusItsCut{"antidNclusITScut", 5, "Minimum number of ITS clusters"};
  Configurable<float> antidNclusTpcCut{"antidNclusTPCcut", 70, "Minimum number of TPC clusters"};
  Configurable<float> antidNsigmaTpcCut{"antidNsigmaTpcCut", 4.f, "TPC PID cut"};
  Configurable<float> antidNsigmaTofCut{"antidNsigmaTofCut", 4.f, "TOF PID cut"};
  Configurable<float> antidDcaCut{"antidDcaCut", 0.1f, "DCA antid to PV"};
  Configurable<float> tpcInnerParamMax{"tpcInnerParamMax", 0.6f, "(temporary) tpc inner param cut"};
  Configurable<float> tofMassMax{"tofMassMax", 0.3f, "(temporary) tof mass cut"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1f, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1f, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5f, "v0radius"};
  Configurable<float> lambdaMassCut{"lambdaMassCut", 0.005f, "maximum deviation from PDG mass"};

  int mRunNumber;
  float d_bz;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  template <class RecV0>
  bool selectLambda(RecV0 const& v0) // TODO: apply ML
  {
    if (std::abs(v0.eta()) > etaMax ||
        v0.v0cosPA() < v0setting_cospa ||
        v0.v0radius() < v0setting_radius) {
      return false;
    }
    auto mLambda = v0.alpha() > 0 ? v0.mLambda() : v0.mAntiLambda();
    if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCut) {
      return false;
    }
    return true;
  }

  template <class T>
  bool selectAntid(T const& track)
  {
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    if (track.sign() > 0.) {
      return false;
    }
    if (track.itsNCls() < antidNclusItsCut ||
        track.tpcNClsFound() < antidNclusTpcCut ||
        track.tpcNClsCrossedRows() < antidNcrossedRows ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcChi2NCl() > 4.f ||
        track.itsChi2NCl() > 36.f) {
      return false;
    }
    return true;
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    mRunNumber = 0;
    d_bz = 0;

    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    histos.add<TH1>("nEv", ";Centrality (%);Entries", {HistType::kTH1D}, {centAxis});

    histos.add<TH1>("q1antid", ";Centrality (%);#it{q}_{1}(#bar{d})", {HistType::kTH1D}, {centAxis});
    histos.add<TH1>("q1sqantid", ";Centrality (%);#it{q}_{1}^{2}(#bar{d})", {HistType::kTH1D}, {centAxis});
    histos.add<TH1>("q2antid", ";Centrality (%);#it{q}_{2}(#bar{d})", {HistType::kTH1D}, {centAxis});

    histos.add<TH1>("q1antiL", ";Centrality (%);#it{q}_{1}(#bar{#Lambda})", {HistType::kTH1D}, {centAxis});
    histos.add<TH1>("q1sqantiL", ";Centrality (%);#it{q}_{1}^{2}(#bar{#Lambda})", {HistType::kTH1D}, {centAxis});
    histos.add<TH1>("q2antiL", ";Centrality (%);#it{q}_{2}(#bar{#Lambda})", {HistType::kTH1D}, {centAxis});

    histos.add<TH1>("q1L", ";Centrality (%);#it{q}_{1}(#Lambda)", {HistType::kTH1D}, {centAxis});
    histos.add<TH1>("q1sqL", ";Centrality (%);#it{q}_{1}^{2}(#Lambda)", {HistType::kTH1D}, {centAxis});
    histos.add<TH1>("q2L", ";Centrality (%);#it{q}_{2}(#Lambda)", {HistType::kTH1D}, {centAxis});

    histos.add<TH1>("q11LantiL", ";Centrality (%);#it{q}_{11}(#Lambda, #bar{#Lambda})", {HistType::kTH1D}, {centAxis});
    histos.add<TH1>("q11Lantid", ";Centrality (%);#it{q}_{11}(#Lambda, #bar{d})", {HistType::kTH1D}, {centAxis});
    histos.add<TH1>("q11antiLantid", ";Centrality (%);#it{q}_{11}(#bar{#Lambda}, #bar{d})", {HistType::kTH1D}, {centAxis});

    // v0 QA
    histos.add<TH1>("massLambda", ";#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", {HistType::kTH1F, {massLambdaAxis}});
    histos.add<TH1>("cosPa", ";cosPa;Entries", {HistType::kTH1F}, {cosPaAxis});
    histos.add<TH1>("radius", ";radius;Entries", {HistType::kTH1F}, {radiusAxis});
    histos.add<TH1>("dcaV0daugh", ";dcaV0daugh;Entries", {HistType::kTH1F}, {dcaV0daughAxis});
    histos.add<TH1>("dcaPosPv", ";dcaPosPv;Entries", {HistType::kTH1F}, {dcaDaughPvAxis});
    histos.add<TH1>("dcaNegPv", ";dcaNegPv;Entries", {HistType::kTH1F}, {dcaDaughPvAxis});

    // antid QA
    histos.add<TH2>("tpcNsigma", ";#it{p}_{TPC} (GeV/#it{c});tpcNsigma", {HistType::kTH2F}, {momAxis, tpcNsigmaAxis});
    histos.add<TH2>("tpcNsigmaGlo", ";#it{p}_{T} (GeV/#it{c});tpcNsigma", {HistType::kTH2F}, {momAxis, tpcNsigmaAxis});
    // histos.add<TH2>("tofNsigma", ";#it{p}_{glo};tofNsigma;Entries", {HistType::kTH2F}, {momAxis, tofNsigmaAxis});
    histos.add<TH2>("tofMass", ";#it{p}_{glo} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", {HistType::kTH2F}, {momAxis, tofMassAxis});
    histos.add<TH2>("tofMassFull", ";#it{p}_{glo} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", {HistType::kTH2F}, {momAxis, tofMassAxis});
    auto hmomCorr = histos.add<TH3>("momCorr", ";#it{p}_{glo} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c});", {HistType::kTH3F}, {momAxisFine, momResAxis, trackingPidAxis});
    histos.add<TH2>("tpcSignal", ";#it{p}_{TPC} (GeV/#it{c});TPC signal (a.u.)", {HistType::kTH2F}, {momAxisFine, tpcAxis});
    histos.add<TH2>("tpcSignalBkg", ";#it{p}_{TPC} (GeV/#it{c});TPC signal (a.u.)", {HistType::kTH2F}, {momAxisFine, tpcAxis});
    auto htpcSignal_glo = histos.add<TH3>("tpcSignal_glo", ";#it{p}_{glo} (GeV/#it{c});TPC signal (a.u.);", {HistType::kTH3F}, {momAxisFine, tpcAxis, trackingPidAxis});
    auto htpcSignalBkg_glo = histos.add<TH3>("tpcSignalBkg_glo", ";#it{p}_{glo} (GeV/#it{c});TPC signal (a.u.);", {HistType::kTH3F}, {momAxisFine, tpcAxis, trackingPidAxis});
    histos.add<TH2>("tofSignal", ";#it{p}_{TPC} (GeV/#it{c});TOF signal (a.u.)", {HistType::kTH2F}, {momAxisFine, tofAxis});
    histos.add<TH2>("tofSignal_glo", ";#it{p}_{T} (GeV/#it{c});TOF signal (a.u.)", {HistType::kTH2F}, {momAxisFine, tofAxis});

    for (int i{1}; i < hmomCorr->GetNbinsZ() + 1; ++i) {
      hmomCorr->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
      htpcSignal_glo->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
      htpcSignalBkg_glo->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();

    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    }
    mRunNumber = bc.runNumber();
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>::iterator const& collision, TracksFull const& tracks, soa::Filtered<aod::V0Datas> const& V0s, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!collision.sel8())
      return;

    if (std::abs(collision.posZ()) > zVtxMax)
      return;

    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

    histos.fill(HIST("zVtx"), collision.posZ());

    double q1antid{0.}, q2antid{0.};
    for (const auto& track : tracks) {
      if (!selectAntid(track)) {
        continue;
      }

      auto trackParCov = getTrackParCov(track);
      gpu::gpustd::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(collVtx, trackParCov, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

      float dcaXYZ = dcaInfo[0];
      if (std::abs(dcaXYZ) > antidDcaCut) {
        continue;
      }

      histos.fill(HIST("tpcSignal"), track.tpcInnerParam(), track.tpcSignal());
      histos.fill(HIST("tpcSignal_glo"), trackParCov.getP(), track.tpcSignal(), track.pidForTracking());

      if (trackParCov.getPt() < antidPtMin || trackParCov.getPt() > antidPtMax) {
        continue;
      }

      double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() / constants::physics::MassDeuteron), cfgBetheBlochParams->get("d", "p0"), cfgBetheBlochParams->get("d", "p1"), cfgBetheBlochParams->get("d", "p2"), cfgBetheBlochParams->get("d", "p3"), cfgBetheBlochParams->get("d", "p4"))};
      double expSigma{expBethe * cfgBetheBlochParams->get("d", "resolution")};
      auto nSigmaTPC = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);

      float beta{track.hasTOF() ? responseBeta.GetBeta(track) : -999.f};
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta));
      float mass{track.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f)};
      bool hasTof = track.hasTOF() && track.tofChi2() < 3;

      if (std::abs(nSigmaTPC) > antidNsigmaTpcCut) {
        continue;
      }
      histos.fill(HIST("tpcNsigma"), track.tpcInnerParam(), nSigmaTPC);
      histos.fill(HIST("momCorr"), trackParCov.getP(), track.tpcInnerParam() - trackParCov.getP(), track.pidForTracking());
      // check contamination
      if (track.tpcInnerParam() < tpcInnerParamMax) {
        histos.fill(HIST("tpcSignalBkg"), track.tpcInnerParam(), track.tpcSignal());
        histos.fill(HIST("tpcSignalBkg_glo"), trackParCov.getP(), track.tpcSignal(), track.pidForTracking());
      }

      if (trackParCov.getPt() > antidPtTof && hasTof) {
        histos.fill(HIST("tofSignal_glo"), trackParCov.getP(), beta);
        histos.fill(HIST("tofSignal"), track.tpcInnerParam(), beta);
      }

      // temporary cut to reject fake matches
      if (track.tpcInnerParam() < tpcInnerParamMax) {
        continue;
      }
      if (trackParCov.getPt() > antidPtTof && !track.hasTOF()) {
        continue;
      }
      histos.fill(HIST("tofMassFull"), trackParCov.getPt(), mass);
      if (trackParCov.getPt() > antidPtTof && track.tofChi2() > 3) {
        continue;
      }
      histos.fill(HIST("tofMass"), trackParCov.getPt(), mass);

      if (trackParCov.getPt() <= antidPtTof || (trackParCov.getPt() > antidPtTof && hasTof && std::abs(mass - o2::constants::physics::MassDeuteron) < tofMassMax)) {
        histos.fill(HIST("tpcNsigmaGlo"), trackParCov.getPt(), nSigmaTPC);
        q1antid += 1.; // TODO: correct for efficiency
        q2antid += 1.;
      }
    }

    double q1L{0.}, q2L{0.}, q1antiL{0.}, q2antiL{0.};
    std::vector<int64_t> trkId;
    for (const auto& v0 : V0s) {
      if (v0.pt() < lambdaPtMin || v0.pt() > lambdaPtMax) {
        continue;
      }

      if (!selectLambda(v0)) {
        continue;
      }

      auto pos = v0.template posTrack_as<TracksFull>();
      auto neg = v0.template negTrack_as<TracksFull>();
      if (std::abs(pos.eta()) > etaMax || std::abs(pos.eta()) > etaMax) {
        continue;
      }

      bool matter = v0.alpha() > 0;

      histos.fill(HIST("massLambda"), matter ? v0.mLambda() : v0.mAntiLambda());
      histos.fill(HIST("cosPa"), v0.v0cosPA());
      histos.fill(HIST("radius"), v0.v0radius());
      histos.fill(HIST("dcaV0daugh"), v0.dcaV0daughters());
      histos.fill(HIST("dcaPosPv"), v0.dcapostopv());
      histos.fill(HIST("dcaNegPv"), v0.dcanegtopv());

      if (matter) {
        q1L += 1.; // TODO: correct for efficiency
        q2L += 1.;
      } else {
        q1antiL += 1.; // TODO: correct for efficiency
        q2antiL += 1.;
      }

      trkId.emplace_back(pos.globalIndex());
      trkId.emplace_back(neg.globalIndex());
    }

    // reject events having multiple v0s from same tracks (TODO: also across collisions?)
    std::sort(trkId.begin(), trkId.end());
    if (std::adjacent_find(trkId.begin(), trkId.end()) != trkId.end()) {
      return;
    }

    histos.fill(HIST("nEv"), collision.centFT0C());

    histos.fill(HIST("q1antid"), collision.centFT0C(), q1antid);
    histos.fill(HIST("q1sqantid"), collision.centFT0C(), std::pow(q1antid, 2));
    histos.fill(HIST("q2antid"), collision.centFT0C(), q2antid);

    histos.fill(HIST("q1L"), collision.centFT0C(), q1L);
    histos.fill(HIST("q1sqL"), collision.centFT0C(), std::pow(q1L, 2));
    histos.fill(HIST("q2L"), collision.centFT0C(), q2L);

    histos.fill(HIST("q1antiL"), collision.centFT0C(), q1antiL);
    histos.fill(HIST("q1sqantiL"), collision.centFT0C(), std::pow(q1antiL, 2));
    histos.fill(HIST("q2antiL"), collision.centFT0C(), q2antiL);

    histos.fill(HIST("q11LantiL"), collision.centFT0C(), q1L * q1antiL);
    histos.fill(HIST("q11Lantid"), collision.centFT0C(), q1L * q1antid);
    histos.fill(HIST("q11antiLantid"), collision.centFT0C(), q1antiL * q1antid);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<antidLambdaEbye>(cfgc)};
}
