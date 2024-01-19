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
#include "CommonConstants/PhysicsConstants.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGLF/DataModel/LFSlimStrangeTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPr, aod::pidTPCPi>;

struct LFStrangeTreeCreator {

  Produces<aod::LambdaTableML> lambdaTableML;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.05f, o2::constants::physics::MassLambda0 + 0.05f}, "binning for the lambda invariant-mass"};
  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};

  ConfigurableAxis cosPaAxis{"cosPaAxis", {1e3, 0.95f, 1.00f}, "binning for the cosPa axis"};
  ConfigurableAxis radiusAxis{"radiusAxis", {1e3, 0.f, 100.f}, "binning for the radius axis"};
  ConfigurableAxis dcaV0daughAxis{"dcaV0daughAxis", {2e2, 0.f, 2.f}, "binning for the dca of V0 daughters"};
  ConfigurableAxis dcaDaughPvAxis{"dcaDaughPvAxis", {1e3, 0.f, 10.f}, "binning for the dca of positive daughter to PV"};

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

  Configurable<float> lambdaPtMin{"lambdaPtMin", 0.5f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 3.0f, "maximum (anti)lambda pT (GeV/c)"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1f, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1f, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5f, "v0radius"};
  Configurable<float> lambdaMassCut{"lambdaMassCut", 0.05f, "maximum deviation from PDG mass"};

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

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    mRunNumber = 0;
    d_bz = 0;

    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    // v0 QA
    histos.add<TH1>("massLambda", ";#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", {HistType::kTH1F, {massLambdaAxis}});
    histos.add<TH1>("cosPa", ";cosPa;Entries", {HistType::kTH1F}, {cosPaAxis});
    histos.add<TH1>("radius", ";radius;Entries", {HistType::kTH1F}, {radiusAxis});
    histos.add<TH1>("dcaV0daugh", ";dcaV0daugh;Entries", {HistType::kTH1F}, {dcaV0daughAxis});
    histos.add<TH1>("dcaPosPv", ";dcaPosPv;Entries", {HistType::kTH1F}, {dcaDaughPvAxis});
    histos.add<TH1>("dcaNegPv", ";dcaNegPv;Entries", {HistType::kTH1F}, {dcaDaughPvAxis});
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

    histos.fill(HIST("zVtx"), collision.posZ());

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
      float ct = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      float tpcNsigmaPi = matter ? neg.tpcNSigmaPi() : pos.tpcNSigmaPi();
      float tpcNsigmaPr = matter ? pos.tpcNSigmaPr() : neg.tpcNSigmaPr();

      histos.fill(HIST("massLambda"), matter ? v0.mLambda() : v0.mAntiLambda());
      histos.fill(HIST("cosPa"), v0.v0cosPA());
      histos.fill(HIST("radius"), v0.v0radius());
      histos.fill(HIST("dcaV0daugh"), v0.dcaV0daughters());
      histos.fill(HIST("dcaPosPv"), v0.dcapostopv());
      histos.fill(HIST("dcaNegPv"), v0.dcanegtopv());

      lambdaTableML(v0.pt(),
                    v0.eta(),
                    collision.centFT0C(),
                    matter,
                    matter ? v0.mLambda() : v0.mAntiLambda(),
                    ct,
                    v0.v0radius(),
                    v0.dcav0topv(),
                    matter ? v0.dcanegtopv() : v0.dcapostopv(),
                    matter ? v0.dcapostopv() : v0.dcanegtopv(),
                    v0.dcaV0daughters(),
                    v0.v0cosPA(),
                    tpcNsigmaPi,
                    tpcNsigmaPr);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LFStrangeTreeCreator>(cfgc)};
}
