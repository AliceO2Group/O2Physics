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

#include <array>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"

#include "PWGLF/DataModel/LFHypernucleiTables.h"
#include "PWGLF/Utils/svCreator.h"

using namespace o2;
using namespace o2::framework;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"He3"};
std::shared_ptr<TH1> hGenPt;
} // namespace

struct hyperTester {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> fitter;

  int dau0pdg{211};
  int dau1pdg{1000020030};
  int motherpdg{1010010030};
  svCreator svfitter{dau0pdg, dau1pdg};
  HistogramRegistry qaRegistry{"QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};
  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", true, "If true, divide tpcInnerParam by the electric charge"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  std::array<float, 6> mBBparamsHe;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (candidate.tpcNClsFound() < 70 ||
        candidate.itsNCls() < 2) {
      return false;
    }

    return true;
  }

  template <typename T>
  float computeNSigmaHe3(const T& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;

    float correctedTPCinnerParam = (heliumPID && cfgCompensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4]);

    double resoTPC{expTPCSignal * mBBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename T>
  bool selectionPIDHe3(const T& candidate)
  {
    auto nSigmaHe3 = computeNSigmaHe3(candidate);
    if (std::abs(nSigmaHe3) < 4) {
      return true;
    }
    return false;
  }

  void init(InitContext&)
  {
    LOG(info) << "Initializing hyperTester";
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    for (int i = 0; i < 5; i++) {
      mBBparamsHe[i] = cfgBetheBlochParams->get("He3", Form("p%i", i));
    }
    mBBparamsHe[5] = cfgBetheBlochParams->get("He3", "resolution");
    LOG(info) << "He3 PID params: " << mBBparamsHe[0] << " " << mBBparamsHe[1] << " " << mBBparamsHe[2] << " " << mBBparamsHe[3] << " " << mBBparamsHe[4] << " " << mBBparamsHe[5];
    svfitter.setFitter(fitter);

    hGenPt = qaRegistry.add<TH1>("hPt", ";pT (GeV/#it{c}); ", HistType::kTH1D, {{60, 0, 10}});
  }

  // loop over data frames
  void process(TracksFull const& tracks, aod::Collisions const& collisions, aod::AmbiguousTracks const& ambiguousTracks, aod::BCsWithTimestamps const& bcs, aod::McTrackLabels const& trackLabels, aod::McParticles const& particlesMC)
  {
    svfitter.clearPools();

    for (auto& track : tracks) {
      // LOG(info) << "Processing track";
      if (!selectionTrack(track)) {
        continue;
      }

      int pdgHypo = selectionPIDHe3(track) ? dau1pdg : dau0pdg;
      svfitter.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }
    LOG(info) << "Processing done!";
    for (int i = 0; i < svfitter.getTrackCandPool().size(); i++) {
      LOG(info) << "Pool " << i << " size: " << svfitter.getTrackCandPool()[i].size();
    }
    auto& svPool = svfitter.getSVCandPool(collisions);
    LOG(info) << "SV pool size: " << svPool.size();

    for (auto& svCand : svPool) {
      auto mcLab0 = trackLabels.rawIteratorAt(svCand.tr0Idx);
      auto mcLab1 = trackLabels.rawIteratorAt(svCand.tr1Idx);

      if (mcLab0.has_mcParticle() && mcLab1.has_mcParticle()) {
        auto mcTrack0 = mcLab0.mcParticle_as<aod::McParticles>();
        auto mcTrack1 = mcLab1.mcParticle_as<aod::McParticles>();
        // LOG(info) << "track0 pdg: " << mcTrack0.pdgCode() << " track1 pdg: " << mcTrack1.pdgCode();
        if (mcTrack0.has_mothers() && mcTrack1.has_mothers()) {
          for (auto& mother0 : mcTrack0.mothers_as<aod::McParticles>()) {
            for (auto& mother1 : mcTrack1.mothers_as<aod::McParticles>()) {
              if (mother1.globalIndex() != mother0.globalIndex())
                continue;
              if (abs(mcTrack0.pdgCode()) != dau0pdg || abs(mcTrack1.pdgCode()) != dau1pdg)
                continue;
              if (std::abs(mother1.pdgCode()) != motherpdg)
                continue;
              hGenPt->Fill(mother1.pt());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hyperTester>(cfgc)};
}