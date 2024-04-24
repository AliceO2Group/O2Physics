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

// Task to produce a SV indices table joinable to the jet tables and 3/2-prong SV for hf jet tagging
//
/// \author Hadi Hassan <hadi.hassan@cern.ch>

#include <type_traits>

#include <TF1.h>
#include <TH1.h>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"

#include "ReconstructionDataFormats/DCA.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// TODO track to collision association

struct SecondaryVertexReconstruction {

  Produces<aod::DataSecondaryVertex3Prongs> sv3prongTableData;
  Produces<aod::DataSecondaryVertex3ProngIndices> sv3prongIndicesTableData;

  Produces<aod::DataSecondaryVertex2Prongs> sv2prongTableData;
  Produces<aod::DataSecondaryVertex2ProngIndices> sv2prongIndicesTableData;

  Produces<aod::MCDSecondaryVertex3Prongs> sv3prongTableMCD;
  Produces<aod::MCDSecondaryVertex3ProngIndices> sv3prongIndicesTableMCD;

  Produces<aod::MCDSecondaryVertex2Prongs> sv2prongTableMCD;
  Produces<aod::MCDSecondaryVertex2ProngIndices> sv2prongIndicesTableMCD;

  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<float> ptMinTrack{"ptMinTrack", -1., "min. track pT"};
  Configurable<float> etaMinTrack{"etaMinTrack", -99999., "min. pseudorapidity"};
  Configurable<float> etaMaxTrack{"etaMaxTrack", 4., "max. pseudorapidity"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "do validation plots"};

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter
  o2::vertexing::DCAFitterN<3> df3; // 3-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;

  HistogramRegistry registry{"registry"};

  int runNumber{0};
  float toMicrometers = 10000.; // from cm to µm
  double bz{0.};

  void init(InitContext const&)
  {
    if (fillHistograms) {
      registry.add("hMass3Prongs", "3-prong candidates;SV inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 30}}});
      registry.add("hLxyS3Prongs", "3-prong Decay length in XY;S#it{L}_{xy};entries", {HistType::kTH1F, {{100, 0., 100.0}}});
      registry.add("hLS3Prongs", "3-prong Decay length in 3D;S#it{L};entries", {HistType::kTH1F, {{100, 0., 100.0}}});
      registry.add("hFe3Prongs", "3-prong Energy fraction carried by the SV from the jet;#it{f}_{E};entries", {HistType::kTH1F, {{100, 0., 1.0}}});
      registry.add("hDcaXY3Prongs", "DCAxy of 3-prong candidate daughters;#it{p}_{T} (GeV/#it{c});#it{d}_{xy} (#mum);entries", {HistType::kTH2F, {{100, 0., 20.}, {200, -500., 500.}}});
      registry.add("hDcaZ3Prongs", "DCAz of 3-prong candidate daughters;#it{p}_{T} (GeV/#it{c});#it{d}_{z} (#mum);entries", {HistType::kTH2F, {{100, 0., 20.}, {200, -500., 500.}}});
      registry.add("hMass2Prongs", "2-prong candidates;SV inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 30}}});
      registry.add("hLxyS2Prongs", "2-prong Decay length in XY;S#it{L}_{xy};entries", {HistType::kTH1F, {{100, 0., 100.0}}});
      registry.add("hLS2Prongs", "2-prong Decay length in 3D;S#it{L};entries", {HistType::kTH1F, {{100, 0., 100.0}}});
      registry.add("hFe2Prongs", "2-prong Energy fraction carried by the SV from the jet;#it{f}_{E};entries", {HistType::kTH1F, {{100, 0., 1.0}}});
      registry.add("hDcaXY2Prongs", "DCAxy of 2-prong candidate daughters;#it{p}_{T} (GeV/#it{c});#it{d}_{xy} (#mum);entries", {HistType::kTH2F, {{100, 0., 20.}, {200, -500., 500.}}});
      registry.add("hDcaZ2Prongs", "DCAz of 2-prong candidate daughters;#it{p}_{T} (GeV/#it{c});#it{d}_{z} (#mum);entries", {HistType::kTH2F, {{100, 0., 20.}, {200, -500., 500.}}});
    }

    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(useAbsDCA);
    df3.setWeightedFinalPCA(useWeightedFinalPCA);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  Filter trackCuts = (aod::jtrack::pt > ptMinTrack && aod::jtrack::eta > etaMinTrack && aod::jtrack::eta < etaMaxTrack);

  using JetCollisionwPIs = soa::Join<JetCollisions, aod::JCollisionPIs>;
  using JetTracksData = soa::Filtered<soa::Join<JetTracks, aod::JTrackPIs>>;
  using JetTracksMCDwPIs = soa::Filtered<soa::Join<JetTracksMCD, aod::JTrackPIs>>;
  using OriginalTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov>;

  template <typename AnyCollision, typename AnyJet, typename AnyParticles>
  void runCreator2Prong(AnyCollision const& collision,
                        AnyJet const& analysisJet,
                        AnyParticles const& listoftracks,
                        std::vector<int>& svIndices)
  {

    const auto& particles = analysisJet.template tracks_as<AnyParticles>();
    for (size_t iprong0 = 0; iprong0 < particles.size(); ++iprong0) {
      const auto& prong0 = particles[iprong0].template track_as<OriginalTracks>();

      if (prong0.pt() < ptMinTrack || prong0.eta() < etaMinTrack || prong0.eta() > etaMaxTrack) {
        continue;
      }

      for (size_t iprong1 = iprong0 + 1; iprong1 < particles.size(); ++iprong1) {
        const auto& prong1 = particles[iprong1].template track_as<OriginalTracks>();

        if (prong1.pt() < ptMinTrack || prong1.eta() < etaMinTrack || prong1.eta() > etaMaxTrack) {
          continue;
        }

        auto trackParVar0 = getTrackParCov(prong0);
        auto trackParVar1 = getTrackParCov(prong1);

        auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
        if (runNumber != bc.runNumber()) {
          initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
          bz = o2::base::Propagator::Instance()->getNominalBz();
        }
        df2.setBz(bz);

        // reconstruct the 2-prong secondary vertex
        if (df2.process(trackParVar0, trackParVar1) == 0) {
          continue;
        }

        const auto& secondaryVertex = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();
        trackParVar0 = df2.getTrack(0);
        trackParVar1 = df2.getTrack(1);

        // get track momenta
        std::array<float, 3> pvec0;
        std::array<float, 3> pvec1;
        trackParVar0.getPxPyPzGlo(pvec0);
        trackParVar1.getPxPyPzGlo(pvec1);

        // get track impact parameters
        // This modifies track momenta!
        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();

        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
        trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
        if (fillHistograms) {
          registry.fill(HIST("hDcaXY2Prongs"), prong0.pt(), impactParameter0.getY() * toMicrometers);
          registry.fill(HIST("hDcaXY2Prongs"), prong1.pt(), impactParameter1.getY() * toMicrometers);
          registry.fill(HIST("hDcaZ2Prongs"), prong0.pt(), impactParameter0.getZ() * toMicrometers);
          registry.fill(HIST("hDcaZ2Prongs"), prong1.pt(), impactParameter1.getZ() * toMicrometers);
        }

        // get uncertainty of the decay length
        double phi, theta;
        getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        // calculate invariant mass
        auto arrayMomenta = std::array{pvec0, pvec1};
        double energySV = prong0.energy(o2::constants::physics::MassPiPlus) + prong1.energy(o2::constants::physics::MassPiPlus);
        double massSV = RecoDecay::m(std::move(arrayMomenta), std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});

        // fill candidate table rows
        if (doprocessData2Prongs) {
          sv2prongTableData(analysisJet.globalIndex(),
                            primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                            secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                            pvec0[0] + pvec1[0],
                            pvec0[1] + pvec1[1],
                            pvec0[2] + pvec1[2],
                            energySV, massSV, chi2PCA, errorDecayLength, errorDecayLengthXY);
          svIndices.push_back(sv2prongTableData.lastIndex());
        } else {
          sv2prongTableMCD(analysisJet.globalIndex(),
                           primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                           secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                           pvec0[0] + pvec1[0],
                           pvec0[1] + pvec1[1],
                           pvec0[2] + pvec1[2],
                           energySV, massSV, chi2PCA, errorDecayLength, errorDecayLengthXY);
          svIndices.push_back(sv2prongTableMCD.lastIndex());
        }

        // fill histograms
        if (fillHistograms) {

          double DecayLengthNormalised = RecoDecay::distance(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, std::array{secondaryVertex[0], secondaryVertex[1], secondaryVertex[2]}) / errorDecayLength;
          double DecayLengthXYNormalised = RecoDecay::distanceXY(std::array{primaryVertex.getX(), primaryVertex.getY()}, std::array{secondaryVertex[0], secondaryVertex[1]}) / errorDecayLengthXY;

          registry.fill(HIST("hMass2Prongs"), massSV);
          registry.fill(HIST("hLxyS2Prongs"), DecayLengthXYNormalised);
          registry.fill(HIST("hLS2Prongs"), DecayLengthNormalised);
          registry.fill(HIST("hFe2Prongs"), energySV / analysisJet.energy() > 1. ? 0.99 : energySV / analysisJet.energy());
        }
      }
    }
  }

  template <typename AnyCollision, typename AnyJet, typename AnyParticles>
  void runCreator3Prong(AnyCollision const& collision,
                        AnyJet const& analysisJet,
                        AnyParticles const& listoftracks,
                        std::vector<int>& svIndices)
  {
    const auto& particles = analysisJet.template tracks_as<AnyParticles>();
    for (size_t iprong0 = 0; iprong0 < particles.size(); ++iprong0) {
      const auto& prong0 = particles[iprong0].template track_as<OriginalTracks>();

      if (prong0.pt() < ptMinTrack || prong0.eta() < etaMinTrack || prong0.eta() > etaMaxTrack) {
        continue;
      }

      for (size_t iprong1 = iprong0 + 1; iprong1 < particles.size(); ++iprong1) {
        const auto& prong1 = particles[iprong1].template track_as<OriginalTracks>();

        if (prong1.pt() < ptMinTrack || prong1.eta() < etaMinTrack || prong1.eta() > etaMaxTrack) {
          continue;
        }

        for (size_t iprong2 = iprong1 + 1; iprong2 < particles.size(); ++iprong2) {
          const auto& prong2 = particles[iprong2].template track_as<OriginalTracks>();

          if (prong2.pt() < ptMinTrack || prong2.eta() < etaMinTrack || prong2.eta() > etaMaxTrack) {
            continue;
          }

          auto trackParVar0 = getTrackParCov(prong0);
          auto trackParVar1 = getTrackParCov(prong1);
          auto trackParVar2 = getTrackParCov(prong2);

          auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
          if (runNumber != bc.runNumber()) {
            initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
            bz = o2::base::Propagator::Instance()->getNominalBz();
          }
          df3.setBz(bz);

          // reconstruct the 3-prong secondary vertex
          if (df3.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
            continue;
          }

          const auto& secondaryVertex = df3.getPCACandidate();
          auto chi2PCA = df3.getChi2AtPCACandidate();
          auto covMatrixPCA = df3.calcPCACovMatrixFlat();
          trackParVar0 = df3.getTrack(0);
          trackParVar1 = df3.getTrack(1);
          trackParVar2 = df3.getTrack(2);

          // get track momenta
          std::array<float, 3> pvec0;
          std::array<float, 3> pvec1;
          std::array<float, 3> pvec2;
          trackParVar0.getPxPyPzGlo(pvec0);
          trackParVar1.getPxPyPzGlo(pvec1);
          trackParVar2.getPxPyPzGlo(pvec2);

          // get track impact parameters
          // This modifies track momenta!
          auto primaryVertex = getPrimaryVertex(collision);
          auto covMatrixPV = primaryVertex.getCov();

          o2::dataformats::DCA impactParameter0;
          o2::dataformats::DCA impactParameter1;
          o2::dataformats::DCA impactParameter2;
          trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
          trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
          trackParVar2.propagateToDCA(primaryVertex, bz, &impactParameter2);
          if (fillHistograms) {
            registry.fill(HIST("hDcaXY3Prongs"), prong0.pt(), impactParameter0.getY() * toMicrometers);
            registry.fill(HIST("hDcaXY3Prongs"), prong1.pt(), impactParameter1.getY() * toMicrometers);
            registry.fill(HIST("hDcaXY3Prongs"), prong2.pt(), impactParameter2.getY() * toMicrometers);
            registry.fill(HIST("hDcaZ3Prongs"), prong0.pt(), impactParameter0.getZ() * toMicrometers);
            registry.fill(HIST("hDcaZ3Prongs"), prong1.pt(), impactParameter1.getZ() * toMicrometers);
            registry.fill(HIST("hDcaZ3Prongs"), prong2.pt(), impactParameter2.getZ() * toMicrometers);
          }

          // get uncertainty of the decay length
          double phi, theta;
          getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          // calculate invariant mass
          auto arrayMomenta = std::array{pvec0, pvec1, pvec2};
          double energySV = prong0.energy(o2::constants::physics::MassPiPlus) + prong1.energy(o2::constants::physics::MassPiPlus) + prong2.energy(o2::constants::physics::MassPiPlus);
          double massSV = RecoDecay::m(std::move(arrayMomenta), std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});

          // fill candidate table rows
          if (doprocessData3Prongs) {
            sv3prongTableData(analysisJet.globalIndex(),
                              primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                              secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                              pvec0[0] + pvec1[0] + pvec2[0],
                              pvec0[1] + pvec1[1] + pvec2[1],
                              pvec0[2] + pvec1[2] + pvec2[2],
                              energySV, massSV, chi2PCA, errorDecayLength, errorDecayLengthXY);
            svIndices.push_back(sv3prongTableData.lastIndex());
          } else {
            sv3prongTableMCD(analysisJet.globalIndex(),
                             primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                             secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                             pvec0[0] + pvec1[0] + pvec2[0],
                             pvec0[1] + pvec1[1] + pvec2[1],
                             pvec0[2] + pvec1[2] + pvec2[2],
                             energySV, massSV, chi2PCA, errorDecayLength, errorDecayLengthXY);
            svIndices.push_back(sv3prongTableMCD.lastIndex());
          }

          // fill histograms
          if (fillHistograms) {
            double DecayLengthNormalised = RecoDecay::distance(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, std::array{secondaryVertex[0], secondaryVertex[1], secondaryVertex[2]}) / errorDecayLength;
            double DecayLengthXYNormalised = RecoDecay::distanceXY(std::array{primaryVertex.getX(), primaryVertex.getY()}, std::array{secondaryVertex[0], secondaryVertex[1]}) / errorDecayLengthXY;

            registry.fill(HIST("hMass3Prongs"), massSV);
            registry.fill(HIST("hLxyS3Prongs"), DecayLengthXYNormalised);
            registry.fill(HIST("hLS3Prongs"), DecayLengthNormalised);
            registry.fill(HIST("hFe3Prongs"), energySV / analysisJet.energy() > 1. ? 0.99 : energySV / analysisJet.energy());
          }
        }
      }
    }
  }

  void processDummy(JetCollisionwPIs::iterator const& collision)
  {
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processDummy, "Dummy process", true);

  void processData3Prongs(JetCollisionwPIs::iterator const& collision, aod::Collisions const& realColl, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, JetTracksData const& jtracks, OriginalTracks const& tracks, aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    for (auto& jet : jets) {
      std::vector<int> svIndices;
      runCreator3Prong(collision.template collision_as<aod::Collisions>(), jet, jtracks, svIndices);
      sv3prongIndicesTableData(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processData3Prongs, "Reconstruct the data 3-prong secondary vertex", false);

  void processData2Prongs(JetCollisionwPIs::iterator const& collision, aod::Collisions const& realColl, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, JetTracksData const& jtracks, OriginalTracks const& tracks, aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    for (auto& jet : jets) {
      std::vector<int> svIndices;
      runCreator2Prong(collision.template collision_as<aod::Collisions>(), jet, jtracks, svIndices);
      sv2prongIndicesTableData(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processData2Prongs, "Reconstruct the data 2-prong secondary vertex", false);

  void processMCD3Prongs(JetCollisionwPIs::iterator const& collision, aod::Collisions const& realColl, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& mcdjets, JetTracksMCDwPIs const& jtracks, OriginalTracks const& tracks, aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    for (auto& jet : mcdjets) {
      std::vector<int> svIndices;
      runCreator3Prong(collision.template collision_as<aod::Collisions>(), jet, jtracks, svIndices);
      sv3prongIndicesTableMCD(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processMCD3Prongs, "Reconstruct the MCD 3-prong secondary vertex", false);

  void processMCD2Prongs(JetCollisionwPIs::iterator const& collision, aod::Collisions const& realColl, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& mcdjets, JetTracksMCDwPIs const& jtracks, OriginalTracks const& tracks, aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    for (auto& jet : mcdjets) {
      std::vector<int> svIndices;
      runCreator2Prong(collision.template collision_as<aod::Collisions>(), jet, jtracks, svIndices);
      sv2prongIndicesTableMCD(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processMCD2Prongs, "Reconstruct the MCD 2-prong secondary vertex", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SecondaryVertexReconstruction>(cfgc, TaskName{"jet-sv-reconstruction-charged"})};
}
