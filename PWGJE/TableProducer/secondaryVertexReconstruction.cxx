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

/// \file secondaryVertexReconstruction.cxx
/// \brief Task to produce a SV indices table joinable to the jet tables and 3/2-prong SV for hf jet tagging
/// \author Hadi Hassan <hadi.hassan@cern.ch>

#include <type_traits>
#include <string>
#include <utility>
#include <vector>

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

  Configurable<float> magneticField{"magneticField", 20.0f, "magnetic field in kG"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> maxRsv{"maxRsv", 999., "max. radius of the reconstruced SV"};
  Configurable<double> maxZsv{"maxZsv", 999., "max. Z coordinates of the reconstruced SV"};
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
      AxisSpec nProngsBins = {6, 2, 8, "nProngs"};
      registry.add("hMassNProngs", "n-prong candidates;SV inv. mass (GeV/#it{c}^{2});nProngs;entries", {HistType::kTH2F, {{500, 0., 30}, nProngsBins}});
      registry.add("hLxySNProngs", "n-prong Decay length in XY;S#it{L}_{xy};nProngs;entries", {HistType::kTH2F, {{100, 0., 100.0}, nProngsBins}});
      registry.add("hLSNProngs", "n-prong Decay length in 3D;S#it{L};nProngs;entries", {HistType::kTH2F, {{100, 0., 100.0}, nProngsBins}});
      registry.add("hFeNProngs", "n-prong Energy fraction carried by the SV from the jet;#it{f}_{E};nProngs;entries", {HistType::kTH2F, {{100, 0., 1.0}, nProngsBins}});
      registry.add("hDcaXYNProngs", "DCAxy of n-prong candidate daughters;#it{p}_{T} (GeV/#it{c});#it{d}_{xy} (#mum);nProngs;entries", {HistType::kTH3F, {{100, 0., 20.}, {200, -500., 500.}, nProngsBins}});
      registry.add("hDcaZNProngs", "DCAz of n-prong candidate daughters;#it{p}_{T} (GeV/#it{c});#it{d}_{z} (#mum);nProngs;entries", {HistType::kTH3F, {{100, 0., 20.}, {200, -500., 500.}, nProngsBins}});
      registry.add("hDispersion", "Vertex dispersion;#sigma_{vtx};nProngs;entries", {HistType::kTH2F, {{200, 0., 1.0}, nProngsBins}});
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

  using JetCollisionwPIs = soa::Join<aod::JetCollisions, aod::JCollisionPIs>;
  using JetTracksData = soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackPIs>>;
  using JetTracksMCDwPIs = soa::Filtered<soa::Join<aod::JetTracksMCD, aod::JTrackPIs>>;
  using OriginalTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov>;

  template <unsigned int numProngs, bool externalMagneticField, typename AnyCollision, typename AnyJet, typename AnyParticles>
  void runCreatorNProng(AnyCollision const& collision,
                        AnyJet const& analysisJet,
                        AnyParticles const& listoftracks,
                        std::vector<int>& svIndices,
                        o2::vertexing::DCAFitterN<numProngs>& df,
                        size_t prongIndex = 0,
                        std::vector<size_t> currentCombination = {})
  {

    const auto& particles = analysisJet.template tracks_as<AnyParticles>();

    if (currentCombination.size() == numProngs) {
      // Create an array of track parameters and covariance matrices for the current combination
      std::array<o2::track::TrackParametrizationWithError<float>, numProngs> trackParVars;
      double energySV = 0.;
      for (unsigned int inum = 0; inum < numProngs; ++inum) {
        const auto& prong = particles[currentCombination[inum]].template track_as<OriginalTracks>();
        energySV += prong.energy(o2::constants::physics::MassPiPlus);
        trackParVars[inum] = getTrackParCov(prong);
      }

      if constexpr (externalMagneticField) {
        bz = magneticField;
      } else {
        auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
        if (runNumber != bc.runNumber()) {
          initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut, false);
          bz = o2::base::Propagator::Instance()->getNominalBz();
        }
      }

      // Use a different fitter depending on the number of prongs
      df.setBz(bz);

      // Reconstruct the secondary vertex
      int processResult = 0;
      try {
        std::apply([&df, &processResult](const auto&... elems) { processResult = df.process(elems...); }, trackParVars);
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        return;
      }
      if (processResult == 0) {
        return;
      }

      const auto& secondaryVertex = df.getPCACandidatePos();
      if (std::sqrt(secondaryVertex[0] * secondaryVertex[0] + secondaryVertex[1] * secondaryVertex[1]) > maxRsv || std::abs(secondaryVertex[2]) > maxZsv) {
        return;
      }

      float dispersion = 0.;
      for (unsigned int inum = 0; inum < numProngs; ++inum) {
        o2::dataformats::VertexBase sv(o2::math_utils::Point3D<float>{secondaryVertex[0], secondaryVertex[1], secondaryVertex[2]}, std::array<float, 6>{0});
        o2::dataformats::DCA dcaSV;
        auto& prong = df.getTrack(inum);
        prong.propagateToDCA(sv, bz, &dcaSV);
        dispersion += (dcaSV.getY() * dcaSV.getY() + dcaSV.getZ() * dcaSV.getZ());
      }
      dispersion = std::sqrt(dispersion / numProngs);

      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();

      registry.fill(HIST("hDispersion"), dispersion, numProngs);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();

      // Get track momenta and impact parameters
      std::array<std::array<float, 3>, numProngs> arrayMomenta;
      std::array<o2::dataformats::DCA, numProngs> impactParameters;
      for (unsigned int inum = 0; inum < numProngs; ++inum) {
        trackParVars[inum].getPxPyPzGlo(arrayMomenta[inum]);
        trackParVars[inum].propagateToDCA(primaryVertex, bz, &impactParameters[inum]);

        if (fillHistograms) {
          const auto& prong = particles[currentCombination[inum]].template track_as<OriginalTracks>();
          registry.fill(HIST("hDcaXYNProngs"), prong.pt(), impactParameters[inum].getY() * toMicrometers, numProngs);
          registry.fill(HIST("hDcaZNProngs"), prong.pt(), impactParameters[inum].getZ() * toMicrometers, numProngs);
        }
      }

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      // calculate invariant mass
      std::array<double, numProngs> massArray;
      std::fill(massArray.begin(), massArray.end(), o2::constants::physics::MassPiPlus);
      double massSV = RecoDecay::m(std::move(arrayMomenta), massArray);

      // fill candidate table rows
      if ((doprocessData3Prongs || doprocessData3ProngsExternalMagneticField) && numProngs == 3) {
        sv3prongTableData(analysisJet.globalIndex(),
                          primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                          secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                          arrayMomenta[0][0] + arrayMomenta[1][0] + arrayMomenta[2][0],
                          arrayMomenta[0][1] + arrayMomenta[1][1] + arrayMomenta[2][1],
                          arrayMomenta[0][2] + arrayMomenta[1][2] + arrayMomenta[2][2],
                          energySV, massSV, chi2PCA, dispersion, errorDecayLength, errorDecayLengthXY);
        svIndices.push_back(sv3prongTableData.lastIndex());
      } else if ((doprocessData2Prongs || doprocessData2ProngsExternalMagneticField) && numProngs == 2) {
        sv2prongTableData(analysisJet.globalIndex(),
                          primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                          secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                          arrayMomenta[0][0] + arrayMomenta[1][0],
                          arrayMomenta[0][1] + arrayMomenta[1][1],
                          arrayMomenta[0][2] + arrayMomenta[1][2],
                          energySV, massSV, chi2PCA, dispersion, errorDecayLength, errorDecayLengthXY);
        svIndices.push_back(sv2prongTableData.lastIndex());
      } else if ((doprocessMCD3Prongs || doprocessMCD3ProngsExternalMagneticField) && numProngs == 3) {
        sv3prongTableMCD(analysisJet.globalIndex(),
                         primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                         secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                         arrayMomenta[0][0] + arrayMomenta[1][0] + arrayMomenta[2][0],
                         arrayMomenta[0][1] + arrayMomenta[1][1] + arrayMomenta[2][1],
                         arrayMomenta[0][2] + arrayMomenta[1][2] + arrayMomenta[2][2],
                         energySV, massSV, chi2PCA, dispersion, errorDecayLength, errorDecayLengthXY);
        svIndices.push_back(sv3prongTableMCD.lastIndex());
      } else if ((doprocessMCD2Prongs || doprocessMCD2ProngsExternalMagneticField) && numProngs == 2) {
        sv2prongTableMCD(analysisJet.globalIndex(),
                         primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                         secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                         arrayMomenta[0][0] + arrayMomenta[1][0],
                         arrayMomenta[0][1] + arrayMomenta[1][1],
                         arrayMomenta[0][2] + arrayMomenta[1][2],
                         energySV, massSV, chi2PCA, dispersion, errorDecayLength, errorDecayLengthXY);
        svIndices.push_back(sv2prongTableMCD.lastIndex());
      } else {
        LOG(error) << "No process specified\n";
      }

      // fill histograms
      if (fillHistograms) {
        double decayLengthNormalised = RecoDecay::distance(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, std::array{secondaryVertex[0], secondaryVertex[1], secondaryVertex[2]}) / errorDecayLength;
        double decayLengthXYNormalised = RecoDecay::distanceXY(std::array{primaryVertex.getX(), primaryVertex.getY()}, std::array{secondaryVertex[0], secondaryVertex[1]}) / errorDecayLengthXY;

        registry.fill(HIST("hMassNProngs"), massSV, numProngs);
        registry.fill(HIST("hLxySNProngs"), decayLengthXYNormalised, numProngs);
        registry.fill(HIST("hLSNProngs"), decayLengthNormalised, numProngs);
        registry.fill(HIST("hFeNProngs"), energySV / analysisJet.energy() > 1. ? 0.99 : energySV / analysisJet.energy(), numProngs);
      }

      return;
    }

    // Recursive call to explore all combinations
    for (size_t iprong = prongIndex; iprong < particles.size(); ++iprong) {

      const auto& testTrack = particles[iprong].template track_as<OriginalTracks>();
      if (testTrack.pt() < ptMinTrack || testTrack.eta() < etaMinTrack || testTrack.eta() > etaMaxTrack) {
        continue;
      }

      currentCombination.push_back(iprong);
      runCreatorNProng<numProngs, externalMagneticField>(
        collision, analysisJet, listoftracks, svIndices, df, iprong + 1, currentCombination);
      currentCombination.pop_back();
    }
  }

  void processDummy(JetCollisionwPIs::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processDummy, "Dummy process", true);

  void processData3Prongs(JetCollisionwPIs::iterator const& collision, aod::Collisions const& /*realColl*/, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, JetTracksData const& tracks, OriginalTracks const& /*tracks*/, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    for (const auto& jet : jets) {
      std::vector<int> svIndices;
      runCreatorNProng<3, false>(collision.template collision_as<aod::Collisions>(), jet, tracks, svIndices, df3);
      sv3prongIndicesTableData(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processData3Prongs, "Reconstruct the data 3-prong secondary vertex", false);

  void processData3ProngsExternalMagneticField(JetCollisionwPIs::iterator const& collision, aod::Collisions const& /*realColl*/, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, JetTracksData const& tracks, OriginalTracks const& /*tracks*/)
  {
    for (const auto& jet : jets) {
      std::vector<int> svIndices;
      runCreatorNProng<3, true>(collision.template collision_as<aod::Collisions>(), jet, tracks, svIndices, df3);
      sv3prongIndicesTableData(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processData3ProngsExternalMagneticField, "Reconstruct the data 3-prong secondary vertex with external magnetic field", false);

  void processData2Prongs(JetCollisionwPIs::iterator const& collision, aod::Collisions const& /*realColl*/, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, JetTracksData const& tracks, OriginalTracks const& /*tracks*/, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    for (const auto& jet : jets) {
      std::vector<int> svIndices;
      runCreatorNProng<2, false>(collision.template collision_as<aod::Collisions>(), jet, tracks, svIndices, df2);
      sv2prongIndicesTableData(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processData2Prongs, "Reconstruct the data 2-prong secondary vertex", false);

  void processData2ProngsExternalMagneticField(JetCollisionwPIs::iterator const& collision, aod::Collisions const& /*realColl*/, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, JetTracksData const& tracks, OriginalTracks const& /*tracks*/)
  {
    for (const auto& jet : jets) {
      std::vector<int> svIndices;
      runCreatorNProng<2, true>(collision.template collision_as<aod::Collisions>(), jet, tracks, svIndices, df2);
      sv2prongIndicesTableData(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processData2ProngsExternalMagneticField, "Reconstruct the data 2-prong secondary vertex with extrernal magnetic field", false);

  void processMCD3Prongs(JetCollisionwPIs::iterator const& collision, aod::Collisions const& /*realColl*/, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& mcdjets, JetTracksMCDwPIs const& tracks, OriginalTracks const& /*tracks*/, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    for (const auto& jet : mcdjets) {
      std::vector<int> svIndices;
      runCreatorNProng<3, false>(collision.template collision_as<aod::Collisions>(), jet, tracks, svIndices, df3);
      sv3prongIndicesTableMCD(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processMCD3Prongs, "Reconstruct the MCD 3-prong secondary vertex", false);

  void processMCD3ProngsExternalMagneticField(JetCollisionwPIs::iterator const& collision, aod::Collisions const& /*realColl*/, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& mcdjets, JetTracksMCDwPIs const& tracks, OriginalTracks const& /*tracks*/)
  {
    for (const auto& jet : mcdjets) {
      std::vector<int> svIndices;
      runCreatorNProng<3, true>(collision.template collision_as<aod::Collisions>(), jet, tracks, svIndices, df3);
      sv3prongIndicesTableMCD(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processMCD3ProngsExternalMagneticField, "Reconstruct the MCD 3-prong secondary vertex with external magnetic field", false);

  void processMCD2Prongs(JetCollisionwPIs::iterator const& collision, aod::Collisions const& /*realColl*/, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& mcdjets, JetTracksMCDwPIs const& tracks, OriginalTracks const& /*tracks*/, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    for (const auto& jet : mcdjets) {
      std::vector<int> svIndices;
      runCreatorNProng<2, false>(collision.template collision_as<aod::Collisions>(), jet, tracks, svIndices, df2);
      sv2prongIndicesTableMCD(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processMCD2Prongs, "Reconstruct the MCD 2-prong secondary vertex", false);

  void processMCD2ProngsExternalMagneticField(JetCollisionwPIs::iterator const& collision, aod::Collisions const& /*realColl*/, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& mcdjets, JetTracksMCDwPIs const& tracks, OriginalTracks const& /*tracks*/)
  {
    for (const auto& jet : mcdjets) {
      std::vector<int> svIndices;
      runCreatorNProng<2, true>(collision.template collision_as<aod::Collisions>(), jet, tracks, svIndices, df2);
      sv2prongIndicesTableMCD(svIndices);
    }
  }
  PROCESS_SWITCH(SecondaryVertexReconstruction, processMCD2ProngsExternalMagneticField, "Reconstruct the MCD 2-prong secondary vertex with external magnetic field", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SecondaryVertexReconstruction>(cfgc, TaskName{"jet-sv-reconstruction-charged"})}; // o2-linter: disable=name/o2-task
}
