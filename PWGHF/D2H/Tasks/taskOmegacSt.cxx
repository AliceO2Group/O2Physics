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

/// \file taskOmegacSt.cxx
/// \brief Task to reconstruct Ωc from strangeness-tracked Ω and pion
///
/// \author Jochen Klein

#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TRandom.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskOmegacSt {
  Configurable<double> bz{"bz", -5., "magnetic field"};
  Configurable<int> materialCorrectionType{"materialCorrectionType", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpMagPath{"grpMagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> matLutPath{"matLutPath", "GLO/Param/MatLUT", "Path of the material LUT"};
  Configurable<bool> propToDCA{"propToDCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int runNumber;
  bool bzOnly = true;

  using TracksExt = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels>;

  HistogramRegistry registry{
    "registry",
    {
      {"hDca", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"hDcaXY", "DCA;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"hDcaZ", "DCA;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"hDcaVsPt", "DCA;DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}}}},
      {"hDcaVsR", "DCA;DCA (cm);R (cm)", {HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}}}},
      {"hDecayLength", "Decay length;L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthId", "Decay length (true #Omega_{c});L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthGen", "Decay length (gen);L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDeltaDecayLength", "#Delta decay length (gen);#Delta L (#mum)", {HistType::kTH1D, {{200, -250., 250.}}}},
      {"hDecayLengthScaled", "Decay length * M/p;L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthScaledId", "Decay length * M/p (true #Omega_{c});L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthScaledGen", "Decay length * M/p (gen);L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hMassOmegac", "inv. mass #Omega + #pi;inv. mass (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 1.5, 3.}}}},
      {"hMassOmegacVsPt", "inv. mass #Omega + #pi;inv. mass (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", {HistType::kTH2D, {{400, 1.5, 3.}, {10, 0., 10.}}}},
      {"hMassOmegacGen", "inv. mass #Omega + #pi (from MC);inv. mass (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 1.5, 3.}}}},
      {"hMassVsPt", "DCA;Mass (GeV/#it{c}^2);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., 10.}, {200, 0., 10.}}}},
      {"hDeltaPtVsPt", "Delta pt;p_{T} (GeV/#it{c});#Delta p_{T} / p_{T}", {HistType::kTH2D, {{200, 0., 10.}, {200, -1., 1.}}}},
    }};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (static_cast<o2::base::Propagator::MatCorrType>(materialCorrectionType.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
      o2::base::Propagator::Instance(true)->setMatLUT(lut);
    }
  }

  void process(aod::Collision const& collision, aod::McCollisions const& mcCollisions,
               aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades,
               aod::V0s const& v0s, TracksExt const& tracks, aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
  {
    const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      runNumber = bc.runNumber();
      auto timestamp = bc.timestamp();

      if (o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, timestamp)) {
        o2::base::Propagator::initFieldFromGRP(grpo);
      } else if (o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpMagPath, timestamp)) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
      } else {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpMagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << timestamp;
      }
    }

    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(bz);
    df2.setPropagateToPCA(propToDCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);

    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(materialCorrectionType.value);
    const auto primaryVertex = getPrimaryVertex(collision);
    o2::dataformats::DCA impactParameterTrk;
    for (const auto& trackedCascade : trackedCascades) {
      const auto trackCasc = trackedCascade.track_as<TracksExt>();
      auto trackParCovTrk = getTrackParCov(trackCasc);
      if (bzOnly) {
        o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovTrk, bz, 2.f, matCorr, &impactParameterTrk);
      } else {
        o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovTrk, 2.f, matCorr, &impactParameterTrk);
      }

      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExt>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.posTrack_as<TracksExt>();
      const auto& ntrack = v0.negTrack_as<TracksExt>();

      registry.fill(HIST("hDca"), std::sqrt(impactParameterTrk.getR2()));
      registry.fill(HIST("hDcaXY"), impactParameterTrk.getY());
      registry.fill(HIST("hDcaZ"), impactParameterTrk.getZ());
      registry.fill(HIST("hDcaVsPt"), impactParameterTrk.getY(), trackCasc.pt());
      registry.fill(HIST("hDcaVsR"), impactParameterTrk.getY(), RecoDecay::sqrtSumOfSquares(trackCasc.x(), trackCasc.y()));
      registry.fill(HIST("hMassVsPt"), trackedCascade.omegaMass(), trackCasc.pt());

      if (!ptrack.has_mcParticle() || !ntrack.has_mcParticle() || !bachelor.has_mcParticle()) {
        continue;
      }

      LOGF(debug, "ptrack (id: %d, pdg: %d) has mother %d", ptrack.mcParticleId(),
           ptrack.mcParticle().pdgCode(), ptrack.mcParticle().has_mothers() ? ptrack.mcParticle().mothersIds()[0] : -1);
      LOGF(debug, "ntrack (id: %d, pdg: %d) has mother %d", ntrack.mcParticleId(),
           ntrack.mcParticle().pdgCode(), ntrack.mcParticle().has_mothers() ? ntrack.mcParticle().mothersIds()[0] : -1);

      LOG(debug) << "bachelor with PDG code: " << bachelor.mcParticle().pdgCode();
      if (ptrack.mcParticle().has_mothers() && ntrack.mcParticle().has_mothers() &&
          ptrack.mcParticle().mothersIds()[0] == ntrack.mcParticle().mothersIds()[0]) {
        const auto v0part = ptrack.mcParticle().mothers_as<aod::McParticles>()[0];
        LOG(debug) << "v0 with PDG code: " << v0part.pdgCode();
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          const auto mother = v0part.mothers_as<aod::McParticles>()[0];
          const auto pdgCode = mother.pdgCode();
          LOG(debug) << "cascade with PDG code: " << pdgCode;
          if (std::abs(pdgCode) == kOmegaMinus) {
            LOG(debug) << "found Omega, looking for pions";
            std::array<double, 2> masses{RecoDecay::getMassPDG(kOmegaMinus), RecoDecay::getMassPDG(kPiPlus)};
            std::array<std::array<float, 3>, 2> momenta;
            std::array<double, 3> primaryVertexPos = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};
            const auto& mcColl = mother.mcCollision();
            std::array<double, 3> primaryVertexPosGen = {mcColl.posX(), mcColl.posY(), mcColl.posZ()};

            for (const auto& track : tracks) {
              if (!track.has_mcParticle()) {
                continue;
              }
              const auto mcpart = track.mcParticle();
              if (mcpart.pdgCode() == (pdgCode > 0 ? kPiPlus : -kPiPlus)) {
                LOGF(debug, "combining Omega with pion %d", track.globalIndex());
                auto trackParCovPion = getTrackParCov(track);
                o2::dataformats::DCA impactParameterPion;
                if (bzOnly) {
                  o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovPion, bz, 2.f, matCorr, &impactParameterPion);
                } else {
                  o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovPion, 2.f, matCorr, &impactParameterPion);
                }

                trackParCovTrk.getPxPyPzGlo(momenta[0]);
                trackParCovPion.getPxPyPzGlo(momenta[1]);
                registry.fill(HIST("hDeltaPtVsPt"), mcpart.pt(), (trackParCovPion.getPt() - mcpart.pt()) / mcpart.pt());
                registry.fill(HIST("hMassOmegac"), RecoDecay::m(momenta, masses));
                registry.fill(HIST("hMassOmegacVsPt"), RecoDecay::m(momenta, masses), RecoDecay::pt(momenta[0], momenta[1]));

                if (df2.process(trackParCovTrk, trackParCovPion)) {
                  const auto& secondaryVertex = df2.getPCACandidate();
                  const auto decayLength = RecoDecay::distance(secondaryVertex, primaryVertexPos);
                  if (std::abs(RecoDecay::m(momenta, masses) - RecoDecay::getMassPDG(analysis::pdg::Code::kOmegaC0)) < 0.02) {
                    registry.fill(HIST("hDecayLength"), decayLength * 1e4);
                    registry.fill(HIST("hDecayLengthScaled"), decayLength * RecoDecay::getMassPDG(analysis::pdg::Code::kOmegaC0) / RecoDecay::p(momenta[0], momenta[1]) * 1e4);
                  }
                  primaryVertexPosGen;
                  if (mother.has_mothers()) {
                    const auto& cand = mother.template mothers_first_as<aod::McParticles>();
                    if (std::abs(cand.pdgCode()) == analysis::pdg::Code::kOmegaC0 && mcpart.has_mothers()) {
                      if (mcpart.mothersIds()[0] == cand.globalIndex()) {
                        std::array<double, 3> secondaryVertexGen = {mother.vx(), mother.vy(), mother.vz()};
                        const auto decayLengthGen = RecoDecay::distance(secondaryVertexGen, primaryVertexPosGen);
                        registry.fill(HIST("hDecayLengthId"), decayLength * 1e4);
                        registry.fill(HIST("hDecayLengthScaledId"), decayLength * RecoDecay::getMassPDG(analysis::pdg::Code::kOmegaC0) / RecoDecay::p(momenta[0], momenta[1]) * 1e4);
                        registry.fill(HIST("hDecayLengthScaledGen"), decayLengthGen * RecoDecay::getMassPDG(analysis::pdg::Code::kOmegaC0) / RecoDecay::p(momenta[0], momenta[1]) * 1e4);
                        registry.fill(HIST("hDecayLengthGen"), decayLengthGen * 1e4);
                        registry.fill(HIST("hDeltaDecayLength"), (decayLength - decayLengthGen) * 1e4);
                      }
                    }
                  }
                }

                // MC-based mass
                float smear = 1.;
                momenta[0] = {smear * mother.px(), smear * mother.py(), smear * mother.pz()};
                // smear = gRandom->Gaus(1., 0.02);
                momenta[1] = {smear * mcpart.px(), smear * mcpart.py(), smear * mcpart.pz()};
                registry.fill(HIST("hMassOmegacGen"), RecoDecay::m(momenta, masses));
              }
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
    adaptAnalysisTask<HfTaskOmegacSt>(cfgc)};
}
