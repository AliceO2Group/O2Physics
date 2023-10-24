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

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/PDG.h"
#include "PWGHF/Core/SelectorCuts.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace stomegac
{

/// definition of tables
namespace hfsttraining
{
DECLARE_SOA_COLUMN(InvMassD0, invMassD0, float);                 //!
DECLARE_SOA_COLUMN(InvMassD0bar, invMassD0bar, float);           //!
DECLARE_SOA_COLUMN(InvMassDplus, invMassDplus, float);           //!
DECLARE_SOA_COLUMN(InvMassDsToKKPi, invMassDsToKKPi, float);     //!
DECLARE_SOA_COLUMN(InvMassDsToPiKK, invMassDsToPiKK, float);     //!
DECLARE_SOA_COLUMN(InvMassLcToPKPi, invMassLcToPKPi, float);     //!
DECLARE_SOA_COLUMN(InvMassLcToPiKP, invMassLcToPiKP, float);     //!
DECLARE_SOA_COLUMN(InvMassXicToPKPi, invMassXicToPKPi, float);   //!
DECLARE_SOA_COLUMN(InvMassXicToPiKP, invMassXicToPiKP, float);   //!
DECLARE_SOA_COLUMN(PT2Prong, pT2Prong, float);                   //!
DECLARE_SOA_COLUMN(PT3Prong, pT3Prong, float);                   //!
DECLARE_SOA_COLUMN(DeltaMassKKFirst, deltaMassKKFirst, float);   //!
DECLARE_SOA_COLUMN(DeltaMassKKSecond, deltaMassKKSecond, float); //!
DECLARE_SOA_COLUMN(PT1, pT1, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY1, dcaPrimXY1, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ1, dcaPrimZ1, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC1, nsigmaPiTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC1, nsigmaKaTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC1, nsigmaPrTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF1, nsigmaPiTOF1, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF1, nsigmaKaTOF1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF1, nsigmaPrTOF1, float);           //!
DECLARE_SOA_COLUMN(PT2, pT2, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY2, dcaPrimXY2, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ2, dcaPrimZ2, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC2, nsigmaPiTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC2, nsigmaKaTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC2, nsigmaPrTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF2, nsigmaPiTOF2, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF2, nsigmaKaTOF2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF2, nsigmaPrTOF2, float);           //!
DECLARE_SOA_COLUMN(FlagOrigin, flagOrigin, int8_t);              //!
DECLARE_SOA_COLUMN(Channel, channel, int8_t);                    //!
DECLARE_SOA_COLUMN(HFSelBit, hfselbit, int8_t);                  //!
DECLARE_SOA_COLUMN(IsInCorrectColl, isInCorrectColl, bool);      //!
} // namespace hfsttraining
}

DECLARE_SOA_TABLE(HfOmegacSt, "AOD", "HFOMEGACST",
                  stomegac::hfsttraining::InvMassD0,
                  stomegac::hfsttraining::InvMassD0bar,
                  stomegac::hfsttraining::PT2Prong,
                  stomegac::hfsttraining::PT1,
                  stomegac::hfsttraining::DCAPrimXY1,
                  stomegac::hfsttraining::DCAPrimZ1,
                  stomegac::hfsttraining::NsigmaPiTPC1,
                  stomegac::hfsttraining::NsigmaKaTPC1,
                  stomegac::hfsttraining::NsigmaPiTOF1,
                  stomegac::hfsttraining::NsigmaKaTOF1,
                  stomegac::hfsttraining::PT2,
                  stomegac::hfsttraining::DCAPrimXY2,
                  stomegac::hfsttraining::DCAPrimZ2,
                  stomegac::hfsttraining::NsigmaPiTPC2,
                  stomegac::hfsttraining::NsigmaKaTPC2,
                  stomegac::hfsttraining::NsigmaPiTOF2,
                  stomegac::hfsttraining::NsigmaKaTOF2,
                  stomegac::hfsttraining::FlagOrigin,
                  stomegac::hfsttraining::IsInCorrectColl);
}

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

  Produces<aod::HfOmegacSt> outputTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> df2;

  bool bzOnly = true;
  int runNumber{0};

  using TracksExt = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
  using TracksExtMc = soa::Join<TracksExt, aod::McTrackLabels>;

  HistogramRegistry registry{
    "registry",
    {
      {"hDca", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"hDcaXY", "DCA;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"hDcaXYVsPt", "DCA;p_{T} (GeV/#it{c};DCA_{xy} (cm)", {HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}}}},
      {"hDcaZ", "DCA;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"hDcaZVsPt", "DCA;p_{T} (GeV/#it{c});DCA_{z} (cm)", {HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}}}},
      {"hDcaVsPt", "DCA;DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}}}},
      {"hDcaVsR", "DCA;DCA (cm);R (cm)", {HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}}}},
      {"hDecayLength", "Decay length;L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthId", "Decay length (true #Omega_{c});L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthGen", "Decay length (gen);L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDeltaDecayLength", "#Delta decay length (gen);#Delta L (#mum)", {HistType::kTH1D, {{200, -250., 250.}}}},
      {"hDecayLengthScaled", "Decay length * M/p;L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthScaledId", "Decay length * M/p (true #Omega_{c});L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthScaledGen", "Decay length * M/p (MC id);L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthScaledMc", "Decay length * M/p (MC);L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hMassOmegac", "inv. mass #Omega + #pi;inv. mass (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 1.5, 3.}}}},
      {"hMassOmegacVsPt", "inv. mass #Omega + #pi;inv. mass (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", {HistType::kTH2D, {{400, 1.5, 3.}, {10, 0., 10.}}}},
      {"hMassOmegacId", "inv. mass #Omega + #pi (MC ID);inv. mass (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 1.5, 3.}}}},
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

    df2.setBz(bz);
    df2.setPropagateToPCA(propToDCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
  }

  void processMc(aod::McCollision const& mcCollision,
                 aod::McParticles const& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() != kOmegaMinus) {
        continue;
      }
      if (mcParticle.has_mothers() && mcParticle.mothers_first_as<aod::McParticles>().pdgCode() == analysis::pdg::Code::kOmegaC0) {
        const auto& mcColl = mcParticle.mcCollision();
        std::array<double, 3> primaryVertexPosGen = {mcColl.posX(), mcColl.posY(), mcColl.posZ()};
        std::array<double, 3> secondaryVertexGen = {mcParticle.vx(), mcParticle.vy(), mcParticle.vz()};
        const auto decayLengthGen = RecoDecay::distance(secondaryVertexGen, primaryVertexPosGen);
        registry.fill(HIST("hDecayLengthScaledMc"), decayLengthGen * o2::analysis::pdg::MassOmegaC0 / mcParticle.mothers_first_as<aod::McParticles>().p() * 1e4);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegacSt, processMc, "Process MC", true);

  void processData(aod::Collision const& collision,
                   aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades,
                   aod::V0s const& v0s, TracksExt const& tracks, aod::BCsWithTimestamps const&)
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
      const auto& v0TrackPos = v0.posTrack_as<TracksExt>();
      const auto& v0TrackNeg = v0.negTrack_as<TracksExt>();

      std::array<double, 3> m{o2::analysis::pdg::MassProton, o2::analysis::pdg::MassPiMinus, o2::analysis::pdg::MassKMinus};
      std::array<std::array<float, 3>, 3> p;
      p[0] = {v0TrackPos.px(), v0TrackPos.py(), v0TrackPos.pz()};
      p[1] = {v0TrackNeg.px(), v0TrackNeg.py(), v0TrackNeg.pz()};
      p[2] = {bachelor.px(), bachelor.py(), bachelor.pz()};
      double massOmega1 = RecoDecay::m(p, m);
      p[0] = {v0TrackNeg.px(), v0TrackNeg.py(), v0TrackNeg.pz()};
      p[1] = {v0TrackPos.px(), v0TrackPos.py(), v0TrackPos.pz()};
      p[2] = {bachelor.px(), bachelor.py(), bachelor.pz()};
      double massOmega2 = RecoDecay::m(p, m);

      registry.fill(HIST("hDca"), std::sqrt(impactParameterTrk.getR2()));
      registry.fill(HIST("hDcaXY"), impactParameterTrk.getY());
      registry.fill(HIST("hDcaXYVsPt"), trackParCovTrk.getPt(), impactParameterTrk.getY());
      registry.fill(HIST("hDcaZ"), impactParameterTrk.getZ());
      registry.fill(HIST("hDcaZVsPt"), trackParCovTrk.getPt(), impactParameterTrk.getZ());
      registry.fill(HIST("hDcaVsPt"), impactParameterTrk.getY(), trackCasc.pt());
      registry.fill(HIST("hDcaVsR"), impactParameterTrk.getY(), RecoDecay::sqrtSumOfSquares(trackCasc.x(), trackCasc.y()));
      registry.fill(HIST("hMassVsPt"), massOmega1, trackCasc.pt());
      registry.fill(HIST("hMassVsPt"), massOmega2, trackCasc.pt());

      if ((std::abs(massOmega1 - o2::analysis::pdg::MassOmegaMinus) < .1) ||
          (std::abs(massOmega2 - o2::analysis::pdg::MassOmegaMinus) < .1)) {
        LOGF(debug, "found candidate in mass range");
        if ((std::abs(bachelor.tpcNSigmaKa()) < 3.) &&
            (((std::abs(v0TrackPos.tpcNSigmaPr()) < 3.) && (std::abs(v0TrackNeg.tpcNSigmaPi()) < 3.)) ||
             ((std::abs(v0TrackPos.tpcNSigmaPi()) < 3.) && (std::abs(v0TrackNeg.tpcNSigmaPr()) < 3.)))) {
          LOGF(debug, ".. species compatible with Omega");
          std::array<double, 2> masses{o2::analysis::pdg::MassOmegaMinus, o2::analysis::pdg::MassPiPlus};
          std::array<std::array<float, 3>, 2> momenta;
          std::array<double, 3> primaryVertexPos = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};

          for (const auto& track : tracks) {
            if (std::abs(track.tpcNSigmaPi()) < 3.) {
              LOGF(debug, "  .. combining with pion candidate %d", track.globalIndex());
              auto trackParCovPion = getTrackParCov(track);
              o2::dataformats::DCA impactParameterPion;
              if (bzOnly) {
                o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovPion, bz, 2.f, matCorr, &impactParameterPion);
              } else {
                o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovPion, 2.f, matCorr, &impactParameterPion);
              }

              trackParCovTrk.getPxPyPzGlo(momenta[0]);
              trackParCovPion.getPxPyPzGlo(momenta[1]);
              registry.fill(HIST("hMassOmegac"), RecoDecay::m(momenta, masses));
              registry.fill(HIST("hMassOmegacVsPt"), RecoDecay::m(momenta, masses), RecoDecay::pt(momenta[0], momenta[1]));

              if (df2.process(trackParCovTrk, trackParCovPion)) {
                const auto& secondaryVertex = df2.getPCACandidate();
                const auto decayLength = RecoDecay::distance(secondaryVertex, primaryVertexPos);
                if (std::abs(RecoDecay::m(momenta, masses) - o2::analysis::pdg::MassOmegaC0) < 0.02) {
                  registry.fill(HIST("hDecayLength"), decayLength * 1e4);
                  registry.fill(HIST("hDecayLengthScaled"), decayLength * o2::analysis::pdg::MassOmegaC0 / RecoDecay::p(momenta[0], momenta[1]) * 1e4);
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegacSt, processData, "Process data", true);

  void processGen(aod::Collision const& collision, aod::McCollisions const& mcCollisions,
                  aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades,
                  aod::V0s const& v0s, TracksExtMc const& tracks, aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
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

    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(materialCorrectionType.value);
    const auto primaryVertex = getPrimaryVertex(collision);
    o2::dataformats::DCA impactParameterTrk;
    for (const auto& trackedCascade : trackedCascades) {
      const auto trackCasc = trackedCascade.track_as<TracksExtMc>();
      auto trackParCovTrk = getTrackParCov(trackCasc);
      if (bzOnly) {
        o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovTrk, bz, 2.f, matCorr, &impactParameterTrk);
      } else {
        o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovTrk, 2.f, matCorr, &impactParameterTrk);
      }

      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.bachelor_as<TracksExtMc>();
      const auto& v0 = casc.v0();
      const auto& v0TrackPos = v0.posTrack_as<TracksExtMc>();
      const auto& v0TrackNeg = v0.negTrack_as<TracksExtMc>();

      if (!v0TrackPos.has_mcParticle() || !v0TrackNeg.has_mcParticle() || !bachelor.has_mcParticle()) {
        continue;
      }

      LOGF(debug, "v0TrackPos (id: %d, pdg: %d) has mother %d", v0TrackPos.mcParticleId(),
           v0TrackPos.mcParticle().pdgCode(), v0TrackPos.mcParticle().has_mothers() ? v0TrackPos.mcParticle().mothersIds()[0] : -1);
      LOGF(debug, "v0TrackNeg (id: %d, pdg: %d) has mother %d", v0TrackNeg.mcParticleId(),
           v0TrackNeg.mcParticle().pdgCode(), v0TrackNeg.mcParticle().has_mothers() ? v0TrackNeg.mcParticle().mothersIds()[0] : -1);

      LOG(debug) << "bachelor with PDG code: " << bachelor.mcParticle().pdgCode();
      if (v0TrackPos.mcParticle().has_mothers() && v0TrackNeg.mcParticle().has_mothers() &&
          v0TrackPos.mcParticle().mothersIds()[0] == v0TrackNeg.mcParticle().mothersIds()[0]) {
        const auto v0part = v0TrackPos.mcParticle().mothers_first_as<aod::McParticles>();
        LOG(debug) << "v0 with PDG code: " << v0part.pdgCode();
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          const auto mother = v0part.mothers_as<aod::McParticles>()[0];
          const auto pdgCode = mother.pdgCode();
          LOG(debug) << "cascade with PDG code: " << pdgCode;
          if (std::abs(pdgCode) == kOmegaMinus) {
            LOG(debug) << "found Omega, looking for pions";
            std::array<double, 2> masses{o2::analysis::pdg::MassOmegaMinus, o2::analysis::pdg::MassPiPlus};
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
                registry.fill(HIST("hMassOmegacId"), RecoDecay::m(momenta, masses));

                if (df2.process(trackParCovTrk, trackParCovPion)) {
                  const auto& secondaryVertex = df2.getPCACandidate();
                  const auto decayLength = RecoDecay::distance(secondaryVertex, primaryVertexPos);
                  if (mother.has_mothers()) {
                    const auto& cand = mother.template mothers_first_as<aod::McParticles>();
                    if (std::abs(cand.pdgCode()) == analysis::pdg::Code::kOmegaC0 && mcpart.has_mothers()) {
                      if (mcpart.mothersIds()[0] == cand.globalIndex()) {
                        registry.fill(HIST("hDecayLengthId"), decayLength * 1e4);
                        registry.fill(HIST("hDecayLengthScaledId"), decayLength * o2::analysis::pdg::MassOmegaC0 / RecoDecay::p(momenta[0], momenta[1]) * 1e4);

                        std::array<double, 3> secondaryVertexGen = {mother.vx(), mother.vy(), mother.vz()};
                        const auto decayLengthGen = RecoDecay::distance(secondaryVertexGen, primaryVertexPosGen);
                        registry.fill(HIST("hDecayLengthGen"), decayLengthGen * 1e4);
                        registry.fill(HIST("hDecayLengthScaledGen"), decayLengthGen * o2::analysis::pdg::MassOmegaC0 / RecoDecay::p(momenta[0], momenta[1]) * 1e4);

                        registry.fill(HIST("hDeltaDecayLength"), (decayLength - decayLengthGen) * 1e4);
                      }
                    }
                  }
                }

                // MC-based mass
                momenta[0] = {mother.px(), mother.py(), mother.pz()};
                momenta[1] = {mcpart.px(), mcpart.py(), mcpart.pz()};
                registry.fill(HIST("hMassOmegacGen"), RecoDecay::m(momenta, masses));
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegacSt, processGen, "Process using MC information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskOmegacSt>(cfgc)};
}
