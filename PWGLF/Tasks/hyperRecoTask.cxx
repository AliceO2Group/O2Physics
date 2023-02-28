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
// Build hypertriton candidates from V0s and tracks
// =====================
//
//
//
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"He3"};
float kHyperPDG = 1010010030;
} // namespace

namespace o2::aod
{
/// FemtoDreamCollision
namespace hyperrec
{
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);      //! bool: true for matter
DECLARE_SOA_COLUMN(Px, px, float);                 //! Momentum of the candidate (x direction)
DECLARE_SOA_COLUMN(Py, py, float);                 //! Momentum of the candidate (y direction)
DECLARE_SOA_COLUMN(Pz, pz, float);                 //! Momentum of the candidate (z direction)
DECLARE_SOA_COLUMN(XDecVtx, xDecVtx, float);       //! Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YDecVtx, yDecVtx, float);       //! Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZDecVtx, zDecVtx, float);       //! Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(DcaV0Daug, dcaV0Daug, float);   //! DCA between daughters
DECLARE_SOA_COLUMN(CosPA, cosPA, double);          //! Cosine of the pointing angle
DECLARE_SOA_COLUMN(NSigmaHe, nSigmaHe, float);     //! Number of sigmas of the He daughter
DECLARE_SOA_COLUMN(NTPCclusHe, nTPCclusHe, int);   //! Number of TPC clusters of the He daughter
DECLARE_SOA_COLUMN(DcaHe, dcaHe, float);           //! DCA between He daughter and V0
DECLARE_SOA_COLUMN(DcaPi, dcaPi, float);           //! DCA between pi daughter and V0
DECLARE_SOA_COLUMN(GenPx, genPx, float);           //! Momentum of the candidate (x direction)
DECLARE_SOA_COLUMN(GenPy, genPy, float);           //! Momentum of the candidate (y direction)
DECLARE_SOA_COLUMN(GenPz, genPz, float);           //! Momentum of the candidate (z direction)
DECLARE_SOA_COLUMN(GenXDecVtx, genXDecVtx, float); //! Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(GenYDecVtx, genYDecVtx, float); //! Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(GenZDecVtx, genZDecVtx, float); //! Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(IsReco, isReco, bool);          //! bool: true for reco
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);      //! bool: true for signal
} // namespace hyperrec

DECLARE_SOA_TABLE(DataHypCands, "AOD", "DATAHYPCANDS",
                  o2::soa::Index<>,
                  o2::aod::hyperrec::Px,
                  o2::aod::hyperrec::Py,
                  o2::aod::hyperrec::Pz,
                  o2::aod::hyperrec::XDecVtx,
                  o2::aod::hyperrec::YDecVtx,
                  o2::aod::hyperrec::ZDecVtx,
                  o2::aod::hyperrec::DcaV0Daug,
                  o2::aod::hyperrec::CosPA,
                  o2::aod::hyperrec::NSigmaHe,
                  o2::aod::hyperrec::NTPCclusHe,
                  o2::aod::hyperrec::DcaHe,
                  o2::aod::hyperrec::DcaPi);

DECLARE_SOA_TABLE(MCHypCands, "AOD", "MCHYPCANDS",
                  o2::soa::Index<>,
                  o2::aod::hyperrec::Px,
                  o2::aod::hyperrec::Py,
                  o2::aod::hyperrec::Pz,
                  o2::aod::hyperrec::XDecVtx,
                  o2::aod::hyperrec::YDecVtx,
                  o2::aod::hyperrec::ZDecVtx,
                  o2::aod::hyperrec::DcaV0Daug,
                  o2::aod::hyperrec::CosPA,
                  o2::aod::hyperrec::NSigmaHe,
                  o2::aod::hyperrec::NTPCclusHe,
                  o2::aod::hyperrec::DcaHe,
                  o2::aod::hyperrec::DcaPi,
                  o2::aod::hyperrec::GenPx,
                  o2::aod::hyperrec::GenPy,
                  o2::aod::hyperrec::GenPz,
                  o2::aod::hyperrec::GenXDecVtx,
                  o2::aod::hyperrec::GenYDecVtx,
                  o2::aod::hyperrec::GenZDecVtx,
                  o2::aod::hyperrec::IsReco,
                  o2::aod::hyperrec::IsSignal);

using DataHypCand = o2::aod::DataHypCands::iterator;
using MCHypCand = o2::aod::MCHypCands::iterator;
} // namespace o2::aod

struct hyperCandidate {
  int posTrackID; // TODO: check whether int is enough
  int negTrackID;
  bool isMatter;
  std::array<float, 3> mom;
  std::array<float, 3> decVtx;
  float dcaV0dau = -1;
  float cosPA;
  float rapidity;
  float nSigmaHe3;
  float nTPCClustersHe3;
  float he3DCAXY;
  float piDCAXY;

  std::array<float, 3> gMom;
  std::array<float, 3> gDecVtx;
  float gRapidity;

  bool isSignal = false; // true MC signal
  bool isReco = false;   // true if the candidate is actually reconstructed
};

struct hyperRecoTask {

  Produces<aod::DataHypCands> outputDataTable;
  Produces<aod::MCHypCands> outputMCTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Selection criteria
  Configurable<double> v0cospa{"hypcospa", 0.95, "V0 CosPA"}; // very open, please
  Configurable<float> dcav0dau{"hypdcaDau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dca3hetopv{"hedcatopv", -1, "DCA helium To PV"};
  Configurable<float> dcapiontopv{"pidcatopv", .1, "DCA pion To PV"};
  Configurable<float> v0radius{"hypradius", 1.0, "hyp radius"};
  Configurable<float> etaMax{"eta", 0.8, "eta daughter"};
  Configurable<float> heliumNsigmaMax{"heliumNsigmaMax", 8, "helium dEdx cut (n sigma)"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;

  // bethe bloch parameters
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};

  o2::track::TrackPar lHyTrack;

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<int> heDauPdg{"heDauDPG", 1010020030, "PDG code of the helium"};

  // std vector of candidates
  std::vector<hyperCandidate> hyperCandidates;

  void init(InitContext const&)
  {

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    // fitter.setWeightedFinalPCA(d_UseWeightedPCA);
    int mat{static_cast<int>(cfgMaterialCorrection)};
    fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));
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
    fitter.setBz(d_bz);
    mRunNumber = bc.runNumber();
  }

  int mRunNumber;
  float d_bz;

  void fillCandidateData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0s const& V0s, TracksFull const& tracks)
  {
    for (auto& v0 : V0s) {

      LOG(info) << "collision: " << collision.globalIndex();
      LOG(info) << "V0 collision ID: " << v0.collisionId();
      auto posTrack = v0.posTrack_as<TracksFull>();
      auto negTrack = v0.negTrack_as<TracksFull>();
      if (abs(posTrack.eta()) > etaMax || abs(negTrack.eta()) > etaMax)
        continue;
      LOG(info) << "pos track: "
                << " pt: " << posTrack.pt();

      LOG(info) << "neg track: "
                << " pt: " << negTrack.pt();

      double expBethePos{tpc::BetheBlochAleph(static_cast<double>(posTrack.tpcInnerParam() * 2 / constants::physics::MassHelium3), cfgBetheBlochParams->get("He3", 0u), cfgBetheBlochParams->get("He3", 1u), cfgBetheBlochParams->get("He3", 2u), cfgBetheBlochParams->get("He3", 3u), cfgBetheBlochParams->get("He3", 4u))};

      double expBetheNeg{tpc::BetheBlochAleph(static_cast<double>(negTrack.tpcInnerParam() * 2 / constants::physics::MassHelium3), cfgBetheBlochParams->get("He3", 0u), cfgBetheBlochParams->get("He3", 1u), cfgBetheBlochParams->get("He3", 2u), cfgBetheBlochParams->get("He3", 3u), cfgBetheBlochParams->get("He3", 4u))};

      double expSigmaPos{expBethePos * cfgBetheBlochParams->get("He3", 5u)};
      double expSigmaNeg{expBetheNeg * cfgBetheBlochParams->get("He3", 5u)};
      auto nSigmaTPCpos = static_cast<float>((posTrack.tpcSignal() - expBethePos) / expSigmaPos);
      auto nSigmaTPCneg = static_cast<float>((negTrack.tpcSignal() - expBetheNeg) / expSigmaNeg);
      if (abs(nSigmaTPCpos) > heliumNsigmaMax || abs(nSigmaTPCneg) > heliumNsigmaMax)
        continue;

      hyperCandidate hypCand;
      hypCand.isMatter = abs(nSigmaTPCpos) < abs(nSigmaTPCneg);

      hypCand.nSigmaHe3 = hypCand.isMatter ? nSigmaTPCpos : nSigmaTPCneg;
      hypCand.nTPCClustersHe3 = hypCand.isMatter ? posTrack.tpcNClsFindable() : negTrack.tpcNClsFindable();

      auto posTrackCov = getTrackParCov(posTrack);
      auto negTrackCov = getTrackParCov(negTrack);

      int nCand = 0;
      try {
        nCand = fitter.process(posTrackCov, negTrackCov);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        continue;
      }
      if (nCand == 0) {
        continue;
      }

      auto& propPosTrack = fitter.getTrack(0);
      auto& propNegTrack = fitter.getTrack(1);

      std::array<float, 3> posTrackP;
      std::array<float, 3> negTrackP;

      propPosTrack.getPxPyPzGlo(posTrackP);
      propNegTrack.getPxPyPzGlo(negTrackP);

      unsigned int posAbsCharge = hypCand.isMatter ? 2 : 1;
      unsigned int negAbsCharge = !hypCand.isMatter ? 2 : 1;

      posTrackP[0] *= posAbsCharge, posTrackP[1] *= posAbsCharge, posTrackP[2] *= posAbsCharge;
      negTrackP[0] *= negAbsCharge, negTrackP[1] *= negAbsCharge, negTrackP[2] *= negAbsCharge;

      auto posPrimVtx = array{collision.posX(), collision.posY(), collision.posZ()};

      // get decay vertex coordinates
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        hypCand.decVtx[i] = vtx[i] - posPrimVtx[i];
        hypCand.mom[i] = posTrackP[i] + negTrackP[i];
      }

      hypCand.dcaV0dau = TMath::Sqrt(fitter.getChi2AtPCACandidate());

      // Apply selections so a skimmed table is created only
      if (hypCand.dcaV0dau > dcav0dau) {
        continue;
      }

      hypCand.cosPA = RecoDecay::cpa(posPrimVtx, array{hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2]}, array{hypCand.mom[0], hypCand.mom[1], hypCand.mom[2]});
      if (hypCand.cosPA < v0cospa) {
        continue;
      }

      // if survived all selections, propagate decay daughters to PV

      gpu::gpustd::array<float, 2> dcaInfo;

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, posTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      hypCand.isMatter ? hypCand.he3DCAXY = dcaInfo[0] : hypCand.piDCAXY = dcaInfo[0];

      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, negTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      hypCand.isMatter ? hypCand.piDCAXY = dcaInfo[0] : hypCand.he3DCAXY = dcaInfo[0];

      // finally, push back the candidate
      hypCand.isReco = true;
      hypCand.posTrackID = posTrack.globalIndex();
      hypCand.negTrackID = negTrack.globalIndex();

      hyperCandidates.push_back(hypCand);
    }
  }

  void fillMCinfo(aod::McTrackLabels const& trackLabels, aod::McParticles const& particlesMC)
  {
    std::vector<unsigned int> filledMothers;
    for (auto& hypCand : hyperCandidates) {

      auto mcLabPos = trackLabels.rawIteratorAt(hypCand.posTrackID);
      auto mcLabNeg = trackLabels.rawIteratorAt(hypCand.negTrackID);

      if (mcLabPos.has_mcParticle() && mcLabNeg.has_mcParticle()) {
        auto mcTrackPos = mcLabPos.mcParticle_as<aod::McParticles>();
        auto mcTrackNeg = mcLabNeg.mcParticle_as<aod::McParticles>();
        if (mcTrackPos.has_mothers() && mcTrackNeg.has_mothers()) {
          for (auto& negMother : mcTrackNeg.mothers_as<aod::McParticles>()) {
            for (auto& posMother : mcTrackPos.mothers_as<aod::McParticles>()) {
              if (posMother.globalIndex() != negMother.globalIndex())
                continue;
              if (!((mcTrackPos.pdgCode() == heDauPdg && mcTrackNeg.pdgCode() == -211) || (mcTrackPos.pdgCode() == 211 && mcTrackNeg.pdgCode() == -1 * heDauPdg))) // TODO: check warning for constant comparison
                continue;
              if (abs(posMother.pdgCode()) != kHyperPDG)
                continue;
              auto posPrimVtx = array{posMother.vx(), posMother.vy(), posMother.vz()};
              auto secVtx = array{mcTrackPos.vx(), mcTrackPos.vy(), mcTrackPos.vz()};
              auto posMom = array{posMother.px(), posMother.py(), posMother.pz()};
              for (int i = 0; i < 3; i++) {
                hypCand.gDecVtx[i] = secVtx[i] - posPrimVtx[i];
                hypCand.gMom[i] = posMom[i];
              }
              hypCand.isSignal = true;
              filledMothers.push_back(posMother.globalIndex());
            }
          }
        }
      }
    }
    for (auto& mcPart : particlesMC) {
      if (abs(mcPart.pdgCode()) != kHyperPDG)
        continue;
      std::array<float, 3> secVtx;
      std::array<float, 3> primVtx = {mcPart.vx(), mcPart.vy(), mcPart.vz()};
      std::array<float, 3> momMother = {mcPart.px(), mcPart.py(), mcPart.pz()};
      bool isHeFound = false;
      for (auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
        if (abs(mcDaught.pdgCode()) == heDauPdg) {
          secVtx = {mcDaught.vx(), mcDaught.vy(), mcDaught.vz()};
          isHeFound = true;
          break;
        }
      }
      if (!isHeFound)
        continue;
      if (std::find(filledMothers.begin(), filledMothers.end(), mcPart.globalIndex()) != std::end(filledMothers))
        continue;
      hyperCandidate hypCand;
      for (int i = 0; i < 3; i++) {
        hypCand.gDecVtx[i] = secVtx[i] - primVtx[i];
        hypCand.gMom[i] = momMother[i];
      }
      hypCand.posTrackID = -1;
      hypCand.negTrackID = -1;
      hypCand.isSignal = true;
    }
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!collision.sel8())
      return;

    if (abs(collision.posZ()) > 10.f) {
      return;

      fillCandidateData(collision, V0s, tracks);
    }

    for (auto& hypCand : hyperCandidates) {
      outputDataTable(hypCand.mom[0], hypCand.mom[1], hypCand.mom[2],
                      hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                      hypCand.dcaV0dau, hypCand.cosPA, hypCand.nSigmaHe3,
                      hypCand.nTPCClustersHe3, hypCand.he3DCAXY, hypCand.piDCAXY);
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processData, "Data analysis", true);

  void processMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0s const& V0s, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!collision.sel8())
      return;

    if (abs(collision.posZ()) > 10.f) {
      return;
      fillCandidateData(collision, V0s, tracks);
      fillMCinfo(trackLabelsMC, particlesMC);
    }

    for (auto& hypCand : hyperCandidates) {
      outputMCTable(hypCand.mom[0], hypCand.mom[1], hypCand.mom[2],
                    hypCand.decVtx[0], hypCand.decVtx[1], hypCand.decVtx[2],
                    hypCand.dcaV0dau, hypCand.cosPA, hypCand.nSigmaHe3,
                    hypCand.nTPCClustersHe3, hypCand.he3DCAXY, hypCand.piDCAXY,
                    hypCand.gMom[0], hypCand.gMom[1], hypCand.gMom[2],
                    hypCand.gDecVtx[0], hypCand.gDecVtx[1], hypCand.gDecVtx[2], hypCand.isReco, hypCand.isSignal);
    }
  }
  PROCESS_SWITCH(hyperRecoTask, processMC, "MC analysis", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hyperRecoTask>(cfgc)};
}
