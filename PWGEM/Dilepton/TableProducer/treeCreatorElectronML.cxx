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
// ========================
//
// This code will create data table for inputs to machine learning for electrons.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector4D.h"

#include <map>
#include <random>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;

using FullTrackExt = FullTracksExt::iterator;

using FullTracksExtMC = soa::Join<FullTracksExt, aod::McTrackLabels>;
using FullTrackExtMC = FullTracksExtMC::iterator;

namespace o2::aod
{

namespace mycollision // reconstructed collision information
{
DECLARE_SOA_COLUMN(Bz, bz, float);         //! in unit of kG.
DECLARE_SOA_COLUMN(MCPosX, mcposX, float); //!
DECLARE_SOA_COLUMN(MCPosY, mcposY, float); //!
DECLARE_SOA_COLUMN(MCPosZ, mcposZ, float); //!
} // namespace mycollision
DECLARE_SOA_TABLE(MyCollisions, "AOD", "MYCOLLISION", //! vertex information of collision
                  o2::soa::Index<>, bc::GlobalBC, bc::RunNumber, collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, evsel::Sel8, mycollision::Bz,
                  mccollision::GeneratorsID, mycollision::MCPosX, mycollision::MCPosY, mycollision::MCPosZ, mult::MultNTracksPV,
                  cent::CentFT0M, cent::CentFT0A, cent::CentFT0C);
using MyCollision = MyCollisions::iterator;

namespace mytrack
{
DECLARE_SOA_INDEX_COLUMN(MyCollision, mycollision);                   //!
DECLARE_SOA_COLUMN(Sign, sign, int);                                  //!
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int);                  //!
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int);      //!
DECLARE_SOA_COLUMN(DCAresXY, dcaresXY, float);                        //!
DECLARE_SOA_COLUMN(DCAresZ, dcaresZ, float);                          //!
DECLARE_SOA_COLUMN(MCVx, mcvx, float);                                //!
DECLARE_SOA_COLUMN(MCVy, mcvy, float);                                //!
DECLARE_SOA_COLUMN(MCVz, mcvz, float);                                //!
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);       //!
DECLARE_SOA_COLUMN(MotherIds, motherIds, std::vector<int>);           //! eta values of the matched tracks
DECLARE_SOA_COLUMN(MotherPdgCodes, motherpdgCodes, std::vector<int>); //! eta values of the matched tracks
} // namespace mytrack

// reconstructed track information
DECLARE_SOA_TABLE(MyTracks, "AOD", "MYTRACK", //!
                  o2::soa::Index<>, mytrack::MyCollisionId, mytrack::Sign,
                  track::Pt, track::Eta, track::Phi, track::Tgl,
                  track::DcaXY, track::DcaZ, mytrack::DCAresXY, mytrack::DCAresZ, track::CZY,
                  track::TPCNClsFindable, mytrack::TPCNClsFound, mytrack::TPCNClsCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, /*pidtpc::TPCNSigmaMu,*/ pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, /*pidtof::TOFNSigmaMu,*/ pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::TOFChi2, track::ITSChi2NCl, track::ITSClusterSizes,
                  mytrack::MCVx, mytrack::MCVy, mytrack::MCVz,
                  mcparticle::PdgCode, mytrack::IsPhysicalPrimary, mytrack::MotherIds, mytrack::MotherPdgCodes);

// iterators
using MyTrack = MyTracks::iterator;

namespace mypair
{
DECLARE_SOA_INDEX_COLUMN(MyCollision, mycollision);                       //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, MyTracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, MyTracks, "_Neg"); //!

DECLARE_SOA_COLUMN(Mass, mass, float);           //!
DECLARE_SOA_COLUMN(Pt, pt, float);               //!
DECLARE_SOA_COLUMN(Eta, eta, float);             //!
DECLARE_SOA_COLUMN(Phi, phi, float);             //!
DECLARE_SOA_COLUMN(PhiV, phiv, float);           //!
DECLARE_SOA_COLUMN(PairDCAxy, pairDCAxy, float); //!
DECLARE_SOA_COLUMN(PairDCAz, pairDCAz, float);   //!
DECLARE_SOA_COLUMN(CosOpAng, cosOpAng, float);   //!
DECLARE_SOA_COLUMN(CosPA, cosPA, float);         //!
DECLARE_SOA_COLUMN(Lxy, lxy, float);             //!
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);     //!
DECLARE_SOA_COLUMN(IsSM, isSM, bool);            //!
DECLARE_SOA_COLUMN(IsHF, isHF, int);             //!
DECLARE_SOA_COLUMN(PairType, pairtype, int);     //!
DECLARE_SOA_COLUMN(IsPrompt, isPrompt, bool);    //!
} // namespace mypair

// reconstructed track information
DECLARE_SOA_TABLE(MyPairs, "AOD", "MYPAIR", //!
                  o2::soa::Index<>, mypair::MyCollisionId, mypair::PosTrackId, mypair::NegTrackId,
                  mypair::Mass, mypair::Pt, mypair::Eta, mypair::Phi, mypair::PhiV, mypair::PairDCAxy, mypair::PairDCAz,
                  mypair::CosOpAng, mypair::CosPA, mypair::Lxy, mypair::Chi2PCA,
                  mypair::IsSM, mypair::IsHF, mypair::PairType, mypair::IsPrompt,
                  mcparticle::PdgCode, mcparticle::Flags,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz,

                  // dynamic column
                  mcparticle::ProducedByGenerator<mcparticle::Flags>,
                  mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                  mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

// iterators
using MyPair = MyPairs::iterator;

} // namespace o2::aod

struct TreeCreatorElectronML {
  enum EM_EEPairType {
    kULS = 0,
    kLSpp = 1,
    kLSnn = 2,
  };

  SliceCache cache;
  Produces<o2::aod::MyCollisions> mycollision;
  Produces<o2::aod::MyPairs> mypair;
  Produces<o2::aod::MyTracks> mytrack;

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{2, 0.f, 2.f}}}},
    },
  };

  // Configurables

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};

  // collision
  Configurable<float> maxVtxZ{"maxVtxZ", 10.0, "max VtxZ [cm]"};
  Configurable<int> maxTrackOccupancy{"maxTrackOccupancy", 999999, "max. track occupancy"};
  Configurable<float> maxFT0Occupancy{"maxFT0Occupancy", 999999., "max. FT0 occupancy"};

  // track
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<int> mintpcclusters{"mintpcclusters", 90, "min. tpc clusters"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 5.0, "max. chi2/NclsITS"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> minpt{"minpt", 0.2, "min. pT"};
  Configurable<float> maxDcaZ{"maxDcaZ", 1.0, "max DCA Z"};
  Configurable<float> maxDcaXY{"maxDcaXY", 1.0, "max DCA XY"};
  Configurable<uint> minITSClusters{"minITSLayers", 5, "min. of ITS clusters"};
  Configurable<uint> minITSClustersIB{"minITSClustersIB", 3, "min. number of ITS clusters in inner barrel"};
  Configurable<float> downSampleEl{"downSampleEl", 1.0, "down scaling factor for electrons"};
  Configurable<float> downSamplePi{"downSamplePi", 1.0, "down scaling factor for pions"};
  Configurable<float> downSampleKa{"downSampleKa", 1.0, "down scaling factor for kaons"};
  Configurable<float> downSamplePr{"downSamplePr", 1.0, "down scaling factor for protons"};
  Configurable<float> downSampleMu{"downSampleMu", 1.0, "down scaling factor for muons"};
  Configurable<float> downSampleNucl{"downSampleNucl", 1.0, "down scaling factor for nuclei"};

  // pid
  Configurable<float> maxNSigmaElTPC{"maxNSigmaElTPC", 999, "max nsigma_el for TPC"};
  Configurable<float> maxNSigmaElTOF{"maxNSigmaElTOF", 999, "max nsgima_el for TOF (if available)"};

  // pair
  Configurable<float> minPairPt{"minPairPt", 0.2, "min. pT,ee"};
  Configurable<float> minMass{"minMass", 0.0, "min. pair invariant mass"};
  Configurable<float> maxMass{"maxMass", 9999.0, "max. pair invariant mass"};
  Configurable<bool> doLS{"doLS", true, "process also LS spectra"};
  Configurable<double> combBgReductionFactor{"combBgReductionFactor", 1.0, "reduction factor for combinatorial background"};

  int mRunNumber;
  float d_bz;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::MatLayerCylSet* lut = nullptr;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

  void init(InitContext&)
  {

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);

    mRunNumber = 0;
    d_bz = 0;
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
    fitter.setWeightedFinalPCA(d_UseWeightedPCA);

    if (useMatCorrType == 1) {
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    } else if (useMatCorrType == 2) {
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    }
    fitter.setMatCorrType(matCorr);
  }

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      fitter.setBz(d_bz);
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter.setBz(d_bz);

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int FindCommonMotherFrom2Prongs(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
  {
    if (!p1.has_mothers())
      return -1;
    if (!p2.has_mothers())
      return -1;

    // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

    int motherid1 = p1.mothersIds()[0];
    auto mother1 = mcparticles.iteratorAt(motherid1);
    int mother1_pdg = mother1.pdgCode();

    int motherid2 = p2.mothersIds()[0];
    auto mother2 = mcparticles.iteratorAt(motherid2);
    int mother2_pdg = mother2.pdgCode();

    // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

    if (motherid1 != motherid2)
      return -1;
    if (mother1_pdg != mother2_pdg)
      return -1;

    if (abs(mother1_pdg) != 22        // photon
        && abs(mother1_pdg) != 111    // pi0
        && abs(mother1_pdg) != 221    // eta
        && abs(mother1_pdg) != 331    // eta'
        && abs(mother1_pdg) != 113    // rho
        && abs(mother1_pdg) != 223    // omega
        && abs(mother1_pdg) != 333    // phi
        && abs(mother1_pdg) != 443    // Jpsi
        && abs(mother1_pdg) != 100443 // psi2S
    ) {
      return -1;
    }
    return motherid1;
  }

  template <typename TTrack>
  bool IsSelected(TTrack const& track)
  {
    if (!track.hasITS()) {
      return false;
    }
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }
    if (track.tpcNClsFound() < mintpcclusters) {
      return false;
    }
    if (track.itsChi2NCl() < -1) { // if tracks are not reconstructed properly, chi2/ITSncls is set to -999;
      return false;
    }
    if (track.itsNCls() < minITSClusters) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < minITSClustersIB) {
      return false;
    }

    if (fabs(track.tpcNSigmaEl()) > maxNSigmaElTPC) {
      return false;
    }
    if (track.hasTOF()) {
      if (fabs(track.tofNSigmaEl()) > maxNSigmaElTOF) {
        return false;
      }
    }
    return true;
  }

  template <typename TMom>
  bool IsSelectedPair(TMom const& v12)
  {
    if (v12.Pt() < minPairPt) {
      return false;
    }
    if (v12.M() < minMass) {
      return false;
    }
    if (v12.M() > maxMass) {
      return false;
    }
    return true;
  }

  bool downSample(int pdg)
  { // pdg=0: combinatorial background for pairs
    float factor = 1.;
    switch (pdg) {
      case 11:
        factor = downSampleEl;
        break;
      case 211:
        factor = downSamplePi;
        break;
      case 321:
        factor = downSampleKa;
        break;
      case 2212:
        factor = downSamplePr;
        break;
      case 13:
        factor = downSampleMu;
        break;
      case 0:
        factor = combBgReductionFactor;
        break;
      default:
        factor = downSampleNucl;
    }
    if (dist01(engine) <= factor) {
      return true;
    }
    return false;
  }

  template <typename TTrack, typename TMCParticles, typename TCollision>
  void doPair(TTrack const& t1, TTrack const& t2, int mc1_id, int mc2_id, int pairtype, TMCParticles const& mctracks, TCollision const& collision, std::map<uint64_t, int>& fNewLabels, std::vector<uint64_t>& fSelected_old_labels, int& fCounter, uint64_t collisionId, std::vector<uint64_t>& collisions_old_labels, int& collisions_counter)
  {
    if (!IsSelected(t1) || !IsSelected(t2)) {
      return;
    }

    auto mc1 = mctracks.iteratorAt(mc1_id);
    auto mc2 = mctracks.iteratorAt(mc2_id);

    if (abs(mc1.pdgCode()) != 11 || abs(mc2.pdgCode()) != 11) {
      return;
    } // charge swap is accepted.

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    if (!IsSelectedPair(v12)) {
      return;
    }

    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), d_bz);
    float pair_dca_xy = sqrt((pow(t1.dcaXY() / sqrt(t1.cYY()), 2) + pow(t2.dcaXY() / sqrt(t2.cYY()), 2)) / 2.);
    float pair_dca_z = sqrt((pow(t1.dcaZ() / sqrt(t1.cZZ()), 2) + pow(t2.dcaZ() / sqrt(t2.cZZ()), 2)) / 2.);

    float cosOpAng = (v1.Px() * v2.Px() + v1.Py() * v2.Py() + v1.Pz() * v2.Pz()) / sqrt((v1.Px() * v1.Px() + v1.Py() * v1.Py() + v1.Pz() * v1.Pz()) * (v2.Px() * v2.Px() + v2.Py() * v2.Py() + v2.Pz() * v2.Pz()));
    float pca = 999.f;
    float lxy = 999.f;
    float cosPA = 999.f;
    o2::aod::pwgem::dilepton::utils::pairutil::isSVFound(fitter, collision, t1, t1, pca, lxy, cosPA);
    float lxy_proper = lxy * v12.M() / v12.Pt();

    bool isSM = false;
    int isHF = o2::aod::pwgem::dilepton::utils::mcutil::IsHF(mc1, mc2, mctracks); // if isHF == true, pdgCode is set to 0, because this pair is correlated HF ee pair decayed from different 2 mothers. Check pdg code of legs.

    bool is_prompt = false;
    int pdgCode = 0;
    uint8_t flags = 0;
    float vx = 0.f;
    float vy = 0.f;
    float vz = 0.f;
    bool is_comb_bg = false;
    if (pairtype == EM_EEPairType::kULS) {
      int common_mother_id = FindCommonMotherFrom2Prongs(mc1, mc2, mctracks);
      if (common_mother_id > 0) {
        isSM = true;
      }
      if (isSM) {
        auto mcpair = mctracks.iteratorAt(common_mother_id);
        is_prompt = true; // only relevant for prompt jpsi
        pdgCode = mcpair.pdgCode();
        flags = mcpair.flags();
        vx = mcpair.vx();
        vy = mcpair.vy();
        vz = mcpair.vz();

        int motherid_pair = mcpair.mothersIds()[0]; // first mother index
        while (motherid_pair > -1) {
          if (motherid_pair < mctracks.size()) { // protect against bad mother indices. why is this needed?
            auto mp = mctracks.iteratorAt(motherid_pair);

            if (std::to_string(mp.pdgCode()).find("5") != std::string::npos) {
              is_prompt = false;
              break;
            }

            if (mp.has_mothers()) {
              motherid_pair = mp.mothersIds()[0];
            } else {
              motherid_pair = -999;
            }
          } else {
            LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid_pair, mctracks.size());
          }
        }
      } else if (isHF == static_cast<int>(o2::aod::pwgem::dilepton::utils::mcutil::EM_HFeeType::kUndef)) {
        is_comb_bg = true;
      }
    } else { // LS
      if (isHF == static_cast<int>(o2::aod::pwgem::dilepton::utils::mcutil::EM_HFeeType::kUndef)) {
        is_comb_bg = true;
      }
    }

    if (!is_comb_bg || downSample(0)) {
      if (!(fNewLabels.find(t1.globalIndex()) != fNewLabels.end())) {
        fNewLabels[t1.globalIndex()] = fCounter;
        fSelected_old_labels.push_back(t1.globalIndex());
        fCounter++;
      }
      if (!(fNewLabels.find(t2.globalIndex()) != fNewLabels.end())) {
        fNewLabels[t2.globalIndex()] = fCounter;
        fSelected_old_labels.push_back(t2.globalIndex());
        fCounter++;
      }
      if (find(collisions_old_labels.begin(), collisions_old_labels.end(), collisionId) == collisions_old_labels.end()) {
        collisions_counter++;
        collisions_old_labels.push_back(collisionId);
      }
      mypair(collisions_counter, fNewLabels[t1.globalIndex()], fNewLabels[t2.globalIndex()],
             v12.M(), v12.Pt(), v12.Eta(), v12.Phi(), phiv, pair_dca_xy, pair_dca_z, cosOpAng, cosPA, lxy_proper, pow(pca, 2),
             isSM, isHF, pairtype, is_prompt, pdgCode, flags,
             vx, vy, vz);
    }
  }

  template <bool isDerived, typename TTrack, typename TMCParticle, typename TMCParticles>
  void doSingleTrack(TTrack& track, TMCParticle& mctrack, TMCParticles& mctracks, uint64_t collisionId, std::vector<uint64_t>& collisions_old_labels, int& collisions_counter, bool use_downsample = true)
  {
    if (!IsSelected(track)) {
      return;
    }

    if ((!use_downsample) || (downSample(abs(mctrack.pdgCode())))) {
      // store all mother relation
      std::vector<int> mothers_id;
      std::vector<int> mothers_pdg;
      if (mctrack.has_mothers()) {
        int motherid = mctrack.mothersIds()[0]; // first mother index
        while (motherid > -1) {
          if (motherid < mctracks.size()) { // protect against bad mother indices. why is this needed?
            auto mp = mctracks.iteratorAt(motherid);
            mothers_id.emplace_back(motherid);
            mothers_pdg.emplace_back(mp.pdgCode());

            if (mp.has_mothers()) {
              motherid = mp.mothersIds()[0];
            } else {
              motherid = -999;
            }
          } else {
            LOGF(info, "Mother label(%d) exceeds the McParticles size(%d)", motherid, mctracks.size());
          }
        }
      }
      if (find(collisions_old_labels.begin(), collisions_old_labels.end(), collisionId) == collisions_old_labels.end()) {
        collisions_counter++;
        collisions_old_labels.push_back(collisionId);
      }

      if constexpr (isDerived) {
        mytrack(collisions_counter,
                track.sign(), track.pt(), track.eta(), track.phi(), track.tgl(), track.dcaXY(), track.dcaZ(), sqrt(track.cYY()), sqrt(track.cZZ()), track.cZY(),
                track.tpcNClsFindable(), track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                track.tpcChi2NCl(), track.tpcInnerParam(),
                track.tpcSignal(), track.tpcNSigmaEl(), /*track.tpcNSigmaMu(),*/ track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                track.beta(), track.tofNSigmaEl(), /*track.tofNSigmaMu(),*/ 0, 0, 0,
                track.tofChi2(), track.itsChi2NCl(), track.itsClusterSizes(),
                mctrack.vx(), mctrack.vy(), mctrack.vz(),
                mctrack.pdgCode(), mctrack.isPhysicalPrimary(), mothers_id, mothers_pdg);
      } else {
        mytrack(collisions_counter,
                track.sign(), track.pt(), track.eta(), track.phi(), track.tgl(), track.dcaXY(), track.dcaZ(), sqrt(track.cYY()), sqrt(track.cZZ()), track.cZY(),
                track.tpcNClsFindable(), track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                track.tpcChi2NCl(), track.tpcInnerParam(),
                track.tpcSignal(), track.tpcNSigmaEl(), /*track.tpcNSigmaMu(),*/ track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                track.beta(), track.tofNSigmaEl(), /*track.tofNSigmaMu(),*/ track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                track.tofChi2(), track.itsChi2NCl(), track.itsClusterSizes(),
                mctrack.vx(), mctrack.vy(), mctrack.vz(),
                mctrack.pdgCode(), mctrack.isPhysicalPrimary(), mothers_id, mothers_pdg);
      }

      mothers_id.shrink_to_fit();
      mothers_pdg.shrink_to_fit();
    }
  }

  template <typename TCollision, typename TMCCollision>
  void doCollision(TCollision& collision, TMCCollision& mccollision, uint64_t globalBC, int runNumber)
  {
    registry.fill(HIST("hEventCounter"), 1.5);
    mycollision(globalBC, runNumber, collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange(), collision.sel8(), d_bz,
                mccollision.generatorsID(), mccollision.posX(), mccollision.posY(), mccollision.posZ(),
                collision.multNTracksPV(), collision.centFT0M(), collision.centFT0A(), collision.centFT0C());
  }

  template <typename TCollision, typename TBC, typename TMCCollision>
  void doCollision(TCollision& collision, TBC& bc, TMCCollision& mccollision)
  {
    doCollision(collision, mccollision, bc.globalBC(), bc.runNumber());
  }

  Filter trackFilter = o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& o2::aod::track::itsChi2NCl < maxchi2its&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& nabs(o2::aod::track::dcaXY) < maxDcaXY&& nabs(o2::aod::track::dcaZ) < maxDcaZ;
  using MyFilteredTracksMC = soa::Filtered<FullTracksExtMC>;
  Preslice<MyFilteredTracksMC> perCollision = aod::track::collisionId;

  Filter collisionFilter = nabs(o2::aod::collision::posZ) < maxVtxZ && o2::aod::evsel::trackOccupancyInTimeRange < maxTrackOccupancy && o2::aod::evsel::ft0cOccupancyInTimeRange < maxFT0Occupancy;
  using MyFilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;

  void processSingleTrack(MyFilteredCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks, aod::McParticles const& mctracks, aod::McCollisions const&)
  {
    int collisions_counter = -1;
    std::vector<uint64_t> collisions_old_labels;

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 0.5);

      auto tracks_coll = tracks.sliceBy(perCollision, collision.globalIndex());
      for (auto& track : tracks_coll) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto mctrack = track.mcParticle_as<aod::McParticles>();
        doSingleTrack<false>(track, mctrack, mctracks, collision.globalIndex(), collisions_old_labels, collisions_counter);
      }
    }

    for (uint64_t collision_label : collisions_old_labels) {
      auto collision = collisions.rawIteratorAt(collision_label);
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      auto mccollision = collision.mcCollision();
      doCollision(collision, bc, mccollision);
    }
  }
  PROCESS_SWITCH(TreeCreatorElectronML, processSingleTrack, "produce ML input for single track level", false);

  using MyFilteredCollisionsSkimmed = soa::Filtered<soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsXY, aod::EMMCEventLabels>>;
  using MyFilteredTracksMCSkimmed = soa::Filtered<soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronsCov, aod::EMPrimaryElectronMCLabels, aod::EMPrimaryElectronEMEventIds>>;
  Preslice<MyFilteredTracksMCSkimmed> perCollisionSkimmed = aod::emprimaryelectron::emeventId;

  void processSingleTrackSkimmed(MyFilteredCollisionsSkimmed const& collisions, MyFilteredTracksMCSkimmed const& tracks, aod::EMMCParticles const& mctracks, aod::EMMCEvents const&)
  {
    int collisions_counter = -1;
    std::vector<uint64_t> collisions_old_labels;

    for (auto& collision : collisions) {
      if (!collision.has_emmcevent()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 0.5);

      auto tracks_coll = tracks.sliceBy(perCollisionSkimmed, collision.globalIndex());
      for (auto& track : tracks_coll) {
        if (!track.has_emmcparticle()) {
          continue;
        }
        auto mctrack = track.emmcparticle_as<aod::EMMCParticles>();
        doSingleTrack<true>(track, mctrack, mctracks, collision.globalIndex(), collisions_old_labels, collisions_counter);
      }
    }

    for (uint64_t collision_label : collisions_old_labels) {
      auto collision = collisions.rawIteratorAt(collision_label);
      auto mccollision = collision.emmcevent();
      doCollision(collision, mccollision, collision.globalBC(), collision.runNumber());
    }
  }
  PROCESS_SWITCH(TreeCreatorElectronML, processSingleTrackSkimmed, "produce ML input for single track level on skimmed data", false);

  Partition<MyFilteredTracksMC> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracks = o2::aod::track::signed1Pt < 0.f;

  void processPair(MyFilteredCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks, aod::McParticles const& mctracks, aod::McCollisions const&)
  {
    std::map<uint64_t, int> fNewLabels;
    int fCounter = 0;

    std::vector<uint64_t> collisions_old_labels;
    int collisions_counter = -1;

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 0.5);
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      std::vector<uint64_t> fSelected_old_labels;

      auto negTracks_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto posTracks_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(negTracks_coll, posTracks_coll))) {
        if (!ele.has_mcParticle() || !pos.has_mcParticle()) {
          continue;
        }
        doPair(pos, ele, pos.mcParticleId(), ele.mcParticleId(), EM_EEPairType::kULS, mctracks, collision, fNewLabels, fSelected_old_labels, fCounter, collision.globalIndex(), collisions_old_labels, collisions_counter);
      }

      if (!doLS) {
        continue;
      }
      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_coll, posTracks_coll))) {
        if (!pos1.has_mcParticle() || !pos2.has_mcParticle()) {
          continue;
        }
        doPair(pos1, pos2, pos1.mcParticleId(), pos2.mcParticleId(), EM_EEPairType::kLSpp, mctracks, collision, fNewLabels, fSelected_old_labels, fCounter, collision.globalIndex(), collisions_old_labels, collisions_counter);
      }

      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_coll, negTracks_coll))) {
        if (!ele1.has_mcParticle() || !ele2.has_mcParticle()) {
          continue;
        }
        doPair(ele1, ele2, ele1.mcParticleId(), ele2.mcParticleId(), EM_EEPairType::kLSnn, mctracks, collision, fNewLabels, fSelected_old_labels, fCounter, collision.globalIndex(), collisions_old_labels, collisions_counter);
      }

      // single tracks, only if selected in at least one pair
      for (uint64_t track_label : fSelected_old_labels) {
        auto track = tracks.rawIteratorAt(track_label);
        auto mctrack = track.mcParticle_as<aod::McParticles>();

        doSingleTrack<false>(track, mctrack, mctracks, collision.globalIndex(), collisions_old_labels, collisions_counter, false);

      } // end of track loop

    } // end of collision loop

    for (uint64_t collision_label : collisions_old_labels) {
      auto collision = collisions.rawIteratorAt(collision_label);
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      auto mccollision = collision.mcCollision();
      doCollision(collision, bc, mccollision);
    }

    fNewLabels.clear();
    fCounter = 0;
  } // end of process
  PROCESS_SWITCH(TreeCreatorElectronML, processPair, "produce ML input for pair level", false);

  Partition<MyFilteredTracksMCSkimmed> posTracksSkimmed = o2::aod::emprimaryelectron::sign > static_cast<int8_t>(0);
  Partition<MyFilteredTracksMCSkimmed> negTracksSkimmed = o2::aod::emprimaryelectron::sign < static_cast<int8_t>(0);

  void processPairSkimmed(MyFilteredCollisionsSkimmed const& collisions, MyFilteredTracksMCSkimmed const& tracks, aod::EMMCParticles const& mctracks, aod::EMMCEvents const&)
  {
    std::map<uint64_t, int> fNewLabels;
    int fCounter = 0;

    std::vector<uint64_t> collisions_old_labels;
    int collisions_counter = -1;
    std::vector<uint64_t> collisions_old_labels_track = {};
    int collisions_counter_track = -1;

    for (auto& collision : collisions) {
      if (!collision.has_emmcevent()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 0.5);

      initCCDB(collision);
      std::vector<uint64_t> fSelected_old_labels;

      auto negTracks_coll = negTracksSkimmed->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
      auto posTracks_coll = posTracksSkimmed->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);

      for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(negTracks_coll, posTracks_coll))) {
        if (!ele.has_emmcparticle() || !pos.has_emmcparticle()) {
          continue;
        }
        doPair(pos, ele, pos.emmcparticleId(), ele.emmcparticleId(), EM_EEPairType::kULS, mctracks, collision, fNewLabels, fSelected_old_labels, fCounter, collision.globalIndex(), collisions_old_labels, collisions_counter);
      }

      if (!doLS) {
        continue;
      }
      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posTracks_coll, posTracks_coll))) {
        if (!pos1.has_emmcparticle() || !pos2.has_emmcparticle()) {
          continue;
        }
        doPair(pos1, pos2, pos1.emmcparticleId(), pos2.emmcparticleId(), EM_EEPairType::kLSpp, mctracks, collision, fNewLabels, fSelected_old_labels, fCounter, collision.globalIndex(), collisions_old_labels, collisions_counter);
      }

      for (auto& [ele1, ele2] : combinations(CombinationsStrictlyUpperIndexPolicy(negTracks_coll, negTracks_coll))) {
        if (!ele1.has_emmcparticle() || !ele2.has_emmcparticle()) {
          continue;
        }
        doPair(ele1, ele2, ele1.emmcparticleId(), ele2.emmcparticleId(), EM_EEPairType::kLSnn, mctracks, collision, fNewLabels, fSelected_old_labels, fCounter, collision.globalIndex(), collisions_old_labels, collisions_counter);
      }

      // single tracks, only if selected in at least one pair
      for (uint64_t track_label : fSelected_old_labels) {
        auto track = tracks.rawIteratorAt(track_label);
        auto mctrack = track.emmcparticle_as<aod::EMMCParticles>();

        doSingleTrack<true>(track, mctrack, mctracks, collision.globalIndex(), collisions_old_labels_track, collisions_counter_track, false);

      } // end of track loop

    } // end of collision loop

    for (uint64_t collision_label : collisions_old_labels) {
      auto collision = collisions.rawIteratorAt(collision_label);
      auto mccollision = collision.emmcevent();
      doCollision(collision, mccollision, collision.globalBC(), collision.runNumber());
    }

    fNewLabels.clear();
    fCounter = 0;
  } // end of process
  PROCESS_SWITCH(TreeCreatorElectronML, processPairSkimmed, "produce ML input for pair level on skimmed data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorElectronML>(cfgc, TaskName{"tree-creator-ele-ml"})};
}
