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
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsVertexing/PVertexer.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "Math/Vector4D.h"
#include "TDatabasePDG.h"
#include "THnSparse.h"
#include "TParticlePDG.h"
#include "TTree.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>
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
struct NPCascCandidate {
  int runNumber;
  int64_t mcParticleId;
  int64_t trackGlobID;
  int64_t trackITSID;
  int64_t collisionID;
  float matchingChi2;
  float deltaPtITS;
  float deltaPt;
  float itsClusSize;
  bool hasReassociatedCluster;
  bool hasFakeReassociation;
  bool isGoodMatch;
  bool isGoodCascade;
  int pdgCodeMom;
  int pdgCodeITStrack;
  bool isFromBeauty;
  bool isFromCharm;
  uint16_t pvContributors;
  uint8_t cascPVContribs;
  float pvTimeResolution;
  float pvX;
  float pvY;
  float pvZ;
  float cascPt;
  float cascEta;
  float cascPhi;
  float protonPt;
  float protonEta;
  float pionPt;
  float pionEta;
  float bachPt;
  float bachEta;
  float cascDCAxy;
  float cascDCAz;
  float protonDCAxy;
  float protonDCAz;
  float pionDCAxy;
  float pionDCAz;
  float bachDCAxy;
  float bachDCAz;
  float casccosPA;
  float v0cosPA;
  double massXi;
  double massOmega;
  double massV0;
  float cascRadius;
  float v0radius;
  float cascLength;
  float v0length;
  int cascNClusITS;
  int protonNClusITS;
  int pionNClusITS;
  int bachNClusITS;
  int protonNClusTPC;
  int pionNClusTPC;
  int bachNClusTPC;
  float protonTPCNSigma;
  float pionTPCNSigma;
  float bachKaonTPCNSigma;
  float bachPionTPCNSigma;
  bool protonHasTOF;
  bool pionHasTOF;
  bool bachHasTOF;
  float protonTOFNSigma;
  float pionTOFNSigma;
  float bachKaonTOFNSigma;
  float bachPionTOFNSigma;
  bool sel8;
  float multFT0C;
  float multFV0A;
  float multFT0M;
  float centFT0C;
  float centFV0A;
  float centFT0M;
  int multNTracksGlobal;
  uint32_t toiMask;
  bool noSameBunchPileup;
};
std::array<bool, 2> isFromHF(auto& particle)
{
  bool fromBeauty = false;
  bool fromCharm = false;
  if (particle.has_mothers()) {
    auto mom = particle.template mothers_as<aod::McParticles>()[0];
    int pdgCodeMom = mom.pdgCode();
    fromBeauty = std::abs(pdgCodeMom) / 5000 == 1 || std::abs(pdgCodeMom) / 500 == 1 || std::abs(pdgCodeMom) == 5;
    fromCharm = std::abs(pdgCodeMom) / 4000 == 1 || std::abs(pdgCodeMom) / 400 == 1 || std::abs(pdgCodeMom) == 4;
    while (mom.has_mothers()) {
      const auto grandma = mom.template mothers_as<aod::McParticles>()[0];
      int pdgCodeGrandma = std::abs(grandma.pdgCode());
      fromBeauty = fromBeauty || (pdgCodeGrandma / 5000 == 1 || pdgCodeGrandma / 500 == 1 || pdgCodeGrandma == 5);
      fromCharm = fromCharm || (pdgCodeGrandma / 4000 == 1 || pdgCodeGrandma / 400 == 1 || pdgCodeGrandma == 4);
      mom = grandma;
    }
  }
  return {fromBeauty, fromCharm};
}

static constexpr int nParticles{4};
static constexpr int nCutsPID{2};
static const std::vector<std::string> matterOrNot{"Matter", "Antimatter"};
static const std::vector<std::string> particlesNames{"K-bachelor", "Pi-bachelor", "Pr", "Pi"};
static const std::vector<std::string> cutsNames{"TPCnSigmaMin", "TPCnSigmaMax"};
static constexpr float cutsPID[nParticles][nCutsPID]{
  {-4.f, +4.f}, /*K bachelor*/
  {-4.f, +4.f}, /*Pi bachelor*/
  {-4.f, +4.f}, /*Pr*/
  {-4.f, +4.f}, /*Pi*/
};

std::vector<NPCascCandidate> gCandidates;
std::vector<NPCascCandidate> gCandidatesNT;

} // namespace

struct NonPromptCascadeTask {

  Produces<o2::aod::NPCascTable> NPCTable;
  Produces<o2::aod::NPCascTableMC> NPCTableMC;
  Produces<o2::aod::NPCascTableNT> NPCTableNT;
  Produces<o2::aod::NPCascTableMCNT> NPCTableMCNT;
  Produces<o2::aod::NPCascTableGen> NPCTableGen;

  using TracksExtData = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullKa, aod::pidTOFFullPi, aod::pidTOFFullPr>;
  using TracksExtMC = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels, aod::pidTPCFullKa, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTOFFullKa, aod::pidTOFFullPi, aod::pidTOFFullPr>;
  using CollisionCandidatesRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0Ms, aod::MultsGlobal>;
  using CollisionCandidatesRun3MC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0Ms, aod::MultsGlobal>;
  using CollisionsWithLabel = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::MultsGlobal>;
  using TracksWithLabel = soa::Join<aod::Tracks, aod::McTrackLabels>;

  Preslice<TracksExtData> perCollision = aod::track::collisionId;
  Preslice<TracksExtMC> perCollisionMC = aod::track::collisionId;

  HistogramRegistry mRegistry;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> cfgPropToPCA{"cfgPropToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> cfgRedoPV{"cfgRedoPV", true, "redo PV"};
  Configurable<bool> cfgUseAbsDCA{"cfgUseAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<double> cfgMaxR{"cfgMaxR", 200., "reject PCA's above this radius"};
  Configurable<double> cfgMaxDZIni{"cfgMaxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> cfgMinParamChange{"cfgMinParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> cfgMinRelChi2Change{"cfgMinRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> cfgSelectOnlyOmegas{"cfgSelectOnlyOmegas", false, "Toggle to select only Omegas"};

  Configurable<int> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};
  Configurable<float> cfgMinCosPA{"cfgMinCosPA", -1.f, "Minimum cosine of pointing angle"};
  Configurable<LabeledArray<float>> cfgCutsPID{"particlesCutsPID", {cutsPID[0], nParticles, nCutsPID, particlesNames, cutsNames}, "Nuclei PID selections"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", true, "Skimmed dataset processing"};

  Configurable<std::string> cfgTriggersOfInterest{"cfgTriggersOfInterest", "fTrackedOmega,fOmegaHighMult,fHighFt0Mult", "Triggers of interest, comma separated for Zorro"};

  Configurable<float> cfgMaxMult{"cfgMaxMult", 8000.f, "Upper range of multiplicty histo"};
  Configurable<float> cfgMaxMultFV0{"cfgMaxMultFV0", 10000.f, "Upper range of multiplicty FV0 histo"};
  Configurable<std::string> cfgPtEdgesdNdeta{"ptEdges", "0,0.2,0.4,0.6,0.8,1,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4,4.5,5,5.5,6,7,8,10", "Pt bin edges (comma-separated)"};

  Zorro mZorro;
  OutputObj<ZorroSummary> mZorroSummary{"ZorroSummary"};
  SliceCache cache;

  Service<o2::ccdb::BasicCCDBManager> mCCDB;
  int mRunNumber = 0;
  float mBz = 0.f;
  o2::vertexing::DCAFitterN<2> mDCAFitter;
  std::array<int, 2> mProcessCounter = {0, 0}; // {Tracked, All}
  std::map<uint64_t, uint32_t> mToiMap;
  //
  HistogramRegistry mRegistryMults{"Multhistos"};
  HistogramRegistry mRegistrydNdeta{"dNdetahistos"};

  //

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

    if (o2::parameters::GRPMagField* grpmag = mCCDB->getForRun<o2::parameters::GRPMagField>(cfgGRPmagPath, mRunNumber)) {
      o2::base::Propagator::initFieldFromGRP(grpmag);
      mBz = static_cast<float>(grpmag->getNominalL3Field());
    }
    mDCAFitter.setBz(mBz);

    if (static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(mCCDB->getForRun<o2::base::MatLayerCylSet>("GLO/Param/MatLUT", mRunNumber));
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  void init(InitContext& context)
  {
    mZorroSummary.setObject(mZorro.getZorroSummary());
    mCCDB->setURL(ccdbUrl);
    mCCDB->setFatalWhenNull(true);
    mCCDB->setCaching(true);
    mCCDB->setLocalObjectValidityChecking();

    mDCAFitter.setPropagateToPCA(cfgPropToPCA);
    mDCAFitter.setMaxR(cfgMaxR);
    mDCAFitter.setMaxDZIni(cfgMaxDZIni);
    mDCAFitter.setMinParamChange(cfgMinParamChange);
    mDCAFitter.setMinRelChi2Change(cfgMinRelChi2Change);
    mDCAFitter.setUseAbsDCA(cfgUseAbsDCA);

    std::vector<double> ptBinning = {0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0};
    // AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    std::array<std::string, 7> cutsNames{"# candidates", "hasTOF", "nClusTPC", "nSigmaTPCbach", "nSigmaTPCprotontrack", "nSigmaTPCpiontrack", "cosPA"};
    auto cutsOmega{std::get<std::shared_ptr<TH2>>(mRegistry.add("h_PIDcutsOmega", ";;Invariant mass (GeV/#it{c}^{2})", HistType::kTH2D, {{cutsNames.size(), -0.5, -0.5 + cutsNames.size()}, {125, 1.650, 1.700}}))};
    auto cutsXi{std::get<std::shared_ptr<TH2>>(mRegistry.add("h_PIDcutsXi", ";;Invariant mass (GeV/#it{c}^{2})", HistType::kTH2D, {{6, -0.5, 5.5}, {125, 1.296, 1.346}}))};
    for (size_t iBin{0}; iBin < cutsNames.size(); ++iBin) {
      cutsOmega->GetYaxis()->SetBinLabel(iBin + 1, cutsNames[iBin].c_str());
      cutsXi->GetYaxis()->SetBinLabel(iBin + 1, cutsNames[iBin].c_str());
    }
    //
    int nBinsMult = cfgMaxMult;
    int nBinsMultFV0 = cfgMaxMultFV0;

    AxisSpec multAxis = {nBinsMult, 0, cfgMaxMult, "Multiplicity FT0M"};
    AxisSpec multAxisFV0 = {nBinsMultFV0, 0, cfgMaxMultFV0, "Multiplicity FT0M"};
    AxisSpec nTracksAxis = {100, 0., 100., "NTracksGlobal"};
    AxisSpec nTracksAxisMC = {100, 0., 100., "NTracksMC"};
    std::vector<double> centBinning;
    std::vector<double> multBinning;
    std::vector<double> trackBinning;
    std::vector<double> runsBinning;
    double cent = -0.0;
    while (cent < 100) {
      if (cent < 0.1) {
        centBinning.push_back(cent);
        cent += 0.01;
      } else if (cent < 1) {
        centBinning.push_back(cent);
        cent += 0.1;
      } else {
        centBinning.push_back(cent);
        cent += 1;
      }
    }
    double ntracks = 0;
    while (ntracks < 100) {
      trackBinning.push_back(ntracks);
      ntracks++;
    }
    double run = 550367 - 0.5;
    for (int i = 0; i < 9000; i++) {
      runsBinning.push_back(run);
      run++;
    }
    AxisSpec centAxis{centBinning, "Centrality (%)"};
    // AxisSpec multAxis{multBinning, "Multiplicity"};
    AxisSpec trackAxisMC{trackBinning, "NTracks MC"};
    AxisSpec trackAxis{trackBinning, "NTracks Global Reco"};
    AxisSpec runsAxis{runsBinning, "Run Number"};

    mRegistryMults.add("hCentMultsRuns", "hCentMultsRuns", HistType::kTHnSparseF, {centAxis, multAxis, centAxis, multAxisFV0, nTracksAxis, runsAxis});
    //
    // dN/deta
    //
    bool runMCdNdeta = context.options().get<bool>("processdNdetaMC");
    // std::cout << "runMCdNdeta: " << runMCdNdeta << std::endl;
    if (runMCdNdeta) {
      std::vector<double> ptBins;
      std::vector<std::string> tokens = o2::utils::Str::tokenize(cfgPtEdgesdNdeta, ',');
      for (auto const& pts : tokens) {
        double pt = 0;
        try {
          pt = std::stof(pts);
        } catch (...) {
          LOG(error) << "Wrong cfgPtEdgesdNdeta string:" << cfgPtEdgesdNdeta << std::endl;
        }
        ptBins.push_back(pt);
      }
      AxisSpec ptAxisMC{ptBins, "pT MC"};
      AxisSpec ptAxisReco{ptBins, "pT Reco"};

      // multMeasured, multMC, ptMeasured, ptMC
      mRegistrydNdeta.add("hdNdetaRM/hdNdetaRM", "hdNdetaRM", HistType::kTHnSparseF, {nTracksAxisMC, nTracksAxis, ptAxisMC, ptAxisReco});
      mRegistrydNdeta.add("hdNdetaRM/hdNdetaRMNotInRecoCol", "hdNdetaRMNotInRecoCol", HistType::kTHnSparseF, {nTracksAxisMC, ptAxisMC});
      mRegistrydNdeta.add("hdNdetaRM/hdNdetaRMNotInRecoTrk", "hdNdetaRMNotInRecoTrk", HistType::kTHnSparseF, {nTracksAxisMC, ptAxisMC});
    }
  }

  template <typename CollisionType, typename TrackType>
  bool recalculatePV(CollisionType const& collision, TrackType const& tracks, int idToRemove, o2::dataformats::VertexBase& primaryVertex)
  {
    // slice tracks by collision
    o2::vertexing::PVertexer vertexer;
    std::vector<o2::track::TrackParCov> pvContributors = {};
    std::vector<bool> pvContributorsMask = {};

    auto tracksInCollision = doprocessTrackedCascadesMC ? tracks.sliceBy(perCollisionMC, collision.globalIndex()) : tracks.sliceBy(perCollision, collision.globalIndex());
    // loop over tracks
    for (auto const& trkInColl : tracksInCollision) { // Loop on tracks
      if (trkInColl.isPVContributor()) {
        pvContributors.push_back(getTrackParCov(trkInColl));
        idToRemove == trkInColl.globalIndex() ? pvContributorsMask.push_back(false) : pvContributorsMask.push_back(true);
      }
    }
    LOG(debug) << "Tracks pushed to the vector: " << pvContributors.size();
    vertexer.init();
    bool canRefit = vertexer.prepareVertexRefit(pvContributors, primaryVertex);
    if (!canRefit) {
      return false;
    }
    // refit the vertex
    auto newPV = vertexer.refitVertex(pvContributorsMask, primaryVertex);
    // set the new vertex to primaryVertex
    primaryVertex.setX(newPV.getX());
    primaryVertex.setY(newPV.getY());
    primaryVertex.setZ(newPV.getZ());
    primaryVertex.setCov(newPV.getCov());
    return true;
  }

  void zorroAccounting(const auto& collisions)
  {
    if (cfgSkimmedProcessing && mProcessCounter[0] != mProcessCounter[1]) {
      mToiMap.clear();
      int runNumber{-1};
      for (const auto& coll : collisions) {
        auto bc = coll.template bc_as<aod::BCsWithTimestamps>();
        if (runNumber != bc.runNumber()) {
          mZorro.initCCDB(mCCDB.service, bc.runNumber(), bc.timestamp(), cfgTriggersOfInterest.value);
          if (mZorro.getNTOIs() > 32) {
            LOG(fatal) << "N TOIs:" << mZorro.getNTOIs() << " Max 32 TOIs possible.";
          }
          mZorro.populateHistRegistry(mRegistry, bc.runNumber());
          runNumber = bc.runNumber();
        }
        bool sel = mZorro.isSelected(bc.globalBC()); /// Just let Zorro do the accounting
        if (sel) {
          std::vector<bool> toivect = mZorro.getTriggerOfInterestResults();
          uint32_t toiMask = 0;
          for (size_t i{0}; i < toivect.size(); i++) {
            toiMask += toivect[i] << i;
          }
          mToiMap[bc.globalBC()] = toiMask;
        }
      }
    }
  }
  template <typename CollisionType>
  void fillMultHistos(CollisionType const& collisions)
  {
    // std::cout << "Filling mult histos" << std::endl;
    for (const auto& coll : collisions) {
      float centFT0M = coll.centFT0M();
      float multFT0M = coll.multFT0M();
      float centFV0A = coll.centFV0A();
      float multFV0A = coll.multFV0A();
      float multNTracks = coll.multNTracksGlobal();
      float run = mRunNumber;
      mRegistryMults.fill(HIST("hCentMultsRuns"), centFT0M, multFT0M, centFV0A, multFV0A, multNTracks, run);
    }
  };

  template <typename TrackType, typename CollisionType>
  void fillCandidatesVector(CollisionType const&, TrackType const& tracks, auto const& cascades, auto& candidates)
  {

    const auto& getCascade = [](auto const& candidate) {
      if constexpr (requires { candidate.cascade(); }) {
        return candidate.cascade();
      } else {
        return candidate;
      }
    };

    candidates.clear();
    for (const auto& candidate : cascades) {

      auto collision = candidate.template collision_as<CollisionType>();
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      auto primaryVertex = getPrimaryVertex(collision);

      const auto& casc = getCascade(candidate);
      const auto& bachelor = casc.template bachelor_as<TrackType>();
      const auto& v0 = casc.v0();
      const auto& ptrack = v0.template posTrack_as<TrackType>();
      const auto& ntrack = v0.template negTrack_as<TrackType>();
      const auto& protonTrack = bachelor.sign() > 0 ? ntrack : ptrack;
      const auto& pionTrack = bachelor.sign() > 0 ? ptrack : ntrack;

      // first bit for the strange track, second for pos v0, third for neg v0, fourth for bachelor
      uint8_t cascPVContribs = 0;
      cascPVContribs |= ptrack.isPVContributor() << 1;
      cascPVContribs |= ntrack.isPVContributor() << 2;
      cascPVContribs |= bachelor.isPVContributor() << 3;

      mRegistry.fill(HIST("h_PIDcutsXi"), 0, 1.322);
      mRegistry.fill(HIST("h_PIDcutsOmega"), 0, 1.675);

      mRegistry.fill(HIST("h_PIDcutsXi"), 1, 1.322);
      mRegistry.fill(HIST("h_PIDcutsOmega"), 1, 1.675);

      if (protonTrack.tpcNClsFound() < cfgCutNclusTPC || pionTrack.tpcNClsFound() < cfgCutNclusTPC || bachelor.tpcNClsFound() < cfgCutNclusTPC) {
        continue;
      }
      mRegistry.fill(HIST("h_PIDcutsXi"), 2, 1.322);
      mRegistry.fill(HIST("h_PIDcutsOmega"), 2, 1.675);

      // QA PID
      float nSigmaTPC[nParticles]{bachelor.tpcNSigmaKa(), bachelor.tpcNSigmaPi(), protonTrack.tpcNSigmaPr(), pionTrack.tpcNSigmaPi()};

      bool isBachelorSurvived = false;
      if (nSigmaTPC[0] > cfgCutsPID->get(0u, 0u) && nSigmaTPC[0] < cfgCutsPID->get(0u, 1u)) {
        mRegistry.fill(HIST("h_PIDcutsOmega"), 3, 1.675);
        isBachelorSurvived = true;
      }

      if (!cfgSelectOnlyOmegas && nSigmaTPC[1] > cfgCutsPID->get(1u, 0u) && nSigmaTPC[1] < cfgCutsPID->get(1u, 1u)) {
        mRegistry.fill(HIST("h_PIDcutsXi"), 3, 1.322);
        isBachelorSurvived = true;
      }

      if (!isBachelorSurvived) {
        continue;
      }

      if (nSigmaTPC[2] < cfgCutsPID->get(2u, 0u) || nSigmaTPC[2] > cfgCutsPID->get(2u, 1u)) {
        continue;
      }

      mRegistry.fill(HIST("h_PIDcutsOmega"), 4, 1.675);
      mRegistry.fill(HIST("h_PIDcutsXi"), 4, 1.322);

      if (nSigmaTPC[3] < cfgCutsPID->get(3u, 0u) || nSigmaTPC[3] > cfgCutsPID->get(3u, 1u)) {
        continue;
      }

      mRegistry.fill(HIST("h_PIDcutsOmega"), 5, 1.675);
      mRegistry.fill(HIST("h_PIDcutsXi"), 5, 1.322);

      auto protonTrkParCov = getTrackParCov(protonTrack);
      auto pionTrkParCov = getTrackParCov(pionTrack);
      auto bachTrkParCov = getTrackParCov(bachelor);

      std::array<std::array<float, 3>, 2> momenta;
      std::array<float, 3> cascadeMomentum;
      o2::math_utils::SVector<double, 3> cascadePos, v0Pos;

      float cascCpa = -1, v0Cpa = -1;
      o2::track::TrackParCov ntCascadeTrack;
      if (mDCAFitter.process(pionTrkParCov, protonTrkParCov)) {
        auto trackParCovV0 = mDCAFitter.createParentTrackParCov(0); // V0 track retrieved from p and pi daughters
        v0Pos = mDCAFitter.getPCACandidate();
        if (mDCAFitter.process(trackParCovV0, bachTrkParCov)) {
          mDCAFitter.getTrackParamAtPCA(0).getPxPyPzGlo(momenta[0]);
          mDCAFitter.getTrackParamAtPCA(1).getPxPyPzGlo(momenta[1]);
          ntCascadeTrack = mDCAFitter.createParentTrackParCov();
          ntCascadeTrack.getPxPyPzGlo(cascadeMomentum);
          std::array<float, 3> pvPos = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};
          cascadePos = mDCAFitter.getPCACandidate();
          cascCpa = RecoDecay::cpa(pvPos, mDCAFitter.getPCACandidate(), cascadeMomentum);
          v0Cpa = RecoDecay::cpa(pvPos, v0Pos, momenta[0]);
        } else {
          continue;
        }
      } else {
        continue;
      }
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> cascadeLvector;
      cascadeLvector.SetPxPyPzE(cascadeMomentum[0], cascadeMomentum[1], cascadeMomentum[2], std::hypot(cascadeMomentum[0], cascadeMomentum[1], cascadeMomentum[2])); /// 0 mass, used only for the momentum

      // Omega
      std::array<double, 2> masses{o2::constants::physics::MassLambda0, o2::constants::physics::MassKPlus};
      const auto massOmega = RecoDecay::m(momenta, masses);

      // Xi
      masses = {o2::constants::physics::MassLambda0, o2::constants::physics::MassPiPlus};
      const auto massXi = RecoDecay::m(momenta, masses);

      // Lambda
      masses = {o2::constants::physics::MassProton, o2::constants::physics::MassPiMinus};
      momenta[0] = {protonTrack.px(), protonTrack.py(), protonTrack.pz()};
      momenta[1] = {pionTrack.px(), pionTrack.py(), pionTrack.pz()};
      const auto v0mass = RecoDecay::m(momenta, masses);

      //// Omega hypohesis -> rejecting Xi, we don't do it in the MC as we can identify the particle with the MC truth
      bool isOmega{std::abs(massXi - constants::physics::MassXiMinus) > 0.005};
      if (cfgSelectOnlyOmegas && !isOmega) {
        continue;
      }

      std::array<bool, 2> fromHF{false, false};
      bool isGoodMatch{false}, isGoodCascade{false};
      int itsTrackPDG{0}, pdgCodeMom{0};
      int64_t mcParticleID{-1};

      if constexpr (TrackType::template contains<aod::McTrackLabels>()) {
        if (protonTrack.mcParticle().has_mothers() && pionTrack.mcParticle().has_mothers() && bachelor.mcParticle().has_mothers()) {
          if (protonTrack.mcParticle().mothersIds()[0] == pionTrack.mcParticle().mothersIds()[0]) {
            const auto v0part = protonTrack.mcParticle().template mothers_first_as<aod::McParticles>();
            if (std::abs(v0part.pdgCode()) == 3122 && v0part.has_mothers()) {
              const auto motherV0 = v0part.template mothers_as<aod::McParticles>()[0];
              if (v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
                if (std::abs(motherV0.pdgCode()) == 3312 || std::abs(motherV0.pdgCode()) == 3334) {
                  isGoodCascade = true;

                  isOmega = (std::abs(motherV0.pdgCode()) == 3334);
                  fromHF = isFromHF(motherV0);
                  mcParticleID = v0part.mothersIds()[0];
                }
              }
            }
          }
        }
      }

      if (cascCpa < cfgMinCosPA) {
        continue;
      }
      if (isOmega) {
        mRegistry.fill(HIST("h_PIDcutsOmega"), 6, massOmega);
      }
      mRegistry.fill(HIST("h_PIDcutsXi"), 6, massXi);

      const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);
      o2::dataformats::DCA motherDCA{-999.f, -999.f}, protonDCA{-999.f, -999.f}, pionDCA{-999.f, -999.f}, bachDCA{-999.f, -999.f};
      o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, protonTrkParCov, mBz, 2.f, matCorr, &protonDCA);
      o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, pionTrkParCov, mBz, 2.f, matCorr, &pionDCA);
      o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, bachTrkParCov, mBz, 2.f, matCorr, &bachDCA);

      float deltaPtITSCascade{-1.e10f}, deltaPtCascade{-1.e10f}, cascITSclsSize{-1.e10f}, matchingChi2{-1.e10f};
      bool hasReassociatedClusters{false}, hasFakeReassociation{false};
      int trackedCascGlobalIndex{-1}, itsTrackGlobalIndex{-1}, cascITSclusters{-1};
      if constexpr (requires { candidate.track(); }) {
        const auto& track = candidate.template track_as<TrackType>();
        const auto& ITStrack = candidate.template itsTrack_as<TrackType>();
        if (cfgRedoPV && ITStrack.isPVContributor()) {
          if (!recalculatePV(collision, tracks, ITStrack.globalIndex(), primaryVertex)) {
            continue;
          }
        }
        cascPVContribs |= ITStrack.isPVContributor() << 0;
        auto trackTrkParCov = getTrackParCov(track);
        o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackTrkParCov, mBz, 2.f, matCorr, &motherDCA);
        hasReassociatedClusters = (track.itsNCls() != ITStrack.itsNCls());
        cascadeLvector.SetCoordinates(track.pt(), track.eta(), track.phi(), 0);
        deltaPtITSCascade = std::hypot(cascadeMomentum[0], cascadeMomentum[1]) - ITStrack.pt();
        deltaPtCascade = std::hypot(cascadeMomentum[0], cascadeMomentum[1]) - track.pt();
        trackedCascGlobalIndex = track.globalIndex();
        itsTrackGlobalIndex = ITStrack.globalIndex();
        cascITSclusters = track.itsNCls();
        cascITSclsSize = candidate.itsClsSize();
        matchingChi2 = candidate.matchingChi2();
        cascadePos = {candidate.decayX(), candidate.decayY(), candidate.decayZ()};
        if constexpr (TrackType::template contains<aod::McTrackLabels>()) {
          isGoodMatch = ((mcParticleID == ITStrack.mcParticleId())) ? true : false;

          if (isGoodMatch) {
            pdgCodeMom = track.mcParticle().has_mothers() ? track.mcParticle().template mothers_as<aod::McParticles>()[0].pdgCode() : 0;
            hasFakeReassociation = track.mcMask() & (1 << 15);
          }
          itsTrackPDG = ITStrack.has_mcParticle() ? ITStrack.mcParticle().pdgCode() : 0;
        }
      } else {
        o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, ntCascadeTrack, mBz, 2.f, matCorr, &motherDCA);
      }
      uint32_t toiMask = 0x0;
      if (mToiMap.count(bc.globalBC())) {
        toiMask = mToiMap[bc.globalBC()];
      }
      candidates.emplace_back(NPCascCandidate{mRunNumber, mcParticleID, trackedCascGlobalIndex, itsTrackGlobalIndex, candidate.collisionId(), matchingChi2, deltaPtITSCascade, deltaPtCascade, cascITSclsSize, hasReassociatedClusters, hasFakeReassociation, isGoodMatch, isGoodCascade, pdgCodeMom, itsTrackPDG, fromHF[0], fromHF[1],
                                              collision.numContrib(), cascPVContribs, collision.collisionTimeRes(), primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                                              cascadeLvector.pt(), cascadeLvector.eta(), cascadeLvector.phi(),
                                              protonTrack.pt(), protonTrack.eta(), pionTrack.pt(), pionTrack.eta(), bachelor.pt(), bachelor.eta(),
                                              motherDCA.getY(), motherDCA.getZ(), protonDCA.getY(), protonDCA.getZ(), pionDCA.getY(), pionDCA.getZ(), bachDCA.getY(), bachDCA.getZ(),
                                              cascCpa, v0Cpa, massXi, massOmega, v0mass,
                                              static_cast<float>(std::hypot(cascadePos[0], cascadePos[1])), static_cast<float>(std::hypot(v0Pos[0], v0Pos[1])), static_cast<float>(std::hypot(cascadePos[0], cascadePos[1], cascadePos[2])), static_cast<float>(std::hypot(v0Pos[0], v0Pos[1], v0Pos[2])),
                                              cascITSclusters, protonTrack.itsNCls(), pionTrack.itsNCls(), bachelor.itsNCls(), protonTrack.tpcNClsFound(), pionTrack.tpcNClsFound(), bachelor.tpcNClsFound(),
                                              protonTrack.tpcNSigmaPr(), pionTrack.tpcNSigmaPi(), bachelor.tpcNSigmaKa(), bachelor.tpcNSigmaPi(),
                                              protonTrack.hasTOF(), pionTrack.hasTOF(), bachelor.hasTOF(),
                                              protonTrack.tofNSigmaPr(), pionTrack.tofNSigmaPi(), bachelor.tofNSigmaKa(), bachelor.tofNSigmaPi(), collision.sel8(), collision.multFT0C(), collision.multFV0A(), collision.multFT0M(), collision.centFT0C(), collision.centFV0A(), collision.centFT0M(), collision.multNTracksGlobal(), toiMask, collision.selection_bit(aod::evsel::kNoSameBunchPileup)});
    }
  }

  template <typename CascadeType>
  void fillDataTable(auto const& candidates)
  {
    for (const auto& c : candidates) {
      getDataTable<CascadeType>()(mRunNumber, c.matchingChi2, c.deltaPtITS, c.deltaPt, c.itsClusSize, c.hasReassociatedCluster,
                                  c.pvContributors, c.cascPVContribs, c.pvTimeResolution, c.pvX, c.pvY, c.pvZ,
                                  c.cascPt, c.cascEta, c.cascPhi,
                                  c.protonPt, c.protonEta, c.pionPt, c.pionEta, c.bachPt, c.bachEta,
                                  c.cascDCAxy, c.cascDCAz, c.protonDCAxy, c.protonDCAz, c.pionDCAxy, c.pionDCAz, c.bachDCAxy, c.bachDCAz,
                                  c.casccosPA, c.v0cosPA,
                                  c.massXi, c.massOmega, c.massV0,
                                  c.cascRadius, c.v0radius, c.cascLength, c.v0length,
                                  c.cascNClusITS, c.protonNClusITS, c.pionNClusITS, c.bachNClusITS, c.protonNClusTPC, c.pionNClusTPC, c.bachNClusTPC,
                                  c.protonTPCNSigma, c.pionTPCNSigma, c.bachKaonTPCNSigma, c.bachPionTPCNSigma,
                                  c.protonHasTOF, c.pionHasTOF, c.bachHasTOF,
                                  c.protonTOFNSigma, c.pionTOFNSigma, c.bachKaonTOFNSigma, c.bachPionTOFNSigma,
                                  c.sel8, c.multFT0C, c.multFV0A, c.multFT0M, c.centFT0C, c.centFV0A, c.centFT0M, c.multNTracksGlobal, c.toiMask, c.noSameBunchPileup);
    }
  }

  template <typename CascadeType>
  void fillMCtable(auto const& mcParticles, auto const& collisions, auto const& candidates)
  {
    for (size_t i = 0; i < candidates.size(); ++i) {
      auto& c = candidates[i];
      if (c.mcParticleId < 0) {
        continue;
      }
      auto particle = mcParticles.iteratorAt(c.mcParticleId);
      int motherDecayDaughters{0};
      if (c.isFromBeauty || c.isFromCharm) {
        auto mom = particle.template mothers_as<aod::McParticles>()[0];
        auto daughters = mom.template daughters_as<aod::McParticles>();
        motherDecayDaughters = daughters.size();
        for (const auto& d : daughters) {
          if (std::abs(d.pdgCode()) == 11 || std::abs(d.pdgCode()) == 13) {
            motherDecayDaughters *= -1;
            break;
          }
        }
      }
      auto mcCollision = particle.template mcCollision_as<aod::McCollisions>();
      auto recCollision = collisions.iteratorAt(c.collisionID);

      getMCtable<CascadeType>()(mRunNumber, c.matchingChi2, c.deltaPtITS, c.deltaPt, c.itsClusSize, c.hasReassociatedCluster, c.isGoodMatch, c.isGoodCascade, c.pdgCodeMom, c.pdgCodeITStrack, c.isFromBeauty, c.isFromCharm,
                                c.pvContributors, c.cascPVContribs, c.pvTimeResolution, c.pvX, c.pvY, c.pvZ, c.cascPt, c.cascEta, c.cascPhi,
                                c.protonPt, c.protonEta, c.pionPt, c.pionEta, c.bachPt, c.bachEta,
                                c.cascDCAxy, c.cascDCAz, c.protonDCAxy, c.protonDCAz, c.pionDCAxy, c.pionDCAz, c.bachDCAxy, c.bachDCAz,
                                c.casccosPA, c.v0cosPA, c.massXi, c.massOmega, c.massV0, c.cascRadius, c.v0radius, c.cascLength, c.v0length,
                                c.cascNClusITS, c.protonNClusITS, c.pionNClusITS, c.bachNClusITS, c.protonNClusTPC, c.pionNClusTPC, c.bachNClusTPC, c.protonTPCNSigma,
                                c.pionTPCNSigma, c.bachKaonTPCNSigma, c.bachPionTPCNSigma, c.protonHasTOF, c.pionHasTOF, c.bachHasTOF,
                                c.protonTOFNSigma, c.pionTOFNSigma, c.bachKaonTOFNSigma, c.bachPionTOFNSigma,
                                c.sel8, c.multFT0C, c.multFV0A, c.multFT0M, c.centFT0C, c.centFV0A, c.centFT0M,
                                particle.pt(), particle.eta(), particle.phi(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
                                particle.pdgCode(), mcCollision.posX() - particle.vx(), mcCollision.posY() - particle.vy(),
                                mcCollision.posZ() - particle.vz(), mcCollision.globalIndex() == recCollision.mcCollisionId(), c.hasFakeReassociation, motherDecayDaughters, c.multNTracksGlobal, c.toiMask, c.noSameBunchPileup);
    }
  }

  template <typename CascadeType>
  auto& getMCtable()
  {
    if constexpr (std::is_same_v<CascadeType, aod::Cascades>) {
      return NPCTableMCNT;
    } else {
      return NPCTableMC;
    }
  }

  template <typename CascadeType>
  auto& getDataTable()
  {
    if constexpr (std::is_same_v<CascadeType, aod::Cascades>) {
      return NPCTableNT;
    } else {
      return NPCTable;
    }
  }

  void processTrackedCascadesMC(CollisionCandidatesRun3MC const& collisions,
                                aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& /*cascades*/,
                                aod::V0s const& /*v0s*/, TracksExtMC const& tracks,
                                aod::McParticles const& mcParticles, aod::McCollisions const&, aod::BCsWithTimestamps const&)
  {
    fillCandidatesVector<TracksExtMC>(collisions, tracks, trackedCascades, gCandidates);
    fillMCtable<aod::AssignedTrackedCascades>(mcParticles, collisions, gCandidates);
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processTrackedCascadesMC, "process cascades from strangeness tracking: MC analysis", true);

  void processCascadesMC(CollisionCandidatesRun3MC const& collisions, aod::Cascades const& cascades,
                         aod::V0s const& /*v0s*/, TracksExtMC const& tracks,
                         aod::McParticles const& mcParticles, aod::McCollisions const&, aod::BCsWithTimestamps const&)
  {
    fillCandidatesVector<TracksExtMC>(collisions, tracks, cascades, gCandidatesNT);
    fillMCtable<aod::Cascades>(mcParticles, collisions, gCandidatesNT);
    fillMultHistos<CollisionCandidatesRun3MC>(collisions);
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processCascadesMC, "process cascades: MC analysis", false);

  void processGenParticles(aod::McParticles const& mcParticles)
  {
    for (const auto& p : mcParticles) {
      auto absCode = std::abs(p.pdgCode());
      if (absCode != 3312 && absCode != 3334) {
        continue;
      }
      auto fromHF = isFromHF(p);
      int pdgCodeMom = p.has_mothers() ? p.template mothers_as<aod::McParticles>()[0].pdgCode() : 0;
      auto mcCollision = p.template mcCollision_as<aod::McCollisions>();
      int motherDecayDaughters{0};
      if (fromHF[0] || fromHF[1]) {
        auto mom = p.template mothers_as<aod::McParticles>()[0];
        auto daughters = mom.template daughters_as<aod::McParticles>();
        motherDecayDaughters = daughters.size();
        for (const auto& d : daughters) {
          if (std::abs(d.pdgCode()) == 11 || std::abs(d.pdgCode()) == 13) {
            motherDecayDaughters *= -1;
            break;
          }
        }
      }
      NPCTableGen(p.pt(), p.eta(), p.phi(), p.pdgCode(), pdgCodeMom, mcCollision.posX() - p.vx(), mcCollision.posY() - p.vy(), mcCollision.posZ() - p.vz(), fromHF[0], fromHF[1], motherDecayDaughters);
    }
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processGenParticles, "process gen cascades: MC analysis", false);

  void processTrackedCascadesData(CollisionCandidatesRun3 const& collisions,
                                  aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& /*cascades*/,
                                  aod::V0s const& /*v0s*/, TracksExtData const& tracks,
                                  aod::BCsWithTimestamps const&)
  {
    mProcessCounter[0]++;
    zorroAccounting(collisions);
    fillCandidatesVector<TracksExtData>(collisions, tracks, trackedCascades, gCandidates);
    fillDataTable<aod::AssignedTrackedCascades>(gCandidates);
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processTrackedCascadesData, "process cascades from strangeness tracking: Data analysis", false);

  void processCascadesData(CollisionCandidatesRun3 const& collisions, aod::Cascades const& cascades,
                           aod::V0s const& /*v0s*/, TracksExtData const& tracks,
                           aod::BCsWithTimestamps const&)
  {
    mProcessCounter[1]++;
    zorroAccounting(collisions);
    fillCandidatesVector<TracksExtData>(collisions, tracks, cascades, gCandidatesNT);
    fillDataTable<aod::Cascades>(gCandidatesNT);
    fillMultHistos<CollisionCandidatesRun3>(collisions);
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processCascadesData, "process cascades: Data analysis", false);

  int getMCMult(aod::McParticles const& mcParticles, int mcCollId, std::vector<float>& ptMCvec)
  {
    int mult = 0;
    for (auto const& mcp : mcParticles) {
      if (mcp.mcCollisionId() == mcCollId) {
        // multiplicity definition:
        bool accept = mcp.isPhysicalPrimary();
        accept = accept && (mcp.eta() < 0.5) && (mcp.eta() > -0.5);
        int q = 0;
        auto pdgEntry = TDatabasePDG::Instance()->GetParticle(mcp.pdgCode());
        if (pdgEntry) {
          q = int(std::round(pdgEntry->Charge() / 3.0));
        } else {
          // LOG(warn) << "No pdg assuming neutral";
        }
        accept = accept && (q != 0);
        if (accept) {
          ++mult;
          ptMCvec.push_back(mcp.pt());
        }
      }
    }
    return mult;
  }
  void processdNdetaMC(CollisionsWithLabel const& colls, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles, TracksWithLabel const& tracks)
  {
    // std::cout << "ProcNegMC" << std::endl;
    // Map: collision index -> reco multiplicity
    std::vector<int> recoMult(colls.size(), 0);
    for (auto const& trk : tracks) {
      int collId = trk.collisionId();
      // Here you can impose same track cuts as for MC (eta, primaries, etc.)
      if (std::abs(trk.eta()) > 0.5f) {
        continue;
      }
      ++recoMult[collId];
    }
    std::vector<int> isReco(mcParticles.size(), 0);
    std::vector<int> isRecoMult(mcParticles.size(), 0);
    std::vector<int> mcReconstructed(mcCollisions.size(), 0);
    // std::cout << " reco cols with mc:" << colls.size() << " tracks:" << tracks.size() << " mccols:" << mcCollisions.size() << "mcParticles:" << mcParticles.size() << std::endl;
    for (auto const& col : colls) {
      int mcCollId = col.mcCollisionId(); // col.template mcCollision_as<aod::McCollisions>();
      int collId = col.globalIndex();
      // auto mc = col.mcCollision();
      // int mcId = mc.globalIndex();
      // std::cout << "globalIndex:" << mcId << " colID:" << mcCollId << std::endl;
      std::vector<float> mcptvec;
      float mult = getMCMult(mcParticles, mcCollId, mcptvec);
      mcReconstructed[mcCollId] = 1;
      for (auto const& trk : tracks) {
        int mcPid = trk.mcParticleId(); // depends on your label table
        if (mcPid >= 0 && mcPid < (int)mcParticles.size()) {
          isReco[mcPid] = 1;
          isRecoMult[mcPid] = mult;
        } else {
          continue;
        }
        if (trk.collisionId() != collId) {
          continue; // different event
        }
        auto mcPar = mcParticles.rawIteratorAt(mcPid);
        // Apply same acceptance as in MC multiplicity
        if (!mcPar.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(mcPar.eta()) > 0.5f) {
          continue;
        }
        int q = 0;
        auto pdgEntry = TDatabasePDG::Instance()->GetParticle(mcPar.pdgCode());
        if (pdgEntry) {
          q = int(std::round(pdgEntry->Charge() / 3.0));
        }
        if (q == 0) {
          continue;
        }
        // float multReco = recoMult[collId];
        float multReco = col.multNTracksGlobal();
        float ptReco = trk.pt();
        float ptMC = mcPar.pt();
        mRegistrydNdeta.fill(HIST("hdNdetaRM/hdNdetaRM"), mult, multReco, ptMC, ptReco);
      }

      // mRegistry.fill(HIST("hNTracksMCVsTracksReco"), mult, col.multNTracksGlobal());
    }
    // count mc particles with no reco tracks
    for (auto const& mcp : mcParticles) {
      int mcPidG = mcp.globalIndex();
      // std::cout << "mcPidG:" << mcPidG << std::endl;
      if (!isReco[mcPidG]) {
        mRegistrydNdeta.fill(HIST("hdNdetaRM/hdNdetaRMNotInRecoTrk"), isRecoMult[mcPidG], mcp.pt());
      }
    }
    // count unreconstructed mc collisions
    for (auto const& mc : mcCollisions) {
      int gindex = mc.globalIndex();
      // std::cout << "mc globalIndex:" << gindex << std::endl;
      if (!mcReconstructed[gindex]) {
        std::vector<float> mcptvec;
        int mult = getMCMult(mcParticles, gindex, mcptvec);
        // std::cout << "===> unreconstructed:" << mult << std::endl;
        for (auto const& pt : mcptvec) {
          mRegistrydNdeta.fill(HIST("hdNdetaRM/hdNdetaRMNotInRecoCol"), mult, pt);
        }
      }
    }
  }
  PROCESS_SWITCH(NonPromptCascadeTask, processdNdetaMC, "process mc dN/deta", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NonPromptCascadeTask>(cfgc)};
}
