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
// Search for hypernuclei kink decay topology
// ==============================================================================

#include <array>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"

#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFHypernucleiTables.h"

static constexpr int POS = 0, NEG = 1;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using VBracket = o2::math_utils::Bracket<int>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCFullTr, aod::pidTOFFullTr>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms>;

namespace
{
constexpr std::array<float, 7> LayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"Triton"};
std::shared_ptr<TH1> hEvents;
std::shared_ptr<TH1> hZvtx;
std::shared_ptr<TH1> hCentFT0A;
std::shared_ptr<TH1> hCentFT0C;
std::shared_ptr<TH1> hCentFT0M;
std::shared_ptr<TH1> hCentFV0A;
std::shared_ptr<TH2> hNsigmaTritSel;
std::shared_ptr<TH2> hDeDxTritSel;
std::shared_ptr<TH2> hDeDxTot;
std::shared_ptr<TH1> hIsMatterGen;
} // namespace

struct kinkCandidate {

  float recoPtHyp() const { return std::hypot(momHyp[0], momHyp[1]); }
  float recoPhiHyp() const { return std::atan2(momHyp[1], momHyp[0]); }
  float recoEtaHyp() const { return std::asinh(momHyp[2] / recoPtHyp()); }

  float recoPtTrit() const { return std::hypot(momTrit[0], momTrit[1]); }
  float recoPhiTrit() const { return std::atan2(momTrit[1], momTrit[0]); }
  float recoEtaTrit() const { return std::asinh(momTrit[2] / recoPtTrit()); }

  float genPt() const { return std::hypot(gMomHyp[0], gMomHyp[1]); }
  float genPtTrit() const { return std::hypot(gMomTrit[0], gMomTrit[1]); }
  float genPhi() const { return std::atan2(gMomHyp[1], gMomHyp[0]); }
  float genEta() const { return std::asinh(gMomHyp[2] / genPt()); }

  int hyperTrackID;
  int tritTrackID;
  bool isMatter = false;

  std::array<float, 3> momHyp = {-999, -999, -999};
  std::array<float, 3> momTrit = {-999, -999, -999};
  std::array<float, 3> primVtx = {-999, -999, -999};
  std::array<float, 3> decVtx = {-999, -999, -999};

  float dcaKinkTopo = -999;
  float nSigmaTPCTrit = -999;
  float nSigmaTOFTrit = -999;
  float tritDCAXY = -999;
  float hypDCAXY = -999;
  float kinkAngle = -999;

  float momTritTPC = -999.f;
  uint16_t tpcSignalTrit = 0u;
  uint8_t nTPCClustersTrit = 0u;
  uint8_t trackingPIDTriton = 0u; // flags for triton PID in tracking

  uint32_t clusterSizeITSHyp = 0u;
  uint32_t clusterSizeITSTrit = 0u;

  std::array<float, 3> gMomHyp;
  std::array<float, 3> gMomTrit;
  std::array<float, 3> gDecVtx;

  bool isSignal = false;        // true MC signal
  bool isReco = false;          // true if the candidate is actually reconstructed
  float itsPt = -999.f;         // pt of the ITS hypertrack even when the topology is not reconstructed, tagged with the MC truth
  uint16_t mcMask = false;      // to study fake its tracks
  bool survEvSelection = false; // true if the corresponding event passed the event selection
  int pdgCode = 0;              // pdg code of the mother particle
};

struct TrackCand {
  int Idxtr;
  int mcParticleIdx = -1;
  VBracket vBracket{};
};

struct hyperKinkRecoTask {

  Produces<aod::DataHypKinkCands> outputDataTable;
  Produces<aod::MCHypKinkCands> outputMCTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Selection criteria
  Configurable<float> invMassLow{"invMassLow", 2.9, "Lower limit for the invariant mass"};
  Configurable<float> invMassHigh{"invMassHigh", 15., "Upper limit for the invariant mass"};
  Configurable<float> maxDCAHypToPV{"maxDCAHypToPV", 0.1, "Max DCA of the hypertriton to the PV"};
  Configurable<float> minDCATritToPV{"minDCATritToPV", 0., "Min DCA of the triton to the PV"};
  Configurable<float> minPtHyp{"minPtHyp", 0.5, "Minimum pT of the hypercandidate"};
  Configurable<float> maxZDiff{"maxZDiff", 20., "Max z difference between the kink daughter and the hypertriton"};
  Configurable<float> maxPhiDiff{"maxPhiDiff", 100, "Max phi difference between the kink daughter and the hypertriton"};
  Configurable<float> timeMarginNS{"timeMarginNS", 600, "Additional time res tolerance in ns"};
  Configurable<float> etaMax{"eta", 1., "eta daughter"};
  Configurable<float> nSigmaTPCCutTrit{"nSigmaTPCTrit", 5, "triton dEdx cut (n sigma)"};
  Configurable<float> nSigmaTOFCutTrit{"nSigmaTOFTrit", 5, "triton TOF cut (n sigma)"};
  Configurable<float> nTPCClusMinTrit{"nTPCClusMinTrit", 80, "triton NTPC clusters cut"};
  Configurable<bool> alwaysAskTOF{"alwaysAskTOF", false, "If true, ask for TOF signal"};
  Configurable<bool> mcSignalOnly{"mcSignalOnly", true, "If true, save only signal in MC"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;
  o2::base::MatLayerCylSet* lut = nullptr;

  // constants
  float tritonMass = o2::constants::physics::MassTriton;
  float pi0Mass = o2::constants::physics::MassPi0;
  float radToDeg = 180. / M_PI;
  Configurable<int> hyperPdg{"hyperPDG", 1010010030, "PDG code of the hyper-mother"};
  Configurable<int> tritDauPdg{"kinkDauPDG", 1000010030, "PDG code of the kink daughter"};

  // bethe bloch parameters
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for Triton"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  // PDG codes

  // histogram axes
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis dedxBins{"dedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {200, -5.f, 5.f}, "Binning for n sigma"};
  ConfigurableAxis zVtxBins{"zVtxBins", {100, -20.f, 20.f}, "Binning for n sigma"};
  ConfigurableAxis centBins{"centBins", {100, 0.f, 100.f}, "Binning for centrality"};

  // std vector of candidates
  std::vector<kinkCandidate> mKinkCandidates;
  // vector to keep track of MC mothers already filled
  std::vector<unsigned int> filledMothers;
  // vector to keep track of the collisions passing the event selection in the MC
  std::vector<bool> isGoodCollision;
  HistogramRegistry qaRegistry{"QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float d_bz;
  std::array<float, 6> mBBparamsTrit;

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
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    int mat{static_cast<int>(cfgMaterialCorrection)};
    fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    const AxisSpec rigidityAxis{rigidityBins, "#it{p}^{TPC}/#it{z}"};
    const AxisSpec dedxAxis{dedxBins, "d#it{E}/d#it{x}"};
    const AxisSpec nSigma3HeAxis{nSigmaBins, "n_{#sigma}({}^{3}He)"};
    const AxisSpec zVtxAxis{zVtxBins, "z_{vtx} (cm)"};
    const AxisSpec centAxis{centBins, "Centrality"};

    hNsigmaTritSel = qaRegistry.add<TH2>("hNsigmaTritSel", "; p_{TPC}/z (GeV/#it{c}); n_{#sigma} ({}^{3}H)", HistType::kTH2F, {rigidityAxis, nSigma3HeAxis});
    hDeDxTritSel = qaRegistry.add<TH2>("hDeDxTritSel", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dedxAxis});
    hDeDxTot = qaRegistry.add<TH2>("hDeDxTot", ";p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dedxAxis});
    hEvents = qaRegistry.add<TH1>("hEvents", ";Events; ", HistType::kTH1D, {{3, -0.5, 2.5}});
    hEvents->GetXaxis()->SetBinLabel(1, "All");
    hEvents->GetXaxis()->SetBinLabel(2, "sel8");
    hEvents->GetXaxis()->SetBinLabel(3, "z vtx");

    if (doprocessMC) {
      hIsMatterGen = qaRegistry.add<TH1>("hIsMatterGen", ";; ", HistType::kTH1D, {{2, -0.5, 1.5}});
      hIsMatterGen->GetXaxis()->SetBinLabel(1, "Matter");
      hIsMatterGen->GetXaxis()->SetBinLabel(2, "Antimatter");
    }
    hZvtx = qaRegistry.add<TH1>("hZvtx", ";z_{vtx} (cm); ", HistType::kTH1D, {{100, -20, 20}});
    hCentFT0M = qaRegistry.add<TH1>("hCentFT0M", ";Centrality; ", HistType::kTH1D, {{100, 0, 100}});
  }

  float angleCutFunction(float x)
  {
    float par1 = 2.99131; // hypertriton mass
    float par2 = 0.07;    // optimized by mdiotti
    float par3 = TMath::Pi();
    return par1 * (par2 / (sqrt((x * x) * (1 - (par2 * par2)) - ((par1 * par1) * (par2 * par2))))) * (180. / par3) + 1;
  }

  template <typename T>
  float computeNSigmaTrit(const T& candidate)
  {
    float nSigmaTrit = -100.f;
    if (mBBparamsTrit[5] > 0) {
      float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() * 2 / constants::physics::MassHelium3), mBBparamsTrit[0], mBBparamsTrit[1], mBBparamsTrit[2], mBBparamsTrit[3], mBBparamsTrit[4]);
      double resoTPC{expTPCSignal * mBBparamsTrit[5]};
      nSigmaTrit = static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
    } else {
      nSigmaTrit = candidate.tpcNSigmaTr();
    }
    return nSigmaTrit;
  }

  template <typename T>
  bool selectHyperTrack(const T& candidate)
  {
    if (candidate.hasITS() && !candidate.hasTPC() && !candidate.hasTOF() && candidate.itsNCls() < 6 &&
        candidate.itsNClsInnerBarrel() == 3 && candidate.itsChi2NCl() < 36 && candidate.pt() > minPtHyp) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectTritTrack(const T& candidate)
  {
    if (!candidate.hasTPC() || !candidate.hasITS()) {
      return false;
    }

    if (alwaysAskTOF && !candidate.hasTOF()) {
      return false;
    }

    float nSigmaTrit = computeNSigmaTrit(candidate);
    float nSigmaTOFTrit = candidate.tofNSigmaTr();

    bool isGoodTPCCand = false;
    if (candidate.itsNClsInnerBarrel() == 0 && candidate.itsNCls() < 4 &&
        candidate.tpcNClsCrossedRows() >= 70 && candidate.tpcChi2NCl() < 4.f &&
        candidate.tpcNClsCrossedRows() > 0.8 * candidate.tpcNClsFindable() && candidate.tpcNClsFound() > 80 && abs(nSigmaTrit) < nSigmaTPCCutTrit) {
      isGoodTPCCand = true;
    }

    if (!isGoodTPCCand) {
      return false;
    }

    if (candidate.hasTOF() && abs(nSigmaTOFTrit) > nSigmaTOFCutTrit) {
      return false;
    }

    hNsigmaTritSel->Fill(candidate.pt(), nSigmaTrit);
    hDeDxTritSel->Fill(candidate.tpcInnerParam(), candidate.tpcSignal());

    return true;
  }

  std::array<std::vector<TrackCand>, 4> makeTracksPool(CollisionsFull const& collisions, TracksFull const& tracks, o2::aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const&)
  {
    std::unordered_map<int, std::pair<int, int>> tmap;
    TrackCand trForpool;

    std::array<std::vector<TrackCand>, 4> pools; // pools of positive and negative seeds sorted in min VtxID
    std::vector<uint8_t> selected(tracks.size(), 0u);
    std::vector<uint64_t> globalBCvector;

    int index{0};
    for (const auto& track : tracks) {
      if (track.has_collision()) {
        if (track.collision_as<CollisionsFull>().has_bc()) {
          globalBCvector.push_back(track.collision_as<CollisionsFull>().bc_as<aod::BCsWithTimestamps>().globalBC());
        }
      } else {
        for (const auto& ambTrack : ambiTracks) {
          if (ambTrack.trackId() == track.globalIndex()) {
            if (!ambTrack.has_bc() || ambTrack.bc_as<aod::BCsWithTimestamps>().size() == 0) {
              globalBCvector.push_back(-1);
              break;
            }
            globalBCvector.push_back(ambTrack.bc_as<aod::BCsWithTimestamps>().begin().globalBC());
            break;
          }
        }
      }

      if (std::abs(track.eta()) < 0.8) {
        if (selectHyperTrack(track)) {
          selected[index] = 1;
        } else if (selectTritTrack(track)) {
          selected[index] = 2;
        }
      }

      index++;
    }
    constexpr auto bOffsetMax = 241; // 6 mus (ITS)

    for (const auto& collision : collisions) {

      hEvents->Fill(0);
      if (!collision.sel8())
        continue;
      hEvents->Fill(1);
      if (std::abs(collision.posZ()) > 10.)
        continue;
      hEvents->Fill(2);
      hZvtx->Fill(collision.posZ());
      hCentFT0M->Fill(collision.centFT0M());

      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc_as<aod::BCsWithTimestamps>().globalBC();
      const auto collIdx = collision.globalIndex();

      index = -1;
      for (const auto& track : tracks) {
        index++;
        if (!selected[index] || !track.has_collision())
          continue;
        const int64_t bcOffset = (int64_t)globalBCvector[track.globalIndex()] - (int64_t)collBC;
        if (std::abs(bcOffset) > bOffsetMax) {
          continue;
        }

        float trackTime{0.};
        float trackTimeRes{0.};
        if (track.isPVContributor()) {
          trackTime = track.collision_as<CollisionsFull>().collisionTime(); // if PV contributor, we assume the time to be the one of the collision
          trackTimeRes = constants::lhc::LHCBunchSpacingNS;                 // 1 BC
        } else {
          trackTime = track.trackTime();
          trackTimeRes = track.trackTimeRes();
        }

        const float deltaTime = trackTime - collTime + bcOffset * constants::lhc::LHCBunchSpacingNS;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;

        float thresholdTime = 0.;
        if (track.isPVContributor()) {
          thresholdTime = trackTimeRes;
        } else if (TESTBIT(track.flags(), o2::aod::track::TrackTimeResIsRange)) {
          thresholdTime = std::sqrt(sigmaTimeRes2); // + timeMargin;
          thresholdTime += timeMarginNS;

        } else {
          // thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2); // + timeMargin;
          thresholdTime = 4. * std::sqrt(sigmaTimeRes2); // + timeMargin;
          thresholdTime += timeMarginNS;
        }

        if (std::abs(deltaTime) > thresholdTime) {
          continue;
        }

        const auto& tref = tmap.find(track.globalIndex());
        if (tref != tmap.end()) {
          LOG(debug) << "Track: " << track.globalIndex() << " already processed with other vertex";
          pools[tref->second.second][tref->second.first].vBracket.setMax(static_cast<int>(collIdx)); // this track was already processed with other vertex, account the latter
          continue;
        }

        int poolIndex = (selected[index] - 1) * 2 + (track.sign() < 0); /// first the two mothers then the two daughters (mom pos 0, mom neg 1, dau pos 2, dau neg 3)
        trForpool.Idxtr = track.globalIndex();
        trForpool.vBracket = {static_cast<int>(collIdx), static_cast<int>(collIdx)};
        pools[poolIndex].emplace_back(trForpool);
        tmap[track.globalIndex()] = {pools[poolIndex].size() - 1, poolIndex};

      } // track Mother loop
    }   // collision loop

    return pools;
  }

  void fillCandidateData(CollisionsFull const& collisions, TracksFull const& tracks, o2::aod::AmbiguousTracks const&, aod::BCsWithTimestamps const&, std::array<std::vector<TrackCand>, 4> trackPool)
  {
    gsl::span<std::vector<TrackCand>> hyperPool{trackPool.data(), 2};
    gsl::span<std::vector<TrackCand>> tritPool{trackPool.data() + 2, 2};

    std::array<std::vector<int>, 2> mVtxTritTrack{}; // 1st pos. and neg. track of the kink pool for each vertex
    for (int i = 0; i < 2; i++) {
      mVtxTritTrack[i].clear();
      mVtxTritTrack[i].resize(collisions.size(), -1);
    }

    for (int pn = 0; pn < 2; pn++) {
      auto& vtxFirstT = mVtxTritTrack[pn];
      const auto& signedTritPool = tritPool[pn];
      for (unsigned i = 0; i < signedTritPool.size(); i++) {
        const auto& t = signedTritPool[i];
        const auto& track = tracks.iteratorAt(t.Idxtr);
        LOG(debug) << "Track with index: " << track.globalIndex() << " min bracket: " << t.vBracket.getMin() << " max bracket: " << t.vBracket.getMax() << " and sign: " << track.sign();
        for (int j{t.vBracket.getMin()}; j <= t.vBracket.getMax(); ++j) {
          if (vtxFirstT[j] == -1) {
            vtxFirstT[j] = i;
          }
        }
      }

      auto& signedHyperPool = hyperPool[pn];
      for (unsigned itp = 0; itp < signedHyperPool.size(); itp++) {
        auto& seedHyper = signedHyperPool[itp];
        LOG(debug) << "Processing hyperseed with global index " << seedHyper.Idxtr << ", bracket: " << seedHyper.vBracket.getMin() << " - " << seedHyper.vBracket.getMax() << " and sign: " << tracks.iteratorAt(seedHyper.Idxtr).sign();
        int firstKinkID = -1;
        for (int j{seedHyper.vBracket.getMin()}; j <= seedHyper.vBracket.getMax(); ++j) {
          LOG(debug) << "Checking vtxFirstT at position " << j << " with value " << vtxFirstT[j];
          if (vtxFirstT[j] != -1) {
            firstKinkID = vtxFirstT[j];
            break;
          }
        }

        if (firstKinkID < 0) {
          continue;
        }
        const auto& trackHyper = tracks.iteratorAt(seedHyper.Idxtr);
        if (!trackHyper.has_collision())
          continue;

        auto const& collision = trackHyper.collision_as<CollisionsFull>();
        o2::dataformats::VertexBase primaryVertex;
        primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
        primaryVertex.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

        o2::track::TrackParCov trackParCovHyper = getTrackParCov(trackHyper);
        o2::base::Propagator::Instance()->PropagateToXBxByBz(trackParCovHyper, LayerRadii[trackHyper.itsNCls() - 1]);

        o2::track::TrackParCov trackParCovHyperPV = getTrackParCov(trackHyper);
        gpu::gpustd::array<float, 2> dcaInfoHyp;
        o2::base::Propagator::Instance()->propagateToDCABxByBz({primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, trackParCovHyperPV, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoHyp);

        if (abs(dcaInfoHyp[0]) > maxDCAHypToPV) {
          continue;
        }

        kinkCandidate kinkCand;
        kinkCand.primVtx = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};
        // now we search for the kink daughter

        for (unsigned itn = firstKinkID; itn < signedTritPool.size(); itn++) {
          auto& seedTrit = signedTritPool[itn];

          if (seedTrit.vBracket.getMin() > seedHyper.vBracket.getMax()) {
            break;
          }

          if (seedTrit.vBracket.isOutside(seedHyper.vBracket)) {
            LOG(debug) << "Brackets do not match";
            continue;
          }

          const auto& trackTrit = tracks.iteratorAt(seedTrit.Idxtr);
          LOG(debug) << "Trying to pair a triton track with global index " << seedTrit.Idxtr;
          o2::track::TrackParCov trackParCovTrit = getTrackParCov(trackTrit);

          // check if the kink daughter is close to the hypertriton
          if (std::abs(trackParCovHyper.getZ() - trackParCovTrit.getZ()) > maxZDiff) {
            continue;
          }
          if ((std::abs(trackParCovHyper.getPhi() - trackParCovTrit.getPhi()) * radToDeg) > maxPhiDiff) {
            continue;
          }

          // propagate to PV
          gpu::gpustd::array<float, 2> dcaInfoTrit;
          o2::base::Propagator::Instance()->propagateToDCABxByBz({primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, trackParCovTrit, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoTrit);
          if (abs(dcaInfoTrit[0]) < minDCATritToPV) {
            continue;
          }

          int nCand = 0;
          try {
            nCand = fitter.process(trackParCovHyper, trackParCovTrit);
          } catch (...) {
            LOG(error) << "Exception caught in DCA fitter process call!";
            continue;
          }
          if (nCand == 0) {
            continue;
          }

          if (!fitter.propagateTracksToVertex()) {
            continue;
          }

          auto propHyperTrack = fitter.getTrack(0);
          auto propTritTrack = fitter.getTrack(1);

          kinkCand.decVtx = fitter.getPCACandidatePos();

          // cut on decay radius to 17 cm
          if (kinkCand.decVtx[0] * kinkCand.decVtx[0] + kinkCand.decVtx[1] * kinkCand.decVtx[1] < 17 * 17) {
            continue;
          }

          propHyperTrack.getPxPyPzGlo(kinkCand.momHyp);
          propTritTrack.getPxPyPzGlo(kinkCand.momTrit);
          float pHyp = propHyperTrack.getP();
          float pTrit = propTritTrack.getP();
          float spKink = kinkCand.momHyp[0] * kinkCand.momTrit[0] + kinkCand.momHyp[1] * kinkCand.momTrit[1] + kinkCand.momHyp[2] * kinkCand.momTrit[2];
          kinkCand.kinkAngle = std::acos(spKink / (pHyp * pTrit));
          float angleCut = angleCutFunction(pHyp);
          if (kinkCand.kinkAngle * radToDeg > angleCut) {
            continue;
          }

          std::array<float, 3> pi0mom{0.f, 0.f, 0.f};
          for (int i = 0; i < 3; i++) {
            pi0mom[i] = kinkCand.momHyp[i] - kinkCand.momTrit[i];
          }
          float pi0E = std::sqrt(pi0mom[0] * pi0mom[0] + pi0mom[1] * pi0mom[1] + pi0mom[2] * pi0mom[2] + pi0Mass * pi0Mass);
          float tritE = std::sqrt(pTrit * pTrit + tritonMass * tritonMass);
          float invMass = std::sqrt((pi0E + tritE) * (pi0E + tritE) - pHyp * pHyp);

          if (invMass < invMassLow || invMass > invMassHigh) {
            continue;
          }

          kinkCand.hyperTrackID = seedHyper.Idxtr;
          kinkCand.hypDCAXY = dcaInfoHyp[0];
          kinkCand.clusterSizeITSHyp = trackHyper.itsClusterSizes();
          kinkCand.isMatter = trackHyper.sign() > 0;
          kinkCand.tritTrackID = seedTrit.Idxtr;
          kinkCand.tritDCAXY = dcaInfoTrit[0];
          kinkCand.clusterSizeITSTrit = trackTrit.itsClusterSizes();
          kinkCand.nSigmaTPCTrit = computeNSigmaTrit(trackTrit);
          kinkCand.nTPCClustersTrit = trackTrit.tpcNClsFound();
          kinkCand.tpcSignalTrit = trackTrit.tpcSignal();
          kinkCand.momTritTPC = trackTrit.tpcInnerParam();
          kinkCand.dcaKinkTopo = std::sqrt(fitter.getChi2AtPCACandidate());
          kinkCand.trackingPIDTriton = trackTrit.pidForTracking();
          kinkCand.isReco = true;
          mKinkCandidates.push_back(kinkCand);
        }
      }
    }
  }

  void fillMCinfo(aod::McTrackLabels const& trackLabels, aod::McParticles const&)
  {

    for (auto& kinkCand : mKinkCandidates) {
      auto mcLabHyper = trackLabels.rawIteratorAt(kinkCand.hyperTrackID);
      auto mcLabTrit = trackLabels.rawIteratorAt(kinkCand.tritTrackID);
      if (mcLabHyper.has_mcParticle() && mcLabTrit.has_mcParticle()) {
        auto mcTrackHyper = mcLabHyper.mcParticle_as<aod::McParticles>();
        auto mcTrackTrit = mcLabTrit.mcParticle_as<aod::McParticles>();

        if (abs(mcTrackHyper.pdgCode()) != hyperPdg || abs(mcTrackTrit.pdgCode()) != tritDauPdg) {
          continue;
        }
        auto tritIdx = mcTrackTrit.globalIndex();
        kinkCand.isSignal = false;
        for (auto& dauMCTracks : mcTrackHyper.daughters_as<aod::McParticles>()) {
          if (abs(dauMCTracks.pdgCode()) == tritDauPdg) {
            if (dauMCTracks.globalIndex() == tritIdx) {
              kinkCand.isSignal = true;
              break;
            }
          }
        }
        auto primVtx = array{mcTrackHyper.vx(), mcTrackHyper.vy(), mcTrackHyper.vz()};
        auto secVtx = array{mcTrackTrit.vx(), mcTrackTrit.vy(), mcTrackTrit.vz()};
        for (int i = 0; i < 3; i++) {
          kinkCand.gDecVtx[i] = secVtx[i] - primVtx[i];
        }
        kinkCand.pdgCode = mcTrackHyper.pdgCode();
        kinkCand.mcMask = mcLabHyper.mcMask();
        filledMothers.push_back(mcTrackHyper.globalIndex());
      }
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
    if (!pidPath.value.empty()) {
      auto tritpid = ccdb->getForTimeStamp<std::array<float, 6>>(pidPath.value + "_Trit", run3grp_timestamp);
      std::copy(tritpid->begin(), tritpid->end(), mBBparamsTrit.begin());
    } else {
      for (int i = 0; i < 5; i++) {
        mBBparamsTrit[i] = cfgBetheBlochParams->get("Triton", Form("p%i", i));
      }
      mBBparamsTrit[5] = cfgBetheBlochParams->get("Triton", "resolution");
    }
    fitter.setBz(d_bz);
    mRunNumber = bc.runNumber();

    o2::base::Propagator::Instance()->setMatLUT(lut);
  }

  void processData(CollisionsFull const& collisions, TracksFull const& tracks, o2::aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcWtmp)
  {

    if (mBBparamsTrit[5] < 0) {
      LOG(info) << "Bethe-Bloch parameters for Triton not set, using default nSigma";
    }

    mKinkCandidates.clear();

    auto firstcollision = collisions.begin();
    auto bc1 = firstcollision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc1);
    LOG(info) << "Processing " << collisions.size() << " collisions";
    auto pools = makeTracksPool(collisions, tracks, ambiTracks, bcWtmp);
    LOG(info) << "Mother track pool created: " << pools[POS].size() << " positive and " << pools[NEG].size() << " negative seeds";
    LOG(info) << "Kink daughter track pool created: " << pools[2 + POS].size() << " positive and " << pools[2 + NEG].size() << " negative seeds";

    fillCandidateData(collisions, tracks, ambiTracks, bcWtmp, pools);
    LOG(info) << "Filled " << mKinkCandidates.size() << " kink candidates";
    for (auto& kinkCand : mKinkCandidates) {
      outputDataTable(kinkCand.primVtx[0], kinkCand.primVtx[1], kinkCand.primVtx[2],
                      kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                      kinkCand.isMatter, kinkCand.recoPtHyp(), kinkCand.recoPhiHyp(), kinkCand.recoEtaHyp(),
                      kinkCand.recoPtTrit(), kinkCand.recoPhiTrit(), kinkCand.recoEtaTrit(),
                      kinkCand.hypDCAXY, kinkCand.tritDCAXY, kinkCand.dcaKinkTopo,
                      kinkCand.clusterSizeITSHyp, kinkCand.clusterSizeITSTrit, kinkCand.trackingPIDTriton,
                      kinkCand.momTritTPC, kinkCand.tpcSignalTrit, kinkCand.nSigmaTPCTrit, kinkCand.nSigmaTOFTrit);
    }
  }
  PROCESS_SWITCH(hyperKinkRecoTask, processData, "Data analysis", true);

  void processMC(CollisionsFull const& collisions, TracksFull const& tracks, o2::aod::AmbiguousTracks const& ambiTracks, aod::BCsWithTimestamps const& bcWtmp, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC)
  {
    filledMothers.clear();
    mKinkCandidates.clear();
    if (mBBparamsTrit[5] < 0) {
      LOG(info) << "Bethe-Bloch parameters for Triton not set, using default nSigma";
    }

    auto firstcollision = collisions.begin();
    auto bc1 = firstcollision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc1);
    LOG(info) << "Processing " << collisions.size() << " collisions";
    auto pools = makeTracksPool(collisions, tracks, ambiTracks, bcWtmp);
    LOG(info) << "Mother track pool created: " << pools[POS].size() << " positive and " << pools[NEG].size() << " negative seeds";
    LOG(info) << "Kink daughter track pool created: " << pools[2 + POS].size() << " positive and " << pools[2 + NEG].size() << " negative seeds";

    fillCandidateData(collisions, tracks, ambiTracks, bcWtmp, pools);
    fillMCinfo(trackLabelsMC, particlesMC);

    // now we fill only the signal candidates that were not reconstructed
    std::vector<int> mcToKinkCandidates;
    mcToKinkCandidates.resize(particlesMC.size(), -1);

    for (auto& mcPart : particlesMC) {

      if (std::abs(mcPart.pdgCode()) != hyperPdg)
        continue;
      std::array<float, 3> secVtx;
      std::array<float, 3> primVtx = {mcPart.vx(), mcPart.vy(), mcPart.vz()};
      std::array<float, 3> momMother = {mcPart.px(), mcPart.py(), mcPart.pz()};
      std::array<float, 3> momTrit;
      bool isTritFound = false;
      int hyperMCIndex = mcPart.globalIndex();
      for (auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(mcDaught.pdgCode()) == tritDauPdg) {
          secVtx = {mcDaught.vx(), mcDaught.vy(), mcDaught.vz()};
          momTrit = {mcDaught.px(), mcDaught.py(), mcDaught.pz()};
          isTritFound = true;
          break;
        }
      }
      if (!isTritFound) {
        continue;
      }

      if (std::find(filledMothers.begin(), filledMothers.end(), mcPart.globalIndex()) != std::end(filledMothers)) {
        continue;
      }
      kinkCandidate kinkCand;
      kinkCand.pdgCode = mcPart.pdgCode();
      for (int i = 0; i < 3; i++) {
        kinkCand.gDecVtx[i] = secVtx[i] - primVtx[i];
        kinkCand.gMomHyp[i] = momMother[i];
        kinkCand.gMomTrit[i] = momTrit[i];
      }
      kinkCand.hyperTrackID = -1;
      kinkCand.tritTrackID = -1;
      kinkCand.isSignal = true;
      mKinkCandidates.push_back(kinkCand);
      mcToKinkCandidates[hyperMCIndex] = mKinkCandidates.size() - 1;
    }

    // look for hypertriton or triton tracks
    for (auto& track : tracks) {
      auto mcLabel = trackLabelsMC.rawIteratorAt(track.globalIndex());
      if (mcLabel.has_mcParticle()) {
        auto mcTrack = mcLabel.mcParticle_as<aod::McParticles>();
        if (mcToKinkCandidates[mcTrack.globalIndex()] < 0 || !track.hasITS()) {
          continue;
        }
        auto& kinkCand = mKinkCandidates[mcToKinkCandidates[mcTrack.globalIndex()]];
        kinkCand.mcMask = mcLabel.mcMask();
        kinkCand.itsPt = track.pt();
      }
    }

    for (auto& kinkCand : mKinkCandidates) {
      if (!kinkCand.isSignal && mcSignalOnly) {
        continue;
      }
      int chargeFactor = -1 + 2 * (kinkCand.pdgCode > 0);
      outputMCTable(kinkCand.primVtx[0], kinkCand.primVtx[1], kinkCand.primVtx[2],
                    kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                    kinkCand.isMatter, kinkCand.recoPtHyp(), kinkCand.recoPhiHyp(), kinkCand.recoEtaHyp(),
                    kinkCand.recoPtTrit(), kinkCand.recoPhiTrit(), kinkCand.recoEtaTrit(),
                    kinkCand.hypDCAXY, kinkCand.tritDCAXY, kinkCand.dcaKinkTopo,
                    kinkCand.clusterSizeITSHyp, kinkCand.clusterSizeITSTrit, kinkCand.trackingPIDTriton,
                    kinkCand.momTritTPC, kinkCand.tpcSignalTrit, kinkCand.nSigmaTPCTrit, kinkCand.nSigmaTOFTrit,
                    kinkCand.gDecVtx[0], kinkCand.gDecVtx[1], kinkCand.gDecVtx[2],
                    kinkCand.genPt() * chargeFactor, kinkCand.genPtTrit(),
                    kinkCand.isReco, kinkCand.isSignal, kinkCand.mcMask, kinkCand.itsPt);
    }
  }
  PROCESS_SWITCH(hyperKinkRecoTask, processMC, "MC analysis", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hyperKinkRecoTask>(cfgc)};
}
