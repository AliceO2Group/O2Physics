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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PHOSBase/Geometry.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>

/// \struct PHOS pi0 analysis
/// \brief Monitoring task for PHOS related quantities
/// \author Dmitri Peresunko, NRC "Kurchatov institute"
/// \since Nov, 2022
///

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phosAlign {

  static constexpr int16_t kCpvX = 7; // grid 13 steps along z and 7 along phi as largest match ellips 20x10 cm
  static constexpr int16_t kCpvZ = 6;
  static constexpr int16_t kCpvCells = 4 * kCpvX * kCpvZ; // 4 modules
  static constexpr float cpvMaxX = 73;                    // max CPV coordinate phi
  static constexpr float cpvMaxZ = 63;                    // max CPV coordinate z

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl>;
  using mcClusters = soa::Join<aod::CaloClusters, aod::PHOSCluLabels>;
  using mcTracks = soa::Join<tracks, aod::McTrackLabels>;

  Configurable<float> mMinE{"mMinCluE", 0.3, "Minimum cluster energy for analysis"};
  Configurable<float> mMinCluTime{"minCluTime", -25.e-9, "Min. cluster time"};
  Configurable<float> mMaxCluTime{"mMaxCluTime", 25.e-9, "Max. cluster time"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  double bz{0.}; // magnetic field
  int runNumber{0};

  class trackMatch
  {
   public:
    trackMatch() = default;
    trackMatch(double x, double z, double p, int32_t l, bool el, bool charge) : pX(x), pZ(z), mom(p), label(l), isEl(el), isPos(charge) {}
    ~trackMatch() = default;

   public:
    double pX = 9999.; // X (phi) track coordinate in PHOS plane
    double pZ = 9999.; // Z (theta) track coordinate in PHOS plane
    double mom = 0.;   // track momentum
    int32_t label = 0;
    bool isEl = false;  // is electron from TPC dEdx
    bool isPos = false; // is positive charge
  };

  HistogramRegistry mHistManager{"phosAlignHistograms"};

  /// \brief Create output histograms
  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS alignment task ...";
    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    const AxisSpec
      axisPi0Pt{100, 0., 20., "p_{T}^{#pi} (GeV/c)", "p_{T}^{#pi} (GeV/c)"},
      axisPi0m{500, 0., 1., "m_{#gamma#gamma} (GeV/c^{2})", "m_{#gamma#gamma} (GeV/c^{2})"},
      axisModComb{10, 0., 10.},
      axisEpEclu{400, 0., 4., "E_{clu} (GeV)", "E_{clu} (GeV)"},
      axisEp{400, 0., 4., "E/p", "E_{clu}/p_{tr}"},
      axisP{400, 0., 4., "p_{tr} (GeV/c)", "p_{tr} (GeV/c)"},
      axisModes{4, 1., 5., "module", "module"},
      axisX{150, -75., 75., "x (cm)", "x (cm)"},
      axisdX{100, -20., 20., "x_{tr}-x_{clu} (cm)", "x_{tr}-x_{clu} (cm)"},
      axisZ{150, -75., 75, "z (cm)", "z (cm)"},
      axisdZ{64, -20., 20., "z_{tr}-z_{clu} (cm)", "z_{tr}-z_{clu} (cm)"};

    // mHistManager.add("mggRe", "inv mass", HistType::kTH3F, {axisPi0m, axisPi0Pt, axisModComb});
    // mHistManager.add("mggMi", "inv mass", HistType::kTH3F, {axisPi0m, axisPi0Pt, axisModComb});
    // mHistManager.add("mggReDisp", "inv mass", HistType::kTH3F, {axisPi0m, axisPi0Pt, axisModComb});
    // mHistManager.add("mggMiDisp", "inv mass", HistType::kTH3F, {axisPi0m, axisPi0Pt, axisModComb});

    mHistManager.add("hEpAll", "E/p ratio, all tracks", HistType::kTH3F, {axisEpEclu, axisEp, axisModes});
    mHistManager.add("hEpEl", "E/p ratio, electrons", HistType::kTH3F, {axisEpEclu, axisEp, axisModes});
    mHistManager.add("hEpAllDisp", "E/p ratio, all tracks", HistType::kTH3F, {axisEpEclu, axisEp, axisModes});
    mHistManager.add("hEpElDisp", "E/p ratio, electrons", HistType::kTH3F, {axisEpEclu, axisEp, axisModes});
    // if(isMC){
    //   mHistManager.add("hMCEpAll", "E/p ratio, all tracks", HistType::kTH3F, {axisEpEclu, axisEpP, axisModes});
    //   mHistManager.add("hMCEpEl", "E/p ratio, electrons", HistType::kTH3F, {axisEpEclu, axisEpP, axisModes});
    //   mHistManager.add("hMCEpAllDisp", "E/p ratio, all tracks", HistType::kTH3F, {axisEpEclu, axisEpP, axisModes});
    //   mHistManager.add("hMCEpElDisp", "E/p ratio, electrons", HistType::kTH3F, {axisEpEclu, axisEpP, axisModes});
    // }
    mHistManager.add("hdZvsZ", "dz(z), all tracks", HistType::kTH3F, {axisdZ, axisZ, axisModes});
    mHistManager.add("hdXvsX", "dx(x), all tracks", HistType::kTH3F, {axisdX, axisX, axisModes});
    mHistManager.add("hdZvsZ_plus", "dz(z), pos tracks", HistType::kTH3F, {axisdZ, axisZ, axisModes});
    mHistManager.add("hdXvsX_plus", "dx(x), pos tracks", HistType::kTH3F, {axisdX, axisX, axisModes});
    mHistManager.add("hdZvsZ_minus", "dz(z), neg tracks", HistType::kTH3F, {axisdZ, axisZ, axisModes});
    mHistManager.add("hdXvsX_minus", "dx(x), neg tracks", HistType::kTH3F, {axisdX, axisX, axisModes});
    mHistManager.add("hdZvsZEl", "dz(z), el tracks", HistType::kTH3F, {axisdZ, axisZ, axisModes});
    mHistManager.add("hdXvsXEl", "dx(x), el tracks", HistType::kTH3F, {axisdX, axisX, axisModes});
    mHistManager.add("hdXdZE", "dx,dz,E_{clu}", HistType::kTH3F, {axisdX, axisdX, axisEpEclu});
    mHistManager.add("hdXdZp", "dx,dz,p_{tr}", HistType::kTH3F, {axisdX, axisdX, axisP});
    mHistManager.add("hXYZ", "xyz", HistType::kTH3F, {{200, -300., 300.}, {100, -500, -250.}, {200, -150., 150.}});

    // mHistManager.add("hdXvsXvsEElM1", "dz(z), el tracks", HistType::kTH3F, {axisdX, axisX, axisEpEclu});
    // mHistManager.add("hdXvsXvsEElM2", "dz(z), el tracks", HistType::kTH3F, {axisdX, axisX, axisEpEclu});
    // mHistManager.add("hdXvsXvsEElM3", "dz(z), el tracks", HistType::kTH3F, {axisdX, axisX, axisEpEclu});
    // mHistManager.add("hdXvsXvsEElM4", "dz(z), el tracks", HistType::kTH3F, {axisdX, axisX, axisEpEclu});

    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");
  }

  /// \brief match tracks and clusters in different PHOS modules
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               mcClusters& clusters,
               mcTracks& tracks,
               aod::BCsWithTimestamps const&)
  {

    // Set the magnetic field from ccdb.
    // The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
    // but this is not true when running on Run2 data/MC already converted into AO2Ds.
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }

    // Nothing to process
    if (clusters.size() == 0) {
      return;
    }
    // calculate coordinates of tracks extrapolated to PHOS and remember track parameters
    std::vector<trackMatch> trackMatchPoints[kCpvCells]; // tracks hit in grid/cell in PHOS
    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        continue;
      }
      if (std::abs(track.collision_as<SelCollisions>().posZ()) > 10.f) {
        continue;
      }

      // calculate coordinate in PHOS plane
      if (std::abs(track.eta()) > 0.3) {
        continue;
      }

      int16_t module;
      float trackX = 999., trackZ = 999.;

      auto trackPar = getTrackPar(track);
      if (!impactOnPHOS(trackPar, track.trackEtaEmcal(), track.trackPhiEmcal(), track.collision_as<SelCollisions>().posZ(), module, trackX, trackZ)) {
        continue;
      }
      int regionIndex = matchIndex(module, trackX, trackZ);

      // PID
      float nsigmaTPCEl = track.tpcNSigmaEl();
      // First rough PID
      const double nSigmaBelowElectronLine = -3.;
      const double nSigmaAboveElectronLine = 5.;
      const double minPnSigmaAbovePionLine = 1.;
      const double maxPnSigmaAbovePionLine = 3.;
      const double nSigmaAbovePionLine = 0;
      bool electron = nsigmaTPCEl > nSigmaBelowElectronLine && nsigmaTPCEl < nSigmaAboveElectronLine;

      if (track.p() > minPnSigmaAbovePionLine && track.p() < maxPnSigmaAbovePionLine) {
        if (track.tpcNSigmaPi() < nSigmaAbovePionLine) {
          electron = false;
        }
      }
      // Strict dEdx
      if (track.p() > minPnSigmaAbovePionLine && track.p() < maxPnSigmaAbovePionLine) {
        if (track.tpcNSigmaPi() < 2.) {
          electron = false;
        }
      }
      const double minPPionRejection = 0.5;
      const double sigmaAroundLine = 1.;
      if (track.p() < minPPionRejection) {
        if (TMath::Abs(track.tpcNSigmaPi()) < sigmaAroundLine) {
          electron = false;
        }
      }
      // Kaon rejection
      const double minPKaonRejection = 1.5;
      if (track.p() < minPKaonRejection) {
        if (TMath::Abs(track.tpcNSigmaKa()) < sigmaAroundLine) {
          electron = false;
        }
      }
      // Proton rejection
      const double minPProtonRejection = 2.;
      if (track.p() < minPProtonRejection) {
        if (TMath::Abs(track.tpcNSigmaPr()) < sigmaAroundLine) {
          electron = false;
        }
      }
      trackMatchPoints[regionIndex].emplace_back(trackX, trackZ, track.p(), track.mcParticleId(), electron, static_cast<bool>(track.sign() > 0));
    }

    for (const auto& clu : clusters) {
      if (clu.e() < mMinE) {
        continue;
      }
      if (clu.time() < mMinCluTime || clu.time() > mMaxCluTime) {
        continue;
      }

      // int label = -1;             // if no MC
      // auto mcList = clu.labels(); // const std::vector<int>
      // if (mcList.size() > 0) {
      //   label = mcList[0];
      // }

      // CPV and track match
      const float cellSizeX = 2 * cpvMaxX / kCpvX;
      const float cellSizeZ = 2 * cpvMaxZ / kCpvZ;
      float posX = clu.x(), posZ = clu.z();
      int module = clu.mod();
      // look 9 CPV regions around PHOS cluster
      int phosIndex = matchIndex(clu.mod(), clu.x(), clu.z());
      std::vector<int> regions;
      regions.push_back(phosIndex);
      if (posX > -cpvMaxX + cellSizeX) {
        if (posZ > -cpvMaxZ + cellSizeZ) { // bottom left
          regions.push_back(phosIndex - kCpvZ - 1);
        }
        regions.push_back(phosIndex - kCpvZ);
        if (posZ < cpvMaxZ - cellSizeZ) { // top left
          regions.push_back(phosIndex - kCpvZ + 1);
        }
      }
      if (posZ > -cpvMaxZ + cellSizeZ) { // bottom
        regions.push_back(phosIndex - 1);
      }
      if (posZ < cpvMaxZ - cellSizeZ) { // top
        regions.push_back(phosIndex + 1);
      }
      if (posX < cpvMaxX - cellSizeX) {
        if (posZ > -cpvMaxZ + cellSizeZ) { // bottom right
          regions.push_back(phosIndex + kCpvZ - 1);
        }
        regions.push_back(phosIndex + kCpvZ);
        if (posZ < cpvMaxZ - cellSizeZ) { // top right
          regions.push_back(phosIndex + kCpvZ + 1);
        }
      }
      float dx = 9999., dz = 9999., trackdist = 9999., trackMom = 0.;
      bool electron = false, posCh = false;
      for (int indx : regions) {
        for (auto pp : trackMatchPoints[indx]) {
          float d = pow(pp.pX - posX, 2) + pow((pp.pZ - posZ), 2);
          if (d < trackdist) {
            trackdist = d;
            dx = pp.pX - posX;
            dz = pp.pZ - posZ;
            trackMom = pp.mom;
            electron = pp.isEl;
            posCh = pp.isPos;
          }
        }
      }
      if (trackMom > 0) { // track found
        bool isDispOK = testLambda(clu.e(), clu.m02(), clu.m20());

        mHistManager.fill(HIST("hEpAll"), clu.e(), trackMom / clu.e(), module);
        if (isDispOK) {
          mHistManager.fill(HIST("hEpAllDisp"), clu.e(), trackMom / clu.e(), module);
          if (electron) {
            mHistManager.fill(HIST("hEpEl"), clu.e(), trackMom / clu.e(), module);
            if (isDispOK) {
              mHistManager.fill(HIST("hEpElDisp"), clu.e(), trackMom / clu.e(), module);
            }
          }
        }
        mHistManager.fill(HIST("hdZvsZ"), dz, posZ, module);
        mHistManager.fill(HIST("hdXvsX"), dx, posX, module);
        mHistManager.fill(HIST("hdXdZE"), dx, dz, clu.e());
        mHistManager.fill(HIST("hdXdZp"), dx, dz, trackMom);
        if (posCh) {
          mHistManager.fill(HIST("hdZvsZ_plus"), dz, posZ, module);
          mHistManager.fill(HIST("hdXvsX_plus"), dx, posX, module);
        } else {
          mHistManager.fill(HIST("hdZvsZ_minus"), dz, posZ, module);
          mHistManager.fill(HIST("hdXvsX_minus"), dx, posX, module);
        }

        if (electron) {
          mHistManager.fill(HIST("hdZvsZEl"), dz, posZ, module);
          mHistManager.fill(HIST("hdXvsXEl"), dx, posX, module);
        }
      }
    }
  }
  int matchIndex(int16_t module, float x, float z)
  {
    // calculate cell index in grid over PHOS detector
    const float cellSizeX = 2 * cpvMaxX / kCpvX;
    const float cellSizeZ = 2 * cpvMaxZ / kCpvZ;
    // in track matching tracks can be beyond CPV surface
    // assign these tracks to the closest cell
    int ix = std::max(0, static_cast<int>((x + cpvMaxX) / cellSizeX));
    int iz = std::max(0, static_cast<int>((z + cpvMaxZ) / cellSizeZ));
    if (ix >= kCpvX) {
      ix = kCpvX - 1;
    }
    if (iz >= kCpvZ) {
      iz = kCpvZ - 1;
    }
    return (module - 1) * kCpvX * kCpvZ + ix * kCpvZ + iz; // modules: 1,2,3,4
  }

  bool impactOnPHOS(o2::track::TrackParametrization<float>& trackPar, float trackEta, float trackPhi, float /*zvtx*/, int16_t& module, float& trackX, float& trackZ)
  {
    // eta,phi was calculated at EMCAL radius.
    // Extrapolate to PHOS assuming zeroB and current vertex
    // Check if direction in PHOS acceptance+20cm and return phos module number and coordinates in PHOS module plane
    const float phiMin = 240. * 0.017453293; // degToRad
    const float phiMax = 323. * 0.017453293; // PHOS+20 cm * degToRad
    const float etaMax = 0.178266;
    if (trackPhi < phiMin || trackPhi > phiMax || abs(trackEta) > etaMax) {
      return false;
    }

    const float dphi = 20. * 0.017453293;
    if (trackPhi < 0.) {
      trackPhi += TMath::TwoPi();
    }
    if (trackPhi > TMath::TwoPi()) {
      trackPhi -= TMath::TwoPi();
    }
    module = 1 + static_cast<int16_t>((trackPhi - phiMin) / dphi);
    if (module < 1) {
      module = 1;
    }
    if (module > 4) {
      module = 4;
    }

    // get PHOS radius
    constexpr float shiftY = -1.26;    // Depth-optimized
    double posL[3] = {0., 0., shiftY}; // local position at the center of module
    double posG[3] = {0};
    geomPHOS->getAlignmentMatrix(module)->LocalToMaster(posL, posG);
    double rPHOS = sqrt(posG[0] * posG[0] + posG[1] * posG[1]);
    double alpha = (230. + 20. * module) * 0.017453293;

    // During main reconstruction track was propagated to radius 460 cm with accounting material
    // now material is not available. Therefore, start from main rec. position and extrapoate to actual radius without material
    // Get track parameters at point where main reconstruction stop
    float xPHOS = 460.f, xtrg = 0.f;
    if (!trackPar.getXatLabR(xPHOS, xtrg, bz, o2::track::DirType::DirOutward)) {
      return false;
    }
    auto prop = o2::base::Propagator::Instance();
    if (!trackPar.rotate(alpha) ||
        !prop->PropagateToXBxByBz(trackPar, xtrg, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    // calculate xyz from old (Phi, eta) and new r
    float r = std::sqrt(trackPar.getX() * trackPar.getX() + trackPar.getY() * trackPar.getY());
    trackPar.setX(r * std::cos(trackPhi - alpha));
    trackPar.setY(r * std::sin(trackPhi - alpha));
    trackPar.setZ(r / std::tan(2. * std::atan(std::exp(-trackEta))));

    if (!prop->PropagateToXBxByBz(trackPar, rPHOS, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    alpha = trackPar.getAlpha();
    double ca = cos(alpha), sa = sin(alpha);
    posG[0] = trackPar.getX() * ca - trackPar.getY() * sa;
    posG[1] = trackPar.getY() * ca + trackPar.getX() * sa;
    posG[2] = trackPar.getZ();

    mHistManager.fill(HIST("hXYZ"), posG[0], posG[1], posG[2]);
    geomPHOS->getAlignmentMatrix(module)->MasterToLocal(posG, posL);
    trackX = posL[0];
    trackZ = posL[1];
    return true;
  }

  //_____________________________________________________________________________
  int moduleCombination(int m1, int m2)
  {
    // enumerates possible module combinations
    // (1,1)=0, (2,2)=1, (3,3)=2, (4,4)=3, (1,2)=(2,1)=4, (2,3)=(3,2)=5, (3,4)=(4,3)=6, (1,3)=(3,1)=7,
    // (2,4)=(4,2)=8, (1,4)=(4,1)=9
    int d = TMath::Abs(m1 - m2);
    if (d == 0) {
      return m1 - 1;
    }
    if (d == 1) {
      return 3 + TMath::Min(m1, m2);
    }
    if (d == 2) {
      return 6 + TMath::Min(m1, m2);
    }
    return 9;
  }
  //_____________________________________________________________________________
  bool testLambda(float pt, float l1, float l2)
  {
    // Parameterization for full dispersion
    // Parameterizatino for full dispersion
    float l2Mean = 1.53126 + 9.50835e+06 / (1. + 1.08728e+07 * pt + 1.73420e+06 * pt * pt);
    float l1Mean = 1.12365 + 0.123770 * TMath::Exp(-pt * 0.246551) + 5.30000e-03 * pt;
    float l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
    float l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
    float c = -0.35 - 0.550 * TMath::Exp(-0.390730 * pt);

    return 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
             0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
             0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma <
           4.;
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosAlign>(cfgc)};
}
