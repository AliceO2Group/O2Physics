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
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  3-body Hyperhelium 4 builder
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
// with contributions from:
// --- navneet.kumar@cern.ch
// --- natasha.sharma@cern.ch
// --- kumar.lokesh@cern.ch
// --- david.dobrigkeit.chinellato@cern.ch

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGLF/DataModel/LFHyperhelium4Tables.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

// For dE/dx association in pre-selection
// using TracksExtraWithPID = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullHe>;

// For MC and dE/dx association
using TracksExtraWithPIDandLabels = soa::Join<aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullHe, aod::McTrackLabels>;

// Pre-selected Decay3Bodys
using TaggedHyHe4Candidates = soa::Join<aod::Decay3Bodys, aod::HyHe4Tags>;

// For MC association in pre-selection
using LabeledTracksExtra = soa::Join<aod::TracksExtra, aod::McTrackLabels>;

struct hyhefourbuilder {
  Produces<aod::StoredHyHe4Datas> hyhe4data;

  Service<o2::ccdb::BasicCCDBManager> ccdb; // <-- for CCDB access

  Configurable<int> tpcrefit{"tpcrefit", 0, "demand TPC refit"};
  Preslice<aod::Decay3Bodys> perCollision = o2::aod::decay3body::collisionId; // used for grouping, here we use the collisionID

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};
  Configurable<int> useMatCorrType{"useMatCorrType", 2, "0: none, 1: TGeo, 2: LUT"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // Basic selection criteria
  Configurable<float> selHyHe4daudca{"selHyHe4daudca", 1, "DCA between HyHe4 daughters"};

  // Filter out HyHe4 that are interesting!
  Filter taggedFilter = aod::hyhe4tag::isInteresting == true;

  // bookkeeping propagation / run number / magnetic field
  int mRunNumber;
  float d_bz;
  float maxSnp;                            // max sine phi for propagation
  float maxStep;                           // max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr; // lut pointer

  // storing output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Define o2 fitter, 3-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<3> fitter;

  void init(InitContext& context)
  {
    const AxisSpec axisMassHyHe4{(int)1000, 3.5f, 4.5f, "Hyperhelium4 Mass Distribution (GeV/c^{2})"};
    const AxisSpec axisNCandidates{(int)100, 0.0f, 100.0f, "Number of 3-body candidates"};

    const AxisSpec axisEventCounter{(int)1, 0.0f, 1.0f, "Number of events"};
    const AxisSpec axisPt{(int)1000, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    // Base counters
    histos.add("hNEvents", "hNEvents", kTH1F, {axisEventCounter});
    histos.add("hNCandidates", "hNCandidates", kTH1F, {axisPt});

    ccdb->setURL(ccdburl);
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

    // initialize O2 3-prong fitter (only once)
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
    fitter.setWeightedFinalPCA(d_UseWeightedPCA);

    // Material correction in the DCA fitter
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    fitter.setMatCorrType(matCorr);
  }

  // Helper struct to pass HyHe4 information
  struct {
    int sign;
    int index3He;
    int indexProton;
    int indexPion;
    float track3HeX;
    float trackProtonX;
    float trackPionX;
    std::array<float, 3> pos;
    std::array<float, 3> mom3He;
    std::array<float, 3> momProton;
    std::array<float, 3> momPion;
    float dcaHyHe4dau;
    float dcaXY3He;
    float dcaXYProton;
    float dcaXYPion;
    float dcaXY;
    float dcaZ;
    float decayRadius;
    float massHyHe4;
    float massAntiHyHe4;
  } hyHe4Candidate;

  // ccdb to fetch the magnetic fieldd
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
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
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
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
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  template <typename TCollision, typename TTrack>
  bool buildHyHe4Candidate(TCollision const& collision, TTrack const& Helium3, TTrack const& proton, TTrack const& pion)
  {
    hyHe4Candidate.index3He = Helium3.globalIndex();
    hyHe4Candidate.indexProton = proton.globalIndex();
    hyHe4Candidate.indexPion = pion.globalIndex();

    //---/---/---/---/---/---/---/---/---/---/---/---/---/
    // Create TrackParCovs for DCA (N.B. propagation irreversibility wrt covmats)
    auto lHelium3TrackForDCA = getTrackParCov(Helium3);
    auto lProtonTrackForDCA = getTrackParCov(proton);
    auto lPionTrackForDCA = getTrackParCov(pion);

    // determine base property of charge for matter vs antimatter
    hyHe4Candidate.sign = +1;
    if (Helium3.sign() < 0)
      hyHe4Candidate.sign = -1;

    //---/---/---/---/---/---/---/---/---/---/---/---/---/
    // Calculate DCA with respect to the collision associated to the V0, not individual tracks
    gpu::gpustd::array<float, 2> dcaInfo;
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lHelium3TrackForDCA, 2.f, fitter.getMatCorrType(), &dcaInfo);
    hyHe4Candidate.dcaXY3He = dcaInfo[0];
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lProtonTrackForDCA, 2.f, fitter.getMatCorrType(), &dcaInfo);
    hyHe4Candidate.dcaXYProton = dcaInfo[0];
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lPionTrackForDCA, 2.f, fitter.getMatCorrType(), &dcaInfo);
    hyHe4Candidate.dcaXYPion = dcaInfo[0];

    //---/---/---/---/---/---/---/---/---/---/---/---/---/
    // Move close to minima
    auto lHelium3Track = getTrackParCov(Helium3);
    auto lProtonTrack = getTrackParCov(proton);
    auto lPionTrack = getTrackParCov(pion);
    int nCand = 0;
    try {
      nCand = fitter.process(lHelium3Track, lProtonTrack, lPionTrack);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return false;
    }
    if (nCand == 0) {
      return false;
    }

    //---/---/---/---/---/---/---/---/---/---/---/---/---/
    // Get base properties
    fitter.getTrack(0).getPxPyPzGlo(hyHe4Candidate.mom3He);
    fitter.getTrack(1).getPxPyPzGlo(hyHe4Candidate.momProton);
    fitter.getTrack(2).getPxPyPzGlo(hyHe4Candidate.momPion);

    hyHe4Candidate.track3HeX = fitter.getTrack(0).getX();
    hyHe4Candidate.trackProtonX = fitter.getTrack(1).getX();
    hyHe4Candidate.trackPionX = fitter.getTrack(2).getX();

    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      hyHe4Candidate.pos[i] = vtx[i];
    }

    // baseline selection
    hyHe4Candidate.dcaHyHe4dau = TMath::Sqrt(fitter.getChi2AtPCACandidate());
    if (hyHe4Candidate.dcaHyHe4dau > selHyHe4daudca)
      return false;

    //---/---/---/---/---/---/---/---/---/---/---/---/---/
    // Calculate DCAxy of the HyHe4 (with bending)
    auto lHyHe4Track = fitter.createParentTrackPar();
    lHyHe4Track.setAbsCharge(2); // to be sure
    dcaInfo[0] = 999;
    dcaInfo[1] = 999;
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lHyHe4Track, 2.f, fitter.getMatCorrType(), &dcaInfo);
    hyHe4Candidate.dcaXY = dcaInfo[0];
    hyHe4Candidate.dcaZ = dcaInfo[1];
    return true;
  }

  template <class TTrackTo, typename T3BodyTable>
  void buildHyHe4Tables(T3BodyTable const& d3Bodys)
  {
    for (const auto& d3body : d3Bodys) {
      auto const& collision = d3body.collision();
      auto const& trackProng0 = d3body.template track0_as<TTrackTo>();
      auto const& trackProng1 = d3body.template track1_as<TTrackTo>();
      auto const& trackProng2 = d3body.template track2_as<TTrackTo>();

      // sign -> if positive, matter; if negative, antimatter
      int sign = trackProng0.sign() + trackProng1.sign() + trackProng2.sign();

      // N.B. TPC signal ordering may not be ideal
      //      consider also momentum ordering (explore kinematics first)

      if (trackProng0.sign() != sign) {                          // prong 0 = pion
        if (trackProng1.tpcSignal() < trackProng2.tpcSignal()) { // 0 = pi, 1 = p, 2 = 3He
          buildHyHe4Candidate(collision, trackProng2, trackProng1, trackProng0);
        }
        if (trackProng1.tpcSignal() > trackProng2.tpcSignal()) { // 0 = pi, 1 = p, 2 = 3He
          buildHyHe4Candidate(collision, trackProng1, trackProng2, trackProng0);
        }
      }
      if (trackProng1.sign() != sign) {                          // prong 1 = pion
        if (trackProng0.tpcSignal() > trackProng2.tpcSignal()) { // 0 = 3He, 1 = pi, 2 = p
          buildHyHe4Candidate(collision, trackProng0, trackProng2, trackProng1);
        }
        if (trackProng0.tpcSignal() < trackProng2.tpcSignal()) { // 0 = p, 1 = pi, 2 = 3He
          buildHyHe4Candidate(collision, trackProng2, trackProng0, trackProng1);
        }
      }
      if (trackProng2.sign() != sign) {                          // prong 2 = pion
        if (trackProng0.tpcSignal() > trackProng1.tpcSignal()) { // 0 = 3He, 1 = p, 2 = pi
          buildHyHe4Candidate(collision, trackProng0, trackProng1, trackProng2);
        }
        if (trackProng0.tpcSignal() < trackProng1.tpcSignal()) { // 0 = p, 1 = 3He, 2 = pi
          buildHyHe4Candidate(collision, trackProng1, trackProng0, trackProng2);
        }
      }

      // Fill the table of candidates for posterior analysis
      hyhe4data(hyHe4Candidate.index3He, hyHe4Candidate.indexProton, hyHe4Candidate.indexPion,
                collision.globalIndex(), d3body.globalIndex(), hyHe4Candidate.sign,
                hyHe4Candidate.track3HeX, hyHe4Candidate.trackProtonX, hyHe4Candidate.trackPionX,
                2.0f * hyHe4Candidate.mom3He[0], 2.0f * hyHe4Candidate.mom3He[1], 2.0f * hyHe4Candidate.mom3He[2],
                hyHe4Candidate.momProton[0], hyHe4Candidate.momProton[1], hyHe4Candidate.momProton[2],
                hyHe4Candidate.momPion[0], hyHe4Candidate.momPion[1], hyHe4Candidate.momPion[2],
                hyHe4Candidate.pos[0], hyHe4Candidate.pos[1], hyHe4Candidate.pos[2],
                hyHe4Candidate.dcaHyHe4dau, hyHe4Candidate.dcaXY3He, hyHe4Candidate.dcaXYProton, hyHe4Candidate.dcaXYPion,
                hyHe4Candidate.dcaXY, hyHe4Candidate.dcaZ);
    }
  }

  void process(aod::Collisions const& collisions, soa::Filtered<TaggedHyHe4Candidates> const& d3bodys, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    long eventCounter = 0;
    for (const auto& collision : collisions) {
      eventCounter++;
      // Fire up CCDB
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = collision.globalIndex();
      auto d3Bodys_thisCollision = d3bodys.sliceBy(perCollision, collIdx); // To access the 3body decays only collisons or I can say the collision with contain the 3body decay

      histos.fill(HIST("hNCandidates"), d3Bodys_thisCollision.size());

      // Do the building
      buildHyHe4Tables<FullTracksExtIU>(d3Bodys_thisCollision);
    }
    histos.fill(HIST("hNEvents"), 0.0, eventCounter);
  }
};

//*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
struct hyHe4Preselector {
  // storing output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Produces<aod::HyHe4Tags> hyhetags; // MC tags
  Preslice<aod::Decay3Bodys> perCollision = o2::aod::decay3body::collisionId;
  //  Configurable<bool> dIfMCgenerateHelium3{"dIfMCgenerateHelium3", true, "if MC, generate MC true Helium3 (yes/no)"};
  Configurable<bool> dIfMCgeneratehyHe4{"dIfMCgeneratehyHe4", true, "if MC, generate MC true hyHe4 (yes/no)"};
  Configurable<bool> dIfMCgenerateAntihyHe4{"dIfMCgenerateAntihyHe4", true, "if MC, generate MC true AntihyHe4 (yes/no)"};

  Configurable<bool> ddEdxPreSelecthyHe4{"ddEdxPreSelhyHe4", true, "pre-select dE/dx compatibility with hyHe4 (yes/no)"};
  Configurable<bool> ddEdxPreSelectAntihyHe4{"ddEdxPreSelectAntihyHe4", true, "pre-select dE/dx compatibility with AntihyHe4 (yes/no)"};

  // dEdx pre-selection compatibility
  Configurable<float> ddEdxPreSelectionWindow{"ddEdxPreSelectionWindow", 7, "Nsigma window for dE/dx preselection"};

  // tpc quality pre-selection
  Configurable<int> dTPCNCrossedRows{"dTPCNCrossedRows", 50, "Minimum TPC crossed rows"};

  void init(InitContext& context)
  {
    const AxisSpec hEventCounter{(int)1, 0.0f, 1.0f, "Number of events"};
    histos.add("hNEvents", "hNEvents", kTH1F, {hEventCounter});
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  /// function to check track quality
  template <class TTrackTo, typename T3Body>
  void checkTrackQuality(T3Body const& d3body, bool& lIsInteresting)
  {
    lIsInteresting = false;
    auto const& trackProng0 = d3body.template track0_as<TTrackTo>();
    auto const& trackProng1 = d3body.template track1_as<TTrackTo>();
    auto const& trackProng2 = d3body.template track2_as<TTrackTo>();
    // cout<<"The number of TPC clusters are "<<trackProng1.tpcNClsCrossedRows() <<endl;
    if ((trackProng0.tpcNClsCrossedRows() >= dTPCNCrossedRows) && (trackProng1.tpcNClsCrossedRows() >= dTPCNCrossedRows) && (trackProng2.tpcNClsCrossedRows() >= dTPCNCrossedRows)) {
      lIsInteresting = true;
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  /// function to check PDG code
  template <class TTrackTo, typename T3Body>
  void checkPDG(T3Body const& d3body, bool& lIsInteresting, bool& lIsTrueHyHelium4, bool& lIsTrueAntiHyHelium4)
  {
    lIsTrueHyHelium4 = false;
    lIsTrueAntiHyHelium4 = false;
    int lPDG = -1;
    lIsInteresting = false;
    auto const& trackProng0 = d3body.template track0_as<TTrackTo>();
    auto const& trackProng1 = d3body.template track1_as<TTrackTo>();
    auto const& trackProng2 = d3body.template track2_as<TTrackTo>();

    int lPDG0 = -1;
    int lPDG1 = -1;
    int lPDG2 = -1;

    // Association check
    // Lets do something to identify the PID of particles
    if (trackProng0.has_mcParticle() && trackProng1.has_mcParticle() && trackProng2.has_mcParticle()) {
      auto lMCtrackProng0 = trackProng0.template mcParticle_as<aod::McParticles>();
      auto lMCtrackProng1 = trackProng1.template mcParticle_as<aod::McParticles>();
      auto lMCtrackProng2 = trackProng2.template mcParticle_as<aod::McParticles>();
      lPDG0 = lMCtrackProng0.pdgCode();
      lPDG1 = lMCtrackProng1.pdgCode();
      lPDG2 = lMCtrackProng2.pdgCode();
      if (lMCtrackProng0.has_mothers() && lMCtrackProng1.has_mothers() && lMCtrackProng2.has_mothers()) {
        for (auto& lProng0Mother : lMCtrackProng0.template mothers_as<aod::McParticles>()) {
          for (auto& lProng1Mother : lMCtrackProng1.template mothers_as<aod::McParticles>()) {
            for (auto& lProng2Mother : lMCtrackProng2.template mothers_as<aod::McParticles>()) {
              if ((lProng0Mother.globalIndex() == lProng1Mother.globalIndex()) && (lProng0Mother.globalIndex() == lProng2Mother.globalIndex())) {
                lPDG = lProng0Mother.pdgCode();
              }
            }
          }
        }
      }
    } // end association check

    if (lPDG == -1010020040 /*&& lPDG0 == -1000020030 && lPDG1 == -2212 && lPDG2 == +211 */ && dIfMCgenerateAntihyHe4) {
      lIsTrueAntiHyHelium4 = true;
      lIsInteresting = true;
    }
    if (lPDG == 1010020040 /*&& lPDG0 == +1000020030 && lPDG1 == +2212 && lPDG2 == -211 */ && dIfMCgeneratehyHe4) {
      lIsTrueHyHelium4 = true;
      lIsInteresting = true;
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  /// This process function ensures that all 3 body candidates are built. It will simply tag everything as true.
  void processBuildAll(aod::Collisions const& collisions, aod::Decay3Bodys const& d3bodys, aod::TracksExtra const&)
  {
    long eventCounter = 0;
    for (const auto& d3body : d3bodys) {
      eventCounter++;
      bool lIsQualityInteresting = false;
      checkTrackQuality<aod::TracksExtra>(d3body, lIsQualityInteresting);
      hyhetags(lIsQualityInteresting, true, true, true, true);
    }
  }
  PROCESS_SWITCH(hyHe4Preselector, processBuildAll, "Switch to build all V0s", false);

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processBuildMCAssociated(aod::Collisions const& collisions, aod::Decay3Bodys const& d3bodys, LabeledTracksExtra const&, aod::McParticles const& particlesMC)
  {
    for (const auto& d3body : d3bodys) {
      bool lIsInteresting = false;
      bool lIsQualityInteresting = false;
      bool lIsTrueHyHelium4 = false;
      bool lIsTrueAntiHyHelium4 = false;
      checkPDG<LabeledTracksExtra>(d3body, lIsInteresting, lIsTrueHyHelium4, lIsTrueAntiHyHelium4);
      checkTrackQuality<LabeledTracksExtra>(d3body, lIsQualityInteresting);
      hyhetags(lIsQualityInteresting * lIsInteresting, lIsTrueHyHelium4, lIsTrueAntiHyHelium4, true, true);
    }
  }
  PROCESS_SWITCH(hyHe4Preselector, processBuildMCAssociated, "Switch to build MC-associated V0s", true);
};

// Extends the v0data table with expression columns
struct hyhe4Initializer {
  Spawns<aod::HyHe4Datas> hyhe4datas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hyhe4Initializer>(cfgc),
    adaptAnalysisTask<hyhefourbuilder>(cfgc),
    adaptAnalysisTask<hyHe4Preselector>(cfgc)};
}
