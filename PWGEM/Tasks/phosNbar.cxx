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
#include <queue>
#include "TRandom.h"

#include "Common/Core/trackUtilities.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"

#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CaloClusters.h"

#include "PHOSBase/Geometry.h"
#include "CommonDataFormat/InteractionRecord.h"

/// \struct phosNbar
/// \brief account Nbar who's clusters appeared within PHOS
/// \author Bakhtin Pavel, Peresunko Dmitri
/// \since Dec. 2023
///
/// This task monitors simply quantities, which allow to identify anti-neutrons,and present some of their propeties
/// - Energy distribution
/// - Time distribution
/// - Count rate in 2D representation

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using PionTracks = soa::Join<aod::pidBayesPi Pions, aod::FullTracks Tracks>; PionTracks const& ptracks

struct phosNbar {

  HistogramRegistry mHistManager{"phosNbaristograms"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using SelCollision = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi>;
  using mcClusters = soa::Join<aod::CaloClusters, aod::PHOSCluLabels>;
  using mcTracks = soa::Join<TrackCandidates, aod::McTrackLabels>;

  // using MatchedClusters = soa::Join<aod::CaloClusters, aod::PHOSMatchedTrack>;

  Configurable<int> mPairingMethod{"mPairingMethod", 0, "0:max CPA, 1: min DCA"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};
  // nbar selection
  Configurable<double> mMinCluE{"minCluE", 0.3, "Minimum cluster energy"};
  Configurable<double> mCpvMinE{"cpvMinE", 200, "Min CPV amplitude"};
  Configurable<int> mNCellMin{"cluNcell", 2, "Min clu multiplicity"};
  Configurable<double> mRvetoSigma{"cluVetoR", 5., "Veto radius in sigma"};
  Configurable<double> mTimeMin{"timeMin", -150.e-9, "Min time cut"};
  Configurable<double> mTimeMax{"timeMax", 150.e-9, "Max time cut"};
  Configurable<double> mDispA{"dispA", -1., "Disp cut, A"};
  Configurable<double> mDispB{"dispB", 4., "Disp cut, B"};
  // track selection
  Configurable<double> mPionDeDxCut{"piondEdx", 3., "pion dE/dx cut in sigma"};
  // Topological cuts
  Configurable<double> mDCAcut{"DCAmax", 0.2, "DCA max cut"};
  Configurable<double> mCPAcut{"CPAmin", 0.96, "CPA min cut"};
  Configurable<uint> mNmix{"nMix", 5, "depth of mixing buffer"};

  static constexpr double c = 29979245800.; // speed of light in cm/sec
  static constexpr double mNbar = 0.939485; // neutron mass
  static constexpr double mpi = 0.13957039; // pion mass
  double mVtxZ{0};                          // primary vertex coordinate in current collision
  double mBz{123456.};                      // Magnetic field to be initialized

  static constexpr int mNZbins = 10; // number of bins of event classification for mixed

  // class to keep nbar candidate parameters
  class nbar
  {
   public:
    nbar() = default;
    nbar(double p, double xPHS, double yPHS, bool zPHS, int mcLabel) : mom(p), x(xPHS), y(yPHS), z(zPHS), label(mcLabel) {}
    ~nbar() = default;
    bool testPIDBit(int ibit) { return (mPIDbits & (1 << ibit)) > 0; }
    void setPIDBit(int ibit) { mPIDbits |= (1 << ibit); }

   public:
    double mom = 0.;  // momentum estimated from Time
    double x = 9999.; // x coordinate in PHOS plane
    double y = 9999.; // y coordinate in PHOS plane
    double z = 9999.; // z coordinate in PHOS plane
    int label = -1;   // label of MC particle
    int mPIDbits = 0; // keep PID cuts
  };

  // class to keep nbar candidate parameters
  class pion
  {
   public:
    pion() = default;
    pion(const pion&) = default;
    pion(const o2::track::TrackParametrization<float>& tr, int mcLabel) : mom(tr), label(mcLabel) {}
    ~pion() = default;
    bool testPIDBit(int ibit) { return (mPIDbits & (1 << ibit)) > 0; }
    void setPIDBit(int ibit) { mPIDbits |= (1 << ibit); }

   public:
    o2::track::TrackParametrization<float> mom; // track parameterization
    int label = -1;                             // label of MC particle
    int mPIDbits = 0;                           // keep PID cuts
  };

  std::vector<pion> piEvent;
  std::vector<nbar> nbarEvent;
  std::array<std::deque<std::vector<pion>>, mNZbins> mixTrackEvts;
  std::array<std::deque<std::vector<nbar>>, mNZbins> mixNbarEvts;

  TH2 *hRePP, *hRePM, *hMiPP, *hMiPM;
  TH3 *hSignalSP, *hSignalSM;
  TH3 *hRePPDCA, *hRePMDCA, *hMiPPDCA, *hMiPMDCA, *hRePPCPA, *hRePMCPA, *hMiPPCPA, *hMiPMCPA;

  /// \brief Create output histograms
  void init(o2::framework::InitContext const&)
  {
    mHistManager.add("evsel", "event selection", HistType::kTH1F, {{10, 0., 10.}});
    mHistManager.add("vtxZ", "Vertex z distribution", HistType::kTH1F, {{100, -20., 20., "z_{vtx} (cm)", "z_{vtx} (cm)"}});
    mHistManager.add("cluCuts", "Spectrum vs cut", HistType::kTH2F, {{100, 0., 5., "E (GeV)", "E (GeV)"}, {10, 0., 10., "cut"}});
    mHistManager.add("nbarCuts", "Spectrum vs cut", HistType::kTH2F, {{100, 0., 5., "E (GeV)", "E (GeV)"}, {10, 0., 10., "cut"}});
    mHistManager.add("trackCuts", "Spectrum vs cut", HistType::kTH2F, {{100, 0., 5., "p_{T} (GeV)", "p_{T} (GeV)"}, {10, 0., 10., "cut"}});

    mHistManager.add("cluTime", "Time vs E clu", HistType::kTH2F, {{200, -100.e-9, 100.e-9, "t (s)", "t (s)"}, {100, 0., 10., "E (GeV)", "E (GeV)"}});
    mHistManager.add("nbarTime", "Time vs E clu nbar", HistType::kTH2F, {{200, -100.e-9, 100.e-9, "t (s)", "t (s)"}, {100, 0., 10., "E (GeV)", "E (GeV)"}});
    mHistManager.add("cluRveto", "CPV radius vs E clu", HistType::kTH2F, {{100, 0., 10., "R (sigmas)", "R (sigmas)"}, {100, 0., 10., "E (GeV)", "E (GeV)"}});
    mHistManager.add("nbarRveto", "CPV radius vs E clu nbar", HistType::kTH2F, {{100, 0., 10., "R (sigmas)", "R (sigmas)"}, {100, 0., 10., "E (GeV)", "E (GeV)"}});

    const AxisSpec
      massAxis{500, 1., 1.5},
      ptAxis{100, 0., 10.},
      dcaAxis{100, 0., 5},
      cpaAxis{100, 0., 1.},
      pidAxis{5, 0., 5.};
    hRePP = (std::get<std::shared_ptr<TH2>>(mHistManager.add("RePiP", "Inv mass", HistType::kTH2F, {massAxis, ptAxis}))).get();
    hRePM = (std::get<std::shared_ptr<TH2>>(mHistManager.add("RePiM", "Inv mass", HistType::kTH2F, {massAxis, ptAxis}))).get();
    hMiPP = (std::get<std::shared_ptr<TH2>>(mHistManager.add("MiPiP", "Inv mass", HistType::kTH2F, {massAxis, ptAxis}))).get();
    hMiPM = (std::get<std::shared_ptr<TH2>>(mHistManager.add("MiPiM", "Inv mass", HistType::kTH2F, {massAxis, ptAxis}))).get();
    hSignalSP = (std::get<std::shared_ptr<TH3>>(mHistManager.add("SignalSP", "Inv mass", HistType::kTH3F, {massAxis, ptAxis, pidAxis}))).get();
    hSignalSM = (std::get<std::shared_ptr<TH3>>(mHistManager.add("SignalSM", "Inv mass", HistType::kTH3F, {massAxis, ptAxis, pidAxis}))).get();

    hRePPDCA = (std::get<std::shared_ptr<TH3>>(mHistManager.add("RePiPDCA", "DCA", HistType::kTH3F, {massAxis, ptAxis, dcaAxis}))).get();
    hRePMDCA = (std::get<std::shared_ptr<TH3>>(mHistManager.add("RePiMDCA", "DCA", HistType::kTH3F, {massAxis, ptAxis, dcaAxis}))).get();
    hMiPPDCA = (std::get<std::shared_ptr<TH3>>(mHistManager.add("MiPiPDCA", "DCA", HistType::kTH3F, {massAxis, ptAxis, dcaAxis}))).get();
    hMiPMDCA = (std::get<std::shared_ptr<TH3>>(mHistManager.add("MiPiMDCA", "DCA", HistType::kTH3F, {massAxis, ptAxis, dcaAxis}))).get();

    hRePPCPA = (std::get<std::shared_ptr<TH3>>(mHistManager.add("RePiPCPA", "CPA", HistType::kTH3F, {massAxis, ptAxis, cpaAxis}))).get();
    hRePMCPA = (std::get<std::shared_ptr<TH3>>(mHistManager.add("RePiMCPA", "CPA", HistType::kTH3F, {massAxis, ptAxis, cpaAxis}))).get();
    hMiPPCPA = (std::get<std::shared_ptr<TH3>>(mHistManager.add("MiPiPCPA", "CPA", HistType::kTH3F, {massAxis, ptAxis, cpaAxis}))).get();
    hMiPMCPA = (std::get<std::shared_ptr<TH3>>(mHistManager.add("MiPiMCPA", "CPA", HistType::kTH3F, {massAxis, ptAxis, cpaAxis}))).get();

    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  /// \brief main function, does the job for data and MC
  template <bool isMC, typename TCollision, typename TTracks, typename TClusters>
  void processAll(TCollision const& collision, TTracks const& tracks, TClusters const& clusters, aod::McParticles const* mcParticles)
  {
    int mixIndex = 0; // index for event classification for mixing
    if (!selectEvent<isMC>(collision, mixIndex)) {
      return;
    }
    selectNbars<isMC>(clusters, mcParticles);
    selectTracks<isMC>(tracks);
    // Fill Real
    double cpa, m, pt;
    math_utils::Point3D<float> vtxV0;
    for (auto tr : piEvent) {
      for (auto nbar : nbarEvent) {
        int cp = 0;
        if constexpr (isMC) { // test parent
          cp = commonParentPDG(tr.label, nbar.label, mcParticles);
        }
        double dca = 999.;
        switch (mPairingMethod) {
          case 0: // maximize CPA
            if (!minimizeCPA(nbar, tr, cpa, vtxV0, m, pt)) {
              continue;
            }
            break;
          case 1: // Minimize DCA
            if (!propagateToDCA(nbar, tr, cpa, dca, vtxV0, m, pt)) {
              continue;
            }
            break;
          default: // do nothing
            continue;
        }
        if (tr.mom.getCharge2Pt() > 0) {
          hRePPDCA->Fill(m, pt, dca);
          hRePPCPA->Fill(m, pt, cpa);
          if ((mPairingMethod == 0 && cpa > mCPAcut) ||
              (mPairingMethod == 1 && dca < mDCAcut && cpa > mCPAcut)) {
            hRePP->Fill(m, pt);
          }
        } else {
          hRePMDCA->Fill(m, pt, dca);
          hRePMCPA->Fill(m, pt, cpa);
          if ((mPairingMethod == 0 && cpa > mCPAcut) ||
              (mPairingMethod == 1 && dca < mDCAcut && cpa > mCPAcut)) {
            hRePM->Fill(m, pt);
          }
        }
        if (cp == -3112) { // Sigma+
          hSignalSP->Fill(m, pt, static_cast<float>(0));
          for (int iPIDbit = 1; iPIDbit < 5; iPIDbit++) {
            if (nbar.testPIDBit(iPIDbit)) {
              hSignalSP->Fill(m, pt, static_cast<float>(iPIDbit));
            }
          }
        }
        if (cp == -3222) { // Sigma-
          hSignalSM->Fill(m, pt, static_cast<float>(0));
          for (int iPIDbit = 1; iPIDbit < 5; iPIDbit++) {
            if (nbar.testPIDBit(iPIDbit)) {
              hSignalSM->Fill(m, pt, static_cast<float>(iPIDbit));
            }
          }
        }
      }
    }
    // FillMixed
    // to avoid too often filling mixed, use events with nbar and tracks
    if (nbarEvent.size() == 0 || piEvent.size() == 0) {
      nbarEvent.clear();
      piEvent.clear();
      return; // do not fill Mixed, do not update stack of events
    }
    for (auto tr : piEvent) {
      for (auto nbarMixEv : mixNbarEvts[mixIndex]) {
        for (auto nbar : nbarMixEv) {
          double dca = 999.;
          switch (mPairingMethod) {
            case 0: // maximize CPA
              if (!minimizeCPA(nbar, tr, cpa, vtxV0, m, pt)) {
                continue;
              }
              break;
            case 1: // Minimize DCA
              if (!propagateToDCA(nbar, tr, cpa, dca, vtxV0, m, pt)) {
                continue;
              }
              break;
            default: // do nothing
              continue;
          }
          if (tr.mom.getCharge2Pt() > 0) {
            hMiPPDCA->Fill(m, pt, dca);
            hMiPPCPA->Fill(m, pt, cpa);
            if (dca < mDCAcut && cpa > mCPAcut) {
              hMiPP->Fill(m, pt);
            }
          } else {
            hMiPMDCA->Fill(m, pt, dca);
            hMiPMCPA->Fill(m, pt, cpa);
            if (dca < mDCAcut && cpa > mCPAcut) {
              hMiPM->Fill(m, pt);
            }
          }
        }
      }
    }
    for (auto trMixEvent : mixTrackEvts[mixIndex]) {
      for (auto tr : trMixEvent) {
        for (auto nbar : nbarEvent) {
          double dca = 999.;
          switch (mPairingMethod) {
            case 0: // maximize CPA
              if (!minimizeCPA(nbar, tr, cpa, vtxV0, m, pt)) {
                continue;
              }
              break;
            case 1: // Minimize DCA
              if (!propagateToDCA(nbar, tr, cpa, dca, vtxV0, m, pt)) {
                continue;
              }
              break;
            default: // do nothing
              continue;
          }
          if (tr.mom.getCharge2Pt() > 0) {
            hMiPPDCA->Fill(m, pt, dca);
            hMiPPCPA->Fill(m, pt, cpa);
            if (dca < mDCAcut && cpa > mCPAcut) {
              hMiPP->Fill(m, pt);
            }
          } else {
            hMiPMDCA->Fill(m, pt, dca);
            hMiPMCPA->Fill(m, pt, cpa);
            if (dca < mDCAcut && cpa > mCPAcut) {
              hMiPM->Fill(m, pt);
            }
          }
        }
      }
    }
    // Fill events to store and remove oldest to keep buffer size
    if (piEvent.size() > 0 && nbarEvent.size() > 0) {
      mixTrackEvts[mixIndex].emplace_back(piEvent);
      if (mixTrackEvts[mixIndex].size() > mNmix) {
        mixTrackEvts[mixIndex].pop_front();
      }
      mixNbarEvts[mixIndex].emplace_back(nbarEvent);
      if (mixNbarEvts[mixIndex].size() > mNmix) {
        mixNbarEvts[mixIndex].pop_front();
      }
    }
    nbarEvent.clear();
    piEvent.clear();
  }

  template <bool isMC, typename TCollision>
  bool selectEvent(TCollision const& col, int& indx)
  {
    mHistManager.fill(HIST("evsel"), 0.);
    mVtxZ = col.posZ();
    mHistManager.fill(HIST("vtxZ"), mVtxZ);
    if (std::abs(mVtxZ) > 10.f) {
      return false;
    }
    mHistManager.fill(HIST("evsel"), 1.);
    if constexpr (!isMC) {
      if (!col.alias_bit(mEvSelTrig)) {
        return false;
      }
    }
    // Remove pileup???
    mHistManager.fill(HIST("evsel"), 2.);

    // so far only binning according to zvtx is implemented
    indx = (mVtxZ + 10.) / 20. * mNZbins;
    if (indx >= mNZbins) {
      indx = mNZbins - 1.;
    }
    return true;
  }
  //----------------------------------------
  template <bool isMC, typename TClusters>
  void selectNbars(TClusters const& clusters, aod::McParticles const* mcParticles)
  {
    // Select clusters produced by nbar and prepare list of nbar candidates
    for (const auto& clu : clusters) {
      bool isNbar = false;
      int label = -1; // if no MC
      if constexpr (isMC) {
        auto mcList = clu.labels(); // const std::vector<int>
        if (mcList.size() > 0) {
          label = mcList[0];
        }
        int iparent = label;
        while (iparent > -1) {
          auto parent = mcParticles->iteratorAt(iparent);
          if (parent.pdgCode() == -2112) {
            isNbar = true;
            break;
          }
          if (parent.mothersIds().size() == 0 || abs(parent.pdgCode()) < 22 || abs(parent.pdgCode()) > 5000) { // no parents, parent not quark/gluon, strings
            break;
          }
          iparent = parent.mothersIds()[0];
        }
      }

      // first simple cuts
      mHistManager.fill(HIST("cluCuts"), clu.e(), 0.);
      if (isNbar) {
        mHistManager.fill(HIST("nbarCuts"), clu.e(), 0.);
      }
      if (clu.e() < mMinCluE) {
        continue;
      }
      mHistManager.fill(HIST("cluCuts"), clu.e(), 1.);
      if (isNbar) {
        mHistManager.fill(HIST("nbarCuts"), clu.e(), 1.);
      }

      if (clu.time() < mTimeMin || clu.time() > mTimeMax) {
        continue;
      }
      mHistManager.fill(HIST("cluCuts"), clu.e(), 2.);
      if (isNbar) {
        mHistManager.fill(HIST("nbarCuts"), clu.e(), 2.);
      }
      mHistManager.fill(HIST("cluRveto"), clu.trackdist(), clu.e());
      if (isNbar) {
        mHistManager.fill(HIST("nbarRveto"), clu.trackdist(), clu.e());
      }
      mHistManager.fill(HIST("cluTime"), clu.time(), clu.e());
      if (isNbar) {
        mHistManager.fill(HIST("nbarTime"), clu.time(), clu.e());
      }
      // estimate momentum
      float t = clu.time();
      float r = sqrt(clu.globalx() * clu.globalx() + clu.globaly() * clu.globaly() + (clu.globalz() - mVtxZ) * (clu.globalz() - mVtxZ));
      float tgamma = r / c;
      if constexpr (isMC) {
        // smear time
        double sigt = 2.e-9; // TODO: realistic resolution RealRes(cluE);
        t = gRandom->Gaus(t, sigt);
      } else {
        // Real data calibrated wrt photon arrival
        t += tgamma;
      }
      // estimate momentum from time
      if (t <= tgamma) // measured time smaller than photon one
      {
        continue;
      }
      double mom = mNbar / std::sqrt(std::pow(t * c / r, 2) - 1.);
      nbarEvent.emplace_back(mom, clu.globalx(), clu.globaly(), clu.globalz(), label);
      // set PID cuts
      if (clu.m02() > 0.2) { // standard exotics cut
        nbarEvent.back().setPIDBit(1);
        mHistManager.fill(HIST("cluCuts"), clu.e(), 3.);
        if (isNbar) {
          mHistManager.fill(HIST("nbarCuts"), clu.e(), 3.);
        }
      }
      if (clu.ncell() >= mNCellMin) {
        nbarEvent.back().setPIDBit(2);
        mHistManager.fill(HIST("cluCuts"), clu.e(), 4.);
        if (isNbar) {
          mHistManager.fill(HIST("nbarCuts"), clu.e(), 4.);
        }
      }
      if (clu.trackdist() > mRvetoSigma) // neutrality
      {
        nbarEvent.back().setPIDBit(3);
        mHistManager.fill(HIST("cluCuts"), clu.e(), 5.);
        if (isNbar) {
          mHistManager.fill(HIST("nbarCuts"), clu.e(), 5.);
        }
      }
      if (clu.m02() > mDispA * clu.m20() + mDispB) {
        nbarEvent.back().setPIDBit(4);
        mHistManager.fill(HIST("cluCuts"), clu.e(), 6.);
        if (isNbar) {
          mHistManager.fill(HIST("nbarCuts"), clu.e(), 6.);
        }
      }
    }
  } // selectNbars

  //----------------------------------------
  template <bool isMC, typename TTracks>
  void selectTracks(TTracks const& tracks)
  {
    // Select pion tracks to pair with nbar

    for (const auto& piontrack : tracks) {
      mHistManager.fill(HIST("trackCuts"), piontrack.pt(), 0.);
      if (piontrack.hasITS() == false) {
        continue;
      }
      mHistManager.fill(HIST("trackCuts"), piontrack.pt(), 1.);
      if (piontrack.hasTPC() == false) {
        continue;
      }
      mHistManager.fill(HIST("trackCuts"), piontrack.pt(), 2.);
      int label = -1;
      if constexpr (isMC) {
        label = piontrack.mcParticleId();
      }
      piEvent.emplace_back(getTrackPar(piontrack), label);

      if (std::abs(piontrack.tpcNSigmaPi()) < mPionDeDxCut) {
        piEvent.back().setPIDBit(0);
        mHistManager.fill(HIST("trackCuts"), piontrack.pt(), 2.);
      }
      // DCA cut???
    }
  }

  //----------------------------------------
  int commonParentPDG(int labPi, int labNbar, aod::McParticles const* mcParticles)
  {
    // Tests if two labels contain common ancestor
    // return 0 if no common parent
    //        PDG code if found
    int iparentPi = labPi;
    while (iparentPi > -1) {
      int iparentN = labNbar;
      while (iparentN > -1) {
        if (iparentPi == iparentN) {
          return mcParticles->iteratorAt(iparentPi).pdgCode();
        }
        auto parentN = mcParticles->iteratorAt(iparentN);
        if (!parentN.has_mothers() || parentN.pdgCode() == 21 || abs(parentN.pdgCode()) < 11 || abs(parentN.pdgCode()) > 5000) { // no parents, parent not quark/gluon, strings
          break;
        }
        iparentN = parentN.mothersIds()[0];
      }
      auto parentPi = mcParticles->iteratorAt(iparentPi);
      if (!parentPi.has_mothers() || parentPi.pdgCode() == 21 || abs(parentPi.pdgCode()) < 11 || abs(parentPi.pdgCode()) > 5000) { // no parents, parent not quark/gluon, strings
        break;
      }
      iparentPi = parentPi.mothersIds()[0];
    }
    return 0; // nothing found
  }

  //----------------------------------------
  bool minimizeCPA(nbar& n, pion& pn, double& cpa, math_utils::Point3D<float>& vtxV0, double& m, double& pt)
  {
    //--------------------------------------------------------------------
    // Function finds optimal decay vertex where CPA (cosine of Pointing Angle) is maximal
    // we assume that nbar is created at vtxV0 and has modulus of momentum measured in PHOS
    // primary vertex
    math_utils::Vector3D<float> vtxPrim(0., 0., mVtxZ);

    // Step along track
    const int npoints = 5;
    const double eps = 0.2; // accuracy of vertex reconstruction
    double xMax = 10., xMin = -10.;
    std::pair<double, double> st[npoints];
    // operate on copy
    o2::track::TrackParametrization<float> tmpT(pn.mom);
    // Maximal number of minimization steps
    int nsteps = 20;
    while (nsteps && xMax - xMin > eps) {
      double dx = (xMax - xMin) / npoints;
      for (int i = 0; i < npoints; i++) {
        double t = xMin + dx * i;
        tmpT.propagateTo(t, mBz);
        vtxV0 = tmpT.getXYZGlo();
        float ptp = tmpT.getPt();
        float cs = cosf(tmpT.getAlpha()), sn = sinf(tmpT.getAlpha());
        float rp = std::sqrt((1.f - tmpT.getSnp()) * (1.f + tmpT.getSnp()));
        math_utils::Vector3D<float> pPi(ptp * (rp * cs - tmpT.getSnp() * sn), ptp * (tmpT.getSnp() * cs + rp * sn), ptp * tmpT.getTgl());

        // recalculate Nbar momentum
        math_utils::Vector3D<float> pNbar(n.x - vtxV0.x(), n.y - vtxV0.y(), n.z - vtxV0.z());
        // assume that total nbar momentum does not change if V0 is shifted wrt primary vertex by ~Sigma lifetime
        pNbar *= n.mom / sqrt(pNbar.Mag2());

        math_utils::Vector3D<float> pSum = pPi + pNbar;
        double denom = sqrt(pSum.Mag2() * (vtxV0 - vtxPrim).Mag2());
        if (denom > 0) {
          cpa = (vtxV0 - vtxPrim).Dot(pSum) / denom;
        } else { // in primary vertex, step off a bit
          cpa = -1.;
        }
        st[i].first = cpa;
        st[i].second = t;
      }
      // Find optimal region
      double maxCPA = -1;
      int iMax = -1;
      for (int i = 0; i < npoints; i++) {
        if (st[i].first > maxCPA) {
          maxCPA = st[i].first;
          iMax = i;
        }
      }
      if (iMax > 0) {
        if (iMax < npoints - 1) {
          xMin = st[iMax - 1].second; // to account step-like functions
          xMax = st[iMax + 1].second;
        } else {
          xMin = st[npoints - 2].second;
          xMax = 2. * st[npoints - 1].second - st[npoints - 2].second;
        }
      } else {
        xMin = 2. * st[0].second - st[1].second;
        xMax = st[1].second;
      }
      nsteps--;
    }
    tmpT.propagateTo(0.5 * (xMin + xMax), mBz);
    // calculate final CPA, Sigma inv mass using find pion and nbar momenta
    // may be re-calculate nbar |p| using reconstructed vertex
    vtxV0 = tmpT.getXYZGlo();
    float ptp = tmpT.getPt();
    float cs = cosf(tmpT.getAlpha()), sn = sinf(tmpT.getAlpha());
    float rp = std::sqrt((1.f - tmpT.getSnp()) * (1.f + tmpT.getSnp()));
    math_utils::Vector3D<float> pPi(ptp * (rp * cs - tmpT.getSnp() * sn), ptp * (tmpT.getSnp() * cs + rp * sn), ptp * tmpT.getTgl());

    // recalculate Nbar momentum
    math_utils::Vector3D<float> pNbar(n.x - vtxV0.x(), n.y - vtxV0.y(), n.z - vtxV0.z());
    // assume that total nbar momentum does not change if V0 is shifted wrt primary vertex by ~Sigma lifetime
    pNbar *= n.mom / sqrt(pNbar.Mag2());

    math_utils::Vector3D<float> pSum = pPi + pNbar;
    double denom = sqrt(pSum.Mag2() * (vtxV0 - vtxPrim).Mag2());
    if (denom > 0) {
      cpa = (vtxV0 - vtxPrim).Dot(pSum) / denom;
    } else { // in primary vertex, step off a bit
      cpa = -1.;
    }
    // m^2 = (E1+E2)^2 - (p1+p2)^2 =
    double Epi2 = pPi.mag2() + mpi * mpi;
    double En = pNbar.mag2() + mNbar * mNbar;
    m = Epi2 + En + 2. * sqrt(Epi2 * En) - pSum.mag2();
    if (m > 0)
      m = sqrt(m);
    pt = pSum.Rho();
    // asymalpha = (pNbar.mag()-pPi.mag())/(pNbar.mag()+pPi.mag());
    return true;
  }

  bool propagateToDCA(nbar& n, pion& pn,
                      double& cpa, double& dca, math_utils::Point3D<float>& vtxV0, double& m, double& pt)
  {
    //--------------------------------------------------------------------
    // This function returns the DCA between the neutron and the track
    // Proparate track to primary vertex,
    // rotate perp. to nbar and propagate to minimal DCA vertex (in fact it will minimize distance in x)
    const double k = (sqrt(5.) - 1.) / 2.; // 0.61803399
    const double minX = -10., maxX = 10.;
    const double eps = 1.e-3;
    double step = 1.; // initial step 1 cm

    o2::math_utils::Point3D<float> vtxPrim(0., 0., mVtxZ); // Primary vertex
    // unit vector in nbar direction
    o2::math_utils::Vector3D<float> pNbar(n.x - vtxPrim.x(), n.y - vtxPrim.y(), n.z - vtxPrim.z());
    pNbar /= sqrt(pNbar.Mag2()); // unit vector

    auto tmpT(pn.mom); // operate on the copy
    // 1: find most probable position
    double phiN = pNbar.phi();
    double xv = -vtxPrim.x() * sin(phiN) - vtxPrim.y() * cos(phiN); // Rotation by phiN+pi/2
    tmpT.rotateParam(phiN + 0.5 * o2::math_utils::pi());            // If rotation failed, it only increases number of iterations
    if (!tmpT.propagateParamTo(xv, mBz)) {                          // Failed to propagate???
      return false;
    }
    vtxV0 = tmpT.getXYZGlo();
    dca = sqrt((vtxV0 - vtxPrim).Cross(pNbar).Mag2());

    // 2: DCA(t) has V-shape. Find range around estimated minimum which has larger DCA
    double x[4] = {xv - k * step, xv + (1. - 2. * k) * step, xv, xv + (1. - k) * step};
    double dcaP[4] = {dca, dca, dca, dca};
    for (int i = 0; i < 4; i++) {
      if (i == 2)
        continue;
      if (!tmpT.propagateParamTo(x[i], mBz)) {
        return false;
      }
      vtxV0 = tmpT.getXYZGlo();
      dcaP[i] = sqrt((vtxV0 - vtxPrim).Cross(pNbar).Mag2());
    }
    while (dcaP[0] < dcaP[1]) { // too narrow range, should be expanded/shifted
      x[2] = x[1];
      dcaP[2] = dcaP[1];
      x[1] = x[0];
      dcaP[1] = dcaP[0];
      x[0] = x[3] - (x[3] - x[0]) / k;
      if (x[0] < minX) {
        return false; // DCA too far from vertex
      }
      if (!tmpT.propagateParamTo(x[0], mBz)) {
        return false;
      }
      vtxV0 = tmpT.getXYZGlo();
      dcaP[0] = sqrt((vtxV0 - vtxPrim).Cross(pNbar).Mag2());
    }
    // same to the right
    while (dcaP[3] < dcaP[2]) {
      x[1] = x[2];
      dcaP[1] = dcaP[2];
      x[2] = x[3];
      dcaP[2] = dcaP[3];
      x[3] = x[0] + (x[3] - x[0]) / k;
      if (x[3] > maxX) {
        return false; // DCA too far from vertex
      }
      if (!tmpT.propagateParamTo(x[3], mBz)) {
        return false;
      }
      vtxV0 = tmpT.getXYZGlo();
      dcaP[3] = sqrt((vtxV0 - vtxPrim).Cross(pNbar).Mag2());
    }

    // 3: minimize using golden section
    math_utils::Point3D<float> vtx[4];
    while (x[2] - x[1] > eps) {
      if (dcaP[1] < dcaP[2]) {
        x[3] = x[2];
        dcaP[3] = dcaP[2];
        x[2] = x[1];
        dcaP[2] = dcaP[1];
        x[1] = x[3] - k * (x[3] - x[0]);
        if (!tmpT.propagateParamTo(x[1], mBz)) {
          return false;
        }
        vtx[1] = tmpT.getXYZGlo();
        dcaP[1] = sqrt((vtx[1] - vtxPrim).Cross(pNbar).Mag2());
      } else {
        x[0] = x[1];
        dcaP[0] = dcaP[1];
        x[1] = x[2];
        dcaP[1] = dcaP[2];
        x[2] = x[0] + k * (x[3] - x[0]);
        if (!tmpT.propagateParamTo(x[2], mBz)) {
          return false;
        }
        vtx[2] = tmpT.getXYZGlo();
        dcaP[2] = sqrt((vtx[2] - vtxPrim).Cross(pNbar).Mag2());
      }
    }
    if (dcaP[1] < dcaP[2]) {
      vtxV0 = vtx[1];
      dca = dcaP[1];
    } else {
      vtxV0 = vtx[2];
      dca = dcaP[2];
    }

    // calculate topological and kinematic parameters
    //  CPA, Sigma inv mass using find pion and nbar momenta
    //  may be re-calculate nbar |p| using reconstructed vertex
    //  tmpT was extrapolated last time to vtx
    float ptp = tmpT.getPt();
    float cs = cosf(tmpT.getAlpha()), sn = sinf(tmpT.getAlpha());
    float rp = std::sqrt((1.f - tmpT.getSnp()) * (1.f + tmpT.getSnp()));
    math_utils::Vector3D<float> pPi(ptp * (rp * cs - tmpT.getSnp() * sn), ptp * (tmpT.getSnp() * cs + rp * sn), ptp * tmpT.getTgl());

    // recalculate Nbar momentum
    pNbar *= n.mom;

    math_utils::Vector3D<float> pSum = pPi + pNbar;
    double denom = sqrt(pSum.Mag2() * (vtxV0 - vtxPrim).Mag2());
    if (denom > 0) {
      cpa = (vtxV0 - vtxPrim).Dot(pSum) / denom;
    } else { // in primary vertex, step off a bit
      cpa = -1.;
    }
    // m^2 = (E1+E2)^2 - (p1+p2)^2 =
    double Epi2 = pPi.mag2() + mpi * mpi;
    double En = pNbar.mag2() + mNbar * mNbar;
    m = Epi2 + En + 2. * sqrt(Epi2 * En) - pSum.mag2();
    if (m > 0)
      m = sqrt(m);
    pt = pSum.Rho();
    // asymalpha = (pNbar.mag()-pPi.mag())/(pNbar.mag()+pPi.mag());
    return true;
  }

  void processData(SelCollision const& coll,
                   aod::BCsWithTimestamps const&, aod::CaloClusters const& clusters, TrackCandidates const& tracks)
  {
    // Initialize B-field
    if (mBz == 123456.) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", coll.bc_as<aod::BCsWithTimestamps>().timestamp());
      base::Propagator::initFieldFromGRP(grpo);
      mBz = base::Propagator::Instance()->getNominalBz();
    }
    processAll<false>(coll, tracks, clusters, static_cast<aod::McParticles const*>(nullptr));
  }
  PROCESS_SWITCH(phosNbar, processData, "process data", false);

  void processMc(SelCollision const& coll,
                 aod::BCsWithTimestamps const&, mcClusters const& clusters, mcTracks const& tracks, aod::McParticles const& mcPart)
  {
    // Initialize B-field
    if (mBz == 123456.) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", coll.bc_as<aod::BCsWithTimestamps>().timestamp());
      base::Propagator::initFieldFromGRP(grpo);
      mBz = base::Propagator::Instance()->getNominalBz();
    }
    processAll<true>(coll, tracks, clusters, &mcPart);
  }
  PROCESS_SWITCH(phosNbar, processMc, "process MC", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosNbar>(cfgc)};
}
