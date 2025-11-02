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

/// \file photonChargedTriggerCorrelation.cxx
/// \author Julius Kinner
/// \brief photon-jet angular correlation analysis
///
/// Analysis for angular correlations between jets and photons via two-particle correlations with charged high-pt triggers
/// Associated hadrons (tracks), pipm, photons (PCM), pi0 (PCM)
/// Also contains checks and monte-carlo (efficiency, purity, mc-true correlation,...)
/// End goal of studying correlations between direct photons and jets

#include "PWGJE/DataModel/PhotonChargedTriggerCorrelation.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/TableHelper.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TMath.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <memory>
#include <random>
#include <string>
#include <vector>

const double absEtaMaxDefault = 0.8;
#define DPHI_SCALE constants::math::TwoPI - constants::math::PIHalf
#define DETA_SCALE 4 * absEtaMaxDefault - 2 * absEtaMaxDefault

using namespace o2;
using namespace o2::framework;

using CorrCollisions = soa::Join<aod::JetCollisions, aod::CollisionsExtraCorr>;
using CorrCollision = CorrCollisions::iterator;
using CorrMcDCollisions = soa::Join<aod::JetCollisionsMCD, aod::CollisionsExtraCorr>;
using CorrMcDCollision = CorrMcDCollisions::iterator;
using CorrMcCollisions = soa::Join<aod::JetMcCollisions, aod::McCollisionsExtraCorr>;
using CorrMcCollision = CorrMcCollisions::iterator;

using BinningZPvMult = ColumnBinningPolicy<aod::jcollision::PosZ, aod::collision_extra_corr::NGlobalTracks>;

// correlation analysis =======================================================================================================================================================================

struct PhotonChargedTriggerCorrelation {
  // configurables

  // general (kenobi)
  Configurable<std::string> pathCcdbEff{"pathCcdbEff", "Users/j/jkinner/efficiency/set_in_config", "base path to the ccdb efficiencies"};
  Configurable<std::string> urlCcdb{"urlCcdb", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> noLaterThanCcdbConfig{"noLaterThanCcdbConfig", -1, "latest acceptable timestamp of creation for the object (-1 for task start time)"};

  // analysis
  Configurable<bool> doEffCorrectionTrigger{"doEffCorrectionTrigger", false, "whether to do on-the-fly mixing correction for triggers"};
  Configurable<bool> doEffCorrectionHadron{"doEffCorrectionHadron", false, "whether to do on-the-fly mixing correction for hadrons"};
  Configurable<bool> doEffCorrectionPipm{"doEffCorrectionPipm", false, "whether to do on-the-fly mixing correction for pipm"};
  Configurable<bool> doEffCorrectionPhotonPCM{"doEffCorrectionPhotonPCM", false, "whether to do on-the-fly mixing correction for photonPCM"};

  Configurable<bool> doTrigEvMixing{"doTrigEvMixing", false, "whether to use trigger events for trigger mixing"};
  Configurable<int> nTriggerSavedForMixing{"nTriggerSavedForMixing", 2048, "number of triggers that are saved for mixing with other events"};
  Configurable<int> nTriggerMixingMcTrue{"nTriggerMixingMcTrue", 8, "number of triggers that are used for mc true mixing"};
  Configurable<int> nTriggerMixingHadron{"nTriggerMixingHadron", 64, "number of triggers that are used for hadron mixing"};
  Configurable<int> nTriggerMixingPipm{"nTriggerMixingPipm", 64, "number of triggers that are used for pipm mixing"};
  Configurable<int> nTriggerMixingPhotonPCM{"nTriggerMixingPhotonPCM", 256, "number of triggers that are saved for photonPCM mixing"};
  Configurable<int> nTriggerMixingH0PCM{"nTriggerMixingH0PCM", 256, "number of triggers that are saved for h0PCM (pi0, eta) mixing"};
  Configurable<int> nNeighboursMixingPhotonPCMPair{"nNeighboursMixingPhotonPCMPair", 32, "number neighbours used for for photonPCM pair mixing"};
  Configurable<std::vector<double>> pi0PCMPeakMassRange{"pi0PCMPeakMassRange", {0.10, 0.15}, "photon-pair mass integration range for pi0PCM"};
  Configurable<std::vector<double>> pi0PCMSideMassRange{"pi0PCMSideMassRange", {0.16, 0.24}, "photon-pair mass integration range outside pi0PCM region"};
  Configurable<std::vector<double>> etaPCMPeakMassRange{"etaPCMPeakMassRange", {0.51, 0.56}, "photon-pair mass integration range for etaPCM"};
  Configurable<std::vector<double>> etaPCMLowSideMassRange{"etaPCMLowSideMassRange", {0.45, 0.50}, "photon-pair mass integration range below etaPCM region"};
  Configurable<std::vector<double>> etaPCMHighSideMassRange{"etaPCMHighSideMassRange", {0.56, 0.65}, "photon-pair mass integration range above etaPCM region"};

  Configurable<bool> doTrigEvEff{"doTrigEvEff", false, "whether to use trigger events for efficiency histograms"};
  Configurable<float> ptCutTrigEvEff{"ptCutTrigEvEff", 4, "pT cut for efficieny calculation in trigger events (to avoid trigger bias)"};
  Configurable<bool> requireSingleCollisionPurity{"requireSingleCollisionPurity", true, "whether particle from single chosen MC-col associated to reco-col (else just type/kin match)"};

  // for histograms
  Configurable<int> nBinsZPv{"nBinsZPv", 100, "number zPv bins in histos for QA"};
  Configurable<int> nBinsZPvSmol{"nBinsZPvSmol", 28, "number zPv bins but smaller"};
  Configurable<int> nBinsMult{"nBinsMult", 200, "number multiplicity bins in histos for QA"};
  Configurable<int> nBinsMultSmol{"nBinsMultSmol", 20, "number multiplicity bins but smaller"};
  Configurable<int> nBinsOccupancy{"nBinsOccupancy", 2000, "number occupancy bins in histos for QA"};

  Configurable<int> nBinsPhi{"nBinsPhi", 72, "number phi bins"};
  Configurable<int> nBinsEta{"nBinsEta", 40, "number eta bins"};
  Configurable<int> nBinsMgg{"nBinsMgg", 160, "number mass-photon-pair bins"};

  Configurable<std::vector<double>> binsPtTrig{"binsPtTrig", {5, 10, 25, 50}, "correlation ptTrig bins"};
  Configurable<std::vector<double>> binsPtAssoc{"binsPtAssoc",
                                                {0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10, 12.5, 15, 17.5, 20, 30, 40},
                                                "correlation ptAssoc bins"};
  Configurable<std::vector<double>> binsDPhi{"binsDPhi",
                                             {0.00 * DPHI_SCALE,
                                              0.04 * DPHI_SCALE, 0.08 * DPHI_SCALE, 0.11 * DPHI_SCALE, 0.14 * DPHI_SCALE,
                                              0.16 * DPHI_SCALE, 0.18 * DPHI_SCALE, 0.20 * DPHI_SCALE, 0.22 * DPHI_SCALE,
                                              0.23 * DPHI_SCALE, 0.24 * DPHI_SCALE, 0.25 * DPHI_SCALE, 0.26 * DPHI_SCALE, 0.27 * DPHI_SCALE, 0.28 * DPHI_SCALE,
                                              0.30 * DPHI_SCALE, 0.32 * DPHI_SCALE, 0.34 * DPHI_SCALE, 0.36 * DPHI_SCALE,
                                              0.39 * DPHI_SCALE, 0.42 * DPHI_SCALE, 0.46 * DPHI_SCALE, 0.50 * DPHI_SCALE,
                                              0.54 * DPHI_SCALE, 0.58 * DPHI_SCALE, 0.61 * DPHI_SCALE, 0.64 * DPHI_SCALE,
                                              0.66 * DPHI_SCALE, 0.68 * DPHI_SCALE, 0.70 * DPHI_SCALE, 0.72 * DPHI_SCALE,
                                              0.74 * DPHI_SCALE, 0.76 * DPHI_SCALE, 0.78 * DPHI_SCALE,
                                              0.80 * DPHI_SCALE, 0.82 * DPHI_SCALE, 0.84 * DPHI_SCALE, 0.86 * DPHI_SCALE,
                                              0.89 * DPHI_SCALE, 0.92 * DPHI_SCALE, 0.96 * DPHI_SCALE, 1.00 * DPHI_SCALE},
                                             "correlation bins DeltaPhi"};
  Configurable<std::vector<double>> binsDEta{"binsDEta",
                                             {0 / 32. * DETA_SCALE,
                                              1 / 32. * DETA_SCALE, 2 / 32. * DETA_SCALE, 3 / 32. * DETA_SCALE, 4 / 32. * DETA_SCALE,
                                              5 / 32. * DETA_SCALE, 6 / 32. * DETA_SCALE, 7 / 32. * DETA_SCALE, 8 / 32. * DETA_SCALE,
                                              9 / 32. * DETA_SCALE, 10 / 32. * DETA_SCALE, 11 / 32. * DETA_SCALE, 12 / 32. * DETA_SCALE, 13 / 32. * DETA_SCALE, 14 / 32. * DETA_SCALE,
                                              59 / 128. * DETA_SCALE, 62 / 128. * DETA_SCALE, 64 / 128. * DETA_SCALE, 66 / 128. * DETA_SCALE, 69 / 128. * DETA_SCALE, 18 / 32. * DETA_SCALE,
                                              19 / 32. * DETA_SCALE, 20 / 32. * DETA_SCALE, 21 / 32. * DETA_SCALE, 22 / 32. * DETA_SCALE, 23 / 32. * DETA_SCALE, 24 / 32. * DETA_SCALE,
                                              25 / 32. * DETA_SCALE, 26 / 32. * DETA_SCALE, 27 / 32. * DETA_SCALE, 28 / 32. * DETA_SCALE,
                                              29 / 32. * DETA_SCALE, 30 / 32. * DETA_SCALE, 31 / 32. * DETA_SCALE, 32 / 32. * DETA_SCALE},
                                             "correlation bins DeltaEta"};
  Configurable<std::vector<double>> binsZPv{"binsZPv",
                                            {-7, -5, -3, -1, 1, 3, 5, 7},
                                            "zPv mixing bins"};
  Configurable<std::vector<double>> binsMult{"binsMult",
                                             {-0.5, 9.5, 14.5, 19.5, 25.5, 32},
                                             "multiplicity mixing bins for mc true"};
  Configurable<std::vector<double>> binsZPvMcTrue{"binsZPvMcTrue",
                                                  {-10000, 10000},
                                                  "zPv mixing bins"};
  Configurable<std::vector<double>> binsMultMcTrue{"binsMultMcTrue",
                                                   {-0.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 40.5, 50.5, 64},
                                                   "multiplicity mixing bins for mc true"};

  // configurables from other tasks

  double etaMax;

  // objects to hold histograms
  HistogramRegistry histos{"histogramRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // ccdb calls
  const int64_t noLaterThanCcdb = noLaterThanCcdbConfig == -1 ? std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() : noLaterThanCcdbConfig;
  Service<ccdb::BasicCCDBManager> ccdb;
  // for mc
  Service<framework::O2DatabasePDG> pdg;

  // random number generation
  static constexpr unsigned int SeedRandomEngine = 12345;
  std::mt19937 randomEngine{SeedRandomEngine};

  // partitions
  SliceCache cache;

  // prepare for per collision slicing
  Preslice<aod::JetTracks> perColTracks = aod::jtrack::collisionId;
  Preslice<aod::Triggers> perColTriggers = aod::corr_particle::jetCollisionId;
  Preslice<aod::Hadrons> perColHadrons = aod::corr_particle::jetCollisionId;
  Preslice<aod::Pipms> perColPipms = aod::corr_particle::jetCollisionId;
  Preslice<aod::PhotonPCMs> perColPhotonPCMs = aod::corr_particle::jetCollisionId;
  Preslice<aod::PhotonPCMPairs> perColPhotonPCMPairs = aod::corr_particle::jetCollisionId;
  Preslice<aod::JetParticles> perColMcParticles = aod::jmcparticle::mcCollisionId;
  Preslice<aod::TriggerParticles> perColTriggerParticles = aod::corr_particle::jetMcCollisionId;

  // combinations binning
  // cumbersome, but still better than having extra configurable or figuring out how to init binningZPvMult later while declaring it here
  std::function<std::vector<double>(std::vector<double> const&, double)> prependValueToVector =
    [](std::vector<double> const& vec, double const value) {
      std::vector<double> resultVec = {value};
      resultVec.insert(resultVec.end(), vec.begin(), vec.end());
      return resultVec;
    };
  BinningZPvMult binningZPvMult{{prependValueToVector(binsZPv.value, VARIABLE_WIDTH), prependValueToVector(binsMult.value, VARIABLE_WIDTH)}, true};

  // declare analysis variables

  // efficiency histograms
  TH1D* h1PtInvEffTrigger;
  TH1D* h1PtInvEffHadron;
  TH1D* h1PtInvEffPipm;
  TH1D* h1PtInvEffPhotonPCM;

  // mixing trigger memory
  struct MixingTrigger {
    float fPt, fPhi, fEta;
    float pt() const { return fPt; }
    float phi() const { return fPhi; }
    float eta() const { return fEta; }
  };
  // class to handle trigger info from previous collisions (beyond single dataframe)
  // organised as zPv- and mult-bin matrix of deque (pt, phi, eta) to save trigger info beyond single dataframe
  // extra bin for mult overflow
  // with adjusted zVtx (see triggerBinValuesZPv in init) and mult overflow -> all events accounted for
  // (possibly replace by some advanced derived data method and O2 event mixing in future?)
  class MixingTriggerMemory
  {
   public:
    // finds bin that value belongs to (assumes ordered bins) (starts at 0; includes underflow (return -1) and overlflow (return bins.size() - 1))
    // should be faster than some std binary search due to small number of bins (zPv, mult)
    static int findIntervalBin(double value, const std::vector<double>& bins)
    {
      const int n = bins.size() - 1;
      if (value < bins[0])
        return -1; // underflow
      for (int i_bin = 0; i_bin < n; i_bin++)
        if (value < bins[i_bin + 1])
          return i_bin;
      return n; // overflow
    }

    MixingTriggerMemory(int const nTriggerSavedForMixingIn, std::vector<double> binsZPv, std::vector<double> binsMult)
    {
      nTriggerSavedForMixing = nTriggerSavedForMixingIn;
      triggerBinValuesZPv = binsZPv;
      triggerBinValuesMult = binsMult;
      // prevent rounding errors in bin finding (multiplicity accounted for by it going to 0 and already considering overflow separately)
      triggerBinValuesZPv.front() *= zPvRoundingErrorAdjust;
      triggerBinValuesZPv.back() *= zPvRoundingErrorAdjust;
      // init correct size of zPv-mult matrix
      savedTriggersZPvMult.resize(binsZPv.size() - 1);
      for (size_t i_zPv = 0; i_zPv < binsZPv.size() - 1; i_zPv++) {
        savedTriggersZPvMult[i_zPv].resize(binsMult.size());
      }
    }

    // save trigger for mixing
    // up to nTriggerSavedForMixing stored (LIFO)
    void saveTrigger(float const pt, float const phi, float const eta, double const zPv, double const mult)
    {
      int const iBinCorrZPv = findIntervalBin(zPv, triggerBinValuesZPv);
      int const iBinCorrMult = findIntervalBin(mult, triggerBinValuesMult);
      // special cases (floating point precision errors, mult overflow) should be taken care of by triggerBinValuesZPv and triggerBinValuesMult
      savedTriggersZPvMult[iBinCorrZPv][iBinCorrMult].push_front(MixingTrigger{pt, phi, eta});
      if (static_cast<int>(savedTriggersZPvMult[iBinCorrZPv][iBinCorrMult].size()) > nTriggerSavedForMixing) {
        savedTriggersZPvMult[iBinCorrZPv][iBinCorrMult].pop_back();
      }
    }

    // return deques of trigger pt, phi, eta in the given zPv/mult bin
    std::deque<MixingTrigger> const& getTriggers(double const zPv, double const mult) const
    {
      int const iBinCorrZPv = findIntervalBin(zPv, triggerBinValuesZPv);
      int const iBinCorrMult = findIntervalBin(mult, triggerBinValuesMult);
      return savedTriggersZPvMult[iBinCorrZPv][iBinCorrMult];
    }

   private:
    double const zPvRoundingErrorAdjust = 1.0001;
    int nTriggerSavedForMixing;
    std::vector<double> triggerBinValuesZPv;
    std::vector<double> triggerBinValuesMult;
    std::vector<std::vector<std::deque<MixingTrigger>>> savedTriggersZPvMult;
  };

  MixingTriggerMemory mixingTriggerMemoryReco{nTriggerSavedForMixing.value, binsZPv.value, binsMult.value};
  MixingTriggerMemory mixingTriggerMemoryTrue{nTriggerSavedForMixing.value, binsZPvMcTrue.value, binsMult.value};
  MixingTriggerMemory mixingTriggerMemoryRecoColTrue{nTriggerSavedForMixing.value, binsZPvMcTrue.value, binsMult.value};

  // functions ================================================================================================================================================================================

  // general (kenobi) /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // get histograms from ccdb
  // save efficiencies from ccdb in histogram registry
  void initCcdbHistograms()
  {
    // trigger
    h1PtInvEffTrigger = nullptr;
    if (doEffCorrectionTrigger) {
      h1PtInvEffTrigger = ccdb->getForTimeStamp<TH1D>(pathCcdbEff.value + "/trigger", noLaterThanCcdb);

      const double* effBinsTrigger = h1PtInvEffTrigger->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtEffTrigger{std::vector<double>(effBinsTrigger, effBinsTrigger + h1PtInvEffTrigger->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedEff/h1_pt_invEff_trigger_ccdb", "h1_pt_invEff_trigger_ccdb", kTH1D, {axisPtEffTrigger}, true);
      for (int iBin = 1; iBin <= h1PtInvEffTrigger->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_trigger_ccdb"))->SetBinContent(iBin, h1PtInvEffTrigger->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_trigger_ccdb"))->SetBinError(iBin, h1PtInvEffTrigger->GetBinError(iBin));
      }
    }
    // hadron
    h1PtInvEffHadron = nullptr;
    if (doEffCorrectionHadron) {
      h1PtInvEffHadron = ccdb->getForTimeStamp<TH1D>(pathCcdbEff.value + "/hadron", noLaterThanCcdb);

      const double* effBinsHadron = h1PtInvEffHadron->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtEffHadron{std::vector<double>(effBinsHadron, effBinsHadron + h1PtInvEffHadron->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedEff/h1_pt_invEff_hadron_ccdb", "h1_pt_invEff_hadron_ccdb", kTH1D, {axisPtEffHadron}, true);
      for (int iBin = 1; iBin <= h1PtInvEffHadron->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_hadron_ccdb"))->SetBinContent(iBin, h1PtInvEffHadron->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_hadron_ccdb"))->SetBinError(iBin, h1PtInvEffHadron->GetBinError(iBin));
      }
    }
    // pipm
    h1PtInvEffPipm = nullptr;
    if (doEffCorrectionPipm) {
      h1PtInvEffPipm = ccdb->getForTimeStamp<TH1D>(pathCcdbEff.value + "/pipm", noLaterThanCcdb);

      const double* effBinsPipm = h1PtInvEffPipm->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtEffPipm{std::vector<double>(effBinsPipm, effBinsPipm + h1PtInvEffPipm->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedEff/h1_pt_invEff_pipm_ccdb", "h1_pt_invEff_pipm_ccdb", kTH1D, {axisPtEffPipm}, true);
      for (int iBin = 1; iBin <= h1PtInvEffPipm->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_pipm_ccdb"))->SetBinContent(iBin, h1PtInvEffPipm->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_pipm_ccdb"))->SetBinError(iBin, h1PtInvEffPipm->GetBinError(iBin));
      }
    }
    // photonPCM
    h1PtInvEffPhotonPCM = nullptr;
    if (doEffCorrectionPhotonPCM) {
      h1PtInvEffPhotonPCM = ccdb->getForTimeStamp<TH1D>(pathCcdbEff.value + "/photonPCM", noLaterThanCcdb);

      const double* effBinsPhotonPCM = h1PtInvEffPhotonPCM->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtEffPhotonPCM{std::vector<double>(effBinsPhotonPCM, effBinsPhotonPCM + h1PtInvEffPhotonPCM->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedEff/h1_pt_invEff_photonPCM_ccdb", "h1_pt_invEff_photonPCM_ccdb", kTH1D, {axisPtEffPhotonPCM}, true);
      for (int iBin = 1; iBin <= h1PtInvEffPhotonPCM->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_photonPCM_ccdb"))->SetBinContent(iBin, h1PtInvEffPhotonPCM->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_photonPCM_ccdb"))->SetBinError(iBin, h1PtInvEffPhotonPCM->GetBinError(iBin));
      }
    }
  }

  // create histograms
  void initHistograms()
  {
    // define axes
    const AxisSpec axisN{1, 0., 1., "#it{N}_{something}"};
    const AxisSpec axisCategories{16, 0., 16., "categories"};

    const AxisSpec axisZPv{nBinsZPv, -10, 10, "#it{z}_{pv}"};
    const AxisSpec axisZPvSmol{nBinsZPvSmol, -7, 7, "#it{z}_{pv}"};
    const AxisSpec axisMult{nBinsMult + 1, -0.5, nBinsMult + 0.5, "multiplicity"};
    const AxisSpec axisMultSmol{nBinsMultSmol + 1, -0.5, nBinsMultSmol + 0.5, "multiplicity"};
    const AxisSpec axisOccupancy{nBinsOccupancy + 1, -0.5, nBinsOccupancy + 0.5, "occupancy"};

    const AxisSpec axisPhi{nBinsPhi, 0, constants::math::TwoPI, "#it{#varphi}"};
    const AxisSpec axisEta{nBinsEta, -etaMax, etaMax, "#it{#eta}"};
    const AxisSpec axisMgg{nBinsMgg, 0, 0.8, "#it{m}_{#gamma#gamma}"};

    const AxisSpec axisPtTrig{binsPtTrig, "#it{p}_{T}^{trig}"};
    const AxisSpec axisPtAssoc{binsPtAssoc, "#it{p}_{T}^{assoc}"};
    const AxisSpec axisDPhi{binsDPhi, "#Delta#it{#varphi}"};
    const AxisSpec axisDEta{binsDEta, "#Delta#it{#eta}"};
    const AxisSpec axisZPvBinning{binsZPv, "#it{z}_{pv} correlation binning"};
    const AxisSpec axisMultBinning{binsMult, "multiplicity correlation binning"};
    const AxisSpec axisZPvBinningMcTrue{binsZPvMcTrue, "#it{z}_{pv} correlation binning for mc true"};
    const AxisSpec axisMultBinningMcTrue{binsMultMcTrue, "multiplicity correlation binning for mc true"};

    // reco info
    histos.add("reco/info/h1_nEvents", "h1_nEvents", kTH1D, {axisCategories});
    histos.get<TH1>(HIST("reco/info/h1_nEvents"))->GetXaxis()->SetBinLabel(1, "#it{N}_{ev}^{sel}");
    histos.get<TH1>(HIST("reco/info/h1_nEvents"))->GetXaxis()->SetBinLabel(2, "#it{N}_{ev}");
    histos.get<TH1>(HIST("reco/info/h1_nEvents"))->GetXaxis()->SetBinLabel(3, "#it{N}_{ev}^{trig}");

    histos.add("reco/info/h2_zPvMult", "h2_zPvMult", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("reco/info/h1_occupancy", "h1_occupancy", kTH1D, {axisOccupancy}, true);
    histos.add("reco/info/h2_zPvMult_trigEv", "h2_zPvMult_trigEv", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("reco/info/h1_occupancy_trigEv", "h1_occupancy_trigEv", kTH1D, {axisOccupancy}, true);

    // reco (correlation) analysis
    histos.add("reco/corr/h3_ptPhiEta_trig", "h3_ptPhiEta_trig", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);

    std::function<void(std::string)> add_corrHists =
      [&](std::string const name_id) {
        histos.add(std::format("reco/corr/h3_ptPhiEta_assoc_{}", name_id).data(), std::format("h3_ptPhiEta_assoc_{}", name_id).data(),
                   kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
        histos.add(std::format("reco/corr/h6_corr_{}", name_id).data(), std::format("h6_corr_{}", name_id).data(),
                   kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
        histos.add(std::format("reco/corr/h6_mix_{}", name_id).data(), std::format("h6_mix_{}", name_id).data(),
                   kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
      };
    // hadron
    histos.add("reco/plain/h3_ptPhiEta_hadron", "h3_ptPhiEta_hadron", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    add_corrHists("hadron");
    // pipm
    histos.add("reco/plain/h3_ptPhiEta_pipm", "h3_ptPhiEta_pipm", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    add_corrHists("pipm");
    // photonPCM
    histos.add("reco/plain/h3_ptPhiEta_photonPCM", "h3_ptPhiEta_photonPCM", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    add_corrHists("photonPCM");
    // photonPCM pairs
    histos.add("reco/plain/h4_ptMggZPvMult_photonPCMPair", "h4_ptMggZPvMult_photonPCMPair", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/plain/h4_ptMggZPvMult_trigEv_photonPCMPair", "h4_ptMggZPvMult_trigEv_photonPCMPair", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h5_ptTrigPtAssocMggZPvMult_assoc_photonPCMPair", "h5_ptTrigPtAssocMggZPvMult_assoc_photonPCMPair", kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    // pi0PCM
    add_corrHists("pi0PCMPeak");
    add_corrHists("pi0PCMSide");
    // etaPCM
    add_corrHists("etaPCMPeak");
    add_corrHists("etaPCMSide");

    // event mixing for photon pairs
    histos.add("reco/plain/h2_zPvMult_photonPCMPair_evMix", "h2_zPvMult_photonPCMPair_evMix", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("reco/plain/h4_ptMggZPvMult_photonPCMPair_evMix", "h4_ptMggZPvMult_photonPCMPair_evMix", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/plain/h4_ptMggZPvMult_trigEv_photonPCMPair_evMix", "h4_ptMggZPvMult_trigEv_photonPCMPair_evMix", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);

    // mc info
    histos.add("mc/info/h1_nEvents_mcTrue", "h1_nEvents_mcTrue", kTH1D, {axisN});
    histos.add("mc/info/h2_zPvMult_mcTrue", "h2_zPvMult_mcTrue", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("mc/info/h1_nTrigEv_mcTrue", "h1_nTrigEv_mcTrue", kTH1D, {axisN});
    histos.add("mc/info/h2_zPvMult_trigEv_mcTrue", "h2_zPvMult_trigEv_mcTrue", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("mc/info/h1_nRecoCol_mcTrue", "h1_nRecoCol_mcTrue", kTH1D, {axisN});
    histos.add("mc/info/h2_zPvMult_recoCol_mcTrue", "h2_zPvMult_recoCol_mcTrue", kTHnSparseD, {axisZPv, axisMult}, true);

    // reco and true collision correlations
    const std::vector<std::string> assocMcCorrHistNames = {"hadron", "pipm", "photon", "pi0", "eta"};
    for (auto const& collision_type : {"true", "recoCol_true"}) {
      histos.add(std::format("mc/{}/corr/h3_ptPhiEta_trig", collision_type).data(), "h3_ptPhiEta_trig", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
      for (auto const& assocName : assocMcCorrHistNames) {
        histos.add(std::format("mc/{}/corr/h6_corr_{}", collision_type, assocName).data(), std::format("h6_corr_{}", assocName).data(),
                   kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinningMcTrue, axisMultBinningMcTrue}, true);
        histos.add(std::format("mc/{}/corr/h6_mix_{}", collision_type, assocName).data(), std::format("h6_mix_{}", assocName).data(),
                   kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinningMcTrue, axisMultBinningMcTrue}, true);
      }
    }

    // mc efficiency/purity
    std::function<void(std::string)> add_effHists =
      [&](std::string const name_id) {
        histos.add(std::format("mc/eff/h3_ptPhiEta_{}", name_id).data(), std::format("h3_ptPhiEta_{}", name_id).data(),
                   kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
        histos.add(std::format("mc/eff/h3_ptZPvMult_{}", name_id).data(), std::format("h3_ptZPvMult_{}", name_id).data(),
                   kTHnSparseD, {axisPtAssoc, axisZPvSmol, axisMultSmol}, true);
      };
    // mc tracks
    add_effHists("mcReco_hadron");
    add_effHists("mcReco_hasCorrectMc_hadron");
    add_effHists("mcTrue_recoCol_hadron");
    // mc pipm PID
    add_effHists("mcReco_pipm");
    add_effHists("mcReco_hasCorrectMc_pipm");
    add_effHists("mcTrue_recoCol_pipm");
    // mc photonPCM
    add_effHists("mcReco_photonPCM");
    add_effHists("mcReco_hasCorrectMc_photonPCM");
    add_effHists("mcTrue_recoCol_photon");
    // mc pi0
    add_effHists("mcTrue_recoCol_pi0");
    // mc eta
    add_effHists("mcTrue_recoCol_eta");

    // test of the test while testing another test. featuring a test
    histos.add("test/h2_mult_comp", "h2_mult_comp", kTH2D, {axisMult, axisMult}, true);
    histos.add("test/h2_tracks_zPvMultDep", "h2_tracks_zPvMultDep", kTH2D, {axisZPv, axisMult}, true);
    histos.add("test/h2_globalTracks_zPvMultDep", "h2_globalTracks_zPvMultDep", kTH2D, {axisZPv, axisMult}, true);
  }

  // selections ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // checks if mcParticle is charged
  template <typename T_mcParticle>
  bool checkChargedMc(T_mcParticle const& mcParticle)
  {
    auto const pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0)
      return false;
    return true;
  }
  // checks if mcParticle should be detected (physicalPrimary, |eta|)
  template <typename T_mcParticle>
  bool checkPrimaryEtaMc(T_mcParticle const& mcParticle)
  {
    if (!mcParticle.isPhysicalPrimary())
      return false;
    if (std::abs(mcParticle.eta()) > etaMax)
      return false;
    return true;
  }
  // checks if mcParticle should be detected as primary track (physicalPrimary, charge, |eta|)
  template <typename T_mcParticle>
  bool checkPrimaryTrackMc(T_mcParticle const& mcParticle)
  {
    if (!checkPrimaryEtaMc(mcParticle))
      return false;
    if (!checkChargedMc(mcParticle))
      return false;
    return true;
  }
  // checks if mcParticle should be detected as 'primary' (|eta| not checked)
  template <typename T_mcParticle>
  bool checkH0Primary(T_mcParticle const& mcParticle, int const pdg)
  {
    if (mcParticle.pdgCode() != pdg)
      return false;
    const auto& h0Daughters = mcParticle.template daughters_as<aod::JetParticles>();
    // identify primary h0 (account for 0 daughters for some reason)
    if (h0Daughters.size() == 0)
      return false;
    for (auto const& h0_daughter : h0Daughters) {
      if (!h0_daughter.isPhysicalPrimary())
        return false;
    }
    return true;
  }
  // checks if mcParticle should be detected as 'primary' pi0->gammagamma (|eta| not checked)
  template <typename T_mcParticle>
  bool checkH0ToGG(T_mcParticle const& mcParticle, int const pdg)
  {
    if (!checkH0Primary(mcParticle, pdg))
      return false;
    // select h0 -> gg
    constexpr int NDaughtersH0ToGG = 2;
    if (mcParticle.template daughters_as<aod::JetParticles>().size() != NDaughtersH0ToGG)
      return false;
    return true;
  }

  // checks if tracks come from photon conversion
  template <typename T_track>
  bool isConversionPhoton(T_track const& posTrack, T_track const& negTrack)
  {
    // check same mother
    auto const& posMothers = posTrack.mcParticle().template mothers_as<aod::JetParticles>();
    auto const& negMothers = negTrack.mcParticle().template mothers_as<aod::JetParticles>();
    if (posMothers.size() != 1 || negMothers.size() != 1)
      return false;
    if (posMothers.begin()->globalIndex() != negMothers.begin()->globalIndex())
      return false;
    // check photon
    if (posMothers.begin()->pdgCode() != PDG_t::kGamma)
      return false;

    return true;
  };
  // checks if tracks come from pi0 double conversion
  template <typename T_track>
  bool isGGFromPi0(T_track const& posTrack1, T_track const& negTrack1, T_track const& posTrack2, T_track const& negTrack2)
  {
    if (!isConversionPhoton(posTrack1, negTrack1) || !isConversionPhoton(posTrack2, negTrack2))
      return false;
    // check same mother
    auto const& mothers1 = (*(posTrack1.mcParticle().template mothers_as<aod::JetParticles>().begin())).template mothers_as<aod::JetParticles>();
    auto const& mothers2 = (*(posTrack2.mcParticle().template mothers_as<aod::JetParticles>().begin())).template mothers_as<aod::JetParticles>();
    constexpr int NMothersPhotonFromPi0 = 2; // for some reason two mothers (same particle) for pi0 decays (contradicts PYTHIA documentation, but whatever)
    if (mothers1.size() != NMothersPhotonFromPi0 || mothers2.size() != NMothersPhotonFromPi0)
      return false;
    if (mothers1.begin()->globalIndex() != mothers2.begin()->globalIndex())
      return false;
    // check pi0
    if (mothers1.begin()->pdgCode() != PDG_t::kPi0)
      return false;

    return true;
  };

  // analysis helpers /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <typename T_h1>
  double getH1ValueAt(T_h1 const* const h1, double const value)
  {
    return h1->GetBinContent(h1->FindFixBin(value));
  }
  // efficiency helpers
  enum class EffParticleType { Trigger,
                               Hadron,
                               Pipm,
                               PhotonPCM };
  // efficiency function
  template <EffParticleType T_effParticleType>
  double getInvEff(double const value)
  {
    if constexpr (T_effParticleType == EffParticleType::Trigger) {
      return doEffCorrectionTrigger ? getH1ValueAt(h1PtInvEffTrigger, value) : 1;
    } else if constexpr (T_effParticleType == EffParticleType::Hadron) {
      return doEffCorrectionHadron ? getH1ValueAt(h1PtInvEffHadron, value) : 1;
    } else if constexpr (T_effParticleType == EffParticleType::Pipm) {
      return doEffCorrectionPipm ? getH1ValueAt(h1PtInvEffPipm, value) : 1;
    } else if constexpr (T_effParticleType == EffParticleType::PhotonPCM) {
      return doEffCorrectionPhotonPCM ? getH1ValueAt(h1PtInvEffPhotonPCM, value) : 1;
    } else {
      return 1;
    }
  }

  // performs 'phi1 - phi2' and pushes it into the interval [-pi/2, 3pi/2]
  inline double getDeltaPhi(double const phi1, double const phi2)
  {
    return RecoDecay::constrainAngle(phi1 - phi2, -1 * constants::math::PIHalf);
  }

  // finds bin that value belongs to (assumes ordered bins) (starts at 0; includes underflow (return -1) and overlflow (return bins.size() - 1))
  // should be faster than some std binary search due to small number of bins (zPv, mult)
  int findIntervalBin(double value, const std::vector<double>& bins)
  {
    const int n = bins.size() - 1;
    if (value < bins[0])
      return -1; // underflow
    for (int i_bin = 0; i_bin < n; i_bin++)
      if (value < bins[i_bin + 1])
        return i_bin;
    return n; // overflow
  }

  // checks that two values belong to the same category (assumes ordered bins)
  // returns -1 for negative result (also for under/overflow values) and bin number (starting at 0) otherwise
  int checkSameBin(double const value1, double const value2, std::vector<double> const& bins)
  {
    // reject underflow
    if (value1 < bins[0])
      return -1;
    // loop over bins
    const int n = bins.size() - 1;
    for (int i_bin = 0; i_bin < n; i_bin++) {
      if (value1 < bins[i_bin + 1]) {
        if (value2 < bins[i_bin + 1] && value2 >= bins[i_bin]) {
          return i_bin;
        }
        return -1;
      }
    }
    // reject overflow
    return -1;
  }

  // check if invariant mass range
  enum class MassRange { pi0PCMPeak,
                         pi0PCMSide,
                         etaPCMPeak,
                         etaPCMSide };

  template <MassRange T_massRange>
  bool checkMassRange(double const mgg)
  {
    if constexpr (T_massRange == MassRange::pi0PCMPeak) {
      return mgg > pi0PCMPeakMassRange.value[0] && mgg < pi0PCMPeakMassRange.value[1];
    } else if constexpr (T_massRange == MassRange::pi0PCMSide) {
      return mgg > pi0PCMSideMassRange.value[0] && mgg < pi0PCMSideMassRange.value[1];
    } else if constexpr (T_massRange == MassRange::etaPCMPeak) {
      return mgg > etaPCMPeakMassRange.value[0] && mgg < etaPCMPeakMassRange.value[1];
    } else if constexpr (T_massRange == MassRange::etaPCMSide) {
      return (mgg > etaPCMLowSideMassRange.value[0] && mgg < etaPCMLowSideMassRange.value[1]) ||
             (mgg > etaPCMHighSideMassRange.value[0] && mgg < etaPCMHighSideMassRange.value[1]);
    }
    return false;
  }

  // analysis /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // generalised correlation functions
  // per collision

  // plain info
  template <typename T_collision, typename T_associatedThisEvent,
            typename T_funcPlain>
  void corrProcessPlain(T_collision const& collision, T_associatedThisEvent const& associatedThisEvent,
                        T_funcPlain&& funcPlain)
  {
    // normal spectra (per event - not per trigger)
    for (auto const& associated : associatedThisEvent) {
      funcPlain(collision, associated);
    }
  }

  // correlation
  template <typename T_collision, typename T_triggersThisEvent, typename T_associatedThisEvent,
            typename T_funcCorrelation>
  void corrProcessCorrelation(T_collision const& collision, T_triggersThisEvent const& triggersThisEvent, T_associatedThisEvent const& associatedThisEvent,
                              T_funcCorrelation&& funcCorrelation)
  {
    // correlation combinations
    for (auto const& [trigger, associated] : soa::combinations(soa::CombinationsFullIndexPolicy(triggersThisEvent, associatedThisEvent))) {
      funcCorrelation(collision, trigger, associated);
    }
  }

  // mixing
  template <typename T_collision, typename T_associatedThisEvent,
            typename T_funcMixing>
  void corrProcessMixing(T_collision const& collision, T_associatedThisEvent const& associatedThisEvent,
                         T_funcMixing&& funcMixing,
                         size_t const nTriggerMixing, size_t const nTriggersThisDataFrame)
  {
    // skip if event does not contain valid trigger
    if (doTrigEvMixing && !collision.trigEv())
      return;

    // mixing loops (more efficient than O2 mixing (for now))
    auto savedTriggers = mixingTriggerMemoryReco.getTriggers(collision.posZ(), collision.nGlobalTracks());
    // number of triggers
    const size_t mixUpToTriggerN = std::min(savedTriggers.size(), nTriggerMixing + nTriggersThisDataFrame);
    const float perTriggerWeight = 1. / (mixUpToTriggerN - nTriggersThisDataFrame); // mixUpToTriggerN <= nTriggersThisDataFrame not problematic since no loop then
    // mixing loops
    for (size_t i_mixingTrigger = nTriggersThisDataFrame; i_mixingTrigger < mixUpToTriggerN; i_mixingTrigger++) {
      for (auto const& associated : associatedThisEvent) {
        funcMixing(collision, savedTriggers[i_mixingTrigger].pt(), savedTriggers[i_mixingTrigger].phi(), savedTriggers[i_mixingTrigger].eta(), associated, perTriggerWeight);
      }
    }
  }

  void init(InitContext& initContext)
  {
    // analysis info
    ccdb->setURL(urlCcdb.value);
    // enabling object caching (otherwise each call goes to CCDB server)
    ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // not later than (avoids replacing objects while a train is running)
    ccdb->setCreatedNotAfter(noLaterThanCcdb);

    // init analysis variables

    // get variables from other tasks
    getTaskOptionValue(initContext, "photon-charged-trigger-producer", "etaMax", etaMax, false);

    // histograms from ccdb
    initCcdbHistograms();

    // create analysis histograms
    initHistograms();
  }

  // reconstructed ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void processInfo(CorrCollision const& collision)
  {
    // all events
    histos.fill(HIST("reco/info/h1_nEvents"), 1.5);

    // event selection
    if (!collision.selEv())
      return;
    histos.fill(HIST("reco/info/h1_nEvents"), 0.5);

    histos.fill(HIST("reco/info/h2_zPvMult"), collision.posZ(), collision.nGlobalTracks());
    histos.fill(HIST("reco/info/h1_occupancy"), collision.trackOccupancyInTimeRange());

    // trigger events
    if (!collision.trigEv())
      return;
    histos.fill(HIST("reco/info/h1_nEvents"), 2.5);

    histos.fill(HIST("reco/info/h2_zPvMult_trigEv"), collision.posZ(), collision.nGlobalTracks());
    histos.fill(HIST("reco/info/h1_occupancy_trigEv"), collision.trackOccupancyInTimeRange());
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processInfo, "process general info on collisions and tracks for analysis and qa", false);

  void processCorrFirst(CorrCollisions const& collisions, aod::Triggers const& triggers)
  {
    // do at beginning of each data frame (before other reco correlation process functions)
    // (PROCESS_SWITCH of this process has to be declared first)

    // [wow, such empty]

    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());

      // trigger loop
      for (auto const& trigger : triggersThisEvent) {
        // trigger info
        histos.fill(HIST("reco/corr/h3_ptPhiEta_trig"), trigger.pt(), trigger.phi(), trigger.eta(),
                    getInvEff<EffParticleType::Trigger>(trigger.pt()));

        // save trigger for mixing
        mixingTriggerMemoryReco.saveTrigger(trigger.pt(), trigger.phi(), trigger.eta(), collision.posZ(), collision.nGlobalTracks());
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrFirst, "process to gather info before correlation processes", false);

  void processCorrHadron(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::Hadrons const& hadrons)
  {
    size_t const nTriggersThisDataFrame = triggers.size();

    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const hadronsThisEvent = hadrons.sliceBy(perColHadrons, collision.globalIndex());

      auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h3_ptPhiEta_hadron"),
                    associated.pt(), associated.phi(), associated.eta(),
                    getInvEff<EffParticleType::Hadron>(associated.pt()));
      };
      corrProcessPlain(collision, hadronsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.jetTrackId())
          return;

        histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_hadron"),
                    associated.pt(), associated.phi(), associated.eta(),
                    getInvEff<EffParticleType::Trigger>(trigger.pt()) * getInvEff<EffParticleType::Hadron>(associated.pt()));
        histos.fill(HIST("reco/corr/h6_corr_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getInvEff<EffParticleType::Trigger>(trigger.pt()) * getInvEff<EffParticleType::Hadron>(associated.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, hadronsThisEvent, funcCorrelation);

      auto const funcMixing = [this](auto const& collision,
                                     float const mixingTriggerPt, float const mixingTriggerPhi, float const mixingTriggerEta, auto const& associated, auto const perTriggerWeight) {
        histos.fill(HIST("reco/corr/h6_mix_hadron"),
                    getDeltaPhi(mixingTriggerPhi, associated.phi()),
                    mixingTriggerEta - associated.eta(),
                    mixingTriggerPt, associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getInvEff<EffParticleType::Trigger>(mixingTriggerPt) * getInvEff<EffParticleType::Hadron>(associated.pt()));
      };
      corrProcessMixing(collision, hadronsThisEvent, funcMixing, nTriggerMixingHadron, nTriggersThisDataFrame);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrHadron, "process standard correlation for associated hardons", false);

  void processCorrPipm(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::Pipms const& pipms)
  {
    size_t const nTriggersThisDataFrame = triggers.size();

    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const pipmsThisEvent = pipms.sliceBy(perColPipms, collision.globalIndex());

      auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h3_ptPhiEta_pipm"),
                    associated.pt(), associated.phi(), associated.eta(),
                    getInvEff<EffParticleType::Pipm>(associated.pt()));
      };
      corrProcessPlain(collision, pipmsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.jetTrackId())
          return;

        histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_pipm"),
                    associated.pt(), associated.phi(), associated.eta(),
                    getInvEff<EffParticleType::Trigger>(trigger.pt()) * getInvEff<EffParticleType::Pipm>(associated.pt()));
        histos.fill(HIST("reco/corr/h6_corr_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getInvEff<EffParticleType::Trigger>(trigger.pt()) * getInvEff<EffParticleType::Pipm>(associated.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, pipmsThisEvent, funcCorrelation);

      auto const funcMixing = [this](auto const& collision,
                                     float const mixingTriggerPt, float const mixingTriggerPhi, float const mixingTriggerEta, auto const& associated, auto const perTriggerWeight) {
        histos.fill(HIST("reco/corr/h6_mix_pipm"),
                    getDeltaPhi(mixingTriggerPhi, associated.phi()),
                    mixingTriggerEta - associated.eta(),
                    mixingTriggerPt, associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getInvEff<EffParticleType::Trigger>(mixingTriggerPt) * getInvEff<EffParticleType::Pipm>(associated.pt()));
      };
      corrProcessMixing(collision, pipmsThisEvent, funcMixing, nTriggerMixingPipm, nTriggersThisDataFrame);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPipm, "process standard correlation for associated pipm", false);

  void processCorrPhotonPCM(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::PhotonPCMs const& photonPCMs)
  {
    size_t const nTriggersThisDataFrame = triggers.size();

    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const photonPCMsThisEvent = photonPCMs.sliceBy(perColPhotonPCMs, collision.globalIndex());

      auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h3_ptPhiEta_photonPCM"),
                    associated.pt(), associated.phi(), associated.eta(),
                    getInvEff<EffParticleType::PhotonPCM>(associated.pt()));
      };
      corrProcessPlain(collision, photonPCMsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.posTrackId() || trigger.jetTrackId() == associated.negTrackId())
          return;

        histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_photonPCM"),
                    associated.pt(), associated.phi(), associated.eta(),
                    getInvEff<EffParticleType::Trigger>(trigger.pt()) * getInvEff<EffParticleType::PhotonPCM>(associated.pt()));
        histos.fill(HIST("reco/corr/h6_corr_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getInvEff<EffParticleType::Trigger>(trigger.pt()) * getInvEff<EffParticleType::PhotonPCM>(associated.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, photonPCMsThisEvent, funcCorrelation);

      auto const funcMixing = [this](auto const& collision,
                                     float const mixingTriggerPt, float const mixingTriggerPhi, float const mixingTriggerEta, auto const& associated, auto const perTriggerWeight) {
        histos.fill(HIST("reco/corr/h6_mix_photonPCM"),
                    getDeltaPhi(mixingTriggerPhi, associated.phi()),
                    mixingTriggerEta - associated.eta(),
                    mixingTriggerPt, associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getInvEff<EffParticleType::Trigger>(mixingTriggerPt) * getInvEff<EffParticleType::PhotonPCM>(associated.pt()));
      };
      corrProcessMixing(collision, photonPCMsThisEvent, funcMixing, nTriggerMixingPhotonPCM, nTriggersThisDataFrame);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPhotonPCM, "process standard correlation for associated photonPCM", false);

  void processCorrPhotonPCMPair(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::PhotonPCMPairs const& photonPCMPairs)
  {
    size_t const nTriggersThisDataFrame = triggers.size();

    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const photonPCMPairsThisEvent = photonPCMPairs.sliceBy(perColPhotonPCMPairs, collision.globalIndex());

      auto const funcPlain = [this](auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h4_ptMggZPvMult_photonPCMPair"), associated.pt(), associated.mgg(), collision.posZ(), collision.nGlobalTracks());
      };
      corrProcessPlain(collision, photonPCMPairsThisEvent, funcPlain);

      auto const funcPlainTrigEv = [this](auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h4_ptMggZPvMult_trigEv_photonPCMPair"), associated.pt(), associated.mgg(), collision.posZ(), collision.nGlobalTracks());
      };
      if (collision.trigEv())
        corrProcessPlain(collision, photonPCMPairsThisEvent, funcPlainTrigEv);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.posTrack1Id() || trigger.jetTrackId() == associated.negTrack1Id() ||
            trigger.jetTrackId() == associated.negTrack2Id() || trigger.jetTrackId() == associated.posTrack2Id())
          return;

        histos.fill(HIST("reco/corr/h5_ptTrigPtAssocMggZPvMult_assoc_photonPCMPair"),
                    trigger.pt(), associated.pt(), associated.mgg(), collision.posZ(), collision.nGlobalTracks(),
                    getInvEff<EffParticleType::Trigger>(trigger.pt()));

        // pi0
        if (checkMassRange<MassRange::pi0PCMPeak>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_pi0PCMPeak"),
                      associated.pt(), associated.phi(), associated.eta(),
                      getInvEff<EffParticleType::Trigger>(trigger.pt()));
          histos.fill(HIST("reco/corr/h6_corr_pi0PCMPeak"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      getInvEff<EffParticleType::Trigger>(trigger.pt()));
          return;
        }
        if (checkMassRange<MassRange::pi0PCMSide>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_pi0PCMSide"),
                      associated.pt(), associated.phi(), associated.eta(),
                      getInvEff<EffParticleType::Trigger>(trigger.pt()));
          histos.fill(HIST("reco/corr/h6_corr_pi0PCMSide"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      getInvEff<EffParticleType::Trigger>(trigger.pt()));
          return;
        }
        // eta
        if (checkMassRange<MassRange::etaPCMPeak>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_etaPCMPeak"),
                      associated.pt(), associated.phi(), associated.eta(),
                      getInvEff<EffParticleType::Trigger>(trigger.pt()));
          histos.fill(HIST("reco/corr/h6_corr_etaPCMPeak"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      getInvEff<EffParticleType::Trigger>(trigger.pt()));
          return;
        }
        if (checkMassRange<MassRange::etaPCMSide>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_etaPCMSide"),
                      associated.pt(), associated.phi(), associated.eta(),
                      getInvEff<EffParticleType::Trigger>(trigger.pt()));
          histos.fill(HIST("reco/corr/h6_corr_etaPCMSide"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      getInvEff<EffParticleType::Trigger>(trigger.pt()));
          return;
        }
      };
      corrProcessCorrelation(collision, triggersThisEvent, photonPCMPairsThisEvent, funcCorrelation);

      auto const funcMixing = [this](auto const& collision,
                                     float const mixingTriggerPt, float const mixingTriggerPhi, float const mixingTriggerEta, auto const& associated, auto const perTriggerWeight) {
        // pi0
        if (checkMassRange<MassRange::pi0PCMPeak>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h6_mix_pi0PCMPeak"),
                      getDeltaPhi(mixingTriggerPhi, associated.phi()),
                      mixingTriggerEta - associated.eta(),
                      mixingTriggerPt, associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      perTriggerWeight * getInvEff<EffParticleType::Trigger>(mixingTriggerPt));
          return;
        }
        if (checkMassRange<MassRange::pi0PCMSide>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h6_mix_pi0PCMSide"),
                      getDeltaPhi(mixingTriggerPhi, associated.phi()),
                      mixingTriggerEta - associated.eta(),
                      mixingTriggerPt, associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      perTriggerWeight * getInvEff<EffParticleType::Trigger>(mixingTriggerPt));
          return;
        }
        // eta
        if (checkMassRange<MassRange::etaPCMPeak>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h6_mix_etaPCMPeak"),
                      getDeltaPhi(mixingTriggerPhi, associated.phi()),
                      mixingTriggerEta - associated.eta(),
                      mixingTriggerPt, associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      perTriggerWeight * getInvEff<EffParticleType::Trigger>(mixingTriggerPt));
          return;
        }
        if (checkMassRange<MassRange::etaPCMSide>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h6_mix_etaPCMSide"),
                      getDeltaPhi(mixingTriggerPhi, associated.phi()),
                      mixingTriggerEta - associated.eta(),
                      mixingTriggerPt, associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      perTriggerWeight * getInvEff<EffParticleType::Trigger>(mixingTriggerPt));
          return;
        }
      };
      corrProcessMixing(collision, photonPCMPairsThisEvent, funcMixing, nTriggerMixingH0PCM, nTriggersThisDataFrame);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPhotonPCMPair, "process standard correlation for associated pi0PCM", false);

  void processCorrPhotonPCMPairMix(CorrCollisions const& collisions, aod::PhotonPCMs const& photonPCMs)
  {
    auto photonPCMsTuple = std::make_tuple(photonPCMs);
    SameKindPair<CorrCollisions, aod::PhotonPCMs, BinningZPvMult> pairs{binningZPvMult, nNeighboursMixingPhotonPCMPair, -1, collisions, photonPCMsTuple, &cache};

    // mixed events
    for (auto pair = pairs.begin(); pair != pairs.end(); pair++) {
      auto const& [collision1, photonPCMs1, collision2, photonPCMs2] = *pair;

      // // check that current und mixing-trigger event are from the same zPv/mult bins
      // if (checkSameBin(collision1.posZ(), collision2.posZ(), binsZPv) == -1) {
      //   std::printf("ERROR: zPv bins do not match\n"); continue;
      // }
      // if (checkSameBin(collision1.nGlobalTracks(), collision2.nGlobalTracks(), binsMult) == -1) {
      //   std::printf("ERROR: multiplicity bins do not match\n"); continue;
      // }

      // event selection
      if (!collision1.selEv() || !collision2.selEv())
        continue;
      // event info
      histos.fill(HIST("reco/plain/h2_zPvMult_photonPCMPair_evMix"), collision1.posZ(), collision1.nGlobalTracks());
      // mixing loop
      for (auto const& [photonPCM1, photonPCM2] : soa::combinations(soa::CombinationsFullIndexPolicy(photonPCMs1, photonPCMs2))) {
        ROOT::Math::PtEtaPhiMVector const p4photonPCM1(photonPCM1.pt(), photonPCM1.eta(), photonPCM1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPCM2(photonPCM2.pt(), photonPCM2.eta(), photonPCM2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPCMPair = p4photonPCM1 + p4photonPCM2;

        histos.fill(HIST("reco/plain/h4_ptMggZPvMult_photonPCMPair_evMix"), p4photonPCMPair.pt(), p4photonPCMPair.M(), collision1.posZ(), collision1.nGlobalTracks());
      }

      // trigger events
      if (!collision1.trigEv() || !collision2.trigEv())
        continue;
      // mixing loop
      for (auto const& [photonPCM1, photonPCM2] : soa::combinations(soa::CombinationsFullIndexPolicy(photonPCMs1, photonPCMs2))) {
        ROOT::Math::PtEtaPhiMVector const p4photonPCM1(photonPCM1.pt(), photonPCM1.eta(), photonPCM1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPCM2(photonPCM2.pt(), photonPCM2.eta(), photonPCM2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPCMPair = p4photonPCM1 + p4photonPCM2;

        histos.fill(HIST("reco/plain/h4_ptMggZPvMult_trigEv_photonPCMPair_evMix"), p4photonPCMPair.pt(), p4photonPCMPair.M(), collision1.posZ(), collision1.nGlobalTracks());
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPhotonPCMPairMix, "process gamma-gamma mixing for photonPCM", false);

  // mc ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void processMcInfo(CorrMcCollisions const& mcCollisions, CorrMcDCollisions const& collisions)
  {
    for (auto const& mcCollision : mcCollisions) {
      // all events
      histos.fill(HIST("mc/info/h1_nEvents_mcTrue"), 0.5);
      histos.fill(HIST("mc/info/h2_zPvMult_mcTrue"), mcCollision.posZ(), mcCollision.nChargedInEtaRange());

      // trigger events
      if (!mcCollision.trigEv())
        continue;
      histos.fill(HIST("mc/info/h1_nTrigEv_mcTrue"), 0.5);
      histos.fill(HIST("mc/info/h2_zPvMult_trigEv_mcTrue"), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
    }
    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;
      histos.fill(HIST("mc/info/h1_nRecoCol_mcTrue"), 0.5);
      histos.fill(HIST("mc/info/h2_zPvMult_recoCol_mcTrue"), collision.mcCollision_as<CorrMcCollisions>().posZ(), collision.mcCollision_as<CorrMcCollisions>().nChargedInEtaRange());
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcInfo, "process general info on mc collisions and tracks for analysis and qa", false);

  // (sad) attempt at reducing code duplication
  enum class McCorrEventType : int { True = 0,
                                     RecoColTrue = 1 };
  enum class McCorrCorrelationType : int { Correlation = 0,
                                           Mixing = 1 };
  enum class McCorrAssociatedType : int { Hadron = 0,
                                          Pipm = 1,
                                          Photon = 2,
                                          Pi0 = 3,
                                          Eta = 4 };
  static constexpr const char* McHistPaths[2][2][5] = {
    {{"mc/true/corr/h6_corr_hadron", "mc/true/corr/h6_corr_pipm", "mc/true/corr/h6_corr_photon",
      "mc/true/corr/h6_corr_pi0", "mc/true/corr/h6_corr_eta"},
     {"mc/true/corr/h6_mix_hadron", "mc/true/corr/h6_mix_pipm", "mc/true/corr/h6_mix_photon",
      "mc/true/corr/h6_mix_pi0", "mc/true/corr/h6_mix_eta"}},
    {{"mc/recoCol_true/corr/h6_corr_hadron", "mc/recoCol_true/corr/h6_corr_pipm", "mc/recoCol_true/corr/h6_corr_photon",
      "mc/recoCol_true/corr/h6_corr_pi0", "mc/recoCol_true/corr/h6_corr_eta"},
     {"mc/recoCol_true/corr/h6_mix_hadron", "mc/recoCol_true/corr/h6_mix_pipm", "mc/recoCol_true/corr/h6_mix_photon",
      "mc/recoCol_true/corr/h6_mix_pi0", "mc/recoCol_true/corr/h6_mix_eta"}}};
  static constexpr const char* getMcHistPath(McCorrEventType eventType, McCorrCorrelationType correlationType, McCorrAssociatedType associatedType)
  {
    return McHistPaths[static_cast<int>(eventType)][static_cast<int>(correlationType)][static_cast<int>(associatedType)];
  }

  // fill mc correaltion histograms based on given associated mc particle
  template <McCorrEventType eventType, McCorrCorrelationType correlationType>
  void fillMcCorrHists(auto const& mcCollision, auto const& trigger, auto const& associated, double const weight)
  {
    // standard particles (marked physical primary)
    if (checkPrimaryEtaMc(associated)) {
      // charged primary ('hadron') selection
      if (checkChargedMc(associated)) {
        histos.fill(HIST(getMcHistPath(eventType, correlationType, McCorrAssociatedType::Hadron)),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), mcCollision.posZ(), mcCollision.nChargedInEtaRange(),
                    weight);
      }
      // pipm selection
      if (std::abs(associated.pdgCode()) == PDG_t::kPiPlus) {
        histos.fill(HIST(getMcHistPath(eventType, correlationType, McCorrAssociatedType::Pipm)),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), mcCollision.posZ(), mcCollision.nChargedInEtaRange(),
                    weight);
        return;
      }
      // photon selection
      if (associated.pdgCode() == PDG_t::kGamma) {
        histos.fill(HIST(getMcHistPath(eventType, correlationType, McCorrAssociatedType::Photon)),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), mcCollision.posZ(), mcCollision.nChargedInEtaRange(),
                    weight);
        return;
      }
      return;
    }
    // decaying particles (not marked physical primary)
    if ((std::abs(associated.eta()) < etaMax)) {
      // pi0 selection
      if (checkH0Primary(associated, PDG_t::kPi0)) {
        histos.fill(HIST(getMcHistPath(eventType, correlationType, McCorrAssociatedType::Pi0)),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), mcCollision.posZ(), mcCollision.nChargedInEtaRange(),
                    weight);
        return;
      }
      // eta selection
      if (checkH0Primary(associated, 221)) {
        histos.fill(HIST(getMcHistPath(eventType, correlationType, McCorrAssociatedType::Eta)),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), mcCollision.posZ(), mcCollision.nChargedInEtaRange(),
                    weight);
        return;
      }
    }
  }

  void processMcTrueCorr(CorrMcCollisions const& mcCollisions, aod::TriggerParticles const& triggerParticles, aod::JetParticles const& mcParticles)
  {
    for (auto const& mcCollision : mcCollisions) {
      // group collision
      auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, mcCollision.globalIndex());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, mcCollision.globalIndex());

      // trigger pairing loop
      for (auto const& trigger : triggerParticlesThisEvent) {
        // trigger info
        histos.fill(HIST("mc/true/corr/h3_ptPhiEta_trig"), trigger.pt(), trigger.phi(), trigger.eta());

        // save trigger for mixing
        mixingTriggerMemoryTrue.saveTrigger(trigger.pt(), trigger.phi(), trigger.eta(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());

        for (auto const& associated : mcParticlesThisEvent) {
          // exclude self correlation
          if (trigger.jetMcParticleId() == associated.globalIndex())
            continue;

          fillMcCorrHists<McCorrEventType::True, McCorrCorrelationType::Correlation>(mcCollision, trigger, associated, 1);
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueCorr, "process mc-true (all collisions) correlation for multiple associated particles", false);

  void processMcTrueMix(CorrMcCollisions const& mcCollisions, aod::TriggerParticles const& triggerParticles, aod::JetParticles const& mcParticles)
  {
    for (auto const& mcCollision : mcCollisions) {
      // group collision
      auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, mcCollision.globalIndex());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, mcCollision.globalIndex());

      const size_t nTriggerParticlesThisDataFrame = triggerParticles.size();
      auto savedTriggers = mixingTriggerMemoryTrue.getTriggers(mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      const size_t mixUpToTriggerN = std::min(savedTriggers.size(), static_cast<size_t>(nTriggerMixingMcTrue) + nTriggerParticlesThisDataFrame);
      const float perTriggerWeight = 1. / (mixUpToTriggerN - nTriggerParticlesThisDataFrame);

      // trigger loop
      for (size_t i_mixingTrigger = nTriggerParticlesThisDataFrame; i_mixingTrigger < mixUpToTriggerN; i_mixingTrigger++) {
        MixingTrigger const& mixingTrigger = savedTriggers[i_mixingTrigger];
        for (auto const& associated : mcParticlesThisEvent) {
          fillMcCorrHists<McCorrEventType::True, McCorrCorrelationType::Mixing>(mcCollision, mixingTrigger, associated, perTriggerWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueMix, "process mc-true (all collisions) correlation mixing for multiple associated particles", false);

  void processMcRecoColTrueCorr(CorrMcDCollisions const& collisions, CorrMcCollisions const&, aod::TriggerParticles const& triggerParticles, aod::JetParticles const& mcParticles)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;

      // group collision
      auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, collision.mcCollisionId());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.mcCollisionId());

      auto const& mcCollision = collision.mcCollision_as<CorrMcCollisions>();

      // trigger pairing loop
      for (auto const& trigger : triggerParticlesThisEvent) {
        // trigger info
        histos.fill(HIST("mc/recoCol_true/corr/h3_ptPhiEta_trig"), trigger.pt(), trigger.phi(), trigger.eta());

        // save trigger for mixing
        mixingTriggerMemoryRecoColTrue.saveTrigger(trigger.pt(), trigger.phi(), trigger.eta(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());

        // hadrons (tracks) and pipm
        for (auto const& associated : mcParticlesThisEvent) {
          // exclude self correlation
          if (trigger.jetMcParticleId() == associated.globalIndex())
            continue;

          fillMcCorrHists<McCorrEventType::RecoColTrue, McCorrCorrelationType::Correlation>(mcCollision, trigger, associated, 1);
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoColTrueCorr, "process mc-true (reco collisions) correlation for multiple associated particles", false);

  void processMcRecoColTrueMix(CorrMcDCollisions const& collisions, CorrMcCollisions const&, aod::TriggerParticles const& triggerParticles, aod::JetParticles const& mcParticles)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;

      // group collision
      auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, collision.mcCollisionId());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.mcCollisionId());

      auto const& mcCollision = collision.mcCollision_as<CorrMcCollisions>();

      const size_t nTriggerParticlesThisDataFrame = triggerParticles.size();
      auto savedTriggers = mixingTriggerMemoryRecoColTrue.getTriggers(mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      const size_t mixUpToTriggerN = std::min(savedTriggers.size(), static_cast<size_t>(nTriggerMixingMcTrue) + nTriggerParticlesThisDataFrame);
      const float perTriggerWeight = 1. / (mixUpToTriggerN - nTriggerParticlesThisDataFrame);

      // trigger loop
      for (size_t i_mixingTrigger = nTriggerParticlesThisDataFrame; i_mixingTrigger < mixUpToTriggerN; i_mixingTrigger++) {
        MixingTrigger const& mixingTrigger = savedTriggers[i_mixingTrigger];
        for (auto const& associated : mcParticlesThisEvent) {
          fillMcCorrHists<McCorrEventType::RecoColTrue, McCorrCorrelationType::Mixing>(mcCollision, mixingTrigger, associated, perTriggerWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoColTrueMix, "process mc-true (reco collisions) correlation mixing for multiple associated particles", false);

  void processMcRecoColEff(CorrMcDCollision const& collision, aod::JetTracksMCD const& tracks,
                           aod::Hadrons const& hadrons, aod::Pipms const& pipms, aod::PhotonPCMs const& photonPCMs,
                           CorrMcCollisions const&, aod::JetParticles const& mcParticles, aod::TriggerParticles const& triggerParticles)
  {
    // event selection
    if (!collision.selEv())
      return;

    auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.mcCollisionId());
    auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, collision.mcCollisionId());

    // hadrons
    for (auto const& hadron : hadrons) {
      if (doTrigEvEff && !collision.trigEv() && hadron.pt() < ptCutTrigEvEff)
        continue;
      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hadron"), hadron.pt(), hadron.phi(), hadron.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hadron"), hadron.pt(), collision.posZ(), collision.nGlobalTracks());
      // purity
      if (!hadron.jetTrack_as<aod::JetTracksMCD>().has_mcParticle())
        continue;
      auto const hadronParticle = hadron.jetTrack_as<aod::JetTracksMCD>().mcParticle();
      if (!checkPrimaryTrackMc(hadronParticle))
        continue;
      if (requireSingleCollisionPurity && hadronParticle.mcCollisionId() != collision.mcCollisionId())
        continue;

      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hasCorrectMc_hadron"), hadron.pt(), hadron.phi(), hadron.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hasCorrectMc_hadron"), hadron.pt(), collision.posZ(), collision.nGlobalTracks());
    }

    // pipm
    for (auto const& pipm : pipms) {
      if (doTrigEvEff && !collision.trigEv() && pipm.pt() < ptCutTrigEvEff)
        continue;
      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_pipm"), pipm.pt(), pipm.phi(), pipm.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_pipm"), pipm.pt(), collision.posZ(), collision.nGlobalTracks());
      // purity
      if (!pipm.jetTrack_as<aod::JetTracksMCD>().has_mcParticle())
        continue;
      auto const pipmParticle = pipm.jetTrack_as<aod::JetTracksMCD>().mcParticle();
      if (std::abs(pipmParticle.pdgCode()) != PDG_t::kPiPlus || !checkPrimaryEtaMc(pipmParticle))
        continue;
      if (requireSingleCollisionPurity && pipmParticle.mcCollisionId() != collision.mcCollisionId())
        continue;

      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hasCorrectMc_pipm"), pipm.pt(), pipm.phi(), pipm.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hasCorrectMc_pipm"), pipm.pt(), collision.posZ(), collision.nGlobalTracks());
    }

    // photonPCM
    for (auto const& photonPCM : photonPCMs) {
      if (doTrigEvEff && !collision.trigEv())
        continue;
      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_photonPCM"), photonPCM.pt(), photonPCM.phi(), photonPCM.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_photonPCM"), photonPCM.pt(), collision.posZ(), collision.nGlobalTracks());

      // purity
      // (V0Legs does not have the tracks reference as index column (just int)??)
      auto const& posTrack = tracks.rawIteratorAt(photonPCM.posTrackId() - tracks.offset());
      auto const& negTrack = tracks.rawIteratorAt(photonPCM.negTrackId() - tracks.offset());
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle())
        continue;
      if (!isConversionPhoton(posTrack, negTrack) || !checkPrimaryEtaMc(*(posTrack.mcParticle().mothers_as<aod::JetParticles>().begin())))
        continue;
      if (requireSingleCollisionPurity && posTrack.mcParticle().mcCollisionId() != collision.mcCollisionId())
        continue;

      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hasCorrectMc_photonPCM"), photonPCM.pt(), photonPCM.phi(), photonPCM.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hasCorrectMc_photonPCM"), photonPCM.pt(), collision.posZ(), collision.nGlobalTracks());
    }

    // mcParticle loop
    for (auto const& mcParticle : mcParticlesThisEvent) {
      bool const countChargedTrigEvEff = !doTrigEvEff || collision.trigEv() || mcParticle.pt() > ptCutTrigEvEff;
      bool const countOtherTrigEvEff = !doTrigEvEff || collision.trigEv();

      // standard particles (marked physical primary)
      if (checkPrimaryEtaMc(mcParticle)) {
        // hadrons
        if (checkChargedMc(mcParticle) && countChargedTrigEvEff) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_hadron"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_hadron"), mcParticle.pt(), collision.mcCollision_as<CorrMcCollisions>().posZ(), collision.nGlobalTracks());
        }
        // pipm
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && countChargedTrigEvEff) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_pipm"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_pipm"), mcParticle.pt(), collision.mcCollision_as<CorrMcCollisions>().posZ(), collision.nGlobalTracks());
        }
        // photons
        if (mcParticle.pdgCode() == PDG_t::kGamma && countOtherTrigEvEff) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_photon"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_photon"), mcParticle.pt(), collision.mcCollision_as<CorrMcCollisions>().posZ(), collision.nGlobalTracks());
        }
      }

      // decaying particles (not marked physical primary)
      if ((std::abs(mcParticle.eta()) < etaMax)) {
        // pi0
        if (checkH0ToGG(mcParticle, PDG_t::kPi0) && countOtherTrigEvEff) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_pi0"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_pi0"), mcParticle.pt(), collision.mcCollision_as<CorrMcCollisions>().posZ(), collision.nGlobalTracks());
        }
        // eta
        if (checkH0ToGG(mcParticle, 221) && countOtherTrigEvEff) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_eta"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_eta"), mcParticle.pt(), collision.mcCollision_as<CorrMcCollisions>().posZ(), collision.nGlobalTracks());
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoColEff, "process MC data to calculate efficiencies and purities", false);

  // test /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void processTest(CorrCollision const& collision,
                   soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::TracksExtra> const&,
                   aod::Hadrons const& hadrons)
  {
    // event selection
    if (!collision.selEv())
      return;

    histos.fill(HIST("test/h2_mult_comp"), collision.nGlobalTracks(), hadrons.size());

    for (auto const& track : tracks) {
      auto const fullTrack = track.track_as<soa::Join<aod::Tracks, aod::TracksExtra>>();

      constexpr float Mincrossedrows = 40;
      constexpr float Maxchi2tpc = 5.0;
      constexpr float Maxchi2its = 6.0;
      constexpr float MaxR = 83.1;
      constexpr float MinPtTrackiu = 0.1;

      if (!fullTrack.hasITS() && !fullTrack.hasTPC())
        continue;
      if (fullTrack.x() * fullTrack.x() + fullTrack.y() * fullTrack.y() > MaxR * MaxR || fullTrack.pt() < MinPtTrackiu)
        continue;
      if (fullTrack.hasTPC()) {
        if (fullTrack.tpcNClsCrossedRows() < Mincrossedrows || fullTrack.tpcChi2NCl() > Maxchi2tpc)
          continue;
      }
      if (fullTrack.hasITS()) {
        if (fullTrack.itsChi2NCl() > Maxchi2its)
          continue;
      }

      histos.fill(HIST("test/h2_tracks_zPvMultDep"), collision.posZ(), collision.nGlobalTracks());
    }

    histos.fill(HIST("test/h2_globalTracks_zPvMultDep"), collision.posZ(), collision.nGlobalTracks(), hadrons.size());
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processTest, "process just to test things", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& configContext)
{
  return WorkflowSpec{adaptAnalysisTask<PhotonChargedTriggerCorrelation>(configContext)};
}
