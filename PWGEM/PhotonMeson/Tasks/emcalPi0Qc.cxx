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

/// \file emcalPi0Qc.cxx
/// \brief Simple monitoring task for EMCal clusters
/// \author Joshua Koenig <joshua.konig@cern.ch>, Goethe University Frankfurt
/// \struct EmcalPi0Qc
/// \since 25.05.2022
///
/// This task is meant to be used for QC for the emcal using properties of the pi0
/// - ...
/// Simple event selection using the flag doEventSel is provided, which selects INT7 events if set to 1
/// For pilot beam data, instead of relying on the event selection, one can veto specific BC IDS using the flag
/// fDoVetoBCID.

#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonDataFormat/InteractionRecord.h>
#include <EMCALBase/Geometry.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TVector3.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

using namespace o2::framework;
using namespace o2::framework::expressions;
using MyCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::EMCALMatchedCollisions>;
using MyBCs = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels, o2::aod::Timestamps>;

struct Photon {
  Photon(float eta_tmp, float phi_tmp, float energy_tmp, int clusid = 0, uint8_t sm_tmp = 0)
  {
    eta = eta_tmp;
    phi = phi_tmp;
    energy = energy_tmp;
    theta = 2 * std::atan2(std::exp(-eta), 1);
    px = energy * std::sin(theta) * std::cos(phi);
    py = energy * std::sin(theta) * std::sin(phi);
    pz = energy * std::cos(theta);
    pt = std::sqrt(px * px + py * py);
    photon.SetPxPyPzE(px, py, pz, energy);
    id = clusid;
    sm = sm_tmp;
    onDCal = (phi < 6 && phi > 4);
  }

  TLorentzVector photon;
  float pt;
  float px;
  float py;
  float pz;
  float eta;
  float phi;
  float energy;
  float theta;
  int id;
  uint8_t sm;
  bool onDCal; // Checks whether photon is in phi region of the DCal, otherwise: EMCal
};

struct Meson {
  Meson(Photon p1, Photon p2) : pgamma1(p1),
                                pgamma2(p2)
  {
    pMeson = p1.photon + p2.photon;
  }
  Photon pgamma1;
  Photon pgamma2;
  TLorentzVector pMeson;

  float getMass() const { return pMeson.M(); }
  float getPt() const { return pMeson.Pt(); }
  float getOpeningAngle() const { return pgamma1.photon.Angle(pgamma2.photon.Vect()); }
};

struct EventMixVec {

  void addEvent(std::vector<Photon> vecGamma)
  {
    if (vecEvtMix.size() < nEVtMixSize) {
      vecEvtMix.push_back(vecGamma);
    } else {
      vecEvtMix.erase(vecEvtMix.begin() + nEVtMixSize - 1);
      vecEvtMix.push_back(vecGamma);
    }
  }

  std::vector<std::vector<Photon>> vecEvtMix;
  unsigned int nEVtMixSize = 20;
};

struct EmcalPi0Qc {
  HistogramRegistry mHistManager{"NeutralMesonHistograms"};
  o2::emcal::Geometry* mGeometry = nullptr;

  Filter emccellfilter = o2::aod::calo::caloType == 1;

  // configurable parameters
  // TODO adapt mDoEventSel switch to also allow selection of other triggers (e.g. EMC7)
  Configurable<bool> mDoEventSel{"mDoEventSel", 0, "demand kINT7"};
  Configurable<bool> mRequireCaloReadout{"mRequireCaloReadout", 0, "require kTVXinEMC"};
  Configurable<bool> mRequireEMCalCells{"mRequireEMCalCells", 0, "require at least one EMC cell in each collision"};
  Configurable<std::string> mVetoBCID{"mVetoBCID", "", "BC ID(s) to be excluded, this should be used as an alternative to the event selection"};
  Configurable<std::string> mSelectBCID{"mSelectBCID", "all", "BC ID(s) to be included, this should be used as an alternative to the event selection"};
  Configurable<double> mVertexCut{"mVertexCut", -1, "apply z-vertex cut with value in cm"};
  Configurable<int> mTimeMin{"mTimeMin", -600, "apply min timing cut (in ns)"};
  Configurable<int> mTimeMax{"mTimeMax", 900, "apply min timing cut (in ns)"};
  Configurable<float> mClusterMinM02Cut{"mClusterMinM02Cut", 0.1, "apply min M02 cut"};
  Configurable<float> mClusterMaxM02Cut{"mClusterMaxM02Cut", 0.7, "apply max M02 cut"};
  Configurable<float> mMinEnergyCut{"mMinEnergyCut", 0.7, "apply min cluster energy cut"};
  Configurable<int> mMinNCellsCut{"mMinNCellsCut", 1, "apply min cluster number of cell cut"};
  Configurable<float> mMinOpenAngleCut{"mMinOpenAngleCut", 0.0202, "apply min opening angle cut"};
  Configurable<std::string> mClusterDefinition{"mClusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<bool> mSplitEMCalDCal{"mSplitEMCalDCal", 0, "Create and fill inv mass histograms for photons on EMCal and DCal individually"};
  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;

  ConfigurableAxis pTBinning{"pTBinning", {500, 0.0f, 50.0f}, "Binning used along pT axis for inv mass histograms"};
  ConfigurableAxis invmassBinning{"invmassBinning", {400, 0.0f, 0.8f}, "Binning used for inv mass axis in inv mass - pT histograms"};

  // define cluster filter. It selects only those clusters which are of the type
  // specified in the string mClusterDefinition,e.g. kV3Default, which is V3 clusterizer with default
  // clusterization parameters
  // o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  // Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);

  // define container for photons
  std::vector<Photon> mPhotons;

  // event mixing class
  EventMixVec evtMix;

  o2::ccdb::CcdbApi ccdbApi;
  int lastRunNumber = -1; // get the runnumber to obtain the SOR of the run to get t - SOR in (s) later
  int64_t tsSOR = -1;
  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::kV3Default;

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    // init ccdb api
    ccdbApi.init("https://alice-ccdb.cern.ch");

    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // create common axes
    LOG(info) << "Creating histograms";
    const AxisSpec bcAxis{3501, -0.5, 3500.5};
    const AxisSpec energyAxis{makeClusterBinning(), "#it{E} (GeV)"};
    const AxisSpec collisionTimeAxis{1440, 0, 1440, "#it{t} - SOR (min)"};
    const AxisSpec clusterTimeAxis{1500, -600, 900, "t_{cl} (ns)"};
    const AxisSpec invmassAxis{invmassBinning, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec ptAxis{pTBinning, "#it{p}_{T} (GeV/#it{c})"};

    if (doprocessCollision) {
      mHistManager.add("events", "events;;#it{count}", HistType::kTH1F, {{7, 0.5, 7.5}});
      auto heventType = mHistManager.get<TH1>(HIST("events"));
      heventType->GetXaxis()->SetBinLabel(1, "All events");
      heventType->GetXaxis()->SetBinLabel(2, "sel8");
      heventType->GetXaxis()->SetBinLabel(3, "EMCal readout");
      heventType->GetXaxis()->SetBinLabel(4, "1+ Contributor");
      heventType->GetXaxis()->SetBinLabel(5, "z<10cm");
      heventType->GetXaxis()->SetBinLabel(6, "unique col");
      heventType->GetXaxis()->SetBinLabel(7, "EMCal cell>0");
      mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", HistType::kTH1F, {{200, -20, 20}});
      mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", HistType::kTH1F, {{200, -20, 20}});
      mHistManager.add("hEventPerTime", "number of events per time", HistType::kTH1F, {collisionTimeAxis});
    }

    if (doprocessAmbiguous) {
      mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", HistType::kTH1F, {bcAxis});
      mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", HistType::kTH1F, {bcAxis});
    }
    // cluster properties
    for (const bool& iBeforeCuts : {false, true}) {
      const char* clusterDirectory = iBeforeCuts ? "ClustersBeforeCuts" : "ClustersAfterCuts";
      mHistManager.add(Form("%s/clusterE", clusterDirectory), "Energy of cluster", HistType::kTH1F, {energyAxis});
      mHistManager.add(Form("%s/clusterE_SimpleBinning", clusterDirectory), "Energy of cluster", HistType::kTH1F, {{400, 0, 100, "#it{E} (GeV)"}});
      mHistManager.add(Form("%s/clusterTime", clusterDirectory), "Time of cluster", HistType::kTH1F, {{500, -250, 250, "#it{t}_{cls} (ns)"}});
      mHistManager.add(Form("%s/clusterEtaPhi", clusterDirectory), "Eta and phi of cluster", HistType::kTH2F, {{100, -1, 1, "#eta"}, {100, 0, o2::constants::math::TwoPI, "#phi"}});
      mHistManager.add(Form("%s/clusterM02", clusterDirectory), "M02 of cluster", HistType::kTH1F, {{400, 0, 5, "#it{M}_{02}"}});
      mHistManager.add(Form("%s/clusterM20", clusterDirectory), "M20 of cluster", HistType::kTH1F, {{400, 0, 2.5, "#it{M}_{20}"}});
      mHistManager.add(Form("%s/clusterNLM", clusterDirectory), "Number of local maxima of cluster", HistType::kTH1I, {{10, 0, 10, "#it{N}_{local maxima}"}});
      mHistManager.add(Form("%s/clusterNCells", clusterDirectory), "Number of cells in cluster", HistType::kTH1I, {{50, 0, 50, "#it{N}_{cells}"}});
      mHistManager.add(Form("%s/clusterDistanceToBadChannel", clusterDirectory), "Distance to bad channel", HistType::kTH1F, {{100, 0, 100, "#it{d}"}});
    }

    // meson related histograms
    mHistManager.add("invMassVsPt", "invariant mass and pT of meson candidates", HistType::kTH2F, {invmassAxis, ptAxis});
    mHistManager.add("invMassVsPtBackground", "invariant mass and pT of background meson candidates", HistType::kTH2F, {invmassAxis, ptAxis});
    mHistManager.add("invMassVsPtMixedBackground", "invariant mass and pT of mixed background meson candidates", HistType::kTH2F, {invmassAxis, ptAxis});

    if (mSplitEMCalDCal) {
      mHistManager.add("invMassVsPt_EMCal", "invariant mass and pT of meson candidates with both clusters on EMCal", HistType::kTH2F, {invmassAxis, ptAxis});
      mHistManager.add("invMassVsPtBackground_EMCal", "invariant mass and pT of background meson candidates with both clusters on EMCal", HistType::kTH2F, {invmassAxis, ptAxis});
      mHistManager.add("invMassVsPtMixedBackground_EMCal", "invariant mass and pT of mixed background meson candidates with both clusters on EMCal", HistType::kTH2F, {invmassAxis, ptAxis});
      mHistManager.add("invMassVsPt_DCal", "invariant mass and pT of meson candidates with both clusters on DCal", HistType::kTH2F, {invmassAxis, ptAxis});
      mHistManager.add("invMassVsPtBackground_DCal", "invariant mass and pT of background meson candidates with both clusters on DCal", HistType::kTH2F, {invmassAxis, ptAxis});
      mHistManager.add("invMassVsPtMixedBackground_DCal", "invariant mass and pT of mixed background meson candidates with both clusters on DCal", HistType::kTH2F, {invmassAxis, ptAxis});
    }

    // add histograms per supermodule
    for (int ism = 0; ism < 20; ++ism) {
      mHistManager.add(Form("clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM%d", ism), Form("Cluster time vs collision timestamp in Supermodule %d", ism), HistType::kTH2F, {clusterTimeAxis, collisionTimeAxis});
      mHistManager.add(Form("clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM%d", ism), Form("Cluster number of cells vs collision timestamp in Supermodule %d", ism), HistType::kTH2F, {{50, 0, 50}, collisionTimeAxis});
      mHistManager.add(Form("clusterM02VsTimeStamp/clusterM02VsTimeStampSM%d", ism), Form("Cluster M02 vs collision timestamp in Supermodule %d", ism), HistType::kTH2F, {{400, 0, 5}, collisionTimeAxis});
      mHistManager.add(Form("mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM%d", ism), Form("invariant mass vs collision timestamp in Supermodule %d", ism), HistType::kTH2F, {invmassAxis, collisionTimeAxis});
    }

    if (mVetoBCID->length()) {
      std::stringstream parser(mVetoBCID.value);
      std::string token;
      int bcid;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Veto BCID " << bcid;
        mVetoBCIDs.push_back(bcid);
      }
    }
    if (mSelectBCID.value != "all") {
      std::stringstream parser(mSelectBCID.value);
      std::string token;
      int bcid;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Select BCID " << bcid;
        mSelectBCIDs.push_back(bcid);
      }
    }
    clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
    LOG(info) << "mDoEventSel = " << mDoEventSel.value;
    LOG(info) << "mRequireCaloReadout = " << mRequireCaloReadout.value;
    LOG(info) << "mRequireEMCalCells = " << mRequireEMCalCells.value;
    LOG(info) << "mSplitEMCalDCal = " << mSplitEMCalDCal.value;
  }

  template <uint8_t supermoduleID>
  void supermoduleHistHelperPhoton(float time, float m02, int NCell, float timeSinceSOR)
  {
    static constexpr std::string_view ClusterTimeHistSM[20] = {"clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM0", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM1", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM2", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM3", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM4", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM5", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM6", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM7", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM8", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM9", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM10", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM11", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM12", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM13", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM14", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM15", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM16", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM17", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM18", "clusterTimeVsTimeStamp/clusterTimeVsTimeStampSM19"};
    static constexpr std::string_view ClusterNcellHistSM[20] = {"clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM0", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM1", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM2", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM3", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM4", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM5", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM6", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM7", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM8", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM9", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM10", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM11", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM12", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM13", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM14", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM15", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM16", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM17", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM18", "clusterNcellVsTimeStamp/clusterNCellVsTimeStampSM19"};
    static constexpr std::string_view ClusterM02HistSM[20] = {"clusterM02VsTimeStamp/clusterM02VsTimeStampSM0", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM1", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM2", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM3", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM4", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM5", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM6", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM7", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM8", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM9", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM10", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM11", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM12", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM13", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM14", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM15", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM16", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM17", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM18", "clusterM02VsTimeStamp/clusterM02VsTimeStampSM19"};
    mHistManager.fill(HIST(ClusterTimeHistSM[supermoduleID]), time, timeSinceSOR);
    mHistManager.fill(HIST(ClusterNcellHistSM[supermoduleID]), NCell, timeSinceSOR);
    mHistManager.fill(HIST(ClusterM02HistSM[supermoduleID]), m02, timeSinceSOR);
  }

  template <uint8_t supermoduleID>
  void supermoduleHistHelperMeson(float minv, float timeSinceSOR)
  {
    static constexpr std::string_view MesonInvMassHistSM[20] = {"mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM0", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM1", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM2", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM3", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM4", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM5", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM6", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM7", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM8", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM9", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM10", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM11", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM12", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM13", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM14", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM15", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM16", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM17", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM18", "mesonInvMassVsTimeStamp/mesonInvMassVsTimeStampSM19"};
    mHistManager.fill(HIST(MesonInvMassHistSM[supermoduleID]), minv, timeSinceSOR);
  }

  void fillSupermoduleHistogramsPhoton(int supermoduleID, float time, float m02, int NCell, float timeSinceSOR)
  {
    // workaround to have the histogram names per supermodule
    // handled as constexpr
    switch (supermoduleID) {
      case 0:
        supermoduleHistHelperPhoton<0>(time, m02, NCell, timeSinceSOR);
        break;
      case 1:
        supermoduleHistHelperPhoton<1>(time, m02, NCell, timeSinceSOR);
        break;
      case 2:
        supermoduleHistHelperPhoton<2>(time, m02, NCell, timeSinceSOR);
        break;
      case 3:
        supermoduleHistHelperPhoton<3>(time, m02, NCell, timeSinceSOR);
        break;
      case 4:
        supermoduleHistHelperPhoton<4>(time, m02, NCell, timeSinceSOR);
        break;
      case 5:
        supermoduleHistHelperPhoton<5>(time, m02, NCell, timeSinceSOR);
        break;
      case 6:
        supermoduleHistHelperPhoton<6>(time, m02, NCell, timeSinceSOR);
        break;
      case 7:
        supermoduleHistHelperPhoton<7>(time, m02, NCell, timeSinceSOR);
        break;
      case 8:
        supermoduleHistHelperPhoton<8>(time, m02, NCell, timeSinceSOR);
        break;
      case 9:
        supermoduleHistHelperPhoton<9>(time, m02, NCell, timeSinceSOR);
        break;
      case 10:
        supermoduleHistHelperPhoton<10>(time, m02, NCell, timeSinceSOR);
        break;
      case 11:
        supermoduleHistHelperPhoton<11>(time, m02, NCell, timeSinceSOR);
        break;
      case 12:
        supermoduleHistHelperPhoton<12>(time, m02, NCell, timeSinceSOR);
        break;
      case 13:
        supermoduleHistHelperPhoton<13>(time, m02, NCell, timeSinceSOR);
        break;
      case 14:
        supermoduleHistHelperPhoton<14>(time, m02, NCell, timeSinceSOR);
        break;
      case 15:
        supermoduleHistHelperPhoton<15>(time, m02, NCell, timeSinceSOR);
        break;
      case 16:
        supermoduleHistHelperPhoton<16>(time, m02, NCell, timeSinceSOR);
        break;
      case 17:
        supermoduleHistHelperPhoton<17>(time, m02, NCell, timeSinceSOR);
        break;
      case 18:
        supermoduleHistHelperPhoton<18>(time, m02, NCell, timeSinceSOR);
        break;
      case 19:
        supermoduleHistHelperPhoton<19>(time, m02, NCell, timeSinceSOR);
        break;
      default:
        break;
    }
  }

  void fillSupermoduleHistogramsMeson(int supermoduleID, float minv, float timeSinceSOR)
  {
    // workaround to have the histogram names per supermodule
    // handled as constexpr
    switch (supermoduleID) {
      case 0:
        supermoduleHistHelperMeson<0>(minv, timeSinceSOR);
        break;
      case 1:
        supermoduleHistHelperMeson<1>(minv, timeSinceSOR);
        break;
      case 2:
        supermoduleHistHelperMeson<2>(minv, timeSinceSOR);
        break;
      case 3:
        supermoduleHistHelperMeson<3>(minv, timeSinceSOR);
        break;
      case 4:
        supermoduleHistHelperMeson<4>(minv, timeSinceSOR);
        break;
      case 5:
        supermoduleHistHelperMeson<5>(minv, timeSinceSOR);
        break;
      case 6:
        supermoduleHistHelperMeson<6>(minv, timeSinceSOR);
        break;
      case 7:
        supermoduleHistHelperMeson<7>(minv, timeSinceSOR);
        break;
      case 8:
        supermoduleHistHelperMeson<8>(minv, timeSinceSOR);
        break;
      case 9:
        supermoduleHistHelperMeson<9>(minv, timeSinceSOR);
        break;
      case 10:
        supermoduleHistHelperMeson<10>(minv, timeSinceSOR);
        break;
      case 11:
        supermoduleHistHelperMeson<11>(minv, timeSinceSOR);
        break;
      case 12:
        supermoduleHistHelperMeson<12>(minv, timeSinceSOR);
        break;
      case 13:
        supermoduleHistHelperMeson<13>(minv, timeSinceSOR);
        break;
      case 14:
        supermoduleHistHelperMeson<14>(minv, timeSinceSOR);
        break;
      case 15:
        supermoduleHistHelperMeson<15>(minv, timeSinceSOR);
        break;
      case 16:
        supermoduleHistHelperMeson<16>(minv, timeSinceSOR);
        break;
      case 17:
        supermoduleHistHelperMeson<17>(minv, timeSinceSOR);
        break;
      case 18:
        supermoduleHistHelperMeson<18>(minv, timeSinceSOR);
        break;
      case 19:
        supermoduleHistHelperMeson<19>(minv, timeSinceSOR);
        break;
      default:
        break;
    }
  }

  PresliceUnsortedOptional<o2::aod::EMCALClusters> perCollision = o2::aod::emcalcluster::collisionId;
  PresliceOptional<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollision(MyBCs const& bcs, MyCollisions const& collisions, o2::aod::EMCALClusters const& clusters, o2::soa::Filtered<o2::aod::Calos> const& cells, o2::aod::EMCALClusterCells const& clusterCells)
  {

    auto cellIter = cells.begin();
    auto bcIter = bcs.begin();
    int runNumber = bcIter.runNumber();
    std::unordered_map<uint64_t, int> cellGlobalBCs;
    // Build map of number of cells for corrected BCs using global BCs
    // used later in the determination whether a BC has EMC cell content (for speed reason)
    for (const auto& cell : cells) {
      cellGlobalBCs[cell.bc_as<MyBCs>().globalBC()]++;
    }

    for (const auto& collision : collisions) {
      mHistManager.fill(HIST("events"), 1); // Fill "All events" bin of event histogram

      if (mDoEventSel.value && (!collision.sel8())) { // Check sel8
        continue;
      }
      mHistManager.fill(HIST("events"), 2);                               // Fill sel8
      if (mRequireCaloReadout.value && !collision.alias_bit(kTVXinEMC)) { // Check whether EMC was read out
        continue;
      }
      mHistManager.fill(HIST("events"), 3); // Fill readout

      if (mDoEventSel.value && collision.numContrib() < 0.5) { // Skip collisions without contributors
        continue;
      }
      mHistManager.fill(HIST("events"), 4); // Fill >1 vtx contr. bin of event histogram

      mHistManager.fill(HIST("eventVertexZAll"), collision.posZ());
      if (mVertexCut > 0 && std::abs(collision.posZ()) > mVertexCut) {
        continue;
      }
      mHistManager.fill(HIST("events"), 5); // Fill z-Vertex selected bin of event histogram
      mHistManager.fill(HIST("eventVertexZSelected"), collision.posZ());

      if (mDoEventSel.value && collision.ambiguous()) { // Skip ambiguous collisions (those that are in BCs including multiple collisions)
        continue;
      }
      mHistManager.fill(HIST("events"), 6); // Fill "One collision in BC" bin of event histogram

      if (mDoEventSel.value) {
        auto found = cellGlobalBCs.find(collision.foundBC_as<MyBCs>().globalBC());
        if (mRequireEMCalCells.value && (found == cellGlobalBCs.end() || found->second == 0)) { // Skip collisions without any readout EMCal cells
          continue;
        }
      }
      mHistManager.fill(HIST("events"), 7); // Fill at least one non0 cell in EMCal of event histogram (Selected)

      // Get BC and run number
      int64_t foundBCId = collision.foundBCId();
      if (foundBCId >= 0) {
        bcIter.setCursor(foundBCId);
      }
      runNumber = bcIter.runNumber();

      // Fetch SOR only when run changes
      if (runNumber != lastRunNumber) {
        std::map<std::string, std::string> headers, metadata;
        headers = ccdbApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", runNumber), metadata, -1);
        tsSOR = atol(headers["SOR"].c_str());
        // LOGP(info, "Run {} | SOR = {} ms", runNumber, tsSOR);
        lastRunNumber = runNumber;
      }

      // Time since SOR in minutes (bc.timestamp() is in ms)
      float timeSinceSORMin = (bcIter.timestamp() - tsSOR) / 1000.0f / 60.f;
      mHistManager.fill(HIST("hEventPerTime"), timeSinceSORMin);

      auto clustersPerColl = clusters.sliceBy(perCollision, collision.globalIndex());
      if (clustersPerColl.size() == 0) {
        continue;
      }
      processClusters(clustersPerColl, clusterCells, cellIter, timeSinceSORMin);
      processMesons(timeSinceSORMin);
    }
  }
  PROCESS_SWITCH(EmcalPi0Qc, processCollision, "Process clusters from collision", false);

  /// \brief Process EMCAL clusters that are not matched to a collision
  /// This is not needed for most users
  void processAmbiguous(o2::aod::BCs::iterator const& bc, o2::aod::EMCALAmbiguousClusters const& clusters)
  {
    LOG(debug) << "processAmbiguous";
    // TODO: remove this loop and put it in separate process function that only takes care of ambiguous clusters
    o2::InteractionRecord eventIR;
    eventIR.setFromLong(bc.globalBC());
    mHistManager.fill(HIST("eventBCAll"), eventIR.bc);
    if (std::find(mVetoBCIDs.begin(), mVetoBCIDs.end(), eventIR.bc) != mVetoBCIDs.end()) {
      LOG(info) << "Event rejected because of veto BCID " << eventIR.bc;
      return;
    }
    if (mSelectBCIDs.size() && (std::find(mSelectBCIDs.begin(), mSelectBCIDs.end(), eventIR.bc) == mSelectBCIDs.end())) {
      return;
    }
    mHistManager.fill(HIST("eventBCSelected"), eventIR.bc);

    processAmbiguousClusters(clusters);
    processMesons();
  }
  PROCESS_SWITCH(EmcalPi0Qc, processAmbiguous, "Process Ambiguous clusters", false);

  /// \brief Process EMCAL clusters that are matched to a collisions
  template <o2::soa::is_table Clusters, o2::soa::is_iterator Cell>
  void processClusters(Clusters const& clusters, o2::aod::EMCALClusterCells const& clusterCells, Cell& cellIter, float timeSinceSOR = 0.f)
  {
    LOG(debug) << "processClusters";
    // clear photon vector
    mPhotons.clear();

    for (const auto& cluster : clusters) {
      if (static_cast<int>(clusDef) != cluster.definition()) {
        continue;
      }

      auto cellsPerCluster = clusterCells.sliceBy(perCluster, cluster.globalIndex());
      auto cellsPerClusterIter = cellsPerCluster.begin();
      cellIter.setCursor(cellsPerClusterIter.caloId());
      auto [supermodule, module, phiInModule, etaInModule] = mGeometry->GetCellIndex(cellIter.cellNumber());

      fillClusterQAHistos<decltype(cluster), 0>(cluster);

      if (clusterRejectedByCut(cluster)) {
        continue;
      }

      fillClusterQAHistos<decltype(cluster), 1>(cluster);

      // put clusters in photon vector
      fillSupermoduleHistogramsPhoton(supermodule, cluster.time(), cluster.m02(), cluster.nCells(), timeSinceSOR);
      mPhotons.push_back(Photon(cluster.eta(), cluster.phi(), cluster.energy(), cluster.id(), supermodule));
    }
  }

  /// \brief Process EMCAL clusters that are not matched to a collisions
  template <o2::soa::is_table Clusters>
  void processAmbiguousClusters(Clusters const& clusters)
  {
    LOG(debug) << "processClusters";
    // clear photon vector
    mPhotons.clear();

    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {

      fillClusterQAHistos<decltype(cluster), 0>(cluster);

      if (clusterRejectedByCut(cluster))
        continue;

      fillClusterQAHistos<decltype(cluster), 1>(cluster);

      // put clusters in photon vector
      mPhotons.push_back(Photon(cluster.eta(), cluster.phi(), cluster.energy(), cluster.id()));
    }
  }

  /// \brief Fills the standard QA histograms for a given cluster
  template <o2::soa::is_iterator Cluster, int BeforeCuts>
  void fillClusterQAHistos(Cluster const& cluster)
  {
    // In this implementation the cluster properties are directly loaded from the flat table,
    // in the future one should consider using the AnalysisCluster object to work with after loading.
    static constexpr std::string_view ClusterQAHistEnergy[2] = {"ClustersBeforeCuts/clusterE", "ClustersAfterCuts/clusterE"};
    static constexpr std::string_view ClusterQAHistEnergySimpleBinning[2] = {"ClustersBeforeCuts/clusterE_SimpleBinning", "ClustersAfterCuts/clusterE_SimpleBinning"};
    static constexpr std::string_view ClusterQAHistTime[2] = {"ClustersBeforeCuts/clusterTime", "ClustersAfterCuts/clusterTime"};
    static constexpr std::string_view ClusterQAHistEtaPhi[2] = {"ClustersBeforeCuts/clusterEtaPhi", "ClustersAfterCuts/clusterEtaPhi"};
    static constexpr std::string_view ClusterQAHistM02[2] = {"ClustersBeforeCuts/clusterM02", "ClustersAfterCuts/clusterM02"};
    static constexpr std::string_view ClusterQAHistM20[2] = {"ClustersBeforeCuts/clusterM20", "ClustersAfterCuts/clusterM20"};
    static constexpr std::string_view ClusterQAHistNLM[2] = {"ClustersBeforeCuts/clusterNLM", "ClustersAfterCuts/clusterNLM"};
    static constexpr std::string_view ClusterQAHistNCells[2] = {"ClustersBeforeCuts/clusterNCells", "ClustersAfterCuts/clusterNCells"};
    static constexpr std::string_view ClusterQAHistDistanceToBadChannel[2] = {"ClustersBeforeCuts/clusterDistanceToBadChannel", "ClustersAfterCuts/clusterDistanceToBadChannel"};
    mHistManager.fill(HIST(ClusterQAHistEnergy[BeforeCuts]), cluster.energy());
    mHistManager.fill(HIST(ClusterQAHistEnergySimpleBinning[BeforeCuts]), cluster.energy());
    mHistManager.fill(HIST(ClusterQAHistTime[BeforeCuts]), cluster.time());
    mHistManager.fill(HIST(ClusterQAHistEtaPhi[BeforeCuts]), cluster.eta(), cluster.phi());
    mHistManager.fill(HIST(ClusterQAHistM02[BeforeCuts]), cluster.m02());
    mHistManager.fill(HIST(ClusterQAHistM20[BeforeCuts]), cluster.m20());
    mHistManager.fill(HIST(ClusterQAHistNLM[BeforeCuts]), cluster.nlm());
    mHistManager.fill(HIST(ClusterQAHistNCells[BeforeCuts]), cluster.nCells());
    mHistManager.fill(HIST(ClusterQAHistDistanceToBadChannel[BeforeCuts]), cluster.distanceToBadChannel());
  }

  /// \brief Return a boolean that states, whether a cluster should be rejected by the applied cluster cuts
  template <o2::soa::is_iterator Cluster>
  bool clusterRejectedByCut(Cluster const& cluster)
  {
    // apply basic cluster cuts
    if (cluster.energy() < mMinEnergyCut) {
      LOG(debug) << "Cluster rejected because of energy cut";
      return true;
    }
    if (cluster.nCells() < mMinNCellsCut) {
      LOG(debug) << "Cluster rejected because of nCells cut";
      return true;
    }
    // Only apply M02 cut when cluster contains more than one cell
    if (cluster.nCells() > 1) {
      if (cluster.m02() < mClusterMinM02Cut || cluster.m02() > mClusterMaxM02Cut) {
        LOG(debug) << "Cluster rejected because of m02 cut";
        return true;
      }
    }
    if (cluster.time() < mTimeMin || cluster.time() > mTimeMax) {
      LOG(debug) << "Cluster rejected because of time cut";
      return true;
    }
    return false;
  }

  /// \brief Process meson candidates, calculate invariant mass and pT and fill histograms
  void processMesons(float timeSinceSOR = 0.f)
  {
    LOG(debug) << "processMesons " << mPhotons.size();

    // if less then 2 clusters are found, skip event
    if (mPhotons.size() < 2)
      return;

    // loop over all photon combinations and build meson candidates
    for (unsigned int ig1 = 0; ig1 < mPhotons.size(); ++ig1) {
      for (unsigned int ig2 = ig1 + 1; ig2 < mPhotons.size(); ++ig2) {

        // build meson from photons
        Meson meson(mPhotons[ig1], mPhotons[ig2]);
        if (meson.getOpeningAngle() > mMinOpenAngleCut) {
          mHistManager.fill(HIST("invMassVsPt"), meson.getMass(), meson.getPt());

          uint8_t sm1 = mPhotons[ig1].sm;
          uint8_t sm2 = mPhotons[ig2].sm;
          if (sm1 == sm2) {
            fillSupermoduleHistogramsMeson(sm1, meson.getMass(), timeSinceSOR);
          }

          if (mSplitEMCalDCal) {
            if (!mPhotons[ig1].onDCal && !mPhotons[ig2].onDCal) {
              mHistManager.fill(HIST("invMassVsPt_EMCal"), meson.getMass(), meson.getPt());
            } else if (mPhotons[ig1].onDCal && mPhotons[ig2].onDCal) {
              mHistManager.fill(HIST("invMassVsPt_DCal"), meson.getMass(), meson.getPt());
            }
          }
        }

        // calculate background candidates (rotation background)
        calculateBackground(meson, ig1, ig2);
      }
      calculateMixedBack(mPhotons[ig1]);
    }

    evtMix.addEvent(mPhotons);
  }

  /// \brief Calculate background (using rotation background method)
  void calculateBackground(const Meson& meson, unsigned int ig1, unsigned int ig2)
  {
    // if less than 3 clusters are present, skip event
    if (mPhotons.size() < 3) {
      return;
    }
    const double rotationAngle = o2::constants::math::PIHalf; // 0.78539816339; // rotaion angle 90°

    TLorentzVector lvRotationPhoton1; // photon candidates which get rotated
    TLorentzVector lvRotationPhoton2; // photon candidates which get rotated
    TVector3 lvRotationPion;          // rotation axis
    for (unsigned int ig3 = 0; ig3 < mPhotons.size(); ++ig3) {
      // continue if photons are identical
      if (ig3 == ig1 || ig3 == ig2) {
        continue;
      }
      // calculate rotation axis
      lvRotationPion = (meson.pMeson).Vect();

      // initialize photons for rotation
      lvRotationPhoton1.SetPxPyPzE(mPhotons[ig1].px, mPhotons[ig1].py, mPhotons[ig1].pz, mPhotons[ig1].energy);
      lvRotationPhoton2.SetPxPyPzE(mPhotons[ig2].px, mPhotons[ig2].py, mPhotons[ig2].pz, mPhotons[ig2].energy);

      // rotate photons around rotation axis
      lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
      lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);

      // initialize Photon objects for rotated photons
      Photon rotPhoton1(lvRotationPhoton1.Eta(), lvRotationPhoton1.Phi(), lvRotationPhoton1.E(), mPhotons[ig1].id);
      Photon rotPhoton2(lvRotationPhoton2.Eta(), lvRotationPhoton2.Phi(), lvRotationPhoton2.E(), mPhotons[ig2].id);

      // build meson from rotated photons
      Meson mesonRotated1(rotPhoton1, mPhotons[ig3]);
      Meson mesonRotated2(rotPhoton2, mPhotons[ig3]);

      // Fill histograms
      if (mesonRotated1.getOpeningAngle() > mMinOpenAngleCut) {
        mHistManager.fill(HIST("invMassVsPtBackground"), mesonRotated1.getMass(), mesonRotated1.getPt());
        if (mSplitEMCalDCal) {
          if (!mPhotons[ig1].onDCal && !mPhotons[ig2].onDCal && !mPhotons[ig3].onDCal) {
            mHistManager.fill(HIST("invMassVsPtBackground_EMCal"), mesonRotated1.getMass(), mesonRotated1.getPt());
          } else if (mPhotons[ig1].onDCal && mPhotons[ig2].onDCal && mPhotons[ig3].onDCal) {
            mHistManager.fill(HIST("invMassVsPtBackground_DCal"), mesonRotated1.getMass(), mesonRotated1.getPt());
          }
        }
      }
      if (mesonRotated2.getOpeningAngle() > mMinOpenAngleCut) {
        mHistManager.fill(HIST("invMassVsPtBackground"), mesonRotated2.getMass(), mesonRotated2.getPt());
        if (mSplitEMCalDCal) {
          if (!mPhotons[ig1].onDCal && !mPhotons[ig2].onDCal && !mPhotons[ig3].onDCal) {
            mHistManager.fill(HIST("invMassVsPtBackground_EMCal"), mesonRotated2.getMass(), mesonRotated2.getPt());
          } else if (mPhotons[ig1].onDCal && mPhotons[ig2].onDCal && mPhotons[ig3].onDCal) {
            mHistManager.fill(HIST("invMassVsPtBackground_DCal"), mesonRotated2.getMass(), mesonRotated2.getPt());
          }
        }
      }
    }
  }

  void calculateMixedBack(Photon gamma)
  {
    for (unsigned int i = 0; i < evtMix.vecEvtMix.size(); ++i) {
      for (unsigned int ig1 = 0; ig1 < evtMix.vecEvtMix[i].size(); ++ig1) {
        Meson meson(gamma, evtMix.vecEvtMix[i][ig1]);
        if (meson.getOpeningAngle() > mMinOpenAngleCut) {
          mHistManager.fill(HIST("invMassVsPtMixedBackground"), meson.getMass(), meson.getPt());
          if (mSplitEMCalDCal) {
            if (!gamma.onDCal && !evtMix.vecEvtMix[i][ig1].onDCal) {
              mHistManager.fill(HIST("invMassVsPtMixedBackground_EMCal"), meson.getMass(), meson.getPt());
            } else if (gamma.onDCal && evtMix.vecEvtMix[i][ig1].onDCal) {
              mHistManager.fill(HIST("invMassVsPtMixedBackground_DCal"), meson.getMass(), meson.getPt());
            }
          }
        }
      }
    }
  }

  /// \brief Create binning for cluster energy/pT axis (variable bin size)
  /// direct port from binning often used in AliPhysics for debugging
  /// \return vector with bin limits
  std::vector<double> makeClusterBinning() const
  {
    std::vector<double> result;
    int nBinsPt = 179;
    double maxPt = 60;
    for (int i = 0; i < nBinsPt + 1; i++) {
      if (i < 100) {
        result.emplace_back(0.10 * i);
      } else if (i < 140) {
        result.emplace_back(10. + 0.25 * (i - 100));
      } else if (i < 180) {
        result.emplace_back(20. + 1.00 * (i - 140));
      } else {
        result.emplace_back(maxPt);
      }
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<EmcalPi0Qc>(cfgc, TaskName{"EmcalPi0QcAssociate"}, SetDefaultProcesses{{{"processCollision", true}, {"processAmbiguous", false}}}),  // o2-linter: disable=name/o2-task (adapted multiple times)
    adaptAnalysisTask<EmcalPi0Qc>(cfgc, TaskName{"EmcalPi0QcAmbiguous"}, SetDefaultProcesses{{{"processCollision", false}, {"processAmbiguous", true}}})}; // o2-linter: disable=name/o2-task (adapted multiple times)
  return workflow;
}
