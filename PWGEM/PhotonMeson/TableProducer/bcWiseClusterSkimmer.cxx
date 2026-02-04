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
///
/// \file bcWiseClusterSkimmer.cxx
///
/// \brief This task creates minimalistic skimmed tables containing EMC clusters and centrality information
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#include "PWGEM/PhotonMeson/DataModel/bcWiseTables.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DetectorsBase/GeometryManager.h"
#include "EMCALBase/Geometry.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod::emdownscaling;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms>;
using MyMCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::McCollisionLabels>;
using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps, aod::BCCentFT0Cs, aod::BCCentFT0Ms, aod::Run3MatchedToBCSparse>;

using SelectedUniqueClusters = soa::Filtered<aod::EMCALClusters>;                                                         // Clusters from collisions with only one collision in the BC
using SelectedUniqueMCClusters = soa::Filtered<soa::Join<aod::EMCALClusters, aod::EMCALMCClusters>>;                      // Clusters from collisions with only one collision in the BC
using SelectedAmbiguousClusters = soa::Filtered<aod::EMCALAmbiguousClusters>;                                             // Clusters from BCs with multiple collisions (no vertex assignment possible)
using SelectedAmbiguousMCClusters = soa::Filtered<soa::Join<aod::EMCALAmbiguousClusters, aod::EMCALAmbiguousMCClusters>>; // Clusters from BCs with multiple collisions (no vertex assignment possible)
using SelectedCells = o2::soa::Filtered<aod::Calos>;

struct bcWiseClusterSkimmer {
  Produces<aod::BCWiseBCs> bcTable;
  Produces<aod::BCWiseClusters> clusterTable;
  Produces<aod::BCWiseCollisions> collisionTable;
  Produces<aod::BCWiseMCPi0s> mcpi0Table;
  Produces<aod::BCWiseMCEtas> mcetaTable;
  Produces<aod::BCWiseMCClusters> mcclusterTable;

  PresliceUnsorted<MyCollisions> perFoundBC = aod::evsel::foundBCId;
  Preslice<SelectedUniqueClusters> perCol = aod::emcalcluster::collisionId;
  Preslice<SelectedAmbiguousClusters> perBC = aod::emcalcluster::bcId;
  Preslice<SelectedCells> cellsPerBC = aod::calo::bcId;

  Configurable<float> cfgMinClusterEnergy{"cfgMinClusterEnergy", 0.5, "Minimum energy of selected clusters (GeV)"};
  Configurable<float> cfgMaxClusterEnergy{"cfgMaxClusterEnergy", 30., "Maximum energy of selected clusters (GeV)"};
  Configurable<float> cfgMinM02{"cfgMinM02", -1., "Minimum M02 of selected clusters"};
  Configurable<float> cfgMaxM02{"cfgMaxM02", 5., "Maximum M02 of selected clusters"};
  Configurable<float> cfgMinTime{"cfgMinTime", -25, "Minimum time of selected clusters (ns)"};
  Configurable<float> cfgMaxTime{"cfgMaxTime", 25, "Maximum time of selected clusters (ns)"};
  Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.8f, "Maximum absolute rapidity of counted generated particles"};
  Configurable<float> cfgMinPtGen{"cfgMinPtGen", 0., "Minimum pT for stored generated mesons (reduce disk space of derived data)"};

  Configurable<bool> cfgRequirekTVXinEMC{"cfgRequirekTVXinEMC", false, "Only store kTVXinEMC triggered BCs"};
  Configurable<bool> cfgRequireGoodRCTQuality{"cfgRequireGoodRCTQuality", false, "Only store BCs with good quality of T0 and EMC in RCT"};
  Configurable<bool> cfgStoreMu{"cfgStoreMu", false, "Calculate and store mu (probablity of a TVX collision in the BC) per BC. Otherwise fill with 0"};
  Configurable<bool> cfgStoreTime{"cfgStoreTime", false, "Calculate and store time since the start of fill. Otherwise fill with 0"};
  Configurable<bool> cfgOnlyCheckFirstTrueContributor{"cfgOnlyCheckFirstTrueContributor", false, "When storing MC cluster info only count as true if the leading contributor is from pi0/eta"};
  ConfigurableAxis cfgMultiplicityBinning{"cfgMultiplicityBinning", {1000, 0, 10000}, "Binning used for the binning of the number of particles in the event"};

  aod::rctsel::RCTFlagsChecker isFT0EMCGoodRCTChecker{aod::rctsel::kFT0Bad, aod::rctsel::kEMCBad};
  parameters::GRPLHCIFData* mLHCIFdata = nullptr;
  int mRunNumber = -1;
  ctpRateFetcher mRateFetcher;

  Filter energyFilter = (aod::emcalcluster::energy > cfgMinClusterEnergy && aod::emcalcluster::energy < cfgMaxClusterEnergy);
  Filter m02Filter = (aod::emcalcluster::nCells == 1 || (aod::emcalcluster::m02 > cfgMinM02 && aod::emcalcluster::m02 < cfgMaxM02));
  Filter timeFilter = (aod::emcalcluster::time > cfgMinTime && aod::emcalcluster::time < cfgMaxTime);
  Filter emccellFilter = aod::calo::caloType == 1;

  HistogramRegistry mHistManager{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  std::map<int32_t, int32_t> fMapPi0Index; // Map to connect the MC index of the pi0 to the one saved in the derived table
  std::map<int32_t, int32_t> fMapEtaIndex; // Map to connect the MC index of the eta to the one saved in the derived table

  void init(o2::framework::InitContext&)
  {
    const int nEventBins = 6;
    mHistManager.add("nBCs", "Number of BCs;;#bf{#it{N}_{BCs}}", HistType::kTH1F, {{nEventBins, -0.5, 5.5}});
    const TString binLabels[nEventBins] = {"All", "FT0", "TVX", "kTVXinEMC", "Cell", "HasMCColl"};
    for (int iBin = 0; iBin < nEventBins; iBin++)
      mHistManager.get<TH1>(HIST("nBCs"))->GetXaxis()->SetBinLabel(iBin + 1, binLabels[iBin]);

    mHistManager.add("CentralityVsGenMultiplicity", "Centrality vs number of generated MC particles;Centrality;#bf{#it{N}_{gen}}", HistType::kTH2F, {{102, 0., 102}, cfgMultiplicityBinning});
    mHistManager.add("CentralityFT0CVsM", "Centrality of FT0C vs FT0M;FT0C Cent;FT0M Cent", HistType::kTH2F, {{102, 0., 102}, {102, 0., 102}});
    mHistManager.add("BCCentVsCollCent", "Centrality of the BC vs Centrality of the collision;BC Centrality;Collision Centrality", HistType::kTH2F, {{102, 0., 102}, {102, 0., 102}});
    mHistManager.add("BCHasFT0FoundFT0", "BCHasFT0FoundFT0;BC has FT0;BC has found FT0", HistType::kTH2F, {{1000, -2., 50000}, {1000, -2., 50000}});
    mHistManager.add("BCFoundFT0Diff", "BCFoundFT0Diff;BC_{ID}^{found FT0} - BC_{ID}^{FT0}", HistType::kTH1F, {{200001, -100000.5, 100000.5}});
    mHistManager.add("BCIDOfMCColl", "BCIDOfMCColl;BCID", HistType::kTH1F, {{4001, -0.5, 4000.5}});
    mHistManager.add("BCIDOfMCCollButNoFT0", "BCIDOfMCCollButNoFT0;BCID", HistType::kTH1F, {{4001, -0.5, 4000.5}});

    mHistManager.add("MCZVtx", "Vertex position of MC collision;#bf{z_{vtx}^{MC} (cm)};#bf{#it{N}_{BCs}}", HistType::kTH1F, {{2000, -100, 100}});
    mHistManager.add("MCZVtx_TVX", "Vertex position of MC collision;#bf{z_{vtx}^{MC} (cm)};#bf{#it{N}_{BCs}}", HistType::kTH1F, {{2000, -100, 100}});

    LOG(info) << "BC wise cluster skimmer cuts:";
    LOG(info) << "------------------------------------";
    LOG(info) << "| Timing cut: " << cfgMinTime << " < t < " << cfgMaxTime;
    LOG(info) << "| Shape cut: " << cfgMinM02 << " < M02 < " << cfgMaxM02;
    LOG(info) << "| Energy cut: " << cfgMinClusterEnergy << " < E < " << cfgMaxClusterEnergy;
    LOG(info) << "| Rapidity cut: |y| < " << cfgRapidityCut;
    LOG(info) << "| Min gen pt: pT > " << cfgMinPtGen;

    o2::emcal::Geometry::GetInstanceFromRunNumber(300000);
    if (cfgRequireGoodRCTQuality)
      isFT0EMCGoodRCTChecker.init({aod::rctsel::kFT0Bad, aod::rctsel::kEMCBad});
  }

  /// \brief Process EMCAL clusters (either ambigous or unique)
  template <typename OutputType, typename InputType>
  OutputType convertForStorage(InputType const& valueIn, Observable observable)
  {
    double valueToBeChecked = valueIn * downscalingFactors[observable];
    if (valueToBeChecked < std::numeric_limits<OutputType>::lowest()) {
      LOG(warning) << "Value " << valueToBeChecked << " of observable " << observable << " below lowest possible value of " << typeid(OutputType).name() << ": " << static_cast<float>(std::numeric_limits<OutputType>::lowest());
      valueToBeChecked = static_cast<float>(std::numeric_limits<OutputType>::lowest());
    }
    if (valueToBeChecked > std::numeric_limits<OutputType>::max()) {
      LOG(warning) << "Value " << valueToBeChecked << " of observable " << observable << " obove highest possible value of " << typeid(OutputType).name() << ": " << static_cast<float>(std::numeric_limits<OutputType>::max());
      valueToBeChecked = static_cast<float>(std::numeric_limits<OutputType>::max());
    }

    return static_cast<OutputType>(valueToBeChecked);
  }

  /// \brief Process EMCAL clusters (either ambigous or unique)
  template <o2::soa::is_table Clusters>
  void processClusters(Clusters const& clusters, const int bcID)
  {
    for (const auto& cluster : clusters) {
      clusterTable(bcID,
                   convertForStorage<int8_t>(cluster.definition(), kDefinition),
                   convertForStorage<int16_t>(cluster.energy(), kEnergy),
                   convertForStorage<int16_t>(cluster.eta(), kEta),
                   convertForStorage<uint16_t>(cluster.phi(), kPhi),
                   convertForStorage<int8_t>(cluster.nCells(), kNCells),
                   convertForStorage<int16_t>(cluster.m02(), kM02),
                   convertForStorage<int16_t>(cluster.time(), kTime),
                   cluster.isExotic());
    }
  }

  template <o2::soa::is_table Clusters>
  void processClusterMCInfo(Clusters const& clusters, const int bcID, aod::McParticles const& mcParticles)
  {
    for (const auto& cluster : clusters) {
      float clusterInducerEnergy = 0.;
      int32_t mesonMCIndex = -1;
      for (size_t i = 0; i < cluster.amplitudeA().size(); i++) {
        int clusterInducerId = cluster.mcParticleIds()[i];
        auto clusterInducer = mcParticles.iteratorAt(clusterInducerId);
        clusterInducerEnergy = clusterInducer.e();
        int daughterId = aod::pwgem::photonmeson::utils::mcutil::FindMotherInChain(clusterInducer, mcParticles, std::vector<int>{111, 221});
        if (daughterId > 0) {
          mesonMCIndex = mcParticles.iteratorAt(daughterId).mothersIds()[0];
          if (mcParticles.iteratorAt(mesonMCIndex).pt() < cfgMinPtGen)
            mesonMCIndex = -1;
        }
        if (mesonMCIndex != -1)
          break;
        if (cfgOnlyCheckFirstTrueContributor)
          break;
      }
      bool isEta = false;
      if (mesonMCIndex >= 0) {
        if (mcParticles.iteratorAt(mesonMCIndex).pdgCode() == 111) {
          if (fMapPi0Index.find(mesonMCIndex) != fMapPi0Index.end()) // Some pi0s might not be found (not gg decay or too large y)
            mesonMCIndex = fMapPi0Index[mesonMCIndex];               // If pi0 was stored in table, change index from the MC index to the pi0 index from this task
          else                                                       // If pi0 was not stored, treat photon as if not from pi0
            mesonMCIndex = -1;
        } else if (mcParticles.iteratorAt(mesonMCIndex).pdgCode() == 221) {
          isEta = true;
          if (fMapEtaIndex.find(mesonMCIndex) != fMapEtaIndex.end()) // Some etas might not be found (not gg decay or too large y)
            mesonMCIndex = fMapEtaIndex[mesonMCIndex];               // If eta was stored in table, change index from the MC index to the eta index from this task
          else                                                       // If eta was not stored, treat photon as if not from eta
            mesonMCIndex = -1;
        } else {
          mesonMCIndex = -1; // Not a pi0 or eta
        }
      }
      mcclusterTable(bcID, mesonMCIndex, isEta, convertForStorage<uint16_t>(clusterInducerEnergy, kEnergy));
    }
  }

  bool isBCSelected(const auto& bc)
  {
    if (cfgRequirekTVXinEMC && !bc.alias_bit(kTVXinEMC))
      return false;
    if (cfgRequireGoodRCTQuality && !isFT0EMCGoodRCTChecker(bc))
      return false;
    return true;
  }

  void setLHCIFData(const auto& bc)
  {
    if (mRunNumber == bc.runNumber())
      return;

    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    uint64_t timeStamp = bc.timestamp();

    std::map<std::string, std::string> metadata;
    mLHCIFdata = ccdbMgr.getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
    if (mLHCIFdata == nullptr)
      LOG(fatal) << "GRPLHCIFData not in database, timestamp:" << timeStamp;
    mRunNumber = bc.runNumber();
    LOG(info) << "LHCIF data fetched for run " << mRunNumber << " and timestamp " << timeStamp;

    return;
  }

  double calculateMu(const auto& bc)
  {
    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    uint64_t timeStamp = bc.timestamp();

    auto bfilling = mLHCIFdata->getBunchFilling();
    double nbc = bfilling.getFilledBCs().size();
    double tvxRate = mRateFetcher.fetch(&ccdbMgr, timeStamp, bc.runNumber(), "T0VTX");
    double nTriggersPerFilledBC = tvxRate / nbc / o2::constants::lhc::LHCRevFreq;
    double mu = -std::log(1 - nTriggersPerFilledBC);

    // LOG(info) << "Time stamp: " << timeStamp << " Run number: " << bc.runNumber() << " Number of filled BCs: " << nbc << " Trigger rate: " << tvxRate << " Mu: " << mu;

    return mu;
  }

  void processEventProperties(const auto& bc, const auto& collisionsInBC, const auto& cellsInBC, const bool hasMCCollision = false)
  {
    bool hasFT0 = bc.has_ft0();
    bool hasFoundFT0 = bc.has_foundFT0();
    int foundFT0ID = hasFoundFT0 ? bc.foundFT0().bcId() : -1;
    int FT0ID = hasFT0 ? bc.ft0().bcId() : -1;
    mHistManager.fill(HIST("BCHasFT0FoundFT0"), static_cast<float>(FT0ID), static_cast<float>(foundFT0ID));
    mHistManager.fill(HIST("BCFoundFT0Diff"), static_cast<float>(foundFT0ID - FT0ID));

    auto bcIDInOrbit = bc.globalBC() % 3564;

    if (hasMCCollision) {
      mHistManager.fill(HIST("BCIDOfMCColl"), static_cast<float>(bcIDInOrbit));
      if (!hasFoundFT0) {
        mHistManager.fill(HIST("BCIDOfMCCollButNoFT0"), static_cast<float>(bcIDInOrbit));
      }
    }

    bool hasTVX = bc.selection_bit(aod::evsel::kIsTriggerTVX);
    bool haskTVXinEMC = bc.alias_bit(kTVXinEMC);
    bool hasEMCCell = cellsInBC.size() > 0;
    mHistManager.fill(HIST("nBCs"), 0);
    if (hasFoundFT0)
      mHistManager.fill(HIST("nBCs"), 1);
    if (hasTVX)
      mHistManager.fill(HIST("nBCs"), 2);
    if (haskTVXinEMC)
      mHistManager.fill(HIST("nBCs"), 3);
    if (hasEMCCell)
      mHistManager.fill(HIST("nBCs"), 4);
    if (hasMCCollision)
      mHistManager.fill(HIST("nBCs"), 5);

    if (cfgStoreMu || cfgStoreTime)
      setLHCIFData(bc);
    double mu = cfgStoreMu ? calculateMu(bc) : 0.;
    float timeSinceSOF = cfgStoreTime ? (bc.timestamp() - mLHCIFdata->getFillNumberTime()) / 1e3 : 0.; // Convert to seconds
    float ft0Amp = hasFoundFT0 ? bc.foundFT0().sumAmpA() + bc.foundFT0().sumAmpC() : 0.;
    double centralityOfCollision = 101.5;
    if (collisionsInBC.size() > 0)
      centralityOfCollision = collisionsInBC.iteratorAt(0).centFT0C();
    double centralityOfBC = bc.centFT0C();

    mHistManager.fill(HIST("BCCentVsCollCent"), centralityOfBC, centralityOfCollision);
    mHistManager.fill(HIST("CentralityFT0CVsM"), centralityOfBC, bc.centFT0M());

    bcTable(hasFoundFT0, hasTVX, haskTVXinEMC, hasEMCCell, convertForStorage<uint8_t>(centralityOfBC, kCent), convertForStorage<uint8_t>(bc.centFT0M(), kCent), convertForStorage<uint16_t>(ft0Amp, kFT0Amp), convertForStorage<uint16_t>(mu, kMu), convertForStorage<uint16_t>(timeSinceSOF, kTimeSinceSOF));

    for (const auto& collision : collisionsInBC)
      collisionTable(bcTable.lastIndex(), convertForStorage<int16_t>(collision.posZ(), kZVtx));
  }

  template <o2::soa::is_iterator TMCParticle, o2::soa::is_table TMCParticles>
  bool isGammaGammaDecay(TMCParticle const& mcParticle, TMCParticles const& mcParticles)
  {
    auto daughtersIds = mcParticle.daughtersIds();
    if (daughtersIds.size() != 2)
      return false;
    for (const auto& daughterId : daughtersIds) {
      if (mcParticles.iteratorAt(daughterId).pdgCode() != 22)
        return false;
    }
    return true;
  }

  template <o2::soa::is_iterator TMCParticle, o2::soa::is_table TMCParticles>
  bool isAccepted(TMCParticle const& mcParticle, TMCParticles const& mcParticles)
  {
    auto daughtersIds = mcParticle.daughtersIds();
    if (daughtersIds.size() != 2)
      return false;
    for (const auto& daughterId : daughtersIds) {
      if (mcParticles.iteratorAt(daughterId).pdgCode() != 22)
        return false;
      int iCellID = -1;
      try {
        iCellID = emcal::Geometry::GetInstance()->GetAbsCellIdFromEtaPhi(mcParticles.iteratorAt(daughterId).eta(), mcParticles.iteratorAt(daughterId).phi());
      } catch (emcal::InvalidPositionException& e) {
        iCellID = -1;
      }
      if (iCellID == -1)
        return false;
    }
    return true;
  }

  void processData(MyBCs const& bcs, MyCollisions const& collisions, aod::FT0s const&, SelectedCells const& cells, SelectedUniqueClusters const& uClusters, SelectedAmbiguousClusters const& aClusters)
  {
    for (const auto& bc : bcs) {
      if (!isBCSelected(bc))
        continue;
      auto collisionsInBC = collisions.sliceBy(perFoundBC, bc.globalIndex());
      auto cellsInBC = cells.sliceBy(cellsPerBC, bc.globalIndex());

      processEventProperties(bc, collisionsInBC, cellsInBC);

      if (collisionsInBC.size() == 1) {
        auto clustersInBC = uClusters.sliceBy(perCol, collisionsInBC.begin().globalIndex());
        processClusters(clustersInBC, bcTable.lastIndex());
      } else {
        auto clustersInBC = aClusters.sliceBy(perBC, bc.globalIndex());
        processClusters(clustersInBC, bcTable.lastIndex());
      }
    }
  }
  PROCESS_SWITCH(bcWiseClusterSkimmer, processData, "Run skimming for data", true);

  Preslice<aod::McCollisions> mcCollperBC = aod::mccollision::bcId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  void processMC(MyBCs const& bcs, MyMCCollisions const& collisions, aod::McCollisions const& mcCollisions, aod::FT0s const&, SelectedCells const& cells, SelectedUniqueMCClusters const& uClusters, SelectedAmbiguousMCClusters const& aClusters, aod::McParticles const& mcParticles)
  {
    for (const auto& bc : bcs) {
      if (!isBCSelected(bc))
        continue;
      auto collisionsInBC = collisions.sliceBy(perFoundBC, bc.globalIndex());
      auto cellsInBC = cells.sliceBy(cellsPerBC, bc.globalIndex());

      auto mcCollisionsBC = mcCollisions.sliceBy(mcCollperBC, bc.globalIndex());

      processEventProperties(bc, collisionsInBC, cellsInBC, mcCollisionsBC.size() > 0);

      for (const auto& mcCollision : mcCollisionsBC) {
        mHistManager.fill(HIST("MCZVtx"), mcCollision.posZ());
        if (bc.selection_bit(aod::evsel::kIsTriggerTVX))
          mHistManager.fill(HIST("MCZVtx_TVX"), mcCollision.posZ());
        auto mcParticlesInColl = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());
        mHistManager.fill(HIST("CentralityVsGenMultiplicity"), bc.centFT0M(), mcParticlesInColl.size());
        for (const auto& mcParticle : mcParticlesInColl) {
          if (std::abs(mcParticle.y()) > cfgRapidityCut || !isGammaGammaDecay(mcParticle, mcParticles) || mcParticle.pt() < cfgMinPtGen)
            continue;
          bool isPrimary = mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator();
          bool isFromWD = (aod::pwgem::photonmeson::utils::mcutil::IsFromWD(mcCollision, mcParticle, mcParticles)) > 0;

          if (mcParticle.pdgCode() == 111) {
            mcpi0Table(bcTable.lastIndex(), convertForStorage<uint16_t>(mcParticle.pt(), kpT), isAccepted(mcParticle, mcParticles), isPrimary, isFromWD);
            fMapPi0Index[mcParticle.globalIndex()] = static_cast<int32_t>(mcpi0Table.lastIndex());
          } else if (mcParticle.pdgCode() == 221) {
            mcetaTable(bcTable.lastIndex(), convertForStorage<uint16_t>(mcParticle.pt(), kpT), isAccepted(mcParticle, mcParticles), isPrimary, isFromWD);
            fMapEtaIndex[mcParticle.globalIndex()] = static_cast<int32_t>(mcetaTable.lastIndex());
          }
        }
      }

      if (collisionsInBC.size() == 1) {
        auto clustersInBC = uClusters.sliceBy(perCol, collisionsInBC.begin().globalIndex());
        processClusters(clustersInBC, bcTable.lastIndex());
        processClusterMCInfo(clustersInBC, bcTable.lastIndex(), mcParticles);
      } else {
        auto clustersInBC = aClusters.sliceBy(perBC, bc.globalIndex());
        processClusters(clustersInBC, bcTable.lastIndex());
        processClusterMCInfo(clustersInBC, bcTable.lastIndex(), mcParticles);
      }
      fMapPi0Index.clear();
      fMapEtaIndex.clear();
    }
  }
  PROCESS_SWITCH(bcWiseClusterSkimmer, processMC, "Run skimming for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<bcWiseClusterSkimmer>(cfgc)}; }
