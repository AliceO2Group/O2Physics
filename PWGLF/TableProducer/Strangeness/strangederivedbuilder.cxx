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
//__________________________________________________
// this task provides general links between collisions
// and strange objects reconstructed in various ways.
// It is meant to help with providing auxiliary information
// when dealing with derived data.

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
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Qvectors.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using TracksWithExtra = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullHe, aod::TOFEvTime, aod::TOFSignal>;
using TracksCompleteIUMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::McTrackLabels>;
using FullTracksExtIUTOF = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFEvTime, aod::TOFSignal>;
using FullCollisions = soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::FT0Mults>;
using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDCollsLabels>;

// simple bit checkers
#define bitset(var, nbit) ((var) |= (1 << (nbit)))
#define bitcheck(var, nbit) ((var) & (1 << (nbit)))

struct strangederivedbuilder {
  SliceCache cache;

  //__________________________________________________
  // fundamental building blocks of derived data
  Produces<aod::StraCollision> strangeColl;        // characterises collisions
  Produces<aod::StraCollLabels> strangeCollLabels; // characterises collisions
  Produces<aod::StraMCCollisions> strangeMCColl;   // characterises collisions / MC
  Produces<aod::StraMCCollMults> strangeMCMults;   // characterises collisions / MC mults
  Produces<aod::StraCents> strangeCents;           // characterises collisions / centrality
  Produces<aod::StraEvSels> strangeEvSels;         // characterises collisions / centrality / sel8 selection
  Produces<aod::StraStamps> strangeStamps;         // provides timestamps, run numbers
  Produces<aod::V0CollRefs> v0collref;             // references collisions from V0s
  Produces<aod::CascCollRefs> casccollref;         // references collisions from cascades
  Produces<aod::KFCascCollRefs> kfcasccollref;     // references collisions from KF cascades
  Produces<aod::TraCascCollRefs> tracasccollref;   // references collisions from tracked cascades

  //__________________________________________________
  // track extra references
  Produces<aod::DauTrackExtras> dauTrackExtras;   // daughter track detector properties
  Produces<aod::DauTrackMCIds> dauTrackMCIds;     // daughter track MC Particle ID
  Produces<aod::DauTrackTPCPIDs> dauTrackTPCPIDs; // daughter track TPC PID
  Produces<aod::DauTrackTOFPIDs> dauTrackTOFPIDs; // daughter track TOF PID
  Produces<aod::V0Extras> v0Extras;               // references DauTracks from V0s
  Produces<aod::CascExtras> cascExtras;           // references DauTracks from cascades
  Produces<aod::StraTrackExtras> straTrackExtras; // references DauTracks from tracked cascades

  //__________________________________________________
  // cascade interlinks
  Produces<aod::CascToTraRefs> cascToTraRefs; // cascades -> tracked
  Produces<aod::CascToKFRefs> cascToKFRefs;   // cascades -> KF
  Produces<aod::TraToCascRefs> traToCascRefs; // tracked -> cascades
  Produces<aod::KFToCascRefs> kfToCascRefs;   // KF -> cascades

  //__________________________________________________
  // mother information
  Produces<aod::V0MCMothers> v0mothers;       // V0 mother references
  Produces<aod::CascMCMothers> cascmothers;   // casc mother references
  Produces<aod::MotherMCParts> motherMCParts; // mc particles for mothers

  //__________________________________________________
  // Q-vectors
  Produces<aod::StraFT0AQVs> StraFT0AQVs;     // FT0A Q-vector
  Produces<aod::StraFT0CQVs> StraFT0CQVs;     // FT0C Q-vector
  Produces<aod::StraFT0MQVs> StraFT0MQVs;     // FT0M Q-vector
  Produces<aod::StraFV0AQVs> StraFV0AQVs;     // FV0A Q-vector
  Produces<aod::StraTPCQVs> StraTPCQVs;       // TPC Q-vector
  Produces<aod::StraFT0CQVsEv> StraFT0CQVsEv; // events used to compute FT0C Q-vector (LF)

  //__________________________________________________
  // Generated binned data
  // this is a hack while the system does not do better
  Produces<aod::GeK0Short> geK0Short;
  Produces<aod::GeLambda> geLambda;
  Produces<aod::GeAntiLambda> geAntiLambda;
  Produces<aod::GeXiMinus> geXiMinus;
  Produces<aod::GeXiPlus> geXiPlus;
  Produces<aod::GeOmegaMinus> geOmegaMinus;
  Produces<aod::GeOmegaPlus> geOmegaPlus;

  //__________________________________________________
  // Found tags for findable exercise
  Produces<aod::V0FoundTags> v0FoundTags;
  Produces<aod::CascFoundTags> cascFoundTags;

  // histogram registry for bookkeeping
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr int nSpecies = 14;
  static constexpr int nParameters = 1;
  static const std::vector<std::string> particleNames;
  static const std::vector<int> particlePDGCodes;
  static const std::vector<std::string> parameterNames;
  static const int defaultParameters[nSpecies][nParameters];
  static constexpr std::string_view particleNamesConstExpr[] = {"Gamma", "K0Short", "Lambda", "AntiLambda",
                                                                "Sigma0", "AntiSigma0", "SigmaPlus", "SigmaMinus",
                                                                "Hypertriton", "AntiHypertriton",
                                                                "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};

  uint32_t enabledBits = 0;

  Configurable<LabeledArray<int>> enableGeneratedInfo{"enableGeneratedInfo",
                                                      {defaultParameters[0], nSpecies,
                                                       nParameters, particleNames, parameterNames},
                                                      "Fill generated particle histograms for each species. 0: no, 1: yes"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisCentrality{"axisCentrality", {100, 0.0f, 100.0f}, "Centrality"};
  ConfigurableAxis axisRawCentrality{"axisRawCentrality", {VARIABLE_WIDTH, 0.000f, 52.320f, 75.400f, 95.719f, 115.364f, 135.211f, 155.791f, 177.504f, 200.686f, 225.641f, 252.645f, 281.906f, 313.850f, 348.302f, 385.732f, 426.307f, 470.146f, 517.555f, 568.899f, 624.177f, 684.021f, 748.734f, 818.078f, 892.577f, 973.087f, 1058.789f, 1150.915f, 1249.319f, 1354.279f, 1465.979f, 1584.790f, 1710.778f, 1844.863f, 1985.746f, 2134.643f, 2291.610f, 2456.943f, 2630.653f, 2813.959f, 3006.631f, 3207.229f, 3417.641f, 3637.318f, 3865.785f, 4104.997f, 4354.938f, 4615.786f, 4885.335f, 5166.555f, 5458.021f, 5762.584f, 6077.881f, 6406.834f, 6746.435f, 7097.958f, 7462.579f, 7839.165f, 8231.629f, 8635.640f, 9052.000f, 9484.268f, 9929.111f, 10389.350f, 10862.059f, 11352.185f, 11856.823f, 12380.371f, 12920.401f, 13476.971f, 14053.087f, 14646.190f, 15258.426f, 15890.617f, 16544.433f, 17218.024f, 17913.465f, 18631.374f, 19374.983f, 20136.700f, 20927.783f, 21746.796f, 22590.880f, 23465.734f, 24372.274f, 25314.351f, 26290.488f, 27300.899f, 28347.512f, 29436.133f, 30567.840f, 31746.818f, 32982.664f, 34276.329f, 35624.859f, 37042.588f, 38546.609f, 40139.742f, 41837.980f, 43679.429f, 45892.130f, 400000.000f}, "raw centrality signal"}; // for QA

  ConfigurableAxis axisNVertices{"axisNVertices", {10, -0.5f, 9.5f}, "N(vertices)"};

  Configurable<bool> fillEmptyCollisions{"fillEmptyCollisions", false, "fill collision entries without candidates"};

  // round Nsigma variables up to a certain level of precision if requested
  // useful to keep derived data sizes under control
  // variables that are rounded include the DCAs but not the CosPA (precision needed)
  Configurable<bool> roundNSigmaVariables{"roundNSigmaVariables", false, "round NSigma variables"};
  Configurable<float> precisionNSigmas{"precisionNSigmas", 0.1f, "precision to keep NSigmas"};

  Configurable<bool> fillRawFT0A{"fillRawFT0A", false, "Fill raw FT0A information for debug"};
  Configurable<bool> fillRawFT0C{"fillRawFT0C", true, "Fill raw FT0C information for debug"};
  Configurable<bool> fillRawFV0A{"fillRawFV0A", false, "Fill raw FV0A information for debug"};
  Configurable<bool> fillRawFDDA{"fillRawFDDA", false, "Fill raw FDDA information for debug"};
  Configurable<bool> fillRawFDDC{"fillRawFDDC", false, "Fill raw FDDC information for debug"};
  Configurable<bool> fillRawZDC{"fillRawZDC", false, "Fill raw ZDC information for debug"};
  Configurable<bool> fillRawNTracksEta1{"fillRawNTracksEta1", true, "Fill raw NTracks |eta|<1 information for debug"};
  Configurable<bool> fillRawNTracksForCorrelation{"fillRawNTracksForCorrelation", true, "Fill raw NTracks for correlation cuts"};
  Configurable<bool> fillTOFInformation{"fillTOFInformation", true, "Fill Daughter Track TOF information"};

  Configurable<bool> qaCentrality{"qaCentrality", false, "qa centrality flag: check base raw values"};
  struct : ConfigurableGroup {
    ConfigurableAxis axisFT0A{"FT0Aamplitude", {100, 0.0f, 2000.0f}, "FT0Aamplitude"};
    ConfigurableAxis axisFT0C{"FT0Camplitude", {100, 0.0f, 2000.0f}, "FT0Camplitude"};
    ConfigurableAxis axisFV0A{"FV0Aamplitude", {100, 0.0f, 2000.0f}, "FV0Aamplitude"};
    ConfigurableAxis axisFDDA{"FDDAamplitude", {100, 0.0f, 2000.0f}, "FDDAamplitude"};
    ConfigurableAxis axisFDDC{"FDDCamplitude", {100, 0.0f, 2000.0f}, "FDDCamplitude"};
    ConfigurableAxis axisZNA{"ZNAamplitude", {100, 0.0f, 250.0f}, "ZNAamplitude"};
    ConfigurableAxis axisZNC{"ZNCamplitude", {100, 0.0f, 250.0f}, "ZNCamplitude"};
  } axisDetectors;

  // For manual sliceBy
  Preslice<aod::V0Datas> V0perCollision = o2::aod::v0data::collisionId;
  Preslice<aod::CascDatas> CascperCollision = o2::aod::cascdata::collisionId;
  Preslice<aod::KFCascDatas> KFCascperCollision = o2::aod::cascdata::collisionId;
  Preslice<aod::TraCascDatas> TraCascperCollision = o2::aod::cascdata::collisionId;
  Preslice<aod::McParticles> mcParticlePerMcCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<UDCollisionsFull> udCollisionsPerCollision = o2::aod::udcollision::collisionId;

  std::vector<uint32_t> genK0Short;
  std::vector<uint32_t> genLambda;
  std::vector<uint32_t> genAntiLambda;
  std::vector<uint32_t> genXiMinus;
  std::vector<uint32_t> genXiPlus;
  std::vector<uint32_t> genOmegaMinus;
  std::vector<uint32_t> genOmegaPlus;

  float roundToPrecision(float number, float step = 0.01)
  {
    // this function rounds a certain number in an axis that is quantized by
    // the variable 'step'; the rounded number is placed halfway between
    // n*step and (n+1)*step such that analysis can be done with absolutely
    // no issue with precision 'step'.
    return step * static_cast<float>(static_cast<int>((number) / step)) + TMath::Sign(1.0f, number) * (0.5f) * step;
  }

  void init(InitContext&)
  {
    // setup map for fast checking if enabled
    static_for<0, nSpecies - 1>([&](auto i) {
      constexpr int index = i.value;
      int f = enableGeneratedInfo->get(particleNames[index].c_str(), "Enable");
      if (f == 1) {
        bitset(enabledBits, index);
      }
    });

    // Creation of histograms: MC generated
    for (Int_t i = 0; i < nSpecies; i++) {
      histos.add(Form("hGenerated%s", particleNames[i].data()), Form("hGenerated%s", particleNames[i].data()), kTH1D, {axisPt});
      histos.add(Form("h2dGenerated%s", particleNames[i].data()), Form("h2dGenerated%s", particleNames[i].data()), kTH2D, {axisCentrality, axisPt});
    }

    histos.add("h2dNVerticesVsCentrality", "h2dNVerticesVsCentrality", kTH2D, {axisCentrality, axisNVertices});

    if (doprocessV0FoundTags || doprocessCascFoundTags) {
      auto h = histos.add<TH1>("hFoundTagsCounters", "hFoundTagsCounters", kTH1D, {{6, -0.5f, 5.5f}});
      h->GetXaxis()->SetBinLabel(1, "Found V0s");
      h->GetXaxis()->SetBinLabel(2, "Findable V0s");
      h->GetXaxis()->SetBinLabel(3, "Findable & found V0s");
      h->GetXaxis()->SetBinLabel(4, "Found Cascades");
      h->GetXaxis()->SetBinLabel(5, "Findable Cascades");
      h->GetXaxis()->SetBinLabel(6, "Findable & found Cascades");
    }

    // for QA and test purposes
    auto hRawCentrality = histos.add<TH1>("hRawCentrality", "hRawCentrality", kTH1F, {axisRawCentrality});

    auto hFT0AMultVsFT0AUD = histos.add<TH2>("hFT0AMultVsFT0AUD", "hFT0AMultVsFT0AUD; FT0-A Mult; FT0-A UD", kTH2F, {axisDetectors.axisFT0A, axisDetectors.axisFT0A});
    auto hFT0CMultVsFT0CUD = histos.add<TH2>("hFT0CMultVsFT0CUD", "hFT0CMultVsFT0CUD; FT0-C Mult; FT0-C UD", kTH2F, {axisDetectors.axisFT0C, axisDetectors.axisFT0C});
    auto hFV0AMultVsFV0AUD = histos.add<TH2>("hFV0AMultVsFV0AUD", "hFV0AMultVsFV0AUD; FV0-A Mult; FV0-A UD", kTH2F, {axisDetectors.axisFV0A, axisDetectors.axisFV0A});
    auto hFDDAMultVsFDDAUD = histos.add<TH2>("hFDDAMultVsFDDAUD", "hFDDAMultVsFDDAUD; FDD-A Mult; FDD-A UD", kTH2F, {axisDetectors.axisFDDA, axisDetectors.axisFDDA});
    auto hFDDCMultVsFDDCUD = histos.add<TH2>("hFDDCMultVsFDDCUD", "hFDDCMultVsFDDCUD; FDD-C Mult; FDD-C UD", kTH2F, {axisDetectors.axisFDDC, axisDetectors.axisFDDC});
    auto hZNAMultVsZNAUD = histos.add<TH2>("hZNAMultVsZNAUD", "hZNAMultVsZNAUD; ZNA Mult; ZNA UD", kTH2F, {axisDetectors.axisZNA, axisDetectors.axisZNA});
    auto hZNCMultVsZNCUD = histos.add<TH2>("hZNCMultVsZNCUD", "hZNCMultVsZNCUD; ZNC Mult; ZNC UD", kTH2F, {axisDetectors.axisZNC, axisDetectors.axisZNC});

    for (int ii = 1; ii < 101; ii++) {
      float value = 100.5f - static_cast<float>(ii);
      hRawCentrality->SetBinContent(ii, value);
    }

    if (doprocessBinnedGenerated) {
      // reserve space for generated vectors if that process enabled
      auto hBinFinder = histos.get<TH2>(HIST("h2dGeneratedK0Short"));
      LOGF(info, "Binned generated processing enabled. Initialising with %i elements...", hBinFinder->GetNcells());
      genK0Short.resize(hBinFinder->GetNcells(), 0);
      genLambda.resize(hBinFinder->GetNcells(), 0);
      genAntiLambda.resize(hBinFinder->GetNcells(), 0);
      genXiMinus.resize(hBinFinder->GetNcells(), 0);
      genXiPlus.resize(hBinFinder->GetNcells(), 0);
      genOmegaMinus.resize(hBinFinder->GetNcells(), 0);
      genOmegaPlus.resize(hBinFinder->GetNcells(), 0);
      LOGF(info, "Binned generated processing: init done.");
    }
  }

  void processCollisionsV0sOnly(soa::Join<aod::Collisions, aod::FT0Mults, aod::FV0Mults, aod::FDDMults, aod::PVMults, aod::ZDCMults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::EvSels, aod::MultsExtra, aod::MultsGlobal> const& collisions, aod::V0Datas const& V0s, aod::BCsWithTimestamps const& /*bcs*/, UDCollisionsFull const& udCollisions)
  {
    for (const auto& collision : collisions) {
      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisColl = V0s.sliceBy(V0perCollision, collIdx);
      bool strange = V0Table_thisColl.size() > 0;

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();

      int gapSide = -1;
      float totalFT0AmplitudeA = -999;
      float totalFT0AmplitudeC = -999;
      float totalFV0AmplitudeA = -999;
      float totalFDDAmplitudeA = -999;
      float totalFDDAmplitudeC = -999;
      float energyCommonZNA = -999;
      float energyCommonZNC = -999;
      if (udCollisions.size() > 0) { // check that the UD collision table is not empty
        auto udCollision = udCollisions.sliceBy(udCollisionsPerCollision, collIdx);
        if (udCollision.size() == 1) { // check that the slicing provide a unique UD collision
          for (auto& udColl : udCollision) {
            gapSide = udColl.gapSide();
            totalFT0AmplitudeA = udColl.totalFT0AmplitudeA();
            totalFT0AmplitudeC = udColl.totalFT0AmplitudeC();
            totalFV0AmplitudeA = udColl.totalFV0AmplitudeA();
            totalFDDAmplitudeA = udColl.totalFDDAmplitudeA();
            totalFDDAmplitudeC = udColl.totalFDDAmplitudeC();
            energyCommonZNA = udColl.energyCommonZNA();
            energyCommonZNC = udColl.energyCommonZNC();
          }
        }
      }

      // casc table sliced
      if (strange || fillEmptyCollisions) {
        strangeColl(collision.posX(), collision.posY(), collision.posZ());
        strangeCents(collision.centFT0M(), collision.centFT0A(),
                     collision.centFT0C(), collision.centFV0A());
        strangeEvSels(collision.sel8(), collision.selection_raw(),
                      collision.multFT0A() * static_cast<float>(fillRawFT0A),
                      collision.multFT0C() * static_cast<float>(fillRawFT0C),
                      collision.multFV0A() * static_cast<float>(fillRawFV0A),
                      collision.multFDDA() * static_cast<float>(fillRawFDDA),
                      collision.multFDDC() * static_cast<float>(fillRawFDDC),
                      collision.multNTracksPVeta1() * static_cast<int>(fillRawNTracksEta1),
                      collision.multPVTotalContributors() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multNTracksGlobal() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multNTracksITSTPC() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multAllTracksTPCOnly() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multAllTracksITSTPC() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multZNA() * static_cast<float>(fillRawZDC),
                      collision.multZNC() * static_cast<float>(fillRawZDC),
                      collision.multZEM1() * static_cast<float>(fillRawZDC),
                      collision.multZEM2() * static_cast<float>(fillRawZDC),
                      collision.multZPA() * static_cast<float>(fillRawZDC),
                      collision.multZPC() * static_cast<float>(fillRawZDC),
                      collision.trackOccupancyInTimeRange(),
                      // UPC info
                      gapSide,
                      totalFT0AmplitudeA, totalFT0AmplitudeC, totalFV0AmplitudeA,
                      totalFDDAmplitudeA, totalFDDAmplitudeC,
                      energyCommonZNA, energyCommonZNC);
        strangeStamps(bc.runNumber(), bc.timestamp());
      }
      for (int i = 0; i < V0Table_thisColl.size(); i++)
        v0collref(strangeColl.lastIndex());
    }
  }

  void processCollisions(soa::Join<aod::Collisions, aod::FT0Mults, aod::FV0Mults, aod::FDDMults, aod::PVMults, aod::ZDCMults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::EvSels, aod::MultsExtra, aod::MultsGlobal> const& collisions, aod::V0Datas const& V0s, aod::CascDatas const& Cascades, aod::KFCascDatas const& KFCascades, aod::TraCascDatas const& TraCascades, aod::BCsWithTimestamps const& /*bcs*/, UDCollisionsFull const& udCollisions)
  {
    // create collision indices beforehand
    std::vector<int> V0CollIndices(V0s.size(), -1);                 // index -1: no collision
    std::vector<int> CascadeCollIndices(Cascades.size(), -1);       // index -1: no collision
    std::vector<int> KFCascadeCollIndices(KFCascades.size(), -1);   // index -1: no collision
    std::vector<int> TraCascadeCollIndices(TraCascades.size(), -1); // index -1: no collision

    for (const auto& collision : collisions) {
      const uint64_t collIdx = collision.globalIndex();

      float centrality = collision.centFT0C();
      if (qaCentrality) {
        auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
        centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(collision.multFT0C()));
      }

      auto V0Table_thisColl = V0s.sliceBy(V0perCollision, collIdx);
      auto CascTable_thisColl = Cascades.sliceBy(CascperCollision, collIdx);
      auto KFCascTable_thisColl = KFCascades.sliceBy(KFCascperCollision, collIdx);
      auto TraCascTable_thisColl = TraCascades.sliceBy(TraCascperCollision, collIdx);
      bool strange = V0Table_thisColl.size() > 0 ||
                     CascTable_thisColl.size() > 0 ||
                     KFCascTable_thisColl.size() > 0 ||
                     TraCascTable_thisColl.size() > 0;

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();

      int gapSide = -1;
      float totalFT0AmplitudeA = -999;
      float totalFT0AmplitudeC = -999;
      float totalFV0AmplitudeA = -999;
      float totalFDDAmplitudeA = -999;
      float totalFDDAmplitudeC = -999;
      float energyCommonZNA = -999;
      float energyCommonZNC = -999;
      if (udCollisions.size() > 0) { // check that the UD collision table is not empty
        auto udCollision = udCollisions.sliceBy(udCollisionsPerCollision, collIdx);
        if (udCollision.size() == 1) { // check that the slicing provide a unique UD collision
          for (auto& udColl : udCollision) {
            gapSide = udColl.gapSide();
            totalFT0AmplitudeA = udColl.totalFT0AmplitudeA();
            totalFT0AmplitudeC = udColl.totalFT0AmplitudeC();
            totalFV0AmplitudeA = udColl.totalFV0AmplitudeA();
            totalFDDAmplitudeA = udColl.totalFDDAmplitudeA();
            totalFDDAmplitudeC = udColl.totalFDDAmplitudeC();
            energyCommonZNA = udColl.energyCommonZNA();
            energyCommonZNC = udColl.energyCommonZNC();

            histos.fill(HIST("hFT0AMultVsFT0AUD"), collision.multFT0A(), udColl.totalFT0AmplitudeA());
            histos.fill(HIST("hFT0CMultVsFT0CUD"), collision.multFT0C(), udColl.totalFT0AmplitudeC());
            histos.fill(HIST("hFV0AMultVsFV0AUD"), collision.multFV0A(), udColl.totalFV0AmplitudeA());
            histos.fill(HIST("hFDDAMultVsFDDAUD"), collision.multFDDA(), udColl.totalFDDAmplitudeA());
            histos.fill(HIST("hFDDCMultVsFDDCUD"), collision.multFDDC(), udColl.totalFDDAmplitudeC());
            histos.fill(HIST("hZNAMultVsZNAUD"), collision.multZNA(), udColl.energyCommonZNA());
            histos.fill(HIST("hZNCMultVsZNCUD"), collision.multZNC(), udColl.energyCommonZNC());
          }
        }
      }

      // casc table sliced
      if (strange || fillEmptyCollisions) {
        strangeColl(collision.posX(), collision.posY(), collision.posZ());
        strangeCents(collision.centFT0M(), collision.centFT0A(),
                     centrality, collision.centFV0A());
        strangeEvSels(collision.sel8(), collision.selection_raw(),
                      collision.multFT0A() * static_cast<float>(fillRawFT0A),
                      collision.multFT0C() * static_cast<float>(fillRawFT0C),
                      collision.multFV0A() * static_cast<float>(fillRawFV0A),
                      collision.multFDDA() * static_cast<float>(fillRawFDDA),
                      collision.multFDDC() * static_cast<float>(fillRawFDDC),
                      collision.multNTracksPVeta1() * static_cast<int>(fillRawNTracksEta1),
                      collision.multPVTotalContributors() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multNTracksGlobal() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multNTracksITSTPC() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multAllTracksTPCOnly() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multAllTracksITSTPC() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multZNA() * static_cast<float>(fillRawZDC),
                      collision.multZNC() * static_cast<float>(fillRawZDC),
                      collision.multZEM1() * static_cast<float>(fillRawZDC),
                      collision.multZEM2() * static_cast<float>(fillRawZDC),
                      collision.multZPA() * static_cast<float>(fillRawZDC),
                      collision.multZPC() * static_cast<float>(fillRawZDC),
                      collision.trackOccupancyInTimeRange(),
                      // UPC info
                      gapSide,
                      totalFT0AmplitudeA, totalFT0AmplitudeC, totalFV0AmplitudeA,
                      totalFDDAmplitudeA, totalFDDAmplitudeC,
                      energyCommonZNA, energyCommonZNC);
        strangeStamps(bc.runNumber(), bc.timestamp());
      }

      for (const auto& v0 : V0Table_thisColl)
        V0CollIndices[v0.globalIndex()] = strangeColl.lastIndex();
      for (const auto& casc : CascTable_thisColl)
        CascadeCollIndices[casc.globalIndex()] = strangeColl.lastIndex();
      for (const auto& casc : KFCascTable_thisColl)
        KFCascadeCollIndices[casc.globalIndex()] = strangeColl.lastIndex();
      for (const auto& casc : TraCascTable_thisColl)
        TraCascadeCollIndices[casc.globalIndex()] = strangeColl.lastIndex();
    }

    // populate references, including those that might not be assigned
    for (const auto& v0 : V0s)
      v0collref(V0CollIndices[v0.globalIndex()]);
    for (const auto& casc : Cascades)
      casccollref(CascadeCollIndices[casc.globalIndex()]);
    for (const auto& casc : KFCascades)
      kfcasccollref(KFCascadeCollIndices[casc.globalIndex()]);
    for (const auto& casc : KFCascades)
      tracasccollref(TraCascadeCollIndices[casc.globalIndex()]);
  }

  void processCollisionsMC(soa::Join<aod::Collisions, aod::FT0Mults, aod::FV0Mults, aod::FDDMults, aod::PVMults, aod::ZDCMults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::EvSels, aod::McCollisionLabels, aod::MultsExtra, aod::MultsGlobal> const& collisions, soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const& /*V0MCCores*/, soa::Join<aod::CascDatas, aod::McCascLabels> const& Cascades, aod::KFCascDatas const& KFCascades, aod::TraCascDatas const& TraCascades, aod::BCsWithTimestamps const& /*bcs*/, UDCollisionsFull const& udCollisions, soa::Join<aod::McCollisions, aod::McCollsExtra, aod::MultsExtraMC> const& mcCollisions, aod::McParticles const&)
  {
    // create collision indices beforehand
    std::vector<int> V0CollIndices(V0s.size(), -1);                 // index -1: no collision
    std::vector<int> V0MCCollIndices(V0s.size(), -1);               // index -1: no collision
    std::vector<int> CascadeCollIndices(Cascades.size(), -1);       // index -1: no collision
    std::vector<int> CascadeMCCollIndices(Cascades.size(), -1);     // index -1: no collision
    std::vector<int> KFCascadeCollIndices(KFCascades.size(), -1);   // index -1: no collision
    std::vector<int> TraCascadeCollIndices(TraCascades.size(), -1); // index -1: no collision

    // ______________________________________________
    // fill all MC collisions, correlate via index later on
    for (const auto& mccollision : mcCollisions) {
      strangeMCColl(mccollision.posX(), mccollision.posY(), mccollision.posZ(), mccollision.impactParameter());
      strangeMCMults(mccollision.multMCFT0A(), mccollision.multMCFT0C(),
                     mccollision.multMCNParticlesEta05(),
                     mccollision.multMCNParticlesEta08(),
                     mccollision.multMCNParticlesEta10());
    }

    // ______________________________________________
    for (const auto& collision : collisions) {
      const uint64_t collIdx = collision.globalIndex();

      float centrality = collision.centFT0C();
      if (qaCentrality) {
        auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
        centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(collision.multFT0C()));
      }

      auto V0Table_thisColl = V0s.sliceBy(V0perCollision, collIdx);
      auto CascTable_thisColl = Cascades.sliceBy(CascperCollision, collIdx);
      auto KFCascTable_thisColl = KFCascades.sliceBy(KFCascperCollision, collIdx);
      auto TraCascTable_thisColl = TraCascades.sliceBy(TraCascperCollision, collIdx);
      bool strange = V0Table_thisColl.size() > 0 ||
                     CascTable_thisColl.size() > 0 ||
                     KFCascTable_thisColl.size() > 0 ||
                     TraCascTable_thisColl.size() > 0;

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();

      int gapSide = -1;
      float totalFT0AmplitudeA = -999;
      float totalFT0AmplitudeC = -999;
      float totalFV0AmplitudeA = -999;
      float totalFDDAmplitudeA = -999;
      float totalFDDAmplitudeC = -999;
      float energyCommonZNA = -999;
      float energyCommonZNC = -999;
      if (udCollisions.size() > 0) { // check that the UD collision table is not empty
        auto udCollision = udCollisions.sliceBy(udCollisionsPerCollision, collIdx);
        if (udCollision.size() == 1) { // check that the slicing provide a unique UD collision
          for (auto& udColl : udCollision) {
            gapSide = udColl.gapSide();
            totalFT0AmplitudeA = udColl.totalFT0AmplitudeA();
            totalFT0AmplitudeC = udColl.totalFT0AmplitudeC();
            totalFV0AmplitudeA = udColl.totalFV0AmplitudeA();
            totalFDDAmplitudeA = udColl.totalFDDAmplitudeA();
            totalFDDAmplitudeC = udColl.totalFDDAmplitudeC();
            energyCommonZNA = udColl.energyCommonZNA();
            energyCommonZNC = udColl.energyCommonZNC();

            histos.fill(HIST("hFT0AMultVsFT0AUD"), collision.multFT0A(), udColl.totalFT0AmplitudeA());
            histos.fill(HIST("hFT0CMultVsFT0CUD"), collision.multFT0C(), udColl.totalFT0AmplitudeC());
            histos.fill(HIST("hFV0AMultVsFV0AUD"), collision.multFV0A(), udColl.totalFV0AmplitudeA());
            histos.fill(HIST("hFDDAMultVsFDDAUD"), collision.multFDDA(), udColl.totalFDDAmplitudeA());
            histos.fill(HIST("hFDDCMultVsFDDCUD"), collision.multFDDC(), udColl.totalFDDAmplitudeC());
            histos.fill(HIST("hZNAMultVsZNAUD"), collision.multZNA(), udColl.energyCommonZNA());
            histos.fill(HIST("hZNCMultVsZNCUD"), collision.multZNC(), udColl.energyCommonZNC());
          }
        }
      }

      // casc table sliced
      if (strange || fillEmptyCollisions) {
        strangeColl(collision.posX(), collision.posY(), collision.posZ());
        strangeCollLabels(collision.mcCollisionId());
        strangeCents(collision.centFT0M(), collision.centFT0A(),
                     centrality, collision.centFV0A());
        strangeEvSels(collision.sel8(), collision.selection_raw(),
                      collision.multFT0A() * static_cast<float>(fillRawFT0A),
                      collision.multFT0C() * static_cast<float>(fillRawFT0C),
                      collision.multFV0A() * static_cast<float>(fillRawFV0A),
                      collision.multFDDA() * static_cast<float>(fillRawFDDA),
                      collision.multFDDC() * static_cast<float>(fillRawFDDC),
                      collision.multNTracksPVeta1() * static_cast<int>(fillRawNTracksEta1),
                      collision.multPVTotalContributors() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multNTracksGlobal() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multNTracksITSTPC() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multAllTracksTPCOnly() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multAllTracksITSTPC() * static_cast<int>(fillRawNTracksForCorrelation),
                      collision.multZNA() * static_cast<float>(fillRawZDC),
                      collision.multZNC() * static_cast<float>(fillRawZDC),
                      collision.multZEM1() * static_cast<float>(fillRawZDC),
                      collision.multZEM2() * static_cast<float>(fillRawZDC),
                      collision.multZPA() * static_cast<float>(fillRawZDC),
                      collision.multZPC() * static_cast<float>(fillRawZDC),
                      collision.trackOccupancyInTimeRange(),
                      // UPC info
                      gapSide,
                      totalFT0AmplitudeA, totalFT0AmplitudeC, totalFV0AmplitudeA,
                      totalFDDAmplitudeA, totalFDDAmplitudeC,
                      energyCommonZNA, energyCommonZNC);
        strangeStamps(bc.runNumber(), bc.timestamp());
      }
      for (const auto& v0 : V0Table_thisColl)
        V0CollIndices[v0.globalIndex()] = strangeColl.lastIndex();
      for (const auto& casc : CascTable_thisColl)
        CascadeCollIndices[casc.globalIndex()] = strangeColl.lastIndex();
      for (const auto& casc : KFCascTable_thisColl)
        KFCascadeCollIndices[casc.globalIndex()] = strangeColl.lastIndex();
      for (const auto& casc : TraCascTable_thisColl)
        TraCascadeCollIndices[casc.globalIndex()] = strangeColl.lastIndex();
    }

    // populate references, including those that might not be assigned
    for (const auto& v0 : V0s) {
      v0collref(V0CollIndices[v0.globalIndex()]);
    }
    for (const auto& casc : Cascades) {
      casccollref(CascadeCollIndices[casc.globalIndex()]);
    }
    for (const auto& casc : KFCascades)
      kfcasccollref(KFCascadeCollIndices[casc.globalIndex()]);
    for (const auto& casc : KFCascades)
      tracasccollref(TraCascadeCollIndices[casc.globalIndex()]);
  }

  void processTrackExtrasV0sOnly(aod::V0Datas const& V0s, TracksWithExtra const& tracksExtra)
  {
    std::vector<int> trackMap(tracksExtra.size(), -1); // index -1: not used

    //__________________________________________________
    // mark tracks that belong to V0s
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.posTrack_as<TracksWithExtra>();
      auto const& negTrack = v0.negTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // Figure out the numbering of the new tracks table
    // assume filling per order
    int nTracks = 0;
    for (int i = 0; i < static_cast<int>(trackMap.size()); i++) {
      if (trackMap[i] >= 0) {
        trackMap[i] = nTracks++;
      }
    }
    //__________________________________________________
    // populate track references
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.posTrack_as<TracksWithExtra>();
      auto const& negTrack = v0.negTrack_as<TracksWithExtra>();
      v0Extras(trackMap[posTrack.globalIndex()],
               trackMap[negTrack.globalIndex()]); // joinable with V0Datas
    }
    //__________________________________________________
    // circle back and populate actual DauTrackExtra table
    for (auto const& tr : tracksExtra) {
      if (trackMap[tr.globalIndex()] >= 0) {
        dauTrackExtras(tr.itsChi2NCl(),
                       tr.detectorMap(), tr.itsClusterSizes(),
                       tr.tpcNClsFound(), tr.tpcNClsCrossedRows());
      }
    }
    // done!
  }

  template <typename V0Datas, typename CascDatas, typename KFCascDatas, typename TraCascDatas, typename tracksWithExtra>
  void fillTrackExtras(V0Datas const& V0s, CascDatas const& Cascades, KFCascDatas const& KFCascades, TraCascDatas const& TraCascades, tracksWithExtra const& tracksExtra)
  {
    std::vector<int> trackMap(tracksExtra.size(), -1); // index -1: not used

    //__________________________________________________
    // mark tracks that belong to V0s
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.template posTrack_as<tracksWithExtra>();
      auto const& negTrack = v0.template negTrack_as<tracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
    }

    //__________________________________________________
    // index tracks that belong to CascDatas
    for (auto const& casc : Cascades) {
      auto bachTrack = casc.template bachelor_as<tracksWithExtra>();
      auto posTrack = casc.template posTrack_as<tracksWithExtra>();
      auto negTrack = casc.template negTrack_as<tracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // index tracks that belong to KFCascDatas
    for (auto const& casc : KFCascades) {
      auto bachTrack = casc.template bachelor_as<tracksWithExtra>();
      auto posTrack = casc.template posTrack_as<tracksWithExtra>();
      auto negTrack = casc.template negTrack_as<tracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // index tracks that belong to TraCascDatas
    for (auto const& casc : TraCascades) {
      auto bachTrack = casc.template bachelor_as<tracksWithExtra>();
      auto posTrack = casc.template posTrack_as<tracksWithExtra>();
      auto negTrack = casc.template negTrack_as<tracksWithExtra>();
      auto strangeTrack = casc.template strangeTrack_as<tracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
      trackMap[strangeTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // Figure out the numbering of the new tracks table
    // assume filling per order
    int nTracks = 0;
    for (int i = 0; i < static_cast<int>(trackMap.size()); i++) {
      if (trackMap[i] >= 0) {
        trackMap[i] = nTracks++;
      }
    }
    //__________________________________________________
    // populate track references
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.template posTrack_as<tracksWithExtra>();
      auto const& negTrack = v0.template negTrack_as<tracksWithExtra>();
      v0Extras(trackMap[posTrack.globalIndex()],
               trackMap[negTrack.globalIndex()]); // joinable with V0Datas
    }
    //__________________________________________________
    // populate track references
    for (auto const& casc : Cascades) {
      auto bachTrack = casc.template bachelor_as<tracksWithExtra>();
      auto posTrack = casc.template posTrack_as<tracksWithExtra>();
      auto negTrack = casc.template negTrack_as<tracksWithExtra>();
      cascExtras(trackMap[posTrack.globalIndex()],
                 trackMap[negTrack.globalIndex()],
                 trackMap[bachTrack.globalIndex()]); // joinable with CascDatas
    }
    //__________________________________________________
    // populate track references
    for (auto const& casc : TraCascades) {
      auto strangeTrack = casc.template strangeTrack_as<tracksWithExtra>();
      straTrackExtras(trackMap[strangeTrack.globalIndex()]); // joinable with TraCascDatas
    }
    //__________________________________________________
    // circle back and populate actual DauTrackExtra table
    for (auto const& tr : tracksExtra) {
      if (trackMap[tr.globalIndex()] >= 0) {
        dauTrackExtras(tr.itsChi2NCl(),
                       tr.detectorMap(), tr.itsClusterSizes(),
                       tr.tpcNClsFound(), tr.tpcNClsCrossedRows());

        // if the table has MC info
        if constexpr (requires { tr.mcParticle(); }) {
          // do your thing with the mcParticleIds only in case the table has the MC info
          dauTrackMCIds(tr.mcParticleId()); // joinable with dauTrackExtras
        }

        // round if requested
        if (roundNSigmaVariables) {
          dauTrackTPCPIDs(tr.tpcSignal(),
                          roundToPrecision(tr.tpcNSigmaEl(), precisionNSigmas),
                          roundToPrecision(tr.tpcNSigmaPi(), precisionNSigmas),
                          roundToPrecision(tr.tpcNSigmaKa(), precisionNSigmas),
                          roundToPrecision(tr.tpcNSigmaPr(), precisionNSigmas),
                          roundToPrecision(tr.tpcNSigmaHe(), precisionNSigmas));
        } else {
          dauTrackTPCPIDs(tr.tpcSignal(), tr.tpcNSigmaEl(),
                          tr.tpcNSigmaPi(), tr.tpcNSigmaKa(),
                          tr.tpcNSigmaPr(), tr.tpcNSigmaHe());
        }

        // populate daughter-level TOF information
        dauTrackTOFPIDs(tr.tofSignal(), tr.tofEvTime(), tr.length());
      }
    }
    // done!
  }

  void processTrackExtras(aod::V0Datas const& V0s, aod::CascDatas const& Cascades, aod::KFCascDatas const& KFCascades, aod::TraCascDatas const& TraCascades, soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullHe, aod::TOFEvTime, aod::TOFSignal> const& tracksExtra, aod::V0s const&)
  {
    fillTrackExtras(V0s, Cascades, KFCascades, TraCascades, tracksExtra);
    // done!
  }

  void processTrackExtrasMC(aod::V0Datas const& V0s, aod::CascDatas const& Cascades, aod::KFCascDatas const& KFCascades, aod::TraCascDatas const& TraCascades, soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullHe, aod::TOFEvTime, aod::TOFSignal> const& tracksExtra, aod::V0s const&)
  {
    fillTrackExtras(V0s, Cascades, KFCascades, TraCascades, tracksExtra);
    // done!
  }

  void processStrangeMothers(soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s, soa::Join<aod::CascDatas, aod::McCascLabels> const& Cascades, aod::McParticles const& mcParticles)
  {
    std::vector<int> motherReference(mcParticles.size(), -1); // index -1: not used / no reference

    //__________________________________________________
    // mark mcParticles for referencing
    for (auto const& v0 : V0s)
      if (v0.has_mcMotherParticle())
        motherReference[v0.mcMotherParticleId()] = 0;
    for (auto const& ca : Cascades)
      if (ca.has_mcMotherParticle())
        motherReference[ca.mcMotherParticleId()] = 0;
    //__________________________________________________
    // Figure out the numbering of the new mcMother table
    // assume filling per order
    int nParticles = 0;
    for (int i = 0; i < static_cast<int>(motherReference.size()); i++) {
      if (motherReference[i] >= 0) {
        motherReference[i] = nParticles++; // count particles of interest
      }
    }
    //__________________________________________________
    // populate track references
    for (auto const& v0 : V0s)
      v0mothers(motherReference[v0.mcMotherParticleId()]); // joinable with V0Datas
    for (auto const& ca : Cascades)
      cascmothers(motherReference[ca.mcMotherParticleId()]); // joinable with CascDatas
    //__________________________________________________
    // populate motherMCParticles
    for (auto const& tr : mcParticles) {
      if (motherReference[tr.globalIndex()] >= 0) {
        motherMCParts(tr.px(), tr.py(), tr.pz(), tr.pdgCode(), tr.isPhysicalPrimary());
      }
    }
  }

  using interlinkedCascades = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink, aod::TraCascDataLink>;

  void processCascadeInterlinkTracked(interlinkedCascades const& /*masterCascades*/, aod::CascIndices const& Cascades, aod::TraCascIndices const& TraCascades)
  {
    // Standard to tracked
    for (auto const& c : Cascades) {
      int indexTracked = -1;
      if (c.has_cascade()) {
        auto cascade = c.cascade_as<interlinkedCascades>();
        indexTracked = cascade.traCascDataId();
      }
      cascToTraRefs(indexTracked);
    }
    // Tracked to standard
    for (auto const& c : TraCascades) {
      int index = -1;
      if (c.has_cascade()) {
        auto cascade = c.cascade_as<interlinkedCascades>();
        index = cascade.cascDataId();
      }
      traToCascRefs(index);
    }
  }

  void processCascadeInterlinkKF(interlinkedCascades const& /*masterCascades*/, aod::CascIndices const& Cascades, aod::KFCascIndices const& KFCascades)
  {
    // Standard to KF
    for (auto const& c : Cascades) {
      int indexKF = -1;
      if (c.has_cascade()) {
        auto cascade = c.cascade_as<interlinkedCascades>();
        indexKF = cascade.kfCascDataId();
      }
      cascToKFRefs(indexKF);
    }
    // KF to standard
    for (auto const& c : KFCascades) {
      int index = -1;
      if (c.has_cascade()) {
        auto cascade = c.cascade_as<interlinkedCascades>();
        index = cascade.cascDataId();
      }
      kfToCascRefs(index);
    }
  }

  void processPureSimulation(aod::McParticles const& mcParticles)
  {
    for (auto& mcp : mcParticles) {
      if (TMath::Abs(mcp.y()) < 0.5) {
        static_for<0, nSpecies - 1>([&](auto i) {
          constexpr int index = i.value;
          if (mcp.pdgCode() == particlePDGCodes[index] && bitcheck(enabledBits, index)) {
            histos.fill(HIST("hGenerated") + HIST(particleNamesConstExpr[index]), mcp.pt());
          }
        });
      }
    }
  }

  void processReconstructedSimulation(aod::McCollision const& /*mcCollision*/, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::FT0Mults>> const& collisions, aod::McParticles const& mcParticles)
  {
    // this process function also checks if a given collision was reconstructed and checks explicitly for splitting, etc
    // identify best-of collision
    int biggestNContribs = -1;
    float bestCentrality = 100.5;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCentrality = collision.centFT0C();
        if (qaCentrality) {
          auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
          bestCentrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(collision.multFT0C()));
        }
      }
    }
    histos.fill(HIST("h2dNVerticesVsCentrality"), bestCentrality, collisions.size());

    for (auto& mcp : mcParticles) {
      if (TMath::Abs(mcp.y()) < 0.5 && mcp.isPhysicalPrimary()) {
        static_for<0, nSpecies - 1>([&](auto i) {
          constexpr int index = i.value;
          if (mcp.pdgCode() == particlePDGCodes[index] && bitcheck(enabledBits, index)) {
            histos.fill(HIST("h2dGenerated") + HIST(particleNamesConstExpr[index]), bestCentrality, mcp.pt());
          }
        });
      }
    }
  }

  void processBinnedGenerated(soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions, aod::McParticles const& mcParticlesEntireTable)
  {
    // set to zero
    std::fill(genK0Short.begin(), genK0Short.end(), 0);
    std::fill(genLambda.begin(), genLambda.end(), 0);
    std::fill(genAntiLambda.begin(), genAntiLambda.end(), 0);
    std::fill(genXiMinus.begin(), genXiMinus.end(), 0);
    std::fill(genXiPlus.begin(), genXiPlus.end(), 0);
    std::fill(genOmegaMinus.begin(), genOmegaMinus.end(), 0);
    std::fill(genOmegaPlus.begin(), genOmegaPlus.end(), 0);

    // this process function also checks if a given collision was reconstructed and checks explicitly for splitting, etc
    for (auto& mcCollision : mcCollisions) {
      const uint64_t mcCollIndex = mcCollision.globalIndex();

      // use one of the generated histograms as the bin finder
      auto hBinFinder = histos.get<TH2>(HIST("h2dGeneratedK0Short"));

      auto mcParticles = mcParticlesEntireTable.sliceBy(mcParticlePerMcCollision, mcCollIndex);
      for (auto& mcp : mcParticles) {
        if (TMath::Abs(mcp.y()) < 0.5 && mcp.isPhysicalPrimary()) {
          auto binNumber = hBinFinder->FindBin(mcCollision.bestCollisionCentFT0C(), mcp.pt()); // caution: pack
          if (mcp.pdgCode() == 310)
            genK0Short[binNumber]++;
          if (mcp.pdgCode() == 3122)
            genLambda[binNumber]++;
          if (mcp.pdgCode() == -3122)
            genAntiLambda[binNumber]++;
          if (mcp.pdgCode() == 3312)
            genXiMinus[binNumber]++;
          if (mcp.pdgCode() == -3312)
            genXiPlus[binNumber]++;
          if (mcp.pdgCode() == 3334)
            genOmegaMinus[binNumber]++;
          if (mcp.pdgCode() == -3334)
            genOmegaPlus[binNumber]++;
        }
      }
    }
    // at end of data frame
    // -> pack information from this DF into a generated histogram, once / DF
    geK0Short(genK0Short);
    geLambda(genLambda);
    geAntiLambda(genAntiLambda);
    geXiMinus(genXiMinus);
    geXiPlus(genXiPlus);
    geOmegaMinus(genOmegaMinus);
    geOmegaPlus(genOmegaPlus);
  }

  void processFT0AQVectors(soa::Join<aod::Collisions, aod::QvectorFT0As>::iterator const& collision)
  {
    StraFT0AQVs(collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.sumAmplFT0A());
  }
  void processFT0CQVectors(soa::Join<aod::Collisions, aod::QvectorFT0Cs>::iterator const& collision)
  {
    StraFT0CQVs(collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.sumAmplFT0C());
  }
  void processFT0CQVectorsLF(soa::Join<aod::Collisions, aod::EPCalibrationTables>::iterator const& collision)
  {
    StraFT0CQVs(collision.qFT0C() * std::cos(2 * collision.psiFT0C()), collision.qFT0C() * std::sin(2 * collision.psiFT0C()), collision.qFT0C());
    StraFT0CQVsEv(collision.triggereventep());
  }
  void processFT0MQVectors(soa::Join<aod::Collisions, aod::QvectorFT0Ms>::iterator const& collision)
  {
    StraFT0MQVs(collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.sumAmplFT0M());
  }
  void processFV0AQVectors(soa::Join<aod::Collisions, aod::QvectorFV0As>::iterator const& collision)
  {
    StraFV0AQVs(collision.qvecFV0ARe(), collision.qvecFV0AIm(), collision.sumAmplFV0A());
  }
  void processTPCQVectors(soa::Join<aod::Collisions, aod::QvectorBPoss, aod::QvectorBNegs>::iterator const& collision)
  {
    StraTPCQVs(collision.qvecBNegRe(), collision.qvecBNegIm(), collision.nTrkBNeg(), collision.qvecBPosRe(), collision.qvecBPosIm(), collision.nTrkBPos());
  }
  void processTPCQVectorsLF(soa::Join<aod::Collisions, aod::EPCalibrationTables>::iterator const& collision)
  {
    StraTPCQVs(collision.qTPCL() * std::cos(2 * collision.psiTPCL()), collision.qTPCL() * std::sin(2 * collision.psiTPCL()), collision.qTPCL(), collision.qTPCR() * std::cos(2 * collision.psiTPCR()), collision.qTPCR() * std::sin(2 * collision.psiTPCR()), collision.qTPCR());
  }

  uint64_t combineProngIndices(uint32_t low, uint32_t high)
  {
    return ((static_cast<uint64_t>(high)) << 32) | (static_cast<uint64_t>(low));
  }

  void processV0FoundTags(aod::V0s const& foundV0s, aod::V0Datas const& findableV0s, aod::FindableV0s const& /* added to avoid troubles */)
  {
    histos.fill(HIST("hFoundTagsCounters"), 0.0f, foundV0s.size());
    histos.fill(HIST("hFoundTagsCounters"), 1.0f, findableV0s.size());

    for (auto const& findableV0 : findableV0s) {
      bool hasBeenFound = false;
      for (auto const& foundV0 : foundV0s) {
        if (foundV0.posTrackId() == findableV0.posTrackId() && foundV0.negTrackId() == findableV0.negTrackId()) {
          hasBeenFound = true;
        }
      }
      v0FoundTags(hasBeenFound);
    }
  }

  using uint128_t = __uint128_t;
  uint128_t combineProngIndices128(uint32_t pos, uint32_t neg, uint32_t bach)
  {
    return ((static_cast<uint128_t>(pos)) << 64) | ((static_cast<uint128_t>(neg)) << 32) | (static_cast<uint128_t>(bach));
  }

  void processCascFoundTags(aod::Cascades const& foundCascades, aod::CascDatas const& findableCascades, aod::V0s const&, aod::FindableCascades const& /* added to avoid troubles */)
  {
    histos.fill(HIST("hFoundTagsCounters"), 3.0f, foundCascades.size());
    histos.fill(HIST("hFoundTagsCounters"), 4.0f, findableCascades.size());

    // pack the found V0s in a long long
    std::vector<uint128_t> foundCascadesPacked;
    foundCascadesPacked.reserve(foundCascades.size());
    for (auto const& foundCascade : foundCascades) {
      auto v0 = foundCascade.v0();
      foundCascadesPacked[foundCascade.globalIndex()] = combineProngIndices128(v0.posTrackId(), v0.negTrackId(), foundCascade.bachelorId());
    }

    bool hasBeenFound = false;
    for (auto const& findableCascade : findableCascades) {
      uint128_t indexPack = combineProngIndices128(findableCascade.posTrackId(), findableCascade.negTrackId(), findableCascade.bachelorId());
      for (uint32_t ic = 0; ic < foundCascades.size(); ic++) {
        if (indexPack == foundCascadesPacked[ic]) {
          hasBeenFound = true;
          histos.fill(HIST("hFoundTagsCounters"), 5.0f);
          break;
        }
      }
      cascFoundTags(hasBeenFound);
    }
  }

  PROCESS_SWITCH(strangederivedbuilder, processCollisionsV0sOnly, "Produce collisions (V0s only)", true);
  PROCESS_SWITCH(strangederivedbuilder, processCollisions, "Produce collisions (V0s + casc)", true);
  PROCESS_SWITCH(strangederivedbuilder, processCollisionsMC, "Produce collisions (V0s + casc)", false);
  PROCESS_SWITCH(strangederivedbuilder, processTrackExtrasV0sOnly, "Produce track extra information (V0s only)", true);
  PROCESS_SWITCH(strangederivedbuilder, processTrackExtras, "Produce track extra information (V0s + casc)", true);
  PROCESS_SWITCH(strangederivedbuilder, processTrackExtrasMC, "Produce track extra information (V0s + casc)", false);
  PROCESS_SWITCH(strangederivedbuilder, processStrangeMothers, "Produce tables with mother info for V0s + casc", true);
  PROCESS_SWITCH(strangederivedbuilder, processCascadeInterlinkTracked, "Produce tables interconnecting cascades", false);
  PROCESS_SWITCH(strangederivedbuilder, processCascadeInterlinkKF, "Produce tables interconnecting cascades", false);
  PROCESS_SWITCH(strangederivedbuilder, processPureSimulation, "Produce pure simulated information", true);
  PROCESS_SWITCH(strangederivedbuilder, processReconstructedSimulation, "Produce reco-ed simulated information", true);
  PROCESS_SWITCH(strangederivedbuilder, processBinnedGenerated, "Produce binned generated information", false);

  // event plane information
  PROCESS_SWITCH(strangederivedbuilder, processFT0AQVectors, "Produce FT0A Q-vectors table", false);
  PROCESS_SWITCH(strangederivedbuilder, processFT0CQVectors, "Produce FT0C Q-vectors table", false);
  PROCESS_SWITCH(strangederivedbuilder, processFT0CQVectorsLF, "Produce FT0C Q-vectors table using LF temporary calibration", false);
  PROCESS_SWITCH(strangederivedbuilder, processFT0MQVectors, "Produce FT0M Q-vectors table", false);
  PROCESS_SWITCH(strangederivedbuilder, processFV0AQVectors, "Produce FV0A Q-vectors table", false);
  PROCESS_SWITCH(strangederivedbuilder, processTPCQVectors, "Produce TPC Q-vectors table", false);
  PROCESS_SWITCH(strangederivedbuilder, processTPCQVectorsLF, "Produce TPC Q-vectors table using LF temporary calibration", false);

  // dedicated findable functionality
  PROCESS_SWITCH(strangederivedbuilder, processV0FoundTags, "Produce FoundV0Tags for findable exercise", false);
  PROCESS_SWITCH(strangederivedbuilder, processCascFoundTags, "Produce FoundCascTags for findable exercise", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangederivedbuilder>(cfgc)};
}

//__________________________________________________
// do not over-populate general namespace, keep scope strangederivedbuilder::
const std::vector<std::string> strangederivedbuilder::particleNames{"Gamma", "K0Short", "Lambda", "AntiLambda",
                                                                    "Sigma0", "AntiSigma0", "SigmaPlus", "SigmaMinus",
                                                                    "Hypertriton", "AntiHypertriton",
                                                                    "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
const std::vector<int> strangederivedbuilder::particlePDGCodes{22, 310, 3122, -3122, 3212, -3212, 3222, 3112,
                                                               1010010030, -1010010030, 3312, -3312, 3334, -3334};
const std::vector<std::string> strangederivedbuilder::parameterNames{"Enable"};

const int strangederivedbuilder::defaultParameters[strangederivedbuilder::nSpecies][strangederivedbuilder::nParameters] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
