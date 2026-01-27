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
// this task provides produces histograms containing
// the number of generated quarkonia per unit of
// percentile and per unit of pT
// It is meant to help with providing auxiliary information
// when dealing with derived data.

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

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
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <map>
#include <utility>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// simple bit checkers
#define bitset(var, nbit) ((var) |= (1 << (nbit)))
#define bitcheck(var, nbit) ((var) & (1 << (nbit)))

struct generatedQuarkoniaMC {
  SliceCache cache;
  //__________________________________________________
  // Generated binned data
  // this is a hack while the system does not do better
  Produces<aod::GeEtaC1S> geEtaC1S;
  Produces<aod::GeJPsi> geJPsi;
  Produces<aod::GeChiC0> geChiC0;
  Produces<aod::GeChiC1> geChiC1;
  Produces<aod::GeHC> geHC;
  Produces<aod::GeChiC2> geChiC2;
  Produces<aod::GeEtaC2S> geEtaC2S;
  Produces<aod::GePsi2S> gePsi2S;

  // histogram registry for bookkeeping
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr int nSpecies = 8;
  static constexpr int nParameters = 1;
  static const std::vector<std::string> particleNames;
  static const std::vector<int> particlePDGCodes;
  static const std::vector<std::string> parameterNames;
  static const int defaultParameters[nSpecies][nParameters];
  static constexpr std::string_view particleNamesConstExpr[] = {"EtaC1S", "JPsi", "ChiC0", "ChiC1",
                                                                "hC", "ChiC2", "EtaC2S", "Psi2S"};

  uint32_t enabledBits = 0;

  Configurable<LabeledArray<int>> enableGeneratedInfo{"enableGeneratedInfo",
                                                      {defaultParameters[0], nSpecies,
                                                       nParameters, particleNames, parameterNames},
                                                      "Fill generated particle histograms for each species. 0: no, 1: yes"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisCentrality{"axisCentrality", {100, 0.0f, 100.0f}, "Centrality"};

  ConfigurableAxis axisNVertices{"axisNVertices", {10, -0.5f, 9.5f}, "N(vertices)"};

  std::vector<uint32_t> genEtaC1S;
  std::vector<uint32_t> genJPsi;
  std::vector<uint32_t> genChiC0;
  std::vector<uint32_t> genChiC1;
  std::vector<uint32_t> genHC;
  std::vector<uint32_t> genChiC2;
  std::vector<uint32_t> genEtaC2S;
  std::vector<uint32_t> genPsi2S;

  // Preslice
  Preslice<aod::McParticles> mcParticlePerMcCollision = o2::aod::mcparticle::mcCollisionId;

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
      histos.add(Form("h2dGenerated%s", particleNames[i].data()), Form("h2dGenerated%s", particleNames[i].data()), kTH2D, {axisCentrality, axisPt});
    }

    histos.add("h2dNVerticesVsCentrality", "h2dNVerticesVsCentrality", kTH2D, {axisCentrality, axisNVertices});

    // reserve space for generated vectors if that process enabled
    auto hBinFinder = histos.get<TH2>(HIST("h2dGeneratedEtaC1S"));
    LOGF(info, "Binned generated processing enabled. Initialising with %i elements...", hBinFinder->GetNcells());
    genEtaC1S.resize(hBinFinder->GetNcells(), 0);
    genJPsi.resize(hBinFinder->GetNcells(), 0);
    genChiC0.resize(hBinFinder->GetNcells(), 0);
    genChiC1.resize(hBinFinder->GetNcells(), 0);
    genHC.resize(hBinFinder->GetNcells(), 0);
    genChiC2.resize(hBinFinder->GetNcells(), 0);
    genEtaC2S.resize(hBinFinder->GetNcells(), 0);
    genPsi2S.resize(hBinFinder->GetNcells(), 0);
    LOGF(info, "Binned generated processing: init done.");
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
        bestCentrality = collision.centFT0M();
      }
    }
    histos.fill(HIST("h2dNVerticesVsCentrality"), bestCentrality, collisions.size());

    for (auto& mcp : mcParticles) {
      if (TMath::Abs(mcp.y()) < 0.5 /* && mcp.isPhysicalPrimary()*/) {
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
    std::fill(genEtaC1S.begin(), genEtaC1S.end(), 0);
    std::fill(genJPsi.begin(), genJPsi.end(), 0);
    std::fill(genChiC0.begin(), genChiC0.end(), 0);
    std::fill(genChiC1.begin(), genChiC1.end(), 0);
    std::fill(genHC.begin(), genHC.end(), 0);
    std::fill(genChiC2.begin(), genChiC2.end(), 0);
    std::fill(genEtaC2S.begin(), genEtaC2S.end(), 0);
    std::fill(genPsi2S.begin(), genPsi2S.end(), 0);

    // this process function also checks if a given collision was reconstructed and checks explicitly for splitting, etc
    for (auto& mcCollision : mcCollisions) {
      const uint64_t mcCollIndex = mcCollision.globalIndex();

      // use one of the generated histograms as the bin finder
      auto hBinFinder = histos.get<TH2>(HIST("h2dGeneratedEtaC1S"));

      auto mcParticles = mcParticlesEntireTable.sliceBy(mcParticlePerMcCollision, mcCollIndex);
      for (auto& mcp : mcParticles) {
        if (TMath::Abs(mcp.y()) < 0.5 /* && mcp.isPhysicalPrimary()*/) {
          auto binNumber = hBinFinder->FindBin(mcCollision.bestCollisionCentFT0C(), mcp.pt()); // caution: pack
          if (mcp.pdgCode() == 441)
            genEtaC1S[binNumber]++;
          if (mcp.pdgCode() == 443) {
            genJPsi[binNumber]++;
            auto daughters = mcp.daughters_as<aod::McParticles>();
            std::cout << mcp.pdgCode() ;
            for(auto dau : daughters) {
              std::cout << " dau: " << dau.pdgCode() << " Primary? " << dau.isPhysicalPrimary();
              auto subdaughters = dau.daughters_as<aod::McParticles>();
              for(auto subdau : subdaughters)
                std::cout << " subdau: " << subdau.pdgCode();
            }
            std::cout << std::endl;
          }
          if (mcp.pdgCode() == 10441)
            genChiC0[binNumber]++;
          if (mcp.pdgCode() == 20443)
            genChiC1[binNumber]++;
          if (mcp.pdgCode() == 10443)
            genHC[binNumber]++;
          if (mcp.pdgCode() == 445)
            genChiC2[binNumber]++;
          if (mcp.pdgCode() == 100441)
            genEtaC2S[binNumber]++;
          if (mcp.pdgCode() == 100443)
            genPsi2S[binNumber]++;
        }
      }
    }
    // at end of data frame
    // -> pack information from this DF into a generated histogram, once / DF
    geEtaC1S(genEtaC1S);
    geJPsi(genJPsi);
    geChiC0(genChiC0);
    geChiC1(genChiC1);
    geHC(genHC);
    geChiC2(genChiC2);
    geEtaC2S(genEtaC2S);
    gePsi2S(genPsi2S);
  }

  PROCESS_SWITCH(generatedQuarkoniaMC, processReconstructedSimulation, "Produce reco-ed simulated information", true);
  PROCESS_SWITCH(generatedQuarkoniaMC, processBinnedGenerated, "Produce binned generated information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<generatedQuarkoniaMC>(cfgc)};
}

//__________________________________________________
// do not over-populate general namespace, keep scope generatedQuarkoniaMC::
const std::vector<std::string> generatedQuarkoniaMC::particleNames{"EtaC1S", "JPsi", "ChiC0", "ChiC1",
                                                                   "hC", "ChiC2", "EtaC2S", "Psi2S"};
const std::vector<int> generatedQuarkoniaMC::particlePDGCodes{441, 443, 10441, 20443, 10443, 445, 100441, 100443};
const std::vector<std::string> generatedQuarkoniaMC::parameterNames{"Enable"};

const int generatedQuarkoniaMC::defaultParameters[generatedQuarkoniaMC::nSpecies][generatedQuarkoniaMC::nParameters] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
