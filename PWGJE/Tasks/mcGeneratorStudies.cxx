// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Task that produces the generated pT spectrum of a given particle for MC studies
//
/// \author Nicolas Strangmann <nicolas.strangmann@cern.ch>, Goethe University Frankfurt / Oak Ridge National Laoratory

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "DetectorsBase/GeometryManager.h"
#include "EMCALBase/Geometry.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyMCCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::EMCALMatchedCollisions>;
using MyBCs = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;

struct MCGeneratorStudies {
  HistogramRegistry mHistManager{"MCGeneratorStudyHistograms"};

  Configurable<float> mVertexCut{"vertexCut", 10.f, "apply z-vertex cut with value in cm"};
  Configurable<float> mRapidityCut{"rapidityCut", 0.9f, "Maximum absolute rapidity of counted generated particles"};
  Configurable<int> mSelectedParticleCode{"particlePDGCode", 111, "PDG code of the particle to be investigated (0 for all)"};
  Configurable<bool> mRequireGammaGammaDecay{"requireGammaGammaDecay", false, "Only count generated particles that decayed into two photons"};
  Configurable<bool> mRequireEMCCellContent{"requireEMCCellContent", false, "Ask forEMCal cell content instead of the kTVXinEMC trigger"};

  Configurable<bool> mRequireTVX{"mRequireTVX", true, "require FT0AND in event cut"};
  Configurable<bool> mRequireSel8{"mRequireSel8", true, "require sel8 in event cut"};
  Configurable<bool> mRequireNoSameBunchPileup{"mRequireNoSameBunchPileup", true, "require no same bunch pileup in event cut"};
  Configurable<bool> mRequireGoodZvtxFT0vsPV{"mRequireGoodZvtxFT0vsPV", true, "require good Zvtx between FT0 vs. PV in event cut"};
  Configurable<bool> mRequireEMCReadoutInMB{"mRequireEMCReadoutInMB", true, "require the EMC to be read out in an MB collision (kTVXinEMC)"};

  void init(InitContext const&)
  {
    AxisSpec pTAxis{250, 0., 25., "#it{p}_{T} (GeV/#it{c})"};

    auto hCollisionCounter = mHistManager.add<TH1>("hCollisionCounter", "Number of collisions after event cuts", HistType::kTH1F, {{7, 0.5, 7.5}});
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "+TVX");         // TVX
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "+|z|<10cm");    // TVX with z < 10cm
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "+Sel8");        // TVX with z < 10cm and Sel8
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "+Good z vtx");  // TVX with z < 10cm and Sel8 and good z xertex
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "+unique");      // TVX with z < 10cm and Sel8 and good z xertex and unique (only collision in the BC)
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "+EMC readout"); // TVX with z < 10cm and Sel8 and good z xertex and unique (only collision in the BC) and kTVXinEMC

    auto hBCCounter = mHistManager.add<TH1>("hBCCounter", "Number of BCs after BC cuts", HistType::kTH1F, {{3, 0.5, 3.5}});
    hBCCounter->GetXaxis()->SetBinLabel(1, "all");
    hBCCounter->GetXaxis()->SetBinLabel(2, "+TVX");
    hBCCounter->GetXaxis()->SetBinLabel(3, "+Collision");

    TString mesonLatexString = (TString)mSelectedParticleCode;
    switch (mSelectedParticleCode) {
      case 0:
        mesonLatexString = "particles";
        break;
      case 111:
        mesonLatexString = "#pi^{0}";
        break;
      case 221:
        mesonLatexString = "#eta";
        break;
    }
    mHistManager.add("Yield", Form("Generated %s in all collisions", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_Accepted", Form("Accepted %s in all collisions", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_T", Form("Generated %s in TVX triggered collisions", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_TZ", Form("Generated %s in TVX collisions with z < 10cm", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_TZS", Form("Generated %s in TVX collisions with z < 10cm and Sel8", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_TZSG", Form("Generated %s in collisions with good z < 10cm and Sel8", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_TZSGU", Form("Generated %s in unique collisions with good z < 10cm and Sel8", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_TZSGUE", Form("Generated %s in unique TVXinEMC collisions with good z < 10cm and Sel8", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_TZSGUE_Accepted", Form("Accepted %s in unique TVXinEMC collisions with good z < 10cm and Sel8", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});

    mHistManager.add("Yield_BC_T", Form("Generated %s in TVX triggered BCs", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("Yield_BC_TC", Form("Generated %s in TVX triggered BCs that have at least one collision", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("NCollisionsMCCollisions", "Number of (MC)Collisions in the BC;#it{N}_(Collisions);#it{N}_(MC Collisions)", kTH2F, {{4, -0.5, 3.5}, {4, -0.5, 3.5}});
    mHistManager.add("NTVXCollisionsMCCollisions", "Number of (MC)Collisions in the TVX triggered BC;#it{N}_(Collisions);#it{N}_(MC Collisions)", kTH2F, {{4, -0.5, 3.5}, {4, -0.5, 3.5}});

    auto hEMCollisionCounter = mHistManager.add<TH1>("hEMCollisionCounter", "collision counter;;Number of events", kTH1F, {{13, 0.5, 13.5}}, false);
    hEMCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(10, "EMC MB Readout");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(11, "EMC L0 Triggered");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(12, "EMC Cell Content");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(13, "accepted");

    o2::emcal::Geometry::GetInstanceFromRunNumber(300000);
  }

  PresliceUnsorted<MyMCCollisions> perFoundBC = aod::evsel::foundBCId;
  Preslice<aod::McCollisions> MCCollperBC = aod::mccollision::bcId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  void process(MyBCs const& bcs, MyMCCollisions const& collisions, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {

    for (const auto& bc : bcs) {

      auto collisionsInFoundBC = collisions.sliceBy(perFoundBC, bc.globalIndex());
      auto MCCollisionsBC = mcCollisions.sliceBy(MCCollperBC, bc.globalIndex());

      mHistManager.fill(HIST("NCollisionsMCCollisions"), collisionsInFoundBC.size(), MCCollisionsBC.size());
      mHistManager.fill(HIST("hBCCounter"), 1);

      if (!mRequireTVX || bc.selection_bit(aod::evsel::kIsTriggerTVX)) { // Count BCs with TVX trigger with and without a collision, as well as the generated particles within

        mHistManager.fill(HIST("NTVXCollisionsMCCollisions"), collisionsInFoundBC.size(), mcCollisions.size());

        mHistManager.fill(HIST("hBCCounter"), 2);

        bool bcHasCollision = collisionsInFoundBC.size() > 0;

        if (bcHasCollision)
          mHistManager.fill(HIST("hBCCounter"), 3);

        for (auto& mcCollision : MCCollisionsBC) {

          auto mcParticles_inColl = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

          for (auto& mcParticle : mcParticles_inColl) {
            if (mSelectedParticleCode != 0 && mcParticle.pdgCode() != mSelectedParticleCode)
              continue;
            if (fabs(mcParticle.y()) > mRapidityCut)
              continue;
            if (!mcParticle.isPhysicalPrimary() && !mcParticle.producedByGenerator())
              continue;
            if (mRequireGammaGammaDecay && !isGammaGammaDecay(mcParticle, mcParticles))
              continue;

            mHistManager.fill(HIST("Yield_BC_T"), mcParticle.pt());

            if (bcHasCollision)
              mHistManager.fill(HIST("Yield_BC_TC"), mcParticle.pt());
          }
        }
      }
    }

    for (auto& collision : collisions) {
      fillEventHistogram(&mHistManager, collision);

      auto mcCollision = collision.mcCollision();
      auto mcParticles_inColl = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

      for (auto& mcParticle : mcParticles_inColl) {
        if (mSelectedParticleCode != 0 && mcParticle.pdgCode() != mSelectedParticleCode)
          continue;
        if (fabs(mcParticle.y()) > mRapidityCut)
          continue;
        if (!mcParticle.isPhysicalPrimary() && !mcParticle.producedByGenerator())
          continue;
        if (mRequireGammaGammaDecay && !isGammaGammaDecay(mcParticle, mcParticles))
          continue;

        mHistManager.fill(HIST("Yield"), mcParticle.pt());
        if (isAccepted(mcParticle, mcParticles))
          mHistManager.fill(HIST("Yield_Accepted"), mcParticle.pt());
        if (!mRequireTVX || collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
          mHistManager.fill(HIST("Yield_T"), mcParticle.pt());
          if (abs(collision.posZ()) < mVertexCut) {
            mHistManager.fill(HIST("Yield_TZ"), mcParticle.pt());
            if (!mRequireSel8 || collision.sel8()) {
              mHistManager.fill(HIST("Yield_TZS"), mcParticle.pt());
              if (!mRequireGoodZvtxFT0vsPV || collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
                mHistManager.fill(HIST("Yield_TZSG"), mcParticle.pt());
                if (!mRequireNoSameBunchPileup || collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
                  mHistManager.fill(HIST("Yield_TZSGU"), mcParticle.pt());
                  if (!mRequireEMCReadoutInMB || (mRequireEMCCellContent ? collision.isemcreadout() : collision.alias_bit(kTVXinEMC))) {
                    mHistManager.fill(HIST("Yield_TZSGUE"), mcParticle.pt());
                    if (isAccepted(mcParticle, mcParticles))
                      mHistManager.fill(HIST("Yield_TZSGUE_Accepted"), mcParticle.pt());
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  template <typename TMCParticle, typename TMCParticles>
  bool isGammaGammaDecay(TMCParticle mcParticle, TMCParticles mcParticles)
  {
    auto daughtersIds = mcParticle.daughtersIds();
    if (daughtersIds.size() != 2)
      return false;
    for (auto& daughterId : daughtersIds) {
      if (mcParticles.iteratorAt(daughterId).pdgCode() != 22)
        return false;
    }
    return true;
  }

  template <typename TMCParticle, typename TMCParticles>
  bool isAccepted(TMCParticle mcParticle, TMCParticles mcParticles)
  {
    auto daughtersIds = mcParticle.daughtersIds();
    if (daughtersIds.size() != 2)
      return false;
    for (auto& daughterId : daughtersIds) {
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

  void fillEventHistogram(HistogramRegistry* fRegistry, MyMCCollisions::iterator const& collision)
  {
    fRegistry->fill(HIST("hEMCollisionCounter"), 1.0);
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      fRegistry->fill(HIST("hEMCollisionCounter"), 2.0);
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      fRegistry->fill(HIST("hEMCollisionCounter"), 3.0);
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      fRegistry->fill(HIST("hEMCollisionCounter"), 4.0);
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      fRegistry->fill(HIST("hEMCollisionCounter"), 5.0);
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      fRegistry->fill(HIST("hEMCollisionCounter"), 6.0);
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX))
      fRegistry->fill(HIST("hEMCollisionCounter"), 7.0);
    if (collision.sel8())
      fRegistry->fill(HIST("hEMCollisionCounter"), 8.0);
    if (abs(collision.posZ()) < 10.0)
      fRegistry->fill(HIST("hEMCollisionCounter"), 9.0);
    if (collision.alias_bit(kTVXinEMC))
      fRegistry->fill(HIST("hEMCollisionCounter"), 10.0);
    if (collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7))
      fRegistry->fill(HIST("hEMCollisionCounter"), 11.0);
    if (collision.isemcreadout())
      fRegistry->fill(HIST("hEMCollisionCounter"), 12.0);
    fRegistry->fill(HIST("hEMCollisionCounter"), 13.0);

    fRegistry->fill(HIST("hCollisionCounter"), 1);
    if (!mRequireTVX || collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      fRegistry->fill(HIST("hCollisionCounter"), 2);
      if (abs(collision.posZ()) < mVertexCut) {
        fRegistry->fill(HIST("hCollisionCounter"), 3);
        if (!mRequireSel8 || collision.sel8()) {
          fRegistry->fill(HIST("hCollisionCounter"), 4);
          if (!mRequireGoodZvtxFT0vsPV || collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
            fRegistry->fill(HIST("hCollisionCounter"), 5);
            if (!mRequireNoSameBunchPileup || collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
              fRegistry->fill(HIST("hCollisionCounter"), 6);
              if (!mRequireEMCReadoutInMB || (mRequireEMCCellContent ? collision.isemcreadout() : collision.alias_bit(kTVXinEMC)))
                fRegistry->fill(HIST("hCollisionCounter"), 7);
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MCGeneratorStudies>(cfgc, TaskName{"mc-generator-studies"})};
}
