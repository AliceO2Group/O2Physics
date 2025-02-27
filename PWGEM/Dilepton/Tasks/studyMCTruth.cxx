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
// ========================
//
// This code is to study MC truth. e.g. S/B
//    Please write to: daiki.sekihata@cern.ch

#include <string>
#include "Math/Vector4D.h"

#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

struct studyMCTruth {

  struct : ConfigurableGroup {
    std::string prefix = "mccut_group";
    Configurable<int> cfgPdgCodeLepton{"cfgPdgCodeLepton", 11, "pdg code for desired lepton"};
    Configurable<float> cfgMinPtGen{"cfgMinPtGen", 0.1, "min. pT of single lepton"};
    Configurable<float> cfgMaxPtGen{"cfgMaxPtGen", 1e+10f, "max. pT of single lepton"};
    Configurable<float> cfgMinEtaGen{"cfgMinEtaGen", -0.8, "min. eta of for single lepton"};
    Configurable<float> cfgMaxEtaGen{"cfgMaxEtaGen", +0.8, "max. eta of for single lepton"};
  } mccuts;

  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<float> cfgMinImpPar{"cfgMinImpPar", -1.f, "min. impact parameter in fm"}; // [0, 4] fm for centrality FT0C 0-10%,  [8, 10] fm for centrality FT0C 30-50% in PbPb
    Configurable<float> cfgMaxImpPar{"cfgMaxImpPar", 999.f, "max. impact parameter in fm"};
  } eventcuts;

  ConfigurableAxis ConfMllBins{"ConfMllBins", {400, 0.f, 4.f}, "mll bins for output histograms"};
  ConfigurableAxis ConfPtllBins{"ConfPtllBins", {100, 0.f, 10.f}, "pTll bins for output histograms"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  float leptonMass = o2::constants::physics::MassElectron;
  void init(o2::framework::InitContext&)
  {
    if (std::abs(mccuts.cfgPdgCodeLepton.value) == 11) {
      leptonMass = o2::constants::physics::MassElectron;
    } else if (std::abs(mccuts.cfgPdgCodeLepton.value) == 13) {
      leptonMass = o2::constants::physics::MassMuon;
    } else {
      LOGF(fatal, "pdg code must be 11 or 13.");
    }
    addHistograms();
  }

  static constexpr std::string_view dileptonSigns[3] = {"uls/", "lspp/", "lsmm/"};
  static constexpr std::string_view evNames[4] = {"allMC/", "selectedMC/", "selectedMC_and_Rec/", "selectedMC_and_selectedRec/"};

  void addHistograms()
  {
    const AxisSpec axis_mll{ConfMllBins, "m_{ll} (GeV/c^{2})"};
    const AxisSpec axis_ptll{ConfPtllBins, "p_{T,ll} (GeV/c)"};

    fRegistry.add("Event/hReccollsPerMCcoll", "Rec. colls per MC coll;Rec. colls per MC coll;Number of MC collisions", kTH1D, {{21, -0.5, 20.5}}, false);
    fRegistry.add("Event/hSelReccollsPerMCcoll", "Selected Rec. colls per MC coll;Selected Rec. colls per MC coll;Number of MC collisions", kTH1D, {{21, -0.5, 20.5}}, false);
    fRegistry.add("Event/hDiffBC", "diffrence in BC;BC_{rec. coll.} - BC_{mc coll.}", kTH1D, {{101, -50.5, +50.5}}, false);
    fRegistry.add("Event/allMC/hZvtx", "MC Zvtx;Z_{vtx} (cm)", kTH1D, {{100, -50, +50}}, false);
    fRegistry.add("Event/allMC/hImpactParameter", "impact parameter;impact parameter b (fm)", kTH1D, {{200, 0, 20}}, false);
    fRegistry.addClone("Event/allMC/", "Event/selectedMC/");
    fRegistry.addClone("Event/allMC/", "Event/selectedMC_and_Rec/");
    fRegistry.addClone("Event/allMC/", "Event/selectedMC_and_selectedRec/");

    fRegistry.add("Pair/allMC/Pi0/uls/hMvsPt", "m_{ll} vs. p_{T,ll}", kTH2D, {axis_mll, axis_ptll}, true);
    fRegistry.addClone("Pair/allMC/Pi0/uls/", "Pair/allMC/Pi0/lspp/");
    fRegistry.addClone("Pair/allMC/Pi0/uls/", "Pair/allMC/Pi0/lsmm/");
    fRegistry.addClone("Pair/allMC/Pi0/", "Pair/allMC/Eta/");
    fRegistry.addClone("Pair/allMC/Pi0/", "Pair/allMC/EtaPrime/");
    fRegistry.addClone("Pair/allMC/Pi0/", "Pair/allMC/Rho/");
    fRegistry.addClone("Pair/allMC/Pi0/", "Pair/allMC/Omega/");
    fRegistry.addClone("Pair/allMC/Pi0/", "Pair/allMC/Phi/");
    fRegistry.addClone("Pair/allMC/Pi0/", "Pair/allMC/JPsi/");
    fRegistry.addClone("Pair/allMC/Pi0/", "Pair/allMC/ccbar/");
    fRegistry.addClone("Pair/allMC/Pi0/", "Pair/allMC/bbbar/");
    fRegistry.addClone("Pair/allMC/", "Pair/selectedMC/");
    fRegistry.addClone("Pair/allMC/", "Pair/selectedMC_and_Rec/");
    fRegistry.addClone("Pair/allMC/", "Pair/selectedMC_and_selectedRec/");

    // allMC = all mc collisions
    // selectedMC = mc collisions selected by your event cuts
    // selectedMC_and_Rec = mc collisions selected by your event cuts and at least 1 reconstructed collision is associated.
    // selectedMC_and_selectedRec = mc collisions selected by your event cuts and at least 1 reconstructed collision is associated, and the associated rec. collision is selected by your cuts.
  }

  template <typename TTrack, typename TMCParticles>
  int FindLF(TTrack const& posmc, TTrack const& negmc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(posmc, negmc, -mccuts.cfgPdgCodeLepton, mccuts.cfgPdgCodeLepton, 111, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -mccuts.cfgPdgCodeLepton, mccuts.cfgPdgCodeLepton, 221, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -mccuts.cfgPdgCodeLepton, mccuts.cfgPdgCodeLepton, 331, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -mccuts.cfgPdgCodeLepton, mccuts.cfgPdgCodeLepton, 113, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -mccuts.cfgPdgCodeLepton, mccuts.cfgPdgCodeLepton, 223, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -mccuts.cfgPdgCodeLepton, mccuts.cfgPdgCodeLepton, 333, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -mccuts.cfgPdgCodeLepton, mccuts.cfgPdgCodeLepton, 443, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <typename TMCParticle>
  bool isSelectedMCParticle(TMCParticle const& mcparticle)
  {
    if (std::abs(mcparticle.pdgCode()) != mccuts.cfgPdgCodeLepton) {
      return false;
    }
    if (!mcparticle.has_mothers()) {
      return false;
    }
    if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TCollision, typename TBC>
  bool isSelectedCollision(TCollision const& collision, TBC const& bc)
  {
    if (collision.posZ() < eventcuts.cfgZvtxMin || eventcuts.cfgZvtxMax < collision.posZ()) {
      return false;
    }

    if (eventcuts.cfgRequireFT0AND && !bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (eventcuts.cfgRequireNoTFB && !bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (eventcuts.cfgRequireNoITSROFB && !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    return true;
  }

  template <int evtype, int signtype, typename TMCLepton, typename TMCParticles>
  void fillTrueInfo(TMCLepton const& t1, TMCLepton const& t2, TMCParticles const& mcParticles)
  {
    if (!isSelectedMCParticle(t1) || !isSelectedMCParticle(t2)) {
      return;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonMass);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonMass);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    if (v12.Rapidity() < mccuts.cfgMinEtaGen || mccuts.cfgMaxEtaGen < v12.Rapidity()) {
      return;
    }

    int mother_id = FindLF(t1, t2, mcParticles);
    int hfll_type = IsHF(t1, t2, mcParticles);

    if (mother_id > 0) { // same mother (photon, LF, Quarkonia)
      const auto& mp = mcParticles.iteratorAt(mother_id);
      if (!(mp.isPhysicalPrimary() || mp.producedByGenerator())) {
        return;
      }
      switch (std::abs(mp.pdgCode())) {
        case 111:
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("Pi0/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 221:
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("Eta/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 331:
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("EtaPrime/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 113:
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("Rho/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 223:
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("Omega/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 333:
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("Phi/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 443:
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("JPsi/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        default:
          break;
      }
    } else if (hfll_type > -1) { // HFll
      switch (hfll_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce): // ULS
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("ccbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case static_cast<int>(EM_HFeeType::kBe_Be): // ULS
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("bbbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_BCe): // ULS
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("bbbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): // ULS
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("bbbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): // LS
          fRegistry.fill(HIST("Pair/") + HIST(evNames[evtype]) + HIST("bbbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
        default:
          break;
      }
    }
  }

  template <typename TMCCollisions, typename TMCParticles, typename TBCs, typename TCollisions>
  void runMC(TMCCollisions const& mcCollisions, TMCParticles const& mcParticles, TBCs const&, TCollisions const& collisions)
  {
    for (const auto& mcCollision : mcCollisions) {
      const auto& bc_from_mcCollision = mcCollision.template bc_as<TBCs>();
      bool isSelectedMC = isSelectedCollision(mcCollision, bc_from_mcCollision);

      auto reccolls_per_mccoll = collisions.sliceBy(recColperMcCollision, mcCollision.globalIndex());
      fRegistry.fill(HIST("Event/hReccollsPerMCcoll"), reccolls_per_mccoll.size());
      int nselreccolls_per_mccoll = 0;
      for (const auto& rec_col : reccolls_per_mccoll) {
        if (isSelectedCollision(rec_col, rec_col.template foundBC_as<TBCs>())) {
          nselreccolls_per_mccoll++;
        }
      } // end of reconstructed collision
      fRegistry.fill(HIST("Event/hSelReccollsPerMCcoll"), nselreccolls_per_mccoll);

      bool isSelectedRec = false;
      bool hasRecCollision = false;
      if (mcCollision.mpemeventId() >= 0) {
        hasRecCollision = true;
        const auto& collision = collisions.rawIteratorAt(mcCollision.mpemeventId()); // most probable reconstructed collision
        const auto& bc_from_collision = collision.template foundBC_as<TBCs>();
        isSelectedRec = isSelectedCollision(collision, bc_from_collision);
        fRegistry.fill(HIST("Event/hDiffBC"), bc_from_collision.globalBC() - bc_from_mcCollision.globalBC());
      }

      fRegistry.fill(HIST("Event/allMC/hZvtx"), mcCollision.posZ());
      fRegistry.fill(HIST("Event/allMC/hImpactParameter"), mcCollision.impactParameter());
      if (isSelectedMC) {
        fRegistry.fill(HIST("Event/selectedMC/hZvtx"), mcCollision.posZ());
        fRegistry.fill(HIST("Event/selectedMC/hImpactParameter"), mcCollision.impactParameter());
        if (hasRecCollision) {
          fRegistry.fill(HIST("Event/selectedMC_and_Rec/hZvtx"), mcCollision.posZ());
          fRegistry.fill(HIST("Event/selectedMC_and_Rec/hImpactParameter"), mcCollision.impactParameter());
          if (isSelectedRec) {
            fRegistry.fill(HIST("Event/selectedMC_and_selectedRec/hZvtx"), mcCollision.posZ());
            fRegistry.fill(HIST("Event/selectedMC_and_selectedRec/hImpactParameter"), mcCollision.impactParameter());
          }
        }
      }

      // store MC true information
      auto posLeptons_per_mccollision = mcPosLeptons.sliceBy(perMcCollision, mcCollision.globalIndex());
      auto negLeptons_per_mccollision = mcNegLeptons.sliceBy(perMcCollision, mcCollision.globalIndex());

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posLeptons_per_mccollision, negLeptons_per_mccollision))) { // ULS
        fillTrueInfo<0, 0>(pos, neg, mcParticles);
        if (isSelectedMC) {
          fillTrueInfo<1, 0>(pos, neg, mcParticles);
          if (hasRecCollision) {
            fillTrueInfo<2, 0>(pos, neg, mcParticles);
            if (isSelectedRec) {
              fillTrueInfo<3, 0>(pos, neg, mcParticles);
            }
          }
        }
      } // end of ULS pair loop

      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posLeptons_per_mccollision, posLeptons_per_mccollision))) { // LS++
        fillTrueInfo<0, 1>(pos1, pos2, mcParticles);
        if (isSelectedMC) {
          fillTrueInfo<1, 1>(pos1, pos2, mcParticles);
          if (hasRecCollision) {
            fillTrueInfo<2, 1>(pos1, pos2, mcParticles);
            if (isSelectedRec) {
              fillTrueInfo<3, 1>(pos1, pos2, mcParticles);
            }
          }
        }
      } // end of LS++ pair loop

      for (auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negLeptons_per_mccollision, negLeptons_per_mccollision))) { // LS--
        fillTrueInfo<0, 2>(neg1, neg2, mcParticles);
        if (isSelectedMC) {
          fillTrueInfo<1, 2>(neg1, neg2, mcParticles);
          if (hasRecCollision) {
            fillTrueInfo<2, 2>(neg1, neg2, mcParticles);
            if (isSelectedRec) {
              fillTrueInfo<3, 2>(neg1, neg2, mcParticles);
            }
          }
        }
      } // end of LS-- pair loop

    } // end of mc collision loop

  } //  end of skimmingMC

  SliceCache cache;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<aod::McCollisionLabels> recColperMcCollision = aod::mccollisionlabel::mcCollisionId;

  using MyMcCollisions = soa::Join<aod::McCollisions, aod::MostProbableEMEventIdsInMC>;

  Filter collisionFilter = eventcuts.cfgMinImpPar < o2::aod::mccollision::impactParameter && o2::aod::mccollision::impactParameter < eventcuts.cfgMaxImpPar;
  using FilteredMyMcCollisions = soa::Filtered<MyMcCollisions>;

  Partition<aod::McParticles> mcPosLeptons = o2::aod::mcparticle::pdgCode == -mccuts.cfgPdgCodeLepton && (mccuts.cfgMinPtGen < o2::aod::mcparticle::pt && o2::aod::mcparticle::pt < mccuts.cfgMaxPtGen) && (mccuts.cfgMinEtaGen < o2::aod::mcparticle::eta && o2::aod::mcparticle::eta < mccuts.cfgMaxEtaGen);
  Partition<aod::McParticles> mcNegLeptons = o2::aod::mcparticle::pdgCode == mccuts.cfgPdgCodeLepton && (mccuts.cfgMinPtGen < o2::aod::mcparticle::pt && o2::aod::mcparticle::pt < mccuts.cfgMaxPtGen) && (mccuts.cfgMinEtaGen < o2::aod::mcparticle::eta && o2::aod::mcparticle::eta < mccuts.cfgMaxEtaGen);

  void processMC(FilteredMyMcCollisions const& mcCollisions, aod::McParticles const& mcParticles, soa::Join<aod::BCsWithTimestamps, aod::BcSels> const& bcs, soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions)
  {
    runMC(mcCollisions, mcParticles, bcs, collisions);
  }
  PROCESS_SWITCH(studyMCTruth, processMC, "process", true);

  void processDummy(FilteredMyMcCollisions const&) {}
  PROCESS_SWITCH(studyMCTruth, processDummy, "process Dummy", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<studyMCTruth>(cfgc, TaskName{"study-mc-truth"})};
}
