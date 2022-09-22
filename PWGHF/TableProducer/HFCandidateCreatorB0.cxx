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

/// \file HFCandidateCreatorB0.cxx
/// \brief Reconstruction of B0 candidates
/// \note Adapted from HFCandidateCreatorXicc
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_cand_b0;	// from HFSecondaryVertex.h
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Perform MC matching."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Reconstruction of B0 candidates
struct HFCandidateCreatorB0 {
  Produces<aod::HfCandB0Base> rowCandidateBase;  // table defined in HFSecondaryVertex.h

  Configurable<double> magneticField{"magneticField", 20., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};

  OutputObj<TH1F> hMassDToPiKPi{TH1F("hMassB0ToPiKPi", "D^{#minus} candidates;inv. mass (p^{#minus} K^{#plus} #pi^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD{TH1F("hPtD", "D^{#minus} candidates;D^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi^{#plus} candidates;#pi^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPAD{TH1F("hCPAD", "D^{#minus} candidates;D^{#minus} cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassB0ToDPi{TH1F("hMassB0ToDPi", "2-prong candidates;inv. mass (B^{0} #rightarrow D^{#minus}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx position (cm^{2});entries", 100, 0., 0.2)};

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massD = RecoDecay::getMassPDG(pdg::Code::kDMinus);
  double massDPi = 0.;

  Configurable<int> d_selectionFlagD{"d_selectionFlagD", 1, "Selection Flag for D"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Filter filterSelectCandidates = (aod::hf_selcandidate_dplus::isSelDplusToPiKPi >= d_selectionFlagD); //FIXME

  void process(aod::Collision const& collision,
               soa::Filtered<soa::Join<
                 aod::HfCandProng3,
                 aod::HFSelDplusToPiKPiCandidate>> const& dCands,
               aod::BigTracks const& tracks)
  {
    // Initialise fitter for B vertex (2-prong vertex filter)
    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(magneticField);
    df2.setPropagateToPCA(b_propdca);
    df2.setMaxR(d_maxr);
    df2.setMaxDZIni(d_maxdzini);
    df2.setMinParamChange(d_minparamchange);
    df2.setMinRelChi2Change(d_minrelchi2change);
    df2.setUseAbsDCA(true);

    // Initial fitter to redo D-vertex to get extrapolated daughter tracks (3-prong vertex filter)
    o2::vertexing::DCAFitterN<3> df3;
    df3.setBz(magneticField);
    df3.setPropagateToPCA(b_propdca);
    df3.setMaxR(d_maxr);
    df3.setMaxDZIni(d_maxdzini);
    df3.setMinParamChange(d_minparamchange);
    df3.setMinRelChi2Change(d_minrelchi2change);
    df3.setUseAbsDCA(true);

    // loop over D candidates
    for (auto& dCand : dCands) {
      if (!TESTBIT(dCand.hfflag(), hf_cand_prong3::DecayType::DPlusToPiKPi)) {
        continue;
      }
      if (dCand.isSelDplusToPiKPi() >= d_selectionFlagD) {
        hMassDToPiKPi->Fill(InvMassDPlus(dCand), dCand.pt());
      }
      hPtD->Fill(dCand.pt());
      hCPAD->Fill(dCand.cpa());

			// track0 <-> pi, track1 <-> K, track2 <-> pi
      auto track0 = dCand.index0_as<aod::BigTracks>();
      auto track1 = dCand.index1_as<aod::BigTracks>();
      auto track2 = dCand.index2_as<aod::BigTracks>();
      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);
      auto collision = track0.collision();

      // reconstruct 3-prong secondary vertex (D±)
      if (df3.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
        continue;
      }
      
      const auto& secondaryVertex = df3.getPCACandidate();
      trackParVar0.propagateTo(secondaryVertex[0], magneticField);
      trackParVar1.propagateTo(secondaryVertex[0], magneticField);
      trackParVar2.propagateTo(secondaryVertex[0], magneticField);

			// D∓ → π∓ K± π∓
      array<float, 3> pvecpiK = {track0.px() + track1.px(), track0.py() + track1.py(), track0.pz() + track1.pz()};
      array<float, 3> pvecD = {pvecpiK[0] + track2.px(), pvecpiK[1] + track2.py(), pvecpiK[2] + track2.pz()};
      auto trackpiK = o2::dataformats::V0(df3.getPCACandidatePos(), pvecpiK, df3.calcPCACovMatrixFlat(),
                                         trackParVar0, trackParVar1, {0, 0}, {0, 0});
      auto trackD = o2::dataformats::V0(df3.getPCACandidatePos(), pvecD, df3.calcPCACovMatrixFlat(),
                                         trackpiK, trackParVar2, {0, 0}, {0, 0});

      int index0D = track0.globalIndex();
      int index1D = track1.globalIndex();
      int index2D = track2.globalIndex();
      //int charge = track0.sign() + track1.sign() + track2.sign();
      
      // loop on D-
      // D- → π- K+ π-
      // we don't have direct access to D sign so we use the sign of the daughters (the pion track0 here)
      if (track0.sign() < 0) {
				// loop over pions      
		    for (auto& trackPion : tracks) {
		    	// minimum pT selection
		      if (trackPion.pt() < ptPionMin) {
		        continue;
		      }
		      // we reject pions that are D daughters
		      if (trackPion.globalIndex() == index0D || trackPion.globalIndex() == index1D || trackPion.globalIndex() == index2D) {
		        continue;
		      }
		      // we only keep pi+ to combine them with Dminus and reconstruct B0
		      if (trackPion.sign() < 0) {
		        continue;
		      }

		      hPtPion->Fill(trackPion.pt());
		      array<float, 3> pvecPion;
		      auto trackParVarPi = getTrackParCov(trackPion);

					// ---------------------------------
		      // reconstruct the 2-prong B0 vertex   
		      if (df2.process(trackD, trackParVarPi) == 0) {
		        continue;
		      }

		      // calculate relevant properties
		      const auto& secondaryVertexB0 = df2.getPCACandidate();
		      auto chi2PCA = df2.getChi2AtPCACandidate();
		      auto covMatrixPCA = df2.calcPCACovMatrixFlat();

		      df2.propagateTracksToVertex();
		      df2.getTrack(0).getPxPyPzGlo(pvecD);
		      df2.getTrack(1).getPxPyPzGlo(pvecPion);

		      auto primaryVertex = getPrimaryVertex(collision);
		      auto covMatrixPV = primaryVertex.getCov();
		      o2::dataformats::DCA impactParameter0;
		      o2::dataformats::DCA impactParameter1;
		      trackD.propagateToDCA(primaryVertex, magneticField, &impactParameter0);
		      trackParVarPi.propagateToDCA(primaryVertex, magneticField, &impactParameter1);

		      hCovSVXX->Fill(covMatrixPCA[0]);
		      hCovPVXX->Fill(covMatrixPV[0]);

		      // get uncertainty of the decay length
		      double phi, theta;
		      getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
		      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
		      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

		      int hfFlag = BIT(hf_cand_b0::DecayType::B0ToDPi);

		      // fill the candidate table for the B0 here:
		      rowCandidateBase(collision.globalIndex(),
		                       collision.posX(), collision.posY(), collision.posZ(),
		                       secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
		                       errorDecayLength, errorDecayLengthXY,
		                       chi2PCA,
		                       pvecD[0], pvecD[1], pvecD[2],
		                       pvecPion[0], pvecPion[1], pvecPion[2],
		                       impactParameter0.getY(), impactParameter1.getY(),
		                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
		                       dCand.globalIndex(), trackPion.globalIndex(),
		                       hfFlag);

		      // calculate invariant mass
		      auto arrayMomenta = array{pvecD, pvecPion};
		      massDPi = RecoDecay::m(std::move(arrayMomenta), array{massD, massPi});
		      if (dCand.isSelDplusToPiKPi() > 0) {
		        hMassB0ToDPi->Fill(massDPi);
		      }
		    } // pi+ loop
      } // if D-

      // loop on D+
      // D+ → π+ K- π+
      // we now loop on the D+
      if (track0.sign() > 0) {
				// loop over pions      
		    for (auto& trackPion : tracks) {
		    	// minimum pT selection
		      if (trackPion.pt() < ptPionMin) {
		        continue;
		      }
		      // we reject pions that are D daughters
		      if (trackPion.globalIndex() == index0D || trackPion.globalIndex() == index1D || trackPion.globalIndex() == index2D) {
		        continue;
		      }
		      // we onnly keep pi- to combine them with D+ and reconstruct B0bar
		      if (trackPion.sign() > 0) {
		        continue;
		      }

		      hPtPion->Fill(trackPion.pt());
		      array<float, 3> pvecPion;
		      auto trackParVarPi = getTrackParCov(trackPion);

					// ---------------------------------
		      // reconstruct the 2-prong B0bar vertex   
		      if (df2.process(trackD, trackParVarPi) == 0) {
		        continue;
		      }

		      // calculate relevant properties
		      const auto& secondaryVertexB0 = df2.getPCACandidate();
		      auto chi2PCA = df2.getChi2AtPCACandidate();
		      auto covMatrixPCA = df2.calcPCACovMatrixFlat();

		      df2.propagateTracksToVertex();
		      df2.getTrack(0).getPxPyPzGlo(pvecD);
		      df2.getTrack(1).getPxPyPzGlo(pvecPion);

		      auto primaryVertex = getPrimaryVertex(collision);
		      auto covMatrixPV = primaryVertex.getCov();
		      o2::dataformats::DCA impactParameter0;
		      o2::dataformats::DCA impactParameter1;
		      trackD.propagateToDCA(primaryVertex, magneticField, &impactParameter0);
		      trackParVarPi.propagateToDCA(primaryVertex, magneticField, &impactParameter1);

		      hCovSVXX->Fill(covMatrixPCA[0]);
		      hCovPVXX->Fill(covMatrixPV[0]);

		      // get uncertainty of the decay length
		      double phi, theta;
		      getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
		      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
		      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

		      int hfFlag = BIT(hf_cand_b0::DecayType::B0ToDPi);

		      // fill the candidate table for the B0 here:
		      rowCandidateBase(collision.globalIndex(),
		                       collision.posX(), collision.posY(), collision.posZ(),
		                       secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
		                       errorDecayLength, errorDecayLengthXY,
		                       chi2PCA,
		                       pvecD[0], pvecD[1], pvecD[2],
		                       pvecPion[0], pvecPion[1], pvecPion[2],
		                       impactParameter0.getY(), impactParameter1.getY(),
		                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
		                       dCand.globalIndex(), trackPion.globalIndex(),
		                       hfFlag);

		      // calculate invariant mass
		      auto arrayMomenta = array{pvecD, pvecPion};
		      massDPi = RecoDecay::m(std::move(arrayMomenta), array{massD, massPi});
		      if (dCand.isSelDplusToPiKPi() > 0) {
		        hMassB0ToDPi->Fill(massDPi);
		      }
		    } // pi- loop
      } // if D+
    }   // D loop
  }     // process
};      // struct

/// Extends the base table with expression columns.
struct HFCandidateCreatorB0Expressions {
  Spawns<aod::HfCandB0Ext> rowCandidateB0;
  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HFCandidateCreatorB0MC {
  Produces<aod::HfCandB0MCRec> rowMCMatchRec; // table defined in HFSecondaryVertex.h
  Produces<aod::HfCandB0MCGen> rowMCMatchGen; // table defined in HFSecondaryVertex.h

  void process(aod::HfCandB0 const& candidates,
               aod::HfCandProng3 const&,
               aod::BigTracksMC const& tracks,
               aod::McParticles const& particlesMC)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      //Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      debug = 0;
      auto dCand = candidate.index0();
      auto arrayDaughters = array{dCand.index0_as<aod::BigTracksMC>(),
                                  dCand.index1_as<aod::BigTracksMC>(),
                                  dCand.index2_as<aod::BigTracksMC>(),
                                  candidate.index1_as<aod::BigTracksMC>()};
      auto arrayDaughtersD = array{dCand.index0_as<aod::BigTracksMC>(),
                                  dCand.index1_as<aod::BigTracksMC>(),
                                  dCand.index2_as<aod::BigTracksMC>()};
      // B0 → D- π+ → (π- K+ π-) π+
      //Printf("Checking B0 → D- π+");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kB0, array{-kPiPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        // B0 → D- π+
        //Printf("Checking B0 → D- π+");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersD, - pdg::Code::kDPlus, array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign, 1);
        if (indexRec > -1) {
          flag = BIT(hf_cand_b0::DecayType::B0ToDPi);
        } else {
          debug = 1;
          LOGF(info, "WARNING: B0 in decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }
      }
      rowMCMatchRec(flag, origin, debug);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      //Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      // B0 → D- π+
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kB0, array{- int(pdg::Code::kDPlus), +kPiPlus}, true)) {
        // Match D- -> π- K+ π-
        auto DCandMC = particlesMC.rawIteratorAt(particle.daughtersIds().front());
        //Printf("Checking D- -> π- K+ π-");
        if (RecoDecay::isMatchedMCGen(particlesMC, DCandMC, - int(pdg::Code::kDPlus), array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign)) {
          flag = sign * BIT(hf_cand_b0::DecayType::B0ToDPi);
        }
      }
      rowMCMatchGen(flag, origin);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HFCandidateCreatorB0>(cfgc, TaskName{"hf-candidate-creator-b0"}),
    adaptAnalysisTask<HFCandidateCreatorB0Expressions>(cfgc, TaskName{"hf-candidate-creator-b0-expressions"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HFCandidateCreatorB0MC>(cfgc, TaskName{"hf-candidate-creator-b0-mc"}));
  }
  return workflow;
}
