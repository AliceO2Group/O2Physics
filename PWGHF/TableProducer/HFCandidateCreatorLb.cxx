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

/// \file HFCandidateCreatorX.cxx
/// \brief Reconstruction of X(3872) candidates
/// \note Adapted from HFCandidateCreator3Prong
///
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef

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
using namespace o2::aod::hf_cand_lb;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Perform MC matching."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Reconstruction of Λb candidates
struct HFCandidateCreatorLb {
  Produces<aod::HfCandLbBase> rowCandidateBase;

  Configurable<double> magneticField{"magneticField", 5., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any Lb is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> ptPionMin{"ptPionMin", 0.15, "minimum pion pT threshold (GeV/c)"};
  Configurable<bool> b_dovalplots{"b_dovalplots", true, "do validation plots"};
  
  OutputObj<TH1F> hMassLcToPKPi{TH1F("hMassLcToPKPi", "Lc+ candidates;inv. mass (#mu^{#plus} #mu^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtLc{TH1F("hPtLc", "J/#psi candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi candidates;candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPALc{TH1F("hCPALc", "J/#psi candidates;cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassLbToLcPi{TH1F("hMassLbToLcPi", "3-prong candidates;inv. mass (Lc+ (#rightarrow pKpi) #pi-) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx position (cm^{2});entries", 100, 0., 0.2)};
  
  double massPi = RecoDecay::getMassPDG(kPiMinus);
  double massLc = RecoDecay::getMassPDG(4122);
  double massLcPi;
  
  Configurable<int> d_selectionFlagLc{"d_selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Filter filterSelectCandidates = (aod::hf_selcandidate_lc::isSelLcpKpi >= d_selectionFlagLc || aod::hf_selcandidate_lc::isSelLcpiKp >= d_selectionFlagLc);

  void process(aod::Collision const& collision,
               soa::Filtered<soa::Join<
               aod::HfCandProng3,
               aod::HFSelLcCandidate>> const& lcCands,
               aod::BigTracks const& tracks)
  {
    // 3-prong vertex fitter (to rebuild Lc vertex)
    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(magneticField);
    df2.setPropagateToPCA(b_propdca);
    df2.setMaxR(d_maxr);
    df2.setMaxDZIni(d_maxdzini);
    df2.setMinParamChange(d_minparamchange);
    df2.setMinRelChi2Change(d_minrelchi2change);
    df2.setUseAbsDCA(true);

    // 3-prong vertex fitter
    o2::vertexing::DCAFitterN<3> df3;
    df3.setBz(magneticField);
    df3.setPropagateToPCA(b_propdca);
    df3.setMaxR(d_maxr);
    df3.setMaxDZIni(d_maxdzini);
    df3.setMinParamChange(d_minparamchange);
    df3.setMinRelChi2Change(d_minrelchi2change);
    df3.setUseAbsDCA(true);

    // loop over Lc candidates
    for (auto& lcCand : lcCands) {
      if (!(lcCand.hfflag() & 1 << hf_cand_prong3::DecayType::LcToPKPi) && !(lcCand.hfflag() & 1 << hf_cand_prong3::DecayType::LcToPKPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YLc(lcCand)) > cutYCandMax) {
        continue;
      }
      if (lcCand.isSelLcpKpi() > 0) {
        hMassLcToPKPi->Fill(InvMassLcpKpi(lcCand));
      }
      //if (jpsiCand.isSelLcToPiKP() > 0) {
        //hMassLcToPKPi->Fill(InvMassLcToMuMu(jpsiCand));
      //}
      hPtLc->Fill(lcCand.pt());
      hCPALc->Fill(lcCand.cpa());
      // create Lc track to pass to DCA fitter; use cand table + rebuild vertex
      const std::array<float, 3> vertexLc = {lcCand.xSecondaryVertex(), lcCand.ySecondaryVertex(), lcCand.zSecondaryVertex()};
      array<float, 3> pvecLc = {lcCand.px(), lcCand.py(), lcCand.pz()};
      auto prong0 = lcCand.index0_as<aod::BigTracks>();
      auto prong1 = lcCand.index1_as<aod::BigTracks>();
      auto prong2 = lcCand.index2_as<aod::BigTracks>();
      auto prong0TrackParCov = getTrackParCov(prong0);
      auto prong1TrackParCov = getTrackParCov(prong1);
      auto prong2TrackParCov = getTrackParCov(prong1);

      if (df3.process(prong0TrackParCov, prong1TrackParCov, prong2TrackParCov) == 0) {
        continue;
      }

      // propogate prong tracks to Lc vertex
      prong0TrackParCov.propagateTo(lcCand.xSecondaryVertex(), magneticField);
      prong1TrackParCov.propagateTo(lcCand.xSecondaryVertex(), magneticField);
      prong2TrackParCov.propagateTo(lcCand.xSecondaryVertex(), magneticField);
      const std::array<float, 6> covLc = df3.calcPCACovMatrixFlat();
      // define the Lc track
      auto trackLc = o2::dataformats::V0(vertexLc, pvecLc, covLc, prong0TrackParCov, prong1TrackParCov, {0, 0}, {0, 0}); //FIXME: also needs covxyz???

      // used to check that prongs used for Lc and X reco are not the same prongs
      int index0Lc = lcCand.index0Id();
      int index1Lc = lcCand.index1Id();
      int index2Lc = lcCand.index2Id();

      // loop over pi- candidates
      for (auto& trackNeg : tracks) {
	if (trackNeg.sign() > 0) { // select only negative tracks - use partitions?
	  continue;
	}
	if (trackNeg.globalIndex() == index2Lc) {
	  continue;
	}
	if (trackNeg.pt() < ptPionMin) {
	  continue;
	}
	hPtPion->Fill(trackNeg.pt());

	auto trackParVarNeg = getTrackParCov(trackNeg);
	array<float, 3> pvecNeg;

	// reconstruct the 3-prong X vertex
	if (df2.process(trackLc, trackParVarNeg) == 0) {
	  continue;
	}

	// calculate relevant properties
	const auto& LbsecondaryVertex = df2.getPCACandidate();
	auto chi2PCA = df2.getChi2AtPCACandidate();
	auto covMatrixPCA = df2.calcPCACovMatrix().Array();
	hCovSVXX->Fill(covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.

	df2.propagateTracksToVertex();          // propagate the pions and Lc to the X vertex
	df2.getTrack(0).getPxPyPzGlo(pvecLc); // update momentum of Lc at the Lb vertex
	df2.getTrack(1).getPxPyPzGlo(pvecNeg);  // momentum of pi- at the Lb vertex

	// get track impact parameters
	// This modifies track momenta!
	auto primaryVertex = getPrimaryVertex(collision);
	auto covMatrixPV = primaryVertex.getCov();
	hCovPVXX->Fill(covMatrixPV[0]);
	o2::dataformats::DCA impactParameter0;
	o2::dataformats::DCA impactParameter1;
	trackParVarNeg.propagateToDCA(primaryVertex, magneticField, &impactParameter1);

	// get uncertainty of the decay length
	double phi, theta;
	getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, LbsecondaryVertex, phi, theta);
	auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
	auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

	int hfFlag = 0;
	if (TESTBIT(lcCand.hfflag(), hf_cand_prong3::DecayType::LcToPKPi)) {
            SETBIT(hfFlag, hf_cand_lb::DecayType::LbToLcPi); 
	}
	/* if (TESTBIT(jpsiCand.hfflag(), hf_cand_prong2::DecayType::LcToEE)) {
            SETBIT(hfFlag, hf_cand_x::DecayType::XToLcToEEPiPi); // dielectron channel
	    }*/

	// fill the candidate table for the Lb here:
	rowCandidateBase(collision.globalIndex(),
			 collision.posX(), collision.posY(), collision.posZ(),
			 LbsecondaryVertex[0], LbsecondaryVertex[1], LbsecondaryVertex[2],
			 errorDecayLength, errorDecayLengthXY,
			 chi2PCA,
			 pvecLc[0], pvecLc[1], pvecLc[2],
			 pvecNeg[0], pvecNeg[1], pvecNeg[2],
			 impactParameter0.getY(), impactParameter1.getY(), 
			 std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()), 
			 lcCand.globalIndex(), trackNeg.globalIndex(),
			 hfFlag);
	
	// calculate invariant mass
	auto arrayMomenta = array{pvecLc, pvecNeg};
	massLcPi = RecoDecay::M(std::move(arrayMomenta), array{massLc, massPi});
	if (lcCand.isSelLcpKpi() > 0) {
	  hMassLbToLcPi->Fill(massLcPi);
	}
      } // pi- loop
    }     // Lc loop
  }       // process
};        // struct

/// Extends the base table with expression columns.
struct HFCandidateCreatorLbExpressions {
  Spawns<aod::HfCandLbExt> rowCandidateLb;
  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HFCandidateCreatorLbMC {
  Produces<aod::HfCandLbMCRec> rowMCMatchRec;
  Produces<aod::HfCandLbMCGen> rowMCMatchGen;

  void process(aod::HfCandLb const& candidates,
               aod::HfCandProng3,
               aod::BigTracksMC const& tracks,
               aod::McParticles const& particlesMC)
  {
    int indexRec = -1;
    int indexRecLc = -1;
    int pdgCodeLb = 4122;
    int pdgCodeLc = 5122;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t channel = 0;

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      //Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      channel = 0;
      //auto lcTrack = candidate.index0();
      auto lcTrack = candidate.index0_as<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>>();
      auto daughter0Lc = lcTrack.index0_as<aod::BigTracksMC>();
      auto daughter1Lc = lcTrack.index1_as<aod::BigTracksMC>();
      auto daughter2Lc = lcTrack.index2_as<aod::BigTracksMC>();
      auto arrayLcDaughters = array{daughter0Lc, daughter1Lc, daughter2Lc};
      auto arrayDaughters = array{candidate.index1_as<aod::BigTracksMC>(),
	                          daughter0Lc,
				  daughter1Lc,
                                  daughter2Lc};
      /*auto arrayDaughters = array{candidate.index0_as<aod::BigTracksMC>(),
                                  candidate.index1_as<aod::BigTracksMC>(),
				  daughter0Lc,
				  daughter1Lc,
                                  daughter2Lc};*/
      Printf("%d",kProton);
      //Fix this part
      // Λb → Λc+ π-
      //Printf("Checking Λb → Λc+ π-");
      indexRecLc = RecoDecay::getMatchedMCRec(particlesMC, arrayLcDaughters, pdg::Code::kLambdaCPlus, array{2212, 321, -211}, true);
      if (indexRecLc > -1) {
      	indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kLambdaB0, array{-211, 2212, 321, -211}, true, &sign, 2);
	if (indexRec > -1) {
          flag = 1 << hf_cand_lb::DecayType::LbToLcPi;
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = particlesMC.iteratorAt(indexRec);
        origin = (RecoDecay::getMother(particlesMC, particle, 5, true) > -1 ? NonPrompt : Prompt);
      }

      rowMCMatchRec(flag, origin, channel);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      //Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      channel = 0;

      // Λb → Λc+ π-
      //Printf("Checking Λb → Λc+ π-");
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kLambdaB0, array{4122, -211}, true)) {
        // Match  Λc+ --> pKπ-
        std::vector<int> arrDaughter;
        RecoDecay::getDaughters(particlesMC, particle, &arrDaughter, array{4122}, 1);
        auto lcCandMC = particlesMC.iteratorAt(arrDaughter[0]);
        if (RecoDecay::isMatchedMCGen(particlesMC, lcCandMC, pdg::Code::kLambdaCPlus, array{2212, 321, -211}, true)) {
          flag = 1 << hf_cand_lb::DecayType::LbToLcPi;
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        origin = (RecoDecay::getMother(particlesMC, particle, 5, true) > -1 ? NonPrompt : Prompt);
      }

      rowMCMatchGen(flag, origin, channel);
    } // candidate loop
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HFCandidateCreatorLb>(cfgc, TaskName{"hf-cand-creator-lb"}),
    adaptAnalysisTask<HFCandidateCreatorLbExpressions>(cfgc, TaskName{"hf-cand-creator-lb-expressions"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HFCandidateCreatorLbMC>(cfgc, TaskName{"hf-cand-creator-lb-mc"}));
  }
  return workflow;
}
