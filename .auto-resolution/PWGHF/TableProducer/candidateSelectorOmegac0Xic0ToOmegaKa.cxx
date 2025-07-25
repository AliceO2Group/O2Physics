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

/// \file candidateSelectorOmegaKa0Xic0ToOmegaKa.cxx
/// \brief OmegaKa0 Xic0 → Omega Ka selection task
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University
/// \author Ruiqi Yin <ruiqi.yin@cern.ch>, Fudan University

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THnSparse.h>

#include <Rtypes.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
// #include "PWGHF/Core/HfMlResponseOmegaKaToOmegaKa.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

enum PidInfoStored {
  PiFromLam = 0,
  PrFromLam,
  KaFromCasc,
  KaFromCharm
};

/// Struct for applying OmegaKa -> Omega Ka selection cuts
struct HfCandidateSelectorToOmegaKa {
  Produces<aod::HfSelToOmegaKaKf> hfSelToOmegaKaKf;
  // Produces<aod::HfMlSelOmegaKaToOmegaKa> hfMlSelToOmegaKa;

  // LF analysis selections
  Configurable<double> radiusCascMin{"radiusCascMin", 0.5, "Min cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 1.1, "Min V0 radius"};
  Configurable<double> cosPAV0Min{"cosPAV0Min", 0.97, "Min valueCosPA V0"};
  Configurable<double> cosPACascMin{"cosPACascMin", 0.97, "Min value CosPA cascade"};
  Configurable<double> dcaCascDauMax{"dcaCascDauMax", 1.0, "Max DCA cascade daughters"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 1.0, "Max DCA V0 daughters"};
  Configurable<float> dcaBachToPvMin{"dcaBachToPvMin", 0.04, "DCA Bach To PV"};
  Configurable<float> v0MassWindow{"v0MassWindow", 0.01, "V0 mass window"};
  Configurable<float> cascadeMassWindow{"cascadeMassWindow", 0.01, "Cascade mass window"};
  Configurable<bool> applyTrkSelLf{"applyTrkSelLf", true, "Apply track selection for LF daughters"};

  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_omegacxic_to_omega_ka::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_omegacxic_to_omega_ka::Cuts[0], hf_cuts_omegacxic_to_omega_ka::NBinsPt, hf_cuts_omegacxic_to_omega_ka::NCutVars, hf_cuts_omegacxic_to_omega_ka::labelsPt, hf_cuts_omegacxic_to_omega_ka::labelsCutVar}, "OmegaC0 candidate selection per pT bin"};

  // limit charm baryon invariant mass spectrum
  Configurable<double> invMassCharmBaryonMin{"invMassCharmBaryonMin", 2.0, "Lower limit invariant mass spectrum charm baryon"}; // Xic0:2.470 Omegac0:2.695
  Configurable<double> invMassCharmBaryonMax{"invMassCharmBaryonMax", 3.1, "Upper limit invariant mass spectrum charm baryon"};

  // kinematic selections
  Configurable<double> etaTrackCharmBachMax{"etaTrackCharmBachMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 0.8, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptCascMin{"ptCascMin", 0.1, "Min pT kaon <- casc"};
  Configurable<double> ptKaFromCharmBaryonMin{"ptKaFromCharmBaryonMin", 0.2, "Min pT Ka <- charm baryon"};

  Configurable<double> impactParameterXYKaFromCharmBaryonMin{"impactParameterXYKaFromCharmBaryonMin", 0., "Min dcaxy pi from charm baryon track to PV"};
  Configurable<double> impactParameterXYKaFromCharmBaryonMax{"impactParameterXYKaFromCharmBaryonMax", 10., "Max dcaxy pi from charm baryon track to PV"};
  Configurable<double> impactParameterXYCascMin{"impactParameterXYCascMin", 0., "Min dcaxy cascade track to PV"};
  Configurable<double> impactParameterXYCascMax{"impactParameterXYCascMax", 10., "Max dcaxy cascade track to PV"};

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};

  Configurable<double> dcaCharmBaryonDauMax{"dcaCharmBaryonDauMax", 2.0, "Max DCA charm baryon daughters"};

  // PID options
  Configurable<bool> usePidTpcOnly{"usePidTpcOnly", false, "Perform PID using only TPC"};
  Configurable<bool> usePidTpcTofCombined{"usePidTpcTofCombined", true, "Perform PID using TPC & TOF"};

  // PID - TPC selections
  Configurable<double> ptPiPidTpcMin{"ptPiPidTpcMin", -1, "Lower bound of track pT for TPC PID for pion selection"};
  Configurable<double> ptPiPidTpcMax{"ptPiPidTpcMax", 9999.9, "Upper bound of track pT for TPC PID for pion selection"};
  Configurable<double> nSigmaTpcPiMax{"nSigmaTpcPiMax", 3., "Nsigma cut on TPC only for pion selection"};
  Configurable<double> nSigmaTpcCombinedPiMax{"nSigmaTpcCombinedPiMax", 0., "Nsigma cut on TPC combined with TOF for pion selection"};

  Configurable<double> ptPrPidTpcMin{"ptPrPidTpcMin", -1, "Lower bound of track pT for TPC PID for proton selection"};
  Configurable<double> ptPrPidTpcMax{"ptPrPidTpcMax", 9999.9, "Upper bound of track pT for TPC PID for proton selection"};
  Configurable<double> nSigmaTpcPrMax{"nSigmaTpcPrMax", 3., "Nsigma cut on TPC only for proton selection"};
  Configurable<double> nSigmaTpcCombinedPrMax{"nSigmaTpcCombinedPrMax", 0., "Nsigma cut on TPC combined with TOF for proton selection"};

  Configurable<double> ptKaPidTpcMin{"ptKaPidTpcMin", -1, "Lower bound of track pT for TPC PID for kaon selection"};
  Configurable<double> ptKaPidTpcMax{"ptKaPidTpcMax", 9999.9, "Upper bound of track pT for TPC PID for kaon selection"};
  Configurable<double> nSigmaTpcKaMax{"nSigmaTpcKaMax", 3., "Nsigma cut on TPC only for kaon selection"};
  Configurable<double> nSigmaTpcCombinedKaMax{"nSigmaTpcCombinedKaMax", 0., "Nsigma cut on TPC combined with TOF for kaon selection"};

  // PID - TOF selections
  Configurable<double> ptPiPidTofMin{"ptPiPidTofMin", -1, "Lower bound of track pT for TOF PID for pion selection"};
  Configurable<double> ptPiPidTofMax{"ptPiPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for pion selection"};
  Configurable<double> nSigmaTofPiMax{"nSigmaTofPiMax", 3., "Nsigma cut on TOF only for pion selection"};
  Configurable<double> nSigmaTofCombinedPiMax{"nSigmaTofCombinedPiMax", 0., "Nsigma cut on TOF combined with TPC for pion selection"};

  Configurable<double> ptPrPidTofMin{"ptPrPidTofMin", -1, "Lower bound of track pT for TOF PID for proton selection"};
  Configurable<double> ptPrPidTofMax{"ptPrPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for proton selection"};
  Configurable<double> nSigmaTofPrMax{"nSigmaTofPrMax", 3., "Nsigma cut on TOF only for proton selection"};
  Configurable<double> nSigmaTofCombinedPrMax{"nSigmaTofCombinedPrMax", 0., "Nsigma cut on TOF combined with TPC for proton selection"};

  Configurable<double> ptKaPidTofMin{"ptKaPidTofMin", -1, "Lower bound of track pT for TOF PID for kaon selection"};
  Configurable<double> ptKaPidTofMax{"ptKaPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for kaon selection"};
  Configurable<double> nSigmaTofKaMax{"nSigmaTofKaMax", 3., "Nsigma cut on TOF only for kaon selection"};
  Configurable<double> nSigmaTofCombinedKaMax{"nSigmaTofCombinedKaMax", 0., "Nsigma cut on TOF combined with TOF for kaon selection"};

  // detector clusters selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 4, "Maximum value of chi2 fit over TPC clusters"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- charm baryon"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- charm baryon"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 36, "Maximum value of chi2 fit over ITS clusters for pi <- charm baryon"};

  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 2.4, 3.1}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0, 20}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {100, 0, 100}, "Centrality bins"};
  ConfigurableAxis thnConfigAxisPtKaon{"thnConfigAxisPtKaon", {100, 0, 10}, "PtPion from Omegac0 bins"};

  struct : ConfigurableGroup {
    //// KF selection
    std::string prefix = "kfSel";
    Configurable<bool> applyCompetingCascRejection{"applyCompetingCascRejection", false, "Apply competing Xi(for OmegaKa) rejection"};
    Configurable<float> cascadeRejMassWindow{"cascadeRejMassWindow", 0.01, "competing Xi(for OmegaKa) rejection mass window"};
    Configurable<float> v0LdlMin{"v0LdlMin", 3., "Minimum value of l/dl of V0"}; // l/dl and Chi2 are to be determined
    Configurable<float> cascLdlMin{"cascLdlMin", 1., "Minimum value of l/dl of casc"};
    Configurable<float> omegaKaLdlMax{"omegaKaLdlMax", 5., "Maximum value of l/dl of OmegaKa"};
    Configurable<float> cTauOmegaKaMax{"cTauOmegaKaMax", 0.4, "lifetime τ of OmegaKa"};
    Configurable<float> v0Chi2OverNdfMax{"v0Chi2OverNdfMax", 100., "Maximum chi2Geo/NDF of V0"};
    Configurable<float> cascChi2OverNdfMax{"cascChi2OverNdfMax", 100., "Maximum chi2Geo/NDF of casc"};
    Configurable<float> omegaKaChi2OverNdfMax{"omegaKaChi2OverNdfMax", 100., "Maximum chi2Geo/NDF of OmegaKa"};
    Configurable<float> chi2TopoV0ToCascMax{"chi2TopoV0ToCascMax", 100., "Maximum chi2Topo/NDF of V0ToCasc"};
    Configurable<float> chi2TopoKaToCascMax{"chi2TopoKaToCascMax", 100., "Maximum chi2Topo/NDF of KaToCasc"};
    Configurable<float> chi2TopoOmegaKaToPvMax{"chi2TopoOmegaKaToPvMax", 100., "Maximum chi2Topo/NDF of OmegaKaToPv"};
    Configurable<float> chi2TopoCascToOmegaKaMax{"chi2TopoCascToOmegaKaMax", 100., "Maximum chi2Topo/NDF of CascToOmegaKa"};
    Configurable<float> chi2TopoKaToOmegaKaMax{"chi2TopoKaToOmegaKaMax", 100., "Maximum chi2Topo/NDF of KaToOmegaKa"};
    Configurable<float> chi2TopoCascToPvMax{"chi2TopoCascToPvMax", 100., "Maximum chi2Topo/NDF of CascToPv"};
    Configurable<float> chi2TopoKaFromOmegaKaToPvMax{"chi2TopoKaFromOmegaKaToPvMax", 100., "Maximum chi2Topo/NDF of CascToPv"};
    Configurable<float> decayLenOmegaKaMax{"decayLenOmegaKaMax", 1.5, "Maximum decay lengthXY of OmegaKa"};
    Configurable<float> decayLenCascMin{"decayLenCascMin", 1., "Minimum decay lengthXY of Cascade"};
    Configurable<float> decayLenLambdaMin{"decayLenLambdaMin", 0., "Minimum decay lengthXY of V0"};
    Configurable<float> cosPaCascToOmegaKaMin{"cosPaCascToOmegaKaMin", 0.995, "Minimum cosPA of cascade<-OmegaKa"};
    Configurable<float> cosPaV0ToCascMin{"cosPaV0ToCascMin", 0.99, "Minimum cosPA of V0<-cascade"};
  } KfconfigurableGroup;

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;
  TrackSelectorKa selectorKaon;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  using TracksSelLf = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;

  HistogramRegistry registry{"registry"}; // for QA of selections

  OutputObj<TH1D> hInvMassCharmBaryon{TH1D("hInvMassCharmBaryon", "Charm baryon invariant mass;inv mass;entries", 500, 2.0, 3.1)};
  OutputObj<TH1D> hPtCharmBaryon{TH1D("hPtCharmBaryon", "Charm baryon transverse momentum before sel;Pt;entries", 3000, 0., 30)};
  OutputObj<TH1D> hPtKaFromCharmBaryon{TH1D("hPtKaFromCharmBaryon", "Ka from charm baryon transverse momentum before sel;Pt;entries", 2000, 0., 20)};

  void init(InitContext const&)
  {
    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#Omega#Ka) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtKaon{thnConfigAxisPtKaon, "Pt of Kaon from Omegac0."};
    std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt, thnAxisPtKaon};
    registry.add("hMassVsPtVsPtKaon", "Thn for Omegac0 or Xic candidates with InvmassOmegaKa&pT&pTKa", HistType::kTHnSparseF, axes);
    registry.get<THnSparse>(HIST("hMassVsPtVsPtKaon"))->Sumw2();

    selectorPion.setRangePtTpc(ptPiPidTpcMin, ptPiPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcPiMax, nSigmaTpcPiMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedPiMax, nSigmaTpcCombinedPiMax);
    selectorPion.setRangePtTof(ptPiPidTofMin, ptPiPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofPiMax, nSigmaTofPiMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedPiMax, nSigmaTofCombinedPiMax);

    selectorProton.setRangePtTpc(ptPrPidTpcMin, ptPrPidTpcMax);
    selectorProton.setRangeNSigmaTpc(-nSigmaTpcPrMax, nSigmaTpcPrMax);
    selectorProton.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedPrMax, nSigmaTpcCombinedPrMax);
    selectorProton.setRangePtTof(ptPrPidTofMin, ptPrPidTofMax);
    selectorProton.setRangeNSigmaTof(-nSigmaTofPrMax, nSigmaTofPrMax);
    selectorProton.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedPrMax, nSigmaTofCombinedPrMax);

    selectorKaon.setRangePtTpc(ptKaPidTpcMin, ptKaPidTpcMax);
    selectorKaon.setRangeNSigmaTpc(-nSigmaTpcKaMax, nSigmaTpcKaMax);
    selectorKaon.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedKaMax, nSigmaTpcCombinedKaMax);
    selectorKaon.setRangePtTof(ptKaPidTofMin, ptKaPidTofMax);
    selectorKaon.setRangeNSigmaTof(-nSigmaTofKaMax, nSigmaTofKaMax);
    selectorKaon.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedKaMax, nSigmaTofCombinedKaMax);

    const AxisSpec axisSel{2, -0.5, 1.5, "status"};

    registry.add("hSelPID", "hSelPID;status;entries", {HistType::kTH1D, {{12, 0., 12.}}});
    registry.add("hStatusCheck", "Check consecutive selections status;status;entries", {HistType::kTH1D, {{12, 0., 12.}}});

    // for QA of the selections (bin 0 -> candidates that did not pass the selection, bin 1 -> candidates that passed the selection)
    registry.add("hSelSignDec", "hSelSignDec;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelEtaPosV0Dau", "hSelEtaPosV0Dau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelEtaNegV0Dau", "hSelEtaNegV0Dau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelEtaKaFromCasc", "hSelEtaKaFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelEtaKaFromCharm", "hSelEtaKaFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelRadCasc", "hSelRadCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelRadV0", "hSelRadV0;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelCosPACasc", "hSelCosPACasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelCosPAV0", "hSelCosPAV0;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCACascDau", "hSelDCACascDau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAV0Dau", "hSelDCAV0Dau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCACharmDau", "hSelDCACharmDau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAXYPrimPi", "hSelDCAXYPrimPi;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAZPrimPi", "hSelDCAZPrimPi;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAXYCasc", "hSelDCAXYCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelKfPtOmega", "hSelKfPtOmega;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelPtKaFromCharm", "hSelPtKaFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityKaFromCharm", "hSelTPCQualityKaFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityPiFromLam", "hSelTPCQualityPiFromLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityPrFromLam", "hSelTPCQualityPrFromLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityKaFromCasc", "hSelTPCQualityKaFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelITSQualityKaFromCharm", "hSelITSQualityKaFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMassLam", "hSelMassLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMassCasc", "hSelMassCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMassCharmBaryon", "hSelMassCharmBaryon;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDcaXYToPvKaFromCasc", "hSelDcaXYToPvKaFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelPtOmegaKa", "hSelPtOmegaKa;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelCompetingCasc", "hSelCompetingCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelV0_Casc_OmegaKaldl", "hSelV0_Casc_OmegaKaldl;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelctauOmegaKa", "hSelctauOmegaKa;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelChi2GeooverNDFV0_Casc_OmegaKa", "hSelChi2GeooverNDFV0_Casc_OmegaKa;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelChi2TopooverNDFV0_Casc_OmegaKa", "hSelChi2TopooverNDFV0_Casc_OmegaKa;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSeldecayLenOmegaKa_Casc_V0", "hSeldecayLenOmegaKa_Casc_V0;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelcosPaCascToOmegaKa_V0ToCasc", "hSelcosPaCascToOmegaKa_V0ToCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hInvMassXiMinus_rej_cut", "hInvMassXiMinus_rej_cut", kTH1D, {{1000, 1.25f, 1.65f}});
  }
  // for pT-dependent cuts (other selections will move into this in futrue)
  // \param hfCandOmegaKa is candidate
  // return true if candidate passes all cuts
  template <typename T1>
  bool selectionTopol(const T1& hfCandOmegaKa)
  {
    auto candpT = hfCandOmegaKa.kfPtOmegaKa();
    auto KaPtFromOmegaKa = hfCandOmegaKa.kfPtKaFromOmegaKa();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT <= ptCandMin || candpT >= ptCandMax) {
      registry.fill(HIST("hSelPtOmegaKa"), 0);
      return false;
    } else {
      registry.fill(HIST("hSelPtOmegaKa"), 1);
    }

    // check that the candidate pT is within the analysis range
    if (KaPtFromOmegaKa < cuts->get(pTBin, "pT ka from OmegaKa")) {
      registry.fill(HIST("hSelPtKaFromCharm"), 0);
      return false;
    } else {
      registry.fill(HIST("hSelPtKaFromCharm"), 1);
    }

    return true;
  } // end template

  void process(aod::HfCandToOmegaKaKf const& candidates,
               TracksSel const& tracks,
               TracksSelLf const& lfTracks)
  {
    // looping over charm baryon candidates
    for (const auto& candidate : candidates) {
      // initializing selection flags
      bool statusPidLambda = false;
      bool statusPidCascade = false;
      bool statusPidCharmBaryon = false;

      bool statusInvMassLambda = false;
      bool statusInvMassCascade = false;
      bool statusInvMassCharmBaryon = false;

      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise

      int infoTpcStored = 0;
      int infoTofStored = 0;

      auto trackV0PosDauId = candidate.posTrackId();                   // positive V0 daughter
      auto trackV0NegDauId = candidate.negTrackId();                   // negative V0 daughter
      auto trackKaFromCascId = candidate.bachelorId();                 // kaon <- cascade
      auto trackKaFromCharmId = candidate.bachelorFromCharmBaryonId(); // pion <- charm baryon
      auto trackV0PosDau = lfTracks.rawIteratorAt(trackV0PosDauId);
      auto trackV0NegDau = lfTracks.rawIteratorAt(trackV0NegDauId);
      auto trackKaFromCasc = lfTracks.rawIteratorAt(trackKaFromCascId);
      auto trackKaFromCharm = tracks.rawIteratorAt(trackKaFromCharmId);

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int8_t signDecay = candidate.signDecay(); // sign of pi <- cascade

      if (signDecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
        registry.fill(HIST("hSelSignDec"), 1); // anti-particle decay
      } else if (signDecay < 0) {
        registry.fill(HIST("hSelSignDec"), 0); // particle decay
      }

      // pt-dependent selection
      if (!selectionTopol(candidate)) {
        resultSelections = false;
      }

      // eta selection
      double etaV0DauPr = candidate.etaV0DauPr();
      double etaV0DauPi = candidate.etaV0DauPi();
      double etaKaFromCasc = candidate.etaBachFromCasc();
      double etaKaFromCharmBaryon = candidate.etaBachFromCharmBaryon();
      if (std::abs(etaV0DauPr) > etaTrackLFDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaPosV0Dau"), 0);
      } else {
        registry.fill(HIST("hSelEtaPosV0Dau"), 1);
      }
      if (std::abs(etaV0DauPi) > etaTrackLFDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaNegV0Dau"), 0);
      } else {
        registry.fill(HIST("hSelEtaNegV0Dau"), 1);
      }
      if (std::abs(etaKaFromCasc) > etaTrackLFDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaKaFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelEtaKaFromCasc"), 1);
      }
      if (std::abs(etaKaFromCharmBaryon) > etaTrackCharmBachMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaKaFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelEtaKaFromCharm"), 1);
      }

      // minimum radius cut (LFcut)
      if (RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxCascadeKf(), candidate.yDecayVtxCascadeKf()) < radiusCascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelRadCasc"), 0);
      } else {
        registry.fill(HIST("hSelRadCasc"), 1);
      }
      if (RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxV0Kf(), candidate.yDecayVtxV0Kf()) < radiusV0Min) {
        resultSelections = false;
        registry.fill(HIST("hSelRadV0"), 0);
      } else {
        registry.fill(HIST("hSelRadV0"), 1);
      }

      // cosPA (LFcut)
      if (candidate.cosPACasc() < cosPACascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelCosPACasc"), 0);
      } else {
        registry.fill(HIST("hSelCosPACasc"), 1);
      }
      if (candidate.cosPAV0() < cosPAV0Min) {
        resultSelections = false;
        registry.fill(HIST("hSelCosPAV0"), 0);
      } else {
        registry.fill(HIST("hSelCosPAV0"), 1);
      }

      // cascade and v0 daughters dca cut (LF cut)
      if (candidate.dcaCascDau() > dcaCascDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelDCACascDau"), 0);
      } else {
        registry.fill(HIST("hSelDCACascDau"), 1);
      }

      if (candidate.dcaV0Dau() > dcaV0DauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAV0Dau"), 0);
      } else {
        registry.fill(HIST("hSelDCAV0Dau"), 1);
      }

      // dca charm baryon daughters cut
      if (candidate.dcaCharmBaryonDau() > dcaCharmBaryonDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelDCACharmDau"), 0);
      } else {
        registry.fill(HIST("hSelDCACharmDau"), 1);
      }

      // cut on charm bachelor Kaon dcaXY
      if ((std::abs(candidate.impactParBachFromCharmBaryonXY()) < impactParameterXYKaFromCharmBaryonMin) || (std::abs(candidate.impactParBachFromCharmBaryonXY()) > impactParameterXYKaFromCharmBaryonMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAXYPrimPi"), 0);
      } else {
        registry.fill(HIST("hSelDCAXYPrimPi"), 1);
      }

      // cut on cascade dcaXY
      if ((std::abs(candidate.impactParCascXY()) < impactParameterXYCascMin) || (std::abs(candidate.impactParCascXY()) > impactParameterXYCascMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAXYCasc"), 0);
      } else {
        registry.fill(HIST("hSelDCAXYCasc"), 1);
      }

      // Charm daughter pT selections
      if (std::abs(candidate.kfPtOmega()) < ptCascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelKfPtOmega"), 0);
      } else {
        registry.fill(HIST("hSelKfPtOmega"), 1);
      }
      if (std::abs(candidate.kfPtKaFromOmegaKa()) < ptKaFromCharmBaryonMin) {
        resultSelections = false;
      }

      //  Competing Ξ rejection(KF)  Try to reject cases in which the candidate has a an inv. mass compatibler to Xi (bachelor pion) instead of Omega (bachelor kaon)
      if (KfconfigurableGroup.applyCompetingCascRejection) {
        if (std::abs(candidate.invMassCascadeRej() - o2::constants::physics::MassXiMinus) < KfconfigurableGroup.cascadeRejMassWindow) {
          resultSelections = false;
          registry.fill(HIST("hSelCompetingCasc"), 0);
        } else {
          registry.fill(HIST("hSelCompetingCasc"), 1);
          registry.fill(HIST("hInvMassXiMinus_rej_cut"), candidate.invMassCascadeRej());
        }
      }

      // v0&Casc&OmegaKa ldl selection
      if ((candidate.v0ldl() < KfconfigurableGroup.v0LdlMin) || (candidate.cascldl() < KfconfigurableGroup.cascLdlMin) || (candidate.omegaKaldl() > KfconfigurableGroup.omegaKaLdlMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelV0_Casc_OmegaKaldl"), 0);
      } else {
        registry.fill(HIST("hSelV0_Casc_OmegaKaldl"), 1);
      }

      // OmegaKa ctau selsection
      if (candidate.cTauOmegaKa() > KfconfigurableGroup.cTauOmegaKaMax) {
        resultSelections = false;
        registry.fill(HIST("hSelctauOmegaKa"), 0);
      } else {
        registry.fill(HIST("hSelctauOmegaKa"), 1);
      }

      // Chi2Geo/NDF V0&Casc&OmegaKa selection
      if ((candidate.chi2GeoV0() > KfconfigurableGroup.v0Chi2OverNdfMax) || (candidate.chi2GeoV0() < 0) || (candidate.chi2GeoCasc() > KfconfigurableGroup.cascChi2OverNdfMax) || (candidate.chi2GeoCasc() < 0) || (candidate.chi2GeoOmegaKa() > KfconfigurableGroup.omegaKaChi2OverNdfMax) || (candidate.chi2GeoOmegaKa() < 0)) {
        resultSelections = false;
        registry.fill(HIST("hSelChi2GeooverNDFV0_Casc_OmegaKa"), 0);
      } else {
        registry.fill(HIST("hSelChi2GeooverNDFV0_Casc_OmegaKa"), 1);
      }

      // Chi2Topo/NDF selection
      if ((candidate.chi2TopoV0ToCasc() > KfconfigurableGroup.chi2TopoV0ToCascMax) || (candidate.chi2TopoV0ToCasc() < 0) || (candidate.chi2TopoKaToCasc() > KfconfigurableGroup.chi2TopoKaToCascMax) || (candidate.chi2TopoKaToCasc() < 0) || (candidate.chi2TopoCascToOmegaKa() > KfconfigurableGroup.chi2TopoCascToOmegaKaMax) || (candidate.chi2TopoCascToOmegaKa() < 0) || (candidate.chi2TopoKaToOmegaKa() > KfconfigurableGroup.chi2TopoKaToOmegaKaMax) || (candidate.chi2TopoKaToOmegaKa() < 0) ||
          (candidate.chi2TopoOmegaKaToPv() > KfconfigurableGroup.chi2TopoOmegaKaToPvMax) || (candidate.chi2TopoOmegaKaToPv() < 0) || (candidate.chi2TopoCascToPv() > KfconfigurableGroup.chi2TopoCascToPvMax) || (candidate.chi2TopoCascToPv() < 0) || (candidate.chi2TopoKaFromOmegaKaToPv() > KfconfigurableGroup.chi2TopoKaFromOmegaKaToPvMax) || (candidate.chi2TopoKaFromOmegaKaToPv() < 0)) {
        resultSelections = false;
        registry.fill(HIST("hSelChi2TopooverNDFV0_Casc_OmegaKa"), 0);
      } else {
        registry.fill(HIST("hSelChi2TopooverNDFV0_Casc_OmegaKa"), 1);
      }

      // DecaylengthXY of OmegaKa&Casc&V0 selection
      if ((std::abs(candidate.decLenCharmBaryon()) > KfconfigurableGroup.decayLenOmegaKaMax) || (std::abs(candidate.decLenCascade()) < KfconfigurableGroup.decayLenCascMin) || (std::abs(candidate.decLenV0()) < KfconfigurableGroup.decayLenLambdaMin)) {
        resultSelections = false;
        registry.fill(HIST("hSeldecayLenOmegaKa_Casc_V0"), 0);
      } else {
        registry.fill(HIST("hSeldecayLenOmegaKa_Casc_V0"), 1);
      }

      // KFPA cut cosPaCascToOmegaKa cosPaV0ToCasc
      if ((candidate.cosPaCascToOmegaKa() < KfconfigurableGroup.cosPaCascToOmegaKaMin) || (candidate.cosPaV0ToCasc() < KfconfigurableGroup.cosPaV0ToCascMin)) {
        resultSelections = false;
        registry.fill(HIST("hSelcosPaCascToOmegaKa_V0ToCasc"), 0);
      } else {
        registry.fill(HIST("hSelcosPaCascToOmegaKa_V0ToCasc"), 1);
      }

      //  TPC clusters selections
      if (applyTrkSelLf) {
        if (!isSelectedTrackTpcQuality(trackPiFromLam, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
          registry.fill(HIST("hSelTPCQualityPiFromLam"), 0);
        } else {
          registry.fill(HIST("hSelTPCQualityPiFromLam"), 1);
        }
        if (!isSelectedTrackTpcQuality(trackPrFromLam, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
          registry.fill(HIST("hSelTPCQualityPrFromLam"), 0);
        } else {
          registry.fill(HIST("hSelTPCQualityPrFromLam"), 1);
        }
        if (!isSelectedTrackTpcQuality(trackKaFromCasc, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
          registry.fill(HIST("hSelTPCQualityKaFromCasc"), 0);
        } else {
          registry.fill(HIST("hSelTPCQualityKaFromCasc"), 1);
        }
      }
      if (!isSelectedTrackTpcQuality(trackKaFromCharm, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelTPCQualityKaFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelTPCQualityKaFromCharm"), 1);
      }

      //  ITS clusters selection
      if (!isSelectedTrackItsQuality(trackKaFromCharm, nClustersItsMin, itsChi2PerClusterMax) || trackKaFromCharm.itsNClsInnerBarrel() < nClustersItsInnBarrMin) {
        resultSelections = false;
        registry.fill(HIST("hSelITSQualityKaFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelITSQualityKaFromCharm"), 1);
      }

      // track-level PID selection

      // for TrackSelectorPID
      int statusPidPrFromLam = -999;
      int statusPidPiFromLam = -999;
      int statusPidKaFromCasc = -999;
      int statusPidKaFromCharmBaryon = -999;

      if (usePidTpcOnly && usePidTpcTofCombined) {
        LOGF(fatal, "Check the PID configurables, usePidTpcOnly and usePidTpcTofCombined can't have the same value");
      } else if (!usePidTpcOnly && !usePidTpcTofCombined) {
        LOGF(fatal, "At least one PID method must be enabled");
      }

      if (trackPiFromLam.hasTPC()) {
        SETBIT(infoTpcStored, PiFromLam);
      }
      if (trackPrFromLam.hasTPC()) {
        SETBIT(infoTpcStored, PrFromLam);
      }
      if (trackKaFromCasc.hasTPC()) {
        SETBIT(infoTpcStored, KaFromCasc);
      }
      if (trackKaFromCharm.hasTPC()) {
        SETBIT(infoTpcStored, KaFromCharm);
      }
      if (trackPiFromLam.hasTOF()) {
        SETBIT(infoTofStored, PiFromLam);
      }
      if (trackPrFromLam.hasTOF()) {
        SETBIT(infoTofStored, PrFromLam);
      }
      if (trackKaFromCasc.hasTOF()) {
        SETBIT(infoTofStored, KaFromCasc);
      }
      if (trackKaFromCharm.hasTOF()) {
        SETBIT(infoTofStored, KaFromCharm);
      }

      if (usePidTpcOnly) {
        statusPidPrFromLam = selectorProton.statusTpc(trackPrFromLam);
        statusPidPiFromLam = selectorPion.statusTpc(trackPiFromLam);
        statusPidKaFromCasc = selectorKaon.statusTpc(trackKaFromCasc);
        statusPidKaFromCharmBaryon = selectorKaon.statusTpc(trackKaFromCharm);
      } else if (usePidTpcTofCombined) {
        statusPidPrFromLam = selectorProton.statusTpcOrTof(trackPrFromLam);
        statusPidPiFromLam = selectorPion.statusTpcOrTof(trackPiFromLam);
        statusPidKaFromCasc = selectorKaon.statusTpcOrTof(trackKaFromCasc);
        statusPidKaFromCharmBaryon = selectorKaon.statusTpcOrTof(trackKaFromCharm);
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted) {
        statusPidLambda = true;
        if (resultSelections) {
          registry.fill(HIST("hStatusCheck"), 0.5);
        }
      } else {
        resultSelections = false;
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted && statusPidKaFromCasc == TrackSelectorPID::Accepted) {
        statusPidCascade = true;
        if (resultSelections) {
          registry.fill(HIST("hStatusCheck"), 1.5);
        }
      } else {
        resultSelections = false;
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted && statusPidKaFromCasc == TrackSelectorPID::Accepted && statusPidKaFromCharmBaryon == TrackSelectorPID::Accepted) {
        statusPidCharmBaryon = true;
        if (resultSelections) {
          registry.fill(HIST("hStatusCheck"), 2.5);
        }
      } else {
        resultSelections = false;
      }

      // invariant mass cuts
      double invMassLambda = candidate.invMassLambda();
      double invMassCascade = candidate.invMassCascade();
      double invMassCharmBaryon = candidate.invMassCharmBaryon();

      if (std::abs(invMassLambda - o2::constants::physics::MassLambda0) < v0MassWindow) {
        statusInvMassLambda = true;
        registry.fill(HIST("hSelMassLam"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 3.5);
        }
      } else {
        registry.fill(HIST("hSelMassLam"), 0);
        resultSelections = false;
      }

      if (std::abs(invMassCascade - o2::constants::physics::MassOmegaMinus) < cascadeMassWindow) {
        statusInvMassCascade = true;
        registry.fill(HIST("hSelMassCasc"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 4.5);
        }
      } else {
        registry.fill(HIST("hSelMassCasc"), 0);
        resultSelections = false;
      }

      if ((invMassCharmBaryon >= invMassCharmBaryonMin) && (invMassCharmBaryon <= invMassCharmBaryonMax)) {
        statusInvMassCharmBaryon = true;
        registry.fill(HIST("hSelMassCharmBaryon"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && statusInvMassCascade && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 5.5);
        }
      } else {
        registry.fill(HIST("hSelMassCharmBaryon"), 0);
        resultSelections = false;
      }

      hfSelToOmegaKaKf(statusPidLambda, statusPidCascade, statusPidCharmBaryon, statusInvMassLambda, statusInvMassCascade, statusInvMassCharmBaryon, resultSelections, infoTpcStored, infoTofStored,
                       trackKaFromCharm.tpcNSigmaKa(), trackKaFromCasc.tpcNSigmaKa(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                       trackKaFromCharm.tofNSigmaKa(), trackKaFromCasc.tofNSigmaKa(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());

      if (resultSelections) {
        if (!statusPidLambda) {
          registry.fill(HIST("hSelPID"), 0.5);
        }
        if (statusPidLambda) {
          registry.fill(HIST("hSelPID"), 1.5);
        }
        if (!statusPidCascade) {
          registry.fill(HIST("hSelPID"), 2.5);
        }
        if (statusPidCascade) {
          registry.fill(HIST("hSelPID"), 3.5);
        }
        if (!statusPidCharmBaryon) {
          registry.fill(HIST("hSelPID"), 4.5);
        }
        if (statusPidCharmBaryon) {
          registry.fill(HIST("hSelPID"), 5.5);
        }
        if (!statusInvMassLambda) {
          registry.fill(HIST("hSelPID"), 6.5);
        }
        if (statusInvMassLambda) {
          registry.fill(HIST("hSelPID"), 7.5);
        }
        if (!statusInvMassCascade) {
          registry.fill(HIST("hSelPID"), 8.5);
        }
        if (statusInvMassCascade) {
          registry.fill(HIST("hSelPID"), 9.5);
        }
        if (!statusInvMassCharmBaryon) {
          registry.fill(HIST("hSelPID"), 10.5);
        }
        if (statusInvMassCharmBaryon) {
          registry.fill(HIST("hSelPID"), 11.5);
        }
      }

      if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && statusInvMassCascade && statusInvMassCharmBaryon && resultSelections) {
        hInvMassCharmBaryon->Fill(invMassCharmBaryon);
        hPtCharmBaryon->Fill(candidate.kfPtOmegaKa());
        hPtKaFromCharmBaryon->Fill(candidate.kfPtKaFromOmegaKa());
        registry.fill(HIST("hMassVsPtVsPtKaon"),
                      candidate.invMassCharmBaryon(),
                      candidate.kfPtOmegaKa(),
                      candidate.kfPtKaFromOmegaKa());
      }
    }
  } // end process
}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorToOmegaKa>(cfgc)};
}
