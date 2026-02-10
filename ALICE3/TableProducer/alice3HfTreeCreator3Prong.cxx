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

/// \file alice3HfTreeCreator3Prong.cxx
/// \brief Writer of 3-prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples or ML training.
///
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Turin Polytechnic University and INFN Turin

#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "ALICE3/DataModel/OTFPIDTrk.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/RICH.h"
#include "ALICE3/Utils/utilsHfAlice3.h"
#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);                     //! Radius of secondary vertex (cm)
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);                                     //! Transverse momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(PProng0, pProng0, float);                                       //! Momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float); //! Normalised impact parameter of prong0
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);                                     //! Transverse momentum of prong1 (GeV/c)
DECLARE_SOA_COLUMN(PProng1, pProng1, float);                                       //! Momentum of prong1 (in GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float); //! Normalised impact parameter of prong1
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);                                     //! Transverse momentum of prong2 (GeV/c)
DECLARE_SOA_COLUMN(PProng2, pProng2, float);                                       //! Momentum of prong2 (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, float); //! Normalised impact parameter of prong2
DECLARE_SOA_COLUMN(M, m, float);                                                   //! Invariant mass of cand (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                 //! Transverse momentum of cand (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                   //! Momentum of cand (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                   //! Rapidity of cand
DECLARE_SOA_COLUMN(Eta, eta, float);                                               //! Pseudorapidity of cand
DECLARE_SOA_COLUMN(Phi, phi, float);                                               //! Azimuth angle of cand
DECLARE_SOA_COLUMN(E, e, float);                                                   //! Energy of cand (GeV)
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                               //! Decay length of cand (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                           //! Transverse decay length of cand (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);           //! Normalised decay length of cand
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);       //! Normalised transverse decay length of cand
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                               //! Cosine pointing angle of cand
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                           //! Cosine pointing angle of cand in transverse plane
// PID columns
// Prong 0
DECLARE_SOA_COLUMN(NSigTrkPi0, nSigTrkPi0, float);       //!
DECLARE_SOA_COLUMN(NSigRichPi0, nSigRichPi0, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofPi0, nSigInnTofPi0, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPi0, nSigOutTofPi0, float); //!
DECLARE_SOA_COLUMN(NSigTrkKa0, nSigTrkKa0, float);       //!
DECLARE_SOA_COLUMN(NSigRichKa0, nSigRichKa0, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofKa0, nSigInnTofKa0, float); //!
DECLARE_SOA_COLUMN(NSigOutTofKa0, nSigOutTofKa0, float); //!
DECLARE_SOA_COLUMN(NSigTrkPr0, nSigTrkPr0, float);       //!
DECLARE_SOA_COLUMN(NSigRichPr0, nSigRichPr0, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofPr0, nSigInnTofPr0, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPr0, nSigOutTofPr0, float); //!
// Prong 1
DECLARE_SOA_COLUMN(NSigTrkPi1, nSigTrkPi1, float);       //!
DECLARE_SOA_COLUMN(NSigRichPi1, nSigRichPi1, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofPi1, nSigInnTofPi1, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPi1, nSigOutTofPi1, float); //!
DECLARE_SOA_COLUMN(NSigTrkKa1, nSigTrkKa1, float);       //!
DECLARE_SOA_COLUMN(NSigRichKa1, nSigRichKa1, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofKa1, nSigInnTofKa1, float); //!
DECLARE_SOA_COLUMN(NSigOutTofKa1, nSigOutTofKa1, float); //!
DECLARE_SOA_COLUMN(NSigTrkPr1, nSigTrkPr1, float);       //!
DECLARE_SOA_COLUMN(NSigRichPr1, nSigRichPr1, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofPr1, nSigInnTofPr1, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPr1, nSigOutTofPr1, float); //!
// Prong 2
DECLARE_SOA_COLUMN(NSigTrkPi2, nSigTrkPi2, float);       //!
DECLARE_SOA_COLUMN(NSigRichPi2, nSigRichPi2, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofPi2, nSigInnTofPi2, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPi2, nSigOutTofPi2, float); //!
DECLARE_SOA_COLUMN(NSigTrkKa2, nSigTrkKa2, float);       //!
DECLARE_SOA_COLUMN(NSigRichKa2, nSigRichKa2, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofKa2, nSigInnTofKa2, float); //!
DECLARE_SOA_COLUMN(NSigOutTofKa2, nSigOutTofKa2, float); //!
DECLARE_SOA_COLUMN(NSigTrkPr2, nSigTrkPr2, float);       //!
DECLARE_SOA_COLUMN(NSigRichPr2, nSigRichPr2, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofPr2, nSigInnTofPr2, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPr2, nSigOutTofPr2, float); //!
// ML scores
DECLARE_SOA_COLUMN(MlScore0, mlScore0, float); //! ML score of the first configured index
DECLARE_SOA_COLUMN(MlScore1, mlScore1, float); //! ML score of the second configured index
DECLARE_SOA_COLUMN(MlScore2, mlScore2, float); //! ML score of the third configured index
} // namespace full

// Topology tables
DECLARE_SOA_TABLE(Alice3CandVtxs, "AOD", "ALICE3CANDVTX", //!
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  a3_hf_cand::XSecondaryVertex,
                  a3_hf_cand::YSecondaryVertex,
                  a3_hf_cand::ZSecondaryVertex,
                  full::RSecondaryVertex);

DECLARE_SOA_TABLE(Alice3CandTopos, "AOD", "ALICE3CANDTOPO", //!
                  a3_hf_cand::Chi2PCA,
                  full::DecayLength,
                  a3_hf_cand::ErrorDecayLength,
                  full::DecayLengthXY,
                  a3_hf_cand::ErrorDecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY);

DECLARE_SOA_TABLE(Alice3DaugTopos, "AOD", "ALICE3DAUGTOPO", //!
                  full::ImpactParameterNormalised0,
                  full::ImpactParameterNormalised1,
                  full::ImpactParameterNormalised2,
                  a3_hf_cand::ImpactParameterY0,
                  a3_hf_cand::ImpactParameterY1,
                  a3_hf_cand::ImpactParameterY2,
                  a3_hf_cand::ErrorImpactParameterY0,
                  a3_hf_cand::ErrorImpactParameterY1,
                  a3_hf_cand::ErrorImpactParameterY2,
                  a3_hf_cand::ImpactParameterZ0,
                  a3_hf_cand::ImpactParameterZ1,
                  a3_hf_cand::ImpactParameterZ2,
                  a3_hf_cand::ErrorImpactParameterZ0,
                  a3_hf_cand::ErrorImpactParameterZ1,
                  a3_hf_cand::ErrorImpactParameterZ2);

// ML tables
DECLARE_SOA_TABLE(Alice3PMls, "AOD", "ALICE3PML",
                  full::MlScore0,
                  full::MlScore1,
                  full::MlScore2)
// Kinematics information table
DECLARE_SOA_TABLE(Alice3Kine3Ps, "AOD", "ALICE3KINE3P",
                  full::PtProng0,
                  full::PtProng1,
                  full::PtProng2,
                  full::Pt,
                  full::P,
                  full::Eta,
                  full::Phi,
                  full::M,
                  full::Y,
                  full::E);

DECLARE_SOA_TABLE(Alice3Cand3PGens, "AOD", "ALICE3CAND3PGEN",
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y);
DECLARE_SOA_TABLE(Alice3DaugMoms, "AOD", "ALICE3DAUGMOM", //!
                  a3_hf_cand::PxProng0, a3_hf_cand::PyProng0, a3_hf_cand::PzProng0, full::PProng0,
                  a3_hf_cand::PxProng1, a3_hf_cand::PyProng1, a3_hf_cand::PzProng1, full::PProng1,
                  a3_hf_cand::PxProng2, a3_hf_cand::PyProng2, a3_hf_cand::PzProng2, full::PProng2);
// MC matching tables
DECLARE_SOA_TABLE(Alice3McRecs, "AOD", "ALICE3MCREC",
                  a3_mc_truth::FlagMcRec,
                  a3_mc_truth::OriginMcRec);

DECLARE_SOA_TABLE(Alice3McGens, "AOD", "ALICE3MCGEN",
                  a3_mc_truth::FlagMcGen,
                  a3_mc_truth::OriginMcGen);
// PID tables
DECLARE_SOA_TABLE(Alice3PidPi0s, "AOD", "ALICE3PIDPI0",
                  full::NSigTrkPi0,
                  full::NSigRichPi0,
                  full::NSigInnTofPi0,
                  full::NSigOutTofPi0)
DECLARE_SOA_TABLE(Alice3PidPi1s, "AOD", "ALICE3PIDPI1",
                  full::NSigTrkPi1,
                  full::NSigRichPi1,
                  full::NSigInnTofPi1,
                  full::NSigOutTofPi1)
DECLARE_SOA_TABLE(Alice3PidPi2s, "AOD", "ALICE3PIDPI2",
                  full::NSigTrkPi2,
                  full::NSigRichPi2,
                  full::NSigInnTofPi2,
                  full::NSigOutTofPi2)
DECLARE_SOA_TABLE(Alice3PidKa0s, "AOD", "ALICE3PIDKA0",
                  full::NSigTrkKa0,
                  full::NSigRichKa0,
                  full::NSigInnTofKa0,
                  full::NSigOutTofKa0)
DECLARE_SOA_TABLE(Alice3PidKa1s, "AOD", "ALICE3PIDKA1",
                  full::NSigTrkKa1,
                  full::NSigRichKa1,
                  full::NSigInnTofKa1,
                  full::NSigOutTofKa1)
DECLARE_SOA_TABLE(Alice3PidKa2s, "AOD", "ALICE3PIDKA2",
                  full::NSigTrkKa2,
                  full::NSigRichKa2,
                  full::NSigInnTofKa2,
                  full::NSigOutTofKa2)
DECLARE_SOA_TABLE(Alice3PidPr0s, "AOD", "ALICE3PIDPR0",
                  full::NSigTrkPr0,
                  full::NSigRichPr0,
                  full::NSigInnTofPr0,
                  full::NSigOutTofPr0)
DECLARE_SOA_TABLE(Alice3PidPr1s, "AOD", "ALICE3PIDPR1",
                  full::NSigTrkPr1,
                  full::NSigRichPr1,
                  full::NSigInnTofPr1,
                  full::NSigOutTofPr1)
DECLARE_SOA_TABLE(Alice3PidPr2s, "AOD", "ALICE3PIDPR2",
                  full::NSigTrkPr2,
                  full::NSigRichPr2,
                  full::NSigInnTofPr2,
                  full::NSigOutTofPr2)

} // namespace o2::aod

/// Writes the full information in an output TTree
struct Alice3HfTreeCreator3Prong {
  Produces<o2::aod::Alice3CandVtxs> rowCandVtxs;
  Produces<o2::aod::Alice3CandTopos> rowCandTopos;
  Produces<o2::aod::Alice3DaugTopos> rowDaugTopos;
  Produces<o2::aod::Alice3PMls> rowCandMls;
  Produces<o2::aod::Alice3Kine3Ps> rowCandKine3Ps;
  Produces<o2::aod::Alice3Cand3PGens> rowCand3PGen;
  Produces<o2::aod::Alice3DaugMoms> rowDaugMoms;
  Produces<o2::aod::Alice3McRecs> rowCand3PMcMatchRec;
  Produces<o2::aod::Alice3McGens> rowCand3PMcMatchGen;
  Produces<o2::aod::Alice3PidPi0s> rowPidPi0;
  Produces<o2::aod::Alice3PidPi1s> rowPidPi1;
  Produces<o2::aod::Alice3PidPi2s> rowPidPi2;
  Produces<o2::aod::Alice3PidKa0s> rowPidKa0;
  Produces<o2::aod::Alice3PidKa1s> rowPidKa1;
  Produces<o2::aod::Alice3PidKa2s> rowPidKa2;
  Produces<o2::aod::Alice3PidPr0s> rowPidPr0;
  Produces<o2::aod::Alice3PidPr1s> rowPidPr1;
  Produces<o2::aod::Alice3PidPr2s> rowPidPr2;

  // Configurables to fill tables
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> fillCandVtxInfo{"fillCandVtxInfo", false, "fill candidate vtx info"};
    Configurable<bool> fillCandTopoInfo{"fillCandTopoInfo", false, "fill candidate topology info"};
    Configurable<bool> fillDaugTopoInfo{"fillDaugTopoInfo", false, "fill daughter topology info"};
    Configurable<bool> fillMlScoreInfo{"fillMlScoreInfo", false, "fill ML scores info"};
    Configurable<bool> fillCandKineInfo{"fillCandKineInfo", false, "fill candidate kinematic info"};
    Configurable<bool> fillCandGenKineInfo{"fillCandGenKineInfo", false, "fill generated candidate kinematic info"};
    Configurable<bool> fillDaugKineInfo{"fillDaugKineInfo", false, "fill daughter kinematic info"};
    Configurable<bool> fillMcMatchRecoInfo{"fillMcMatchRecoInfo", false, "fill MC match reconstruction info"};
    Configurable<bool> fillMcMatchGenInfo{"fillMcMatchGenInfo", false, "fill MC match generation info"};
    Configurable<bool> fillPid{"fillPid", false, "fill PID info"};
  } fillTables;
  // parameters for production of training samples
  Configurable<bool> fillSwapMassHypo{"fillSwapMassHypo", false, "Flag to fill derived tables with swapped mass hypothesis"};
  Configurable<bool> fillOnlySignal{"fillOnlySignal", true, "Flag to fill derived tables with signal"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background"};
  Configurable<float> downSampleFactor{"downSampleFactor", 1., "Fraction of cands to keep"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelperAlice3 hfHelper;

  using CandsLcRec = soa::Filtered<soa::Join<aod::Alice3Cand3Ps, aod::Alice3Sel3Ps, aod::Alice3PidLcs, aod::Alice3McRecFlags>>;
  using CandsLcRecWMl = soa::Filtered<soa::Join<aod::Alice3Cand3Ps, aod::Alice3Sel3Ps, aod::Alice3PidLcs, aod::Alice3Ml3Ps, aod::Alice3McRecFlags>>;
  using CandsMcGen = soa::Join<aod::McParticles, aod::Alice3McGenFlags>;

  Filter filterSelectCandidates = (aod::a3_hf_sel_3prong::isSelMassHypo0 == true || aod::a3_hf_sel_3prong::isSelMassHypo1 == true);
  Filter filterSelectGenCands = nabs(aod::a3_mc_truth::flagMcGen) == static_cast<int>(CharmHadAlice3::Lc);

  Partition<CandsLcRec> recoLcCandSig = nabs(o2::aod::a3_mc_truth::flagMcRec) == static_cast<int>(CharmHadAlice3::Lc);
  Partition<CandsLcRec> recoLcCandBkg = nabs(o2::aod::a3_mc_truth::flagMcRec) == 0;
  Partition<CandsLcRecWMl> recoLcCandSigWMl = nabs(o2::aod::a3_mc_truth::flagMcRec) == static_cast<int>(CharmHadAlice3::Lc);
  Partition<CandsLcRecWMl> recoLcCandBkgWMl = nabs(o2::aod::a3_mc_truth::flagMcRec) == 0;

  void init(InitContext const&)
  {
    const std::array<bool, 2> doprocess{doprocessLc, doprocessLcWMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }
  }

  /// Reserve space in the output tables
  /// \tparam TCand Type of candidate
  /// \param nCands Number of candidates to reserve space for
  /// \param cand Candidate to be used to check which PID tables to reserve
  template <typename TCand>
  void reserveTables(size_t nCands, const TCand& cand)
  {
    if (fillTables.fillCandVtxInfo)
      rowCandVtxs.reserve(nCands);
    if (fillTables.fillCandTopoInfo)
      rowCandTopos.reserve(nCands);
    if (fillTables.fillDaugTopoInfo)
      rowDaugTopos.reserve(nCands);
    if (fillTables.fillMlScoreInfo)
      rowCandMls.reserve(nCands);
    if (fillTables.fillCandKineInfo)
      rowCandKine3Ps.reserve(nCands);
    if (fillTables.fillCandGenKineInfo)
      rowCand3PGen.reserve(nCands);
    if (fillTables.fillDaugKineInfo)
      rowDaugMoms.reserve(nCands);
    if (fillTables.fillMcMatchRecoInfo)
      rowCand3PMcMatchRec.reserve(nCands);
    if (fillTables.fillMcMatchGenInfo)
      rowCand3PMcMatchGen.reserve(nCands);
    // PID tables
    if (fillTables.fillPid) {
      if constexpr (requires { cand.nSigTrkPi0(); })
        rowPidPi0.reserve(nCands);
      if constexpr (requires { cand.nSigTrkPi1(); })
        rowPidPi1.reserve(nCands);
      if constexpr (requires { cand.nSigTrkPi2(); })
        rowPidPi2.reserve(nCands);
      if constexpr (requires { cand.nSigTrkKa0(); })
        rowPidKa0.reserve(nCands);
      if constexpr (requires { cand.nSigTrkKa1(); })
        rowPidKa1.reserve(nCands);
      if constexpr (requires { cand.nSigTrkKa2(); })
        rowPidKa2.reserve(nCands);
      if constexpr (requires { cand.nSigTrkPr0(); })
        rowPidPr0.reserve(nCands);
      if constexpr (requires { cand.nSigTrkPr1(); })
        rowPidPr1.reserve(nCands);
      if constexpr (requires { cand.nSigTrkPr2(); })
        rowPidPr2.reserve(nCands);
    }
  }

  /// Fill reconstructed candidate tables
  /// \tparam CharmHadAlice3: charm hadron type
  /// \tparam IsSwapMassHypo: whether to swap mass hypothesis or not
  /// \tparam T: candidate type
  /// \param cand: candidate to be used to fill the tables
  template <CharmHadAlice3 CharmHad, bool IsSwapMassHypo, typename T>
  void fillRecoTables(const T& cand)
  {
    if (fillTables.fillCandVtxInfo) {
      rowCandVtxs(
        cand.posX(),
        cand.posY(),
        cand.posZ(),
        cand.xSecondaryVertex(),
        cand.ySecondaryVertex(),
        cand.zSecondaryVertex(),
        cand.rSecondaryVertex());
    }
    if (fillTables.fillCandTopoInfo) {
      rowCandTopos(
        cand.chi2PCA(),
        cand.decayLength(),
        cand.errorDecayLength(),
        cand.decayLengthXY(),
        cand.errorDecayLengthXY(),
        cand.decayLengthNormalised(),
        cand.decayLengthXYNormalised(),
        cand.cpa(),
        cand.cpaXY());
    }
    if (fillTables.fillDaugTopoInfo) {
      rowDaugTopos(
        cand.impactParameterNormalised0(),
        cand.impactParameterNormalised1(),
        cand.impactParameterNormalised2(),
        cand.impactParameterY0(),
        cand.impactParameterY1(),
        cand.impactParameterY2(),
        cand.errorImpactParameterY0(),
        cand.errorImpactParameterY1(),
        cand.errorImpactParameterY2(),
        cand.impactParameterZ0(),
        cand.impactParameterZ1(),
        cand.impactParameterZ2(),
        cand.errorImpactParameterZ0(),
        cand.errorImpactParameterZ1(),
        cand.errorImpactParameterZ2());
    }
    if constexpr (requires { cand.mlScore0(); }) {
      if (fillTables.fillMlScoreInfo) {
        rowCandMls(cand.mlScore0(), cand.mlScore1(), cand.mlScore2());
      }
    }
    if (fillTables.fillCandKineInfo) {
      rowCandKine3Ps(
        cand.ptProng0(), cand.ptProng1(), cand.ptProng2(),
        cand.pt(), RecoDecay::p(cand.px(), cand.py(), cand.pz()), cand.eta(), cand.phi(),
        hfHelper.getCandMass<CharmHad, IsSwapMassHypo>(cand), hfHelper.getCandY<CharmHad>(cand), hfHelper.getCandEnergy<CharmHad>(cand));
    }
    if (fillTables.fillDaugKineInfo) {
      rowDaugMoms(
        cand.pxProng0(), cand.pyProng0(), cand.pzProng0(), RecoDecay::p(cand.pxProng0(), cand.pyProng0(), cand.pzProng0()),
        cand.pxProng1(), cand.pyProng1(), cand.pzProng1(), RecoDecay::p(cand.pxProng1(), cand.pyProng1(), cand.pzProng1()),
        cand.pxProng2(), cand.pyProng2(), cand.pzProng2(), RecoDecay::p(cand.pxProng2(), cand.pyProng2(), cand.pzProng2()));
    }
    if (fillTables.fillMcMatchRecoInfo)
      rowCand3PMcMatchRec(cand.flagMcRec(), cand.originMcRec());

    // Fill PID tables
    if (fillTables.fillPid) {
      if constexpr (requires { cand.nSigTrkPi0(); })
        rowPidPi0(cand.nSigTrkPi0(), cand.nSigRichPi0(), cand.nSigInnTofPi0(), cand.nSigOutTofPi0());
      if constexpr (requires { cand.nSigTrkPi1(); })
        rowPidPi1(cand.nSigTrkPi1(), cand.nSigRichPi1(), cand.nSigInnTofPi1(), cand.nSigOutTofPi1());
      if constexpr (requires { cand.nSigTrkPi2(); })
        rowPidPi2(cand.nSigTrkPi2(), cand.nSigRichPi2(), cand.nSigInnTofPi2(), cand.nSigOutTofPi2());
      if constexpr (requires { cand.nSigTrkKa0(); })
        rowPidKa0(cand.nSigTrkKa0(), cand.nSigRichKa0(), cand.nSigInnTofKa0(), cand.nSigOutTofKa0());
      if constexpr (requires { cand.nSigTrkKa1(); })
        rowPidKa1(cand.nSigTrkKa1(), cand.nSigRichKa1(), cand.nSigInnTofKa1(), cand.nSigOutTofKa1());
      if constexpr (requires { cand.nSigTrkKa2(); })
        rowPidKa2(cand.nSigTrkKa2(), cand.nSigRichKa2(), cand.nSigInnTofKa2(), cand.nSigOutTofKa2());
      if constexpr (requires { cand.nSigTrkPr0(); })
        rowPidPr0(cand.nSigTrkPr0(), cand.nSigRichPr0(), cand.nSigInnTofPr0(), cand.nSigOutTofPr0());
      if constexpr (requires { cand.nSigTrkPr1(); })
        rowPidPr1(cand.nSigTrkPr1(), cand.nSigRichPr1(), cand.nSigInnTofPr1(), cand.nSigOutTofPr1());
      if constexpr (requires { cand.nSigTrkPr2(); })
        rowPidPr2(cand.nSigTrkPr2(), cand.nSigRichPr2(), cand.nSigInnTofPr2(), cand.nSigOutTofPr2());
    }
  }

  /// Function to fill generated tables
  /// \tparam CharmHad Type of 3prong particle
  /// \tparam T Type of generated candidates collection
  /// \param parts Generated candidates collection
  template <CharmHadAlice3 CharmHad, typename T>
  void fillGenTables(const T& parts)
  {
    for (const auto& part : parts) {
      if (fillTables.fillCandGenKineInfo) {
        rowCand3PGen(part.pt(), part.eta(), part.phi(), hfHelper.getCandY<CharmHad>(part));
      }
      if (fillTables.fillMcMatchGenInfo) {
        rowCand3PMcMatchGen(part.flagMcGen(), part.originMcGen());
      }
    }
  }

  /// Function to fill both reco and gen tables
  /// from any candidate collection
  /// \tparam CharmHad Type of 3prong particle
  /// \tparam TCandsRec Type of reconstructed candidates collection
  /// \tparam TCandsGen Type of generated candidates collection
  /// \param candsRec Reconstructed candidates collection
  /// \param candsGen Generated candidates collection
  template <CharmHadAlice3 CharmHad, typename TCandsRec, typename TCandsGen>
  void fillRecoGenTables(const TCandsRec& candsRec,
                         const TCandsGen& candsGen)
  {
    reserveTables(candsRec.size(), *candsRec.begin());
    for (const auto& cand : candsRec) {
      if (downSampleFactor < 1.) {
        float const pseudoRndm = cand.ptProng0() * 1000. - static_cast<int64_t>(cand.ptProng0() * 1000);
        if (cand.pt() < ptMaxForDownSample && pseudoRndm >= downSampleFactor) {
          continue;
        }
      }
      if (cand.isSelMassHypo0()) {
        fillRecoTables<CharmHad, false>(cand);
      }
      if (fillSwapMassHypo) {
        if (cand.isSelMassHypo1()) {
          fillRecoTables<CharmHad, true>(cand);
        }
      }
    }
    fillGenTables<CharmHad>(candsGen);
  }

  void processLc(CandsLcRec const& selCandsLcRec,
                 CandsMcGen const& parts,
                 aod::McCollisions const&)
  {
    if (fillOnlySignal) {
      fillRecoGenTables<CharmHadAlice3::Lc>(recoLcCandSig, parts);
    } else if (fillOnlyBackground) {
      fillRecoGenTables<CharmHadAlice3::Lc>(recoLcCandBkg, parts);
    } else {
      fillRecoGenTables<CharmHadAlice3::Lc>(selCandsLcRec, parts);
    }
  }
  PROCESS_SWITCH(Alice3HfTreeCreator3Prong, processLc, "Process Lc", true);

  void processLcWMl(CandsLcRecWMl const& selCandsLcRec,
                    CandsMcGen const& parts,
                    aod::McCollisions const&)
  {
    if (fillOnlySignal) {
      fillRecoGenTables<CharmHadAlice3::Lc>(recoLcCandSig, parts);
    } else if (fillOnlyBackground) {
      fillRecoGenTables<CharmHadAlice3::Lc>(recoLcCandBkg, parts);
    } else {
      fillRecoGenTables<CharmHadAlice3::Lc>(selCandsLcRec, parts);
    }
  }
  PROCESS_SWITCH(Alice3HfTreeCreator3Prong, processLcWMl, "Process Lc with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3HfTreeCreator3Prong>(cfgc)};
}
