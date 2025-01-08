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

/// \file taskPidStudies.cxx
/// \brief task for studies of PID performance
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università and INFN Torino
/// \author Stefano Politanò <stefano.politano@cern.ch>, INFN Torino
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino
/// \author Luca Aglietta <luca.aglietta@unito.it>, Università and INFN Torino

#include "TPDGCode.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace pid_studies
{
enum Particle { NotMatched = 0,
                K0s,
                Lambda,
                Omega };
// V0s
DECLARE_SOA_COLUMN(MassK0, massK0, float);                 //! Candidate mass
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);         //! Candidate mass
DECLARE_SOA_COLUMN(MassAntiLambda, massAntiLambda, float); //! Candidate mass
DECLARE_SOA_COLUMN(Pt, pt, float);                         //! Transverse momentum of the candidate (GeV/c)
DECLARE_SOA_COLUMN(PtPos, ptPos, float);                   //! Transverse momentum of positive track (GeV/c)
DECLARE_SOA_COLUMN(PtNeg, ptNeg, float);                   //! Transverse momentum of negative track (GeV/c)
DECLARE_SOA_COLUMN(PtPosTpc, ptPosTpc, float);             //! Transverse Momentum of positive track at inner wall of TPC (GeV/c)
DECLARE_SOA_COLUMN(PtNegTpc, ptNegTpc, float);             //! Transverse Momentum of negative track at inner wall of TPC (GeV/c)
DECLARE_SOA_COLUMN(Radius, radius, float);                 //! Radius
DECLARE_SOA_COLUMN(Cpa, cpa, float);                       //! Cosine of pointing angle
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float); //! DCA between V0 daughters
DECLARE_SOA_COLUMN(DcaV0ToPv, dcaV0ToPv, float);           //! DCA V0 to PV
DECLARE_SOA_COLUMN(NSigmaTpcPosPi, nSigmaTpcPosPi, float); //! nSigmaTPC of positive track with pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcNegPi, nSigmaTpcNegPi, float); //! nSigmaTPC of negative track with pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcPosPr, nSigmaTpcPosPr, float); //! nSigmaTPC of positive track with proton hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcNegPr, nSigmaTpcNegPr, float); //! nSigmaTPC of negative track with proton hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPosPi, nSigmaTofPosPi, float); //! nSigmaTOF of positive track with pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTofNegPi, nSigmaTofNegPi, float); //! nSigmaTOF of negative track with pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPosPr, nSigmaTofPosPr, float); //! nSigmaTOF of positive track with proton hypothesis
DECLARE_SOA_COLUMN(NSigmaTofNegPr, nSigmaTofNegPr, float); //! nSigmaTOF of negative track with proton hypothesis
DECLARE_SOA_COLUMN(AlphaArm, alphaArm, float);             //! Armenteros alpha
DECLARE_SOA_COLUMN(QtArm, qtArm, float);                   //! Armenteros Qt

// Cascades
DECLARE_SOA_COLUMN(MassOmega, massOmega, float);             //! Candidate mass
DECLARE_SOA_COLUMN(MassXi, massXi, float);                   //! Candidate mass
DECLARE_SOA_COLUMN(PtBach, ptBach, float);                   //! Transverse momentum of the bachelor (GeV/c)
DECLARE_SOA_COLUMN(PtBachTpc, ptBachTpc, float);             //! Transverse momentum of the bachelor at inner wall of TPC (GeV/c)
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                 //! Daughter lambda mass (GeV/c^2)
DECLARE_SOA_COLUMN(V0cosPA, v0cosPA, float);                 //! V0 CPA
DECLARE_SOA_COLUMN(CascCosPa, cascCosPa, float);             //! Cascade CPA
DECLARE_SOA_COLUMN(NSigmaTpcBachKa, nSigmaTpcBachKa, float); //! nSigmaTPC of bachelor with kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTofBachKa, nSigmaTofBachKa, float); //! nSigmaTOF of bachelor with kaon hypothesis

// Common columns
DECLARE_SOA_COLUMN(OccupancyFt0c, occupancyFt0c, float);   //! Occupancy from FT0C
DECLARE_SOA_COLUMN(OccupancyIts, occupancyIts, float);     //! Occupancy from ITS
DECLARE_SOA_COLUMN(CentralityFT0C, centralityFT0C, float); //! Centrality from FT0C
DECLARE_SOA_COLUMN(CentralityFT0M, centralityFT0M, float); //! Centrality from FT0M
DECLARE_SOA_COLUMN(CandFlag, candFlag, int);               //! Flag for MC matching
} // namespace pid_studies

DECLARE_SOA_TABLE(PidV0s, "AOD", "PIDV0S", //! Table with PID information
                  pid_studies::MassK0,
                  pid_studies::MassLambda,
                  pid_studies::MassAntiLambda,
                  pid_studies::Pt,
                  pid_studies::PtPos,
                  pid_studies::PtNeg,
                  pid_studies::PtPosTpc,
                  pid_studies::PtNegTpc,
                  pid_studies::Radius,
                  pid_studies::Cpa,
                  pid_studies::DcaV0Daughters,
                  pid_studies::DcaV0ToPv,
                  pid_studies::NSigmaTpcPosPi,
                  pid_studies::NSigmaTpcNegPi,
                  pid_studies::NSigmaTpcPosPr,
                  pid_studies::NSigmaTpcNegPr,
                  pid_studies::NSigmaTofPosPi,
                  pid_studies::NSigmaTofNegPi,
                  pid_studies::NSigmaTofPosPr,
                  pid_studies::NSigmaTofNegPr,
                  pid_studies::AlphaArm,
                  pid_studies::QtArm,
                  pid_studies::OccupancyFt0c,
                  pid_studies::OccupancyIts,
                  pid_studies::CentralityFT0C,
                  pid_studies::CentralityFT0M,
                  pid_studies::CandFlag);

DECLARE_SOA_TABLE(PidCascades, "AOD", "PIDCASCADES", //! Table with PID information
                  pid_studies::MassOmega,
                  pid_studies::Pt,
                  pid_studies::PtBach,
                  pid_studies::PtBachTpc,
                  pid_studies::Radius,
                  pid_studies::MLambda,
                  pid_studies::V0cosPA,
                  pid_studies::MassXi,
                  pid_studies::CascCosPa,
                  pid_studies::DcaV0Daughters,
                  pid_studies::DcaV0ToPv,
                  pid_studies::NSigmaTpcBachKa,
                  pid_studies::NSigmaTofBachKa,
                  pid_studies::OccupancyFt0c,
                  pid_studies::OccupancyIts,
                  pid_studies::CentralityFT0C,
                  pid_studies::CentralityFT0M,
                  pid_studies::CandFlag);
} // namespace o2::aod

struct HfTaskPidStudies {
  Produces<o2::aod::PidV0s> pidV0;
  Produces<o2::aod::PidCascades> pidCascade;

  Configurable<float> massK0Min{"massK0Min", 0.4, "Minimum mass for K0"};
  Configurable<float> massK0Max{"massK0Max", 0.6, "Maximum mass for K0"};
  Configurable<float> massLambdaMin{"massLambdaMin", 1.0, "Minimum mass for lambda"};
  Configurable<float> massLambdaMax{"massLambdaMax", 1.3, "Maximum mass for lambda"};
  Configurable<float> massOmegaMin{"massOmegaMin", 1.5, "Minimum mass for omega"};
  Configurable<float> massOmegaMax{"massOmegaMax", 1.8, "Maximum mass for omega"};
  Configurable<float> radiusMax{"radiusMax", 2.3, "Maximum decay radius (cm)"};
  Configurable<float> cosPaMin{"cosPaMin", 0.98, "Minimum cosine of pointing angle"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.2, "Maximum DCA among the V0 daughters (cm)"};
  Configurable<float> dcaV0ToPvMax{"dcaV0ToPvMax", 0.2, "Maximum DCA of the V0 from the primary vertex (cm)"};
  Configurable<float> cosPaV0Min{"cosPaV0Min", 0.95, "Minimum cosine of pointing angle for V0 stemming from cascade decays"};
  Configurable<float> qtArmenterosMinForK0{"qtArmenterosMinForK0", 0.12, "Minimum Armenteros' qt for K0"};
  Configurable<float> qtArmenterosMaxForLambda{"qtArmenterosMaxForLambda", 0.12, "Minimum Armenteros' qt for (anti)Lambda"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of candidates to keep"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  using PidTracks = soa::Join<aod::Tracks, aod::TracksExtra,
                              aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                              aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using CollSels = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms>;
  using V0sMcRec = soa::Join<aod::V0Datas, aod::V0CoreMCLabels>;
  using CascsMcRec = soa::Join<aod::CascDatas, aod::CascCoreMCLabels>;

  void init(InitContext&)
  {
    if ((doprocessV0Mc && doprocessV0Data) || (doprocessCascMc && doprocessCascData)) {
      LOGP(fatal, "Both data and MC process functions were enabled! Please check your configuration!");
    }
  }

  template <bool isV0, typename Cand>
  void fillTree(Cand const& candidate, const int flag)
  {
    float pseudoRndm = candidate.pt() * 1000. - static_cast<int64_t>(candidate.pt() * 1000);
    if (candidate.pt() < ptMaxForDownSample && pseudoRndm > downSampleBkgFactor) {
      return;
    }

    const auto& coll = candidate.template collision_as<CollSels>();
    if constexpr (isV0) {
      const auto& posTrack = candidate.template posTrack_as<PidTracks>();
      const auto& negTrack = candidate.template negTrack_as<PidTracks>();
      pidV0(
        candidate.mK0Short(),
        candidate.mLambda(),
        candidate.mAntiLambda(),
        candidate.pt(),
        posTrack.pt(),
        negTrack.pt(),
        posTrack.tpcInnerParam() / std::cosh(candidate.positiveeta()),
        negTrack.tpcInnerParam() / std::cosh(candidate.negativeeta()),
        candidate.v0radius(),
        candidate.v0cosPA(),
        candidate.dcaV0daughters(),
        candidate.dcav0topv(),
        posTrack.tpcNSigmaPi(),
        negTrack.tpcNSigmaPi(),
        posTrack.tpcNSigmaPr(),
        negTrack.tpcNSigmaPr(),
        posTrack.tofNSigmaPi(),
        negTrack.tofNSigmaPi(),
        posTrack.tofNSigmaPr(),
        negTrack.tofNSigmaPr(),
        candidate.alpha(),
        candidate.qtarm(),
        coll.ft0cOccupancyInTimeRange(),
        coll.trackOccupancyInTimeRange(),
        coll.centFT0C(),
        coll.centFT0M(),
        flag);
    } else {
      const auto& bachTrack = candidate.template bachelor_as<PidTracks>();
      pidCascade(
        candidate.mOmega(),
        candidate.pt(),
        candidate.bachelorpt(),
        bachTrack.tpcInnerParam() / std::cosh(candidate.bacheloreta()),
        candidate.cascradius(),
        candidate.mLambda(),
        candidate.v0cosPA(coll.posX(), coll.posY(), coll.posZ()),
        candidate.mXi(),
        candidate.casccosPA(coll.posX(), coll.posY(), coll.posZ()),
        candidate.dcaV0daughters(),
        candidate.dcav0topv(coll.posX(), coll.posY(), coll.posZ()),
        bachTrack.tpcNSigmaKa(),
        bachTrack.tofNSigmaKa(),
        coll.ft0cOccupancyInTimeRange(),
        coll.trackOccupancyInTimeRange(),
        coll.centFT0C(),
        coll.centFT0M(),
        flag);
    }
  }

  template <typename T1>
  int isMatched(const T1& cand)
  {
    if constexpr (std::is_same<T1, V0sMcRec::iterator>::value) {
      if (!cand.has_v0MCCore()) {
        return aod::pid_studies::Particle::NotMatched;
      }
      auto v0MC = cand.template v0MCCore_as<aod::V0MCCores>();
      if (v0MC.pdgCode() == kK0Short && v0MC.pdgCodeNegative() == -kPiPlus && v0MC.pdgCodePositive() == kPiPlus) {
        return aod::pid_studies::Particle::K0s;
      }
      if (v0MC.pdgCode() == kLambda0 && v0MC.pdgCodeNegative() == -kPiPlus && v0MC.pdgCodePositive() == kProton) {
        return aod::pid_studies::Particle::Lambda;
      }
      if (v0MC.pdgCode() == -kLambda0 && v0MC.pdgCodeNegative() == -kProton && v0MC.pdgCodePositive() == kPiPlus) {
        return -aod::pid_studies::Particle::Lambda;
      }
    }
    if constexpr (std::is_same<T1, CascsMcRec::iterator>::value) {
      if (!cand.has_cascMCCore()) {
        return aod::pid_studies::Particle::NotMatched;
      }
      auto cascMC = cand.template cascMCCore_as<aod::CascMCCores>();
      if (cascMC.pdgCode() == kOmegaMinus &&
          cascMC.pdgCodeBachelor() == -kKPlus &&
          cascMC.pdgCodeV0() == kLambda0 &&
          cascMC.pdgCodePositive() == kProton &&
          cascMC.pdgCodeNegative() == -kPiPlus) {
        return aod::pid_studies::Particle::Omega;
      }
      if (cascMC.pdgCode() == -kOmegaMinus &&
          cascMC.pdgCodeBachelor() == kKPlus &&
          cascMC.pdgCodeV0() == -kLambda0 &&
          cascMC.pdgCodePositive() == kPiPlus &&
          cascMC.pdgCodeNegative() == -kProton) {
        return -aod::pid_studies::Particle::Omega;
      }
    }
    return aod::pid_studies::Particle::NotMatched;
  }

  template <typename V0Cand>
  bool isSelectedV0AsK0s(const V0Cand& v0)
  {
    if (v0.mK0Short() < massK0Min || v0.mK0Short() > massK0Max) {
      return false;
    }
    if (v0.qtarm() < qtArmenterosMinForK0) {
      return false;
    }
    if (v0.v0radius() > radiusMax) {
      return false;
    }
    if (v0.v0cosPA() < cosPaMin) {
      return false;
    }
    if (v0.dcaV0daughters() > dcaV0DaughtersMax) {
      return false;
    }
    if (v0.dcav0topv() > dcaV0ToPvMax) {
      return false;
    }
    return true;
  }

  template <typename V0Cand>
  bool isSelectedV0AsLambda(const V0Cand& v0)
  {
    if ((v0.mLambda() < massLambdaMin || v0.mLambda() > massLambdaMax) &&
        (v0.mAntiLambda() < massLambdaMin || v0.mAntiLambda() > massLambdaMax)) {
      return false;
    }
    if (v0.qtarm() > qtArmenterosMaxForLambda) {
      return false;
    }
    if (v0.v0radius() > radiusMax) {
      return false;
    }
    if (v0.v0cosPA() < cosPaMin) {
      return false;
    }
    if (v0.dcaV0daughters() > dcaV0DaughtersMax) {
      return false;
    }
    if (v0.dcav0topv() > dcaV0ToPvMax) {
      return false;
    }
    return true;
  }

  template <typename CascCand>
  bool isSelectedCascAsOmega(const CascCand& casc)
  {
    if (casc.mOmega() < massOmegaMin || casc.mOmega() > massOmegaMax) {
      return false;
    }
    if (casc.mLambda() < massLambdaMin || casc.mLambda() > massLambdaMax) {
      return false;
    }
    if (casc.cascradius() > radiusMax) {
      return false;
    }
    const auto& coll = casc.template collision_as<CollSels>();
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < cosPaMin) {
      return false;
    }
    if (casc.dcaV0daughters() > dcaV0DaughtersMax) {
      return false;
    }
    if (casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()) > dcaV0ToPvMax) {
      return false;
    }
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < cosPaV0Min) {
      return false;
    }
    return true;
  }

  void processV0Mc(V0sMcRec const& V0s,
                   aod::V0MCCores const&,
                   CollSels const&,
                   PidTracks const&)
  {
    for (const auto& v0 : V0s) {
      if (isSelectedV0AsK0s(v0) || isSelectedV0AsLambda(v0)) {
        int matched = isMatched(v0);
        if (matched != aod::pid_studies::Particle::NotMatched) {
          fillTree<true>(v0, matched);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskPidStudies, processV0Mc, "Process MC", true);

  void processV0Data(aod::V0Datas const& V0s,
                     CollSels const&,
                     PidTracks const&)
  {
    for (const auto& v0 : V0s) {
      if (isSelectedV0AsK0s(v0) || isSelectedV0AsLambda(v0)) {
        fillTree<true>(v0, aod::pid_studies::Particle::NotMatched);
      }
    }
  }
  PROCESS_SWITCH(HfTaskPidStudies, processV0Data, "Process data", false);

  void processCascMc(CascsMcRec const& cascades,
                     aod::CascMCCores const&,
                     CollSels const&,
                     PidTracks const&)
  {
    for (const auto& casc : cascades) {
      if (isSelectedCascAsOmega(casc)) {
        int matched = isMatched(casc);
        if (matched != aod::pid_studies::Particle::NotMatched) {
          fillTree<false>(casc, matched);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskPidStudies, processCascMc, "Process MC", true);

  void processCascData(aod::CascDatas const& cascades,
                       CollSels const&,
                       PidTracks const&)
  {
    for (const auto& casc : cascades) {
      if (isSelectedCascAsOmega(casc)) {
        fillTree<false>(casc, aod::pid_studies::Particle::NotMatched);
      }
    }
  }
  PROCESS_SWITCH(HfTaskPidStudies, processCascData, "Process data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskPidStudies>(cfgc)};
}
