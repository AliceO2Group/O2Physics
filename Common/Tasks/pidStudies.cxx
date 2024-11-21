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

/// \file pidStudies.cxx
/// \brief task for studies of PID performance
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università and INFN Torino
/// \author Stefano Politanò <stefano.politano@cern.ch>, INFN Torino
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino
/// \author Marcello Di Costanzo <luca.aglietta@unito.it>, Università and INFN Torino

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
// #include "/home/mdicosta/alice/O2Physics/Common/TableProducer/Converters/mcCollisionConverter.cxx"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


namespace o2::aod
{
namespace pid_studies
{
DECLARE_SOA_COLUMN(MassK0, massK0, float);                  //! Candidate mass
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);          //! Candidate mass
DECLARE_SOA_COLUMN(MassOmega, massOmega, float);            //! Candidate mass
DECLARE_SOA_COLUMN(PtPos, ptPos, float);                    //! Transverse momentum of positive track (GeV/c)
DECLARE_SOA_COLUMN(PtNeg, ptNeg, float);                    //! Transverse momentum of negative track (GeV/c)
DECLARE_SOA_COLUMN(Radius, radius, float);                  //! Radius
DECLARE_SOA_COLUMN(Cpa, cpa, float);                        //! Cosine of pointing angle
DECLARE_SOA_COLUMN(NSigmaTpcPosPi, nSigmaTpcPosPi, float);  //! nSigmaTPC of positive track with pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcNegPi, nSigmaTpcNegPi, float);  //! nSigmaTPC of negative track with pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcPosKa, nSigmaTpcPosKa, float);  //! nSigmaTPC of positive track with kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcNegKa, nSigmaTpcNegKa, float);  //! nSigmaTPC of negative track with kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcPosPr, nSigmaTpcPosPr, float);  //! nSigmaTPC of positive track with proton hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcNegPr, nSigmaTpcNegPr, float);  //! nSigmaTPC of negative track with proton hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPosPi, nSigmaTofPosPi, float);  //! nSigmaTOF of positive track with pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTofNegPi, nSigmaTofNegPi, float);  //! nSigmaTOF of negative track with pion hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPosKa, nSigmaTofPosKa, float);  //! nSigmaTOF of positive track with kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTofNegKa, nSigmaTofNegKa, float);  //! nSigmaTOF of negative track with kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPosPr, nSigmaTofPosPr, float);  //! nSigmaTOF of positive track with proton hypothesis
DECLARE_SOA_COLUMN(NSigmaTofNegPr, nSigmaTofNegPr, float);  //! nSigmaTOF of negative track with proton hypothesis
DECLARE_SOA_COLUMN(AlphaArm, alphaArm, float);              //! Armenteros alpha
DECLARE_SOA_COLUMN(QtArm, qtArm, float);                    //! Armenteros Qt
DECLARE_SOA_COLUMN(OccupancyFt0c, occupancyFt0c, float);    //! Occupancy from FT0C
DECLARE_SOA_COLUMN(OccupancyIts, occupancyIts, float);      //! Occupancy from ITS
DECLARE_SOA_COLUMN(CentralityFT0C, centralityFT0C, float);  //! Centrality from FT0C
DECLARE_SOA_COLUMN(CentralityFT0M, centralityFT0M, float);  //! Centrality from FT0M
DECLARE_SOA_COLUMN(CandFlag, candFlag, int);                //! Flag for MC matching
} // namespace pid_studies

DECLARE_SOA_TABLE(pidV0s, "AOD", "PIDV0S", //! Table with PID information
                  pid_studies::MassK0,
                  pid_studies::MassLambda,
                  pid_studies::PtPos,
                  pid_studies::PtNeg,
                  pid_studies::Radius,
                  pid_studies::Cpa,
                  pid_studies::NSigmaTpcPosPi,
                  pid_studies::NSigmaTpcNegPi,
                  pid_studies::NSigmaTpcPosKa,
                  pid_studies::NSigmaTpcNegKa,
                  pid_studies::NSigmaTpcPosPr,
                  pid_studies::NSigmaTpcNegPr,
                  pid_studies::NSigmaTofPosPi,
                  pid_studies::NSigmaTofNegPi,
                  pid_studies::NSigmaTofPosKa,
                  pid_studies::NSigmaTofNegKa,
                  pid_studies::NSigmaTofPosPr,
                  pid_studies::NSigmaTofNegPr,
                  pid_studies::AlphaArm,
                  pid_studies::QtArm,
                  pid_studies::OccupancyFt0c,
                  pid_studies::OccupancyIts,
                  pid_studies::CentralityFT0C,                  
                  pid_studies::CentralityFT0M,                  
                  pid_studies::CandFlag                  
                  );

DECLARE_SOA_TABLE(pidCascades, "AOD", "PIDCASCADES", //! Table with PID information
                  // pid_studies::MassOmega,
                  // pid_studies::PtPos,
                  // pid_studies::PtNeg,
                  // pid_studies::Radius,
                  // pid_studies::Cpa,
                  // pid_studies::NSigmaTpcPosKa,
                  // pid_studies::NSigmaTpcNegKa,
                  pid_studies::OccupancyFt0c,
                  pid_studies::OccupancyIts,
                  pid_studies::CentralityFT0C,                  
                  pid_studies::CentralityFT0M,                  
                  pid_studies::CandFlag                  
                  );
} // namespace o2::aod


struct pidStudies {
  Produces <o2::aod::pidV0s> pidV0;
  Produces <o2::aod::pidCascades> pidCascade;
  HistogramRegistry registry{"registry", {}};

  using PIDTracks = soa::Join<aod::Tracks, aod::TracksExtra,
                            aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                            aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using CollSels = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms>;
  using V0sMCRec = soa::Join<aod::V0Datas, aod::V0CoreMCLabels>;
  using CascsMCRec = soa::Join<aod::CascDatas, aod::CascCoreMCLabels>;

  Configurable<float> massK0Min{"massK0Min", 0.4, "Minimum mass for K0"};
  Configurable<float> massK0Max{"massK0Max", 0.6, "Maximum mass for K0"};
  Configurable<float> massLambdaMin{"massLambdaMin", 1.0, "Minimum mass for lambda"};
  Configurable<float> massLambdaMax{"massLambdaMax", 1.3, "Maximum mass for lambda"};
  Configurable<float> massOmegaMin{"massOmegaMin", 1.5, "Minimum mass for omega"};
  Configurable<float> massOmegaMax{"massOmegaMax", 1.8, "Maximum mass for omega"};

  void init(InitContext&)
  {
  }

  template <typename Cand>
  void fillTree(Cand const& candidate, const int& flag)
  {
    if constexpr (std::is_same<Cand, V0sMCRec::iterator>::value) {
      const auto& posTrack = candidate.template posTrack_as<PIDTracks>();
      const auto& negTrack = candidate.template negTrack_as<PIDTracks>();
      pidV0(
        candidate.mK0Short(),
        candidate.mLambda(),
        posTrack.pt(),
        negTrack.pt(),
        candidate.v0radius(),
        candidate.v0cosPA(),
        posTrack.tofNSigmaPi(),
        negTrack.tofNSigmaPi(),
        posTrack.tofNSigmaKa(),
        negTrack.tofNSigmaKa(),
        posTrack.tofNSigmaPr(),
        negTrack.tofNSigmaPr(),
        posTrack.tpcNSigmaPi(),
        negTrack.tpcNSigmaPi(),
        posTrack.tpcNSigmaKa(),
        negTrack.tpcNSigmaKa(),
        posTrack.tpcNSigmaPr(),
        negTrack.tpcNSigmaPr(),
        candidate.alpha(),
        candidate.qtarm(),
        candidate.template collision_as<CollSels>().ft0cOccupancyInTimeRange(),
        candidate.template collision_as<CollSels>().trackOccupancyInTimeRange(),
        candidate.template collision_as<CollSels>().centFT0C(),
        candidate.template collision_as<CollSels>().centFT0M(),
        flag
      );
    } else {
      // const auto& posTrack = candidate.template posTrack_as<PIDTracks>();
      // const auto& negTrack = candidate.template negTrack_as<PIDTracks>();
      pidCascades(
        // candidate.mOmega(),
        // posTrack.pt(),
        // negTrack.pt(),
        // candidate.v0radius(),
        // candidate.v0cosPA(),
        // posTrack.tpcNSigmaKa(),
        // negTrack.tpcNSigmaKa(),
        candidate.template collision_as<CollSels>().ft0cOccupancyInTimeRange(),
        candidate.template collision_as<CollSels>().trackOccupancyInTimeRange(),
        candidate.template collision_as<CollSels>().centFT0C(),
        candidate.template collision_as<CollSels>().centFT0M(),
        flag
      );      
    }
  } 

  template <typename T1>
  int isMatched(const T1& cand) {
    int matched{0};
    if constexpr (std::is_same<T1, V0sMCRec::iterator>::value) {
      if (!cand.has_v0MCCore())
        return matched;
      auto v0MC = cand.template v0MCCore_as<aod::V0MCCores>();     
      // LOG(info) << "v0MC.pdgCode() " << v0MC.pdgCode();
      // LOG(info) << "v0MC.pdgCodeNegative() " << v0MC.pdgCodeNegative();
      // LOG(info) << "v0MC.pdgCodePositive() " << v0MC.pdgCodePositive();
      if (v0MC.pdgCode() == 3122 && v0MC.pdgCodeNegative() == -211 
          && v0MC.pdgCodePositive() == 2212) {
          matched = 1;
      }
      if (v0MC.pdgCode() == -3122 && v0MC.pdgCodeNegative() == -2212 
          && v0MC.pdgCodePositive() == 211) {
          matched = -1;
      }
      if (v0MC.pdgCode() == 310 && v0MC.pdgCodeNegative() == -211 
          && v0MC.pdgCodePositive() == 211) {
          matched = 2;
      } 
      if (v0MC.pdgCode() == -310 && v0MC.pdgCodeNegative() == -211 
          && v0MC.pdgCodePositive() == 211) {
          matched = -2;
      }
    }
    if constexpr (std::is_same<T1, CascsMCRec::iterator>::value) {
      if (!cand.has_cascMCCore())
        return matched;
      LOG(info) << "Inside checker";
      auto cascMC = cand.template cascMCCore_as<aod::CascMCCores>();     
      // LOG(info) << "cascMC.pdgCode() " << cascMC.pdgCode();
      // LOG(info) << "cascMC.pdgCodeBachelor() " << cascMC.pdgCodeBachelor();
      // LOG(info) << "cascMC.pdgCodeV0() " << cascMC.pdgCodeV0();
      // LOG(info) << "cascMC.pdgCodePositive() " << cascMC.pdgCodePositive();
      // LOG(info) << "cascMC.pdgCodeNegative() " << cascMC.pdgCodeNegative();
      if (cascMC.pdgCode() > 0) {
        if (cascMC.pdgCode() == 3334 && 
            cascMC.pdgCodeBachelor() == -321 && 
            cascMC.pdgCodeV0() == 3122 &&
            cascMC.pdgCodePositive() == 2212 &&
            cascMC.pdgCodeNegative() == -211) {
            matched = 3;
        }
      } else {
        if (cascMC.pdgCode() == -3334 && 
            cascMC.pdgCodeBachelor() == 321 && 
            cascMC.pdgCodeV0() == -3122 &&
            cascMC.pdgCodePositive() == 211 &&
            cascMC.pdgCodeNegative() == -2212) {
            matched = -3;
        }
      }
    }
    return matched;   
  }

  void processMC(V0sMCRec const& V0s, aod::V0MCCores const& V0sMC, CascsMCRec const& cascs, 
                 aod::CascMCCores const& cascsMC, CollSels const&, PIDTracks const&)  {
    for (const auto& v0 : V0s) {
      if (v0.mK0Short() > massK0Min && v0.mK0Short() < massK0Max ||
          v0.mLambda() > massLambdaMin && v0.mLambda() < massLambdaMax ||
          v0.mAntiLambda() > massLambdaMin && v0.mAntiLambda() < massLambdaMax) {
        int matched = isMatched(v0); 
        if(matched != 0) {
          fillTree(v0, matched);
        }
      }
    }
    for (const auto& casc : cascs) {
      LOG(info) << "checking omega";
      if (casc.mOmega() > massOmegaMin && casc.mOmega() < massOmegaMax) {
        LOG(info) << "mass passed";
        int matched = isMatched(casc); 
        if(matched != 0) {
          LOG(info) << "Matched omega";
          // fillTree(casc, matched);
        }
      }
    }
    
  }
  PROCESS_SWITCH(pidStudies, processMC, "process MC", true);

  void processData(aod::V0Datas const& V0s, aod::Cascades const& cascades, CollSels const&, PIDTracks const&)
  {
    for (const auto& v0 : V0s) {
      if (v0.mK0Short() > massK0Min && v0.mK0Short() < massK0Max ||
          v0.mLambda() > massLambdaMin && v0.mLambda() < massLambdaMax ||
          v0.mAntiLambda() > massLambdaMin && v0.mAntiLambda() < massLambdaMax) {
        fillTree(v0, 0);
      }
    }
  }
  PROCESS_SWITCH(pidStudies, processData, "process data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pidStudies>(cfgc)};
}
