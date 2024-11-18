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

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


namespace o2::aod
{
namespace pid_studies
{
DECLARE_SOA_COLUMN(MassK0, massK0, float);                  //! Candidate mass
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);          //! Candidate mass
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
DECLARE_SOA_COLUMN(OccupancyFt0c, occupancyFt0c, float);    //! Occupancy of FT0C
DECLARE_SOA_COLUMN(OccupancyIts, occupancyIts, float);      //! Occupancy of ITS
DECLARE_SOA_COLUMN(CentralityFT0C, centralityFT0C, float);  //! Centrality from FT0C
DECLARE_SOA_COLUMN(CentralityFT0M, centralityFT0M, float);  //! Centrality from FT0M
} // namespace pid_studies

DECLARE_SOA_TABLE(pidInformation, "AOD", "PIDSTUDIES", //! Table with PID information
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
                  pid_studies::CentralityFT0M                  
                  );
} // namespace o2::aod


struct pidStudies {
  Produces <o2::aod::pidInformation> pidInformation;
  HistogramRegistry registry{"registry", {}};

  using PIDTracks = soa::Join<aod::Tracks, aod::TracksExtra,
                            aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                            aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using CollSels = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms>;

  void init(InitContext&)
  {
  }

  template <bool isMc, typename Cand>
  void fillTree(Cand const& candidate)
  {
    const auto& posTrack = candidate.template posTrack_as<PIDTracks>();
    const auto& negTrack = candidate.template negTrack_as<PIDTracks>();
    pidInformation(
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
      candidate.template collision_as<CollSels>().centFT0M()
    );
  }

  void processData(aod::V0Datas const& V0s, aod::Cascades const& cascades, CollSels const&, PIDTracks const&)
  {
    for (const auto& v0 : V0s) {
      fillTree<false>(v0);
    }
  }
  PROCESS_SWITCH(pidStudies, processData, "process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pidStudies>(cfgc)};
}
