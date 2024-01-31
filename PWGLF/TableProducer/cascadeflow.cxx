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
///
/// \brief Task to create derived data
///

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/cascqaanalysis.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "TRandom3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TrkPidInfo = soa::Join<aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa>;
using DauTracks = soa::Join<aod::TracksIU, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, TrkPidInfo>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

namespace cascadev2 {
enum species { Xi = 0,
               Omega = 1 };
constexpr double massSigmaParameters[4][2]{
  {4.9736e-3, 0.006815},
  {-2.39594, -2.257},
  {1.8064e-3, 0.00138},
  {1.03468e-1, 0.1898}};
static const std::vector<std::string> massSigmaParameterNames{"p0", "p1", "p2", "p3"};
static const std::vector<std::string> speciesNames{"Xi", "Omega"};
}

struct cascadeFlow {

  // Tables to produce
  Produces<aod::CascTraining> trainingSample;

  Configurable<double> sideBandStart{"sideBandStart", 5, "Start of the sideband region in number of sigmas"};
  Configurable<double> sideBandEnd{"sideBandEnd", 7, "End of the sideband region in number of sigmas"};
  Configurable<double> downsample{"downsample", 1., "Downsample training output tree"};
  Configurable<LabeledArray<double>> parSigmaMass{
    "parSigmaMass",
    {cascadev2::massSigmaParameters[0], 4, 2,
     cascadev2::massSigmaParameterNames, cascadev2::speciesNames},
    "Mass resolution parameters: [0]*exp([1]*x)+[2]*exp([3]*x)"};
  
  float getNsigmaMass(const cascadev2::species s, const float pt, const float nsigma = 6)
  {
    const auto sigma = parSigmaMass->get(0u, s) * exp(parSigmaMass->get(1, s) * pt) + parSigmaMass->get(2, s) * exp(parSigmaMass->get(3, s) * pt);
    return nsigma * sigma;
  }

  template<class collision_t, class cascade_t>
  void fillTrainingTable(collision_t coll, cascade_t casc, int pdgCode) {
    trainingSample(coll.centFT0M(), 
        casc.sign(), 
        casc.pt(), 
        casc.eta(), 
        casc.mXi(), 
        casc.mOmega(), 
        casc.mLambda(), 
        casc.cascradius(), 
        casc.v0radius(),
        casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()),
        casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()),
        casc.dcapostopv(),
        casc.dcanegtopv(),
        casc.dcabachtopv(),
        casc.dcacascdaughters(),
        casc.dcaV0daughters(),
        casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()),
        casc.bachBaryonCosPA(),
        casc.bachBaryonDCAxyToPV(),
        pdgCode); 
  }

  //  void processTrainingBackground(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs> const& Cascades, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&) {
  void processTrainingBackground(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs> const& Cascades, DauTracks const&) {
    
    if (!coll.sel8() || std::abs(coll.posZ()) > 10.) {
      return;
    }

    for (auto& casc : Cascades) {
      if (gRandom->Uniform() > downsample) {
        continue;
      }
      
      float sigmaRangeXi[2]{getNsigmaMass(cascadev2::Xi, casc.pt(), sideBandStart), getNsigmaMass(cascadev2::Xi, casc.pt(), sideBandEnd)};
      float sigmaRangeOmega[2]{getNsigmaMass(cascadev2::Omega, casc.pt(), sideBandStart), getNsigmaMass(cascadev2::Omega, casc.pt(), sideBandEnd)};

      if ((std::abs(casc.mXi() - constants::physics::MassXiMinus) < sigmaRangeXi[0] || 
	   std::abs(casc.mXi() - constants::physics::MassXiMinus) > sigmaRangeXi[1]) &&
	  (std::abs(casc.mOmega() - constants::physics::MassOmegaMinus) < sigmaRangeOmega[0] || 
	   std::abs(casc.mOmega() - constants::physics::MassOmegaMinus) > sigmaRangeOmega[1])) {
	continue;
      }

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      fillTrainingTable(coll, casc, 0); 
    }
  }
  
  PROCESS_SWITCH(cascadeFlow, processTrainingBackground, "Process to create the training dataset for the background", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeFlow>(cfgc, TaskName{"lf-cascade-flow"})};
}
