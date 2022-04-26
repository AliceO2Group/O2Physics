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
#include <cmath>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MathUtils/Utils.h"
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;

struct vertexingfwd {

  // Configurable<int> rangeBC{"rangeBC", 10, "Range for collision BCId and BCglobalIndex correspondance"};
  const int rangeBC = 10;

  HistogramRegistry registry{
    "registry",
    {
      {"EventsNtrkZvtx", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}}, //
      {"EventSelection", "; status; events", {HistType::kTH1F, {{4, 0.5, 4.5}}}},                                    //
      {"TrackDCAxy", "; DCA_{xy}; counts", {HistType::kTH1F, {{100, -10, 10}}}},                                     //
      {"TrackDCAx", "; DCA_{x}; counts", {HistType::kTH1F, {{100, -10, 10}}}},                                       //
      {"TrackDCAy", "; DCA_{y}; counts", {HistType::kTH1F, {{100, -10, 10}}}}                                        //
    }                                                                                                                //
  };

  //     std::vector<double> propagateToZlinear(double zEnd, aod::AmbiguousTrack const& track)
  //     {
  //       // Track parameters and their covariances linearly extrapolated to the plane at "zEnd".
  //
  //       //Calculate the jacobian related to the track parameters extrapolated to "zEnd"
  //       auto dZ = (zEnd - track.z());
  //       auto phi0 = track.phi();
  //       auto tanl0 = track.tgl();
  //       auto invtanl0 = 1.0 / tanl0;
  //       auto [sinphi0, cosphi0] = o2::math_utils::sincosd(phi0);
  //       auto n = dZ * invtanl0;
  //
  //       std::vector<double> T ={track.x(),track.y(),zEnd};
  //       // Extrapolate track parameters to "zEnd"
  //       T[0] += n * cosphi0;
  //       T[1] += n * sinphi0;
  //
  //       //Update to do : Add the jacobian and turn it to a method if possible
  //       return(T);
  ////     }

  void process(aod::AmbiguousTracks const& ambitracks, aod::BCs const&, aod::Collisions const& collisions, aod::Tracks const& tracks) // AmbiguousMFTTracks and fwd doesn't work yet
  {
    for (auto& ambitrack : ambitracks) {
      LOGF(info, "------------------------------------ We look at ambitrack %d which has %d possible BCs", ambitrack.globalIndex(), ambitrack.bc().size());
      auto extAmbiTrack = tracks.iteratorAt(ambitrack.globalIndex());
      printf("phi = %f\n", extAmbiTrack.phi());

      for (auto& bc : ambitrack.bc()) {
        // LOGF(info, "  BC %d with global BC %lld", bc.globalIndex(), bc.globalBC());

        for (auto& collision : collisions) {
          registry.fill(HIST("EventSelection"), 1.);

          // printf("collision BC ID = %lld, bc.globalIndex() %lld\n", collision.bcId(), bc.globalIndex());
          if ((collision.bcId() > (bc.globalIndex() - rangeBC)) && (collision.bcId() < (bc.globalIndex() + rangeBC))) {
            registry.fill(HIST("EventsNtrkZvtx"), ambitracks.size(), collision.posZ());
            registry.fill(HIST("EventSelection"), 2.);
            if (collision.bcId() == bc.globalIndex())
              registry.fill(HIST("EventSelection"), 3.);

            printf("collision BC ID = %d, bc.globalIndex() %lld\n", collision.bcId(), bc.globalIndex());
            printf("collision pos Z %f\n", collision.posZ());

            // T = propagateToZlinear(collision.posZ(), track);
            // printf("T[0] %f, T[1] %f\n", T[0], T[1]);
            /*
            SMatrix5 tpars(extAmbiTrack.x(), extAmbiTrack.y(), extAmbiTrack.phi(), extAmbiTrack.tgl(), extAmbiTrack.signed1Pt());

            std::vector<double> v1{extAmbiTrack.cXX(), extAmbiTrack.cXY(), extAmbiTrack.cYY(), extAmbiTrack.cPhiX(), extAmbiTrack.cPhiY(),
                extAmbiTrack.cPhiPhi(), extAmbiTrack.cTglX(), extAmbiTrack.cTglY(), extAmbiTrack.cTglPhi(), extAmbiTrack.cTglTgl(),
                extAmbiTrack.c1PtX(), extAmbiTrack.c1PtY(), extAmbiTrack.c1PtPhi(), extAmbiTrack.c1PtTgl(), extAmbiTrack.c1Pt21Pt2()};

            SMatrix55 tcovs(v1.begin(), v1.end());

            o2::track::TrackParCovFwd pars1{extAmbiTrack.z(), tpars, tcovs, extAmbiTrack.chi2()};

            pars1.propagateToZlinear(collision.posZ());

            //auto dca = std::sqrt((pars1.getX()-collision.posX())*(pars1.getX()-collision.posX())+(pars1.getY()-collision.posY())*(pars1.getY()-collision.posY()));

            const auto dcaX(pars1.getX()-collision.posX());
            const auto dcaY(pars1.getY()-collision.posY());
            auto dcaXY = std::sqrt(dcaX*dcaX + dcaY*dcaY);

            registry.fill(HIST("TrackDCAxy"), dcaXY);
            registry.fill(HIST("TrackDCAx"), dcaX);
            registry.fill(HIST("TrackDCAy"), dcaY);
*/
          }
        }
      }
    }
  }

  PROCESS_SWITCH(vertexingfwd, process, "Process ambiguous track DCA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexingfwd>(cfgc)};
}
