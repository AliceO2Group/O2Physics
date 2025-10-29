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

/// \file   match-mft-ft0.cxx
/// \author Sarah Herrmann <sarah.herrmann@cern.ch>
///
/// \brief This code loops over every MFT tracks (except orphan tracks) and propagates
///        them to the FT0-C, matching the signals in some BC to reduce track ambiguity
///        It produces a table containing for each MFT track a list of BCs with an FT0C match
///        called aod::BCofMFT
/// \date 03/09/24
/// \note https://github.com/AliceO2Group/AliceO2/blob/dev/DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackFwd.h

#include "Common/DataModel/MatchMFTFT0.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/TrackFwd.h> //for propagate

#include <Math/MatrixFunctions.h>
#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <TGeoGlobalMagField.h>

#include <RtypesCore.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using namespace o2;
using namespace o2::framework;

// Creating a table BC to FT0 and filling it
struct bctoft0c {
  Produces<aod::MatchedToFT0> mf;
  struct {
    std::vector<int> ft0ids;
  } filler;
  void process(aod::BCs::iterator const& bc, soa::SmallGroups<aod::FT0s> const& ft0s)
  {
    filler.ft0ids.clear();
    for (auto const& ft0 : ft0s) {
      filler.ft0ids.emplace_back(ft0.globalIndex());
    }
    mf(bc.globalIndex(), filler.ft0ids);
  }
};

using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedToFT0>;

template <typename T>
T getCompatibleBCs(aod::AmbiguousMFTTrack const& atrack, aod::Collision const& collOrig, T const& bcs, int deltaBC)
{
  // this method is unused for now, MFTtracks with no collisions (orphan tracks) are not considered for the matching with FT0C
  auto compBCs = atrack.bc_as<T>(); // BC + info on FT0
  auto bcIter = compBCs.begin();    // first element of compBC
  uint64_t firstBC = bcIter.globalBC();

  bcIter.moveToEnd();                  // does it move to the end or the next one after the end ?
  --bcIter;                            // to avoid a seg fault
  uint64_t lastBC = bcIter.globalBC(); // gives the last COMPATIBLE BC in compBCs

  auto bcIt = collOrig.bc_as<T>();

  int64_t minBCId = bcIt.globalIndex();
  auto minGlobalBC = bcIt.globalBC();

  if (bcIt.globalBC() < firstBC + deltaBC) {
    while (bcIt != bcs.end() && bcIt.globalBC() < firstBC + deltaBC) {
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();

      ++bcIt;
    }
    if (bcIt == bcs.end()) {
      --bcIt;
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();
    }
  } else {
    // here bcIt.globalBC() >= firstBC + deltaBC

    while (bcIt != bcs.begin() && bcIt.globalBC() > firstBC + deltaBC) {
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();

      --bcIt;
    }
  }

  int64_t maxBCId = bcIt.globalIndex();

  while (bcIt != bcs.end() && bcIt.globalBC() < lastBC + deltaBC) {
    maxBCId = bcIt.globalIndex();

    ++bcIt;
  }

  if (bcIt != bcs.end() && maxBCId >= minBCId) {
    T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};
    bcs.copyIndexBindings(slice);
    return slice;
  } else {
    T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId)}, (uint64_t)minBCId};
    bcs.copyIndexBindings(slice);
    return slice;
  }
}

template <typename T>
T getCompatibleBCs(aod::MFTTracks::iterator const& track, aod::Collision const& collOrig, T const& bcs, int deltaBC)
{

  // define firstBC and lastBC (globalBC of beginning and end of the range, when no shift is applied)

  auto bcIt = collOrig.bc_as<T>();
  // auto timstp = bcIt.timestamp();

  int64_t firstBC = bcIt.globalBC() + (track.trackTime() - track.trackTimeRes()) / o2::constants::lhc::LHCBunchSpacingNS;
  int64_t lastBC = firstBC + 2 * track.trackTimeRes() / o2::constants::lhc::LHCBunchSpacingNS + 1; // to have a delta = 198 BC

  // printf(">>>>>>>>>>>>>>>>>>>>>>>>>>> last-first %lld\n", lastBC-firstBC);

  // int collTimeResInBC = collOrig.collisionTimeRes()/o2::constants::lhc::LHCBunchSpacingNS;

  // int64_t collFirstBC = bcIt.globalBC() + (collOrig.collisionTime() - collOrig.collisionTimeRes())/o2::constants::lhc::LHCBunchSpacingNS;
  // int64_t collLastBC = collFirstBC + 2*collOrig.collisionTimeRes()/o2::constants::lhc::LHCBunchSpacingNS +1;

  int64_t minBCId = bcIt.globalIndex();

  if ((int64_t)bcIt.globalBC() < firstBC + deltaBC) {
    while (bcIt != bcs.end() && (int64_t)bcIt.globalBC() < firstBC + deltaBC) {
      minBCId = bcIt.globalIndex();

      ++bcIt;
    }
    if (bcIt == bcs.end()) {
      --bcIt;
      // allows to avoid bcIt==bcs.end() in the following
    }
    // minGlobalBC needs to be >= to firstBC+deltaBC
    minBCId = bcIt.globalIndex();

  } else {
    // here bcIt.globalBC() >= firstBC + deltaBC

    while (bcIt != bcs.begin() && (int64_t)bcIt.globalBC() >= (int64_t)firstBC + deltaBC) {
      minBCId = bcIt.globalIndex();
      --bcIt;
    }
    if (bcIt == bcs.begin() && (int64_t)bcIt.globalBC() >= (int64_t)firstBC + deltaBC) {
      minBCId = bcIt.globalIndex();
    }
    ++bcIt; // retrieve the pointer which gave minBCId and minGlobalBC
    if (bcIt == bcs.end()) {
      --bcIt; // go back if we got to the end of the list
    }
  }

  int64_t maxBCId = bcIt.globalIndex();

  if ((int64_t)bcIt.globalBC() > (int64_t)lastBC + deltaBC) {
    // the previous minimum is actually bigger than the right boundary

    if (bcIt != bcs.begin()) {
      --bcIt;                                                    // let's check the previous element in the BC list
      if ((int64_t)bcIt.globalBC() < (int64_t)firstBC + deltaBC) // if this previous element is smaller than the left boundary
      {
        // means that the slice of compatible BCs is empty

        T slice{{bcs.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
        // bcs.copyIndexBindings(slice); REMOVED IT BECAUSE I DON'T KNOW WHAT IT DOES HERE
        return slice; // returns an empty slice
      }
    }
  }

  if ((int64_t)bcIt.globalBC() < (int64_t)firstBC + deltaBC) {
    // the previous minimum is actually smaller than the right boundary
    ++bcIt;

    if (bcIt != bcs.end() && ((int64_t)bcIt.globalBC() > (int64_t)lastBC + deltaBC)) {
      // check the following element

      T slice{{bcs.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
      // bcs.copyIndexBindings(slice); REMOVED IT BECAUSE I DON'T KNOW WHAT IT DOES HERE
      return slice; // returns an empty slice
    }
  }

  while (bcIt != bcs.end() && (int64_t)bcIt.globalBC() <= (int64_t)lastBC + deltaBC) {
    maxBCId = bcIt.globalIndex();

    ++bcIt;
  }

  if (maxBCId < minBCId) {
    if (bcIt == bcs.end()) {
      printf("at the end of the bcs iterator %d\n", 1);
    }
    T slice{{bcs.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
    // bcs.copyIndexBindings(slice); REMOVED IT BECAUSE I DON'T KNOW WHAT IT DOES HERE
    return slice; // returns an empty slice
  }

  T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};
  bcs.copyIndexBindings(slice);
  return slice;
}

struct matchmftft0 {
  Produces<aod::BCofMFT> BcMft;
  struct {
    std::vector<int> BCids;
  } filler;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber = -1;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  int count = 0;
  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<bool> strictBCSel{"strictBCSel", false, "force the BC of the match to have FT0A&C signals"};
  Configurable<int> shiftBC{"shiftBC", 0, "shift in BC wrt normal"}; // should be kept at zero except if the time-alignment MFT-FT0C must be redone

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  std::vector<std::vector<float>> channelCoord = {{103.2, 17.8, -813.1}, {76.9, 17.8, -815.9}, {103.1, 44.2, -812.1}, {76.8, 44.2, -814.9}, {103.2, 78.7, -810}, {76.8, 79, -812.9}, {103.2, 105, -807.1}, {76.8, 105.3, -810}, {43.2, 78.8, -815}, {43.2, 105.1, -812.1}, {16.8, 78.9, -815.9}, {16.8, 105.2, -813}, {-16.8, 105.2, -813}, {-16.8, 78.9, -815.9}, {-43.2, 105.1, -812.1}, {-43.2, 78.8, -815}, {-76.8, 105.3, -810}, {-76.8, 79, -812.9}, {-103.2, 105, -807.1}, {-103.2, 78.7, -810}, {-76.8, 44.2, -814.9}, {-103.1, 44.2, -812.1}, {-76.9, 17.8, -815.9}, {-103.2, 17.8, -813.1}, {-103.2, -17.8, -813.1}, {-76.9, -17.8, -815.9}, {-103.1, -44.2, -812.1}, {-76.8, -44.2, -814.9}, {-103.2, -78.7, -810}, {-76.8, -79, -812.9}, {-103.2, -105, -807.1}, {-76.8, -105.3, -810}, {-43.2, -78.8, -815}, {-43.2, -105.1, -812.1}, {-16.8, -78.9, -815.9}, {-16.8, -105.2, -813}, {16.8, -105.2, -813}, {16.8, -78.9, -815.9}, {43.2, -105.1, -812.1}, {43.2, -78.8, -815}, {76.8, -105.3, -810}, {76.8, -79, -812.9}, {103.2, -105, -807.1}, {103.2, -78.7, -810}, {76.8, -44.2, -814.9}, {103.1, -44.2, -812.1}, {76.9, -17.8, -815.9}, {103.2, -17.8, -813.1}, {163, 18.7, -804.1}, {137, 18.9, -808.9}, {163, 45.2, -803.1}, {137, 45.3, -807.9}, {163, 78.6, -800.1}, {137, 79.1, -804.9}, {163, 104.9, -797.2}, {137, 105.4, -801.9}, {103.4, 138, -802}, {102.9, 164, -797.2}, {77.1, 138, -804.9}, {76.6, 164, -800}, {43.3, 139, -807}, {43.2, 165, -802.1}, {16.9, 139, -807.9}, {16.7, 165, -803}, {-16.7, 165, -803}, {-16.9, 139, -807.9}, {-43.2, 165, -802.1}, {-43.3, 139, -807}, {-76.6, 164, -800}, {-77.1, 138, -804.9}, {-102.9, 164, -797.2}, {-103.4, 138, -802}, {-137, 105.4, -801.9}, {-163, 104.9, -797.2}, {-137, 79.1, -804.9}, {-163, 78.6, -800.1}, {-137, 45.3, -807.9}, {-163, 45.2, -803.1}, {-137, 18.9, -808.9}, {-163, 18.7, -804.1}, {-163, -18.7, -804.1}, {-137, -18.9, -808.9}, {-163, -45.2, -803.1}, {-137, -45.3, -807.9}, {-163, -78.6, -800.1}, {-137, -79.1, -804.9}, {-163, -104.9, -797.2}, {-137, -105.4, -801.9}, {-103.4, -138, -802}, {-102.9, -164, -797.2}, {-77.1, -138, -804.9}, {-76.6, -164, -800}, {-43.3, -139, -807}, {-43.2, -165, -802.1}, {-16.9, -139, -807.9}, {-16.7, -165, -803}, {16.7, -165, -803}, {16.9, -139, -807.9}, {43.2, -165, -802.1}, {43.3, -139, -807}, {76.6, -164, -800}, {77.1, -138, -804.9}, {102.9, -164, -797.2}, {103.4, -138, -802}, {137, -105.4, -801.9}, {163, -104.9, -797.2}, {137, -79.1, -804.9}, {163, -78.6, -800.1}, {137, -45.3, -807.9}, {163, -45.2, -803.1}, {137, -18.9, -808.9}, {163, -18.7, -804.1}};

  HistogramRegistry registry{
    "registry",
    {{"UnMatchedTracksXY", "; #it{x} (cm); #it{y} (cm);", {HistType::kTH2F, {{701, -35.05, 35.05}, {701, -35.05, 35.05}}}},
     {"MatchedTracksXY", "; #it{x} (cm); #it{y} (cm);", {HistType::kTH2F, {{701, -35.05, 35.05}, {701, -35.05, 35.05}}}},
     {"AllTracksXY", "; #it{x} (cm); #it{y} (cm);", {HistType::kTH2F, {{701, -35.05, 35.05}, {701, -35.05, 35.05}}}},
     {"DistChannelToProp", "; D (cm); #count", {HistType::kTH1D, {{101, 0, 100}}}},
     {"NchannelsPerBC", "; N_{channelC}; #count", {HistType::kTH1D, {{101, 0, 100}}}},
     {"NgoodBCperTrack", "; N_{goodBC}; #count", {HistType::kTH1D, {{11, 0, 10}}}},
     {"NgoodBCperTrackINDIV", "; N_{goodBC}; #count", {HistType::kTH1D, {{11, 0, 10}}}},
     {"NCompBCwFT0C", "; N_{compBC}; #count", {HistType::kTH1D, {{21, -0.5, 20.5}}}},
     {"NCompBCwFT0s", "; N_{compBC}; #count", {HistType::kTH1D, {{21, -0.5, 20.5}}}},
     {"DiffInBCINDIV", "; indivBC-firstBC (globalBC); #count", {HistType::kTH1I, {{199, 0, 199}}}},
     {"DiffInBC", "; goodBC-firstBC (globalBC); #count", {HistType::kTH1I, {{199, 0, 199}}}}}};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {

    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current()
              << " A for run " << bc.runNumber()
              << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag); // for some reason this is necessary for the next next line
    runNumber = bc.runNumber();

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());

    Bz = field->getBz(centerMFT); // gives error if the propagator is not initFielded
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz;
  }

  bool isInFT0Acc(double x, double y)
  {
    // returns true if the propagated x and y positions are in an active zone of the FT0-C, false if they in a dead zone

    if ((abs(x) < 6.365) && (abs(y) < 6.555)) {
      // track outside the FT0-C acceptance (in the central hole)
      return false;
    }

    if (((x > -12.75) && (x < -11.85)) || ((x > -6.55) && (x < -5.75)) || ((x > -0.35) && (x < 0.45)) || ((x > 5.75) && (x < 6.55)) || ((x > 11.85) && (x < 12.75))) {
      // track outside the FT0-C acceptance (in the vertical line holes)
      return false;
    }

    if (((y > -12.95) && (y < -11.95)) || ((y > -6.65) && (y < -5.85)) || ((y > -0.55) && (y < 0.45)) || ((y > 5.75) && (y < 6.65)) || ((y > 11.95) && (y < 12.85))) {
      // track outside the FT0-C acceptance (in the horizontal line holes)
      return false;
    }

    return true;
  }

  void processMFT(aod::MFTTracks const& mfttracks,
                  aod::Collisions const&, ExtBCs const& bcs,
                  aod::FT0s const&)
  {
    initCCDB(bcs.begin());

    int i = 0; // counts the number of channels having non-zero amplitude
    // for a particular BC
    double D = 0.0; // distance between (xe,ye,ze) and (xc,yc,zc)
    double minD;
    double globalMinD;

    for (auto& track : mfttracks) {
      filler.BCids.clear();
      globalMinD = 999.;          // minimum D for all BC
      ExtBCs::iterator closestBC; // compatible BC with the D the smallest
      // beware: there could be several BC with the same smallest D
      // not a very useful variable

      if (!track.has_collision()) {
        BcMft(track.globalIndex(), filler.BCids); // empty
        continue;
      }
      auto collOrig = track.collision();

      auto bcSlice = getCompatibleBCs(track, collOrig, bcs, shiftBC);

      // firstBC= global BC of the beginning of the ROF (shifted by shiftBC)
      int64_t firstBC = collOrig.bc_as<ExtBCs>().globalBC() + (track.trackTime() - track.trackTimeRes()) / o2::constants::lhc::LHCBunchSpacingNS + shiftBC;

      bool rofHasBoth = false; // ROF with both FT0C and FT0A signal in the same BC

      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());

      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      // we propagate the MFT track to the mean z position of FT0-C
      // getTrackPar() doesn't work because mft tracks don't have alpha
      trackPar.propagateToZhelix(-82.6, Bz); // z in cm

      if (!isInFT0Acc(trackPar.getX(), trackPar.getY())) {
        // track outside the FT0-C acceptance
        BcMft(track.globalIndex(), filler.BCids); // empty
        continue;
      }

      std::vector<ExtBCs::iterator> goodBC; // contains the BCs matched with the current MFT track
      int nCompBCwft0C = 0;
      int nCompBCwft0s = 0; // Number of compatible BCs with FT0A AND C signals

      bool hasft0A = false;
      bool hasft0C = false;
      for (auto& bc : bcSlice) {
        hasft0C = false;
        hasft0A = false;
        // printf("----------bcId %lld\n", bc.globalIndex());
        if (!bc.has_ft0s()) {
          continue;
        }

        auto ft0s = bc.ft0s();
        i = 0; // reinitialise
        D = 0.0;
        minD = 999.9;
        for (auto const& ft0 : ft0s) {
          // printf("---------ft0.bcId %d\n", ft0.bcId());
          if (ft0.channelA().size() > 0) {
            // BC with signals in FT0A
            hasft0A = true;
          }

          if (ft0.channelC().size() > 0) {
            hasft0C = true;
          }
          for (auto channelId : ft0.channelC()) {

            std::vector<float> Xc = channelCoord[channelId]; //(xc,yc,zc) coordinates
            // D in cm
            D = sqrt(pow(Xc[0] * 0.1 - trackPar.getX(), 2) + pow(Xc[1] * 0.1 - trackPar.getY(), 2) + pow(Xc[2] * 0.1 - 1.87 - trackPar.getZ(), 2));
            // printf("----channelId %u, D %f, n %d\n", channelId, D, n);//should be between 96 and 207
            if (D < minD) {
              minD = D;
            }

            registry.fill(HIST("DistChannelToProp"), D);
          }

          i += ft0.channelC().size();
        }

        registry.fill(HIST("NchannelsPerBC"), i);
        if (hasft0C) {
          nCompBCwft0C++; // number of compatible BCs that have ft0-C signal
        }

        //----------------------- BC selection here ------------------------
        // if strictBCSel true we are only considering BC with signals from both FT0A and FT0C
        if (!(hasft0A && hasft0C) && strictBCSel) {
          continue;
          // we go to the next BC
        }
        nCompBCwft0s++;
        if (hasft0A && hasft0C) {
          rofHasBoth = true;
        }

        //----------------------- end of BC selection ------------------------

        if (minD < 2) // 20 mm
        {
          goodBC.push_back(bc); // goodBC is a vector of bc
          filler.BCids.emplace_back(bc.globalIndex());
        }
        if (minD < globalMinD) {
          globalMinD = minD;
          closestBC = bc;
        }
      }

      if (!rofHasBoth) {
        // there isn't a coincidence of FT0A and C inside the considered MFT ROF
        // MFT track is probably noise, we don't select it
        filler.BCids.clear();
        BcMft(track.globalIndex(), filler.BCids); // empty
        continue;
      }
      registry.fill(HIST("NgoodBCperTrack"), goodBC.size());
      if (goodBC.size() == 0) {
        registry.fill(HIST("UnMatchedTracksXY"), trackPar.getX(), trackPar.getY());
      }
      if (goodBC.size() > 0) {
        registry.fill(HIST("MatchedTracksXY"), trackPar.getX(), trackPar.getY());
        int64_t diff = goodBC[0].globalBC() - firstBC;
        registry.fill(HIST("DiffInBC"), diff);
      }
      registry.fill(HIST("AllTracksXY"), trackPar.getX(), trackPar.getY());
      registry.fill(HIST("NCompBCwFT0C"), nCompBCwft0C);
      registry.fill(HIST("NCompBCwFT0s"), nCompBCwft0s);

      if (nCompBCwft0s == 1) {
        registry.fill(HIST("NgoodBCperTrackINDIV"), goodBC.size());

        // position of the goodBC in the ROF for isolated colliding BCs
        if (goodBC.size() > 0) {
          int64_t diff = goodBC[0].globalBC() - firstBC;
          registry.fill(HIST("DiffInBCINDIV"), diff);
        }
      }

      BcMft(track.globalIndex(), filler.BCids);
    } // loop of mfttracks
  }
  PROCESS_SWITCH(matchmftft0, processMFT, "Process MFT tracks with collisions", true);
};

struct checkmatchinmc {
  // checks if the matching works as expected in MC
  // only doprocessMFTMCcheck==true if you are analysing MC

  HistogramRegistry registryMC{
    "registryMC",
    {}};

  void init(InitContext const&)
  {
    if (doprocessMFTMCcheck) {
      registryMC.add({"TrackIsMatched", "; isMFTTrackMatched; #count", {HistType::kTH1I, {{2, 0, 2}}}});
      registryMC.add({"DiffInBCTrue", "; goodBC-trueBC (globalBC); #count", {HistType::kTH1I, {{800, -400, 400}}}});
      registryMC.add({"TrueBCAmongMatched", "; isTrueBCAmongMatchedOnes; #count", {HistType::kTH1D, {{2, 0, 2}}}});
    }
  }

  using MFTTracksLabeledWithFT0 = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels, aod::BCofMFT>;

  void processMFTMCcheck(MFTTracksLabeledWithFT0 const& mfttracks,
                         aod::McCollisions const&, ExtBCs const&, aod::McParticles const&)
  {

    for (auto& mfttrack : mfttracks) {

      if (!mfttrack.has_bcs()) // mft tracks having a match in FT0-C
      {
        registryMC.fill(HIST("TrackIsMatched"), 0);
        continue;
      }

      registryMC.fill(HIST("TrackIsMatched"), 1); // around 50% of all MFT tracks are matched with FT0-C in data
      // around 90% of MFT tracks falling in the active FT0-C regions are matched
      if (!mfttrack.has_mcParticle()) {
        continue;
      }

      o2::aod::McParticle particle = mfttrack.mcParticle();
      int64_t trueMFTBC = particle.mcCollision().bc_as<ExtBCs>().globalBC();
      ;
      bool isTrueBCAmongMatchedOnes = false;

      for (auto& bc : mfttrack.bcs_as<ExtBCs>()) {      //
        int64_t bcDiffTrue = bc.globalBC() - trueMFTBC; // difference between the muon's BC and the MFT track's BC
        registryMC.fill(HIST("DiffInBCTrue"), bcDiffTrue);
        if (bcDiffTrue == 0) {
          isTrueBCAmongMatchedOnes = true;
        }
      }
      if (isTrueBCAmongMatchedOnes) {
        registryMC.fill(HIST("TrueBCAmongMatched"), 1);
      } else {
        registryMC.fill(HIST("TrueBCAmongMatched"), 0);
      }
    }
  }
  PROCESS_SWITCH(checkmatchinmc, processMFTMCcheck, "Process MFT tracks and check matching with MC information", false);

  void processDummy(aod::Collisions const&)
  {
    // do nothing
  }
  PROCESS_SWITCH(checkmatchinmc, processDummy, "Do nothing if not MC", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<bctoft0c>(cfgc),
                        adaptAnalysisTask<matchmftft0>(cfgc),
                        adaptAnalysisTask<checkmatchinmc>(cfgc)};
  return workflow;
}
