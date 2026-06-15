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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
// Class to define and fill histograms
//

#ifndef PWGDQ_CORE_MIXINGHANDLER_H_
#define PWGDQ_CORE_MIXINGHANDLER_H_

#include "PWGDQ/Core/VarManager.h"

#include <TArrayF.h>
#include <TNamed.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <vector>

class MixingHandler : public TNamed
{

 public:
  // Struct to define track properties relevant for mixing and few utility functions
  struct MixingTrack {
    float pt;
    float eta;
    float phi;
    uint32_t filteringFlags;
    // Clear a bit once the track was used in mixing for that bit for the required pool depth.
    void ClearBit(uint32_t mask) { filteringFlags &= ~mask; }
    void Print() const
    {
      std::cout << "pt: " << pt << ", eta: " << eta << ", phi: " << phi << ", filteringFlags: " << filteringFlags << std::endl;
    }
  };

  // Struct to define events used in mixing and few utility functions.
  // An event is defined as two vectors of tracks (typically the legs of a two-body
  // decay or the two-particles in a correlation analysis)
  struct MixingEvent {
    std::vector<MixingTrack> tracks1;
    std::vector<MixingTrack> tracks2;
    // bit map for active filtering bits of all the tracks
    uint32_t filteringMask = 0;
    // counters to keep track of how many times the event was used for mixing (for each track cut separately)
    std::array<short, 64> counters = {0};
    // add a track to the event and update the filtering mask accordingly
    void AddTrack1(const MixingTrack& track)
    {
      tracks1.push_back(track);
      filteringMask |= track.filteringFlags;
    }
    void AddTrack2(const MixingTrack& track)
    {
      tracks2.push_back(track);
      filteringMask |= track.filteringFlags;
    }
    // Clear bits in the filtering mask.
    void ClearFilteringMask(uint32_t mask) { filteringMask &= ~mask; }
    // 1) increment the counters for a given track cut bit mask and if the counters reached the pool depth,
    // 2) clear the corresponding bit in the tracks filtering flags to exclude them from further mixing
    // 3) for each track, if there are no more active bits in the filtering mask, then remove the track from the event
    void IncrementCounters(uint32_t mask, short poolDepth)
    {
      for (int i = 0; i < 32; i++) {
        uint32_t bitMask = (static_cast<uint32_t>(1) << i);
        if (mask & bitMask) {
          counters[i]++;
          if (counters[i] >= poolDepth) {
            auto clearExpiredBit = [bitMask](auto& tracks) {
              for (auto& track : tracks) {
                track.ClearBit(bitMask);
              }
              tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                                          [](const MixingTrack& track) { return track.filteringFlags == 0; }),
                           tracks.end());
            };
            clearExpiredBit(tracks1);
            clearExpiredBit(tracks2);
            ClearFilteringMask(bitMask);
          }
        }
      }
    }
    void Print() const
    {
      std::cout << "Event filtering mask: ";
      for (int i = 0; i < 32; i++) {
        if (filteringMask & (1ULL << i)) {
          std::cout << "1";
        } else {
          std::cout << "0";
        }
      }
      std::cout << std::endl;
      for (int i = 0; i < 32; i++) {
        if (filteringMask & (1ULL << i)) {
          std::cout << "Counter " << i << ": " << counters[i] << std::endl;
        }
      }
      std::cout << "Tracks 1: " << std::endl;
      for (const auto& track : tracks1) {
        track.Print();
      }
      std::cout << "Tracks 2: " << std::endl;
      for (const auto& track : tracks2) {
        track.Print();
      }
    }
  };

  struct MixingPool {
    std::vector<MixingEvent> events;

    // check which events in the pool are empty (i.e. no active tracks for mixing) and remove them from the pool
    void CleanPool()
    {
      events.erase(std::remove_if(events.begin(), events.end(),
                                  [](const MixingEvent& event) { return event.tracks1.empty() && event.tracks2.empty(); }),
                   events.end());
    }
    // The function that performs the mixing is called outside this class, but the pool provides the events and tracks to be mixed and takes care of updating the events after mixing
    // (e.g. incrementing the counters and removing the tracks that reached the pool depth for a given cut)
    void UpdatePool(const MixingEvent& event, short poolDepth)
    {
      for (auto& event : events) {
        event.IncrementCounters(event.filteringMask, poolDepth);
      }
      CleanPool();
      events.push_back(event);
    }
    // getter for the events in the pool
    const std::vector<MixingEvent>& GetEvents() const { return events; }

    void Print() const
    {
      std::cout << "Mixing pool with " << events.size() << " events:" << std::endl;
      for (const auto& event : events) {
        event.Print();
      }
    }
  };

  MixingHandler();
  MixingHandler(const char* name, const char* title);
  virtual ~MixingHandler();

  // setters
  void AddMixingVariable(int var, std::vector<float> binLims);
  void SetPoolDepth(short depth) { fPoolDepth = depth; }

  // getters
  // int GetNMixingVariables() const { return fVariables.size(); }
  // int GetMixingVariable(VarManager::Variables var); // returns the position in the internal varible list of the handler. Useful for checks, mostly
  // std::vector<float> GetMixingVariableLimits(VarManager::Variables var);
  MixingPool& GetPool(int category) { return fPools[category]; }
  short GetPoolDepth() const { return fPoolDepth; }

  void Init();
  int FindEventCategory(float* values);
  int GetBinFromCategory(VarManager::Variables var, int category) const;

 private:
  MixingHandler(const MixingHandler& handler);
  MixingHandler& operator=(const MixingHandler& handler);

  // User options
  bool fIsInitialized; // check if the mixing handler is initialized

  // bin limits for the variables used for mixing, the number of vectors corresponds to the number of variables and the content of each vector corresponds to the bin limits for that variable
  std::vector<std::vector<float>> fVariableLimits;
  std::map<int, int> fVariables; // key: variable, value: position in the internal variable list of the handler (used to map the variables to the values passed to FindEventCategory)

  short fPoolDepth;                 // number of events to be kept in each pool
  std::map<int, MixingPool> fPools; // key: category, value: pool of events corresponding to that category

  ClassDef(MixingHandler, 2);
};

#endif // PWGDQ_CORE_MIXINGHANDLER_H_
