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
#ifndef PWGCF_CORE_ANALYSISCONFIGURABLECUTS_H_
#define PWGCF_CORE_ANALYSISCONFIGURABLECUTS_H_

#include <string>
#include <vector>
#include <Rtypes.h>
#include <TObject.h>
#include <TNamed.h>
#include <TMath.h>

namespace o2
{
namespace analysis
{
/// \class EventSelectionCuts
/// \brief Class which implements configurable event selection cuts
///
class EventSelectionCuts
{
 public:
  int mOfflinetrigger = 1;                     ///< offline trigger, default MB = 1
  std::string mCentmultestmator = "V0M";       ///< centrality / multiplicity estimation, default V0M
  int mRemovepileupcode = 1;                   ///< Procedure for pile-up removal, default V0M vs TPCout tracks = 1
  std::string mRemovepileupfn = "-2500+5.0*x"; ///< function for pile-up removal, procedure dependent, defaul V0M vs TPCout tracks for LHC15o HIR

 private:
  ClassDefNV(EventSelectionCuts, 1);
};

/// \class DptDptBinningCuts
/// \brief Class which implements configurable acceptance cuts
///
class DptDptBinningCuts
{
 public:
  int mZVtxbins = 28;       ///< the number of z_vtx bins default 28
  float mZVtxmin = -7.0;    ///< the minimum z_vtx value, default -7.0 cm
  float mZVtxmax = 7.0;     ///< the maximum z_vtx value, default 7.0 cm
  int mPTbins = 18;         ///< the number of pT bins, default 18
  float mPTmin = 0.2;       ///< the minimum pT value, default 0.2 GeV
  float mPTmax = 2.0;       ///< the maximum pT value, default 2.0 GeV
  int mEtabins = 16;        ///< the number of eta bins default 16
  float mEtamin = -0.8;     ///< the minimum eta value, default -0.8
  float mEtamax = 0.8;      ///< the maximum eta value, default 0.8
  int mPhibins = 72;        ///< the number of phi bins, default 72
  float mPhibinshift = 0.5; ///< the shift in the azimuthal origen, defoult 0.5, i.e half a bin

 private:
  ClassDefNV(DptDptBinningCuts, 1);
};

/// \brief Simple class for a generic check within a concrete range of a magnitude
class CheckRangeCfg
{
 public:
  bool mDoIt = false;    ///< do the actual check
  float mLowValue = 0.0; ///< range lowest value
  float mUpValue = 0.0;  ///< range upper value
 private:
  ClassDefNV(CheckRangeCfg, 1);
};

/// \brief Simple class for configuring a track selection object
class TrackSelectionCfg
{
 public:
  bool mUseIt = false;          ///< use this track selection configuration
  bool mOnGen = false;          ///< apply it to generator level also
  int mTPCclusters = 0;         ///< minimum number of TPC clusters
  int mTPCxRows = 70;           ///< minimum number of TPC crossed rows
  float mTPCXRoFClusters = 0.8; ///< minimu value of the TPC ratio no of crossed rows over findable clusters
  float mDCAxy = 2.4;           ///< maximum DCA on xy plane
  float mDCAz = 3.2;            ///< maximum DCA on z axis

 private:
  ClassDefNV(TrackSelectionCfg, 1);
};

/// \brief Simple class for configuring the fine tuning a track selection object
/// The tune should change the track selection objects and probably do further checks
/// after the actual track selection has ben performed
class TrackSelectionTuneCfg
{
 public:
  bool mUseIt = false;                        ///< use this track selection tuning configuration
  int mTPCclusters = 0;                       ///< minimum number of TPC clusters
  bool mUseTPCclusters = false;               ///< use or not the number of TPC clusters
  int mTPCxRows = 70;                         ///< minimum number of TPC crossed rows
  bool mUseTPCxRows = false;                  ///< use or not the number of TPC crossed rows
  float mTPCXRoFClusters = 0.8;               ///< minimum value of the TPC ratio no of crossed rows over findable clusters
  bool mUseTPCXRoFClusters = false;           ///< use or not the TPC ratio of no of crossed rows over findable clusters
  float mFractionTpcSharedClusters = 0.4;     ///< the maximum fraction of TPC shared clusters
  bool mUseFractionTpcSharedClusters = false; ///< use or not the fraction of TPC share clusters
  float mDCAxy = 2.4;                         ///< maximum DCA on xy plane
  bool mUseDCAxy = false;                     ///< use or not the maximum DCA on the xy plane
  float mDCAz = 3.2;                          ///< maximum DCA on z axis
  bool mUseDCAz = false;                      ///< use or not the maximum DCA on z asis

 private:
  ClassDefNV(TrackSelectionTuneCfg, 2);
};

/// \brief Simple class to configure a selection based on PID
/// The nsigmas information is for closeness to the line of interest
/// and for separation to the other lines
class TrackSelectionPIDCfg
{
 public:
  bool mUseIt = false;
  std::vector<float> mMinNSigmasTPC = {0.0f, 0.0f, -3.0f, -3.0f, -3.0f}; ///< nsigmas TPC lower limit for e, mu, pi, Ka, and p
  std::vector<float> mMaxNSigmasTPC = {0.0f, 0.0f, 3.0f, 3.0f, 3.0f};    ///< nsigmas TPC upper limit for e, mu, pi, Ka, and p
  float mPThreshold = 0.0;                                               ///< momentum threshold for considering TOF information
  bool mRequireTOF = true;                                               ///< require or not the presence of TOF when the momentum threshold is passed
  std::vector<float> mMinNSigmasTOF = {0.0f, 0.0f, -3.0f, -3.0f, -3.0f}; ///< nsigmas TOF lower limit for e, mu, pi, Ka, and p
  std::vector<float> mMaxNSigmasTOF = {0.0f, 0.0f, 3.0f, 3.0f, 3.0f};    ///< nsigmas TOF upper limit for e, mu, pi, Ka, and p
  bool m2Dcut = true;                                                    ///< use an elliptic cut using TPC and TOF nsigmas
  bool mExclude = false;                                                 ///< should the identified track be excluded for analysis?
  float mPtMin = 0.2;                                                    ///< increase the lower pT limit for this species
  float mPtMax = 2.0;                                                    ///< decrease the upper pT limit for this species
 private:
  ClassDefNV(TrackSelectionPIDCfg, 2);
};

class SimpleInclusiveCut : public TNamed
{
 public:
  int mX = 1;
  float mY = 2.f;
  SimpleInclusiveCut();
  SimpleInclusiveCut(const char*, int, float);
  SimpleInclusiveCut(const SimpleInclusiveCut&) = default;
  ~SimpleInclusiveCut() override = default;

  SimpleInclusiveCut& operator=(const SimpleInclusiveCut&);

 private:
  ClassDef(SimpleInclusiveCut, 1);
};

} // namespace analysis
} // namespace o2
#endif // PWGCF_CORE_ANALYSISCONFIGURABLECUTS_H_
