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

#ifndef PWGUD_CORE_UPCCUTPARHOLDER_H_
#define PWGUD_CORE_UPCCUTPARHOLDER_H_

#include <Rtypes.h>

// object to hold customizable cut values
class UPCCutparHolder
{
 public:
  // constructor
  UPCCutparHolder(bool useFwdCuts = true,
                  int trackType = 3,
                  float fwdPtLow = 0.5,
                  float fwdPtHigh = 4.,
                  float fwdEtaLow = -4.0,
                  float fwdEtaHigh = -2.5,
                  float muonRAtAbsorberEndLow = 17.6,
                  float muonRAtAbsorberEndHigh = 89.5,
                  float muonPDcaHighFirst = 594.0,
                  float muonPDcaHighSecond = 324.0,
                  float fwdChi2Low = 0.0,
                  float fwdChi2High = 10000.0,
                  bool useBarCuts = true,
                  float barPtLow = 0.,
                  float barPtHigh = 1000.,
                  float barEtaLow = -0.9,
                  float barEtaHigh = 0.9,
                  int ITSNClusLow = 4,
                  int ITSNClusHigh = 9,
                  float ITSChi2Low = 0.,
                  float ITSChi2High = 5.,
                  int TPCNClusCRLow = 70,
                  int TPCNClusCRHigh = 161,
                  float TPCChi2Low = 0.,
                  float TPCChi2High = 4.,
                  bool checkMaxDcaXY = true,
                  float dcaZLow = -3.,
                  float dcaZHigh = 3.,
                  bool requireTOF = false,
                  bool requireITSTPC = false,
                  int maxNContrib = 2,
                  int ambigSwitch = 0)
    : fUseFwdCuts{useFwdCuts},
      fTrackType{trackType},
      fFwdPtLow{fwdPtLow},
      fFwdPtHigh{fwdPtHigh},
      fFwdEtaLow{fwdEtaLow},
      fFwdEtaHigh{fwdEtaHigh},
      fMuonRAtAbsorberEndLow{muonRAtAbsorberEndLow},
      fMuonRAtAbsorberEndHigh{muonRAtAbsorberEndHigh},
      fMuonPDcaHighFirst{muonPDcaHighFirst},
      fMuonPDcaHighSecond{muonPDcaHighSecond},
      fFwdChi2Low{fwdChi2Low},
      fFwdChi2High{fwdChi2High},
      fUseBarCuts{useBarCuts},
      fBarPtLow{barPtLow},
      fBarPtHigh{barPtHigh},
      fBarEtaLow{barEtaLow},
      fBarEtaHigh{barEtaHigh},
      fITSNClusLow{ITSNClusLow},
      fITSNClusHigh{ITSNClusHigh},
      fITSChi2Low{ITSChi2Low},
      fITSChi2High{ITSChi2High},
      fTPCNClusCRLow{TPCNClusCRLow},
      fTPCNClusCRHigh{TPCNClusCRHigh},
      fTPCChi2Low{TPCChi2Low},
      fTPCChi2High{TPCChi2High},
      fCheckMaxDcaXY{checkMaxDcaXY},
      fDcaZLow{dcaZLow},
      fDcaZHigh{dcaZHigh},
      fRequireTOF{requireTOF},
      fRequireITSTPC{requireITSTPC},
      fMaxNContrib{maxNContrib},
      fAmbigSwitch{ambigSwitch} {}

 public:
  // setters
  void setUseFwdCuts(bool useFwdCuts);
  void setTrackType(int trackType);
  void setFwdPtLow(float fwdPtLow);
  void setFwdPtHigh(float fwdPtHigh);
  void setFwdEtaLow(float fwdEtaLow);
  void setFwdEtaHigh(float fwdEtaHigh);
  void setMuonRAtAbsorberEndLow(float muonRAtAbsorberEndLow);
  void setMuonRAtAbsorberEndHigh(float muonRAtAbsorberEndHigh);
  void setMuonPDcaHighFirst(float muonPDcaHighFirst);
  void setMuonPDcaHighSecond(float muonPDcaHighSecond);
  void setFwdChi2Low(float fwdChi2Low);
  void setFwdChi2High(float fwdChi2High);
  void setUseBarCuts(bool useBarCuts);
  void setBarPtLow(float barPtLow);
  void setBarPtHigh(float barPtHigh);
  void setBarEtaLow(float barEtaLow);
  void setBarEtaHigh(float barEtaHigh);
  void setITSNClusLow(int ITSNClusLow);
  void setITSNClusHigh(int ITSNClusHigh);
  void setITSChi2Low(float ITSChi2Low);
  void setITSChi2High(float ITSChi2High);
  void setTPCNClusCRLow(int TPCNClusCRLow);
  void setTPCNClusCRHigh(int TPCNClusCRHigh);
  void setTPCChi2Low(float TPCChi2Low);
  void setTPCChi2High(float TPCChi2High);
  void setCheckMaxDcaXY(bool checkMaxDcaXY);
  void setDcaZLow(float dcaZLow);
  void setDcaZHigh(float dcaZHigh);
  void setRequireTOF(bool requireTOF);
  void setRequireITSTPC(bool requireITSTPC);
  void setProduceITSITS(bool produceITSITS);
  void setMaxNContrib(int maxNContrib);
  void setAmbigSwitch(int ambigSwitch);

  // getters
  bool getUseFwdCuts() const;
  int getTrackType() const;
  float getFwdPtLow() const;
  float getFwdPtHigh() const;
  float getFwdEtaLow() const;
  float getFwdEtaHigh() const;
  float getMuonRAtAbsorberEndLow() const;
  float getMuonRAtAbsorberEndHigh() const;
  float getMuonPDcaHighFirst() const;
  float getMuonPDcaHighSecond() const;
  float getFwdChi2Low() const;
  float getFwdChi2High() const;
  bool getUseBarCuts() const;
  float getBarPtLow() const;
  float getBarPtHigh() const;
  float getBarEtaLow() const;
  float getBarEtaHigh() const;
  int getITSNClusLow() const;
  int getITSNClusHigh() const;
  float getITSChi2Low() const;
  float getITSChi2High() const;
  int getTPCNClusCRLow() const;
  int getTPCNClusCRHigh() const;
  float getTPCChi2Low() const;
  float getTPCChi2High() const;
  bool getCheckMaxDcaXY() const;
  float getDcaZLow() const;
  float getDcaZHigh() const;
  bool getRequireTOF() const;
  bool getRequireITSTPC() const;
  bool getProduceITSITS() const;
  int getMaxNContrib() const;
  int getAmbigSwitch() const;

 private:
  bool fUseFwdCuts{true}; // Use cuts for forward tracks

  // Filter by Fwd. track type:
  //   -1 -> no filter,
  //   0  -> MFT-MCH-MID,
  //   2  -> MFT-MCH,
  //   3  -> MCH-MID,
  //   4  -> MCH
  //   5  -> MCH-MID and MCH
  // See ForwardTrackTypeEnum
  int fTrackType{3};

  // basic
  float fFwdPtLow{0.5};    // Minimal Pt for forward tracks
  float fFwdPtHigh{4.};    // Maximal Pt for forward tracks
  float fFwdEtaLow{-4.0};  // Maximal Eta for forward tracks
  float fFwdEtaHigh{-2.5}; // Maximal Eta for forward tracks

  // quality
  float fMuonRAtAbsorberEndLow{17.6};  // "Minimal muon R at absorber end
  float fMuonRAtAbsorberEndHigh{89.5}; // "Maximal muon R at absorber end
  float fMuonPDcaHighFirst{594.0};     // "Primary PDCA cut: Maximal value for R < 26.5
  float fMuonPDcaHighSecond{324.0};    // "Additional PDCA cut: Maximal value for R >= 26.5
  float fFwdChi2Low{0.0};              // "Minimal Chi2 for forward tracks
  float fFwdChi2High{10000.0};         // "Maximal Chi2 for forward tracks

  // cuts for central-barrel tracks
  bool fUseBarCuts{true}; // Use cuts for barrel tracks
  // basic
  float fBarPtLow{0.};     // Minimal Pt for barrel tracks
  float fBarPtHigh{1000.}; // Maximal Pt for barrel tracks
  float fBarEtaLow{-0.9};  // Maximal Eta for barrel tracks
  float fBarEtaHigh{0.9};  // Maximal Eta for barrel tracks
  // quality: ITS
  int fITSNClusLow{4};    // Minimal number of ITS clusters
  int fITSNClusHigh{9};   // Maximal number of ITS clusters
  float fITSChi2Low{0.};  // Minimal Chi2 in ITS per cluster
  float fITSChi2High{5.}; // Maximal Chi2 in ITS per cluster
  // quality: TPC
  int fTPCNClusCRLow{70};   // Minimal number of TPC clusters (crossed rows)
  int fTPCNClusCRHigh{161}; // Maximal number of TPC clusters (crossed rows)
  float fTPCChi2Low{0.};    // Minimal Chi2 in TPC per cluster
  float fTPCChi2High{4.};   // Maximal Chi2 in TPC per cluster
  // quality: DCA
  bool fCheckMaxDcaXY{true}; // Apply cut on maximal DCA_xy
  float fDcaZLow{-3.};       // Minimal DCA_z for barrel tracks
  float fDcaZHigh{3.};       // Maximal DCA_z for barrel tracks
  // quality: matching
  bool fRequireTOF{false};    // Require all tracks in event candidates to have TOF matches
  bool fRequireITSTPC{false}; // Require all tracks in event candidates to have ITS-TPC matches
  bool fProduceITSITS{false}; // Produce candidates using only ITS-TPC tracks as well

  // tracks from collisions: consider only tracks from collisions with N tracks less or equal than fMaxNContrib
  int fMaxNContrib{2}; // Central barrel: consider tracks from collisions with N contributors <= maxNContrib
  int fAmbigSwitch{0}; // Central barrel: 0 -- loop over all tracks, 1 -- loop only over tracks with vertices

  ClassDefNV(UPCCutparHolder, 1);
};

#endif // PWGUD_CORE_UPCCUTPARHOLDER_H_
