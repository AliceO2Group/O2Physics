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

/// \file DhCorrelationExtraction.h
/// \brief Class for D-h correlation extraction
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#ifndef PWGHF_HFC_MACROS_DHCORRELATIONEXTRACTION_H_
#define PWGHF_HFC_MACROS_DHCORRELATIONEXTRACTION_H_

#include <Rtypes.h>
#include <RtypesCore.h>
#include <TAttMarker.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TObject.h>
#include <TString.h>

class DhCorrelationExtraction : public TObject
{

 public:
  enum DmesonSpecie { kD0toKpi = 0,
                      kDplusKpipi,
                      kDsToKKPi,
                      kDStarD0pi };
  enum selectAnalysisType { kSE,
                            kME };
  enum selectInvMassRegion { kSign,
                             kSideb };

  DhCorrelationExtraction(); // default constructor
  DhCorrelationExtraction(const DhCorrelationExtraction& source);
  virtual ~DhCorrelationExtraction();

  /// Methods to set the input configuration
  // Input files, directories and histograms
  Bool_t SetDmesonSpecie(DmesonSpecie k);
  void SetInputFilenameMass(TString filenameMass) { fFileNameMass = filenameMass; }
  void SetInputFilenameSE(TString filenameSE) { fFileNameSE = filenameSE; }
  void SetInputFilenameME(TString filenameME) { fFileNameME = filenameME; }
  void SetDirNameSE(TString dirNameSE) { fDirNameSE = dirNameSE; }
  void SetDirNameME(TString dirNameME) { fDirNameME = dirNameME; }
  void SetMassHistoNameSgn(TString massHistoNameSgn) { fMassHistoNameSgn = massHistoNameSgn; }
  void SetMassHistoNameBkg(TString massHistoNameBkg) { fMassHistoNameBkg = massHistoNameBkg; }
  void SetMassHistoNameSBs(TString massHistoNameSBs) { fMassHistoNameSBs = massHistoNameSBs; }
  void SetSECorrelHistoSignalName(TString correlNameSigSE) { fSECorrelSignalRegionName = correlNameSigSE; }
  void SetSECorrelHistoSidebandName(TString correlNameSbSE) { fSECorrelSidebandsName = correlNameSbSE; }
  void SetMECorrelHistoSignalName(TString correlNameSigME) { fMECorrelSignalRegionName = correlNameSigME; }
  void SetMECorrelHistoSidebandName(TString correlNameSbME) { fMECorrelSidebandsName = correlNameSbME; }

  // Input conditions: PtCand, PtHad, PoolBins
  void SetNpools(Int_t npools) { fNpools = npools; }
  void SetCorrectPoolsSeparately(Bool_t usePools) { fCorrectPoolsSeparately = usePools; }
  void SetDeltaEtaRange(Double_t etaLow = -1., Double_t etaHigh = 1)
  {
    fDeltaEtaMin = etaLow;
    fDeltaEtaMax = etaHigh;
  }
  void SetSubtractSoftPiInMEdistr(Bool_t subtractSoftPiME) { fSubtractSoftPiME = subtractSoftPiME; }
  void SetBkgScaleFactor(Double_t scaleFactor) { fBkgScaleFactor = scaleFactor; }
  void SetSignalYieldforNorm(Double_t sgnYield) { fSgnYieldNorm = sgnYield; }
  void SetRebin2DcorrelHisto(Int_t rebinDeltaEta, Int_t rebinDeltaPhi)
  {
    fRebin2Dhisto = kTRUE;
    fRebinAxisDeltaEta = rebinDeltaEta;
    fRebinAxisDeltaPhi = rebinDeltaPhi;
  }
  void GetSignalAndBackgroundForNorm(Double_t PtCandMin, Double_t PtCandMax);
  void NormalizeMEplot(TH2D*& histoME, TH2D*& histoMEsoftPi);
  void SetDebugLevel(Int_t debug) { fDebug = debug; }

  /// Analysis methods
  TH2D* GetCorrelHisto(Int_t SEorME, Int_t SorSB, Int_t pool, Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax);
  Bool_t ReadInputSEandME();
  Bool_t ReadInputInvMass();
  Bool_t ExtractCorrelations(Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax, TString codeName);
  TH1D* GetCorrectedCorrHisto() { return fCorrectedCorrHisto; }

  /// Histogram style
  void SetTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, Style_t markerStyle = kFullCircle, Color_t markerColor = kRed + 1, Double_t markerSize = 1.4, Color_t lineColor = kRed + 1, Int_t lineWidth = 3, Float_t hTitleXaxisOffset = 1.0, Float_t hTitleYaxisOffset = 1.0, Float_t hTitleXaxisSize = 0.060, Float_t hTitleYaxisSize = 0.060, Float_t hLabelXaxisSize = 0.060, Float_t hLabelYaxisSize = 0.060, Bool_t centerXaxisTitle = false, Bool_t centerYaxisTitle = false);
  void SetTH2HistoStyle(TH2D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, TString hZaxisTitle, Float_t hTitleXaxisOffset = 1.8, Float_t hTitleYaxisOffset = 1.8, Float_t hTitleZaxisOffset = 1.2, Float_t hTitleXaxisSize = 0.060, Float_t hTitleYaxisSize = 0.060, Float_t hTitleZaxisSize = 0.060, Float_t hLabelXaxisSize = 0.060, Float_t hLabelYaxisSize = 0.060, Float_t hLabelZaxisSize = 0.060, Bool_t centerXaxisTitle = true, Bool_t centerYaxisTitle = true);

 private:
  TFile* fFileMass; // file containing the mass histograms
  TFile* fFileSE;   // file containing the Same Event (SE) output
  TFile* fFileME;   // file containing the Mixed Event (ME) output

  TDirectoryFile* fDirMass; // TDirectory for mass histos
  TDirectoryFile* fDirSE;   // TDirectory for SE info
  TDirectoryFile* fDirME;   // TDirectory for ME info

  TH1D* fCorrectedCorrHisto; // Corrected correlation histogram

  DmesonSpecie fDmesonSpecies;       // D meson specie
  TString fDmesonLabel;              // D meson label
  TString fFileNameMass;             // File cntaining inv. mass histograms
  TString fFileNameSE;               // File contaning Same Event (SE) output
  TString fFileNameME;               // File contaning Mixed Event (ME) output
  TString fDirNameSE;                // Directory in the file containing SE output
  TString fDirNameME;                // Directory in the file containing ME output
  TString fMassHistoNameSgn;         // Inv. mass histo name signal yield
  TString fMassHistoNameBkg;         // Inv. mass histo name background yield
  TString fMassHistoNameSBs;         // Inv. mass histo name sideband yield
  TString fSECorrelSignalRegionName; // THnSparse name containing SE output for signal region
  TString fSECorrelSidebandsName;    // THnSparse name containing SE output for sideband region
  TString fMECorrelSignalRegionName; // THnSparse name containing ME output for signal region
  TString fMECorrelSidebandsName;    // THnSparse name containing ME output for sideband region

  Int_t fNpools;            // number of pools used for the ME correction
  Int_t fRebinAxisDeltaEta; // rebin deltaEta axis
  Int_t fRebinAxisDeltaPhi; // rebin deltaPhi axis
  Int_t fDebug;             // debug level

  Double_t fDeltaEtaMin;    // deltaEta min value
  Double_t fDeltaEtaMax;    // deltaEta max value
  Double_t fBkgScaleFactor; // Bkg/SB factor to scale correlation plots obtained in the sideband region
  Double_t fSgnYieldNorm;   // Signal yield (used for normalize correlation plots after bkg subtraction)

  Bool_t fCorrectPoolsSeparately; // Possibility to do the ME correction pool-by-pool (kTRUE) or merging all pools (kFALSE)
  Bool_t fSubtractSoftPiME;       // Soft pion subtraction (for D0 case)
  Bool_t fRebin2Dhisto;           // Flag to rebin the 2D correlation plots
};

#endif // PWGHF_HFC_MACROS_DHCORRELATIONEXTRACTION_H_
