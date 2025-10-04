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

#include <TAttMarker.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TObject.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

class DhCorrelationExtraction : public TObject
{

 public:
  enum DmesonSpecie { kD0toKpi = 0,
                      kDplusKpipi,
                      kDsToKKPi,
                      kDStarD0pi };
  enum SelectAnalysisType { kSE,
                            kME };
  enum SelectInvMassRegion { kSign,
                             kSideb };
  enum SelectDmesonOrigin { kPrompt,
                            kFD };
  enum SelectParticleType { kPrimaryPart,
                            kAllPart };

  DhCorrelationExtraction(); // default constructor
  DhCorrelationExtraction(const DhCorrelationExtraction& source);
  ~DhCorrelationExtraction() override;

  /// Methods to set the input configuration
  // Input files, directories and histograms
  Bool_t setDmesonSpecie(DmesonSpecie k);
  void setInputFilenameMass(TString filenameMass) { fFileNameMass = filenameMass; }
  void setInputFilenameSe(TString filenameSE) { fFileNameSE = filenameSE; }
  void setInputFilenameMe(TString filenameME) { fFileNameME = filenameME; }
  void setInputFilenameSecPart(TString filenameSecPart) { fFileSecPartName = filenameSecPart; }
  void setInputFilenameBiasBtoD(TString filenamePromptMcRec, TString filenameNonPromptMcRec)
  {
    fFilePromptMcRecName = filenamePromptMcRec;
    fFileNonPromptMcRecName = filenameNonPromptMcRec;
  }
  void setDirNameSe(TString dirNameSE) { fDirNameSE = dirNameSE; }
  void setDirNameMe(TString dirNameME) { fDirNameME = dirNameME; }
  void setDirNameSecPart(TString dirNameSecPart) { fDirSecPartName = dirNameSecPart; }
  void setMassHistoNameSgn(TString massHistoNameSgn) { fMassHistoNameSgn = massHistoNameSgn; }
  void setMassHistoNameBkg(TString massHistoNameBkg) { fMassHistoNameBkg = massHistoNameBkg; }
  void setMassHistoNameSBs(TString massHistoNameSBs) { fMassHistoNameSBs = massHistoNameSBs; }
  void setSeCorrelHistoSignalName(TString correlNameSigSE) { fSECorrelSignalRegionName = correlNameSigSE; }
  void setSeCorrelHistoSidebandName(TString correlNameSbSE) { fSECorrelSidebandsName = correlNameSbSE; }
  void setSeCorrelHistoSidebandLeftName(TString correlNameSbSE) { fSECorrelSidebandLeftName = correlNameSbSE; }
  void setSeCorrelHistoSidebandRightName(TString correlNameSbSE) { fSECorrelSidebandRightName = correlNameSbSE; }
  void setMeCorrelHistoSignalName(TString correlNameSigME) { fMECorrelSignalRegionName = correlNameSigME; }
  void setMeCorrelHistoSidebandName(TString correlNameSbME) { fMECorrelSidebandsName = correlNameSbME; }
  void setMeCorrelHistoSidebandLeftName(TString correlNameSbME) { fMECorrelSidebandLeftName = correlNameSbME; }
  void setMeCorrelHistoSidebandRightName(TString correlNameSbME) { fMECorrelSidebandRightName = correlNameSbME; }
  void setHistoSecPartName(TString histoPrimaryPartName, TString histoAllPartName)
  {
    fHistoPrimaryPartName = histoPrimaryPartName;
    fHistoAllPartName = histoAllPartName;
  }
  void setInputFilenameFdTemplate(TString filenameFDTemplate) { fFileFDTemplateName = filenameFDTemplate; }
  void setInputFilenameFdPromptFrac(TString filenameFDPromptFrac) { fFileFDPromptFracName = filenameFDPromptFrac; }
  void setInputHistoNameFdTemplatePrompt(TString hNameFDTemplatePrompt) { fHistoFDTemplatePromptName = hNameFDTemplatePrompt; }
  void setInputHistoNameFdTemplateNonPrompt(TString hNameFDTemplateNonPrompt) { fHistoFDTemplateNonPromptName = hNameFDTemplateNonPrompt; }
  void setInputHistoNameFdPromptFrac(TString hNameFDPromptFrac) { fHistoFDPromptFracName = hNameFDPromptFrac; }

  // Input conditions: PtCand, PtHad, PoolBins
  void setNpools(Int_t npools) { fNpools = npools; }
  void setCorrectPoolsSeparately(Bool_t usePools) { fCorrectPoolsSeparately = usePools; }
  void setDeltaEtaRange(Double_t etaLow = -1., Double_t etaHigh = 1)
  {
    fDeltaEtaMin = etaLow;
    fDeltaEtaMax = etaHigh;
  }
  void setSubtractSoftPiInMEdistr(Bool_t subtractSoftPiME) { fSubtractSoftPiME = subtractSoftPiME; }
  void setBkgScaleFactor(Double_t scaleFactor) { fBkgScaleFactor = scaleFactor; }
  void setSignalYieldforNorm(Double_t sgnYield) { fSgnYieldNorm = sgnYield; }
  void setBkgYield(Double_t bkgYield) { fBkgYield = bkgYield; }
  void setSbYield(Double_t sbYield) { fSBYield = sbYield; }
  void setRebin2DcorrelHisto(Int_t rebinDeltaEta, Int_t rebinDeltaPhi)
  {
    fRebinAxisDeltaEta = rebinDeltaEta;
    fRebinAxisDeltaPhi = rebinDeltaPhi;
  }
  void setRebinOptions(Bool_t rebinAngCorr, Bool_t rebinFDCorr, Bool_t rebinSecPart)
  {
    fRebinAngCorr = rebinAngCorr;
    fRebinFDCorr = rebinFDCorr;
    fRebinSecPart = rebinSecPart;
  }
  void getSignalAndBackgroundForNorm(Double_t ptCandMin, Double_t ptCandMax);
  void normalizeMePlot(TH2D*& histoME, TH2D*& histoMEsoftPi) const;
  void setDebugLevel(Int_t debug) { fDebug = debug; }
  void setDividedSidebands(Bool_t dividedSideb, Bool_t useSidebLeft, Bool_t useSidebRight)
  {
    fSidebandDivided = dividedSideb;
    fUseSidebLeft = useSidebLeft;
    fUseSidebRight = useSidebRight;
  }
  void setFdSubtraction(Bool_t subtractFD) { fFDsubtraction = subtractFD; }
  void setSecPartContamination(Bool_t secPartContamination) { fSecPartContamination = secPartContamination; }
  void setCorrBiasBtoD(Bool_t corrbiasBtoD) { fCorrBiasBtoD = corrbiasBtoD; }
  void setBinCandAndHad(Int_t binCand, Int_t binHad)
  {
    fBinPtCand = binCand;
    fBinPtHad = binHad;
  }

  /// Analysis methods
  TH2D* getCorrelHisto(Int_t sEorMe, Int_t sorSb, Int_t pool, Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax);
  TH2D* getFdTemplateHisto(Int_t promptOrFd, Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax);
  TH1D* getCorrelHistoSecondaryPart(Int_t primaryPart, Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax);
  TH1D* reflectCorrHistogram(TH1D*& histo);
  TH1D* reflectHistoRun2(TH1D* h, Double_t scale);
  TH1D* evaluateMcClosModulations(Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax);
  Double_t getFdPromptFrac(Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax);
  Double_t calculateBaseline(TH1D*& histo, Bool_t totalRange = kTRUE, Bool_t reflected = kFALSE);
  Double_t calculateBaselineError(TH1D*& histo, Bool_t totalRange = kTRUE, Bool_t reflected = kFALSE);
  Bool_t readInputSeAndMe();
  Bool_t readInputInvMass();
  Bool_t readInputFdSubtr();
  Bool_t readInputSecondaryPartContamination();
  Bool_t extractCorrelations(Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax, TString codeName);
  TH1D* getCorrectedCorrHisto() { return fCorrectedCorrHisto; }
  TH1D* getCorrectedCorrHistoBaselineSubtr() { return fCorrectedCorrHistoBaselineSubtr; }
  TH1D* getCorrectedCorrHistoReflected() { return fCorrectedCorrHistoReflected; }
  TH1D* getCorrectedCorrHistoReflectedBaselineSubtr() { return fCorrectedCorrHistoReflectedBaselineSubtr; }

  /// Histogram style
  void setTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, Style_t markerStyle = kFullCircle, Color_t markerColor = kRed + 1, Double_t markerSize = 1.4, Color_t lineColor = kRed + 1, Int_t lineWidth = 3, Float_t hTitleXaxisOffset = 1.0, Float_t hTitleYaxisOffset = 1.0, Float_t hTitleXaxisSize = 0.060, Float_t hTitleYaxisSize = 0.060, Float_t hLabelXaxisSize = 0.060, Float_t hLabelYaxisSize = 0.060, Bool_t centerXaxisTitle = false, Bool_t centerYaxisTitle = false);
  void setTH2HistoStyle(TH2D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, TString hZaxisTitle, Float_t hTitleXaxisOffset = 1.8, Float_t hTitleYaxisOffset = 1.8, Float_t hTitleZaxisOffset = 1.2, Float_t hTitleXaxisSize = 0.060, Float_t hTitleYaxisSize = 0.060, Float_t hTitleZaxisSize = 0.060, Float_t hLabelXaxisSize = 0.060, Float_t hLabelYaxisSize = 0.060, Float_t hLabelZaxisSize = 0.060, Bool_t centerXaxisTitle = true, Bool_t centerYaxisTitle = true);

 private:
  TFile* fFileMass;         // File containing the mass histograms
  TFile* fFileSE;           // File containing the Same Event (SE) output
  TFile* fFileME;           // File containing the Mixed Event (ME) output
  TFile* fFileFDTemplate;   // File containing FD angular correlation templates
  TFile* fFileFDPromptFrac; // File containing prompt fraction (used fo FD subtraction)
  TFile* fFileSecPart;      // File containing secondary particle contaminaion teplates
  TFile* fFilePromptMc;     // File containing prompt ratio taken from MC Closure test study (to use for B to D bias correction)
  TFile* fFileNonPromptMc;  // File containing non-prompt ratio taken from MC Closure test study (to use for B to D bias correction)

  TDirectoryFile* fDirMass;    // TDirectory for mass histos
  TDirectoryFile* fDirSE;      // TDirectory for SE info
  TDirectoryFile* fDirME;      // TDirectory for ME info
  TDirectoryFile* fDirSecPart; // TDirectory for seondary particle correction

  TH1D* fCorrectedCorrHisto;                       // Corrected correlation histogram
  TH1D* fCorrectedCorrHistoBaselineSubtr;          // Corrected correlation histogram with baseline subtracion
  TH1D* fCorrectedCorrHistoReflected;              // Corrected correlation histogram relected in azimuth
  TH1D* fCorrectedCorrHistoReflectedBaselineSubtr; // Corrected correlation histogram reflected in azimuth with baseline subtraction

  DmesonSpecie fDmesonSpecies;           // D meson specie
  TString fDmesonLabel;                  // D meson label
  TString fFileNameMass;                 // File name containing inv. mass histograms
  TString fFileNameSE;                   // File name contaning Same Event (SE) output
  TString fFileNameME;                   // File name contaning Mixed Event (ME) output
  TString fFileSecPartName;              // File name contaning secondary particle correction output
  TString fFileFDTemplateName;           // File name contaning FD angular correlation templates
  TString fFileFDPromptFracName;         // File name contaning prompt fraction (used for FD subtraction)
  TString fFilePromptMcRecName;          // File name contaning prompt angular correlation (used for B to d bias correction)
  TString fFileNonPromptMcRecName;       // File name contaning non-prompt angular correlation (used for B to d bias correction)
  TString fDirNameSE;                    // Directory in the file containing SE output
  TString fDirNameME;                    // Directory in the file containing ME output
  TString fDirSecPartName;               // Directory in the file containing secondary particle correction output
  TString fMassHistoNameSgn;             // Inv. mass histo name signal yield
  TString fMassHistoNameBkg;             // Inv. mass histo name background yield
  TString fMassHistoNameSBs;             // Inv. mass histo name sideband yield
  TString fSECorrelSignalRegionName;     // THnSparse name containing SE output for signal region
  TString fSECorrelSidebandsName;        // THnSparse name containing SE output for sideband region
  TString fSECorrelSidebandLeftName;     // THnSparse name containing SE output for sideband left region
  TString fSECorrelSidebandRightName;    // THnSparse name containing SE output for sideband right region
  TString fMECorrelSignalRegionName;     // THnSparse name containing ME output for signal region
  TString fMECorrelSidebandsName;        // THnSparse name containing ME output for sideband regions
  TString fMECorrelSidebandLeftName;     // THnSparse name containing ME output for sideband left region
  TString fMECorrelSidebandRightName;    // THnSparse name containing ME output for sideband right region
  TString fHistoFDTemplatePromptName;    // Prompt angular correlation histogram name
  TString fHistoFDTemplateNonPromptName; // FD angular correlation histogram name
  TString fHistoFDPromptFracName;        // Prompt fraction histogram name
  TString fHistoPrimaryPartName;         // Primary particle histogram (to be used for secondary particle contamination correction)
  TString fHistoAllPartName;             // All particle histogram (to be used for secondary particle contamination correction)

  Int_t fNpools;            // Number of pools used for the ME correction
  Int_t fRebinAxisDeltaEta; // Rebin deltaEta axis value
  Int_t fRebinAxisDeltaPhi; // Rebin deltaPhi axis value
  Int_t fDebug;             // Debug level
  Int_t fBinPtCand;         // Pt bin of the candidate
  Int_t fBinPtHad;          // Pt bin of the hadron

  Double_t fDeltaEtaMin;    // DeltaEta min value
  Double_t fDeltaEtaMax;    // DeltaEta max value
  Double_t fBkgScaleFactor; // Bkg/SB factor to scale correlation plots obtained in the sideband region
  Double_t fSgnYieldNorm;   // Signal yield (used for normalize correlation plots after bkg subtraction)
  Double_t fBkgYield;       // Bkg yield under signal peak region
  Double_t fSBYield;        // Sideband yield

  Bool_t fCorrectPoolsSeparately; // Possibility to do the ME correction pool-by-pool (kTRUE) or merging all pools (kFALSE)
  Bool_t fSubtractSoftPiME;       // Soft pion subtraction (for D0 case)
  Bool_t fRebinAngCorr;           // Rebin angular correlaion distributons (SE and ME)
  Bool_t fRebinFDCorr;            // Rebin angular correlaion distributon templates used for FD correction (theory driven)
  Bool_t fRebinSecPart;           // Rebin angular correlaion distributon templates used for secodary particle contamination correction
  Bool_t fSidebandDivided;        // To be set to TRUE if two sideband corrlaion histograms are passed inteh config file
  Bool_t fUseSidebLeft;           // To be set to TRUE if only sideband left is used for the bkg correction
  Bool_t fUseSidebRight;          // To be set to TRUE if only sideband right is used for the bkg correction
  Bool_t fFDsubtraction;          // Enable feed-down (FD) correction
  Bool_t fSecPartContamination;   // Enable seconday particle contamination correction
  Bool_t fCorrBiasBtoD;           // Enable bias B to D correction (perfrmed with angular correlaion templates taken from MC simulaions)
};

#endif // PWGHF_HFC_MACROS_DHCORRELATIONEXTRACTION_H_
