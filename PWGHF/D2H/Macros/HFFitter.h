#ifndef HFFITTER_H
#define HFFITTER_H


#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <RooWorkspace.h>


#include <TNamed.h>

using namespace RooFit;
using namespace std;

class TF1;
class TH1F;

class HFFitter : public TNamed {
public:
  enum ETypeOfBkg{ kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5};
  enum ETypeOfSgn{ kGaus=0, k2Gaus=1, k2GausSigmaRatioPar=2, kGausSec=3 };
  enum ETypeOfRfl{ kGaus1=0, kGaus2=1, kPol3=2, kPol6=3};
  HFFitter();
  HFFitter(const TH1F* histoToFit, Double_t minvalue, Double_t maxvalue, Int_t fittypeb=kExpo, Int_t fittypes=kGaus);
  ~HFFitter();
  void     SetHistogramFit(const TH1F* histoToFit){
    if (fHistoInvMass) delete fHistoInvMass;
    fHistoInvMass=(TH1F*)histoToFit->Clone("fHistoInvMass");
    fHistoInvMass->SetDirectory(0);
  }
  void     SetRangeFit(Double_t minvalue, Double_t maxvalue){
    fMinMass=minvalue; 
    fMaxMass=maxvalue;
  }
  void     SetFitFunctions(Int_t fittypeb, Int_t fittypes){
    fTypeOfFit4Bkg=fittypeb; 
    fTypeOfFit4Sgn=fittypes;
  }
  void		MassFitter(Bool_t draw=kTRUE);
  void		SetNSigma4SideBands(Double_t ns=4.) {
  	fNSigma4SideBands=ns;
  }
  void 		SetInitialReflOverS(Double_t rovers) {fRflOverSig = rovers;}

  void		SetFixReflOverS(Double_t rovers) {
	SetInitialReflOverS(rovers);
	fFixRflOverSig = kTRUE;
  }
  
  void 		SetTemplateReflections(const TH1 *h, Int_t fittyper=kGaus2) {
	if(!h) {
	  fReflections = kFALSE;
	}
	fHistoTemplRfl = (TH1F*)h->Clone("fHistoTemplRfl");

  }

  Double_t 	GetChiSquareOverNDF() const{return fChi2;}
  
  
  Double_t 	GetRawYield()const {return (int)fRawYield;}
  Double_t 	GetRawYieldError()const {return (int)fRawYield;}
  Double_t 	GetMean() const {return fMeanSignal->getVal();}
  Double_t 	GetMeanUncertainty() const {return fMeanSignal->getError();}
  Double_t 	GetSigma()const {return fSigmaSignal->getVal();}
  Double_t 	GetSigmaUncertainty()const { return fSigmaSignal->getError();}
  Double_t 	GetReflOverSig()const{
    if(fRflFunc) return fRflOverSig;
    else return 0;
  }
  void		Signal(Double_t &signal, Double_t &errsignal) const;
  void		Background(Double_t &bkg, Double_t &errbkg) const;
  void 		Significance(Double_t &significance, Double_t &errsignificance) const;
  void 		CheckForSignal(Double_t &estimSignal);

  void		DrawHere(TVirtualPad* c,Int_t writeFitInfo = 2);
  void		DrawPull(TVirtualPad* c);

private:
  HFFitter(const HFFitter &source);
  HFFitter& operator=(const HFFitter& source);
  //void		fillWorkspace(RooWorkspace &w);
  void		fillWorkspace(RooWorkspace &w);

  TH1F*			fHistoInvMass;			/// histogram to fit
  Double_t		fMinMass;				/// lower mass limit
  Double_t		fMaxMass;				/// upper mass limit
  Int_t			fTypeOfFit4Bkg;			/// background fit function
  Double_t		fMassParticle;			///	pdg value of particle mass
  Int_t			fTypeOfFit4Sgn;			/// signal fit function
  Int_t			fTypeOfFit4Rfl;			/// reflection fit function
  Double_t		fMass;					/// signal gaussian mean value
  Double_t		fSecMass;				/// Second peak mean value
  Double_t		fMassErr;				/// uncertainty on signal gaussian mean value
  Double_t		fSigmaSgn;				/// signal gaussian sigma
  Double_t		fSecSigma;				/// Second peak gaussian sigma
  Int_t			fNSigma4SideBands;		///	number of sigmas to veto the signal peak
  Int_t			fNSigma4Sgn;			/// number of sigmas to veto the signal peak
  Double_t		fSigmaSgnErr;			/// uncertainty on signal gaussian sigma
  Double_t		fSigmaSgn2Gaus;			/// signal 2gaussian sigma
  Double_t		fFixedMean;				/// switch for fix mean of gaussian
  Bool_t		fBoundMean;				/// 
  RooRealVar*	fMeanSignal;
  RooRealVar*	fSigmaSignal;
  Double_t		fMassLowerLim;
  Double_t		fMassUpperLim;
  Bool_t		fFixedSigma;
  Bool_t		fBoundSigma;
  Double_t		fSigmaVar;
  Double_t		fParSig;
  Double_t 		fRflOverSig;
  Bool_t 		fReflections;
  Double_t		fRawYield;
  Double_t		fRawYieldErr;
  Double_t		fBkgYield;
  Double_t		fBkgYieldErr;
  Double_t		fSignificance;
  Double_t		fSignificanceErr;
  Double_t		fChi2;
  Double_t 		fFixRflOverSig;
  RooAbsPdf* 	fSigFunc;
  RooAbsPdf* 	fBkgFunc;
  RooAbsPdf* 	fRflFunc;
  Double_t		fIntegralHisto;
  Double_t		fIntegralBkg;
  Double_t		fIntegralSig;
  RooRealVar*	fNSig;
  RooRealVar*	fNBkg;
  RooAbsPdf*	fTotFunc;
  RooPlot*		mframe;
  RooPlot*		rframe;
  RooPlot*		residualframe;
  RooPlot*		pullframe;
  RooRealVar*	mass;
  RooWorkspace* w1;
  TH1F*			fHistoTemplRfl;

  ClassDef(HFFitter,1);
};

#endif
