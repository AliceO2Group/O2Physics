// \author Zhen Zhang <zhenz@cern.ch>
// \author Mingyu Zhang <mingyu.zang@cern.ch>
// \author Xinye Peng  <xinye.peng@cern.ch>
// \author Biao Zhang <biao.zhang@cern.ch>

// RooFit includes
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooDataHist.h"

#include "HFFitter.h"
using namespace RooFit;
using namespace std;


ClassImp(HFFitter);

HFFitter::HFFitter() :
  TNamed(),
  fHistoInvMass(0x0),
  fMinMass(0),
  fMaxMass(5),
  fTypeOfFit4Bkg(kExpo),
  fMassParticle(1.864),
  fTypeOfFit4Sgn(kGaus),
  fTypeOfFit4Rfl(1),
  fMass(1.865),
  fSecMass(1.969),
  fMassErr(0.),
  fSigmaSgn(0.012),
  fSecSigma(0.006),
  fNSigma4SideBands(4.),
  fNSigma4Sgn(3.),
  fSigmaSgnErr(0.),
  fSigmaSgn2Gaus(0.012),
  fFixedMean(kFALSE),
  fBoundMean(kFALSE),
  fMeanSignal(0x0),
  fSigmaSignal(0x0),
  fMassLowerLim(0),
  fMassUpperLim(0),
  fFixedSigma(kFALSE),
  fBoundSigma(kFALSE),
  fSigmaVar(0.012),
  fParSig(0.1),
  fRflOverSig(0),
  fReflections(kFALSE),
  fRawYield(0),
  fRawYieldErr(0),
  fBkgYield(0),
  fBkgYieldErr(0),
  fSignificance(0),
  fSignificanceErr(0),
  fChi2(0),
  fSigFunc(0x0),
  fBkgFunc(0x0),
  fRflFunc(0x0),
  fIntegralHisto(0),
  fIntegralBkg(0),
  fIntegralSig(0),
  fNSig(0x0),
  fNBkg(0x0),
  fTotFunc(0x0),
  mframe(0x0),
  rframe(0x0),
  residualframe(0x0),
  pullframe(0x0),
  w1(0x0),
  fHistoTemplRfl(0x0)
{ 
      // default constructor
}


HFFitter::HFFitter(const TH1F *histoToFit, Double_t minvalue, Double_t maxvalue, Int_t fittypeb, Int_t fittypes) :
  TNamed(),
  fHistoInvMass(0x0),
  fMinMass(minvalue),
  fMaxMass(maxvalue),
  fTypeOfFit4Bkg(fittypeb),
  fMassParticle(1.864),
  fTypeOfFit4Sgn(fittypes),
  fTypeOfFit4Rfl(1),
  fMass(1.865),
  fSecMass(1.969),
  fMassErr(0.),
  fSigmaSgn(0.012),
  fSecSigma(0.006),
  fNSigma4SideBands(3.),
  fNSigma4Sgn(3.),
  fSigmaSgnErr(0.),
  fSigmaSgn2Gaus(0.012),
  fFixedMean(kFALSE),
  fBoundMean(kFALSE),
  fMeanSignal(0x0),
  fSigmaSignal(0x0),
  fMassLowerLim(0),
  fMassUpperLim(0),
  fFixedSigma(kFALSE),
  fBoundSigma(kFALSE),
  fSigmaVar(0.012),
  fParSig(0.1),
  fRflOverSig(0),
  fReflections(kFALSE),
  fRawYield(0),
  fRawYieldErr(0),
  fBkgYield(0),
  fBkgYieldErr(0),
  fSignificance(0),
  fSignificanceErr(0),
  fChi2(0),
  fSigFunc(0x0),
  fBkgFunc(0x0),
  fRflFunc(0x0),
  fIntegralHisto(0),
  fIntegralBkg(0),
  fIntegralSig(0),
  fNSig(0x0),
  fNBkg(0x0),
  fTotFunc(0x0),
  mframe(0x0),
  rframe(0x0),
  residualframe(0x0),
  pullframe(0x0),
  w1(0x0),
  fHistoTemplRfl(0x0)
{
      // standard constructor
  fHistoInvMass=(TH1F*)histoToFit->Clone("fHistoInvMass");
  fHistoInvMass->SetDirectory(0);
}


HFFitter::~HFFitter()
{
  
  ///destructor
 
  delete fHistoInvMass;
  delete fHistoTemplRfl;
  delete fMeanSignal;
  delete fSigmaSignal;
  delete fSigFunc;
  delete fBkgFunc;
  delete fRflFunc;
  delete fTotFunc;
  delete fNSig;
  delete fNBkg;
  delete mframe;
  delete rframe;
  delete residualframe;
  delete pullframe;
  delete w1;
}


void HFFitter::MassFitter(Bool_t draw){
    fIntegralHisto = fHistoInvMass->Integral();


	w1 = new RooWorkspace("w", kTRUE);

	fillWorkspace(*w1);

	RooRealVar *mass = w1->var("mass");
    RooDataHist data("data","data",*mass,Import(*fHistoInvMass));    // Binned dataset from histogram to fit


	if (fTypeOfFit4Sgn==3) {									// Check Second peak or not
      mass->setRange("SBL", fMinMass, fMass - fNSigma4SideBands * fSigmaSgn);  // left sideband
      mass->setRange("SBR", fMass + fNSigma4SideBands * fSigmaSgn, fSecMass - fNSigma4SideBands * fSecSigma);  // Right sideband
      mass->setRange("SEC", fSecMass + fNSigma4SideBands * fSecSigma, fMaxMass);  // Right sideband
      mass->setRange("signal", fSecMass - fNSigma4SideBands *fSecSigma, fSecMass + fNSigma4SideBands * fSecSigma);  // signal interval

	}
	else {
      mass->setRange("SBL", fMinMass, fMass - fNSigma4SideBands * fSigmaSgn);  // left sideband
      mass->setRange("SBR", fMass + fNSigma4SideBands * fSigmaSgn, fMaxMass);  // Right sideband
      mass->setRange("signal", fMass - fNSigma4SideBands *fSigmaSgn, fMass + fNSigma4SideBands * fSigmaSgn);  // signal interval

	}
    mass->setRange("FULL",fMinMass, fMaxMass);    // full range
	
    mframe = mass->frame(Title("Invariant mass"),Range("FULL"));    // define an frame to plot

    data.plotOn(mframe,Name("data_c"));      // plot the dataset on the frame


	// --------------------- background fit functon type ---------------------------------- //
	
	RooAbsPdf *fbkgFunc;		// background fit PDF(probability density function)

	switch (fTypeOfFit4Bkg) {
	case 0:         // exponential fuction
	  { 
		fbkgFunc = w1->pdf("fbkgFunc0");
	  }
	  break;
	case 1:         // linear fuction
	  {
		fbkgFunc = w1->pdf("fbkgFunc1");
	  }
	  break;
	case 2:
	  {
		fbkgFunc = w1->pdf("fbkgFunc2");
	  }
	  break;
	case 3:
	  {
		fbkgFunc = w1->pdf("fbkgFunc3");
	  }
	  break;
	case 4:
  	  {
		fbkgFunc = w1->pdf("fbkgFunc4");
 	  }
	  break;
	case 5:
	  {
		fbkgFunc = w1->pdf("fbkgFunc5");
	  }
	  break;
	case 6:
	  {
		fbkgFunc = w1->pdf("fbkgFunc6");
	  }
	  break;
	default:
	  break;
	}
    
	// ----------------------------------------------------------------------------------- //

	// ----------------------signal fit function type ------------------------------------ //

	RooAbsPdf *fsigFunc;		// signal fit PDF

    switch (fTypeOfFit4Sgn) {
    case 0:
	  {
		fsigFunc = w1->pdf("fsigFunc0");
		fSigmaSignal= w1->var("sigma_sig0");
		fMeanSignal = w1->var("mean_sig");
	  }
      break;
    case 1:
	  {
		fsigFunc = w1->pdf("fsigFunc1");
		fSigmaSignal= w1->var("sigma2_sig1");
		fMeanSignal = w1->var("mean_sig");
	  }
      break;
    case 2:
	  {
		fsigFunc = w1->pdf("fsigFunc2");
		fSigmaSignal= w1->var("sigma2_sig2");
		fMeanSignal = w1->var("mean_sig");
	  }
      break;
    case 3:
	  {
		fsigFunc = w1->pdf("fsigFunc3");
		fMeanSignal = w1->var("mean_sig2");
		fSigmaSignal= w1->var("sigma2_sig3");
	  }
      break;
    default:
      break;
    }
    
	// ---------------------------------------------------------------------------------- //

	// -------------------------- reflection fit function type -------------------------- //

	RooAbsPdf *frflFunc;		// reflection fit PDF

    if (fHistoTemplRfl) {

      switch (fTypeOfFit4Rfl) {
      case 0:
	    {
		  frflFunc = w1->pdf("frflFunc0");
	    }
        break;
      case 1:
	    {
		  frflFunc = w1->pdf("frflFunc1");
	    }
        break;
      case 2:
	    {
		  frflFunc = w1->pdf("frflFunc2");
	    }
	    break;
      case 3:
	    {
		  frflFunc = w1->pdf("frflFunc3");
	    }
	    break;
      default:
        break;
        }
	}
	
	// ---------------------------------------------------------------------------------- //

    printf("\n---First fit with only backround on the sidebands - Exclusion region = %d sigma ---\n",fNSigma4SideBands);

	// define number of background and background fit function 
    fNBkg = new RooRealVar("fNBkg","number of background", 0.5 * fIntegralHisto, 0, fIntegralHisto);     // num of background
    fBkgFunc = new RooAddPdf("fBkgFunc", "background fit function", RooArgList(*fbkgFunc), RooArgList(*fNBkg));
	
	if (fTypeOfFit4Sgn==3) {	
      fBkgFunc->fitTo(data, Range("SBL,SBR,SEC"),Save());
      fBkgFunc->plotOn(mframe,Range("FULL"),LineColor(kRed-4),Name("Bkg_c"));
	}
	else {
      fBkgFunc->fitTo(data, Range("SBL,SBR"),Save());
      fBkgFunc->plotOn(mframe,Range("FULL"),LineColor(kRed-4),Name("Bkg_c"));
	}


	// define reflection template
    if (fHistoTemplRfl) {
	  RooRealVar *mass_r = w1->var("mass_r");
      RooDataHist reflhist("reflhist", "refl for Fit", *mass_r,Import(*fHistoTemplRfl));
      rframe = mass_r->frame(Range("FULL"));

	  RooRealVar nrfl("nrfl","number of reflection template from simulation", 0.5 * fHistoTemplRfl->Integral(), 0, fHistoTemplRfl->Integral());
	  RooAddPdf fRflFuncTemp("fRflFuncTemp", "template reflection fit function", RooArgList(*frflFunc), RooArgList(nrfl));
	  fRflFuncTemp.fitTo(reflhist);
	}

	RooAbsReal *bk = fBkgFunc->createIntegral(*mass,NormSet(*mass),Range("signal"));	//bkg integral
	fIntegralBkg = bk->getValV();

	Double_t estimSignal;
	CheckForSignal(estimSignal);
	Background(fBkgYield, fBkgYieldErr);

    if (fHistoTemplRfl) {
	RooRealVar nRfl("nRfl", "reflection scale to data", fRflOverSig * estimSignal);

	fRflFunc = new RooAddPdf("fRflFunc", "reflection fit function", RooArgList(*frflFunc), RooArgList(nRfl));

	fRflFunc->plotOn(rframe, Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(kGreen));
	}

	printf("\n---Fit with background + signal ---\n");
    fNSig = new RooRealVar("fNSig","number of signal", 0.5 * fIntegralHisto, 0, fIntegralHisto);    // num of signal
    fSigFunc = new RooAddPdf("fSigFunc", "signal fit function", RooArgList(*fsigFunc),RooArgList(*fNSig));
	fTotFunc = new RooAddPdf("fBkgPlusSigFunc", "background + signal fit function", RooArgList(*fbkgFunc,*fsigFunc), RooArgList(*fNBkg,*fNSig));
	fTotFunc->fitTo(data,Range("FULL"));		// fit the dataset with bkg + sig
    fSigFunc->plotOn(mframe,Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(kGreen));
    fTotFunc->plotOn(mframe,Range("FULL"),Name("Tot_c"));

	fChi2 = mframe->chiSquare("Tot_c","data_c");		// calculate chi2

	RooAbsReal *sg = fSigFunc->createIntegral(*mass,NormSet(*mass),Range("signal"));	//bkg integral
	fIntegralSig = sg->getValV();

	Signal(fRawYield,fRawYieldErr);
	Significance(fSignificance, fSignificanceErr);

	RooHist* hresid = mframe->residHist("data_c","Bkg_c");
	RooHist* hpull = mframe->pullHist("data_c","Tot_c");
	residualframe =  mass->frame(Title("Residual Distribution"));
	pullframe = mass->frame(Title("Pull Distribution"));
	
	residualframe->addPlotable(hresid,"p");
	fSigFunc->plotOn(residualframe,Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(kBlue));
    pullframe->addPlotable(hpull, "P");


}

// ---------------------------- Fit function ----------------------------------------------------- //
void HFFitter::fillWorkspace(RooWorkspace &w)
{
	//Declare observable variable
    RooRealVar mass("mass", "mass", fMinMass, fMaxMass, "GeV/c");


	//	-------------------------background fit function--------------------------------------- //

	// exp
	RooRealVar tau("tau", "tau", 1, -10., 10.);
    RooAbsPdf* fbkgFunc0 = new RooExponential("fbkgFunc0", "background fit function", mass, tau);
	w.import(*fbkgFunc0);

	// liner
    RooRealVar a1("a1", "a1", -10, 10);
    RooRealVar b1("b1", "b1", -10, 10);
    RooAbsPdf* fbkgFunc1 = new RooGenericPdf("fbkgFunc1", "a1 + b1 * mass",RooArgSet(mass, a1, b1));
	w.import(*fbkgFunc1);

	// poly2
    RooRealVar p02("p02","p02",0.5,-1.,1.);        // poly2
    RooRealVar p12("p12","p12",0.2,-1.,1.);
    RooRealVar p22("p22","p22",0.2,-1.,1.);
    RooAbsPdf* fbkgFunc2 = new RooPolynomial("fbkgFunc2","background fit function", mass, RooArgSet(p02,p12,p22));      // background function
	w.import(*fbkgFunc2);

	// no background
    RooRealVar a3("a3", "a3",0);           // no background
    RooAbsPdf* fbkgFunc3 = new RooGenericPdf("fbkgFunc3", "a3", RooArgSet(a3));
	w.import(*fbkgFunc3);

	// pow
    RooRealVar a4("a4","a4", 0., 100000.);       // pow function
    RooRealVar b4("b4", "b4", 1, 0, 10);
    RooAbsPdf* fbkgFunc4 = new RooGenericPdf("fbkgFunc4", "(mass-a4)^b4",RooArgSet(mass,a4,b4));
	w.import(*fbkgFunc4);

	// pow*exp
    RooRealVar a5("a5","a5", 0., 100000.);       // pow*exp
    RooRealVar b5("b5", "b5", 1, 0, 10);
    RooAbsPdf* fbkgFunc5 = new RooGenericPdf("fbkgFunc5", "sqrt(mass-a5)*exp(-b5*(mass-a5))",RooArgSet(mass,a5,b5));
	w.import(*fbkgFunc5);

	// poly3
    RooRealVar p06("p06","p06",0.5,-1.,1.);        // poly3
    RooRealVar p16("p16","p16",0.2,-1.,1.);
    RooRealVar p26("p26","p26",0.2,-1.,1.);
    RooRealVar p36("p36","p36",0.2,-1.,1.);
    RooAbsPdf* fbkgFunc6 = new RooPolynomial("fbkgFunc6","background PDF", mass, RooArgSet(p06,p16,p26,p36));
	w.import(*fbkgFunc6);

	// --------------------------------------------------------------------------------------- //


	//	----------------------------signal fit function--------------------------------------- //

    RooRealVar mean_sig("mean_sig", "mean for signal fit", fMass, fMass - 0.03, fMass + 0.03);

	// Guassian
    RooRealVar sigma_sig0("sigma_sig0", "sigma for signal", fSigmaSgn, fSigmaSgn - 0.01, fSigmaSgn + 0.01);
    RooAbsPdf* fsigFunc0 = new RooGaussian("fsigFunc0", "signal PDF", mass, mean_sig, sigma_sig0);
	w.import(*fsigFunc0);

	// Two Gaussian
    RooRealVar sigma1_sig1("sigma1_sig1", "sigma1_sig", fSigmaSgn, fSigmaSgn - 0.01, fSigmaSgn + 0.01);
    RooRealVar sigma2_sig1("sigma2_sig1", "sigma2_sig", fSigmaSgn, fSigmaSgn - 0.01, fSigmaSgn + 0.01);
    RooGaussian g11("g11", "g1", mass, mean_sig, sigma1_sig1);
    RooGaussian g21("g21", "g2", mass, mean_sig, sigma2_sig1);
    RooRealVar frac_sig1("frac_sig1", "frac of two gauss", 0.5, 0, 1.);
    RooAbsPdf* fsigFunc1 = new RooAddPdf("fsigFunc1", "signal PDF", RooArgList(g11, g21), frac_sig1);
	w.import(*fsigFunc1);

	// Gaussian ratio
    RooRealVar sigma1_sig2("sigma1_sig2", "sigma1_sig", fSigmaSgn, fSigmaSgn - 0.1, fSigmaSgn + 0.1);
    RooRealVar ratio2("ratio2", "ratio of sigma12", 0, 10);
    RooRealVar sigma2_sig2("sigma2_sig2", "sigma2_sig", sigma1_sig2.getVal()*ratio2.getVal());
    RooGaussian g12("g12", "g12", mass, mean_sig, sigma1_sig2);
    RooGaussian g22("g22", "g22", mass, mean_sig, sigma2_sig2);
    RooRealVar frac_sig2("frac_sig2", "frac of two gauss", 0.5, 0, 1.);
    RooAbsPdf* fsigFunc2 = new RooAddPdf("fsigFunc2", "signal PDF", RooArgList(g12, g22), frac_sig2);
	w.import(*fsigFunc2);

	// two peak
    RooRealVar mean_sig1("mean_sig1", "mean for signal1 fit", 1.85, fMinMass, fMaxMass);
    RooRealVar mean_sig2("mean_sig2", "mean for signal2 fit", 1.95, fMinMass, fMaxMass);
    RooRealVar sigma1_sig3("sigma1_sig3", "sigma1_sig", fSigmaSgn, fSigmaSgn - 0.01, fSigmaSgn + 0.01);
    RooRealVar sigma2_sig3("sigma2_sig3", "sigma2_sig", fSigmaSgn, fSigmaSgn - 0.01, fSigmaSgn + 0.01);
    RooGaussian g13("g13", "g1", mass, mean_sig1, sigma1_sig3);
    RooGaussian g23("g23", "g2", mass, mean_sig2, sigma2_sig3);
    RooRealVar frac_sig3("frac_sig3", "frac of two gauss", 0.5, 0, 1.);
    RooAbsPdf* fsigFunc3 = new RooAddPdf("fsigFunc3", "signal PDF", RooArgList(g13, g23), frac_sig3);
	w.import(*fsigFunc3);

	// --------------------------------------------------------------------------------------- //




	//	----------------------------reflection fit function----------------------------------- //

   	RooRealVar mass_r("mass_r","mass_r", fMinMass, fMaxMass,"GeV/c^{2}");
	
	RooRealVar mean_rfl0("mean_rfl0", "mean for reflection", fMass, fMass - 0.03, fMass + 0.03);
	RooRealVar sigma_rfl0("sigma_rfl0", "sigma for reflection", fSigmaSgn, fSigmaSgn - 0.1, fSigmaSgn + 0.1);
	RooAbsPdf* frflFunc0 = new RooGaussian("frflFunc0", "reflection PDF", mass_r, mean_rfl0, sigma_rfl0);
	w.import(*frflFunc0);

	RooRealVar mean_rfl11("mean_rfl11", "mean", fMass, fMinMass, fMaxMass );
	RooRealVar mean_rfl21("mean_rfl21", "mean", fMass, fMinMass, fMaxMass);
	RooRealVar sigma1_rfl1("sigma1_rfl1", "sigma1_rfl", fSigmaSgn, fSigmaSgn - 0.1, fSigmaSgn + 0.1);
	RooRealVar sigma2_rfl1("sigma2_rfl1", "sigma2_rfl", fSigmaSgn, fSigmaSgn - 0.1, fSigmaSgn + 0.1);
	RooGaussian g1_rfl1("g1_rfl1", "g1", mass_r, mean_rfl11, sigma1_rfl1);
	RooGaussian g2_rfl1("g2_rfl1", "g2", mass_r, mean_rfl21, sigma2_rfl1);
	RooRealVar frac_rfl1("frac_rfl", "frac of two gauss", 0.5, 0, 1.);
	RooAbsPdf* frflFunc1 = new RooAddPdf("frflFunc1", "reflection PDF", RooArgList(g1_rfl1, g2_rfl1), frac_rfl1);
	w.import(*frflFunc1);

	RooRealVar p0_rfl2("p0_rfl2","p0",0.5,-1.,1.);        // poly3
	RooRealVar p1_rfl2("p1_rfl2","p1",0.2,-1.,1.);
	RooRealVar p2_rfl2("p2_rfl2","p2",0.2,-1.,1.);
	RooRealVar p3_rfl2("p3_rfl2","p3",0.2,-1.,1.);
	RooAbsPdf* frflFunc2 = new RooPolynomial("frflFunc2","reflection PDF", mass_r, RooArgSet(p0_rfl2,p1_rfl2,p2_rfl2,p3_rfl2));
	w.import(*frflFunc2);

	RooRealVar p0_rfl3("p0_rfl3","p0",0.5,-1.,1.);        // poly6
	RooRealVar p1_rfl3("p1_rfl3","p1",0.2,-1.,1.);
	RooRealVar p2_rfl3("p2_rfl3","p2",0.2,-1.,1.);
	RooRealVar p3_rfl3("p3_rfl3","p3",0.2,-1.,1.);
	RooRealVar p4_rfl3("p4_rfl3","p3",0.2,-1.,1.);
	RooRealVar p5_rfl3("p5_rfl3","p3",0.2,-1.,1.);
	RooRealVar p6_rfl3("p6_rfl3","p3",0.2,-1.,1.);
	RooAbsPdf* frflFunc3 = new RooPolynomial("frflFunc3","reflection PDF", mass_r, RooArgSet(p0_rfl3,p1_rfl3,p2_rfl3,p3_rfl3,p4_rfl3,p5_rfl3,p6_rfl3));
	w.import(*frflFunc3);

	// --------------------------------------------------------------------------------------- //

}

// ---------------------------------------------------------------------------------------------- //


// --------------------------------- draw the fit output -----------------------------------------//

void HFFitter::DrawHere(TVirtualPad* c,Int_t writeFitInfo){
  /// Core method to draw the fit output
  ///

    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);
    c->cd();

    if(writeFitInfo > 0){
	  TPaveText *pinfos=new TPaveText(0.12,0.65,0.47,0.89,"NDC");
	  TPaveText *pinfom=new TPaveText(0.6,0.7,1.,.87,"NDC");
	  pinfos->SetBorderSize(0);
	  pinfos->SetFillStyle(0);
	  pinfom->SetBorderSize(0);
	  pinfom->SetFillStyle(0);
	  pinfom->SetTextColor(kBlue);
	
	  pinfos->AddText(Form("S = %.0f #pm %.0f ",fRawYield,fRawYieldErr));
	  pinfos->AddText(Form("B (%d#sigma) = %.0f #pm %.0f",fNSigma4SideBands,fBkgYield,fBkgYieldErr));
	  pinfos->AddText(Form("S/B (%d#sigma) = %.4g ",fNSigma4SideBands,fRawYield/fBkgYield));
	  if(fRflFunc)  pinfos->AddText(Form("Refl/Sig =  %.3f #pm %.3f ",fRflOverSig,0.0));
	  pinfos->AddText(Form("Signif (%d#sigma) = %.1f #pm %.1f ",fNSigma4SideBands,fSignificance,fSignificanceErr));
	  pinfos->AddText(Form("Chi2/NDF (%d#sigma) = %.3f ",fNSigma4SideBands,mframe->chiSquare()));
	
	  pinfom->AddText(Form("mean = %.3f #pm %.3f",fMeanSignal->getVal(), fMeanSignal->getError()));
	  pinfom->AddText(Form("sigma = %.3f #pm %.3f",fSigmaSignal->getVal(), fSigmaSignal->getError()));
	
	  mframe->addObject(pinfos);
	  mframe->addObject(pinfom);

	  mframe->Draw();
    if (fHistoTemplRfl) {
	  rframe->Draw("same");
	}
  }
  c->Update();
  return;
}

// ---------------------------------------------------------------------------------------------- //


// --------------------------------- draw the fit output -----------------------------------------//

void HFFitter::DrawPull(TVirtualPad* c){
	
	//c->Divide(2);
	c->cd();
	residualframe->GetYaxis()->SetTitle("");
	residualframe->Draw();
	//c->cd(2);
	//pullframe->Draw();
}



// ---------------------------------------------------------------------------------------------- //


// -------------------------------------- calculate signal -------------- //
void HFFitter::Signal(Double_t &signal, Double_t &errsignal) const {
	

	signal 		   = fNSig->getVal() * fIntegralSig;
	errsignal      = fNSig->getError() * fIntegralSig;

	return;
}

// ---------------------------------------------------------------------------------------------- //


// -------------------------------------- calculate background -------------- //
void HFFitter::Background(Double_t &bkg, Double_t &errbkg) const {
	
	bkg 		   = fNBkg->getVal() * fIntegralBkg;
	errbkg      = fNBkg->getError() * fIntegralBkg;
	
	return;
}

// ---------------------------------------------------------------------------------------------- //

// -------------------------------------- calculate background -------------- //
void HFFitter::Significance(Double_t &significance, Double_t &errsignificance) const {
	
	Double_t signal, errsignal;
	Signal(signal,errsignal);

	Double_t bkg, errbkg;
	Background(bkg, errbkg);
	
	Double_t m_sigerrSq = errsignal * errsignal;
	Double_t m_bkgerrSq = errbkg	* errbkg;
	Double_t m_tot		= signal	+ bkg;

	significance = signal/sqrt(signal + bkg);
	errsignificance = significance * sqrt((m_sigerrSq + m_bkgerrSq)/(fNSigma4SideBands * m_tot * m_tot) + (bkg/m_tot) * (m_sigerrSq/signal/signal)); 

	return;
}

// ------------------------------------- Check Signnal ------------------------------------------- //

void HFFitter::CheckForSignal(Double_t &estimSignal) {
	
	Double_t minForSig = fMass - fNSigma4SideBands * fSigmaSgn;
	Double_t maxForSig = fMass + fNSigma4SideBands * fSigmaSgn;
	
	Int_t binForMinSig = fHistoInvMass->FindBin(minForSig);
	Int_t binForMaxSig = fHistoInvMass->FindBin(maxForSig);

	Double_t sum = 0;

	for (Int_t i = binForMinSig; i <= binForMaxSig; i++) {
	  sum += fHistoInvMass->GetBinContent(i);
	}
	Double_t bkg,errbkg;
	Background(bkg,errbkg);

	estimSignal = sum - bkg;

}

// ---------------------------------------------------------------------------------------------- //


