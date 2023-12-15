
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include <TGraphErrors.h>

using namespace o2;
using namespace o2::framework;

struct trackTuner {

  TGraphErrors* grDcaXYResVsPtPionCurrent;
  TGraphErrors* grDcaXYResVsPtPionUpgrded;

  TGraphErrors* grDcaZResVsPtPionCurrent;
  TGraphErrors* grDcaZResVsPtPionUpgrded;

  TGraphErrors* grDcaXYMeanVsPtPionCurrent;
  TGraphErrors* grDcaXYMeanVsPtPionUpgrded;

  TGraphErrors* grDcaZMeanVsPtPionCurrent;
  TGraphErrors* grDcaZMeanVsPtPionUpgrded;

  TGraphErrors* grOneOverPtPionCurrent;
  TGraphErrors* grOneOverPtPionUpgrded;

  TGraphErrors* grDcaXYPullVsPtPionCurrent;
  TGraphErrors* grDcaXYPullVsPtPionUpgrded;

  TGraphErrors* grDcaZPullVsPtPionCurrent;
  TGraphErrors* grDcaZPullVsPtPionUpgrded;

  void getDcaGraphs(std::string pathCurrFileDcaXY, std::string pathCurrFileDcaZ, std::string pathUpgrFileDcaXY, std::string pathUpgrFileDcaZ)
  {

    const TString fnameDCAxyFileCurr = pathCurrFileDcaXY.data();
    const TString fnameDCAxyFileUpgr = pathUpgrFileDcaXY.data();

    const TString fnameDCAzFileCurr = pathCurrFileDcaZ.data();
    const TString fnameDCAzFileUpgr = pathUpgrFileDcaZ.data();

    if ((fnameDCAxyFileCurr == "") | (fnameDCAzFileCurr == "") | (fnameDCAxyFileUpgr == "") | (fnameDCAzFileUpgr == "")) {
      LOG(fatal) << "No Correction DCA files provided!";
      return;
    }

    std::unique_ptr<TFile> fileCurrDcaXY(TFile::Open(fnameDCAxyFileCurr.Data(), "READ"));
    std::unique_ptr<TFile> fileUpgrDcaXY(TFile::Open(fnameDCAxyFileUpgr.Data(), "READ"));

    std::unique_ptr<TFile> fileCurrDcaZ(TFile::Open(fnameDCAzFileCurr.Data(), "READ"));
    std::unique_ptr<TFile> fileUpgrDcaZ(TFile::Open(fnameDCAzFileUpgr.Data(), "READ"));

    TString grDcaResName = "tge_DCA_res_withoutPVrefit_all";
    TString grDcaMeanName = "tge_DCA_mean_withoutPVrefit_all";
    TString grDcaPullName = "tge_DCAPulls_res_withoutPVrefit_all";

    grDcaXYResVsPtPionCurrent = (TGraphErrors*)fileCurrDcaXY->Get(grDcaResName.Data());
    grDcaXYResVsPtPionUpgrded = (TGraphErrors*)fileUpgrDcaXY->Get(grDcaResName.Data());

    grDcaXYMeanVsPtPionCurrent = (TGraphErrors*)fileCurrDcaXY->Get(grDcaMeanName.Data());
    grDcaXYMeanVsPtPionUpgrded = (TGraphErrors*)fileUpgrDcaXY->Get(grDcaMeanName.Data());

    grDcaXYPullVsPtPionCurrent = (TGraphErrors*)fileCurrDcaXY->Get(grDcaPullName.Data());
    grDcaXYPullVsPtPionUpgrded = (TGraphErrors*)fileUpgrDcaXY->Get(grDcaPullName.Data());

    grDcaZResVsPtPionCurrent = (TGraphErrors*)fileCurrDcaZ->Get(grDcaResName.Data());
    grDcaZResVsPtPionUpgrded = (TGraphErrors*)fileUpgrDcaZ->Get(grDcaResName.Data());

    grDcaZMeanVsPtPionCurrent = (TGraphErrors*)fileCurrDcaZ->Get(grDcaMeanName.Data());
    grDcaZMeanVsPtPionUpgrded = (TGraphErrors*)fileUpgrDcaZ->Get(grDcaMeanName.Data());

    grDcaZPullVsPtPionCurrent = (TGraphErrors*)fileCurrDcaZ->Get(grDcaPullName.Data());
    grDcaZPullVsPtPionUpgrded = (TGraphErrors*)fileUpgrDcaZ->Get(grDcaPullName.Data());
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void tuneTrackParams(T1 mcparticle, T2& trackParCov, T3 matCorr, T4 dcaInfoCov, T5 isUpdateTrackCovMat, T6 isUpdateCurvature, T7 isUpdatePulls)
  {

    Double_t ptMC = TMath::Abs(mcparticle.pt());

    Double_t DcaXYResCurrent = 0.0; // sd0rpo=0.;
    Double_t DcaZResCurrent = 0.0;  // sd0zo =0.;

    Double_t DcaXYResUpgrded = 0.0; // sd0rpn=0.;
    Double_t DcaZResUpgrded = 0.0;  // sd0zn =0.;

    Double_t OneOverPtCurrent = 0.0; // spt1o =0.;
    Double_t OneOverPtUpgrded = 0.0; // spt1n =0.;

    Double_t DcaXYMeanCurrent = 0.0; // sd0mrpo=0.;
    Double_t DcaXYMeanUpgrded = 0.0; // sd0mrpn=0.;

    Double_t DcaXYPullCurrent = 1.0;
    Double_t DcaXYPullUpgrded = 1.0;

    Double_t DcaZPullCurrent = 1.0;
    Double_t DcaZPullUpgrded = 1.0;

    DcaXYResCurrent = EvalGraph(ptMC, grDcaXYResVsPtPionCurrent);
    DcaXYResUpgrded = EvalGraph(ptMC, grDcaXYResVsPtPionUpgrded);

    // DcaXYResCurrent = 1.0;
    // DcaXYResUpgrded = 1.5;

    DcaZResCurrent = EvalGraph(ptMC, grDcaZResVsPtPionCurrent);
    DcaZResUpgrded = EvalGraph(ptMC, grDcaZResVsPtPionUpgrded);

    // DcaZResCurrent = 1.0;
    // DcaZResUpgrded = 1.0;

    // OneOverPtCurrent = EvalGraph(ptMC, grOneOverPtPionCurrent );
    // OneOverPtUpgrded = EvalGraph(ptMC, grOneOverPtPionUpgrded );

    // OneOverPtCurrent = 1.0;
    // OneOverPtUpgrded = 2.0;

    DcaXYMeanCurrent = EvalGraph(ptMC, grDcaXYMeanVsPtPionCurrent);
    DcaXYMeanUpgrded = EvalGraph(ptMC, grDcaXYMeanVsPtPionUpgrded);

    // DcaXYMeanCurrent = 0.0;
    // DcaXYMeanUpgrded = 0.0;

    DcaXYPullCurrent = EvalGraph(ptMC, grDcaXYPullVsPtPionCurrent);
    DcaXYPullUpgrded = EvalGraph(ptMC, grDcaXYPullVsPtPionUpgrded);

    DcaZPullCurrent = EvalGraph(ptMC, grDcaZPullVsPtPionCurrent);
    DcaZPullUpgrded = EvalGraph(ptMC, grDcaZPullVsPtPionUpgrded);

    //  Unit conversion, is it required ??
    DcaXYResCurrent *= 1.e-4;
    DcaZResCurrent *= 1.e-4;

    DcaXYResUpgrded *= 1.e-4;
    DcaZResUpgrded *= 1.e-4;

    DcaXYMeanCurrent *= 1.e-4;
    DcaXYMeanUpgrded *= 1.e-4;

    // Apply the smearing
    // ---------------------------------------------
    // Double_t pt1o  =param  [4];
    Double_t trackParOneOverPtCurrent = trackParCov.getQ2Pt();
    int sign = trackParCov.getQ2Pt() / TMath::Abs(trackParCov.getQ2Pt());

    // Double_t pt1mc =parammc[4];
    Double_t trackParOneOverPtMC = sign / mcparticle.pt();

    std::cout << std::scientific << std::setprecision(5);
    o2::dataformats::VertexBase vtxMC;
    vtxMC.setPos({mcparticle.vx(), mcparticle.vy(), mcparticle.vz()});
    vtxMC.setCov(0, 0, 0, 0, 0, 0); // ??? or All ZEROs // == 1 cm2? wrt prop point

    // std::cout << " sign " << sign << std::endl;
    // std::cout << " trackParCov.getQ2Pt() " << trackParOneOverPtCurrent << " " << trackParCov.getQ2Pt() << std::endl;
    // std::cout << " sign/mcparticle.pt() " << trackParOneOverPtMC << std::endl;
    // std::cout << " (curvReco-curvMC)/curvMC " << (trackParOneOverPtCurrent - trackParOneOverPtMC)/trackParOneOverPtMC * 100.0 << "%" << std::endl;
    // std::cout << " trackParCov.getPtInv() " << trackParCov.getPtInv() <<  std::endl;
    // std::cout << " 1/trackParCov.getPtInv() " << 1./trackParCov.getPtInv() << " & mcparticle.pt() " << mcparticle.pt() << std::endl;

    std::cout << mcparticle.pt() << " " << 1 / trackParCov.getPtInv() << " " << (trackParOneOverPtCurrent - trackParOneOverPtMC) / trackParOneOverPtMC * 100.0 << std::endl;

    // std::cout << "Before Propagation to Production Point -> alpha: " << trackParCov.getAlpha() << ", DCAxy: " << trackParCov.getY() << ", DCAz: " << trackParCov.getZ() << std::endl;

    // propagate to DCA with respect to the Production point
    o2::base::Propagator::Instance()->propagateToDCABxByBz(vtxMC, trackParCov, 2.f, matCorr, dcaInfoCov);

    // std::cout << "After Propagation to Production Point -> alpha: " << trackParCov.getAlpha() << ", DCAxy: " << trackParCov.getY() << ", DCAz: " << trackParCov.getZ() << std::endl;

    //////////////////////////////  DCAs modifications Start  /////////////////////////////////

    std::cout << "track.y(): " << trackParCov.getY() << std::endl;

    // Double_t d0zo  =param  [1];
    Double_t trackParDcaZCurrent = trackParCov.getZ();

    // Double_t d0rpo =param  [0];
    Double_t trackParDcaXYCurrent = trackParCov.getY();

    float mcVxRotated = mcparticle.vx() * TMath::Cos(trackParCov.getAlpha()) + mcparticle.vy() * TMath::Sin(trackParCov.getAlpha()); // invert
    float mcVyRotated = mcparticle.vy() * TMath::Cos(trackParCov.getAlpha()) - mcparticle.vx() * TMath::Sin(trackParCov.getAlpha());

    // std::cout << "mcVy " <<  mcparticle.vy()   <<  std::endl;
    std::cout << "mcVxRotated " << mcVxRotated << std::endl;
    std::cout << "mcVyRotated " << mcVyRotated << std::endl;

    // std::array<float, 3>  arrayXYZ = { mcVxRotated, mcVyRotated , mcparticle.vz()};
    // std::array<float, 3>  arrayPxPyPz = {mcparticle.px(), mcparticle.py(), mcparticle.pz()};

    // const int matSize = 21;
    // std::array<float, matSize> arrayCovMatMC = {0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0. };
    // o2::track::TrackParametrizationWithError<float> trackParCovMC{arrayXYZ, arrayPxPyPz,  arrayCovMatMC, trackParCov.getSign()};

    // Double_t d0rpmc=parammc[0];
    // Double_t trackParDcaXYMC = trackParCovMC.getY(); // here

    Double_t trackParDcaXYMC = mcVyRotated; // here
    // std::cout << "trackParCovMC.getY() " <<  trackParCovMC.getY()   <<  std::endl;

    // Double_t d0zmc =parammc[1];
    Double_t trackParDcaZMC = mcparticle.vz();

    // Double_t dd0zo =d0zo-d0zmc;
    Double_t diffDcaZFromMCCurent = trackParDcaZCurrent - trackParDcaZMC;

    // Double_t dd0zn =dd0zo *(sd0zo >0. ? (sd0zn /sd0zo ) : 1.);
    Double_t diffDcaZFromMCUpgrded = diffDcaZFromMCCurent * (DcaZResCurrent > 0. ? (DcaZResUpgrded / DcaZResCurrent) : 1.);

    // Double_t d0zn  =d0zmc+dd0zn;
    Double_t trackParDcaZUpgrded = trackParDcaZMC + diffDcaZFromMCUpgrded;

    // Double_t dd0rpo=d0rpo-d0rpmc;
    Double_t diffDcaXYFromMCCurent = trackParDcaXYCurrent - trackParDcaXYMC;

    // Double_t dd0rpn=dd0rpo*(sd0rpo>0. ? (sd0rpn/sd0rpo) : 1.);
    Double_t diffDcaXYFromMCUpgrded = diffDcaXYFromMCCurent * (DcaXYResCurrent > 0. ? (DcaXYResUpgrded / DcaXYResCurrent) : 1.);

    // Double_t dd0mrpn=TMath::Abs(sd0mrpn)-TMath::Abs(sd0mrpo);
    // Double_t diffDcaXYMeanUpgMinusCur = TMath::Abs(DcaXYMeanUpgrded) - TMath::Abs(DcaXYMeanCurrent) ;
    Double_t diffDcaXYMeanUpgMinusCur = DcaXYMeanUpgrded - DcaXYMeanCurrent;

    // Double_t d0rpn =d0rpmc+dd0rpn-dd0mrpn;
    Double_t trackParDcaXYUpgrded = trackParDcaXYMC + diffDcaXYFromMCUpgrded - diffDcaXYMeanUpgMinusCur;

    std::cout << DcaZResCurrent << ", " << DcaZResUpgrded << ", diff(DcaZ - DcaZMC): " << diffDcaZFromMCCurent << ", diff upgraded: " << diffDcaZFromMCUpgrded << ", DcaZ Upgrded : " << trackParDcaZUpgrded << std::endl;
    std::cout << DcaXYResCurrent << ", " << DcaXYResUpgrded << ", diff(DcaY - DcaYMC): " << diffDcaXYFromMCCurent << ", diff upgraded: " << diffDcaXYFromMCUpgrded << ", DcaY Upgrded :" << trackParDcaXYUpgrded << std::endl;

    // option mimic data
    // ----------------------
    // if(fMimicData){
    //     // dd0mrpn=sd0mrpn-sd0mrpo;
    //     diffDcaXYMeanUpgMinusCur = DcaXYMeanUpgrded - DcaXYMeanCurrent;
    //     // d0rpn = d0rpmc+dd0rpn+dd0mrpn;
    //     trackParDcaXYUpgrded =  diffDcaXYFromMCCurent + diffDcaXYFromMCUpgrded + diffDcaXYMeanUpgMinusCur;
    // }

    // setting updated track parameters
    // --------------------------------
    // param[0]=d0rpn;
    // Double_t oldDCAxyValue = trackParCov.getY();
    trackParCov.setY(trackParDcaXYUpgrded);
    // trackParCov.setY(oldDCAxyValue);
    // param[1]=d0zn ;
    trackParCov.setZ(trackParDcaZUpgrded);

    if (isUpdateCurvature) {
      // --------------------------------------
      // Double_t dpt1o =pt1o-pt1mc;
      Double_t diffOneOverPtFromMCCurent = trackParOneOverPtCurrent - trackParOneOverPtMC;

      // Double_t dpt1n =dpt1o *(spt1o >0. ? (spt1n /spt1o ) : 1.);
      Double_t diffOneOverPtFromMCUpgrded = diffOneOverPtFromMCCurent * (OneOverPtCurrent > 0. ? (OneOverPtUpgrded / OneOverPtCurrent) : 1.);

      // Double_t pt1n  = pt1mc+dpt1n;
      Double_t trackParOneOverPtUpgrded = trackParOneOverPtMC + diffOneOverPtFromMCUpgrded;

      // param[4]=pt1n ;
      trackParCov.setQ2Pt(trackParOneOverPtUpgrded);
    }
    // std::cout << "Inside tuneTrackParams() before modifying trackParCov.getY(): " << trackParCov.getY() << " trackParOneOverPtMC = " << trackParOneOverPtMC << " diffOneOverPtFromMCUpgrded = "  << diffOneOverPtFromMCUpgrded <<  std::endl;

    // Updating Single Track Covariance matrices

    double sigmaY2 = 0.0;
    double sigmaZY = 0.0;
    double sigmaZ2 = 0.0;
    double sigmaSnpY = 0.0;
    double sigmaSnpZ = 0.0;
    double sigmaTglY = 0.0;
    double sigmaTglZ = 0.0;
    double sigma1PtY = 0.0;
    double sigma1PtZ = 0.0;
    double sigma1PtSnp = 0.0;
    double sigma1PtTgl = 0.0;
    double sigma1Pt2 = 0.0;

    if (isUpdateTrackCovMat) {
      //       if(sd0rpo>0.)            covar[0]*=(sd0rpn/sd0rpo)*(sd0rpn/sd0rpo);//yy
      sigmaY2 = trackParCov.getSigmaY2();
      if (DcaXYResCurrent > 0.)
        sigmaY2 *= ((DcaXYResUpgrded / DcaXYResCurrent) * (DcaXYResUpgrded / DcaXYResCurrent));
      trackParCov.setCov(sigmaY2, 0);

      //       if(sd0zo>0. && sd0rpo>0.)covar[1]*=(sd0rpn/sd0rpo)*(sd0zn/sd0zo);//yz
      sigmaZY = trackParCov.getSigmaZY();
      if (DcaZResCurrent > 0. && DcaXYResCurrent > 0.)
        sigmaZY *= ((DcaXYResUpgrded / DcaXYResCurrent) * (DcaXYResUpgrded / DcaXYResCurrent));
      trackParCov.setCov(sigmaZY, 1);

      //       if(sd0zo>0.)             covar[2]*=(sd0zn/sd0zo)*(sd0zn/sd0zo);//zz
      sigmaZ2 = trackParCov.getSigmaZ2();
      if (DcaZResCurrent > 0.)
        sigmaZ2 *= ((DcaZResUpgrded / DcaZResCurrent) * (DcaZResUpgrded / DcaZResCurrent));
      trackParCov.setCov(sigmaZ2, 2);

      //       if(sd0rpo>0.)            covar[3]*=(sd0rpn/sd0rpo);//yl
      sigmaSnpY = trackParCov.getSigmaSnpY();
      if (DcaXYResCurrent > 0.)
        sigmaSnpY *= ((DcaXYResUpgrded / DcaXYResCurrent));
      trackParCov.setCov(sigmaSnpY, 3);

      //       if(sd0zo>0.)             covar[4]*=(sd0zn/sd0zo);//zl
      sigmaSnpZ = trackParCov.getSigmaSnpZ();
      if (DcaZResCurrent > 0.)
        sigmaSnpZ *= ((DcaZResUpgrded / DcaZResCurrent));
      trackParCov.setCov(sigmaSnpZ, 4);

      //       if(sd0rpo>0.)            covar[6]*=(sd0rpn/sd0rpo);//ysenT
      sigmaTglY = trackParCov.getSigmaTglY();
      if (DcaXYResCurrent > 0.)
        sigmaTglY *= ((DcaXYResUpgrded / DcaXYResCurrent));
      trackParCov.setCov(sigmaTglY, 6);

      //       if(sd0zo>0.)             covar[7]*=(sd0zn/sd0zo);//zsenT
      sigmaTglZ = trackParCov.getSigmaTglZ();
      if (DcaZResCurrent > 0.)
        sigmaTglZ *= ((DcaZResUpgrded / DcaZResCurrent));
      trackParCov.setCov(sigmaTglZ, 7);

      //       if(sd0rpo>0. && spt1o>0.)covar[10]*=(sd0rpn/sd0rpo)*(spt1n/spt1o);//ypt
      sigma1PtY = trackParCov.getSigma1PtY();
      if (DcaXYResCurrent > 0. && OneOverPtCurrent > 0.)
        sigma1PtY *= ((DcaXYResUpgrded / DcaXYResCurrent) * (OneOverPtUpgrded / OneOverPtCurrent));
      trackParCov.setCov(sigma1PtY, 10);

      //       if(sd0zo>0. && spt1o>0.) covar[11]*=(sd0zn/sd0zo)*(spt1n/spt1o);//zpt
      sigma1PtZ = trackParCov.getSigma1PtZ();
      if (DcaZResCurrent > 0. && OneOverPtCurrent > 0.)
        sigma1PtZ *= ((DcaZResUpgrded / DcaZResCurrent) * (OneOverPtUpgrded / OneOverPtCurrent));
      trackParCov.setCov(sigma1PtZ, 11);

      //       if(spt1o>0.)             covar[12]*=(spt1n/spt1o);//sinPhipt
      sigma1PtSnp = trackParCov.getSigma1PtSnp();
      if (OneOverPtCurrent > 0.)
        sigma1PtSnp *= (OneOverPtUpgrded / OneOverPtCurrent);
      trackParCov.setCov(sigma1PtSnp, 12);

      //       if(spt1o>0.)             covar[13]*=(spt1n/spt1o);//tanTpt
      sigma1PtTgl = trackParCov.getSigma1PtTgl();
      if (OneOverPtCurrent > 0.)
        sigma1PtTgl *= (OneOverPtUpgrded / OneOverPtCurrent);
      trackParCov.setCov(sigma1PtTgl, 13);

      //       if(spt1o>0.)             covar[14]*=(spt1n/spt1o)*(spt1n/spt1o);//ptpt
      sigma1Pt2 = trackParCov.getSigma1Pt2();
      if (OneOverPtCurrent > 0.)
        sigma1Pt2 *= (OneOverPtUpgrded / OneOverPtCurrent);
      trackParCov.setCov(sigma1Pt2, 14);
    }

    if (isUpdatePulls) {
      double ratioDCAxyPulls = DcaXYPullCurrent / DcaXYPullUpgrded;
      double ratioDCAzPulls = DcaZPullCurrent / DcaZPullUpgrded;

      // covar[0]*=pullcorr*pullcorr;//yy

      sigmaY2 *= (ratioDCAxyPulls * ratioDCAxyPulls);
      trackParCov.setCov(sigmaY2, 0);

      // covar[1]*=pullcorr;//yz
      sigmaZY *= ratioDCAxyPulls;
      trackParCov.setCov(sigmaZY, 1);

      sigmaZ2 *= (ratioDCAzPulls * ratioDCAzPulls);
      trackParCov.setCov(sigmaZ2, 2);

      // covar[3]*=pullcorr;//yl
      sigmaSnpY *= ratioDCAxyPulls;
      trackParCov.setCov(sigmaSnpY, 3);

      sigmaSnpZ *= ratioDCAzPulls;
      trackParCov.setCov(sigmaSnpZ, 4);

      // covar[6]*=pullcorr;//ysenT
      sigmaTglY *= ratioDCAxyPulls;
      trackParCov.setCov(sigmaTglY, 6);

      sigmaTglZ *= ratioDCAzPulls;
      trackParCov.setCov(sigmaTglZ, 7);

      // covar[10]*=pullcorr;//ypt
      sigma1PtY *= ratioDCAxyPulls;
      trackParCov.setCov(sigma1PtY, 10);

      sigma1PtZ *= ratioDCAzPulls;
      trackParCov.setCov(sigma1PtZ, 11);
    }
  }

  // to be declared
  // ---------------
  Int_t getPhiBin(Double_t phi) const
  {
    Double_t pi = TMath::Pi();
    if (phi > 2. * pi || phi < 0.)
      return -1;
    if ((phi <= (pi / 4.)) || (phi > 7. * (pi / 4.)))
      return 0;
    if ((phi > (pi / 4.)) && (phi <= 3. * (pi / 4.)))
      return 1;
    if ((phi > 3. * (pi / 4.)) && (phi <= 5. * (pi / 4.)))
      return 2;
    if ((phi > (5. * pi / 4.)) && (phi <= 7. * (pi / 4.)))
      return 3;

    return -1;
  }

  Double_t EvalGraph(Double_t x, const TGraphErrors* graph) const
  {

    if (!graph) {
      printf("\tEvalGraph fails !\n");
      return 0.;
    }
    Int_t nPoints = graph->GetN();
    Double_t xMin = graph->GetX()[0];
    Double_t xMax = graph->GetX()[nPoints - 1];
    if (x > xMax)
      return graph->Eval(xMax);
    if (x < xMin)
      return graph->Eval(xMin);
    return graph->Eval(x);
  }
};
