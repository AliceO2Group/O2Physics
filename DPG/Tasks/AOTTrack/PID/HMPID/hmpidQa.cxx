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

#include "tableHMPID.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <TF1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>

#include <HMPIDBase/Param.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

double expectedSignal(double mass, float mom)
{
  // expected theoretical values

  double chAngleTh = -999.;

  const double nMean = 1.288; // radiator mean refraction index

  double cosChAngleTh = (TMath::Sqrt(mass * mass + mom * mom)) / (nMean * mom);
  if (cosChAngleTh > 1)
    return chAngleTh;

  chAngleTh = TMath::ACos(cosChAngleTh);

  return chAngleTh;
}

double expectedSigma(int iPart, float mom)
{
  double sigmaRing = 0;
  const int nMaxSigmas = 6;

  // fit parameters from sigma extrapolation
  double fitSigmaExtraPions[nMaxSigmas] = {42.8961, -49.8723, 26.2311, -6.59093, 0.754578, -0.0286546};
  double fitSigmaExtraKaons[nMaxSigmas] = {76.9786, -103.655, 61.7533, -18.7436, 2.87855, -0.178318};
  double fitSigmaExtraProtons[nMaxSigmas] = {299.466, -383.277, 201.127, -52.2554, 6.67285, -0.334106};

  // create sigma vs p functions
  TF1* sigmaVsMomPions = new TF1("sigmaVsMomPions", "pol5", 0., 6.);
  TF1* sigmaVsMomKaons = new TF1("sigmaVsMomKaons", "pol5", 0., 6.);
  TF1* sigmaVsMomProtons = new TF1("sigmaVsMomProtons", "pol5", 0., 6.);

  for (int i = 0; i < nMaxSigmas; i++) {
    sigmaVsMomPions->SetParameter(i, fitSigmaExtraPions[i]);
    sigmaVsMomKaons->SetParameter(i, fitSigmaExtraKaons[i]);
    sigmaVsMomProtons->SetParameter(i, fitSigmaExtraProtons[i]);
  }

  const int idPions = 0, idKaons = 1, idProtons = 2;

  if (iPart == idPions) {
    sigmaRing = gRandom->Gaus(sigmaVsMomPions->Eval(mom) / 1000., 0.1 * sigmaVsMomPions->Eval(mom) / 1000.);
  }
  if (iPart == idKaons) {
    sigmaRing = 0.8 * gRandom->Gaus(sigmaVsMomKaons->Eval(mom) / 1000., 0.15 * sigmaVsMomKaons->Eval(mom) / 1000.);
  }
  if (iPart == idProtons) {
    sigmaRing = 0.6 * gRandom->Gaus(sigmaVsMomProtons->Eval(mom) / 1000., 0.1 * sigmaVsMomProtons->Eval(mom) / 1000.);
  }

  delete sigmaVsMomPions;
  delete sigmaVsMomKaons;
  delete sigmaVsMomProtons;

  return sigmaRing;
}

void getProbability(float hmpidSignal, float hmpidMomentum, double* probs)
{
  // Calculates probability to be a pion-kaon-proton with the "amplitude" method
  // from the given Cerenkov angle and momentum assuming no initial particle composition (class taken by AliROOT)

  int nSpecies = 3;
  float angleZeroHeight = 900.;

  if (hmpidSignal <= 0) {
    // HMPID does not find anything reasonable for this track, assign 0.33 for all species
    for (int iPart = 0; iPart < nSpecies; iPart++)
      probs[iPart] = 1.0 / nSpecies;
    return;
  }

  // assign mass in GeV/c^2
  double mass[] = {o2::constants::physics::MassPionCharged, o2::constants::physics::MassKaonCharged, o2::constants::physics::MassProton};

  double hTot = 0;                  // Initialize the total height of the amplitude method
  double* h = new double[nSpecies]; // number of charged particles to be considered

  bool desert = kTRUE; // Flag to evaluate if ThetaC is far ("desert") from the given Gaussians

  for (int iPart = 0; iPart < nSpecies; iPart++) { // for each particle

    h[iPart] = 0;                                                   // reset the height
    double thetaCerTh = expectedSignal(mass[iPart], hmpidMomentum); // theoretical Theta Cherenkov
    if (thetaCerTh > angleZeroHeight)
      continue; // no light emitted, zero height
    double sigmaRing = expectedSigma(iPart, hmpidMomentum);
    float maxSigmaRing = 4 * sigmaRing;

    if (sigmaRing == 0)
      continue;

    if (TMath::Abs(hmpidSignal - thetaCerTh) < maxSigmaRing)
      desert = kFALSE;
    h[iPart] = TMath::Gaus(thetaCerTh, hmpidSignal, sigmaRing, kTRUE);
    hTot += h[iPart]; // total height of all theoretical heights for normalization

  } // species loop

  for (int iPart = 0; iPart < nSpecies; iPart++) { // species loop to assign probabilities

    if (!desert)
      probs[iPart] = h[iPart] / hTot;
    else
      probs[iPart] = 1.0 / nSpecies; // all theoretical values are far away from experemental one
  }

  delete[] h;
}

struct HmpidQa {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsP{"nBinsP", 1200, "Number of momentum bins"};
  Configurable<float> minP{"minP", -20.f, "Minimum momentum plotted (GeV/c)"};
  Configurable<float> maxP{"maxP", 20.f, "Maximum momentum plotted (GeV/c)"};

  Configurable<int> nBinsCh{"nBinsCh", 800, "Number of ch angle bins"};
  Configurable<float> minCh{"minCh", 0.f, "Minimum ch angle plotted (rad)"};
  Configurable<float> maxCh{"maxCh", 0.8f, "Maximum ch angle plotted (rad)"};

  /// filters configurables for primary tracks
  Configurable<float> nSigmaTpcMin{"nSigmaTpcMin", -3.0, "nSigmaTpcMin"};
  Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", +3.0, "nSigmaTpcMax"};
  Configurable<float> nSigmaTofMin{"nSigmaTofMin", -3.0, "nSigmaTofMin"};
  Configurable<float> nSigmaTofMax{"nSigmaTofMax", +3.5, "nSigmaTofMax"};
  Configurable<float> minReqClusterIts{"minReqClusterIts", 4.0, "min number of clusters required in ITS"};
  Configurable<float> minTpcNClsFound{"minTpcNClsFound", 50.0f, "minTpcNClsFound"};
  Configurable<float> minNCrossedRowsTpc{"minNCrossedRowsTpc", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> maxChi2Its{"maxChi2Its", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2Tpc{"maxChi2Tpc", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAXY{"maxDCAXY", 0.5f, "maxDCAXY"};
  Configurable<float> maxDCAZ{"maxDCAZ", 0.5f, "maxDCAZ"};

  // QA filters
  Configurable<float> cutDistanceMipTrack{"cutDistanceMipTrack", 3.0f, "cut distance between MIP and track"};
  Configurable<float> cutQmip{"cutQmip", 120.0f, "cut on Q MIP"};
  Configurable<float> cutMinMomGlobalTrack{"cutMinMomGlobalTrack", 1.5f, "minimum momentum of global track"};
  Configurable<float> minProbParticle{"minProbParticle", 0.7f, "minimum particle probability"};
  Configurable<float> maxDistanceForProb{"maxDistanceForProb", 1.5f, "maximum distance for probability calculation"};
  Configurable<float> maxBoxHit{"maxBoxHit", 100.0f, "maximum box hit position"};
  Configurable<float> minBoxHit{"minBoxHit", 40.0f, "minimum box hit position"};

  // variables for chamber_number and HVs/PCs
  const int rich0 = 0, rich1 = 1, rich2 = 2, rich3 = 3, rich4 = 4, rich5 = 5, rich6 = 6;
  const int hv0 = 0, hv1 = 1, hv2 = 2, hv3 = 3, hv4 = 4, hv5 = 5;
  // total number of chambers and HVs/PCs
  static const int nCh = 7, nSec = 6, nPc = 6;

  ////////////////////////////////////
  ////////////////////////////////////
  // load geometry
  o2::hmpid::Param* fParam = o2::hmpid::Param::instanceNoGeo();

  void init(InitContext const&)
  {
    AxisSpec momAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    AxisSpec cherenkAxis{nBinsCh, minCh, maxCh, "#theta_{Ch} (rad)"};

    histos.add("nPhotons_vs_sin2Ch", "nPhotons_vs_sin2Ch", kTProfile, {{40, 0.0, 0.5}});

    histos.add("ChAngle_LowPt", "ChAngle_LowPt", kTH1F, {cherenkAxis});
    histos.add("ChAngle_HighPt", "ChAngle_HighPt", kTH1F, {cherenkAxis});

    histos.add("hmpidSignal", "hmpidSignal", kTH1F, {cherenkAxis});

    // th2f for spectra
    histos.add("pTvsChAngle", "pTvsChAngle", kTH2F, {{500, 0, 10., "#it{p}_{T} (GeV/#it{c})"}, {cherenkAxis}});

    // charge identification
    histos.add("pTvsChAnglePos", "pTvsChAnglePos", kTH2F, {{500, 0, 10., "#it{p}_{T} (GeV/#it{c})"}, {cherenkAxis}});
    histos.add("pTvsChAngleNeg", "pTvsChAngleNeg", kTH2F, {{500, 0, 10., "#it{p}_{T} (GeV/#it{c})"}, {cherenkAxis}});

    histos.add("hmpidMomvsTrackMom", "hmpidMomvsTrackMom", kTH2F, {{1200, 0, 30, "Track #it{p} (GeV/#it{c})"}, {1200, 0, 30, "HMPID #it{p} (GeV/#it{c})"}});
    histos.add("hmpidCkovvsMom", "hmpidCkovvsMom", kTH2F, {{1000, 0, 10, "#it{p} (GeV/#it{c})"}, cherenkAxis});
    histos.add("TrackMom", "TrackMom", kTH1F, {momAxis});
    histos.add("hmpidMom", "hmpidMom", kTH1F, {momAxis});

    histos.add("hmpidNPhotons", "hmpidNPhotons", kTH1F, {{50, 2, 50, "Number of photons"}});

    histos.add("hmpidCkovvsMom_nocut", "hmpidCkovvsMom_nocut", kTH2F, {{1000, 0, 10, "#it{p} (GeV/#it{c})"}, cherenkAxis});

    histos.add("hmpidPhotsCharge", "hmpidPhotsCharge", kTH1F, {{180, 4, 210}});
    histos.add("hmpidQMip", "hmpidQMip", kTH1F, {{1000, 200, 2200, "Charge (ADC)"}});

    // information on particle position
    histos.add("hmpidXTrack", "hmpidXTrack", kTH1F, {{270, 0, 135, "X track (cm)"}});
    histos.add("hmpidYTrack", "hmpidYTrack", kTH1F, {{270, 0, 135, "Y track (cm)"}});
    histos.add("hmpidXMip", "hmpidXMip", kTH1F, {{270, 0, 135, "X mip (cm)"}});
    histos.add("hmpidYMip", "hmpidYMip", kTH1F, {{270, 0, 135, "X mip (cm)"}});
    histos.add("hmpidXResiduals", "hmpidXResiduals", kTH1F, {{400, -20, 20, "X Residuals (cm)"}});
    histos.add("hmpidYResiduals", "hmpidYResiduals", kTH1F, {{400, -20, 20, "Y Residuals (cm)"}});
    // 2D map for the mip and the track
    histos.add("hmpidXYTrack", "hmpidXYTrack", kTH2F, {{270, 0, 135, "X track (cm)"}, {270, 0, 135, "Y track (cm)"}});
    histos.add("hmpidXYMip", "hmpidXYMip", kTH2F, {{270, 0, 135, "X mip (cm)"}, {270, 0, 135, "Y mip (cm)"}});

    // histos per chamber
    for (int iCh = 0; iCh < nCh; iCh++) {
      histos.add(Form("hmpidXTrack%i", iCh), Form("hmpidXTrack%i", iCh), kTH1F, {{270, 0, 135, "X track (cm)"}});
      histos.add(Form("hmpidYTrack%i", iCh), Form("hmpidYTrack%i", iCh), kTH1F, {{270, 0, 135, "Y track (cm)"}});
      histos.add(Form("hmpidXMip%i", iCh), Form("hmpidXMip%i", iCh), kTH1F, {{270, 0, 135, "X mip (cm)"}});
      histos.add(Form("hmpidYMip%i", iCh), Form("hmpidYMip%i", iCh), kTH1F, {{270, 0, 135, "X mip (cm)"}});
      histos.add(Form("hmpidXResiduals%i", iCh), Form("hmpidXResiduals%i", iCh), kTH1F, {{400, -20, 20, "X Residuals (cm)"}});
      histos.add(Form("hmpidYResiduals%i", iCh), Form("hmpidYResiduals%i", iCh), kTH1F, {{400, -20, 20, "Y Residuals (cm)"}});

      // residuals discriminated for charge sign
      histos.add(Form("hmpidXResidualsPos%i", iCh), Form("hmpidXResidualsPos%i", iCh), kTH1F, {{400, -20, 20, "X Residuals (cm)"}});
      histos.add(Form("hmpidYResidualsPos%i", iCh), Form("hmpidYResidualsPos%i", iCh), kTH1F, {{400, -20, 20, "Y Residuals (cm)"}});

      histos.add(Form("hmpidXResidualsNeg%i", iCh), Form("hmpidXResidualsNeg%i", iCh), kTH1F, {{400, -20, 20, "X Residuals (cm)"}});
      histos.add(Form("hmpidYResidualsNeg%i", iCh), Form("hmpidYResidualsNeg%i", iCh), kTH1F, {{400, -20, 20, "Y Residuals (cm)"}});

      histos.add(Form("hmpidNPhotons%i", iCh), Form("hmpidNPhotons%i", iCh), kTH1F, {{50, 2, 50, "Number of photons"}});

      histos.add(Form("hmpidQMip%i", iCh), Form("hmpidQMip%i", iCh), kTH1F, {{1000, 200, 2200, "Charge (ADC)"}});
      histos.add(Form("hmpidClusSize%i", iCh), Form("hmpidClusSize%i", iCh), kTH1F, {{15, 0, 15, "MIP Cluster size"}});

      histos.add(Form("TrackMom%i", iCh), Form("TrackMom%i", iCh), kTH1F, {momAxis});
      histos.add(Form("hmpidMom%i", iCh), Form("hmpidMom%i", iCh), kTH1F, {momAxis});

      histos.add(Form("hmpidPhotsCharge%i", iCh), Form("hmpidPhotsCharge%i", iCh), kTH1F, {{180, 4, 210}});
      histos.add(Form("hmpidXYMip%i", iCh), Form("hmpidXYMip%i", iCh), kTH2F, {{270, 0, 135, "X mip (cm)"}, {270, 0, 135, "Y mip (cm)"}});

      histos.add(Form("nPhotons_vs_sin2Ch%i", iCh), Form("N. of Photons vs sin^{2}(#theta_{Ch}) - chamber%i", iCh), kTProfile, {{40, 0.0, 0.5}});

      // histos per HV sector
      for (int iSec = 0; iSec < nSec; iSec++) {
        histos.add(Form("hmpidQMip_RICH%i_HV%i", iCh, iSec), Form("hmpidQMip_RICH%i_HV%i", iCh, iSec), kTH1F, {{2000, 200, 2200, "Charge (ADC)"}});
        histos.add(Form("hmpidNPhotons_RICH%i_HV%i", iCh, iSec), Form("hmpidNPhotons_RICH%i_HV%i", iCh, iSec), kTH1F, {{50, 2, 50, "Number of photons"}});
        histos.add(Form("hmpidPhotsCharge_RICH%i_HV%i", iCh, iSec), Form("hmpidPhotsCharge_RICH%i_HV%i", iCh, iSec), kTH1F, {{180, 4, 210}});
      }

      // plot n_ph vs sin2Ch per PC
      for (int iPc = 0; iPc < nPc; iPc++) {
        histos.add(Form("nPhotons_vs_sin2Ch_RICH%i_PC%i", iCh, iPc), Form("N. of Photons vs sin^{2}(#theta_{Ch}) - chamber%i, photocathode%i", iCh, iPc), kTProfile, {{20, 0.0, 0.4}});
      }
    }
  }

  void process(aod::HmpidAnalysis const& hmpidtable)
  {
    // photocathods limits
    static float xMinPc[nPc];
    static float yMinPc[nPc];
    static float xMaxPc[nPc];
    static float yMaxPc[nPc];

    for (int iPc = 0; iPc < nPc; iPc++) {
      xMaxPc[iPc] = (fParam->maxPcX(iPc)) - 10.;
      yMaxPc[iPc] = (fParam->maxPcY(iPc)) - 10.;
      xMinPc[iPc] = (fParam->minPcX(iPc)) + 10.;
      yMinPc[iPc] = (fParam->minPcY(iPc)) + 10.;
    }

    for (const auto& hmpid : hmpidtable) // loop on tracks contained in the table
    {

      // filters on primary tracks
      if (hmpid.itsNCluster() < minReqClusterIts)
        continue;
      if (hmpid.tpcNCluster() < minTpcNClsFound)
        continue;
      if (hmpid.tpcNClsCrossedRows() < minNCrossedRowsTpc)
        continue;
      if (hmpid.tpcChi2() > maxChi2Tpc)
        continue;
      if (hmpid.itsChi2() > maxChi2Its)
        continue;
      if (TMath::Abs(hmpid.dcaXY()) > maxDCAXY)
        continue;
      if (TMath::Abs(hmpid.dcaZ()) > maxDCAZ)
        continue;

      // evaluate distance mip-track
      const float distanceMipToTrack = std::hypot(hmpid.xTrack() - hmpid.xMip(), hmpid.yTrack() - hmpid.yMip());

      // quality conditions to check
      const bool physicalChAngle = (hmpid.chAngle() > 0);
      const bool mipChargeCondition = (hmpid.chargeMip() > cutQmip);
      const bool distanceCondition = (distanceMipToTrack < cutDistanceMipTrack);

      // fill histograms
      histos.fill(HIST("hmpidMomvsTrackMom"), std::fabs(hmpid.momentumTrack()), std::fabs(hmpid.momentumHmpid()));
      histos.fill(HIST("TrackMom"), std::fabs(hmpid.momentumTrack()));
      histos.fill(HIST("hmpidMom"), std::fabs(hmpid.momentumHmpid()));

      histos.fill(HIST("hmpidSignal"), hmpid.chAngle());

      if (physicalChAngle && distanceCondition && mipChargeCondition) {
        double pT = static_cast<double>(hmpid.momentumTrack() / TMath::CosH(hmpid.etaTrack()));
        histos.fill(HIST("pTvsChAngle"), pT, hmpid.chAngle());
        if (hmpid.momentumHmpid() > 0) {
          histos.fill(HIST("pTvsChAnglePos"), pT, hmpid.chAngle());
        }
        if (hmpid.momentumHmpid() < 0) {
          histos.fill(HIST("pTvsChAngleNeg"), pT, hmpid.chAngle());
        }
      }

      float sin2changle = 0.;

      if (distanceCondition && mipChargeCondition) {

        histos.fill(HIST("hmpidNPhotons"), hmpid.nPhotons());

        sin2changle = static_cast<float>(TMath::Power(TMath::Sin(hmpid.chAngle()), 2));
        if (hmpid.xMip() <= maxBoxHit && hmpid.xMip() >= minBoxHit && hmpid.yMip() <= maxBoxHit && hmpid.yMip() >= minBoxHit) {
          histos.fill(HIST("nPhotons_vs_sin2Ch"), sin2changle, hmpid.nPhotons());
        }

        for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
          if (hmpid.photonsCharge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge"), hmpid.photonsCharge()[i]);
        }
      }

      if (distanceCondition) {
        histos.fill(HIST("hmpidQMip"), hmpid.chargeMip());
      }

      histos.fill(HIST("hmpidXTrack"), hmpid.xTrack());
      histos.fill(HIST("hmpidYTrack"), hmpid.yTrack());
      histos.fill(HIST("hmpidXMip"), hmpid.xMip());
      histos.fill(HIST("hmpidYMip"), hmpid.yMip());
      if (hmpid.momentumTrack() > cutMinMomGlobalTrack) {
        histos.fill(HIST("hmpidXResiduals"), hmpid.xMip() - hmpid.xTrack());
        histos.fill(HIST("hmpidYResiduals"), hmpid.yMip() - hmpid.yTrack());
      }
      histos.fill(HIST("hmpidXYTrack"), hmpid.xTrack(), hmpid.yTrack());
      histos.fill(HIST("hmpidXYMip"), hmpid.xMip(), hmpid.yMip());

      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      // fill histograms per chamber
      if (hmpid.chamber() == rich0) {
        histos.fill(HIST("hmpidXTrack0"), hmpid.xTrack());
        histos.fill(HIST("hmpidYTrack0"), hmpid.yTrack());
        histos.fill(HIST("hmpidXMip0"), hmpid.xMip());
        histos.fill(HIST("hmpidYMip0"), hmpid.yMip());
        histos.fill(HIST("hmpidXResiduals0"), hmpid.xMip() - hmpid.xTrack());
        histos.fill(HIST("hmpidYResiduals0"), hmpid.yMip() - hmpid.yTrack());

        if (hmpid.momentumTrack() > cutMinMomGlobalTrack) {
          if (hmpid.momentumHmpid() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos0"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsPos0"), hmpid.yMip() - hmpid.yTrack());
          }

          if (hmpid.momentumHmpid() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg0"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsNeg0"), hmpid.yMip() - hmpid.yTrack());
          }
        }

        if (distanceCondition) {
          histos.fill(HIST("hmpidQMip0"), hmpid.chargeMip());
        }
        histos.fill(HIST("hmpidClusSize0"), hmpid.clusterSize());
        histos.fill(HIST("TrackMom0"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom0"), std::fabs(hmpid.momentumHmpid()));
        histos.fill(HIST("hmpidXYMip0"), hmpid.xMip(), hmpid.yMip());

        if (distanceCondition && mipChargeCondition) {
          histos.fill(HIST("hmpidNPhotons0"), hmpid.nPhotons());
          sin2changle = static_cast<float>(TMath::Power(TMath::Sin(hmpid.chAngle()), 2));
          if (hmpid.xMip() <= maxBoxHit && hmpid.xMip() >= minBoxHit && hmpid.yMip() <= maxBoxHit && hmpid.yMip() >= minBoxHit) {
            histos.fill(HIST("nPhotons_vs_sin2Ch0"), sin2changle, hmpid.nPhotons());
          }

          for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
            if (hmpid.photonsCharge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge0"), hmpid.photonsCharge()[i]);
          }
        }

        //////////////////////////////////////////////////////////////////
        // plot per HV sector
        if (fParam->inHVSector(hmpid.yMip()) == hv0) {

          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH0_HV0"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV0"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV0"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv1) {

          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH0_HV1"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV1"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV1"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv2) {

          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH0_HV2"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV2"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV2"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv3) {

          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH0_HV3"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV3"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV3"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv4) {

          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH0_HV4"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV4"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV4"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv5) {

          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH0_HV5"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV5"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV5"), hmpid.nPhotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distanceCondition && mipChargeCondition && physicalChAngle) // condizione da verificare a priori
        {
          if (hmpid.xMip() >= xMinPc[0] && hmpid.xMip() <= xMaxPc[0] && hmpid.yMip() >= yMinPc[0] && hmpid.yMip() <= yMaxPc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC0"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[1] && hmpid.xMip() <= xMaxPc[1] && hmpid.yMip() >= yMinPc[1] && hmpid.yMip() <= yMaxPc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC1"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[2] && hmpid.xMip() <= xMaxPc[2] && hmpid.yMip() >= yMinPc[2] && hmpid.yMip() <= yMaxPc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC2"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[3] && hmpid.xMip() <= xMaxPc[3] && hmpid.yMip() >= yMinPc[3] && hmpid.yMip() <= yMaxPc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC3"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[4] && hmpid.xMip() <= xMaxPc[4] && hmpid.yMip() >= yMinPc[4] && hmpid.yMip() <= yMaxPc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC4"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[5] && hmpid.xMip() <= xMaxPc[5] && hmpid.yMip() >= yMinPc[5] && hmpid.yMip() <= yMaxPc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC5"), sin2changle, hmpid.nPhotons());
          }
        }
      }

      if (hmpid.chamber() == rich1) {
        histos.fill(HIST("hmpidXTrack1"), hmpid.xTrack());
        histos.fill(HIST("hmpidYTrack1"), hmpid.yTrack());
        histos.fill(HIST("hmpidXMip1"), hmpid.xMip());
        histos.fill(HIST("hmpidYMip1"), hmpid.yMip());
        histos.fill(HIST("hmpidXResiduals1"), hmpid.xMip() - hmpid.xTrack());
        histos.fill(HIST("hmpidYResiduals1"), hmpid.yMip() - hmpid.yTrack());

        if (hmpid.momentumTrack() > cutMinMomGlobalTrack) {
          if (hmpid.momentumHmpid() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos1"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsPos1"), hmpid.yMip() - hmpid.yTrack());
          }

          if (hmpid.momentumHmpid() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg1"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsNeg1"), hmpid.yMip() - hmpid.yTrack());
          }
        }

        if (distanceCondition) {
          histos.fill(HIST("hmpidQMip1"), hmpid.chargeMip());
        }
        histos.fill(HIST("hmpidClusSize1"), hmpid.clusterSize());
        histos.fill(HIST("TrackMom1"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom1"), std::fabs(hmpid.momentumHmpid()));
        histos.fill(HIST("hmpidXYMip1"), hmpid.xMip(), hmpid.yMip());

        if (distanceCondition && mipChargeCondition) {
          histos.fill(HIST("hmpidNPhotons1"), hmpid.nPhotons());
          sin2changle = static_cast<float>(TMath::Power(TMath::Sin(hmpid.chAngle()), 2));
          if (hmpid.xMip() <= maxBoxHit && hmpid.xMip() >= minBoxHit && hmpid.yMip() <= maxBoxHit && hmpid.yMip() >= minBoxHit) {
            histos.fill(HIST("nPhotons_vs_sin2Ch1"), sin2changle, hmpid.nPhotons());
          }

          for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
            if (hmpid.photonsCharge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge1"), hmpid.photonsCharge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.yMip()) == hv0) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH1_HV0"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV0"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV0"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv1) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH1_HV1"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV1"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV1"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv2) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH1_HV2"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV2"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV2"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv3) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH1_HV3"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV3"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV3"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv4) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH1_HV4"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV4"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV4"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv5) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH1_HV5"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV5"), hmpid.photonsCharge()[i]);
              } // hmpidPhotsCharge_RICH%i_HV%i
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV5"), hmpid.nPhotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distanceCondition && mipChargeCondition && physicalChAngle) // condizione da verificare a priori
        {
          if (hmpid.xMip() >= xMinPc[0] && hmpid.xMip() <= xMaxPc[0] && hmpid.yMip() >= yMinPc[0] && hmpid.yMip() <= yMaxPc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC0"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[1] && hmpid.xMip() <= xMaxPc[1] && hmpid.yMip() >= yMinPc[1] && hmpid.yMip() <= yMaxPc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC1"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[2] && hmpid.xMip() <= xMaxPc[2] && hmpid.yMip() >= yMinPc[2] && hmpid.yMip() <= yMaxPc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC2"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[3] && hmpid.xMip() <= xMaxPc[3] && hmpid.yMip() >= yMinPc[3] && hmpid.yMip() <= yMaxPc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC3"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[4] && hmpid.xMip() <= xMaxPc[4] && hmpid.yMip() >= yMinPc[4] && hmpid.yMip() <= yMaxPc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC4"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[5] && hmpid.xMip() <= xMaxPc[5] && hmpid.yMip() >= yMinPc[5] && hmpid.yMip() <= yMaxPc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC5"), sin2changle, hmpid.nPhotons());
          }
        }
      }

      if (hmpid.chamber() == rich2) {
        histos.fill(HIST("hmpidXTrack2"), hmpid.xTrack());
        histos.fill(HIST("hmpidYTrack2"), hmpid.yTrack());
        histos.fill(HIST("hmpidXMip2"), hmpid.xMip());
        histos.fill(HIST("hmpidYMip2"), hmpid.yMip());
        histos.fill(HIST("hmpidXResiduals2"), hmpid.xMip() - hmpid.xTrack());
        histos.fill(HIST("hmpidYResiduals2"), hmpid.yMip() - hmpid.yTrack());

        if (hmpid.momentumTrack() > cutMinMomGlobalTrack) {
          if (hmpid.momentumHmpid() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos2"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsPos2"), hmpid.yMip() - hmpid.yTrack());
          }

          if (hmpid.momentumHmpid() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg2"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsNeg2"), hmpid.yMip() - hmpid.yTrack());
          }
        }

        if (distanceCondition) {
          histos.fill(HIST("hmpidQMip2"), hmpid.chargeMip());
        }
        histos.fill(HIST("hmpidClusSize2"), hmpid.clusterSize());
        histos.fill(HIST("TrackMom2"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom2"), std::fabs(hmpid.momentumHmpid()));
        histos.fill(HIST("hmpidXYMip2"), hmpid.xMip(), hmpid.yMip());

        if (distanceCondition && mipChargeCondition) {
          histos.fill(HIST("hmpidNPhotons2"), hmpid.nPhotons());
          sin2changle = static_cast<float>(TMath::Power(TMath::Sin(hmpid.chAngle()), 2));
          if (hmpid.xMip() <= maxBoxHit && hmpid.xMip() >= minBoxHit && hmpid.yMip() <= maxBoxHit && hmpid.yMip() >= minBoxHit) {
            histos.fill(HIST("nPhotons_vs_sin2Ch2"), sin2changle, hmpid.nPhotons());
          }

          for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
            if (hmpid.photonsCharge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge2"), hmpid.photonsCharge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.yMip()) == hv0) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH2_HV0"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV0"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV0"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv1) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH2_HV1"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV1"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV1"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv2) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH2_HV2"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV2"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV2"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv3) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH2_HV3"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV3"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV3"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv4) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH2_HV4"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV4"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV4"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv5) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH2_HV5"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV5"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV5"), hmpid.nPhotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distanceCondition && mipChargeCondition && physicalChAngle) // condizione da verificare a priori
        {
          if (hmpid.xMip() >= xMinPc[0] && hmpid.xMip() <= xMaxPc[0] && hmpid.yMip() >= yMinPc[0] && hmpid.yMip() <= yMaxPc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC0"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[1] && hmpid.xMip() <= xMaxPc[1] && hmpid.yMip() >= yMinPc[1] && hmpid.yMip() <= yMaxPc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC1"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[2] && hmpid.xMip() <= xMaxPc[2] && hmpid.yMip() >= yMinPc[2] && hmpid.yMip() <= yMaxPc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC2"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[3] && hmpid.xMip() <= xMaxPc[3] && hmpid.yMip() >= yMinPc[3] && hmpid.yMip() <= yMaxPc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC3"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[4] && hmpid.xMip() <= xMaxPc[4] && hmpid.yMip() >= yMinPc[4] && hmpid.yMip() <= yMaxPc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC4"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[5] && hmpid.xMip() <= xMaxPc[5] && hmpid.yMip() >= yMinPc[5] && hmpid.yMip() <= yMaxPc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC5"), sin2changle, hmpid.nPhotons());
          }
        }
      }

      if (hmpid.chamber() == rich3) {
        histos.fill(HIST("hmpidXTrack3"), hmpid.xTrack());
        histos.fill(HIST("hmpidYTrack3"), hmpid.yTrack());
        histos.fill(HIST("hmpidXMip3"), hmpid.xMip());
        histos.fill(HIST("hmpidYMip3"), hmpid.yMip());
        histos.fill(HIST("hmpidXResiduals3"), hmpid.xMip() - hmpid.xTrack());
        histos.fill(HIST("hmpidYResiduals3"), hmpid.yMip() - hmpid.yTrack());

        if (hmpid.momentumTrack() > cutMinMomGlobalTrack) {
          if (hmpid.momentumHmpid() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos3"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsPos3"), hmpid.yMip() - hmpid.yTrack());
          }

          if (hmpid.momentumHmpid() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg3"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsNeg3"), hmpid.yMip() - hmpid.yTrack());
          }
        }

        if (distanceCondition) {
          histos.fill(HIST("hmpidQMip3"), hmpid.chargeMip());
        }
        histos.fill(HIST("hmpidClusSize3"), hmpid.clusterSize());
        histos.fill(HIST("TrackMom3"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom3"), std::fabs(hmpid.momentumHmpid()));
        histos.fill(HIST("hmpidXYMip3"), hmpid.xMip(), hmpid.yMip());

        if (distanceCondition && mipChargeCondition) {
          histos.fill(HIST("hmpidNPhotons3"), hmpid.nPhotons());
          sin2changle = static_cast<float>(TMath::Power(TMath::Sin(hmpid.chAngle()), 2));
          if (hmpid.xMip() <= maxBoxHit && hmpid.xMip() >= minBoxHit && hmpid.yMip() <= maxBoxHit && hmpid.yMip() >= minBoxHit) {
            histos.fill(HIST("nPhotons_vs_sin2Ch3"), sin2changle, hmpid.nPhotons());
          }

          for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
            if (hmpid.photonsCharge()[i] > 0) {
              histos.fill(HIST("hmpidPhotsCharge3"), hmpid.photonsCharge()[i]);
            }
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.yMip()) == hv0) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH3_HV0"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV0"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV0"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv1) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH3_HV1"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV1"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV1"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv2) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH3_HV2"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV2"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV2"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv3) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH3_HV3"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV3"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV3"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv4) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH3_HV4"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV4"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV4"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv5) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH3_HV5"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV5"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV5"), hmpid.nPhotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distanceCondition && mipChargeCondition && physicalChAngle) // condizione da verificare a priori
        {
          if (hmpid.xMip() >= xMinPc[0] && hmpid.xMip() <= xMaxPc[0] && hmpid.yMip() >= yMinPc[0] && hmpid.yMip() <= yMaxPc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC0"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[1] && hmpid.xMip() <= xMaxPc[1] && hmpid.yMip() >= yMinPc[1] && hmpid.yMip() <= yMaxPc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC1"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[2] && hmpid.xMip() <= xMaxPc[2] && hmpid.yMip() >= yMinPc[2] && hmpid.yMip() <= yMaxPc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC2"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[3] && hmpid.xMip() <= xMaxPc[3] && hmpid.yMip() >= yMinPc[3] && hmpid.yMip() <= yMaxPc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC3"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[4] && hmpid.xMip() <= xMaxPc[4] && hmpid.yMip() >= yMinPc[4] && hmpid.yMip() <= yMaxPc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC4"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[5] && hmpid.xMip() <= xMaxPc[5] && hmpid.yMip() >= yMinPc[5] && hmpid.yMip() <= yMaxPc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC5"), sin2changle, hmpid.nPhotons());
          }
        }
      }

      if (hmpid.chamber() == rich4) {
        histos.fill(HIST("hmpidXTrack4"), hmpid.xTrack());
        histos.fill(HIST("hmpidYTrack4"), hmpid.yTrack());
        histos.fill(HIST("hmpidXMip4"), hmpid.xMip());
        histos.fill(HIST("hmpidYMip4"), hmpid.yMip());
        histos.fill(HIST("hmpidXResiduals4"), hmpid.xMip() - hmpid.xTrack());
        histos.fill(HIST("hmpidYResiduals4"), hmpid.yMip() - hmpid.yTrack());

        if (hmpid.momentumTrack() > cutMinMomGlobalTrack) {
          if (hmpid.momentumHmpid() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos4"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsPos4"), hmpid.yMip() - hmpid.yTrack());
          }

          if (hmpid.momentumHmpid() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg4"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsNeg4"), hmpid.yMip() - hmpid.yTrack());
          }
        }

        if (distanceCondition) {
          histos.fill(HIST("hmpidQMip4"), hmpid.chargeMip());
        }
        histos.fill(HIST("hmpidClusSize4"), hmpid.clusterSize());
        histos.fill(HIST("TrackMom4"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom4"), std::fabs(hmpid.momentumHmpid()));
        histos.fill(HIST("hmpidXYMip4"), hmpid.xMip(), hmpid.yMip());

        if (distanceCondition && mipChargeCondition) {
          histos.fill(HIST("hmpidNPhotons4"), hmpid.nPhotons());
          sin2changle = static_cast<float>(TMath::Power(TMath::Sin(hmpid.chAngle()), 2));
          if (hmpid.xMip() <= maxBoxHit && hmpid.xMip() >= minBoxHit && hmpid.yMip() <= maxBoxHit && hmpid.yMip() >= minBoxHit) {
            histos.fill(HIST("nPhotons_vs_sin2Ch4"), sin2changle, hmpid.nPhotons());
          }

          for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
            if (hmpid.photonsCharge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge4"), hmpid.photonsCharge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.yMip()) == hv0) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH4_HV0"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV0"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV0"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv1) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH4_HV1"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV1"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV1"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv2) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH4_HV2"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV2"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV2"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv3) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH4_HV3"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV3"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV3"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv4) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH4_HV4"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV4"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV4"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv5) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH4_HV5"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV5"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV5"), hmpid.nPhotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distanceCondition && mipChargeCondition && physicalChAngle) // condizione da verificare a priori
        {
          if (hmpid.xMip() >= xMinPc[0] && hmpid.xMip() <= xMaxPc[0] && hmpid.yMip() >= yMinPc[0] && hmpid.yMip() <= yMaxPc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC0"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[1] && hmpid.xMip() <= xMaxPc[1] && hmpid.yMip() >= yMinPc[1] && hmpid.yMip() <= yMaxPc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC1"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[2] && hmpid.xMip() <= xMaxPc[2] && hmpid.yMip() >= yMinPc[2] && hmpid.yMip() <= yMaxPc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC2"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[3] && hmpid.xMip() <= xMaxPc[3] && hmpid.yMip() >= yMinPc[3] && hmpid.yMip() <= yMaxPc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC3"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[4] && hmpid.xMip() <= xMaxPc[4] && hmpid.yMip() >= yMinPc[4] && hmpid.yMip() <= yMaxPc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC4"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[5] && hmpid.xMip() <= xMaxPc[5] && hmpid.yMip() >= yMinPc[5] && hmpid.yMip() <= yMaxPc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC5"), sin2changle, hmpid.nPhotons());
          }
        }
      }

      if (hmpid.chamber() == rich5) {
        histos.fill(HIST("hmpidXTrack5"), hmpid.xTrack());
        histos.fill(HIST("hmpidYTrack5"), hmpid.yTrack());
        histos.fill(HIST("hmpidXMip5"), hmpid.xMip());
        histos.fill(HIST("hmpidYMip5"), hmpid.yMip());
        histos.fill(HIST("hmpidXResiduals5"), hmpid.xMip() - hmpid.xTrack());
        histos.fill(HIST("hmpidYResiduals5"), hmpid.yMip() - hmpid.yTrack());

        if (hmpid.momentumTrack() > cutMinMomGlobalTrack) {
          if (hmpid.momentumHmpid() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos5"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsPos5"), hmpid.yMip() - hmpid.yTrack());
          }

          if (hmpid.momentumHmpid() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg5"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsNeg5"), hmpid.yMip() - hmpid.yTrack());
          }
        }

        if (distanceCondition) {
          histos.fill(HIST("hmpidQMip5"), hmpid.chargeMip());
        }
        histos.fill(HIST("hmpidClusSize5"), hmpid.clusterSize());
        histos.fill(HIST("TrackMom5"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom5"), std::fabs(hmpid.momentumHmpid()));
        histos.fill(HIST("hmpidXYMip5"), hmpid.xMip(), hmpid.yMip());

        if (distanceCondition && mipChargeCondition) {
          histos.fill(HIST("hmpidNPhotons5"), hmpid.nPhotons());
          sin2changle = static_cast<float>(TMath::Power(TMath::Sin(hmpid.chAngle()), 2));
          if (hmpid.xMip() <= maxBoxHit && hmpid.xMip() >= minBoxHit && hmpid.yMip() <= maxBoxHit && hmpid.yMip() >= minBoxHit) {
            histos.fill(HIST("nPhotons_vs_sin2Ch5"), sin2changle, hmpid.nPhotons());
          }

          for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
            if (hmpid.photonsCharge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge5"), hmpid.photonsCharge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.yMip()) == hv0) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH5_HV0"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV0"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV0"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv1) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH5_HV1"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV1"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV1"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv2) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH5_HV2"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV2"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV2"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv3) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH5_HV3"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV3"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV3"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv4) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH5_HV4"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV4"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV4"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv5) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH5_HV5"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV5"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV5"), hmpid.nPhotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distanceCondition && mipChargeCondition && physicalChAngle) // condizione da verificare a priori
        {
          if (hmpid.xMip() >= xMinPc[0] && hmpid.xMip() <= xMaxPc[0] && hmpid.yMip() >= yMinPc[0] && hmpid.yMip() <= yMaxPc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC0"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[1] && hmpid.xMip() <= xMaxPc[1] && hmpid.yMip() >= yMinPc[1] && hmpid.yMip() <= yMaxPc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC1"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[2] && hmpid.xMip() <= xMaxPc[2] && hmpid.yMip() >= yMinPc[2] && hmpid.yMip() <= yMaxPc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC2"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[3] && hmpid.xMip() <= xMaxPc[3] && hmpid.yMip() >= yMinPc[3] && hmpid.yMip() <= yMaxPc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC3"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[4] && hmpid.xMip() <= xMaxPc[4] && hmpid.yMip() >= yMinPc[4] && hmpid.yMip() <= yMaxPc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC4"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[5] && hmpid.xMip() <= xMaxPc[5] && hmpid.yMip() >= yMinPc[5] && hmpid.yMip() <= yMaxPc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC5"), sin2changle, hmpid.nPhotons());
          }
        }
      }

      if (hmpid.chamber() == rich6) {
        histos.fill(HIST("hmpidXTrack6"), hmpid.xTrack());
        histos.fill(HIST("hmpidYTrack6"), hmpid.yTrack());
        histos.fill(HIST("hmpidXMip6"), hmpid.xMip());
        histos.fill(HIST("hmpidYMip6"), hmpid.yMip());
        histos.fill(HIST("hmpidXResiduals6"), hmpid.xMip() - hmpid.xTrack());
        histos.fill(HIST("hmpidYResiduals6"), hmpid.yMip() - hmpid.yTrack());

        if (hmpid.momentumTrack() > cutMinMomGlobalTrack) {
          if (hmpid.momentumHmpid() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos6"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsPos6"), hmpid.yMip() - hmpid.yTrack());
          }

          if (hmpid.momentumHmpid() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg6"), hmpid.xMip() - hmpid.xTrack());
            histos.fill(HIST("hmpidYResidualsNeg6"), hmpid.yMip() - hmpid.yTrack());
          }
        }

        if (distanceCondition) {
          histos.fill(HIST("hmpidQMip6"), hmpid.chargeMip());
        }
        histos.fill(HIST("hmpidClusSize6"), hmpid.clusterSize());
        histos.fill(HIST("TrackMom6"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom6"), std::fabs(hmpid.momentumHmpid()));
        histos.fill(HIST("hmpidXYMip6"), hmpid.xMip(), hmpid.yMip());

        if (distanceCondition && mipChargeCondition) {
          histos.fill(HIST("hmpidNPhotons6"), hmpid.nPhotons());
          sin2changle = static_cast<float>(TMath::Power(TMath::Sin(hmpid.chAngle()), 2));
          if (hmpid.xMip() <= maxBoxHit && hmpid.xMip() >= minBoxHit && hmpid.yMip() <= maxBoxHit && hmpid.yMip() >= minBoxHit) {
            histos.fill(HIST("nPhotons_vs_sin2Ch6"), sin2changle, hmpid.nPhotons());
          }

          for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
            if (hmpid.photonsCharge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge6"), hmpid.photonsCharge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.yMip()) == hv0) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH6_HV0"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV0"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV0"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv1) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH6_HV1"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV1"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV1"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv2) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH6_HV2"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV2"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV2"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv3) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH6_HV3"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV3"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV3"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv4) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH6_HV4"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV4"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV4"), hmpid.nPhotons());
          }
        }

        if (fParam->inHVSector(hmpid.yMip()) == hv5) {
          if (distanceCondition) {
            histos.fill(HIST("hmpidQMip_RICH6_HV5"), hmpid.chargeMip());
          }
          if (distanceCondition && mipChargeCondition && physicalChAngle) {
            for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
              if (hmpid.photonsCharge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV5"), hmpid.photonsCharge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV5"), hmpid.nPhotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distanceCondition && mipChargeCondition && physicalChAngle) // condizione da verificare a priori
        {
          if (hmpid.xMip() >= xMinPc[0] && hmpid.xMip() <= xMaxPc[0] && hmpid.yMip() >= yMinPc[0] && hmpid.yMip() <= yMaxPc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC0"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[1] && hmpid.xMip() <= xMaxPc[1] && hmpid.yMip() >= yMinPc[1] && hmpid.yMip() <= yMaxPc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC1"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[2] && hmpid.xMip() <= xMaxPc[2] && hmpid.yMip() >= yMinPc[2] && hmpid.yMip() <= yMaxPc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC2"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[3] && hmpid.xMip() <= xMaxPc[3] && hmpid.yMip() >= yMinPc[3] && hmpid.yMip() <= yMaxPc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC3"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[4] && hmpid.xMip() <= xMaxPc[4] && hmpid.yMip() >= yMinPc[4] && hmpid.yMip() <= yMaxPc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC4"), sin2changle, hmpid.nPhotons());
          }

          if (hmpid.xMip() >= xMinPc[5] && hmpid.xMip() <= xMaxPc[5] && hmpid.yMip() >= yMinPc[5] && hmpid.yMip() <= yMaxPc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC5"), sin2changle, hmpid.nPhotons());
          }
        }
      }

      double probsHMP[3];

      getProbability(hmpid.chAngle(), std::fabs(hmpid.momentumHmpid()), probsHMP);

      if (distanceMipToTrack > maxDistanceForProb || hmpid.chargeMip() < cutQmip)
        continue;

      histos.fill(HIST("hmpidCkovvsMom_nocut"), std::fabs(hmpid.momentumHmpid()), hmpid.chAngle());

      if (probsHMP[0] < minProbParticle && probsHMP[1] < minProbParticle && probsHMP[2] < minProbParticle)
        continue;
      // if(hmpid.momentumTrack()<0.75 && hmpid.nPhotons()<7 && hmpid.chAngle()<0.52) continue;

      histos.fill(HIST("hmpidCkovvsMom"), std::fabs(hmpid.momentumHmpid()), hmpid.chAngle());

    } // close loop on tracks

  } // close process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<HmpidQa>(cfg)}; }
