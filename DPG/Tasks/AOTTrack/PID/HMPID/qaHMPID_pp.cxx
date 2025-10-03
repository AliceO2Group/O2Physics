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

// O2 includes
#include "tableHMPID.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

#include "TF1.h"
#include <TRandom.h>
#include <TTree.h>

#include "HMPIDBase/Param.h"

//////////////////////////////////////////
// task per leggere la table creata dal task HMPIDanalysis per performance e stablity --dataset PbPb and pp
/// !!! TABLE CREATA DAL TASK CHE TIENE CONTO DEL PHOTSCHARGE COME VETTORE

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// distance 2D between MIP and Track extrapolated
float distance2D(float x1, float y1, float x2, float y2)
{
  float d = TMath::Sqrt(TMath::Power(x2 - x1, 2) + TMath::Power(y2 - y1, 2));
  return d;
}

Double_t expectedSignal(double mass, float mom)
{
  // expected theoretical values

  double chAngle_th = -999.;

  const Double_t nmean = 1.288; // radiator mean refraction index

  Double_t cos_chAngle_th = (TMath::Sqrt(mass * mass + mom * mom)) / (nmean * mom);
  if (cos_chAngle_th > 1)
    return chAngle_th;

  chAngle_th = TMath::ACos(cos_chAngle_th);

  return chAngle_th;
}

Double_t expectedSigma(int iPart, float mom)
{
  Double_t sigma_ring = 0;

  // fit parameters from sigma extrapolation
  Double_t fit_sigmaExtra_pions[6] = {42.8961, -49.8723, 26.2311, -6.59093, 0.754578, -0.0286546};
  Double_t fit_sigmaExtra_kaons[6] = {76.9786, -103.655, 61.7533, -18.7436, 2.87855, -0.178318};
  Double_t fit_sigmaExtra_protons[6] = {299.466, -383.277, 201.127, -52.2554, 6.67285, -0.334106};

  // create sigma vs p functions
  TF1* sigma_vs_p_pions = new TF1("sigma_vs_p_pions", "pol5", 0., 6.);
  TF1* sigma_vs_p_kaons = new TF1("sigma_vs_p_kaons", "pol5", 0., 6.);
  TF1* sigma_vs_p_protons = new TF1("sigma_vs_p_protons", "pol5", 0., 6.);

  for (int i = 0; i < 6; i++) {
    sigma_vs_p_pions->SetParameter(i, fit_sigmaExtra_pions[i]);
    sigma_vs_p_kaons->SetParameter(i, fit_sigmaExtra_kaons[i]);
    sigma_vs_p_protons->SetParameter(i, fit_sigmaExtra_protons[i]);
  }

  if (iPart == 0) {
    sigma_ring = gRandom->Gaus(sigma_vs_p_pions->Eval(mom) / 1000., 0.1 * sigma_vs_p_pions->Eval(mom) / 1000.);
  }
  if (iPart == 1) {
    sigma_ring = 0.8 * gRandom->Gaus(sigma_vs_p_kaons->Eval(mom) / 1000., 0.15 * sigma_vs_p_kaons->Eval(mom) / 1000.);
  }
  if (iPart == 2) {
    sigma_ring = 0.6 * gRandom->Gaus(sigma_vs_p_protons->Eval(mom) / 1000., 0.1 * sigma_vs_p_protons->Eval(mom) / 1000.);
  }

  delete sigma_vs_p_pions;
  delete sigma_vs_p_kaons;
  delete sigma_vs_p_protons;

  return sigma_ring;
}

void getProbability(float hmpidSignal, float hmpidMomentum, double* probs)
{
  // Calculates probability to be a pion-kaon-proton with the "amplitude" method
  // from the given Cerenkov angle and momentum assuming no initial particle composition (class taken by AliROOT)

  Double_t probabilityParticle;

  Int_t nSpecies = 3;

  if (hmpidSignal <= 0) {
    // HMPID does not find anything reasonable for this track, assign 0.33 for all species
    for (Int_t iPart = 0; iPart < nSpecies; iPart++)
      probs[iPart] = 1.0 / (Float_t)nSpecies;
    return;
  }

  // assign mass in GeV/c^2
  Double_t mass[] = {0.139570, 0.493677, 0.938272};

  /*
  Double_t chAngle_th = expectedSignal(mass, hmpidMomentum);

  //check
  if (chAngle_th > 900.)
  {
    return 0.0; // No light emitted, probability is zero
  }

  /////////////////////////////////////////////////////////////

  Double_t sigma_ring = expectedSigma(hmpidMomentum, particleType);

  if(sigma_ring==0)
  {
    return 0.0; //no ring
  }

  /////////////////////////////////////////////////////////////

  // Check if within the "desert" range
  if (TMath::Abs(hmpidSignal - chAngle_th) >= 4 * sigma_ring)
  {
    return 1.0/((Float_t) nSpecies); // Assume equal probability in "desert" case
  }

  // Calculate Gaussian height for the selected species
  Double_t hSpecies = TMath::Gaus(chAngle_th, hmpidSignal, sigma_ring, kTRUE);

  // To normalize the probability, calculate the total height (hTot) across all species
  Double_t hTot = 0;
  for (int i = 0; i < nSpecies; i++)
  {
    Double_t chAngle_th_i = expectedSignal(mass_species[i], hmpidMomentum);
    Double_t sigma_ring_i = expectedSigma(hmpidMomentum, particles_species[i]);
    if (sigma_ring_i > 0 && chAngle_th_i <= 900)
    {
      hTot += TMath::Gaus(chAngle_th_i, hmpidSignal, sigma_ring_i, kTRUE);
    }


  }

  probabilityParticle = hSpecies / hTot;

  // Return the normalized probability for the specific species
  return (hTot > 0) ? probabilityParticle : 0.0;
*/

  Double_t hTot = 0;                    // Initialize the total height of the amplitude method
  Double_t* h = new Double_t[nSpecies]; // number of charged particles to be considered

  Bool_t desert = kTRUE; // Flag to evaluate if ThetaC is far ("desert") from the given Gaussians

  for (Int_t iPart = 0; iPart < nSpecies; iPart++) { // for each particle

    h[iPart] = 0; // reset the height
    Double_t thetaCerTh = expectedSignal(mass[iPart], hmpidMomentum);
    ; // theoretical Theta Cherenkov
    if (thetaCerTh > 900.)
      continue; // no light emitted, zero height
    Double_t sigmaRing = expectedSigma(iPart, hmpidMomentum);

    if (sigmaRing == 0)
      continue;

    if (TMath::Abs(hmpidSignal - thetaCerTh) < 4 * sigmaRing)
      desert = kFALSE;
    h[iPart] = TMath::Gaus(thetaCerTh, hmpidSignal, sigmaRing, kTRUE);
    hTot += h[iPart]; // total height of all theoretical heights for normalization

  } // species loop

  for (Int_t iPart = 0; iPart < nSpecies; iPart++) { // species loop to assign probabilities

    if (!desert)
      probs[iPart] = h[iPart] / hTot;
    else
      probs[iPart] = 1.0 / (Float_t)nSpecies; // all theoretical values are far away from experemental one
  }

  delete[] h;
}

struct pidHmpidQapp {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsP{"nBinsP", 1200, "Number of momentum bins"};
  Configurable<float> minP{"minP", -20.f, "Minimum momentum plotted (GeV/c)"};
  Configurable<float> maxP{"maxP", 20.f, "Maximum momentum plotted (GeV/c)"};

  Configurable<int> nBinsCh{"nBinsCh", 800, "Number of ch angle bins"};
  Configurable<float> minCh{"minCh", 0.f, "Minimum ch angle plotted (rad)"};
  Configurable<float> maxCh{"maxCh", 0.8f, "Maximum ch angle plotted (rad)"};

  /// valori per filtri su tracce primarie
  Configurable<float> nsigmaTPCMin{"nsigmaTPCMin", -3.0, "nsigmaTPCMin"};
  Configurable<float> nsigmaTPCMax{"nsigmaTPCMax", +3.0, "nsigmaTPCMax"};
  Configurable<float> nsigmaTOFMin{"nsigmaTOFMin", -3.0, "nsigmaTOFMin"};
  Configurable<float> nsigmaTOFMax{"nsigmaTOFMax", +3.5, "nsigmaTOFMax"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 4.0, "min number of clusters required in ITS"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 50.0f, "minTPCnClsFound"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.5f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.5f, "maxDCAz"};

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
    for (int iCh = 0; iCh < 7; iCh++) {
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
      for (int iSec = 0; iSec < 6; iSec++) {
        histos.add(Form("hmpidQMip_RICH%i_HV%i", iCh, iSec), Form("hmpidQMip_RICH%i_HV%i", iCh, iSec), kTH1F, {{2000, 200, 2200, "Charge (ADC)"}});
        histos.add(Form("hmpidNPhotons_RICH%i_HV%i", iCh, iSec), Form("hmpidNPhotons_RICH%i_HV%i", iCh, iSec), kTH1F, {{50, 2, 50, "Number of photons"}});
        histos.add(Form("hmpidPhotsCharge_RICH%i_HV%i", iCh, iSec), Form("hmpidPhotsCharge_RICH%i_HV%i", iCh, iSec), kTH1F, {{180, 4, 210}});
      }

      // plot n_ph vs sin2Ch per PC
      for (int pc = 0; pc < 6; pc++) {
        histos.add(Form("nPhotons_vs_sin2Ch_RICH%i_PC%i", iCh, pc), Form("N. of Photons vs sin^{2}(#theta_{Ch}) - chamber%i, photocathode%i", iCh, pc), kTProfile, {{20, 0.0, 0.4}});
      }
    }
  }

  void process(aod::HMPID_analysis const& hmpidtable)
  {
    // photocathods limits
    const int n_photocathodes = 6;
    static float x_max_pc[n_photocathodes];
    static float y_max_pc[n_photocathodes];
    static float x_min_pc[n_photocathodes];
    static float y_min_pc[n_photocathodes];

    for (int ipc = 0; ipc < 6; ipc++) {
      x_max_pc[ipc] = (fParam->maxPcX(ipc)) - 10.;
      y_max_pc[ipc] = (fParam->maxPcY(ipc)) - 10.;
      x_min_pc[ipc] = (fParam->minPcX(ipc)) + 10.;
      y_min_pc[ipc] = (fParam->minPcY(ipc)) + 10.;
    }

    Double_t cut_d_mip_track = 3.0; // cut for the distnace mip-track

    for (auto& hmpid : hmpidtable) // loop on tracks contained in the table
    {

      // filters on primary tracks
      if (hmpid.itsNcluster() < minReqClusterITS)
        continue;
      if (hmpid.tpcNcluster() < minTPCnClsFound)
        continue;
      if (hmpid.tpcNClsCrossedRows() < minNCrossedRowsTPC)
        continue;
      if (hmpid.tpcChi2() > maxChi2TPC)
        continue;
      if (hmpid.itsChi2() > maxChi2ITS)
        continue;
      if (TMath::Abs(hmpid.dcaxy()) > maxDCAxy)
        continue;
      if (TMath::Abs(hmpid.dcaz()) > maxDCAz)
        continue;

      // fill histograms
      histos.fill(HIST("hmpidMomvsTrackMom"), fabs(hmpid.momentumTrack()), fabs(hmpid.momentumHMPID()));
      histos.fill(HIST("TrackMom"), fabs(hmpid.momentumTrack()));
      histos.fill(HIST("hmpidMom"), fabs(hmpid.momentumHMPID()));

      histos.fill(HIST("hmpidSignal"), hmpid.chAngle());

      if (hmpid.chAngle() > 0 && distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {
        Double_t pT = (Double_t)(hmpid.momentumTrack() / TMath::CosH(hmpid.etatrack()));
        histos.fill(HIST("pTvsChAngle"), pT, hmpid.chAngle());
        if (hmpid.momentumHMPID() > 0) {
          histos.fill(HIST("pTvsChAnglePos"), pT, hmpid.chAngle());
        }
        if (hmpid.momentumHMPID() < 0) {
          histos.fill(HIST("pTvsChAngleNeg"), pT, hmpid.chAngle());
        }
      }

      float sin2changle;

      if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {

        histos.fill(HIST("hmpidNPhotons"), hmpid.nphotons());

        sin2changle = (float)TMath::Power(TMath::Sin(hmpid.chAngle()), 2);
        if (hmpid.xmip() <= 100. && hmpid.xmip() >= 40. && hmpid.ymip() <= 100. && hmpid.ymip() >= 40.) {
          histos.fill(HIST("nPhotons_vs_sin2Ch"), sin2changle, hmpid.nphotons());
        }

        for (int i = 0; i < 10; i++) {
          if (hmpid.photons_charge()[i] > 0)
            histos.fill(HIST("hmpidPhotsCharge"), hmpid.photons_charge()[i]);
        }
      }

      if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
        histos.fill(HIST("hmpidQMip"), hmpid.chargeMIP());
      }

      histos.fill(HIST("hmpidXTrack"), hmpid.xtrack());
      histos.fill(HIST("hmpidYTrack"), hmpid.ytrack());
      histos.fill(HIST("hmpidXMip"), hmpid.xmip());
      histos.fill(HIST("hmpidYMip"), hmpid.ymip());
      if (hmpid.momentumTrack() > 1.5) {
        histos.fill(HIST("hmpidXResiduals"), hmpid.xmip() - hmpid.xtrack());
        histos.fill(HIST("hmpidYResiduals"), hmpid.ymip() - hmpid.ytrack());
      }
      histos.fill(HIST("hmpidXYTrack"), hmpid.xtrack(), hmpid.ytrack());
      histos.fill(HIST("hmpidXYMip"), hmpid.xmip(), hmpid.ymip());

      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////
      // fill histograms per chamber
      if (hmpid.chamber() == 0) {
        histos.fill(HIST("hmpidXTrack0"), hmpid.xtrack());
        histos.fill(HIST("hmpidYTrack0"), hmpid.ytrack());
        histos.fill(HIST("hmpidXMip0"), hmpid.xmip());
        histos.fill(HIST("hmpidYMip0"), hmpid.ymip());
        histos.fill(HIST("hmpidXResiduals0"), hmpid.xmip() - hmpid.xtrack());
        histos.fill(HIST("hmpidYResiduals0"), hmpid.ymip() - hmpid.ytrack());

        if (hmpid.momentumTrack() > 1.5) {
          if (hmpid.momentumHMPID() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos0"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsPos0"), hmpid.ymip() - hmpid.ytrack());
          }

          if (hmpid.momentumHMPID() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg0"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsNeg0"), hmpid.ymip() - hmpid.ytrack());
          }
        }

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
          histos.fill(HIST("hmpidQMip0"), hmpid.chargeMIP());
        }
        histos.fill(HIST("hmpidClusSize0"), hmpid.clustersize());
        histos.fill(HIST("TrackMom0"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom0"), fabs(hmpid.momentumHMPID()));
        histos.fill(HIST("hmpidXYMip0"), hmpid.xmip(), hmpid.ymip());

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {
          histos.fill(HIST("hmpidNPhotons0"), hmpid.nphotons());
          sin2changle = (float)TMath::Power(TMath::Sin(hmpid.chAngle()), 2);
          if (hmpid.xmip() <= 100. && hmpid.xmip() >= 40. && hmpid.ymip() <= 100. && hmpid.ymip() >= 40.) {
            histos.fill(HIST("nPhotons_vs_sin2Ch0"), sin2changle, hmpid.nphotons());
          }

          for (int i = 0; i < 10; i++) {
            if (hmpid.photons_charge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge0"), hmpid.photons_charge()[i]);
          }
        }

        //////////////////////////////////////////////////////////////////
        // plot per HV sector
        if (fParam->inHVSector(hmpid.ymip()) == 0) {

          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH0_HV0"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV0"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV0"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 1) {

          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH0_HV1"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV1"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV1"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 2) {

          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH0_HV2"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV2"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV2"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 3) {

          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH0_HV3"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV3"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV3"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 4) {

          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH0_HV4"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV4"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV4"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 5) {

          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH0_HV5"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH0_HV5"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH0_HV5"), hmpid.nphotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) // condizione da verificare a priori
        {
          if (hmpid.xmip() >= x_min_pc[0] && hmpid.xmip() <= x_max_pc[0] && hmpid.ymip() >= y_min_pc[0] && hmpid.ymip() <= y_max_pc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC0"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[1] && hmpid.xmip() <= x_max_pc[1] && hmpid.ymip() >= y_min_pc[1] && hmpid.ymip() <= y_max_pc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC1"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[2] && hmpid.xmip() <= x_max_pc[2] && hmpid.ymip() >= y_min_pc[2] && hmpid.ymip() <= y_max_pc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC2"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[3] && hmpid.xmip() <= x_max_pc[3] && hmpid.ymip() >= y_min_pc[3] && hmpid.ymip() <= y_max_pc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC3"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[4] && hmpid.xmip() <= x_max_pc[4] && hmpid.ymip() >= y_min_pc[4] && hmpid.ymip() <= y_max_pc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC4"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[5] && hmpid.xmip() <= x_max_pc[5] && hmpid.ymip() >= y_min_pc[5] && hmpid.ymip() <= y_max_pc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH0_PC5"), sin2changle, hmpid.nphotons());
          }
        }
      }

      if (hmpid.chamber() == 1) {
        histos.fill(HIST("hmpidXTrack1"), hmpid.xtrack());
        histos.fill(HIST("hmpidYTrack1"), hmpid.ytrack());
        histos.fill(HIST("hmpidXMip1"), hmpid.xmip());
        histos.fill(HIST("hmpidYMip1"), hmpid.ymip());
        histos.fill(HIST("hmpidXResiduals1"), hmpid.xmip() - hmpid.xtrack());
        histos.fill(HIST("hmpidYResiduals1"), hmpid.ymip() - hmpid.ytrack());

        if (hmpid.momentumTrack() > 1.5) {
          if (hmpid.momentumHMPID() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos1"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsPos1"), hmpid.ymip() - hmpid.ytrack());
          }

          if (hmpid.momentumHMPID() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg1"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsNeg1"), hmpid.ymip() - hmpid.ytrack());
          }
        }

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
          histos.fill(HIST("hmpidQMip1"), hmpid.chargeMIP());
        }
        histos.fill(HIST("hmpidClusSize1"), hmpid.clustersize());
        histos.fill(HIST("TrackMom1"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom1"), fabs(hmpid.momentumHMPID()));
        histos.fill(HIST("hmpidXYMip1"), hmpid.xmip(), hmpid.ymip());

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {
          histos.fill(HIST("hmpidNPhotons1"), hmpid.nphotons());
          sin2changle = (float)TMath::Power(TMath::Sin(hmpid.chAngle()), 2);
          if (hmpid.xmip() <= 100. && hmpid.xmip() >= 40. && hmpid.ymip() <= 100. && hmpid.ymip() >= 40.) {
            histos.fill(HIST("nPhotons_vs_sin2Ch1"), sin2changle, hmpid.nphotons());
          }

          for (int i = 0; i < 10; i++) {
            if (hmpid.photons_charge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge1"), hmpid.photons_charge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.ymip()) == 0) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH1_HV0"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV0"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV0"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 1) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH1_HV1"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV1"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV1"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 2) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH1_HV2"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV2"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV2"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 3) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH1_HV3"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV3"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV3"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 4) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH1_HV4"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV4"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV4"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 5) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH1_HV5"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH1_HV5"), hmpid.photons_charge()[i]);
              } // hmpidPhotsCharge_RICH%i_HV%i
            }
            histos.fill(HIST("hmpidNPhotons_RICH1_HV5"), hmpid.nphotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) // condizione da verificare a priori
        {
          if (hmpid.xmip() >= x_min_pc[0] && hmpid.xmip() <= x_max_pc[0] && hmpid.ymip() >= y_min_pc[0] && hmpid.ymip() <= y_max_pc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC0"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[1] && hmpid.xmip() <= x_max_pc[1] && hmpid.ymip() >= y_min_pc[1] && hmpid.ymip() <= y_max_pc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC1"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[2] && hmpid.xmip() <= x_max_pc[2] && hmpid.ymip() >= y_min_pc[2] && hmpid.ymip() <= y_max_pc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC2"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[3] && hmpid.xmip() <= x_max_pc[3] && hmpid.ymip() >= y_min_pc[3] && hmpid.ymip() <= y_max_pc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC3"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[4] && hmpid.xmip() <= x_max_pc[4] && hmpid.ymip() >= y_min_pc[4] && hmpid.ymip() <= y_max_pc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC4"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[5] && hmpid.xmip() <= x_max_pc[5] && hmpid.ymip() >= y_min_pc[5] && hmpid.ymip() <= y_max_pc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH1_PC5"), sin2changle, hmpid.nphotons());
          }
        }
      }

      if (hmpid.chamber() == 2) {
        histos.fill(HIST("hmpidXTrack2"), hmpid.xtrack());
        histos.fill(HIST("hmpidYTrack2"), hmpid.ytrack());
        histos.fill(HIST("hmpidXMip2"), hmpid.xmip());
        histos.fill(HIST("hmpidYMip2"), hmpid.ymip());
        histos.fill(HIST("hmpidXResiduals2"), hmpid.xmip() - hmpid.xtrack());
        histos.fill(HIST("hmpidYResiduals2"), hmpid.ymip() - hmpid.ytrack());

        if (hmpid.momentumTrack() > 1.5) {
          if (hmpid.momentumHMPID() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos2"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsPos2"), hmpid.ymip() - hmpid.ytrack());
          }

          if (hmpid.momentumHMPID() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg2"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsNeg2"), hmpid.ymip() - hmpid.ytrack());
          }
        }

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
          histos.fill(HIST("hmpidQMip2"), hmpid.chargeMIP());
        }
        histos.fill(HIST("hmpidClusSize2"), hmpid.clustersize());
        histos.fill(HIST("TrackMom2"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom2"), fabs(hmpid.momentumHMPID()));
        histos.fill(HIST("hmpidXYMip2"), hmpid.xmip(), hmpid.ymip());

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {
          histos.fill(HIST("hmpidNPhotons2"), hmpid.nphotons());
          sin2changle = (float)TMath::Power(TMath::Sin(hmpid.chAngle()), 2);
          if (hmpid.xmip() <= 100. && hmpid.xmip() >= 40. && hmpid.ymip() <= 100. && hmpid.ymip() >= 40.) {
            histos.fill(HIST("nPhotons_vs_sin2Ch2"), sin2changle, hmpid.nphotons());
          }

          for (int i = 0; i < 10; i++) {
            if (hmpid.photons_charge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge2"), hmpid.photons_charge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.ymip()) == 0) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH2_HV0"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV0"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV0"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 1) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH2_HV1"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV1"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV1"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 2) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH2_HV2"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV2"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV2"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 3) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH2_HV3"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV3"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV3"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 4) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH2_HV4"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV4"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV4"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 5) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH2_HV5"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH2_HV5"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH2_HV5"), hmpid.nphotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) // condizione da verificare a priori
        {
          if (hmpid.xmip() >= x_min_pc[0] && hmpid.xmip() <= x_max_pc[0] && hmpid.ymip() >= y_min_pc[0] && hmpid.ymip() <= y_max_pc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC0"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[1] && hmpid.xmip() <= x_max_pc[1] && hmpid.ymip() >= y_min_pc[1] && hmpid.ymip() <= y_max_pc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC1"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[2] && hmpid.xmip() <= x_max_pc[2] && hmpid.ymip() >= y_min_pc[2] && hmpid.ymip() <= y_max_pc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC2"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[3] && hmpid.xmip() <= x_max_pc[3] && hmpid.ymip() >= y_min_pc[3] && hmpid.ymip() <= y_max_pc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC3"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[4] && hmpid.xmip() <= x_max_pc[4] && hmpid.ymip() >= y_min_pc[4] && hmpid.ymip() <= y_max_pc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC4"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[5] && hmpid.xmip() <= x_max_pc[5] && hmpid.ymip() >= y_min_pc[5] && hmpid.ymip() <= y_max_pc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH2_PC5"), sin2changle, hmpid.nphotons());
          }
        }
      }

      if (hmpid.chamber() == 3) {
        histos.fill(HIST("hmpidXTrack3"), hmpid.xtrack());
        histos.fill(HIST("hmpidYTrack3"), hmpid.ytrack());
        histos.fill(HIST("hmpidXMip3"), hmpid.xmip());
        histos.fill(HIST("hmpidYMip3"), hmpid.ymip());
        histos.fill(HIST("hmpidXResiduals3"), hmpid.xmip() - hmpid.xtrack());
        histos.fill(HIST("hmpidYResiduals3"), hmpid.ymip() - hmpid.ytrack());

        if (hmpid.momentumTrack() > 1.5) {
          if (hmpid.momentumHMPID() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos3"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsPos3"), hmpid.ymip() - hmpid.ytrack());
          }

          if (hmpid.momentumHMPID() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg3"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsNeg3"), hmpid.ymip() - hmpid.ytrack());
          }
        }

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
          histos.fill(HIST("hmpidQMip3"), hmpid.chargeMIP());
        }
        histos.fill(HIST("hmpidClusSize3"), hmpid.clustersize());
        histos.fill(HIST("TrackMom3"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom3"), fabs(hmpid.momentumHMPID()));
        histos.fill(HIST("hmpidXYMip3"), hmpid.xmip(), hmpid.ymip());

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {
          histos.fill(HIST("hmpidNPhotons3"), hmpid.nphotons());
          sin2changle = (float)TMath::Power(TMath::Sin(hmpid.chAngle()), 2);
          if (hmpid.xmip() <= 100. && hmpid.xmip() >= 40. && hmpid.ymip() <= 100. && hmpid.ymip() >= 40.) {
            histos.fill(HIST("nPhotons_vs_sin2Ch3"), sin2changle, hmpid.nphotons());
          }

          for (int i = 0; i < 10; i++) {
            if (hmpid.photons_charge()[i] > 0) {
              histos.fill(HIST("hmpidPhotsCharge3"), hmpid.photons_charge()[i]);
            }
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.ymip()) == 0) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH3_HV0"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV0"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV0"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 1) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH3_HV1"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV1"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV1"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 2) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH3_HV2"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV2"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV2"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 3) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH3_HV3"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV3"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV3"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 4) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH3_HV4"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV4"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV4"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 5) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH3_HV5"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH3_HV5"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH3_HV5"), hmpid.nphotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) // condizione da verificare a priori
        {
          if (hmpid.xmip() >= x_min_pc[0] && hmpid.xmip() <= x_max_pc[0] && hmpid.ymip() >= y_min_pc[0] && hmpid.ymip() <= y_max_pc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC0"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[1] && hmpid.xmip() <= x_max_pc[1] && hmpid.ymip() >= y_min_pc[1] && hmpid.ymip() <= y_max_pc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC1"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[2] && hmpid.xmip() <= x_max_pc[2] && hmpid.ymip() >= y_min_pc[2] && hmpid.ymip() <= y_max_pc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC2"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[3] && hmpid.xmip() <= x_max_pc[3] && hmpid.ymip() >= y_min_pc[3] && hmpid.ymip() <= y_max_pc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC3"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[4] && hmpid.xmip() <= x_max_pc[4] && hmpid.ymip() >= y_min_pc[4] && hmpid.ymip() <= y_max_pc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC4"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[5] && hmpid.xmip() <= x_max_pc[5] && hmpid.ymip() >= y_min_pc[5] && hmpid.ymip() <= y_max_pc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH3_PC5"), sin2changle, hmpid.nphotons());
          }
        }
      }

      if (hmpid.chamber() == 4) {
        histos.fill(HIST("hmpidXTrack4"), hmpid.xtrack());
        histos.fill(HIST("hmpidYTrack4"), hmpid.ytrack());
        histos.fill(HIST("hmpidXMip4"), hmpid.xmip());
        histos.fill(HIST("hmpidYMip4"), hmpid.ymip());
        histos.fill(HIST("hmpidXResiduals4"), hmpid.xmip() - hmpid.xtrack());
        histos.fill(HIST("hmpidYResiduals4"), hmpid.ymip() - hmpid.ytrack());

        if (hmpid.momentumTrack() > 1.5) {
          if (hmpid.momentumHMPID() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos4"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsPos4"), hmpid.ymip() - hmpid.ytrack());
          }

          if (hmpid.momentumHMPID() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg4"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsNeg4"), hmpid.ymip() - hmpid.ytrack());
          }
        }

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
          histos.fill(HIST("hmpidQMip4"), hmpid.chargeMIP());
        }
        histos.fill(HIST("hmpidClusSize4"), hmpid.clustersize());
        histos.fill(HIST("TrackMom4"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom4"), fabs(hmpid.momentumHMPID()));
        histos.fill(HIST("hmpidXYMip4"), hmpid.xmip(), hmpid.ymip());

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {
          histos.fill(HIST("hmpidNPhotons4"), hmpid.nphotons());
          sin2changle = (float)TMath::Power(TMath::Sin(hmpid.chAngle()), 2);
          if (hmpid.xmip() <= 100. && hmpid.xmip() >= 40. && hmpid.ymip() <= 100. && hmpid.ymip() >= 40.) {
            histos.fill(HIST("nPhotons_vs_sin2Ch4"), sin2changle, hmpid.nphotons());
          }

          for (int i = 0; i < 10; i++) {
            if (hmpid.photons_charge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge4"), hmpid.photons_charge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.ymip()) == 0) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH4_HV0"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV0"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV0"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 1) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH4_HV1"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV1"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV1"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 2) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH4_HV2"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV2"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV2"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 3) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH4_HV3"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV3"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV3"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 4) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH4_HV4"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV4"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV4"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 5) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH4_HV5"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH4_HV5"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH4_HV5"), hmpid.nphotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) // condizione da verificare a priori
        {
          if (hmpid.xmip() >= x_min_pc[0] && hmpid.xmip() <= x_max_pc[0] && hmpid.ymip() >= y_min_pc[0] && hmpid.ymip() <= y_max_pc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC0"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[1] && hmpid.xmip() <= x_max_pc[1] && hmpid.ymip() >= y_min_pc[1] && hmpid.ymip() <= y_max_pc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC1"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[2] && hmpid.xmip() <= x_max_pc[2] && hmpid.ymip() >= y_min_pc[2] && hmpid.ymip() <= y_max_pc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC2"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[3] && hmpid.xmip() <= x_max_pc[3] && hmpid.ymip() >= y_min_pc[3] && hmpid.ymip() <= y_max_pc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC3"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[4] && hmpid.xmip() <= x_max_pc[4] && hmpid.ymip() >= y_min_pc[4] && hmpid.ymip() <= y_max_pc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC4"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[5] && hmpid.xmip() <= x_max_pc[5] && hmpid.ymip() >= y_min_pc[5] && hmpid.ymip() <= y_max_pc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH4_PC5"), sin2changle, hmpid.nphotons());
          }
        }
      }

      if (hmpid.chamber() == 5) {
        histos.fill(HIST("hmpidXTrack5"), hmpid.xtrack());
        histos.fill(HIST("hmpidYTrack5"), hmpid.ytrack());
        histos.fill(HIST("hmpidXMip5"), hmpid.xmip());
        histos.fill(HIST("hmpidYMip5"), hmpid.ymip());
        histos.fill(HIST("hmpidXResiduals5"), hmpid.xmip() - hmpid.xtrack());
        histos.fill(HIST("hmpidYResiduals5"), hmpid.ymip() - hmpid.ytrack());

        if (hmpid.momentumTrack() > 1.5) {
          if (hmpid.momentumHMPID() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos5"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsPos5"), hmpid.ymip() - hmpid.ytrack());
          }

          if (hmpid.momentumHMPID() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg5"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsNeg5"), hmpid.ymip() - hmpid.ytrack());
          }
        }

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
          histos.fill(HIST("hmpidQMip5"), hmpid.chargeMIP());
        }
        histos.fill(HIST("hmpidClusSize5"), hmpid.clustersize());
        histos.fill(HIST("TrackMom5"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom5"), fabs(hmpid.momentumHMPID()));
        histos.fill(HIST("hmpidXYMip5"), hmpid.xmip(), hmpid.ymip());

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {
          histos.fill(HIST("hmpidNPhotons5"), hmpid.nphotons());
          sin2changle = (float)TMath::Power(TMath::Sin(hmpid.chAngle()), 2);
          if (hmpid.xmip() <= 100. && hmpid.xmip() >= 40. && hmpid.ymip() <= 100. && hmpid.ymip() >= 40.) {
            histos.fill(HIST("nPhotons_vs_sin2Ch5"), sin2changle, hmpid.nphotons());
          }

          for (int i = 0; i < 10; i++) {
            if (hmpid.photons_charge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge5"), hmpid.photons_charge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.ymip()) == 0) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH5_HV0"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV0"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV0"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 1) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH5_HV1"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV1"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV1"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 2) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH5_HV2"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV2"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV2"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 3) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH5_HV3"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV3"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV3"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 4) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH5_HV4"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV4"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV4"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 5) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH5_HV5"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH5_HV5"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH5_HV5"), hmpid.nphotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) // condizione da verificare a priori
        {
          if (hmpid.xmip() >= x_min_pc[0] && hmpid.xmip() <= x_max_pc[0] && hmpid.ymip() >= y_min_pc[0] && hmpid.ymip() <= y_max_pc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC0"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[1] && hmpid.xmip() <= x_max_pc[1] && hmpid.ymip() >= y_min_pc[1] && hmpid.ymip() <= y_max_pc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC1"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[2] && hmpid.xmip() <= x_max_pc[2] && hmpid.ymip() >= y_min_pc[2] && hmpid.ymip() <= y_max_pc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC2"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[3] && hmpid.xmip() <= x_max_pc[3] && hmpid.ymip() >= y_min_pc[3] && hmpid.ymip() <= y_max_pc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC3"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[4] && hmpid.xmip() <= x_max_pc[4] && hmpid.ymip() >= y_min_pc[4] && hmpid.ymip() <= y_max_pc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC4"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[5] && hmpid.xmip() <= x_max_pc[5] && hmpid.ymip() >= y_min_pc[5] && hmpid.ymip() <= y_max_pc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH5_PC5"), sin2changle, hmpid.nphotons());
          }
        }
      }

      if (hmpid.chamber() == 6) {
        histos.fill(HIST("hmpidXTrack6"), hmpid.xtrack());
        histos.fill(HIST("hmpidYTrack6"), hmpid.ytrack());
        histos.fill(HIST("hmpidXMip6"), hmpid.xmip());
        histos.fill(HIST("hmpidYMip6"), hmpid.ymip());
        histos.fill(HIST("hmpidXResiduals6"), hmpid.xmip() - hmpid.xtrack());
        histos.fill(HIST("hmpidYResiduals6"), hmpid.ymip() - hmpid.ytrack());

        if (hmpid.momentumTrack() > 1.5) {
          if (hmpid.momentumHMPID() > 0) {
            // fill residual histos for positive charges
            histos.fill(HIST("hmpidXResidualsPos6"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsPos6"), hmpid.ymip() - hmpid.ytrack());
          }

          if (hmpid.momentumHMPID() < 0) {
            // fill residual histos for negative charges
            histos.fill(HIST("hmpidXResidualsNeg6"), hmpid.xmip() - hmpid.xtrack());
            histos.fill(HIST("hmpidYResidualsNeg6"), hmpid.ymip() - hmpid.ytrack());
          }
        }

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
          histos.fill(HIST("hmpidQMip6"), hmpid.chargeMIP());
        }
        histos.fill(HIST("hmpidClusSize6"), hmpid.clustersize());
        histos.fill(HIST("TrackMom6"), hmpid.momentumTrack());
        histos.fill(HIST("hmpidMom6"), fabs(hmpid.momentumHMPID()));
        histos.fill(HIST("hmpidXYMip6"), hmpid.xmip(), hmpid.ymip());

        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120.) {
          histos.fill(HIST("hmpidNPhotons6"), hmpid.nphotons());
          sin2changle = (float)TMath::Power(TMath::Sin(hmpid.chAngle()), 2);
          if (hmpid.xmip() <= 100. && hmpid.xmip() >= 40. && hmpid.ymip() <= 100. && hmpid.ymip() >= 40.) {
            histos.fill(HIST("nPhotons_vs_sin2Ch6"), sin2changle, hmpid.nphotons());
          }

          for (int i = 0; i < 10; i++) {
            if (hmpid.photons_charge()[i] > 0)
              histos.fill(HIST("hmpidPhotsCharge6"), hmpid.photons_charge()[i]);
          }
        }

        // plot per HV sector
        if (fParam->inHVSector(hmpid.ymip()) == 0) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH6_HV0"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV0"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV0"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 1) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH6_HV1"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV1"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV1"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 2) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH6_HV2"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV2"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV2"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 3) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH6_HV3"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV3"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV3"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 4) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH6_HV4"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV4"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV4"), hmpid.nphotons());
          }
        }

        if (fParam->inHVSector(hmpid.ymip()) == 5) {
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track) {
            histos.fill(HIST("hmpidQMip_RICH6_HV5"), hmpid.chargeMIP());
          }
          if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) {
            for (int i = 0; i < 10; i++) {
              if (hmpid.photons_charge()[i] > 0) {
                histos.fill(HIST("hmpidPhotsCharge_RICH6_HV5"), hmpid.photons_charge()[i]);
              }
            }
            histos.fill(HIST("hmpidNPhotons_RICH6_HV5"), hmpid.nphotons());
          }
        }

        //////////////////////////////////////////////////////////////////
        // fill plot photocathode
        if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) < cut_d_mip_track && hmpid.chargeMIP() > 120. && hmpid.chAngle() > 0) // condizione da verificare a priori
        {
          if (hmpid.xmip() >= x_min_pc[0] && hmpid.xmip() <= x_max_pc[0] && hmpid.ymip() >= y_min_pc[0] && hmpid.ymip() <= y_max_pc[0]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC0"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[1] && hmpid.xmip() <= x_max_pc[1] && hmpid.ymip() >= y_min_pc[1] && hmpid.ymip() <= y_max_pc[1]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC1"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[2] && hmpid.xmip() <= x_max_pc[2] && hmpid.ymip() >= y_min_pc[2] && hmpid.ymip() <= y_max_pc[2]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC2"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[3] && hmpid.xmip() <= x_max_pc[3] && hmpid.ymip() >= y_min_pc[3] && hmpid.ymip() <= y_max_pc[3]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC3"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[4] && hmpid.xmip() <= x_max_pc[4] && hmpid.ymip() >= y_min_pc[4] && hmpid.ymip() <= y_max_pc[4]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC4"), sin2changle, hmpid.nphotons());
          }

          if (hmpid.xmip() >= x_min_pc[5] && hmpid.xmip() <= x_max_pc[5] && hmpid.ymip() >= y_min_pc[5] && hmpid.ymip() <= y_max_pc[5]) {
            histos.fill(HIST("nPhotons_vs_sin2Ch_RICH6_PC5"), sin2changle, hmpid.nphotons());
          }
        }
      }

      // EVALUATE PROBABILITES --- syntax: float hmpidSignal, float hmpidMomentum, const std::string& particleType
      /* Double_t prob_pion = getProbability(hmpid.chAngle(), fabs(hmpid.momentumHMPID()), "pion");
       Double_t prob_kaon = getProbability(hmpid.chAngle(), fabs(hmpid.momentumHMPID()), "kaon");
       Double_t prob_proton = getProbability(hmpid.chAngle(), fabs(hmpid.momentumHMPID()), "proton");*/

      Double_t probsHMP[3];

      getProbability(hmpid.chAngle(), fabs(hmpid.momentumHMPID()), probsHMP);

      // Printf("INIZIO CALCOLO PROBABILITIES");
      /* if(hmpid.chAngle() > 0.49 && hmpid.chAngle() < 0.54 && fabs(hmpid.momentumHMPID()) > 0.5 && fabs(hmpid.momentumHMPID()) < 0.66)
        {
         Printf("prob pion = %f", probsHMP[0]);
         Printf("prob kaon = %f", probsHMP[1]);
         Printf("prob proton = %f", probsHMP[2]);
        }
 */

      if (distance2D(hmpid.xtrack(), hmpid.ytrack(), hmpid.xmip(), hmpid.ymip()) > 1.5 || hmpid.chargeMIP() < 120.)
        continue;

      histos.fill(HIST("hmpidCkovvsMom_nocut"), fabs(hmpid.momentumHMPID()), hmpid.chAngle());

      if (probsHMP[0] < 0.7 && probsHMP[1] < 0.7 && probsHMP[2] < 0.7)
        continue;
      // if(hmpid.momentumTrack()<0.75 && hmpid.nphotons()<7 && hmpid.chAngle()<0.52) continue;

      histos.fill(HIST("hmpidCkovvsMom"), fabs(hmpid.momentumHMPID()), hmpid.chAngle());

      // pTree->Fill();
      // photCharge.clear();

    } // close loop on tracks

  } // close process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<pidHmpidQapp>(cfg)}; }
