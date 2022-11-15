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

/// \file HFCandidateCreatorOmegac.cxx
/// \brief Reconstruction of Omegac0 candidates
/// \note Extended from HFTrackIndexSkimsCreator.cxx and lambdakzerobuilder.cxx
/// \author Federica Zanone <federica.zanone@cern.ch>, HEIDELBERG UNIVERSITY & GSI

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_casc;
using namespace o2::aod::hf_cand_omegac;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Perform MC matching."}}; //to be set in .json file
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

// Reconstruction of omegac candidates
struct HfCandidateCreatorOmegac {
  Produces<aod::HfCandOmegacBase> rowCandidateBase; //produced table

  // - - - - - - - - - - CONFIGURABLES - - - - - - - - - -
  Configurable<double> magneticField{"d_bz", 5., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> b_dovalplots{"b_dovalplots", true, "do validation plots"};
  Configurable<bool> useabsdca{"useabsdca", true, "use absolute DCA for vertexing"};
  Configurable<bool> b_refitwithmatcorr{"b_refitwithmatcorr", true, "when doing propagateTracksToVertex, propagate tracks to vtx with material corrections and rerun minimization"};
  /*
  UseAbsDCA:
  TRUE -> simple average of tracks position propagated to respective X_dca parameters and rotated to the lab. frame
  FALSE -> weighted (by tracks covariances) average of tracks position propagated to respective X_dca parameters and rotated to the lab. frame
  */

  /*ABOUT DCA FITTER: To get the most precise results one can request ft.setRefitWithMatCorr(true): in this case when propagateTracksToVertex() is 
  called, the tracks will be propagated to the V0 with requested material corrections, one new V0 minimization will be done and only after that 
  the final propagation to final V0 position will be done. Since this is CPU consiming, it is reccomended to disable propagation to V0 by default
  (ft.setPropagateToPCA(false)) and call separately ft.propagateTracksToVertex() after preliminary checks on the V0 candidate  
  (see AliceO2/Detectors/Vertexing/README.md )*/

  // - - - - - - - - - - OUTPUTOBJECTS - - - - - - - - - -
  //These istograms will be added to the file AnalysisResults.root
  OutputObj<TH1F> hPtPrimaryPi{TH1F("hPtPrimaryPi", "p_T primary #pi;p_T (GeV/#it{c});entries", 500, 0, 20)};
  OutputObj<TH1F> hxVertexOmegac{TH1F("hxVertexOmegac", "x Omegac vertex;xVtx;entries", 500, -10, 10)};
  OutputObj<TH1F> hEtaCascade{TH1F("hEtaCascade", "Cascade pseudorapidity;eta;entries", 200, -1., 1.)};
  OutputObj<TH1F> hInvMassCascade{TH1F("hInvMassCascade", "Cascade invariant mass;inv mass;entries", 500, 0.9, 1.7)};
  OutputObj<TH1F> hCosPAV0{TH1F("hCosPAV0", "CosPA V0;cosPA;entries", 200, -1., 1.)};
  OutputObj<TH1F> hCTauCascade{TH1F("hCTauCascade", "Cascade ctau;ctau;entries", 500, 0., 10.)};

  OutputObj<TH1F> hMassLambdaLF{TH1F("hMassLambdaLF", "LF Lambda mass;invmass;entries", 1000, 0.8, 1.3)};

  /*Using flags in filters (for example d_selectionFlagXic): first define the filter with the flag, then specify in the tables taken as arguments in process: process(soa::Filtered<nometabella>), 
  by doing so I will automatically take from the table <nometabella> only the elements that satisfy the condition specified in the filter. Example:
  Configurable<int> d_selectionFlagXic{"d_selectionFlagXic", 1, "Selection Flag for Xic"}; --> set flag value in a configurable
  Filter filterSelectCandidates = (aod::hf_selcandidate_xic::isSelXicToPKPi >= d_selectionFlagXic || aod::hf_selcandidate_xic::isSelXicToPiKP >= d_selectionFlagXic);  --> define the filter */

  //Configurable<double> cutPtPionMin{"cutPtPionMin", 1., "min. pt pion track"};

  // - - - - - - - - - - PROCESS - - - - - - - - - -
  void process(aod::Collision const& collision,
               aod::HfCascades const& cascades,
               aod::BigTracks const& tracks,
               aod::V0Datas const&,
               aod::V0sLinked const&)
               
  {
  
    // 2-prong vertex fitter to build the omegac vertex
    o2::vertexing::DCAFitterN<2> df;
    df.setBz(magneticField);
    df.setPropagateToPCA(b_propdca);
    df.setMaxR(d_maxr);
    df.setMaxDZIni(d_maxdzini);
    df.setMinParamChange(d_minparamchange);
    df.setMinRelChi2Change(d_minrelchi2change);
    df.setUseAbsDCA(useabsdca);
    df.setRefitWithMatCorr(b_refitwithmatcorr);

    // 2-prong vertex fitter to build the cascade vertex
    o2::vertexing::DCAFitterN<2> dfc;
    dfc.setBz(magneticField);
    dfc.setPropagateToPCA(b_propdca);
    dfc.setMaxR(d_maxr);
    dfc.setMaxDZIni(d_maxdzini);
    dfc.setMinParamChange(d_minparamchange);
    dfc.setMinRelChi2Change(d_minrelchi2change);
    dfc.setUseAbsDCA(useabsdca);
    dfc.setRefitWithMatCorr(b_refitwithmatcorr);

    // 2-prong vertex fitter to build the V0 vertex
    o2::vertexing::DCAFitterN<2> dfv;
    dfv.setBz(magneticField);
    dfv.setPropagateToPCA(b_propdca);
    dfv.setMaxR(d_maxr);
    dfv.setMaxDZIni(d_maxdzini);
    dfv.setMinParamChange(d_minparamchange);
    dfv.setMinRelChi2Change(d_minrelchi2change);
    dfv.setUseAbsDCA(useabsdca);
    dfv.setRefitWithMatCorr(b_refitwithmatcorr);

    //INFO aboout DCA fitter: AliceO2/Detectors/Vertexing/include/DetectorsVertexing/DCAFitterN.h 

    //loop over cascades reconstructed by HFTrackIndexSkimsCreator.cxx (struct HfTrackIndexSkimsCreatorCascades)
    for (auto& casc : cascades) {

//--------------------------accessing particles in the decay chain------------

      //cascade daughter - charged particle
      int index_trackxidaucharged = casc.index0Id();                 //pion <- xi index from cascade table (not used)
      auto trackxidaucharged = casc.index0_as<aod::BigTracks>();     //pion <- xi track from BigTracks table
      //cascade daughter - V0
      int indexv0 = casc.v0Id();                                     //VO index from cascades table
      if (!casc.v0_as<aod::V0sLinked>().has_v0Data()){               //check that V0 data are stored 
        continue;
      }
      auto v0element = casc.v0_as<aod::V0sLinked>().v0Data();  //V0 element from LF table containing V0 info
      //V0 positive daughter
      auto trackv0dau0 = v0element.posTrack_as<aod::BigTracks>();    //p <- V0 track (positive track) from BigTracks table
      //V0 negative daughter
      auto trackv0dau1 = v0element.negTrack_as<aod::BigTracks>();    //pion <- V0 track (negative track) from BigTracks table

      /*
      LOGF(info, "- - - - - - -   Particles processed   - - - - - - -");
      LOGF(info, "Process cascade %d with V0 daughter %d and charged daughter %d ", casc.globalIndex(), casc.v0Id(), casc.index0Id());
      LOGF(info, "Process V0 (V0Id %d) daughters %d (pos) and %d (neg) ", casc.v0Id(), trackv0dau0.globalIndex(), trackv0dau1.globalIndex());
      */
//--------------------------reconstruct V0--------------------------------------

      //pseudorapidity
      double pseudorap_v0posdau = trackv0dau0.eta();
      double pseudorap_v0negdau = trackv0dau1.eta();

      //pion & p <- V0 track to be processed with DCAfitter
      auto trackParCovV0Dau0 = getTrackParCov(trackv0dau0);
      auto trackParCovV0Dau1 = getTrackParCov(trackv0dau1);

      //info from LF table
      array<float, 3> pvecV0LFtable = {v0element.px(), v0element.py(), v0element.pz()};
      array<float, 3> vertexV0 = {v0element.x(), v0element.y(), v0element.z()};
      //LOGF(info, "Momentum V0 from LF table: %f (px) %f (py) %f (pz)", v0element.px(), v0element.py(), v0element.pz());
      //LOGF(info, "Vertex V0 from LF table: %f (x) %f (y) %f (z)", v0element.x(), v0element.y(), v0element.z());

      //reconstruct V0 with DCAFitter
      int nfv = dfv.process(trackParCovV0Dau0, trackParCovV0Dau1); //process restituisce un int che corrisponde al numero di vertici trovati
      //LOGF(info, "- - - - - - -   V0 DCA fitter   - - - - - - -");
      //LOGF(info, "Process value (V0): %d", nfv);
      if (nfv == 0) {
        continue;
      }
      auto vertexV0FromFitter = dfv.getPCACandidate(); //After finding the vertex tracks are propagated to the PCA if setPropagateToPCA is set to true
      auto chi2PCA_v0 = dfv.getChi2AtPCACandidate();
      //NOT NEEDED, PropagateTracksToVertex() IS CALLED
      //trackParCovV0Dau0.propagateTo(vertexV0FromFitter[0],magneticField);
      //trackParCovV0Dau1.propagateTo(vertexV0FromFitter[0],magneticField);
      array<float, 3> pvecV0Dau0;
      array<float, 3> pvecV0Dau1;
      dfv.propagateTracksToVertex();
      dfv.getTrack(0).getPxPyPzGlo(pvecV0Dau0); // update V0daughter0 momentum at V0 decay vertex (track number follows the order of tracks used in the function df.process)
      dfv.getTrack(1).getPxPyPzGlo(pvecV0Dau1); // update V0daughter1 momentum at V0 decay vertex
      array<float, 3> pvecV0_m = {pvecV0Dau0[0]+pvecV0Dau1[0], pvecV0Dau0[1]+pvecV0Dau1[1], pvecV0Dau0[2]+pvecV0Dau1[2]};
    
      std::array<float, 3> coordvtx_v0 = dfv.getPCACandidatePos();
      std::array<float, 3> momvtx_v0 = pvecV0_m;
      std::array<float, 6> covvtx_v0 = dfv.calcPCACovMatrixFlat();
      //LOGF(info,"Decay vertex V0 from DCAFitter: %f (x), %f (y), %f (z)", coordvtx_v0[0], coordvtx_v0[1], coordvtx_v0[2]);
      //LOGF(info,"Momentum V0 from DCAFitter: %f (px), %f (py), %f (pz)", momvtx_v0[0], momvtx_v0[1], momvtx_v0[2]);
      //LOGF(info,"Covariance (decay vertex) V0 from DCAFitter: %f (x), %f (y), %f (z), %f, %f, %f", covvtx_v0[0], covvtx_v0[1], covvtx_v0[2], covvtx_v0[3], covvtx_v0[4], covvtx_v0[5]);

      //create V0 track
      auto trackV0 = o2::dataformats::V0(dfv.getPCACandidatePos(), pvecV0_m, dfv.calcPCACovMatrixFlat(), trackParCovV0Dau0, trackParCovV0Dau1, {0, 0}, {0, 0});
      //see  AliceO2/DataFormats/Reconstruction/include/ReconstructionDataFormats/V0.h 

//-----------------------------reconstruct cascade------------------------------

      //pseudorapidity
      double pseudorap_pifromcas = trackxidaucharged.eta();

      //pion <- casc track to be processed with DCAfitter
      auto trackParVar_xidaucharged = getTrackParCov(trackxidaucharged);

      //reconstruct cascade with DCAFitter
      int nfc = dfc.process(trackV0, trackParVar_xidaucharged);
      //LOGF(info, "- - - - - - -   Cascade DCA fitter   - - - - - - -");
      //LOGF(info, "Process value (cascade): %d", nfc);
      if (nfc == 0) {
        continue;
      }
      auto vertexcascFromFitter = dfc.getPCACandidate(); //After finding the vertex tracks are propagated to the DCA
      auto chi2PCA_cascade = dfc.getChi2AtPCACandidate();
      //NOT NEEDED, PropagateTracksToVertex() IS CALLED
      //trackV0.propagateTo(vertexcascFromFitter[0],magneticField);
      //trackParVar_xidaucharged.propagateTo(vertexcascFromFitter[0],magneticField);
      array<float, 3> pvecV0_d;
      array<float, 3> pvecpionfromcasc;
      dfc.propagateTracksToVertex();
      dfc.getTrack(0).getPxPyPzGlo(pvecV0_d); // update V0 momentum at cascade decay vertex
      dfc.getTrack(1).getPxPyPzGlo(pvecpionfromcasc); // update pi <- cascade momentum at cascade decay vertex
      array<float, 3> pveccasc_m = {pvecV0_d[0]+pvecpionfromcasc[0], pvecV0_d[1]+pvecpionfromcasc[1], pvecV0_d[2]+pvecpionfromcasc[2]};

      std::array<float, 3> coordvtx_casc = dfc.getPCACandidatePos();
      std::array<float, 3> momvtx_casc = pveccasc_m;
      std::array<float, 6> covvtx_casc = dfc.calcPCACovMatrixFlat();
      //LOGF(info,"Decay vertex cascade from DCAFitter: %f (x), %f (y), %f (z)", coordvtx_casc[0], coordvtx_casc[1], coordvtx_casc[2]);
      //LOGF(info,"Momentum cascade from DCAFitter: %f (px), %f (py), %f (pz)", momvtx_casc[0], momvtx_casc[1], momvtx_casc[2]);
      //LOGF(info,"Covariance (decay vertex) cascade from DCAFitter: %f (x), %f (y), %f (z), %f, %f, %f", covvtx_casc[0], covvtx_casc[1], covvtx_casc[2], covvtx_casc[3], covvtx_casc[4], covvtx_casc[5]);

      //create cascade track
      auto trackcasc = o2::dataformats::V0(dfc.getPCACandidatePos(), pveccasc_m, dfc.calcPCACovMatrixFlat(), trackV0, trackParVar_xidaucharged, {0, 0}, {0, 0});
      auto trackcasc_copy = trackcasc;

//-------------------combining cascade and pion tracks--------------------------

      for (auto& trackpion : tracks) {

        /*cut on primary pion pT
        if (trackpion.pt() < cutPtPionMin) {  //primary pion pT cut
          continue;
        }*/

        //ask for opposite sign daughters (omegac daughters)
        if (trackpion.sign() * trackxidaucharged.sign() >= 0) {
          continue;
        }

        //check not to take the same particle twice in the decay chain
        if (trackpion.globalIndex() == trackxidaucharged.globalIndex() || trackpion.globalIndex() == trackv0dau0.globalIndex() || trackpion.globalIndex() == trackv0dau1.globalIndex() || trackpion.globalIndex() == casc.globalIndex()) {
          continue;
        }

        //pseudirapidity
        double pseudorap_pifromome = trackpion.eta();

        //primary pion track to be processed with DCAFitter
        auto trackParVarPi = getTrackParCov(trackpion);
        auto trackParVarPi_copy = trackParVarPi;

        //reconstruct omegac with DCAFitter
        int nfo = df.process(trackcasc, trackParVarPi);
        //LOGF(info, "- - - - - - -   Omegac DCA fitter   - - - - - - -");
        //LOGF(info, "Process value (omegac): %d", nfo);
        // reconstruct the 2-prong omegac
        if (nfo == 0) {
          continue;
        }
        auto vertexomegacFromFitter = df.getPCACandidate();
        auto chi2PCA_omegac = df.getChi2AtPCACandidate();
        auto covMatrixPCA = df.calcPCACovMatrixFlat();
        //NOT NEEDED, PropagateTracksToVertex() IS CALLED
        //trackcasc.propagateTo(vertexomegacFromFitter[0],magneticField);
        //trackParVarPi.propagateTo(vertexomegacFromFitter[0],magneticField);
        array<float, 3> pveccasc_d;
        array<float, 3> pvecpionfromomegac;
        df.propagateTracksToVertex();
        df.getTrack(0).getPxPyPzGlo(pveccasc_d); // update cascade momentum at omegac decay vertex
        df.getTrack(1).getPxPyPzGlo(pvecpionfromomegac); // update pi <- omegac momentum at omegac decay vertex
        array<float, 3> pvecomegac = {pveccasc_d[0]+pvecpionfromomegac[0], pveccasc_d[1]+pvecpionfromomegac[1], pveccasc_d[2]+pvecpionfromomegac[2]};  //omegac momentum
        
        std::array<float, 3> coordvtx_omegac = df.getPCACandidatePos();
        std::array<float, 3> momvtx_omegac = pvecomegac;
        std::array<float, 6> covvtx_omegac = df.calcPCACovMatrixFlat();
        //LOGF(info,"Decay vertex omegac from DCAFitter: %f (x), %f (y), %f (z)", coordvtx_omegac[0], coordvtx_omegac[1], coordvtx_omegac[2]);
        //LOGF(info,"Momentum omegac from DCAFitter: %f (px), %f (py), %f (pz)", momvtx_omegac[0], momvtx_omegac[1], momvtx_omegac[2]);
        //LOGF(info,"Covariance (decay vertex) omegac from DCAFitter: %f (x), %f (y), %f (z), %f, %f, %f", covvtx_omegac[0], covvtx_omegac[1], covvtx_omegac[2], covvtx_omegac[3], covvtx_omegac[4], covvtx_omegac[5]);

        //create omegac track
        auto trackomegac = o2::dataformats::V0(df.getPCACandidatePos(), pvecomegac, df.calcPCACovMatrixFlat(), trackcasc, trackParVarPi, {0, 0}, {0, 0});

        //impact parameter omegac
        o2::dataformats::DCA impactParameterOmegac;
        auto primaryVertex = getPrimaryVertex(collision);
        trackomegac.propagateToDCA(primaryVertex, magneticField, &impactParameterOmegac);

        //impact parameter omegac daughters
        auto covMatrixPV = primaryVertex.getCov();
        o2::dataformats::DCA impactParameterCasc; //impact parameter cascade
        o2::dataformats::DCA impactParameterPrimaryPi; //impact parameter pion
        trackcasc_copy.propagateToDCA(primaryVertex, magneticField, &impactParameterCasc);
        trackParVarPi_copy.propagateToDCA(primaryVertex, magneticField, &impactParameterPrimaryPi);

        //get uncertainty of the decay length
        double phi, theta;
        getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, vertexomegacFromFitter, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        //set hfFlag
        int hfFlag = 1 << DecayType::OmegacToXiPi;

        hxVertexOmegac->Fill(vertexomegacFromFitter[0]);

        //primary pi  pT spectrum
        double ptprimarypi = sqrt((pvecpionfromomegac[0]*pvecpionfromomegac[0])+(pvecpionfromomegac[1]*pvecpionfromomegac[1]));
        hPtPrimaryPi->Fill(ptprimarypi);

        //computing invariant mass under the hypothesis of particles ID corresponding to the decay chain
        double m_lambda = v0element.mLambda(); //from LF table, V0 mass under lambda hypothesis
        double m_antilambda = v0element.mAntiLambda(); //from LF table, V0 mass under anti-lambda hypothesis

        const std::array<float, 2> arrMas_cascade = {1.11568, 0.13957};
        double m_cascade = RecoDecay::m(array{pvecV0_d, pvecpionfromcasc}, arrMas_cascade);
        /*double m_cascade = 0.;
        if(m2_cascade > 0){
          m_cascade = sqrt(m2_cascade);
        }*/

        const std::array<float, 2> arrMas_omegac = {1.32171, 0.13957};
        double m_omegac = RecoDecay::m(array{pveccasc_d, pvecpionfromomegac}, arrMas_omegac);
        /*double m_omegac = 0.;
        if(m2_omegac > 0){
          m_omegac = sqrt(m2_omegac);
        }*/

        //computing cosPA for neutral particles
        double cpa_V0 = RecoDecay::cpa(coordvtx_casc, coordvtx_v0, pvecV0_m);
        double cpa_omegac = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, coordvtx_omegac, pvecomegac);
        double cpaxy_V0 = RecoDecay::cpaXY(coordvtx_casc, coordvtx_v0, pvecV0_m);
        double cpaxy_omegac = RecoDecay::cpaXY(array{collision.posX(), collision.posY(), collision.posZ()}, coordvtx_omegac, pvecomegac);

        //computing decay length and ctau
        double declen_omegac = RecoDecay::distance(array{collision.posX(), collision.posY(), collision.posZ()}, coordvtx_omegac);
        double declen_cascade = RecoDecay::distance(coordvtx_omegac, coordvtx_casc);
        double declen_V0 = RecoDecay::distance(coordvtx_casc, coordvtx_v0);
        double ct_omegac = RecoDecay::ct(pvecomegac, declen_omegac, 2.6952);
        double ct_cascade = RecoDecay::ct(pveccasc_d, declen_cascade, 1.32171);
        double ct_V0 = RecoDecay::ct(pvecV0_m, declen_V0, 1.11568);

        //computing eta
        double pseudorap_omegac = RecoDecay::eta(pvecomegac);
        double pseudorap_cascade = RecoDecay::eta(pveccasc_d);
        double pseudorap_v0 = RecoDecay::eta(pvecV0_m);

        //fill test histograms
        hEtaCascade->Fill(pseudorap_cascade);
        hInvMassCascade->Fill(m_cascade);
        hCosPAV0->Fill(cpa_V0);
        hCTauCascade->Fill(ct_cascade);
        if(trackxidaucharged.sign()>0){
            hMassLambdaLF->Fill(m_antilambda);
        } else if(trackxidaucharged.sign()<0){
            hMassLambdaLF->Fill(m_lambda);
        }
        

        //fill the table
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         vertexomegacFromFitter[0], vertexomegacFromFitter[1], vertexomegacFromFitter[2],
                         vertexcascFromFitter[0], vertexcascFromFitter[1], vertexcascFromFitter[2],
                         vertexV0FromFitter[0], vertexV0FromFitter[1], vertexV0FromFitter[2],                         
                         //errorDecayLength, errorDecayLengthXY,
                         trackxidaucharged.sign(),
                         chi2PCA_omegac, chi2PCA_v0, chi2PCA_cascade,
                         pvecomegac[0], pvecomegac[1], pvecomegac[2],
                         pveccasc_d[0], pveccasc_d[1], pveccasc_d[2],
                         pvecpionfromomegac[0], pvecpionfromomegac[1], pvecpionfromomegac[2],
                         pvecV0_d[0], pvecV0_d[1], pvecV0_d[2],
                         pvecpionfromcasc[0], pvecpionfromcasc[1], pvecpionfromcasc[2],
                         pvecV0Dau0[0], pvecV0Dau0[1], pvecV0Dau0[2],
                         pvecV0Dau1[0], pvecV0Dau1[1], pvecV0Dau1[2],
                         impactParameterCasc.getY(), impactParameterPrimaryPi.getY(),  //questo dovrebbe essere nel piano (x,y), non solo nella direzione y, vero? si
                         impactParameterCasc.getZ(), impactParameterPrimaryPi.getZ(),
                         std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterPrimaryPi.getSigmaY2()),
                         v0element.globalIndex(),
                         v0element.posTrackId(), v0element.negTrackId(),
                         casc.globalIndex(), 
                         trackpion.globalIndex(), //index pi <- omegac
                         trackxidaucharged.globalIndex(), //index pi <- cascade
                         impactParameterOmegac.getY(), impactParameterOmegac.getZ(),
                         hfFlag, ptprimarypi,
                         m_lambda, m_antilambda, m_cascade, m_omegac,
                         cpa_V0, cpa_omegac, cpaxy_V0, cpaxy_omegac,
                         ct_omegac, ct_cascade, ct_V0,
                         pseudorap_v0posdau, pseudorap_v0negdau, pseudorap_pifromcas, pseudorap_pifromome,
                         pseudorap_omegac, pseudorap_cascade, pseudorap_v0
                         ); 

      } // loop over pions
    }   // loop over candidates
  }     // end of process
};      //end of struct

/*
//Extends the base table with expression columns
struct HfCandidateCreatorOmegacExpressions {
  Spawns<aod::HfCandOmegacExt> rowCandidateOmegac;
  void init(InitContext const&) {}
};*/
/*Spawns<ExtendedTable> handle; User-defined extended tables need to be requested by adding
Spawns<> template to the task. The table is created before the task is run and is
accessible through the handle variable that acts as a pointer.
Looking at HFSecondaryVertex.h seems like you need to define the extended table to make the momenta 
components available to compute other quantities that have to be put in the table.
*/


/// Performs MC matching.
struct HfCandidateCreatorOmegacMc {
  Produces<aod::HfCandOmegacMCRec> rowMCMatchRec;
  Produces<aod::HfCandOmegacMCGen> rowMCMatchGen;

  void process(aod::HfCandOmegacBase const& candidates,
               aod::BigTracksMC const& tracks,
               aod::McParticles const& particlesMC) //this table contains basic MC particle properties
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    //int8_t origin = 0; //to be used for prompt/non prompt
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      //Printf("New rec. candidate");
      flag = 0;
      //origin = 0;
      debug = 0;
      //here I use the getters I implemented for the elements of the table HfCandOmegacBase (the one defined above for data analysis -> for the getters see SecondaryVertex.h)
      auto arrayDaughters = array{candidate.index1_as<aod::BigTracksMC>(), //pi <- omegac (sintax taken from CandidateCreatorLambdaB)
                                  candidate.index2_as<aod::BigTracksMC>(), //pi <- cascade
                                  candidate.posTrack_as<aod::BigTracksMC>(), //p <- lambda
                                  candidate.negTrack_as<aod::BigTracksMC>()}; //pi <- lambda
      auto arrayDaughtersCasc = array{candidate.index2_as<aod::BigTracksMC>(),
                                      candidate.posTrack_as<aod::BigTracksMC>(),
                                      candidate.negTrack_as<aod::BigTracksMC>()};
      auto arrayDaughtersV0 = array{candidate.posTrack_as<aod::BigTracksMC>(),
                                    candidate.negTrack_as<aod::BigTracksMC>()};

      // Omegac → p pi pi pi
      //Printf("Checking Omegac → p π π π");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, 4332, array{211, -211, 2212, -211}, true, &sign, 3); //pdg::Code::kOmegac0=4332 - proton 2212 - pi+ 211
      if (indexRec == -1){
        debug=1;
      }
      if (indexRec > -1) {
        // cascade → lambda pi
        //Printf("Checking cascade → lambda pi");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersCasc, 3312, array{-211, 2212, -211}, true, &sign, 2); //pdg::Code::kXiMinus=3312
        if (indexRec == -1){
          debug=2;
        }
        if(indexRec > -1){
          // v0 → p pi
          //Printf("Checking v0 → p pi");
          indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersV0, 3122, array{2212, -211}, true, &sign, 1); //pdg::Code::kLambda=3122
          if (indexRec == -1){
            debug=3;
          }
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::OmegacToXiPi);
        }
       }
      }
      if(debug == 2 || debug == 3){
        LOGF(info, "WARNING: Omegac0 decays in the expected final state but the condition on the intermediate states are not fulfilled");
      }
      rowMCMatchRec(flag, debug);
      //debug=0 -> tutti i decadimenti sono ok
      //debug=1 -> decadimento omegac in stato finale non corretto
      //debug=2 -> decadimento omegac in stato finale corretto ma decadimento cascade in stato finale non corretto
      //debug=3 -> decadimento omegac in stato finale corretto, decadimento cascade in stato finale corretto ma decadimento lambda in stato finale non corretto

    } //close loop over candidates

    // Match generated particles.
    for (auto& particle : particlesMC) {
      //Printf("New gen. candidate");
      flag = 0;
      //origin = 0;
      // Omegac → Xi pi
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, 4332, array{3312, 211}, true)) {
        // Match Xi -> lambda pi
        auto cascMC = particlesMC.rawIteratorAt(particle.daughtersIds().front());
        //Printf("Checking cascade → lambda pi");
        if (RecoDecay::isMatchedMCGen(particlesMC, cascMC, 3312, array{3122, -211}, true)) {
          //lambda -> p pi
          auto v0MC = particlesMC.rawIteratorAt(cascMC.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(particlesMC, v0MC, 3122, array{2212, -211}, true, &sign)){
            flag = sign * (1 << DecayType::OmegacToXiPi);
          }
        }
      }
      //rowMCMatchGen(flag, origin);
      rowMCMatchGen(flag);
    }

  } //close process
}; //close struct



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorOmegac>(cfgc)};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfCandidateCreatorOmegacMc>(cfgc));
  }
  return workflow;
}

/*
NOTE:
- Acceptin/rejecting collision/track --> done in HFTrackIndexSkimsCreator (with a filter within the cascades struct)

- description table BigTracks:
  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullEl, aod::pidTOFFullEl,
  aod::pidTPCFullMu, aod::pidTOFFullMu, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr,
  aod::pidTOFFullPr, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels>>;
  (all the candidatecreators in HF use the table BigTracks to make the combinatorics with a primary pion)

- for function computing momentum, cos(PA), ctau... see O2Physics/Common/Core/RecoDecay.h & trackUtilities.h
*/

//TO CHECK:
//- getPCACandidatePos restituisce coordinate cilindriche o cartesiane?
//- quale momento usare per calcolare ctau della cascade (e' carica!)