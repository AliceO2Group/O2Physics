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
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \since Sep 2022

#ifndef PWGCF_JCORRAN_CORE_JFFLUCANALYSIS_H_
#define PWGCF_JCORRAN_CORE_JFFLUCANALYSIS_H_

#include <experimental/type_traits>
#include "JHistManager.h"
#include "JQVectors.h"
#include <TComplex.h>
#include <tuple>

class JFFlucAnalysis
{
 public:
  JFFlucAnalysis();
  explicit JFFlucAnalysis(const char* name);
  explicit JFFlucAnalysis(const JFFlucAnalysis& a);    // not implemented
  JFFlucAnalysis& operator=(const JFFlucAnalysis& ap); // not implemented

  ~JFFlucAnalysis();
  void UserCreateOutputObjects();
  void Init();
  TComplex Q(int n, int p);
  TComplex Two(int n1, int n2);
  TComplex Four(int n1, int n2, int n3, int n4);
  void UserExec(Option_t* option);
  void Terminate(Option_t*);

  inline void SetEventCentralityAndBin(float cent, UInt_t cbin)
  {
    fCent = cent;
    fCBin = cbin;
  }
  inline float GetEventCentrality() const { return fCent; }
  inline void SetEventImpactParameter(float ip) { fImpactParameter = ip; }
  inline void SetEventVertex(const Double_t* vtx) { fVertex = vtx; }
  inline void SetEtaRange(Double_t eta_min, Double_t eta_max)
  {
    fEta_min = eta_min;
    fEta_max = eta_max;
  }
  inline void SetEventTracksQA(unsigned int tpc, unsigned int glb)
  {
    fTPCtrks = static_cast<float>(tpc);
    fGlbtrks = static_cast<float>(glb);
  }
  inline void SetEventFB32TracksQA(unsigned int fb32, unsigned int fb32tof)
  {
    fFB32trks = static_cast<float>(fb32);
    fFB32TOFtrks = static_cast<float>(fb32tof);
  }
  enum SubEvent {
    kSubEvent_A = 0x1,
    kSubEvent_B = 0x2
  };
  inline void SelectSubevents(UInt_t _subeventMask)
  {
    subeventMask = _subeventMask;
  }
  // set the number of bins before initialization (UserCreateOutputObjects)
  inline void SetNumBins(UInt_t _numBins)
  {
    numBins = _numBins;
  }
  enum {
    kFlucPhiCorrection = 0x2,
    kFlucEbEWeighting = 0x4
  };
  inline void AddFlags(UInt_t _flags)
  {
    flags |= _flags;
  }

  enum { kH0,
         kH1,
         kH2,
         kH3,
         kH4,
         kH5,
         kH6,
         kH7,
         kH8,
         kH9,
         kH10,
         kH11,
         kH12,
         kNH }; // harmonics
  enum { kK0,
         kK1,
         kK2,
         kK3,
         kK4,
         nKL }; // order
  using JQVectorsT = JQVectors<TComplex, kNH, nKL, true>;
  inline void SetJQVectors(const JQVectorsT* _pqvecs) { pqvecs = _pqvecs; }

  template <class T>
  using hasWeightNUA = decltype(std::declval<T&>().weightNUA());
  template <class T>
  using hasWeightEff = decltype(std::declval<T&>().weightEff());

  template <class JInputClass>
  inline void FillQA(JInputClass& inputInst)
  {
    fh_ntracks[fCBin]->Fill(inputInst.size());
    fh_ImpactParameter->Fill(fImpactParameter);
    fh_cent->Fill(fCent);

    fh_TrkQA_TPCvsCent->Fill(fCent, fTPCtrks);
    fh_TrkQA_TPCvsGlob->Fill(fGlbtrks, fTPCtrks);
    fh_TrkQA_FB32_vs_FB32TOF->Fill(fFB32trks, fFB32TOFtrks);

    for (auto& track : inputInst) {
      if (!(flags & kFlucPhiCorrection)) {
        fh_phieta[fCBin]->Fill(track.phi(), track.eta());
        fh_phietaz[fCBin]->Fill(track.phi(), track.eta(), fVertex[2]);
      }

      if (TMath::Abs(track.eta()) < fEta_min || TMath::Abs(track.eta()) > fEta_max)
        continue;

      Double_t corrInv = 1.0;
      using JInputClassIter = typename JInputClass::iterator;
      if constexpr (std::experimental::is_detected<hasWeightEff, const JInputClassIter>::value)
        corrInv /= track.weightEff();
      fh_eta[fCBin]->Fill(track.eta(), corrInv);
      fh_pt[fCBin]->Fill(track.pt(), corrInv);
      if constexpr (std::experimental::is_detected<hasWeightNUA, const JInputClassIter>::value)
        corrInv /= track.weightNUA();
      fh_phi[fCBin][(UInt_t)(track.eta() > 0.0)]->Fill(track.phi(), corrInv);
    }

    for (UInt_t iaxis = 0; iaxis < 3; iaxis++)
      fh_vertex[iaxis]->Fill(fVertex[iaxis]);
  }

  static Double_t pttJacek[74];
  static UInt_t NpttJacek;

#define kcNH kH6 // max second dimension + 1
 private:
  const Double_t* fVertex; //!
  Float_t fCent;
  Float_t fImpactParameter;
  UInt_t fCBin;
  float fTPCtrks;
  float fGlbtrks;
  float fFB32trks;
  float fFB32TOFtrks;
  UInt_t subeventMask;
  UInt_t numBins; // total number of bins
  UInt_t flags;

  Double_t fEta_min;
  Double_t fEta_max;

  const JQVectorsT* pqvecs;

  JHistManager* fHMG; //!

  JBin fBin_Subset;  //!
  JBin fBin_h;       //!
  JBin fBin_k;       //!
  JBin fBin_hh;      //!
  JBin fBin_kk;      //!
  JBin fHistCentBin; //!
  JBin fVertexBin;   //! // x, y, z
  JBin fCorrBin;     //!

  JTH1D fh_cent;            //! // for cent dist
  JTH1D fh_ImpactParameter; //! // for impact parameter for mc
  JTH1D fh_vertex;          //!
  JTH1D fh_pt;              //! // for pt dist of tracks
  JTH1D fh_eta;             //! // for eta dist of tracks
  JTH1D fh_phi;             //! // for phi dist [ic][isub]
  JTH2D fh_phieta;          //!
  JTH3D fh_phietaz;         //!

  JTH1D fh_psi_n;       //!
  JTH1D fh_cos_n_phi;   //!
  JTH1D fh_sin_n_phi;   //!
  JTH1D fh_cos_n_psi_n; //!
  JTH1D fh_sin_n_psi_n; //!

  JTH1D fh_ntracks; //! // for number of tracks dist
  JTH1D fh_vn;      //!  // single vn^k  array [ih][ik][iCent]
  JTH1D fh_vna;     //! // single vn^k with autocorrelation removed (up to a limited order)
  JTH1D fh_vn_vn;   //! // combination for <vn*vn> [ih][ik][ihh][ikk][iCent]

  JTH1D fh_correlator;            //! // some more complex correlators
  JTH2D fh_TrkQA_TPCvsGlob;       //! // QA histos
  JTH2D fh_TrkQA_TPCvsCent;       //! // QA histos
  JTH2D fh_TrkQA_FB32_vs_FB32TOF; //!

  // additional variables for ptbins(Standard Candles only)
  enum { kPt0,
         kPt1,
         kPt2,
         kPt3,
         kPt4,
         kPt5,
         kPt6,
         kPt7,
         N_ptbins };
  JBin fBin_Nptbins;             //!
  JTH1D fh_SC_ptdep_4corr;       //! // for < vn^2 vm^2 >
  JTH1D fh_SC_ptdep_2corr;       //!  // for < vn^2 >
  JTH1D fh_SC_with_QC_4corr;     //! // for <vn^2 vm^2>
  JTH1D fh_SC_with_QC_2corr;     //! // for <vn^2>
  JTH1D fh_SC_with_QC_2corr_gap; //!
  // JTH1D fh_evt_SP_QC_ratio_2p;     //! // check SP QC evt by evt ratio
  // JTH1D fh_evt_SP_QC_ratio_4p;     //! // check SP QC evt by evt ratio
};

#endif // PWGCF_JCORRAN_CORE_JFFLUCANALYSIS_H_
