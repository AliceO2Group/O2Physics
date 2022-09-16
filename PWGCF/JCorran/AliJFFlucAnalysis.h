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

#ifndef AliJFFlucAnalysis_cxx
#define AliJFFlucAnalysis_cxx

#include "AliJHistManager.h"
#include <Math/Vector4D.h>
#include <Math/LorentzVector.h>
#include <TComplex.h>
#include <iterator>
#include <memory>

using namespace ROOT;
using namespace ROOT::Math;

#ifndef JFLUC_USE_INPUT_LIST
template <class ValueType>
class TrackIterInterface
{
 public:
  virtual ValueType& deref() = 0;
  virtual void increment() = 0;
  virtual bool equals(const TrackIterInterface<ValueType>*) = 0;
};

template <class ValueType>
class TrackIterBase : public std::iterator<std::input_iterator_tag, ValueType>
{
 public:
  // TrackIterBase(TrackIterInterface<ValueType> *_pm) : pm(_pm){};
  TrackIterBase(std::unique_ptr<TrackIterInterface<ValueType>> _pm) : pm(std::move(_pm)){};
  ~TrackIterBase(){};
  // TrackIterInterface<ValueType> *pm;
  std::unique_ptr<TrackIterInterface<ValueType>> pm;
  TrackIterBase& operator++()
  {
    pm->increment();
    return *this;
  }
  ValueType& operator*()
  {
    return pm->deref();
  }
  bool operator==(const TrackIterBase<ValueType>& m) const
  {
    return pm->equals(m.pm.get());
  }
  bool operator!=(const TrackIterBase<ValueType>& m) const
  {
    return !(*this == m);
  }
};

typedef TrackIterBase<PtEtaPhiEVector> JTrackIter;
typedef TrackIterInterface<PtEtaPhiEVector> JTrackIterInterface;
class TracksBase
{
 public:
  TracksBase() {}
  virtual ~TracksBase(){};
  virtual JTrackIter begin() = 0;
  virtual JTrackIter end() = 0;
  virtual size_t size() const = 0;
};
#endif

class AliJFFlucAnalysis
{
 public:
  AliJFFlucAnalysis();
  AliJFFlucAnalysis(const char* name);
  AliJFFlucAnalysis(const AliJFFlucAnalysis& a);             // not implemented
  AliJFFlucAnalysis& operator=(const AliJFFlucAnalysis& ap); // not implemented

  ~AliJFFlucAnalysis();
  void UserCreateOutputObjects();
  void Init();
  void UserExec(Option_t* option);
  void Terminate(Option_t*);

#ifndef JFLUC_USE_INPUT_LIST
  // this generic class can be used by standalone toyMC/hydro/etc.
  void SetInputList(TracksBase* _fInputList) { fInputList = _fInputList; }
#else
  void SetInputList(std::vector<PtEtaPhiEVector>* _fInputList)
  {
    fInputList = _fInputList;
  }
#endif
  void SetEventCentralityAndBin(float cent, UInt_t cbin)
  {
    fCent = cent;
    fCBin = cbin;
  }
  float GetEventCentrality() const { return fCent; }
  void SetEventImpactParameter(float ip) { fImpactParameter = ip; }
  void SetEventVertex(const double* vtx) { fVertex = vtx; }

  void SetEtaRange(double eta_min, double eta_max)
  {
    fEta_min = eta_min;
    fEta_max = eta_max;
  }
  void SetEventTracksQA(unsigned int tpc, unsigned int glb)
  {
    fTPCtrks = (float)tpc;
    fGlbtrks = (float)glb;
  }
  void SetEventFB32TracksQA(unsigned int fb32, unsigned int fb32tof)
  {
    fFB32trks = (float)fb32;
    fFB32TOFtrks = (float)fb32tof;
  }

  void Fill_QA_plot(double eta1, double eta2);

  // new function for QC method //
  void CalculateQvectorsQC(double, double);
  TComplex Q(int n, int p);
  TComplex Two(int n1, int n2);
  TComplex Four(int n1, int n2, int n3, int n4);

  enum SUBEVENT {
    SUBEVENT_A = 0x1,
    SUBEVENT_B = 0x2
  };
  void SelectSubevents(UInt_t _subeventMask)
  {
    subeventMask = _subeventMask;
  }
  // set the number of bins before initialization (UserCreateOutputObjects)
  void SetNumBins(UInt_t _numBins)
  {
    numBins = _numBins;
  }
  enum {
    FLUC_PHI_CORRECTION = 0x2,
    FLUC_EBE_WEIGHTING = 0x4
  };
  void AddFlags(UInt_t _flags)
  {
    flags |= _flags;
  }

  static Double_t pttJacek[74];
  static UInt_t NpttJacek;

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
         nKL };  // order
#define kcNH kH6 // max second dimension + 1
 private:
#ifndef JFLUC_USE_INPUT_LIST
  TracksBase* fInputList; // no-copy interface to tracks (for O2)
#else
  std::vector<PtEtaPhiEVector>* fInputList; // simple list for simulation codes
#endif
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

  double fEta_min;
  double fEta_max;

  TComplex QvectorQC[kNH][nKL];
  TComplex QvectorQCeta10[2][kNH][nKL]; // ksub

  AliJHistManager* fHMG; //!

  AliJBin fBin_Subset;  //!
  AliJBin fBin_h;       //!
  AliJBin fBin_k;       //!
  AliJBin fBin_hh;      //!
  AliJBin fBin_kk;      //!
  AliJBin fHistCentBin; //!
  AliJBin fVertexBin;   //! // x, y, z
  AliJBin fCorrBin;     //!

  AliJTH1D fh_cent;            //! // for cent dist
  AliJTH1D fh_ImpactParameter; //! // for impact parameter for mc
  AliJTH1D fh_vertex;          //!
  AliJTH1D fh_pt;              //! // for pt dist of tracks
  AliJTH1D fh_eta;             //! // for eta dist of tracks
  AliJTH1D fh_phi;             //! // for phi dist [ic][isub]
  AliJTH2D fh_phieta;          //!
  AliJTH3D fh_phietaz;         //!

  AliJTH1D fh_psi_n;       //!
  AliJTH1D fh_cos_n_phi;   //!
  AliJTH1D fh_sin_n_phi;   //!
  AliJTH1D fh_cos_n_psi_n; //!
  AliJTH1D fh_sin_n_psi_n; //!

  AliJTH1D fh_ntracks; //! // for number of tracks dist
  AliJTH1D fh_vn;      //!  // single vn^k  array [ih][ik][iCent]
  AliJTH1D fh_vna;     //! // single vn^k with autocorrelation removed (up to a limited order)
  AliJTH1D fh_vn_vn;   //! // combination for <vn*vn> [ih][ik][ihh][ikk][iCent]
  /*AliJTH1D fh_cn_4c;//!  // QC
  AliJTH1D fh_cn_2c;//!  // QC
  AliJTH1D fh_cn_cn_2c;//! // QC
  AliJTH1D fh_cn_2c_eta10;//!  // QC
  AliJTH1D fh_cn_cn_2c_eta10;//! // QC*/

  AliJTH1D fh_correlator;            //! // some more complex correlators
  AliJTH2D fh_TrkQA_TPCvsGlob;       //! // QA histos
  AliJTH2D fh_TrkQA_TPCvsCent;       //! // QA histos
  AliJTH2D fh_TrkQA_FB32_vs_FB32TOF; //!

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
  double NSubTracks_pt[2][N_ptbins];
  AliJBin fBin_Nptbins;       //!
  AliJTH1D fh_SC_ptdep_4corr; //! // for < vn^2 vm^2 >
  AliJTH1D fh_SC_ptdep_2corr; //!  // for < vn^2 >
  // additinal variables for SC with QC
  AliJTH1D fh_SC_with_QC_4corr;       //! // for <vn^2 vm^2>
  AliJTH1D fh_SC_with_QC_2corr;       //! // for <vn^2>
  AliJTH1D fh_SC_with_QC_2corr_eta10; //!
  AliJTH1D fh_evt_SP_QC_ratio_2p;     //! // check SP QC evt by evt ratio
  AliJTH1D fh_evt_SP_QC_ratio_4p;     //! // check SP QC evt by evt ratio
};

#endif
