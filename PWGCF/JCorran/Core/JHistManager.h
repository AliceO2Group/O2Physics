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
/// \since Sep 2022

#ifndef UTILS_ALIJARRAY_H
#define UTILS_ALIJARRAY_H

#include <vector>
#include <TString.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TPRegexp.h>
#include <TVector.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TObjString.h>
#include <TClass.h>
#include "Framework/Logger.h"

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define JERROR(...) LOGF(error, "JERROR: " __FILE__ ":" STR(__LINE__) ": " __VA_ARGS__)

class JArrayBase;
class JBin;
class JArrayAlgorithm;
class JArrayAlgorithmSimple;
class JTH1;
class JHistManager;
template <typename t>
class JTH1Derived;
template <typename t>
class JTH1DerivedPlayer;

//////////////////////////////////////////////////////
//  Utils
//////////////////////////////////////////////////////
std::vector<TString> Tokenize(TString s, TString d, int quote = 1);
TString Join(std::vector<TString>& ss, TString del = " ");
bool OutOf(int, int, int);

typedef std::vector<int> ArrayInt;
typedef std::vector<double> ArrayDouble;
// typedef JArray<void*> ArrayVoid*;

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// JNamed                                                              //
//                                                                      //
//
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class JNamed
{
 public:
  JNamed(TString name, TString title, TString opt, int mode);
  virtual ~JNamed();
  TString GetName() { return fName; }
  TString GetTitle() { return fTitle; }
  TString GetOption() { return fOption; }
  TString GetOption(TString key);
  bool HasOption(TString key) { return GetOption(key) != UndefinedOption(); }
  bool HasOption(TString key, TString val) { return GetOption(key) == val; } // TODO sensitive?
  void SetName(const char* s) { fName = s; }
  void SetTitle(const char* s) { fTitle = s; }
  void SetNameTitle(TString n, TString t)
  {
    SetName(n);
    SetTitle(t);
  }
  void SetFullOption(const char* s) { fOption = s; }
  void SetOption(TString key, TString value = "");
  void RemoveOption(TString key); // TODO
  // void SetOptionWithString( TString s );// TODO
  static TString UndefinedOption();

 protected:
  TString fName;
  TString fTitle;
  TString fOption;
  int fMode;
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// JBin                                                              //
//                                                                      //
//
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class JBin : public JNamed
{
 public:
  enum { kSingle,
         kRange,
         kString,
         kNMode };
  JBin();
  JBin(TString config, JHistManager* hmg);
  JBin(const JBin& obj);
  JBin& operator=(const JBin& obj);
  JBin& Set(TString name, TString iname, TString Title, int mode = kRange);
  void AddToManager(JHistManager* hmg);
  JBin& SetBin(const int n, const float* v);
  JBin& SetBin(const int n, const double* v);
  JBin& SetBin(TVector* v);
  JBin& SetBin(const TString v);
  JBin& SetBin(const int n);

  int GetBin(double x);

  double GetMin() { return fBinD[0]; }
  double GetMax() { return fBinD[RawSize() - 1]; }

  TString BuildTitle(int i);
  int RawSize() { return fBinD.size(); }
  int Size() { return fMode == kRange ? fBinD.size() - 1 : fBinD.size(); }
  double At(int i) { return fBinD[i]; }
  TString GetIndexName() { return fIndexName; }

  void Print();
  TString GetString();

  operator int() { return Size(); }

  static TString GetModeString(int i);
  static int GetMode(TString mode);

 private:
  void AddBin(const TString& v);
  void AddBin(float v);
  virtual void FixBin();

  std::vector<double> fBinD;
  std::vector<TString> fBinStr;
  bool fIsFixedBin;
  TString fIndexName;
  JHistManager* fHMG;
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// JArrayBase                                                        //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//________________________________________________________________________
class JArrayBase : public JNamed
{
 public:
  enum { kNormal,
         kSingle };
  virtual ~JArrayBase();
  JArrayBase& operator=(const JArrayBase& obj);

  int AddDim(int i)
  {
    fDim.push_back(i);
    return fDim.size();
  }
  int Dimension() { return fDim.size(); }
  int GetEntries() { return fArraySize; }
  int SizeOf(int i) { return fDim.at(i); }

  ArrayInt& Index() { return fIndex; }
  int Index(int d);
  void SetIndex(int i, int d);
  void ClearIndex()
  {
    fIndex.clear();
    fIndex.resize(Dimension(), 0);
  }

  void* GetItem();
  void* GetSingleItem();

  /// void LockBin(bool is=true){}//TODO
  // bool IsBinLocked(){ return fIsBinLocked; }

  virtual void FixBin();
  bool IsBinFixed() { return fIsBinFixed; }

  bool OutOfDim(int d) { return OutOf(d, 0, Dimension() - 1); }
  bool OutOfSize(int i, int d) { return OutOfDim(d) || OutOf(i, 0, SizeOf(d) - 1); }

  // Virtual
  virtual void* BuildItem() = 0;
  virtual TString BuildName() = 0;
  virtual TString BuildTitle() = 0;
  virtual void Print() = 0;
  virtual TString GetString() = 0;

  // int Resize( int size, int dim=-1 ); // NextStep
  void InitIterator();
  bool Next(void*& item);

 protected:
  JArrayBase(); // Prevent direct creation of JArrayBase
  JArrayBase(const JArrayBase& obj);

  ArrayInt fDim;   // Comment test
  ArrayInt fIndex; /// Comment test
  int fArraySize;  /// Comment test3
  int fNGenerated;
  bool fIsBinFixed;
  bool fIsBinLocked;
  JArrayAlgorithm* fAlg;
  friend class JArrayAlgorithm;
};

//________________________________________________________________________
class JArrayAlgorithm
{
 public:
  JArrayAlgorithm(JArrayBase* cmd); // TODO Move to private
  JArrayAlgorithm(const JArrayAlgorithm& obj);
  JArrayAlgorithm& operator=(const JArrayAlgorithm& obj);
  virtual ~JArrayAlgorithm();
  int Dimension() { return fCMD->Dimension(); }
  int SizeOf(int i) { return fCMD->SizeOf(i); }
  int GetEntries() { return fCMD->GetEntries(); }
  int Index(int i) { return fCMD->Index(i); }
  virtual int BuildArray() = 0;
  virtual void* GetItem() = 0;
  virtual void SetItem(void* item) = 0;
  virtual void InitIterator() = 0;
  virtual bool Next(void*& item) = 0;
  virtual void** GetRawItem() = 0;
  virtual void* GetPosition() = 0;
  virtual bool IsCurrentPosition(void* pos) = 0;
  virtual void SetPosition(void* pos) = 0;
  virtual void DeletePosition(void* pos) = 0;

 protected:
  JArrayBase* fCMD;
};

//________________________________________________________________________
class JArrayAlgorithmSimple : public JArrayAlgorithm
{
 public:
  JArrayAlgorithmSimple(JArrayBase* cmd);
  JArrayAlgorithmSimple(const JArrayAlgorithmSimple& obj);
  JArrayAlgorithmSimple& operator=(const JArrayAlgorithmSimple& obj);
  virtual ~JArrayAlgorithmSimple();
  virtual int BuildArray();
  int GlobalIndex();
  void ReverseIndex(int iG);
  virtual void* GetItem();
  virtual void SetItem(void* item);
  virtual void InitIterator() { fPos = 0; }
  virtual void** GetRawItem() { return &fArray[GlobalIndex()]; }
  virtual bool Next(void*& item)
  {
    item = fPos < GetEntries() ? (void*)fArray[fPos] : NULL;
    if (fPos < GetEntries())
      ReverseIndex(fPos);
    return fPos++ < GetEntries();
  }
  virtual void* GetPosition() { return static_cast<void*>(new int(fPos)); }
  virtual bool IsCurrentPosition(void* pos) { return *static_cast<int*>(pos) == fPos; }
  virtual void SetPosition(void* pos)
  {
    fPos = *static_cast<int*>(pos);
    ReverseIndex(fPos);
  }
  virtual void DeletePosition(void* pos) { delete static_cast<int*>(pos); }

 private:
  ArrayInt fDimFactor;
  void** fArray;
  int fPos;
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// JTH1                                                              //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//________________________________________________________________________
class JTH1 : public JArrayBase
{
 public:
  JTH1();
  JTH1(TString config, JHistManager* hmg);
  JTH1(const JTH1& obj);
  JTH1& operator=(const JTH1& obj);
  virtual ~JTH1();

  int AddDim(int i)
  {
    fDim.push_back(i);
    return fDim.size();
  }
  int AddDim(JBin* bin);
  int AddDim(TString v);
  void AddToManager(JHistManager* hmg);
  JBin* GetBinPtr(int i) { return fBins.at(i); }

  // Virtual from JArrayBase
  virtual void* BuildItem();
  virtual TString GetString();
  virtual void Print();
  virtual void FixBin();
  // Virtual from this
  virtual Int_t Write();
  // virtual Int_t WriteAll();
  virtual const char* ClassName() { return "JTH1"; }

  // Not Virtual
  virtual TString BuildName();
  virtual TString BuildTitle();
  bool IsLoadMode();
  void SetTemplate(TH1* h);
  TH1* GetTemplatePtr() { return fTemplate; }

 protected:
  TDirectory* fDirectory;
  TDirectory* fSubDirectory;
  JHistManager* fHMG;
  TH1* fTemplate;
  std::vector<JBin*> fBins;
};
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// JTH1Derived                                                       //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
template <typename T>
class JTH1Derived : public JTH1
{
 protected:
 public:
  JTH1Derived();
  JTH1Derived(TString config, JHistManager* hmg) : JTH1(config, hmg), fPlayer(this) {}
  virtual ~JTH1Derived();

  JTH1DerivedPlayer<T>& operator[](int i)
  {
    fPlayer.Init();
    fPlayer[i];
    return fPlayer;
  }
  T* operator->() { return static_cast<T*>(GetSingleItem()); }
  operator T*() { return static_cast<T*>(GetSingleItem()); }
  // Virtual from JArrayBase

  // Virtual from JTH1
  virtual const char* ClassName() { return Form("J%s", T::Class()->GetName()); }

  JTH1Derived<T>& operator<<(int i)
  {
    AddDim(i);
    return *this;
  }
  JTH1Derived<T>& operator<<(JBin& v)
  {
    AddDim(&v);
    return *this;
  }
  JTH1Derived<T>& operator<<(TString v)
  {
    AddDim(v);
    return *this;
  }
  JTH1Derived<T>& operator<<(T v)
  {
    SetTemplate(&v);
    return *this;
  }
  void SetWith(JTH1Derived<T>& v, TString name, TString title = "")
  {
    SetTemplate(v.GetTemplatePtr());
    fName = name;
    fTitle = title;
    GetTemplatePtr()->SetName(name);
    GetTemplatePtr()->SetTitle(title);
    for (int i = 0; i < v.Dimension(); i++)
      AddDim(v.GetBinPtr(i));
    AddDim("END");
  }
  void SetWith(JTH1Derived<T>& v, T tem)
  {
    SetTemplate(&tem);
    for (int i = 0; i < v.Dimension(); i++)
      AddDim(v.GetBinPtr(i));
    AddDim("END");
  }

 protected:
  JTH1DerivedPlayer<T> fPlayer;
};

//////////////////////////////////////////////////////////////////////////
// JTH1DerivedPlayer                                                 //
//////////////////////////////////////////////////////////////////////////
template <typename T>
class JTH1DerivedPlayer
{
 public:
  JTH1DerivedPlayer(JTH1Derived<T>* cmd) : fLevel(0), fCMD(cmd){};
  JTH1DerivedPlayer<T>& operator[](int i)
  {
    if (fLevel > fCMD->Dimension()) {
      JERROR("Exceed Dimension");
    }
    if (OutOf(i, 0, fCMD->SizeOf(fLevel) - 1)) {
      JERROR("wrong Index %d of %dth in %s", i, fLevel, fCMD->GetName().Data());
    }
    fCMD->SetIndex(i, fLevel++);
    return *this;
  }
  void Init()
  {
    fLevel = 0;
    fCMD->ClearIndex();
  }
  T* operator->() { return static_cast<T*>(fCMD->GetItem()); }
  operator T*() { return static_cast<T*>(fCMD->GetItem()); }
  operator TObject*() { return static_cast<TObject*>(fCMD->GetItem()); }
  operator TH1*() { return static_cast<TH1*>(fCMD->GetItem()); }

 private:
  int fLevel;
  JTH1Derived<T>* fCMD;
};

typedef JTH1Derived<TH1D> JTH1D;
typedef JTH1Derived<TH2D> JTH2D;
typedef JTH1Derived<TH3D> JTH3D;
typedef JTH1Derived<TProfile> JTProfile;

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// JHistManager                                                       //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class JHistManager : public JNamed
{
 public:
  JHistManager(TString name, TString dirname = "");
  JHistManager(const JHistManager& obj);
  JHistManager& operator=(const JHistManager& obj);
  void Add(JBin* o);
  void Add(JTH1* o);

  int GetNBin() { return fBin.size() > fBinNames.size() ? fBin.size() : fBinNames.size(); }      // TODO
  int GetNHist() { return fHist.size() > fHistNames.size() ? fHist.size() : fHistNames.size(); } // TODO
  void Print();
  int LoadConfig();
  TDirectory* GetDirectory() { return fDirectory; }
  void SetDirectory(TDirectory* d) { fDirectory = d; }
  static JHistManager* GlobalManager();
  static JHistManager* CurrentManager(JHistManager* hmg = NULL);
  JHistManager* cd() { return JHistManager::CurrentManager(this); }
  void SetLoadMode(bool b = true) { fIsLoadMode = b; }
  bool IsLoadMode() { return fIsLoadMode; }
  TString GetString()
  {
    TString st;
    for (int i = 0; i < GetNBin(); i++)
      st += fBin[i]->GetString() + "\n";
    for (int i = 0; i < GetNHist(); i++) {
      st += fHist[i]->GetString() + "\n";
    }
    return st;
  }
  void Write();
  void WriteConfig();

  JBin* GetBin(TString name);
  JBin* GetBuiltBin(TString name);
  JTH1* GetTH1(TString name);
  JTH1* GetBuiltTH1(TString name);
  JTProfile& GetTProfile(TString name) { return dynamic_cast<JTProfile&>(*GetTH1(name)); }
  JTH1D& GetTH1D(TString name) { return dynamic_cast<JTH1D&>(*GetTH1(name)); }
  JTH2D& GetTH2D(TString name) { return dynamic_cast<JTH2D&>(*GetTH1(name)); }
  JTH3D& GetTH3D(TString name) { return dynamic_cast<JTH3D&>(*GetTH1(name)); }
  bool fIsLoadMode;

  TString GetHistName(int i) { return fHistNames[i]; }

  JTH1* GetJTH1(int i) { return GetTH1(fHistNames[i]); }
  int GetNJTH1() { return fHistNames.size(); }
  bool HistogramExists(TString name);

 private:
  TDirectory* fDirectory;
  TString fConfigStr;
  std::vector<JBin*> fBin;
  std::vector<JTH1*> fHist;
  std::vector<JHistManager*> fManager;
  std::vector<TString> fBinNames;
  std::vector<TString> fBinConfigs;
  std::vector<TString> fHistNames;
  std::vector<TString> fHistConfigs;
};

#endif
