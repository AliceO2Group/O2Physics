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

#include "JHistManager.h"
#include <TMath.h>
using namespace std;
//////////////////////////////////////////////////////
//  JBin
//////////////////////////////////////////////////////

JNamed::JNamed(TString name, TString title, TString opt, int mode) : fName(name),
                                                                     fTitle(title),
                                                                     fOption(opt),
                                                                     fMode(mode)
{
  // constructor
}

JNamed::~JNamed()
{
  // virtual destructor for base class
}

TString JNamed::GetOption(TString key)
{
  TPMERegexp a("&" + key + "=?([^&]*)", "i");
  int nMatch = a.Match(fOption);
  if (nMatch < 2)
    return UndefinedOption();
  return a[1];
}
void JNamed::SetOption(TString key, TString value)
{
  TPMERegexp a("&" + key + "=?[^&]*", "i");
  int nMatch = a.Match(fOption);
  TString newOpt = "&" + key + (value.Length() ? "=" + value : "");
  if (value == UndefinedOption())
    newOpt = "";
  if (nMatch < 1)
    fOption += newOpt;
  else
    fOption.ReplaceAll(a[0], newOpt);
}
void JNamed::RemoveOption(TString key)
{
  SetOption(key, UndefinedOption());
}
TString JNamed::UndefinedOption()
{
  // static TString undefinedOption = "Undefined";
  // return undefinedOption;
  return "Undefined";
}

//////////////////////////////////////////////////////
//  JBin
//////////////////////////////////////////////////////

//_____________________________________________________
JBin::JBin() : JNamed("JBin", "%.2f-%2.f", "&Mode=Range", kRange),
               fBinD(0),
               fBinStr(0),
               fIsFixedBin(false),
               fIndexName("H"),
               fHMG(NULL)
{
  ;
}
//_____________________________________________________
JBin::JBin(TString config, JHistManager* hmg) : JNamed("JBin", "%.2f-%2.f", "&Mode=Range", kRange),
                                                fBinD(0),
                                                fBinStr(0),
                                                fIsFixedBin(false),
                                                fIndexName("H"),
                                                fHMG(NULL)
{
  // cout<< config<<endl;
  LOGF(info, "JBin: %s", config.Data());
  std::vector<TString> t = Tokenize(config, " \t,");
  TString type = t[0];
  SetName(t[1]);
  fIndexName = t[2];
  SetTitle(t[3]);
  fTitle.ReplaceAll("\"", "");
  SetFullOption(t[4]);
  fMode = GetMode(GetOption("mode"));
  AddToManager(hmg);
  TString s;
  for (UInt_t i = 5; i < t.size(); i++)
    s += " " + t[i];
  SetBin(s);
}

//_____________________________________________________
JBin::JBin(const JBin& obj) : JNamed(obj.fName, obj.fTitle, obj.fOption, obj.fMode),
                              fBinD(obj.fBinD),
                              fBinStr(obj.fBinStr),
                              fIsFixedBin(obj.fIsFixedBin),
                              fIndexName(obj.fIndexName),
                              fHMG(obj.fHMG)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
JBin& JBin::operator=(const JBin& obj)
{
  // assignment operator
  if (this != &obj) {
    // TODO: proper implementation
  }
  return *this;
}

//_____________________________________________________
void JBin::FixBin()
{
  if (fIsFixedBin)
    return;
  fIsFixedBin = true;
  if (!fHMG)
    AddToManager(JHistManager::CurrentManager());
}
//_____________________________________________________
void JBin::AddToManager(JHistManager* hmg)
{
  hmg->Add(this);
}
//_____________________________________________________
JBin& JBin::Set(TString name, TString iname, TString title, int mode)
{
  SetNameTitle(name, title);
  fIndexName = iname;
  fMode = mode;
  SetOption("mode", GetModeString(mode));
  return *this;
}
//_____________________________________________________
TString JBin::GetModeString(int i)
{
  static TString mode[] = {"Single", "Range", "String"};
  if (i < 0 || i > 2)
    return "";
  return mode[i];
}
int JBin::GetMode(TString mode)
{
  for (int i = 0; i < kNMode; i++)
    if (mode == GetModeString(i))
      return i;
  return -1;
}
//_____________________________________________________
JBin& JBin::SetBin(const int n, const float* v)
{
  for (int i = 0; i < n; i++) {
    AddBin(v[i]);
  }
  FixBin();
  return *this;
}
//_____________________________________________________
JBin& JBin::SetBin(const int n, const double* v)
{
  for (int i = 0; i < n; i++) {
    AddBin(v[i]);
  }
  FixBin();
  return *this;
}
JBin& JBin::SetBin(TVector* v)
{
  for (int i = 0; i < v->GetNrows(); i++) {
    AddBin((v->GetMatrixArray())[i]);
  }
  FixBin();
  return *this;
}
//_____________________________________________________
JBin& JBin::SetBin(const TString v)
{
  std::vector<TString> ar = Tokenize(v, "\t ,");
  for (UInt_t i = 0; i < ar.size(); i++) {
    AddBin(ar[i]);
  }
  FixBin();
  return *this;
}

//_____________________________________________________
JBin& JBin::SetBin(const int n)
{
  for (UInt_t i = 0; i < UInt_t(n); i++) {
    AddBin(i);
  }
  FixBin();
  return *this;
}
//_____________________________________________________
void JBin::AddBin(const TString& v)
{
  if (fIsFixedBin) {
    JERROR("You can't Add Bini %s", GetName().Data());
  }
  fBinStr.push_back((v == "_") ? "" : v);
  fBinD.push_back(v.Atof());
}
//_____________________________________________________
void JBin::AddBin(float v)
{
  if (fIsFixedBin) {
    JERROR("You can't Add Bin %s", GetName().Data());
  }
  fBinD.push_back(v);
  fBinStr.push_back(Form("%f", v));
}
//_____________________________________________________
TString JBin::BuildTitle(int i)
{
  if (i < 0 || i > Size())
    return "";
  if (fMode == kSingle)
    return TString(Form(fTitle.Data(), fBinD[i]));
  if (fMode == kRange)
    return TString(Form(fTitle.Data(), fBinD[i], fBinD[i + 1]));
  if (fMode == kString)
    return TString(Form(fTitle.Data(), fBinStr[i].Data()));
  JERROR("Bad Mode of JBin type %d in %s", fMode, fName.Data());
  return "";
}
//_____________________________________________________
TString JBin::GetString()
{
  SetOption("mode", GetModeString(fMode));
  return "JBin\t" + fName + "\t" + fIndexName + "\t\"" + fTitle + "\"" + "\t" + fOption + "\t" + Join(fBinStr, " ");
}
//_____________________________________________________
void JBin::Print()
{
  // std::cout<<"*"+GetString()<<std::endl;
  LOGF(info, "*%s", GetString().Data());
}

int JBin::GetBin(double x)
{
  auto i = TMath::BinarySearch(fBinD.size(), &fBinD[0], x);
  if (fMode == kRange && i + 1 >= fBinD.size())
    return -1;
  return i;
}

//////////////////////////////////////////////////////
// JArrayBase
//////////////////////////////////////////////////////

//_____________________________________________________
JArrayBase::JArrayBase() : JNamed("JArayBase", "", "&Dir=default&Lazy", 0),
                           // JNamed("JArayBase","","&Dir=default&LessLazy",0),
                           fDim(0),
                           fIndex(0),
                           fArraySize(0),
                           fNGenerated(0),
                           fIsBinFixed(false),
                           fIsBinLocked(false),
                           fAlg(NULL)
{
  // constrctor
}
//_____________________________________________________
JArrayBase::~JArrayBase()
{
  // destructor
  if (fAlg)
    delete fAlg;
}

//_____________________________________________________
JArrayBase::JArrayBase(const JArrayBase& obj) : JNamed(obj.fName, obj.fTitle, obj.fOption, obj.fMode),
                                                fDim(obj.fDim),
                                                fIndex(obj.fIndex),
                                                fArraySize(obj.fArraySize),
                                                fNGenerated(obj.fNGenerated),
                                                fIsBinFixed(obj.fIsBinFixed),
                                                fIsBinLocked(obj.fIsBinLocked),
                                                fAlg(obj.fAlg)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
JArrayBase& JArrayBase::operator=(const JArrayBase& obj)
{
  // assignment operator
  if (this != &obj) {
    // TODO: proper implementation
  }
  return *this;
}

//_____________________________________________________
void* JArrayBase::GetItem()
{
  void* item = fAlg->GetItem();
  if (!item) {
    BuildItem();
    item = fAlg->GetItem();
  }
  return item;
}
//_____________________________________________________
void* JArrayBase::GetSingleItem()
{
  if (fMode == kSingle)
    return GetItem();
  JERROR("This is not single array");
  return NULL;
}
//_____________________________________________________
void JArrayBase::FixBin()
{
  if (Dimension() == 0) {
    AddDim(1);
    SetOption("Single");
    fMode = kSingle;
    if (HasOption("dir", "default"))
      RemoveOption("dir");
  }
  ClearIndex();
  fAlg = new JArrayAlgorithmSimple(this);
  fArraySize = fAlg->BuildArray();
}
//_____________________________________________________
int JArrayBase::Index(int d)
{
  if (OutOfDim(d))
    JERROR("Wrong Dim");
  return fIndex[d];
}
void JArrayBase::SetIndex(int i, int d)
{
  if (OutOfSize(i, d))
    JERROR("Wrong Index");
  fIndex[d] = i;
}

void JArrayBase::InitIterator() { fAlg->InitIterator(); }
bool JArrayBase::Next(void*& item) { return fAlg->Next(item); }

//////////////////////////////////////////////////////
//  JArrayAlgorithm
//////////////////////////////////////////////////////

//_____________________________________________________
JArrayAlgorithm::JArrayAlgorithm(JArrayBase* cmd) : fCMD(cmd)
{
  // constructor
}
//_____________________________________________________
JArrayAlgorithm::~JArrayAlgorithm()
{
  // destructor
}

//_____________________________________________________
JArrayAlgorithm::JArrayAlgorithm(const JArrayAlgorithm& obj) : fCMD(obj.fCMD)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
JArrayAlgorithm& JArrayAlgorithm::operator=(const JArrayAlgorithm& obj)
{
  // assignment operator
  if (this != &obj) {
    *fCMD = *(obj.fCMD);
  }
  return *this;
}

//////////////////////////////////////////////////////
//  JArrayAlgorithmSimple
//////////////////////////////////////////////////////

//_____________________________________________________
JArrayAlgorithmSimple::JArrayAlgorithmSimple(JArrayBase* cmd) : JArrayAlgorithm(cmd),
                                                                fDimFactor(0),
                                                                fArray(NULL),
                                                                fPos(0)
{
  // constructor
}
//_____________________________________________________
JArrayAlgorithmSimple::~JArrayAlgorithmSimple()
{
  // Dimension, GetEntries, SizeOf
  if (fArray)
    delete[] fArray;
}

//_____________________________________________________
JArrayAlgorithmSimple::JArrayAlgorithmSimple(const JArrayAlgorithmSimple& obj) : JArrayAlgorithm(obj.fCMD),
                                                                                 fDimFactor(obj.fDimFactor),
                                                                                 fArray(obj.fArray),
                                                                                 fPos(obj.fPos)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
JArrayAlgorithmSimple& JArrayAlgorithmSimple::operator=(const JArrayAlgorithmSimple& obj)
{
  // assignment operator TODO: proper implementation
  if (this != &obj) {
    *fCMD = *(obj.fCMD);
  }
  return *this;
}
//_____________________________________________________
int JArrayAlgorithmSimple::BuildArray()
{
  fDimFactor.resize(Dimension(), 1);
  for (int i = Dimension() - 2; i >= 0; i--) {
    fDimFactor[i] = fDimFactor[i + 1] * SizeOf(i + 1);
  } // TODO split to BuildArray and lazyArray in GetItem
  int arraySize = fDimFactor[0] * SizeOf(0);
  fArray = new void*[arraySize];
  for (int i = 0; i < arraySize; i++)
    fArray[i] = NULL;
  return arraySize;
}
//_____________________________________________________
int JArrayAlgorithmSimple::GlobalIndex()
{
  int iG = 0;
  for (int i = 0; i < Dimension(); i++) // Index is checked by fCMD
    iG += Index(i) * fDimFactor[i];
  // TODO check iG
  return iG;
}
void JArrayAlgorithmSimple::ReverseIndex(int iG)
{
  int n = iG;
  for (int i = 0; i < Dimension(); i++) {
    int n1 = static_cast<int>(n / fDimFactor[i]);
    fCMD->SetIndex(n1, i);
    n -= n1 * fDimFactor[i];
  }
}
void* JArrayAlgorithmSimple::GetItem()
{
  return fArray[GlobalIndex()];
}
void JArrayAlgorithmSimple::SetItem(void* item)
{
  fArray[GlobalIndex()] = item;
}

//////////////////////////////////////////////////////
//  JTH1
//////////////////////////////////////////////////////
//_____________________________________________________
JTH1::JTH1() : fDirectory(NULL),
               fSubDirectory(NULL),
               fHMG(NULL),
               fTemplate(NULL),
               fBins(0)
{
  // default constructor
  fName = "JTH1";
}

//_____________________________________________________
JTH1::JTH1(TString config, JHistManager* hmg) : fDirectory(NULL),
                                                fSubDirectory(NULL),
                                                fHMG(NULL),
                                                fTemplate(NULL),
                                                fBins(0)
{
  // constructor
  std::vector<TString> t = Tokenize(config, " \t,");
  TString type = t[0];
  SetName(t[1]);
  SetTitle(t[2]);
  fTitle.ReplaceAll("\"", "");
  SetFullOption(t[3]);
  fMode = HasOption("mode", "Single") ? kSingle : kNormal;
  AddToManager(hmg);
  TString s;
  for (UInt_t i = 4; i < t.size(); i++)
    s += " " + t[i];
  AddDim(s);
  FixBin();
}
//_____________________________________________________
JTH1::~JTH1()
{
  // destructor
  if (fNGenerated == 0 && fTemplate)
    delete fTemplate;
}

//_____________________________________________________
JTH1::JTH1(const JTH1& obj) : JArrayBase(),
                              fDirectory(obj.fDirectory),
                              fSubDirectory(obj.fSubDirectory),
                              fHMG(obj.fHMG),
                              fTemplate(obj.fTemplate),
                              fBins(obj.fBins)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
JTH1& JTH1::operator=(const JTH1& obj)
{
  // assignment operator
  if (this != &obj) {
    // TODO: proper implementation
  }
  return *this;
}

//_____________________________________________________
int JTH1::AddDim(JBin* bin)
{
  int ndim = this->JArrayBase::AddDim(bin->Size());
  fBins.resize(ndim, NULL);
  fBins[ndim - 1] = bin;
  return ndim;
}

int JTH1::AddDim(TString v)
{
  if (v == "END") {
    FixBin();
  } else {
    std::vector<TString> o = Tokenize(v, "\t ,");
    for (UInt_t i = 0; i < o.size(); i++) {
      TString& s = o[i];
      if (s.Length() == 0)
        continue;
      if (s.IsFloat()) { // TODO IsInt? IsDigit?
        AddDim(s.Atoi());
        continue;
      }
      JBin* b = NULL;
      if (fHMG)
        b = fHMG->GetBin(s);
      if (b)
        this->AddDim(b);
      else
        JERROR("Wrong terminator of Array : \"%s\" in %s", s.Data(), fName.Data());
    }
  }
  return Dimension();
}
//_____________________________________________________
Int_t JTH1::Write()
{
  TDirectory* owd = gDirectory;
  InitIterator();
  void* item;
  if (fSubDirectory)
    fSubDirectory->cd();
  // else fDirectory->cd();
  while (Next(item)) {
    if (!item)
      continue;
    TH1* obj = static_cast<TH1*>(item);
    obj->Write();
    // obj->Write( 0, TObject::kOverwrite );
  }
  if (owd != gDirectory)
    owd->cd();
  return 0;
}
//_____________________________________________________
TString JTH1::GetString()
{
  TString s = Form("%s\t%s\t\"%s\"\t%s\t",
                   ClassName(), fName.Data(), fTitle.Data(), fOption.Data());
  for (int i = 0; i < Dimension(); i++) {
    if (fBins.size() > i && fBins[i] != NULL) {
      s += " " + fBins[i]->GetName();
    } else {
      s += TString(" ") + Form("%d", SizeOf(i));
    }
  }
  return s;
}
//_____________________________________________________
void JTH1::FixBin()
{
  this->JArrayBase::FixBin();

  if (!fHMG) {
    AddToManager(JHistManager::CurrentManager());
  }
  if (!fDirectory)
    fDirectory = fHMG->GetDirectory();
}

//_____________________________________________________
void JTH1::AddToManager(JHistManager* hmg)
{
  if (fHMG)
    return; // TODO handle error
  fHMG = hmg;
  hmg->Add(this);
}
//_____________________________________________________
void JTH1::Print()
{
  // std::cout<<"*"<<GetString()<<std::endl;
  LOGF(info, "*%s", GetString().Data());
  // TODO more details.
}
//_____________________________________________________
void JTH1::SetTemplate(TH1* h)
{
  if (fTemplate)
    return; /// TDOO give error
  fTemplate = static_cast<TH1*>(h->Clone());
  fTemplate->Sumw2();
  fTemplate->SetDirectory(0);
  fName = h->GetName();
  fTitle = h->GetTitle();
}
//_____________________________________________________
TString JTH1::BuildName()
{
  TString name = fName;
  if (!HasOption("Single")) {
    for (UInt_t i = 0; i < Dimension(); i++) {
      name += ((fBins.size() > i && fBins[i] != NULL) ? fBins[i]->GetIndexName() : "H") + Form("%02d", Index(i));
    }
  }
  return name;
}
//_____________________________________________________
TString JTH1::BuildTitle()
{
  TString title = fTitle;
  for (int i = 0; i < Dimension(); i++)
    title += ((static_cast<int>(fBins.size()) > i && fBins[i] != NULL) ? " " + fBins[i]->BuildTitle(Index(i)) : "") + Form("%02d", Index(i));
  return title;
}
//_____________________________________________________
void* JTH1::BuildItem()
{
  TDirectory* owd = gDirectory;
  gROOT->cd();
  TString name = BuildName();
  TH1* item = NULL;
  if (!fSubDirectory) {
    if (!HasOption("dir")) {
      fSubDirectory = fDirectory;
    } else {
      fSubDirectory = fDirectory->GetDirectory(fName);
      if (!fSubDirectory && !IsLoadMode()) {
        fSubDirectory = fDirectory->mkdir(fName);
      }
    }
  }
  if (IsLoadMode()) {
    // if( fSubDirectory ) JDEBUG(2, fSubDirectory->GetName() );
    if (fSubDirectory)
      item = dynamic_cast<TH1*>(fSubDirectory->Get(name));
    if (!item) {
      void** rawItem = fAlg->GetRawItem();
      InitIterator();
      void* tmp;
      while (Next(tmp)) {
        item = static_cast<TH1*>(fSubDirectory->Get(BuildName()));
        if (item)
          break;
      }
      if (item) {
        item = dynamic_cast<TH1*>((static_cast<TH1*>(item))->Clone(name));
        if (!item) {
          JERROR("Any of %s doesn't exists. I need at least one", fName.Data());
          return NULL;
        }
        item->Reset();
        item->SetTitle(BuildTitle());
        item->SetDirectory(0);
        *rawItem = item;
      }
    }
    if (!item) {
      JERROR("Any of %s doesn't exists. I need at least one", fName.Data());
      return NULL;
    }
  } else { //  Gen Mode
    TH1* titem = NULL;
    if (fNGenerated == 0)
      titem = fTemplate;
    else
      titem = static_cast<TH1*>(fTemplate->Clone());
    titem->SetDirectory(fSubDirectory);
    titem->Reset();
    titem->SetName(BuildName());
    titem->SetTitle(BuildTitle());
    fNGenerated++;
    item = titem;
  }
  if (item)
    fAlg->SetItem(item);
  owd->cd();
  return item;
}
//_____________________________________________________
bool JTH1::IsLoadMode()
{
  return fHMG->IsLoadMode();
}

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// JTH1Derived                                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
template <typename T>
JTH1Derived<T>::JTH1Derived() : JTH1(), fPlayer(this)
{
}
template <typename T>
JTH1Derived<T>::~JTH1Derived()
{
}

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// JHistManager                                                       //
//                                                                      //
// Array Base Class                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
JHistManager::JHistManager(TString name, TString dirname) : JNamed(name, "", "", 0),
                                                            fIsLoadMode(false),
                                                            fDirectory(gDirectory),
                                                            fConfigStr(),
                                                            fBin(0),
                                                            fHist(0),
                                                            fManager(0),
                                                            fBinNames(0),
                                                            fBinConfigs(0),
                                                            fHistNames(0),
                                                            fHistConfigs(0)
{
  // constructor
  if (dirname.Length() == 0)
    dirname = name;
  if (dirname.Length() > 0) {
    fDirectory = static_cast<TDirectory*>(gDirectory->Get(dirname));
    if (fDirectory) {
      LOGF(warning, "Hist directory %s exists", dirname.Data());
      // gSystem->Exit(1);  // We might actually want the directory to exist, so no exit
    }
    if (!fDirectory) {
      fDirectory = gDirectory->mkdir(dirname);
    }
  }
  if (!fDirectory) {
    JERROR("Fail to generate Hist directory %s", dirname.Data());
    // gSystem->Exit(1);
  }
  this->cd();
}

//_____________________________________________________
JHistManager::JHistManager(const JHistManager& obj) : JNamed(obj.fName, obj.fTitle, obj.fOption, obj.fMode),
                                                      fIsLoadMode(obj.fIsLoadMode),
                                                      fDirectory(obj.fDirectory),
                                                      fConfigStr(obj.fConfigStr),
                                                      fBin(obj.fBin),
                                                      fHist(obj.fHist),
                                                      fManager(obj.fManager),
                                                      fBinNames(obj.fBinNames),
                                                      fBinConfigs(obj.fBinConfigs),
                                                      fHistNames(obj.fHistNames),
                                                      fHistConfigs(obj.fHistConfigs)
{
  // copy constructor TODO: proper handling of pointer data members
}

//_____________________________________________________
JHistManager& JHistManager::operator=(const JHistManager& obj)
{
  // assignment operator
  if (this != &obj) {
    // TODO: proper implementation
  }
  return *this;
}

JHistManager* JHistManager::GlobalManager()
{
  static JHistManager* singleton = new JHistManager("GlobalHistManager");
  return singleton;
}

JHistManager* JHistManager::CurrentManager(JHistManager* hmg)
{
  static JHistManager* currentManager = NULL; //;JHistManager::GlobalManager();
  if (hmg)
    currentManager = hmg;
  return currentManager;
}

JBin* JHistManager::GetBuiltBin(TString s)
{
  for (UInt_t i = 0; i < fBin.size(); i++)
    if (fBin[i]->GetName() == s)
      return fBin[i];
  return NULL;
}
JBin* JHistManager::GetBin(TString s)
{
  JBin* h = GetBuiltBin(s);
  if (h)
    return h;
  for (UInt_t i = 0; i < GetNBin(); i++)
    if (fBinNames[i] == s) {
      return new JBin(fBinConfigs[i], this);
    }
  return NULL;
}
JTH1* JHistManager::GetBuiltTH1(TString s)
{
  for (UInt_t i = 0; i < fHist.size(); i++)
    if (fHist[i]->GetName() == s)
      return fHist[i];
  return NULL;
}
// Note: Returning NULL crashes the code, something should be done about this.
// The error given by compiler is: non-const lvalue reference to type 'JTH1D' (aka 'JTH1Derived<TH1D>') cannot bind to a temporary of type 'long'
// The reoson for crash is that NULL cannot be used in dynamic_cast<JTH1D&>
JTH1* JHistManager::GetTH1(TString s)
{
  JTH1* h = GetBuiltTH1(s);
  if (h)
    return h;
  for (UInt_t i = 0; i < GetNHist(); i++)
    if (fHistNames[i] == s) {
      if (fHistConfigs[i].BeginsWith("JTH1D"))
        return new JTH1D(fHistConfigs[i], this);
      if (fHistConfigs[i].BeginsWith("JTH2D"))
        return new JTH2D(fHistConfigs[i], this);
      if (fHistConfigs[i].BeginsWith("JTH3D"))
        return new JTH3D(fHistConfigs[i], this);
      if (fHistConfigs[i].BeginsWith("JTProfile"))
        return new JTProfile(fHistConfigs[i], this);
    }
  return NULL;
}
void JHistManager::Add(JBin* o)
{
  if (!o)
    return;
  if (GetBuiltBin(o->GetName()))
    return; // TODO error handle
  fBin.push_back(o);
}
void JHistManager::Add(JTH1* o)
{
  if (!o)
    return;
  if (GetBuiltTH1(o->GetName()))
    return; // TODO error handle
  fHist.push_back(o);
}
void JHistManager::Print()
{
  if (IsLoadMode()) {
    // cout<<fConfigStr<<endl;
    LOGF(info, "%s", fConfigStr.Data());
    return;
  }
  LOGF(info, "============ JHistManager : %s ===================\n", fName.Data());
  LOGF(info, "\n---- JBin ----\n");
  for (UInt_t i = 0; i < GetNBin(); i++) {
    fBin[i]->Print();
  }
  LOGF(info, "\n---- JTH1 ----\n");
  for (UInt_t i = 0; i < GetNHist(); i++) {
    fHist[i]->Print();
  }
}
void JHistManager::Write()
{
  for (UInt_t i = 0; i < GetNHist(); i++)
    fHist[i]->Write();
}

void JHistManager::WriteConfig()
{
  TDirectory* owd = fDirectory;
  // cout<<"DEBUG_T1: "<<fDirectory<<endl;
  // cout<<"DEBUG_T2: "<<GetName()<<"\t"<<fDirectory->GetName()<<endl;
  // exit(1);
  //  TODO 1.Error Check 2.Duplicaition check
  TDirectory* fHistConfigDir = fDirectory->mkdir("HistManager");
  fHistConfigDir->cd();
  TObjString* config = new TObjString(GetString().Data());
  config->Write("Config");
  owd->cd();
}

int JHistManager::LoadConfig()
{
  SetLoadMode(true);
  TObjString* strobj = static_cast<TObjString*>(fDirectory->Get("HistManager/Config"));
  if (!strobj)
    return 0; // TODO
  TString config = strobj->String();
  fConfigStr = config;
  vector<TString> lines = Tokenize(config, "\n");
  LOGF(info, "Read Config.%d objects found", (int)lines.size());
  for (UInt_t i = 0; i < lines.size(); i++) {
    TString line = lines.at(i);
    std::vector<TString> t = Tokenize(line, " \t,");
    if (line.BeginsWith("JBin")) {
      fBinNames.push_back(t[1]);
      fBinConfigs.push_back(line);
    } else if (line.BeginsWith("J")) {
      fHistNames.push_back(t[1]);
      fHistConfigs.push_back(line);
    }
  }
  return 1;
}

bool JHistManager::HistogramExists(TString name)
{
  for (UInt_t i = 0; i < fHistNames.size(); i++) {
    if (fHistNames[i] == name)
      return true;
  }
  return false;
}

//////////////////////////////////////////////////////
//  Utils
//////////////////////////////////////////////////////
vector<TString> Tokenize(TString s, TString d, int quote)
{
  // int nd = d.Length();
  bool flagBeforeToken = 0;
  bool inQuote = 0;
  TString tok = "";
  vector<TString> toks;
  s += d[0];
  for (int i = 0; i < s.Length(); i++) {
    if (quote == 1 && s[i] == '\"') {
      inQuote = !inQuote;
    }
    if (d.First(s[i]) != kNPOS && !inQuote) {
      if (flagBeforeToken == 0 && tok.Length() > 0) {
        toks.push_back(tok);
        tok.Clear();
        flagBeforeToken = 1;
      }
    } else {
      tok += s[i];
      flagBeforeToken = 0;
    }
  }
  return toks;
}

TString Join(vector<TString>& ss, TString del)
{
  if (ss.size() < 1)
    return "";
  TString s = ss[0];
  for (UInt_t i = 1; i < ss.size(); i++)
    s += del + ss[i];
  return s;
}

template class JTH1Derived<TH1D>;
template class JTH1Derived<TH2D>;
template class JTH1Derived<TH3D>;
template class JTH1Derived<TProfile>;

bool OutOf(int i, int x, int y) { return (i < x || i > y); }

#include <TFile.h>

void ttestJArray()
{
  JHistManager* fHMG;
  JBin fCentBin;
  JBin fVtxBin;
  JBin fPTtBin;
  JBin fPTaBin;
  JBin fXEBin;
  JBin fKLongBin;
  JBin fRGapBin;
  JBin fEtaGapBin;
  JBin fPhiGapBin;
  JBin fMassBin;
  JBin fTypBin;
  JBin fTypBin3;
  JBin fPairPtBin;
  JTH1D fhTriggPtBin;
  JTH1D fhTriggMult;
  JTH1D fhIphiTrigg;
  JTH1D fhIetaTrigg;
  JTH2D test1;
  JTProfile test2;

  TFile* f = new TFile("test.root", "RECREATE");
  fHMG = JHistManager::GlobalManager();
  fCentBin.Set("Cent", "C", "C %2.0f-%2.0f%%").SetBin("0 100");
  fVtxBin.Set("Vtx", "V", "Vtx %.0f-%.0f").SetBin("-10 10");
  fPTtBin.Set("PTt", "T", "p_{Tt} %.1f-%.1f").SetBin("3 5 8 10 15 20");
  fPTaBin.Set("PTa", "A", "p_{Tt} %.1f-%.1f").SetBin("3 5 8 10 15 20");

  fhTriggMult
    << TH1D("hTriggMult", "", 100, -0.5, 99.5)
    << fCentBin << fPTtBin << "END";
  fhIphiTrigg
    << TH1D("fhIphiTrigg", "", 3, -0.1, 0.1)
    << fCentBin << fPTtBin << "END";
  fhIetaTrigg
    << TH1D("hIetaTrigg", "", 80, -5, 5)
    << fCentBin << fPTtBin << "END"; // inclusive eta
  fhTriggPtBin
    << TH1D("hTriggPtBin", "", 10, 0, 10)
    << fCentBin << fVtxBin << fPTtBin << "END";

  fhTriggMult[0][0]->Fill(1);
  fHMG->Print();

  f->Write();
  fHMG->Write();
  fHMG->WriteConfig();
}
