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

// NOLINTBEGIN(readability/casting, readability/braces)
#ifndef PWGCF_GENERICFRAMEWORK_CORE_MYCANVAS_H_
#define PWGCF_GENERICFRAMEWORK_CORE_MYCANVAS_H_

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

// Helper to detect if T has SetFillStyle(int)
template <typename T, typename = void>
struct has_GetNbinsX : std::false_type {
};

template <typename T>
struct has_GetNbinsX<T, std::void_t<decltype(std::declval<T>().GetNbinsX())>> : std::true_type {
};

using namespace std;
class MyCanvas
{
 public:
  using ObjType = std::variant<TH1*, TGraphErrors*, TH1D*, TGraphAsymmErrors*>;
  MyCanvas();
  explicit MyCanvas(const char* name, double w = 500, double h = 500, double lm = 0.05, double bm = 0.08, double tm = 0.04, double rm = 0.04);
  ~MyCanvas();
  TCanvas* BuildSingleCanvas();
  TCanvas* BuildRatioCanvas(double pr);
  TCanvas* BuildMultiRatioCanvas(int N, float pr, bool separatePad = false);
  TCanvas* BuildNbyMCanvas(int N, int M, bool separatePad = false);
  TCanvas* BuildNbyMRatioCanvas(int N, int M, double pr = 0.33, bool separatePad = false, bool alternateAxisSide = false);
  void cd(int Npad = 1);
  TPad* GetPad(int Npad = 1);
  void cdCanvas() { mCanvas->cd(); };
  void SetLogy(int Npad = 1)
  {
    TPad* tmp = (TPad*)mCanvas->GetPad(Npad);
    tmp->SetLogy();
  }
  void SetLogx(int Npad = 1)
  {
    TPad* tmp = (TPad*)mCanvas->GetPad(Npad);
    tmp->SetLogx();
  }
  void SetLogz(int Npad = 1)
  {
    TPad* tmp = (TPad*)mCanvas->GetPad(Npad);
    tmp->SetLogz();
  }
  void SaveAs(const char* name, const char* opt);
  void Print(const char* name) { mCanvas->Print(name); }
  void Write() { mCanvas->Write(); }
  void AddLegEntriesAbove(bool add) { addAbove = add; }
  template <typename T>
  void AddLegEntry(T entry, TString name, const char* drawopt, int legInd)
  {
    BuildLegend(0.5, 0.5, 0.6, 0.6, legInd, 1);
    mLeg[legInd]->AddEntry(entry, name, drawopt);
    return;
  }
  void AddObject(ObjType stat, int Npad = 0, int color = -1, int style = -1, ObjType sys = {});
  void AddModel(ObjType model, int Npad, int color, int fillstyle, float alpha = 1.);
  void AddLine(double x1, double y1, double x2, double y2, int nPad = 0);
  void AddText(double x, double y, const char* text, bool NDC = true, int textsize = 0, int textcolor = kBlack);
  void Draw(const char* opt = "");
  void DrawLegend(double x1, double y1, int legInd = 0, int nCol = 1);
  void DrawLegend(double x1, double y1, double x2, double y2, int legInd = 0, int nCol = 1);
  void DrawLegend(TString loc, int legInd = 0, int nCol = 1);
  void Update() { mCanvas->Update(); }
  void UseSameYscale(bool newval) { mUseSameYscale = newval; }
  int GetTextSize() { return mTextSize; }
  void SetYScale(double newval, double ratioval = 0)
  {
    mYscale = newval;
    if (ratioval)
      mRatioYScale = ratioval;
    else
      mRatioYScale = newval;
  } // Controls how large a percentage above/below the max/min values of the histograms should be drawn
  void SetYLimits(double ymin = -999, double ymax = -999)
  {
    mYmin = ymin;
    mYmax = ymax;
  }
  void SetYLimits(double ymin, double ymax, int padNo) { padlimits[padNo] = std::make_pair(ymin, ymax); }
  void SetYRatioLimits(double yRatioMin, double yRatioMax)
  {
    mYRatioMin = yRatioMin;
    mYRatioMax = yRatioMax;
  }
  void SetTextSize(int newval) { mTextSize = newval; }
  void SetTitleOffset(double newval, TString axis = "Y")
  {
    if (axis.EqualTo("Y"))
      mOffset[1] = newval;
    else
      mOffset[0] = newval;
  }
  void SetTitle(TString newval) { mTitle = newval; }
  void SetLabelOffset(double newval) { mLabelOffset = newval; }
  void ForceXLimits(double xmin, double xmax)
  {
    mXmin = xmin;
    mXmax = xmax;
  };
  void SetModelDrawOption(std::string option) { fModelDrawOpt = option; }
  void SetAdjustTicks(bool newval) { mAdjustTicks = newval; }
  void DrawGraphWithErrorOutline(ObjType gvar, int lineWidth, int lineColor);
  bool isInitialised() { return initialised; }
  void SetLineWidth(int w) { mLineWidth = w; }

  template <typename T>
  TGraphErrors* HistToGraph(T h); // public, don't want to reimplement elsewhere
  TGraphErrors* HistToGraph(TGraphErrors* gr)
  {
    gr->SetFillColorAlpha(gr->GetMarkerColor(), 0.4);
    return gr;
  }
  struct DrawObject {
    ObjType stat;
    ObjType sys;
    int color;
    int style;
    double min;
    double max;
  };
  struct ModelObject {
    ObjType model;
    int color;
    int fillstyle;
    float alpha;
    double min;
    double max;
  };

 private:
  bool initialised;
  const char* mName;
  double mWidth;
  double mHeight;
  double mTopMargin;
  double mBottomMargin;
  double mLeftMargin;
  double mRightMargin;
  int mM;
  int mN;
  bool addAbove;
  TCanvas* mCanvas;
  vector<TLegend*> mLeg;
  int mTextSize;
  vector<double> mOffset;
  TString mTitle;
  double mLabelOffset;
  double mYmin;
  double mYmax;
  double mXmin;
  double mXmax;
  double mYRatioMin;
  double mYRatioMax;
  bool mAddCentLabel;
  double mYscale;
  double mRatioYScale;
  bool mUseSameYscale;
  bool mSeparatePad;
  bool mAlternateAxisSide;
  bool fRatioCanvas;
  bool mAdjustTicks;
  bool mLineWidth;
  int mLegIndex;
  std::string fModelDrawOpt;
  vector<int> mStyles;
  vector<int> mColors;
  map<int, double> stylemap;
  vector<int> mFstyles;
  void BuildLegend(double x1, double y1, double x2, double y2, int legInd = 0, int nCol = 0);
  void ApplyFormat();
  void SetLimits();
  float SafeAdd(float x, float eps);
  float SafeSub(float x, float eps);
  void ApplyPadStyle(TPad* pad);
  static bool isGraphType(const ObjType& obj)
  {
    return std::visit(
      [](auto* ptr) {
        using T = std::decay_t<decltype(ptr)>;
        return std::is_same_v<T, TGraphErrors*> || std::is_same_v<T, TGraphAsymmErrors*>;
      },
      obj);
  }
  vector<vector<DrawObject>> objVector;
  vector<vector<ModelObject>> modelVector;
  std::map<int, std::pair<double, double>> padlimits; //!

  ClassDef(MyCanvas, 1);
};
#endif // PWGCF_GENERICFRAMEWORK_CORE_MYCANVAS_H_
// NOLINTEND(readability/casting, readability/braces)
