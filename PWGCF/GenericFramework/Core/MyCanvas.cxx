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

// NOLINTBEGIN(readability/casting, readability/braces, whitespace/parens, whitespace/newline, whitespace/comma, whitespace/operators)
#include "MyCanvas.h"

#include <TROOT.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
ClassImp(MyCanvas)

  MyCanvas::MyCanvas() : initialised(0),
                         mName("c"),
                         mWidth(800),
                         mHeight(600),
                         mTopMargin(0.0),
                         mBottomMargin(0.0),
                         mLeftMargin(0.0),
                         mRightMargin(0.0),
                         mM(1),
                         mN(1),
                         addAbove(false),
                         mCanvas(0),
                         mLeg(0),
                         mTextSize(30),
                         mOffset{1, 1.4},
                         mTitle(""),
                         mLabelOffset(0),
                         mYmin(-999),
                         mYmax(-999),
                         mXmin(-999),
                         mXmax(-999),
                         mYRatioMin(-999),
                         mYRatioMax(-999),
                         mYscale(0.25),
                         mRatioYScale(0.25),
                         mUseSameYscale(true),
                         mSeparatePad(false),
                         mAlternateAxisSide(false),
                         mAdjustTicks(false),
                         mLineWidth(1),
                         mLegIndex(0),
                         fModelDrawOpt("zsame3"),
                         mStyles{kFullCircle, kFullSquare, kFullDiamond, kFullStar, kFullCross, kFullCrossX, kFullDoubleDiamond, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenStar, kOpenCross, kOpenCrossX, kOpenDoubleDiamond},
                         // mColors{kBlack,kP10Blue,kP10Yellow,kP10Red,kP10Gray,kP10Violet,kP10Brown,kP10Orange,kP10Green,kP10Ash,kP10Cyan},
                         stylemap{{kFullCircle, 1}, {kFullSquare, 1}, {kFullDiamond, 1.5}, {kFullStar, 1.5}, {kFullCross, 1.5}, {kFullCrossX, 1.5}, {kFullDoubleDiamond, 2}, {kFullTriangleUp, 1.5}, {kOpenCircle, 1}, {kOpenSquare, 1}, {kOpenDiamond, 1.5}, {kOpenStar, 1.5}, {kOpenCross, 1.5}, {kOpenCrossX, 1.5}, {kOpenDoubleDiamond, 2}, {kOpenTriangleUp, 1.5}},
                         mFstyles{1001, -1, 3002, 3008, 3014}
{
  gStyle->SetErrorX(0.0);
}
MyCanvas::MyCanvas(const char* name, double w, double h, double lm, double bm, double tm, double rm) : initialised(0),
                                                                                                       mName(name),
                                                                                                       mWidth(w),
                                                                                                       mHeight(h),
                                                                                                       mTopMargin(tm),
                                                                                                       mBottomMargin(bm),
                                                                                                       mLeftMargin(lm),
                                                                                                       mRightMargin(rm),
                                                                                                       mM(1),
                                                                                                       mN(1),
                                                                                                       addAbove(false),
                                                                                                       mCanvas(0),
                                                                                                       mLeg(0),
                                                                                                       mTextSize(30),
                                                                                                       mOffset{1, 1.4},
                                                                                                       mTitle(""),
                                                                                                       mLabelOffset(0),
                                                                                                       mYmin(-999),
                                                                                                       mYmax(-999),
                                                                                                       mXmin(-999),
                                                                                                       mXmax(-999),
                                                                                                       mYRatioMin(-999),
                                                                                                       mYRatioMax(-999),
                                                                                                       mYscale(0.25),
                                                                                                       mRatioYScale(0.25),
                                                                                                       mUseSameYscale(true),
                                                                                                       mSeparatePad(false),
                                                                                                       mAlternateAxisSide(false),
                                                                                                       mAdjustTicks(false),
                                                                                                       mLineWidth(1),
                                                                                                       mLegIndex(0),
                                                                                                       fModelDrawOpt("zsame3"),
                                                                                                       mStyles{kFullCircle, kFullSquare, kFullDiamond, kFullStar, kFullCross, kFullCrossX, kFullDoubleDiamond, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenStar, kOpenCross, kOpenCrossX, kOpenDoubleDiamond},
                                                                                                       // mColors{kBlack,kP10Blue,kP10Yellow,kP10Red,kP10Gray,kP10Violet,kP10Brown,kP10Orange,kP10Green,kP10Ash,kP10Cyan},
                                                                                                       stylemap{{kFullCircle, 1}, {kFullSquare, 1}, {kFullDiamond, 1.5}, {kFullStar, 1.5}, {kFullCross, 1.5}, {kFullCrossX, 1.5}, {kFullDoubleDiamond, 2}, {kFullTriangleUp, 1.5}, {kOpenCircle, 1}, {kOpenSquare, 1}, {kOpenDiamond, 1.5}, {kOpenStar, 1.5}, {kOpenCross, 1.5}, {kOpenCrossX, 1.5}, {kOpenDoubleDiamond, 2}, {kOpenTriangleUp, 1.5}},
                                                                                                       mFstyles{1001, -1, 3002, 3008, 3014}
{
  gStyle->SetErrorX(0.0);
}
MyCanvas::~MyCanvas()
{
  delete mCanvas;
  for (auto* legend : mLeg)
    delete legend;
}
TCanvas* MyCanvas::BuildSingleCanvas()
{
  // if(rm<0) mRightMargin = mTopMargin;
  mCanvas = new TCanvas(mName, mName, mWidth, mHeight);
  TPad* pad = new TPad("pad", "pad", 0, 0, 1, 1);
  pad->SetNumber(1);
  pad->SetMargin(mLeftMargin, mRightMargin, mBottomMargin, mTopMargin);
  mCanvas->cd();
  pad->Draw();
  mCanvas->cd();
  objVector.resize(1);
  modelVector.resize(1);
  initialised = true;
  return mCanvas;
}
TCanvas* MyCanvas::BuildRatioCanvas(double pr)
{
  fRatioCanvas = true;
  mUseSameYscale = false;
  mCanvas = new TCanvas(mName, mName, mWidth, mHeight);
  double th = mHeight * (1 - pr) - mTopMargin - mBottomMargin;
  double bh = mHeight * pr - mTopMargin - mBottomMargin;
  double bhrel = bh / mHeight;
  mN = 2;

  TPad* top = new TPad("top", "top", 0, bhrel, 1, 1);
  top->SetNumber(1);
  top->SetMargin(mLeftMargin, mTopMargin, 0, (mTopMargin * mHeight) / th);
  mCanvas->cd();
  top->Draw();
  mCanvas->cd();
  TPad* bot = new TPad("bot", "bot", 0, 0, 1, bhrel);
  bot->SetNumber(2);
  bot->SetMargin(mLeftMargin, mTopMargin, (mBottomMargin * mHeight) / bh, 0);
  bot->Draw();
  mCanvas->cd();
  objVector.resize(2);
  modelVector.resize(2);
  initialised = true;
  return mCanvas;
}
TCanvas* MyCanvas::BuildMultiRatioCanvas(int N, float pr, bool separatePad)
{
  mSeparatePad = separatePad;
  mCanvas = new TCanvas(mName, mName, mWidth, mHeight);

  // Calculate height for equal visible area per pad
  float H = (pr - mBottomMargin) / N;

  float currentTop = pr;   // Start from top of small pad area
  const float epsilon = 0; // Small overlap to avoid white gaps
  for (int i = 0; i <= N; ++i) {
    float vposu, vposd, vfactor;
    float vmard = 0.0, vmaru = 0.0;

    float hposl = 0.0, hposr = 1.0;
    float hfactor = hposr - hposl;
    float hmarl = mLeftMargin / hfactor;
    float hmarr = mRightMargin / hfactor;

    if (i == 0) {
      // Main pad
      vposd = pr;
      vposu = 1.0;
      vfactor = vposu - vposd;
      vmaru = mTopMargin / vfactor;
    } else {
      // Small pad
      if (i == N) {
        // Last pad
        vposd = SafeSub(currentTop - (H + mBottomMargin), epsilon);
        vposu = SafeAdd(currentTop, epsilon);
        vfactor = vposu - vposd;
        vmard = mBottomMargin / vfactor;
      } else {
        vposd = SafeSub(currentTop - H, epsilon);
        vposu = SafeAdd(currentTop, epsilon);
        vfactor = vposu - vposd;
      }
      currentTop = vposd; // no epsilon subtraction here to avoid accumulating underflow
    }

    TString padName = (i == 0) ? "padMain" : Form("pad%i", i);
    TPad* pad = new TPad(padName, padName, hposl, vposd, hposr, vposu);
    pad->SetNumber(i + 1);

    pad->SetLeftMargin(hmarl);
    pad->SetRightMargin(hmarr);
    pad->SetTopMargin(vmaru);
    pad->SetBottomMargin(vmard);

    pad->SetFrameBorderMode(0);
    pad->SetBorderMode(0);
    pad->SetBorderSize(0);
    // ApplyPadStyle(pad);
    pad->SetTicks(1, 1);
    mCanvas->cd();
    pad->Draw();
  }

  mCanvas->cd();
  objVector.resize(N + 1);
  modelVector.resize(N + 1);
  initialised = true;
  return mCanvas;
}

TCanvas* MyCanvas::BuildNbyMCanvas(int N, int M, bool separatePad)
{
  mSeparatePad = separatePad;
  mCanvas = new TCanvas(mName, mName, mWidth, mHeight);
  // Setup Pad layout:
  // Float_t vSpacing = (mSeparatePad)?0.05:0.0;
  int Nvmar = (mSeparatePad) ? M : 1;
  Float_t vStep = (1. - Nvmar * mBottomMargin - Nvmar * mTopMargin) / M; //- (M-1) * vSpacing

  // Float_t hSpacing = (mSeparatePad)?0.02:0.0;
  int Nmar = (mSeparatePad) ? N : 1;
  Float_t hStep = (1. - Nmar * mLeftMargin - Nmar * mRightMargin) / N; //- (N-1) * hSpacing

  Float_t vposd, vposu, vmard, vmaru, vfactor;
  Float_t hposl, hposr, hmarl, hmarr, hfactor;

  // printf("hStep = %.3f\nvStep = %.3f\n",hStep,vStep);
  for (Int_t i = 0; i < N; i++) {

    if (i == 0) {
      hposl = 0.0;
      hposr = mLeftMargin + hStep;
      if (i == N - 1 || mSeparatePad)
        hposr += mRightMargin;
      hfactor = hposr - hposl;
      hmarl = mLeftMargin / hfactor;
      hmarr = (i == N - 1 || mSeparatePad) ? mRightMargin / hfactor : 0.0;
    } else if (i == N - 1) {
      hposl = hposr;
      hposr = hposl + hStep + mRightMargin;
      if (mSeparatePad)
        hposr += mLeftMargin;
      hfactor = hposr - hposl;
      hmarl = (mSeparatePad) ? mLeftMargin / hfactor : 0.0;
      hmarr = mRightMargin / hfactor;
    } else {
      hposl = hposr;
      hposr = hposl + hStep;
      if (mSeparatePad)
        hposr += mLeftMargin + mRightMargin;
      hfactor = hposr - hposl;
      hmarl = (mSeparatePad) ? mLeftMargin / hfactor : 0.0;
      hmarr = (mSeparatePad) ? mRightMargin / hfactor : 0.0;
    }

    // printf("i = %i\n hposl = %.3f, hposr = %.3f, hfactor = %.3f\n lmar = %.3f, rmar = %.3f\n",i,hposl,hposr,hfactor,hmarl,hmarr);
    for (Int_t j = 0; j < M; j++) {
      if (j == 0) {
        vposu = 1.;
        vposd = vposu - vStep - mTopMargin;
        if (mSeparatePad)
          vposd -= mBottomMargin;
        if (j == M - 1)
          vposd = 0.0;
        vfactor = vposu - vposd;
        vmard = (mSeparatePad || j == M - 1) ? mBottomMargin / vfactor : 0.0;
        vmaru = mTopMargin / vfactor; // 0.0
      } else if (j == M - 1) {
        vposu = vposd;
        vposd = 0.0;
        vfactor = vposu - vposd;
        vmard = mBottomMargin / vfactor; // 0.0;
        vmaru = (mSeparatePad) ? mTopMargin / vfactor : 0.0;
      } else {
        vposu = vposd;
        vposd = vposu - vStep;
        if (mSeparatePad)
          vposd -= (mBottomMargin + mTopMargin);
        vfactor = vposu - vposd;
        vmard = (mSeparatePad) ? mBottomMargin / vfactor : 0.0; // 0.0;
        vmaru = (mSeparatePad) ? mTopMargin / vfactor : 0.0;
      }
      // printf("i = %i and j = %i\n vposd = %.3f, vposu = %.3f, vfactor = %.3f\n dmar = %.3f, umar = %.3f\n",i,j,vposd,vposu,vfactor,vmard,vmaru);
      TPad* nPad = new TPad(Form("pad%i%i", i, j), Form("pad%i%i", i, j), hposl, vposd, hposr, vposu);
      nPad->SetNumber((i + 1) + N * j);
      nPad->SetLeftMargin(hmarl);
      nPad->SetRightMargin(hmarr);
      nPad->SetBottomMargin(vmard);
      nPad->SetTopMargin(vmaru);

      nPad->SetFrameBorderMode(0);
      nPad->SetBorderMode(0);
      nPad->SetBorderSize(0);
      mCanvas->cd();
      nPad->Draw();
    }
  }
  mCanvas->cd();
  objVector.resize(N * M);
  modelVector.resize(N * M);
  initialised = true;
  return mCanvas;
}
/* TCanvas* MyCanvas::BuildNbyMRatioCanvas(
    int N,
    int M,
    double pr,
    bool separatePad,
    bool alternateAxisSide)
{
    mSeparatePad = separatePad;
    mAlternateAxisSide = alternateAxisSide;
    fRatioCanvas = true;
    mUseSameYscale = false;

    mCanvas = new TCanvas(mName, mName, mWidth, mHeight);

    if (N <= 0 || M <= 0) {
        std::cout << "BuildNbyMRatioCanvas: N and M must be positive" << std::endl;
        return mCanvas;
    }

    if (pr <= 0.0 || pr >= 1.0) {
        std::cout << "BuildNbyMRatioCanvas: pr must be between 0 and 1" << std::endl;
        return mCanvas;
    }

    const double totalHorizontalMargins =
        separatePad
            ? N * (mLeftMargin + mRightMargin)
            : mLeftMargin + mRightMargin;

    const double totalVerticalMargins =
        separatePad
            ? M * (mTopMargin + mBottomMargin)
            : mTopMargin + mBottomMargin;

    const double hStep =
        (1.0 - totalHorizontalMargins) / static_cast<double>(N);

    const double vStep =
        (1.0 - totalVerticalMargins) / static_cast<double>(M);

    if (hStep <= 0.0 || vStep <= 0.0) {
        std::cout << "BuildNbyMRatioCanvas: Margins leave no drawable canvas area" << std::endl;
        return mCanvas;
    }

    const double mainFrameHeight  = (1.0 - pr) * vStep;
    const double ratioFrameHeight = pr * vStep;

    double xLeft = 0.0;

    for (int i = 0; i < N; ++i) {
        const double leftExtra =
            (separatePad || i == 0) ? mLeftMargin : 0.0;

        const double rightExtra =
            (separatePad || i == N - 1) ? mRightMargin : 0.0;

        const double padWidth = leftExtra + hStep + rightExtra;
        const double xRight = xLeft + padWidth;

        // ROOT margins are fractions of the complete TPad dimensions.
        const double leftPadMargin  = leftExtra / padWidth;
        const double rightPadMargin = rightExtra / padWidth;

        double yTop = 1.0;

        for (int row = 0; row < M; ++row) {
            const int mainIndex  = 2 * row;
            const int ratioIndex = mainIndex + 1;

            const double topExtra =
                (separatePad || row == 0) ? mTopMargin : 0.0;

            const double bottomExtra =
                (separatePad || row == M - 1)
                    ? mBottomMargin
                    : 0.0;

            const double mainPadHeight = topExtra + mainFrameHeight;
            const double mainBottom = yTop - mainPadHeight;

            TPad* mainPad = new TPad(
                Form("pad%i%i", i, mainIndex),
                Form("pad%i%i", i, mainIndex),
                xLeft,
                mainBottom,
                xRight,
                yTop
            );

            mainPad->SetNumber(
                2 * N * row + 2 * i + 1
            );

            mainPad->SetLeftMargin(leftPadMargin);
            mainPad->SetRightMargin(rightPadMargin);
            mainPad->SetTopMargin(topExtra / mainPadHeight);
            mainPad->SetBottomMargin(0.0);

            mainPad->SetFrameBorderMode(0);
            mainPad->SetBorderMode(0);
            mainPad->SetBorderSize(0);

            mCanvas->cd();
            mainPad->Draw();

            const double ratioPadHeight =
                ratioFrameHeight + bottomExtra;

            const double ratioBottom =
                mainBottom - ratioPadHeight;

            TPad* ratioPad = new TPad(
                Form("pad%i%i", i, ratioIndex),
                Form("pad%i%i", i, ratioIndex),
                xLeft,
                ratioBottom,
                xRight,
                mainBottom
            );

            ratioPad->SetNumber(
                2 * N * row + 2 * i + 2
            );

            ratioPad->SetLeftMargin(leftPadMargin);
            ratioPad->SetRightMargin(rightPadMargin);
            ratioPad->SetTopMargin(0.0);
            ratioPad->SetBottomMargin(
                bottomExtra / ratioPadHeight
            );

            ratioPad->SetFrameBorderMode(0);
            ratioPad->SetBorderMode(0);
            ratioPad->SetBorderSize(0);

            mCanvas->cd();
            ratioPad->Draw();

            yTop = ratioBottom;
        }

        xLeft = xRight;
    }

    mCanvas->cd();

    objVector.resize(2 * N * M);
    modelVector.resize(2 * N * M);

    initialised = true;
    return mCanvas;
} */
TCanvas* MyCanvas::BuildNbyMRatioCanvas(int N, int M, double pr, bool separatePad, bool alternateAxisSide)
{
  mSeparatePad = separatePad;
  mAlternateAxisSide = alternateAxisSide;
  fRatioCanvas = true;
  mUseSameYscale = false;
  mCanvas = new TCanvas(mName, mName, mWidth, mHeight);
  // Setup Pad layout:
  // Float_t vSpacing = (mSeparatePad)?0.05:0.0;
  int Nvmar = (mSeparatePad) ? M : 1;
  Float_t vStep = (1. - Nvmar * mBottomMargin - Nvmar * mTopMargin) / M; //- (M-1) * vSpacing
  // Float_t hSpacing = (mSeparatePad)?0.02:0.0;
  int Nmar = (mSeparatePad) ? N : 1;
  Float_t hStep = (1. - Nmar * mLeftMargin - Nmar * mRightMargin) / N; //- (N-1) * hSpacing
  Float_t vposd, vposu, vmard, vmaru, vfactor;
  Float_t hposl, hposr, hmarl, hmarr, hfactor;
  for (Int_t i = 0; i < N; i++) {
    if (i == 0) {
      hposl = 0.0;
      hposr = mLeftMargin + hStep;
      if (i == N - 1 || mSeparatePad)
        hposr += mRightMargin;
      hfactor = hposr - hposl;
      hmarl = mLeftMargin / hfactor;
      hmarr = (i == N - 1 || mSeparatePad) ? mRightMargin / hfactor : 0.0;
    } else if (i == N - 1) {
      hposl = hposr;
      hposr = hposl + hStep + mRightMargin;
      if (mSeparatePad)
        hposr += mLeftMargin;
      hfactor = hposr - hposl;
      hmarl = (mSeparatePad) ? mLeftMargin / hfactor : 0.0;
      hmarr = mRightMargin / hfactor;
    } else {
      hposl = hposr;
      hposr = hposl + hStep;
      if (mSeparatePad)
        hposr += mLeftMargin + mRightMargin;
      hfactor = hposr - hposl;
      hmarl = (mSeparatePad) ? mLeftMargin / hfactor : 0.0;
      hmarr = (mSeparatePad) ? mRightMargin / hfactor : 0.0;
    }
    // printf("i = %i\n hposl = %.3f, hposr = %.3f, hfactor = %.3f\n lmar = %.3f, rmar = %.3f\n",i,hposl,hposr,hfactor,hmarl,hmarr);
    for (Int_t j = 0; j < 2 * M; j++) {
      if (j == 0) {
        vposu = 1.;
        vposd = vposu - (1 - pr) * vStep - mTopMargin;
        if (mSeparatePad)
          vposd -= mBottomMargin;
        if (j == 2 * M - 1)
          vposd = 0.0;
        vfactor = vposu - vposd;
        // vmard = (mSeparatePad || j == 2*M-1)?mBottomMargin / vfactor:0.0;
        vmard = (j == 2 * M - 1) ? mBottomMargin / vfactor : 0.0;
        vmaru = mTopMargin / vfactor; // 0.0
      } else if (j == 2 * M - 1) {
        vposu = vposd;
        vposd = 0.;
        vfactor = vposu - vposd;
        vmard = mBottomMargin / vfactor; // 0.0;
        vmaru = 0.0;
      } else {
        if (j % 2) {
          vposu = vposd;
          vposd = vposu - pr * vStep; // - mBottomMargin;
          vfactor = vposu - vposd;
          vmard = (mSeparatePad) ? mBottomMargin / vfactor : 0.0; // 0.0;
          vmaru = 0.0;
        } else {
          vposu = vposd;
          vposd = vposu - (1 - pr) * vStep; // - mTopMargin;
          vfactor = vposu - vposd;
          vmard = 0.0; // 0.0;
          vmaru = (mSeparatePad) ? mTopMargin / vfactor : 0.0;
        }
      }
      TPad* nPad = new TPad(Form("pad%i%i", i, j), Form("pad%i%i", i, j), hposl, vposd, hposr, vposu);
      if (j % 2)
        nPad->SetNumber(2 * (i + 1) + N * (j - 1));
      else
        nPad->SetNumber((2 * i + 1) + N * j);
      nPad->SetLeftMargin(hmarl);
      nPad->SetRightMargin(hmarr);
      nPad->SetBottomMargin(vmard);
      nPad->SetTopMargin(vmaru);

      nPad->SetFrameBorderMode(0);
      nPad->SetBorderMode(0);
      nPad->SetBorderSize(0);
      mCanvas->cd();
      nPad->Draw();
    }
  }
  mCanvas->cd();
  objVector.resize(N * 2 * M);
  modelVector.resize(N * 2 * M);
  initialised = true;
  return mCanvas;
}

void MyCanvas::AddObject(ObjType stat, int Npad, int color, int style, ObjType sys)
{
  if (Npad > int(objVector.size())) {
    printf("Pad no. exceeds number of pads in canvas\n");
    return;
  }
  DrawObject tmp_obj = {stat, sys, color, style, std::visit([](auto&& arg) { return arg->GetMinimum(); }, stat), std::visit([](auto&& arg) { return arg->GetMaximum(); }, stat)};
  objVector[Npad - 1].push_back(tmp_obj);
  return;
}
void MyCanvas::AddModel(ObjType model, int Npad, int color, int fillstyle, float alpha)
{
  if (Npad > int(modelVector.size())) {
    printf("Pad no. exceeds number of pads in canvas\n");
    return;
  }
  ModelObject tmp_obj = {model, color, fillstyle, alpha, std::visit([](auto&& arg) { return arg->GetMinimum(); }, model), std::visit([](auto&& arg) { return arg->GetMaximum(); }, model)};
  modelVector[Npad - 1].push_back(tmp_obj);
  return;
}
void MyCanvas::Draw(const char* opt)
{
  SetLimits();
  ApplyFormat();
  int padcounter = 0;
  for (auto& pad : objVector) {
    ++padcounter;
    int objcounter = 0;
    mCanvas->cd(padcounter);
    /*         std::visit([this,opt,objcounter,padcounter](auto&& arg){
                    arg->Draw("A");
            }, objVector[padcounter][0].stat);  */

    if (!modelVector[padcounter - 1].empty()) {
      int modelcounter = 0;
      for (auto& model : modelVector[padcounter - 1]) {
        ++modelcounter;
        std::visit([model, this](auto&& arg) { arg->SetFillColorAlpha(model.color, model.alpha); }, model.model);
        std::visit([model, this](auto&& arg) { arg->SetLineColorAlpha(model.color, model.alpha); }, model.model);
        if (fModelDrawOpt.find("l") == std::string::npos && fModelDrawOpt.find("L") == std::string::npos)
          std::visit([model, this](auto&& arg) { arg->SetLineWidth(2); }, model.model);
        else
          std::visit([model, this](auto&& arg) { arg->SetLineWidth(2); }, model.model);
        std::visit([model, this](auto&& arg) { arg->SetFillStyle(model.fillstyle); }, model.model);
        std::visit([this, modelcounter, padcounter](auto&& arg) {
          std::string optstr = fModelDrawOpt;
          if (mAlternateAxisSide && ((padcounter - 1) / 2) % 2)
            optstr.insert(0, "Y+");
          if (modelcounter == 1)
            optstr.insert(0, "A");
          arg->Draw(optstr.c_str());
        },
                   model.model);
        if (model.fillstyle == -1) {
          DrawGraphWithErrorOutline(model.model, 2, model.color);
          std::visit([model, this](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, TGraphErrors*> || std::is_same_v<T, TGraphAsymmErrors*>)
              DrawGraphWithErrorOutline(arg, 2, model.color);
          },
                     model.model);
        }
      }
    }
    bool onlyGraphs = true;
    for (auto& obj : pad) {
      std::visit([&onlyGraphs](auto&& arg) {
                if(!isGraphType(arg))
                onlyGraphs = false; },
                 obj.stat);
    }
    for (auto& obj : pad) {
      ++objcounter;
      std::visit([obj, this](auto&& arg) { arg->SetMarkerStyle(obj.style); }, obj.stat);
      // std::visit([obj,this](auto&& arg){ arg->SetMarkerSize(stylemap[obj.style]); }, obj.stat);
      std::visit([obj, this](auto&& arg) { if(stylemap.find(obj.style) == stylemap.end()) arg->SetMarkerStyle(1); else arg->SetMarkerSize(stylemap[obj.style]); }, obj.stat);
      std::visit([obj, this](auto&& arg) { if(arg) arg->SetMarkerStyle(obj.style); }, obj.sys);
      std::visit([obj, this](auto&& arg) { if(arg) { if(stylemap.find(obj.style) == stylemap.end()) arg->SetMarkerStyle(1); else arg->SetMarkerSize(stylemap[obj.style]); } }, obj.sys);

      std::visit([obj, this](auto&& arg) { arg->SetMarkerColor(obj.color); }, obj.stat);
      std::visit([obj, this](auto&& arg) { arg->SetLineColor(obj.color); }, obj.stat);
      std::visit([obj, this](auto&& arg) { if(arg) arg->SetMarkerColor(obj.color); }, obj.sys);
      std::visit([obj, this](auto&& arg) { if(arg) arg->SetLineColor(obj.color); }, obj.sys);

      std::visit([this](auto&& arg) { arg->SetLineWidth(mLineWidth); }, obj.stat);
      std::visit([this](auto&& arg) { if(arg) arg->SetLineWidth(mLineWidth); }, obj.sys);

      std::visit([this](auto&& arg) { if (!mTitle.IsNull() && mTitle.Length() > 0) arg->SetTitle(Form("%s",mTitle.Data())); }, obj.stat);
      std::visit([this, opt, objcounter, padcounter, onlyGraphs](auto&& arg) {
        std::string optstr = opt;
        if (mAlternateAxisSide && ((padcounter - 1) / 2) % 2)
          optstr.insert(0, "Y+");
        // if(objcounter==1) optstr.insert(0,"a");
        if (onlyGraphs)
          optstr.insert(0, "a");
        arg->Draw(optstr.c_str());
      },
                 obj.stat);

      std::visit([this](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, TH1D*>) {
          // Safe to treat arg as TH1D*
          TGraphErrors* gr = HistToGraph(arg);
          gr->SetFillStyle(0);
          gr->Draw("pzsame2");
        } else {
          // Do nothing or handle other types if desired
        }
      },
                 obj.sys);
    }
    gPad->RedrawAxis();
  }
}
void MyCanvas::AddLine(double x1, double y1, double x2, double y2, int nPad)
{
  auto drawLine = [=]() {
    TLine* line = new TLine(x1, y1, x2, y2);
    line->SetBit(kCanDelete);
    line->SetLineStyle(2);
    line->SetLineColorAlpha(kBlack, 0.4);
    line->Draw("same");
  };
  if (nPad) {
    mCanvas->cd(nPad);
    drawLine();
  } else {
    for (int n(0); n < int(objVector.size()); ++n) {
      mCanvas->cd(n + 1);
      drawLine();
    }
  }
  return;
}
void MyCanvas::AddText(double x, double y, const char* text, bool NDC, int textsize, int textcolor)
{
  TLatex l;
  l.SetTextFont(43);
  (!textsize) ? l.SetTextSize(mTextSize) : l.SetTextSize(textsize);
  l.SetTextColor(textcolor);
  if (NDC)
    l.DrawLatexNDC(x, y, text);
  else
    l.DrawLatex(x, y, text);
  return;
}
void MyCanvas::ApplyFormat()
{
  gStyle->SetTickLength(0.9, "x");
  int padCounter = 0;
  for (auto& pad : modelVector) {
    ++padCounter;
    for (auto& model : pad) {
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTickLength(0.9); }, model.model);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTitleFont(43); }, model.model);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTitleSize(mTextSize); }, model.model);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTitleOffset(mOffset[0]); }, model.model);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetLabelFont(43); }, model.model);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetLabelSize(mTextSize); }, model.model);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetTitleFont(43); }, model.model);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetTitleSize(mTextSize); }, model.model);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetTitleOffset(mOffset[1]); }, model.model);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetLabelFont(43); }, model.model);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetLabelSize(mTextSize); }, model.model);
      if (mLabelOffset)
        std::visit([this](auto&& arg) { arg->GetYaxis()->SetLabelOffset(mLabelOffset); }, model.model);
      if (mAdjustTicks) {
        TPad* current_pad = (TPad*)mCanvas->GetPad(padCounter);
        double tickx = 0.01 * mHeight / (current_pad->GetWh() * current_pad->GetAbsHNDC());
        // double ticky = 0.03*mCanvas->GetWw()/(current_pad->GetWw()*current_pad->GetAbsWNDC());
        double padW = current_pad->GetWw() * current_pad->GetAbsWNDC();
        double tickScaleY = (current_pad->GetUymax() - current_pad->GetUymin()) / (current_pad->GetY2() - current_pad->GetY1() - current_pad->GetBottomMargin()) * padW;
        double ticky = 0.00003 * tickScaleY;
        printf("tickx for pad %i = %g\n", padCounter, tickx);
        printf("ticky for pad %i = %g\n", padCounter, ticky);
        std::visit([tickx, this](auto&& arg) { arg->GetXaxis()->SetTickLength(tickx); }, model.model);
        std::visit([ticky, this](auto&& arg) { arg->GetYaxis()->SetTickLength(ticky); }, model.model);
      }
    }
  }
  for (auto& pad : objVector) {
    for (auto& obj : pad) {
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTickLength(0.9); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTitleFont(43); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTitleFont(43); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTitleSize(mTextSize); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTitleOffset(mOffset[0]); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetLabelFont(43); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetLabelSize(mTextSize); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetTitleFont(43); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetTitleSize(mTextSize); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetTitleOffset(mOffset[1]); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetLabelFont(43); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetLabelSize(mTextSize); }, obj.stat);
      if (mLabelOffset)
        std::visit([this](auto&& arg) { arg->GetYaxis()->SetLabelOffset(mLabelOffset); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetXaxis()->SetTickLength(0.02); }, obj.stat);
      std::visit([this](auto&& arg) { arg->GetYaxis()->SetTickLength(0.02); }, obj.stat);
      if (mAdjustTicks) {
        TPad* current_pad = (TPad*)mCanvas->GetPad(padCounter);
        double tickx = 0.01 * mHeight / (current_pad->GetWh() * current_pad->GetAbsHNDC());
        double padW = current_pad->GetWw() * current_pad->GetAbsWNDC();
        double tickScaleY = (current_pad->GetUymax() - current_pad->GetUymin()) / (current_pad->GetY2() - current_pad->GetY1() - current_pad->GetBottomMargin()) * padW;
        double ticky = 0.00003 * tickScaleY;
        std::visit([tickx, this](auto&& arg) { arg->GetXaxis()->SetTickLength(tickx); }, obj.stat);
        std::visit([ticky, this](auto&& arg) { arg->GetYaxis()->SetTickLength(ticky); }, obj.stat);
      }
    }
  }
  return;
}
void MyCanvas::SetLimits()
{
  if (mXmin != -999 && mXmax != -999) {
    for (auto& pad : objVector) {
      for (auto& obj : pad) {
        std::visit([this](auto&& arg) { arg->GetXaxis()->SetRangeUser(mXmin, mXmax); }, obj.stat);
      }
    }
  }
  bool setMin = 1;
  bool setMax = 1;
  if (mYmin != -999) { // If limits are manually set, then set all histograms with those limits
    printf("Setting manual limits!\n");
    setMin = 0;
    for (auto v : objVector) {
      for (auto obj : v) {
        std::visit([this](auto&& arg) { arg->SetMinimum(mYmin); }, obj.stat);
      }
    }
  }
  if (mYmax != -999) {
    setMax = 0;
    for (auto v : objVector) {
      for (auto obj : v) {
        std::visit([this](auto&& arg) { arg->SetMaximum(mYmax); }, obj.stat);
      }
    }
  }
  if (!padlimits.empty()) {
    for (auto& el : padlimits) {
      printf("Setting manual limits (%g,%g) for pad %i!\n", el.second.first, el.second.second, el.first);
      for (auto& model : modelVector[el.first - 1]) {
        std::visit([&el, this](auto&& arg) { arg->SetMinimum(el.second.first); }, model.model);
        std::visit([&el, this](auto&& arg) { arg->SetMaximum(el.second.second); }, model.model);
      }
      for (auto& obj : objVector[el.first - 1]) {
        std::visit([&el, this](auto&& arg) { arg->SetMinimum(el.second.first); }, obj.stat);
        std::visit([&el, this](auto&& arg) { arg->SetMaximum(el.second.second); }, obj.stat);
      }
    }
    setMin = 0;
    setMax = 0;
  }
  if (setMin || setMax) { // Otherwise find max and min for each canvas and set histograms accordingly
    vector<double> max, min;
    for (auto& pad : objVector) {
      double tmp_max = std::visit([](auto&& arg) { return arg->GetMaximum(); }, pad[0].stat);
      double tmp_min = std::visit([](auto&& arg) { return arg->GetMinimum(); }, pad[0].stat);
      for (auto& obj : pad) {
        tmp_max = std::visit([tmp_max](auto&& arg) { if(arg->GetMaximum()>tmp_max) return arg->GetMaximum(); else return tmp_max; }, obj.stat);
        tmp_min = std::visit([tmp_min](auto&& arg) { if(arg->GetMinimum()<tmp_min) return arg->GetMinimum(); else return tmp_min; }, obj.stat);
      }
      max.push_back(tmp_max);
      min.push_back(tmp_min);
    }
    double global_max = *std::max_element(max.begin(), max.end()); // If pads should share common y scale, then find the global max and min
    double global_min = *std::min_element(min.begin(), min.end());
    int padcounter = 0;
    for (auto& pad : objVector) {
      ++padcounter;
      double scale = (padcounter % 2) ? mRatioYScale : mYscale;
      for (auto& obj : pad) {
        if (mUseSameYscale) {
          double diff = global_max - global_min;
          if (setMax)
            std::visit([scale, diff, global_max](auto&& arg) { arg->SetMaximum(global_max + scale * diff); }, obj.stat);
          if (setMin)
            std::visit([scale, diff, global_min](auto&& arg) { arg->SetMinimum(global_min - scale * diff); }, obj.stat);
        } else {
          double diff = max[padcounter - 1] - min[padcounter - 1];
          if (setMax)
            std::visit([scale, diff, min, max, padcounter](auto&& arg) { arg->SetMaximum(max[padcounter - 1] + scale * diff); }, obj.stat);
          if (setMin)
            std::visit([scale, diff, min, max, padcounter](auto&& arg) { arg->SetMinimum(min[padcounter - 1] - scale * diff); }, obj.stat);
        }
      }
    }
  }
  if (fRatioCanvas && mYRatioMax != -999 && mYRatioMin != -999) {
    int padcounter = 0;
    for (auto& pad : objVector) {
      ++padcounter;
      if (!(padcounter % 2)) {
        for (auto& obj : pad) {
          std::visit([this](auto&& arg) { arg->SetMinimum(mYRatioMin); }, obj.stat);
          std::visit([this](auto&& arg) { arg->SetMaximum(mYRatioMax); }, obj.stat);
        }
      }
    }
  }
  return;
}
void MyCanvas::cd(int Npad)
{
  (Npad == 0) ? mCanvas->cd() : mCanvas->cd(Npad);
  return;
}
TPad* MyCanvas::GetPad(int Npad)
{
  if (!mCanvas) {
    return nullptr;
  }
  if (Npad <= 0) {
    return dynamic_cast<TPad*>(mCanvas);
  }
  return dynamic_cast<TPad*>(mCanvas->GetPad(Npad));
}
void MyCanvas::SaveAs(const char* name, const char* opt)
{
  mCanvas->SaveAs(name, opt);
  return;
}
void MyCanvas::BuildLegend(double x1, double y1, double x2, double y2, int legInd, int nCol)
{
  if (static_cast<size_t>(legInd + 1) > mLeg.size())
    mLeg.push_back(new TLegend(x1, y1, x2, y2));
  mLeg[legInd]->SetBorderSize(0);
  mLeg[legInd]->SetFillStyle(0);
  mLeg[legInd]->SetTextFont(43);
  mLeg[legInd]->SetTextSize(mTextSize);
  if (nCol > 1)
    mLeg[legInd]->SetNColumns(nCol);
  mLeg[legInd]->SetTextSize(mTextSize);
  return;
}
void MyCanvas::DrawLegend(double x1, double y1, int legInd, int nCol)
{
  mLeg[legInd]->Draw();
  gPad->Update();
  mLeg[legInd]->SetX1NDC(x1);
  mLeg[legInd]->SetY1NDC(y1);
  mLeg[legInd]->SetX2NDC(x1 + 0.2);
  mLeg[legInd]->SetY2NDC(y1 + 0.3);
  if (nCol > 1)
    mLeg[legInd]->SetNColumns(nCol);
  mLeg[legInd]->SetTextSize(mTextSize);
  gPad->Modified();
  return;
}
void MyCanvas::DrawLegend(double x1, double y1, double x2, double y2, int legInd, int nCol)
{
  mLeg[legInd]->Draw();
  gPad->Update();
  mLeg[legInd]->SetX1NDC(x1);
  mLeg[legInd]->SetY1NDC(y1);
  mLeg[legInd]->SetX2NDC(x2);
  mLeg[legInd]->SetY2NDC(y2);
  if (nCol > 1)
    mLeg[legInd]->SetNColumns(nCol);
  mLeg[legInd]->SetTextSize(mTextSize);
  gPad->Modified();
  return;
}
void MyCanvas::DrawLegend(TString loc, int legInd, int nCol)
{
  mLeg[legInd]->Draw();
  double x1 = 0.0;
  double y1 = 0.0;
  if (loc.Contains("top"))
    y1 = (addAbove) ? 0.95 - gPad->GetTopMargin() - 0.05 : 0.95 - gPad->GetTopMargin();
  if (loc.Contains("left"))
    x1 = gPad->GetLeftMargin() + 0.05;
  if (loc.Contains("right"))
    x1 = 0.7;
  if (loc.Contains("bottom"))
    y1 = (addAbove) ? gPad->GetBottomMargin() + 0.05 : gPad->GetBottomMargin() + 0.1;
  if (loc.Contains("mid"))
    y1 = 0.5;
  if (loc.Contains("center"))
    x1 = 0.5;
  if (x1 == 0.0 || y1 == 0.0) {
    printf("Specify location by combining (bottom,mid,top) and (left,center,right)\n");
    return;
  }
  int nrows = mLeg[legInd]->GetNRows();
  gPad->Update();
  mLeg[legInd]->SetX1NDC(x1);
  mLeg[legInd]->SetY1NDC((addAbove) ? y1 : y1 - 0.05 * nrows);
  mLeg[legInd]->SetX2NDC(x1 + 0.2);
  mLeg[legInd]->SetY2NDC((addAbove) ? y1 + 0.05 * nrows : y1);
  mLeg[legInd]->SetNColumns(nCol);
  mLeg[legInd]->SetTextSize(mTextSize);
  gPad->Modified();
  return;
}
template <typename T>
TGraphErrors* MyCanvas::HistToGraph(T h)
{
  int nBins = h->GetNbinsX();
  vector<double> binEdges;
  for (auto i(0); i <= nBins; ++i)
    binEdges.push_back(h->GetBinLowEdge(i + 1));
  vector<double> val;
  vector<double> ey;
  vector<double> ex;
  vector<double> x;
  TGraphErrors* gr;
  for (Int_t bin(0); bin < nBins; ++bin) {
    val.push_back(h->GetBinContent(bin + 1));
    ey.push_back(h->GetBinError(bin + 1));
    // ex.push_back(h->GetBinWidth(bin+1)*0.18);
    ex.push_back(0.05);
    x.push_back((binEdges[bin] + binEdges[bin + 1]) / 2);
  }
  gr = new TGraphErrors(nBins, &x[0], &val[0], &ex[0], &ey[0]);
  gr->SetFillStyle(1001);
  gr->SetLineWidth(h->GetLineWidth());
  gr->SetFillColorAlpha(h->GetMarkerColor(), 0.4);
  gr->SetLineColor(h->GetLineColor());
  gr->SetMarkerColor(h->GetMarkerColor());
  gr->SetMarkerStyle(h->GetMarkerStyle());
  gr->SetMarkerSize(h->GetMarkerSize());
  gr->SetMinimum(h->GetMinimum());
  gr->SetMinimum(h->GetMaximum());
  return gr;
}
float MyCanvas::SafeAdd(float x, float eps)
{
  return std::min(1.0f, x + eps);
}

float MyCanvas::SafeSub(float x, float eps)
{
  return std::max(0.0f, x - eps);
}
void MyCanvas::ApplyPadStyle(TPad* pad)
{
  pad->SetFillColor(gStyle->GetPadColor());
  // pad->SetTickx(gStyle->GetPadTickX());
  // pad->SetTicky(gStyle->GetPadTickY());
  pad->SetFrameLineWidth(gStyle->GetFrameLineWidth());
  pad->SetFrameFillColor(gStyle->GetFrameFillColor());
  pad->SetFrameBorderMode(gStyle->GetFrameBorderMode());
  pad->SetFrameBorderSize(gStyle->GetFrameBorderSize());
}
void MyCanvas::DrawGraphWithErrorOutline(ObjType gvar, int lineWidth = 2, int lineColor = kBlack)
{
  std::visit([lineWidth, lineColor](auto* g) {
    using T = std::decay_t<decltype(g)>;
    if constexpr (std::is_same_v<T, TGraphErrors*> || std::is_same_v<T, TGraphAsymmErrors*>) {
      if (!g)
        return;

      // Draw the error band with transparent fill
      g->SetFillStyle(-1); // No fill
      // g->SetFillColor(0);
      g->Draw("3");

      // Prepare vectors
      int n = g->GetN();
      std::vector<double> x(n), yUpper(n), yLower(n);

      for (int i = 0; i < n; ++i) {
        double xi, yi;
        g->GetPoint(i, xi, yi);
        x[i] = xi;

        if constexpr (std::is_same_v<std::decay_t<decltype(*g)>, TGraphErrors>) {
          double ey = g->GetErrorY(i);
          yUpper[i] = yi + ey;
          yLower[i] = yi - ey;
        } else if constexpr (std::is_same_v<std::decay_t<decltype(*g)>, TGraphAsymmErrors>) {
          double eyh = g->GetErrorYhigh(i);
          double eyl = g->GetErrorYlow(i);
          yUpper[i] = yi + eyh;
          yLower[i] = yi - eyl;
        }
      }

      // Create and draw upper and lower outlines
      TGraph* upper = new TGraph(n, x.data(), yUpper.data());
      TGraph* lower = new TGraph(n, x.data(), yLower.data());
      upper->SetLineWidth(lineWidth);
      lower->SetLineWidth(lineWidth);
      upper->SetLineColor(lineColor);
      lower->SetLineColor(lineColor);
      upper->SetLineStyle(1); // dashed
      lower->SetLineStyle(1);
      upper->Draw("L SAME");
      lower->Draw("L SAME");
    }
  },
             gvar);
}
// NOLINTEND(readability/casting, readability/braces, whitespace/parens, whitespace/newline, whitespace/comma, whitespace/operators)
