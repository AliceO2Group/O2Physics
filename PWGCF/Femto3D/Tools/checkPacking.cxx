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

///
/// \file   checkPacking.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  exec to check the packing and depacking of femto tables
/// \since  03/05/2024
///

#include "PWGCF/Femto3D/DataModel/singletrackselector.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom.h"

using namespace o2;

template <typename T>
bool process(std::string outputName, int nevents = 100000)
{
  class Container
  {
   public:
    Container() {}
    void operator()(const float toPack) { mPacked = aod::singletrackselector::packInTable<T>(toPack); }
    void test(const float toPack)
    {
      auto bin = aod::singletrackselector::packInTable<T>(toPack);
      LOG(info) << toPack << " goes to " << aod::singletrackselector::unPack<T>(bin) << " bin " << static_cast<int>(bin);
    }
    T::binned_t mPacked = 0;
    float unpack() { return aod::singletrackselector::unPack<T>(mPacked); }
  } container;

  T::print();
  const float min = T::binned_min;
  const float max = T::binned_max;
  container.test(0);
  container.test(0 + T::bin_width);
  container.test(0 - T::bin_width);
  container.test(0 - T::bin_width * 0.5);
  const int nbins = (max - min) / T::bin_width;
  std::vector<float> xbins;
  for (int i = 0; i <= nbins; i++) {
    const float x = min + i * T::bin_width;
    const auto ix = aod::singletrackselector::packInTable<T>(x);
    const float u = aod::singletrackselector::unPack<T>(ix);
    LOG(info) << "Bin " << i << "/" << xbins.size() << " " << x << " => " << static_cast<int>(ix) << " " << u;
    if (i > 1) {
      if (ix == aod::singletrackselector::packInTable<T>(xbins.back())) {
        continue;
      }
    }
    xbins.push_back(u);
  }
  LOG(info) << "Min = " << min << " Max = " << max;
  TH1F* hgaus = new TH1F("hgaus", "", nbins, min+T::bin_width*0.5, max+0.5*T::bin_width);
  hgaus->Print();
  LOG(info) << "Bin width = " << T::bin_width << " vs histo " << hgaus->GetXaxis()->GetBinWidth(1);
  // for (int i = 1; i <= hgaus->GetNbinsX(); i++) {
  //   // Check the packing of the bins center
  //   const float x = hgaus->GetXaxis()->GetBinCenter(i);
  //   const auto ix = aod::singletrackselector::packInTable<T>(x);
  //   if (i != ix) {
  //     LOG(info) << "Bin " << i << " " << hgaus->GetBinLowEdge(i) << " - " << hgaus->GetBinLowEdge(i + 1) << " is mapped wrongly to " << ix;
  //   }
  // }
  hgaus->SetLineColor(2);
  hgaus->SetLineStyle(1);
  TH1F* hgausPacked = static_cast<TH1F*>(hgaus->Clone("hgausPacked"));
  hgausPacked->SetLineColor(4);
  hgausPacked->SetLineStyle(2);

  TH1F* huniform = static_cast<TH1F*>(hgaus->Clone("huniform"));
  huniform->SetLineColor(2);
  huniform->SetLineStyle(1);
  TH1F* huniformPacked = static_cast<TH1F*>(hgaus->Clone("huniformPacked"));
  huniformPacked->SetLineColor(4);
  huniformPacked->SetLineStyle(2);

  for (int i = 0; i < nevents; i++) {
    float nsigma = gRandom->Gaus(0, 1);
    hgaus->Fill(nsigma);
    container(nsigma);
    hgausPacked->Fill(container.unpack());

    nsigma = gRandom->Uniform(-10, 10);
    huniform->Fill(nsigma);
    container(nsigma);
    huniformPacked->Fill(container.unpack());
  }

  TCanvas* can = new TCanvas("can");
  hgaus->Draw();
  hgausPacked->Draw("same");
  outputName = "/tmp/" + outputName + ".pdf";
  can->SaveAs("/tmp/asd.root");
  can->SaveAs(Form("%s[", outputName.c_str()));
  can->SaveAs(outputName.c_str());

  huniform->Draw();
  huniformPacked->Draw("same");
  can->SaveAs(outputName.c_str());
  can->SaveAs(Form("%s]", outputName.c_str()));
  const bool gausOk = (hgaus->GetBinContent(hgaus->FindBin(0)) == hgausPacked->GetBinContent(hgausPacked->FindBin(0)));
  if (!gausOk) {
    LOG(info) << "Gaus packing/unpacking failed";
  }
  const bool uniformOk = (huniform->GetBinContent(huniform->FindBin(0)) == huniformPacked->GetBinContent(huniformPacked->FindBin(0)));
  if (!uniformOk) {
    LOG(info) << "Uniform packing/unpacking failed";
  }
  return gausOk && uniformOk;
}

int main(int /*argc*/, char* /*argv*/[])
{

  LOG(info) << "Checking the packing and unpacking of PID signals (nsigmas) in the Femto PID response.";
  if (process<o2::aod::singletrackselector::binning::nsigma>("Nsigma", 100000)) {
    LOG(info) << "Packing and unpacking of PID signals (nsigmas) in the Femto PID response is correct.";
  } else {
    LOG(fatal) << "Packing and unpacking of PID signals (nsigmas) in the Femto PID response is incorrect.";
  }

} // main
