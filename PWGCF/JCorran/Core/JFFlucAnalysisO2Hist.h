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
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \since Sep 2024

#ifndef PWGCF_JCORRAN_CORE_JFFLUCANALYSISO2HIST_H_
#define PWGCF_JCORRAN_CORE_JFFLUCANALYSISO2HIST_H_

#include "JFFlucAnalysis.h"

#include "Framework/HistogramRegistry.h"

class JFFlucAnalysisO2Hist : public JFFlucAnalysis
{
 public:
  JFFlucAnalysisO2Hist(o2::framework::HistogramRegistry&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, o2::framework::AxisSpec&, uint16_t, const TString&);
  ~JFFlucAnalysisO2Hist();
};

#endif // PWGCF_JCORRAN_CORE_JFFLUCANALYSISO2HIST_H_
