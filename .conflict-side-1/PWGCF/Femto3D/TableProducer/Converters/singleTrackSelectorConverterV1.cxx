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
//
/// \brief Converter for the different versions of the singletrackselector tables
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 03 May 2024

#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <fairlogger/Logger.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::aod;
//::singletrackselector; // the namespace defined in .h

struct singleTrackSelectorConverter {
  Produces<o2::aod::SingleTrackSels_v2> tableRow;

  void init(InitContext&) {}

  void process(o2::aod::SingleTrackSels_v1 const& tracks)
  {
    tableRow.reserve(tracks.size());
    for (auto const& track : tracks) {
      tableRow(track.singleCollSelId(),
               track.p(),
               track.eta(),
               track.phi(),
               track.sign(),
               track.tpcNClsFound(),
               track.tpcNClsShared(),
               track.itsClsMap(),
               track.itsClusterSizes(),
               singletrackselector::packSymmetric<singletrackselector::binning::dca>(track.dcaXY()),
               singletrackselector::packSymmetric<singletrackselector::binning::dca>(track.dcaZ()),
               singletrackselector::packInTable<singletrackselector::binning::chi2>(track.tpcChi2NCl()),
               singletrackselector::packInTable<singletrackselector::binning::chi2>(track.itsChi2NCl()),
               singletrackselector::packInTable<singletrackselector::binning::rowsOverFindable>(track.tpcCrossedRowsOverFindableCls()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaPi()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaPi()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaKa()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaKa()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaPr()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaPr()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaDe()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaDe()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaHe()),
               singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaHe()));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<singleTrackSelectorConverter>(cfgc)};
}
