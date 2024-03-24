// Copyright 2019-2099 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Rivet/Analysis.hh"
#include "Rivet/Config/RivetConfig.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet
{
/** @brief Example analysis */
class ALICE_YYYY_I1234567 : public Analysis
{
 public:
  typedef ALICE::PrimaryParticles PrimProj;

  /** Constructor */
  ALICE_YYYY_I1234567() : Analysis("ALICE_YYYY_I1234567") {}

  /**
   * @name Analysis methods
   * @{
   */
  /** Book histograms and initialise projections before the run */
  void init()
  {
    // Initialise and register projections
    declare(PrimProj(Cuts::abseta < 5 && Cuts::abscharge3 > 0), "APRIM");
    // book histograms
    book(_h, 1, 1, 1);
  }
  /** Perform the per-event analysis */
  void analyze(const Event& event)
  {
    for (auto p : apply<PrimProj>(event, "APRIM").particles())
      _h->fill(p.eta());
  }
  /** Normalise histograms etc., after the run */
  void finalize()
  {
    std::cout << "Scaling histogram by " << numEvents()
              << " " << sumW() << std::endl;
    scale(_h, 1. / sumW());
  }
  /** @} */

  /**
   * @name Histograms
   * @{
   */
  Histo1DPtr _h;
  /* @} */
};
RIVET_DECLARE_PLUGIN(ALICE_YYYY_I1234567);
} // namespace Rivet

//
// EOF
//
