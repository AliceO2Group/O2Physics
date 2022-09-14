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
/// \brief
/// \author Roman Lavicka, roman.lavicka@cern.ch
/// \since  12.07.2022

#ifndef O2_ANALYSIS_UPCFILTERCENTRALBARREL_H
#define O2_ANALYSIS_UPCFILTERCENTRALBARREL_H

enum MyParticle { P_ELECTRON = 0, P_MUON = 1, P_PION = 2, P_KAON = 3, P_PROTON = 4};

namespace o2::aod
{

	namespace upccentralbarrel
	{
		DECLARE_SOA_COLUMN(IsTwoTracks, isTwoTracks, bool); //! Exactly two tracks were found
		DECLARE_SOA_COLUMN(IsFourTracks, isFourTracks, bool); //! Exactly four tracks were found
		DECLARE_SOA_COLUMN(FlagWhichParticle, flagWhichParticle, int); //! Flag the particle according to TPC and TOF
		DECLARE_SOA_COLUMN(ReachedTOF, reachedTOF, bool);  //! Mark, if track hits TOF
	} // namespace upccentralbarrel

	DECLARE_SOA_TABLE(UPCTrackCandidates, "AOD", "UPCCANDIDATES", //! Table with UPC track candidates
										o2::soa::Index<>,
										collision::PosX, collision::PosY, collision::PosZ, collision::NumContrib,
										track::Pt, track::P,
										upccentralbarrel::IsTwoTracks, upccentralbarrel::IsFourTracks,
										upccentralbarrel::FlagWhichParticle, upccentralbarrel::ReachedTOF);
	using UPCTrackCandidate = UPCTrackCandidates::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_UPCFILTERCENTRALBARREL_H
