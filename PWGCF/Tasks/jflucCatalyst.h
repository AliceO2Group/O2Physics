#ifndef JFLUC_CATALYST_H
#define JFLUC_CATALYST_H

namespace o2::aod{
namespace particleTrack{
DECLARE_SOA_INDEX_COLUMN(Collision,collision);
DECLARE_SOA_COLUMN(Pt,pt,float);
DECLARE_SOA_COLUMN(Eta,eta,float);
DECLARE_SOA_COLUMN(Phi,phi,float);
DECLARE_SOA_COLUMN(WeightNUA,weightNUA,float);
DECLARE_SOA_COLUMN(WeightEff,weightEff,float);
}

namespace collisionData{
DECLARE_SOA_INDEX_COLUMN(Collision,collision);
DECLARE_SOA_COLUMN(Cent,cent,float);
}

DECLARE_SOA_TABLE(ParticleTrack,"AOD","PARTICLETRACK",
	particleTrack::CollisionId,
	particleTrack::Pt,particleTrack::Eta,particleTrack::Phi,
	particleTrack::WeightNUA,particleTrack::WeightEff);
DECLARE_SOA_TABLE(CollisionData,"AOD","COLLISIONDATA",
	collisionData::CollisionId,
	collisionData::Cent);
}

#endif

