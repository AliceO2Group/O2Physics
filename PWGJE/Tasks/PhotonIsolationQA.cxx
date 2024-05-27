#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "CommonDataFormat/InteractionRecord.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//using tracksPID = o2::soa::Join<o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::FullTracks, o2::aod::TrackSelection>;
using myGlobTracks = o2::soa::Join<o2::aod::FullTracks, o2::aod::TrackSelection>;
using collisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>; //monte carlo change
using MCcollisions = o2::aod::McCollisions;// mc
using mcCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>; //mc
using filteredMCCells = o2::soa::Filtered<mcCells>; //mc
using MCClusters = o2::soa::Join<o2::aod::EMCALMCClusters, o2::aod::EMCALClusters>; //mc

struct Photon_Isolation_QA{
    HistogramRegistry Cluster_Info{"Cluster Info", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
    HistogramRegistry Matched_Track_Info{"Track Info", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
    HistogramRegistry Photon_Isolation{"Photon_Isolation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
    HistogramRegistry MC_Info{"MC_Info", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
   
    using o2HistType = HistType;
    using o2Axis = AxisSpec;

    Configurable<float> maxPosZ{"maxPosZ", 10.0f, "maximum z position of collision in cm"};
    Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
    Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
    Configurable<int> mClusterDefinition{"clusterDefinition", 0, "cluster definition to be selected, e.g. 10=kV3Default, 0 = kV1Default"};
    Configurable<float> minTime{"minTime", -30., "Minimum cluster time for time cut"};
    Configurable<float> maxTime{"maxTime", +35., "Maximum cluster time for time cut"};
    Configurable<float> minClusterEnergy{"minClusterEnergy", 0.7f, "Minimal cluster energy"};
    Configurable<int> minNCells{"minNCelss", 2, "Minimal amount of cells per cluster"};
    Configurable<int> maxNLM{"maxNLM", 2, "Maximal amount of local maxima per cluster"};
    Configurable<bool> ExoticContribution{"ExoticContribution", false, "Exotic cluster in the data"};
    Configurable<float> TrackMatchingRadius{"TrackMatchingRadius", 0.05, "maximum radius between track and cluster in eta phi plane"};
    Configurable<float> minDPhi{"minDPhi", 0.01, "Minimum dPhi between track and cluster"};

    Filter MCPosZFilter = nabs(aod::mccollision::posZ) < maxPosZ;  //monte change
    Filter PosZFilter = nabs(aod::collision::posZ) < maxPosZ; 
    Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::energy > minClusterEnergy) && (o2::aod::emcalcluster::nCells >= minNCells) && (o2::aod::emcalcluster::nlm <= maxNLM) && (o2::aod::emcalcluster::isExotic == ExoticContribution);
    Filter emccellfilter = aod::calo::caloType == 1; //mc emcal cell

    using selectedCollisions = soa::Filtered<collisionEvSelIt>;
    using selectedMCcollisions = soa::Filtered<MCcollisions>;
    using selectedClusters = soa::Filtered<aod::EMCALClusters>;
    using selectedMCClusters = soa::Filtered<MCClusters>;

    Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;

    int eventSelection = -1;
    int trackSelection = -1;
    void init(o2::framework::InitContext&){
        eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
        trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

        const o2Axis PosZ_Axis{100, -15, 15, "Z Positions (cm)"};
        const o2Axis Num_Cluster_Axis{10, 0, 10, "N Clusters"};
        const o2Axis Shower_Shape_Long_Axis{100, 0, 4, "#sigma^{2}_{long}"};
        const o2Axis Shower_Shape_Short_Axis{100, 0, 4, "#sigma^{2}_{short}"};
        const o2Axis Phi_Axis{100, 1, 6, "#phi"};
        const o2Axis Eta_Axis{100, -0.9, 0.9, "#eta"};
        const o2Axis Energy_Axis{100, 0, 100, "E (GeV/c)"};
        const o2Axis NCells_Axis{10,0, 10, "N Cells"};
        const o2Axis Pt_Axis{100, 0, 100, "P_{T} (GeV/c)"};
        const o2Axis P_Axis{100, 0, 100, "P (GeV/c)"};
        const o2Axis RatioP_Axis{50, 0, 10, "E/P"};
        const o2Axis Num_Track_Axis{10, 0, 10, "N Tracks"};
        const o2Axis PtIso_Axis{100, -10, 15, "P_{T, Iso} (GeV/c)"};
        const o2Axis Rho_Axis{100, 0,100, "#rho (#frac{GeV/c}{A})"};

        // Define histograms for Z-possition collision and cluster information
        Cluster_Info.add("hPosZ", "Z Position of collision", o2HistType::kTH1F, {PosZ_Axis});
        Cluster_Info.add("hNumClusters", "Number of cluster per collision", o2HistType::kTH1F, {Num_Cluster_Axis});
        Cluster_Info.add("hClusterLocation", "Location of shower in eta phi plane", o2HistType::kTH2F, {{Eta_Axis},{Phi_Axis}});
        Cluster_Info.add("hEnergy", "Energy of the cluster", o2HistType::kTH1F, {Energy_Axis});
        Cluster_Info.add("hEnergy_ShowerShapeLong", "Energy vs Shower shape long axis", o2HistType::kTH2F, {{Energy_Axis},{Shower_Shape_Long_Axis}});
        Cluster_Info.add("hEnergy_ShowerShapeShort", "Energy vs Shower shape short axis", o2HistType::kTH2F, {{Energy_Axis},{Shower_Shape_Short_Axis}});
        Cluster_Info.add("hEnergy_NCells", "Energy vs Number of cells in cluster", o2HistType::kTH2F, {{Energy_Axis},{NCells_Axis}});
        Cluster_Info.add("hEnergyvsSigmaLong", "Energy cluster vs Long shower shape", o2HistType::kTH2F, {{Energy_Axis},{Shower_Shape_Long_Axis}});

        //Define histograms for tracks matched to cluster
        Matched_Track_Info.add("hNumTracks", "Number of matched tracks per cluster", o2HistType::kTH1F, {Num_Track_Axis});
        Matched_Track_Info.add("hEvsNumTracks", "Energy of cluster vs matched tracks", o2HistType::kTH2F, {{Energy_Axis}, {Num_Track_Axis}});
        Matched_Track_Info.add("hTrackPt", "Pt of matched track", o2HistType::kTH1F, {Pt_Axis});
        Matched_Track_Info.add("hTrackEta", "Eta of matched track", o2HistType::kTH1F, {Eta_Axis});
        Matched_Track_Info.add("hTrackPhi", "Phi of matched track", o2HistType::kTH1F, {Phi_Axis});
        Matched_Track_Info.add("hRatioClusterETrackP", "Ratio between energy of cluster and P of tracks", o2HistType::kTH1F, {RatioP_Axis});
        Matched_Track_Info.add("hRatioClusterETrackPrackPt", "Ratio between E and P vs Pt of matched tracks", o2HistType::kTH2F, {{RatioP_Axis}, {Pt_Axis}});

        //Define histograms for isolated photon candidates
        Photon_Isolation.add("hPtIso", "Pt_Iso", o2HistType::kTH1F, {PtIso_Axis});
        Photon_Isolation.add("hRho_Perpen_Cone", "Density of perpendicular cone", o2HistType::kTH1F, {Rho_Axis});
        Photon_Isolation.add("hShowerShape", "Shower shape", o2HistType::kTH2F, {{Shower_Shape_Long_Axis},{Shower_Shape_Short_Axis}});
        Photon_Isolation.add("hPtTracksOfCluster", "Tracks matched to the cluster", o2HistType::kTH1F, {Pt_Axis});
        Photon_Isolation.add("hSigmaLongvsPtIso", "Long shower shape vs Pt_Iso", o2HistType::kTH2F, {{Shower_Shape_Long_Axis},{PtIso_Axis}});
        Photon_Isolation.add("hABCDControlRegion", "Yield Control Regions", o2HistType::kTH1F, {{4, 0.0, 4.0}}); 
        Photon_Isolation.add("hABCDControlRegionE_10", "Yield Control Regions for E cluster > 10 GeV", o2HistType::kTH1F, {{4, 0.0, 4.0}}); 
        Photon_Isolation.add("hABCDControlRegionE_20", "Yield Control Regions for E cluster > 20 GeV", o2HistType::kTH1F, {{4, 0.0, 4.0}}); 
        Photon_Isolation.add("hABCDControlRegionE_30", "Yield Control Regions for E cluster > 30 GeV", o2HistType::kTH1F, {{4, 0.0, 4.0}}); 
        Photon_Isolation.add("hABCDControlRegionE_40", "Yield Control Regions for E cluster > 40 GeV", o2HistType::kTH1F, {{4, 0.0, 4.0}}); 
        Photon_Isolation.add("hABCDControlRegionE_50", "Yield Control Regions for E cluster > 50 GeV", o2HistType::kTH1F, {{4, 0.0, 4.0}}); 
        
        std::vector<std::string> bin_names = {"A", "B", "C", "D"};
        for (int i = 0; i < bin_names.size(); i++){
            Photon_Isolation.get<TH1>(HIST("hABCDControlRegion"))->GetXaxis()->SetBinLabel(i + 1, bin_names[i].c_str());
            Photon_Isolation.get<TH1>(HIST("hABCDControlRegionE_10"))->GetXaxis()->SetBinLabel(i + 1, bin_names[i].c_str());
            Photon_Isolation.get<TH1>(HIST("hABCDControlRegionE_20"))->GetXaxis()->SetBinLabel(i + 1, bin_names[i].c_str());
            Photon_Isolation.get<TH1>(HIST("hABCDControlRegionE_30"))->GetXaxis()->SetBinLabel(i + 1, bin_names[i].c_str());
            Photon_Isolation.get<TH1>(HIST("hABCDControlRegionE_40"))->GetXaxis()->SetBinLabel(i + 1, bin_names[i].c_str());
            Photon_Isolation.get<TH1>(HIST("hABCDControlRegionE_50"))->GetXaxis()->SetBinLabel(i + 1, bin_names[i].c_str());
        }

        MC_Info.add("hPosZ", "Z Position of collision", o2HistType::kTH1F, {PosZ_Axis});
        MC_Info.add("hMcParticlesToCluster", "Energy of particles linked to mc cluster", o2HistType::kTH1F, {{100, 0, 20}});
        MC_Info.add("hMcParticleStatusCode", "status code of mc particle linked to mc cluster", o2HistType::kTH1F, {{100, 0, 100}});
        MC_Info.add("hECluster", "Energy of MC cluster", o2HistType::kTH1F, {Energy_Axis});
        MC_Info.add("hSigmaLongvsPtIso", "Long shower shape vs Pt_Iso", o2HistType::kTH2F, {{Shower_Shape_Long_Axis},{PtIso_Axis}});
        MC_Info.add("hPtIso", "Pt_Iso", o2HistType::kTH1F, {PtIso_Axis});
        MC_Info.add("hRho_Perpen_Cone", "Density of perpendicular cone", o2HistType::kTH1F, {Rho_Axis});
        MC_Info.add("hShowerShape", "Shower shape", o2HistType::kTH2F, {{Shower_Shape_Long_Axis},{Shower_Shape_Short_Axis}});

        MC_Info.add("hABCDControlRegion", "Yield Control Regions", o2HistType::kTH1F, {{5, 0.0, 5.0}}); 
        MC_Info.add("hABCDControlRegionE_10", "Yield Control Regions for E cluster > 10 GeV", o2HistType::kTH1F, {{5, 0.0, 5.0}}); 
        MC_Info.add("hABCDControlRegionE_20", "Yield Control Regions for E cluster > 20 GeV", o2HistType::kTH1F, {{5, 0.0, 5.0}}); 
        MC_Info.add("hABCDControlRegionE_30", "Yield Control Regions for E cluster > 30 GeV", o2HistType::kTH1F, {{5, 0.0, 5.0}}); 
        MC_Info.add("hABCDControlRegionE_40", "Yield Control Regions for E cluster > 40 GeV", o2HistType::kTH1F, {{5, 0.0, 5.0}}); 
        MC_Info.add("hABCDControlRegionE_50", "Yield Control Regions for E cluster > 50 GeV", o2HistType::kTH1F, {{5, 0.0, 5.0}}); 
        
        std::vector<std::string> bin_names_mc = {"A", "B", "C", "D", "True Bckgr A"};
        for (int i = 0; i < bin_names_mc.size(); i++){
            MC_Info.get<TH1>(HIST("hABCDControlRegion"))->GetXaxis()->SetBinLabel(i + 1, bin_names_mc[i].c_str());
            MC_Info.get<TH1>(HIST("hABCDControlRegionE_10"))->GetXaxis()->SetBinLabel(i + 1, bin_names_mc[i].c_str());
            MC_Info.get<TH1>(HIST("hABCDControlRegionE_20"))->GetXaxis()->SetBinLabel(i + 1, bin_names_mc[i].c_str());
            MC_Info.get<TH1>(HIST("hABCDControlRegionE_30"))->GetXaxis()->SetBinLabel(i + 1, bin_names_mc[i].c_str());
            MC_Info.get<TH1>(HIST("hABCDControlRegionE_40"))->GetXaxis()->SetBinLabel(i + 1, bin_names_mc[i].c_str());
            MC_Info.get<TH1>(HIST("hABCDControlRegionE_50"))->GetXaxis()->SetBinLabel(i + 1, bin_names_mc[i].c_str());
        }
    }
        
    // boolian returns true if a track is matched to the cluster with E_Cluster/P_track < 1.75. Otherwise returns false
    bool track_matching(const auto& cluster, o2::aod::EMCALMatchedTracks const& matched_tracks){
        for (const auto& match : matched_tracks){
            double dEta = abs(match.track_as<myGlobTracks>().trackEtaEmcal() - cluster.eta()); 
            double dPhi = abs(match.track_as<myGlobTracks>().trackPhiEmcal() - cluster.phi());
            if (dPhi > M_PI) dPhi = 2.*M_PI - dPhi;
            double distance = sqrt(pow(dEta, 2) + pow(dPhi,2));
            if (distance <= TrackMatchingRadius) { 
                double abs_p = abs(match.track_as<myGlobTracks>().p());
                if ((cluster.energy()/abs_p) < 1.75){
                    return true;
                    } 
                }
            }
        return false;
    }

    // sums the pt of all tracks around the cluster within a radius of R = 0.4
    double sum_Pt_tracks_in_cone(const auto& cluster, myGlobTracks const& tracks, double Cone_Radius = 0.4){
        double sum_Pt = 0.;
        for (const auto& track : tracks){
            double dphi = abs(cluster.phi() - track.phi());
            if (dphi > M_PI) dphi = 2.*M_PI - dphi;
            double distance = sqrt(pow((cluster.eta() - track.eta()), 2) + pow(dphi,2));
            if (distance < Cone_Radius){
                sum_Pt += track.pt();
            }
        }
        return sum_Pt;
    }
    
    //Calculates the Pt density of tracks of the UE with the perpendicular cone method
    double Rho_Perpendicular_Cone(const auto& cluster, myGlobTracks const& tracks, double Perpendicular_Cone_Radius = 0.4){
        double sum_Pt = 0.;
        double Perpendicular_Phi1 = cluster.phi() + (1./2.)*M_PI;
        double Perpendicular_Phi2 = cluster.phi() - (1./2.)*M_PI;
        if (Perpendicular_Phi1 > 2*M_PI) Perpendicular_Phi1 = Perpendicular_Phi1 - 2*M_PI;
        if (Perpendicular_Phi2 < 0.) Perpendicular_Phi2 = 2*M_PI + Perpendicular_Phi2;
        for (const auto& track : tracks){
            double dphi1 = abs(Perpendicular_Phi1 - track.phi());
            double dphi2 = abs(Perpendicular_Phi2 - track.phi());
            if (dphi1 > M_PI) dphi1 = 2.*M_PI - dphi1;
            if (dphi2 > M_PI) dphi2 = 2.*M_PI - dphi2;
            double distance1 = sqrt(pow((cluster.eta() - track.eta()), 2) + pow(dphi1,2));
            double distance2 = sqrt(pow((cluster.eta() - track.eta()), 2) + pow(dphi2,2));
            if (distance1 < Perpendicular_Cone_Radius) sum_Pt += track.pt();
            if (distance2 < Perpendicular_Cone_Radius) sum_Pt += track.pt();
        }
        double Rho = sum_Pt/(2*M_PI*pow(Perpendicular_Cone_Radius, 2));
        return Rho;
    }

    //Calculates the Pt_Isolation. (Check if the photon candidate is isolated)
    double Pt_Iso(double sum_Pt_tracks_in_cone, double rho, double Cone_Radius = 0.4){
        double Pt_Iso = sum_Pt_tracks_in_cone - rho*M_PI*pow(Cone_Radius, 2);
        return Pt_Iso;
    }
    
    //fill ABCD histograms
    void fillABCDHistograms(HistogramRegistry registry, double regionValue, const auto& cluster) {
        registry.fill(HIST("hABCDControlRegion"), regionValue);
            if (cluster.energy() > 10.) registry.fill(HIST("hABCDControlRegionE_10"), regionValue);
            if (cluster.energy() > 20.) registry.fill(HIST("hABCDControlRegionE_20"), regionValue);
            if (cluster.energy() > 30.) registry.fill(HIST("hABCDControlRegionE_30"), regionValue);
            if (cluster.energy() > 40.) registry.fill(HIST("hABCDControlRegionE_40"), regionValue);
            if (cluster.energy() > 50.) registry.fill(HIST("hABCDControlRegionE_50"), regionValue);
    }

    // process monte carlo data
    void processMC(selectedCollisions::iterator const& theCollision, selectedMCClusters const& mcclusters, aod::StoredMcParticles_001 const&, myGlobTracks const& tracks, o2::aod::EMCALClusterCells const& emccluscells, o2::aod::EMCALMatchedTracks const& matchedtracks){
        MC_Info.fill(HIST("hPosZ"), theCollision.posZ());
        for (auto& mccluster : mcclusters){
            MC_Info.fill(HIST("hECluster"), mccluster.energy());
            
            auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, mccluster.globalIndex());

            if (!track_matching(mccluster, tracksofcluster)){ //no track with significant momentum is matched to cluster
                double Pt_Cone = sum_Pt_tracks_in_cone(mccluster, tracks);
                double Rho_Perpen_Cone = Rho_Perpendicular_Cone(mccluster, tracks);
                double Pt_iso = Pt_Iso(Pt_Cone, Rho_Perpen_Cone);
                MC_Info.fill(HIST("hSigmaLongvsPtIso"), mccluster.m02(), Pt_iso);
                MC_Info.fill(HIST("hPtIso"), Pt_iso);
                MC_Info.fill(HIST("hRho_Perpen_Cone"), Rho_Perpen_Cone);
                MC_Info.fill(HIST("hShowerShape"), mccluster.m02(), mccluster.m20());
                // fill abcd histograms
                if ((Pt_iso < 1.5) && (mccluster.m02() < 0.3) && (mccluster.m02() > 0.1)) {
                    fillABCDHistograms(MC_Info, 0.5, mccluster);
                }
                if ((Pt_iso > 4.0) && (mccluster.m02() < 0.3) && (mccluster.m02() > 0.1)) {
                    fillABCDHistograms(MC_Info, 1.5, mccluster);
                }
                if ((Pt_iso < 1.5) && (mccluster.m02() < 2.0) && (mccluster.m02() > 0.4)) {
                    fillABCDHistograms(MC_Info, 2.5, mccluster);
                }
                if ((Pt_iso > 4.0) && (mccluster.m02() < 2.0) && (mccluster.m02() > 0.4)) {
                    fillABCDHistograms(MC_Info, 3.5, mccluster);
                }
            
                //acces mc true info
                auto ClusterParticles = mccluster.mcParticle_as<aod::StoredMcParticles_001>();
                bool background = true;
                for (auto& clusterparticle : ClusterParticles){
                    if (clusterparticle.pdgCode() == 22){
                        MC_Info.fill(HIST("hMcParticlesToCluster"), clusterparticle.e());
                        MC_Info.fill(HIST("hMcParticleStatusCode"), clusterparticle.getGenStatusCode());
                        if ((clusterparticle.getGenStatusCode() >= 21) && (clusterparticle.getGenStatusCode() <= 29)) background = false;//Particles from the hardest subprocess 
                    }
                }
                if (background){
                    if ((Pt_iso < 1.5) && (mccluster.m02() < 0.3) && (mccluster.m02() > 0.1)) {
                        fillABCDHistograms(MC_Info, 4.5, mccluster);
                    }
                }
            }
        }
    }

    PROCESS_SWITCH(Photon_Isolation_QA, processMC, "proces MC data", true);

    void processData(selectedCollisions::iterator const& theCollision, selectedClusters const& clusters, o2::aod::EMCALClusterCells const& emccluscells, o2::aod::EMCALMatchedTracks const& matchedtracks, myGlobTracks const& alltracks){
        Cluster_Info.fill(HIST("hPosZ"), theCollision.posZ());

        if (clusters.size() > 0){
            Cluster_Info.fill(HIST("hNumClusters"), clusters.size());
        }
        
        for (const auto& cluster : clusters){ 
            Cluster_Info.fill(HIST("hClusterLocation"), cluster.eta(), cluster.phi());
            Cluster_Info.fill(HIST("hEnergy"), cluster.energy());
            Cluster_Info.fill(HIST("hEnergy_ShowerShapeLong"), cluster.energy(), cluster.m02());
            Cluster_Info.fill(HIST("hEnergy_ShowerShapeShort"), cluster.energy(), cluster.m20());
            Cluster_Info.fill(HIST("hEnergy_NCells"), cluster.energy(), cluster.nCells());
            Cluster_Info.fill(HIST("hEnergyvsSigmaLong"), cluster.energy(), cluster.m02());
            
            auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
            Matched_Track_Info.fill(HIST("hNumTracks"), tracksofcluster.size());
            Matched_Track_Info.fill(HIST("hEvsNumTracks"), cluster.energy(), tracksofcluster.size());
            
            if (!track_matching(cluster, tracksofcluster)){ //no track with significant momentum is matched to cluster
                double Pt_Cone = sum_Pt_tracks_in_cone(cluster, alltracks);
                double Rho_Perpen_Cone = Rho_Perpendicular_Cone(cluster, alltracks);
                double Pt_iso = Pt_Iso(Pt_Cone, Rho_Perpen_Cone);
                Photon_Isolation.fill(HIST("hPtIso"), Pt_iso);
                Photon_Isolation.fill(HIST("hRho_Perpen_Cone"), Rho_Perpen_Cone);
                Photon_Isolation.fill(HIST("hShowerShape"), cluster.m02(), cluster.m20());
                Photon_Isolation.fill(HIST("hSigmaLongvsPtIso"), cluster.m02(), Pt_iso);
                for (const auto& match : tracksofcluster){
                    Photon_Isolation.fill(HIST("hPtTracksOfCluster"), match.track_as<myGlobTracks>().pt());
                }
                // fill abcd histograms
                if ((Pt_iso < 1.5) && (cluster.m02() < 0.3) && (cluster.m02() > 0.1)) {
                    fillABCDHistograms(Photon_Isolation, 0.5, cluster);
                }
                if ((Pt_iso > 4.0) && (cluster.m02() < 0.3) && (cluster.m02() > 0.1)) {
                    fillABCDHistograms(Photon_Isolation, 1.5, cluster);
                }
                if ((Pt_iso < 1.5) && (cluster.m02() < 2.0) && (cluster.m02() > 0.4)) {
                    fillABCDHistograms(Photon_Isolation, 2.5, cluster);
                }
                if ((Pt_iso > 4.0) && (cluster.m02() < 2.0) && (cluster.m02() > 0.4)) {
                    fillABCDHistograms(Photon_Isolation, 3.5, cluster);
                }
            }
            
            for (const auto& match : tracksofcluster){
                Matched_Track_Info.fill(HIST("hTrackPt"), match.track_as<myGlobTracks>().pt());
                Matched_Track_Info.fill(HIST("hTrackEta"), match.track_as<myGlobTracks>().eta());
                Matched_Track_Info.fill(HIST("hTrackPhi"), match.track_as<myGlobTracks>().phi());
                Matched_Track_Info.fill(HIST("hRatioClusterETrackP"), cluster.energy()/abs(match.track_as<myGlobTracks>().p()));
                Matched_Track_Info.fill(HIST("hRatioClusterETrackPrackPt"), cluster.energy()/abs(match.track_as<myGlobTracks>().p()), match.track_as<myGlobTracks>().pt());
            }
        }
    }

    PROCESS_SWITCH(Photon_Isolation_QA, processData, "proces data", true);
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Photon_Isolation_QA>(cfgc, TaskName{"photon-isolation-qa"})};
};
