////////////////////////////////////////////////////////////////////////
// Class:       VisualiseSlice
// Plugin Type: analyzer (art v3_01_02)
// File:        VisualiseSlice_module.cc
//
// Generated at Mon Sep  4 08:24:35 2023 by Isobel Mawby using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT
#include "TTree.h"

// Services
#include "art/Framework/Services/Optional/TFileService.h"

// Hyperon
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"

// larpandora
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

// lardataobj
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

// nusimdata
#include "nusimdata/SimulationBase/MCParticle.h"

// searchingfornues
#include "ubana/searchingfornues/Selection/CommonDefs/PIDFuncs.h"
#include "ubana/searchingfornues/Selection/CommonDefs/LLR_PID.h"
#include "ubana/searchingfornues/Selection/CommonDefs/LLRPID_proton_muon_lookup.h"

namespace pandora {
  class VisualiseSlice;
}

double DEFAULT_INT = -999;
double DEFAULT_DOUBLE = -999.0;
double DEFAULT_BOOL = false;

class pandora::VisualiseSlice : public art::EDAnalyzer {
public:
  explicit VisualiseSlice(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VisualiseSlice(VisualiseSlice const&) = delete;
  VisualiseSlice(VisualiseSlice&&) = delete;
  VisualiseSlice& operator=(VisualiseSlice const&) = delete;
  VisualiseSlice& operator=(VisualiseSlice&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void Reset();
  void ResetPfo();
  void InitialiseTrees(); 
  void FillPandoraMaps(art::Event const& evt);
  void FillMCParticleMaps(art::Event const& evt);
  bool IsEM(const art::Ptr<simb::MCParticle> &mcParticle);
  int GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle);
  void FillMCSliceInfo(art::Event const& evt);
  void FillNuVertexInfo(art::Event const& evt);

  void FillPFPInformation(art::Event const& evt);
  bool FillPFPMatchInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp, 
      int &matchedTrackID);
  void CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
      std::vector<art::Ptr<recob::Hit>> &hits);
  void FillPFPVisualisationInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPHitInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillNuVertexInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPMiscInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPShowerInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPEnergyInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  double GetMedianValue(const std::vector<float> &inputVector);
  void FillPIDInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
private:

    bool m_applySCECorrections;

    // Product labels
    std::string m_MCParticleModuleLabel;
    std::string m_HitModuleLabel;
    std::string m_PandoraModuleLabel;
    std::string m_FlashMatchModuleLabel;
    std::string m_BacktrackModuleLabel;
    std::string m_TrackModuleLabel;
    std::string m_ShowerModuleLabel;
    std::string m_CalorimetryModuleLabel;
    std::string m_PIDModuleLabel;
    int m_sliceMode;

    // BDT trees
    TTree * m_tree;

    // BDT tree variables
    int m_event;
    int m_run;
    int m_subrun;
    int m_truePDG;
    double m_recoNuVertexX;
    double m_recoNuVertexY;
    double m_recoNuVertexZ;

    std::vector<double> m_spacePointsX;
    std::vector<double> m_spacePointsY;
    std::vector<double> m_spacePointsZ;
    double m_completeness;
    double m_purity;
    int m_generation;
    int m_pandoraPFPCode;
    double m_trackShowerScore;
    double m_nuVertexSeparation;
    double m_dca;
    double m_nuVertexChargeDistribution;
    double m_trackParentSeparation;
    double m_showerOpeningAngle;
    double m_showerLength;
    int m_nHits2D;
    int m_nHits3D;
    int m_nHitsU;
    int m_nHitsV;
    int m_nHitsW;
    double m_totalEnergy;
    double m_initialdEdx;
    double m_chiPIDProton;
    double m_chiPIDMuon;
    double m_chiPIDPion;
    double m_braggPIDProton;
    double m_braggPIDMuon;
    double m_braggPIDPion;
    double m_LLRPIDReduced;

    double m_sliceCompleteness;
    double m_slicePurity;
    bool m_trueNuSliceMode;
    bool m_pandoraNuSliceMode;
    bool m_flashMatchNuSliceMode;

    // Analyzer Variables
    int m_trueNuSliceID;
    int m_pandoraNuSliceID;
    int m_flashMatchNuSliceID;

    std::map<int, int> m_hitToTrackID;
    std::map<int, std::vector<int>> m_trackIDToHits;

    // Linking TrackID -> MCParticle
    lar_pandora::MCParticleMap m_mcParticleMap;

    // Linking Self() -> PFParticle
    lar_pandora::PFParticleMap m_pfpMap;

    // LLR_PID
    searchingfornues::LLRPID m_LLRPIDCalculator;
    searchingfornues::ProtonMuonLookUpParameters m_ProtonMuonParams;

    // Get print statements?
    bool m_debug;
};

pandora::VisualiseSlice::VisualiseSlice(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset},
    m_applySCECorrections(pset.get<bool>("ApplySCECorrections")),
    m_MCParticleModuleLabel(pset.get<std::string>("MCParticleModuleLabel")),
    m_HitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_PandoraModuleLabel(pset.get<std::string>("PandoraModuleLabel")),
    m_FlashMatchModuleLabel(pset.get<std::string>("FlashMatchModuleLabel")),
    m_BacktrackModuleLabel(pset.get<std::string>("BacktrackModuleLabel")),
    m_TrackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_ShowerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_CalorimetryModuleLabel(pset.get<std::string>(m_applySCECorrections ? "CalorimetryModuleLabelSCE" : "CalorimetryModuleLabel")),
    m_PIDModuleLabel(pset.get<std::string>(m_applySCECorrections ? "PIDModuleLabelSCE" : "PIDModuleLabel")),
    m_sliceMode(pset.get<int>("SliceMode")),
    m_debug(pset.get<bool>("Debug"))
{
    m_trueNuSliceMode = DEFAULT_BOOL;
    m_pandoraNuSliceMode = DEFAULT_BOOL;
    m_flashMatchNuSliceMode= DEFAULT_BOOL;

    Reset();
    InitialiseTrees();

    if (m_sliceMode == 0)
    {
        m_trueNuSliceMode = true;

        if (m_debug) std::cout << "Running in true Nu Slice mode" << std::endl;
    }
    else if (m_sliceMode == 1)
    {
        m_pandoraNuSliceMode = true;

        if (m_debug) std::cout << "Running in pandora Nu Slice mode" << std::endl;
    }
    else if (m_sliceMode == 2)
    {
        m_flashMatchNuSliceMode = true;

        if (m_debug) std::cout << "Running in flash match Nu Slice mode" << std::endl;
    }

    m_LLRPIDCalculator.set_dedx_binning(0, m_ProtonMuonParams.dedx_edges_pl_0);
    m_LLRPIDCalculator.set_par_binning(0, m_ProtonMuonParams.parameters_edges_pl_0);
    m_LLRPIDCalculator.set_lookup_tables(0, m_ProtonMuonParams.dedx_pdf_pl_0);

    m_LLRPIDCalculator.set_dedx_binning(1, m_ProtonMuonParams.dedx_edges_pl_1);
    m_LLRPIDCalculator.set_par_binning(1, m_ProtonMuonParams.parameters_edges_pl_1);
    m_LLRPIDCalculator.set_lookup_tables(1, m_ProtonMuonParams.dedx_pdf_pl_1);

    m_LLRPIDCalculator.set_dedx_binning(2, m_ProtonMuonParams.dedx_edges_pl_2);
    m_LLRPIDCalculator.set_par_binning(2, m_ProtonMuonParams.parameters_edges_pl_2);
    m_LLRPIDCalculator.set_lookup_tables(2, m_ProtonMuonParams.dedx_pdf_pl_2);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::Reset() 
{
    m_event = DEFAULT_INT;
    m_run = DEFAULT_INT;
    m_subrun = DEFAULT_INT;

    ResetPfo();

    m_trueNuSliceID = DEFAULT_INT;
    m_pandoraNuSliceID = DEFAULT_INT;
    m_flashMatchNuSliceID = DEFAULT_INT;

    m_recoNuVertexX = DEFAULT_DOUBLE;
    m_recoNuVertexY = DEFAULT_DOUBLE;
    m_recoNuVertexZ = DEFAULT_DOUBLE;

    m_sliceCompleteness = DEFAULT_DOUBLE;
    m_slicePurity = DEFAULT_DOUBLE;

    m_hitToTrackID.clear();
    m_trackIDToHits.clear();
    m_mcParticleMap.clear();
    m_pfpMap.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::ResetPfo() 
{
    m_truePDG = DEFAULT_INT;

    m_spacePointsX.clear();
    m_spacePointsY.clear();
    m_spacePointsZ.clear();

    m_completeness = DEFAULT_DOUBLE;
    m_purity = DEFAULT_DOUBLE;
    m_generation = DEFAULT_INT;
    m_pandoraPFPCode = DEFAULT_INT;
    m_nHits2D = DEFAULT_INT;
    m_nHits3D = DEFAULT_INT;
    m_nHitsU = DEFAULT_INT;
    m_nHitsV = DEFAULT_INT;
    m_nHitsW = DEFAULT_INT;
    m_nuVertexSeparation = DEFAULT_DOUBLE;
    m_dca = DEFAULT_DOUBLE;
    m_nuVertexChargeDistribution = DEFAULT_DOUBLE;
    m_totalEnergy = DEFAULT_DOUBLE;
    m_initialdEdx = DEFAULT_DOUBLE;
    m_trackShowerScore = DEFAULT_DOUBLE;
    m_showerOpeningAngle = DEFAULT_DOUBLE;
    m_showerLength = DEFAULT_DOUBLE;
    m_trackParentSeparation = DEFAULT_DOUBLE;
    m_chiPIDProton = DEFAULT_DOUBLE;
    m_chiPIDMuon = DEFAULT_DOUBLE;
    m_chiPIDPion = DEFAULT_DOUBLE;
    m_braggPIDProton = DEFAULT_DOUBLE;
    m_braggPIDMuon = DEFAULT_DOUBLE;
    m_braggPIDPion = DEFAULT_DOUBLE;
    m_LLRPIDReduced = DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::InitialiseTrees()
{
    art::ServiceHandle<art::TFileService> tfs;

    m_tree = tfs->make<TTree>("VisualisationTree", "VisualisationTree");

    m_tree->Branch("Event", &m_event);
    m_tree->Branch("Run", &m_run);
    m_tree->Branch("Subrun", &m_subrun);
    m_tree->Branch("TrueNuSliceID", &m_trueNuSliceID);
    m_tree->Branch("PandoraNuSliceID", &m_pandoraNuSliceID);
    m_tree->Branch("FlashMatchNuSliceID", &m_flashMatchNuSliceID);
    m_tree->Branch("TrueNuSliceMode", &m_trueNuSliceMode);
    m_tree->Branch("PandoraNuSliceMode", &m_pandoraNuSliceMode);
    m_tree->Branch("FlashMatchNuSliceMode", &m_flashMatchNuSliceMode);
    m_tree->Branch("RecoNuVertexX", &m_recoNuVertexX);
    m_tree->Branch("RecoNuVertexY", &m_recoNuVertexY);
    m_tree->Branch("RecoNuVertexZ", &m_recoNuVertexZ);

    m_tree->Branch("SliceCompleteness", &m_sliceCompleteness);
    m_tree->Branch("SlicePurity", &m_slicePurity);

    m_tree->Branch("TruePDG", &m_truePDG);
    m_tree->Branch("SpacePointsX", &m_spacePointsX);
    m_tree->Branch("SpacePointsY", &m_spacePointsY);
    m_tree->Branch("SpacePointsZ", &m_spacePointsZ);

    m_tree->Branch("Completeness", &m_completeness);
    m_tree->Branch("Purity", &m_purity);
    m_tree->Branch("Generation", &m_generation);
    m_tree->Branch("PandoraPFPCode", &m_pandoraPFPCode);
    m_tree->Branch("NHits3D", &m_nHits3D);
    m_tree->Branch("NHits2D", &m_nHits2D);
    m_tree->Branch("NHitsU", &m_nHitsU);
    m_tree->Branch("NHitsV", &m_nHitsV);
    m_tree->Branch("NHitsW", &m_nHitsW);
    m_tree->Branch("NuVertexSeparation", &m_nuVertexSeparation);
    m_tree->Branch("DCA", &m_dca);
    m_tree->Branch("NuVertexChargeDistribution", &m_nuVertexChargeDistribution);
    m_tree->Branch("TotalEnergy", &m_totalEnergy);
    m_tree->Branch("InitialdEdx", &m_initialdEdx);
    m_tree->Branch("TrackShowerScore", &m_trackShowerScore);
    m_tree->Branch("ShowerOpeningAngle", &m_showerOpeningAngle);
    m_tree->Branch("ShowerLength", &m_showerLength);
    m_tree->Branch("TrackParentSeparation", &m_trackParentSeparation);
    m_tree->Branch("ChiPIDProton", &m_chiPIDProton);
    m_tree->Branch("ChiPIDMuon", &m_chiPIDMuon);
    m_tree->Branch("ChiPIDPion", &m_chiPIDPion);
    m_tree->Branch("BraggPIDProton", &m_braggPIDProton);
    m_tree->Branch("BraggPIDMuon", &m_braggPIDMuon);
    m_tree->Branch("BraggPIDPion", &m_braggPIDPion);
    m_tree->Branch("LLRPIDReduced", &m_LLRPIDReduced);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::analyze(art::Event const& evt)
{
    Reset();

    std::cout << "m_PandoraModuleLabel: " << m_PandoraModuleLabel << std::endl;

    m_event = evt.event();
    m_run = evt.run();
    m_subrun = evt.subRun();

    if (m_debug) std::cout << "Filling Pandora Maps..." << std::endl;
    FillPandoraMaps(evt);
    if (m_debug) std::cout << "Filling MCParticle Maps..." << std::endl;
    FillMCParticleMaps(evt);
    if (m_debug) std::cout << "Filling MC Slice Info..." << std::endl;
    FillMCSliceInfo(evt);
    if (m_debug) std::cout << "Filling NuVertex Info..." << std::endl;
    FillNuVertexInfo(evt);
    if (m_debug) std::cout << "Filling PFParticle Information..." << std::endl;
    FillPFPInformation(evt);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPandoraMaps(art::Event const& evt)
{
    // MCParticle map
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

    if (!evt.getByLabel(m_MCParticleModuleLabel, mcParticleHandle))
        throw cet::exception("VisualiseSlice") << "No MCParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(mcParticleVector, mcParticleHandle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, m_mcParticleMap);

    // PFParticle map
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector, pfpHandle);

    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector, m_pfpMap);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillMCParticleMaps(art::Event const& evt)
{
    // Get event hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVector;

    if (!evt.getByLabel(m_HitModuleLabel, hitHandle))
        throw cet::exception("VisualiseSlice") << "No Hit Data Products Found!" << std::endl;

    art::fill_ptr_vector(hitVector, hitHandle);

    // Get backtracker info
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assocMCPart = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitHandle, evt, m_BacktrackModuleLabel);

    // Truth match
    for (unsigned int hitIndex = 0; hitIndex < hitVector.size(); hitIndex++)
    {
        const art::Ptr<recob::Hit> &hit = hitVector[hitIndex];
        const std::vector<art::Ptr<simb::MCParticle>> &matchedMCParticles = assocMCPart.at(hit.key());
        auto matchedDatas = assocMCPart.data(hit.key());

        for (unsigned int mcParticleIndex = 0; mcParticleIndex < matchedMCParticles.size(); ++mcParticleIndex)
        {
            const art::Ptr<simb::MCParticle> &matchedMCParticle = matchedMCParticles.at(mcParticleIndex);
            auto matchedData = matchedDatas.at(mcParticleIndex);

            if (matchedData->isMaxIDE != 1)
                continue;

            // Get hit view
            const int trackID = IsEM(matchedMCParticle) ? GetLeadEMTrackID(matchedMCParticle) : matchedMCParticle->TrackId();

            m_hitToTrackID[hit.key()] = trackID;
            m_trackIDToHits[trackID].push_back(hit.key());

        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool pandora::VisualiseSlice::IsEM(const art::Ptr<simb::MCParticle> &mcParticle)
{
    return ((std::abs(mcParticle->PdgCode()) == 11) || (mcParticle->PdgCode() == 22));
}

///////////////////////////////////////////////////////////////////////////////////////////

// If its an EM particle, we have to move up the EM hierarchy
int pandora::VisualiseSlice::GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle)
{
    int trackID = mcParticle->TrackId();
    art::Ptr<simb::MCParticle> motherMCParticle = mcParticle;

    do
    {
        trackID = motherMCParticle->TrackId();
        const int motherID = motherMCParticle->Mother();

        if (m_mcParticleMap.find(motherID) == m_mcParticleMap.end())
            break;

        motherMCParticle = m_mcParticleMap.at(motherID);
    }
    while (IsEM(motherMCParticle));

    return trackID;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillMCSliceInfo(art::Event const& evt)
{
    // Get slice information
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel(m_PandoraModuleLabel, sliceHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::Hit> hitAssoc = art::FindManyP<recob::Hit>(sliceHandle, evt, m_PandoraModuleLabel);
    art::fill_ptr_vector(sliceVector, sliceHandle);

    /////////////////////////////////////////////
    // Find the true nu slice ID 
    /////////////////////////////////////////////
    int highestHitNumber(-1);
    std::map<int, int> sliceSignalHitMap;
    int totalTrueHits(0);

    for (art::Ptr<recob::Slice> &slice : sliceVector)
    {
        sliceSignalHitMap[slice->ID()] = 0;

        const std::vector<art::Ptr<recob::Hit>> &sliceHits(hitAssoc.at(slice.key()));

        for (const art::Ptr<recob::Hit> &sliceHit : sliceHits)
        {
            if (m_hitToTrackID.find(sliceHit.key()) == m_hitToTrackID.end())
                continue;

            ++sliceSignalHitMap[slice->ID()];
            ++totalTrueHits;
        }

        if ((sliceSignalHitMap[slice->ID()] > highestHitNumber) && (sliceSignalHitMap[slice->ID()] > 0))
        {
            highestHitNumber = sliceSignalHitMap[slice->ID()];
            m_trueNuSliceID = slice->ID();
        }
    }

    /////////////////////////////////////////////
    // Find the pandora nu slice ID 
    /////////////////////////////////////////////
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc = art::FindManyP<recob::PFParticle>(sliceHandle, evt, m_PandoraModuleLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfpHandle, evt, m_PandoraModuleLabel);

    double bestTopologicalScore(-std::numeric_limits<double>::max());

    for (const art::Ptr<recob::Slice> &slice : sliceVector)
    {
        const std::vector<art::Ptr<recob::PFParticle>> slicePFPs = pfpAssoc.at(slice.key());

        for (const art::Ptr<recob::PFParticle> &pfp : slicePFPs)
        {
            // only topological score for the primary pfp
            if (!pfp->IsPrimary())
                continue;

            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMeta = metadataAssn.at(pfp.key());

            if (pfpMeta.empty())
                continue;

            const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfpMeta.at(0)->GetPropertiesMap());

            if (!pfParticlePropertiesMap.empty() && (pfParticlePropertiesMap.find("NuScore") != pfParticlePropertiesMap.end()))
            {
                const double thisTopologicalScore = pfParticlePropertiesMap.at("NuScore");

                if (thisTopologicalScore > bestTopologicalScore)
                {
                    bestTopologicalScore = thisTopologicalScore;
                    m_pandoraNuSliceID = slice->ID();
                }
            }
        }
    }

    /////////////////////////////////////////////
    // Find the flash match nu slice ID 
    /////////////////////////////////////////////
    art::Handle<std::vector<recob::PFParticle>> flashMatchPFPHandle;
    std::vector<art::Ptr<recob::PFParticle>> flashMatchPFPVector;

    if (!evt.getByLabel(m_FlashMatchModuleLabel, flashMatchPFPHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(flashMatchPFPVector, flashMatchPFPHandle);

    std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(flashMatchPFPVector, neutrinoPFPs);

    if (neutrinoPFPs.size() > 1)
    {
        throw cet::exception("SigmaRecoAnalyser") << "Too many neutrinos found!" << std::endl;
    }
    else if (neutrinoPFPs.size() == 1)
    {
        bool found = false;

        art::FindManyP<recob::Slice> flashMatchSliceAssoc = art::FindManyP<recob::Slice>(flashMatchPFPHandle, evt, m_FlashMatchModuleLabel);
        const std::vector<art::Ptr<recob::Slice>> &flashMatchSlices = flashMatchSliceAssoc.at(neutrinoPFPs[0].key());

        if (!flashMatchSlices.empty())
        {
            // Get hits, and then work out if they also live in the pandoraPatRec reco?
            const art::Ptr<recob::Slice> &flashMatchSlice(flashMatchSlices.at(0));

            art::Handle<std::vector<recob::Slice>> flashMatchSliceHandle;

            if (!evt.getByLabel(m_FlashMatchModuleLabel, flashMatchSliceHandle))
                throw cet::exception("PhotonBDTModuleLabel") << "No Flash Match Slice Data Products Found!" << std::endl;

            art::FindManyP<recob::Hit> flashMatchHitAssoc = art::FindManyP<recob::Hit>(flashMatchSliceHandle, evt, m_FlashMatchModuleLabel);
            const std::vector<art::Ptr<recob::Hit>> &flashMatchSliceHits(flashMatchHitAssoc.at(flashMatchSlice.key()));

            // Loop through pandoraPatRecSlices
            for (art::Ptr<recob::Slice> &slice : sliceVector)
            {
                const std::vector<art::Ptr<recob::Hit>> &sliceHits(hitAssoc.at(slice.key()));

                for (const art::Ptr<recob::Hit> &sliceHit : sliceHits)
                {
                    for (const art::Ptr<recob::Hit> &flashMatchSliceHit : flashMatchSliceHits)
                    {
                        if (sliceHit.key() == flashMatchSliceHit.key())
                        {
                            found = true;
                            m_flashMatchNuSliceID = slice.key();
                            break;
                        }
                    }

                    if (found)
                        break;
                }

                if (found)
                    break;
            }
        }
    }

    /////////////////////////////////////////////
    // Completeness/Purity (this is so inefficienct, i'm sorry)
    /////////////////////////////////////////////
    for (art::Ptr<recob::Slice> &slice : sliceVector)
    {
        if ((m_trueNuSliceMode && (slice->ID() == m_trueNuSliceID)) || (m_pandoraNuSliceMode && (slice->ID() == m_pandoraNuSliceID)) || (m_flashMatchNuSliceMode && (slice->ID() == m_flashMatchNuSliceID)))
        {
            const std::vector<art::Ptr<recob::Hit>> &sliceHits(hitAssoc.at(slice.key()));

            const int nSliceHits = sliceHits.size();
            const int nSliceTrueHits = sliceSignalHitMap.find(slice->ID()) == sliceSignalHitMap.end() ? 0 : sliceSignalHitMap.at(slice->ID());

            m_sliceCompleteness = totalTrueHits == 0 ? 0.0 : static_cast<float>(nSliceTrueHits) / static_cast<float>(totalTrueHits);
            m_slicePurity = nSliceHits == 0 ? 0.0 : static_cast<float>(nSliceTrueHits) / static_cast<float>(nSliceHits);

            break;
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillNuVertexInfo(art::Event const& evt)
{
    // Get slice information
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel(m_PandoraModuleLabel, sliceHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc = art::FindManyP<recob::PFParticle>(sliceHandle, evt, m_PandoraModuleLabel);
    art::fill_ptr_vector(sliceVector, sliceHandle);

    // Get PFP handle
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    for (art::Ptr<recob::Slice> &slice : sliceVector)
    {
        if ((m_trueNuSliceMode && (slice->ID() == m_trueNuSliceID)) || (m_pandoraNuSliceMode && (slice->ID() == m_pandoraNuSliceID)) || 
            (m_flashMatchNuSliceMode && (slice->ID() == m_flashMatchNuSliceID)))
        {
            const std::vector<art::Ptr<recob::PFParticle>> &slicePFPs(pfpAssoc.at(slice.key()));

            std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
            lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(slicePFPs, neutrinoPFPs);

            if (neutrinoPFPs.size() != 1)
                return;

            art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfpHandle, evt, m_PandoraModuleLabel);

            const std::vector<art::Ptr<recob::Vertex>> &nuVertex(vertexAssoc.at(neutrinoPFPs.at(0).key()));

            if (nuVertex.empty())
                return;

            m_recoNuVertexX = nuVertex.at(0)->position().X();
            m_recoNuVertexY = nuVertex.at(0)->position().Y();
            m_recoNuVertexZ = nuVertex.at(0)->position().Z();

            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPInformation(art::Event const& evt)
{
    // Get slice information
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel(m_PandoraModuleLabel, sliceHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No Slice Data Products Found!" << std::endl;

    art::fill_ptr_vector(sliceVector, sliceHandle);

    art::FindManyP<recob::PFParticle> pfpAssoc = art::FindManyP<recob::PFParticle>(sliceHandle, evt, m_PandoraModuleLabel);

    for (const art::Ptr<recob::Slice> &slice : sliceVector)
    {
        if (slice->ID() != m_trueNuSliceID)
            continue;

        if (m_debug) std::cout << "Found Slice, Filling PFParticle Information..." << std::endl;

        const std::vector<art::Ptr<recob::PFParticle>> &pfps(pfpAssoc.at(slice.key()));

        for (const art::Ptr<recob::PFParticle> &pfp : pfps)
        {
            // Is the PFP in the neutrino hierarchy??
            const art::Ptr<recob::PFParticle> parentPFP = lar_pandora::LArPandoraHelper::GetParentPFParticle(m_pfpMap, pfp);

            if (!lar_pandora::LArPandoraHelper::IsNeutrino(parentPFP))
                continue;

            if (m_debug) std::cout << "New PFParticle!" << std::endl;

            ResetPfo();

            int matchedTrackID(-1);

            // ATTN: NEED TO NOT DO THIS
            if (!FillPFPMatchInfo(evt, pfp, matchedTrackID))
                continue;

            if (m_debug) std::cout << "Fill PFP Visualisation Information..." << std::endl;
            FillPFPVisualisationInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Hit Information..." << std::endl;
            FillPFPHitInfo(evt, pfp);
            if (m_debug) std::cout << "Fill Nu Vertex Information..." << std::endl;
            FillNuVertexInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Misc Information..." << std::endl;
            FillPFPMiscInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Shower Information..." << std::endl;
            FillPFPShowerInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PFP Energy Information..." << std::endl;
            FillPFPEnergyInfo(evt, pfp);
            if (m_debug) std::cout << "Fill PID Information..." << std::endl;
            FillPIDInfo(evt, pfp);
            if (m_debug) std::cout << "Fill Tree..." << std::endl;

             m_tree->Fill();
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool pandora::VisualiseSlice::FillPFPMatchInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp, 
    int &matchedTrackID)
{
    std::vector<art::Ptr<recob::Hit>> pfpHits;
    CollectHitsFromClusters(evt, pfp, pfpHits);

    std::map<int, int> countingMap;
    for (const art::Ptr<recob::Hit> pfpHit : pfpHits)
    {
        if (m_hitToTrackID.find(pfpHit.key()) == m_hitToTrackID.end())
            continue;

        const int trackID(m_hitToTrackID.at(pfpHit.key()));

        if (countingMap.find(trackID) == countingMap.end())
            countingMap[trackID] = 1;
        else
            ++countingMap[trackID];
    }

    if (countingMap.empty())
        return false;

    int maxHits = -1;

    for (auto &entry : countingMap)
    {
        if ((entry.second > maxHits) || ((entry.second == maxHits) && (entry.first > matchedTrackID)))
        {
            maxHits = entry.second;
            matchedTrackID = entry.first;
        }
    }

    if (m_mcParticleMap.find(matchedTrackID) == m_mcParticleMap.end())
        return false;

    const art::Ptr<simb::MCParticle> &matchedMCParticle(m_mcParticleMap.at(matchedTrackID)); 

    m_truePDG = matchedMCParticle->PdgCode();

    const int nTrueHits = m_trackIDToHits.at(matchedTrackID).size();

    m_completeness = static_cast<double>(maxHits) / static_cast<double>(nTrueHits);
    m_purity = static_cast<double>(maxHits) / pfpHits.size();

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
   std::vector<art::Ptr<recob::Hit>> &hits)
{
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

   if (!evt.getByLabel(m_PandoraModuleLabel, pfparticleHandle))
       throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

   art::Handle<std::vector<recob::Cluster>> clusterHandle;

   if (!evt.getByLabel(m_PandoraModuleLabel, clusterHandle)) 
       throw cet::exception("VisualiseSlice") << "No Cluster Data Products Found!" << std::endl;

   art::FindManyP<recob::Cluster> pfparticleClustersAssoc = art::FindManyP<recob::Cluster>(pfparticleHandle, evt, m_PandoraModuleLabel);
   art::FindManyP<recob::Hit> clusterHitAssoc = art::FindManyP<recob::Hit>(clusterHandle, evt, m_PandoraModuleLabel);

   std::vector<art::Ptr<recob::Cluster>> clusters = pfparticleClustersAssoc.at(pfparticle.key());

   for (const art::Ptr<recob::Cluster> cluster : clusters)
   {
       std::vector<art::Ptr<recob::Hit>> clusterHits = clusterHitAssoc.at(cluster.key());
       hits.insert(hits.end(), clusterHits.begin(), clusterHits.end());
   }
}


///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPVisualisationInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::SpacePoint> spacePointAssoc = art::FindManyP<recob::SpacePoint>(pfpHandle, evt, m_PandoraModuleLabel);
    const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacePoints = spacePointAssoc.at(pfp.key());

    for (const art::Ptr<recob::SpacePoint> &spacePoint : pfpSpacePoints)
    {
        m_spacePointsX.push_back(spacePoint->XYZ()[0]);
        m_spacePointsY.push_back(spacePoint->XYZ()[1]);
        m_spacePointsZ.push_back(spacePoint->XYZ()[2]);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPHitInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::SpacePoint> spacePointAssoc = art::FindManyP<recob::SpacePoint>(pfpHandle, evt, m_PandoraModuleLabel);
    const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacePoints = spacePointAssoc.at(pfp.key());

    m_nHits3D = pfpSpacePoints.size();

    std::vector<art::Ptr<recob::Hit>> pfpHits;
    CollectHitsFromClusters(evt, pfp, pfpHits);

    int uCount(0), vCount(0), wCount(0);

    for (const art::Ptr<recob::Hit> &hit : pfpHits)
    {
        // Get hit view
        const geo::WireID hit_WireID(hit->WireID());
        const geo::View_t hit_View(hit->View());
        const geo::View_t pandoraView(lar_pandora::LArPandoraGeometry::GetGlobalView(hit_WireID.Cryostat, hit_WireID.TPC, hit_View));

        if (pandoraView == geo::kW || pandoraView == geo::kY)
            ++wCount;
        else if (pandoraView == geo::kU)
            ++uCount;
        else if (pandoraView == geo::kV)
            ++vCount;
    }

    m_nHits2D = pfpHits.size();
    m_nHitsU = uCount;
    m_nHitsV = vCount;
    m_nHitsW = wCount;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillNuVertexInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    ////////////////////////
    // NuVertex Separation
    ////////////////////////
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfpHandle, evt, m_PandoraModuleLabel);

    const std::vector<art::Ptr<recob::Vertex>> &pfpVertex(vertexAssoc.at(pfp.key()));

    if (pfpVertex.empty())
    {
        std::cout << "pfp has no vertex!?" << std::endl;
        throw;
        return;
    }

    const double dX(m_recoNuVertexX - pfpVertex.at(0)->position().X());
    const double dY(m_recoNuVertexY - pfpVertex.at(0)->position().Y());
    const double dZ(m_recoNuVertexZ - pfpVertex.at(0)->position().Z());

    m_nuVertexSeparation = std::sqrt((dX * dX) + (dY * dY) + (dZ * dZ));

    ////////////////////////
    // NuVertex Charge Distribution (in 3D)
    ////////////////////////
    const geo::Vector_t geoNuVertex(m_recoNuVertexX, m_recoNuVertexY, m_recoNuVertexZ);
    const geo::Vector_t geoPFPVertex(pfpVertex.at(0)->position().X(), pfpVertex.at(0)->position().Y(), pfpVertex.at(0)->position().Z());
    const geo::Vector_t nuVertexAxis((geoPFPVertex - geoNuVertex).Unit());

    // Need SpacePoints
    art::FindManyP<recob::SpacePoint> spacePointAssoc = art::FindManyP<recob::SpacePoint>(pfpHandle, evt, m_PandoraModuleLabel);
    const std::vector<art::Ptr<recob::SpacePoint>> &spacePoints(spacePointAssoc.at(pfp.key()));

    if (!spacePoints.empty())
    {
        // Need charge info
        art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;

        if (!evt.getByLabel(m_PandoraModuleLabel, spacePointHandle))
            throw cet::exception("PhotonBDTNtuple") << "No SpacePoint Data Products Found!" << std::endl;

        art::FindManyP<recob::Hit> hitAssoc = art::FindManyP<recob::Hit>(spacePointHandle, evt, m_PandoraModuleLabel);

        // Nu vertex charge distribution
        double totalCharge(0.f);
        double nuVertexChargeDistribution = 0.f;

        for (const art::Ptr<recob::SpacePoint> &spacePoint : spacePoints)
        {
            const std::vector<art::Ptr<recob::Hit>> &hit(hitAssoc.at(spacePoint.key()));

            if (hit.empty())
                continue;

            const double hitCharge(std::fabs(hit.at(0)->Integral()));
            totalCharge += hitCharge;

            const geo::Vector_t geoSpacePoint(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]);
            const geo::Vector_t displacement(geoSpacePoint - geoNuVertex);
            const double thisT(std::sqrt(nuVertexAxis.Cross(displacement).Mag2()));

            nuVertexChargeDistribution += (thisT * hitCharge);
        }

        if (totalCharge > std::numeric_limits<double>::epsilon())
            nuVertexChargeDistribution /= totalCharge;

        m_nuVertexChargeDistribution = nuVertexChargeDistribution;
    }

    ////////////////////////
    // Distance of closest approach
    ////////////////////////
    art::FindManyP<recob::Shower> showerAssoc = art::FindManyP<recob::Shower>(pfpHandle, evt, m_ShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower>> &showers = showerAssoc.at(pfp.key());

    if (!showers.empty())
    {
        const art::Ptr<recob::Shower> &shower = showers.at(0);

        const TVector3 nuVertexPosition(m_recoNuVertexX, m_recoNuVertexY, m_recoNuVertexZ);
        const double alpha = std::fabs((shower->ShowerStart() - nuVertexPosition).Dot(shower->Direction()));
        const TVector3 r = shower->ShowerStart() - (alpha * shower->Direction());
        m_dca = (r - nuVertexPosition).Mag();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPMiscInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    ////////////////////////
    // Generation
    ////////////////////////
    m_generation = lar_pandora::LArPandoraHelper::GetGeneration(m_pfpMap, pfp);
    m_pandoraPFPCode = pfp->PdgCode();

    ////////////////////////
    // TrackShower Score
    ////////////////////////
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfpHandle, evt, m_PandoraModuleLabel);
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(pfp.key());

    if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("TrackScore") != pfpMetadata[0]->GetPropertiesMap().end()))
        m_trackShowerScore = pfpMetadata[0]->GetPropertiesMap().at("TrackScore");

    ////////////////////////
    // Parent Separation
    ////////////////////////
    if (m_generation > 2)
    {
        const int parentID(pfp->Parent());
        art::FindManyP<recob::SpacePoint> spacePointAssoc = art::FindManyP<recob::SpacePoint>(pfpHandle, evt, m_PandoraModuleLabel);

        if (m_pfpMap.find(parentID) != m_pfpMap.end())
        {
            const art::Ptr<recob::PFParticle> parentPFP(m_pfpMap.at(parentID));
            const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacePoints(spacePointAssoc.at(pfp.key()));
            const std::vector<art::Ptr<recob::SpacePoint>> &parentSpacePoints(spacePointAssoc.at(parentPFP.key()));

            double closestDistanceSq(std::numeric_limits<double>::max());

            for (const art::Ptr<recob::SpacePoint> &pfpSpacePoint : pfpSpacePoints)
            {
                for (const art::Ptr<recob::SpacePoint> &parentSpacePoint : parentSpacePoints)
                {
                    const double dX = pfpSpacePoint->XYZ()[0] - parentSpacePoint->XYZ()[0];
                    const double dY = pfpSpacePoint->XYZ()[1] - parentSpacePoint->XYZ()[1];
                    const double dZ = pfpSpacePoint->XYZ()[2] - parentSpacePoint->XYZ()[2];
                    const double thisSep((dX * dX) + (dY * dY) + (dZ * dZ));

                    if (thisSep < closestDistanceSq)
                        closestDistanceSq = thisSep;
                }
            }

            m_trackParentSeparation = std::sqrt(closestDistanceSq);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPShowerInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Shower> showerAssoc = art::FindManyP<recob::Shower>(pfpHandle, evt, m_ShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower>> &shower = showerAssoc.at(pfp.key());

    if (shower.empty())
        return;

    m_showerOpeningAngle = shower.at(0)->OpenAngle() * 180.0 / 3.14;
    m_showerLength = shower.at(0)->Length();
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPFPEnergyInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Track> trackAssoc = art::FindManyP<recob::Track>(pfpHandle, evt, m_TrackModuleLabel);

    // Get the track
    const std::vector<art::Ptr<recob::Track>> &track = trackAssoc.at(pfp.key());

    if (track.empty())
        return;

    art::Handle<std::vector<recob::Track>> trackHandle;

    if (!evt.getByLabel(m_TrackModuleLabel, trackHandle))
        throw cet::exception("VisualiseSlice") << "No Track Data Products Found!" << std::endl;

    art::FindManyP<anab::Calorimetry> caloAssoc = art::FindManyP<anab::Calorimetry>(trackHandle, evt, m_CalorimetryModuleLabel);

    const std::vector<art::Ptr<anab::Calorimetry>> &trackCalo = caloAssoc.at(track.at(0).key());

    if (trackCalo.empty())
        return;

    // Get the correct plane (collection) because it's not always the same                                                                                                        
    bool foundCorrectPlane = false;
    size_t index = 0;

    for (size_t i = 0; i < trackCalo.size(); ++i)
    {
        if (trackCalo.at(i)->PlaneID().Plane == 2)
        {
            foundCorrectPlane = true;
            index = i;
            break;
        }
    }

    if (!foundCorrectPlane)
        return;

    std::vector<float> calibrated_dEdX = trackCalo.at(index)->dEdx();
    const std::vector<float> residualRange = trackCalo.at(index)->ResidualRange();
    const std::vector<geo::Point_t> theXYZPositions = trackCalo.at(index)->XYZ();

    float minRange = std::numeric_limits<float>::max();
    float maxRange = -std::numeric_limits<float>::max();

    for (size_t i = 0; i < calibrated_dEdX.size(); ++i)
    {
        // Sometimes I think the dEdx fails? sometimes, everybody cries...
        // I think you can pick this out with a magnitude of zero
        const geo::Point_t xyzPosition = theXYZPositions[i];
        double magnitude = sqrt((xyzPosition.X() * xyzPosition.X()) + (xyzPosition.Y() * xyzPosition.Y()) + (xyzPosition.Z() * xyzPosition.Z()));

         if (magnitude < std::numeric_limits<double>::epsilon())
             continue;

         if (residualRange[i] < minRange)
             minRange = residualRange[i];

         if (residualRange[i] > maxRange)
             maxRange = residualRange[i];
    }

    // Median energy in first 4cm
    std::vector<float> initialSegmentdEdx;

    for (size_t i = 0; i < calibrated_dEdX.size(); ++i)
    {
        if (((maxRange - residualRange[i]) > 0) && ((maxRange - residualRange[i]) < 4.0))
            initialSegmentdEdx.push_back(calibrated_dEdX[i]);
    }

    m_initialdEdx = GetMedianValue(initialSegmentdEdx);


    // SUM THE ENERGY OF THE HITS!!!!! - I have no idea how uboone do this..
    m_totalEnergy = trackCalo.at(index)->KineticEnergy();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double pandora::VisualiseSlice::GetMedianValue(const std::vector<float> &inputVector)
{
    if (inputVector.empty())
        return -999.0;

    // if odd
    if (inputVector.size() % 2 != 0)
    {
        const int index = std::floor(static_cast<double>(inputVector.size()) / 2.0);
        return inputVector.at(index);
    }

    // if even
    const int firstIndex = inputVector.size() / static_cast<int>(2);
    const int secondIndex = firstIndex - 1;

    return (inputVector.at(firstIndex) + inputVector.at(secondIndex)) / 2.0;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::FillPIDInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("VisualiseSlice") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Track> trackAssoc = art::FindManyP<recob::Track>(pfpHandle, evt, m_TrackModuleLabel);

    // Get the track
    const std::vector<art::Ptr<recob::Track>> &track = trackAssoc.at(pfp.key());

    if (track.empty())
        return;

    // Now get the PID
    art::Handle<std::vector<recob::Track>> trackHandle;

    if (!evt.getByLabel(m_TrackModuleLabel, trackHandle))
        throw cet::exception("VisualiseSlice") << "No Track Data Products Found!" << std::endl;

    art::FindManyP<anab::ParticleID> pidAssoc = art::FindManyP<anab::ParticleID>(trackHandle, evt, m_PIDModuleLabel);

    const std::vector<art::Ptr<anab::ParticleID>> &pids = pidAssoc.at(track.at(0).key());

    if (pids.empty())
        return;

    const art::Ptr<anab::ParticleID> &pid = pids.at(0);

    // Look at the chi2 PID for induction plane
    m_chiPIDProton = searchingfornues::PID(pid, "Chi2", anab::kGOF, anab::kForward, 2212, 2);
    m_chiPIDMuon = searchingfornues::PID(pid, "Chi2", anab::kGOF, anab::kForward, 13, 2);
    m_chiPIDPion = searchingfornues::PID(pid, "Chi2", anab::kGOF, anab::kForward, 211, 2);

    // Look at Bragg peak for induction plane
    m_braggPIDProton = std::max(searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2),
                                searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2));
    m_braggPIDMuon = std::max(searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2),
                              searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2));
    m_braggPIDPion = std::max(searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 211, 2),
                              searchingfornues::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 211, 2));

    // PID LLR calculator
    art::FindManyP<anab::Calorimetry> caloAssoc = art::FindManyP<anab::Calorimetry>(trackHandle, evt, m_CalorimetryModuleLabel);
    const std::vector<art::Ptr<anab::Calorimetry>> &trackCalo = caloAssoc.at(track.at(0).key());

    if (trackCalo.empty())
        return;

    float LLRPID = 0;

    for (const art::Ptr<anab::Calorimetry> &calo : trackCalo)
    {
        if (calo->ResidualRange().size() == 0) continue;

        auto const &plane = calo->PlaneID().Plane;
        auto const &dEdxValues = calo->dEdx();
        auto const &resRange = calo->ResidualRange();
        auto const &pitch = calo->TrkPitchVec();
        std::vector<std::vector<float>> paramValues;
        paramValues.push_back(resRange);
        paramValues.push_back(pitch);

        float caloEnergy = 0;

        for (unsigned int i = 0; i < dEdxValues.size(); i++)
            caloEnergy += dEdxValues[i] * pitch[i];

        const float thisLLRPID = m_LLRPIDCalculator.LLR_many_hits_one_plane(dEdxValues, paramValues, plane);

        LLRPID += thisLLRPID;
    }

    m_LLRPIDReduced = atan(LLRPID / 100.f) * 2.f / 3.14159266;  
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::beginJob()
{
    if (m_debug) std::cout << "Starting..." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////

void pandora::VisualiseSlice::endJob()
{
    if (m_debug) std::cout << "Finished..." << std::endl;
}

DEFINE_ART_MODULE(pandora::VisualiseSlice)
