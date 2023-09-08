////////////////////////////////////////////////////////////////////////
// Class:       PhotonBDTNtuple
// Plugin Type: analyzer (art v3_01_02)
// File:        PhotonBDTNtuple_module.cc
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

namespace hyperon {
  class PhotonBDTNtuple;
}

double DEFAULT_INT = -999;
double DEFAULT_DOUBLE = -999.0;
double DEFAULT_BOOL = false;

class hyperon::PhotonBDTNtuple : public art::EDAnalyzer {
public:
  explicit PhotonBDTNtuple(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonBDTNtuple(PhotonBDTNtuple const&) = delete;
  PhotonBDTNtuple(PhotonBDTNtuple&&) = delete;
  PhotonBDTNtuple& operator=(PhotonBDTNtuple const&) = delete;
  PhotonBDTNtuple& operator=(PhotonBDTNtuple&&) = delete;

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
  void FillMCSigmaInfo(art::Event const& evt);
  bool IsSignalEvent(const GeneratorTruth & genTruth, const G4Truth &G4Truth);
  void IdentifySignalPhoton(const G4Truth &G4Truth);
  void FillPFPInformation(art::Event const& evt);
  bool FillPFPMatchInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp, 
      int &matchedTrackID);
  void CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
      std::vector<art::Ptr<recob::Hit>> &hits);
  void FillPFPHitInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillNuVertexInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPMiscInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPShowerInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  void FillPFPEnergyInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
  double GetMedianValue(const std::vector<float> &inputVector);
private:

    // Product labels
    fhicl::ParameterSet m_G4Labels;
    fhicl::ParameterSet m_generatorLabels;
    std::string m_MCParticleModuleLabel;
    std::string m_HitModuleLabel;
    std::string m_PandoraModuleLabel;
    std::string m_BacktrackModuleLabel;
    std::string m_TrackModuleLabel;
    std::string m_ShowerModuleLabel;
    std::string m_CalorimetryModuleLabel;

    // BDT trees
    TTree * m_electronTree;
    TTree * m_photonTree;
    TTree * m_sigmaSignalTree;
    TTree * m_sigmaBackgroundTree;

    // BDT tree variables
    int m_event;
    int m_run;
    int m_subrun;
    int m_truePDG;
    int m_isTrueSigmaPhoton;
    double m_completeness;
    double m_purity;
    int m_generation;
    int m_pandoraPFPCode;
    double m_trackShowerScore;
    double m_nuVertexSeparation;
    double m_trackParentSeparation;
    double m_showerOpeningAngle;
    double m_showerLength;
    int m_nHits2D;
    int m_nHits3D;
    int m_nHitsU;
    int m_nHitsV;
    int m_nHitsW;
    double m_totalEnergy;
    double m_nuVertexChargeDistribution;
    double m_initialdEdx;
    std::vector<double> m_sigmaSystemEnergy;
    std::vector<double> m_sigmaSystemOpeningAngle;
    std::vector<double> m_sigmaSystemVertexSeparation;
    std::vector<double> m_sigmaSystemTrack1;
    std::vector<double> m_sigmaSystemTrack2;

    // Analyzer Variables
    int m_trueNuSliceID;
    bool m_isTrueSigmaEvent;
    int m_mcTruthIndex;
    int m_truePhotonID;

    std::map<int, int> m_hitToTrackID;
    std::map<int, std::vector<int>> m_trackIDToHits;

    // Linking TrackID -> MCParticle
    lar_pandora::MCParticleMap m_mcParticleMap;

    // Linking Self() -> PFParticle
    lar_pandora::PFParticleMap m_pfpMap;

    bool m_debug;
};

hyperon::PhotonBDTNtuple::PhotonBDTNtuple(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset},
    m_G4Labels(pset.get<fhicl::ParameterSet>("Geant4")),
    m_generatorLabels(pset.get<fhicl::ParameterSet>("Generator")),
    m_MCParticleModuleLabel(pset.get<std::string>("MCParticleModuleLabel")),
    m_HitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_PandoraModuleLabel(pset.get<std::string>("PandoraModuleLabel")),
    m_BacktrackModuleLabel(pset.get<std::string>("BacktrackModuleLabel")),
    m_TrackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_ShowerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_CalorimetryModuleLabel(pset.get<std::string>("CalorimetryModuleLabel")),
    m_debug(pset.get<bool>("Debug"))
{
    Reset();
    InitialiseTrees();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::Reset() 
{
    m_event = DEFAULT_INT;
    m_run = DEFAULT_INT;
    m_subrun = DEFAULT_INT;

    ResetPfo();

    m_trueNuSliceID = DEFAULT_INT;
    m_isTrueSigmaEvent = DEFAULT_BOOL;
    m_mcTruthIndex = DEFAULT_INT;
    m_truePhotonID = DEFAULT_INT;
    m_hitToTrackID.clear();
    m_trackIDToHits.clear();
    m_mcParticleMap.clear();
    m_pfpMap.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::ResetPfo() 
{
    m_truePDG = DEFAULT_INT;
    m_isTrueSigmaPhoton = DEFAULT_INT;
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
    m_totalEnergy = DEFAULT_DOUBLE;
    m_nuVertexChargeDistribution = DEFAULT_DOUBLE;
    m_initialdEdx = DEFAULT_DOUBLE;
    m_trackShowerScore = DEFAULT_DOUBLE;
    m_showerOpeningAngle = DEFAULT_DOUBLE;
    m_showerLength = DEFAULT_DOUBLE;
    m_trackParentSeparation = DEFAULT_DOUBLE;
    m_sigmaSystemEnergy.clear();
    m_sigmaSystemOpeningAngle.clear();
    m_sigmaSystemVertexSeparation.clear();
    m_sigmaSystemTrack1.clear();
    m_sigmaSystemTrack2.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::InitialiseTrees()
{
    art::ServiceHandle<art::TFileService> tfs;

    m_electronTree = tfs->make<TTree>("ElectronTree", "ElectronTree");
    m_photonTree = tfs->make<TTree>("PhotonTree", "PhotonTree");
    m_sigmaSignalTree = tfs->make<TTree>("SigmaSignalTree", "SigmaSignalTree");
    m_sigmaBackgroundTree = tfs->make<TTree>("SigmaBackgroundTree", "SigmaBackgroundTree");

    for (TTree* tree : {m_electronTree, m_photonTree, m_sigmaSignalTree, m_sigmaBackgroundTree})
    {
        tree->Branch("Event", &m_event);
        tree->Branch("Run", &m_run);
        tree->Branch("Subrun", &m_subrun);
        tree->Branch("TruePDG", &m_truePDG);
        tree->Branch("IsTrueSigmaPhoton", &m_isTrueSigmaPhoton);
        tree->Branch("Completeness", &m_completeness);
        tree->Branch("Purity", &m_purity);
        tree->Branch("Generation", &m_generation);
        tree->Branch("PandoraPFPCode", &m_pandoraPFPCode);
        tree->Branch("NHits3D", &m_nHits3D);
        tree->Branch("NHits2D", &m_nHits2D);
        tree->Branch("NHitsU", &m_nHitsU);
        tree->Branch("NHitsV", &m_nHitsV);
        tree->Branch("NHitsW", &m_nHitsW);
        tree->Branch("NuVertexSeparation", &m_nuVertexSeparation);
        tree->Branch("TotalEnergy", &m_totalEnergy);
        tree->Branch("NuVertexChargeDistribution", &m_nuVertexChargeDistribution);
        tree->Branch("InitialdEdx", &m_initialdEdx);
        tree->Branch("TrackShowerScore", &m_trackShowerScore);
        tree->Branch("ShowerOpeningAngle", &m_showerOpeningAngle);
        tree->Branch("ShowerLength", &m_showerLength);
        tree->Branch("TrackParentSeparation", &m_trackParentSeparation);
        tree->Branch("SigmaSystemEnergy", &m_sigmaSystemEnergy);
        tree->Branch("SigmaSystemOpeningAngle", &m_sigmaSystemOpeningAngle);
        tree->Branch("SigmaSystemVertexSeparation", &m_sigmaSystemVertexSeparation);
        tree->Branch("SigmaSystemTrack1", &m_sigmaSystemTrack1);
        tree->Branch("SigmaSystemTrack2", &m_sigmaSystemTrack2);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::analyze(art::Event const& evt)
{
    Reset();

    std::cout << "m_PandoraModuleLabel: " << m_PandoraModuleLabel << std::endl;

    m_event = evt.event();
    m_run = evt.run();
    m_subrun = evt.subRun();

    // MCParticle map
    art::Handle<std::vector<recob::Shower>> showerHandle;
    std::vector<art::Ptr<recob::Shower>> showerVector;

    if (!evt.getByLabel(m_ShowerModuleLabel, showerHandle))
        throw cet::exception("PhotonBDTNtuple") << "No Shower Data Products Found!" << std::endl;

    art::fill_ptr_vector(showerVector, showerHandle);

    std::cout << "showerVector.size(): " << showerVector.size() << std::endl;

    //////
    art::Handle<std::vector<recob::Track>> trackHandle;
    std::vector<art::Ptr<recob::Track>> trackVector;

    if (!evt.getByLabel(m_TrackModuleLabel, trackHandle))
        throw cet::exception("PhotonBDTNtuple") << "No Track Data Products Found!" << std::endl;

    art::fill_ptr_vector(trackVector, trackHandle);

    std::cout << "trackVector.size(): " << trackVector.size() << std::endl;


    if (m_debug) std::cout << "Filling Pandora Maps..." << std::endl;
    FillPandoraMaps(evt);
    if (m_debug) std::cout << "Filling MCParticle Maps..." << std::endl;
    FillMCParticleMaps(evt);
    if (m_debug) std::cout << "Filling MC Slice Info..." << std::endl;
    FillMCSliceInfo(evt);
    if (m_debug) std::cout << "Filling MC Sigma Info..." << std::endl;
    FillMCSigmaInfo(evt);
    if (m_debug) std::cout << "Filling PFParticle Information..." << std::endl;
    FillPFPInformation(evt);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::FillPandoraMaps(art::Event const& evt)
{
    // MCParticle map
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

    if (!evt.getByLabel(m_MCParticleModuleLabel, mcParticleHandle))
        throw cet::exception("PhotonBDTNtuple") << "No MCParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(mcParticleVector, mcParticleHandle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, m_mcParticleMap);

    // PFParticle map
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("PhotonBDTNtuple") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector, pfpHandle);

    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector, m_pfpMap);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::FillMCParticleMaps(art::Event const& evt)
{
    // Get event hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVector;

    if (!evt.getByLabel(m_HitModuleLabel, hitHandle))
        throw cet::exception("PhotonBDTNtuple") << "No Hit Data Products Found!" << std::endl;

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

bool hyperon::PhotonBDTNtuple::IsEM(const art::Ptr<simb::MCParticle> &mcParticle)
{
    return ((std::abs(mcParticle->PdgCode()) == 11) || (mcParticle->PdgCode() == 22));
}

///////////////////////////////////////////////////////////////////////////////////////////

// If its an EM particle, we have to move up the EM hierarchy
int hyperon::PhotonBDTNtuple::GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle)
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

void hyperon::PhotonBDTNtuple::FillMCSliceInfo(art::Event const& evt)
{
    // Get slice information
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel(m_PandoraModuleLabel, sliceHandle))
        throw cet::exception("PhotonBDTModuleLabel") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::Hit> hitAssoc = art::FindManyP<recob::Hit>(sliceHandle, evt, m_PandoraModuleLabel);
    art::fill_ptr_vector(sliceVector, sliceHandle);

    // Find slice that contains the nu hierarchy hits
    int highestHitNumber(-1);
    std::map<int, int> sliceSignalHitMap;

    for (art::Ptr<recob::Slice> &slice : sliceVector)
    {
        sliceSignalHitMap[slice->ID()] = 0;

        const std::vector<art::Ptr<recob::Hit>> &sliceHits(hitAssoc.at(slice.key()));

        for (const art::Ptr<recob::Hit> &sliceHit : sliceHits)
        {
            if (m_hitToTrackID.find(sliceHit.key()) == m_hitToTrackID.end())
                continue;

            ++sliceSignalHitMap[slice->ID()];
        }

        if ((sliceSignalHitMap[slice->ID()] > highestHitNumber) && (sliceSignalHitMap[slice->ID()] > 0))
        {
            highestHitNumber = sliceSignalHitMap[slice->ID()];
            m_trueNuSliceID = slice->ID();
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::FillMCSigmaInfo(art::Event const& evt)
{
    SubModuleGeneratorTruth* generatorSubModule = new SubModuleGeneratorTruth(evt, m_generatorLabels);
    GeneratorTruth genTruth = generatorSubModule->GetGeneratorTruth();

    SubModuleG4Truth* G4SubModule = new SubModuleG4Truth(evt, m_G4Labels);
    G4Truth G4Truth = G4SubModule->GetG4Info();

    m_isTrueSigmaEvent = IsSignalEvent(genTruth, G4Truth);

    if (m_isTrueSigmaEvent)
        IdentifySignalPhoton(G4Truth);
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::PhotonBDTNtuple::IsSignalEvent(const GeneratorTruth & genTruth, const G4Truth &G4Truth)
{
    bool isSignal = false;

    for (int i = 0; i < genTruth.NMCTruths; i++)
    {
        if (genTruth.Mode.at(i) != "HYP")
            continue;

        if (!G4Truth.InActiveTPC.at(i))
            continue;

        if (genTruth.Neutrino.at(i).PDG != -14)
            continue;

        if (!G4Truth.IsSigmaZeroCharged.at(i))
            continue;

        if (G4Truth.IsAssociatedHyperon.at(i))
            continue;

        isSignal = true;
        m_mcTruthIndex = i;

        break;
    }

    return isSignal;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::IdentifySignalPhoton(const G4Truth &G4Truth)
{
    bool found = false;

    for (const SimParticle &simParticle : G4Truth.SigmaZeroDecayPhoton)
    {
        if (simParticle.MCTruthIndex != m_mcTruthIndex)
            continue;

        if (simParticle.PDG != 22)
            continue;

        if (m_mcParticleMap.find(simParticle.ArtID) == m_mcParticleMap.end())
            continue;

        found = true;
        m_truePhotonID = simParticle.ArtID;

        break;
    }

    if (!found)
        throw cet::exception("PhotonBDTNtuple") << "No MCParticle Found in Truth!" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::FillPFPInformation(art::Event const& evt)
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
            if (m_debug) std::cout << "New PFParticle!" << std::endl;

            ResetPfo();

            int matchedTrackID(-1);

            if (!FillPFPMatchInfo(evt, pfp, matchedTrackID))
                continue;

            /*
            if ((m_truePDG != 22) && (std::abs(m_truePDG) != 11))
                continue;

            if (m_debug) std::cout << "PFParticle Is Photon/Electron" << std::endl;
            */

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
            if (m_debug) std::cout << "Fill Tree..." << std::endl;

            // Set limits for easy viewing...
            m_trackShowerScore = std::max(m_trackShowerScore, -2.0);
            m_nuVertexSeparation = std::max(m_nuVertexSeparation, -2.0);
            m_trackParentSeparation = std::max(m_trackParentSeparation, -2.0);
            m_showerOpeningAngle  = std::max(m_showerOpeningAngle, -2.0);
            m_showerLength = std::min(std::max(m_showerLength, -2.0), 1000.0);
            m_nHits3D = std::min(std::max(m_nHits3D, -2), 1000);
            m_nHits2D = std::min(std::max(m_nHits2D, -2), 3000);
            m_nuVertexChargeDistribution = std::min(std::max(m_nuVertexChargeDistribution, -2.0), 200.0);
            m_initialdEdx = std::min(std::max(m_initialdEdx, -2.0), 10.0);

            // Fill tree...
            if (m_truePDG == 22)
            {
                m_photonTree->Fill();

                if (m_isTrueSigmaPhoton)
                    m_sigmaSignalTree->Fill();
                else
                    m_sigmaBackgroundTree->Fill();
            }
            else
            {
                m_electronTree->Fill();
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::PhotonBDTNtuple::FillPFPMatchInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp, 
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
    m_isTrueSigmaPhoton = m_isTrueSigmaEvent && (matchedTrackID == m_truePhotonID);

    const int nTrueHits = m_trackIDToHits.at(matchedTrackID).size();

    m_completeness = static_cast<double>(maxHits) / static_cast<double>(nTrueHits);
    m_purity = static_cast<double>(maxHits) / pfpHits.size();

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
   std::vector<art::Ptr<recob::Hit>> &hits)
{
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

   if (!evt.getByLabel(m_PandoraModuleLabel, pfparticleHandle))
       throw cet::exception("PhotonBDTNtuple") << "No PFParticle Data Products Found!" << std::endl;

   art::Handle<std::vector<recob::Cluster>> clusterHandle;

   if (!evt.getByLabel(m_PandoraModuleLabel, clusterHandle)) 
       throw cet::exception("PhotonBDTNtuple") << "No Cluster Data Products Found!" << std::endl;

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

void hyperon::PhotonBDTNtuple::FillPFPHitInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("PhotonBDTNtuple") << "No PFParticle Data Products Found!" << std::endl;

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

        if (pandoraView == geo::kW)
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

void hyperon::PhotonBDTNtuple::FillNuVertexInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    ////////////////////////
    // NuVertex Separation
    ////////////////////////
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("PhotonBDTNtuple") << "No PFParticle Data Products Found!" << std::endl;

    // Get Slice in which pfp lives..
    art::FindManyP<recob::Slice> sliceAssoc = art::FindManyP<recob::Slice>(pfpHandle, evt, m_PandoraModuleLabel);
    const std::vector<art::Ptr<recob::Slice>> &slice(sliceAssoc.at(pfp.key()));

    // Now find neutrino
    art::Handle<std::vector<recob::Slice>> sliceHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, sliceHandle))
        throw cet::exception("PhotonBDTNtuple") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc = art::FindManyP<recob::PFParticle>(sliceHandle, evt, m_PandoraModuleLabel);
    const std::vector<art::Ptr<recob::PFParticle>> &slicePFPs(pfpAssoc.at(slice.at(0).key()));

    std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(slicePFPs, neutrinoPFPs);

    if (neutrinoPFPs.size() != 1)
        return;

    art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfpHandle, evt, m_PandoraModuleLabel);

    const std::vector<art::Ptr<recob::Vertex>> &nuVertex(vertexAssoc.at(neutrinoPFPs.at(0).key()));
    const std::vector<art::Ptr<recob::Vertex>> &pfpVertex(vertexAssoc.at(pfp.key()));

    if (nuVertex.empty() || pfpVertex.empty())
        return;

    const double dX(nuVertex.at(0)->position().X() - pfpVertex.at(0)->position().X());
    const double dY(nuVertex.at(0)->position().Y() - pfpVertex.at(0)->position().Y());
    const double dZ(nuVertex.at(0)->position().Z() - pfpVertex.at(0)->position().Z());

    m_nuVertexSeparation = std::sqrt((dX * dX) + (dY * dY) + (dZ * dZ));

    ////////////////////////
    // NuVertex Charge Distribution (in 3D)
    ////////////////////////
    const geo::Vector_t geoNuVertex(nuVertex.at(0)->position().X(), nuVertex.at(0)->position().Y(), nuVertex.at(0)->position().Z());
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
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::FillPFPMiscInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
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
        throw cet::exception("PhotonBDTNtuple") << "No PFParticle Data Products Found!" << std::endl;

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

void hyperon::PhotonBDTNtuple::FillPFPShowerInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("PhotonBDTNtuple") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Shower> showerAssoc = art::FindManyP<recob::Shower>(pfpHandle, evt, m_ShowerModuleLabel);
    const std::vector<art::Ptr<recob::Shower>> &shower = showerAssoc.at(pfp.key());

    if (shower.empty())
        return;

    m_showerOpeningAngle = shower.at(0)->OpenAngle() * 180.0 / 3.14;
    m_showerLength = shower.at(0)->Length();
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::FillPFPEnergyInfo(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;

    if (!evt.getByLabel(m_PandoraModuleLabel, pfpHandle))
        throw cet::exception("PhotonBDTNtuple") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Track> trackAssoc = art::FindManyP<recob::Track>(pfpHandle, evt, m_TrackModuleLabel);

    // Get the track
    const std::vector<art::Ptr<recob::Track>> &track = trackAssoc.at(pfp.key());

    if (track.empty())
        return;

    art::Handle<std::vector<recob::Track>> trackHandle;

    if (!evt.getByLabel(m_TrackModuleLabel, trackHandle))
        throw cet::exception("PhotonBDTNtuple") << "No Track Data Products Found!" << std::endl;

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
    m_totalEnergy = trackCalo.at(index)->KineticEnergy();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double hyperon::PhotonBDTNtuple::GetMedianValue(const std::vector<float> &inputVector)
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

void hyperon::PhotonBDTNtuple::beginJob()
{
    if (m_debug) std::cout << "Starting..." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::PhotonBDTNtuple::endJob()
{
    if (m_debug) std::cout << "Finished..." << std::endl;
}

DEFINE_ART_MODULE(hyperon::PhotonBDTNtuple)
