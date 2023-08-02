////////////////////////////////////////////////////////////////////////
// Class:       SigmaRecoAnalyser
// Plugin Type: analyzer (Unknown Unknown)
// File:        SigmaRecoAnalyser_module.cc
//
// Generated at Tue Jul 25 07:26:20 2023 by Isobel Mawby using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT
#include "TTree.h"
#include "TVector3.h"

// Services
#include "art_root_io/TFileService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// Hyperon
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"

// larpandora
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// lardataobj
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"

const int DEFAULT_INT = -999;
const double DEFAULT_DOUBLE = -999.0;
const std::string DEFAULT_STRING = "DEFAULT"; 
const bool DEFAULT_BOOL = false;
const int MUON_INDEX = 0;
const int PROTON_INDEX = 1;
const int PION_INDEX = 2;
const int GAMMA_INDEX = 3;
const double COMPLETENESS_THRESHOLD = 0.1;
const double PURITY_THRESHOLD = 0.5;

namespace hyperon {
    class SigmaRecoAnalyser;
}

class hyperon::SigmaRecoAnalyser : public art::EDAnalyzer {
public:
    explicit SigmaRecoAnalyser(fhicl::ParameterSet const& p);

    // Plugins should not be copied or assigned.
    SigmaRecoAnalyser(SigmaRecoAnalyser const&) = delete;
    SigmaRecoAnalyser(SigmaRecoAnalyser&&) = delete;
    SigmaRecoAnalyser& operator=(SigmaRecoAnalyser const&) = delete;
    SigmaRecoAnalyser& operator=(SigmaRecoAnalyser&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

    void Reset();
    void InvestigateSlices(art::Event const& evt);
    void CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        std::vector<art::Ptr<recob::Hit>> &hits, const std::string &moduleName, const std::string &instanceName);
    bool IsSignalEvent(const GeneratorTruth & genTruth, const G4Truth &G4Truth);
    void FillPandoraMaps(art::Event const& evt);
    bool IdentifySignalParticles(const G4Truth &G4Truth);
    bool IdentifySignalParticle(const G4Truth &G4Truth, const int index);
    void FillMCParticleTopologyInfo(const GeneratorTruth & genTruth, const G4Truth &G4Truth);
    void FillMCParticleHitInfo(art::Event const& evt);
    bool IsEM(const art::Ptr<simb::MCParticle> &mcParticle);
    int GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle);
    void FillMCSliceInfo(art::Event const& evt);
    void PerformMatching(art::Event const& evt);
    void FindMCParticleMatches(art::Event const& evt, const int particleIndex);
    void CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        std::vector<art::Ptr<recob::Hit>> &hits);
    void FillMatchingInfo(art::Event const& evt);
    double GetSliceScore(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp);
    void FillEventRecoInfo(art::Event const& evt);

private:
    typedef std::map<int, std::vector<std::pair<int,int>>> MatchingMap; // [Particle Type Index, [PFP ID, nHits]]

    // Product labels
    fhicl::ParameterSet m_G4Labels;
    fhicl::ParameterSet m_generatorLabels;
    std::string m_MCParticleLabel;
    std::string m_PFParticleLabel;
    std::string m_PFParticleInstance;
    std::string m_HitLabel;
    std::string m_BacktrackLabel;
    std::string m_ClusterLabel;
    std::string m_ClusterInstance;
    //std::string m_SliceLabel;

    // Debug Info?
    bool m_debugMode;

    // ID info 
    int m_mcTruthIndex;
    std::vector<int> m_trueParticleID;
    std::vector<int> m_bestMatchedID;
    std::vector<std::vector<int>> m_matchIDs;
    int m_trueNuSliceID;
    int m_trueNuSliceNHits;

    // Analyser tree
    TTree * m_tree;

    // Event ID
    unsigned int m_eventID;
    int m_run;
    int m_subRun;
    int m_event;

    // Truth stuff (event level) 
    bool m_foundTrueNuVertexSCE;
    double m_trueNuVertexSCEX;
    double m_trueNuVertexSCEY;
    double m_trueNuVertexSCEZ;
    double m_trueProtonPiOpeningAngle;
    double m_trueGammaLambdaOpeningAngle;
    double m_trueGammaLambdaVertexSep;

    // Truth stuff (particle level)
    std::vector<int> m_truePDG;
    std::vector<int> m_nTrueHits;
    std::vector<double> m_trueNuVertexSep;

    // Reco stuff (event level)
    double m_recoNuSliceScore;
    int m_recoNuSliceID;
    bool m_foundRecoNuVertex;
    double m_recoNuVertexX;
    double m_recoNuVertexY;
    double m_recoNuVertexZ;
    bool m_protonPionMerged;
    bool m_isPionContaminatedByProton;
    bool m_isProtonContaminatedByPion;

    // Reco stuff (particle level)
    std::vector<int> m_nMatches;
    std::vector<double> m_bestMatchCompleteness;
    std::vector<double> m_bestMatchPurity;
    std::vector<double> m_bestMatchTrackScore;
    std::vector<int> m_bestMatchGeneration;
    std::vector<double> m_bestMatchSliceScore;
    std::vector<double> m_bestMatchRecoVertexX;
    std::vector<double> m_bestMatchRecoVertexY;
    std::vector<double> m_bestMatchRecoVertexZ;
    std::vector<double> m_bestMatchNuVertexSep;

    std::vector<int> m_matchFoundInSlice;
    std::vector<int> m_matchFoundInOtherSlice;

    int m_nParticlesFoundInSlice;
    int m_nParticlesFoundInOtherSlice;
    int m_nParticlesFoundOverall;

    // Linking particleTypeIndex -> list of hit.key()
    typedef std::map<int, int> TrueHitMap;
    TrueHitMap m_trueHitMap;

    // Linking TrackID -> MCParticle
    lar_pandora::MCParticleMap m_mcParticleMap;

    // Linking Self() -> PFParticle
    lar_pandora::PFParticleMap m_pfpMap;
    lar_pandora::PFParticleMap m_pfpMap_AO;

    // Matching map
    MatchingMap m_nuSliceMatchingMap;
    MatchingMap m_otherSliceMatchingMap;
};

///////////////////////////////////////////////////////////////////////////////////////////

hyperon::SigmaRecoAnalyser::SigmaRecoAnalyser(fhicl::ParameterSet const& pset)
    : EDAnalyzer{pset},
    m_G4Labels(pset.get<fhicl::ParameterSet>("Geant4")),
    m_generatorLabels(pset.get<fhicl::ParameterSet>("Generator")),
    m_MCParticleLabel(pset.get<std::string>("MCParticleLabel")),
    m_PFParticleLabel(pset.get<std::string>("PFParticleLabel")),
    m_PFParticleInstance(pset.get<std::string>("PFParticleInstance", "")),
    m_HitLabel(pset.get<std::string>("HitLabel")),
    m_BacktrackLabel(pset.get<std::string>("BacktrackLabel")),
    m_ClusterLabel(pset.get<std::string>("ClusterLabel")),
    m_ClusterInstance(pset.get<std::string>("ClusterInstance", "")),
    //m_SliceLabel(pset.get<std::string>("SliceLabel")),
    m_debugMode(pset.get<bool>("DebugMode"))
{
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::analyze(art::Event const& evt)
{
    // Reset tree
    Reset();

    if (m_debugMode) std::cout << "//////////////////" << std::endl;
    if (m_debugMode) std::cout << "New Event!" << std::endl;
    if (m_debugMode) std::cout << "//////////////////" << std::endl;

    /////////////////////////
    // Event info
    /////////////////////////
    m_eventID = evt.id().event();
    m_run = evt.run();
    m_subRun = evt.subRun();
    m_event = evt.event();

    m_truePDG[MUON_INDEX] = -13;
    m_truePDG[PROTON_INDEX] = 2212;
    m_truePDG[PION_INDEX] = -211;
    m_truePDG[GAMMA_INDEX] = 22;


    /////////////////////////
    // Investigating slices
    /////////////////////////
    InvestigateSlices(evt);


    /////////////////////////
    // Truth info
    /////////////////////////
    if (m_debugMode) std::cout << "Getting G4 info..." << std::endl;

    SubModuleGeneratorTruth* generatorSubModule = new SubModuleGeneratorTruth(evt, m_generatorLabels);
    GeneratorTruth genTruth = generatorSubModule->GetGeneratorTruth();

    if (m_debugMode) std::cout << "Getting generator info..." << std::endl;

    SubModuleG4Truth* G4SubModule = new SubModuleG4Truth(evt, m_G4Labels);
    G4Truth G4Truth = G4SubModule->GetG4Info();

    if (!IsSignalEvent(genTruth, G4Truth))
    {
        if (m_debugMode) std::cout << "Event is not signal, returning..." << std::endl;
        return;
    }

    if (m_debugMode) std::cout << "Have signal event!" << std::endl;

    FillPandoraMaps(evt);

    if (!IdentifySignalParticles(G4Truth))
    {
        if (m_debugMode) std::cout << "Couldn't identify the signal MCParticles, returning..." << std::endl;
        return;
    }

    if (m_debugMode) std::cout << "Fill MC info..." << std::endl;
    FillMCParticleTopologyInfo(genTruth, G4Truth);
    FillMCParticleHitInfo(evt);
    FillMCSliceInfo(evt);

    /////////////////////////
    // Reco info
    /////////////////////////
    FillEventRecoInfo(evt);

    if (m_debugMode) std::cout << "Perform matching..." << std::endl;
    PerformMatching(evt);

    if (m_debugMode) std::cout << "Fill reco info..." << std::endl;
    FillMatchingInfo(evt);

    m_tree->Fill();
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::InvestigateSlices(art::Event const& evt)
{
    // Slice vector - OG
    art::Handle<std::vector<recob::Slice>> sliceHandle_OG;
    std::vector<art::Ptr<recob::Slice>> sliceVector_OG;

    if (!evt.getByLabel("pandora", sliceHandle_OG))
        throw cet::exception("SigmaRecoAnalyser") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc_OG = art::FindManyP<recob::PFParticle>(sliceHandle_OG, evt, "pandora");
    art::fill_ptr_vector(sliceVector_OG, sliceHandle_OG);

    // Fill slice ID -> PFParticle map (OG)
    std::map<int, std::vector<art::Ptr<recob::PFParticle>>> sliceMap_OG;

    for (art::Ptr<recob::Slice> &slice_OG : sliceVector_OG)
        sliceMap_OG[slice_OG->ID()] = pfpAssoc_OG.at(slice_OG.key());

    // Slice vector - AO
    art::InputTag sliceInputTag("pandoraPatRec", "allOutcomes");
    art::Handle<std::vector<recob::Slice>> sliceHandle_AO;
    std::vector<art::Ptr<recob::Slice>> sliceVector_AO;

    if (!evt.getByLabel(sliceInputTag, sliceHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc_AO = art::FindManyP<recob::PFParticle>(sliceHandle_AO, evt, sliceInputTag);
    art::fill_ptr_vector(sliceVector_AO, sliceHandle_AO);

    // Fill slice ID -> PFParticle map (AO)
    std::map<int, std::vector<art::Ptr<recob::PFParticle>>> sliceMap_AO;

    for (art::Ptr<recob::Slice> &slice_AO : sliceVector_AO)
        sliceMap_AO[slice_AO->ID()] = pfpAssoc_AO.at(slice_AO.key());

    //////////////////////////////////////////////////////
    // PFParticle map
    art::Handle<std::vector<recob::PFParticle>> pfpHandle_AO;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector_AO;

    art::InputTag pfpInputTag_AO("pandoraPatRec", "allOutcomes");

    if (!evt.getByLabel(pfpInputTag_AO, pfpHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector_AO, pfpHandle_AO);

    lar_pandora::PFParticleMap pfpMap_AO;
    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector_AO, pfpMap_AO);

    for (auto &entry : sliceMap_AO)
    {
        std::cout << "entry.first: " << entry.first << std::endl;
        std::cout << "------------" << std::endl;
        for (art::Ptr<recob::PFParticle> &fu : sliceMap_AO.at(entry.first))
        {
            std::cout << "Is parent a neutrino?: " << (lar_pandora::LArPandoraHelper::IsNeutrino(lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpMap_AO, fu))) << std::endl;
            std::cout << "generation: " << lar_pandora::LArPandoraHelper::GetGeneration(pfpMap_AO, fu) << std::endl;
        }
    }
    //////////////////////////////////////////////////////

    // Look at first slice
    std::cout << "pfps_OG.size(): " << sliceMap_OG.at(1).size() << std::endl;
    std::cout << "pfps_AO.size(): " << sliceMap_AO.at(1).size() << std::endl;

    std::vector<art::Ptr<recob::Hit>> allHits_OG;
    std::vector<art::Ptr<recob::Hit>> allHits_AO;

    for (art::Ptr<recob::PFParticle> pfp : sliceMap_OG.at(1))
    {
        std::vector<art::Ptr<recob::Hit>> thisHits;
        CollectHitsFromClusters(evt, pfp, thisHits, "pandora", "");
        allHits_OG.insert(allHits_OG.begin(), thisHits.begin(), thisHits.end());
    }

    for (art::Ptr<recob::PFParticle> pfp : sliceMap_AO.at(1))
    {
        if(!lar_pandora::LArPandoraHelper::IsNeutrino(lar_pandora::LArPandoraHelper::GetParentPFParticle(pfpMap_AO, pfp)))
            continue;

        std::vector<art::Ptr<recob::Hit>> thisHits;
        CollectHitsFromClusters(evt, pfp, thisHits, "pandoraPatRec", "allOutcomes");
        allHits_AO.insert(allHits_AO.begin(), thisHits.begin(), thisHits.end());
    }

    std::cout << "allHits_OG.size(): " << allHits_OG.size() << std::endl;
    std::cout << "allHits_AO.size(): " << allHits_AO.size() << std::endl;

    for (art::Ptr<recob::Hit> &hit : allHits_OG)
    {
        bool found = false;

        for (const art::Ptr<recob::Hit> &otherHit : allHits_AO)
        {
            if (otherHit.key() == hit.key())
                found = true;
        }
        
        if (!found)
            std::cout << "DIDN'T FIND HIT!!!" << std::endl;
    }

    for (art::Ptr<recob::Hit> &hit : allHits_AO)
    {
        bool found = false;

        for (const art::Ptr<recob::Hit> &otherHit : allHits_OG)
        {
            if (otherHit.key() == hit.key())
                found = true;
        }
        
        if (!found)
            std::cout << "DIDN'T FIND HIT!!!" << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::SigmaRecoAnalyser::IsSignalEvent(const GeneratorTruth & genTruth, const G4Truth &G4Truth)
{
    bool isSignalSigmaZero = false;

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

        isSignalSigmaZero = true;
        m_mcTruthIndex = i;

        break;
    }

    return isSignalSigmaZero;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillPandoraMaps(art::Event const& evt)
{
    // MCParticle map
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

    if (!evt.getByLabel(m_MCParticleLabel, mcParticleHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No MCParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(mcParticleVector, mcParticleHandle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, m_mcParticleMap);

    if (m_debugMode) std::cout << "m_mcParticleMap.size(): " << m_mcParticleMap.size() << std::endl;

    // PFParticle map
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    art::InputTag pfpInputTag(m_PFParticleLabel, m_PFParticleInstance);

    if (!evt.getByLabel(pfpInputTag, pfpHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector, pfpHandle);

    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector, m_pfpMap);

    if (m_debugMode) std::cout << "m_pfpMap.size(): " << m_pfpMap.size() << std::endl;

    // PFParticle map - all outcomes
    art::InputTag pfpInputTag_AO("pandoraPatRec", "allOutcomes");
    art::Handle<std::vector<recob::PFParticle>> pfpHandle_AO;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector_AO;

    if (!evt.getByLabel(pfpInputTag_AO, pfpHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector_AO, pfpHandle_AO);

    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector_AO, m_pfpMap_AO);

    if (m_debugMode) std::cout << "m_pfpMap_AO.size(): " << m_pfpMap_AO.size() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::SigmaRecoAnalyser::IdentifySignalParticles(const G4Truth &G4Truth)
{
    if (!IdentifySignalParticle(G4Truth, MUON_INDEX))
    {
        if (m_debugMode) std::cout << "Couldn't find MC muon :(" << std::endl;
        return false;
    }

    if (!IdentifySignalParticle(G4Truth, PROTON_INDEX))
    {
        if (m_debugMode) std::cout << "Couldn't find MC proton :(" << std::endl;
        return false;
    }

    if (!IdentifySignalParticle(G4Truth, PION_INDEX))
    {
        if (m_debugMode) std::cout << "Couldn't find MC pion :(" << std::endl;
        return false;
    }

    if (!IdentifySignalParticle(G4Truth, GAMMA_INDEX))
    {
        if (m_debugMode) std::cout << "Couldn't find MC gamma :(" << std::endl;
        return false;
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::SigmaRecoAnalyser::IdentifySignalParticle(const G4Truth &G4Truth, const int particleTypeIndex)
{
    const std::vector<SimParticle> &simParticles(particleTypeIndex == MUON_INDEX ? G4Truth.Lepton : 
        particleTypeIndex == PROTON_INDEX ? G4Truth.Decay : particleTypeIndex == PION_INDEX ? G4Truth.Decay : G4Truth.SigmaZeroDecayPhoton);

    const int pdg(particleTypeIndex == MUON_INDEX ? -13 : particleTypeIndex == PROTON_INDEX ? 2212 : particleTypeIndex == PION_INDEX ? -211 : 22);

    bool found = false;

    for (const SimParticle &simParticle : simParticles)
    {
        if (simParticle.MCTruthIndex != m_mcTruthIndex)
            continue;

        if (simParticle.PDG != pdg)
            continue;

        if (m_mcParticleMap.find(simParticle.ArtID) == m_mcParticleMap.end())
            continue;

        found = true;
        m_trueParticleID[particleTypeIndex] = simParticle.ArtID;
        break;
    }

    return found;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillMCParticleTopologyInfo(const GeneratorTruth & genTruth, const G4Truth &G4Truth)
{
    const art::Ptr<simb::MCParticle> &muonMCParticle(m_mcParticleMap.at(m_trueParticleID[MUON_INDEX]));
    const art::Ptr<simb::MCParticle> &protonMCParticle(m_mcParticleMap.at(m_trueParticleID[PROTON_INDEX]));
    const art::Ptr<simb::MCParticle> &pionMCParticle(m_mcParticleMap.at(m_trueParticleID[PION_INDEX]));
    const art::Ptr<simb::MCParticle> &gammaMCParticle(m_mcParticleMap.at(m_trueParticleID[GAMMA_INDEX]));

    const TVector3 &protonMom(protonMCParticle->Momentum().Vect());
    const TVector3 &pionMom(pionMCParticle->Momentum().Vect());
    const TVector3 &lambdaMom(protonMom + pionMom);
    const TVector3 &gammaMom(gammaMCParticle->Momentum().Vect());

    m_trueProtonPiOpeningAngle = protonMom.Angle(pionMom);
    m_trueGammaLambdaOpeningAngle = gammaMom.Angle(lambdaMom);

    const TVector3 nuVertex(genTruth.TruePrimaryVertex_X.at(m_mcTruthIndex), genTruth.TruePrimaryVertex_Y.at(m_mcTruthIndex), 
        genTruth.TruePrimaryVertex_Z.at(m_mcTruthIndex));


    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableSimSpatialSCE() == true)
    {
        m_foundTrueNuVertexSCE = true;

        auto offset = SCE->GetPosOffsets(geo::Point_t(nuVertex.X(), nuVertex.Y(), nuVertex.Z()));

        m_trueNuVertexSCEX = nuVertex.X() - offset.X() + 0.6;
        m_trueNuVertexSCEY = nuVertex.Y() + offset.Y();
        m_trueNuVertexSCEZ = nuVertex.Z() + offset.Z();
    }

    const TVector3 muonVertex(muonMCParticle->Vx(), muonMCParticle->Vy(), muonMCParticle->Vz());
    const TVector3 protonVertex(protonMCParticle->Vx(), protonMCParticle->Vy(), protonMCParticle->Vz());
    const TVector3 pionVertex(pionMCParticle->Vx(), pionMCParticle->Vy(), pionMCParticle->Vz());

    // Gamma, end point is first energy deposit (who designs this?)
    const TVector3 gammaVertex({gammaMCParticle->EndX(), gammaMCParticle->EndY(), gammaMCParticle->EndZ()});

    m_trueGammaLambdaVertexSep = (protonVertex - gammaVertex).Mag();
    m_trueNuVertexSep[MUON_INDEX] = (muonVertex - nuVertex).Mag();
    m_trueNuVertexSep[PROTON_INDEX] = (protonVertex - nuVertex).Mag();
    m_trueNuVertexSep[PION_INDEX] = (pionVertex - nuVertex).Mag();
    m_trueNuVertexSep[GAMMA_INDEX] = (gammaVertex - nuVertex).Mag();
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillMCParticleHitInfo(art::Event const& evt)
{
    // Get event hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVector;

    if (!evt.getByLabel(m_HitLabel, hitHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No Hit Data Products Found!" << std::endl;

    art::fill_ptr_vector(hitVector, hitHandle);

    // Get backtracker info
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assocMCPart = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitHandle, evt, m_BacktrackLabel);

    m_nTrueHits[MUON_INDEX] = 0;
    m_nTrueHits[PROTON_INDEX] = 0;
    m_nTrueHits[PION_INDEX] = 0;
    m_nTrueHits[GAMMA_INDEX] = 0;

    // Truth match to found IDs
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

            int trackID = IsEM(matchedMCParticle) ? GetLeadEMTrackID(matchedMCParticle) : matchedMCParticle->TrackId();

            m_trueHitMap[hit.key()] = trackID;

            if (trackID == m_trueParticleID[MUON_INDEX])
                ++m_nTrueHits[MUON_INDEX];
            else if (trackID == m_trueParticleID[PROTON_INDEX])
                ++m_nTrueHits[PROTON_INDEX];
            else if (trackID == m_trueParticleID[PION_INDEX])
                ++m_nTrueHits[PION_INDEX];
            else if (trackID == m_trueParticleID[GAMMA_INDEX])
                ++m_nTrueHits[GAMMA_INDEX];
        }
    }

    if (m_debugMode)
    {
        for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        {
            std::cout << "particleTypeIndex: " << particleTypeIndex << std::endl;
            std::cout << "nTrueHits: " << m_nTrueHits[particleTypeIndex] << std::endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::SigmaRecoAnalyser::IsEM(const art::Ptr<simb::MCParticle> &mcParticle)
{
    return ((std::abs(mcParticle->PdgCode()) == 11) || (mcParticle->PdgCode() == 22));
}

///////////////////////////////////////////////////////////////////////////////////////////
// If its an EM particle, we have to move up the EM hierarchy
int hyperon::SigmaRecoAnalyser::GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle)
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

void hyperon::SigmaRecoAnalyser::FillMCSliceInfo(art::Event const& evt)
{
    // Get slice information
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel("pandora", sliceHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::Hit> hitAssoc = art::FindManyP<recob::Hit>(sliceHandle, evt, "pandora");
    art::fill_ptr_vector(sliceVector, sliceHandle);

    // Find slice that contains the nu hierarchy hits
    std::map<int, int> sliceSignalHitMap;

    int highestHitNumber(-1);
    int highestHitSliceID(-1);

    for (art::Ptr<recob::Slice> &slice : sliceVector)
    {
        sliceSignalHitMap[slice->ID()] = 0;

        const std::vector<art::Ptr<recob::Hit>> &sliceHits(hitAssoc.at(slice.key()));

        for (const art::Ptr<recob::Hit> &sliceHit : sliceHits)
        {
            if (m_trueHitMap.find(sliceHit.key()) == m_trueHitMap.end())
                continue;

            int trackID = m_trueHitMap.at(sliceHit.key());

            if ((trackID == m_trueParticleID[MUON_INDEX]) || (trackID == m_trueParticleID[PROTON_INDEX]) || 
                (trackID == m_trueParticleID[PION_INDEX]) || (trackID == m_trueParticleID[GAMMA_INDEX]))
            {
                ++sliceSignalHitMap[slice->ID()]; 
            }
        }

        if ((sliceSignalHitMap[slice->ID()] > highestHitNumber) && (sliceSignalHitMap[slice->ID()] > 0))
        {
            highestHitNumber = sliceSignalHitMap[slice->ID()];
            highestHitSliceID = slice->ID();
        }
    }

    if (highestHitSliceID >= 0)
    {
        m_trueNuSliceID = highestHitSliceID;
        m_trueNuSliceNHits = highestHitNumber;
    }

    if (m_debugMode) std::cout << "True nu slice has ID: " <<m_trueNuSliceID << ", and contains: " << m_trueNuSliceNHits << " hits" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::PerformMatching(art::Event const& evt)
{
    // Fill matching map
    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        FindMCParticleMatches(evt, particleTypeIndex);

    // Order the matching map (in order of highest -> lowest shared hits)
    for (auto &entry : m_nuSliceMatchingMap)
    {
        std::sort(entry.second.begin(), entry.second.end(), 
            [](const std::pair<int, int> &a, std::pair<int, int> &b) -> bool 
            {
                // tie-breaker
                if (a.second == b.second)
                    return a.first > b.first;

                return a.second > b.second;
            }
        );
    }

    // Order the matching map (in order of highest -> lowest shared hits)
    for (auto &entry : m_otherSliceMatchingMap)
    {
        std::sort(entry.second.begin(), entry.second.end(), 
            [](const std::pair<int, int> &a, std::pair<int, int> &b) -> bool 
            {
                // tie-breaker
                if (a.second == b.second)
                    return a.first > b.first;

                return a.second > b.second;
            }
        );
    }

    if (m_debugMode)
    {
        std::cout << "Neutrino slice matching map: " << std::endl;

        for (auto &entry : m_nuSliceMatchingMap)
        {
            std::cout << "index: " << entry.first << std::endl;
            std::cout << entry.second.size() << " matches" << std::endl;

            for (auto &pair : entry.second)
                std::cout << "(PFParticle ID: " << pair.first << ", shared hits: " << pair.second << std::endl;
        }

        std::cout << "Other slice matching map: " << std::endl;

        for (auto &entry : m_otherSliceMatchingMap)
        {
            std::cout << "index: " << entry.first << std::endl;
            std::cout << entry.second.size() << " matches" << std::endl;

            for (auto &pair : entry.second)
                std::cout << "(PFParticle ID: " << pair.first << ", shared hits: " << pair.second << std::endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FindMCParticleMatches(art::Event const& evt, const int particleIndex)
{
    // Slice vector
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel("pandora", sliceHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc = art::FindManyP<recob::PFParticle>(sliceHandle, evt, "pandora");
    art::fill_ptr_vector(sliceVector, sliceHandle);

    // Search in the slice vector first...
    bool foundInNuSlice = false;
    for (const art::Ptr<recob::Slice> &slice : sliceVector)
    {
        if (slice->ID() != m_recoNuSliceID)
            continue;

        const std::vector<art::Ptr<recob::PFParticle>> &pfps(pfpAssoc.at(slice.key()));

        for (const art::Ptr<recob::PFParticle> &pfp : pfps)
        {
            std::map<int, int> trackIDToHitMap;

            std::vector<art::Ptr<recob::Hit>> pfpHits;
            CollectHitsFromClusters(evt, pfp, pfpHits, "pandora", "");

            for (const art::Ptr<recob::Hit> pfpHit : pfpHits)
            {
                if (m_trueHitMap.find(pfpHit.key()) == m_trueHitMap.end())
                    continue;

                const int trackID(m_trueHitMap.at(pfpHit.key()));

                if (trackIDToHitMap.find(trackID) == trackIDToHitMap.end())
                    trackIDToHitMap[trackID] = 1;
                else
                    ++trackIDToHitMap[trackID];
            }

            if (trackIDToHitMap.find(m_trueParticleID[particleIndex]) == trackIDToHitMap.end())
                continue;

            int maxHits = -1;
            int maxOwnerID = -1;

            for (auto &entry : trackIDToHitMap)
            {
                if ((entry.second > maxHits) || ((entry.second == maxHits) && (entry.first > maxOwnerID)))
                {
                    maxHits = entry.second;
                    maxOwnerID = entry.first;
                }
            }

            if (maxOwnerID == m_trueParticleID[particleIndex])
            {
                m_nuSliceMatchingMap[particleIndex].push_back(std::pair<int, int>(pfp->Self(), trackIDToHitMap[m_trueParticleID[particleIndex]]));
                foundInNuSlice = true;
            }

            if (particleIndex == PROTON_INDEX)
            {
                if ((maxOwnerID != m_trueParticleID[PROTON_INDEX]) && (maxOwnerID == m_trueParticleID[PION_INDEX]))
                {
                    m_isPionContaminatedByProton = true;
                    foundInNuSlice = true;
                }
            }

            if (particleIndex == PION_INDEX)
            {
                if ((maxOwnerID != m_trueParticleID[PION_INDEX]) && (maxOwnerID == m_trueParticleID[PROTON_INDEX]))
                {
                    m_isProtonContaminatedByPion = true;
                    foundInNuSlice = true;
                }
            }
        }
    }

    if (foundInNuSlice)
        return;

    //////////////////////////////
    // Now search in other slice!
    //////////////////////////////

    // Slice vector - All outcomes
    art::InputTag sliceInputTag_AO("pandoraPatRec", "allOutcomes");
    art::Handle<std::vector<recob::Slice>> sliceHandle_AO;
    std::vector<art::Ptr<recob::Slice>> sliceVector_AO;

    if (!evt.getByLabel(sliceInputTag_AO, sliceHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc_AO = art::FindManyP<recob::PFParticle>(sliceHandle_AO, evt, sliceInputTag_AO);
    art::fill_ptr_vector(sliceVector_AO, sliceHandle_AO);

    for (const art::Ptr<recob::Slice> &slice_AO : sliceVector_AO)
    {
        if (slice_AO->ID() == m_recoNuSliceID)
            continue;

        int maxHits = -1;
        int maxOwnerID = -1;

        const std::vector<art::Ptr<recob::PFParticle>> &pfps_AO(pfpAssoc_AO.at(slice_AO.key()));

        for (const art::Ptr<recob::PFParticle> &pfp_AO : pfps_AO)
        {
            // If not neutrino...
            if (!lar_pandora::LArPandoraHelper::IsNeutrino(lar_pandora::LArPandoraHelper::GetParentPFParticle(m_pfpMap_AO, pfp_AO)))
                continue;

            std::map<int, int> trackIDToHitMap;

            std::vector<art::Ptr<recob::Hit>> pfpHits_AO;
            CollectHitsFromClusters(evt, pfp_AO, pfpHits_AO, "pandoraPatRec", "allOutcomes");

            for (const art::Ptr<recob::Hit> pfpHit_AO : pfpHits_AO)
            {
                if (m_trueHitMap.find(pfpHit_AO.key()) == m_trueHitMap.end())
                    continue;

                const int trackID(m_trueHitMap.at(pfpHit_AO.key()));

                if (trackIDToHitMap.find(trackID) == trackIDToHitMap.end())
                    trackIDToHitMap[trackID] = 1;
                else
                    ++trackIDToHitMap[trackID];
            }

            if (trackIDToHitMap.find(m_trueParticleID[particleIndex]) == trackIDToHitMap.end())
                continue;

            for (auto &entry : trackIDToHitMap)
            {
                if ((entry.second > maxHits) || ((entry.second == maxHits) && (entry.first > maxOwnerID)))
                {
                    maxHits = entry.second;
                    maxOwnerID = entry.first;
                }
            }

            if (maxOwnerID == m_trueParticleID[particleIndex])
                m_otherSliceMatchingMap[particleIndex].push_back(std::pair<int, int>(pfp_AO->Self(), trackIDToHitMap[m_trueParticleID[particleIndex]]));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
   std::vector<art::Ptr<recob::Hit>> &hits)
{
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

    art::InputTag pfpInputTag(m_PFParticleLabel, m_PFParticleInstance);

   if (!evt.getByLabel(pfpInputTag, pfparticleHandle))
       throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

   art::Handle<std::vector<recob::Cluster>> clusterHandle;

    art::InputTag clusterInputTag(m_ClusterLabel, m_ClusterInstance);

   if (!evt.getByLabel(clusterInputTag, clusterHandle)) 
       throw cet::exception("SigmaRecoAnalyser") << "No Cluster Data Products Found!" << std::endl;

   art::FindManyP<recob::Cluster> pfparticleClustersAssoc = art::FindManyP<recob::Cluster>(pfparticleHandle, evt, pfpInputTag);
   art::FindManyP<recob::Hit> clusterHitAssoc = art::FindManyP<recob::Hit>(clusterHandle, evt, clusterInputTag);

   std::vector<art::Ptr<recob::Cluster>> clusters = pfparticleClustersAssoc.at(pfparticle.key());

   for (const art::Ptr<recob::Cluster> cluster : clusters)
   {
       std::vector<art::Ptr<recob::Hit>> clusterHits = clusterHitAssoc.at(cluster.key());
       hits.insert(hits.end(), clusterHits.begin(), clusterHits.end());
   }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
   std::vector<art::Ptr<recob::Hit>> &hits, const std::string &moduleName, const std::string &instanceName)
{
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

   art::InputTag pfpInputTag(moduleName, instanceName);

   if (!evt.getByLabel(pfpInputTag, pfparticleHandle))
       throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

   art::Handle<std::vector<recob::Cluster>> clusterHandle;

   art::InputTag clusterInputTag(moduleName, instanceName);

   if (!evt.getByLabel(clusterInputTag, clusterHandle)) 
       throw cet::exception("SigmaRecoAnalyser") << "No Cluster Data Products Found!" << std::endl;

   art::FindManyP<recob::Cluster> pfparticleClustersAssoc = art::FindManyP<recob::Cluster>(pfparticleHandle, evt, pfpInputTag);
   art::FindManyP<recob::Hit> clusterHitAssoc = art::FindManyP<recob::Hit>(clusterHandle, evt, clusterInputTag);

   std::vector<art::Ptr<recob::Cluster>> clusters = pfparticleClustersAssoc.at(pfparticle.key());

   for (const art::Ptr<recob::Cluster> cluster : clusters)
   {
       std::vector<art::Ptr<recob::Hit>> clusterHits = clusterHitAssoc.at(cluster.key());
       hits.insert(hits.end(), clusterHits.begin(), clusterHits.end());
   }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillMatchingInfo(art::Event const& evt)
{
    // Standard reco output
    art::InputTag pfpInputTag(m_PFParticleLabel, m_PFParticleInstance);
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

    if (!evt.getByLabel(pfpInputTag, pfparticleHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfparticleHandle, evt, pfpInputTag);
    art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfparticleHandle, evt, pfpInputTag);

    // All outcomes output
    art::InputTag pfpInputTag_AO("pandoraPatRec", "allOutcomes");
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle_AO;

    if (!evt.getByLabel(pfpInputTag_AO, pfparticleHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn_AO = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfparticleHandle_AO, evt, pfpInputTag_AO);
    art::FindManyP<recob::Vertex> vertexAssoc_AO = art::FindManyP<recob::Vertex>(pfparticleHandle_AO, evt, pfpInputTag_AO);

    m_nParticlesFoundInSlice = 0;
    m_nParticlesFoundInOtherSlice = 0;

    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        if (m_debugMode) std::cout << "--------------- particleTypeIndex: " << particleTypeIndex << std::endl;

        int nuSliceMatches = 0;
        int otherSliceMatches = 0;

        bool foundMatchInSlice = false;
        bool foundMatchInOtherSlice = false;

        const std::vector<std::pair<int, int>> &nuSliceMatchPairs(m_nuSliceMatchingMap[particleTypeIndex]);
        const std::vector<std::pair<int, int>> &otherSliceMatchPairs(m_otherSliceMatchingMap[particleTypeIndex]);

        nuSliceMatches = nuSliceMatchPairs.size();
        otherSliceMatches = otherSliceMatchPairs.size();

        if ((nuSliceMatches == 0) && (otherSliceMatches == 0))
            continue;

        foundMatchInSlice = (nuSliceMatches != 0);
        foundMatchInOtherSlice = (otherSliceMatches != 0);

        const std::pair<int, int> &matchPair(foundMatchInSlice ? nuSliceMatchPairs.front() : otherSliceMatchPairs.front());
        const art::Ptr<recob::PFParticle> pfp(foundMatchInSlice ? m_pfpMap.at(matchPair.first) : m_pfpMap_AO.at(matchPair.first));

        // Now fill out bestMatchPfo info
        const double completeness(static_cast<double>(matchPair.second) / static_cast<double>(m_nTrueHits[particleTypeIndex]));

        std::vector<art::Ptr<recob::Hit>> pfpHits;
        CollectHitsFromClusters(evt, pfp, pfpHits, (foundMatchInSlice ? "pandora" : "pandoraPatRec"), (foundMatchInSlice ? "" : "allOutcomes"));

        const double purity(static_cast<double>(matchPair.second) / static_cast<double>(pfpHits.size()));

        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = (foundMatchInSlice ? metadataAssn.at(pfp.key()) : metadataAssn_AO.at(pfp.key()));

        double trackScore(DEFAULT_DOUBLE);

        if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("TrackScore") != pfpMetadata[0]->GetPropertiesMap().end()))
            trackScore = pfpMetadata[0]->GetPropertiesMap().at("TrackScore");

        //double sliceScore(GetSliceScore(evt, pfp));
        int generation(lar_pandora::LArPandoraHelper::GetGeneration(foundMatchInSlice ? m_pfpMap : m_pfpMap_AO, pfp));

        std::vector<art::Ptr<recob::Vertex>> vertex = foundMatchInSlice ? vertexAssoc.at(pfp.key()) : vertexAssoc_AO.at(pfp.key());

        double vertexX(DEFAULT_DOUBLE), vertexY(DEFAULT_DOUBLE), vertexZ(DEFAULT_DOUBLE);
        double nuVertexSep(DEFAULT_DOUBLE);

        if (!vertex.empty() && m_foundRecoNuVertex)
        {
            vertexX = vertex[0]->position().X();
            vertexY = vertex[0]->position().Y();
            vertexZ = vertex[0]->position().Z();

            double deltaX(vertexX - m_recoNuVertexX);
            double deltaY(vertexY - m_recoNuVertexY);
            double deltaZ(vertexZ - m_recoNuVertexZ);

            nuVertexSep = std::sqrt((deltaX * deltaX) + (deltaY * deltaY) + (deltaZ *  deltaZ));
        }

        m_bestMatchCompleteness[particleTypeIndex] = completeness;
        m_bestMatchPurity[particleTypeIndex] = purity;
        m_bestMatchTrackScore[particleTypeIndex] = trackScore;
        m_bestMatchGeneration[particleTypeIndex] = generation;
        //m_bestMatchSliceScore[particleTypeIndex] = sliceScore;
        m_bestMatchRecoVertexX[particleTypeIndex] = vertexX;
        m_bestMatchRecoVertexY[particleTypeIndex] = vertexY;
        m_bestMatchRecoVertexZ[particleTypeIndex] = vertexZ;
        m_bestMatchNuVertexSep[particleTypeIndex] = nuVertexSep;

        if (m_debugMode)
        {
            std::cout << "bestMatchCompleteness: " << m_bestMatchCompleteness[particleTypeIndex] << std::endl;
            std::cout << "bestMatchPurity: " << m_bestMatchPurity[particleTypeIndex] << std::endl;
            std::cout << "bestMatchTrackScore: " << m_bestMatchTrackScore[particleTypeIndex] << std::endl;
            std::cout << "bestMatchGeneration: " << m_bestMatchGeneration[particleTypeIndex] << std::endl;
            std::cout << "bestMatchSliceScore: " << m_bestMatchSliceScore[particleTypeIndex] << std::endl;
            std::cout << "bestMatchNuVertexSep: " << m_bestMatchNuVertexSep[particleTypeIndex] << std::endl;
        }

        if (foundMatchInSlice)
            ++m_nParticlesFoundInSlice;

        if (!foundMatchInSlice && foundMatchInOtherSlice)
            ++m_nParticlesFoundInOtherSlice;

        m_nMatches[particleTypeIndex] = foundMatchInSlice ? nuSliceMatches : otherSliceMatches;

        m_matchFoundInSlice[particleTypeIndex] = foundMatchInSlice;
        m_matchFoundInOtherSlice[particleTypeIndex] = !foundMatchInSlice && foundMatchInOtherSlice;

        if (m_debugMode)
        {
            std::cout << "nMatches: " << m_nMatches[particleTypeIndex] << std::endl;
            std::cout << "matchFoundInSlice: " << m_matchFoundInSlice[particleTypeIndex] << std::endl;
            std::cout << "matchFoundInOtherSlice " << m_matchFoundInOtherSlice[particleTypeIndex] << std::endl;
        }
    }

    if (m_matchFoundInSlice[PION_INDEX] && !m_matchFoundInSlice[PROTON_INDEX] && m_isPionContaminatedByProton)
        m_protonPionMerged = true;

    if (m_matchFoundInSlice[PROTON_INDEX] && !m_matchFoundInSlice[PION_INDEX] && m_isProtonContaminatedByPion)
        m_protonPionMerged = true;

    if (m_debugMode)
    {
        std::cout << "Overall matching status: " << std::endl;
        std::cout << "nParticlesFoundInSlice: " << m_nParticlesFoundInSlice << std::endl;
        std::cout << "nParticlesFoundInOtherSlice: " << m_nParticlesFoundInOtherSlice << std::endl;
        std::cout << "protonPionMerged? " << m_protonPionMerged << std::endl;
    }

    m_nParticlesFoundOverall = m_nParticlesFoundInSlice + m_nParticlesFoundInOtherSlice;
}

///////////////////////////////////////////////////////////////////////////////////////////

double hyperon::SigmaRecoAnalyser::GetSliceScore(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfp)
{
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

    art::InputTag pfpInputTag(m_PFParticleLabel, m_PFParticleInstance);

    if (!evt.getByLabel(pfpInputTag, pfparticleHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    double sliceScore = DEFAULT_DOUBLE;
    const art::Ptr<recob::PFParticle> &parentPFParticle(lar_pandora::LArPandoraHelper::GetParentPFParticle(m_pfpMap, pfp));
    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfparticleHandle, evt, pfpInputTag);
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(parentPFParticle.key());

    if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("NuScore") != pfpMetadata[0]->GetPropertiesMap().end()))
        sliceScore = pfpMetadata[0]->GetPropertiesMap().at("NuScore");

    if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("IsClearCosmic") != pfpMetadata[0]->GetPropertiesMap().end()))
        std::cout << "IsClearCosmic: " << pfpMetadata[0]->GetPropertiesMap().at("IsClearCosmic") << std::endl;

    return sliceScore;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillEventRecoInfo(art::Event const& evt)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    art::InputTag pfpInputTag(m_PFParticleLabel, m_PFParticleInstance);

    if (!evt.getByLabel(pfpInputTag, pfpHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector, pfpHandle);

    std::vector<art::Ptr<recob::PFParticle>> neutrinoPFPs;
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(pfpVector, neutrinoPFPs);

    if (neutrinoPFPs.size() > 1)
    {
        throw cet::exception("SigmaRecoAnalyser") << "Too many neutrinos found!" << std::endl;
    }
    else if (neutrinoPFPs.size() == 0)
    {
        if (m_debugMode) std::cout << "No neutrino found!! :(" << std::endl;
    }
    else
    {
        // Slice score
        art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfpHandle, evt, pfpInputTag);
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(neutrinoPFPs[0].key());

        if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("NuScore") != pfpMetadata[0]->GetPropertiesMap().end()))
            m_recoNuSliceScore = pfpMetadata[0]->GetPropertiesMap().at("NuScore");

        // Slice ID
        art::FindManyP<recob::Slice> sliceAssoc = art::FindManyP<recob::Slice>(pfpHandle, evt, "pandora");
        std::vector<art::Ptr<recob::Slice>> slice = sliceAssoc.at(neutrinoPFPs[0].key());
        m_recoNuSliceID = slice[0]->ID();

        // Nu vertex
        art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfpHandle, evt, pfpInputTag);
        std::vector<art::Ptr<recob::Vertex>> vertex = vertexAssoc.at(neutrinoPFPs[0].key());

        if (!vertex.empty())
        {
            m_foundRecoNuVertex = true;
            m_recoNuVertexX = vertex[0]->position().X();
            m_recoNuVertexY = vertex[0]->position().Y();
            m_recoNuVertexZ = vertex[0]->position().Z();
        }

        if (m_debugMode)
        {
            std::cout << "m_recoNuSliceScore: " << m_recoNuSliceScore << std::endl;
            std::cout << "m_recoNuSliceID: " << m_recoNuSliceID << std::endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_tree = tfs->make<TTree>("sigmaAna", "");

    //EventID
    m_tree->Branch("EventID", &m_eventID);
    m_tree->Branch("Run", &m_run);
    m_tree->Branch("SubRun", &m_subRun);
    m_tree->Branch("Event", &m_event);

    // Truth stuff (event level) 
    m_tree->Branch("FoundTrueNuVertexSCE", &m_foundTrueNuVertexSCE);
    m_tree->Branch("TrueNuVertexSCEX", &m_trueNuVertexSCEX);
    m_tree->Branch("TrueNuVertexSCEY", &m_trueNuVertexSCEY);
    m_tree->Branch("TrueNuVertexSCEZ", &m_trueNuVertexSCEZ);
    m_tree->Branch("TrueProtonPiOpeningAngle", &m_trueProtonPiOpeningAngle);
    m_tree->Branch("TrueGammaLambdaOpeningAngle", &m_trueGammaLambdaOpeningAngle);
    m_tree->Branch("TrueGammaLambdaVertexSep", &m_trueGammaLambdaVertexSep);

    // Truth matching stuff
    m_tree->Branch("MCTruthIndex", &m_mcTruthIndex);    
    m_tree->Branch("TrueParticleID", &m_trueParticleID);
    m_tree->Branch("BestMatchedID", &m_bestMatchedID);
    m_tree->Branch("MatchIDs", &m_matchIDs);
    m_tree->Branch("TrueNuSliceID", &m_trueNuSliceID);
    m_tree->Branch("TrueNuSliceNHits", &m_trueNuSliceNHits);

    // Truth stuff (particle level)
    m_tree->Branch("TruePDG", &m_truePDG);
    m_tree->Branch("NTrueHits", &m_nTrueHits);
    m_tree->Branch("TrueNuVertexSep", &m_trueNuVertexSep);

    // Reco stuff (event level)
    m_tree->Branch("RecoNuSliceScore", &m_recoNuSliceScore);
    m_tree->Branch("RecoNuSliceID", &m_recoNuSliceID);
    m_tree->Branch("FoundRecoNuVertex", &m_foundRecoNuVertex);
    m_tree->Branch("RecoNuVertexX", &m_recoNuVertexX);
    m_tree->Branch("RecoNuVertexY", &m_recoNuVertexY);
    m_tree->Branch("RecoNuVertexZ", &m_recoNuVertexZ);
    m_tree->Branch("ProtonPionMerged", &m_protonPionMerged);

    // Reco stuff (particle level)
    m_tree->Branch("NMatches", &m_nMatches);
    m_tree->Branch("MatchFoundInSlice", &m_matchFoundInSlice);
    m_tree->Branch("BestMatchCompleteness", &m_bestMatchCompleteness);
    m_tree->Branch("BestMatchPurity", &m_bestMatchPurity);
    m_tree->Branch("BestMatchTrackScore", &m_bestMatchTrackScore);
    m_tree->Branch("BestMatchGeneration", &m_bestMatchGeneration);
    m_tree->Branch("BestMatchSliceScore", &m_bestMatchSliceScore);
    m_tree->Branch("BestMatchRecoVertexX", &m_bestMatchRecoVertexX);
    m_tree->Branch("BestMatchRecoVertexY", &m_bestMatchRecoVertexY);
    m_tree->Branch("BestMatchRecoVertexZ", &m_bestMatchRecoVertexZ);
    m_tree->Branch("BestMatchNuVertexSep", &m_bestMatchNuVertexSep);
    m_tree->Branch("MatchFoundInOtherSlice", &m_matchFoundInOtherSlice);
    m_tree->Branch("NParticlesFoundInSlice", &m_nParticlesFoundInSlice);
    m_tree->Branch("NParticlesFoundInOtherSlice", &m_nParticlesFoundInOtherSlice);
    m_tree->Branch("NParticlesFoundOverall", &m_nParticlesFoundOverall);
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::Reset()
{
    m_trueHitMap.clear();
    m_pfpMap.clear();
    m_mcParticleMap.clear();
    m_nuSliceMatchingMap.clear();
    m_otherSliceMatchingMap.clear();

    // ID info    
    m_mcTruthIndex = DEFAULT_INT;
    m_trueParticleID.clear();
    m_trueParticleID.resize(4, DEFAULT_INT);
    m_bestMatchedID.clear();
    m_bestMatchedID.resize(4, DEFAULT_INT);
    m_matchIDs.clear();
    m_matchIDs.resize(4, std::vector<int>());
    m_trueNuSliceID = DEFAULT_INT;
    m_trueNuSliceNHits = DEFAULT_INT;

    //EventID
    m_eventID = DEFAULT_INT;
    m_run = DEFAULT_INT;
    m_subRun = DEFAULT_INT;
    m_event = DEFAULT_INT;

    // Truth stuff (event level) 
    m_foundTrueNuVertexSCE = DEFAULT_BOOL;
    m_trueNuVertexSCEX = DEFAULT_DOUBLE;
    m_trueNuVertexSCEY = DEFAULT_DOUBLE;
    m_trueNuVertexSCEZ = DEFAULT_DOUBLE;
    m_trueProtonPiOpeningAngle = DEFAULT_DOUBLE;
    m_trueGammaLambdaOpeningAngle = DEFAULT_DOUBLE;
    m_trueGammaLambdaVertexSep = DEFAULT_DOUBLE;

    // Truth stuff (particle level)
    m_truePDG.clear();
    m_truePDG.resize(4, DEFAULT_INT);
    m_nTrueHits.clear();
    m_nTrueHits.resize(4, DEFAULT_INT);
    m_trueNuVertexSep.clear();
    m_trueNuVertexSep.resize(4, DEFAULT_DOUBLE);

    // Reco stuff (event level)
    m_recoNuSliceScore = DEFAULT_DOUBLE;
    m_recoNuSliceID = DEFAULT_DOUBLE;
    m_foundRecoNuVertex = DEFAULT_BOOL;
    m_recoNuVertexX = DEFAULT_DOUBLE;
    m_recoNuVertexY = DEFAULT_DOUBLE;
    m_recoNuVertexZ = DEFAULT_DOUBLE;
    m_protonPionMerged = DEFAULT_BOOL;
    m_isPionContaminatedByProton = DEFAULT_BOOL;
    m_isProtonContaminatedByPion = DEFAULT_BOOL;

    // Reco stuff (particle level)
    m_nMatches.clear();
    m_nMatches.resize(4, DEFAULT_INT);
    m_matchFoundInSlice.clear();
    m_matchFoundInSlice.resize(4, DEFAULT_INT);
    m_bestMatchCompleteness.clear();
    m_bestMatchCompleteness.resize(4, DEFAULT_DOUBLE);
    m_bestMatchPurity.clear();
    m_bestMatchPurity.resize(4, DEFAULT_DOUBLE);
    m_bestMatchTrackScore.clear();
    m_bestMatchTrackScore.resize(4, DEFAULT_DOUBLE);
    m_bestMatchGeneration.clear();
    m_bestMatchGeneration.resize(4, DEFAULT_INT);
    m_bestMatchSliceScore.clear();
    m_bestMatchSliceScore.resize(4, DEFAULT_DOUBLE);
    m_bestMatchRecoVertexX.clear();
    m_bestMatchRecoVertexX.resize(4, DEFAULT_DOUBLE);
    m_bestMatchRecoVertexY.clear();
    m_bestMatchRecoVertexY.resize(4, DEFAULT_DOUBLE);
    m_bestMatchRecoVertexZ.clear();
    m_bestMatchRecoVertexZ.resize(4, DEFAULT_DOUBLE);
    m_bestMatchNuVertexSep.clear();
    m_bestMatchNuVertexSep.resize(4, DEFAULT_DOUBLE);
    m_matchFoundInOtherSlice.clear();
    m_matchFoundInOtherSlice.resize(4, DEFAULT_INT);
    m_nParticlesFoundInSlice = DEFAULT_INT;
    m_nParticlesFoundInOtherSlice = DEFAULT_INT;
    m_nParticlesFoundOverall = DEFAULT_INT;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(hyperon::SigmaRecoAnalyser)
