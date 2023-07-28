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
    bool IsSignalEvent(const GeneratorTruth & genTruth, const G4Truth &G4Truth);
    void FillPandoraMaps(art::Event const& evt);
    bool IdentifySignalParticles(const G4Truth &G4Truth);
    bool IdentifySignalParticle(const G4Truth &G4Truth, const int index);
    void FillMCParticleTopologyInfo(const GeneratorTruth & genTruth, const G4Truth &G4Truth);
    void FillMCParticleHitInfo(art::Event const& evt);
    bool IsEM(const art::Ptr<simb::MCParticle> &mcParticle);
    int GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle);
    void PerformMatching(art::Event const& evt);
    void FindMCParticleMatches(art::Event const& evt, const int particleIndex,
        const std::vector<art::Ptr<recob::PFParticle>> &pfpVector);
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

    // Debug Info?
    bool m_debugMode;

    // ID info 
    int m_mcTruthIndex;
    std::vector<int> m_trueParticleID;
    std::vector<int> m_bestMatchedID;
    std::vector<std::vector<int>> m_matchIDs;

    // Analyser tree
    TTree * m_tree;

    // Event ID
    unsigned int m_eventID;
    int m_run;
    int m_subRun;
    int m_event;

    // Truth stuff (event level) 
    double m_trueProtonPiOpeningAngle;
    double m_trueGammaLambdaOpeningAngle;
    double m_trueGammaLambdaVertexSep;

    // Truth stuff (particle level)
    std::vector<int> m_truePDG;
    std::vector<int> m_nTrueHits;
    std::vector<double> m_trueNuVertexSep;

    // Reco stuff (event level)
    double m_sliceNuScore;
    bool m_foundRecoNuVertex;
    double m_recoNuVertexX;
    double m_recoNuVertexY;
    double m_recoNuVertexZ;
    bool m_protonPionMerged;

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

    // Linking TrackID -> MCParticle
    lar_pandora::MCParticleMap m_mcParticleMap;

    // Linking Self() -> PFParticle
    lar_pandora::PFParticleMap m_pfpMap;

    // Matching map
    MatchingMap m_matchingMap;
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

            int trackID = matchedMCParticle->TrackId();

            if (trackID == m_trueParticleID[MUON_INDEX])
                m_nTrueHits[MUON_INDEX] = (m_nTrueHits[MUON_INDEX] == DEFAULT_INT ? 1 : ++m_nTrueHits[MUON_INDEX]);
            else if (trackID == m_trueParticleID[PROTON_INDEX])
                m_nTrueHits[PROTON_INDEX] = (m_nTrueHits[PROTON_INDEX] == DEFAULT_INT ? 1 : ++m_nTrueHits[PROTON_INDEX]);
            else if (trackID == m_trueParticleID[PION_INDEX])
                m_nTrueHits[PION_INDEX] = (m_nTrueHits[PION_INDEX] == DEFAULT_INT ? 1 : ++m_nTrueHits[PION_INDEX]);
            else if (IsEM(matchedMCParticle) && (GetLeadEMTrackID(matchedMCParticle) == m_trueParticleID[GAMMA_INDEX]))
                m_nTrueHits[GAMMA_INDEX] = (m_nTrueHits[GAMMA_INDEX] == DEFAULT_INT ? 1 : ++m_nTrueHits[GAMMA_INDEX]);
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

void hyperon::SigmaRecoAnalyser::PerformMatching(art::Event const& evt)
{
    // Get all reco particles
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    art::InputTag pfpInputTag(m_PFParticleLabel, m_PFParticleInstance);

    if (!evt.getByLabel(pfpInputTag, pfpHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector, pfpHandle);

    // Fill matching map
    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        FindMCParticleMatches(evt, particleTypeIndex, pfpVector);

    // Order the matching map (in order of highest -> lowest shared hits)
    for (auto &entry : m_matchingMap)
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
        std::cout << "Matching map: " << std::endl;

        for (auto &entry : m_matchingMap)
        {
            std::cout << "index: " << entry.first << std::endl;
            std::cout << entry.second.size() << " matches" << std::endl;

            for (auto &pair : entry.second)
                std::cout << "(PFParticle ID: " << pair.first << ", shared hits: " << pair.second << std::endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FindMCParticleMatches(art::Event const& evt, const int particleIndex, 
    const std::vector<art::Ptr<recob::PFParticle>> &pfpVector)
{
    // Get event hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVector;

    if (!evt.getByLabel(m_HitLabel, hitHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No Hit Data Products Found!" << std::endl;

    art::fill_ptr_vector(hitVector, hitHandle);

    // Get backtracker info
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assocMCPart = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitHandle, evt, m_BacktrackLabel);

    for (const art::Ptr<recob::PFParticle> &pfp : pfpVector)
    {
        std::map<int, int> trackIDToHitMap;

        std::vector<art::Ptr<recob::Hit>> pfpHits;
        CollectHitsFromClusters(evt, pfp, pfpHits);

        for (const art::Ptr<recob::Hit> pfpHit : pfpHits)
        {
            const std::vector<art::Ptr<simb::MCParticle>> &matchedMCParticles = assocMCPart.at(pfpHit.key());
            auto matchedDatas = assocMCPart.data(pfpHit.key());

            for (unsigned int mcParticleIndex = 0; mcParticleIndex < matchedMCParticles.size(); ++mcParticleIndex)
            {
                const art::Ptr<simb::MCParticle> &matchedMCParticle = matchedMCParticles.at(mcParticleIndex);
                auto matchedData = matchedDatas.at(mcParticleIndex);

                if (matchedData->isMaxIDE != 1)
                    continue;

                const int matchedTrackID = IsEM(matchedMCParticle) ? GetLeadEMTrackID(matchedMCParticle) : matchedMCParticle->TrackId();

                trackIDToHitMap[matchedTrackID]++;
            }
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
            m_matchingMap[particleIndex].push_back(std::pair<int, int>(pfp->Self(), trackIDToHitMap[m_trueParticleID[particleIndex]]));

        if (particleIndex == PROTON_INDEX)
        {
            if ((maxOwnerID != PROTON_INDEX) && (maxOwnerID == PION_INDEX))
                m_protonPionMerged = true;
        }

        if (particleIndex == PION_INDEX)
        {
            if ((maxOwnerID != PION_INDEX) && (maxOwnerID == PROTON_INDEX))
                m_protonPionMerged = true;
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

void hyperon::SigmaRecoAnalyser::FillMatchingInfo(art::Event const& evt)
{
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

    art::InputTag pfpInputTag(m_PFParticleLabel, m_PFParticleInstance);

    if (!evt.getByLabel(pfpInputTag, pfparticleHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfparticleHandle, evt, pfpInputTag);
    art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfparticleHandle, evt, pfpInputTag);

    m_nParticlesFoundInSlice = 0;
    m_nParticlesFoundInOtherSlice = 0;

    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        if (m_debugMode) std::cout << "--------------- particleTypeIndex: " << particleTypeIndex << std::endl;

        int nuSliceMatches = 0;
        int otherSliceMatches = 0;

        bool foundMatchInSlice = false;
        bool foundMatchInOtherSlice = false;

        const std::vector<std::pair<int, int>> &matchPairs(m_matchingMap[particleTypeIndex]);

        for (const std::pair<int, int> &matchPair : matchPairs)
        {
            const art::Ptr<recob::PFParticle> pfp(m_pfpMap.at(matchPair.first));
            const int parentPDG(lar_pandora::LArPandoraHelper::GetParentNeutrino(m_pfpMap, pfp));

            // Best match should come from neutrino slice (if it exists in nu slice)
            // Come from another slice if it exists elsewhere)
            if (std::abs(parentPDG) == 14)
                ++nuSliceMatches;
            else
                ++otherSliceMatches;

            const bool setBestSliceMatch(!foundMatchInSlice && (nuSliceMatches == 1));
            const bool setOtherSliceMatch(!foundMatchInOtherSlice && (nuSliceMatches == 0));

            if (!setBestSliceMatch && !setOtherSliceMatch)
                continue;

            // Now work out the best match variables 
            const double completeness(static_cast<double>(matchPair.second) / static_cast<double>(m_nTrueHits[particleTypeIndex]));

            std::vector<art::Ptr<recob::Hit>> pfpHits;
            CollectHitsFromClusters(evt, pfp, pfpHits);

            const double purity(static_cast<double>(matchPair.second) / static_cast<double>(pfpHits.size()));

            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(pfp.key());

            double trackScore(DEFAULT_DOUBLE);

            if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("TrackScore") != pfpMetadata[0]->GetPropertiesMap().end()))
                trackScore = pfpMetadata[0]->GetPropertiesMap().at("TrackScore");

            double sliceScore(GetSliceScore(evt, pfp));
            int generation(lar_pandora::LArPandoraHelper::GetGeneration(m_pfpMap, pfp));

            std::vector<art::Ptr<recob::Vertex>> vertex = vertexAssoc.at(pfp.key());

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

            if (setBestSliceMatch)
                foundMatchInSlice = true;
            else
                foundMatchInOtherSlice = true;

            m_bestMatchCompleteness[particleTypeIndex] = completeness;
            m_bestMatchPurity[particleTypeIndex] = purity;
            m_bestMatchTrackScore[particleTypeIndex] = trackScore;
            m_bestMatchGeneration[particleTypeIndex] = generation;
            m_bestMatchSliceScore[particleTypeIndex] = sliceScore;
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
        std::cout << "NO NEUTRINO FOUND" << std::endl;
        m_sliceNuScore = -1;
    }
    else
    {
        art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfpHandle, evt, pfpInputTag);
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(neutrinoPFPs[0].key());

        if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("NuScore") != pfpMetadata[0]->GetPropertiesMap().end()))
            m_sliceNuScore = pfpMetadata[0]->GetPropertiesMap().at("NuScore");

        std::cout << "m_sliceNuScore: " << m_sliceNuScore << std::endl;

        // Vertices
        art::FindManyP<recob::Vertex> vertexAssoc = art::FindManyP<recob::Vertex>(pfpHandle, evt, pfpInputTag);
        std::vector<art::Ptr<recob::Vertex>> vertex = vertexAssoc.at(neutrinoPFPs[0].key());

        if (!vertex.empty())
        {
            m_foundRecoNuVertex = true;
            m_recoNuVertexX = vertex[0]->position().X();
            m_recoNuVertexY = vertex[0]->position().Y();
            m_recoNuVertexZ = vertex[0]->position().Z();
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
    m_tree->Branch("TrueProtonPiOpeningAngle", &m_trueProtonPiOpeningAngle);
    m_tree->Branch("TrueGammaLambdaOpeningAngle", &m_trueGammaLambdaOpeningAngle);
    m_tree->Branch("TrueGammaLambdaVertexSep", &m_trueGammaLambdaVertexSep);

    // Truth matching stuff
    m_tree->Branch("MCTruthIndex", &m_mcTruthIndex);    
    m_tree->Branch("TrueParticleID", &m_trueParticleID);
    m_tree->Branch("BestMatchedID", &m_bestMatchedID);
    m_tree->Branch("MatchIDs", &m_matchIDs);

    // Truth stuff (particle level)
    m_tree->Branch("TruePDG", &m_truePDG);
    m_tree->Branch("NTrueHits", &m_nTrueHits);
    m_tree->Branch("TrueNuVertexSep", &m_trueNuVertexSep);

    // Reco stuff (event level)
    m_tree->Branch("SliceNuScore", &m_sliceNuScore);
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
    m_mcParticleMap.clear();
    m_matchingMap.clear();

    
    m_mcTruthIndex = DEFAULT_INT;
    m_trueParticleID.clear();
    m_trueParticleID.resize(4, DEFAULT_INT);
    m_bestMatchedID.clear();
    m_bestMatchedID.resize(4, DEFAULT_INT);
    m_matchIDs.clear();
    m_matchIDs.resize(4, std::vector<int>());

    //EventID
    m_eventID = DEFAULT_INT;
    m_run = DEFAULT_INT;
    m_subRun = DEFAULT_INT;
    m_event = DEFAULT_INT;

    // Truth stuff (event level) 
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
    m_sliceNuScore = DEFAULT_DOUBLE;
    m_foundRecoNuVertex = DEFAULT_BOOL;
    m_recoNuVertexX = DEFAULT_DOUBLE;
    m_recoNuVertexY = DEFAULT_DOUBLE;
    m_recoNuVertexZ = DEFAULT_DOUBLE;
    m_protonPionMerged = DEFAULT_BOOL;

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
