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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

const int DEFAULT_INT = -999;
const double DEFAULT_DOUBLE = -999.0;
const std::string DEFAULT_STRING = "DEFAULT"; 
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
    void FillEventRecoInfo(art::Event const& evt);

private:
    typedef std::map<int, std::vector<std::pair<int,int>>> MatchingMap; // [Particle Type Index, [PFP ID, nHits]]

    // Product labels
    fhicl::ParameterSet m_G4Labels;
    fhicl::ParameterSet m_generatorLabels;
    std::string m_MCParticleLabel;
    std::string m_PFParticleLabel;
    std::string m_HitLabel;
    std::string m_BacktrackLabel;
    std::string m_ClusterLabel;

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
    std::string m_mode;
    std::string m_CCNC;
    double m_trueProtonPiOpeningAngle;
    double m_trueGammaLambdaOpeningAngle;
    double m_trueGammaLambdaVertexSep;

    // Truth stuff (particle level)
    std::vector<int> m_truePDG;
    std::vector<int> m_nTrueHits;
    std::vector<double> m_trueNuVertexSep;

    // Reco stuff (event level)
    double m_sliceNuScore;

    // Reco stuff (particle level)
    std::vector<int> m_nMatches;
    std::vector<int> m_matchFoundInSlice;
    std::vector<double> m_bestMatchCompleteness;
    std::vector<double> m_bestMatchPurity;
    std::vector<double> m_bestMatchTrackScore;
    std::vector<int> m_bestMatchGeneration;
    std::vector<double> m_bestMatchSliceScore;

    std::vector<int> m_matchFoundInOtherSlice;
    std::vector<double> m_otherSliceCompleteness;
    std::vector<double> m_otherSlicePurity;
    std::vector<double> m_otherSliceTrackScore;
    std::vector<int> m_otherSliceGeneration;
    std::vector<double> m_otherSliceSliceScore;

    // Linking MCParticle -> TrackID
    lar_pandora::MCParticleMap m_mcParticleMap;
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
    m_HitLabel(pset.get<std::string>("HitLabel")),
    m_BacktrackLabel(pset.get<std::string>("BacktrackLabel")),
    m_ClusterLabel(pset.get<std::string>("ClusterLabel"))
{
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::analyze(art::Event const& evt)
{
    // Reset tree
    Reset();

    if (m_debugMode) std::cout << "New Event" << std::endl;

    m_eventID = evt.id().event();
    m_run = evt.run();
    m_subRun = evt.subRun();
    m_event = evt.event();

    // Get truth info
    if (m_debugMode) std::cout << "Getting G4 Info" << std::endl;

    SubModuleGeneratorTruth* generatorSubModule = new SubModuleGeneratorTruth(evt, m_generatorLabels);
    GeneratorTruth genTruth = generatorSubModule->GetGeneratorTruth();

    SubModuleG4Truth* G4SubModule = new SubModuleG4Truth(evt, m_G4Labels);
    G4Truth G4Truth = G4SubModule->GetG4Info();

    if (!IsSignalEvent(genTruth, G4Truth))
    {
        if (m_debugMode) std::cout << "Event is not signal, moving on..." << std::endl;
        return;
    }

    //////////////////////////////////////////
    // Identify truth IDs - should be its own function
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

    if (!evt.getByLabel(m_MCParticleLabel, mcParticleHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No MCParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(mcParticleVector, mcParticleHandle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, m_mcParticleMap);


    // Get all reco particles
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    if (!evt.getByLabel(m_PFParticleLabel, pfpHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector, pfpHandle);

    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector, m_pfpMap);
    //////////////////////////////////////////

    if (!IdentifySignalParticles(G4Truth))
        return;

    FillMCParticleTopologyInfo(genTruth, G4Truth);

    FillMCParticleHitInfo(evt);

    // Make true -> reco matches
    PerformMatching(evt);

    FillMatchingInfo(evt);

    FillEventRecoInfo(evt);

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

bool hyperon::SigmaRecoAnalyser::IdentifySignalParticles(const G4Truth &G4Truth)
{
    if (!IdentifySignalParticle(G4Truth, MUON_INDEX))
        return false;

    if (!IdentifySignalParticle(G4Truth, PROTON_INDEX))
        return false;

    if (!IdentifySignalParticle(G4Truth, PION_INDEX))
        return false;

    if (!IdentifySignalParticle(G4Truth, GAMMA_INDEX))
        return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::SigmaRecoAnalyser::IdentifySignalParticle(const G4Truth &G4Truth, const int index)
{
    const std::vector<SimParticle> &simParticles(index == MUON_INDEX ? G4Truth.Lepton : 
        index == PROTON_INDEX ? G4Truth.Decay : index == PION_INDEX ? G4Truth.Decay : G4Truth.SigmaZeroDecayPhoton);

    const int pdg(index == MUON_INDEX ? -13 : index == PROTON_INDEX ? 2212 : index == PION_INDEX ? -211 : 22);

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
        m_trueParticleID[index] = simParticle.ArtID;
        break;
    }

    return found;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillMCParticleTopologyInfo(const GeneratorTruth & genTruth, const G4Truth &G4Truth)
{
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
    const TVector3 protonVertex(protonMCParticle->Vx(), protonMCParticle->Vy(), protonMCParticle->Vz());
    const TVector3 pionVertex(pionMCParticle->Vx(), pionMCParticle->Vy(), pionMCParticle->Vz());

    // Gamma, end point is first energy deposit (who designs this?)
    const TVector3 gammaVertex({gammaMCParticle->EndX(), gammaMCParticle->EndY(), gammaMCParticle->EndZ()});

    m_trueGammaLambdaVertexSep = (protonVertex - gammaVertex).Mag();
    m_trueNuVertexSep[MUON_INDEX] = 0.0;
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
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::SigmaRecoAnalyser::IsEM(const art::Ptr<simb::MCParticle> &mcParticle)
{
    return ((std::abs(mcParticle->PdgCode()) == 11) || (mcParticle->PdgCode() == 22));
}

///////////////////////////////////////////////////////////////////////////////////////////
// If its an EM particle, we have to move up the EM hierarchy (who designs this, honestly)
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

    if (!evt.getByLabel(m_PFParticleLabel, pfpHandle))
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
        std::vector<art::Ptr<recob::Hit>> pfpHits;
        CollectHitsFromClusters(evt, pfp, pfpHits);

        int sharedHits = 0;

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

                const int matchedTrackID = matchedMCParticle->TrackId();

                if (matchedTrackID == m_trueParticleID[particleIndex])
                    sharedHits++;
            }
        }

        const double completeness = static_cast<double>(sharedHits) / static_cast<double>(m_nTrueHits[particleIndex]);
        const double purity = static_cast<double>(sharedHits) / static_cast<double>(pfpHits.size());

        if ((completeness > COMPLETENESS_THRESHOLD) && (purity > PURITY_THRESHOLD))
            m_matchingMap[particleIndex].push_back(std::pair<int, int>(pfp->Self(), sharedHits));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
   std::vector<art::Ptr<recob::Hit>> &hits)
{
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

   if (!evt.getByLabel(m_PFParticleLabel, pfparticleHandle))
       throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

   art::Handle<std::vector<recob::Cluster>> clusterHandle;

   if (!evt.getByLabel(m_ClusterLabel, clusterHandle)) 
       throw cet::exception("SigmaRecoAnalyser") << "No Cluster Data Products Found!" << std::endl;

   art::FindManyP<recob::Cluster> pfparticleClustersAssoc = art::FindManyP<recob::Cluster>(pfparticleHandle, evt, m_PFParticleLabel);
   art::FindManyP<recob::Hit> clusterHitAssoc = art::FindManyP<recob::Hit>(clusterHandle, evt, m_ClusterLabel);

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

   if (!evt.getByLabel(m_PFParticleLabel, pfparticleHandle))
       throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

   art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfparticleHandle, evt, m_PFParticleLabel);

    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        bool foundMatchInSlice = false;
        bool foundMatchInOtherSlice = false;

        const std::vector<std::pair<int, int>> &matchPairs(m_matchingMap[particleTypeIndex]);

        for (const std::pair<int, int> &matchPair : matchPairs)
        {
            const art::Ptr<recob::PFParticle> pfp(m_pfpMap.at(matchPair.first));
            const int parentPDG(lar_pandora::LArPandoraHelper::GetParentNeutrino(m_pfpMap, pfp));

            // If in neutrino slice...
            if (std::abs(parentPDG) == 14)
                ++m_nMatches[particleTypeIndex];

            const bool setBestSliceMatch(!foundMatchInSlice && (m_nMatches[particleTypeIndex] == 1));
            const bool setOtherSliceMatch(!foundMatchInOtherSlice && (m_nMatches[particleTypeIndex] == 0));

            if (!setBestSliceMatch && !setOtherSliceMatch)
                continue;

            const double completeness(static_cast<double>(matchPair.second) / static_cast<double>(m_nTrueHits[particleTypeIndex]));

            std::vector<art::Ptr<recob::Hit>> pfpHits;
            CollectHitsFromClusters(evt, pfp, pfpHits);

            const double purity(static_cast<double>(matchPair.second) / static_cast<double>(pfpHits.size()));

            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(pfp.key());

            double trackScore(DEFAULT_DOUBLE);

            if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("TrackScore") != pfpMetadata[0]->GetPropertiesMap().end()))
                trackScore = pfpMetadata[0]->GetPropertiesMap().at("TrackScore");

            double sliceScore(DEFAULT_DOUBLE);

            if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("NuScore") != pfpMetadata[0]->GetPropertiesMap().end()))
                sliceScore = pfpMetadata[0]->GetPropertiesMap().at("NuScore");

            int generation(lar_pandora::LArPandoraHelper::GetGeneration(m_pfpMap, pfp));

            if (setBestSliceMatch)
            {
                m_bestMatchCompleteness[particleTypeIndex] = completeness;
                m_bestMatchPurity[particleTypeIndex] = purity;
                m_bestMatchTrackScore[particleTypeIndex] = trackScore;
                m_bestMatchGeneration[particleTypeIndex] = generation;
                m_bestMatchSliceScore[particleTypeIndex] = sliceScore;
                foundMatchInSlice = true;
            }
            else
            {
                m_otherSliceCompleteness[particleTypeIndex] = completeness;
                m_otherSlicePurity[particleTypeIndex] = purity;
                m_otherSliceTrackScore[particleTypeIndex] = trackScore;
                m_otherSliceGeneration[particleTypeIndex] = generation;
                m_otherSliceSliceScore[particleTypeIndex] = sliceScore;
                foundMatchInOtherSlice = true;
            }
        }

        m_matchFoundInSlice[particleTypeIndex] = foundMatchInSlice;
        m_matchFoundInOtherSlice[particleTypeIndex] = foundMatchInOtherSlice;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillEventRecoInfo(art::Event const& evt)
{
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector;

    if (!evt.getByLabel(m_PFParticleLabel, pfpHandle))
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
        m_sliceNuScore = -1;
    }
    else
    {
        art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfpHandle, evt, m_PFParticleLabel);
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn.at(neutrinoPFPs[0].key());

        if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("NuScore") != pfpMetadata[0]->GetPropertiesMap().end()))
            m_sliceNuScore = pfpMetadata[0]->GetPropertiesMap().at("NuScore");
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
    m_tree->Branch("Mode", &m_mode);
    m_tree->Branch("CCNC", &m_CCNC);
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

    // Reco stuff (particle level)
    m_tree->Branch("NMatches", &m_nMatches);
    m_tree->Branch("MatchFoundInSlice", &m_matchFoundInSlice);
    m_tree->Branch("BestMatchCompleteness", &m_bestMatchCompleteness);
    m_tree->Branch("BestMatchPurity", &m_bestMatchPurity);
    m_tree->Branch("BestMatchTrackScore", &m_bestMatchTrackScore);
    m_tree->Branch("BestMatchGeneration", &m_bestMatchGeneration);
    m_tree->Branch("BestMatchSliceScore", &m_bestMatchSliceScore);
    m_tree->Branch("MatchFoundInOtherSlice", &m_matchFoundInOtherSlice);
    m_tree->Branch("OtherSliceCompleteness", &m_otherSliceCompleteness);
    m_tree->Branch("OtherSlicePurity", &m_otherSlicePurity);
    m_tree->Branch("OtherSliceTrackScore", &m_otherSliceTrackScore);
    m_tree->Branch("OtherSliceGeneration", &m_otherSliceGeneration);
    m_tree->Branch("OtherSliceSliceScore", &m_otherSliceSliceScore);
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
    m_mode = DEFAULT_STRING;
    m_CCNC = DEFAULT_STRING;
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
    m_matchFoundInOtherSlice.clear();
    m_matchFoundInOtherSlice.resize(4, DEFAULT_INT);
    m_otherSliceCompleteness.clear();
    m_otherSliceCompleteness.resize(4, DEFAULT_DOUBLE);
    m_otherSlicePurity.clear();
    m_otherSlicePurity.resize(4, DEFAULT_DOUBLE);
    m_otherSliceTrackScore.clear();
    m_otherSliceTrackScore.resize(4, DEFAULT_DOUBLE);
    m_otherSliceGeneration.clear();
    m_otherSliceGeneration.resize(4, DEFAULT_INT);
    m_otherSliceSliceScore.clear();
    m_otherSliceSliceScore.resize(4, DEFAULT_DOUBLE);
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(hyperon::SigmaRecoAnalyser)
