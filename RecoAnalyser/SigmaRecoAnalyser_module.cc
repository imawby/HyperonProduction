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
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

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

    typedef std::map<int, std::vector<std::pair<int,int>>> MatchingMap; // [Particle Type Index, [PFP ID, nHits]]

    void Reset();
    //void InvestigateSlices(art::Event const& evt);
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
    void FindMCParticleMatches(art::Event const& evt, const int particleIndex, const int nuSliceID,
        MatchingMap &nuSliceMatchingMap, MatchingMap &nuSliceMatchingMap_Repass, MatchingMap &otherSliceMatchingMap);
    void CollectHitsFromClusters(art::Event const& evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        std::vector<art::Ptr<recob::Hit>> &hits);
    void OrderMatchingMap(MatchingMap &matchingMap);
    void FillMatchingInfo(art::Event const& evt, const bool isTrueNuSlice);
    void FillContaminationInfo(art::Event const& evt, const bool isTrueNuSlice);
    void FillEventRecoInfo(art::Event const& evt);

private:
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
    std::string m_ModifiedPFParticleLabel;
    std::string m_RepassPFParticleLabel;

    // Does input include modified reco?
    bool m_modifiedRecoInput;

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
    std::vector<int> m_nTrueHitsU;
    std::vector<int> m_nTrueHitsV;
    std::vector<int> m_nTrueHitsW;
    std::vector<double> m_trueNuVertexSep;

    // Reco stuff (event level)
    double m_recoNuSliceScore;
    int m_recoNuSliceID;
    bool m_foundCorrectNuSlice;
    bool m_foundRecoNuVertex;
    double m_recoNuVertexX;
    double m_recoNuVertexY;
    double m_recoNuVertexZ;

    // Reco stuff - match found?
    std::vector<int> m_matchFoundInTrueNuSlice;
    std::vector<int> m_matchFoundInTrueOtherSlice;
    std::vector<int> m_matchFoundInRecoNuSlice;
    std::vector<int> m_matchFoundInRecoOtherSlice;

    int m_nParticlesFoundInTrueNuSlice;
    int m_nParticlesFoundInTrueOtherSlice;
    int m_nParticlesFoundTrueOverall;

    int m_nParticlesFoundInRecoNuSlice;
    int m_nParticlesFoundInRecoOtherSlice;
    int m_nParticlesFoundRecoOverall;

    bool m_pionMergedTrueSlice;
    int m_pionContaminantHitsTrueSlice;
    bool m_protonMergedTrueSlice;
    int m_protonContaminantHitsTrueSlice;

    bool m_pionMergedRecoSlice;
    int m_pionContaminantHitsRecoSlice;
    bool m_protonMergedRecoSlice;
    int m_protonContaminantHitsRecoSlice;

    // Reco stuff - best match stuff
    std::vector<bool> m_trueBestMatchRecoRepass;
    std::vector<double> m_trueBestMatchCompleteness;
    std::vector<double> m_trueBestMatchPurity;
    std::vector<double> m_trueBestMatchTrackScore;
    std::vector<int> m_trueBestMatchGeneration;
    std::vector<double> m_trueBestMatchSliceID;
    std::vector<double> m_trueBestMatchRecoVertexX;
    std::vector<double> m_trueBestMatchRecoVertexY;
    std::vector<double> m_trueBestMatchRecoVertexZ;
    std::vector<double> m_trueBestMatchNuVertexSep;

    std::vector<bool> m_recoBestMatchRecoRepass;
    std::vector<double> m_recoBestMatchCompleteness;
    std::vector<double> m_recoBestMatchPurity;
    std::vector<double> m_recoBestMatchTrackScore;
    std::vector<int> m_recoBestMatchGeneration;
    std::vector<double> m_recoBestMatchSliceID;
    std::vector<double> m_recoBestMatchRecoVertexX;
    std::vector<double> m_recoBestMatchRecoVertexY;
    std::vector<double> m_recoBestMatchRecoVertexZ;
    std::vector<double> m_recoBestMatchNuVertexSep;

    // Linking particleTypeIndex -> list of hit.key()
    typedef std::map<int, int> TrueHitMap;
    TrueHitMap m_trueHitMap;

    // Linking TrackID -> MCParticle
    lar_pandora::MCParticleMap m_mcParticleMap;

    // Linking Self() -> PFParticle
    lar_pandora::PFParticleMap m_pfpMap;
    lar_pandora::PFParticleMap m_pfpMap_AO;
    lar_pandora::PFParticleMap m_pfpMap_Repass;

    // Linking slice ID -> slice key
    std::map<int, int> m_sliceMap_AO;
    std::map<int, int> m_sliceMap_Repass;

    // Matching map
    MatchingMap m_trueNuSliceMatchingMap;
    MatchingMap m_trueNuSliceMatchingMap_Repass;
    MatchingMap m_trueOtherSliceMatchingMap;
    MatchingMap m_recoNuSliceMatchingMap;
    MatchingMap m_recoNuSliceMatchingMap_Repass;
    MatchingMap m_recoOtherSliceMatchingMap;
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
    m_ModifiedPFParticleLabel(pset.get<std::string>("ModifiedPFParticleLabel", "")),
    m_RepassPFParticleLabel(pset.get<std::string>("RepassPFParticleLabel", "")),
    m_modifiedRecoInput(pset.get<bool>("ModifiedRecoInput")),
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
    //InvestigateSlices(evt);

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
    FillMatchingInfo(evt, true);
    FillMatchingInfo(evt, false);

    FillContaminationInfo(evt, true);
    FillContaminationInfo(evt, false);

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

    // PFParticle map - all outcomes
    art::InputTag pfpInputTag_AO("pandoraPatRec", "allOutcomes");
    art::Handle<std::vector<recob::PFParticle>> pfpHandle_AO;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector_AO;

    if (!evt.getByLabel(pfpInputTag_AO, pfpHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(pfpVector_AO, pfpHandle_AO);

    lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector_AO, m_pfpMap_AO);

    if (m_debugMode) std::cout << "m_pfpMap_AO.size(): " << m_pfpMap_AO.size() << std::endl;

    // Slice map - all outcomes
    art::InputTag sliceInputTag_AO("pandoraPatRec", "allOutcomes");
    art::Handle<std::vector<recob::Slice>> sliceHandle_AO;
    std::vector<art::Ptr<recob::Slice>> sliceVector_AO;

    if (!evt.getByLabel(sliceInputTag_AO, sliceHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No Slice Data Products Found!" << std::endl;

    art::fill_ptr_vector(sliceVector_AO, sliceHandle_AO);

    for (art::Ptr<recob::Slice> &slice_AO : sliceVector_AO)
        m_sliceMap_AO[slice_AO->ID()] = slice_AO.key();

    // Repass PFP Map
    if (m_modifiedRecoInput)
    {
        art::Handle<std::vector<recob::PFParticle>> pfpHandle_Repass;
        std::vector<art::Ptr<recob::PFParticle>> pfpVector_Repass;

        art::InputTag pfpInputTag_Repass("pandoraJam", "");

        if (!evt.getByLabel(pfpInputTag_Repass, pfpHandle_Repass))
            throw cet::exception("SigmaRecoAnalyser") << "No Repass PFParticle Data Products Found!" << std::endl;

        art::fill_ptr_vector(pfpVector_Repass, pfpHandle_Repass);

        lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfpVector_Repass, m_pfpMap_Repass);

        if (m_debugMode) std::cout << "m_pfpMap_Repass.size(): " << m_pfpMap_Repass.size() << std::endl;
    }
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

    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        m_nTrueHits[particleTypeIndex] = 0;
        m_nTrueHitsW[particleTypeIndex] = 0;
        m_nTrueHitsU[particleTypeIndex] = 0;
        m_nTrueHitsV[particleTypeIndex] = 0;
    }

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

            // Get hit view
            const geo::WireID hit_WireID(hit->WireID());
            const geo::View_t hit_View(hit->View());
            const geo::View_t pandoraView(lar_pandora::LArPandoraGeometry::GetGlobalView(hit_WireID.Cryostat, hit_WireID.TPC, hit_View));
            const int trackID = IsEM(matchedMCParticle) ? GetLeadEMTrackID(matchedMCParticle) : matchedMCParticle->TrackId();

            m_trueHitMap[hit.key()] = trackID;

            for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
            {
                if (trackID == m_trueParticleID[particleTypeIndex])
                {
                    ++m_nTrueHits[particleTypeIndex];

                    if (pandoraView == geo::kW)
                        ++m_nTrueHitsW[particleTypeIndex];
                    else if (pandoraView == geo::kU)
                        ++m_nTrueHitsU[particleTypeIndex];
                    else if (pandoraView == geo::kV)
                        ++m_nTrueHitsV[particleTypeIndex];
                }
            }
        }
    }

    if (m_debugMode)
    {
        for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        {
            std::cout << "particleTypeIndex: " << particleTypeIndex << std::endl;
            std::cout << "nTrueHits: " << m_nTrueHits[particleTypeIndex] << std::endl;
            std::cout << "nTrueHitsU: " << m_nTrueHitsU[particleTypeIndex] << std::endl;
            std::cout << "nTrueHitsV: " << m_nTrueHitsV[particleTypeIndex] << std::endl;
            std::cout << "nTrueHitsW: " << m_nTrueHitsW[particleTypeIndex] << std::endl;
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
    {
        FindMCParticleMatches(evt, particleTypeIndex, m_trueNuSliceID, m_trueNuSliceMatchingMap, m_trueNuSliceMatchingMap_Repass, m_trueOtherSliceMatchingMap);
        FindMCParticleMatches(evt, particleTypeIndex, m_recoNuSliceID, m_recoNuSliceMatchingMap, m_recoNuSliceMatchingMap_Repass, m_recoOtherSliceMatchingMap);
    }

    OrderMatchingMap(m_trueNuSliceMatchingMap);
    OrderMatchingMap(m_trueNuSliceMatchingMap_Repass);
    OrderMatchingMap(m_trueOtherSliceMatchingMap);
    OrderMatchingMap(m_recoNuSliceMatchingMap);
    OrderMatchingMap(m_recoNuSliceMatchingMap_Repass);
    OrderMatchingMap(m_recoOtherSliceMatchingMap);
    
    if (m_debugMode)
    {
        std::cout << "Reco neutrino slice matching map: " << std::endl;

        for (auto &entry : m_recoNuSliceMatchingMap)
        {
            std::cout << "index: " << entry.first << std::endl;
            std::cout << entry.second.size() << " matches" << std::endl;

            for (auto &pair : entry.second)
                std::cout << "(PFParticle ID: " << pair.first << ", shared hits: " << pair.second << std::endl;
        }

        std::cout << "Reco neutrino slice matching map (Reco repass): " << std::endl;

        for (auto &entry : m_recoNuSliceMatchingMap_Repass)
        {
            std::cout << "index: " << entry.first << std::endl;
            std::cout << entry.second.size() << " matches" << std::endl;

            for (auto &pair : entry.second)
                std::cout << "(PFParticle ID: " << pair.first << ", shared hits: " << pair.second << std::endl;
        }

        std::cout << "Reco other slice matching map: " << std::endl;

        for (auto &entry : m_recoOtherSliceMatchingMap)
        {
            std::cout << "index: " << entry.first << std::endl;
            std::cout << entry.second.size() << " matches" << std::endl;

            for (auto &pair : entry.second)
                std::cout << "(PFParticle ID: " << pair.first << ", shared hits: " << pair.second << std::endl;
        }
    }
    
}


///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FindMCParticleMatches(art::Event const& evt, const int particleIndex, const int nuSliceID, 
    MatchingMap &nuSliceMatchingMap, MatchingMap &nuSliceMatchingMap_Repass, MatchingMap &otherSliceMatchingMap)
{
    // Slice vector - All outcomes
    art::InputTag sliceInputTag_AO("pandoraPatRec", "allOutcomes");
    art::Handle<std::vector<recob::Slice>> sliceHandle_AO;
    std::vector<art::Ptr<recob::Slice>> sliceVector_AO;

    if (!evt.getByLabel(sliceInputTag_AO, sliceHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::PFParticle> pfpAssoc_AO = art::FindManyP<recob::PFParticle>(sliceHandle_AO, evt, sliceInputTag_AO);
    art::fill_ptr_vector(sliceVector_AO, sliceHandle_AO);

    // Modifed reconstruction
    art::Handle<std::vector<int>> pfpHandle_Modified;

    if (m_modifiedRecoInput)
    {
        if (!evt.getByLabel("pandoraJam", pfpHandle_Modified))
            throw cet::exception("SigmaRecoAnalyser") << "No Modified PFParticle Data Products Found!" << std::endl;
    }

    // Repass Reconstruction 
    art::Handle<std::vector<recob::PFParticle>> pfpHandle_Repass;
    std::vector<art::Ptr<recob::PFParticle>> pfpVector_Repass;

    if (m_modifiedRecoInput)
    {
        art::InputTag pfpInputTag_Repass("pandoraJam", "");

        if (!evt.getByLabel(pfpInputTag_Repass, pfpHandle_Repass))
            throw cet::exception("SigmaRecoAnalyser") << "No Repass PFParticle Data Products Found!" << std::endl;

        art::fill_ptr_vector(pfpVector_Repass, pfpHandle_Repass);
    }

    // Search in the neutrino slice first...
    bool foundInNuSlice = false;

    const std::vector<art::Ptr<recob::PFParticle>> &nuSlicePFPs_AO(pfpAssoc_AO.at(m_sliceMap_AO[nuSliceID]));

    for (const art::Ptr<recob::PFParticle> &pfp_AO : nuSlicePFPs_AO)
    {
        // If not neutrino...
        if (!lar_pandora::LArPandoraHelper::IsNeutrino(lar_pandora::LArPandoraHelper::GetParentPFParticle(m_pfpMap_AO, pfp_AO)))
            continue;

        if (m_modifiedRecoInput)
        {
            bool modified = false;

            for (long unsigned int i = 0; i < pfpHandle_Modified.product()->size(); ++i)
            {
                if (static_cast<long unsigned int>(pfpHandle_Modified.product()->at(i)) == pfp_AO.key())
                {
                    modified = true;
                    break;
                }
            }

            if (modified)
                continue;
        }
        
        std::map<int, int> trackIDToHitMap;

        std::vector<art::Ptr<recob::Hit>> pfpHits;
        CollectHitsFromClusters(evt, pfp_AO, pfpHits, "pandoraPatRec", "allOutcomes");

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
            nuSliceMatchingMap[particleIndex].push_back(std::pair<int, int>(pfp_AO->Self(), trackIDToHitMap[m_trueParticleID[particleIndex]]));
            foundInNuSlice = true;
        }
    }

    // Search in the repass PFParticles...
    for (const art::Ptr<recob::PFParticle> &pfp_Repass : pfpVector_Repass)
    {
        std::map<int, int> trackIDToHitMap;

        std::vector<art::Ptr<recob::Hit>> pfpHits;
        CollectHitsFromClusters(evt, pfp_Repass, pfpHits, "pandoraJam", "");

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
            nuSliceMatchingMap_Repass[particleIndex].push_back(std::pair<int, int>(pfp_Repass->Self(), trackIDToHitMap[m_trueParticleID[particleIndex]]));
            foundInNuSlice = true;
        }
    }

    if (foundInNuSlice)
        return;

    //////////////////////////////
    // Now search in other slices!
    //////////////////////////////
    for (const art::Ptr<recob::Slice> &slice_AO : sliceVector_AO)
    {
        if (slice_AO->ID() == nuSliceID)
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
                otherSliceMatchingMap[particleIndex].push_back(std::pair<int, int>(pfp_AO->Self(), trackIDToHitMap[m_trueParticleID[particleIndex]]));
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

void hyperon::SigmaRecoAnalyser::OrderMatchingMap(MatchingMap &matchingMap)
{
    // Order the matching map (in order of highest -> lowest shared hits)
    for (auto &entry : matchingMap)
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

void hyperon::SigmaRecoAnalyser::FillMatchingInfo(art::Event const& evt, const bool isTrueNuSlice)
{
    // All outcomes output
    art::InputTag pfpInputTag_AO("pandoraPatRec", "allOutcomes");
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle_AO;

    if (!evt.getByLabel(pfpInputTag_AO, pfparticleHandle_AO))
        throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn_AO = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfparticleHandle_AO, evt, pfpInputTag_AO);
    art::FindManyP<recob::Vertex> vertexAssoc_AO = art::FindManyP<recob::Vertex>(pfparticleHandle_AO, evt, pfpInputTag_AO);
    art::FindManyP<recob::Slice> sliceAssoc_AO = art::FindManyP<recob::Slice>(pfparticleHandle_AO, evt, pfpInputTag_AO);

    // Repass output
    art::InputTag pfpInputTag_Repass("pandoraJam", "");
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle_Repass;

    if (m_modifiedRecoInput)
    {
        if (!evt.getByLabel(pfpInputTag_Repass, pfparticleHandle_Repass))
            throw cet::exception("SigmaRecoAnalyser") << "No PFParticle Data Products Found!" << std::endl;
    }

    int nParticlesFoundInNuSlice = 0;
    int nParticlesFoundInOtherSlice = 0;

    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        if (m_debugMode) std::cout << "--------------- particleTypeIndex: " << particleTypeIndex << std::endl;

        const std::vector<std::pair<int, int>> &nuSliceMatchPairs(isTrueNuSlice ? m_trueNuSliceMatchingMap[particleTypeIndex] : m_recoNuSliceMatchingMap[particleTypeIndex]);
        const std::vector<std::pair<int, int>> &nuSliceMatchPairs_Repass(isTrueNuSlice ? m_trueNuSliceMatchingMap_Repass[particleTypeIndex] : m_recoNuSliceMatchingMap_Repass[particleTypeIndex]);
        const std::vector<std::pair<int, int>> &otherSliceMatchPairs(isTrueNuSlice ? m_trueOtherSliceMatchingMap[particleTypeIndex] : m_recoOtherSliceMatchingMap[particleTypeIndex]);

        const int nuSliceMatches = nuSliceMatchPairs.size() + nuSliceMatchPairs_Repass.size();
        const int otherSliceMatches = otherSliceMatchPairs.size();

        bool foundMatchInNuSlice = (nuSliceMatches != 0);
        bool foundMatchInOtherSlice = (otherSliceMatches != 0);

        if (isTrueNuSlice)
        {
            m_matchFoundInTrueNuSlice[particleTypeIndex] = foundMatchInNuSlice;
            m_matchFoundInTrueOtherSlice[particleTypeIndex] = !foundMatchInNuSlice && foundMatchInOtherSlice;
        }
        else
        {
            m_matchFoundInRecoNuSlice[particleTypeIndex] = foundMatchInNuSlice;
            m_matchFoundInRecoOtherSlice[particleTypeIndex] = !foundMatchInNuSlice && foundMatchInOtherSlice;
        }

        if (m_debugMode)
        {
            std::cout << (isTrueNuSlice ? "matchFoundInTrueNuSlice: " : "matchFoundInRecoNuSlice: ") << 
                (isTrueNuSlice ? m_matchFoundInTrueNuSlice[particleTypeIndex] : m_matchFoundInRecoNuSlice[particleTypeIndex]) << std::endl;
            std::cout << (isTrueNuSlice ? "matchFoundInTrueOtherSlice: " : "matchFoundInRecoOtherSlice: ") << 
                (isTrueNuSlice ? m_matchFoundInTrueOtherSlice[particleTypeIndex] : m_matchFoundInRecoOtherSlice[particleTypeIndex]) << std::endl;
        }

        if ((nuSliceMatches == 0) && (otherSliceMatches == 0))
            continue;

        if (foundMatchInNuSlice)
            ++nParticlesFoundInNuSlice;

        if (!foundMatchInNuSlice && foundMatchInOtherSlice)
            ++nParticlesFoundInOtherSlice;

        const bool isRepassOutput((nuSliceMatchPairs_Repass.size() != 0) && ((nuSliceMatchPairs.size() == 0) || (nuSliceMatchPairs.front().second < nuSliceMatchPairs_Repass.front().second)));

        const std::pair<int, int> &matchPair(foundMatchInNuSlice ? (isRepassOutput ? nuSliceMatchPairs_Repass.front() : nuSliceMatchPairs.front()) : otherSliceMatchPairs.front());
        const art::Ptr<recob::PFParticle> pfp(isRepassOutput ? m_pfpMap_Repass.at(matchPair.first) : m_pfpMap_AO.at(matchPair.first));

        //////////////////////////////////////
        // Now fill out bestMatchPfo info
        //////////////////////////////////////
        // Completeness & Purity
        std::vector<art::Ptr<recob::Hit>> pfpHits;
        CollectHitsFromClusters(evt, pfp, pfpHits, isRepassOutput ? "pandoraJam" : "pandoraPatRec", isRepassOutput ? "" : "allOutcomes");

        const double completeness(static_cast<double>(matchPair.second) / static_cast<double>(m_nTrueHits[particleTypeIndex]));
        const double purity(static_cast<double>(matchPair.second) / static_cast<double>(pfpHits.size()));

        // Track score
        double trackScore(DEFAULT_DOUBLE);

        if (isRepassOutput)
        {
            art::FindManyP<larpandoraobj::PFParticleMetadata> metadataAssn_Repass = art::FindManyP<larpandoraobj::PFParticleMetadata>(pfparticleHandle_Repass, evt, pfpInputTag_Repass);

            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn_Repass.at(pfp.key());

            if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("TrackScore") != pfpMetadata[0]->GetPropertiesMap().end()))
                trackScore = pfpMetadata[0]->GetPropertiesMap().at("TrackScore");
        }
        else
        {
            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetadata = metadataAssn_AO.at(pfp.key());

            if (!pfpMetadata.empty() && (pfpMetadata[0]->GetPropertiesMap().find("TrackScore") != pfpMetadata[0]->GetPropertiesMap().end()))
                trackScore = pfpMetadata[0]->GetPropertiesMap().at("TrackScore");
        }

        // Slice ID
        int sliceID(DEFAULT_INT);

        if (isRepassOutput)
        {
            art::FindManyP<recob::Slice> sliceAssoc_Repass = art::FindManyP<recob::Slice>(pfparticleHandle_Repass, evt, pfpInputTag_Repass);
            std::vector<art::Ptr<recob::Slice>> slices = sliceAssoc_Repass.at(pfp.key());
            if (!slices.empty()) sliceID = slices[0]->ID();
        }
        else
        {
            std::vector<art::Ptr<recob::Slice>> slices = sliceAssoc_AO.at(pfp.key());
            if (!slices.empty()) sliceID = slices[0]->ID();
        }

        // Generation
        int generation(lar_pandora::LArPandoraHelper::GetGeneration(isRepassOutput ? m_pfpMap_Repass : m_pfpMap_AO, pfp));

        // Vertex info
        std::vector<art::Ptr<recob::Vertex>> vertex;

        if (isRepassOutput)
        {
            art::FindManyP<recob::Vertex> vertexAssoc_Repass = art::FindManyP<recob::Vertex>(pfparticleHandle_Repass, evt, pfpInputTag_Repass);
            vertex = vertexAssoc_Repass.at(pfp.key());
        }
        else
        {
            std::vector<art::Ptr<recob::Vertex>> vertex = vertexAssoc_AO.at(pfp.key());
        }

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

        if (isTrueNuSlice)
        {
            m_trueBestMatchRecoRepass[particleTypeIndex] = isRepassOutput;
            m_trueBestMatchCompleteness[particleTypeIndex] = completeness;
            m_trueBestMatchPurity[particleTypeIndex] = purity;
            m_trueBestMatchTrackScore[particleTypeIndex] = trackScore;
            m_trueBestMatchGeneration[particleTypeIndex] = generation;
            m_trueBestMatchSliceID[particleTypeIndex] = sliceID;
            m_trueBestMatchRecoVertexX[particleTypeIndex] = vertexX;
            m_trueBestMatchRecoVertexY[particleTypeIndex] = vertexY;
            m_trueBestMatchRecoVertexZ[particleTypeIndex] = vertexZ;
            m_trueBestMatchNuVertexSep[particleTypeIndex] = nuVertexSep;
        }
        else
        {
            m_recoBestMatchRecoRepass[particleTypeIndex] = isRepassOutput;
            m_recoBestMatchCompleteness[particleTypeIndex] = completeness;
            m_recoBestMatchPurity[particleTypeIndex] = purity;
            m_recoBestMatchTrackScore[particleTypeIndex] = trackScore;
            m_recoBestMatchGeneration[particleTypeIndex] = generation;
            m_recoBestMatchSliceID[particleTypeIndex] = sliceID;
            m_recoBestMatchRecoVertexX[particleTypeIndex] = vertexX;
            m_recoBestMatchRecoVertexY[particleTypeIndex] = vertexY;
            m_recoBestMatchRecoVertexZ[particleTypeIndex] = vertexZ;
            m_recoBestMatchNuVertexSep[particleTypeIndex] = nuVertexSep;
        }

        if (m_debugMode)
        {
            std::cout << (isTrueNuSlice ? "trueBestMatchRecoRepass: " : "recoBestMatchRecoRepass: ") << (isRepassOutput? "yes" : "no") << std::endl;
            std::cout << (isTrueNuSlice ? "trueBestMatchCompleteness: " : "recoBestMatchCompleteness: ") << completeness << std::endl;
            std::cout << (isTrueNuSlice ? "trueBestMatchPurity: " : "recoBestMatchPurity: ") << purity << std::endl;
            std::cout << (isTrueNuSlice ? "trueBestMatchTrackScore: " : "recoBestMatchTrackScore: ") << trackScore << std::endl;
            std::cout << (isTrueNuSlice ? "trueBestMatchGeneration: " : "recoBestMatchGeneration: ") << generation << std::endl;
            std::cout << (isTrueNuSlice ? "trueBestMatchSliceID: " : "recoBestMatchSliceID: ") << sliceID << std::endl;
            std::cout << (isTrueNuSlice ? "trueBestMatchNuVertexSep: " : "recoBestMatchNuVertexSep: ") << nuVertexSep << std::endl;
        }
    }

    const int nParticlesFoundOverall = nParticlesFoundInNuSlice + nParticlesFoundInOtherSlice;

    if (isTrueNuSlice)
    {
        m_nParticlesFoundInTrueNuSlice = nParticlesFoundInNuSlice;
        m_nParticlesFoundInTrueOtherSlice = nParticlesFoundInOtherSlice;
        m_nParticlesFoundTrueOverall = nParticlesFoundOverall;
    }
    else
    {
        m_nParticlesFoundInRecoNuSlice = nParticlesFoundInNuSlice;
        m_nParticlesFoundInRecoOtherSlice = nParticlesFoundInOtherSlice;
        m_nParticlesFoundRecoOverall = nParticlesFoundOverall;
    }

    if (m_debugMode)
    {
        std::cout << "Overall matching status: " << std::endl;
        std::cout << (isTrueNuSlice ? "nParticlesFoundInTrueNuSlice: " : "nParticlesFoundInRecoNuSlice: ") << nParticlesFoundInNuSlice << std::endl;
        std::cout << (isTrueNuSlice ? "nParticlesFoundInTrueOtherSlice: " : "nParticlesFoundInRecoOtherSlice: ") << nParticlesFoundInOtherSlice << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillContaminationInfo(art::Event const& evt, const bool isTrueNuSlice)
{
    const bool foundProtonMatch = (isTrueNuSlice ? m_matchFoundInTrueNuSlice[PROTON_INDEX] == 1 : m_matchFoundInRecoNuSlice[PROTON_INDEX] == 1);
    const bool foundPionMatch = (isTrueNuSlice ? m_matchFoundInTrueNuSlice[PION_INDEX] == 1 : m_matchFoundInRecoNuSlice[PION_INDEX] == 1);

    if (foundProtonMatch && foundPionMatch)
        return;

    if (!foundProtonMatch && !foundPionMatch)
        return;

    const int foundIndex = foundProtonMatch ? PROTON_INDEX : PION_INDEX;
    const int soughtIndex = foundProtonMatch ? PION_INDEX : PROTON_INDEX;
    const bool isRepassReco(isTrueNuSlice ? m_trueBestMatchRecoRepass[foundIndex] : m_recoBestMatchRecoRepass[foundIndex]);
    const MatchingMap &matchingMap(isTrueNuSlice ? (isRepassReco ? m_trueNuSliceMatchingMap_Repass : m_trueNuSliceMatchingMap) : 
        (isRepassReco ? m_recoNuSliceMatchingMap_Repass : m_recoNuSliceMatchingMap));
    const lar_pandora::PFParticleMap &pfpMap(isRepassReco ? m_pfpMap_Repass : m_pfpMap_AO);
    const art::Ptr<recob::PFParticle> &foundPFP(pfpMap.at(matchingMap.at(foundIndex).front().first));
    const int contaminantTrackID(m_trueParticleID.at(soughtIndex));
    bool foundContaminationEvidence = false;
    int nContaminantHits = 0;

    std::vector<art::Ptr<recob::Hit>> pfpHits;
    CollectHitsFromClusters(evt, foundPFP, pfpHits, isRepassReco ? "pandoraJam" : "pandoraPatRec", isRepassReco ? "pandoraJam" : "allOutcomes");

    for (const art::Ptr<recob::Hit> pfpHit : pfpHits)
    {
        if (m_trueHitMap.find(pfpHit.key()) == m_trueHitMap.end())
            continue;

        const int trackID(m_trueHitMap.at(pfpHit.key()));

        if (trackID == contaminantTrackID)
        {
            foundContaminationEvidence = true;
            ++nContaminantHits;
        }
    }

    if (foundContaminationEvidence)
    {
        if (isTrueNuSlice)
        {
            foundProtonMatch ? m_pionMergedTrueSlice = foundContaminationEvidence : m_protonMergedTrueSlice = foundContaminationEvidence;
            foundProtonMatch ? m_pionContaminantHitsTrueSlice = nContaminantHits : m_protonContaminantHitsTrueSlice = nContaminantHits;

            if (m_matchFoundInTrueOtherSlice[soughtIndex] == 1)
            {
                m_matchFoundInTrueOtherSlice[soughtIndex] = 0;
                m_nParticlesFoundInTrueOtherSlice -= 1;
                m_nParticlesFoundTrueOverall -= 1;
            }
        }
        else
        {
            foundProtonMatch ? m_pionMergedRecoSlice = foundContaminationEvidence : m_protonMergedRecoSlice = foundContaminationEvidence;
            foundProtonMatch ? m_pionContaminantHitsRecoSlice = nContaminantHits : m_protonContaminantHitsRecoSlice = nContaminantHits;

            if (m_matchFoundInRecoOtherSlice[soughtIndex] == 1)
            {
                m_matchFoundInRecoOtherSlice[soughtIndex] = 0;
                m_nParticlesFoundInRecoOtherSlice -= 1;
                m_nParticlesFoundRecoOverall -= 1;
            }
        }
    }

    if (m_debugMode)
    {
        std::cout << "Proton/Pion merging information: " << std::endl;

        if (isTrueNuSlice)
        {
            std::cout << "protonMergedTrueSlice: " << m_protonMergedTrueSlice << std::endl;
            std::cout << "protonContaminantHitsTrueSlice: " << m_protonContaminantHitsTrueSlice << std::endl; 
            std::cout << "pionMergedTrueSlice: " << m_pionMergedTrueSlice << std::endl;
            std::cout << "pionContaminantHitsTrueSlice: " << m_pionContaminantHitsTrueSlice << std::endl; 
        }
        else
        {
            std::cout << "protonMergedRecoSlice: " << m_protonMergedRecoSlice << std::endl;
            std::cout << "protonContaminantHitsRecoSlice: " << m_protonContaminantHitsRecoSlice << std::endl; 
            std::cout << "pionMergedRecoSlice: " << m_pionMergedRecoSlice << std::endl;
            std::cout << "pionContaminantHitsRecoSlice: " << m_pionContaminantHitsRecoSlice << std::endl; 
        }
    }
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

        m_foundCorrectNuSlice = (m_recoNuSliceID == m_trueNuSliceID);

        if (m_debugMode)
        {
            std::cout << "m_recoNuSliceScore: " << m_recoNuSliceScore << std::endl;
            std::cout << "m_recoNuSliceID: " << m_recoNuSliceID << std::endl;
            std::cout << "m_foundCorrectNuSlice: " << m_foundCorrectNuSlice << std::endl;
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
    m_tree->Branch("NTrueHitsU", &m_nTrueHitsU);
    m_tree->Branch("NTrueHitsV", &m_nTrueHitsV);
    m_tree->Branch("NTrueHitsW", &m_nTrueHitsW);
    m_tree->Branch("TrueNuVertexSep", &m_trueNuVertexSep);

    // Reco stuff (event level)
    m_tree->Branch("RecoNuSliceScore", &m_recoNuSliceScore);
    m_tree->Branch("RecoNuSliceID", &m_recoNuSliceID);
    m_tree->Branch("FoundCorrectNuSlice", &m_foundCorrectNuSlice);
    m_tree->Branch("FoundRecoNuVertex", &m_foundRecoNuVertex);
    m_tree->Branch("RecoNuVertexX", &m_recoNuVertexX);
    m_tree->Branch("RecoNuVertexY", &m_recoNuVertexY);
    m_tree->Branch("RecoNuVertexZ", &m_recoNuVertexZ);

    // Reco stuff - match found?
    m_tree->Branch("MatchFoundInTrueNuSlice", &m_matchFoundInTrueNuSlice);
    m_tree->Branch("MatchFoundInTrueOtherSlice", &m_matchFoundInTrueOtherSlice);
    m_tree->Branch("MatchFoundInRecoNuSlice", &m_matchFoundInRecoNuSlice);
    m_tree->Branch("MatchFoundInRecoOtherSlice", &m_matchFoundInRecoOtherSlice);

    m_tree->Branch("NParticlesFoundInTrueNuSlice", &m_nParticlesFoundInTrueNuSlice);
    m_tree->Branch("NParticlesFoundInTrueOtherSlice", &m_nParticlesFoundInTrueOtherSlice);
    m_tree->Branch("NParticlesFoundTrueOverall", &m_nParticlesFoundTrueOverall);

    m_tree->Branch("NParticlesFoundInRecoNuSlice", &m_nParticlesFoundInRecoNuSlice);
    m_tree->Branch("NParticlesFoundInRecoOtherSlice", &m_nParticlesFoundInRecoOtherSlice);
    m_tree->Branch("NParticlesFoundRecoOverall", &m_nParticlesFoundRecoOverall);

    m_tree->Branch("ProtonMergedTrueSlice", &m_protonMergedTrueSlice);
    m_tree->Branch("ProtonContaminantHitsTrueSlice", &m_protonContaminantHitsTrueSlice);
    m_tree->Branch("PionMergedTrueSlice", &m_pionMergedTrueSlice);
    m_tree->Branch("PionContaminantHitsTrueSlice", &m_pionContaminantHitsTrueSlice);

    m_tree->Branch("ProtonMergedRecoSlice", &m_protonMergedRecoSlice);
    m_tree->Branch("ProtonContaminantHitsRecoSlice", &m_protonContaminantHitsRecoSlice);
    m_tree->Branch("PionMergedRecoSlice", &m_pionMergedRecoSlice);
    m_tree->Branch("PionContaminantHitsRecoSlice", &m_pionContaminantHitsRecoSlice);

    // Reco stuff - best match stuff  
    m_tree->Branch("TrueBestMatchRecoRepass", &m_trueBestMatchRecoRepass);
    m_tree->Branch("TrueBestMatchCompleteness", &m_trueBestMatchCompleteness);
    m_tree->Branch("TrueBestMatchPurity", &m_trueBestMatchPurity);
    m_tree->Branch("TrueBestMatchTrackScore", &m_trueBestMatchTrackScore);
    m_tree->Branch("TrueBestMatchGeneration", &m_trueBestMatchGeneration);
    m_tree->Branch("TrueBestMatchSliceID", &m_trueBestMatchSliceID);
    m_tree->Branch("TrueBestMatchRecoVertexX", &m_trueBestMatchRecoVertexX);
    m_tree->Branch("TrueBestMatchRecoVertexY", &m_trueBestMatchRecoVertexY);
    m_tree->Branch("TrueBestMatchRecoVertexZ", &m_trueBestMatchRecoVertexZ);
    m_tree->Branch("TrueBestMatchNuVertexSep", &m_trueBestMatchNuVertexSep);

    m_tree->Branch("RecoBestMatchRecoRepass", &m_recoBestMatchRecoRepass);
    m_tree->Branch("RecoBestMatchCompleteness", &m_recoBestMatchCompleteness);
    m_tree->Branch("RecoBestMatchPurity", &m_recoBestMatchPurity);
    m_tree->Branch("RecoBestMatchTrackScore", &m_recoBestMatchTrackScore);
    m_tree->Branch("RecoBestMatchGeneration", &m_recoBestMatchGeneration);
    m_tree->Branch("RecoBestMatchSliceID", &m_recoBestMatchSliceID);
    m_tree->Branch("RecoBestMatchRecoVertexX", &m_recoBestMatchRecoVertexX);
    m_tree->Branch("RecoBestMatchRecoVertexY", &m_recoBestMatchRecoVertexY);
    m_tree->Branch("RecoBestMatchRecoVertexZ", &m_recoBestMatchRecoVertexZ);
    m_tree->Branch("RecoBestMatchNuVertexSep", &m_recoBestMatchNuVertexSep);
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::Reset()
{
    m_trueHitMap.clear();
    m_pfpMap.clear();
    m_pfpMap_AO.clear();
    m_pfpMap_Repass.clear();
    m_sliceMap_AO.clear();
    m_mcParticleMap.clear();
    m_trueNuSliceMatchingMap.clear();
    m_trueNuSliceMatchingMap_Repass.clear();
    m_trueOtherSliceMatchingMap.clear();
    m_recoNuSliceMatchingMap.clear();
    m_recoNuSliceMatchingMap_Repass.clear();
    m_recoOtherSliceMatchingMap.clear();

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
    m_nTrueHitsU.clear();
    m_nTrueHitsU.resize(4, DEFAULT_INT);
    m_nTrueHitsV.clear();
    m_nTrueHitsV.resize(4, DEFAULT_INT);
    m_nTrueHitsW.clear();
    m_nTrueHitsW.resize(4, DEFAULT_INT);
    m_trueNuVertexSep.clear();
    m_trueNuVertexSep.resize(4, DEFAULT_DOUBLE);

    // Reco stuff (event level)
    m_recoNuSliceScore = DEFAULT_DOUBLE;
    m_recoNuSliceID = DEFAULT_DOUBLE;
    m_foundCorrectNuSlice = DEFAULT_BOOL;
    m_foundRecoNuVertex = DEFAULT_BOOL;
    m_recoNuVertexX = DEFAULT_DOUBLE;
    m_recoNuVertexY = DEFAULT_DOUBLE;
    m_recoNuVertexZ = DEFAULT_DOUBLE;

    // Reco stuff - match found?
    m_matchFoundInTrueNuSlice.clear();
    m_matchFoundInTrueNuSlice.resize(4, DEFAULT_INT);
    m_matchFoundInTrueOtherSlice.clear();
    m_matchFoundInTrueOtherSlice.resize(4, DEFAULT_INT);
    m_matchFoundInRecoNuSlice.clear();
    m_matchFoundInRecoNuSlice.resize(4, DEFAULT_INT);
    m_matchFoundInRecoOtherSlice.clear();
    m_matchFoundInRecoOtherSlice.resize(4, DEFAULT_INT);

    m_nParticlesFoundInTrueNuSlice = DEFAULT_INT;
    m_nParticlesFoundInTrueOtherSlice = DEFAULT_INT;
    m_nParticlesFoundTrueOverall = DEFAULT_INT;
    m_nParticlesFoundInRecoNuSlice = DEFAULT_INT;
    m_nParticlesFoundInRecoOtherSlice = DEFAULT_INT;
    m_nParticlesFoundRecoOverall = DEFAULT_INT;
    m_pionMergedTrueSlice = DEFAULT_BOOL;
    m_pionContaminantHitsTrueSlice = DEFAULT_INT;
    m_protonMergedTrueSlice = DEFAULT_BOOL;
    m_protonContaminantHitsTrueSlice = DEFAULT_INT;
    m_pionMergedRecoSlice = DEFAULT_BOOL;
    m_pionContaminantHitsRecoSlice = DEFAULT_INT;
    m_protonMergedRecoSlice = DEFAULT_BOOL;
    m_protonContaminantHitsRecoSlice = DEFAULT_INT;

    // Reco stuff - best match stuff   
    m_trueBestMatchRecoRepass.clear();
    m_trueBestMatchRecoRepass.resize(4, DEFAULT_BOOL);
    m_trueBestMatchCompleteness.clear();
    m_trueBestMatchCompleteness.resize(4, DEFAULT_DOUBLE);
    m_trueBestMatchPurity.clear();
    m_trueBestMatchPurity.resize(4, DEFAULT_DOUBLE);
    m_trueBestMatchTrackScore.clear();
    m_trueBestMatchTrackScore.resize(4, DEFAULT_DOUBLE);
    m_trueBestMatchGeneration.clear();
    m_trueBestMatchGeneration.resize(4, DEFAULT_INT);
    m_trueBestMatchSliceID.clear();
    m_trueBestMatchSliceID.resize(4, DEFAULT_DOUBLE);
    m_trueBestMatchRecoVertexX.clear();
    m_trueBestMatchRecoVertexX.resize(4, DEFAULT_DOUBLE);
    m_trueBestMatchRecoVertexY.clear();
    m_trueBestMatchRecoVertexY.resize(4, DEFAULT_DOUBLE);
    m_trueBestMatchRecoVertexZ.clear();
    m_trueBestMatchRecoVertexZ.resize(4, DEFAULT_DOUBLE);
    m_trueBestMatchNuVertexSep.clear();
    m_trueBestMatchNuVertexSep.resize(4, DEFAULT_DOUBLE);

    m_recoBestMatchRecoRepass.clear();
    m_recoBestMatchRecoRepass.resize(4, DEFAULT_BOOL);
    m_recoBestMatchCompleteness.clear();
    m_recoBestMatchCompleteness.resize(4, DEFAULT_DOUBLE);
    m_recoBestMatchPurity.clear();
    m_recoBestMatchPurity.resize(4, DEFAULT_DOUBLE);
    m_recoBestMatchTrackScore.clear();
    m_recoBestMatchTrackScore.resize(4, DEFAULT_DOUBLE);
    m_recoBestMatchGeneration.clear();
    m_recoBestMatchGeneration.resize(4, DEFAULT_INT);
    m_recoBestMatchSliceID.clear();
    m_recoBestMatchSliceID.resize(4, DEFAULT_DOUBLE);
    m_recoBestMatchRecoVertexX.clear();
    m_recoBestMatchRecoVertexX.resize(4, DEFAULT_DOUBLE);
    m_recoBestMatchRecoVertexY.clear();
    m_recoBestMatchRecoVertexY.resize(4, DEFAULT_DOUBLE);
    m_recoBestMatchRecoVertexZ.clear();
    m_recoBestMatchRecoVertexZ.resize(4, DEFAULT_DOUBLE);
    m_recoBestMatchNuVertexSep.clear();
    m_recoBestMatchNuVertexSep.resize(4, DEFAULT_DOUBLE);
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(hyperon::SigmaRecoAnalyser)


/*
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


        if (particleIndex == PROTON_INDEX)
        {
            if ((maxOwnerID != m_trueParticleID[PROTON_INDEX]) && (maxOwnerID == m_trueParticleID[PION_INDEX]))
            {
                m_isPionContaminatedByProton = true;
                foundInTrueNuSlice = true;
            }
        }

        if (particleIndex == PION_INDEX)
        {
            if ((maxOwnerID != m_trueParticleID[PION_INDEX]) && (maxOwnerID == m_trueParticleID[PROTON_INDEX]))
            {
                m_isProtonContaminatedByPion = true;
                foundInTrueNuSlice = true;
            }
        }
*/
