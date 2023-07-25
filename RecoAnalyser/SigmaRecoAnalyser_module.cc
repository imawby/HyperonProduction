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

// Services
#include "art_root_io/TFileService.h"

// Hyperon
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"

// larpandora
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

const int DEFAULT_INT = -999;
const double DEFAULT_DOUBLE = -999.0;
const std::string DEFAULT_STRING = "DEFAULT"; 
const int MUON_INDEX = 0;
const int PROTON_INDEX = 1;
const int PION_INDEX = 2;
const int GAMMA_INDEX = 3;


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
    bool IdentifySignalParticle(const std::vector<SimParticle> &simParticles,
        const int pdg, int &signalParticleID);

private:

    // Product labels
    fhicl::ParameterSet m_G4Labels;
    fhicl::ParameterSet m_generatorLabels;
    std::string m_MCParticleLabel;

    // Debug Info?
    bool m_debugMode;

    // ID info 
    int m_mcTruthIndex;
    int m_trueMuonID;
    int m_trueProtonID;
    int m_truePionID;
    int m_trueGammaID;
    int m_bestMatchedMuonID;
    int m_bestMatchedProtonID;
    int m_bestMatchedPionID;
    int m_bestMatchedGammaID;
    std::vector<int> m_muonMatchIDs;
    std::vector<int> m_protonMatchIDs;
    std::vector<int> m_pionMatchIDs;
    std::vector<int> m_gammaMatchIDs;

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
    int m_gammaFoundInSlice;
    int m_gammaFoundInOtherSlice;
    double m_sliceNuScore;

    // Reco stuff (particle level)
    std::vector<int> m_nMatches;
    std::vector<int> m_bestMatchFound;
    std::vector<double> m_bestMatchCompleteness;
    std::vector<double> m_bestMatchPurity;
    std::vector<double> m_bestMatchTrackScore;
    std::vector<int> m_bestMatchGeneration;

    // Linking MCParticle -> TrackID
    lar_pandora::MCParticleMap m_mcParticleMap;
};

///////////////////////////////////////////////////////////////////////////////////////////

hyperon::SigmaRecoAnalyser::SigmaRecoAnalyser(fhicl::ParameterSet const& pset)
    : EDAnalyzer{pset},
    m_G4Labels(pset.get<fhicl::ParameterSet>("Geant4"))
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

    // Identify truth IDs
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

    if (!evt.getByLabel(m_MCParticleLabel, mcParticleHandle))
        throw cet::exception("SigmaRecoAnalyser") << "No MCParticle Data Products Found!" << std::endl;

    art::fill_ptr_vector(mcParticleVector, mcParticleHandle);
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, m_mcParticleMap);

    if (!IdentifySignalParticles(G4Truth))
        return;


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
    if (!IdentifySignalParticle(G4Truth.Lepton, -13, m_trueMuonID))
        return false;

    if (!IdentifySignalParticle(G4Truth.Decay, 2212, m_trueProtonID))
        return false;

    if (!IdentifySignalParticle(G4Truth.Lepton, -211, m_truePionID))
        return false;

    if (!IdentifySignalParticle(G4Truth.SigmaZeroDecayPhoton, 22, m_trueGammaID))
        return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::SigmaRecoAnalyser::IdentifySignalParticle(const std::vector<SimParticle> &simParticles,
    const int pdg, int &signalParticleID)
{
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
        signalParticleID = simParticle.ArtID;
        break;
    }

    return found;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::FillTruthEventInfo(const GeneratorTruth & genTruth, const G4Truth &G4Truth)
{
    const art::Ptr<simb::MCParticle> &protonMCParticle(m_mcParticleMap.at(m_trueProtonID));
    const art::Ptr<simb::MCParticle> &pionMCParticle(m_mcParticleMap.at(m_truePionID));
    const art::Ptr<simb::MCParticle> &gammaMCParticle(m_mcParticleMap.at(m_trueGammaID));

    const TVector3 &protonMom(protonMCParticle->Momentum().Vect());
    const TVector3 &pionMom(pionMCParticle->Momentum().Vect());
    const TVector3 &lambdaMom(protonMom + pionMom);
    const TVector3 &gammaMom(gammaMCParticle->Momentum().Vect());

    m_trueProtonPiOpeningAngle = protonMom.Angle(pionMom);
    m_trueGammaLambdaOpeningAngle = gammaMom.Angle(lambdaMom);

    // NEED NU VERTEX AGGGHHH

m_trueGammaLambdaVertexSep 

    m_tree->Branch("TrueProtonPiOpeningAngle", &m_trueProtonPiOpeningAngle);
    m_tree->Branch("TrueGammaLambdaOpeningAngle", &m_trueGammaLambdaOpeningAngle);
    m_tree->Branch("TrueGammaLambdaVertexSep", &m_trueGammaLambdaVertexSep);

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
    m_tree->Branch("TrueMuonID", &m_trueMuonID);
    m_tree->Branch("TrueProtonID", &m_trueProtonID);
    m_tree->Branch("TruePionID", &m_truePionID);
    m_tree->Branch("TrueGammaID", &m_trueGammaID);
    m_tree->Branch("BestMatchedMuonID", &m_bestMatchedMuonID);
    m_tree->Branch("BestMatchedProtonID", &m_bestMatchedProtonID);
    m_tree->Branch("BestMatchedPionID", &m_bestMatchedPionID);
    m_tree->Branch("BestMatchedGammaID", &m_bestMatchedGammaID);
    m_tree->Branch("MuonMatchIDs", &m_muonMatchIDs);
    m_tree->Branch("ProtonMatchIDs", &m_protonMatchIDs);
    m_tree->Branch("PionMatchIDs", &m_pionMatchIDs);
    m_tree->Branch("GammaMatchIDs", &m_gammaMatchIDs);

    // Truth stuff (particle level)
    m_tree->Branch("TruePDG", &m_truePDG);
    m_tree->Branch("NTrueHits", &m_nTrueHits);
    m_tree->Branch("TrueNuVertexSep", &m_trueNuVertexSep);

    // Reco stuff (event level)
    m_tree->Branch("GammaFoundInSlice", &m_gammaFoundInSlice);
    m_tree->Branch("GammaFoundInOtherSlice", &m_gammaFoundInOtherSlice);
    m_tree->Branch("SliceNuScore", &m_sliceNuScore);

    // Reco stuff (particle level)
    m_tree->Branch("NMatches", &m_nMatches);
    m_tree->Branch("BestMatchFound", &m_bestMatchFound);
    m_tree->Branch("BestMatchCompleteness", &m_bestMatchCompleteness);
    m_tree->Branch("BestMatchPurity", &m_bestMatchPurity);
    m_tree->Branch("BestMatchTrackScore", &m_bestMatchTrackScore);
    m_tree->Branch("BestMatchGeneration", &m_bestMatchGeneration);
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::Reset()
{
    m_mcParticleMap.clear();

    m_mcTruthIndex = DEFAULT_INT;
    m_trueMuonID = DEFAULT_INT;
    m_trueProtonID = DEFAULT_INT;
    m_truePionID = DEFAULT_INT;
    m_trueGammaID = DEFAULT_INT;
    m_bestMatchedMuonID = DEFAULT_INT;
    m_bestMatchedProtonID = DEFAULT_INT;
    m_bestMatchedPionID = DEFAULT_INT;
    m_bestMatchedGammaID = DEFAULT_INT;
    m_muonMatchIDs.clear();
    m_protonMatchIDs.clear();
    m_pionMatchIDs.clear();
    m_gammaMatchIDs.clear();

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
    m_truePDG.resize(3, DEFAULT_INT);
    m_nTrueHits.clear();
    m_nTrueHits.resize(3, DEFAULT_INT);
    m_trueNuVertexSep.clear();
    m_trueNuVertexSep.resize(3, DEFAULT_DOUBLE);

    // Reco stuff (event level)
    m_gammaFoundInSlice = DEFAULT_INT;
    m_gammaFoundInOtherSlice = DEFAULT_INT;
    m_sliceNuScore = DEFAULT_DOUBLE;

    // Reco stuff (particle level)
    m_nMatches.clear();
    m_nMatches.resize(3, DEFAULT_INT);
    m_bestMatchFound.clear();
    m_bestMatchFound.resize(3, DEFAULT_INT);
    m_bestMatchCompleteness.clear();
    m_bestMatchCompleteness.resize(3, DEFAULT_DOUBLE);
    m_bestMatchPurity.clear();
    m_bestMatchPurity.resize(3, DEFAULT_DOUBLE);
    m_bestMatchTrackScore.clear();
    m_bestMatchTrackScore.resize(3, DEFAULT_DOUBLE);
    m_bestMatchGeneration.clear();
    m_bestMatchGeneration.resize(3, DEFAULT_INT);
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::SigmaRecoAnalyser::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(hyperon::SigmaRecoAnalyser)
