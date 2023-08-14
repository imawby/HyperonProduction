////////////////////////////////////////////////////////////////////////
// Class:       LambdaVertexProducer
// Plugin Type: producer (art v3_01_02)
// File:        LambdaVertexProducer_module.cc
//
// Generated at Wed Aug 21 17:07:38 2019 by Giuseppe Cerati using cetskelgen
// Modified by C Thorpe Sept 2022.
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "ubana/HyperonProduction/Modules/SubModules/SubModuleReco.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <memory>

const int MUON_INDEX = 0;
const int PROTON_INDEX = 1;
const int PION_INDEX = 2;
const int GAMMA_INDEX = 3;

namespace hyperon {
   class LambdaVertexProducer;
}

class hyperon::LambdaVertexProducer : public art::EDProducer {
   public:
      explicit LambdaVertexProducer(fhicl::ParameterSet const& p);
      LambdaVertexProducer(LambdaVertexProducer const&) = delete;
      LambdaVertexProducer(LambdaVertexProducer&&) = delete;
      LambdaVertexProducer& operator=(LambdaVertexProducer const&) = delete;
      LambdaVertexProducer& operator=(LambdaVertexProducer&&) = delete;

      void produce(art::Event& e) override;
      bool IdentifySignalParticles(const G4Truth &G4Truth, const int mcTruthIndex, const lar_pandora::MCParticleMap &mcParticleMap, 
          int &muonTrackID, int &protonTrackID, int &pionTrackID, int &gammaTrackID);
      bool IdentifySignalParticle(const G4Truth &G4Truth, const int mcTruthIndex, const lar_pandora::MCParticleMap &mcParticleMap, 
          const int particleTypeIndex, int &particleTrackID);
      void FillMCParticleHitInfo(art::Event const& evt, const lar_pandora::MCParticleMap &mcParticleMap,
          std::map<int, int> &trueHitMap);
      bool IsEM(const art::Ptr<simb::MCParticle> &mcParticle);
      int GetLeadEMTrackID(const lar_pandora::MCParticleMap &mcParticleMap, const art::Ptr<simb::MCParticle> &mcParticle);
      void FillSliceMap(art::Event& e, std::map<int, int> &sliceMap);
      int GetTrueNuSliceID(art::Event const& evt, const std::map<int, int> &trueHitMap, 
          const int muonTrackID, const int protonTrackID, const int pionTrackID, const int gammaTrackID);
      void CollectChildren(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, lar_pandora::PFParticleMap &pfparticleMap,
          std::vector<art::Ptr<recob::PFParticle>> &collectedParticles, std::vector<art::Ptr<recob::Hit>> &collectedHits);
      void CollectHits(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, std::vector<art::Ptr<recob::Hit>> &hits);
      bool GetMatchedTrackID(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, 
          const std::map<int, int> &trueHitMap, int &matchedTrackID);
      bool GetDecayVertex(art::Event& e, std::vector<art::Ptr<recob::PFParticle>> &protonPFParticles, 
          std::vector<art::Ptr<recob::PFParticle>> &pionPFParticles, double &decayVertexX, double &decayVertexY, double &decayVertexZ);
   private:
      std::string f_MCParticleLabel;
      std::string f_PFParticleLabel;
      std::string f_HitLabel;
      std::string f_ClusterLabel;
      std::string f_TruthMatchingLabel;
      std::string f_PFPInstanceName;

      fhicl::ParameterSet f_Generator;
      fhicl::ParameterSet f_G4;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////

hyperon::LambdaVertexProducer::LambdaVertexProducer(fhicl::ParameterSet const& p)
   : EDProducer{p},
   f_MCParticleLabel(p.get<std::string>("MCParticleLabel", "largeant")),
   f_PFParticleLabel(p.get<std::string>("PFParticleLabel")),
   f_HitLabel(p.get<std::string>("HitLabel", "gaushit")),
   f_ClusterLabel(p.get<std::string>("ClusterLabel")),
   f_TruthMatchingLabel(p.get<std::string>("TruthMatchingLabel", "gaushitTruthMatch")),
   f_PFPInstanceName(p.get<std::string>("PFPInstanceName", "")),
   f_Generator(p.get<fhicl::ParameterSet>("Generator")),
   f_G4(p.get<fhicl::ParameterSet>("Geant4"))
{
   produces<std::vector<recob::Vertex>>();
   produces<std::vector<recob::Hit>>();
   produces<std::vector<int>>();
   produces<art::Assns<recob::Slice,recob::Vertex,void>>();
   produces<art::Assns<recob::Slice,recob::Hit,void>>();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::LambdaVertexProducer::produce(art::Event& e)
{
   std::unique_ptr< std::vector<recob::Vertex> > outputVertex(new std::vector<recob::Vertex>);
   std::unique_ptr< std::vector<recob::Hit> > outputHitVector(new std::vector<recob::Hit>);
   std::unique_ptr< std::vector<int> > outputPFParticleVector(new std::vector<int>);
   std::unique_ptr< art::Assns<recob::Slice, recob::Vertex> > outputSliceVertexAssoc(new art::Assns<recob::Slice, recob::Vertex>);
   std::unique_ptr< art::Assns<recob::Slice, recob::Hit> > outputSliceHitAssoc(new art::Assns<recob::Slice, recob::Hit>);

   ///////////////////////////////////////////
   // First check that we have a signal event
   ///////////////////////////////////////////

   SubModuleGeneratorTruth* Generator_SM = new SubModuleGeneratorTruth(e, f_Generator);
   GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();

   SubModuleG4Truth* G4_SM = new SubModuleG4Truth(e, f_G4);
   G4Truth G4T = G4_SM->GetG4Info();

   bool isSignalSigmaZero = false;
   int mcTruthIndex = -1;

   for (int i = 0; i < GenT.NMCTruths; i++)
   {
       if (GenT.Mode.at(i) != "HYP")
           continue;

       if (!G4T.InActiveTPC.at(i))
           continue;

       if (GenT.Neutrino.at(i).PDG != -14)
           continue;

       if (!G4T.IsSigmaZeroCharged.at(i))
           continue;

       if (G4T.IsAssociatedHyperon.at(i))
           continue;

       isSignalSigmaZero = true;
       mcTruthIndex = i;
       break;
   }
   
   std::cout << "isSignalSigmaZero: " << (isSignalSigmaZero ? "yes" : "no") << std::endl;

   // Put empty products into file if not signal   
   if (!isSignalSigmaZero)
   {
       e.put(std::move(outputVertex));
       e.put(std::move(outputHitVector));
       e.put(std::move(outputPFParticleVector));
       e.put(std::move(outputSliceVertexAssoc));
       e.put(std::move(outputSliceHitAssoc));
       return;
   }
   
   ///////////////////////////////////////////
   // Now check whether we need to modify the reco
   ///////////////////////////////////////////

   // For truth matching
   art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
   std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

   if (!e.getByLabel(f_MCParticleLabel, mcParticleHandle)) 
       throw cet::exception("LambdaVertexProducer") << "No MCParticle Data Products Found!" << std::endl;

   art::fill_ptr_vector(mcParticleVector, mcParticleHandle);

   lar_pandora::MCParticleMap mcParticleMap;
   lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, mcParticleMap);

   int muonTrackID(-1), protonTrackID(-1), pionTrackID(-1), gammaTrackID(-1);

   if (!IdentifySignalParticles(G4T, mcTruthIndex, mcParticleMap, muonTrackID, protonTrackID, pionTrackID, gammaTrackID))
   {
       e.put(std::move(outputVertex));
       e.put(std::move(outputHitVector));
       e.put(std::move(outputPFParticleVector));
       e.put(std::move(outputSliceVertexAssoc));
       e.put(std::move(outputSliceHitAssoc));
       return;
   }

   std::cout << "muonTrackID: " << muonTrackID << std::endl;
   std::cout << "protonTrackID: " << protonTrackID << std::endl;
   std::cout << "pionTrackID: " << pionTrackID << std::endl;
   std::cout << "gammaTrackID: " << gammaTrackID << std::endl;


   std::map<int, int> trueHitMap;
   FillMCParticleHitInfo(e, mcParticleMap, trueHitMap);

   int trueNuSliceID(GetTrueNuSliceID(e, trueHitMap, muonTrackID, protonTrackID, pionTrackID, gammaTrackID));

   if (trueNuSliceID == -1)
   {
       e.put(std::move(outputVertex));
       e.put(std::move(outputHitVector));
       e.put(std::move(outputPFParticleVector));
       e.put(std::move(outputSliceVertexAssoc));
       e.put(std::move(outputSliceHitAssoc));
       return;
   }

   // Get event PFParticles
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;
   std::vector<art::Ptr<recob::PFParticle>> pfparticleVector;

   art::InputTag pfpInputTag(f_PFParticleLabel, f_PFPInstanceName);
   if (!e.getByLabel(pfpInputTag, pfparticleHandle))
       throw cet::exception("LambdaVertexProducer") << "No PFParticle Data Products Found!" << std::endl;

   art::fill_ptr_vector(pfparticleVector, pfparticleHandle);

   // Now fill the Pandora PFParticle map
   lar_pandora::PFParticleMap pfparticleMap;
   lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleVector, pfparticleMap);

   // Now fill the slice map
   std::map<int, int> sliceMap;
   FillSliceMap(e, sliceMap);

   // Now search for the PFParticle matches (look only in the true neutrino slice!)
   std::vector<art::Ptr<recob::PFParticle>> pionPFParticles;
   std::vector<art::Ptr<recob::PFParticle>> protonPFParticles;

   // Get true nu slice pfps...
   art::InputTag sliceInputTag(f_PFParticleLabel, f_PFPInstanceName);
   art::Handle<std::vector<recob::Slice>> sliceHandle;

   if (!e.getByLabel(sliceInputTag, sliceHandle))
       throw cet::exception("SigmaRecoAnalyser") << "No Slice Data Products Found!" << std::endl;

   art::FindManyP<recob::PFParticle> pfpAssoc = art::FindManyP<recob::PFParticle>(sliceHandle, e, sliceInputTag);
   std::vector<art::Ptr<recob::PFParticle>> nuSlicePFPs = pfpAssoc.at(sliceMap[trueNuSliceID]);

   for (const art::Ptr<recob::PFParticle> pfparticle : nuSlicePFPs)
   {
       // Make sure that particle is in the neutrino hierarchy
       if (lar_pandora::LArPandoraHelper::GetParentNeutrino(pfparticleMap, pfparticle) == 0)
           continue;

       int matchedTrackID(-1);
       if (!GetMatchedTrackID(e, pfparticle, trueHitMap, matchedTrackID))
           continue;

       std::cout << "matchedTrackID: " << matchedTrackID << std::endl;
       std::cout << "matched PDG: " << mcParticleMap.at(matchedTrackID)->PdgCode() << std::endl;

       if (matchedTrackID == protonTrackID)
           protonPFParticles.push_back(pfparticle);
       else if(matchedTrackID == pionTrackID)
           pionPFParticles.push_back(pfparticle);
   }

   const bool foundPion = (pionPFParticles.size() != 0);
   const bool foundProton = (protonPFParticles.size() != 0);
   const bool modifyReco = foundPion || foundProton;

   if (!modifyReco)
   {
       e.put(std::move(outputVertex));
       e.put(std::move(outputHitVector));
       e.put(std::move(outputPFParticleVector));
       e.put(std::move(outputSliceVertexAssoc));
       e.put(std::move(outputSliceHitAssoc));

       return;
   }

   // Now we collect the PFParticles in the hierarchy and their hits
   std::vector<art::Ptr<recob::PFParticle>> lambdaHierarachyPFParticles;
   std::vector<art::Ptr<recob::Hit>> lambdaHierarachyHits;

   for (art::Ptr<recob::PFParticle> pfparticle : protonPFParticles)
       this->CollectChildren(e, pfparticle, pfparticleMap, lambdaHierarachyPFParticles, lambdaHierarachyHits);

   for (art::Ptr<recob::PFParticle> pfparticle : pionPFParticles)
       this->CollectChildren(e, pfparticle, pfparticleMap, lambdaHierarachyPFParticles, lambdaHierarachyHits);

    double decayVertexX(-999.0), decayVertexY(-999.0), decayVertexZ(-999.0);
    if (!this->GetDecayVertex(e, protonPFParticles, pionPFParticles, decayVertexX, decayVertexY, decayVertexZ))
    {
        std::cout << "could not find decay vertex..." << std::endl;

        e.put(std::move(outputVertex));
        e.put(std::move(outputHitVector));
        e.put(std::move(outputPFParticleVector));
        e.put(std::move(outputSliceVertexAssoc));
        e.put(std::move(outputSliceHitAssoc));

        return;
    }

    double xyz[3] = {decayVertexX, decayVertexY, decayVertexZ};
    recob::Vertex decayVertex(xyz, -1);

    outputVertex->push_back(decayVertex);

    // Make an association to the slice
   art::FindManyP<recob::Slice> pfparticleSliceAssoc = art::FindManyP<recob::Slice>(pfparticleHandle, e, pfpInputTag);
   std::vector<art::Ptr<recob::Slice>> slices = pfparticleSliceAssoc.at(foundProton ? protonPFParticles.front().key() : pionPFParticles.front().key());
   util::CreateAssn(*this, e, *outputVertex, slices[0], *outputSliceVertexAssoc);

   for (const art::Ptr<recob::Hit> hit : lambdaHierarachyHits)
   {
       outputHitVector->push_back(*hit.get());
       outputSliceHitAssoc->addSingle(slices[0], hit);
   }

   for (const art::Ptr<recob::PFParticle> pfparticle : lambdaHierarachyPFParticles)
       outputPFParticleVector->push_back(pfparticle.key());

   e.put(std::move(outputVertex));
   e.put(std::move(outputHitVector));
   e.put(std::move(outputPFParticleVector));
   e.put(std::move(outputSliceVertexAssoc));
   e.put(std::move(outputSliceHitAssoc));

   return;

}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::LambdaVertexProducer::IdentifySignalParticles(const G4Truth &G4Truth, const int mcTruthIndex, 
    const lar_pandora::MCParticleMap &mcParticleMap, int &muonTrackID, int &protonTrackID, int &pionTrackID, 
    int &gammaTrackID)
{
    if (!IdentifySignalParticle(G4Truth, mcTruthIndex, mcParticleMap, MUON_INDEX, muonTrackID))
    {
        std::cout << "Couldn't find MC muon :(" << std::endl;
        return false;
    }

    if (!IdentifySignalParticle(G4Truth, mcTruthIndex, mcParticleMap, PROTON_INDEX, protonTrackID))
    {
        std::cout << "Couldn't find MC proton :(" << std::endl;
        return false;
    }

    if (!IdentifySignalParticle(G4Truth, mcTruthIndex, mcParticleMap, PION_INDEX, pionTrackID))
    {
        std::cout << "Couldn't find MC pion :(" << std::endl;
        return false;
    }

    if (!IdentifySignalParticle(G4Truth, mcTruthIndex, mcParticleMap, GAMMA_INDEX, gammaTrackID))
    {
        std::cout << "Couldn't find MC gamma :(" << std::endl;
        return false;
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::LambdaVertexProducer::IdentifySignalParticle(const G4Truth &G4Truth, const int mcTruthIndex, 
    const lar_pandora::MCParticleMap &mcParticleMap, const int particleTypeIndex, int &particleTrackID)
{
    const std::vector<SimParticle> &simParticles(particleTypeIndex == MUON_INDEX ? G4Truth.Lepton : 
        particleTypeIndex == PROTON_INDEX ? G4Truth.Decay : particleTypeIndex == PION_INDEX ? G4Truth.Decay : G4Truth.SigmaZeroDecayPhoton);

    const int pdg(particleTypeIndex == MUON_INDEX ? -13 : particleTypeIndex == PROTON_INDEX ? 2212 : particleTypeIndex == PION_INDEX ? -211 : 22);

    bool found = false;

    for (const SimParticle &simParticle : simParticles)
    {
        if (simParticle.MCTruthIndex != mcTruthIndex)
            continue;

        if (simParticle.PDG != pdg)
            continue;

        if (mcParticleMap.find(simParticle.ArtID) == mcParticleMap.end())
            continue;

        found = true;
        particleTrackID = simParticle.ArtID;
        break;
    }

    return found;
}

///////////////////////////////////////////////////////////////////////////////////////////

void hyperon::LambdaVertexProducer::FillMCParticleHitInfo(art::Event const& evt, const lar_pandora::MCParticleMap &mcParticleMap, 
    std::map<int, int> &trueHitMap)
{
    // Get event hits
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVector;

    if (!evt.getByLabel(f_HitLabel, hitHandle))
        throw cet::exception("LambdaVertexProducer") << "No Hit Data Products Found!" << std::endl;

    art::fill_ptr_vector(hitVector, hitHandle);

    // Get backtracker info
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> assocMCPart = art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitHandle, evt, f_TruthMatchingLabel);

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

            const int trackID = IsEM(matchedMCParticle) ? GetLeadEMTrackID(mcParticleMap, matchedMCParticle) : matchedMCParticle->TrackId();

            trueHitMap[hit.key()] = trackID;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::LambdaVertexProducer::IsEM(const art::Ptr<simb::MCParticle> &mcParticle)
{
    return ((std::abs(mcParticle->PdgCode()) == 11) || (mcParticle->PdgCode() == 22));
}

///////////////////////////////////////////////////////////////////////////////////////////
// If its an EM particle, we have to move up the EM hierarchy
int hyperon::LambdaVertexProducer::GetLeadEMTrackID(const lar_pandora::MCParticleMap &mcParticleMap, const art::Ptr<simb::MCParticle> &mcParticle)
{
    int trackID = mcParticle->TrackId();
    art::Ptr<simb::MCParticle> motherMCParticle = mcParticle;

    do
    {
        trackID = motherMCParticle->TrackId();
        const int motherID = motherMCParticle->Mother();

        if (mcParticleMap.find(motherID) == mcParticleMap.end())
            break;

        motherMCParticle = mcParticleMap.at(motherID);
    }
    while (IsEM(motherMCParticle));

    return trackID;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::LambdaVertexProducer::FillSliceMap(art::Event& e, std::map<int, int> &sliceMap)
{
   art::InputTag sliceInputTag(f_PFParticleLabel, f_PFPInstanceName);
   art::Handle<std::vector<recob::Slice>> sliceHandle;
   std::vector<art::Ptr<recob::Slice>> sliceVector;

   if (!e.getByLabel(sliceInputTag, sliceHandle))
       throw cet::exception("LambdaVertexProducer") << "No Slice Data Products Found!" << std::endl;

   art::fill_ptr_vector(sliceVector, sliceHandle);

   for (art::Ptr<recob::Slice> &slice : sliceVector)
        sliceMap[slice->ID()] = slice.key();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

int hyperon::LambdaVertexProducer::GetTrueNuSliceID(art::Event const& evt, const std::map<int, int> &trueHitMap, 
    const int muonTrackID, const int protonTrackID, const int pionTrackID, const int gammaTrackID)
{
    // Get slice information
    art::InputTag sliceInputTag(f_PFParticleLabel, f_PFPInstanceName);
    art::Handle<std::vector<recob::Slice>> sliceHandle;
    std::vector<art::Ptr<recob::Slice>> sliceVector;

    if (!evt.getByLabel(sliceInputTag, sliceHandle))
        throw cet::exception("LambdaVertexProducer") << "No Slice Data Products Found!" << std::endl;

    art::FindManyP<recob::Hit> hitAssoc = art::FindManyP<recob::Hit>(sliceHandle, evt, sliceInputTag);
    art::fill_ptr_vector(sliceVector, sliceHandle);

    // Find slice that contains the nu hierarchy hits
    std::map<int, int> sliceSignalHitMap;

    int trueNuSliceID(-1);
    int highestHitNumber(-1);
    int highestHitSliceID(-1);

    for (art::Ptr<recob::Slice> &slice : sliceVector)
    {
        sliceSignalHitMap[slice->ID()] = 0;

        const std::vector<art::Ptr<recob::Hit>> &sliceHits(hitAssoc.at(slice.key()));

        for (const art::Ptr<recob::Hit> &sliceHit : sliceHits)
        {
            if (trueHitMap.find(sliceHit.key()) == trueHitMap.end())
                continue;

            int trackID = trueHitMap.at(sliceHit.key());

            if ((trackID == muonTrackID) || (trackID == protonTrackID) || 
                (trackID == pionTrackID) || (trackID == gammaTrackID))
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
        trueNuSliceID = highestHitSliceID;

    return trueNuSliceID;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::LambdaVertexProducer::CollectChildren(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, 
    lar_pandora::PFParticleMap &pfparticleMap, std::vector<art::Ptr<recob::PFParticle>> &collectedParticles, 
    std::vector<art::Ptr<recob::Hit>> &collectedHits)
{
    bool found = false;

    for (const art::Ptr<recob::PFParticle> &collectedPFP : collectedParticles)
    {
        if (collectedPFP->Self() == pfparticle->Self())
            found = true;
    }

    if (!found)
    {
        collectedParticles.push_back(pfparticle);

        std::vector<art::Ptr<recob::Hit>> hits;
        this->CollectHits(e, pfparticle, hits);

        collectedHits.insert(collectedHits.end(), hits.begin(), hits.end());
    }

    const int nChildren = pfparticle->NumDaughters();

    for (int i = 0; i < nChildren; ++i)
    {
        const int childIndex = pfparticle->Daughter(i);
        const art::Ptr<recob::PFParticle> childPFParticle = pfparticleMap.at(childIndex);

        this->CollectChildren(e, childPFParticle, pfparticleMap, collectedParticles, collectedHits);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::LambdaVertexProducer::CollectHits(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, 
   std::vector<art::Ptr<recob::Hit>> &hits)
{
   art::InputTag pfpInputTag(f_PFParticleLabel, f_PFPInstanceName);
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

   if (!e.getByLabel(pfpInputTag, pfparticleHandle))
       throw cet::exception("LambdaVertexProducer") << "No PFParticle Data Products Found!" << std::endl;

   art::InputTag clusterInputTag(f_ClusterLabel, f_PFPInstanceName);
   art::Handle<std::vector<recob::Cluster>> clusterHandle;

   if (!e.getByLabel(clusterInputTag, clusterHandle)) 
       throw cet::exception("LambdaVertexProducer") << "No Cluster Data Products Found!" << std::endl;

   art::FindManyP<recob::Cluster> pfparticleClustersAssoc = art::FindManyP<recob::Cluster>(pfparticleHandle, e, pfpInputTag);
   art::FindManyP<recob::Hit> clusterHitAssoc = art::FindManyP<recob::Hit>(clusterHandle, e, clusterInputTag);

   std::vector<art::Ptr<recob::Cluster>> clusters = pfparticleClustersAssoc.at(pfparticle.key());

   for (const art::Ptr<recob::Cluster> cluster : clusters)
   {
       std::vector<art::Ptr<recob::Hit>> clusterHits = clusterHitAssoc.at(cluster.key());


       std::cout << "clusterHits.size(): " << clusterHits.size() << std::endl;

       hits.insert(hits.end(), clusterHits.begin(), clusterHits.end());
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::LambdaVertexProducer::GetMatchedTrackID(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, 
    const std::map<int, int> &trueHitMap, int &matchedTrackID)   
{
    int nMaxHits = -1;
    bool matchFound = false;
    std::unordered_map<int, int> mcParticleToHitsMap;

    std::vector<art::Ptr<recob::Hit>> hits;
    this->CollectHits(e, pfparticle, hits);

    for (const art::Ptr<recob::Hit> &hit : hits)
    {
        if (trueHitMap.find(hit.key()) == trueHitMap.end())
            continue;

        int thisTrackID = trueHitMap.at(hit.key());

        mcParticleToHitsMap[thisTrackID]++;

        if (mcParticleToHitsMap[thisTrackID] > nMaxHits)
        {
            nMaxHits = mcParticleToHitsMap[thisTrackID];
            matchedTrackID = thisTrackID;
            matchFound = true;
        }
    }

    return matchFound;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool hyperon::LambdaVertexProducer::GetDecayVertex(art::Event& e, std::vector<art::Ptr<recob::PFParticle>> &protonPFParticles, 
    std::vector<art::Ptr<recob::PFParticle>> &pionPFParticles, double &decayVertexX, double &decayVertexY, double &decayVertexZ)
{
    // Get SCE shifted true neutrino vertex
    art::InputTag GenInputTag("generator");
    art::ValidHandle<std::vector<simb::MCTruth>> inputMCTruth = e.getValidHandle<std::vector<simb::MCTruth>>(GenInputTag);
    std::vector<art::Ptr<simb::MCTruth>> inputMCVector;

    art::fill_ptr_vector(inputMCVector, inputMCTruth);

    if (inputMCVector.size() < 1) 
        return false;

    const art::Ptr<simb::MCTruth> mcTruth = inputMCVector.at(0);

    double nuVertexX = mcTruth->GetNeutrino().Nu().Vx();
    double nuVertexY = mcTruth->GetNeutrino().Nu().Vy();
    double nuVertexZ = mcTruth->GetNeutrino().Nu().Vz();

    // Need to apply the SCE correction
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>(); 
    auto scecorr = SCE->GetPosOffsets({nuVertexX, nuVertexY, nuVertexZ});

    double xOffset = -scecorr.X()+0.6;
    double yOffset = scecorr.Y();
    double zOffset = scecorr.Z();
    recob::tracking::Point_t shiftedVertex(nuVertexX + xOffset, nuVertexY + yOffset, nuVertexZ + zOffset);

    // Now find best position
    double closestSeparationSq(std::numeric_limits<double>::max());
    bool found = false;

    art::InputTag pfpInputTag(f_PFParticleLabel, f_PFPInstanceName);
    art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;

    if (!e.getByLabel(pfpInputTag, pfparticleHandle))
        throw cet::exception("LambdaVertexProducer") << "No PFParticle Data Products Found!" << std::endl;

    art::FindManyP<recob::Vertex> pfparticleVertexAssoc = art::FindManyP<recob::Vertex>(pfparticleHandle, e, pfpInputTag);

    for (art::Ptr<recob::PFParticle> pfparticle : protonPFParticles)
    {
        std::vector<art::Ptr<recob::Vertex>> vertexVector(pfparticleVertexAssoc.at(pfparticle.key()));

        if (vertexVector.size() < 1)
            continue;

        art::Ptr<recob::Vertex> pfparticleVertex = vertexVector.at(0);

        double x = pfparticleVertex->position().X();
        double y = pfparticleVertex->position().Y();
        double z = pfparticleVertex->position().Z();
        double dx = x - nuVertexX;
        double dy = y - nuVertexY;
        double dz = z - nuVertexZ;
        double separationSq = (dx * dx) + (dy * dy) + (dz * dz);

        if (separationSq < closestSeparationSq)
        {
            decayVertexX = x;
            decayVertexY = y;
            decayVertexZ = z;
            closestSeparationSq = separationSq;
            found = true;
        }
    }

    for (art::Ptr<recob::PFParticle> pfparticle : pionPFParticles)
    {
        std::vector<art::Ptr<recob::Vertex>> vertexVector(pfparticleVertexAssoc.at(pfparticle.key()));

        if (vertexVector.size() < 1)
            continue;

        art::Ptr<recob::Vertex> pfparticleVertex = vertexVector.at(0);

        double x = pfparticleVertex->position().X();
        double y = pfparticleVertex->position().Y();
        double z = pfparticleVertex->position().Z();
        double dx = x - nuVertexX;
        double dy = y - nuVertexY;
        double dz = z - nuVertexZ;
        double separationSq = (dx * dx) + (dy * dy) + (dz * dz);

        if (separationSq < closestSeparationSq)
        {
            decayVertexX = x;
            decayVertexY = y;
            decayVertexZ = z;
            closestSeparationSq = separationSq;
            found = true;
        }
    }


   return found;
}


DEFINE_ART_MODULE(hyperon::LambdaVertexProducer)
