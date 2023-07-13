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

#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"

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
      void CollectChildren(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, lar_pandora::PFParticleMap &pfparticleMap,
          std::vector<art::Ptr<recob::PFParticle>> &collectedParticles, std::vector<art::Ptr<recob::Hit>> &collectedHits);
      void CollectHits(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, std::vector<art::Ptr<recob::Hit>> &hits);
      void GetTruthMatch(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, const std::vector<art::Ptr<simb::MCParticle>> &mcParticleVector, 
          simb::MCParticle const* &matchedParticle);
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
   f_PFParticleLabel(p.get<std::string>("PFParticleLabel", "pandoraPatRec")),
   f_HitLabel(p.get<std::string>("HitLabel", "gaushit")),
   f_ClusterLabel(p.get<std::string>("ClusterLabel", "pandoraPatRec")),
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
   }

   
   std::cout << "isSignalSigmaZero: " << (isSignalSigmaZero ? "yes" : "no") << std::endl;
   /*
   if (!isSignalSigmaZero)
   {
       e.put(std::move(outputVertex));
       e.put(std::move(outputHitVector));
       e.put(std::move(outputPFParticleVector));
       e.put(std::move(outputSliceVertexAssoc));
       e.put(std::move(outputSliceHitAssoc));
       return;
   }
   */
   ///////////////////////////////////////////
   // Now check whether we need to modify the reco
   ///////////////////////////////////////////

   // Get event PFParticles
   art::Handle<std::vector<recob::PFParticle>> pfparticleHandle;
   std::vector<art::Ptr<recob::PFParticle>> pfparticleVector;

   art::InputTag pfpInputTag(f_PFParticleLabel, f_PFPInstanceName);
   if (!e.getByLabel(pfpInputTag, pfparticleHandle))
       throw cet::exception("LambdaVertexProducer") << "No PFParticle Data Products Found!" << std::endl;

   art::fill_ptr_vector(pfparticleVector, pfparticleHandle);

   // For truth matching
   art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
   std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

   if (!e.getByLabel(f_MCParticleLabel, mcParticleHandle)) 
       throw cet::exception("LambdaVertexProducer") << "No MCParticle Data Products Found!" << std::endl;

   art::fill_ptr_vector(mcParticleVector, mcParticleHandle);

   lar_pandora::MCParticleMap mcParticleMap;
   lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, mcParticleMap);

   // Now fill the Pandora PFParticle map
   lar_pandora::PFParticleMap pfparticleMap;
   lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfparticleVector, pfparticleMap);

   // Now search for the PFParticle matches
   std::vector<art::Ptr<recob::PFParticle>> pionPFParticles;
   std::vector<art::Ptr<recob::PFParticle>> protonPFParticles;

   for (const art::Ptr<recob::PFParticle> pfparticle : pfparticleVector)
   {
       // Make sure that particle is in the neutrino hierarchy (there should only be one)
       if (lar_pandora::LArPandoraHelper::GetParentNeutrino(pfparticleMap, pfparticle) == 0)
           continue;

       // Truth matching 
       simb::MCParticle const* matchedParticle = NULL;
       this->GetTruthMatch(e, pfparticle, mcParticleVector, matchedParticle);

       // The neutrino will have zero hits..
       if (!matchedParticle)
           continue;

       // If its an EM particle, we have to move up the EM hierarchy (who designs this, honestly)
       art::Ptr<simb::MCParticle> motherMCParticle = art::Ptr(art::ProductID(10000), matchedParticle, art::Ptr<simb::MCParticle>::key_type(1000));

       bool found = true;
       do
       {
           const int motherID = motherMCParticle->Mother();

           if (mcParticleMap.find(motherID) == mcParticleMap.end())
           {
               found = false;
               break;
           }

           motherMCParticle = mcParticleMap.at(motherID);
       }
       while ((std::abs(motherMCParticle->PdgCode()) == 11) || (motherMCParticle->PdgCode() == 22));

       if (!found)
       {
           std::cout << "couldn't find mother for: " << matchedParticle->PdgCode() << std::endl;
           continue;
       }


       int motherPDG = motherMCParticle->PdgCode();

       std::cout << "matchedParticle->PdgCode(): " << matchedParticle->PdgCode() << std::endl;
       std::cout << "matchedParticle->Mother(): " << motherPDG << std::endl;


       // Now check if it is one of our lambda0 decay products...
       if (motherPDG != 3122)
           continue;

       if (matchedParticle->PdgCode() == 2212)
           protonPFParticles.push_back(pfparticle);
       else if(matchedParticle->PdgCode() == -211)
           pionPFParticles.push_back(pfparticle);
   }

   const bool foundPion = (pionPFParticles.size() != 0);
   const bool foundProton = (protonPFParticles.size() != 0);
   //const bool modifyReco = (foundPion && !foundProton) || (!foundPion && foundProton);
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

   // Get true vertex
   /*
   simb::MCParticle const* matchedDecayPrimary = NULL;

   if (foundProton)
       this->GetTruthMatch(e, protonPFParticles.front(), mcParticleVector, matchedDecayPrimary);

   if (foundPion)
       this->GetTruthMatch(e, pionPFParticles.front(), mcParticleVector, matchedDecayPrimary);

   std::cout << "evt.isRealData(): " << e.isRealData() << std::endl;
    std::cout << "matchedDecayPrimary->Vx(): " << matchedDecayPrimary->Vx()<< std::endl;
    std::cout << "matchedDecayPrimary->Vy(): " << matchedDecayPrimary->Vy()<< std::endl;
    std::cout << "matchedDecayPrimary->Vz(): " << matchedDecayPrimary->Vz()<< std::endl;
   */

    double decayVertexX(-999.0), decayVertexY(-999.0), decayVertexZ(-999.0);
    if (!this->GetDecayVertex(e, protonPFParticles, pionPFParticles, decayVertexX, decayVertexY, decayVertexZ))
    {
        std::cout << "could not find decay vertex..." << std::endl;
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
   {
       outputPFParticleVector->push_back(pfparticle.key());
       std::cout << "Generation: " << lar_pandora::LArPandoraHelper::GetGeneration(pfparticleMap, pfparticle) << std::endl;
       std::cout << "pfparticle: " << pfparticle << std::endl;
       std::cout << "pfparticle.key(): " << pfparticle.key() << std::endl;
   }

   // Now we need the true vertex.. 
   std::cout << "protonPFParticles.size(): " << protonPFParticles.size() << std::endl;
   std::cout << "pionPFParticles.size(): " << pionPFParticles.size() << std::endl;
   std::cout << "lambdaHierarachyPFParticles.size(): " << lambdaHierarachyPFParticles.size() << std::endl;
   std::cout << "lambdaHierarachyHits.size(): " << lambdaHierarachyHits.size() << std::endl;
   std::cout << "outputHitVector.size(): " << outputHitVector->size() << std::endl;

   e.put(std::move(outputVertex));
   e.put(std::move(outputHitVector));
   e.put(std::move(outputPFParticleVector));
   e.put(std::move(outputSliceVertexAssoc));
   e.put(std::move(outputSliceHitAssoc));

   return;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::LambdaVertexProducer::CollectChildren(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, 
    lar_pandora::PFParticleMap &pfparticleMap, std::vector<art::Ptr<recob::PFParticle>> &collectedParticles, 
    std::vector<art::Ptr<recob::Hit>> &collectedHits)
{
    if (std::find(collectedParticles.begin(), collectedParticles.end(), pfparticle) == collectedParticles.end())
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
       hits.insert(hits.end(), clusterHits.begin(), clusterHits.end());
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hyperon::LambdaVertexProducer::GetTruthMatch(art::Event& e, const art::Ptr<recob::PFParticle> &pfparticle, 
   const std::vector<art::Ptr<simb::MCParticle>> &mcParticleVector, simb::MCParticle const* &matchedParticle)
{
   art::Handle<std::vector<recob::Hit>> hitHandle;

   if (!e.getByLabel(f_HitLabel, hitHandle)) 
       throw cet::exception("LambdaVertexProducer") << "No Hit Data Products Found!" << std::endl;

   art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> truthMatching = art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitHandle, e, f_TruthMatchingLabel);

    int nMaxHits = -1;
    std::unordered_map<int, int> mcParticleToHitsMap;

    std::vector<art::Ptr<recob::Hit>> hits;
    this->CollectHits(e, pfparticle, hits);

    for (const art::Ptr<recob::Hit> &hit : hits)
    {
        std::vector<simb::MCParticle const*> mcParticleVector;
        std::vector<anab::BackTrackerHitMatchingData const*> matchInfoVector;

        truthMatching.get(hit.key(), mcParticleVector, matchInfoVector);

        for (simb::MCParticle const* pMCParticle : mcParticleVector)
        {
            mcParticleToHitsMap[pMCParticle->TrackId()]++;

            if (mcParticleToHitsMap[pMCParticle->TrackId()] > nMaxHits)
            {
                nMaxHits = mcParticleToHitsMap[pMCParticle->TrackId()];
                matchedParticle = pMCParticle;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
