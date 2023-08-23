#ifndef _SubModuleValidation_h_
#define _SubModuleValidation_h_

#include <string>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"

// Pandora is the best
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace hyperon {

struct ValidationData {

    int m_trueNuSliceID;
    int m_mcTruthIndex;
    bool m_isSignalSigmaZero;
    int m_trueMuonTrackID;
    int m_trueProtonTrackID;
    int m_truePionTrackID;
    int m_trueGammaTrackID;
};

class SubModuleValidation {

   public:
      SubModuleValidation(art::Event const& evt, bool isData, fhicl::ParameterSet pset);
      void PrepareInfo(SubModuleGeneratorTruth *const gen_SM, SubModuleG4Truth *const g4_SM, 
          std::vector<std::vector<RecoParticle>> &trackPrimaries, std::vector<std::vector<RecoParticle>> &showerPrimaries);
      ValidationData GetInfo();

   private:
      void FillMaps();
      void GetSigmaZeroTrackIDs(SubModuleGeneratorTruth *const gen_SM, SubModuleG4Truth *const g4_SM);
      void SetSigmaZeroIsSignal(SubModuleGeneratorTruth *const gen_SM, SubModuleG4Truth *const g4_SM);
      int IdentifySignalParticle(SubModuleG4Truth *const g4_SM, const int particleTypeIndex);
      void FillMCParticleHitInfo();
      bool IsEM(const art::Ptr<simb::MCParticle> &mcParticle);
      int GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle);
      void FillMCSliceInfo();
      bool IsInNuHierarchy(const int trackID, const std::vector<int> &nuChildrenIDs);
      void FindMCParticleMatch(const SubModuleG4Truth* G4T, RecoParticle &recoParticle);
      void CollectHitsFromClusters(const art::Ptr<recob::PFParticle> &pfparticle, 
          std::vector<art::Ptr<recob::Hit>> &hits);

      art::Handle<std::vector<simb::MCTruth>> m_mcTruthHandle;
      art::Handle<std::vector<simb::MCParticle>> m_mcParticleHandle;
      art::Handle<std::vector<recob::Hit>> m_hitHandle;
      art::Handle<std::vector<recob::Slice>> m_sliceHandle;
      art::Handle<std::vector<recob::PFParticle>> m_pfpHandle;
      art::Handle<std::vector<recob::Cluster>> m_clusterHandle;

      std::vector<art::Ptr<simb::MCTruth>> m_mcTruthVector;
      std::vector<art::Ptr<simb::MCParticle>> m_mcParticleVector;
      std::vector<art::Ptr<recob::Hit>> m_hitVector;
      std::vector<art::Ptr<recob::Slice>> m_sliceVector;
      std::vector<art::Ptr<recob::PFParticle>> m_pfpVector;

      art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>* m_assocMCPartHit;
      art::FindManyP<recob::Hit>* m_assocHitSlice;
      art::FindManyP<recob::Cluster>* m_assocClusterPFP;
      art::FindManyP<recob::Hit>* m_assocHitCluster;

      // Not Configurable
      ValidationData m_validationData;
      lar_pandora::MCParticleMap m_mcParticleMap;
      std::map<int, int> m_hitToTrackIDMap;
      std::map<int, std::vector<int>> m_trackIDToHitMap;

      // Configurables
      bool m_isData;
      std::string m_MCTruthLabel;
      std::string m_MCParticleLabel;
      std::string m_HitLabel;
      std::string m_PandoraModuleLabel;
      std::string m_PandoraInstanceLabel;
      std::string m_BacktrackLabel;
};

}

#endif
