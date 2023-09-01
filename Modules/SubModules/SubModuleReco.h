#ifndef _SubModuleReco_h_
#define _SubModuleReco_h_

#include <string>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "ubana/HyperonProduction/Headers/ParticleTypes.h"
#include "ubana/HyperonProduction/Headers/LLR_PID.h"
#include "ubana/HyperonProduction/Headers/LLRPID_proton_muon_lookup.h"
#include "ubana/HyperonProduction/Headers/LLR_PID_K.h"
#include "ubana/HyperonProduction/Headers/LLRPID_kaon_proton_lookup.h"
#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"
#include "ubana/HyperonProduction/Alg/PIDManager.h"
#include "ubana/HyperonProduction/Alg/BDTHandle.h"

// Pandora is the best
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TVector3.h"

using std::string;

namespace hyperon {

struct RecoData {

   // Information for each slice 
   int FlashMatchedNuSliceID;
   int PandoraNuSliceID;
   std::vector<int> SliceID;
   std::vector<TVector3> RecoPrimaryVertex;
   std::vector<int> NPrimaryDaughters; 
   std::vector<int> NPrimaryTrackDaughters;
   std::vector<int> NPrimaryShowerDaughters;
   std::vector<std::vector<RecoParticle>> TrackPrimaryDaughters;
   std::vector<std::vector<RecoParticle>> ShowerPrimaryDaughters;
   std::vector<std::vector<TVector3>> TrackStarts;
};

class SubModuleReco {

   public:
      SubModuleReco(art::Event const& e,bool isdata, fhicl::ParameterSet pset);

      void PrepareInfo(); 
      bool GetPrimaryVertex(const art::Ptr<recob::Slice> &slice, TVector3 &nuVertex3D);
      void FillPrimaryInfo(const art::Ptr<recob::Slice> &slice, const bool useRepassLabels, const TVector3 &nuVertex3D, 
          std::vector<RecoParticle> &trackPrimaries, std::vector<RecoParticle> &showerPrimaries, std::vector<TVector3> &trackStarts);
      void SetIndices(std::vector<bool> IsSignal,std::vector<bool> IsSignalSigmaZero);
      RecoData GetInfo();
      void SetResRangeCutoff(double cutoff){ m_resRangeCutoff = cutoff; }

   private:
      art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle_SingleOutcome;
      std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle_SingleOutcome;

      art::Handle<std::vector<recob::Slice>> Handle_Slice;
      std::vector<art::Ptr<recob::Slice>> Vect_Slice;
      art::Handle<std::vector<recob::Slice>> Handle_Slice_Reprocessed;
      std::vector<art::Ptr<recob::Slice>> Vect_Slice_Reprocessed;

      art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
      std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;
      art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle_Reprocessed;
      std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle_Reprocessed;
      art::Handle<std::vector<int>> Handle_PFParticle_Modified;
      std::vector<int> Vect_PFParticle_Modified;

      art::Handle<std::vector<recob::Track>> Handle_Track;
      std::vector<art::Ptr<recob::Track>> Vect_Track;
      art::Handle<std::vector<recob::Track>> Handle_Track_Reprocessed;
      std::vector<art::Ptr<recob::Track>> Vect_Track_Reprocessed;

      art::Handle<std::vector<recob::Shower>> Handle_Shower;
      std::vector<art::Ptr<recob::Shower>> Vect_Shower;
      art::Handle<std::vector<recob::Shower>> Handle_Shower_Reprocessed;
      std::vector<art::Ptr<recob::Shower>> Vect_Shower_Reprocessed;

      art::Handle<std::vector<recob::Hit>> Handle_Hit;
      std::vector<art::Ptr<recob::Hit>> Vect_Hit;
      art::Handle<std::vector<recob::Hit>> Handle_Hit_Reprocessed;
      std::vector<art::Ptr<recob::Hit>> Vect_Hit_Reprocessed;

      art::FindManyP<recob::Slice>* Assoc_PFParticleSlice_SingleOutcome; 

      art::FindManyP<recob::PFParticle>* Assoc_SlicePFParticle;
      art::FindManyP<recob::Vertex>* Assoc_PFParticleVertex;
      art::FindManyP<recob::Track>* Assoc_PFParticleTrack;
      art::FindManyP<recob::Shower>* Assoc_PFParticleShower;
      art::FindManyP<larpandoraobj::PFParticleMetadata>* Assoc_PFParticleMetadata;
      art::FindManyP<recob::Hit>* Assoc_TrackHit;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* Assoc_MCParticleBacktracker;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit;
      art::FindManyP<anab::Calorimetry>* Assoc_TrackCalo;
      art::FindManyP<anab::ParticleID>* Assoc_TrackPID;

      art::FindManyP<recob::PFParticle>* Assoc_SlicePFParticle_Reprocessed;
      art::FindManyP<recob::Vertex>* Assoc_PFParticleVertex_Reprocessed;
      art::FindManyP<recob::Track>* Assoc_PFParticleTrack_Reprocessed;
      art::FindManyP<recob::Shower>* Assoc_PFParticleShower_Reprocessed;
      art::FindManyP<larpandoraobj::PFParticleMetadata>* Assoc_PFParticleMetadata_Reprocessed;
      art::FindManyP<recob::Hit>* Assoc_TrackHit_Reprocessed;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit_Reprocessed;
      art::FindManyP<anab::Calorimetry>* Assoc_TrackCalo_Reprocessed;
      art::FindManyP<anab::ParticleID>* Assoc_TrackPID_Reprocessed;

      searchingfornues::LLRPID llr_pid_calculator;
      searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;

      searchingfornuesk::LLRPIDK llr_pid_calculator_kaon;
      searchingfornuesk::KaonProtonLookUpParameters kaonproton_parameters;

      RecoData theData;

      PIDManager m_PIDCalc;
      bool m_isData;
      std::string m_GeneratorModuleLabel;
      std::string m_G4ModuleLabel;
      std::string m_PandoraSingleOutcomeModuleLabel;
      std::string m_PandoraModuleLabel;
      std::string m_PandoraInstanceLabel;
      std::string m_TrackModuleLabel;
      std::string m_ShowerModuleLabel;
      std::string m_PIDModuleLabel;
      std::string m_CaloModuleLabel;
      std::string m_HitModuleLabel;
      bool m_modifiedReco;
      std::string m_PandoraModuleLabelReprocessed;
      std::string m_PandoraInstanceLabelReprocessed;
      std::string m_TrackModuleLabelReprocessed;
      std::string m_ShowerModuleLabelReprocessed;
      std::string m_PIDModuleLabelReprocessed;
      std::string m_CaloModuleLabelReprocessed;
      std::string m_HitModuleLabelReprocessed;
      std::string m_ModifiedPFPListLabel;
      bool m_doGetPIDs;
      double m_resRangeCutoff; 
      bool m_includeCosmics;
      bool m_reducedFileSize;

      int GetFlashMatchedNuSliceID();
      int GetPandoraNuSliceID();
      RecoParticle MakeRecoParticle(const art::Ptr<recob::PFParticle> &pfp, const bool useRepassLabels,
          const TVector3 &nuVertex3D);
      void GetPFPMetadata(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P, const bool useRepassLabels);
      void GetTrackData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P, const bool useRepassLabels);
      void GetPIDs(const art::Ptr<recob::Track> &trk,RecoParticle &P, const bool useRepassLabels);
      void GetVertexData(const art::Ptr<recob::PFParticle> &pfp, RecoParticle &P, const bool useRepassLabels,
          const TVector3 &nuVertex3D);
      void GetVertexData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P, const bool useRepassLabels);
};

}

#endif


