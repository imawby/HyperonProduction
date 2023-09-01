#ifndef _RecoVariables_cxx_
#define _RecoVariables_cxx_

#include "SubModuleReco.h"

using namespace hyperon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleReco::SubModuleReco(art::Event const& e, bool isdata, fhicl::ParameterSet pset) :
              m_isData(isdata),
              m_GeneratorModuleLabel(pset.get<std::string>("GeneratorModuleLabel")),
              m_G4ModuleLabel(pset.get<std::string>("G4ModuleLabel")),
              m_PandoraSingleOutcomeModuleLabel(pset.get<std::string>("PandoraSingleOutcomeModuleLabel")),
              m_PandoraModuleLabel(pset.get<std::string>("PandoraModuleLabel")),
              m_PandoraInstanceLabel(pset.get<std::string>("PandoraInstanceLabel", "")),
              m_TrackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
              m_ShowerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
              m_PIDModuleLabel(pset.get<std::string>("PIDModuleLabel")),
              m_CaloModuleLabel(pset.get<std::string>("CaloModuleLabel")),
              m_HitModuleLabel(pset.get<std::string>("HitModuleLabel")),
              m_modifiedReco(pset.get<bool>("ModifiedReco", false)),
              m_PandoraModuleLabelReprocessed(m_modifiedReco ? pset.get<std::string>("PandoraModuleLabelReprocessed") : ""),
              m_PandoraInstanceLabelReprocessed(pset.get<std::string>("PandoraInstanceLabelReprocessed", "")),
              m_TrackModuleLabelReprocessed(m_modifiedReco ? pset.get<std::string>("TrackModuleLabelReprocessed") : ""),
              m_ShowerModuleLabelReprocessed(m_modifiedReco ? pset.get<std::string>("ShowerModuleLabelReprocessed") : ""),
              m_PIDModuleLabelReprocessed(m_modifiedReco ? pset.get<std::string>("PIDModuleLabelReprocessed") : ""),
              m_CaloModuleLabelReprocessed(m_modifiedReco ? pset.get<std::string>("CaloModuleLabelReprocessed") : ""),
              m_HitModuleLabelReprocessed(m_modifiedReco ? pset.get<std::string>("HitModuleLabelReprocessed") : ""),
              m_ModifiedPFPListLabel(m_modifiedReco ? pset.get<std::string>("ModifiedPFPListLabel") : ""),
              m_doGetPIDs(pset.get<bool>("DoGetPIDs", true)),
              m_resRangeCutoff(pset.get<double>("ResRangeCutoff", 5.0)),
              m_reducedFileSize(pset.get<bool>("ReducedFileSize", false))
{
   if (!e.getByLabel(m_PandoraSingleOutcomeModuleLabel, Handle_PFParticle_SingleOutcome))
       throw cet::exception("SubModuleReco") << "No Single Outcome PFParticle Data Products Found!" << std::endl;

   const art::InputTag pandoraTag(m_PandoraModuleLabel, m_PandoraInstanceLabel);

   if(!e.getByLabel(pandoraTag, Handle_PFParticle)) 
      throw cet::exception("SubModuleReco") << "No PFParticle Data Products Found!" << std::endl;

   if(!e.getByLabel(pandoraTag, Handle_Slice)) 
      throw cet::exception("SubModuleReco") << "No Slice Data Products Found!" << std::endl;

   if(!e.getByLabel(m_TrackModuleLabel, Handle_Track)) 
      throw cet::exception("SubModuleReco") << "No Track Data Products Found!" << std::endl;

   if(!e.getByLabel(m_ShowerModuleLabel, Handle_Shower)) 
      throw cet::exception("SubModuleReco") << "No Shower Data Products Found!" << std::endl;

   if(!e.getByLabel(m_HitModuleLabel, Handle_Hit)) 
      throw cet::exception("SubModuleReco") << "No Hit Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_PFParticle_SingleOutcome, Handle_PFParticle_SingleOutcome);
   art::fill_ptr_vector(Vect_PFParticle, Handle_PFParticle);
   art::fill_ptr_vector(Vect_Slice, Handle_Slice);
   art::fill_ptr_vector(Vect_Track, Handle_Track);
   art::fill_ptr_vector(Vect_Shower, Handle_Shower);
   art::fill_ptr_vector(Vect_Hit, Handle_Hit);

   Assoc_PFParticleSlice_SingleOutcome = new art::FindManyP<recob::Slice>(Vect_PFParticle_SingleOutcome, e, m_PandoraSingleOutcomeModuleLabel);

   Assoc_SlicePFParticle = new art::FindManyP<recob::PFParticle>(Vect_Slice, e, pandoraTag);
   Assoc_PFParticleVertex = new art::FindManyP<recob::Vertex>(Vect_PFParticle, e, pandoraTag);
   Assoc_PFParticleTrack = new art::FindManyP<recob::Track>(Vect_PFParticle, e, m_TrackModuleLabel);    
   Assoc_PFParticleShower = new art::FindManyP<recob::Shower>(Vect_PFParticle, e, m_ShowerModuleLabel);    
   Assoc_PFParticleMetadata = new art::FindManyP<larpandoraobj::PFParticleMetadata>(Vect_PFParticle, e, pandoraTag);   
   Assoc_TrackHit = new  art::FindManyP<recob::Hit>(Vect_Track, e, m_TrackModuleLabel);

   if(m_doGetPIDs)
   {
      Assoc_TrackCalo = new art::FindManyP<anab::Calorimetry>(Vect_Track, e, m_CaloModuleLabel);
      Assoc_TrackPID = new art::FindManyP<anab::ParticleID>(Vect_Track, e, m_PIDModuleLabel);
   }

   if (m_modifiedReco)
   {
       const art::InputTag pandoraTagReprocessed(m_PandoraModuleLabelReprocessed, m_PandoraInstanceLabelReprocessed);

       e.getByLabel(m_ModifiedPFPListLabel, Handle_PFParticle_Modified);
       e.getByLabel(pandoraTagReprocessed, Handle_PFParticle_Reprocessed);
       e.getByLabel(pandoraTagReprocessed, Handle_Slice_Reprocessed);
       e.getByLabel(m_TrackModuleLabelReprocessed, Handle_Track_Reprocessed);
       e.getByLabel(m_ShowerModuleLabelReprocessed, Handle_Shower_Reprocessed);
       e.getByLabel(m_HitModuleLabelReprocessed, Handle_Hit_Reprocessed);

       art::fill_ptr_vector(Vect_PFParticle_Reprocessed, Handle_PFParticle_Reprocessed);
       art::fill_ptr_vector(Vect_Slice_Reprocessed, Handle_Slice_Reprocessed);
       art::fill_ptr_vector(Vect_Track_Reprocessed, Handle_Track_Reprocessed);
       art::fill_ptr_vector(Vect_Shower_Reprocessed, Handle_Shower_Reprocessed);
       art::fill_ptr_vector(Vect_Hit_Reprocessed, Handle_Hit_Reprocessed);

       Assoc_SlicePFParticle_Reprocessed = new art::FindManyP<recob::PFParticle>(Vect_Slice_Reprocessed, e, pandoraTagReprocessed);
       Assoc_PFParticleVertex_Reprocessed = new art::FindManyP<recob::Vertex>(Vect_PFParticle_Reprocessed, e, pandoraTagReprocessed);
       Assoc_PFParticleTrack_Reprocessed = new art::FindManyP<recob::Track>(Vect_PFParticle_Reprocessed, e, m_TrackModuleLabelReprocessed);
       Assoc_PFParticleShower_Reprocessed = new art::FindManyP<recob::Shower>(Vect_PFParticle_Reprocessed, e, m_ShowerModuleLabelReprocessed);   
       Assoc_PFParticleMetadata_Reprocessed = new art::FindManyP<larpandoraobj::PFParticleMetadata>(Vect_PFParticle_Reprocessed, e, pandoraTagReprocessed);
       Assoc_TrackHit_Reprocessed = new  art::FindManyP<recob::Hit>(Vect_Track_Reprocessed, e, m_TrackModuleLabelReprocessed);
    
       if(m_doGetPIDs)
       {
           Assoc_TrackCalo_Reprocessed = new art::FindManyP<anab::Calorimetry>(Vect_Track_Reprocessed, e, m_CaloModuleLabelReprocessed);
           Assoc_TrackPID_Reprocessed = new art::FindManyP<anab::ParticleID>(Vect_Track_Reprocessed, e, m_PIDModuleLabelReprocessed);
       }
   }

   llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
   llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
   llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

   llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
   llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
   llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

   llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
   llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
   llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

   llr_pid_calculator_kaon.set_dedx_binning(0, kaonproton_parameters.dedx_edges_pl_0);
   llr_pid_calculator_kaon.set_par_binning(0, kaonproton_parameters.parameters_edges_pl_0);
   llr_pid_calculator_kaon.set_lookup_tables(0, kaonproton_parameters.dedx_pdf_pl_0);

   llr_pid_calculator_kaon.set_dedx_binning(1, kaonproton_parameters.dedx_edges_pl_1);
   llr_pid_calculator_kaon.set_par_binning(1, kaonproton_parameters.parameters_edges_pl_1);
   llr_pid_calculator_kaon.set_lookup_tables(1, kaonproton_parameters.dedx_pdf_pl_1);

   llr_pid_calculator_kaon.set_dedx_binning(2, kaonproton_parameters.dedx_edges_pl_2);
   llr_pid_calculator_kaon.set_par_binning(2, kaonproton_parameters.parameters_edges_pl_2);
   llr_pid_calculator_kaon.set_lookup_tables(2, kaonproton_parameters.dedx_pdf_pl_2);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::PrepareInfo()
{
   // Get chosen nu slice ID
   theData.FlashMatchedNuSliceID = GetFlashMatchedNuSliceID();
   theData.PandoraNuSliceID = GetPandoraNuSliceID();

   // Go through slices... 
   for (const art::Ptr<recob::Slice> &slice : Vect_Slice)
   {
       TVector3 nuVertex3D(-999.0, -999.0, -999.0);
       if (!GetPrimaryVertex(slice, nuVertex3D))
           continue;

       std::vector<RecoParticle> trackPrimaries, showerPrimaries;
       std::vector<TVector3> trackStarts;

       this->FillPrimaryInfo(slice, false, nuVertex3D, trackPrimaries, showerPrimaries, trackStarts);

       // Look at reprocessed reco
       if (m_modifiedReco)
           this->FillPrimaryInfo(slice, true, nuVertex3D, trackPrimaries, showerPrimaries, trackStarts);

       int nTrackPrimaries = trackPrimaries.size();
       int nShowerPrimaries = showerPrimaries.size();

       if ((nTrackPrimaries == 0) && (nShowerPrimaries == 0))
           continue;

       if (m_reducedFileSize)
       {
           if ((nTrackPrimaries < 3) || (nShowerPrimaries < 1))
           {
               /*
               nTrackPrimaries = 0;
               nShowerPrimaries = 0;
               trackPrimaries.clear();
               showerPrimaries.clear();
               trackStarts.clear();
               */

               continue;
           }
       }

       theData.SliceID.push_back(slice->ID());
       theData.RecoPrimaryVertex.push_back(nuVertex3D);
       theData.TrackPrimaryDaughters.push_back(trackPrimaries);
       theData.ShowerPrimaryDaughters.push_back(showerPrimaries);
       theData.NPrimaryDaughters.push_back(nTrackPrimaries + nShowerPrimaries);
       theData.NPrimaryTrackDaughters.push_back(nTrackPrimaries);
       theData.NPrimaryShowerDaughters.push_back(nShowerPrimaries);
       theData.TrackStarts.push_back(trackStarts);
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

int SubModuleReco::GetFlashMatchedNuSliceID()
{
    int nNeutrinos = 0;
    int nuSliceID = -1;

    for (const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle_SingleOutcome)
    {
        if (pfp->IsPrimary() && isNeutrino(pfp->PdgCode()))
        {
            ++nNeutrinos;

            const std::vector<art::Ptr<recob::Slice>> nuSlice = Assoc_PFParticleSlice_SingleOutcome->at(pfp.key());

            if (nuSlice.size() != 1)
                continue;

            nuSliceID = nuSlice.at(0)->ID();
        }
    }

    if (nNeutrinos > 1)
        throw cet::exception("SubModuleReco") << "Too Many Reconstructed Neutrinos In Event!" << std::endl;

    return nuSliceID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

int SubModuleReco::GetPandoraNuSliceID()
{
    double bestTopologicalScore(-std::numeric_limits<double>::max());
    int nuSliceID = -1;

    for (const art::Ptr<recob::Slice> &slice : Vect_Slice)
    {
        const std::vector<art::Ptr<recob::PFParticle>> slicePFPs = Assoc_SlicePFParticle->at(slice.key());

        for (const art::Ptr<recob::PFParticle> &pfp : slicePFPs)
        {
            // only topological score for the primary pfp
            if (!pfp->IsPrimary())
                continue;

            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMeta = Assoc_PFParticleMetadata->at(pfp.key());

            if (pfpMeta.empty())
                continue;

            const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfpMeta.at(0)->GetPropertiesMap());

            if (!pfParticlePropertiesMap.empty() && (pfParticlePropertiesMap.find("NuScore") != pfParticlePropertiesMap.end()))
            {
                const double thisTopologicalScore = pfParticlePropertiesMap.at("NuScore");

                if (thisTopologicalScore > bestTopologicalScore)
                {
                    bestTopologicalScore = thisTopologicalScore;
                    nuSliceID = slice->ID();
                    break;
                }
            }
        }
    }

    return nuSliceID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleReco::GetPrimaryVertex(const art::Ptr<recob::Slice> &slice, TVector3 &nuVertex3D)
{
   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
   const std::vector<art::Ptr<recob::PFParticle>> slicePFPs = Assoc_SlicePFParticle->at(slice.key());

   for (const art::Ptr<recob::PFParticle> &pfp : slicePFPs)
   {
      if (pfp->IsPrimary() && isNeutrino(pfp->PdgCode()))
      {
         std::vector<art::Ptr<recob::Vertex>> pfpVertex = Assoc_PFParticleVertex->at(pfp.key());

         if (pfpVertex.size() == 0)
             return false;

         const art::Ptr<recob::Vertex> vtx = pfpVertex.at(0);

         geo::Point_t point = {vtx->position().X(), vtx->position().Y(), vtx->position().Z()};
         geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

         nuVertex3D = TVector3(vtx->position().X() + sce_corr.X(), vtx->position().Y() - sce_corr.Y(), vtx->position().Z() - sce_corr.Z());

         return true;
      }
   }

   return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::FillPrimaryInfo(const art::Ptr<recob::Slice> &slice, const bool useRepassLabels, const TVector3 &nuVertex3D,
   std::vector<RecoParticle> &trackPrimaries, std::vector<RecoParticle> &showerPrimaries, std::vector<TVector3> &trackStarts)
{   
   const std::vector<art::Ptr<recob::PFParticle>> slicePFPs = useRepassLabels ? Assoc_SlicePFParticle_Reprocessed->at(slice.key()) : 
       Assoc_SlicePFParticle->at(slice.key());

   lar_pandora::PFParticleMap pfparticleMap;
   lar_pandora::LArPandoraHelper::BuildPFParticleMap(slicePFPs, pfparticleMap);

   int trackCounter = 0;

   for (const art::Ptr<recob::PFParticle> &pfp : slicePFPs)
   {
       if (m_modifiedReco && !useRepassLabels)
       {
           bool modified = false;

           for (long unsigned int i = 0; i < Handle_PFParticle_Modified.product()->size(); ++i)
           {
               if (static_cast<long unsigned int>(Handle_PFParticle_Modified.product()->at(i)) == pfp.key())
               {
                   modified = true;
                   break;
               }
           }

           if (modified)
               continue;
       }

       const int generation = lar_pandora::LArPandoraHelper::GetGeneration(pfparticleMap, pfp);
       const int parentNeutrinoPDG = lar_pandora::LArPandoraHelper::GetParentNeutrino(pfparticleMap, pfp);

       // Do not care about CR
       if (parentNeutrinoPDG == 0)
           continue;

       // Only care about primaries
       if (generation != 2)
           continue;

       RecoParticle P = MakeRecoParticle(pfp, useRepassLabels, nuVertex3D);

       P.Parentage = 1;
       P.InNuSlice = true;
       P.Key = pfp.key();   
       P.Self = pfp->Self();

       if (P.PDG == 13)
       {
           ++trackCounter;
           P.TrackVectorIndex = trackCounter - 1;

           std::vector<art::Ptr<recob::Track>> pfpTracks = useRepassLabels ? Assoc_PFParticleTrack_Reprocessed->at(pfp.key()) : 
               Assoc_PFParticleTrack->at(pfp.key());

           art::Ptr<recob::Track> trk = pfpTracks.at(0);
           trackStarts.push_back(TVector3(trk->Start().X(), trk->Start().Y(), trk->Start().Z()));
           trackPrimaries.push_back(P);
       }
       else if (P.PDG == 11) 
           showerPrimaries.push_back(P);
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoParticle SubModuleReco::MakeRecoParticle(const art::Ptr<recob::PFParticle> &pfp, const bool useRepassLabels,
   const TVector3 &nuVertex3D)
{
   RecoParticle P;

   std::vector<art::Ptr<recob::Track>> pfpTracks = useRepassLabels ? Assoc_PFParticleTrack_Reprocessed->at(pfp.key()) : Assoc_PFParticleTrack->at(pfp.key());
   std::vector<art::Ptr<recob::Shower>> pfpShowers = useRepassLabels ? Assoc_PFParticleShower_Reprocessed->at(pfp.key()) : Assoc_PFParticleShower->at(pfp.key());

   // Make sure it has an associated fitted object
   if ((pfp->PdgCode() == 13 && pfpTracks.size() != 1) || (pfp->PdgCode() == 11 && pfpShowers.size() != 1))
   {
       P.PDG = 0;
       return P;
   }

   P.PDG = pfp->PdgCode();

   GetPFPMetadata(pfp, P, useRepassLabels);

   if (P.PDG == 13)
   {
      GetTrackData(pfp, P, useRepassLabels);
      GetVertexData(pfp, P, useRepassLabels, nuVertex3D);
   }

   return P;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetPFPMetadata(const art::Ptr<recob::PFParticle> &pfp, RecoParticle &P, const bool useRepassLabels){

    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMeta = useRepassLabels ? Assoc_PFParticleMetadata_Reprocessed->at(pfp.key()) : Assoc_PFParticleMetadata->at(pfp.key());

   for(const art::Ptr<larpandoraobj::PFParticleMetadata> &meta : pfpMeta){

      const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(meta->GetPropertiesMap());

      if (!pfParticlePropertiesMap.empty()){
         for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it){
            if(it->first == "TrackScore"){
               P.TrackShowerScore = it->second;
            }
         }
      }	
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetTrackData(const art::Ptr<recob::PFParticle> &pfp, RecoParticle &P, const bool useRepassLabels)
{
   std::vector<art::Ptr<recob::Track>> pfpTracks = useRepassLabels ? Assoc_PFParticleTrack_Reprocessed->at(pfp.key()) : Assoc_PFParticleTrack->at(pfp.key());

   if(pfpTracks.size() != 1) return;

   art::Ptr<recob::Track> trk = pfpTracks.at(0);

   // Sets track length/position related variables
   SetTrackVariables(P, trk);

   if (m_doGetPIDs) GetPIDs(trk, P, useRepassLabels);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetPIDs(const art::Ptr<recob::Track> &trk,RecoParticle &P, const bool useRepassLabels){

   std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = useRepassLabels ? Assoc_TrackCalo_Reprocessed->at(trk.key()) : Assoc_TrackCalo->at(trk.key());
   std::vector<art::Ptr<anab::ParticleID>> trackPID = useRepassLabels ? Assoc_TrackPID_Reprocessed->at(trk.key()) : Assoc_TrackPID->at(trk.key());
   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

   PIDStore store = m_PIDCalc.GetPIDs(trk,caloFromTrack,AlgScoresVec);

   P.Track_LLR_PID = store.LLR;
   P.Track_LLR_PID_Kaon = store.LLR_Kaon;
   P.Track_LLR_PID_Kaon_Partial = store.LLR_Kaon_Partial;
   P.MeandEdX_Plane0 = store.MeandEdX_Plane0;
   P.MeandEdX_Plane1 = store.MeandEdX_Plane1;
   P.MeandEdX_Plane2 = store.MeandEdX_Plane2;
   P.MeandEdX_ThreePlane = store.MeandEdX_3Plane;
   P.Track_Bragg_PID_Kaon = store.Bragg_Kaon_3Plane;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetVertexData(const art::Ptr<recob::PFParticle> &pfp, RecoParticle &P, const bool useRepassLabels,
   const TVector3 &nuVertex3D)
{
   std::vector<art::Ptr<recob::Vertex>> pfpVertex = useRepassLabels ? Assoc_PFParticleVertex_Reprocessed->at(pfp.key()) : Assoc_PFParticleVertex->at(pfp.key());

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   if (pfpVertex.size() == 0)
       return;

   const art::Ptr<recob::Vertex> vtx = pfpVertex.at(0);

   geo::Point_t point = {vtx->position().X(),vtx->position().Y(),vtx->position().Z()};
   geo::Vector_t sce_corr = SCE->GetPosOffsets(point);
   TVector3 pos(vtx->position().X()+sce_corr.X(),vtx->position().Y()-sce_corr.Y(),vtx->position().Z()-sce_corr.Z());

   P.SetVertex(pos);
   P.Displacement = (pos - nuVertex3D).Mag();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoData SubModuleReco::GetInfo(){
   return theData;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
