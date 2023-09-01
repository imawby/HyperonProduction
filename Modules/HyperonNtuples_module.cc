////////////////////////////////////////////////////////////////////////
// Class:       HyperonNtuples
// Plugin Type: analyzer (art v3_03_01)
// File:        HyperonNtuples_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

//root includes
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

//local includes

//objects and helpers
#include "ubana/HyperonProduction/Objects/SimParticle.h"
#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"
#include "ubana/HyperonProduction/Objects/LinkDef.h"

//algorithms
#include "ubana/HyperonProduction/Alg/ConnectednessHelper.h"

//submodules
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleReco.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleValidation.h"

namespace hyperon {
   class HyperonNtuples;
}


class hyperon::HyperonNtuples : public art::EDAnalyzer {
   public:
      explicit HyperonNtuples(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      HyperonNtuples(HyperonNtuples const&) = delete;
      HyperonNtuples(HyperonNtuples&&) = delete;
      HyperonNtuples& operator=(HyperonNtuples const&) = delete;
      HyperonNtuples& operator=(HyperonNtuples&&) = delete;

      // Required functions.
      void analyze(art::Event const& e) override;

      // Selected optional functions.
      void beginJob() override;
      void endJob() override;

      void FinishEvent();

      //check if event contains a reco'd muon, proton and pion from Lambda decay
      //records their positions in track vector if they exist
      //void StoreTrackTruth();

      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

   private:

      // Output trees
      TTree * OutputTree;
      TTree * MetaTree;

      // Basic event info
      unsigned int t_EventID;
      int t_run,t_subrun,t_event;

      double t_Weight=1.0;

      // Generator/Geant4 truth info

      int t_NMCTruths=0;	
      int t_NMCTruthsInTPC=0;	

      std::vector<std::string> t_Mode;
      std::vector<std::string> t_CCNC;

      // Flags applying to the entire event
      bool t_EventHasNeutronScatter;
      bool t_EventHasHyperon;

      // Flags applying to each MCTruth
      std::vector<bool> t_InActiveTPC;
      std::vector<bool> t_IsHyperon;
      std::vector<bool> t_IsLambda;
      std::vector<bool> t_IsLambdaCharged;
      std::vector<bool> t_IsSigmaZero; 		
      std::vector<bool> t_IsSigmaZeroCharged; 		
      std::vector<bool> t_IsAssociatedHyperon;
      std::vector<bool> t_IsSignal;
      std::vector<bool> t_IsSignalSigmaZero;

      bool t_EventHasFinalStateNeutron;

      std::vector<SimParticle> t_Neutrino;
      std::vector<SimParticle> t_Lepton;
      std::vector<SimParticle> t_Hyperon;
      std::vector<SimParticle> t_PrimaryNucleon;
      std::vector<SimParticle> t_PrimaryPion;
      std::vector<SimParticle> t_PrimaryKaon; 
      std::vector<SimParticle> t_Decay; 
      std::vector<SimParticle> t_SigmaZeroDecayPhoton;
      std::vector<SimParticle> t_SigmaZeroDecayLambda;
      std::vector<SimParticle> t_KaonDecay;

      std::vector<double> t_TruePrimaryVertex_X;
      std::vector<double> t_TruePrimaryVertex_Y;
      std::vector<double> t_TruePrimaryVertex_Z;
      std::vector<double> t_DecayVertex_X;
      std::vector<double> t_DecayVertex_Y;
      std::vector<double> t_DecayVertex_Z;

      ////////////////////////////
      //   Output for each slice
      ////////////////////////////

      int t_FlashMatchedNuSliceID;
      int t_PandoraNuSliceID;
      std::vector<int> t_SliceID;

      std::vector<bool> t_GoodReco;
      std::vector<int> t_NPrimaryDaughters;
      std::vector<int> t_NPrimaryTrackDaughters;
      std::vector<int> t_NPrimaryShowerDaughters;

      std::vector<std::vector<RecoParticle>> t_TrackPrimaryDaughters;
      std::vector<std::vector<RecoParticle>> t_ShowerPrimaryDaughters;   
      //std::vector<TVector3> t_RecoPrimaryVertex;
      std::vector<double> t_RecoPrimaryVertexX;
      std::vector<double> t_RecoPrimaryVertexY;
      std::vector<double> t_RecoPrimaryVertexZ;

      std::vector<std::vector<std::vector<int>>> t_Conn_SeedIndexes_Plane0;
      std::vector<std::vector<std::vector<int>>> t_Conn_OutputIndexes_Plane0;
      std::vector<std::vector<std::vector<int>>> t_Conn_OutputSizes_Plane0;
      std::vector<std::vector<std::vector<int>>> t_Conn_SeedChannels_Plane0;
      std::vector<std::vector<std::vector<int>>> t_Conn_SeedTicks_Plane0;

      std::vector<std::vector<std::vector<int>>> t_Conn_SeedIndexes_Plane1;
      std::vector<std::vector<std::vector<int>>> t_Conn_OutputIndexes_Plane1;
      std::vector<std::vector<std::vector<int>>> t_Conn_OutputSizes_Plane1;
      std::vector<std::vector<std::vector<int>>> t_Conn_SeedChannels_Plane1;
      std::vector<std::vector<std::vector<int>>> t_Conn_SeedTicks_Plane1;

      std::vector<std::vector<std::vector<int>>> t_Conn_SeedIndexes_Plane2;
      std::vector<std::vector<std::vector<int>>> t_Conn_OutputIndexes_Plane2;
      std::vector<std::vector<std::vector<int>>> t_Conn_OutputSizes_Plane2;
      std::vector<std::vector<std::vector<int>>> t_Conn_SeedChannels_Plane2;
      std::vector<std::vector<std::vector<int>>> t_Conn_SeedTicks_Plane2;

      int t_trueNuSliceID;
      int t_trueMuonTrackID;
      int t_trueProtonTrackID;
      int t_truePionTrackID;
      int t_trueGammaTrackID;

      std::vector<std::string> t_SysDials;
      std::vector<std::vector<double>> t_SysWeights;

      /////////////////////////
      // Metadata for sample //
      /////////////////////////

      int m_NEvents;
      int m_NHyperons;
      int m_NSignal;      
      int m_NSignalSigmaZero;      
      int m_NGoodReco;

      double m_POT = 0; //total POT of the sample

      //////////////////////////
      //   FHICL PARAMETERS   //
      //////////////////////////

      bool f_GetGeneratorInfo;
      bool f_GetG4Info;
      bool f_GetRecoInfo;
      bool f_GetConnInfo;
      bool f_GetValidationInfo;

      fhicl::ParameterSet f_Generator;
      fhicl::ParameterSet f_G4;
      fhicl::ParameterSet f_Reco;
      fhicl::ParameterSet f_Validation;
      std::string f_WireLabel;
      std::vector<art::InputTag> f_WeightLabels;
      std::string f_POTSummaryLabel;

      bool f_IsData;
      bool f_Debug = false;

      ///////////////////////
      //      Objects      //
      ///////////////////////

      ConnectednessHelper Conn_Helper;
};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

hyperon::HyperonNtuples::HyperonNtuples(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   f_GetGeneratorInfo(p.get<bool>("GetGeneratorInfo", true)),   
   f_GetG4Info(p.get<bool>("GetG4Info", true)),   
   f_GetRecoInfo(p.get<bool>("GetRecoInfo", true)),   
   f_GetConnInfo(p.get<bool>("GetConnInfo", true)),   
   f_GetValidationInfo(p.get<bool>("GetValidationInfo", true)),
   f_Generator(p.get<fhicl::ParameterSet>("Generator")),
   f_G4(p.get<fhicl::ParameterSet>("Geant4")),
   f_Reco(p.get<fhicl::ParameterSet>("Reco")),
   f_Validation(p.get<fhicl::ParameterSet>("Validation")),
   f_WireLabel(p.get<std::string>("WireLabel")),
   f_WeightLabels(p.get<std::vector<art::InputTag>>("WeightCalculators",{})),
   f_POTSummaryLabel(p.get<std::string>("POTSummaryLabel")),
   f_IsData(p.get<bool>("IsData")),
   f_Debug(p.get<bool>("Debug",false)),
   Conn_Helper(p.get<bool>("DrawConnectedness",false))
{
   if(f_WeightLabels.size()){
      std::cout << "Getting weights from data products with tags:" << std::endl;
      for(size_t i=0;i<f_WeightLabels.size();i++) std::cout << f_WeightLabels.at(i) << std::endl;
   }

   if(f_Reco.get<bool>("IncludeCosmics",false) && f_GetConnInfo)
      std::cout << std::endl << "HyperonNTuples WARNING: Requesting connectedness information with cosmics included. This will take a very long time." << std::endl
                << "Set GetConnInfo fhicl parameter to false disable connectedness" << std::endl << std::endl;
      
}

void hyperon::HyperonNtuples::analyze(art::Event const& e)
{
   if(f_Debug) std::cout << "New Event" << std::endl;

   //begin by resetting everything

   t_Weight = 1.0;
   //t_Mode = "NONE";
   t_Mode.clear();
   t_CCNC.clear();
   t_NMCTruths = 0;
   t_NMCTruthsInTPC = 0;
   t_EventHasFinalStateNeutron = false;

   t_InActiveTPC.clear();
   t_IsHyperon.clear();
   t_IsLambda.clear();
   t_IsLambdaCharged.clear();
   t_IsSigmaZero.clear();
   t_IsSigmaZeroCharged.clear();
   t_IsAssociatedHyperon.clear();
   t_IsSignal.clear();	
   t_IsSignalSigmaZero.clear();	
   t_EventHasNeutronScatter = false;
   t_EventHasHyperon = false;

   t_Neutrino.clear();
   t_Lepton.clear();
   t_Hyperon.clear();
   t_PrimaryNucleon.clear();
   t_PrimaryPion.clear();
   t_PrimaryKaon.clear();
   t_Decay.clear();
   t_KaonDecay.clear();
   t_SigmaZeroDecayPhoton.clear();
   t_SigmaZeroDecayLambda.clear();

   t_TruePrimaryVertex_X.clear();
   t_TruePrimaryVertex_Y.clear();
   t_TruePrimaryVertex_Z.clear();
   t_DecayVertex_X.clear();
   t_DecayVertex_Y.clear();
   t_DecayVertex_Z.clear();

   t_FlashMatchedNuSliceID = -1;
   t_PandoraNuSliceID = -1;
   t_SliceID.clear();

   t_GoodReco.clear();
   t_NPrimaryDaughters.clear(); //number of primary daughters
   t_NPrimaryTrackDaughters.clear(); //num of track like primary daughters
   t_NPrimaryShowerDaughters.clear(); //num of shower like primary daughters

   t_TrackPrimaryDaughters.clear();
   t_ShowerPrimaryDaughters.clear();

   //t_RecoPrimaryVertex.clear(); //.SetXYZ(-1000,-1000,-1000); //position of reco'd primary vertex
   t_RecoPrimaryVertexX.clear();
   t_RecoPrimaryVertexY.clear();
   t_RecoPrimaryVertexZ.clear();

   t_Conn_SeedIndexes_Plane0.clear();
   t_Conn_OutputIndexes_Plane0.clear();
   t_Conn_OutputSizes_Plane0.clear();
   t_Conn_SeedChannels_Plane0.clear();
   t_Conn_SeedTicks_Plane0.clear();

   t_Conn_SeedIndexes_Plane1.clear();
   t_Conn_OutputIndexes_Plane1.clear();
   t_Conn_OutputSizes_Plane1.clear();
   t_Conn_SeedChannels_Plane1.clear();
   t_Conn_SeedTicks_Plane1.clear();

   t_Conn_SeedIndexes_Plane2.clear();
   t_Conn_OutputIndexes_Plane2.clear();
   t_Conn_OutputSizes_Plane2.clear();
   t_Conn_SeedChannels_Plane2.clear();
   t_Conn_SeedTicks_Plane2.clear();

   t_trueNuSliceID = -1;
   t_trueMuonTrackID = -1;
   t_trueProtonTrackID = -1;
   t_truePionTrackID = -1;
   t_trueGammaTrackID = -1;

   t_SysDials.clear();
   t_SysWeights.clear();

   // General Event Info

   t_EventID = e.id().event();
   t_run = e.run();
   t_subrun = e.subRun();
   t_event = e.event();

   //////////////////////////////
   // Event Generator Info
   //////////////////////////////
   SubModuleGeneratorTruth* Generator_SM(nullptr);

   if (!f_IsData && f_GetGeneratorInfo)
   {
      if(f_Debug) std::cout << "Getting EG Info" << std::endl;

      Generator_SM = new SubModuleGeneratorTruth(e,f_Generator);
      GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();

      t_Weight *= GenT.Weight;
      t_CCNC = GenT.CCNC;
      t_Mode = GenT.Mode;
      t_NMCTruths = GenT.NMCTruths;
      t_NMCTruthsInTPC = GenT.NMCTruthsInTPC;
      t_Neutrino = GenT.Neutrino;
      t_TruePrimaryVertex_X = GenT.TruePrimaryVertex_X;
      t_TruePrimaryVertex_Y = GenT.TruePrimaryVertex_Y;
      t_TruePrimaryVertex_Z = GenT.TruePrimaryVertex_Z;

      t_EventHasFinalStateNeutron = GenT.EventHasFinalStateNeutron;
   }

   //////////////////////////////
   // G4 Info
   //////////////////////////////
   SubModuleG4Truth* G4_SM(nullptr);

   if (!f_IsData && f_GetG4Info)
   {
      if(f_Debug) std::cout << "Getting G4 Info" << std::endl;

      G4_SM = new SubModuleG4Truth(e, f_G4);
      G4Truth G4T = G4_SM->GetG4Info();

      t_InActiveTPC = G4T.InActiveTPC;
      t_IsHyperon = G4T.IsHyperon;
      t_IsLambda = G4T.IsLambda;
      t_IsLambdaCharged = G4T.IsLambdaCharged;
      t_IsSigmaZero = G4T.IsSigmaZero;
      t_IsSigmaZeroCharged = G4T.IsSigmaZeroCharged;
      t_IsAssociatedHyperon = G4T.IsAssociatedHyperon;
      t_EventHasNeutronScatter = G4T.EventHasNeutronScatter;       
      t_EventHasHyperon = G4T.EventHasHyperon;       
      t_Weight *= G4T.Weight;
      t_Lepton = G4T.Lepton;
      t_Hyperon = G4T.Hyperon;
      t_PrimaryNucleon = G4T.PrimaryNucleon;
      t_PrimaryPion = G4T.PrimaryPion;
      t_PrimaryKaon = G4T.PrimaryKaon;
      t_Decay = G4T.Decay;
      t_KaonDecay = G4T.KaonDecay;
      t_SigmaZeroDecayPhoton = G4T.SigmaZeroDecayPhoton;
      t_SigmaZeroDecayLambda = G4T.SigmaZeroDecayLambda;
      t_DecayVertex_X = G4T.DecayVertex_X;
      t_DecayVertex_Y = G4T.DecayVertex_Y;
      t_DecayVertex_Z = G4T.DecayVertex_Z;

      t_IsSignal.resize(t_NMCTruths);
      t_IsSignalSigmaZero.resize(t_NMCTruths);

      for(int i_t=0;i_t<t_NMCTruths;i_t++){
         t_IsSignal[i_t] = t_Mode.at(i_t) == "HYP" && t_InActiveTPC.at(i_t) && t_Neutrino.at(i_t).PDG == -14 && t_IsLambdaCharged.at(i_t) && !t_IsAssociatedHyperon.at(i_t);
         t_IsSignalSigmaZero[i_t] = t_Mode.at(i_t) == "HYP" && t_InActiveTPC.at(i_t) && t_Neutrino.at(i_t).PDG == -14 && t_IsSigmaZeroCharged.at(i_t) && !t_IsAssociatedHyperon.at(i_t);
      }
   }

   //////////////////////////////
   // Reconstructed Info
   //////////////////////////////
   SubModuleReco* Reco_SM(nullptr);
   SubModuleValidation* Validation_SM(nullptr);
   if (f_GetRecoInfo)
   {
      if(f_Debug) std::cout << "Getting Reconstructed Info" << std::endl;

      Reco_SM = new SubModuleReco(e,f_IsData,f_Reco);
      Reco_SM->PrepareInfo();
      //Reco_SM->SetIndices(t_IsSignal,t_IsSignalSigmaZero);
      RecoData RecoD =  Reco_SM->GetInfo();   

      t_FlashMatchedNuSliceID = RecoD.FlashMatchedNuSliceID;
      t_PandoraNuSliceID = RecoD.PandoraNuSliceID;
      t_SliceID = RecoD.SliceID;
      t_NPrimaryDaughters = RecoD.NPrimaryDaughters;
      t_NPrimaryTrackDaughters = RecoD.NPrimaryTrackDaughters;
      t_NPrimaryShowerDaughters = RecoD.NPrimaryShowerDaughters;
      t_TrackPrimaryDaughters = RecoD.TrackPrimaryDaughters;
      t_ShowerPrimaryDaughters = RecoD.ShowerPrimaryDaughters;

      for (const TVector3 &primaryVertex : RecoD.RecoPrimaryVertex)
      {
          t_RecoPrimaryVertexX.push_back(primaryVertex.X());
          t_RecoPrimaryVertexY.push_back(primaryVertex.Y());
          t_RecoPrimaryVertexZ.push_back(primaryVertex.Z());
      }

      // Results of connectedness test on different combinations of tracks
      if(f_Debug) std::cout << "Performing Connectedness Tests" << std::endl;

      for (unsigned int iSlice = 0; iSlice < RecoD.SliceID.size(); ++iSlice)
      {
         if(f_GetConnInfo)
         {
             CTOutcome ConnData = Conn_Helper.PrepareAndTestEvent(e, f_WireLabel, RecoD.TrackStarts.at(iSlice));

             t_Conn_SeedIndexes_Plane0.push_back(ConnData.SeedIndexes_Plane0);
             t_Conn_OutputIndexes_Plane0.push_back(ConnData.OutputIndexes_Plane0);
             t_Conn_OutputSizes_Plane0.push_back(ConnData.OutputSizes_Plane0);
             t_Conn_SeedChannels_Plane0.push_back(ConnData.SeedChannels_Plane0);
             t_Conn_SeedTicks_Plane0.push_back(ConnData.SeedTicks_Plane0);
             t_Conn_SeedIndexes_Plane1.push_back(ConnData.SeedIndexes_Plane1);
             t_Conn_OutputIndexes_Plane1.push_back(ConnData.OutputIndexes_Plane1);
             t_Conn_OutputSizes_Plane1.push_back(ConnData.OutputSizes_Plane1);
             t_Conn_SeedChannels_Plane1.push_back(ConnData.SeedChannels_Plane1);
             t_Conn_SeedTicks_Plane1.push_back(ConnData.SeedTicks_Plane1);
             t_Conn_SeedIndexes_Plane2.push_back(ConnData.SeedIndexes_Plane2);
             t_Conn_OutputIndexes_Plane2.push_back(ConnData.OutputIndexes_Plane2);
             t_Conn_OutputSizes_Plane2.push_back(ConnData.OutputSizes_Plane2);
             t_Conn_SeedChannels_Plane2.push_back(ConnData.SeedChannels_Plane2);
             t_Conn_SeedTicks_Plane2.push_back(ConnData.SeedTicks_Plane2);
         }
      }

      //////////////////////////////
      // Validation Info
      //////////////////////////////
      if (!f_IsData && G4_SM && f_GetValidationInfo)
      {
         if(f_Debug) std::cout << "Getting Validation Info" << std::endl;

         Validation_SM = new SubModuleValidation(e, f_IsData, f_Validation);
         Validation_SM->PrepareInfo(Generator_SM, G4_SM, t_TrackPrimaryDaughters, t_ShowerPrimaryDaughters);

         ValidationData ValidationD = Validation_SM->GetInfo();

         t_trueNuSliceID = ValidationD.m_trueNuSliceID;
         t_trueMuonTrackID = ValidationD.m_trueMuonTrackID;
         t_trueProtonTrackID = ValidationD.m_trueProtonTrackID;
         t_truePionTrackID = ValidationD.m_truePionTrackID;
         t_trueGammaTrackID = ValidationD.m_trueGammaTrackID;


         std::cout << "t_FlashMatchedNuSliceID: " << t_FlashMatchedNuSliceID << std::endl;
         std::cout << "t_PandoraNuSliceID: " << t_PandoraNuSliceID << std::endl;
         std::cout << "t_trueNuSliceID: " << t_trueNuSliceID << std::endl;

      }
   }

    // delete everything.. 
    delete Generator_SM;
    delete G4_SM;
    delete Reco_SM;
    delete Validation_SM;


   if(!f_IsData){

      std::vector<std::map<std::string,std::vector<double>>> theweightmap(t_NMCTruths); 

      for(size_t i_w=0;i_w<f_WeightLabels.size();i_w++){

         std::cout << "Getting new weight products with label " << f_WeightLabels.at(i_w) << std::endl;

         art::Handle<std::vector<evwgh::MCEventWeight>> Handle_EventWeight;
         std::vector<art::Ptr<evwgh::MCEventWeight>> Vect_EventWeight;

         if(!e.getByLabel(f_WeightLabels.at(i_w),Handle_EventWeight)) 
            throw cet::exception("HyperonNtuples") << "No EventWeight Found!" << std::endl;

         art::fill_ptr_vector(Vect_EventWeight,Handle_EventWeight);

         if(!Vect_EventWeight.size())
            throw cet::exception("HyperonNtuples") << "Weight vector empty!" << std::endl;

         if(Vect_EventWeight.size() != (size_t)t_NMCTruths)
            throw cet::exception("HyperonNtuples") << "Weight vector size != NMCTruths" << std::endl;

         for(size_t i_tr=0;i_tr<Vect_EventWeight.size();i_tr++){       

            std::cout << "Getting weights for truth " << i_tr << std::endl;

            std::map<std::string,std::vector<double>> theWeights = Vect_EventWeight.at(i_tr)->fWeight;
            std::map<std::string,std::vector<double>>::iterator it;

            for(it = theWeights.begin();it != theWeights.end();it++){

               if(it->first ==  "empty") continue;

               bool dial_found=false;

               std::map<std::string,std::vector<double>>::iterator it2;
               for(it2 = theweightmap.at(i_tr).begin();it2 != theweightmap.at(i_tr).end();it2++){
                  if(it->first == it2->first){
                     dial_found = true;
                     theweightmap.at(i_tr)[it->first].insert(theweightmap.at(i_tr)[it->first].end(),it->second.begin(),it->second.end());
                  }
               }

               if(!dial_found)
                  theweightmap.at(i_tr)[it->first] = it->second;

            }
         } // i_tr
      }

      // Organise the weights
      if(theweightmap.size()){
         std::map<std::string,std::vector<double>>::iterator it;
         for(it = theweightmap.at(0).begin();it != theweightmap.at(0).end();it++){
            std::cout << "Organising weights for dial " << it->first << std::endl;

            t_SysDials.push_back(it->first);
            t_SysWeights.push_back(it->second);

            for(size_t i_tr=1;i_tr<theweightmap.size();i_tr++){
               if(theweightmap.at(i_tr).find(it->first) == theweightmap.at(i_tr).end()) 
                  throw cet::exception("HyperonNtuples") << "Dial " << it->first << " not found in weights for MC truth " << i_tr << std::endl;
               if(theweightmap.at(i_tr)[it->first].size() != t_SysWeights.back().size())
                  throw cet::exception("HyperonNtuples") << "Dial " << it->first << " weight vector mismatch" << std::endl;                        
               for(size_t i_w=0;i_w<t_SysWeights.back().size();i_w++)
                  t_SysWeights.back().at(i_w) *= theweightmap.at(i_tr)[it->first].at(i_w);
            }                  
         }
      }
   }

   FinishEvent();
}

///////////////////////////////////////////////////////////////	
// Finished processing event - update Metadata and fill tree //
///////////////////////////////////////////////////////////////

void hyperon::HyperonNtuples::FinishEvent(){

   if(f_Debug) std::cout << "Finishing Event" << std::endl;

   OutputTree->Fill();

   m_NEvents++;

   if(std::find(t_IsHyperon.begin(), t_IsHyperon.end(), true) != t_IsHyperon.end()) m_NHyperons++;
   if(std::find(t_IsSignal.begin(), t_IsSignal.end(), true) != t_IsSignal.end()) m_NSignal++;
   if(std::find(t_IsSignalSigmaZero.begin(), t_IsSignalSigmaZero.end(), true) != t_IsSignalSigmaZero.end()) m_NSignalSigmaZero++;
   //if(t_GoodReco) m_NGoodReco++;

   if(f_Debug) std::cout << "Finished event" << std::endl;

}

///////////////////////////////////////////////////////////////	

void hyperon::HyperonNtuples::beginJob(){

   if(f_Debug) std::cout << "Begin job" << std::endl;

   art::ServiceHandle<art::TFileService> tfs;

   //////////////////////////////////////////
   //             Output Tree	           //
   //////////////////////////////////////////

   OutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");

   OutputTree->Branch("IsData",&f_IsData);
   OutputTree->Branch("EventID",&t_EventID);
   OutputTree->Branch("run",&t_run);
   OutputTree->Branch("subrun",&t_subrun);
   OutputTree->Branch("event",&t_event);

   OutputTree->Branch("Weight",&t_Weight);
   OutputTree->Branch("Mode","vector<string>",&t_Mode);
   OutputTree->Branch("CCNC","vector<string>",&t_CCNC);
   OutputTree->Branch("NMCTruths",&t_NMCTruths);
   OutputTree->Branch("NMCTruthsInTPC",&t_NMCTruthsInTPC);
   OutputTree->Branch("InActiveTPC","vector<bool>",&t_InActiveTPC);
   OutputTree->Branch("IsHyperon","vector<bool>",&t_IsHyperon);
   OutputTree->Branch("IsLambda","vector<bool>",&t_IsLambda);
   OutputTree->Branch("IsLambdaCharged","vector<bool>",&t_IsLambdaCharged);
   OutputTree->Branch("IsSigmaZero","vector<bool>",&t_IsSigmaZero);
   OutputTree->Branch("IsSigmaZeroCharged","vector<bool>",&t_IsSigmaZeroCharged);
   OutputTree->Branch("IsAssociatedHyperon","vector<bool>",&t_IsAssociatedHyperon);
   OutputTree->Branch("IsSignal","vector<bool>",&t_IsSignal);
   OutputTree->Branch("IsSignalSigmaZero","vector<bool>",&t_IsSignalSigmaZero);
   OutputTree->Branch("GoodReco",&t_GoodReco);
   OutputTree->Branch("EventHasNeutronScatter",&t_EventHasNeutronScatter);
   OutputTree->Branch("EventHasHyperon",&t_EventHasHyperon);
   OutputTree->Branch("EventHasFinalStateNeutron",&t_EventHasFinalStateNeutron);

   OutputTree->Branch("Neutrino","vector<SimParticle>",&t_Neutrino);
   OutputTree->Branch("Lepton","vector<SimParticle>",&t_Lepton);
   OutputTree->Branch("Hyperon","vector<SimParticle>",&t_Hyperon);
   OutputTree->Branch("PrimaryNucleon","vector<SimParticle>",&t_PrimaryNucleon);
   OutputTree->Branch("PrimaryPion","vector<SimParticle>",&t_PrimaryPion);
   OutputTree->Branch("PrimaryKaon","vector<SimParticle>",&t_PrimaryKaon);
   OutputTree->Branch("Decay","vector<SimParticle>",&t_Decay);
   OutputTree->Branch("SigmaZeroDecayPhoton","vector<SimParticle>",&t_SigmaZeroDecayPhoton);
   OutputTree->Branch("SigmaZeroDecayLambda","vector<SimParticle>",&t_SigmaZeroDecayLambda);
   OutputTree->Branch("KaonDecay","vector<SimParticle>",&t_KaonDecay);
   OutputTree->Branch("TruePrimaryVertex_X",&t_TruePrimaryVertex_X);
   OutputTree->Branch("TruePrimaryVertex_Y",&t_TruePrimaryVertex_Y);
   OutputTree->Branch("TruePrimaryVertex_Z",&t_TruePrimaryVertex_Z);

   OutputTree->Branch("DecayVertex_X",&t_DecayVertex_X);
   OutputTree->Branch("DecayVertex_Y",&t_DecayVertex_Y);
   OutputTree->Branch("DecayVertex_Z",&t_DecayVertex_Z);

   OutputTree->Branch("FlashMatchedNuSliceID", &t_FlashMatchedNuSliceID);
   OutputTree->Branch("PandoraNuSliceID", &t_PandoraNuSliceID);
   OutputTree->Branch("SliceID", &t_SliceID);
   OutputTree->Branch("RecoPrimaryVertexX", &t_RecoPrimaryVertexX);
   OutputTree->Branch("RecoPrimaryVertexY", &t_RecoPrimaryVertexY);
   OutputTree->Branch("RecoPrimaryVertexZ", &t_RecoPrimaryVertexZ);
   OutputTree->Branch("NPrimaryTrackDaughters", "std::vector<int>", &t_NPrimaryTrackDaughters);
   OutputTree->Branch("NPrimaryShowerDaughters", "std::vector<int>", &t_NPrimaryShowerDaughters);
   OutputTree->Branch("TracklikePrimaryDaughters", "std::vector<std::vector<RecoParticle>>",&t_TrackPrimaryDaughters);
   OutputTree->Branch("ShowerlikePrimaryDaughters", "std::vector<std::vector<RecoParticle>>",&t_ShowerPrimaryDaughters);

   OutputTree->Branch("ConnSeedIndexes_Plane0",&t_Conn_SeedIndexes_Plane0);
   OutputTree->Branch("ConnOutputIndexes_Plane0",&t_Conn_OutputIndexes_Plane0);
   OutputTree->Branch("ConnOutputSizes_Plane0",&t_Conn_OutputSizes_Plane0);
   OutputTree->Branch("ConnSeedChannels_Plane0",&t_Conn_SeedChannels_Plane0);
   OutputTree->Branch("ConnSeedTicks_Plane0",&t_Conn_SeedTicks_Plane0);
   OutputTree->Branch("ConnSeedIndexes_Plane1",&t_Conn_SeedIndexes_Plane1);
   OutputTree->Branch("ConnOutputIndexes_Plane1",&t_Conn_OutputIndexes_Plane1);
   OutputTree->Branch("ConnOutputSizes_Plane1",&t_Conn_OutputSizes_Plane1);
   OutputTree->Branch("ConnSeedChannels_Plane1",&t_Conn_SeedChannels_Plane1);
   OutputTree->Branch("ConnSeedTicks_Plane1",&t_Conn_SeedTicks_Plane1);
   OutputTree->Branch("ConnSeedIndexes_Plane2",&t_Conn_SeedIndexes_Plane2);
   OutputTree->Branch("ConnOutputIndexes_Plane2",&t_Conn_OutputIndexes_Plane2);
   OutputTree->Branch("ConnOutputSizes_Plane2",&t_Conn_OutputSizes_Plane2);
   OutputTree->Branch("ConnSeedChannels_Plane2",&t_Conn_SeedChannels_Plane2);
   OutputTree->Branch("ConnSeedTicks_Plane2",&t_Conn_SeedTicks_Plane2);

   OutputTree->Branch("TrueNuSliceID", &t_trueNuSliceID);
   OutputTree->Branch("TrueMuonTrackID", &t_trueMuonTrackID);
   OutputTree->Branch("TrueProtonTrackID", &t_trueProtonTrackID);
   OutputTree->Branch("TruePionTrackID", &t_truePionTrackID);
   OutputTree->Branch("TrueGammaTrackID", &t_trueGammaTrackID);

   //OutputTree->Branch("SysDials",&t_SysDials);
   //OutputTree->Branch("SysWeights","vector<vector<vector<double>>>",&t_SysWeights);

   OutputTree->Branch("SysDials",&t_SysDials);
   OutputTree->Branch("SysWeights","vector<vector<double>>",&t_SysWeights);

   //////////////////////////////////////////
   //             Metadata Tree	           //
   //////////////////////////////////////////

   m_NEvents=0;
   m_NHyperons=0;
   m_NSignal=0;
   m_NSignalSigmaZero=0;
   m_NGoodReco=0;
   m_POT=0;

   MetaTree=tfs->make<TTree>("MetaTree","Metadata Info Tree");

   MetaTree->Branch("NEvents",&m_NEvents);
   MetaTree->Branch("NHyperons",&m_NHyperons);
   MetaTree->Branch("NSignal",&m_NSignal);
   MetaTree->Branch("NSignalSigmaZero",&m_NSignalSigmaZero);
   MetaTree->Branch("NGoodReco",&m_NGoodReco);

   MetaTree->Branch("POT",&m_POT);

   if(f_Debug) std::cout << "Finished begin job" << std::endl;

}

void hyperon::HyperonNtuples::endJob()
{
   MetaTree->Fill();
}

void hyperon::HyperonNtuples::beginSubRun(const art::SubRun& sr)
{
   if(f_Debug) std::cout << "Getting Subrun POT Info" << std::endl;

   art::Handle<sumdata::POTSummary> POTHandle;
   if(sr.getByLabel(f_POTSummaryLabel,POTHandle)) m_POT += POTHandle->totpot;	
}

void hyperon::HyperonNtuples::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(hyperon::HyperonNtuples)
