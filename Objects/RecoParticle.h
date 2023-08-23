#ifndef _RecoParticle_h_
#define _RecoParticle_h_

#include <iostream>

#include "TLorentzVector.h"
#include "TVector3.h"

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

#ifdef __MAKE_ROOT_DICT__
class RecoParticle : public TObject{
#else
class RecoParticle {
#endif

public:

RecoParticle(){}
~RecoParticle(){}

int Self = -1;
bool InNuSlice = false;
int Key = -1;
int TrackVectorIndex = -1;

// General reco info
int PDG = -1; // Pandora PDG code (11 or 13)
int Parentage = -1; // 1 - neutrino daughter, 2 - neutrino granddaughter, 3 - other
int ParentIndex=-1; // -1 - neutrino candidate or no parent
double TrackShowerScore = -999.0;
double X = -999.0;
double Y = -999.0;
double Z = -999.0;
double Displacement = -999.0; // Distance from PV

// Track variables
 double TrackLength = -999.0;
 double TrackDirectionX = -999.0;
 double TrackDirectionY = -999.0;
 double TrackDirectionZ = -999.0;
 double TrackStartX = -999.0;
 double TrackStartY = -999.0;
 double TrackStartZ = -999.0;
 double TrackEndX = -999.0;
 double TrackEndY = -999.0;
 double TrackEndZ = -999.0;
 double TrackPID = -999.0; // 3 plane PID score
 double MeandEdX_Plane0 = -999.0;
 double MeandEdX_Plane1 = -999.0;
 double MeandEdX_Plane2 = -999.0;
 double MeandEdX_ThreePlane = -999.0;
 double Track_LLR_PID = -999.0; // LLR PID
 double Track_LLR_PID_Kaon = -999.0; // LLR PID with Kaon hypothesis
 double Track_LLR_PID_Kaon_Partial = -999.0; // LLR PID with Kaon hypothesis using last 5cm of track
 double Track_Bragg_PID_Kaon = -999.0;
 double ProtonMomentum = -999.0;
 double MuonMomentum = -999.0;
 double KaonMomentum = -999.0; // Track kinematics
 double TrackWiggliness = -999.0;

 // Truth info
 bool HasTruth = false; // False if reco particle has no corresponding MC particle
 int MCTruthIndex = -1;
 int TrackTruePDG = -1;
 double TrackTrueE = -999.0;
 double TrackTruePx = -999.0;
 double TrackTruePy = -999.0;
 double TrackTruePz = -999.0;
 double TrackTrueEndE = -999.0;
 double TrackTrueEndPx = -999.0;
 double TrackTrueEndPy = -999.0;
 double TrackTrueEndPz = -999.0;
 double TrackTrueModMomentum  = -999.0;
 double TrackTrueEndModMomentum = -999.0;
 double TrackTrueKE = -999.0;
 double TrackTrueEndKE = -999.0;
 double TrackTrueLength = -999.0;
 int TrackTrueOrigin = - 1; // 1 - primary , 2 - hyperon decay, 3 - other, 4 - kaon decay, 5 - Sigma0 decay
 double TrackTruthPurity  = -999.0;
 double TrackTruthCompleteness = -999.0;
 int NMatchedHits = -1;

inline void SetVertex(TVector3 V);
inline void SetTrackPositions(TVector3 Start,TVector3 End);
inline void Print();

#ifdef __MAKE_ROOT_DICT__
ClassDef(RecoParticle,1);
#endif

};


inline void RecoParticle::SetVertex(TVector3 V){

   X = V.X();
   Y = V.Y();
   Z = V.Z();

}

inline void RecoParticle::SetTrackPositions(TVector3 Start,TVector3 End){

   TrackStartX = Start.X();
   TrackStartY = Start.Y();
   TrackStartZ = Start.Z();

   TrackEndX = End.X();
   TrackEndY = End.Y();
   TrackEndZ = End.Z();

}

inline void RecoParticle::Print(){

   std::cout << "Reco Info:" << std::endl;
   std::cout << "Self: " << Self << std::endl;
   std::cout << "InNuSlice: " << InNuSlice << std::endl;
   std::cout << "Key: " << Key << std::endl;
   std::cout << "TrackVectorIndex: " << TrackVectorIndex << std::endl;
   std::cout << "PDG Code: " << PDG << "  Track/Shower score: " << TrackShowerScore << std::endl;
   std::cout << "Parentage: " << Parentage << std::endl;
   std::cout << "X: " << X << ", Y: " << Y << ", Z: " << Z << std::endl;
   std::cout << "Displacement: " << Displacement << std::endl;
   std::cout << "Track length: " << TrackLength << "  PID score: " << TrackPID <<  std::endl;
   std::cout << "TrackDirectionX: " << TrackDirectionX << ", TrackDirectionY: " << TrackDirectionY << ", TrackDirectionZ: " << TrackDirectionZ << std::endl;
   std::cout << "MeandEdX_Plane0: " << MeandEdX_Plane0 << ", MeandEdX_Plane1: " << MeandEdX_Plane1 << ", MeandEdX_Plane2: " << MeandEdX_Plane2 << std::endl;
   std::cout << "MeandEdX_ThreePlane: " << MeandEdX_ThreePlane << std::endl;
   std::cout << "Track_LLR_PID: " << Track_LLR_PID << std::endl;
   std::cout << "Track_LLR_PID_Kaon: " << Track_LLR_PID_Kaon << std::endl;
   std::cout << "Track_LLR_PID_Kaon_Partial: " << Track_LLR_PID_Kaon_Partial << std::endl;
   std::cout << "Track_Bragg_PID_Kaon: " << Track_Bragg_PID_Kaon << std::endl;
   std::cout << "ProtonMomentum: " << ProtonMomentum << std::endl;
   std::cout << "MuonMomentum: " << MuonMomentum << std::endl;
   std::cout << "KaonMomentum: " << KaonMomentum << std::endl;
   std::cout << "TrackWiggliness: " << TrackWiggliness << std::endl;
   std::cout << "Truth Info:" << std::endl;
   std::cout << "HasTruth: " << HasTruth << std::endl;
   std::cout << "MCTruthIndex: " << MCTruthIndex << std::endl;
   std::cout << "PDG: " << TrackTruePDG << "  Origin: " << TrackTrueOrigin << std::endl;
   std::cout << "TrackTrueE: " << TrackTrueE << ", TrackTruePx: " << TrackTruePx << ", TrackTruePy: " << TrackTruePy << ", TrackTruePz: " << TrackTruePz << std::endl;
   std::cout << "TrackTrueEndE: " << TrackTrueEndE << ", TrackTrueEndPx: " << TrackTrueEndPx << ", TrackTrueEndPy: " << TrackTrueEndPy << ", TrackTrueEndPz: " << TrackTrueEndPz << std::endl;
   std::cout << "TrackTrueModMomentum: " << TrackTrueModMomentum << std::endl;
   std::cout << "TrackTrueKE: " << TrackTrueKE << ", TrackTrueEndKE: " << TrackTrueEndKE << std::endl;
   std::cout << "TrackTrueLength: " << TrackTrueLength << std::endl;
   std::cout << "TrackTruthPurity: " << TrackTruthPurity << std::endl;
   std::cout << "TrackTruthCompleteness: " << TrackTruthCompleteness << std::endl;
   std::cout << "NMatchedHits: " << NMatchedHits << std::endl;
}

#endif
