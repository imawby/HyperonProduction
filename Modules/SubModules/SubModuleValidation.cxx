#ifndef _ValidationVariables_cxx_
#define _ValidationVariables_cxx_

#include "SubModuleValidation.h"

using namespace hyperon;

int MUON_INDEX = 0;
int PROTON_INDEX = 1;
int PION_INDEX = 2;
int GAMMA_INDEX = 3;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleValidation::SubModuleValidation(art::Event const& evt, bool isData, fhicl::ParameterSet pset) :
    m_isData(isData),
    m_MCParticleLabel(pset.get<std::string>("MCParticleModuleLabel")),
    m_HitLabel(pset.get<std::string>("HitModuleLabel")),
    m_PandoraModuleLabel(pset.get<std::string>("PandoraModuleLabel")),
    m_PandoraInstanceLabel(pset.get<std::string>("PandoraInstanceLabel")),
    m_BacktrackLabel(pset.get<std::string>("BacktrackModuleLabel"))
{

    const art::InputTag pandoraTag(m_PandoraModuleLabel, m_PandoraInstanceLabel);

    if (!evt.getByLabel(m_MCParticleLabel, m_mcParticleHandle))
        throw cet::exception("SubModuleValidation") << "No MCParticle Data Products Found!" << std::endl;

    if (!evt.getByLabel(m_HitLabel, m_hitHandle))
        throw cet::exception("SubModuleValidation") << "No Hit Data Products Found!" << std::endl;

    if (!evt.getByLabel(pandoraTag, m_sliceHandle))
        throw cet::exception("SubModuleValidation") << "No Slice Data Products Found!" << std::endl;

    if (!evt.getByLabel(pandoraTag, m_pfpHandle))
        throw cet::exception("SubModuleValidation") << "No PFParticle Data Products Found!" << std::endl;

    if (!evt.getByLabel(pandoraTag, m_clusterHandle))
        throw cet::exception("SubModuleValidation") << "No Cluster Data Products Found!" << std::endl;

    // Fill those vectors!
    art::fill_ptr_vector(m_mcParticleVector, m_mcParticleHandle);
    art::fill_ptr_vector(m_hitVector, m_hitHandle);
    art::fill_ptr_vector(m_sliceVector, m_sliceHandle);
    art::fill_ptr_vector(m_pfpVector, m_pfpHandle);

    // Associations
    m_assocMCPartHit = new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(m_hitHandle, evt, m_BacktrackLabel);
    m_assocHitSlice = new art::FindManyP<recob::Hit>(m_sliceHandle, evt, pandoraTag);
    m_assocClusterPFP = new art::FindManyP<recob::Cluster>(m_pfpHandle, evt, pandoraTag);
    m_assocHitCluster = new art::FindManyP<recob::Hit>(m_clusterHandle, evt, pandoraTag);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleValidation::PrepareInfo(SubModuleGeneratorTruth *const gen_SM, SubModuleG4Truth *const g4_SM, 
    std::vector<std::vector<RecoParticle>> &trackPrimaries, std::vector<std::vector<RecoParticle>> &showerPrimaries)
{
    FillMaps();
    GetSigmaZeroTrackIDs(gen_SM, g4_SM);
    FillMCSliceInfo();

    for (unsigned int iSlice = 0; iSlice < trackPrimaries.size(); ++iSlice)
    {
        std::vector<RecoParticle> &sliceTrackPrimaries(trackPrimaries.at(iSlice));

        for (RecoParticle &recoParticle : sliceTrackPrimaries)
            FindMCParticleMatch(g4_SM, recoParticle);

        std::vector<RecoParticle> &sliceShowerPrimaries(showerPrimaries.at(iSlice));

        for (RecoParticle &recoParticle : sliceShowerPrimaries)
            FindMCParticleMatch(g4_SM, recoParticle);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////// 

void SubModuleValidation::FillMaps()
{
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(m_mcParticleVector, m_mcParticleMap);

    FillMCParticleHitInfo();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleValidation::GetSigmaZeroTrackIDs(SubModuleGeneratorTruth *const gen_SM, SubModuleG4Truth *const g4_SM)
{
    SetSigmaZeroIsSignal(gen_SM, g4_SM);

    if (!m_validationData.m_isSignalSigmaZero)
        return;

    m_validationData.m_trueMuonTrackID = IdentifySignalParticle(g4_SM, MUON_INDEX);
    m_validationData.m_trueProtonTrackID = IdentifySignalParticle(g4_SM, PROTON_INDEX);
    m_validationData.m_truePionTrackID = IdentifySignalParticle(g4_SM, PION_INDEX);
    m_validationData.m_trueGammaTrackID = IdentifySignalParticle(g4_SM, GAMMA_INDEX);
}

///////////////////////////////////////////////////////////////////////////////////////////

void SubModuleValidation::SetSigmaZeroIsSignal(SubModuleGeneratorTruth *const gen_SM, SubModuleG4Truth *const g4_SM)
{
    bool isSignalSigmaZero = false;

    for (int i = 0; i < gen_SM->GetGeneratorTruth().NMCTruths; i++)
    {
        if (gen_SM->GetGeneratorTruth().Mode.at(i) != "HYP")
            continue;

        if (!g4_SM->GetG4Info().InActiveTPC.at(i))
            continue;

        if (gen_SM->GetGeneratorTruth().Neutrino.at(i).PDG != -14)
            continue;

        if (!g4_SM->GetG4Info().IsSigmaZeroCharged.at(i))
            continue;

        if (g4_SM->GetG4Info().IsAssociatedHyperon.at(i))
            continue;

        isSignalSigmaZero = true;
        m_validationData.m_mcTruthIndex = i;

        break;
    }

    m_validationData.m_isSignalSigmaZero = isSignalSigmaZero;
}

///////////////////////////////////////////////////////////////////////////////////////////

int SubModuleValidation::IdentifySignalParticle(SubModuleG4Truth *const g4_SM, const int particleTypeIndex)
{
    const std::vector<SimParticle> &simParticles(particleTypeIndex == MUON_INDEX ? g4_SM->GetG4Info().Lepton : 
        particleTypeIndex == PROTON_INDEX ? g4_SM->GetG4Info().Decay : particleTypeIndex == PION_INDEX ? g4_SM->GetG4Info().Decay : 
        g4_SM->GetG4Info().SigmaZeroDecayPhoton);

    const int pdg(particleTypeIndex == MUON_INDEX ? -13 : particleTypeIndex == PROTON_INDEX ? 2212 : particleTypeIndex == PION_INDEX ? -211 : 22);

    int trueTrackID = -1;

    for (const SimParticle &simParticle : simParticles)
    {
        if (simParticle.MCTruthIndex != m_validationData.m_mcTruthIndex)
            continue;

        if (simParticle.PDG != pdg)
            continue;

        if (m_mcParticleMap.find(simParticle.ArtID) == m_mcParticleMap.end())
            continue;

        trueTrackID = simParticle.ArtID;

        break;
    }

    return trueTrackID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleValidation::FillMCParticleHitInfo()
{
    // Truth match to found IDs
    for (unsigned int hitIndex = 0; hitIndex < m_hitVector.size(); hitIndex++)
    {
        const art::Ptr<recob::Hit> &hit = m_hitVector[hitIndex];
        const std::vector<art::Ptr<simb::MCParticle>> &matchedMCParticles = m_assocMCPartHit->at(hit.key());
        auto matchedDatas = m_assocMCPartHit->data(hit.key());

        for (unsigned int mcParticleIndex = 0; mcParticleIndex < matchedMCParticles.size(); ++mcParticleIndex)
        {
            const art::Ptr<simb::MCParticle> &matchedMCParticle = matchedMCParticles.at(mcParticleIndex);
            auto matchedData = matchedDatas.at(mcParticleIndex);

            if (matchedData->isMaxIDE != 1)
                continue;

            const int trackID = IsEM(matchedMCParticle) ? GetLeadEMTrackID(matchedMCParticle) : matchedMCParticle->TrackId();

            m_hitToTrackIDMap[hit.key()] = trackID;
            m_trackIDToHitMap[trackID].push_back(hit.key());
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleValidation::IsEM(const art::Ptr<simb::MCParticle> &mcParticle)
{
    return ((std::abs(mcParticle->PdgCode()) == 11) || (mcParticle->PdgCode() == 22));
}

///////////////////////////////////////////////////////////////////////////////////////////

int SubModuleValidation::GetLeadEMTrackID(const art::Ptr<simb::MCParticle> &mcParticle)
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

void SubModuleValidation::FillMCSliceInfo()
{
    int highestHitNumber(-1);
    int highestHitSliceID(-1);

    for (art::Ptr<recob::Slice> &slice : m_sliceVector)
    {
        int nuHits = 0;

        const std::vector<art::Ptr<recob::Hit>> &sliceHits(m_assocHitSlice->at(slice.key()));

        for (const art::Ptr<recob::Hit> &sliceHit : sliceHits)
        {
            if (m_hitToTrackIDMap.find(sliceHit.key()) == m_hitToTrackIDMap.end())
                continue;

            const int trackID = m_hitToTrackIDMap.at(sliceHit.key());

            if (IsInNuHierarchy(trackID))
                ++nuHits;
        }

        if ((nuHits > highestHitNumber) && (nuHits > 0))
        {
            highestHitNumber = nuHits;
            highestHitSliceID = slice->ID();
        }
    }

    if (highestHitSliceID >= 0)
        m_validationData.m_trueNuSliceID = highestHitSliceID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleValidation::IsInNuHierarchy(const int trackID)
{
    int motherID = trackID;

    while (m_mcParticleMap.find(motherID) != m_mcParticleMap.end())
    {
        const art::Ptr<simb::MCParticle> motherMCParticle = m_mcParticleMap.at(motherID);
        const int pdg(std::abs(motherMCParticle->PdgCode()));

        if ((pdg == 12) || (pdg == 14) || (pdg == 16))
            return true;

        motherID = motherMCParticle->Mother();
    }

    return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleValidation::FindMCParticleMatch(const SubModuleG4Truth* g4_SM, RecoParticle &recoParticle)
{
    if (recoParticle.Key < 0)
        return;

    const art::Ptr<recob::PFParticle> pfp(m_pfpVector.at(recoParticle.Key));

    std::vector<art::Ptr<recob::Hit>> pfpHits;
    CollectHitsFromClusters(pfp, pfpHits);

    std::map<int, int> trackIDToCountMap;

    for (const art::Ptr<recob::Hit> pfpHit : pfpHits)
    {
        if (m_hitToTrackIDMap.find(pfpHit.key()) == m_hitToTrackIDMap.end())
            continue;

        const int trackID(m_hitToTrackIDMap.at(pfpHit.key()));

        if (trackIDToCountMap.find(trackID) == trackIDToCountMap.end())
            trackIDToCountMap[trackID] = 1;
        else
            ++trackIDToCountMap[trackID];
    }

    int maxHits = -1;
    int maxOwnerID = -1;

    for (auto &entry : trackIDToCountMap)
    {
        if ((entry.second > maxHits) || ((entry.second == maxHits) && (entry.first > maxOwnerID)))
        {
            maxHits = entry.second;
            maxOwnerID = entry.first;
        }
    }

    if (maxOwnerID < 0)
    {
        recoParticle.HasTruth = false;
        return;
    }

    const int totalMCHits = m_trackIDToHitMap.at(maxOwnerID).size();
    const double completeness = static_cast<double>(maxHits) / totalMCHits;
    const double purity = static_cast<double>(maxHits) / pfpHits.size();

    const art::Ptr<simb::MCParticle> &mcParticle(m_mcParticleMap.at(maxOwnerID));

    recoParticle.HasTruth = true;
    recoParticle.MCTruthIndex = maxOwnerID;
    recoParticle.TrackTruePDG = mcParticle->PdgCode(); 
    recoParticle.TrackTrueE = mcParticle->Momentum().E();
    recoParticle.TrackTruePx = mcParticle->Momentum().X();
    recoParticle.TrackTruePy = mcParticle->Momentum().Y();
    recoParticle.TrackTruePz = mcParticle->Momentum().Z();
    recoParticle.TrackTrueModMomentum = sqrt((recoParticle.TrackTruePx * recoParticle.TrackTruePx) +
                                             (recoParticle.TrackTruePy * recoParticle.TrackTruePy) +
                                             (recoParticle.TrackTruePz * recoParticle.TrackTruePz));
    recoParticle.TrackTrueKE = recoParticle.TrackTrueE - mcParticle->Mass();


    recoParticle.TrackTrueEndE = mcParticle->EndMomentum().E();
    recoParticle.TrackTrueEndPx = mcParticle->EndMomentum().X();
    recoParticle.TrackTrueEndPy = mcParticle->EndMomentum().Y();
    recoParticle.TrackTrueEndPz = mcParticle->EndMomentum().Z();
    recoParticle.TrackTrueEndModMomentum = sqrt((recoParticle.TrackTrueEndPx * recoParticle.TrackTrueEndPx) +
                                                (recoParticle.TrackTrueEndPy * recoParticle.TrackTrueEndPy) +
                                                (recoParticle.TrackTrueEndPz * recoParticle.TrackTrueEndPz));
    recoParticle.TrackTrueEndKE = recoParticle.TrackTrueEndE - mcParticle->Mass();

    recoParticle.TrackTrueOrigin = g4_SM->GetOrigin(mcParticle->TrackId());

    const double trueLengthDeltaX(mcParticle->Position().X() - mcParticle->EndPosition().X());
    const double trueLengthDeltaY(mcParticle->Position().Y() - mcParticle->EndPosition().Y());
    const double trueLengthDeltaZ(mcParticle->Position().Z() - mcParticle->EndPosition().Z());

    recoParticle.TrackTrueLength = sqrt((trueLengthDeltaX * trueLengthDeltaX) + (trueLengthDeltaY * trueLengthDeltaY) + 
                                        (trueLengthDeltaZ * trueLengthDeltaZ));
    recoParticle.TrackTruthPurity = purity;
    recoParticle.TrackTruthCompleteness = completeness;
    recoParticle.NMatchedHits = maxHits;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleValidation::CollectHitsFromClusters(const art::Ptr<recob::PFParticle> &pfparticle, 
   std::vector<art::Ptr<recob::Hit>> &hits)
{
   std::vector<art::Ptr<recob::Cluster>> clusters = m_assocClusterPFP->at(pfparticle.key());

   for (const art::Ptr<recob::Cluster> cluster : clusters)
   {
       std::vector<art::Ptr<recob::Hit>> clusterHits = m_assocHitCluster->at(cluster.key());
       hits.insert(hits.end(), clusterHits.begin(), clusterHits.end());
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

ValidationData SubModuleValidation::GetInfo()
{
    return m_validationData;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
