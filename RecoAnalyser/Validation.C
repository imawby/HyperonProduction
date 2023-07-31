const int MUON_INDEX = 0;
const int PROTON_INDEX = 1;
const int PION_INDEX = 2;
const int GAMMA_INDEX = 3;

const double EPSILON = 0.0000000001;
const std::vector<int> LINE_COLOR = {kBlue, kRed, kSpring+3, kBlack};
const std::vector LEGEND_STRING = {"Muon", "Proton", "Pion", "Gamma"};

void DrawHistogramCollection(const std::vector<TH1D*> &histogramCollection, const std::string &title,
    const std::string &xAxis, const std::string &yAxis);
void ObtainEfficiencyHistogramCollection(const std::vector<TH1D*> &weightedHistogramCollection, 
    const std::vector<TH1D*> &fullHistogramCollection);
void NormaliseHistogramCollection(const std::vector<TH1D*> &histogramCollection);
void DrawSliceEfficiency(const std::vector<TH1D*> &fullHistogramCollection, const std::vector<TH1D*> &thisSliceEffCollection,
    const std::vector<TH1D*> &otherSliceEffCollection, const int particleTypeIndex);

void Validation()
{
    gStyle->SetOptStat(0000);


    TFile * file = new TFile("sigmaRecoAna.root", "READ");
    TTree * tree = (TTree*)file->Get("sigmaRecoAnalyser/sigmaAna");

    std::vector<int> *t_nTrueHits(nullptr);
    std::vector<int> *t_matchFoundInSlice(nullptr);
    std::vector<int> *t_matchFoundInOtherSlice(nullptr);
    std::vector<double> *t_trueNuVertexSep(nullptr);
    std::vector<double> *t_trackScore(nullptr);
    std::vector<double> *t_completeness(nullptr);
    std::vector<double> *t_purity(nullptr);
    std::vector<int> *t_generation(nullptr);
    int t_nParticlesFoundInSlice;
    int t_protonPionMerged;

    tree->SetBranchAddress("NTrueHits", &t_nTrueHits);
    tree->SetBranchAddress("MatchFoundInSlice", &t_matchFoundInSlice);
    tree->SetBranchAddress("MatchFoundInOtherSlice", &t_matchFoundInOtherSlice);
    tree->SetBranchAddress("TrueNuVertexSep", &t_trueNuVertexSep);
    tree->SetBranchAddress("BestMatchTrackScore", &t_trackScore);
    tree->SetBranchAddress("BestMatchCompleteness", &t_completeness);
    tree->SetBranchAddress("BestMatchPurity", &t_purity);
    tree->SetBranchAddress("BestMatchGeneration", &t_generation);
    tree->SetBranchAddress("NParticlesFoundInSlice", &t_nParticlesFoundInSlice);
    tree->SetBranchAddress("ProtonPionMerged", &t_protonPionMerged);

    /////////////////////
    // True histograms
    /////////////////////
    // True hits
    double nTrueHitsMin(-1.5), nTrueHitsMax(199.5), nTrueHitsBin(40);
    TH1D * nTrueHitsDist_0 = new TH1D("nTrueHitsDist_0", "nTrueHitsDist_0", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsDist_1 = new TH1D("nTrueHitsDist_1", "nTrueHitsDist_1", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsDist_2 = new TH1D("nTrueHitsDist_2", "nTrueHitsDist_2", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsDist_3 = new TH1D("nTrueHitsDist_3", "nTrueHitsDist_3", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    std::vector<TH1D*> nTrueHitsDists = {nTrueHitsDist_0, nTrueHitsDist_1, nTrueHitsDist_2, nTrueHitsDist_3};

    TH1D * nTrueHitsThisSliceEff_0 = new TH1D("nTrueHitsThisSliceEff_0", "nTrueHitsThisSliceEff_0", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsThisSliceEff_1 = new TH1D("nTrueHitsThisSliceEff_1", "nTrueHitsThisSliceEff_1", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsThisSliceEff_2 = new TH1D("nTrueHitsThisSliceEff_2", "nTrueHitsThisSliceEff_2", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsThisSliceEff_3 = new TH1D("nTrueHitsThisSliceEff_3", "nTrueHitsThisSliceEff_3", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    std::vector<TH1D*> nTrueHitsThisSliceEff = {nTrueHitsThisSliceEff_0, nTrueHitsThisSliceEff_1, nTrueHitsThisSliceEff_2, nTrueHitsThisSliceEff_3};

    TH1D * nTrueHitsOtherSliceEff_0 = new TH1D("nTrueHitsOtherSliceEff_0", "nTrueHitsOtherSliceEff_0", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsOtherSliceEff_1 = new TH1D("nTrueHitsOtherSliceEff_1", "nTrueHitsOtherSliceEff_1", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsOtherSliceEff_2 = new TH1D("nTrueHitsOtherSliceEff_2", "nTrueHitsOtherSliceEff_2", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    TH1D * nTrueHitsOtherSliceEff_3 = new TH1D("nTrueHitsOtherSliceEff_3", "nTrueHitsOtherSliceEff_3", nTrueHitsBin, nTrueHitsMin, nTrueHitsMax);
    std::vector<TH1D*> nTrueHitsOtherSliceEff = {nTrueHitsOtherSliceEff_0, nTrueHitsOtherSliceEff_1, nTrueHitsOtherSliceEff_2, nTrueHitsOtherSliceEff_3};

    // Distance from nu vertex
    double nTrueDistanceMin(-1.5), nTrueDistanceMax(100.5), nTrueDistanceBin(50);
    TH1D * nTrueDistanceDist_0 = new TH1D("nTrueDistanceDist_0", "nTrueDistanceDist_0", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceDist_1 = new TH1D("nTrueDistanceDist_1", "nTrueDistanceDist_1", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceDist_2 = new TH1D("nTrueDistanceDist_2", "nTrueDistanceDist_2", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceDist_3 = new TH1D("nTrueDistanceDist_3", "nTrueDistanceDist_3", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    std::vector<TH1D*> nTrueDistanceDists = {nTrueDistanceDist_0, nTrueDistanceDist_1, nTrueDistanceDist_2, nTrueDistanceDist_3};

    TH1D * nTrueDistanceThisSliceEff_0 = new TH1D("nTrueDistanceThisSliceEff_0", "nTrueDistanceThisSliceEff_0", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceThisSliceEff_1 = new TH1D("nTrueDistanceThisSliceEff_1", "nTrueDistanceThisSliceEff_1", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceThisSliceEff_2 = new TH1D("nTrueDistanceThisSliceEff_2", "nTrueDistanceThisSliceEff_2", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceThisSliceEff_3 = new TH1D("nTrueDistanceThisSliceEff_3", "nTrueDistanceThisSliceEff_3", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    std::vector<TH1D*> nTrueDistanceThisSliceEff = {nTrueDistanceThisSliceEff_0, nTrueDistanceThisSliceEff_1, nTrueDistanceThisSliceEff_2, nTrueDistanceThisSliceEff_3};

    TH1D * nTrueDistanceOtherSliceEff_0 = new TH1D("nTrueDistanceOtherSliceEff_0", "nTrueDistanceOtherSliceEff_0", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceOtherSliceEff_1 = new TH1D("nTrueDistanceOtherSliceEff_1", "nTrueDistanceOtherSliceEff_1", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceOtherSliceEff_2 = new TH1D("nTrueDistanceOtherSliceEff_2", "nTrueDistanceOtherSliceEff_2", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    TH1D * nTrueDistanceOtherSliceEff_3 = new TH1D("nTrueDistanceOtherSliceEff_3", "nTrueDistanceOtherSliceEff_3", nTrueDistanceBin, nTrueDistanceMin, nTrueDistanceMax);
    std::vector<TH1D*> nTrueDistanceOtherSliceEff = {nTrueDistanceOtherSliceEff_0, nTrueDistanceOtherSliceEff_1, nTrueDistanceOtherSliceEff_2, nTrueDistanceOtherSliceEff_3};

    /////////////////////
    // Reco histograms - best match if in slice...
    /////////////////////
    // Track Score
    double trackScoreMin(-1.05), trackScoreMax(1.05), trackScoreBin(22);
    TH1D * trackScoreDist_0 = new TH1D("trackScoreDist_0", "trackScoreDist_0", trackScoreBin, trackScoreMin, trackScoreMax);
    TH1D * trackScoreDist_1 = new TH1D("trackScoreDist_1", "trackScoreDist_1", trackScoreBin, trackScoreMin, trackScoreMax);
    TH1D * trackScoreDist_2 = new TH1D("trackScoreDist_2", "trackScoreDist_2", trackScoreBin, trackScoreMin, trackScoreMax);
    TH1D * trackScoreDist_3 = new TH1D("trackScoreDist_3", "trackScoreDist_3", trackScoreBin, trackScoreMin, trackScoreMax);
    std::vector<TH1D*> trackScoreDists = {trackScoreDist_0, trackScoreDist_1, trackScoreDist_2, trackScoreDist_3};

    // Completeness
    double completenessMin(-1.05), completenessMax(1.05), completenessBin(22);
    TH1D * completenessDist_0 = new TH1D("completenessDist_0", "completenessDist_0", completenessBin, completenessMin, completenessMax);
    TH1D * completenessDist_1 = new TH1D("completenessDist_1", "completenessDist_1", completenessBin, completenessMin, completenessMax);
    TH1D * completenessDist_2 = new TH1D("completenessDist_2", "completenessDist_2", completenessBin, completenessMin, completenessMax);
    TH1D * completenessDist_3 = new TH1D("completenessDist_3", "completenessDist_3", completenessBin, completenessMin, completenessMax);
    std::vector<TH1D*> completenessDists = {completenessDist_0, completenessDist_1, completenessDist_2, completenessDist_3};

    // Purity
    double purityMin(-1.05), purityMax(1.05), purityBin(22);
    TH1D * purityDist_0 = new TH1D("purityDist_0", "purityDist_0", purityBin, purityMin, purityMax);
    TH1D * purityDist_1 = new TH1D("purityDist_1", "purityDist_1", purityBin, purityMin, purityMax);
    TH1D * purityDist_2 = new TH1D("purityDist_2", "purityDist_2", purityBin, purityMin, purityMax);
    TH1D * purityDist_3 = new TH1D("purityDist_3", "purityDist_3", purityBin, purityMin, purityMax);
    std::vector<TH1D*> purityDists = {purityDist_0, purityDist_1, purityDist_2, purityDist_3};

    double nEvents(tree->GetEntries());
    double nEvents0Shower0Track(0.0);
    double nEvents0Shower1Track(0.0);
    double nEvents0Shower2Track(0.0);
    double nEvents0Shower3Track(0.0);
    double nEvents0Shower4PlusTrack(0.0);
    double nEvents1Shower0Track(0.0);
    double nEvents1Shower1Track(0.0);
    double nEvents1Shower2Track(0.0);
    double nEvents1Shower3Track(0.0);
    double nEvents1Shower4PlusTrack(0.0);
    double nEvents2PlusShower0Track(0.0);
    double nEvents2PlusShower1Track(0.0);
    double nEvents2PlusShower2Track(0.0);
    double nEvents2PlusShower3Track(0.0);
    double nEvents2PlusShower4PlusTrack(0.0);
    std::vector<double> foundThisSlice(4, 0.0);
    std::vector<double> foundOtherSlice(4, 0.0);
    std::vector<double> protonPionMerged(4, 0.0);
    std::vector<double> noSlice(4, 0.0);
    std::vector<double> totalParticle(4, 0.0);
    std::vector<double> thisSliceTrackID(4, 0.0);
    std::vector<double> thisSliceShowerID(4, 0.0);
    std::vector<double> thisSliceGen1(4, 0.0);
    std::vector<double> thisSliceGen2(4, 0.0);
    std::vector<double> thisSliceGen3(4, 0.0);
    std::vector<double> thisSliceGen4Plus(4, 0.0);
    double allThisSlice(0.0), allThisSlicePlusGen(0.0), allThisSlicePlusID(0.0), allThisSlicePlusGenPlusID(0.0);

    // Events where all reco particles have > 15 hits
    double nEventsAboveThreshold(0.0);
    std::vector<double> foundThisSlice_15(4, 0.0);
    std::vector<double> foundOtherSlice_15(4, 0.0);
    std::vector<double> protonPionMerged_15(4, 0.0);
    std::vector<double> noSlice_15(4, 0.0);
    std::vector<double> thisSliceTrackID_15(4, 0.0);
    std::vector<double> thisSliceShowerID_15(4, 0.0);
    std::vector<double> thisSliceGen2_15(4, 0.0);
    std::vector<double> thisSliceGen3_15(4, 0.0);
    std::vector<double> thisSliceGen4Plus_15(4, 0.0);
    double allThisSlice_15(0.0), allThisSlicePlusGen_15(0.0), allThisSlicePlusID_15(0.0), allThisSlicePlusGenPlusID_15(0.0);

    /////////////////////
    // Loop over tree, fill counts and histograms
    /////////////////////
    for (unsigned int treeIndex = 0; treeIndex < tree->GetEntries(); ++treeIndex)
    {
        tree->GetEntry(treeIndex);

        const bool isEventAboveThreshold = (t_nTrueHits->at(MUON_INDEX) > 14) && (t_nTrueHits->at(PROTON_INDEX) > 14) && (t_nTrueHits->at(PION_INDEX) > 14) &&
            (t_nTrueHits->at(GAMMA_INDEX) > 14);

        if (isEventAboveThreshold)
            ++nEventsAboveThreshold;

        int nShowers(0);
        int nTracks(0);

        for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        {
            // Fill histograms
            const double boundedHits = std::min(std::max(nTrueHitsMin + EPSILON, static_cast<double>(t_nTrueHits->at(particleTypeIndex))), nTrueHitsMax - EPSILON);
            const double boundedDistance = std::min(std::max(nTrueDistanceMin + EPSILON, static_cast<double>(t_trueNuVertexSep->at(particleTypeIndex))), nTrueDistanceMax - EPSILON);
            const double boundedTrackScore = std::min(std::max(trackScoreMin + EPSILON, static_cast<double>(t_trackScore->at(particleTypeIndex))), trackScoreMax - EPSILON);

            nTrueHitsDists[particleTypeIndex]->Fill(boundedHits);
            nTrueDistanceDists[particleTypeIndex]->Fill(boundedDistance);

            if (t_matchFoundInSlice->at(particleTypeIndex))
            {
                trackScoreDists[particleTypeIndex]->Fill(boundedTrackScore);
                completenessDists[particleTypeIndex]->Fill(t_completeness->at(particleTypeIndex));
                purityDists[particleTypeIndex]->Fill(t_purity->at(particleTypeIndex));
            }

            nTrueHitsThisSliceEff[particleTypeIndex]->Fill(boundedHits, t_matchFoundInSlice->at(particleTypeIndex));
            nTrueDistanceThisSliceEff[particleTypeIndex]->Fill(boundedDistance, t_matchFoundInSlice->at(particleTypeIndex));

            nTrueHitsOtherSliceEff[particleTypeIndex]->Fill(boundedHits, t_matchFoundInOtherSlice->at(particleTypeIndex));
            nTrueDistanceOtherSliceEff[particleTypeIndex]->Fill(boundedDistance, t_matchFoundInOtherSlice->at(particleTypeIndex));

            // Fill counts
            ++totalParticle[particleTypeIndex];

            if (t_matchFoundInSlice->at(particleTypeIndex))
            {
                const bool isTrackID = (!(t_trackScore->at(particleTypeIndex) < 0.0) && (t_trackScore->at(particleTypeIndex) > 0.5));
                const bool isShowerID = (!(t_trackScore->at(particleTypeIndex) < 0.0) && !(t_trackScore->at(particleTypeIndex) > 0.5));

                if (isTrackID) ++nTracks;
                if (isShowerID) ++nShowers;

                ++foundThisSlice[particleTypeIndex];

                if (isEventAboveThreshold) ++foundThisSlice_15[particleTypeIndex];

                // Track score counts
                if (isTrackID)
                {
                    ++thisSliceTrackID[particleTypeIndex];

                    if (isEventAboveThreshold) ++thisSliceTrackID_15[particleTypeIndex];
                }
                else if (isShowerID)
                {
                    ++thisSliceShowerID[particleTypeIndex];

                    if (isEventAboveThreshold) ++thisSliceShowerID_15[particleTypeIndex];
                }

                // Generation counts
                if (t_generation->at(particleTypeIndex) == 2)
                {
                    ++thisSliceGen2[particleTypeIndex];

                    if (isEventAboveThreshold) ++thisSliceGen2_15[particleTypeIndex];

                }
                else if (t_generation->at(particleTypeIndex) == 3)
                {
                    ++thisSliceGen3[particleTypeIndex];

                    if (isEventAboveThreshold) ++thisSliceGen3_15[particleTypeIndex];
                }
                else if (t_generation->at(particleTypeIndex) > 3)
                {
                    ++thisSliceGen4Plus[particleTypeIndex];

                    if (isEventAboveThreshold) ++thisSliceGen4Plus_15[particleTypeIndex];
                }
            }
            else if (((particleTypeIndex == PROTON_INDEX) || (particleTypeIndex == PION_INDEX)) && t_protonPionMerged && (t_matchFoundInSlice->at(PROTON_INDEX) || t_matchFoundInSlice->at(PION_INDEX)))
            {
                //std::cout << "HERE";
                ++protonPionMerged[particleTypeIndex];

                if (isEventAboveThreshold) ++protonPionMerged_15[particleTypeIndex];
            }
            else if (t_matchFoundInOtherSlice->at(particleTypeIndex))
            {
                ++foundOtherSlice[particleTypeIndex];

                if (isEventAboveThreshold) ++foundOtherSlice_15[particleTypeIndex];
            }
            else
            {
                ++noSlice[particleTypeIndex];

                if (isEventAboveThreshold) ++noSlice_15[particleTypeIndex];
            }
        } //end of particle index loop

        // Track/shower multiplicity counts
        if (isEventAboveThreshold)
            {
        if (nShowers == 0)
        {
            if (nTracks == 0) ++nEvents0Shower0Track;
            if (nTracks == 1) ++nEvents0Shower1Track;
            if (nTracks == 2) ++nEvents0Shower2Track;
            if (nTracks == 3) ++nEvents0Shower3Track;
            if (nTracks > 3) ++nEvents0Shower4PlusTrack;
        }
        if (nShowers == 1)
        {
            if (nTracks == 0) ++nEvents1Shower0Track;
            if (nTracks == 1) ++nEvents1Shower1Track;
            if (nTracks == 2) ++nEvents1Shower2Track;
            if (nTracks == 3) ++nEvents1Shower3Track;
            if (nTracks > 3) ++nEvents1Shower4PlusTrack;
        }
        if (nShowers > 1)
        {
            if (nTracks == 0) ++nEvents2PlusShower0Track;
            if (nTracks == 1) ++nEvents2PlusShower1Track;
            if (nTracks == 2) ++nEvents2PlusShower2Track;
            if (nTracks == 3) ++nEvents2PlusShower3Track;
            if (nTracks > 3) ++nEvents2PlusShower4PlusTrack;
        }
            }
        // Event level correctness counts
        if (t_nParticlesFoundInSlice == 4)
        {
            allThisSlice += 1.0;

            if (isEventAboveThreshold) ++allThisSlice_15;

            const bool correctGen = (t_generation->at(MUON_INDEX) == 2) && (t_generation->at(PROTON_INDEX) == 2) && (t_generation->at(PION_INDEX) == 2) &&
                (t_generation->at(GAMMA_INDEX) == 2);

            if (correctGen)
            {
                ++allThisSlicePlusGen;

                if (isEventAboveThreshold) ++allThisSlicePlusGen_15;
            }

            bool correctID = true;

            //std::cout << "---------------------------------------" << std::endl;

            for (int particleType : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
            {
                const bool trackID = (!(t_trackScore->at(particleType) < 0.0) && (t_trackScore->at(particleType) > 0.5));
                const bool showerID = (!(t_trackScore->at(particleType) < 0.0) && !(t_trackScore->at(particleType) > 0.5));

                //std::cout << "particleType: " << particleType << ", trackID: " << trackID << ", showerID: " << showerID << std::endl;

                if ((particleType == MUON_INDEX) && showerID)
                {
                    correctID = false;
                }
                else if ((particleType == PROTON_INDEX) && showerID)
                {
                    correctID = false;
                }
                else if ((particleType == PION_INDEX) && showerID)
                {
                    correctID = false;
                }
                else if ((particleType == GAMMA_INDEX) && trackID)
                {
                    correctID = false;
                }
            }

            if (correctID)
            {
                ++allThisSlicePlusID;

                if (isEventAboveThreshold) ++allThisSlicePlusID_15;
            }

            if (correctGen && correctID)
            {
                ++allThisSlicePlusGenPlusID;

                if (isEventAboveThreshold) ++allThisSlicePlusGenPlusID_15;
            }
        }
    } // end of tree loop
    std::cout << std::endl;

    /////////////////////
    // Normalise count information
    /////////////////////
    nEvents0Shower0Track = nEvents0Shower0Track * 100.0 / nEventsAboveThreshold;
    nEvents0Shower1Track = nEvents0Shower1Track * 100.0 / nEventsAboveThreshold;
    nEvents0Shower2Track = nEvents0Shower2Track * 100.0 / nEventsAboveThreshold;
    nEvents0Shower3Track = nEvents0Shower3Track * 100.0 / nEventsAboveThreshold;
    nEvents0Shower4PlusTrack = nEvents0Shower4PlusTrack * 100.0 / nEventsAboveThreshold;
    double nEvents0Shower = nEvents0Shower0Track + nEvents0Shower1Track + nEvents0Shower2Track + nEvents0Shower3Track + nEvents0Shower4PlusTrack;

    nEvents1Shower0Track = nEvents1Shower0Track * 100.0 / nEventsAboveThreshold;
    nEvents1Shower1Track = nEvents1Shower1Track * 100.0 / nEventsAboveThreshold;
    nEvents1Shower2Track = nEvents1Shower2Track * 100.0 / nEventsAboveThreshold;
    nEvents1Shower3Track = nEvents1Shower3Track * 100.0 / nEventsAboveThreshold;
    nEvents1Shower4PlusTrack = nEvents1Shower4PlusTrack * 100.0 / nEventsAboveThreshold;
    double nEvents1Shower = nEvents1Shower0Track + nEvents1Shower1Track + nEvents1Shower2Track + nEvents1Shower3Track + nEvents1Shower4PlusTrack;

    nEvents2PlusShower0Track = nEvents2PlusShower0Track * 100.0 / nEventsAboveThreshold;
    nEvents2PlusShower1Track = nEvents2PlusShower1Track * 100.0 / nEventsAboveThreshold;
    nEvents2PlusShower2Track = nEvents2PlusShower2Track * 100.0 / nEventsAboveThreshold;
    nEvents2PlusShower3Track = nEvents2PlusShower3Track * 100.0 / nEventsAboveThreshold;
    nEvents2PlusShower4PlusTrack = nEvents2PlusShower4PlusTrack * 100.0 / nEventsAboveThreshold;
    double nEvents2PlusShower = nEvents2PlusShower0Track + nEvents2PlusShower1Track + nEvents2PlusShower2Track + nEvents2PlusShower3Track + nEvents2PlusShower4PlusTrack;

    double nEvents0Track = nEvents0Shower0Track + nEvents1Shower0Track + nEvents2PlusShower0Track;
    double nEvents1Track = nEvents0Shower1Track + nEvents1Shower1Track + nEvents2PlusShower1Track;
    double nEvents2Track = nEvents0Shower2Track + nEvents1Shower2Track + nEvents2PlusShower2Track;
    double nEvents3Track = nEvents0Shower3Track + nEvents1Shower3Track + nEvents2PlusShower3Track;
    double nEvents4PlusTrack = nEvents0Shower4PlusTrack + nEvents1Shower4PlusTrack + nEvents2PlusShower4PlusTrack;

    for (int particleType : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        thisSliceTrackID[particleType] = thisSliceTrackID[particleType] * 100.0 / foundThisSlice[particleType];
        thisSliceShowerID[particleType] = thisSliceShowerID[particleType] * 100.0 / foundThisSlice[particleType];
        thisSliceGen1[particleType] = thisSliceGen1[particleType] * 100.0 / foundThisSlice[particleType];
        thisSliceGen2[particleType] = thisSliceGen2[particleType] * 100.0 / foundThisSlice[particleType];
        thisSliceGen3[particleType] = thisSliceGen3[particleType] * 100.0 / foundThisSlice[particleType];
        thisSliceGen4Plus[particleType] = thisSliceGen4Plus[particleType] * 100.0 / foundThisSlice[particleType];

        foundThisSlice[particleType] = foundThisSlice[particleType] * 100.0 / nEvents;
        protonPionMerged[particleType] = protonPionMerged[particleType] * 100.0 / nEvents;
        foundOtherSlice[particleType] = foundOtherSlice[particleType] * 100.0 / nEvents;
        noSlice[particleType] = noSlice[particleType] * 100.0 / nEvents;
    }

    allThisSlice = allThisSlice * 100.0 / nEvents;
    allThisSlicePlusGen = allThisSlicePlusGen * 100.0 / nEvents;
    allThisSlicePlusID = allThisSlicePlusID * 100.0 / nEvents;
    allThisSlicePlusGenPlusID = allThisSlicePlusGenPlusID * 100.0 / nEvents;

    for (int particleType : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        thisSliceTrackID_15[particleType] = thisSliceTrackID_15[particleType] * 100.0 / foundThisSlice_15[particleType];
        thisSliceShowerID_15[particleType] = thisSliceShowerID_15[particleType] * 100.0 / foundThisSlice_15[particleType];
        thisSliceGen2_15[particleType] = thisSliceGen2_15[particleType] * 100.0 / foundThisSlice_15[particleType];
        thisSliceGen3_15[particleType] = thisSliceGen3_15[particleType] * 100.0 / foundThisSlice_15[particleType];
        thisSliceGen4Plus_15[particleType] = thisSliceGen4Plus_15[particleType] * 100.0 / foundThisSlice_15[particleType];

        foundThisSlice_15[particleType] = foundThisSlice_15[particleType] * 100.0 / nEventsAboveThreshold;
        protonPionMerged_15[particleType] = protonPionMerged_15[particleType] * 100.0 / nEventsAboveThreshold;
        foundOtherSlice_15[particleType] = foundOtherSlice_15[particleType] * 100.0 / nEventsAboveThreshold;
        noSlice_15[particleType] = noSlice_15[particleType] * 100.0 / nEventsAboveThreshold;
    }

    allThisSlice_15 = allThisSlice_15 * 100.0 / nEventsAboveThreshold;
    allThisSlicePlusGen_15 = allThisSlicePlusGen_15 * 100.0 / nEventsAboveThreshold;
    allThisSlicePlusID_15 = allThisSlicePlusID_15 * 100.0 / nEventsAboveThreshold;
    allThisSlicePlusGenPlusID_15 = allThisSlicePlusGenPlusID_15 * 100.0 / nEventsAboveThreshold;

    std::cout << std::endl;
    std::cout << "Existence information..." << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "| Particle | Nu slice | Other slice |  Merged  |    None    |" << std::endl;
    std::cout << "|  Muon    | " << foundThisSlice[MUON_INDEX] << "% | " << foundOtherSlice[MUON_INDEX] << "% | " <<
        "  N/A   | " << noSlice[MUON_INDEX] << "% |" << std::endl;
    std::cout << "|  Proton  | " << foundThisSlice[PROTON_INDEX] << "% | " << foundOtherSlice[PROTON_INDEX] << "% | " << 
        protonPionMerged[PROTON_INDEX] << "% | " << noSlice[PROTON_INDEX] << "% |" << std::endl;
    std::cout << "|  Pion    | " << foundThisSlice[PION_INDEX] << "% | " << foundOtherSlice[PION_INDEX] << "% | " << 
        protonPionMerged[PION_INDEX] << "% | " << noSlice[PION_INDEX] << "% |" << std::endl;
    std::cout << "|  Gamma   | " << foundThisSlice[GAMMA_INDEX] << "% | " << foundOtherSlice[GAMMA_INDEX] << "% | " << 
        "  N/A   | " << noSlice[GAMMA_INDEX] << "% |" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "Track/shower ID for nu slice reco... " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "| Particle | Track | Shower |" << std::endl;
    std::cout << "|  Muon    | " << thisSliceTrackID[MUON_INDEX] << "% | " << thisSliceShowerID[MUON_INDEX] << "% | " << std::endl;
    std::cout << "|  Proton  | " << thisSliceTrackID[PROTON_INDEX] << "% | " << thisSliceShowerID[PROTON_INDEX] << "% | " << std::endl;
    std::cout << "|  Pion    | " << thisSliceTrackID[PION_INDEX] << "% | " << thisSliceShowerID[PION_INDEX] << "% | " << std::endl;
    std::cout << "|  Gamma   | " << thisSliceTrackID[GAMMA_INDEX] << "% | " << thisSliceShowerID[GAMMA_INDEX] << "% | " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "Track/shower ID for nu slice reco... " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "| Particle |  Gen 2  |  Gen 3  |  Gen 4+  |" << std::endl;
    std::cout << "|  Muon    | " << thisSliceGen2[MUON_INDEX] << "% | " << thisSliceGen3[MUON_INDEX] << "% | " << 
        thisSliceGen4Plus[MUON_INDEX] << "% | " << std::endl;
    std::cout << "|  Proton  | " << thisSliceGen2[PROTON_INDEX] << "% | " << thisSliceGen3[PROTON_INDEX] << "% | " << 
        thisSliceGen4Plus[PROTON_INDEX] << "% | " << std::endl;
    std::cout << "|  Pion    | " << thisSliceGen2[PION_INDEX] << "% | " << thisSliceGen3[PION_INDEX] << "% | " << 
        thisSliceGen4Plus[PION_INDEX] << "% | " << std::endl;
    std::cout << "|  Gamma   | " << thisSliceGen2[GAMMA_INDEX] << "% | " << thisSliceGen3[GAMMA_INDEX] << "% | " << 
        thisSliceGen4Plus[GAMMA_INDEX] << "% | " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "Reco particle matrix for nu slice..." << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "|         | 0 Shower | 1 Shower | 2 Shower |           | " << std::endl;
    std::cout << "| 0 Track | " << nEvents0Shower0Track << "% | " << nEvents1Shower0Track << "% |" << nEvents2PlusShower0Track << "% |" << nEvents0Track << "% | " << std::endl; 
    std::cout << "| 1 Track | " << nEvents0Shower1Track << "% | " << nEvents1Shower1Track << "% |" << nEvents2PlusShower1Track << "% |" << nEvents1Track << "% | " << std::endl; 
    std::cout << "| 2 Track | " << nEvents0Shower2Track << "% | " << nEvents1Shower2Track << "% |" << nEvents2PlusShower2Track << "% |" << nEvents2Track << "% | " << std::endl; 
    std::cout << "| 3 Track | " << nEvents0Shower3Track << "% | " << nEvents1Shower3Track << "% |" << nEvents2PlusShower3Track << "% |" << nEvents3Track << "% | " << std::endl; 
    std::cout << "| 4 Track | " << nEvents0Shower4PlusTrack << "% | " << nEvents1Shower4PlusTrack << "% |" << nEvents2PlusShower4PlusTrack << "% | " << nEvents4PlusTrack << "% |" << std::endl; 
    std::cout << "|         | " << nEvents0Shower << "% | " << nEvents1Shower << "% | " << nEvents2PlusShower << "% |         | " << std::endl; 
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "How did we do?" << std::endl;
    std::cout << allThisSlice << "% of events had all particles found in nu slice, where" << std::endl;
    std::cout << allThisSlicePlusGen << "% also had correct hierarchy tiers" << std::endl;   
    std::cout << allThisSlicePlusID << "% also had correct track/shower IDs" << std::endl;   
    std::cout << allThisSlicePlusGenPlusID << "% also had correct hierachy tier & track/shower IDs" << std::endl;   

    std::cout << std::endl;
    std::cout << "Okay, but how do we do when above reco thresholds (each particle has >= 15 hits" << std::endl;
    std::cout << nEventsAboveThreshold * 100.0 / nEvents << "% of events above threshold" << std::endl;
    std::cout << "Existence information..." << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "| Particle | Nu slice | Other slice |  Merged  |    None    |" << std::endl;
    std::cout << "|  Muon    | " << foundThisSlice_15[MUON_INDEX] << "% | " << foundOtherSlice_15[MUON_INDEX] << "% | " <<
        "  N/A   | " << noSlice_15[MUON_INDEX] << "% |" << std::endl;
    std::cout << "|  Proton  | " << foundThisSlice_15[PROTON_INDEX] << "% | " << foundOtherSlice_15[PROTON_INDEX] << "% | " << 
        protonPionMerged_15[PROTON_INDEX] << "% | " << noSlice_15[PROTON_INDEX] << "% |" << std::endl;
    std::cout << "|  Pion    | " << foundThisSlice_15[PION_INDEX] << "% | " << foundOtherSlice_15[PION_INDEX] << "% | " << 
        protonPionMerged_15[PION_INDEX] << "% | " << noSlice_15[PION_INDEX] << "% |" << std::endl;
    std::cout << "|  Gamma   | " << foundThisSlice_15[GAMMA_INDEX] << "% | " << foundOtherSlice_15[GAMMA_INDEX] << "% | " << 
        "  N/A   | " << noSlice_15[GAMMA_INDEX] << "% |" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "Track/shower ID for nu slice reco... " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "| Particle | Track | Shower |" << std::endl;
    std::cout << "|  Muon    | " << thisSliceTrackID_15[MUON_INDEX] << "% | " << thisSliceShowerID_15[MUON_INDEX] << "% | " << std::endl;
    std::cout << "|  Proton  | " << thisSliceTrackID_15[PROTON_INDEX] << "% | " << thisSliceShowerID_15[PROTON_INDEX] << "% | " << std::endl;
    std::cout << "|  Pion    | " << thisSliceTrackID_15[PION_INDEX] << "% | " << thisSliceShowerID_15[PION_INDEX] << "% | " << std::endl;
    std::cout << "|  Gamma   | " << thisSliceTrackID_15[GAMMA_INDEX] << "% | " << thisSliceShowerID_15[GAMMA_INDEX] << "% | " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "Track/shower ID for nu slice reco... " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "| Particle |  Gen 2  |  Gen 3  |  Gen 4+  |" << std::endl;
    std::cout << "|  Muon    | " << thisSliceGen2_15[MUON_INDEX] << "% | " << thisSliceGen3_15[MUON_INDEX] << "% | " << 
        thisSliceGen4Plus_15[MUON_INDEX] << "% | " << std::endl;
    std::cout << "|  Proton  | " << thisSliceGen2_15[PROTON_INDEX] << "% | " << thisSliceGen3_15[PROTON_INDEX] << "% | " << 
        thisSliceGen4Plus_15[PROTON_INDEX] << "% | " << std::endl;
    std::cout << "|  Pion    | " << thisSliceGen2[PION_INDEX] << "% | " << thisSliceGen3[PION_INDEX] << "% | " << 
        thisSliceGen4Plus_15[PION_INDEX] << "% | " << std::endl;
    std::cout << "|  Gamma   | " << thisSliceGen2_15[GAMMA_INDEX] << "% | " << thisSliceGen3_15[GAMMA_INDEX] << "% | " << 
        thisSliceGen4Plus_15[GAMMA_INDEX] << "% | " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "How did we do?" << std::endl;
    std::cout << allThisSlice_15 << "% of events had all particles found in nu slice, where" << std::endl;
    std::cout << allThisSlicePlusGen_15 << "% also had correct hierarchy tiers" << std::endl;   
    std::cout << allThisSlicePlusID_15 << "% also had correct track/shower IDs" << std::endl;   
    std::cout << allThisSlicePlusGenPlusID_15 << "% also had correct hierachy tier & track/shower IDs" << std::endl;

    /////////////////////
    // Divide efficiency histograms
    /////////////////////
    
    // True hits
    ObtainEfficiencyHistogramCollection(nTrueHitsThisSliceEff, nTrueHitsDists);
    ObtainEfficiencyHistogramCollection(nTrueHitsOtherSliceEff, nTrueHitsDists);

    // True distance
    //ObtainEfficiencyHistogramCollection(nTrueDistanceThisSliceEff, nTrueDistanceDists);
    //ObtainEfficiencyHistogramCollection(nTrueDistanceOtherSliceEff, nTrueDistanceDists);    
    

    /////////////////////
    // Normalise histograms
    /////////////////////
    //NormaliseHistogramCollection(nTrueHitsDists);
    //NormaliseHistogramCollection(nTrueDistanceDists);
    //NormaliseHistogramCollection(trackScoreDists);
    //NormaliseHistogramCollection(completenessDists);
    NormaliseHistogramCollection(purityDists);

    /////////////////////
    // Draw histograms
    /////////////////////
    TCanvas * c1 = new TCanvas("c1", "c1");

    // Distributions
    //DrawHistogramCollection(nTrueHitsDists, "", "nTrueHits", "Proportion of Events");
    //DrawHistogramCollection(nTrueDistanceDists, "", "TrueNuVertexSep", "Proportion of Events");
    //DrawHistogramCollection(trackScoreDists, "", "Track Score", "Proportion of Events");
    //DrawHistogramCollection(completenessDists, "", "Completeness", "Proportion of Events");
    DrawHistogramCollection(purityDists, "", "Purity", "Proportion of Events");

    // Each particle type efficiency
    //TCanvas * cSliceEfficiency_MUON = new TCanvas("cSliceEfficiency_MUON", "cSliceEfficiency_MUON");
    //DrawSliceEfficiency(nTrueHitsDists, nTrueHitsThisSliceEff, nTrueHitsOtherSliceEff, MUON_INDEX);
    //DrawSliceEfficiency(nTrueDistanceDists, nTrueDistanceThisSliceEff, nTrueDistanceOtherSliceEff, MUON_INDEX);

    //TCanvas * cSliceEfficiency_PROTON = new TCanvas("cSliceEfficiency_PROTON", "cSliceEfficiency_PROTON");
    //DrawSliceEfficiency(nTrueHitsDists, nTrueHitsThisSliceEff, nTrueHitsOtherSliceEff, PROTON_INDEX);
    //DrawSliceEfficiency(nTrueDistanceDists, nTrueDistanceThisSliceEff, nTrueDistanceOtherSliceEff, PROTON_INDEX);

    //TCanvas * cSliceEfficiency_PION = new TCanvas("cSliceEfficiency_PION", "cSliceEfficiency_PION");
    //DrawSliceEfficiency(nTrueHitsDists, nTrueHitsThisSliceEff, nTrueHitsOtherSliceEff, PION_INDEX);
    //DrawSliceEfficiency(nTrueDistanceDists, nTrueDistanceThisSliceEff, nTrueDistanceOtherSliceEff, PION_INDEX);

    //TCanvas * cSliceEfficiency_GAMMA = new TCanvas("cSliceEfficiency_GAMMA", "cSliceEfficiency_GAMMA");
    //DrawSliceEfficiency(nTrueHitsDists, nTrueHitsThisSliceEff, nTrueHitsOtherSliceEff, GAMMA_INDEX);
    //DrawSliceEfficiency(nTrueDistanceDists, nTrueDistanceThisSliceEff, nTrueDistanceOtherSliceEff, GAMMA_INDEX);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void ObtainEfficiencyHistogramCollection(const std::vector<TH1D*> &weightedHistogramCollection, 
    const std::vector<TH1D*> &fullHistogramCollection)
{

    for (unsigned int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        for (unsigned int bin = 1; bin <= weightedHistogramCollection[particleTypeIndex]->GetNbinsX(); ++bin)
        {
            const double weighted(weightedHistogramCollection[particleTypeIndex]->GetBinContent(bin));
            const double full(fullHistogramCollection[particleTypeIndex]->GetBinContent(bin));
            const double newBinContent(full < std::numeric_limits<double>::epsilon() ? 0.0 : weighted / full);

            weightedHistogramCollection[particleTypeIndex]->SetBinContent(bin, newBinContent);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void NormaliseHistogramCollection(const std::vector<TH1D*> &histogramCollection)
{
    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        if (histogramCollection[particleTypeIndex]->Integral() > std::numeric_limits<float>::epsilon())
            histogramCollection[particleTypeIndex]->Scale(1.f / histogramCollection[particleTypeIndex]->Integral());
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void DrawHistogramCollection(const std::vector<TH1D*> &histogramCollection, const std::string &title,
    const std::string &xAxis, const std::string &yAxis)
{
    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
    {
        std::string jam = title + ";" + xAxis + ";" + yAxis;
        histogramCollection[particleTypeIndex]->SetTitle(jam.c_str());
        histogramCollection[particleTypeIndex]->GetYaxis()->SetRangeUser(0.0, 1.0);
        histogramCollection[particleTypeIndex]->SetLineColor(LINE_COLOR[particleTypeIndex]);
        histogramCollection[particleTypeIndex]->SetLineWidth(2);
    }

    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        histogramCollection[particleTypeIndex]->Draw("hist same");

    TLegend * legend = new TLegend(0.1,0.7,0.48,0.9);

    for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        legend->AddEntry(histogramCollection[particleTypeIndex], LEGEND_STRING[particleTypeIndex], "l");

    legend->Draw("same");

    //delete legend;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void DrawSliceEfficiency(const std::vector<TH1D*> &fullHistogramCollection, const std::vector<TH1D*> &thisSliceEffCollection,
    const std::vector<TH1D*> &otherSliceEffCollection, const int particleTypeIndex)
{
    TH1D * totalEff = (TH1D*)thisSliceEffCollection[particleTypeIndex]->Clone();
    totalEff->Add(otherSliceEffCollection[particleTypeIndex]);

    thisSliceEffCollection[particleTypeIndex]->SetLineStyle(9);
    otherSliceEffCollection[particleTypeIndex]->SetLineStyle(2);
    totalEff->SetLineStyle(1);

    for (TH1D * hist : {thisSliceEffCollection[particleTypeIndex], otherSliceEffCollection[particleTypeIndex], totalEff})
    {
        hist->GetYaxis()->SetRangeUser(0.0, 1.0);
        hist->SetTitle(";nTrueHits;Efficiency");
        hist->SetLineWidth(2);
        hist->SetLineColor(LINE_COLOR[particleTypeIndex]);
    }

    //fullHistogramCollection[particleTypeIndex]->SetFillStyle(1001);
    fullHistogramCollection[particleTypeIndex]->SetFillColorAlpha(LINE_COLOR[particleTypeIndex], 0.2);
    fullHistogramCollection[particleTypeIndex]->GetYaxis()->SetRangeUser(0.0, 1.0);
    fullHistogramCollection[particleTypeIndex]->SetTitle(";nTrueHits;Efficiency");

    if (fullHistogramCollection[particleTypeIndex]->Integral() > std::numeric_limits<float>::epsilon())
        fullHistogramCollection[particleTypeIndex]->Scale(1.f / fullHistogramCollection[particleTypeIndex]->Integral());


    for (TH1D * hist : {thisSliceEffCollection[particleTypeIndex], otherSliceEffCollection[particleTypeIndex], fullHistogramCollection[particleTypeIndex]})
        hist->Draw("hist same");

    TLegend * legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(thisSliceEffCollection[particleTypeIndex], "Neutrino Slice", "l");
    legend->AddEntry(otherSliceEffCollection[particleTypeIndex], "Other Slice", "l");
    //legend->AddEntry(totalEff, "Total", "l");
    legend->AddEntry(fullHistogramCollection[particleTypeIndex], "Normalised True Distribution", "f");
    legend->Draw("same");

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
