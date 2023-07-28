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

    tree->SetBranchAddress("NTrueHits", &t_nTrueHits);
    tree->SetBranchAddress("MatchFoundInSlice", &t_matchFoundInSlice);
    tree->SetBranchAddress("MatchFoundInOtherSlice", &t_matchFoundInOtherSlice);
    tree->SetBranchAddress("TrueNuVertexSep", &t_trueNuVertexSep);
    tree->SetBranchAddress("BestMatchTrackScore", &t_trackScore);

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

    // Track Score
    double trackScoreMin(-1.05), trackScoreMax(1.05), trackScoreBin(22);
    TH1D * trackScoreDist_0 = new TH1D("trackScoreDist_0", "trackScoreDist_0", trackScoreBin, trackScoreMin, trackScoreMax);
    TH1D * trackScoreDist_1 = new TH1D("trackScoreDist_1", "trackScoreDist_1", trackScoreBin, trackScoreMin, trackScoreMax);
    TH1D * trackScoreDist_2 = new TH1D("trackScoreDist_2", "trackScoreDist_2", trackScoreBin, trackScoreMin, trackScoreMax);
    TH1D * trackScoreDist_3 = new TH1D("trackScoreDist_3", "trackScoreDist_3", trackScoreBin, trackScoreMin, trackScoreMax);
    std::vector<TH1D*> trackScoreDists = {trackScoreDist_0, trackScoreDist_1, trackScoreDist_2, trackScoreDist_3};

    /////////////////////
    // Fill histograms
    /////////////////////
    for (unsigned int treeIndex = 0; treeIndex < tree->GetEntries(); ++treeIndex)
    {
        tree->GetEntry(treeIndex);

        for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        {
            const double boundedHits = std::min(std::max(nTrueHitsMin + EPSILON, static_cast<double>(t_nTrueHits->at(particleTypeIndex))), nTrueHitsMax - EPSILON);

            nTrueHitsDists[particleTypeIndex]->Fill(boundedHits);
            nTrueHitsThisSliceEff[particleTypeIndex]->Fill(boundedHits, t_matchFoundInSlice->at(particleTypeIndex));
            nTrueHitsOtherSliceEff[particleTypeIndex]->Fill(boundedHits, t_matchFoundInOtherSlice->at(particleTypeIndex));
        }

        for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        {
            const double boundedDistance = std::min(std::max(nTrueDistanceMin + EPSILON, static_cast<double>(t_trueNuVertexSep->at(particleTypeIndex))), nTrueDistanceMax - EPSILON);

            nTrueDistanceDists[particleTypeIndex]->Fill(boundedDistance);
            nTrueDistanceThisSliceEff[particleTypeIndex]->Fill(boundedDistance, t_matchFoundInSlice->at(particleTypeIndex));
            nTrueDistanceOtherSliceEff[particleTypeIndex]->Fill(boundedDistance, t_matchFoundInOtherSlice->at(particleTypeIndex));
        }

        for (int particleTypeIndex : {MUON_INDEX, PROTON_INDEX, PION_INDEX, GAMMA_INDEX})
        {
            const double boundedTrackScore = std::min(std::max(trackScoreMin + EPSILON, static_cast<double>(t_trackScore->at(particleTypeIndex))), trackScoreMax - EPSILON);

            if (t_matchFoundInSlice->at(particleTypeIndex))
                trackScoreDists[particleTypeIndex]->Fill(boundedTrackScore);
        }
    }

    /////////////////////
    // Divide efficiency histograms
    /////////////////////
    /*
    // True hits
    ObtainEfficiencyHistogramCollection(nTrueHitsThisSliceEff, nTrueHitsDists);
    ObtainEfficiencyHistogramCollection(nTrueHitsOtherSliceEff, nTrueHitsDists);

    // True distance
    ObtainEfficiencyHistogramCollection(nTrueDistanceThisSliceEff, nTrueDistanceDists);
    ObtainEfficiencyHistogramCollection(nTrueDistanceOtherSliceEff, nTrueDistanceDists);    
    */

    /////////////////////
    // Normalise histograms
    /////////////////////
    NormaliseHistogramCollection(nTrueHitsDists);
    NormaliseHistogramCollection(nTrueDistanceDists);
    NormaliseHistogramCollection(trackScoreDists);

    /////////////////////
    // Draw histograms
    /////////////////////
    TCanvas * c1 = new TCanvas("c1", "c1");

    // Distributions
    //DrawHistogramCollection(nTrueHitsDists, "", "nTrueHits", "Prop. of events");
    //DrawHistogramCollection(nTrueDistanceDists, "", "TrueNuVertexSep", "Prop. of events");
    DrawHistogramCollection(trackScoreDists, "IN NU SLICE", "Track Score", "Prop. of events");

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

    for (TH1D * hist : {thisSliceEffCollection[particleTypeIndex], otherSliceEffCollection[particleTypeIndex], totalEff, fullHistogramCollection[particleTypeIndex]})
    {
        hist->GetYaxis()->SetRangeUser(0.0, 1.0);
        hist->SetTitle(";nTrueHits;Efficiency");
        hist->SetLineWidth(2);
        hist->SetLineColor(LINE_COLOR[particleTypeIndex]);
    }

    fullHistogramCollection[particleTypeIndex]->SetFillStyle(1001);
    fullHistogramCollection[particleTypeIndex]->SetFillColor(LINE_COLOR[particleTypeIndex]);

    if (fullHistogramCollection[particleTypeIndex]->Integral() > std::numeric_limits<float>::epsilon())
        fullHistogramCollection[particleTypeIndex]->Scale(1.f / fullHistogramCollection[particleTypeIndex]->Integral());


    for (TH1D * hist : {thisSliceEffCollection[particleTypeIndex], otherSliceEffCollection[particleTypeIndex], fullHistogramCollection[particleTypeIndex]})
        hist->Draw("hist same");

    TLegend * legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(thisSliceEffCollection[particleTypeIndex], "Neutrino Slice", "l");
    legend->AddEntry(otherSliceEffCollection[particleTypeIndex], "Other Slice", "l");
    legend->AddEntry(totalEff, "Total", "l");
    legend->AddEntry(fullHistogramCollection[particleTypeIndex], "Distribution", "f");
    legend->Draw("same");

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
