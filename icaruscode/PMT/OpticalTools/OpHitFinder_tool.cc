////////////////////////////////////////////////////////////////////////
/// \file   OpHitFinder.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "icaruscode/PMT/OpticalTools/IOpHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"

#include <cmath>
#include <fstream>

namespace light
{

class OpHitFinder : public IOpHitFinder
{
public:
    explicit OpHitFinder(const fhicl::ParameterSet& pset);
    
    ~OpHitFinder();
    
    void configure(const fhicl::ParameterSet& pset)                 override;
    void outputHistograms(art::TFileDirectory&)               const override;
    
    void FindOpHits(const raw::OpDetWaveform&, OpHitVec&)     const override;
    
private:
    // fhicl parameters
    float          fSPEArea;         //< conversion between phe and Adc*ns
    float          fSaturationCut;   //< Saturation occurs at this point
    int            fMaxSatChannels;  //< maximum saturated channels
    mutable size_t fEventCount;      //< Keep track of the number of events processed

    float getBaseline(const raw::OpDetWaveform&) const;
    
    std::unique_ptr<reco_tool::ICandidateHitFinder> fHitFinderTool;  ///< For finding candidate hits
};
    
//----------------------------------------------------------------------
// Constructor.
OpHitFinder::OpHitFinder(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
OpHitFinder::~OpHitFinder()
{
}
    
void OpHitFinder::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fSPEArea        = pset.get< float >("SPEArea"              );
    fSaturationCut  = pset.get< float >("SaturationCut",  3000.);
    fMaxSatChannels = pset.get< int   >("MaxSatChannels",    10);
    
    fEventCount     = 0;

    fHitFinderTool  = art::make_tool<reco_tool::ICandidateHitFinder>(pset.get<fhicl::ParameterSet>("CandidateHits"));

    return;
}

    
void OpHitFinder::FindOpHits(const raw::OpDetWaveform& opDetWaveform,
                             OpHitVec&                 opHitVec) const
{
    // The plan here:
    // 1) Get the baseline
    // 2) Copy to a local vector doing baseline subtraction and inversion
    // 3) Set up and call the standard gaushit finder tools for finding peaks
    // 4) Return the parameters for an ophit
    float baseline = getBaseline(opDetWaveform);

    std::vector<float> locWaveform;
    
    locWaveform.resize(opDetWaveform.size());
    
    // The aim here is to baseline correct AND invert the waveform
    std::transform(opDetWaveform.begin(),opDetWaveform.end(),locWaveform.begin(),[baseline](const auto& val){return baseline - val;});
    
    std::pair<std::vector<float>::iterator,std::vector<float>::iterator> minMaxItr = std::minmax_element(locWaveform.begin(),locWaveform.end());
    
    reco_tool::ICandidateHitFinder::HitCandidateVec      hitCandidateVec;
    reco_tool::ICandidateHitFinder::MergeHitCandidateVec mergedCandidateHitVec;
    
    bool  notSaturated(true);
    int   numSatChans(0);
    float maxAdc = *minMaxItr.second;
    float minAdc = *minMaxItr.first;
    
    if (maxAdc - minAdc > fSaturationCut)
    {
        float unSatAdc = 0.5 * (maxAdc - minAdc) + minAdc;
        std::vector<float>::iterator unSatItr = std::find_if(minMaxItr.second,locWaveform.end(),[unSatAdc](const auto& elem){return elem < unSatAdc;});
        
        numSatChans = std::distance(minMaxItr.second,unSatItr);
        
        if (numSatChans > fMaxSatChannels) notSaturated = false;
    }
        
    if (notSaturated)
    {
        std::vector<float> tempVec = locWaveform;
        recob::Wire::RegionsOfInterest_t::datarange_t rangeData(size_t(0),std::move(tempVec));
        
        fHitFinderTool->findHitCandidates(rangeData, 0, 0, fEventCount, hitCandidateVec);
        fHitFinderTool->MergeHitCandidates(rangeData, hitCandidateVec, mergedCandidateHitVec);
    }
    else
    {
        reco_tool::ICandidateHitFinder::HitCandidate hitCandidate;
        
        hitCandidate.startTick     = 0;
        hitCandidate.stopTick      = 0;
        hitCandidate.maxTick       = 0;
        hitCandidate.minTick       = 0;
        hitCandidate.maxDerivative = 0;
        hitCandidate.minDerivative = 0;
        hitCandidate.hitCenter     = std::distance(locWaveform.begin(),minMaxItr.second);
        hitCandidate.hitSigma      = std::min(40,numSatChans);
        hitCandidate.hitHeight     = maxAdc - minAdc;
        
        std::cout << "***>> Saturated PMT, numSatChans: " << numSatChans << ", range: " << maxAdc-minAdc << std::endl;
        
        hitCandidateVec.push_back(hitCandidate);
        mergedCandidateHitVec.push_back(hitCandidateVec);
    }
    
    // Recover the channel number
    raw::Channel_t chNumber = opDetWaveform.ChannelNumber();

    // Go through the hit candidates and convert to ophits
    // Note that the "merged candidates" represent lists of candidate hits that
    // are in one pulse train... so we need a double loop
    for(const auto& mergedCands : mergedCandidateHitVec)
    {
        for(const auto& candidateHit : mergedCands)
        {
            float peakMean    = candidateHit.hitCenter;
            float peakSigma   = candidateHit.hitSigma;
            float amplitude   = candidateHit.hitHeight;
            float peakArea    = amplitude * peakSigma * 2.5066;  // * sqrt(2pi)
            float nPhotoElec  = peakArea / fSPEArea;
            
            float peakTimeAbs = peakMean;   // NOTE: these times will need to be corrected
            int   frame       = 1;          //       also this needs to be the clock frame
            float fastTotal   = 0.;         //       not sure what this is
            
            opHitVec.emplace_back(chNumber, peakMean, peakTimeAbs, frame, 2.35 * peakSigma, peakArea, amplitude, nPhotoElec, fastTotal);//including hit info
        }
    }
    
    fEventCount++;

    return;
}

float OpHitFinder::getBaseline(const raw::OpDetWaveform& locWaveform) const
{
    // Fill a map to determine the most probable value
    std::map<raw::ADC_Count_t,int> adcFrequencyMap;
    
    raw::ADC_Count_t maxBin(0);
    int              maxCount(0);
    
    for(const auto& adc : locWaveform)
    {
        int& adcFrequency = adcFrequencyMap[adc];
        
        if (++adcFrequency > maxCount)
        {
            maxBin   = adc;
            maxCount = adcFrequency;
        }
    }
    
    float mostProbableBaseline(0.);
    int   mostProbableCount(0);
    
    for(raw::ADC_Count_t adcBin = maxBin - 3; adcBin <= maxBin + 3; adcBin++)
    {
        try{
            mostProbableBaseline += adcFrequencyMap.at(adcBin) * float(adcBin);
            mostProbableCount    += adcFrequencyMap.at(adcBin);
        }
        catch(...) {}
    }
    
    mostProbableBaseline /= mostProbableCount;
   
    return mostProbableBaseline;
}

    
void OpHitFinder::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "OpHitFinderPlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fOpHitFinderVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "OpHitFinderPlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "OpHitFinder;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fOpHitFinderVec.at(bin).Re());
    }
*/
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(OpHitFinder)
}
