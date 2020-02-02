////////////////////////////////////////////////////////////////////////
//
// Class:       FilterNoiseICARUS
// Module Type: producer
// File:        FilterNoiseICARUS_module.cc
//
//              The intent of this module is to both "decode" artdaq fragments
//              and convert to RawDigits and also to do the initial noise 
//              filtering, specifically the coherent noise. 
//
// Configuration parameters:
//
// DigitModuleLabel      - the source of the RawDigit collection
//
//
// Modeled after example from Mike Wang (mwang@fnal.gov)
// Copied/Modified by Tracy Usher (usher@slac.stanford.edu) on January 27, 2020
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/ReplicatedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RawData/RawDigit.h"

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/Decode/DecoderTools/IDecoderFilter.h"

#include "icarussigproc/ICARUSSigProcDefs.h"

namespace daq 
{

class FilterNoiseICARUS : public art::ReplicatedProducer
{
public:

    // Copnstructors, destructor.
    explicit FilterNoiseICARUS(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame);
    virtual ~FilterNoiseICARUS();

    // Overrides.
    virtual void configure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e, art::ProcessingFrame const& frame);
    virtual void beginJob(art::ProcessingFrame const& frame);
    virtual void endJob(art::ProcessingFrame const& frame);

private:
    // Define the RawDigit collection
    using RawDigitCollection    = std::vector<raw::RawDigit>;
    using RawDigitCollectionPtr = std::unique_ptr<RawDigitCollection>;

    // Function to save our RawDigits
    void saveRawDigits(const icarussigproc::ArrayFloat&, 
                       const icarussigproc::VectorFloat&, 
                       RawDigitCollectionPtr&, size_t);

    // Tools for decoding fragments depending on type
    std::vector<std::unique_ptr<IDecoderFilter>> fDecoderToolVec;      ///< Decoder tools

    // Fcl parameters.
    art::InputTag                                fFragmentsLabel;      ///< The input artdaq fragment label
    bool                                         fOutputPedestalCor;   ///< Should we output pedestal corrected (not noise filtered)?
    std::string                                  fOutputPedCorPath;    ///< Path to assign to the output if asked for

    // Statistics.
    int                                          fNumEvent;             ///< Number of events seen.

    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const*                     fGeometry;             ///< pointer to Geometry service
    detinfo::DetectorProperties const*           fDetectorProperties;   ///< Detector properties service
};

DEFINE_ART_MODULE(FilterNoiseICARUS)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
FilterNoiseICARUS::FilterNoiseICARUS(fhicl::ParameterSet const & pset, art::ProcessingFrame const& frame) :
                      art::ReplicatedProducer(pset, frame),
                      fNumEvent(0)
{

    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    configure(pset);

    // Recover the vector of fhicl parameters for the ROI tools
    const fhicl::ParameterSet& decoderTools = pset.get<fhicl::ParameterSet>("DecoderToolVec");
    
    fDecoderToolVec.resize(decoderTools.get_pset_names().size());
    
    for(const std::string& decoderTool : decoderTools.get_pset_names())
    {
        const fhicl::ParameterSet& decoderToolParamSet = decoderTools.get<fhicl::ParameterSet>(decoderTool);
        
        // Get instance of tool
        fDecoderToolVec.push_back(art::make_tool<IDecoderFilter>(decoderToolParamSet));
    }

    produces<std::vector<raw::RawDigit>>();

    if (fOutputPedestalCor)
        produces<std::vector<raw::RawDigit>>(fOutputPedCorPath);

    // Report.
    mf::LogInfo("FilterNoiseICARUS") << "FilterNoiseICARUS configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
FilterNoiseICARUS::~FilterNoiseICARUS()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void FilterNoiseICARUS::configure(fhicl::ParameterSet const & pset)
{
    fFragmentsLabel    = pset.get<art::InputTag>("FragmentsLabel",    "daq:PHYSCRATEDATA");
    fOutputPedestalCor = pset.get<bool         >("OutputPedestalCor",               false);
    fOutputPedCorPath  = pset.get<std::string  >("OutputPedCorPath",                "RAW");
}

//----------------------------------------------------------------------------
/// Begin job method.
void FilterNoiseICARUS::beginJob(art::ProcessingFrame const&)
{ 
    return;
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method.
///
void FilterNoiseICARUS::produce(art::Event & event, art::ProcessingFrame const&)
{
    ++fNumEvent;

    // Allocate the output collection of RawDigits
    RawDigitCollectionPtr rawDigitCollection(new std::vector<raw::RawDigit>);

    // Are we outputting the raw (pedestal corrected only) RawDigits? 
    RawDigitCollectionPtr rawRawDigitCollection(new std::vector<raw::RawDigit>);

    art::Handle<artdaq::Fragments> daq_handle;
    event.getByLabel(fFragmentsLabel, daq_handle);

    std::cout << "**** Processing raw data fragments ****" << std::endl;

    // Get the right tool for the job
    daq::IDecoderFilter* decoderTool = fDecoderToolVec.back().get();

    // Loop over daq fragments
    for (auto const &fragment: *daq_handle) 
    {
        std::cout << "--> Processing fragment ID: " << fragment.fragmentID() << std::endl;

        //process_fragment(event, rawfrag, product_collection, header_collection);
        decoderTool->process_fragment(fragment);

        // Useful numerology
        // convert fragment to Nevis fragment
        icarus::PhysCrateFragment physCrateFragment(fragment);

        size_t nBoardsPerFragment = physCrateFragment.nBoards();
        size_t nChannelsPerBoard  = physCrateFragment.nChannelsPerBoard();

        // Set base channel for both the board and the board/fragment
        size_t boardFragOffset    = nChannelsPerBoard * nBoardsPerFragment * fragment.fragmentID();

        // Save the filtered RawDigits
        saveRawDigits(decoderTool->getWaveLessCoherent(),decoderTool->getTruncRMSVals(),rawDigitCollection,boardFragOffset);

        // Optionally, save the pedestal corrected RawDigits
        if (fOutputPedestalCor)
            saveRawDigits(decoderTool->getPedSubtractedWaveforms(),decoderTool->getFullRMSVals(),rawRawDigitCollection,boardFragOffset);
    }

    // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
    std::sort(rawDigitCollection->begin(),rawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});

    // Now transfer ownership to the event store
    event.put(std::move(rawDigitCollection));

    if (fOutputPedestalCor)
    {
        // Want the RawDigits to be sorted in channel order... has to be done somewhere so why not now?
        std::sort(rawRawDigitCollection->begin(),rawRawDigitCollection->end(),[](const auto& left,const auto&right){return left.Channel() < right.Channel();});

        // Now transfer ownership to the event store
        event.put(std::move(rawRawDigitCollection),fOutputPedCorPath);
    }

    return;
}

void FilterNoiseICARUS::saveRawDigits(const icarussigproc::ArrayFloat&  dataArray, 
                                      const icarussigproc::VectorFloat& rmsVec,
                                      RawDigitCollectionPtr&                  rawDigitCol, 
                                      size_t                                  channel)
{
        raw::RawDigit::ADCvector_t wvfm(dataArray[0].size());

        // Loop over the channels to recover the RawDigits after filtering
        for(size_t chanIdx = 0; chanIdx != dataArray.size(); chanIdx++)
        {
            const icarussigproc::VectorFloat& dataVec = dataArray[chanIdx];

            // Need to convert from float to short int
            std::transform(dataVec.begin(),dataVec.end(),wvfm.begin(),[](const auto& val){return short(std::round(val));});

            rawDigitCol->emplace_back(channel++,wvfm.size(),wvfm); 
            rawDigitCol->back().SetPedestal(0.,rmsVec[chanIdx]);
        }//loop over channel indices

    return;
}

//----------------------------------------------------------------------------
/// End job method.
void FilterNoiseICARUS::endJob(art::ProcessingFrame const&)
{
    mf::LogInfo("FilterNoiseICARUS") << "Looked at " << fNumEvent << " events" << std::endl;
}

} // end of namespace