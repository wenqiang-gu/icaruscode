/**
 *  @file   FakeParticle_tool.cc
 *
 *  @brief  This tool creates a fake particle and overlays on input data
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ICARUS package includes
#include "icaruscode/Utilities/SignalShapingServiceICARUS.h"
#include "icaruscode/Decode/DecoderTools/IFakeParticle.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 *  @brief  FakeParticle class definiton
 */
class FakeParticle : virtual public IFakeParticle
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit FakeParticle(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~FakeParticle();

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) override;

    /**
     *  @brief Creates a fake particle and overlays on the input fragment
     *
     *  @param waveforms  The waveform container to place fake particle on
     */
    virtual void overlayFakeParticle(ArrayFloat& waveforms) override; 

private:

    // fhicl variables
    std::vector<size_t>                fWireRange;             //< The range of wires our particle extends over
    std::vector<size_t>                fTickRange;             //< The range of ticks our particle extends over
    int                                fNumElectronsPerMM;     //< The number of electrons/mm to deposit
    size_t                             fPlaneToSimulate;       //< The plane to simulate

    // Some useful variables
    float                              fMMPerTick;             //< Convert ticks in us to mm
    float                              fMMPerWire;             //< Convert wire pitch to mm
    float                              fTanThetaTW;            //< tangent of angle in time/wire space
    float                              fTanTheta;              //< tangent in euclidean space
    float                              fSinTheta;              //< sine of this angle
    float                              fCosTheta;              //< save some time by also doing cosine

    const geo::Geometry*               fGeometry;              //< pointer to the Geometry service
    const detinfo::DetectorProperties* fDetector;              //< Pointer to the detector properties
    util::SignalShapingServiceICARUS*  fSignalShapeService;    //< Access to the response functions
};

FakeParticle::FakeParticle(fhicl::ParameterSet const &pset)
{
    this->configure(pset);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

FakeParticle::~FakeParticle()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void FakeParticle::configure(fhicl::ParameterSet const &pset)
{
    fWireRange         = pset.get<std::vector<size_t>>("WireRange",             std::vector<size_t>() = {0,575});
    fTickRange         = pset.get<std::vector<size_t>>("TickRange",         std::vector<size_t>() = {1000,4000});
    fNumElectronsPerMM = pset.get<int                >("NumElectronsPerMM",                                6000);
    fPlaneToSimulate   = pset.get<size_t             >("PlaneToSimulate",                                     2);

    fGeometry           = art::ServiceHandle<geo::Geometry const>{}.get();
    fDetector           = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fSignalShapeService = art::ServiceHandle<util::SignalShapingServiceICARUS>{}.get();

    // Convert ticks in us to mm by taking the drift velocity and multiplying by the tick period
    double driftVelocity = fDetector->DriftVelocity() * 10.;   // should be mm/us
    double samplingRate  = fDetector->SamplingRate() / 1000.;  // sampling rate returned in ns

    fMMPerTick = driftVelocity * samplingRate;

    fMMPerWire = fGeometry->WirePitch() * 10.;  // wire pitch returned in cm, want mm

    // Get the sin/tan of the angle in wire space
    fTanThetaTW = float(fTickRange[1] - fTickRange[0]) / float(fWireRange[1] - fWireRange[0]);
    fTanTheta   = fTanThetaTW * fMMPerTick / fMMPerWire;
    fCosTheta   = 1. / sqrt(1. + fTanTheta * fTanTheta);
    fSinTheta   = fTanTheta * fCosTheta;

    return;
}

void FakeParticle::overlayFakeParticle(ArrayFloat& waveforms)
{
    // We have assumed the input waveform array will have 576 wires... 
    // Our "range" must be contained within that. By assumption we start at wire 0, so really just need
    // to set the max range 
    size_t maxWire = std::min(waveforms.size(),fWireRange[1]);

    // Create a temporary waveform to handle the input charge
    std::vector<float> tempWaveform(waveforms[0].size(),0.);

    // Also recover the gain
    float asicGain = fSignalShapeService->GetASICGain(0) * fDetector->SamplingRate() / 1000.;  // something like 67.4

    // Get a base channel number for the plane we want
    raw::ChannelID_t channel = fGeometry->PlaneWireToChannel(fPlaneToSimulate, 0);

    // Loop over the wire range
    for(size_t wireIdx = fWireRange[0]; wireIdx < maxWire; wireIdx++)
    {
        // As the track angle becomes more parallel to the electron drift direction there will be more charge
        // deposited per mm that will impact a given wire. So we need to try to accommodate this here.
        // Start by computing the starting tick (based on how far we have gone so far)
        size_t startTick = size_t(fTanThetaTW * float(wireIdx - fWireRange[0])) + fTickRange[0];

        // If this has gone outside the maximum tick then we are done
        if (!(startTick < fTickRange[1] && startTick < tempWaveform.size())) break;

        // Begin by computing the number of ticks for this given wire, converted to tick units and 
        // always at least one tick
        size_t deltaTicks = size_t(fMMPerWire * std::fabs(fTanTheta)) + 1;
        size_t endTick    = startTick + deltaTicks;

        // Trim back if we look to step outside of the max range 
        endTick = std::min(endTick,fTickRange[1]);
        endTick = std::min(endTick,tempWaveform.size());

        // Ok, now loop through the ticks and deposit charge. 
        // This will be the number of electrons per mm (input) 
        // x the total arc length of track for this step
        float arcLenPerWire    = fMMPerWire / fCosTheta;
        float arcLenPerTick    = arcLenPerWire / float(deltaTicks);
        float numElectronsWire = fNumElectronsPerMM * arcLenPerWire;
        float fracElecPerTick  = numElectronsWire * arcLenPerTick / arcLenPerWire;

        // Make sure the waveform is zeroed
        std::fill(tempWaveform.begin(),tempWaveform.end(),0.);

        // Now loop through the ticks and deposit charge
        for(size_t tickIdx = startTick; tickIdx < endTick; tickIdx++)
            tempWaveform[tickIdx] = fracElecPerTick / asicGain;

        // Convolute with the response functions
        fSignalShapeService->Convolute(channel, tempWaveform);

        // And now add this to the waveform in question 
        VectorFloat& waveform = waveforms[wireIdx];

        std::transform(waveform.begin(),waveform.end(),tempWaveform.begin(),waveform.begin(),std::plus<float>());
    }

    return;
}

DEFINE_ART_CLASS_TOOL(FakeParticle)
} // namespace daq
