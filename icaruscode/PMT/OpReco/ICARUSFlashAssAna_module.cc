////////////////////////////////////////////////////////////////////////
// Class:       ICARUSFlashAssAna
// Plugin Type: analyzer (art v3_06_03)
// File:        ICARUSFlashAssAna_module.cc
//
// Generated at Tue Jun 29 13:43:54 2021 by Andrea Scarpelli using cetskelgen
// from cetlib version v3_11_01.
//
// Module that dumps the assciation between Flashes and OpHit
//
// mailto:ascarpel@bnl.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "art_root_io/TFileService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

#include "TTree.h"

#include <vector>
#include <map>
#include <numeric> // std::accumulate



namespace opana {
  class ICARUSFlashAssAna;
}


class opana::ICARUSFlashAssAna : public art::EDAnalyzer {

  public:

    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> TriggerLabel {
        Name("TriggerLabel"),
        Comment("Label for the Trigger fragment label")
      };

      fhicl::Atom<bool> DumpWaveformsInfo {
        Name("DumpWaveformsInfo"),
        Comment("Set the option to save some aggregated waveform information")
      };

      fhicl::Sequence<art::InputTag> OpDetWaveformLabels {
        Name("OpDetWaveformLabels"),
        Comment("Tags for the raw::OpDetWaveform data products")
      };

      fhicl::Sequence<art::InputTag> OpHitLabels {
        Name("OpHitLabels"),
        Comment("Tags for the recob::OpHit data products")
      };

      fhicl::Sequence<art::InputTag> FlashLabels {
        Name("FlashLabels"),
        Comment("Tags for the recob::Flashe data products")
      };

      fhicl::Atom<float> PEOpHitThreshold {
        Name("PEOpHitThreshold"),
        Comment("Threshold in PE for an OpHit to be considered in the information calculated for a flash")
      };

      fhicl::Atom<bool> Debug {
        Name("Debug"),
        Comment("Be more verbose")
      };

    };

    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit ICARUSFlashAssAna(Parameters const& config);

    ICARUSFlashAssAna(ICARUSFlashAssAna const&) = delete;
    ICARUSFlashAssAna(ICARUSFlashAssAna&&) = delete;
    ICARUSFlashAssAna& operator=(ICARUSFlashAssAna const&) = delete;
    ICARUSFlashAssAna& operator=(ICARUSFlashAssAna&&) = delete;

    void analyze(art::Event const& e) override;

    void beginJob() override;
    void endJob() override;

    template<typename T> T Median( std::vector<T> data ) const;

    int getCryostatByChannel( const int channel );

    int getSideByChannel( const int channel );

    void processOpHits( art::Event const& e, int cryo );

    void processOpHitsFlash( std::vector<art::Ptr<recob::OpHit>> const &ophits, 
                        int &multiplicity_left, int &multiplicity_right, 
                        float &sum_pe_left, float &sum_pe_right, 
                        float *xyz, TTree *ophittree   ); 

    static std::string_view firstLine(std::string const& s, const char* endl = "\r");

  private:

    art::InputTag fTriggerLabel;
    bool fSaveWaveformInfo;
    std::vector<art::InputTag> fOpDetWaveformLabels;
    std::vector<art::InputTag> fOpHitLabels;
    std::vector<art::InputTag> fFlashLabels;
    float fPEOpHitThreshold;
    bool fDebug;


    TTree *fEventTree;
    std::vector<TTree*> fOpDetWaveformTrees;
    std::vector<TTree*> fOpFlashTrees;
    std::vector<TTree*> fOpHitTrees;
    std::vector<TTree*> fOpHitFlashTrees;

    int m_run;
    int m_event;
    int m_timestamp;
    //int m_nflashes;
    //int m_nophit;
    short m_baseline;
    short m_chargesum;
    int m_nticks;
    float m_beam_gate_start=-99999;
    float m_beam_gate_width=-99999;
    int m_beam_type=-1;
    int m_flash_id;
    int m_multiplicity;
    int m_multiplicity_left;
    int m_multiplicity_right;
    float m_sum_pe;
    float m_sum_pe_left;
    float m_sum_pe_right;
    float m_flash_time;
    //float m_flash_x;
    //float m_flash_width_x;
    float m_flash_y;
    float m_flash_width_y;
    float m_flash_z;
    float m_flash_width_z;

    int m_channel_id;
    float m_integral; // in ADC x tick
    float m_amplitude; // in ADC
    float m_start_time;
    float m_width;
    float m_abs_start_time;
    float m_pe;
    float m_fast_to_total;

    std::vector<float> m_pmt_x;
    std::vector<float> m_pmt_y;
    std::vector<float> m_pmt_z;

};


opana::ICARUSFlashAssAna::ICARUSFlashAssAna(Parameters const& config)
  : EDAnalyzer(config)
  , fTriggerLabel( config().TriggerLabel() )
  , fSaveWaveformInfo( config().DumpWaveformsInfo() )
  , fOpDetWaveformLabels( config().OpDetWaveformLabels() )
  , fOpHitLabels( config().OpHitLabels() )
  , fFlashLabels( config().FlashLabels() )
  , fPEOpHitThreshold( config().PEOpHitThreshold() )
  , fDebug( config().Debug() )
{ }


void opana::ICARUSFlashAssAna::beginJob() {

  art::ServiceHandle<art::TFileService const> tfs;

  TTree* fGeoTree = tfs->make<TTree>("geotree", "geometry information" );
  fGeoTree->Branch("pmt_x",&m_pmt_x);
  fGeoTree->Branch("pmt_y",&m_pmt_y);
  fGeoTree->Branch("pmt_z",&m_pmt_z);
  
  auto const geop = lar::providerFrom<geo::Geometry>();

  double PMTxyz[3];
  for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {

    geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);

    //std::cout << PMTxyz[0] << " " << PMTxyz[1] << " " << PMTxyz[2] << std::endl;

    m_pmt_x.push_back(PMTxyz[0]);
    m_pmt_y.push_back(PMTxyz[1]);
    m_pmt_z.push_back(PMTxyz[2]);

  }

  fGeoTree->Fill();

  fEventTree = tfs->make<TTree>("eventstree", "higher level information on the event" );
  fEventTree->Branch("run", &m_run, "run/I");
  fEventTree->Branch("event", &m_event, "event/I");
  fEventTree->Branch("timestamp", &m_timestamp, "timestamp/I");
  //fEventTree->Branch("nflashes", &m_nflashes, "nflashes/I");
  //fEventTree->Branch("nophits", &m_nophit, "nophits/I");
  fEventTree->Branch("beam_gate_start", &m_beam_gate_start, "beam_gate_start/F");
  fEventTree->Branch("beam_gate_width", &m_beam_gate_width, "beam_gate_width/F");
  fEventTree->Branch("beam_type", &m_beam_type, "beam_type/I");


  // This tree will hold some aggregated optical waveform information
  // The flag must be enabled to have the information saved
  if( !fOpDetWaveformLabels.empty() && fSaveWaveformInfo ) { 
  
    for( auto const & label : fOpDetWaveformLabels ) {

      std::string name = label.label()+"wfttree";
      std::string info = "TTree with aggregated optical waveform information with label: " + label.label();

      TTree* ttree = tfs->make<TTree>(name.c_str(), info.c_str());
      ttree->Branch("run", &m_run, "run/I");
      ttree->Branch("event", &m_event, "event/I");
      ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
      ttree->Branch("channel_id", &m_channel_id, "channel_id/I");
      ttree->Branch("baseline", &m_baseline, "baseline/s");
      ttree->Branch("chargesum", &m_chargesum, "chargesum/s");
      ttree->Branch("nticks", &m_nticks, "nticks/I");
   
      fOpDetWaveformTrees.push_back(ttree);

    }

  }


  // This ttree will hold the ophit information when a flash is not found in the event
  // NB: information of the optical hits in events where flashes are present are lost
  
  if ( !fOpHitLabels.empty() ) { 
    
    for( auto const & label : fOpHitLabels ) { 

      std::string name = label.label()+"_ttree"; 
      std::string info = "TTree for the recob::OpHit objects with label " + label.label() + " in events without flashes.";

      TTree* ttree = tfs->make<TTree>(name.c_str(), info.c_str());
      ttree->Branch("run", &m_run, "run/I");
      ttree->Branch("event", &m_event, "event/I");
      ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
      ttree->Branch("channel_id", &m_channel_id, "channel_id/I");
      ttree->Branch("integral", &m_integral, "integral/F");
      ttree->Branch("amplitude", &m_amplitude, "amplitude/F");
      ttree->Branch("start_time", &m_start_time, "start_time/F");
      ttree->Branch("abs_start_time", &m_abs_start_time, "abs_start_time/F");
      ttree->Branch("pe", &m_pe, "pe/F");
      ttree->Branch("width", &m_width, "width/F");
      ttree->Branch("fast_to_total", &m_fast_to_total, "fast_to_total/F");

      fOpHitTrees.push_back(ttree);

    }
  }



  if ( !fFlashLabels.empty() ) {

    for( auto const & label : fFlashLabels ) {

        // TTree for the flash in a given cryostat
        std::string name = label.label()+"_flashtree";
        std::string info = "TTree for the recob::Flashes with label "+label.label();

        TTree* ttree = tfs->make<TTree>(name.c_str(), info.c_str() );
        ttree->Branch("run", &m_run, "run/I");
        ttree->Branch("event", &m_event, "event/I");
        ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
        ttree->Branch("flash_id", &m_flash_id, "flash_id/I");
        ttree->Branch("multiplicity", &m_multiplicity, "multiplicity/I");
        ttree->Branch("multiplicity_right", &m_multiplicity_right, "multiplicity_right/I" );
        ttree->Branch("multiplicity_left", &m_multiplicity_left, "multiplicity_left/I" );
        ttree->Branch("sum_pe", &m_sum_pe, "sum_pe/F");
        ttree->Branch("sum_pe_right", &m_sum_pe_right, "sum_pe_right/F");
        ttree->Branch("sum_pe_left", &m_sum_pe_left, "sum_pe_left/F");
        ttree->Branch("flash_time", &m_flash_time, "flash_time/F");
        //ttree->Branch("flash_x", &m_flash_x, "flash_x/F");
        //ttree->Branch("flash_width_x", &m_flash_width_x, "flash_width_x/F");
        ttree->Branch("flash_y", &m_flash_y, "flash_y/F");
        ttree->Branch("flash_width_y", &m_flash_width_y, "flash_width_y/F");
        ttree->Branch("flash_z", &m_flash_z, "flash_z/F");
        ttree->Branch("flash_width_z", &m_flash_width_z, "flash_width_z/F");

        fOpFlashTrees.push_back( ttree );

        // Now the ttree for the OpHit associated in the flash
        name = label.label()+"_ophittree";
        info = "Three for the recob::OpHit associated with an OpHitFlash"+label.label();

        TTree* ophittree = tfs->make<TTree>(name.c_str(), info.c_str() );
        ophittree->Branch("run", &m_run, "run/I");
        ophittree->Branch("event", &m_event, "event/I");
        ophittree->Branch("timestamp", &m_timestamp, "timestamp/I");
        ophittree->Branch("flash_id", &m_flash_id, "flash_id/I");
        ophittree->Branch("channel_id", &m_channel_id, "channel_id/I");
        ophittree->Branch("integral", &m_integral, "integral/F");
        ophittree->Branch("amplitude", &m_amplitude, "amplitude/F");
        ophittree->Branch("start_time", &m_start_time, "start_time/F");
        ophittree->Branch("abs_start_time", &m_abs_start_time, "abs_start_time/F");
        ophittree->Branch("pe", &m_pe, "pe/F");
        ophittree->Branch("width", &m_width, "width/F");
        ophittree->Branch("fast_to_total", &m_fast_to_total, "fast_to_total/F");

        fOpHitFlashTrees.push_back( ophittree );

    }
  }

}


    
template<typename T>
  T opana::ICARUSFlashAssAna::Median( std::vector<T> data ) const {

    std::nth_element( data.begin(), data.begin() + data.size()/2, data.end() );
    
    return data[ data.size()/2 ];

  }



int opana::ICARUSFlashAssAna::getCryostatByChannel( const int channel ) {


  /*
  - [0:180]: East cryostat
  - [180:360]: West cryostat
  */


  int cryo = channel / 180; // always round down

  return cryo;

}


int opana::ICARUSFlashAssAna::getSideByChannel( const int channel ) {

  /* 
  Channels are numbered from east to west, from North (cryo side) to South (beam side)
  We look in the opposide direction wrt to the beam direction South->North: 

  - Left is the east wall of each cryostat;

  - Right is the west side of each cryostat;

  - [ 0:89 ] and [180:269] are on the left, 
    the return value of the function is 0;

  - [ 90-179 ] and [ 270:359 ] are on the right,
    the return value of the function is 1;
  */


  int side = channel / 90; // always round down

  return side % 2;
}


void opana::ICARUSFlashAssAna::processOpHits( art::Event const& e, int cryo ) {


  if( !fOpHitLabels.empty() ){  

    for( size_t i=0; i<fOpHitLabels.size(); i++ ) { 

      auto const label = fOpHitLabels[i];

      art::Handle<std::vector<recob::OpHit>> ophit_handle;
      e.getByLabel( label, ophit_handle );


      // We want our flashes to be valid and not empty
      if( ophit_handle.isValid() && !ophit_handle->empty() ) {


        for( auto const & ophit : *ophit_handle ) {

          //auto const & ophit = (*ophit_handle)[idx];

          const int channel_id = ophit.OpChannel(); 

          if( getCryostatByChannel(channel_id) != cryo ){ continue; }

          m_channel_id = channel_id;
          m_integral = ophit.Area(); // in ADC x tick
          m_amplitude = ophit.Amplitude(); // in ADC
          m_start_time = ophit.PeakTime();
          m_width = ophit.Width();
          m_abs_start_time = ophit.PeakTimeAbs();
          m_pe = ophit.PE();
          m_fast_to_total = ophit.FastToTotal();

          fOpHitTrees[i]->Fill();

        }
      }

      else { 
        if( fDebug ){ mf::LogError("ICARUSFlashAssAna")
              << "Invalid recob::OpHit with label"+label.label()+"\n"; }
      }
    }
  }

  else {
    if( fDebug ){ mf::LogError("ICARUSFlashAssAna")
                << "No recob::OpHit labels selected\n"; }

  }

  return;

} 


void opana::ICARUSFlashAssAna::processOpHitsFlash( std::vector<art::Ptr<recob::OpHit>> const &ophits, 
                                              int &multiplicity_left, int &multiplicity_right, 
                                              float &sum_pe_left, float &sum_pe_right, 
                                              float *xyz, TTree *ophittree  ) {


  std::unordered_map<int, float > sumpe_map;

  // We caluclate the total charge clustered in the flash per channel taking part to the flash
  for( auto const ophit : ophits ) {

    if ( ophit->PE() < fPEOpHitThreshold ) { continue; }

    const int channel_id = ophit->OpChannel(); 

    sumpe_map[ channel_id ]+=ophit->PE() ;

    //xyz[0] += m_pmt_x[channel_id]*ophit->PE();
    //xyz[1] += m_pmt_y[channel_id]*ophit->PE();
    //xyz[2] += m_pmt_z[channel_id]*ophit->PE();

    m_channel_id = channel_id;
    m_integral = ophit->Area(); // in ADC x tick
    m_amplitude = ophit->Amplitude(); // in ADC
    m_start_time = ophit->PeakTime();
    m_width = ophit->Width();
    m_abs_start_time = ophit->PeakTimeAbs();
    m_pe = ophit->PE();
    m_fast_to_total = ophit->FastToTotal();

    ophittree->Fill();

  }

  m_multiplicity_left = std::accumulate( sumpe_map.begin(), sumpe_map.end(), 0,
           [&](int value, const std::map<int, float>::value_type& p) {
             return getSideByChannel(p.first)==0 ? ++value : value ;
           });

  m_multiplicity_right =std::accumulate( sumpe_map.begin(), sumpe_map.end(), 0,
           [&](int value, const std::map<int, float>::value_type& p) {
             return getSideByChannel(p.first)==1 ? ++value : value ;
           });

  m_sum_pe_left = std::accumulate( sumpe_map.begin(), sumpe_map.end(), 0.0, 
           [&](float value, const std::map<int, float>::value_type& p) {
             return getSideByChannel(p.first)==0 ? value+p.second : value ; 
           });
  
  m_sum_pe_right = std::accumulate( sumpe_map.begin(), sumpe_map.end(), 0.0,
            [&](float value, const std::map<int, float>::value_type& p) {
              return getSideByChannel(p.first)==1 ? value+p.second : value ;
            }); 

  //for( int i=0; i<3; i++ ){ xyz[i] /= (m_sum_pe_left+ m_sum_pe_right); }
  
}


void opana::ICARUSFlashAssAna::endJob() {

}


void opana::ICARUSFlashAssAna::analyze(art::Event const& e) {


  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second 


  /*
  This part is for the trigger information
  */


  // We work out the trigger information here 
  if( !fTriggerLabel.empty() ) { 

      art::Handle<std::vector<sim::BeamGateInfo>> beamgate_handle;
      e.getByLabel( fTriggerLabel, beamgate_handle );

      if( beamgate_handle.isValid() ) {

        for( auto const & beamgate : *beamgate_handle ) {

          m_beam_gate_start = beamgate.Start(); 
          m_beam_gate_width = beamgate.Width(); 
          m_beam_type = beamgate.BeamType() ;

        }

      }

      else {
        if( fDebug) { std::cout << "Invalid Trigger Data product " << fTriggerLabel.label() << "\n" ; }
      }

  }

  else {
     if( fDebug ){ std::cout << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ; }
  }


  /* 
  Now we work on the waveforms if we are allowed to
  */

    if( !fOpDetWaveformLabels.empty() && fSaveWaveformInfo ) { 
    
      for( size_t i=0; i<fOpDetWaveformLabels.size(); i++ ) {

        auto const label = fOpDetWaveformLabels[i];

        art::Handle<std::vector<raw::OpDetWaveform>> wfm_handle;
        e.getByLabel( label, wfm_handle );

        if( wfm_handle.isValid() && !wfm_handle->empty() ) {

          for( auto const & wave : *wfm_handle ){

            m_channel_id = wave.ChannelNumber();
            m_nticks = wave.Waveform().size();
            m_baseline = Median( wave.Waveform() ); 
            m_chargesum = std::accumulate( wave.Waveform().begin(), wave.Waveform().end(), 0, 
             [ & ](short x, short y){ return ((m_baseline-x) + (m_baseline-y)) ; } );

            fOpDetWaveformTrees[i]->Fill();
            
          }
        }
      }
    }


  /*
    Now we take care of the flashes: we separate the case where we have a flash and the case where we have not a flash
  */

  if ( !fFlashLabels.empty() ) {

    for ( size_t i=0; i<fFlashLabels.size(); i++  ) {

      auto const label = fFlashLabels[i];

      art::Handle<std::vector<recob::OpFlash>> flash_handle;
      e.getByLabel( label, flash_handle );


      // We want our flashes to be valid and not empty
      if( flash_handle.isValid() && !flash_handle->empty() && flash_handle->size() > 0 ) {

        art::FindManyP<recob::OpHit> ophitsPtr( flash_handle, e, label );


        for ( size_t idx=0; idx<flash_handle->size(); idx++ ) {

          m_flash_id = idx;
          auto const & flash = (*flash_handle)[idx];

          m_flash_time = flash.Time();
          m_sum_pe = flash.TotalPE();
            
          auto const & ophits = ophitsPtr.at(idx);

          // Get the multiplicity, the position and the number of PE per Side
          float xyz[3] = {0.0, 0.0, 0.0};
          processOpHitsFlash( ophits, 
                              m_multiplicity_left, m_multiplicity_right, 
                                m_sum_pe_left, m_sum_pe_right, xyz, fOpHitFlashTrees[i] );

          /*
          std::cout << "\tflash id: " << idx << ", time: " << m_flash_time;
          std::cout << ", multiplicity left: " << m_multiplicity_left << ", multiplicity right: " << m_multiplicity_right;
          std::cout << ", sum pe left: " << m_sum_pe_left << ", sum pe right: " << m_sum_pe_right;
          std::cout << " coor: [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "]";
          std::cout << " coor: [" << 0.0 << ", " << flash.YCenter() << ", " << flash.ZCenter() << "]";
          std::cout << " coor: [" << 0.0 << ", " << flash.YWidth() << ", " << flash.ZWidth() << "]";
          std::cout  << "\n";
          */

          m_multiplicity = m_multiplicity_left+m_multiplicity_right;

          //m_flash_x = 0.0;
          //m_flash_width_x = 0.0;
          m_flash_y = flash.YCenter(); 
          m_flash_width_y = flash.YWidth();
          m_flash_z = flash.ZCenter();
          m_flash_width_z = flash.ZWidth();

          fOpFlashTrees[i]->Fill();

        }

      }

      else {

        if( fDebug) { mf::LogError("ICARUSFlashAssAna")
           << "Invalid recob::OpFlash with label"+label.label()+"\n"; }


        // We save the ophits anyways in absence of flashes
        processOpHits(e, i);

      }

    } 
  }

  else {
    
    if( fDebug ) { mf::LogError("ICARUSFlashAssAna")
          << "No recob::OpFlash labels selected\n"; }

    // We save the ophits anyways even in absence of flashes
    processOpHits(e, 0);
    processOpHits(e, 1);

  }


  fEventTree->Fill();

}


DEFINE_ART_MODULE(opana::ICARUSFlashAssAna)