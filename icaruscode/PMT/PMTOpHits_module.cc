////////////////////////////////////////////////////////////////////////
// Class:       PMTOpHits
// Plugin Type: analyzer (art v2_07_03)
// File:        PMTOpHits_module.cc
//
// Generated at Mon Sep 25 15:10:31 2017 by Andrea Falcone using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larana/OpticalDetector/SimPhotonCounter.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "nug4/ParticleNavigation/ParticleList.h"

#include "lardataobj/Simulation/sim.h"

#include  "lardataobj/RecoBase/OpHit.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

////////////////////////////////// Define some constant variable //////////////////////////////
const int nPMTs = 360;
//const int PMTs_per_TPC = 90;
const int MaxPhotons = 10000;

namespace icarus {
  class PMTOpHits;
}


class icarus::PMTOpHits : public art::EDAnalyzer 
{
public:
explicit PMTOpHits(fhicl::ParameterSet const & p);

// Plugins should not be copied or assigned.
PMTOpHits(PMTOpHits const &) = delete;
PMTOpHits(PMTOpHits &&) = delete;
PMTOpHits & operator = (PMTOpHits const &) = delete;
PMTOpHits & operator = (PMTOpHits &&) = delete;

// Required functions.
void analyze(art::Event const & e) override;

// Selected optional functions.
void beginJob() override;
//void reconfigure(fhicl::ParameterSet const & p);

private:

//TRandom* Ran;
 
TTree* opTree;
TTree* mcTree;

////////////////////////////////// Variable in tree//////////////////////////////

int event;

int event_type;

int is_Neutrino;
int Neutrino_Interaction;

////////////////////////////////// Variables in OpHit Tree////////////////////////////////// 

int opturned_PMT;

double opPMTx[nPMTs];
double opPMTy[nPMTs];
double opPMTz[nPMTs];

int opCryostat[nPMTs];	
int opTPC[nPMTs];

float opphotons_collected[nPMTs];
float opQE_photons_collected[nPMTs];

float opphoton_time[nPMTs][MaxPhotons];
float opfirstphoton_time[nPMTs];

double vertex_x;
double vertex_y;
double vertex_z;
double vertex_t;

float optotal_coll_photons;

////////////////////////////////// Variables in mcOpHit Tree////////////////////////////////// 
int mcturned_PMT;

double mcPMTx[nPMTs];
double mcPMTy[nPMTs];
double mcPMTz[nPMTs];

int mcCryostat[nPMTs];	
int mcTPC[nPMTs];

float mcphotons_collected[nPMTs];
float mcQE_photons_collected[nPMTs];


float mcphoton_time[nPMTs][MaxPhotons];
float mcfirstphoton_time[nPMTs];

float mctotal_coll_photons;
  
art::InputTag photonLabel;
art::InputTag chargeLabel;
art::InputTag ophitLabel;
art::InputTag mcophitLabel;
art::InputTag typoLabel;

};


icarus::PMTOpHits::PMTOpHits(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  ophitLabel  (p.get<art::InputTag>("ophit","ophit")),
  mcophitLabel(p.get<art::InputTag>("mcophit","mcophit")),
  typoLabel   (p.get<art::InputTag>("typo","generator"))
 // More initializers here.
{
    std::cout << " PMT OpHits constructor " << std::endl;
}

void icarus::PMTOpHits::analyze(art::Event const & evt)
{
	////////////////////////////////// Create the LArsoft services and service handle//////////////////////////////

	art::ServiceHandle<geo::Geometry> geom;

	std::vector<recob::OpHit> const& ophits        = *(evt.getValidHandle<std::vector<recob::OpHit>>(ophitLabel));
	std::vector<recob::OpHit> const& mchits      = *(evt.getValidHandle<std::vector<recob::OpHit>>(mcophitLabel));

	////////////////////////////////// Event number//////////////////////////////

	event = evt.id().event();

	std::vector< art::Handle< std::vector<simb::MCTruth> > > type;
	evt.getManyByType(type);

	for(size_t mcl = 0; mcl < type.size(); ++mcl)
	{	
		art::Handle< std::vector<simb::MCTruth> > mclistHandle = type[mcl];
	
		for(size_t m = 0; m < mclistHandle->size(); ++m)
		{
			art::Ptr<simb::MCTruth> mct(mclistHandle, m);	

			event_type=mct->GetParticle(0).PdgCode();	

			vertex_x=mct->GetParticle(0).Vx();	
			vertex_y=mct->GetParticle(0).Vy();
	        vertex_z=mct->GetParticle(0).Vz();
			vertex_t=mct->GetParticle(0).T();

			if (event_type==12||event_type==-12||event_type==14||event_type==-14||event_type==16||event_type==-16)
			{
				is_Neutrino=1;
				Neutrino_Interaction=mct->GetNeutrino().InteractionType();
			}	
			else
			{
				is_Neutrino=0;
				Neutrino_Interaction=-9999;			
			}		
		}
	}

	opturned_PMT=0;
	optotal_coll_photons=0;

	mcturned_PMT=0;
	mctotal_coll_photons=0;

	for (std::size_t channel = 0; channel < 360; channel++) {
		
		////////////////////////////////// Putting at 0 all the variables//////////////////////////////
		opphotons_collected[channel] = 0;
		mcphotons_collected[channel] = 0;
	
		int op_i = 0;
		int mc_i = 0;

		////////////////////////////////// OpHits part //////////////////////////////////////////////////
		for (recob::OpHit ophit : ophits) {

			std::size_t pmt = ophit.OpChannel();

			if (channel != pmt) {continue;}

			op_i++;

			opphotons_collected[channel] = opphotons_collected[channel] + ophit.PE();

			double xyz[3];
			geom->OpDetGeoFromOpChannel(channel).GetCenter(xyz);

			opPMTx[channel] = xyz[0];
			opPMTy[channel] = xyz[1];
			opPMTz[channel] = xyz[2]; 

			opfirstphoton_time[channel] = 100000000;

			// Get the earliest photon time for each channel
			if (opphotons_collected[channel]>0) {
				opphoton_time[channel][op_i] = ophit.PeakTime();
				if (ophit.PeakTime() < opfirstphoton_time[channel]) {
					opfirstphoton_time[channel] = ophit.PeakTime();
				}
			}
		}

		if (opphotons_collected[channel]>0){
			opturned_PMT++;
			optotal_coll_photons= optotal_coll_photons + opphotons_collected[channel];
		}

		////////////////////////////////// mcOpHits part //////////////////////////////////////////////////
		for (recob::OpHit mchit : mchits) {

			std::size_t pmt = mchit.OpChannel();

			if (channel != pmt) {continue;}

			mcphotons_collected[channel] = mcphotons_collected[channel] + mchit.PE();

			double xyz[3];
			geom->OpDetGeoFromOpChannel(channel).GetCenter(xyz);

			mcPMTx[channel] = xyz[0];
			mcPMTy[channel] = xyz[1];
			mcPMTz[channel] = xyz[2]; 

			mcfirstphoton_time[channel] = 100000000;

			// Get the earliest photon time for each channel
			if (mcphotons_collected[channel]>0) {
				mcphoton_time[channel][mc_i] = mchit.PeakTime();
				if (mchit.PeakTime() < mcfirstphoton_time[channel]) {
					mcfirstphoton_time[channel] = mchit.PeakTime();
				}
			}
		}

		if (mcphotons_collected[channel]>0){
			mcturned_PMT++;
			mctotal_coll_photons= mctotal_coll_photons + mcphotons_collected[channel];
		}

		// Get PMT position information for each channel
		double xyz[3];

		geom->OpDetGeoFromOpChannel(channel).GetCenter(xyz);

		opPMTx[channel] = xyz[0];
		opPMTy[channel] = xyz[1];
		opPMTz[channel] = xyz[2]; 

		mcPMTx[channel] = xyz[0];
		mcPMTy[channel] = xyz[1];
		mcPMTz[channel] = xyz[2]; 

		if (opPMTx[channel]<0){opCryostat[channel]=0;}
		if (opPMTx[channel]>0){opCryostat[channel]=1;}

		if (opPMTx[channel]<-200){opTPC[channel]=0;}
		if (opPMTx[channel]>-200 && opPMTx[channel]<0){opTPC[channel]=1;}
		if (opPMTx[channel]<200 && opPMTx[channel]>0){opTPC[channel]=2;}
		if (opPMTx[channel]>200){opTPC[channel]=3;}

		if (mcPMTx[channel]<0){mcCryostat[channel]=0;}
		if (mcPMTx[channel]>0){mcCryostat[channel]=1;}

		if (mcPMTx[channel]<-200){mcTPC[channel]=0;}
		if (mcPMTx[channel]>-200 && mcPMTx[channel]<0){mcTPC[channel]=1;}
		if (mcPMTx[channel]<200 && mcPMTx[channel]>0){mcTPC[channel]=2;}
		if (mcPMTx[channel]>200){mcTPC[channel]=3;}
	}

	opTree->Fill();
	std::cout << " finished filling " << "opTree" << std::endl;
	mcTree->Fill();
	std::cout << " finished filling " << "mcTree" << std::endl;
}

void icarus::PMTOpHits::beginJob()
{
    std::cout << " PMTOpHits beginjob " << std::endl;

	art::ServiceHandle<art::TFileService> tfs;

	opTree = tfs->make<TTree>("ophittree","tree for the ophit response");

	opTree->Branch("event",&event,"event/I");
	opTree->Branch("event_type",&event_type,"event_type/I");
	opTree->Branch("is_Neutrino",&is_Neutrino,"is_Neutrino/I");
	opTree->Branch("Neutrino_Interaction",&Neutrino_Interaction,"Neutrino_Interaction/I");
	opTree->Branch("Cryostat",&opCryostat,("Cryostat[" + std::to_string(nPMTs) + "]/I").c_str());
	opTree->Branch("TPC",&opTPC,("TPC[" + std::to_string(nPMTs) + "]/I").c_str());
	opTree->Branch("PMTx",&opPMTx,("PMTx[" + std::to_string(nPMTs) + "]/D").c_str());
	opTree->Branch("PMTy",&opPMTy,("PMTy[" + std::to_string(nPMTs) + "]/D").c_str());
	opTree->Branch("PMTz",&opPMTz,("PMTz[" + std::to_string(nPMTs) + "]/D").c_str());
	opTree->Branch("turned_PMT",&opturned_PMT,"turned_PMT/I");
	opTree->Branch("total_coll_photons",&optotal_coll_photons,"total_coll_photons/F");
	opTree->Branch("photons_collected",&opphotons_collected,("photons_collected[" + std::to_string(nPMTs) + "]/F").c_str());
	opTree->Branch("firstphoton_time",&opfirstphoton_time,("firstphoton_time[" + std::to_string(nPMTs) + "]/F").c_str());
	opTree->Branch("photon_time",&opphoton_time,"photon_time[360][10000]/F");
	opTree->Branch("vertex_x",&vertex_x,"vertex_x/D");
	opTree->Branch("vertex_y",&vertex_y,"vertex_y/D");
	opTree->Branch("vertex_z",&vertex_z,"vertex_z/D");
	opTree->Branch("vertex_t",&vertex_t,"vertex_t/D");

	mcTree = tfs->make<TTree>("mcophittree","tree for the mcophit response");

	mcTree->Branch("event",&event,"event/I");
	mcTree->Branch("event_type",&event_type,"event_type/I");
	mcTree->Branch("is_Neutrino",&is_Neutrino,"is_Neutrino/I");
	mcTree->Branch("Neutrino_Interaction",&Neutrino_Interaction,"Neutrino_Interaction/I");
	mcTree->Branch("Cryostat",&mcCryostat,("Cryostat[" + std::to_string(nPMTs) + "]/I").c_str());
	mcTree->Branch("TPC",&mcTPC,("TPC[" + std::to_string(nPMTs) + "]/I").c_str());
	mcTree->Branch("PMTx",&mcPMTx,("PMTx[" + std::to_string(nPMTs) + "]/D").c_str());
	mcTree->Branch("PMTy",&mcPMTy,("PMTy[" + std::to_string(nPMTs) + "]/D").c_str());
	mcTree->Branch("PMTz",&mcPMTz,("PMTz[" + std::to_string(nPMTs) + "]/D").c_str());
	mcTree->Branch("turned_PMT",&mcturned_PMT,"turned_PMT/I");
	mcTree->Branch("total_coll_photons",&mctotal_coll_photons,"total_coll_photons/F");
	mcTree->Branch("photons_collected",&mcphotons_collected,("photons_collected[" + std::to_string(nPMTs) + "]/F").c_str());
	mcTree->Branch("firstphoton_time",&mcfirstphoton_time,("firstphoton_time[" + std::to_string(nPMTs) + "]/F").c_str());
	mcTree->Branch("photon_time",&mcphoton_time,"photon_time[360][10000]/F");
	mcTree->Branch("vertex_x",&vertex_x,"vertex_x/D");
	mcTree->Branch("vertex_y",&vertex_y,"vertex_y/D");
	mcTree->Branch("vertex_z",&vertex_z,"vertex_z/D");
	mcTree->Branch("vertex_t",&vertex_t,"vertex_t/D");
}

DEFINE_ART_MODULE(icarus::PMTOpHits)
