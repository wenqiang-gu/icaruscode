////////////////////////////////////////////////////////////////////////
// Class:       ICARUSTPCCompress
// Plugin Type: producer (art v3_03_01)
// File:        ICARUSTPCCompress_module.cc
//
// Generated at Fri Jan 10 13:58:09 2020 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"

// #include <zstd.h>
// #include "zlib.h"

class ICARUSTPCCompress;


class ICARUSTPCCompress : public art::EDProducer {
public:
  explicit ICARUSTPCCompress(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSTPCCompress(ICARUSTPCCompress const&) = delete;
  ICARUSTPCCompress(ICARUSTPCCompress&&) = delete;
  ICARUSTPCCompress& operator=(ICARUSTPCCompress const&) = delete;
  ICARUSTPCCompress& operator=(ICARUSTPCCompress&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

};


ICARUSTPCCompress::ICARUSTPCCompress(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  produces<std::vector<raw::RawDigit>>();
}

void ICARUSTPCCompress::produce(art::Event& e)
{

  std::unique_ptr<std::vector<raw::RawDigit>> product_collection(new std::vector<raw::RawDigit>);

  auto const& rawdigits_handle = e.getValidHandle<std::vector<raw::RawDigit>>("daq");
  std::cout << "rawdigits size=" << rawdigits_handle->size() << std::endl;
  for (const raw::RawDigit& rd : *rawdigits_handle) {

    // compress using Huffman method is larsoft
    raw::RawDigit::ADCvector_t cADCs(rd.ADCs());
    raw::Compress(cADCs,raw::kHuffman);
    product_collection->emplace_back(rd.Channel(),rd.NADC(),cADCs,raw::kHuffman);

    // test uncompression
    raw::RawDigit::ADCvector_t uADCs(rd.NADC());
    raw::Uncompress(cADCs,uADCs,raw::kHuffman);
    std::cout << "channel=" << rd.Channel() << " size in=" << rd.NADC() << " out=" << cADCs.size() << " back=" << uADCs.size() << std::endl;
    for (size_t i = 0; i<rd.NADC(); i++) {
      if (rd.ADC(i)-uADCs[i] != 0) {	
	std::cout << rd.ADC(i) << " " << uADCs[i] << " - diff=" << rd.ADC(i)-uADCs[i] << std::endl;
	exit(1);
      }
    }

    // size_t size_of_data = sizeof(short)*rd.NADC();
    // size_t estimated_size_of_compressed_data = ZSTD_compressBound(size_of_data);
    // void* compressed_data = malloc(estimated_size_of_compressed_data);
    // size_t actual_size_of_compressed_data = ZSTD_compress(compressed_data, estimated_size_of_compressed_data, &ADCs, size_of_data, 22);
    // //size_t estimated_size_of_decompressed_data = ZSTD_decompressBound(actual_size_of_compressed_data);
    // size_t estimated_size_of_decompressed_data = size_of_data;
    // void* decompressed_data = malloc(estimated_size_of_decompressed_data);

    // size_t const dSize = ZSTD_decompress(decompressed_data, estimated_size_of_decompressed_data, compressed_data, actual_size_of_compressed_data);

    // std::cout << "zstd estimated_size_of_compressed_data=" << estimated_size_of_compressed_data 
    // 	      << " actual_size_of_compressed_data=" << actual_size_of_compressed_data 
    // 	      << " dSize=" << dSize
    // 	      << std::endl;

    // short* decompressed_data_s = static_cast<short *>(decompressed_data);

    // for (size_t i = 0; i<rd.NADC(); i++) {
    //   if (rd.ADC(i)-decompressed_data_s[i] != 0) {	
    // 	std::cout << rd.ADC(i) << " " << decompressed_data_s[i] << " - diff=" << rd.ADC(i)-decompressed_data_s[i] << std::endl;
    // 	exit(1);
    //   }
    // }

  }

  e.put(std::move(product_collection));
}

DEFINE_ART_MODULE(ICARUSTPCCompress)
