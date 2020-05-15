#pragma once

#include <cinttypes>
#include <cassert>
#include <vector>
#include "TageParams.hpp"


namespace dabble {

namespace bdp {
    struct BdpParams : public tage::TageParams
    {
       BdpParams(const uint32_t x_use_alt_ctr_num_bits,
                 const uint32_t x_use_alt_num_ctr,
                 const uint32_t x_num_to_allocate_on_mispredict,
                 const uint32_t x_useful_reset_threshold,
                 const std::vector<uint32_t> &x_tagged_table_size,
                 const std::vector<uint32_t> &x_tag_num_bits,
                 const std::vector<uint32_t> &x_geometric_len,
                 const uint32_t x_pred_hys_num_bits,
                 const uint32_t x_useful_ctr_num_bits,
                 const uint32_t x_bimodal_tbl_size,
                 const uint32_t x_bimodal_hys_num_bits,
                 const uint32_t x_bimodal_hys_frac,
                 const uint32_t x_addr_num_bits) :
           tage::TageParams( x_use_alt_ctr_num_bits,
                             x_use_alt_num_ctr,
                             x_num_to_allocate_on_mispredict,
                             x_useful_reset_threshold,
                             x_tagged_table_size,
                             x_tag_num_bits,
                             x_geometric_len,
                             x_pred_hys_num_bits,
                             x_useful_ctr_num_bits,
                             x_bimodal_tbl_size,
                             x_bimodal_hys_num_bits,
                             x_bimodal_hys_frac,
                             x_addr_num_bits)
        {}

    }; // class BdpParams

}; // namespace bdp
}; // namespace dabble
