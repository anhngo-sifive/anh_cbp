#pragma once

#include <cinttypes>
#include <cassert>
#include <vector>

namespace dabble {
namespace tage {
    struct TageParams
    {
        TageParams(const uint32_t x_use_alt_ctr_num_bits,
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
            use_alt_ctr_num_bits(x_use_alt_ctr_num_bits),
            use_alt_num_ctr(x_use_alt_num_ctr),
            num_to_allocate_on_mispredict(x_num_to_allocate_on_mispredict),
            useful_reset_threshold(x_useful_reset_threshold),
            tagged_table_size(x_tagged_table_size),
            tag_num_bits(x_tag_num_bits),
            geometric_len(x_geometric_len),
            tagged_pred_hys_num_bits(x_pred_hys_num_bits),
            useful_ctr_num_bits(x_useful_ctr_num_bits),
            bimodal_tbl_size(x_bimodal_tbl_size),
            bimodal_hys_num_bits(x_bimodal_hys_num_bits),
            bimodal_hys_frac(x_bimodal_hys_frac),
            addr_num_bits(x_addr_num_bits)
            {}
        const uint32_t use_alt_ctr_num_bits;          ///< Use-Alt counter number of bits
        const uint32_t use_alt_num_ctr;               ///< Number of Use-Alt counters
        const uint32_t num_to_allocate_on_mispredict; ///< Number of tagged entries to alloc on mispredict
        const uint32_t useful_reset_threshold;        ///< Threshold to decrement useful bit
        const std::vector<uint32_t> tagged_table_size;///< Size of Tagged table (per Tagged table)
        const std::vector<uint32_t> tag_num_bits;     ///< Size of the tag in Tagged table (per Tagged table)
        const std::vector<uint32_t> geometric_len;    ///< The geometric history length (per Tagged table)
        const uint32_t tagged_pred_hys_num_bits;      ///< Tagged table counter's hysteresis' size
        const uint32_t useful_ctr_num_bits;           ///< Tagged table's Useful Counter size
        const uint32_t bimodal_tbl_size;              ///< Bimodal table size (in entries)
        const uint32_t bimodal_hys_num_bits;          ///< Bimodal counter hysteresis size
        const uint32_t bimodal_hys_frac;              ///< Number of bimodal entries sharing one hysteresis
        const uint32_t addr_num_bits;                 ///< Number of address bit (PC or tgt) to used in hashing

    }; // class TageParams
}; // namespace tage
}; // namespace dabble
