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
                   const uint32_t x_pc_num_bits,
                   const uint32_t x_pos_num_bits) :
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
            pc_num_bits(x_pc_num_bits),
            pos_num_bits(x_pos_num_bits)
            {}
        const uint32_t use_alt_ctr_num_bits;          ///< Use-Alt counter number of bits
        const uint32_t use_alt_num_ctr;               ///< Number of Use-Alt counters
        const uint32_t num_to_allocate_on_mispredict; ///< Num of tagged entries to alloc on mispred
        const uint32_t useful_reset_threshold;        ///< Threshold to decrement useful bit
        const std::vector<uint32_t> tagged_table_size;///< Size of Tagged table (per Tagged table)
        const std::vector<uint32_t> tag_num_bits;     ///< Size of the tag in Tagged table (per table)
        const std::vector<uint32_t> geometric_len;    ///< The geometric history length (per table)
        const uint32_t tagged_pred_hys_num_bits;      ///< Tagged table counter's hysteresis' size
        const uint32_t useful_ctr_num_bits;           ///< Tagged table's Useful Counter size
        const uint32_t bimodal_tbl_size;              ///< Bimodal table size (in entries)
        const uint32_t bimodal_hys_num_bits;          ///< Bimodal counter hysteresis size
        const uint32_t bimodal_hys_frac;              ///< Num bimodal entries per one hysteresis
        const uint32_t pc_num_bits;                   ///< Number of PC bits to used in idx/tag calc
        const uint32_t pos_num_bits;                  ///< Number of bits for 'pos'

    }; // class TageParams
}; // namespace tage
}; // namespace dabble
