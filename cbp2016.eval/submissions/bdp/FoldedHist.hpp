#pragma once

#include <cinttypes>
#include "GHR.hpp"

namespace dabble {
namespace tage {

    class FoldedHist
    {
    public:
        void initialize(uint32_t original_length, uint32_t folded_length)
        {
            assert(folded_length < (sizeof(folded_hist_)*8));
            folded_length_       = folded_length;
            original_length_     = original_length;
            evict_bit_y_XOR_pos_ = original_length % folded_length;
            final_mask_          =  (1ULL << folded_length) - 1;
        }

        // Update folded_history after the GHR has been pushed a new history bit
        void updateWithGHR(const GHR &ghr)
        {

            // See PPM pager
            // Also see "A Comprehensive Study of Dynamic Global History
            //           Branch Prediction", figure 6

            const uint64_t evict_bit_y  = ghr[original_length_]; // oldest bit
            const uint64_t insert_bit_x = ghr[0];                // newest bit

            // Push in insert_bit_x
            folded_hist_ =  (folded_hist_ << 1) + insert_bit_x;

            // XOR in the evict_bit_folded
            const uint64_t evict_bit_folded  = (folded_hist_ >> folded_length_);
            folded_hist_ ^= evict_bit_folded;

            // XOR in the evict_bit_y
            folded_hist_ ^= (evict_bit_y << evict_bit_y_XOR_pos_);

            // final mask off
            folded_hist_ &= final_mask_;
        }

        const uint64_t &getValue() const  { return folded_hist_; }
    private:
        uint32_t folded_length_=0;      // length after folding
        uint32_t original_length_=0;    // length before folding
        uint32_t evict_bit_y_XOR_pos_=0; // see PPM paper
        uint64_t final_mask_=0;
        uint64_t folded_hist_=0;

    }; // class FoldedHist

}; // namespace tage
}; // namespace dabble
