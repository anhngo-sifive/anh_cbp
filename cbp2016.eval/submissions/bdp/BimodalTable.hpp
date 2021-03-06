#pragma once

#include <cinttypes>
#include <cassert>
#include <vector>
//#include "sparta/utils/SpartaAssert.hpp"
#include "PredCounter.hpp"

namespace dabble {
namespace tage {
    template<typename PredT>
    class BimodalTable
    {

    public:
        BimodalTable(uint32_t tbl_sz,
                     uint32_t hys_numbits,
                     uint32_t hys_frac)  // # of entries that share same hysteresis counter
        {
            //sparta_assert(((tbl_sz&(tbl_sz-1))==0),
            //              "Expected power-of-2, tbl_sz=" << tbl_sz);
            //sparta_assert(((hys_frac&(hys_frac-1))==0),
            //              "Expected power-of-2, hys_frac=" << hys_frac);
            tbl_size_ = tbl_sz;
            tbl_size_mask_ = tbl_size_-1;
            hys_frac_mask_ = ~(hys_frac-1);
            bimodal_tbl_.resize(tbl_size_, PredCounter<PredT>(hys_numbits));
        }

        void initialize(PredT initial_pred_val) {
            for (uint32_t i=0; i<tbl_size_; ++i) {
                bimodal_tbl_[i].reset(initial_pred_val);
            }
        }

        void updatePredCounter(const int32_t idx, const PredT actual_val) {
            // Read and use hysteresis from the shared entry
            uint32_t hysteresis = bimodal_tbl_.at(idx&hys_frac_mask_).getHysteresis();
            bimodal_tbl_.at(idx).setHysteresis(hysteresis);
            bimodal_tbl_[idx].update(actual_val);

            // Update the shared entry with new hysteresis
            hysteresis = bimodal_tbl_[idx].getHysteresis();
            bimodal_tbl_[idx&hys_frac_mask_].setHysteresis(hysteresis);
            std::cout << " Base ctrupdate idx=" << std::dec << (idx&hys_frac_mask_) << std::endl;
        }

        void lookupPrediction(const uint64_t PC_shifted,
                              PredT &pred,
                              int32_t &idx) const {
            idx = calcIdx_(PC_shifted);
            pred =  bimodal_tbl_.at(idx).getPred();
        }

        BimodalTable(const BimodalTable&) = delete;
        BimodalTable& operator=(const BimodalTable&) = delete;

    private:

        int32_t calcIdx_(uint64_t PC_shifted) const { return (PC_shifted & tbl_size_mask_); }
        uint32_t tbl_size_=0;
        uint32_t tbl_size_mask_=0;
        uint32_t hys_frac_mask_=0; ///< mask to determine shared hysteresis entry
        std::vector<PredCounter<PredT>> bimodal_tbl_;
    }; // class BimodalTable
}; // namespace tage
}; // namespace dabble
