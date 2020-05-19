#pragma once

#include <cinttypes>
#include "sparta/utils/SpartaAssert.hpp"

namespace dabble {
namespace tage {

    // Note that this class is used in branch prediction,
    // and not related to sparta::Counter
    template <typename PredT>
    class PredCounter
    {
    public:
        explicit PredCounter(uint32_t hysteresis_numbits) :
            MAX_((1<<hysteresis_numbits)-1),
            hysteresis_(0)
        {
            sparta_assert(hysteresis_numbits>0,
                          "hysteresis_numbits must be greater than 0");
            sparta_assert(hysteresis_numbits<(sizeof(hysteresis_) * 8),
                          "Cannot support hysteresis_numbits=" << hysteresis_numbits);
        }

        void reset(PredT init_pred_val) {
            pred_ = init_pred_val;
            hysteresis_ = 0;
        }

        void update(PredT actual_val) {
            if (pred_ == actual_val) {
                if (hysteresis_ < MAX_) {
                    ++ hysteresis_;
                }
            }
            else  {
                if (hysteresis_ > 0) {
                    -- hysteresis_;
                }
                else {
                    pred_ = actual_val;
                }
            }
        }

        void setHysteresis(const uint32_t hys) { hysteresis_ = hys; }
        uint32_t  getHysteresis() const { return hysteresis_; }
        bool isWeak() const { return (hysteresis_==0); }
        PredT getPred() const { return pred_; }
    private:
        uint32_t  MAX_=0;
        uint32_t hysteresis_ = 0;
        PredT    pred_;  // Predicted value
    }; // class PredCounter

}; // namespace tage
}; // namespace dabble
