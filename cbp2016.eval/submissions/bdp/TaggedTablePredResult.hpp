#pragma once

#include <cinttypes>
#include <cassert>

namespace dabble {
namespace tage {

    // PredResult for TaggedTable
    template<typename PredT>
    class TaggedTablePredResult
    {
    public:
        void setPred(PredT pred, bool is_weak_pred,
                     int32_t computed_idx, uint64_t computed_tag)
        {
            pred_ = pred;
            is_weak_pred_ = is_weak_pred;
            computed_idx_ = computed_idx;
            computed_tag_ = computed_tag;
        }

        PredT getPred() const { return pred_; }
        bool getWeakPred() const { return is_weak_pred_; }
        int32_t getComputedIdx() const { return computed_idx_; }
        uint64_t getComputedTag() const { return computed_tag_; }
    private:
        PredT    pred_ = false;         // prediction value
        bool     is_weak_pred_ = false; // whether the prediction is weak
        int32_t  computed_idx_ = -1;    // computed idx used for lookup
        uint64_t computed_tag_ = 0;     // computed tag used for lookup

    }; // class TaggedTablePredResult
}; // namespace tage
}; // namespace dabble
