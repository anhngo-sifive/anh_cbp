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
        void setPred(PredT pred,
                     bool is_weak_pred,
                     int32_t idx,
                     uint64_t tag,
                     uint32_t pos)
        {
            pred_ = pred;
            is_weak_pred_ = is_weak_pred;
            idx_ = idx;
            tag_ = tag;
            pos_ = pos;
        }

        PredT getPred() const { return pred_; }
        bool getWeakPred() const { return is_weak_pred_; }
        int32_t getIdx() const { return idx_; }
        uint64_t getTag() const { return tag_; }
        uint32_t getPos() const { return pos_; }
    private:
        PredT    pred_ = false;         // prediction value
        bool     is_weak_pred_ = false; // whether the prediction is weak
        int32_t  idx_ = -1;    //  idx used for lookup
        uint64_t tag_ = 0;     //  tag used for lookup
        uint32_t pos_ = 0;     //  pos used for lookup

    }; // class TaggedTablePredResult
}; // namespace tage
}; // namespace dabble
