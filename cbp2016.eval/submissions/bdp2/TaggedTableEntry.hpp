#pragma once

#include <cinttypes>
#include <iostream>
#include "UnsignedCounter.hpp"
#include "PredCounter.hpp"

namespace dabble {
namespace tage {

    template<typename PredT>
    class TaggedTableEntry {
    public:
        TaggedTableEntry(uint32_t ctr_num_bits,
                         uint32_t useful_num_bits) :
            tag_(0),
            pred_ctr_(ctr_num_bits),
            useful_ctr_(useful_num_bits)
        {
        }

        void reset(const uint64_t tag,
                   const uint32_t pos,
                   const PredT &init_pred_val) {
            tag_ = tag;
            pos_ = pos;
            pred_ctr_.reset(init_pred_val);
            useful_ctr_.resetToZero();
        }

        void resetUsefulCounter() { useful_ctr_.resetToZero(); }
        void updateUsefulCounter(bool dir) { useful_ctr_.update(dir); }
        void updatePredCounter(PredT actual_val) { pred_ctr_.update(actual_val); }

        PredT getPred() const { return pred_ctr_.getPred(); }
        bool isPredWeak() const { return pred_ctr_.isWeak(); }
        uint64_t getTag() const { return tag_; }
        uint32_t getPos() const { return pos_; }
        const UnsignedCounter &getUsefulCounter() const { return useful_ctr_; }

        friend std::ostream& operator<<(std::ostream& os, const TaggedTableEntry& entry)
        {
            os << "["
               <<  std::hex << entry.tag_ << std::dec
               << "|" << entry.pred_ctr_.getValue()
               << "|" << entry.useful_ctr_.getValue()
               << "]";
            return os;
        }

    private:
        uint64_t            tag_=0;
        uint32_t            pos_=0;      // position within fetch group
        PredCounter<PredT>  pred_ctr_;   // prediction counter
        UnsignedCounter     useful_ctr_; // useful counter
    }; // TaggedTableEntry
}; // namespace tage
}; // namespace dabble
