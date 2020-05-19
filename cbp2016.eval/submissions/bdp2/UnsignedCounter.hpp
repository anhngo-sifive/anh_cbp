
#pragma once

#include <cinttypes>
#include <cassert>

namespace dabble {
namespace tage {
    // Unsigned Prediction Counter
    // Note that this class is used in branch prediction,
    // and not related to sparta::Counter
    class UnsignedCounter
    {
    public:
        explicit UnsignedCounter(uint32_t numbits) :
            MAX_((1<<numbits)-1),
            counter_(0)
        {
            assert(numbits>0);
            assert(numbits<(sizeof(counter_) * 8));
        }

        void resetToZero() { counter_ = 0; }


        // Sat-increment if taken, sat-decrement if !taken
        void update(bool direction)
        {
            if (direction) {
                if (counter_ < MAX_) {
                    ++ counter_;
                }
            }
            else {
                if (counter_ > 0) {
                    -- counter_;
                }
            }
        }

        uint32_t getValue() const { return counter_; }
    private:
        uint32_t MAX_=0;
        uint32_t counter_=0;
    }; // class UnsignedCounter
}; // namespace tage
}; // namespace dabble
