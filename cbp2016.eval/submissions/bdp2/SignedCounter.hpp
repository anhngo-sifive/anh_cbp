#pragma once

#include <cinttypes>
#include <cassert>

namespace dabble {
namespace tage {

    // Signed Prediction Counter
    // Note that this class is used in branch prediction,
    // and not related to sparta::Counter
    class SignedCounter
    {
    public:
        explicit SignedCounter(uint32_t numbits) :
            MAX_((1<<(numbits-1))-1),
            MIN_(-(1<<(numbits-1))),
            counter_(0) // init to weak positive
        {
            assert(numbits<(sizeof(counter_) * 8));
        }

        void reset(bool taken) {
            counter_ = taken ? 0 : -1;
        }

        void update(bool dir) {
            if (dir) {
                if (counter_ < MAX_) {
                    ++ counter_;
                }
            }
            else  {
                if (counter_ > MIN_) {
                    -- counter_;
                }
            }
        }

        bool isGTEZero() const { return (counter_ >=0) ; }
        bool isWeak() const { return (counter_==0) || (counter_==-1); }
        int32_t getValue() const { return counter_; }
    private:
        int32_t MAX_=0;
        int32_t MIN_=0;
        int32_t counter_ = 0;
    }; // class PredCounter
}; // namespace tage
}; // namespace dabble
