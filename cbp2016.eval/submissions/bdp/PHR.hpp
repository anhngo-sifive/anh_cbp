#pragma once

#include <cinttypes>
#include <cassert>

namespace dabble {
namespace tage {
    class PHR
    {
    public:
        explicit PHR(uint32_t length) :
            length_(length),
            length_mask_((1ULL << length) -1),
            phr_(0)   {
            assert(length < (sizeof(phr_)*8));
        }

        void addHistory(uint64_t path, uint32_t num_bits)
        {
            phr_ <<= num_bits;
            phr_ += (path & ((1 <<num_bits) - 1));
            phr_ &= length_mask_;
        }

        void addHistoryBit(bool new_bit)
        {
            phr_ <<= 1;
            phr_ += new_bit;
            phr_ &= length_mask_;
        }

        const uint64_t &getHistory() const { return phr_; }
        uint32_t getLength() const { return length_; }
    private:
        uint32_t length_=0;
        uint64_t length_mask_=0;
        uint64_t phr_=0;
    }; // class PHR
}; // namespace tage
}; // namespace dabble
