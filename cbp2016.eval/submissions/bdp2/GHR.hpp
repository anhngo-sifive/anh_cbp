#pragma once

#include <cinttypes>
#include <bitset>
#include <cassert>

namespace dabble {
namespace tage {

    // A simple GHR class, leveraging std::bitset
    class GHR
    {
    public:
        static const uint32_t MAXLEN=512;
        using BitVector = std::bitset<MAXLEN>;

        GHR(const uint32_t pc_bits_per_grain,
            const uint32_t tgt_bits_per_grain,
            const uint32_t hist_bits_per_br) :
            pc_bits_per_grain_(pc_bits_per_grain),
            tgt_bits_per_grain_(tgt_bits_per_grain),
            hist_bits_per_br_(hist_bits_per_br),
            total_bits_to_hash_(pc_bits_per_grain_+tgt_bits_per_grain_),
            pc_mask_((0x1ULL << pc_bits_per_grain)-1),
            tgt_mask_((0x1ULL << tgt_bits_per_grain)-1),
            sig_mask_((0x1ULL << hist_bits_per_br)-1)
        {
            sparta_assert(total_bits_to_hash_ <= sizeof(sig_mask_)*8,
                          "Too many bits from PC and target to compute signature");
        }

        // Add a new history bit, 1 for taken, 0 for not taken
        void updateWithTaken(bool taken) {
            ghr_ = (ghr_<< 1);
            ghr_.set(0,taken);
        }

        void updateWithPcAndTarget(const uint64_t pc,
                                   const uint64_t target)
        {

            // Updating GHR with PC and target (see MAL-583)
            // 1. create a signature by folding down PC+target
            // 2. shift GHR by 1, and XOR in signature
            //
            // NOTE that the PC and target should be pre-shifted to account
            //      for either 2B or 4B inst width

            uint64_t pc_plus_tgt = ((pc & pc_mask_)<<tgt_bits_per_grain_) | (target & tgt_mask_);

            // Fold pc_plus_tgt down to hist_bits_per_br
            uint64_t signature =  fold_(pc_plus_tgt, total_bits_to_hash_, hist_bits_per_br_);

            // shift GHR by 1 and XOR in the signature
            ghr_ <<= 1;
            ghr_ ^= signature;    // XXX (PERF) a bitset is created for XORing
        }

        // Accessor
        bool operator[](uint32_t index) const {
            assert(index<MAXLEN);
            return ghr_[index];
        }

        const BitVector &getBitVector() const { return ghr_; }
        BitVector &getBitVector() { return ghr_; }

        // Convert of a string of 0's and 1's
        std::string to_string() const { return ghr_.to_string(); }

        // Convert to uint64_t.  Will throw exception if overflow.
        uint64_t to_ullong() const { return ghr_.to_ullong(); }
    private:
        // We always store MAXLEN number of bits, the
        // actual number of bits used is per user
        BitVector ghr_;

        constexpr uint64_t fold_(uint64_t val, uint32_t val_len, uint32_t hash_len) {
            val &= ((0x1ULL<<val_len) - 1);
            uint64_t hash = 0;
            uint32_t bits_consumed = 0;
            while (bits_consumed<val_len) {
                hash    ^= val;
                val     >>= hash_len;
                bits_consumed  += hash_len;
            }
            hash &= ((0x1ULL<<hash_len)-1);
            return hash;
        }

        const uint32_t pc_bits_per_grain_;    ///< Num of PC bits to use
        const uint32_t tgt_bits_per_grain_;   ///< Num of target bits to use
        const uint32_t hist_bits_per_br_;     ///< Num of signature bits per branch
        const uint32_t total_bits_to_hash_;   ///< Total num of pc and target bits
        const uint64_t pc_mask_;
        const uint64_t tgt_mask_;
        const uint64_t sig_mask_;
    }; // class GHR
}; // namespace tage
}; // namespace dabble
