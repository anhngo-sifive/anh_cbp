#pragma once

#include <cinttypes>
#include <cassert>
#include <vector>
#include <math.h>
#include <algorithm>
//#include "sparta/utils/SpartaAssert.hpp"
#include "FoldedHist.hpp"
#include "TaggedTableEntry.hpp"
#include "TaggedTablePredResult.hpp"
#include "GHR.hpp"
#include "PHR.hpp"

namespace dabble {
namespace tage {

    template<typename PredT>
    class TaggedTable
    {
    public:
        TaggedTable(uint32_t bank_num,
                    uint32_t tbl_num_entries,
                    uint32_t tag_num_bits,
                    uint32_t pred_ctr_num_bits,
                    uint32_t useful_ctr_num_bits,
                    uint32_t geometric_ghr_len,
                    uint32_t addr_num_bits);

        void initialize(uint64_t initial_tag, PredT initial_pred_val);
        bool lookupPrediction(const uint64_t PC,
                              const GHR &ghr,
                              const PHR &phr,
                              const bool is_indirect,
                              const bool ghr_enable_mallard_hash,
                              TaggedTablePredResult<PredT> &presult) const;
        void printIdxTag(const uint64_t PC, const PHR &phr) const    {
            std::cout << std::dec << "   " << (bank_num_+1)
                      << std::setw(4) << " idx=" <<  calcIdx_cbp_(PC, phr)
                      << std::hex << std::setw(8) << " t=" << calcTag_classic_(PC)
                      << std::endl;

        }
        void updateFoldedHist(const GHR &ghr);
        TaggedTableEntry<PredT> &getEntry(const uint32_t idx) { return tagged_tbl_.at(idx); }
        const TaggedTableEntry<PredT> &getEntry(const uint32_t idx) const;
        void gracefulDecrementAllUsefulCounters();
        void setDbgOstream(std::ostream &os) { dbg_ostream_ = &os; }

        TaggedTable(const TaggedTable&) = delete;
        TaggedTable& operator=(const TaggedTable&) = delete;
        void printFoldedHist() const {
            std::cout << "  " << bank_num_+1 << ".  c_i=" << std::hex <<std::setw(8) << folded_ghr_for_idx_.getValue() << std::endl;
            std::cout << "  " << bank_num_+1 << ". c_t0=" << std::hex <<std::setw(8) << folded_ghr_for_tag_0_.getValue() << std::endl;
            std::cout << "  " << bank_num_+1 << ". c_t1=" << std::hex <<std::setw(8) << folded_ghr_for_tag_1_.getValue() << std::endl;
      }

    private:
        /**
         * \brief Mallard hash of PC and GHR (approximation)
         * \param  hash_len length of hash produced
         * \param  rotate_amount number of bits to rotate for each grain
         */
        uint64_t hashPcAndGhr_(uint64_t  PC,
                               const GHR &ghr,
                               uint32_t hash_len,
                               uint32_t rotate_amount) const;
        uint32_t calcIdx_mallard_(uint64_t PC,
                                  const GHR &ghr,
                                  uint32_t rotate_amount) const;
        uint64_t calcTag_mallard_(uint64_t PC,
                                  const GHR &ghr,
                                  uint32_t rotate_amount) const;
        uint32_t calcIdx_classic_(uint64_t PC) const;
        uint64_t calcTag_classic_(uint64_t PC) const;

        uint32_t calcIdx_cbp_(uint64_t PC, const PHR &phr)  const;

        uint64_t getLsb_(const GHR::BitVector &bvect, uint32_t num_bits) const  {
            // XXX (PERF) extracting one bit a a time
            uint64_t lsb_bits=0;
            while (num_bits>0) {
                -- num_bits;
                lsb_bits = (lsb_bits << 1) + bvect[num_bits];
            }
            return lsb_bits;
        }
        void rotateRight_(uint64_t &val, uint32_t val_len, uint32_t rot_len) const  {
            rot_len = rot_len % val_len;  // Make sure rot_len is doesn't wrap
            const uint64_t mask = (0x1ULL<<val_len)-1;
            val &= mask;
            val = (val >> rot_len) | (val << (val_len-rot_len));
            val &= mask;
        }

        // the index functions for the tagged tables uses path history as in the OGEHL predictor
        //F serves to mix path history: not very important impact
        int F (long long A, int size, int bank) const
        {
            int A1, A2;
            A = A & ((1 << size) - 1);
            A1 = (A & ((1 << tbl_idx_num_bits_) - 1));
            A2 = (A >> tbl_idx_num_bits_);
            A2 = ((A2 << bank) & ((1 << tbl_idx_num_bits_) - 1)) + (A2 >> (tbl_idx_num_bits_ - bank));
            A = A1 ^ A2;
            A = ((A << bank) & ((1 << tbl_idx_num_bits_) - 1)) + (A >> (tbl_idx_num_bits_ - bank));
            return (A);
	}


        std::ostream *dbg_ostream_=nullptr; // Ostream for debugging/logging

        const uint32_t bank_num_;
        const uint32_t tbl_num_entries_;    // number of entries, must be power-of-2
        const uint32_t tbl_idx_num_bits_;   // log2(tlb_size)
        const uint32_t tbl_idx_mask_;       // mask to apply to enforce table size
        const uint32_t tbl_tag_num_bits_;   // number of tag bits
        const uint64_t tbl_tag_mask_;       // mask to apply to enforce tag size
        const uint32_t addr_num_bits_;      // Num address bits used in hashing
        const uint64_t addr_mask_ = (0x1ULL << addr_num_bits_) - 1;
        using EntryType = TaggedTableEntry<PredT>;
        std::vector<EntryType> tagged_tbl_;

        uint32_t geometric_ghr_len_=0; // table's assigned geometric GHR length (in bits)
        uint64_t tag_mask_=0;          // mask to apply to enforce tag size

        // 3 folded GHRs: 1 for index computation, 2 for tags
        FoldedHist folded_ghr_for_idx_;
        FoldedHist folded_ghr_for_tag_0_;  // For computing tags (need 2 per PPM paper)
        FoldedHist folded_ghr_for_tag_1_;
    }; // class TaggedTable

    ///////////////////////////////////////////

    template<typename PredT>
    TaggedTable<PredT>::TaggedTable(uint32_t bank_num,
                                    uint32_t tbl_num_entries,
                                    uint32_t tag_num_bits,
                                    uint32_t pred_ctr_num_bits,
                                    uint32_t useful_ctr_num_bits,
                                    uint32_t geometric_ghr_len,
                                    uint32_t addr_num_bits) :
        bank_num_(bank_num),
        tbl_num_entries_(tbl_num_entries),
        tbl_idx_num_bits_(log2(tbl_num_entries)),
        tbl_idx_mask_(tbl_num_entries - 1),
        tbl_tag_num_bits_(tag_num_bits),
        tbl_tag_mask_((0x1ULL<<tag_num_bits)-1),
        addr_num_bits_(addr_num_bits)

    {
        //sparta_assert((tbl_num_entries&(tbl_num_entries-1))==0,
        //              "Number of entries must be power-of-2");
        geometric_ghr_len_ = geometric_ghr_len;
        tag_mask_          = (1ULL<<tag_num_bits) - 1;
        tagged_tbl_.resize(tbl_num_entries_,
                           TaggedTableEntry<PredT>(pred_ctr_num_bits, useful_ctr_num_bits));

        folded_ghr_for_idx_.initialize(  geometric_ghr_len, tbl_idx_num_bits_);
        folded_ghr_for_tag_0_.initialize(geometric_ghr_len, tag_num_bits);
        folded_ghr_for_tag_1_.initialize(geometric_ghr_len, tag_num_bits-1);
    }

    template<typename PredT>
    void TaggedTable<PredT>::initialize(uint64_t initial_tag, PredT initial_pred_value)
    {
        for (uint32_t idx=0; idx<tbl_num_entries_; ++idx) {
            getEntry(idx).reset(initial_tag, initial_pred_value);
        }
    }
    template<typename PredT>
    bool TaggedTable<PredT>::lookupPrediction(const uint64_t PC,
                                              const GHR &ghr,
                                              const PHR &phr,
                                              const bool is_indirect,
                                              const bool enable_mallard_hash,
                                              TaggedTablePredResult<PredT> &presult) const
    {
        (void) phr;

        /* Tag and index calculations use 4 history vectors:
         * 1. ghr, globally maintained and is passed in
         * 2. folded_ghr for index,  tagged-table specific
         * 3. 2 folded ghr for tags, tagged-table specific
         *
         * Calculated tag and idx are saved for later use in update
         */
        const uint64_t tag = enable_mallard_hash ?
            calcTag_mallard_(PC, ghr, bank_num_*2+is_indirect) :
            calcTag_classic_(PC);
        const int32_t  idx = enable_mallard_hash ?
            calcIdx_mallard_(PC, ghr, bank_num_*2+is_indirect) :
            calcIdx_cbp_(PC, phr);

        // Access table and determine if tag matched
        const auto &tentry = getEntry(idx);
        const PredT pred = tentry.getPred();
        const bool is_weak_pred = tentry.isPredWeak();
        const bool tag_matched = (tentry.getTag() == tag);

        presult.setPred(pred, is_weak_pred, idx, tag);

        if (dbg_ostream_) {
            *dbg_ostream_ << "TaggedTable" << bank_num_
                          << " idx=" << idx
                          << " tag=0x" << std::hex << tag << "-" << tentry.getTag() << std::dec
                          << " pred=" << pred
                          << " matched=" << tag_matched
                          << std::endl;
        }

        return tag_matched;
    }

    template<typename PredT>
    void TaggedTable<PredT>::updateFoldedHist(const GHR &ghr)
    {
        folded_ghr_for_idx_.updateWithGHR(ghr);
        folded_ghr_for_tag_0_.updateWithGHR(ghr);
        folded_ghr_for_tag_1_.updateWithGHR(ghr);
    }

    template<typename PredT>
    const TaggedTableEntry<PredT> &TaggedTable<PredT>::getEntry(const uint32_t idx) const
    {
        //sparta_assert(idx < tbl_num_entries_, "Invalid idx=" << idx);
        return tagged_tbl_.at(idx);
    }

    template<typename PredT>
    void TaggedTable<PredT>::gracefulDecrementAllUsefulCounters()
    {
        for (uint32_t idx=0; idx<tbl_num_entries_; ++idx) {
            getEntry(idx).updateUsefulCounter(false); // decrement useful counter
        }
    }

    template<typename PredT>
    uint64_t TaggedTable<PredT>::hashPcAndGhr_(uint64_t PC,
                                               const GHR  &ghr,
                                               uint32_t hash_len,
                                               uint32_t rotate_amount)  const
    {
        /*Hash of PC and GHR (Mallard implementation, see BDP.scala):
         *  1. Concatenate the GHR & PC
         *  2. Fold down to hash_len by shifting right and XORing
         *  3. To make sure each table produces a unique hash,
         *     rotate_amount rotates each grain before XORing.
         *     rotate_amount is based on bank_num & bdp_vs_ijtp
         *     (a grain is the pre-hash hash_len num of bits)
         * NOTE The PC and target should be pre-shifted to account
         *      for either 2B or 4B inst width
         * NOTE In Mallard, there's a limit to the number of XORs that can be
         *      used to compute the hash.  To stay in limit, bit sampling may
         *      be used.  Dabble always uses all bits.
         */

        // 1. concatenate GHR+PC
        const uint64_t hash_mask = (0x1ULL<<hash_len)-1;
        auto bitvector = ghr.getBitVector();
        bitvector <<= addr_num_bits_;
        bitvector |= ( PC & addr_mask_);
        const uint32_t bits_to_hash = geometric_ghr_len_ + addr_num_bits_;

        // 2. Fold down, rotating each grain as we go
        uint64_t hash = 0;
        uint32_t bits_consumed = 0;
        uint32_t grain_id = 0;
        while (bits_consumed < bits_to_hash) {
            uint64_t grain = getLsb_(bitvector, hash_len);
            rotateRight_(grain, hash_len, rotate_amount + grain_id);
            bits_consumed += hash_len;
            ++ grain_id;
            bitvector >>= hash_len;
            if (bits_consumed > bits_to_hash) {
                //sparta_assert((bits_to_hash + hash_len) > bits_consumed,
                //              "Too many bits consumed=" <<  bits_consumed);
                // mask off the extra MSBs for the last grain
                grain &= ((0x1ULL << (bits_to_hash + hash_len - bits_consumed)) - 1);
            }
            hash ^= grain;
        }

        hash &= hash_mask;
        return hash;
    }

    template<typename PredT>
    uint32_t TaggedTable<PredT>::calcIdx_mallard_(uint64_t PC,
                                                  const GHR  &ghr,
                                                  uint32_t rotate_amount)  const
    {
        return hashPcAndGhr_(PC, ghr, tbl_idx_num_bits_, rotate_amount);
    }

    template<typename PredT>
    uint64_t TaggedTable<PredT>::calcTag_mallard_(uint64_t PC,
                                                  const GHR &ghr,
                                                  uint32_t rotate_amount) const
    {
        return hashPcAndGhr_(PC, ghr, tbl_tag_num_bits_, rotate_amount);
    }

    // Index computation--ful hash of PC, ghist and phist
    template<typename PredT>
    uint32_t TaggedTable<PredT>::calcIdx_classic_(uint64_t PC)  const
    {
        const uint64_t hashed_PC  = PC ^ (PC >> (abs (int32_t(tbl_idx_num_bits_) - int32_t(bank_num_)) + 1));
        const uint32_t index      = hashed_PC ^ folded_ghr_for_idx_.getValue();

        return (index & tbl_idx_mask_);
    }

    // Index computation--ful hash of PC, ghist and phist
    template<typename PredT>
    uint32_t TaggedTable<PredT>::calcIdx_cbp_(uint64_t PC, const PHR &phr)  const
    {

	    int index;
	    int M = (geometric_ghr_len_ > phr.getLength()) ? phr.getLength() : geometric_ghr_len_;
        index =
            PC ^ (PC >> (abs (int32_t(tbl_idx_num_bits_) - int32_t(bank_num_+1)) + 1))
            ^ folded_ghr_for_idx_.getValue() ^ F (phr.getHistory(), M, bank_num_+1);

        return (index & tbl_idx_mask_);
    }

    // Tag computation
    // Two slightly different folded GHRs are used to avoid periodic pattern
    // in the global history (see PPM paper)
    template<typename PredT>
    uint64_t TaggedTable<PredT>::calcTag_classic_(uint64_t PC) const
    {
        const uint64_t tag = PC ^ folded_ghr_for_tag_0_.getValue() ^ (folded_ghr_for_tag_1_.getValue() << 1);
        return (tag & tag_mask_);
    }
}; // namespace tage
}; // namespace dabble
