#pragma once

#include <cinttypes>
#include <memory>
//#include "sparta/utils/SpartaAssert.hpp"
//#include "sparta/utils/SpartaSharedPointer.hpp"
#include "SignedCounter.hpp"
#include "GHR.hpp"
#include "TageParams.hpp"
#include "BimodalTable.hpp"
#include "TaggedTable.hpp"
#include "TagePredResult.hpp"

namespace dabble {
namespace tage {

// Branch Direction Predictor
template <typename PredT>
class Tage
{
public:
    Tage(const TageParams &params,
         PredT init_pred);
    bool lookupPrediction(const uint64_t pc,
                          const GHR &ghr,
                          const PHR &phr,
                          const bool is_indirect,
                          const bool enable_mallard_hash,
                          PredResult<PredT> &presult);
    void updatePredictor(const uint64_t pc,
                         const PredResult<PredT> &presult,
                         const PredT actual_val);

    uint32_t getNumTaggedTables() const { return num_tagged_tables_; }

    void setDbgOstream(std::ostream &os);
    void updateFoldedHist(const tage::GHR &ghr);
    void printFoldedHist() const {
        for(uint32_t bank=0; bank<num_tagged_tables_; ++bank) {
            tagged_tbl_.at(bank)->printFoldedHist();
        }
    }
    Tage(const Tage&) = delete;
    Tage& operator=(const Tage&) = delete;

protected:
    int32_t calcUseAltIdx_(uint64_t pc);

    void updateUseAlt_(uint64_t pc,
                       PredT longest_match_pred,
                       PredT alt_pred,
                       PredT actual_val);

    uint32_t seed_=1;
    uint32_t myrand_();

    bool determineNeedToAllocate_(const uint64_t pc,
                                  const PredResult<PredT> &presult,
                                  const PredT actual_val);

    // Allocate up to N entries per mispredict
    void handleAllocate_(const uint64_t pc, const PredResult<PredT> &presult, const PredT actual_val);

    void handleUpdatePrediction_(const PredResult<PredT> &presult, const PredT actual_val);

    std::ostream *dbg_ostream_=nullptr; // Ostream for logging

    // BDP parameters
    const uint32_t num_tagged_tables_;
    const uint32_t num_to_allocate_on_mispredict_;
    const uint32_t useful_reset_threshold_;  // threshold of when to decrement useful counter
    const bool     enable_use_alt_;
    const uint32_t use_alt_num_ctr_mask_;

    using BimodalTablePtr = std::unique_ptr<BimodalTable<PredT>>;
    using TaggedTablePtr  = std::unique_ptr<TaggedTable<PredT>>;
    BimodalTablePtr             bimodal_tbl_;
    std::vector<TaggedTablePtr> tagged_tbl_;
    std::vector<SignedCounter>  use_alt_on_na_; /// Counters for use alt_pred on newly allocate

    int32_t cumulative_no_allocate_ = 0;       /// Running count of no-allocates

}; // class Tage


//////////////////////////////////////////////////////////////////////////////////




template <typename PredT>
Tage<PredT>::Tage(const TageParams &params,
                  const PredT init_pred) :
    num_tagged_tables_(params.tagged_table_size.size()),
    num_to_allocate_on_mispredict_(params.num_to_allocate_on_mispredict),
    useful_reset_threshold_(params.useful_reset_threshold),
    enable_use_alt_(params.use_alt_num_ctr>0),
    use_alt_num_ctr_mask_(params.use_alt_num_ctr-1)
{
    //sparta_assert(num_tagged_tables_ == params.tagged_table_size.size(), "Number of tables and tabe_size array must match");
    //sparta_assert(num_tagged_tables_ == params.tag_num_bits.size(), "Number of tables and tag_num_bits array must match");
    //sparta_assert(num_tagged_tables_ == params.geometric_len.size(), "Number of tables and geometric_len array must match");

    use_alt_on_na_.resize(params.use_alt_num_ctr,
                          SignedCounter(params.use_alt_ctr_num_bits));

    bimodal_tbl_.reset(new BimodalTable<PredT>(params.bimodal_tbl_size,
                                               params.bimodal_hys_num_bits,
                                               params.bimodal_hys_frac));
    bimodal_tbl_->initialize(init_pred);

    const uint64_t init_tag=0;
    tagged_tbl_.resize(num_tagged_tables_);
    for (uint32_t bank=0; bank<num_tagged_tables_; ++bank) {
        tagged_tbl_.at(bank).reset(new TaggedTable<PredT>(bank,
                                                          params.tagged_table_size.at(bank),
                                                          params.tag_num_bits.at(bank),
                                                          params.tagged_pred_hys_num_bits,
                                                          params.useful_ctr_num_bits,
                                                          params.geometric_len.at(bank),
                                                          params.addr_num_bits));
        tagged_tbl_.at(bank)->initialize(init_tag, init_pred);
    }
}

template <typename PredT>
bool Tage<PredT>::lookupPrediction(const uint64_t pc,
                                   const GHR &ghr,
                                   const PHR &phr,
                                   const bool is_indirect,
                                   const bool enable_mallard_hash,
                                   PredResult<PredT> &presult)
{

    TaggedTablePredResult<PredT> tagged_presult;
    int32_t longest_match_bank = -1;

#if 1
    for (int32_t bank=0; bank<num_tagged_tables_ ; ++bank) {
        tagged_tbl_[bank]->printIdxTag(pc, phr);
    }
#endif

    // Longest history tagged-table look-up
    for (int32_t bank=(num_tagged_tables_ - 1); bank>=0; --bank) {
        const bool tag_matched  = tagged_tbl_[bank]->lookupPrediction(pc,
                                                                      ghr,
                                                                      phr,
                                                                      is_indirect,
                                                                      enable_mallard_hash,
                                                                      tagged_presult);
        presult.saveIdxAndTag(bank,
                              tagged_presult.getComputedIdx(),
                              tagged_presult.getComputedTag());  // save the calculated tag and idx
        if (tag_matched){
            longest_match_bank = bank;
            presult.setLongestMatchPred(bank,
                                        tagged_presult.getPred(),
                                        tagged_presult.getWeakPred());
            break;
        }
    }

    // Alternate tagged-table lookup
    for (int32_t bank=(longest_match_bank - 1); bank>=0; --bank) {
        const bool tag_matched  = tagged_tbl_[bank]->lookupPrediction(pc,
                                                                      ghr,
                                                                      phr,
                                                                      is_indirect,
                                                                      enable_mallard_hash,
                                                                      tagged_presult);
        presult.saveIdxAndTag(bank,
                              tagged_presult.getComputedIdx(),
                              tagged_presult.getComputedTag());  // save the calculated tag and idx
        if (tag_matched) {
            presult.setAltPred(bank,
                               tagged_presult.getPred());
            break;
        }
    }

    // Bimodal lookup
    int32_t  bimodal_idx=-1;
    PredT    bimodal_pred=0;
    bimodal_tbl_->lookupPrediction(pc, bimodal_pred, bimodal_idx);
    presult.setBimodalPred(bimodal_idx, bimodal_pred);

    // Longest-match, alt, and bimodal preds have been saved in presult.
    // Now summarize it.
    // NOTE:  use-alt counter is a heuristic to indicate if
    //        alt is better than longest-match
    bool use_alt = false;
    if (enable_use_alt_) {
        const int32_t use_alt_idx = calcUseAltIdx_(pc);
        use_alt = (use_alt_on_na_.at(use_alt_idx).isGTEZero());
    }
    presult.summarizePred(use_alt);

    return true; // Tage lookup is always valid
}


template <typename PredT>
void Tage<PredT>::updatePredictor(const uint64_t pc,
                                  const PredResult<PredT> &presult,
                                  PredT actual_val)
{
    bool need_allocate = determineNeedToAllocate_(pc, presult, actual_val);

    if (need_allocate) {
        handleAllocate_(pc, presult, actual_val);
    }

    handleUpdatePrediction_(presult, actual_val);

}

template <typename PredT>
void Tage<PredT>::updateUseAlt_(uint64_t pc,
                                PredT longest_match_pred,
                                PredT alt_pred,
                                PredT actual_val)
{
    if (longest_match_pred != alt_pred) {
        const int32_t use_alt_idx = calcUseAltIdx_(pc);
        const bool alt_is_correct = (alt_pred == actual_val);
        use_alt_on_na_.at(use_alt_idx).update(alt_is_correct);
    }
}


template <typename PredT>
uint32_t Tage<PredT>::myrand_()
{
    // Placeholder to see if/how Mallard implements it
    // A simple random generator from the web.
    seed_ = (seed_ * 32719 + 3) % 32749;
    return seed_;
}

template <typename PredT>
int32_t Tage<PredT>::calcUseAltIdx_(uint64_t pc)
{
    return ((pc^(pc>>2)) & use_alt_num_ctr_mask_);
}

template <typename PredT>
bool Tage<PredT>::determineNeedToAllocate_(const uint64_t pc,
                                           const PredResult<PredT> &presult,
                                           const PredT actual_val)
{
    const bool  mispredicted       = (presult.getSummarizedPred() != actual_val);
    const bool  longest_match_valid= presult.isLongestMatchValid();
    const auto  longest_match_bank = presult.getLongestMatchBank();
    const PredT longest_match_pred = presult.getLongestMatchPred();
    const PredT alt_pred           = presult.getAltPred();
    const bool  max_bank_used      = longest_match_bank >= (int32_t(num_tagged_tables_) - 1);
    bool need_allocate = mispredicted && !max_bank_used;



    // update Use-Alt
    if (longest_match_valid) {
        // Manage the selection between longest matching and alternate matching
        //    for "pseudo"-newly allocated longest matching entry.
        // An entry is considered as newly allocated if its prediction counter is weak.
        const bool pseudo_new_alloc = presult.isLongestMatchWeak();

        if ( pseudo_new_alloc ) {
            //If longest-match was delivering the correct prediction, no need to allocate a new entry
            //even if the overall prediction was false.
            if ( longest_match_pred == actual_val ) {
                need_allocate = false;
            }

            if (enable_use_alt_) {
                updateUseAlt_(pc, longest_match_pred, alt_pred, actual_val);
            }
        }
    }

    return need_allocate;
}

// Allocate up to N entries per mispredict
template <typename PredT>
void Tage<PredT>::handleAllocate_(const uint64_t pc, const PredResult<PredT> &presult, const PredT actual_val)
{
    const auto longest_match_bank = presult.getLongestMatchBank();
    uint32_t num_no_allocate = 0;
    uint32_t num_allocate    = 0;
    // TAGE can randomly skip the table after longest-match when allocating
    // NOTE: removing this feature dropped MPKI for a few CBP eval traces.
    // XXX check to see how Mallard does this
    const uint32_t start_bank = (longest_match_bank + 1);
    for (uint32_t bank=start_bank; bank<num_tagged_tables_; ++bank) {
        const auto saved_idx = presult.getSavedIdx(bank);
        auto &tentry = tagged_tbl_.at(bank)->getEntry(saved_idx);
        if (tentry.getUsefulCounter().getValue() == 0) {
            tentry.reset(presult.getSavedTag(bank), actual_val);
            std::cout << "Allocate, pc=" << std::hex << (pc<<1) << std::dec
                      << " bank=" << (bank+1)
                      << " idx=" << saved_idx
                      << " tag=" << std::hex << presult.getSavedTag(bank)
                      << " actual=" << actual_val
                      << std::endl;
            ++num_allocate;
            if (num_allocate >= num_to_allocate_on_mispredict_) {
                break;
            }
            ++bank; // don't allocate on two adjacent banks
        }
        else {
            ++num_no_allocate;
        }
    }

    // Gracefully decrement useful counters based on cumulative_no_allocate
    cumulative_no_allocate_ += (num_no_allocate - num_allocate );
    if (cumulative_no_allocate_ < 0) {
        cumulative_no_allocate_ = 0;
    }
    if (cumulative_no_allocate_ > int32_t(useful_reset_threshold_)) {
        for (uint32_t bank=0; bank<num_tagged_tables_; ++bank) {
            tagged_tbl_.at(bank)->gracefulDecrementAllUsefulCounters();
        }
        cumulative_no_allocate_ = 0;
    }
}

template <typename PredT>
void Tage<PredT>::handleUpdatePrediction_(const PredResult<PredT> &presult, const PredT actual_val)
{
    const bool  longest_match_valid= presult.isLongestMatchValid();
    const auto  longest_match_bank = presult.getLongestMatchBank();
    const PredT longest_match_pred = presult.getLongestMatchPred();
    const bool  alt_valid = presult.isAltValid();
    const PredT alt_pred  = presult.getAltPred();
    const auto  alt_bank  = presult.getAltBank();

    //update prediction counter
    if (longest_match_valid) {

        // 1. Update alternate:  if longest-match is wrong and weak
        const bool update_alt =  presult.isLongestMatchWeak() &&
                                 (longest_match_pred != actual_val);
        if ( update_alt ) {
            if (alt_valid) {
                const auto alt_idx   = presult.getSavedIdx(alt_bank);
                auto &alt_tentry = tagged_tbl_.at(alt_bank)->getEntry(alt_idx);
                alt_tentry.updatePredCounter(actual_val);
                std::cout << " Alt ctrupdate bank=" << std::dec << (alt_bank+1) << " idx=" << alt_idx << std::endl;
            }
            else {
                bimodal_tbl_->updatePredCounter(presult.getBimodalIdx(), actual_val);
            }
        }

        // 2. Update longest-matching
        const auto longest_match_idx = presult.getSavedIdx(longest_match_bank);
        auto &longest_match_tentry = tagged_tbl_.at(longest_match_bank)->getEntry(longest_match_idx);
        longest_match_tentry.updatePredCounter(actual_val);
        std::cout << " Longest ctrupdate bank="<<  std::dec << (longest_match_bank+1) << " idx=" << longest_match_idx << std::endl;

        // 3. Reset 'useful' if hysteresis becomes weak
        if (  longest_match_tentry.isPredWeak() ) {
            longest_match_tentry.resetUsefulCounter();
        }

    }
    else {
        bimodal_tbl_->updatePredCounter(presult.getBimodalIdx(),
                                        actual_val);
    }

    // Update 'useful': if longest-matching is correct and alt is incorrect
    // NOTE: no need to make sure longest_match_pred and alt_pred are valid
    //       since they default to bimodal
    if ( (longest_match_pred != alt_pred) &&
         (longest_match_pred == actual_val) ) {
        const auto longest_match_idx = presult.getSavedIdx(longest_match_bank);
        auto &longest_match_tentry = tagged_tbl_.at(longest_match_bank)->getEntry(longest_match_idx);
        longest_match_tentry.updateUsefulCounter(true); // increment usefulness
    }
}

template <typename PredT>
void Tage<PredT>::updateFoldedHist(const tage::GHR &ghr)
{
    for(uint32_t bank=0; bank<num_tagged_tables_; ++bank) {
        tagged_tbl_.at(bank)->updateFoldedHist(ghr);
    }
}

template <typename PredT>
void Tage<PredT>::setDbgOstream(std::ostream &os) {
    dbg_ostream_ = &os;
    for (uint32_t bank=0; bank<num_tagged_tables_; ++bank) {
        tagged_tbl_[bank]->setDbgOstream(os);
    }
}

}; // namespace tage
}; // namespace dabble
