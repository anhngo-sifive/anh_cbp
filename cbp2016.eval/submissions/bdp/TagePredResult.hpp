#pragma once
#include <type_traits>
#include <cinttypes>
#include <cassert>
#include <iostream>
#include <vector>

namespace dabble {
namespace tage {
    class IdxTagInfo {
    public:
        void setIdxAndTag(int32_t idx, uint64_t tag) { idx_=idx; tag_=tag; }
        bool isValid() const { return (idx_>=0); }
        int32_t getIdx() const { assert(isValid()); return idx_; }
        uint64_t getTag() const { assert(isValid()); return tag_; }
    private:
        int32_t idx_  = -1;
        uint64_t tag_ = 0;

    };

    template<typename PredT>
    class PredResult
    {
    public:
        explicit PredResult(const uint32_t num_tagged_tables);
        void saveIdxAndTag(int32_t bank, int32_t idx, uint64_t tag);
        void setLongestMatchPred(int32_t bank, PredT pred,
                                 bool is_weak);
        void setAltPred(int32_t bank, PredT pred);
        void setBimodalPred(int32_t idx, PredT pred);
        void summarizePred(bool use_alt);

        int32_t  getSavedIdx(int32_t bank) const { return saved_idx_tag_.at(bank).getIdx(); }
        uint64_t getSavedTag(int32_t bank) const { return saved_idx_tag_.at(bank).getTag(); }

        bool     isLongestMatchValid() const { return longest_match_bank_>=0; }
        int32_t  getLongestMatchBank() const { return longest_match_bank_; }
        PredT    getLongestMatchPred() const { return longest_match_pred_; }
        bool     isLongestMatchWeak()  const { return longest_match_weak_; }

        bool     isAltValid() const { return alt_bank_>=0; }
        int32_t  getAltBank() const { return alt_bank_; }
        PredT    getAltPred() const { return alt_pred_; }

        bool     isValid() const { return isBimodalValid(); }    // does PredResult contain valid prediction
        bool     isBimodalValid() const { return bimodal_idx_>=0; }
        PredT    getBimodalPred() const { return bimodal_pred_; }
        int32_t  getBimodalIdx() const { return bimodal_idx_; }

        // Bank==-1 means prediction came from bimodal
        PredT    getSummarizedPred() const { return summarized_pred_; }
        int32_t  getSummarizedBank() const { return summarized_bank_; }

        uint64_t getPC() const { assert(0); return 0;}

        friend std::ostream& operator<<(std::ostream& os, const PredResult& presult) {
            os << "[" << presult.longest_match_bank_ << "," << presult.longest_match_pred_
               << "|" << presult.alt_bank_ << "," << presult.alt_pred_
               << "|" << presult.bimodal_idx_ << "," << presult.bimodal_pred_
               << "]";
            return os;
        }

    private:
        void setSummarizedPred_(int32_t bank, PredT pred);
        std::vector<IdxTagInfo> saved_idx_tag_;  // saved idx and tag; per table (...the ones that were looked up)

        int32_t  longest_match_bank_ = -1;
        PredT    longest_match_pred_ = 0;
        bool     longest_match_weak_ = false;
        int32_t  alt_bank_ = -1;
        PredT    alt_pred_ = 0;
        int32_t  bimodal_idx_  = -1;
        PredT    bimodal_pred_ = 0;

        // Summarized pred
        int32_t  summarized_bank_ = -1;
        PredT    summarized_pred_ = 0;
    }; // class PredResult

    ////////////////////////////////////////////////////

template<typename PredT>
PredResult<PredT>::PredResult(const uint32_t num_tagged_tables)
{
    saved_idx_tag_.resize(num_tagged_tables);
}

template<typename PredT>
void PredResult<PredT>::saveIdxAndTag(int32_t bank, int32_t idx, uint64_t tag)
{
    assert(!saved_idx_tag_.at(bank).isValid());
    saved_idx_tag_.at(bank).setIdxAndTag(idx, tag);
}


template<typename PredT>
void PredResult<PredT>::setLongestMatchPred(int32_t bank, PredT pred, bool is_weak)
{
    assert(bank>=0);
    longest_match_bank_ = bank;
    longest_match_pred_ = pred;
    longest_match_weak_ = is_weak;
}

template<typename PredT>
void PredResult<PredT>::setAltPred(int32_t bank, PredT pred)
{
    assert(bank>=0);
    alt_bank_ = bank;
    alt_pred_ = pred;
}

template<typename PredT>
void PredResult<PredT>::setBimodalPred(int32_t idx, PredT pred)
{
    assert(idx>=0);
    bimodal_idx_ = idx;
    bimodal_pred_= pred;
}

template<typename PredT>
void PredResult<PredT>::summarizePred(bool use_alt)
{

    /* Prediction Computation:
       1. Use the component with the longest matching history
       2. if prediction_counter==weak && use_alt>=0, use alternate pred
    */
    if (isLongestMatchValid()) {
        if (isLongestMatchWeak() && use_alt) {
            if (isAltValid()) {
                setSummarizedPred_(alt_bank_, alt_pred_);
            }
            else {
                setSummarizedPred_(-1, bimodal_pred_);
            }
        }
        else {
            setSummarizedPred_(longest_match_bank_, longest_match_pred_);
        }
    }
    else {
        if (isBimodalValid()) {
            setSummarizedPred_(-1, bimodal_pred_);
        }
    }

    // default to bimodal if pred is invalid
    if (!isLongestMatchValid()) {
        longest_match_pred_ = bimodal_pred_;
    }
    if (!isAltValid()) {
        alt_pred_ = bimodal_pred_;
    }
}

template<typename PredT>
void PredResult<PredT>::setSummarizedPred_(int32_t bank, PredT pred)
{
    summarized_bank_ = bank;
    summarized_pred_ = pred;
}

}; // namespace tage
}; // namespace dabble
