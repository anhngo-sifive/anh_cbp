#pragma once

#include "Tage.hpp"
#include "IjtpParams.hpp"

namespace dabble {

namespace ijtp
{
    class Ijtp : public tage::Tage<uint64_t>
    {
    public:
        explicit Ijtp(const IjtpParams &ijtp_params) :
            tage::Tage<uint64_t>(ijtp_params, initial_pred_)
        {}

    private:
        static const uint64_t initial_pred_ = 0;
    }; // class Ijtp

}; // namespace ijtp
}; // namespace dabble
