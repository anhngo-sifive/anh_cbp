#pragma once

#include "Tage.hpp"
#include "BdpParams.hpp"

namespace dabble {

namespace bdp
{
    class Bdp : public tage::Tage<bool>
    {
    public:
        explicit Bdp(const BdpParams &bdp_params) :
            tage::Tage<bool>(bdp_params, initial_pred_)
        {}

    private:
        static const bool initial_pred_ = true;
    }; // class Bdp

}; // namespace bdp
}; // namespace dabble
