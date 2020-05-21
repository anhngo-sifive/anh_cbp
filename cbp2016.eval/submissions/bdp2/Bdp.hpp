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
            tage::Tage<bool>(bdp_params)
        {
            const bool     bimodal_init_pred = false;
            const uint32_t bimodal_init_hys = 1;
            initializeBimodalTable(bimodal_init_pred, bimodal_init_hys);

            const bool     tagged_init_pred = true;
            const uint32_t tagged_init_hys = 0;
            initializeTaggedTables(tagged_init_pred, tagged_init_hys);
        }
    }; // class Bdp

}; // namespace bdp
}; // namespace dabble
