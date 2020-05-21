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
            tage::Tage<uint64_t>(ijtp_params)
        {
            const uint64_t bimodal_init_pred = 0;
            const uint32_t bimodal_init_hys = 0;
            initializeBimodalTable(bimodal_init_pred, bimodal_init_hys);

            const uint64_t tagged_init_pred = 0;
            const uint32_t tagged_init_hys = 0;
            initializeTaggedTables(tagged_init_pred, tagged_init_hys);
        }
    }; // class Ijtp

}; // namespace ijtp
}; // namespace dabble
