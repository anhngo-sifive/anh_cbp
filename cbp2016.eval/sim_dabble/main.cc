
///////////////////////////////////////////////////////////////////////
//  Copyright 2015 Samsung Austin Semiconductor, LLC.                //
///////////////////////////////////////////////////////////////////////

//Description : Main file for CBP2016

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <map>
using namespace std;

#include "utils.h"
//#include "bt9.h"
#include "bt9_reader.h"
//#include "predictor.cc"
#include "predictor.h"
#include "../submissions/bdp/Bdp.hpp"
#include "../submissions/bdp/GHR.hpp"
#include "../submissions/bdp/BdpPredResult.hpp"
using namespace dabble;


#define COUNTER     unsigned long long


using vuint32 = std::vector<uint32_t>;
uint32_t bdp_use_alt_ctr_num_bits =                               4;
uint32_t bdp_use_alt_num_ctr =                                    0;
uint32_t bdp_num_to_allocate_on_mispredict =                      1;
uint32_t bdp_useful_reset_threshold =                          1023;
uint32_t bdp_tagged_pred_hys_num_bits =                           1;
uint32_t bdp_useful_ctr_num_bits =                                1;
vuint32  bdp_tagged_table_size =       vuint32({512,1024,1024,512});
vuint32  bdp_tag_num_bits =                  vuint32({14,14,14,14});
vuint32  bdp_geometric_len =                 vuint32({23,47,71,97});
uint32_t bdp_bimodal_tbl_size =                                4096;
uint32_t bdp_bimodal_hys_num_bits =                               1;
uint32_t bdp_bimodal_hys_frac =                                   4;
uint32_t addr_num_bits =                                         20;

uint32_t ghr_pc_bits_per_grain =                                 10;
uint32_t ghr_tgt_bits_per_grain =                                10;
uint32_t ghr_hist_bits_per_br =                                   7;
bool     ghr_enable_mallard_hash =                             false;
bool     ghr_dont_update_on_NT =                               false;

uint32_t bdp_num_tagged_tables = bdp_tagged_table_size.size();
uint32_t addr_shift_amount = 1;

void CheckHeartBeat(UINT64 numIter, UINT64 numMispred)
{
  UINT64 dotInterval=1000000;
  UINT64 lineInterval=30*dotInterval;

 UINT64 d1K   =1000;
 UINT64 d10K  =10000;
 UINT64 d100K =100000;
 UINT64 d1M   =1000000;
 UINT64 d10M  =10000000;
 UINT64 d30M  =30000000;
 UINT64 d60M  =60000000;
 UINT64 d100M =100000000;
 UINT64 d300M =300000000;
 UINT64 d600M =600000000;
 UINT64 d1B   =1000000000;
 UINT64 d10B  =10000000000;


//  if(numIter % lineInterval == 0){ //prints line every 30 million branches
//    printf("\n");
//    fflush(stdout);
//  }
  if(numIter == d1K){ //prints MPKI after 100K branches
    printf("  MPKBr_1K         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d10K){ //prints MPKI after 100K branches
    printf("  MPKBr_10K         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d100K){ //prints MPKI after 100K branches
    printf("  MPKBr_100K         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }
  if(numIter == d1M){
    printf("  MPKBr_1M         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d10M){ //prints MPKI after 100K branches
    printf("  MPKBr_10M         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d30M){ //prints MPKI after 100K branches
    printf("  MPKBr_30M         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d60M){ //prints MPKI after 100K branches
    printf("  MPKBr_60M         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d100M){ //prints MPKI after 100K branches
    printf("  MPKBr_100M         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d300M){ //prints MPKI after 100K branches
    printf("  MPKBr_300M         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d600M){ //prints MPKI after 100K branches
    printf("  MPKBr_600M         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d1B){ //prints MPKI after 100K branches
    printf("  MPKBr_1B         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

  if(numIter == d10B){ //prints MPKI after 100K branches
    printf("  MPKBr_10B         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(numIter));
    fflush(stdout);
  }

}//void CheckHeartBeat

OpType BrClass2OpType(bt9::BrClass br_class)
{
          OpType opType = OPTYPE_ERROR;

          if (br_class.type == bt9::BrClass::Type::RET) {
            if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
              opType = OPTYPE_RET_COND;
            else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
              opType = OPTYPE_RET_UNCOND;
            else {
              opType = OPTYPE_ERROR;
            }
          }
          else if (br_class.directness == bt9::BrClass::Directness::INDIRECT) {
            if (br_class.type == bt9::BrClass::Type::CALL) {
              if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
                opType = OPTYPE_CALL_INDIRECT_COND;
              else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
                opType = OPTYPE_CALL_INDIRECT_UNCOND;
              else {
                opType = OPTYPE_ERROR;
              }
            }
            else if (br_class.type == bt9::BrClass::Type::JMP) {
              if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
                opType = OPTYPE_JMP_INDIRECT_COND;
              else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
                opType = OPTYPE_JMP_INDIRECT_UNCOND;
              else {
                opType = OPTYPE_ERROR;
              }
            }
            else {
              opType = OPTYPE_ERROR;
            }
          }
          else if (br_class.directness == bt9::BrClass::Directness::DIRECT) {
            if (br_class.type == bt9::BrClass::Type::CALL) {
              if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL) {
                opType = OPTYPE_CALL_DIRECT_COND;
              }
              else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL) {
                opType = OPTYPE_CALL_DIRECT_UNCOND;
              }
              else {
                opType = OPTYPE_ERROR;
              }
            }
            else if (br_class.type == bt9::BrClass::Type::JMP) {
              if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL) {
                opType = OPTYPE_JMP_DIRECT_COND;
              }
              else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL) {
                opType = OPTYPE_JMP_DIRECT_UNCOND;
              }
              else {
                opType = OPTYPE_ERROR;
              }
            }
            else {
              opType = OPTYPE_ERROR;
            }
          }
          else {
            opType = OPTYPE_ERROR;
          }
          return opType;
}


void updateGhrForConditional(const uint64_t pc,
                             const uint64_t tgt,
                             const bool     taken,
                             tage::GHR      &ghr,
                             bdp::Bdp       &bdp)

{
    // NOTE: Mallard GHR update scheme only updates on a taken branch/jump

    if(!taken && ghr_dont_update_on_NT) {
        return;
    }

    if (ghr_enable_mallard_hash) {
        ghr.updateWithPcAndTarget(pc,tgt);
    }
    else {
        ghr.updateWithTaken(taken);
        bdp.updateFoldedHist(ghr);
#if 1
        static int cnt=0;
        ++cnt;
        std::cout << std::dec << cnt << ": pc=" << std::hex << (pc<<1) << " ghist=" <<  ghr.to_string() << std::endl;
        bdp.printFoldedHist();
#endif
    }
}


// usage: predictor <trace>

int main(int argc, char* argv[]){

  if (argc != 2) {
    printf("usage: %s <trace>\n", argv[0]);
    exit(-1);
  }

  ///////////////////////////////////////////////
  // Init variables
  ///////////////////////////////////////////////

  tage::GHR ghr(ghr_pc_bits_per_grain,
                ghr_tgt_bits_per_grain,
                ghr_hist_bits_per_br);
  bdp::BdpParams bdp_params(bdp_use_alt_ctr_num_bits,
                            bdp_use_alt_num_ctr,
                            bdp_num_to_allocate_on_mispredict,
                            bdp_useful_reset_threshold,
                            bdp_tagged_table_size,
                            bdp_tag_num_bits,
                            bdp_geometric_len,
                            bdp_tagged_pred_hys_num_bits,
                            bdp_useful_ctr_num_bits,
                            bdp_bimodal_tbl_size,
                            bdp_bimodal_hys_num_bits,
                            bdp_bimodal_hys_frac,
                            addr_num_bits);
  bdp::Bdp bdp(bdp_params);


    PREDICTOR  *brpred = new PREDICTOR();  // this instantiates the predictor code
  ///////////////////////////////////////////////
  // read each trace recrod, simulate until done
  ///////////////////////////////////////////////

    std::string trace_path;
    trace_path = argv[1];
    bt9::BT9Reader bt9_reader(trace_path);

    std::string key = "total_instruction_count:";
    std::string value;
    bt9_reader.header.getFieldValueStr(key, value);
    UINT64     total_instruction_counter = std::stoull(value, nullptr, 0);
    UINT64 current_instruction_counter = 0;
    key = "branch_instruction_count:";
    bt9_reader.header.getFieldValueStr(key, value);
    UINT64     branch_instruction_counter = std::stoull(value, nullptr, 0);
    UINT64     numMispred =0;
    UINT64 cond_branch_instruction_counter=0;
    UINT64 uncond_branch_instruction_counter=0;
    UINT64 pred_taken=0;
    UINT64 pred_ntaken=0;

  ///////////////////////////////////////////////
  // read each trace record, simulate until done
  ///////////////////////////////////////////////

      OpType opType= OPTYPE_ERROR;
      UINT64 pc;
      bool actual_taken;
      UINT64 actual_target;
      UINT64 numIter = 0;

      for (auto it = bt9_reader.begin(); it != bt9_reader.end(); ++it) {
          if (it->getSrcNode()->brNodeIndex() == 0) {
              // skip first node
              continue;
          }

          // CheckHeartBeat(++numIter, numMispred); //Here numIter will be equal to number of branches read
        try {
          bt9::BrClass br_class = it->getSrcNode()->brClass();
          opType = BrClass2OpType(br_class);
          assert(opType != OPTYPE_ERROR);

          pc = it->getSrcNode()->brVirtualAddr();

          actual_taken = it->getEdge()->isTakenPath();
          actual_target = it->getEdge()->brVirtualTarget();

          if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL) {
              bool predDir = false;
#if 1
              bdp::PredResult bdp_presult{bdp_num_tagged_tables};
              const bool is_indirect = false;
              const auto shifted_pc  = (pc >> addr_shift_amount);
              bdp.lookupPrediction(pc,
                                   ghr,
                                   is_indirect,
                                   ghr_enable_mallard_hash,
                                   bdp_presult);   // <-- bdp.1:  lookup

              updateGhrForConditional(shifted_pc,
                                      (actual_target>>addr_shift_amount),
                                      actual_taken,
                                      ghr,
                                      bdp); // bdp.2 update GHR

              predDir = bdp_presult.getSummarizedPred();
              bdp.updatePredictor(shifted_pc,
                                  bdp_presult,
                                  actual_taken); // <-- bdp.3:  train tables
#else
              predDir = brpred->GetPrediction(pc);
              brpred->UpdatePredictor(pc, opType, actual_taken, predDir, actual_target);
#endif

              if(predDir != actual_taken){
                  numMispred++; // update mispred stats
              }
              cond_branch_instruction_counter++;

              if (predDir) {
                  ++ pred_taken;
              }
              else {
                  ++ pred_ntaken;
              }
          }
          else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL) { // for predictors that want to track unconditional branches
              uncond_branch_instruction_counter++;
              //brpred->TrackOtherInst(pc, opType, actual_taken, actual_target);
          }
          else {
              fprintf(stderr, "CONDITIONALITY ERROR\n");
              printf("CONDITIONALITY ERROR\n");
              exit(-1); //this should never happen, if it does please email CBP org chair.
          }

/************************************************************************************************************/
        }
        catch (const std::out_of_range & ex) {
          std::cout << ex.what() << '\n';
          break;
        }

      } //for (auto it = bt9_reader.begin(); it != bt9_reader.end(); ++it)


    ///////////////////////////////////////////
    //print_stats
    ///////////////////////////////////////////

    //NOTE: competitors are judged solely on MISPRED_PER_1K_INST. The additional stats are just for tuning your predictors.

      printf("  TRACE \t : %s\n" , trace_path.c_str());
      printf("  NUM_INSTRUCTIONS            \t : %10llu\n",   total_instruction_counter);
      printf("  NUM_BR                      \t : %10llu\n",   branch_instruction_counter-1); //JD2_2_2016 NOTE there is a dummy branch at the beginning of the trace...
      printf("  NUM_PRED_TAKEN              \t : %10llu\n",   pred_taken);
      printf("  NUM_PRED_NTAKEN             \t : %10llu\n",   pred_ntaken);
      printf("  NUM_UNCOND_BR               \t : %10llu\n",   uncond_branch_instruction_counter);
      printf("  NUM_CONDITIONAL_BR          \t : %10llu\n",   cond_branch_instruction_counter);
      printf("  NUM_MISPREDICTIONS          \t : %10llu\n",   numMispred);
      printf("  MISPRED_PER_1K_INST         \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(total_instruction_counter));
      printf("  MPKB                        \t : %10.4f\n",   1000.0*(double)(numMispred)/(double)(cond_branch_instruction_counter));
      printf("\n");
}
