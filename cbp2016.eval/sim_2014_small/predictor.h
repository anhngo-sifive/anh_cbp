
#ifndef _PREDICTOR_H_
#define _PREDICTOR_H_

//To get the predictor storage budget on stderr  uncomment the next line

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <inttypes.h>
#include <cmath>
#include "utils.h"
#include "bt9.h"
#include "bt9_reader.h"

#define NNN 0			// number of extra entries allocated on a TAGE misprediction


#define HYSTSHIFT 0 // bimodal hysteresis shared by 4 entries
#define LOGB 14 // log of number of entries in bimodal predictor

#define PHISTWIDTH 27		// width of the path history used in TAGE

#define UWIDTH 1  // u counter width on TAGE
#define CWIDTH 2  // predictor counter width on the TAGE tagged tables

#define HISTBUFFERLENGTH 4096	// we use a 4K entries history buffer to store the branch history

//the counter(s) to chose between longest match and alternate prediction on TAGE when weak counters
#define LOGSIZEUSEALT 0
#define SIZEUSEALT  (1<<(LOGSIZEUSEALT))
#define INDUSEALT (PC & (SIZEUSEALT -1))
std::vector<std::vector<int8_t>>  use_alt_on_na;
bool enable_use_alt=false;

int8_t BIM=0;

std::string getGhistString(const std::vector<uint8_t> &ghist, int ptr) {
    ptr += 511;
    std::stringstream tmp;
    for (uint32_t i=0; i<512; ++i) {
        ptr = ptr % ghist.size();
        bool hbit = ghist.at(ptr);
        tmp << hbit;
        --ptr;
    }
    return tmp.str();
}

// utility class for index computation
// this is the cyclic shift register for folding
// a long global history into a smaller number of bits; see P. Michaud's PPM-like predictor at CBP-1
class folded_history
{
public:
    unsigned comp;
    int CLENGTH;
    int OLENGTH;
    int OUTPOINT;

    folded_history () {}

    void init (int original_length, int compressed_length, int N)
	{
	    comp = 0;
	    OLENGTH = original_length;
	    CLENGTH = compressed_length;
	    OUTPOINT = OLENGTH % CLENGTH;

	}

    void update (std::vector<uint8_t> &hist, int PT)
	{
	    comp = (comp << 1) ^ hist.at(PT & (HISTBUFFERLENGTH - 1));
	    comp ^= hist.at((PT + OLENGTH) & (HISTBUFFERLENGTH - 1)) << OUTPOINT;
	    comp ^= (comp >> CLENGTH);
	    comp = (comp) & ((1 << CLENGTH) - 1);
	}

};

class bentry			// TAGE bimodal table entry
{
public:
#if 0
    int8_t pred = 0;
    int8_t hyst = 1;
#else
    int8_t pred = 1;
    int8_t hyst = 0;
#endif
};

class gentry			// TAGE global table entry
{
public:
    int8_t ctr = 0;
    uint   tag = 0;
    int8_t u   = 0;
};



#define DBG

int TICK;// for the reset of the u counter

// number of tagged tables
#define NHIST 12

std::vector<uint8_t> ghist(HISTBUFFERLENGTH,0);
int ptghist = 0;
long long phist=0;		//path history
std::vector<folded_history> ch_i(NHIST + 1);	//utility for computing TAGE indices
std::vector<std::vector<folded_history>> ch_t(2,std::vector<folded_history>(NHIST + 1));	//utility for computing TAGE tags


// number of tagged tables
#if 1
  #define LOGB 12 // log of number of entries in bimodal predictor
  #define NHIST 4
  std::vector<int>  m{0,23,47,71,97}; // history lengths
  std::vector<int>  TB{0,14,14,14,14};		// tag width for the different tagged tables
  std::vector<int>  logg={0,9,10,10,9};		// log of number entries of the different tagged tables
#else
  #define LOGB 14 // log of number of entries in bimodal predictor
  #define NHIST 12
  std::vector<int>  m{0,8,12,18,27,40,60,90,135,203,305,459,690}; // history lengths
  std::vector<int>  TB{0,8, 9, 9,10,10,11,11,12,12,13,13,14};		// tag width for the different tagged tables
  std::vector<int>  logg={0,10,10,10,10,10,10,10,10,10,10,10,10};		// log of number entries of the different tagged tables
#endif

//For the TAGE predictor
std::vector<bentry> btable(1 << LOGB);			//bimodal TAGE table
std::vector<std::vector<gentry>> gtable(NHIST + 1);	// tagged TAGE tables
std::vector<int>  GI(NHIST + 1);		// indexes to the different tables are computed only once
std::vector<uint> GTAG(NHIST + 1);	// tags for the different tables are computed only once
int  BI;				// index of the bimodal table
bool pred_taken;		// prediction
bool alttaken;			// alternate  TAGEprediction
bool tage_pred;			// TAGE prediction
bool LongestMatchPred;
int  HitBank;			// longest matching bank
int  AltBank;			// alternate matching bank
int  Seed;		        // for the pseudo-random number generator



class PREDICTOR
{
public:
    PREDICTOR (void)
	{
	    reinit ();
	}


    void reinit ()
	{
	    for (int i = 1; i <= NHIST; i++)  {
            gtable.at(i).resize(1 << (logg.at(i)));
	    }

	    for (int i = 1; i <= NHIST; i++)  {
            ch_i[i].init (m.at(i), (logg.at(i)), i - 1);
            ch_t[0][i].init (ch_i[i].OLENGTH, TB.at(i), i);
            ch_t[1][i].init (ch_i[i].OLENGTH, TB.at(i) - 1, i + 2);
	    }
	    Seed = 0;

	    TICK = 0;
	    phist = 0;
	    Seed = 0;

        use_alt_on_na.resize(SIZEUSEALT);
	    for (int i = 0; i < SIZEUSEALT; i++)  {
            use_alt_on_na.at(i).resize(2);
            use_alt_on_na.at(i).at(0) = 0;
            use_alt_on_na.at(i).at(1) = 0;
	    }

	    TICK = 0;
	    phist = 0;

	}

    // index function for the bimodal table
    int bindex (uint32_t PC)
	{
	    return ((PC) & ((1 << (LOGB)) - 1));
	}


    // the index functions for the tagged tables uses path history as in the OGEHL predictor
    //F serves to mix path history: not very important impact
    int F (long long A, int size, int bank)
	{
	    int A1, A2;
	    A = A & ((1 << size) - 1);
	    A1 = (A & ((1 << logg.at(bank)) - 1));
	    A2 = (A >> logg.at(bank));
	    A2 = ((A2 << bank) & ((1 << logg.at(bank)) - 1)) + (A2 >> (logg.at(bank) - bank));
	    A = A1 ^ A2;
	    A = ((A << bank) & ((1 << logg.at(bank)) - 1)) + (A >> (logg.at(bank) - bank));
	    return (A);
	}

    // gindex computes a full hash of PC, ghist and phist
    int gindex (unsigned int PC, int bank, long long hist,
                std::vector<folded_history> &ch_i)
	{
	    int index;
	    int M = (m.at(bank) > PHISTWIDTH) ? PHISTWIDTH : m.at(bank);
	    index =
            PC ^ (PC >> (abs (logg.at(bank) - bank) + 1))
            ^ ch_i.at(bank).comp ^ F (hist, M, bank);
	    return (index & ((1 << (logg.at(bank))) - 1));
	}

    //  tag computation
    uint16_t gtag (unsigned int PC, int bank, std::vector<folded_history> &ch0,
                   std::vector<folded_history> &ch1)
	{
	    int tag = PC ^ ch0.at(bank).comp ^ (ch1.at(bank).comp << 1);
	    return (tag & ((1 << TB.at(bank)) - 1));
	}

    // up-down saturating counter
    void ctrupdate (int8_t &ctr, bool taken, int nbits)
	{
	    if (taken)  {
            if (ctr < ((1 << (nbits - 1)) - 1))
                ctr++;
	    }
	    else  {
            if (ctr > -(1 << (nbits - 1)))
                ctr--;
	    }
	}


    bool getbim ()
	{
	    BIM = (btable.at(BI).pred << 1) + (btable.at(BI >> HYSTSHIFT).hyst);
	    return (btable.at(BI).pred > 0);
	}

    void baseupdate (bool Taken)
	{
	    int inter = BIM;
	    if (Taken)  {
            if (inter < 3)
                inter += 1;
	    }
	    else if (inter > 0)
            inter--;
#ifdef DBG
        std::cout << " Base ctrupdate idx=" << std::dec << BI << std::endl;
#endif
	    btable.at(BI).pred = inter >> 1;
	    btable.at(BI >> HYSTSHIFT).hyst = (inter & 1);
	};

    //just a simple pseudo random number generator: use available information
    // to allocate entries  in the loop predictor
    int MYRANDOM ()
	{
	    Seed++;
	    Seed ^= phist;
	    Seed = (Seed >> 21) + (Seed << 11);
	    return (Seed);
	};


    //  TAGE PREDICTION: same code at fetch or retire time but the index and tags must recomputed
    void Tagepred (UINT64 PC)
	{
	    HitBank = 0;
	    AltBank = 0;
	    for (int i = 1; i <= NHIST; i++) {
            GI.at(i) = gindex (PC, i, phist, ch_i);
            GTAG.at(i) = gtag (PC, i, ch_t.at(0), ch_t.at(1));
#ifdef DBG
            std::cout << std::dec << " -Pred bank=" << i
                      << " idx=" <<  GI.at(i)
                      << std::hex << " t=" << GTAG.at(i)
                      << std::dec << " pred=" << (gtable.at(i).at(GI.at(i)).ctr>=0)
                      << std::endl;
#endif
	    }

	    BI = PC & ((1 << LOGB) - 1);


        //Look for the bank with longest matching history
	    for (int i = NHIST; i > 0; i--)  {
            if (gtable.at(i).at(GI.at(i)).tag == GTAG.at(i)) {
                HitBank = i;
                LongestMatchPred = (gtable.at(HitBank).at(GI.at(HitBank)).ctr >= 0);
                break;
            }
	    }

        //Look for the alternate bank
	    for (int i = HitBank - 1; i > 0; i--) {
            if (gtable.at(i).at(GI.at(i)).tag == GTAG.at(i)) {
                AltBank = i;
                break;
            }
	    }

        //computes the prediction and the alternate prediction
	    if (HitBank > 0)  {
            if (AltBank > 0)
                alttaken = (gtable.at(AltBank).at(GI.at(AltBank)).ctr >= 0);
            else
                alttaken = getbim ();

            //if the entry is recognized as a newly allocated entry and
            //USE_ALT_ON_NA is positive  use the alternate prediction
            int index = (INDUSEALT ^ LongestMatchPred) & (SIZEUSEALT -1);
            bool Huse_alt_on_na =
                (use_alt_on_na.at(index).at(HitBank > (NHIST / 3)) >= 0);

            if (!enable_use_alt) {
                tage_pred = LongestMatchPred;
            }
            else {
                if ((!Huse_alt_on_na)
                    || (abs (2 * gtable.at(HitBank).at(GI.at(HitBank)).ctr + 1) > 1))
                    tage_pred = LongestMatchPred;
                else
                    tage_pred = alttaken;
            }

	    }
	    else  {
            alttaken = getbim ();
            tage_pred = alttaken;
            LongestMatchPred = alttaken;
	    }



	}
    //compute the prediction

    bool GetPrediction (UINT64 PC)
	{
        // computes the TAGE table addresses and the partial tags
	    Tagepred (PC);
	    pred_taken = tage_pred;
	    return pred_taken;
	}

    void HistoryUpdate (uint64_t PC,
                        uint16_t brtype,
                        bool taken,
                        uint32_t target,
                        long long &phr,
                        int &Y,
                        std::vector<folded_history> &H,
                        std::vector<folded_history> &G,
                        std::vector<folded_history> &J)
	{
        //special treatment for unconditional branchs;
	    int maxt;
	    if (brtype == OPTYPE_RET_COND || brtype == OPTYPE_JMP_DIRECT_COND || brtype == OPTYPE_JMP_INDIRECT_COND || brtype == OPTYPE_CALL_DIRECT_COND || brtype == OPTYPE_CALL_INDIRECT_COND)
            maxt = 1;
	    else
            maxt = 4;

	    int T = ((PC) << 1) + taken;
	    int PATH = PC;

	    for (int t = 0; t < maxt; t++)   {
            bool DIR = (T & 1);
            T >>= 1;
            int PATHBIT = (PATH & 127);
            PATH >>= 1;
            //update  history
            Y--;
            ghist.at(Y & (HISTBUFFERLENGTH - 1)) = DIR;
            phr = (phr << 1) ^ PATHBIT;
            phr &= ((1ULL << PHISTWIDTH) - 1);
#ifdef DBG
            std::cout << " updateHistory: pc=" << std::hex << PC
                      << " phist=" << std::setw(8) << phist
                      << " ghist=" <<  getGhistString(ghist, ptghist)
                      << std::endl;
#endif
            for (int i = 1; i <= NHIST; i++) {
                H.at(i).update (ghist, Y);
                G.at(i).update (ghist, Y);
                J.at(i).update (ghist, Y);
#ifdef DBG
                std::cout << "  -FoldedHist bank=" << std::dec << i
                          << " c_i="  << std::hex << H.at(i).comp
                          << " c_t0=" << G.at(i).comp
                          << " c_t1=" <<  J.at(i).comp
                          << std::endl;

#endif
            }
	    }

        //END UPDATE  HISTORIES
	}

    // PREDICTOR UPDATE
    void UpdatePredictor (UINT64 PC, bool resolveDir, bool predDir,
                          UINT64 branchTarget)
	{
        //TAGE UPDATE
	    bool ALLOC = ((tage_pred != resolveDir) & (HitBank < NHIST));
	    //do not allocate too often if the overall prediction is correct

        // Use-alt update
	    if (HitBank > 0)  {
            // Manage the selection between longest matching and alternate matching
            // for "pseudo"-newly allocated longest matching entry
            bool PseudoNewAlloc =
                (abs (2 * gtable.at(HitBank).at(GI.at(HitBank)).ctr + 1) <= 1);
            // an entry is considered as newly allocated if its prediction counter is weak
            if (PseudoNewAlloc) 	{
                if (LongestMatchPred == resolveDir)
                    ALLOC = false;
                // if it was delivering the correct prediction, no need to allocate a new entry
                //even if the overall prediction was false
                if (LongestMatchPred != alttaken) {
                    int index = (INDUSEALT ^ LongestMatchPred) & (SIZEUSEALT -1);
                    ctrupdate (use_alt_on_na.at(index).at(HitBank > (NHIST / 3)),
                               (alttaken == resolveDir), 4);
                }
            }
	    }


	    if (ALLOC) {
            int T = NNN;
            int A = 1;
            int Penalty = 0;
            int NA = 0;
            for (int i = HitBank + A; i <= NHIST; i += 1)  {
                if (gtable.at(i).at(GI.at(i)).u == 0) {
#ifdef DBG
                    std::cout << " Tagged Allocate, pc=" << std::hex << PC << std::dec
                              << " bank=" << i
                              << " idx=" << GI.at(i)
                              << " tag=" << std::hex <<  GTAG.at(i)
                              << " actual=" << resolveDir
                              << std::endl;
#endif
                    gtable.at(i).at(GI.at(i)).tag = GTAG.at(i);
                    gtable.at(i).at(GI.at(i)).ctr = (resolveDir) ? 0 : -1;
                    NA++;
                    if (T <= 0)  {
                        break;
                    }
                    i +=  1;
                    T -= 1;
                }
                else  {
                    Penalty++;
                }
            }
            TICK += (Penalty -NA );
            //just the best formula for the Championship
            if (TICK < 0)
                TICK = 0;
            if (TICK > 1023) {
                for (int i = 1; i <= NHIST; i++) {
                    for (int j = 0; j <= (1 << logg.at(i)) - 1; j++) {
                        // substracting 1 to a whole array is not that realistic
                        gtable.at(i).at(j).u >>= 1;
                    }
                }
                TICK = 0;
            }


	    }

        //update predictions
	    if (HitBank > 0)   {
            // 1. update alt
            if (abs (2 * gtable.at(HitBank).at(GI.at(HitBank)).ctr + 1) == 1)
                if (LongestMatchPred != resolveDir)  {
                    // acts as a protection
                    if (AltBank > 0) {
#ifdef DBG
                        std::cout << " Alt ctrupdate bank=" << std::dec << AltBank << " idx=" << GI.at(AltBank) << std::endl;
#endif
                        ctrupdate (gtable.at(AltBank).at(GI.at(AltBank)).ctr,
                                   resolveDir, CWIDTH);

                    }
                    if (AltBank == 0)
                        baseupdate (resolveDir);
                }

            // 2. update longest-match
#ifdef DBG
            std::cout << " Longest ctrupdate bank="<<  std::dec << HitBank << " idx=" << GI.at(HitBank) << std::endl;
#endif
            ctrupdate (gtable.at(HitBank).at(GI.at(HitBank)).ctr, resolveDir, CWIDTH);

            // 3. sign changes: no way it can have been useful
            if (abs (2 * gtable.at(HitBank).at(GI.at(HitBank)).ctr + 1) == 1)
                gtable.at(HitBank).at(GI.at(HitBank)).u = 0;

	    }
	    else
            baseupdate (resolveDir);

	    if (LongestMatchPred != alttaken)
            if (LongestMatchPred == resolveDir)  {
                if (gtable.at(HitBank).at(GI.at(HitBank)).u < (1 << UWIDTH) - 1)
                    gtable.at(HitBank).at(GI.at(HitBank)).u++;
            }
        //END TAGE UPDATE

	    HistoryUpdate (PC,
                       OPTYPE_JMP_DIRECT_COND,
                       resolveDir,
                       branchTarget,
                       phist,
                       ptghist,
                       ch_i,
                       ch_t.at(0),
                       ch_t.at(1));

        //END PREDICTOR UPDATE


	}

    void TrackOtherInst (UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget)
	{

	    bool taken = true;

	    switch (opType)
	    {
	    case    OPTYPE_RET_UNCOND:
	    case    OPTYPE_JMP_DIRECT_UNCOND:
	    case    OPTYPE_JMP_INDIRECT_UNCOND:
	    case    OPTYPE_CALL_DIRECT_UNCOND:
	    case    OPTYPE_CALL_INDIRECT_UNCOND:

		HistoryUpdate (PC,
                       opType,
                       taken,
                       branchTarget,
                       phist,
                       ptghist,
                       ch_i,
                       ch_t.at(0),
                       ch_t.at(1));
		break;


	    default:;
	    }


	}

};





#endif
