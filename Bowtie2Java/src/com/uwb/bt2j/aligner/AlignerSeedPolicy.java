/* #ifndef ALIGNER_SEED_POLICY_H_
#define ALIGNER_SEED_POLICY_H_

#include "scoring.h"
#include "simple_func.h"

#define DEFAULT_SEEDMMS 0
#define DEFAULT_SEEDLEN 22
#define DEFAULT_LOCAL_SEEDLEN 20

#define DEFAULT_IVAL SIMPLE_FUNC_SQRT
#define DEFAULT_IVAL_A 1.15f
#define DEFAULT_IVAL_B 0.0f

#define DEFAULT_UNGAPPED_HITS 6 */
package com.uwb.bt2j.aligner;
/**
 * Encapsulates the set of all parameters that affect what the
 * SeedAligner does with reads.
 */
class SeedAlignmentPolicy {
	/**
	 * Parse alignment policy when provided in this format:
	 * <lab>=<val>;<lab>=<val>;<lab>=<val>...
	 *
	 * And label=value possibilities are:
	 *
	 * Bonus for a match
	 * -----------------
	 *
	 * MA=xx (default: MA=0, or MA=2 if --local is set)
	 *
	 *    xx = Each position where equal read and reference characters match up
	 *         in the alignment contriubtes this amount to the total score.
	 *
	 * Penalty for a mismatch
	 * ----------------------
	 *
	 * MMP={Cxx|Q|RQ} (default: MMP=C6)
	 *
	 *   Cxx = Each mismatch costs xx.  If MMP=Cxx is specified, quality
	 *         values are ignored when assessing penalities for mismatches.
	 *   Q   = Each mismatch incurs a penalty equal to the mismatched base's
	 *         value.
	 *   R   = Each mismatch incurs a penalty equal to the mismatched base's
	 *         rounded quality value.  Qualities are rounded off to the
	 *         nearest 10, and qualities greater than 30 are rounded to 30.
	 *
	 * Penalty for position with N (in either read or reference)
	 * ---------------------------------------------------------
	 *
	 * NP={Cxx|Q|RQ} (default: NP=C1)
	 *
	 *   Cxx = Each alignment position with an N in either the read or the
	 *         reference costs xx.  If NP=Cxx is specified, quality values are
	 *         ignored when assessing penalities for Ns.
	 *   Q   = Each alignment position with an N in either the read or the
	 *         reference incurs a penalty equal to the read base's quality
	 *         value.
	 *   R   = Each alignment position with an N in either the read or the
	 *         reference incurs a penalty equal to the read base's rounded
	 *         quality value.  Qualities are rounded off to the nearest 10,
	 *         and qualities greater than 30 are rounded to 30.
	 *
	 * Penalty for a read gap
	 * ----------------------
	 *
	 * RDG=xx,yy (default: RDG=5,3)
	 *
	 *   xx    = Read gap open penalty.
	 *   yy    = Read gap extension penalty.
	 *
	 * Total cost incurred by a read gap = xx + (yy * gap length)
	 *
	 * Penalty for a reference gap
	 * ---------------------------
	 *
	 * RFG=xx,yy (default: RFG=5,3)
	 *
	 *   xx    = Reference gap open penalty.
	 *   yy    = Reference gap extension penalty.
	 *
	 * Total cost incurred by a reference gap = xx + (yy * gap length)
	 *
	 * Minimum score for valid alignment
	 * ---------------------------------
	 *
	 * MIN=xx,yy (defaults: MIN=-0.6,-0.6, or MIN=0.0,0.66 if --local is set)
	 *
	 *   xx,yy = For a read of length N, the total score must be at least
	 *           xx + (read length * yy) for the alignment to be valid.  The
	 *           total score is the sum of all negative penalties (from
	 *           mismatches and gaps) and all positive bonuses.  The minimum
	 *           can be negative (and is by default in global alignment mode).
	 *
	 * N ceiling
	 * ---------
	 *
	 * NCEIL=xx,yy (default: NCEIL=0.0,0.15)
	 *
	 *   xx,yy = For a read of length N, the number of alignment
	 *           positions with an N in either the read or the
	 *           reference cannot exceed
	 *           ceiling = xx + (read length * yy).  If the ceiling is
	 *           exceeded, the alignment is considered invalid.
	 *
	 * Seeds
	 * -----
	 *
	 * SEED=mm,len,ival (default: SEED=0,22)
	 *
	 *   mm   = Maximum number of mismatches allowed within a seed.
	 *          Must be >= 0 and <= 2.  Note that 2-mismatch mode is
	 *          not fully sensitive; i.e. some 2-mismatch seed
	 *          alignments may be missed.
	 *   len  = Length of seed.
	 *   ival = Interval between seeds.  If not specified, seed
	 *          interval is determined by IVAL.
	 *
	 * Seed interval
	 * -------------
	 *
	 * IVAL={L|S|C},xx,yy (default: IVAL=S,1.0,0.0)
	 *
	 *   L  = let interval between seeds be a linear function of the
	 *        read length.  xx and yy are the ant and linear
	 *        coefficients respectively.  In other words, the interval
	 *        equals a * len + b, where len is the read length.
	 *        Intervals less than 1 are rounded up to 1.
	 *   S  = let interval between seeds be a function of the sqaure
	 *        root of the  read length.  xx and yy are the
	 *        coefficients.  In other words, the interval equals
	 *        a * sqrt(len) + b, where len is the read length.
	 *        Intervals less than 1 are rounded up to 1.
	 *   C  = Like S but uses cube root of length instead of square
	 *        root.
	 *
	 * Example 1:
	 *
	 *  SEED=1,10,5 and read sequence is TGCTATCGTACGATCGTAC:
	 *
	 *  The following seeds are extracted from the forward
	 *  representation of the read and aligned to the reference
	 *  allowing up to 1 mismatch:
	 *
	 *  Read:    TGCTATCGTACGATCGTACA
	 *
	 *  Seed 1+: TGCTATCGTA
	 *  Seed 2+:      TCGTACGATC
	 *  Seed 3+:           CGATCGTACA
	 *
	 *  ...and the following are extracted from the reverse-complement
	 *  representation of the read and align to the reference allowing
	 *  up to 1 mismatch:
	 *
	 *  Seed 1-: TACGATAGCA
	 *  Seed 2-:      GATCGTACGA
	 *  Seed 3-:           TGTACGATCG
	 *
	 * Example 2:
	 *
	 *  SEED=1,20,20 and read sequence is TGCTATCGTACGATC.  The seed
	 *  length is 20 but the read is only 15 characters long.  In this
	 *  case, Bowtie2 automatically shrinks the seed length to be equal
	 *  to the read length.
	 *
	 *  Read:    TGCTATCGTACGATC
	 *
	 *  Seed 1+: TGCTATCGTACGATC
	 *  Seed 1-: GATCGTACGATAGCA
	 *
	 * Example 3:
	 *
	 *  SEED=1,10,10 and read sequence is TGCTATCGTACGATC.  Only one seed
	 *  fits on the read; a second seed would overhang the end of the read
	 *  by 5 positions.  In this case, Bowtie2 extracts one seed.
	 *
	 *  Read:    TGCTATCGTACGATC
	 *
	 *  Seed 1+: TGCTATCGTA
	 *  Seed 1-: TACGATAGCA
	 */
	static void parseString(
		String s,
		Boolean        local,
		Boolean        noisyHpolymer,
		Boolean        ignoreQuals,
		int        bonusMatchType,
		int        bonusMatch,
		int        penMmcType,
		int        penMmcMax,
		int        penMmcMin,
		int        penNType,
		int        penN,
		int        penRdEx,
		int        penRfEx,
		int        penRdExLinear,
		int        penRfExLinear,
		SimpleFunc costMin,
		SimpleFunc nCeil,
		Boolean      nCatPair,
		int        multiseedMms,
		int        multiseedLen,
		SimpleFunc multiseedIval,
		int    failStreak,
		int     seedRounds)
	{
		bonusMatchType    = local ? DEFAULT_MATCH_BONUS_TYPE_LOCAL : DEFAULT_MATCH_BONUS_TYPE;
		bonusMatch        = local ? DEFAULT_MATCH_BONUS_LOCAL : DEFAULT_MATCH_BONUS;
		penMmcType        = ignoreQuals ? DEFAULT_MM_PENALTY_TYPE_IGNORE_QUALS :
										  DEFAULT_MM_PENALTY_TYPE;
		penMmcMax         = DEFAULT_MM_PENALTY_MAX;
		penMmcMin         = DEFAULT_MM_PENALTY_MIN;
		penNType          = DEFAULT_N_PENALTY_TYPE;
		penN              = DEFAULT_N_PENALTY;
		
		 double DMAX = Double.MAX_VALUE;
		costMin.init(
			local ? SIMPLE_FUNC_LOG : SIMPLE_FUNC_LINEAR,
			local ? DEFAULT_MIN__LOCAL  : DEFAULT_MIN_,
			local ? DEFAULT_MIN_LINEAR_LOCAL : DEFAULT_MIN_LINEAR);
		nCeil.init(
			SIMPLE_FUNC_LINEAR, 0.0f, DMAX,
			DEFAULT_N_CEIL_, DEFAULT_N_CEIL_LINEAR);
		multiseedIval.init(
			DEFAULT_IVAL, 1.0f, DMAX,
			DEFAULT_IVAL_B, DEFAULT_IVAL_A);
		nCatPair          = DEFAULT_N_CAT_PAIR;

		if(!noisyHpolymer) {
			penRdEx  = DEFAULT_READ_GAP_;
			penRdExLinear = DEFAULT_READ_GAP_LINEAR;
			penRfEx  = DEFAULT_REF_GAP_;
			penRfExLinear = DEFAULT_REF_GAP_LINEAR;
		} else {
			penRdEx  = DEFAULT_READ_GAP__BADHPOLY;
			penRdExLinear = DEFAULT_READ_GAP_LINEAR_BADHPOLY;
			penRfEx  = DEFAULT_REF_GAP__BADHPOLY;
			penRfExLinear = DEFAULT_REF_GAP_LINEAR_BADHPOLY;
		}
		
		multiseedMms      = DEFAULT_SEEDMMS;
		multiseedLen      = gDefaultSeedLen;
		
		EList<String> toks(MISC_CAT);
		String tok;
		Scanner scan = new Scanner(s);
		scan.useDelimiter(";");
		int setting = 0;
		// Get each ;-separated token
		// Old: getline(ss, tok, ';')
		while(scan.hasNext()) {
			tok = scan.next();
			setting++;
			EList<String> etoks(MISC_CAT);
			String etok;
			// Divide into tokens on either side of =
			Scanner ess = new Scanner(tok);
			ess.useDelimiter("=");
			while(ess.hasNext()) {
				etok = ess.next();
				etoks.push_back(etok);
			}
			// Must be exactly 1 =
			if(etoks.size() != 2) {
				System.err.println("Error parsing alignment policy setting " + str(setting)
					 + "; must be bisected by = sign\n" +
					 "Policy: " + str(s));
				assert(false); throw 1;
			}
			// LHS is tag, RHS value
			String tag = etoks[0], val = etoks[1];
			// Separate value into comma-separated tokens
			EList<String> ctoks(MISC_CAT);
			String ctok;
			Scanner css = new Scanner(val);
			while(css.hasNext()) {
				ctok = css.next();
				ctoks.push_back(ctok);
			}
			if(ctoks.size() == 0) {
				System.err.println("Error parsing alignment policy setting " + str(setting)
					 + "; RHS must have at least 1 token\n" +
					 "Policy: " + str(s));
				assert(false); throw 1;
			}
			for(int i = 0; i < ctoks.size(); i++) {
				if(ctoks[i].length() == 0) {
					System.err.println("Error parsing alignment policy setting " + str(setting)
						 + "; token " + str(i+1) + " on RHS had length=0\n"
						 + "Policy: " + str(s));
					assert(false); throw 1;
				}
			}
			// Bonus for a match
			// MA=xx (default: MA=0, or MA=10 if --local is set)
			if(tag == "MA") {
				if(ctoks.size() != 1) {
					System.err.println("Error parsing alignment policy setting " + str(setting)
						 + "; RHS must have 1 token\n" +
						 + "Policy: " + str(s));
					assert(false); throw 1;
				}
				String tmp = ctoks[0];
				Scanner tmpss = new Scanner(tmp);
				bonusMatch = tmpss.nextInt();
			}
			// Scoring for mismatches
			// MMP={Cxx|Q|RQ}
			//        Cxx = ant, where ant is integer xx
			//        Qxx = equal to quality, scaled
			//        R   = equal to maq-rounded quality value (rounded to nearest
			//              10, can't be greater than 30)
			else if(tag == "MMP") {
				if(ctoks.size() > 3) {
					System.err.println("Error parsing alignment policy setting "
						 + "'" + str(tag) + "'"
						 + "; RHS must have at most 3 tokens\n"
						 + "Policy: '" + str(s) + "'");
					assert(false); throw 1;
				}
				if(ctoks[0][0] == 'C') {
					String tmp = ctoks[0].substr(1);
					// Parse ant penalty
					Scanner tmpss = new Scanner(tmp);
					penMmcMax = tmpss.nextInt();
					penMmcMin = penMmcMax;
					// Parse ant penalty
					penMmcType = COST_MODEL_ANT;
				} else if(ctoks[0][0] == 'Q') {
					if(ctoks.size() >= 2) {
						String tmp = ctoks[1];
						Scanner tmpss = new Scanner(tmp);
						penMmcMax = tmpss.nextInt();
					} else {
						penMmcMax = DEFAULT_MM_PENALTY_MAX;
					}
					if(ctoks.size() >= 3) {
						String tmp = ctoks[2];
						Scanner tmpss = new Scanner(tmp);
						penMmcMin = tmpss.nextInt();
					} else {
						penMmcMin = DEFAULT_MM_PENALTY_MIN;
					}
					if(penMmcMin > penMmcMax) {
						System.err.println("Error: Maximum mismatch penalty (" + str(penMmcMax)
							 + ") is less than minimum penalty (" + penMmcMin);
						throw 1;
					}
					// Set type to =quality
					penMmcType = COST_MODEL_QUAL;
				} else if(ctoks[0][0] == 'R') {
					// Set type to=Maq-quality
					penMmcType = COST_MODEL_ROUNDED_QUAL;
				} else {
					System.err.println("Error parsing alignment policy setting "
						 + "'" + str(tag) + "'"
						 + "; RHS must start with C, Q or R\n" +
						 "Policy: '" + str(s) + "'");
					assert(false); throw 1;
				}
			}
			// Scoring for mismatches where read char=N
			// NP={Cxx|Q|RQ}
			//        Cxx = ant, where ant is integer xx
			//        Q   = equal to quality
			//        R   = equal to maq-rounded quality value (rounded to nearest
			//              10, can't be greater than 30)
			else if(tag == "NP") {
				if(ctoks.size() != 1) {
					System.err.println("Error parsing alignment policy setting "
						 + "'" + str(tag) + "'"
						 + "; RHS must have 1 token\n" +
						 "Policy: '" + str(s) + "'");
					assert(false); throw 1;
				}
				if(ctoks[0][0] == 'C') {
					String tmp = ctoks[0].substr(1);
					// Parse ant penalty
					Scanner tmpss = new Scanner(tmp);
					penN = tmpss.nextInt();
					// Parse ant penalty
					penNType = COST_MODEL_ANT;
				} else if(ctoks[0][0] == 'Q') {
					// Set type to =quality
					penNType = COST_MODEL_QUAL;
				} else if(ctoks[0][0] == 'R') {
					// Set type to=Maq-quality
					penNType = COST_MODEL_ROUNDED_QUAL;
				} else {
					System.err.println("Error parsing alignment policy setting "
						 + "'" + str(tag) + "'"
						 + "; RHS must start with C, Q or R\n"
						 + "Policy: '" + s + "'");
					assert(false); throw 1;
				}
			}
			// Scoring for read gaps
			// RDG=xx,yy,zz
			//        xx = read gap open penalty
			//        yy = read gap extension penalty ant coefficient
			//             (defaults to open penalty)
			//        zz = read gap extension penalty linear coefficient
			//             (defaults to 0)
			else if(tag == "RDG") {
				if(ctoks.size() >= 1) {
					Scanner tmpss = new Scanner(ctoks[0]);
					penRdEx = tmpss.nextInt();
				} else {
					penRdEx = noisyHpolymer ?
						DEFAULT_READ_GAP__BADHPOLY :
						DEFAULT_READ_GAP_;
				}
				if(ctoks.size() >= 2) {
					Scanner tmpss = new Scanner(ctoks[1]);
					penRdExLinear = tmpss.nextInt();
				} else {
					penRdExLinear = noisyHpolymer ?
						DEFAULT_READ_GAP_LINEAR_BADHPOLY :
						DEFAULT_READ_GAP_LINEAR;
				}
			}
			// Scoring for reference gaps
			// RFG=xx,yy,zz
			//        xx = ref gap open penalty
			//        yy = ref gap extension penalty ant coefficient
			//             (defaults to open penalty)
			//        zz = ref gap extension penalty linear coefficient
			//             (defaults to 0)
			else if(tag == "RFG") {
				if(ctoks.size() >= 1) {
					Scanner tmpss = new Scanner(ctoks[0]);
					penRfEx = tmpss.nextInt();
				} else {
					penRfEx = noisyHpolymer ?
						DEFAULT_REF_GAP__BADHPOLY :
						DEFAULT_REF_GAP_;
				}
				if(ctoks.size() >= 2) {
					Scanner tmpss = new Scanner(ctoks[1]);
					penRfExLinear = tmpss.nextInt();
				} else {
					penRfExLinear = noisyHpolymer ?
						DEFAULT_REF_GAP_LINEAR_BADHPOLY :
						DEFAULT_REF_GAP_LINEAR;
				}
			}
			// Minimum score as a function of read length
			// MIN=xx,yy
			//        xx = ant coefficient
			//        yy = linear coefficient
			else if(tag == "MIN") {
				PARSE_FUNC(costMin);
			}
			// Per-read N ceiling as a function of read length
			// NCEIL=xx,yy
			//        xx = N ceiling ant coefficient
			//        yy = N ceiling linear coefficient (set to 0 if unspecified)
			else if(tag == "NCEIL") {
				PARSE_FUNC(nCeil);
			}
			/*
			 * Seeds
			 * -----
			 *
			 * SEED=mm,len,ival (default: SEED=0,22)
			 *
			 *   mm   = Maximum number of mismatches allowed within a seed.
			 *          Must be >= 0 and <= 2.  Note that 2-mismatch mode is
			 *          not fully sensitive; i.e. some 2-mismatch seed
			 *          alignments may be missed.
			 *   len  = Length of seed.
			 *   ival = Interval between seeds.  If not specified, seed
			 *          interval is determined by IVAL.
			 */
			else if(tag == "SEED") {
				if(ctoks.size() > 1) {
					System.err.println("Error parsing alignment policy setting "
						 + "'" + str(tag) + "'; RHS must have 1 token, "
						 + "had " + str(ctoks.size()) + ".  "
						 + "Policy: '" + str(s) + "'");
					assert(false); throw 1;
				}
				if(ctoks.size() >= 1) {
					Scanner tmpss = new Scanner(ctoks[0]);
					multiseedMms = tmps.nextInt();
					if(multiseedMms > 1) {
						System.err.println("Error: -N was set to " + str(multiseedMms) + ", but cannot be set greater than 1");
						throw 1;
					}
					if(multiseedMms < 0) {
						System.err.println("Error: -N was set to a number less than 0 (" + multiseedMms + ")");
						throw 1;
					}
				}
			}
			else if(tag == "SEEDLEN") {
				if(ctoks.size() > 1) {
					System.err.println("Error parsing alignment policy setting "
						 + "'" + str(tag) + "'; RHS must have 1 token, "
						 + "had " + str(ctoks.size()) + ".  "
						 + "Policy: '" + s + "'");
					assert(false); throw 1;
				}
				if(ctoks.size() >= 1) {
					Scanner tmpss = new Scanner(ctoks[0]);
					multiseedLen = tmpss.nextInt();
				}
			}
			else if(tag == "DPS") {
				if(ctoks.size() > 1) {
					System.err.println("Error parsing alignment policy setting "
						 + "'" + str(tag) + "'; RHS must have 1 token, "
						 + "had " + str(ctoks.size()) + ".  "
						 + "Policy: '" + str(s) + "'");
					assert(false); throw 1;
				}
				if(ctoks.size() >= 1) {
					Scanner tmpss = new Scanner(ctoks[0]);
					failStreak = tmpss.nextInt();
				}
			}
			else if(tag == "ROUNDS") {
				if(ctoks.size() > 1) {
					System.err.println("Error parsing alignment policy setting "
						 + "'" + str(tag) + "'; RHS must have 1 token, "
						 + "had " + str(ctoks.size()) + ".  "
						 + "Policy: '" + str(s) + "'");
					assert(false); throw 1;
				}
				if(ctoks.size() >= 1) {
					Scanner tmpss = new Scanner(ctoks[0]);
					seedRounds = tmpss.nextInt();
				}
			}
			/*
			 * Seed interval
			 * -------------
			 *
			 * IVAL={L|S|C},a,b (default: IVAL=S,1.0,0.0)
			 *
			 *   L  = let interval between seeds be a linear function of the
			 *        read length.  xx and yy are the ant and linear
			 *        coefficients respectively.  In other words, the interval
			 *        equals a * len + b, where len is the read length.
			 *        Intervals less than 1 are rounded up to 1.
			 *   S  = let interval between seeds be a function of the sqaure
			 *        root of the  read length.  xx and yy are the
			 *        coefficients.  In other words, the interval equals
			 *        a * sqrt(len) + b, where len is the read length.
			 *        Intervals less than 1 are rounded up to 1.
			 *   C  = Like S but uses cube root of length instead of square
			 *        root.
			 */
			else if(tag == "IVAL") {
				PARSE_FUNC(multiseedIval);
			}
			else {
				// Unknown tag
				System.err.println("Unexpected alignment policy setting "
					 + "'" + str(tag) + "'" + endl
					 + "Policy: '" + str(s) + "'");
				assert(false); throw 1;
			}
		}
	}
}

//#endif /*ndef ALIGNER_SEED_POLICY_H_*/
