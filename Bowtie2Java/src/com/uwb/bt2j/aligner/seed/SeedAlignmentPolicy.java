package com.uwb.bt2j.aligner.seed;

import com.uwb.bt2j.aligner.Scoring;
import com.uwb.bt2j.aligner.Scoring.CostModel;
import com.uwb.bt2j.indexer.EList;
import com.uwb.bt2j.aligner.SimpleFunc;

public class SeedAlignmentPolicy {
	
	public static final int DEFAULT_SEEDMMS = 0;
	public static final int DEFAULT_SEEDLEN = 22;
	public static final int DEFAULT_LOCAL_SEEDLEN = 20;
	public static final byte DEFAULT_IVAL = SimpleFunc.SIMPLE_FUNC_SQRT;
	public static final float DEFAULT_IVAL_A = 1.15f;
	public static final float DEFAULT_IVAL_B = 0.0f;
	public static final int DEFAULT_UNGAPPED_HITS = 6;
	public static int gDefaultSeedLen;
	
	public static void parseString(
			String s,
			boolean        local,
			boolean        noisyHpolymer,
			boolean        ignoreQuals,
			int        bonusMatchType,
			int        bonusMatch,
			int        penMmcType,
			int        penMmcMax,
			int        penMmcMin,
			CostModel        penNType,
			int        penN,
			int        penRdExConst,
			int        penRfExConst,
			int        penRdExLinear,
			int        penRfExLinear,
			SimpleFunc costMin,
			SimpleFunc nCeil,
			boolean       nCatPair,
			int        multiseedMms,
			int        multiseedLen,
			SimpleFunc multiseedIval,
			double     failStreak,
			double     seedRounds
			) {
		bonusMatchType    = local ? Scoring.DEFAULT_MATCH_BONUS_TYPE_LOCAL : Scoring.DEFAULT_MATCH_BONUS_TYPE;
		bonusMatch        = local ? Scoring.DEFAULT_MATCH_BONUS_LOCAL : Scoring.DEFAULT_MATCH_BONUS;
		penMmcType        = ignoreQuals ? Scoring.DEFAULT_MM_PENALTY_TYPE_IGNORE_QUALS :
			Scoring.DEFAULT_MM_PENALTY_TYPE;
		penMmcMax         = Scoring.DEFAULT_MM_PENALTY_MAX;
		penMmcMin         = Scoring.DEFAULT_MM_PENALTY_MIN;
		penNType          = Scoring.DEFAULT_N_PENALTY_TYPE;
		penN              = Scoring.DEFAULT_N_PENALTY;
		
		double DMAX = Double.MAX_VALUE;
		costMin.init(
			local ? SimpleFunc.SIMPLE_FUNC_LOG : SimpleFunc.SIMPLE_FUNC_LINEAR,
			local ? Scoring.DEFAULT_MIN_CONST_LOCAL  : Scoring.DEFAULT_MIN_CONST,
			local ? Scoring.DEFAULT_MIN_LINEAR_LOCAL : Scoring.DEFAULT_MIN_LINEAR);
		nCeil.init(
			SimpleFunc.SIMPLE_FUNC_LINEAR, 0.0f, DMAX,
			Scoring.DEFAULT_N_CEIL_CONST, Scoring.DEFAULT_N_CEIL_LINEAR);
		multiseedIval.init(
			DEFAULT_IVAL, 1.0f, DMAX,
			DEFAULT_IVAL_B, DEFAULT_IVAL_A);
		nCatPair          = Scoring.DEFAULT_N_CAT_PAIR;

		if(!noisyHpolymer) {
			penRdExConst  = Scoring.DEFAULT_READ_GAP_CONST;
			penRdExLinear = Scoring.DEFAULT_READ_GAP_LINEAR;
			penRfExConst  = Scoring.DEFAULT_REF_GAP_CONST;
			penRfExLinear = Scoring.DEFAULT_REF_GAP_LINEAR;
		} else {
			penRdExConst  = Scoring.DEFAULT_READ_GAP_CONST_BADHPOLY;
			penRdExLinear = Scoring.DEFAULT_READ_GAP_LINEAR_BADHPOLY;
			penRfExConst  = Scoring.DEFAULT_REF_GAP_CONST_BADHPOLY;
			penRfExLinear = Scoring.DEFAULT_REF_GAP_LINEAR_BADHPOLY;
		}
		
		multiseedMms      = DEFAULT_SEEDMMS;
		multiseedLen      = gDefaultSeedLen;
		
		EList<String> toks(MISC_CAT);
		String tok;
		istringstream ss(s);
		int setting = 0;
		// Get each ;-separated token
		while(getline(ss, tok, ';')) {
			setting++;
			EList<string> etoks(MISC_CAT);
			string etok;
			// Divide into tokens on either side of =
			istringstream ess(tok);
			while(getline(ess, etok, '=')) {
				etoks.push_back(etok);
			}
			// Must be exactly 1 =
			if(etoks.size() != 2) {
				cerr << "Error parsing alignment policy setting " << setting
				     << "; must be bisected by = sign" << endl
					 << "Policy: " << s.c_str() << endl;
				assert(false); throw 1;
			}
			// LHS is tag, RHS value
			string tag = etoks[0], val = etoks[1];
			// Separate value into comma-separated tokens
			EList<string> ctoks(MISC_CAT);
			string ctok;
			istringstream css(val);
			while(getline(css, ctok, ',')) {
				ctoks.push_back(ctok);
			}
			if(ctoks.size() == 0) {
				cerr << "Error parsing alignment policy setting " << setting
				     << "; RHS must have at least 1 token" << endl
					 << "Policy: " << s.c_str() << endl;
				assert(false); throw 1;
			}
			for(size_t i = 0; i < ctoks.size(); i++) {
				if(ctoks[i].length() == 0) {
					cerr << "Error parsing alignment policy setting " << setting
					     << "; token " << i+1 << " on RHS had length=0" << endl
						 << "Policy: " << s.c_str() << endl;
					assert(false); throw 1;
				}
			}
			// Bonus for a match
			// MA=xx (default: MA=0, or MA=10 if --local is set)
			if(tag == "MA") {
				if(ctoks.size() != 1) {
					cerr << "Error parsing alignment policy setting " << setting
					     << "; RHS must have 1 token" << endl
						 << "Policy: " << s.c_str() << endl;
					assert(false); throw 1;
				}
				string tmp = ctoks[0];
				istringstream tmpss(tmp);
				tmpss >> bonusMatch;
			}
			// Scoring for mismatches
			// MMP={Cxx|Q|RQ}
			//        Cxx = constant, where constant is integer xx
			//        Qxx = equal to quality, scaled
			//        R   = equal to maq-rounded quality value (rounded to nearest
			//              10, can't be greater than 30)
			else if(tag == "MMP") {
				if(ctoks.size() > 3) {
					cerr << "Error parsing alignment policy setting "
					     << "'" << tag.c_str() << "'"
					     << "; RHS must have at most 3 tokens" << endl
						 << "Policy: '" << s.c_str() << "'" << endl;
					assert(false); throw 1;
				}
				if(ctoks[0][0] == 'C') {
					string tmp = ctoks[0].substr(1);
					// Parse constant penalty
					istringstream tmpss(tmp);
					tmpss >> penMmcMax;
					penMmcMin = penMmcMax;
					// Parse constant penalty
					penMmcType = COST_MODEL_CONSTANT;
				} else if(ctoks[0][0] == 'Q') {
					if(ctoks.size() >= 2) {
						string tmp = ctoks[1];
						istringstream tmpss(tmp);
						tmpss >> penMmcMax;
					} else {
						penMmcMax = DEFAULT_MM_PENALTY_MAX;
					}
					if(ctoks.size() >= 3) {
						string tmp = ctoks[2];
						istringstream tmpss(tmp);
						tmpss >> penMmcMin;
					} else {
						penMmcMin = DEFAULT_MM_PENALTY_MIN;
					}
					if(penMmcMin > penMmcMax) {
						cerr << "Error: Maximum mismatch penalty (" << penMmcMax
						     << ") is less than minimum penalty (" << penMmcMin
							 << endl;
						throw 1;
					}
					// Set type to =quality
					penMmcType = COST_MODEL_QUAL;
				} else if(ctoks[0][0] == 'R') {
					// Set type to=Maq-quality
					penMmcType = COST_MODEL_ROUNDED_QUAL;
				} else {
					cerr << "Error parsing alignment policy setting "
					     << "'" << tag.c_str() << "'"
					     << "; RHS must start with C, Q or R" << endl
						 << "Policy: '" << s.c_str() << "'" << endl;
					assert(false); throw 1;
				}
			}
			// Scoring for mismatches where read char=N
			// NP={Cxx|Q|RQ}
			//        Cxx = constant, where constant is integer xx
			//        Q   = equal to quality
			//        R   = equal to maq-rounded quality value (rounded to nearest
			//              10, can't be greater than 30)
			else if(tag == "NP") {
				if(ctoks.size() != 1) {
					cerr << "Error parsing alignment policy setting "
					     << "'" << tag.c_str() << "'"
					     << "; RHS must have 1 token" << endl
						 << "Policy: '" << s.c_str() << "'" << endl;
					assert(false); throw 1;
				}
				if(ctoks[0][0] == 'C') {
					string tmp = ctoks[0].substr(1);
					// Parse constant penalty
					istringstream tmpss(tmp);
					tmpss >> penN;
					// Parse constant penalty
					penNType = COST_MODEL_CONSTANT;
				} else if(ctoks[0][0] == 'Q') {
					// Set type to =quality
					penNType = COST_MODEL_QUAL;
				} else if(ctoks[0][0] == 'R') {
					// Set type to=Maq-quality
					penNType = COST_MODEL_ROUNDED_QUAL;
				} else {
					cerr << "Error parsing alignment policy setting "
					     << "'" << tag.c_str() << "'"
					     << "; RHS must start with C, Q or R" << endl
						 << "Policy: '" << s.c_str() << "'" << endl;
					assert(false); throw 1;
				}
			}
			// Scoring for read gaps
			// RDG=xx,yy,zz
			//        xx = read gap open penalty
			//        yy = read gap extension penalty constant coefficient
			//             (defaults to open penalty)
			//        zz = read gap extension penalty linear coefficient
			//             (defaults to 0)
			else if(tag == "RDG") {
				if(ctoks.size() >= 1) {
					istringstream tmpss(ctoks[0]);
					tmpss >> penRdExConst;
				} else {
					penRdExConst = noisyHpolymer ?
						DEFAULT_READ_GAP_CONST_BADHPOLY :
						DEFAULT_READ_GAP_CONST;
				}
				if(ctoks.size() >= 2) {
					istringstream tmpss(ctoks[1]);
					tmpss >> penRdExLinear;
				} else {
					penRdExLinear = noisyHpolymer ?
						DEFAULT_READ_GAP_LINEAR_BADHPOLY :
						DEFAULT_READ_GAP_LINEAR;
				}
			}
			// Scoring for reference gaps
			// RFG=xx,yy,zz
			//        xx = ref gap open penalty
			//        yy = ref gap extension penalty constant coefficient
			//             (defaults to open penalty)
			//        zz = ref gap extension penalty linear coefficient
			//             (defaults to 0)
			else if(tag == "RFG") {
				if(ctoks.size() >= 1) {
					istringstream tmpss(ctoks[0]);
					tmpss >> penRfExConst;
				} else {
					penRfExConst = noisyHpolymer ?
						DEFAULT_REF_GAP_CONST_BADHPOLY :
						DEFAULT_REF_GAP_CONST;
				}
				if(ctoks.size() >= 2) {
					istringstream tmpss(ctoks[1]);
					tmpss >> penRfExLinear;
				} else {
					penRfExLinear = noisyHpolymer ?
						DEFAULT_REF_GAP_LINEAR_BADHPOLY :
						DEFAULT_REF_GAP_LINEAR;
				}
			}
			// Minimum score as a function of read length
			// MIN=xx,yy
			//        xx = constant coefficient
			//        yy = linear coefficient
			else if(tag == "MIN") {
				PARSE_FUNC(costMin);
			}
			// Per-read N ceiling as a function of read length
			// NCEIL=xx,yy
			//        xx = N ceiling constant coefficient
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
					cerr << "Error parsing alignment policy setting "
					     << "'" << tag.c_str() << "'; RHS must have 1 token, "
						 << "had " << ctoks.size() << ".  "
						 << "Policy: '" << s.c_str() << "'" << endl;
					assert(false); throw 1;
				}
				if(ctoks.size() >= 1) {
					istringstream tmpss(ctoks[0]);
					tmpss >> multiseedMms;
					if(multiseedMms > 1) {
						cerr << "Error: -N was set to " << multiseedMms << ", but cannot be set greater than 1" << endl;
						throw 1;
					}
					if(multiseedMms < 0) {
						cerr << "Error: -N was set to a number less than 0 (" << multiseedMms << ")" << endl;
						throw 1;
					}
				}
			}
			else if(tag == "SEEDLEN") {
				if(ctoks.size() > 1) {
					cerr << "Error parsing alignment policy setting "
					     << "'" << tag.c_str() << "'; RHS must have 1 token, "
						 << "had " << ctoks.size() << ".  "
						 << "Policy: '" << s.c_str() << "'" << endl;
					assert(false); throw 1;
				}
				if(ctoks.size() >= 1) {
					istringstream tmpss(ctoks[0]);
					tmpss >> multiseedLen;
				}
			}
			else if(tag == "DPS") {
				if(ctoks.size() > 1) {
					cerr << "Error parsing alignment policy setting "
					     << "'" << tag.c_str() << "'; RHS must have 1 token, "
						 << "had " << ctoks.size() << ".  "
						 << "Policy: '" << s.c_str() << "'" << endl;
					assert(false); throw 1;
				}
				if(ctoks.size() >= 1) {
					istringstream tmpss(ctoks[0]);
					tmpss >> failStreak;
				}
			}
			else if(tag == "ROUNDS") {
				if(ctoks.size() > 1) {
					cerr << "Error parsing alignment policy setting "
					     << "'" << tag.c_str() << "'; RHS must have 1 token, "
						 << "had " << ctoks.size() << ".  "
						 << "Policy: '" << s.c_str() << "'" << endl;
					assert(false); throw 1;
				}
				if(ctoks.size() >= 1) {
					istringstream tmpss(ctoks[0]);
					tmpss >> seedRounds;
				}
			}
			/*
			 * Seed interval
			 * -------------
			 *
			 * IVAL={L|S|C},a,b (default: IVAL=S,1.0,0.0)
			 *
			 *   L  = let interval between seeds be a linear function of the
			 *        read length.  xx and yy are the constant and linear
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
				cerr << "Unexpected alignment policy setting "
					 << "'" << tag.c_str() << "'" << endl
					 << "Policy: '" << s.c_str() << "'" << endl;
				assert(false); throw 1;
			}
		}
	}
	
	public static int parseFuncType(String otype) {
		String type = otype;
		if(type == "C" || type == "Constant") {
			return SimpleFunc.SIMPLE_FUNC_CONST;
		} else if(type == "L" || type == "Linear") {
			return SimpleFunc.SIMPLE_FUNC_LINEAR;
		} else if(type == "S" || type == "Sqrt") {
			return SimpleFunc.SIMPLE_FUNC_SQRT;
		} else if(type == "G" || type == "Log") {
			return SimpleFunc.SIMPLE_FUNC_LOG;
		}
		System.err.println("Error: Bad function type '" + otype
		          + "'.  Should be C (constant), L (linear), "
		          + "S (square root) or G (natural log).");
		return 1;
	}
}
