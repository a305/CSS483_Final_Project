package com.uwb.bt2j.aligner;
public class Scoring<T> {

 public enum CostModel {
    COST_MODEL_ROUNDED_QUAL(1),
    COST_MODEL_QUAL(2),
	COST_MODEL_CONSTANT(3);
	  
	 private int x;
	  CostModel(int y){x = y;}
  };

  public static final CostModel DEFAULT_MATCH_BONUS_TYPE = CostModel.COST_MODEL_CONSTANT;
  public static final int DEFAULT_MATCH_BONUS = 0;
  public static final CostModel DEFAULT_MATCH_BONUS_TYPE_LOCAL = CostModel.COST_MODEL_CONSTANT;
  public static final int DEFAULT_MATCH_BONUS_LOCAL = 2;
  public static final CostModel DEFAULT_MM_PENALTY_TYPE = CostModel.COST_MODEL_QUAL;
  public static final CostModel DEFAULT_MM_PENALTY_TYPE_IGNORE_QUALS = CostModel.COST_MODEL_CONSTANT;
  public static final int DEFAULT_MM_PENALTY_MAX = 6;
  public static final int DEFAULT_MM_PENALTY_MIN = 2;
  public static final CostModel DEFAULT_N_PENALTY_TYPE = CostModel.COST_MODEL_CONSTANT;
  public static final int DEFAULT_N_PENALTY = 1;
  public static final float DEFAULT_MIN_CONST = -0.6f;
  public static final float DEFAULT_MIN_LINEAR = -0.6f;
  public static final float DEFAULT_MIN_CONST_LOCAL = 20.0f;
  public static final float DEFAULT_MIN_LINEAR_LOCAL = 8.0f;
  public static final float DEFAULT_N_CEIL_CONST = 0.0F;
  public static final float DEFAULT_n_ceil_linear = 0.15f;
  public static final Boolean DEFAULT_N_CAT_PAIR = false;
  public static final int DEFAULT_READ_GAP_CONST = 5;
  public static final int DEFAULT_READ_GAP_LINEAR = 3;
  public static final int DEFAULT_READ_GAP_CONST_BADHPOLY = 3;
  public static final int DEFAULT_READ_GAP_LINEAR_BADHPOLY = 1;
  public static final int DEFAULT_REF_GAP_CONST = 5;
  public static final int DEFAULT_REF_GAP_LINEAR = 3;
  public static final int DEFAULT_REF_GAP_CONST_BADHPOLY = 3;
  public static final int DEFAULT_REF_GAP_LINEAR_BADHPOLY = 1;
  
  public CostModel matchType;
  public int matchConst;
  public CostModel mmcostType;
  public int mmpMax;
  public int mmpMin;
  public SimpleFunc scoreMin;
  public SimpleFunc nCeil;
  public int npenType;
  public int npen;
  public Boolean ncatpair;
  public int rdGapConst;
  public int rfGapConst;
  public int rdGapLinear;
  public int rfGapLinear;
  public int gapbar;
  public Boolean monotone;
  public float matchBonuses[];
  public int mmpens[];
  public int npens[];
  
  protected Boolean qualsMatter_;
  
  public Scoring(int   mat,          // reward for a match
			int   mmcType,      // how to penalize mismatches
		    int   mmpMax_,      // maximum mismatch penalty
		    int   mmpMin_,      // minimum mismatch penalty
		    SimpleFunc scoreMin_,   // minimum score for valid alignment; const coeff
			SimpleFunc nCeil_,      // max # ref Ns allowed in alignment; const coeff
		    int   nType,        // how to penalize Ns in the read
		    int   n,            // constant if N pelanty is a constant
			Boolean  ncat,         // whether to concatenate mates before N filtering
		    int   rdGpConst,    // constant coeff for cost of gap in the read
		    int   rfGpConst,    // constant coeff for cost of gap in the ref
		    int   rdGpLinear,   // coeff of linear term for cost of gap in read
		    int   rfGpLinear,   // coeff of linear term for cost of gap in ref
			int     gapbar_)    // # rows at top/bot can only be entered diagonally
  {
	  matchType    = CostModel.COST_MODEL_CONSTANT;
		matchConst   = mat;
		mmcostType   = mmcType;
		mmpMax       = mmpMax_;
		mmpMin       = mmpMin_;
		scoreMin     = scoreMin_;
		nCeil        = nCeil_;
		npenType     = nType;
		npen         = n;
		ncatpair     = ncat;
		rdGapConst   = rdGpConst;
		rfGapConst   = rfGpConst;
		rdGapLinear  = rdGpLinear;
		rfGapLinear  = rfGpLinear;
		qualsMatter_ = mmcostType != CostModel.COST_MODEL_CONSTANT;
		gapbar       = gapbar_;
		monotone     = matchType == CostModel.COST_MODEL_CONSTANT && matchConst == 0;
		initPens<Integer>(mmpens, mmcostType, mmpMin_, mmpMax_);
		initPens<Integer>(npens, npenType, npen, npen);
		initPens<float>(matchBonuses, matchType, matchConst, matchConst);
  }
  
  public Boolean scoreFilter(long minsc, double rdlen) {
	  long sc = (long)(rdlen * match(30));
		return sc >= minsc;
  }
  
  public int maxReadGaps(long minsc, double rdlen) {
		// Score if all characters match.  TODO: remove assumption that match bonus
		// is independent of quality value.
		long sc = (long)(rdlen * match(30));
		// Now convert matches to read gaps until sc calls below minsc
		Boolean first = true;
		int num = 0;
		while(sc >= minsc) {
			if(first) {
				first = false;
				// Subtract both penalties
				sc -= readGapOpen();
			} else {
				// Subtract just the extension penalty
				sc -= readGapExtend();
			}
			num++;
		}
		return num-1;
  }
  
  public int maxRefGaps(long minsc, double rdlen) {
	// Score if all characters match.  TODO: remove assumption that match bonus
		// is independent of quality value.
		long sc = (long)(rdlen * match(30));
		// Now convert matches to read gaps until sc calls below minsc
		Boolean first = true;
		int num = 0;
		while(sc >= minsc) {
			sc -= match(30);
			if(first) {
				first = false;
				// Subtract both penalties
				sc -= refGapOpen();
			} else {
				// Subtract just the extension penalty
				sc -= refGapExtend();
			}
			num++;
		}
		return num-1;
  }
  
  public void setMatchBonus(int bonus) {
		matchType  = CostModel.COST_MODEL_CONSTANT;
		matchConst = bonus;
		initPens<float>(matchBonuses, matchType, matchConst, matchConst);
  }
  
  public void setMmPen(CostModel mmType_, int mmpMax_, int mmpMin_) {
		mmcostType = mmType_;
		mmpMax     = mmpMax_;
		mmpMin     = mmpMin_;
		initPens<int>(mmpens, mmcostType, mmpMin, mmpMax);
  }
  
  public void setNPen(int nType, int n) {
	  npenType     = nType;
		npen         = n;
		initPens<int>(npens, npenType, npen, npen);
  }
  
  public static float linearFunc(long x, float cnst, float lin) {
	  return (float)((double)cnst + ((double)lin * x));
  }
  
  public int mm(int rdc, int refm, int q) {
	  return (rdc > 3 || refm > 15) ? npens[q] : mmpens[q];
  }
  
  public int mm(int rdc, int q) {
	  return (rdc > 3) ? npens[q] : mmpens[q];
  }
  
  public int mm(int q) {
	  return q < 255 ? mmpens[q] : mmpens[255];
  }
  
  public int score(int rdc, int refm, int q) {
	  if(rdc > 3 || refm > 15) {
			return -npens[q];
		}
		if((refm & (1 << rdc)) != 0) {
			return (int)matchBonuses[q];
		} else {
			return -mmpens[q];
		}
  }
  
  public int score(int rdc, int refm, int q, int& ns) {
	  if(rdc > 3 || refm > 15) {
			ns++;
			return -npens[q];
		}
		if((refm & (1 << rdc)) != 0) {
			return (int)matchBonuses[q];
		} else {
			return -mmpens[q];
		}
  }
  
  public final long match(int q) {
		return (long)((q < 255 ? matchBonuses[q] : matchBonuses[255]) + 0.5f);
  }
  
  public final long match() {
    return match(30);
  }
  
  public final long perfectScore(long rdlen) {
		if(monotone) {
			return 0;
		} else {
			return rdlen * match(30);
		}
  }
  
  public final Boolean qualitiesMatter() {
	  return qualsMatter_;
  }
  
  public final int n(int q) {
	  return q < 255 ? npens[q] : npens[255];
  }
  
  public final int ins(int ext) {
	  if(ext == 0) return readGapOpen();
		return readGapExtend();
  }
  
  public final int del(int ext) {
	  if(ext == 0) return refGapOpen();
		return refGapExtend();
  }
  
  public final Boolean scoreFilter() {
    
  }
  
  public final Boolean nFilter(BTDnaString rd, double ns) {
	  double rdlen = rd.length();
		double maxns = nCeil.f<double>((double)rdlen);
		for(double i = 0; i < rdlen; i++) {
			if(rd[i] == 4) {
				ns++;
				if(ns > maxns) {
					return false; // doesn't pass
				}
			}
		}
		return true; // passes
  }
  
  public void nFilterPair(
		  BTDnaString rd1, // mate 1
			BTDnaString rd2, // mate 2
			double ns1,            // # Ns in mate 1
			double ns2,            // # Ns in mate 2
			Boolean filt1,            // true -> mate 1 rejected by filter
			Boolean filt2) {
		// Both fail to pass by default
		filt1 = filt2 = false;
		if(rd1 != null && rd2 != null && ncatpair) {
			double rdlen1 = rd1.length();
			double rdlen2 = rd2.length();
			double maxns = nCeil.f<double>((double)(rdlen1 + rdlen2));
			for(double i = 0; i < rdlen1; i++) {
				if((rd1)[i] == 4) ns1++;
				if(ns1 > maxns) {
					// doesn't pass
					return;
				}
			}
			for(double i = 0; i < rdlen2; i++) {
				if((*rd2)[i] == 4) ns2++;
				if(ns2 > maxns) {
					// doesn't pass
					return;
				}
			}
			// Both pass
			filt1 = filt2 = true;
		} else {
			if(rd1 != null) filt1 = nFilter(rd1, ns1);
			if(rd2 != null) filt2 = nFilter(rd2, ns2);
		}
  }
  
  public final int readGapOpen() {
	  return rdGapConst + rdGapLinear;
  }
  
  public final int refGapOpen() {
	  return rfGapConst + rfGapLinear;
  }
  
  public final int readGapExtend() {
	  return rdGapLinear;
  }
  
  public final int refGapExtend() {
	  return rfGapLinear;
  }
  
  public void initPens(T pens, CostModel type, int consMin, int consMax) {
		if(type == CostModel.COST_MODEL_ROUNDED_QUAL) {
			for(int i = 0; i < 256; i++) {
				pens[i] = (T)qualRounds[i];
			}
		} else if(type == CostModel.COST_MODEL_QUAL) {
			for(int i = 0; i < 256; i++) {
				int ii = Math.min(i, 40); // TODO: Bit hacky, this
				float frac = (float)ii / 40.0f;
				pens[i] = consMin + (T)(frac * (consMax-consMin));
			}
		} else if(type == CostModel.COST_MODEL_CONSTANT) {
			for(int i = 0; i < 256; i++) {
				pens[i] = (T)consMax;
			}
		}
  }
  
  public static Scoring base1() {
	  double DMAX = Double.MAX_VALUE;
		SimpleFunc scoreMin(CostModel.SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 37.0f, 0.3f);
		SimpleFunc nCeil(CostModel.SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 2.0f, 0.1f);
		return Scoring(
			1,                       // reward for a match
			CostModel.COST_MODEL_CONSTANT,     // how to penalize mismatches
			3,                       // max mismatch pelanty
			3,                       // min mismatch pelanty
			scoreMin,                // score min: 37 + 0.3x
			nCeil,                   // n ceiling: 2 + 0.1x
			CostModel.COST_MODEL_CONSTANT,     // how to penalize Ns in the read
			3,                       // constant if N pelanty is a constant
			false,                   // concatenate mates before N filtering?
			11,                      // constant coeff for gap in read
			11,                      // constant coeff for gap in ref
			4,                       // linear coeff for gap in read
			4,                       // linear coeff for gap in ref
			5);                      // 5 rows @ top/bot diagonal-entrance-only
  }
}
