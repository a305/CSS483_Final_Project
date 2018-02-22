class Scoring {

  enum CostModel {
    COST_MODEL_ROUNDED_QUAL = 1,
    COST_MODEL_QUAL,
	  COST_MODEL_CONSTANT
  };

  public static final int DEFAULT_MATCH_BONUS_TYPE = COST_MODEL_CONSTANT;
  public static final int DEFAULT_MATCH_BONUS = 0;
  public static final int DEFAULT_MATCH_BONUS_TYPE_LOCAL = COST_MODEL_CONSTANT;
  public static final int DEFAULT_MATCH_BONUS_LOCAL = 2;
  public static final int DEFAULT_MM_PENALTY_TYPE = COST_MODEL_QUAL;
  public static final int DEFAULT_MM_PENALTY_TYPE_IGNORE_QUALS = COST_MODEL_CONSTANT;
  public static final int DEFAULT_MM_PENALTY_MAX = 6;
  public static final int DEFAULT_MM_PENALTY_MIN = 2;
  public static final int DEFAULT_N_PENALTY_TYPE = COST_MODEL_CONSTANT;
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
  
  public int matchType;
  public int matchConst;
  public int mmcostType;
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
  public float matchBonuses[256];
  public int mmpens[256];
  public int npens[256];
  
  protected Boolean qualsMatter_;
  
  public Scoring() {
    
  }
  
  public void setMatchBonus() {
    
  }
  
  public void setMmPen() {
    
  }
  
  public void setNPen() {
    
  }
  
  public static float linearFunc() {
    
  }
  
  public int mm() {
    
  }
  
  public int score() {
    
  }
  
  public final long match() {
    
  }
  
  public final long perfectScore() {
    
  }
  
  public final Boolean qualitiesMatter() {
    
  }
  
  public final int n() {
    
  }
  
  public final int ins() {
    
  }
  
  public final int del() {
    
  }
  
  public final Boolean scoreFilter() {
    
  }
  
  public final int maxReadGaps() {
    
  }
  
  public final int maxRefGaps() {
    
  }
  
  public final Boolean nFilter() {
    
  }
  
  public void nFilterPair() {
    
  }
  
  public final int readGapOpen() {
    
  }
  
  public final int refGapOpen() {
    
  }
  
  public final int readGapExtend() {
    
  }
  
  public final int refGapExtend() {
  
  }
  
  public void initPens() {
    
  }
  
  public static Scoring base1() {
  
  }
}
