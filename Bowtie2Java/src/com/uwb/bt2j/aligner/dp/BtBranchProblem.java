package com.uwb.bt2j.aligner.dp;

import com.uwb.bt2j.aligner.Scoring;
import com.uwb.bt2j.util.Checkpointer;

class BtBranchProblem {

  protected String qry_;
  protected String qual_;
  protected double qrylen_;
  protected String ref_;
  protected double reflen_;
  protected double treflen_;
  protected double refid_;
  protected double refoff_;
  protected Boolean fw_;
  protected DPRect rect_;
  protected double row_;
  protected double col_;
  protected double targ_;
  protected Checkpointer cper_;
  protected Boolean fill_;
  protected Boolean usecp_;
  protected Scoring sc_;
  protected double nceil_;
  
  public BtBranchProblem() {
    reset();
  }
  
  public void initRef(
		 String qry,    // query string (along rows)
		 String qual,   // query quality string (along rows)
		 int  qrylen, // query string (along rows) length
		 String ref,    // reference string (along columns)
		 long reflen, // in-rectangle reference string length
		 long treflen,// total reference string length
	 	 long refid,  // reference id
		 long refoff, // reference offset
		 boolean fw,     // orientation of problem
		 DPRect rect,   // dynamic programming rectangle filled out
		 Checkpointer cper,   // checkpointer
		 Scoring sc,     // scoring scheme
		 int nceil)  // max # Ns allowed in alignment
  {
		qry_     = qry;
		qual_    = qual;
		qrylen_  = qrylen;
		ref_     = ref;
		reflen_  = reflen;
		treflen_ = treflen;
		refid_   = refid;
		refoff_  = refoff;
		fw_      = fw;
		rect_    = rect;
		cper_    = cper;
		sc_      = sc;
		nceil_   = nceil;
  }
  
  public void initBt(
		  double   row,   // row
			double   col,   // column
			boolean     fill,  // use a filling rather than a backtracking strategy
			boolean     usecp, // use checkpoints to short-circuit while backtracking
			long targ)  // target score
  {
	  	row_    = row;
		col_    = col;
		targ_   = targ;
		fill_   = fill;
		usecp_  = usecp;
  }
  
  public void reset() {
	  	qry_ = qual_ = ref_ = null;
		cper_ = null;
		rect_ = null;
		sc_ = null;
		qrylen_ = reflen_ = treflen_ = refid_ = refoff_ = row_ = col_ = targ_ = nceil_ = 0;
		fill_ = fw_ = usecp_ = false;
  }
  
  public boolean inited() {
	  return qry_ != null;
  }
  
  public double reflen() {
	  return reflen_;
  }
  
  public double treflen() {
	  return treflen_;
  }
  
  
}