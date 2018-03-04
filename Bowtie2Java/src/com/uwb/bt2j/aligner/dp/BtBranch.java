package com.uwb.bt2j.aligner.dp;

import com.uwb.bt2j.aligner.Edit;

public class BtBranch {

  protected double parentId_;
  protected double penalty_;
  protected double score_st_;
  protected double score_en_;
  protected double len_;
  protected long row_;
  protected long col_;
  protected Edit e_;
  
  protected boolean root_;
  protected boolean curtailed_;
  
  public enum BTStatus {
	  BT_NOT_FOUND(1),      // could not obtain the backtrace because it
              // overlapped a previous solution
	  BT_FOUND(2),              // obtained a valid backtrace
	  BT_REJECTED_N(3),         // backtrace rejected because it had too many Ns
	  BT_REJECTED_CORE_DIAG(4);  // backtrace rejected because it failed to overlap a
              // core diagonal
	  private int x;
	  BTStatus(int y){x = y;}
  }
  
  public BtBranch() {
	  reset();
  }
  
  public BtBranch(
		  BtBranchProblem prob,
			int parentId,
			long penalty,
			long score_en,
			long row,
			long col,
			Edit e,
			int hef,
			boolean root,
			boolean extend){
	  init(prob, parentId, penalty, score_en, row, col, e, hef, root, extend);
  }
  
  public void reset() {
	  	parentId_ = 0;
		score_st_ = score_en_ = len_ = row_ = col_ = 0;
		curtailed_ = false;
		e_.reset();
  }
  
  public void init(
		  BtBranchProblem prob,
			int parentId,
			long penalty,
			long score_en,
			long row,
			long col,
			Edit e,
			int hef,
			boolean root,
			boolean extend)
  {
	  score_en_ = score_en;
		penalty_ = penalty;
		score_st_ = score_en_;
		row_ = row;
		col_ = col;
		parentId_ = parentId;
		e_ = e;
		root_ = root;
		assert(!root_ || parentId == 0);
		assert_lt(row, (int64_t)prob.qrylen_);
		assert_lt(col, (int64_t)prob.reflen_);
		// First match to check is diagonally above and to the left of the cell
		// where the edit occurs
		int64_t rowc = row;
		int64_t colc = col;
		len_ = 0;
		if(e.inited() && e.isMismatch()) {
			rowc--; colc--;
			len_ = 1;
		}
		int64_t match = prob.sc_->match();
		bool cp = prob.usecp_;
		size_t iters = 0;
		curtailed_ = false;
		if(extend) {
			while(rowc >= 0 && colc >= 0) {
				int rfm = prob.ref_[colc];
				assert_range(0, 16, rfm);
				int rdc = prob.qry_[rowc];
				bool matches = (rfm & (1 << rdc)) != 0;
				if(!matches) {
					// What's the mismatch penalty?
					break;
				}
				// Get score from checkpointer
				score_st_ += match;
				if(cp && rowc - 1 >= 0 && colc - 1 >= 0 &&
				   prob.cper_->isCheckpointed(rowc - 1, colc - 1))
				{
					// Possibly prune
					int16_t cpsc;
					cpsc = prob.cper_->scoreTriangle(rowc - 1, colc - 1, hef);
					if(cpsc + score_st_ < prob.targ_) {
						curtailed_ = true;
						break;
					}
				}
				iters++;
				rowc--; colc--;
			}
		}
		assert_geq(rowc, -1);
		assert_geq(colc, -1);
		len_ = (int64_t)row - rowc;
  }
  
  public boolean isSolution(BtBranchProblem prob) {
	  	boolean end2end = prob.sc_.monotone;
		return score_st_ == prob.targ_ && (!end2end || endsInFirstRow());
  }
  
  public boolean isValid(BtBranchProblem prob) {
	  long scoreFloor = prob.sc_.monotone ? Long.MIN_VALUE : 0;
		if(score_st_ < scoreFloor) {
			// Dipped below the score floor
			return false;
		}
		if(isSolution(prob)) {
			// It's a solution, so it's also valid
			return true;
		}
		if((long)len_ > row_) {
			// Went all the way to the top row
			//assert_leq(score_st_, prob.targ_);
			return score_st_ == prob.targ_;
		} else {
			long match = prob.sc_.match();
			double bonusLeft = (row_ + 1 - len_) * match;
			return score_st_ + bonusLeft >= prob.targ_;
		}
  }
  
  public boolean overlap(BtBranchProblem prob, BtBranch bt) {
	// Calculate this branch's diagonal
			double fromend = prob.qrylen_ - row_ - 1;
			double diag = fromend + col_;
			double lo = 0;
			double hi = row_ + 1;
			if(len_ == 0) {
				lo = row_;
			} else {
				lo = row_ - (len_ - 1);
			}
			// Calculate other branch's diagonal
			double ofromend = prob.qrylen_ - bt.row_ - 1;
			double odiag = ofromend + bt.col_;
			if(diag != odiag) {
				return false;
			}
			double olo = 0;
			long ohi = bt.row_ + 1;
			if(bt.len_ == 0) {
				olo = bt.row_;
			} else {
				olo = bt.row_ - (bt.len_ - 1);
			}
			double losm = olo;
			double hism = ohi;
			if(hi - lo < ohi - olo) {
				double tmp = lo;
				lo = losm;
				losm = tmp;
				tmp = hi;
				hi = hism;
				hism = tmp;
			}
			if((lo <= losm && hi > losm) || (lo <  hism && hi >= hism)) {
				return true;
			}
			return false;
  }
  
  public boolean endsInFirstRow() {
	  	return (long)len_ == row_+1;
  }
  
  public double uppermostRow() {
	  	return row_ + 1 - (long)len_;
  }
  
  public long leftmostCol() {
	  return col_ + 1 - (long)len_;
  }
}
