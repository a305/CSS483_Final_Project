package com.uwb.bt2j.util;

import com.uwb.bt2j.indexer.EList;
import com.uwb.bt2j.indexer.Quad;

public class Checkpointer {
	public int   perpow2_;   // 1 << perpow2_ - 2 is the # of uncheckpointed
    // anti-diags between checkpointed anti-diag pairs
	public int   per_;       // 1 << perpow2_
	public int   lomask_;    // mask for extracting low bits
	public int   nrow_;      // # rows in current rectangle
	public int   ncol_;      // # cols in current rectangle
	public long  perf_;      // perfect score
	public boolean     local_;     // local alignment?

	public int   ndiag_;     // # of double-diags

	public int   locol_;     // leftmost column committed
	public int   hicol_;     // rightmost column committed

	// Map for committing scores from vector columns to checkpointed diagonals
	public EList<Integer> commitMap_;
	public boolean firstCommit_;

	public EList<Quad> qdiag1s_; // checkpoint H/E/F values for diagonal 1
	public EList<Quad> qdiag2s_; // checkpoint H/E/F values for diagonal 2

	EList<Quad> qrows_;   // checkpoint H/E/F values for rows

	// We store columns in this way to reduce overhead of populating them
	public boolean          is8_;     // true -> fill used 8-bit cells
	public int        niter_;   // # __m128i words per column
	EList_m128i   qcols_;   // checkpoint E/F/H values for select columns

	public boolean          debug_;   // get debug checkpoints? (i.e. fill qcolsD_?)
	EList_m128i   qcolsD_;  // checkpoint E/F/H values for all columns (debug)
	
	public Checkpointer() {
		reset();
	}
	
	public void init(
			int nrow,          // # of rows
			int ncol,          // # of columns
			int perpow2,       // checkpoint every 1 << perpow2 diags (& next)
			long perfectScore, // what is a perfect score?  for sanity checks
			boolean is8,             // 8-bit?
			boolean doTri,           // triangle shaped?
			boolean local,           // is alignment local?  for sanity checks
			boolean debug)           // gather debug checkpoints?
	{
		nrow_ = nrow;
		ncol_ = ncol;
		perpow2_ = perpow2;
		per_ = 1 << perpow2;
		lomask_ = ~(0xffffffff << perpow2);
		perf_ = perfectScore;
		local_ = local;
		ndiag_ = (ncol + nrow - 1 + 1) / per_;
		locol_ = Integer.MAX_VALUE;
		hicol_ = Integer.MIN_VALUE;
//		debug_ = debug;
		debug_ = true;
		commitMap_.clear();
		firstCommit_ = true;
		int perword = (is8 ? 16 : 8);
		is8_ = is8;
		niter_ = ((nrow_ + perword - 1) / perword);
		if(doTri) {
			// Save a pair of anti-diagonals every per_ anti-diagonals for
			// backtrace purposes
			qdiag1s_.resize(ndiag_ * nrow_);
			qdiag2s_.resize(ndiag_ * nrow_);
		} else {
			// Save every per_ columns and rows for backtrace purposes
			qrows_.resize((nrow_ / per_) * ncol_);
			qcols_.resize((ncol_ / per_) * (niter_ << 2));
		}
		if(debug_) {
			// Save all columns for debug purposes
			qcolsD_.resize(ncol_ * (niter_ << 2));
		}
	}
	
	public boolean debug() {
		return debug_;
	}
	
	public long debugCell(int row, int col, int hef) {
		__m128i ptr = qcolsD_.ptr() + hef;
		// Fast forward to appropriate column
		ptr += ((col * niter_) << 2);
		int mod = row % niter_; // which m128i
		int div = row / niter_; // offset into m128i
		// Fast forward to appropriate word
		ptr += (mod << 2);
		// Extract score
		short sc = (is8_ ? ((byte)ptr)[div] : ((short)ptr)[div]);
		long asc = Long.MIN_VALUE;
		// Convert score
		if(is8_) {
			if(local_) {
				asc = sc;
			} else {
				if(sc == 0) asc = Long.MIN_VALUE;
				else asc = sc - 0xff;
			}
		} else {
			if(local_) {
				asc = sc + 0x8000;
			} else {
				if(sc != Short.MIN_VALUE) asc = sc - 0x7fff;
			}
		}
		return asc;
	}
	
	public boolean isCheckpointed(int row, int col) {
		int mod = (row + col) & lomask_;
		return mod >= per_ - 2;
	}
	
	public long scoreTriangle(int row, int col, int ref) {
		boolean diag1 = ((row + col) & lomask_) == per_ - 2;
		int off = (row + col) >> perpow2_;
		if(diag1) {
			if(qdiag1s_[off * nrow_ + row].sc[hef] == Short.MIN_VALUE) {
				return Long.MIN_VALUE;
			} else {
				return qdiag1s_[off * nrow_ + row].sc[hef];
			}
		} else {
			if(qdiag2s_[off * nrow_ + row].sc[hef] == Short.MIN_VALUE) {
				return Long.MIN_VALUE;
			} else {
				return qdiag2s_[off * nrow_ + row].sc[hef];
			}
		}
	}
	
	public long scoreSquare(int row, int col, int ref) {
		// Is it in a checkpointed row?  Note that checkpointed rows don't
				// necessarily have the horizontal contributions calculated, so we want
				// to use the column info in that case.
				if((row & lomask_) == lomask_ && hef != 1) {
					long sc = qrows_[(row >> perpow2_) * ncol_ + col].sc[hef];
					if(sc == Short.MIN_VALUE) return Long.MIN_VALUE;
					return sc;
				}
				hef--;
				if(hef == -1) hef = 2;
				// Fast forward to appropriate column
				__m128i* ptr = qcols_.ptr() + hef;
				ptr += (((col >> perpow2_) * niter_) << 2);
				int mod = row % niter_; // which m128i
				int div = row / niter_; // offset into m128i
				// Fast forward to appropriate word
				ptr += (mod << 2);
				// Extract score
				short sc = (is8_ ? ((byte)ptr)[div] : ((short)ptr)[div]);
				long asc = Long.MIN_VALUE;
				// Convert score
				if(is8_) {
					if(local_) {
						asc = sc;
					} else {
						if(sc == 0) asc = Long.MIN_VALUE;
						else asc = sc - 0xff;
					}
				} else {
					if(local_) {
						asc = sc + 0x8000;
					} else {
						if(sc != Short.MIN_VALUE) asc = sc - 0x7fff;
					}
				}
				return asc;
	}
	
	public void commitCol(__m128i pvH, __m128i pvE, __m128i pvF, int coli) {
		
	}
	
	public boolean inited() {
		return nrow_ > 0;
	}
	
	public void reset() {
		perpow2_ = per_ = lomask_ = nrow_ = ncol_ = 0;
		local_ = false;
		niter_ = ndiag_ = locol_ = hicol_ = 0;
		perf_ = 0;
		firstCommit_ = true;
		is8_ = debug_ = false;
	}
	
	public int per()     { return per_;     }
	public int perpow2() { return perpow2_; }
	public int lomask()  { return lomask_;  }
	public int locol()   { return locol_;   }
	public int hicol()   { return hicol_;   }
	public int nrow()    { return nrow_;    }
	public int ncol()    { return ncol_;    }
	
	Quad qdiag1sPtr() { return qdiag1s_.ptr(); }
	Quad qdiag2sPtr() { return qdiag2s_.ptr(); }
}
