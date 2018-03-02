package com.uwb.bt2j.aligner;

import com.uwb.bt2j.util.EList;

public class SSEMatrix {
	public static final double E = 0;
	public static final double F = 1;
	public static final double H = 2;
	public static final double TMP = 3;
	public boolean             inited_;      // initialized?
	public double           nrow_;        // # rows
	public double           ncol_;        // # columns
	public double           nvecrow_;     // # vector rows (<= nrow_)
	public double           nveccol_;     // # vector columns (<= ncol_)
	public double           wperv_;       // # words per vector
	public double           vecshift_;    // # bits to shift to divide by words per vec
	public double           nvecPerCol_;  // # vectors per column
	public double           nvecPerCell_; // # vectors per matrix cell (4)
	public double           colstride_;   // # vectors b/t adjacent cells in same row
	public double           rowstride_;   // # vectors b/t adjacent cells in same col
	public EList_m128i      matbuf_;      // buffer for holding vectors
	public ELList<Short> masks_;       // buffer for masks/backtracking flags
	public EList<Boolean>      reset_;       // true iff row in masks_ has been reset
	
	public SSEMatrix(int cat) {
		nvecPerCell_ = 4;
		matbuf_ = cat;
	}
	
	public void init(double nrow,double ncol,double wperv) {
		nrow_ = nrow;
		ncol_ = ncol;
		wperv_ = wperv;
		nvecPerCol_ = (nrow + (wperv-1)) / wperv;
		// The +1 is so that we don't have to special-case the final column;
		// instead, we just write off the end of the useful part of the table
		// with pvEStore.
		try {
			matbuf_.resizeNoCopy((ncol+1) * nvecPerCell_ * nvecPerCol_);
		} catch(Exception e) {
			System.err.println( "Tried to allocate DP matrix with " + (ncol+1)
			     + " columns, " + nvecPerCol_
				 + " vectors per column, and and " + nvecPerCell_
				 + " vectors per cell" );
		}
		vecshift_ = (wperv_ == 8) ? 3 : 4;
		nvecrow_ = (nrow + (wperv_-1)) >> vecshift_;
		nveccol_ = ncol;
		colstride_ = nvecPerCol_ * nvecPerCell_;
		rowstride_ = nvecPerCell_;
		inited_ = true;
	}
	
	public double colstride() {
		return colstride_;
	}
	
	public double rowstride(){
		return rowstride_;
	}
	
	public int eltSlow(double row, double col, double mat) {
		// Move to beginning of column/row
		double rowelt = row / nvecrow_;
		double rowvec = row % nvecrow_;
		double eltvec = (col * colstride_) + (rowvec * rowstride_) + mat;
		if(wperv_ == 16) {
			return (int)((matbuf_.ptr() + eltvec))[rowelt];
		} else {
			return (int)((matbuf_.ptr() + eltvec))[rowelt];
		}
	}
	
	public int elt(double row, double col, double mat) {
		// Move to beginning of column/row
				double rowelt = row / nvecrow_;
				double rowvec = row % nvecrow_;
				double eltvec = (col * colstride_) + (rowvec * rowstride_) + mat;
				if(wperv_ == 16) {
					return (int)((matbuf_.ptr() + eltvec))[rowelt];
				} else {
					return (int)((matbuf_.ptr() + eltvec))[rowelt];
				}
	}
	
	public int eelt(double row, double col) {
		return elt(row, col, E);
	}
	
	public int felt(double row, double col) {
		return elt(row, col, F);
	}
	
	public int helt(double row, double col) {
		return elt(row, col, H);
	}
	
	public boolean reportedThrough(double row, double col) {
		return (masks_[row][col] & (1 << 0)) != 0;
	}
	
	public void setReportedThrough(double row, double col) {
		masks_[row][col] |= (1 << 0);
	}
	
	public boolean isHMaskSet(double row, double col) {
		
	}
	
	public boolean hMaskSet(double row, double col, int mask) {
		
	}
	
	public boolean isEMaskSet(double row, double col) {
		
	}
	
	public boolean eMaskSet(double row, double col, int mask) {
		
	}
	
	public boolean isFMaskSet(double row, double col) {
		
	}
	
	public boolean fMaskSet(double row, double col, int mask) {
		
	}
	
	public void analyzeCell(
			double row,          // current row
			double col,          // current column
			double ct,           // current cell type: E/F/H
			int refc,
			int readc,
			int readq,
			Scoring sc,   // scoring scheme
			long offsetsc,    // offset to add to each score
			RandomSource rand,  // rand gen for choosing among equal options
			boolean empty,         // out: =true iff no way to backtrace
			int cur,            // out: =type of transition
			boolean branch,        // out: =true iff we chose among >1 options
			boolean canMoveThru,   // out: =true iff ...
			boolean reportedThru) // out: =true iff ...
	{
		
	}
	
	public void initMasks() {
		masks_.resize(nrow_);
		reset_.resizeNoCopy(nrow_);
		reset_.fill(false);
	}
	
	public double nrow() {
		return nrow_;
	}
	
	public double ncol() {
		return ncol_;
	}
	
	public void resetRow(double i) {
		masks_[i].resizeNoCopy(ncol_);
		masks_[i].fillZero();
		reset_[i] = true;
	}
}
