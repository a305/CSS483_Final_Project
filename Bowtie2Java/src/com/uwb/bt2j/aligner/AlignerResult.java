package com.uwb.bt2j/aligner;

public class AlignerResult {
  public void reset() {
    ned_.clear();
	  aed_.clear();
	  score_.invalidate();
	  refcoord_.reset();
	  refival_.reset();
	  shapeSet_     = false;
	  rdlen_        = 0;
	  reflen_       = 0;
	  rdrows_       = 0;
	  rdextent_     = 0;
	  rdexrows_     = 0;
	  rfextent_     = 0;
	  refns_        = 0;
	  type_         = ALN_RES_TYPE_UNPAIRED;
	  fraglen_      = -1;
	  trimSoft_     = false;
	  trim5p_       = 0;
	  trim3p_       = 0;
	  pretrimSoft_  = true;
	  pretrim5p_    = 0;
	  pretrim3p_    = 0;
	  seedmms_      = 0; // number of mismatches allowed in seed
	  seedlen_      = 0; // length of seed
	  seedival_     = 0; // interval between seeds
	  minsc_        = 0; // minimum score
	  nuc5p_        = 0;
	  nuc3p_        = 0;
	  fraglenSet_   = false;
  }
  
  public void setShape(long id,
  long off,
  long reflen,
  Boolean fw,
  double rdlen,
  Boolean pretrimSoft,
  double pretrim5p,
  double pretrim3p) {
    rdlen_       = rdlen;
	rdrows_      = rdlen;
	refcoord_.init(id, off, fw);
	pretrimSoft_  = pretrimSoft;
	pretrim5p_    = pretrim5p;
	pretrim3p_    = pretrim3p;
	trimSoft_     = trimSoft;
	trim5p_       = trim5p;
	trim3p_       = trim3p;
	// Propagate trimming to the edits.  We assume that the pos fields of the
	// edits are set w/r/t to the rows of the dynamic programming table, and
	// haven't taken trimming into account yet.
	//
	// TODO: The division of labor between the aligner and the AlnRes is not
	// clean.  Perhaps the trimming and *all* of its side-effects should be
	// handled by the aligner.
	double trimBeg = fw ? trim5p : trim3p;
	if(trimBeg > 0) {
		for(double i = 0; i < ned_.size(); i++) {
			// Shift by trim5p, since edits are w/r/t 5p end
			ned_[i].pos -= (double)trimBeg;
		}
	}
	// Length after all soft trimming and any hard trimming that occurred
	// during alignment
	rdextent_ = rdlen;
	if(pretrimSoft_) {
		rdextent_ -= (pretrim5p + pretrim3p); // soft trim
	}
	rdextent_ -= (trim5p + trim3p); // soft or hard trim from alignment
	rdexrows_ = rdextent_;
	calcRefExtent();
	refival_.init(id, off, fw, rfextent_);
	reflen_ = reflen;
	shapeSet_ = true;
  }
  
  public void init(
  double rdLen,
  AlnScore score,
  EList<Edit> ned,
  double ned_i,
  double ned_n,
  EList<Edit> aed,
  double aed_i,
  double aed_n,
  Coord, refcoord,
  long reflen,
  int seedmms,
  int seedlen,
  int seedival,
  long minsc,
  int nuc5p,
  int nuc3p,
  Boolean pretrimSoft,
  double pretrim5p,
  double pretrim3p,
  Boolean trimSoft,
  double trim5p,
  double trim3p) {
    	rdlen_  = rdlen;
	rdrows_ = rdlen;
	score_  = score;
	ned_.clear();
	aed_.clear();
	if(ned != null) {
		for(double i = ned_i; i < ned_i + ned_n; i++) {
			ned_.push_back(ned[i]);
		}
	}
	if(aed != null) {
		for(double i = aed_i; i < aed_i + aed_n; i++) {
			aed_.push_back(aed[i]);
		}
	}
	refcoord_     = refcoord;
	reflen_       = reflen;
	seedmms_      = seedmms;
	seedlen_      = seedlen;
	seedival_     = seedival;
	minsc_        = minsc;
	nuc5p_        = nuc5p;
	nuc3p_        = nuc3p;
	pretrimSoft_  = pretrimSoft;
	pretrim5p_    = pretrim5p;
	pretrim3p_    = pretrim3p;
	trimSoft_     = trimSoft;
	trim5p_       = trim5p;
	trim3p_       = trim3p;
	rdextent_     = rdlen;      // # read characters after any hard trimming
	if(pretrimSoft) {
		rdextent_ -= (pretrim5p + pretrim3p);
	}
	if(trimSoft) {
		rdextent_ -= (trim5p + trim3p);
	}
	rdexrows_ = rdextent_;
	calcRefExtent();
	setShape(
		refcoord.ref(), // id of reference aligned to
		refcoord.off(), // offset of first aligned char into ref seq
		reflen,         // length of reference sequence aligned to
		refcoord.fw(),  // aligned to Watson strand?
		rdlen,          // length of read after hard trimming, before soft
		pretrimSoft,    // whether trimming prior to alignment was soft
		pretrim5p,      // # poss trimmed form 5p end before alignment
		pretrim3p,      // # poss trimmed form 3p end before alignment
		trimSoft,       // whether local-alignment trimming was soft
		trim5p,         // # poss trimmed form 5p end during alignment
		trim3p);        // # poss trimmed form 3p end during alignment
	shapeSet_ = true;
  }
}
