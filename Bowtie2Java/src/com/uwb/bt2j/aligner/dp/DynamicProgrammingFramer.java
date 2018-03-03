package com.uwb.bt2j.aligner.dp;
class DynamicProgrammingFramer {
  protected boolean trimToRef_;
  
  public DynamicProgrammingFramer(boolean trimToRef) {
	  trimToRef_ = trimToRef;
  }
  
  public boolean frameSeedExtensionRect(
		  long  off,      // ref offset implied by seed hit assuming no gaps
			double   rdlen,    // length of read sequence used in DP table
			long  reflen,   // length of reference sequence aligned to
			double   maxrdgap, // max # of read gaps permitted in opp mate alignment
			double   maxrfgap, // max # of ref gaps permitted in opp mate alignment
			long  maxns,    // # Ns permitted
			double   maxhalf,  // max width in either direction
			DPRect  rect) {
	// Set N, the maximum number of reference or read gaps permitted, whichever
		// is larger.  Also, enforce ceiling: can't be larger than 'maxhalf'.
		double maxgap = Math.max(maxrdgap, maxrfgap);
		maxgap = Math.min(maxgap, maxhalf);
		// Leave room for "LHS gap" and "LHS extra" diagonals
		double refl = off - 2 * maxgap;               // inclusive
		// Leave room for "RHS gap" and "RHS extra" diagonals
		double refr = off + (rdlen - 1) + 2 * maxgap; // inclusive
		double triml = 0, trimr = 0;
		// Check if we have to trim to fit the extents of the reference
		if(trimToRef_) {
			maxns = 0; // no leeway
		} else if(maxns == (long)rdlen) {
			maxns--;
		}
		// Trim from RHS of rectangle
		if(refr >= reflen + maxns) {
			trimr = (double)(refr - (reflen + maxns - 1));
		}
		// Trim from LHS of rectangle
		if(refl < -maxns) {
			triml = (double)(-refl) - (double)maxns;
		}
		rect.refl_pretrim = refl;
		rect.refr_pretrim = refr;
		rect.refl  = refl + triml;
		rect.refr  = refr - trimr;
		rect.triml = triml;
		rect.trimr = trimr;
		rect.maxgap = maxgap;
		// Remember which diagonals are "core" as offsets from the LHS of the
		// untrimmed rectangle
		rect.corel = maxgap;
		rect.corer = rect.corel + 2 * maxgap; // inclusive
		return !rect.entirelyTrimmed();
  }
  
  public boolean frameFindMateRect(
		  boolean anchorLeft,  // true iff anchor alignment is to the left
			long ll,       // leftmost Watson off for LHS of opp alignment
			long lr,       // rightmost Watson off for LHS of opp alignment
			long rl,       // leftmost Watson off for RHS of opp alignment
			long rr,       // rightmost Watson off for RHS of opp alignment
			double  rdlen,    // length of opposite mate
			long reflen,   // length of reference sequence aligned to
			double  maxrdgap, // max # of read gaps permitted in opp mate alignment
			double  maxrfgap, // max # of ref gaps permitted in opp mate alignment
			long maxns,    // max # Ns permitted
			double  maxhalf,  // max width in either direction
			DPRect rect  
		  ) {
		if(anchorLeft) {
			return frameFindMateAnchorLeftRect(
				ll,
				lr,
				rl,
				rr,
				rdlen,
				reflen,
				maxrdgap,
				maxrfgap,
				maxns,
				maxhalf,
				rect);
		} else {
			return frameFindMateAnchorRightRect(
				ll,
				lr,
				rl,
				rr,
				rdlen,
				reflen,
				maxrdgap,
				maxrfgap,
				maxns,
				maxhalf,
				rect);
		}
  }
  
  public Boolean frameFindMateAnchorLeftRect(
		  long ll,       // leftmost Watson off for LHS of opp alignment
			long lr,       // rightmost Watson off for LHS of opp alignment
			long rl,       // leftmost Watson off for RHS of opp alignment
			long rr,       // rightmost Watson off for RHS of opp alignment
			double  rdlen,    // length of opposite mate
			long reflen,   // length of reference sequence aligned to
			double  maxrdgap, // max # of read gaps permitted in opp mate alignment
			double  maxrfgap, // max # of ref gaps permitted in opp mate alignment
			long maxns,    // max # Ns permitted in alignment
			double  maxhalf,  // max width in either direction
			DPRect rect)     // out: DP rectangle
  {
	  double triml = 0, trimr = 0;
		double maxgap = Math.max(maxrdgap, maxrfgap);
		maxgap = Math.max(maxgap, maxhalf);
		// Amount of padding we have to add to account for the fact that alignments
		// ending between en_left/en_right might start in various columns in the
		// first row
		double pad_left = maxgap;
		double pad_right = maxgap;
		long en_left  = rl;
		long en_right = rr;
		double st_left  = en_left - (rdlen-1);
		double en_right_pad = en_right + pad_right;
		double st_left_pad  = st_left  - pad_left;
		double refl = st_left_pad;
		double refr = en_right_pad;
		if(trimToRef_) {
			maxns = 0;
		} else if(maxns == (long)rdlen) {
			maxns--;
		}
		// Trim from the RHS of the rectangle?
		if(refr >= reflen + maxns) {
			trimr = (double)(refr - (reflen + maxns - 1));
		}
		// Trim from the LHS of the rectangle?
		if(refl < -maxns) {
			triml = (double)(-refl) - (double)maxns;
		}
		double width = (double)(refr - refl + 1);
		rect.refl_pretrim = refl;
		rect.refr_pretrim = refr;
		rect.refl  = refl + triml;
		rect.refr  = refr - trimr;
		rect.triml = triml;
		rect.trimr = trimr;
		rect.maxgap = maxgap;
		rect.corel = maxgap;
		rect.corer = width - maxgap - 1; // inclusive
		return !rect.entirelyTrimmed();
  }
  
  public Boolean frameFindMateAnchorRightRect(
		  long ll,       // leftmost Watson off for LHS of opp alignment
			long lr,       // rightmost Watson off for LHS of opp alignment
			long rl,       // leftmost Watson off for RHS of opp alignment
			long rr,       // rightmost Watson off for RHS of opp alignment
			double  rdlen,    // length of opposite mate
			long reflen,   // length of reference sequence aligned to
			double  maxrdgap, // max # of read gaps permitted in opp mate alignment
			double  maxrfgap, // max # of ref gaps permitted in opp mate alignment
			long maxns,    // max # Ns permitted in alignment
			double  maxhalf,  // max width in either direction
			DPRect rect)     // out: DP rectangle
  {
	  double triml = 0, trimr = 0;
		double maxgap = Math.max(maxrdgap, maxrfgap);
		maxgap = Math.max(maxgap, maxhalf);
		double pad_left = maxgap;
		double pad_right = maxgap;
		long st_left = ll;
		long st_right = lr;
		double en_right = st_right + (rdlen-1);
		double en_right_pad = en_right + pad_right;
		double st_left_pad  = st_left  - pad_left;
		// We have enough info to deduce where the boundaries of our rectangle
		// should be.  Finalize the boundaries, ignoring reference trimming for now
		double refl = st_left_pad;
		double refr = en_right_pad;
		if(trimToRef_) {
			maxns = 0;
		} else if(maxns == (long)rdlen) {
			maxns--;
		}
		// Trim from the RHS of the rectangle?
		if(refr >= reflen + maxns) {
			trimr = (double)(refr - (reflen + maxns - 1));
		}
		// Trim from the LHS of the rectangle?
		if(refl < -maxns) {
			triml = (double)(-refl) - (double)maxns;
		}
		double width = (double)(refr - refl + 1);
		rect.refl_pretrim = refl;
		rect.refr_pretrim = refr;
		rect.refl  = refl + triml;
		rect.refr  = refr - trimr;
		rect.triml = triml;
		rect.trimr = trimr;
		rect.maxgap = maxgap;
		rect.corel = maxgap;
		rect.corer = width - maxgap - 1; // inclusive
		return !rect.entirelyTrimmed();
  }
  
  protected void trimToRef(
			double   reflen,  // in: length of reference sequence aligned to
			long refl,    // in/out: ref pos of upper LHS of parallelogram
			long refr,    // in/out: ref pos of lower RHS of parallelogram
			double  trimup,  // out: number of bases trimmed from upstream end
			double  trimdn)  // out: number of bases trimmed from downstream end
  {
		if(refl < 0) {
			trimup = (double)(-refl);
			//refl = 0;
		}
		if(refr >= (long)reflen) {
			trimdn = (double)(refr - reflen + 1);
			//refr = (long)reflen-1;
		}
  }
}
