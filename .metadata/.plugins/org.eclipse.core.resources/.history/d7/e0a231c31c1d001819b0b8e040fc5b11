package com.uwb.bt2j/aligner;

public class AlignerResult {
	protected Boolean shapeSet_;
	protected Boolean fraglenSet_;
	protected Boolean pretrimSoft_;
	protected Boolean trimSoft_;
	protected double pretrim5p_;
	protected double pretrim3p_;
	protected double trim5p_;
	protected double trim3p_;
	protected double rdlen_;
	protected double rdrows_;
	protected double rdextent_;
	protected double rdexrows_;
	protected double rfextent_;
	protected double refns_;
	protected int seedmms_;
	protected int seedlen_;
	protected int seedival_;
	protected int nuc5p_;
	protected int nuc3p_;
	protected int type_;
	protected AlignmentScore score_;
	protected AlignmentScore oscore_;
	protected EList<Edit> ned_;
	protected EList<Edit> aed_;
	protected Coord refcoord_;
	protected long reflen_;
	protected Interval refival_;
	protected long minsc_;
	protected long fraglen_;
	
public enum ResultType {
	ALN_RES_TYPE_UNPAIRED(1),   // unpaired alignment
	ALN_RES_TYPE_UNPAIRED_MATE1(2), // mate #1 in pair, aligned unpaired
	ALN_RES_TYPE_UNPAIRED_MATE2(3), // mate #2 in pair, aligned unpaired
	ALN_RES_TYPE_MATE1(4),          // mate #1 in paired-end alignment
	ALN_RES_TYPE_MATE2(5);           // mate #2 in paired-end alignment
	private int x;
	ResultType(int y) {x = y;}
}

	public AlignerResult() {
		ned_ = 7;
		aed_ = 7;
		reset();	
	}
	
	public void reverseEdits() {
		ned_.reverse();
		aed_.reverse();
	}
	
	public void invertEdits() {
		Edit.invertPoss(ned_, rdexrows_, false);
		Edit.invertPoss(aed_, rdexrows_, false);
	}
	
	public Boolean empty() {
		if(!VALID_AL_SCORE(score_)) {
			return true;
		} else {
			return false;
		}	
	}
	
	public long refid() {
		return refcoord_.ref();	
	}
	
	public int orient() {
		return refcoord_.orient();	
	}
	
	public long refoff() {
		return refcoord_.off();	
	}
	
	public void getCoords(Coord st, Coord en) {
		st.init(refcoord_);
		en.init(refcoord_);
		en.adjustOff(refExtent() - 1);
	}
	
	public void getExtendedCoords(Coord st, Coord en, AlignerFlags flags) {
		getCoords(st, en);
		// Take trimming into account
		if (!flags.scUnMapped()) {
			long trim_st  = (fw() ? trim5p_ : trim3p_);
			long trim_en  = (fw() ? trim3p_ : trim5p_);
			trim_st += (fw() ? pretrim5p_ : pretrim3p_);
			trim_en += (fw() ? pretrim3p_ : pretrim5p_);
			st.adjustOff(-trim_st);
			en.adjustOff( trim_en);
		}	
	}
	
	public Boolean within(long id, long off, Boolean fw, double extent) {
		if(refcoord_.ref() == id &&
		   refcoord_.off() >= off &&
		   refcoord_.off() + refExtent() <= off + extent &&
		   refcoord_.fw() == fw)
		{
			return true;
		}
		return false;
	}
	
	public void setScore(AlignmentScore score) {
		score_ = score;
	}
	
	public void setNucs(Boolean fw, int nup, int ndn) {
		nuc5p_ = fw ? nup : ndn;
		nuc3p_ = fw ? ndn : nup;
	}
	
	public Coord refcoords() {
		return refcoord_;	
	}
	
	public Interval refival() {
		return refival_;	
	}
	
	public Boolean fw() {
		return refcoord_.fw();	
	}
	
	public AlignmentScore score() {
		return score_;	
	}
	
	public AlignmentScore oscore() {
		return oscore_;	
	}
	
	public EList<Edit> ned() {
		return ned_;	
	}
	
	public EList<Edit> aed() {
		return aed_;	
	}
	
	public double readExtent() {
		return rdextent_;	
	}
	
	public double readExtentRow() {
		return rdextrows_;	
	}
	
	public double readLength() {
		return rdlen_;	
	}
	
	public double refExtent() {
		return rfextent_;	
	}
	
	public long reflen() {
		return reflen_;
	}
	
	public double refNucExtent() {
		return rfextent_;	
	}
	
	public void printStacked(Read rd, OutputStream o) {
		printStacked(refcoord_.fw() ? rd.patFw : rd.patRc, o);
	}
	
	public void printStacked(BTDnaString seq, OutputStream o) {
		Edit.printQAlign(o, seq, ned_);
		// Print reference offset below reference string
		o.write("^" + "\n");
		o.write("(" + refcoord_.ref() + "," + refcoord_.off() + ")" + "\n");
	}
	
	public void setParams(int seedmms, int seedlen, int seedival) {
		seedmms_ = seedmms;
		seedlen_ = seedlen;
		seedival_ = seedival;
		minsc_ = minsc;
	}
	
	public int seedmms() {
		return seedmms_;	
	}
	
	public int seedlen() {
		return seedlen_;	
	}
	
	public int seedival() {
		return seedival;	
	}
	
	public long minScore() {
		return minsc_;	
	}
	
	public Boolean trimmedRow5p(double i) {
		return i < trim5p_ || rdrows_ - i - 1 < trim3p_;
	}
	
	public Boolean trimmedPos5p(double i) {
		return i < trim5p_ || rdlen_ - i - 1 < trim3p_;
	}
	
	public Boolean alignedRow5p(double i) {
		return !trimmedRow5p(i);
	}
	
	public Boolean alignedPos5p(double i) {
		return !trimmedPos5p(i);
	}
	
	public Boolean readUnpaired() {
		return type_ == ALN_RES_TYPE_UNPAIRED;
	}
	
	public Boolean alignedPaired() {
		return type_ == ALN_RES_TYPE_MATE1 ||
		       type_ == ALN_RES_TYPE_MATE2;	
	}
	
	public Boolean readMate1() {
		return type_ == ALN_RES_TYPE_MATE1 ||
		       type_ == ALN_RES_TYPE_UNPAIRED_MATE1;
	}
	
	public Boolean alignedMate1() {
		return type_ == ALN_RES_TYPE_MATE1;	
	}
	
	public Boolean readMate2() {
		return type_ == ALN_RES_TYPE_MATE2 ||
		       type_ == ALN_RES_TYPE_UNPAIRED_MATE2;
	}
	
	public Boolean alignedMate2() {
		return type_ == ALN_RES_TYPE_MATE2;	
	}
	
	public void setMateParams(int type, AlignmentResult omate, AlignerFlags flags) {
		type_ = type;
		fraglen_ = 0;
		if(omate != null) {
			oscore_ = omate.score_;
			// When should we calculate a fragment length here?  There are a
			// couple reasonable ideas:
			// 1. When mates align concordantly
			// 2. When both mates align to the same reference string
			// BWA seems to do 2., so that's what we'll do here.
			Boolean sameChr = true;
			if((sameChr && refcoord_.ref() == omate.refcoord_.ref()) ||
			   flags.alignedConcordant())
			{
				setFragmentLength(omate, flags);
			}
		}
	}
	
	public long setFragmentLength(AlignmentResult omate, AlignerFlags flags) {
		Coord st, en;
		Coord ost, oen;
		getExtendedCoords(st, en, flags);
		omate.getExtendedCoords(ost, oen, flags);
		Boolean imUpstream = st.off() < ost.off();
		long up = Math.min(st.off(), ost.off());
		long dn = Math.max(en.off(), oen.off());
		fraglen_ = 1 + dn - up;
		if(!imUpstream) {
			fraglen_ = -fraglen_;
		}
		fraglenSet_ = true;
		return fraglen_;
	}
	
	public long fragmentLength() {
		return fraglen_;	
	}
	
	public double trimmed5p(Boolean soft) {
		double trim = 0;
		if(pretrimSoft_ == soft) trim += pretrim5p_;
		if(trimSoft_ == soft) trim += trim5p_;
		return trim;
	}
	
	public double trimmed3p(Boolean soft) {
		double trim = 0;
		if(pretrimSoft_ == soft) trim += pretrim3p_;
		if(trimSoft_ == soft) trim += trim3p_;
		return trim;
	}
	
	public double trimmedLeft(Boolean soft) {
		return fw() ? trimmed5p(soft) : trimmed3p(soft);	
	}
	
	public double trimmedRight(Boolean soft) {
		return fw() ? trimmed3p(soft) : trimmed3p(soft);	
	}
	
	public void setRefNs(double refns) {
		refns_ = refns;
	}
	
	public double refNs() {
		return refns_;	
	}
	
	public void initStacked(Read rd, StackedAln st) {
		double trimLS = trimmed5p(true);
		double trimLH = trimmed5p(false);
		double trimRS = trimmed3p(true);
		double trimRH = trimmed3p(false);
		double len_trimmed = rd.length() - trimLS - trimRS;
		if(!fw()) {
			Edit.invertPoss(const_cast<EList<Edit>&>(ned_), len_trimmed, false);
			swap(trimLS, trimRS);
			swap(trimLH, trimRH);
		}
		st.init(
			fw() ? rd.patFw : rd.patRc,
			ned_, trimLS, trimLH, trimRS, trimRH);
		if(!fw()) {
			Edit.invertPoss(const_cast<EList<Edit>&>(ned_), len_trimmed, false);
		}
	}
	
	protected void calcRefExtent() {
		rfextent_ = rdextent_;
		for(double i = 0; i < ned_.size(); i++) {
			if(ned_[i].isRefGap()) rfextent_--;
			if(ned_[i].isReadGap()) rfextent_++;
		}	
	}
	
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
	
	public void clipLeft(double rd_amt, double rf_amt) {
		if(fw()) {
		trim5p_ += rd_amt;
		Edit.clipLo(ned_, rdexrows_, rd_amt);
		Edit.clipLo(aed_, rdexrows_, rd_amt);
	} else {
		trim3p_ += rd_amt;
		Edit.clipHi(ned_, rdexrows_, rd_amt);
		Edit.clipHi(aed_, rdexrows_, rd_amt);
	}
	rdexrows_ -= rd_amt;
	rdextent_ -= rd_amt;
	rfextent_ -= rf_amt;
	refcoord_.adjustOff(rf_amt);
	refival_.adjustOff(rf_amt);	
	}
	
	public void clipRight(double rd_amt, double rf_amt) {
		if(fw()) {
		trim3p_ += rd_amt;
		Edit.clipHi(ned_, rdexrows_, rd_amt);
		Edit.clipHi(aed_, rdexrows_, rd_amt);
	} else {
		trim5p_ += rd_amt;
		Edit.clipLo(ned_, rdexrows_, rd_amt);
		Edit.clipLo(aed_, rdexrows_, rd_amt);
	}
	rdexrows_ -= rd_amt;
	rdextent_ -= rd_amt;
	rfextent_ -= rf_amt;
	}
	
	public void clipOutside(Boolean soft, long refi, long reff) {
			// Overhang on LHS
	long left = refcoord_.off();
	if(left < refi) {
		double rf_amt = (double)(refi - left);
		double rf_i = rf_amt;
		double nedsz = ned_.size();
		if(!fw()) {
			Edit.invertPoss(ned_, rdexrows_, false);
		}
		for(double i = 0; i < nedsz; i++) {
			if(ned_[i].pos > rf_i) break;
			if(ned_[i].isRefGap()) rf_i++;
		}
		if(!fw()) {
			Edit.invertPoss(ned_, rdexrows_, false);
		}
		clipLeft(rf_i, rf_amt);
	}
	// Overhang on RHS
	long right = refcoord_.off() + refNucExtent();
	if(right > reff) {
		double rf_amt = (double)(right - reff);
		double rf_i = rf_amt;
		double nedsz = ned_.size();
		if(fw()) {
			Edit.invertPoss(ned_, rdexrows_, false);
		}
		for(double i = 0; i < nedsz; i++) {
			if(ned_[i].pos > rf_i) break;
			if(ned_[i].isRefGap()) rf_i++;
		}
		if(fw()) {
			Edit.invertPoss(ned_, rdexrows_, false);
		}
		clipRight(rf_i, rf_amt);
	}
	}
	
	public Boolean overlap(AlignerResult res) {
		if(fw() != res.fw() || refid() != res.refid()) {
		// Must be same reference and same strand in order to overlap
		return false;
	}
	long my_left     = refoff();     // my leftmost aligned char
	long other_left  = res.refoff(); // other leftmost aligned char
	long my_right    = my_left    + refExtent();
	long other_right = other_left + res.refExtent();
	if(my_right < other_left || other_right < my_left) {
		// The rectangular hulls of the two alignments don't overlap, so
		// they can't overlap at any cell
		return false;
	}
	// Reference and strand are the same and hulls overlap.  Now go read
	// position by read position testing if any align identically with the
	// reference.
	
	// Edits are ordered and indexed from 5' to 3' to start with.  We
	// reorder them to go from left to right along the Watson strand.
	if(!fw()) {
		invertEdits();
	}
	if(!res.fw()) {
		res.invertEdits();
	}
	double nedidx = 0, onedidx = 0;
	Boolean olap = false;
	// For each row, going left to right along Watson reference strand...
	for(double i = 0; i < rdexrows_; i++) {
		double diff = 1;  // amount to shift to right for next round
		double odiff = 1; // amount to shift to right for next round
		// Unless there are insertions before the next position, we say
		// that there is one cell in this row involved in the alignment
		my_right = my_left + 1;
		other_right = other_left + 1;
		while(nedidx < ned_.size() && ned_[nedidx].pos == i) {
			if(ned_[nedidx].isRefGap()) {
				// Next my_left will be in same column as this round
				diff = 0;
			}
			nedidx++;
		}
		while(onedidx < res.ned_.size() && res.ned_[onedidx].pos == i) {
			if(res.ned_[onedidx].isRefGap()) {
				// Next my_left will be in same column as this round
				odiff = 0;
			}
			onedidx++;
		}
		if(i < rdexrows_ - 1) {
			// See how many inserts there are before the next read
			// character
			double nedidx_next  = nedidx;
			double onedidx_next = onedidx;
			while(nedidx_next < ned_.size() &&
				  ned_[nedidx_next].pos == i+1)
			{
				if(ned_[nedidx_next].isReadGap()) {
					my_right++;
				}
				nedidx_next++;
			}
			while(onedidx_next < res.ned_.size() &&
				  res.ned_[onedidx_next].pos == i+1)
			{
				if(res.ned_[onedidx_next].isReadGap()) {
					other_right++;
				}
				onedidx_next++;
			}
		}
		// Contained?
		olap =
			(my_left >= other_left && my_right <= other_right) ||
			(other_left >= my_left && other_right <= my_right);
		// Overlapping but not contained?
		if(!olap) {
			olap =
				(my_left <= other_left && my_right > other_left) ||
				(other_left <= my_left && other_right > my_left);
		}
		if(olap) {
			break;
		}
		// How to do adjust my_left and my_right
		my_left = my_right + diff - 1;
		other_left = other_right + odiff - 1;
	}
	if(!fw()) {
		invertEdits();
	}
	if(!res.fw()) {
		res.invertEdits();
	}
	return olap;
	}
	
	public Boolean matchesRef(
	Read rd,
	BitPairReference ref,
	BTDnaString rf,
	BTDnaString rdseq,
	BTString qseq,
	SStringExpandable<char> raw_refbuf,
	SStringExpandable<double> destU32,
	EList<Boolean> matches) {
		Boolean fw = refcoord_.fw();
	// Adjust reference string length according to edits
	double refallen = refNucExtent();
	raw_refbuf.resize(refallen + 16);
	raw_refbuf.clear();
	int nsOnLeft = 0;
	if(refcoord_.off() < 0) {
		nsOnLeft = -((int)refcoord_.off());
	}
	int off = ref.getStretch(
		reinterpret_cast<double>(raw_refbuf.wbuf()),
		(double)refcoord_.ref(),
		(double)max<TRefOff>(refcoord_.off(), 0),
		refallen,
		destU32);
	String refbuf = raw_refbuf.wbuf() + off;
	double trim5 = 0, trim3 = 0;
	if(trimSoft_) {
		trim5 += trim5p_;
		trim3 += trim3p_;
	}
	if(pretrimSoft_) {
		trim5 += pretrim5p_;
		trim3 += pretrim3p_;
	}
	rf.clear();
	rdseq.clear();
	rdseq = rd.patFw;
	if(!fw) {
		rdseq.reverseComp();
	}
	// rdseq is the nucleotide sequence from upstream to downstream on the
	// Watson strand.  ned_ are the nucleotide edits from upstream to
	// downstream.  rf contains the reference characters.
	Edit.toRef(rdseq, ned_, rf, fw, trim5, trim3);
	matches.clear();
	Boolean matchesOverall = true;
	matches.resize(refallen);
	matches.fill(true);
	for(double i = 0; i < refallen; i++) {
		if((int)i < nsOnLeft) {
			if((int)rf[i] != 4) {
				matches[i] = false;
				matchesOverall = false;
			}
		} else {
			if((int)rf[i] != (int)refbuf[i-nsOnLeft]) {
				matches[i] = false;
				matchesOverall = false;
			}
		}
	}
	if(!matchesOverall) {
		// Print a friendly message showing the difference between the
		// reference sequence obtained with Edit::toRef and the actual
		// reference sequence
		System.err.println();
		Edit.printQAlignNoCheck(
			cerr,
			"    ",
			rdseq,
			ned_);
		System.err.print( "    ");
		for(double i = 0; i < refallen; i++) {
			System.err.print( (matches[i] ? " " : "*"));
		}
		System.err.println();
		System.err.print( "    ");
		for(double i = 0; i < refallen-nsOnLeft; i++) {
			System.err.print( "ACGTN"[(int)refbuf[i]]);
		}
		System.err.println();
		Edit.printQAlign(
			cerr,
			"    ",
			rdseq,
			ned_);
		System.err.println();
	}
	return matchesOverall;
	}
	
	public void printSeq(Read rd, BTDnaString dns, BTString o) {
		double len = dns.length();
	double st = 0;
	double en = len;
	for(double i = st; i < en; i++) {
		int c = dns.get(i);
		o.append("ACGT"[c]);
	}
	}
	
	public void printQuals(Read rd, BTString dqs, BTString 0) {
		double len = dqs.length();
	// Print decoded qualities from upstream to downstream Watson
	for(double i = 1; i < len-1; i++) {
		o.append(dqs.get(i));
	}
	}
}
