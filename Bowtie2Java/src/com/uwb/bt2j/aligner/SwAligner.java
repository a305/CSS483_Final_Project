package com.uwb.bt2j.aligner;

import java.io.OutputStream;

import com.uwb.bt2j.aligner.dp.BtBranchTracer;
import com.uwb.bt2j.util.strings.BTDnaString;
import com.uwb.bt2j.util.strings.BTString;

public class SwAligner {
	protected BTDnaString  rd_;     // read sequence
	protected BTString     qu_;     // read qualities
	protected BTDnaString  rdfw_;   // read sequence for fw read
	protected BTDnaString  rdrc_;   // read sequence for rc read
	protected BTString     qufw_;   // read qualities for fw read
	protected BTString     qurc_;   // read qualities for rc read
	protected long            rdi_;    // offset of first read char to align
	protected long            rdf_;    // offset of last read char to align
	protected boolean                fw_;     // true iff read sequence is original fw read
	protected long              refidx_; // id of reference aligned against
	protected long             reflen_; // length of entire reference sequence
	protected DPRect       rect_;   // DP rectangle
	protected char               rf_;     // reference sequence
	protected long             rfi_;    // offset of first ref char to align to
	protected long             rff_;    // offset of last ref char to align to (excl)
	protected double              rdgap_;  // max # gaps in read
	protected double              rfgap_;  // max # gaps in reference
	protected boolean                enable8_;// enable 8-bit sse
	protected boolean                extend_; // true iff this is a seed-extend problem
	protected Scoring      sc_;     // penalties for edit types
	protected long            minsc_;  // penalty ceiling for valid alignments
	protected int                 nceil_;  // max # Ns allowed in ref portion of aln

	protected boolean                sse8succ_;  // whether 8-bit worked
	protected boolean                sse16succ_; // whether 16-bit worked
	protected SSEData             sseU8fw_;   // buf for fw query, 8-bit score
	protected SSEData             sseU8rc_;   // buf for rc query, 8-bit score
	protected SSEData             sseI16fw_;  // buf for fw query, 16-bit score
	protected SSEData             sseI16rc_;  // buf for rc query, 16-bit score
	protected boolean                sseU8fwBuilt_;   // built fw query profile, 8-bit score
	protected boolean                sseU8rcBuilt_;   // built rc query profile, 8-bit score
	protected boolean                sseI16fwBuilt_;  // built fw query profile, 16-bit score
	protected boolean                sseI16rcBuilt_;  // built rc query profile, 16-bit score

	protected SSEMetrics			sseU8ExtendMet_;
	protected SSEMetrics			sseU8MateMet_;
	protected SSEMetrics			sseI16ExtendMet_;
	protected SSEMetrics			sseI16MateMet_;

	protected int                 state_;        // state
	protected boolean                initedRead_;   // true iff initialized with initRead
	protected boolean                readSse16_;    // true . sse16 from now on for read
	protected boolean                initedRef_;    // true iff initialized with initRef
	protected EList<Integer>     rfwbuf_;       // buffer for wordized ref stretches
	
	protected EList<DpNucFrame>    btnstack_;    // backtrace stack for nucleotides
	protected EList<SizeTPair>     btcells_;     // cells involved in current backtrace

	protected NBest<DpBtCandidate> btdiag_;      // per-diagonal backtrace candidates
	protected EList<DpBtCandidate> btncand_;     // cells we might backtrace from
	protected EList<DpBtCandidate> btncanddone_; // candidates that we investigated
	protected double              btncanddoneSucc_; // # investigated and succeeded
	protected double              btncanddoneFail_; // # investigated and failed
	
	protected BtBranchTracer       bter_;        // backtracer
	
	protected Checkpointer         cper_;        // structure for saving checkpoint cells
	protected double               cperMinlen_;  // minimum length for using checkpointer
	protected double               cperPerPow2_; // checkpoint every 1 << perpow2 diags (& next)
	protected boolean                 cperEf_;      // store E and F in addition to H?
	protected boolean                 cperTri_;     // checkpoint for triangular mini-fills?
	
	protected double              colstop_;      // bailed on DP loop after this many cols
	protected double              lastsolcol_;   // last DP col with valid cell
	protected double              cural_;        // index of next alignment to be given
	
	protected long nbtfiltst_; // # candidates filtered b/c starting cell was seen
	protected long nbtfiltsc_; // # candidates filtered b/c score uninteresting
	protected long nbtfiltdo_; // # candidates filtered b/c dominated by other cell
	
	protected OutputStream dpLog_;
	protected boolean firstRead_;
	protected double dpRows() {
		return rdf_ - rdi_;
	}
	protected long alignNucleotidesEnd2EndSseU8(  // unsigned 8-bit elements
			int flag, booleanean debug);
	protected long alignNucleotidesLocalSseU8(    // unsigned 8-bit elements
			int flag, boolean debug);
	protected long alignNucleotidesEnd2EndSseI16( // signed 16-bit elements
			int flag, boolean debug);
	protected long alignNucleotidesLocalSseI16(   // signed 16-bit elements
			int flag, boolean debug);
	protected long alignGatherEE8(                // unsigned 8-bit elements
			int flag, boolean debug);
	protected long alignGatherLoc8(               // unsigned 8-bit elements
			int flag, boolean debug);
	protected long alignGatherEE16(               // signed 16-bit elements
			int flag, boolean debug);
	protected long alignGatherLoc16(              // signed 16-bit elements
			int flag, boolean debug);
	protected void buildQueryProfileEnd2EndSseU8(boolean fw);
	protected void buildQueryProfileLocalSseU8(boolean fw);
	protected void buildQueryProfileEnd2EndSseI16(boolean fw);
	protected void buildQueryProfileLocalSseI16(boolean fw);
	
	protected boolean gatherCellsNucleotidesLocalSseU8(long best);
	protected boolean gatherCellsNucleotidesEnd2EndSseU8(long best);

	protected boolean gatherCellsNucleotidesLocalSseI16(long best);
	protected boolean gatherCellsNucleotidesEnd2EndSseI16(long best);
	
	protected boolean backtraceNucleotidesLocalSseU8(
			long       escore, // in: expected score
			SwResult      res,    // out: store results (edits and scores) here
			double        off,    // out: store diagonal projection of origin
			double        nbts,   // out: # backtracks
			double         row,    // start in this rectangle row
			double         col,    // start in this rectangle column
			RandomSource  rand);  // random gen, to choose among equal paths

	protected boolean backtraceNucleotidesLocalSseI16(
			long       escore, // in: expected score
			SwResult      res,    // out: store results (edits and scores) here
			double        off,    // out: store diagonal projection of origin
			double       nbts,   // out: # backtracks
			double         row,    // start in this rectangle row
			double         col,    // start in this rectangle column
			RandomSource  rand);  // random gen, to choose among equal paths

	protected boolean backtraceNucleotidesEnd2EndSseU8(
			long       escore, // in: expected score
			SwResult      res,    // out: store results (edits and scores) here
			double        off,    // out: store diagonal projection of origin
			double        nbts,   // out: # backtracks
			double         row,    // start in this rectangle row
			double         col,    // start in this rectangle column
			RandomSource  rand);  // random gen, to choose among equal paths

	protected boolean backtraceNucleotidesEnd2EndSseI16(
			long       escore, // in: expected score
			SwResult      res,    // out: store results (edits and scores) here
			double        off,    // out: store diagonal projection of origin
			double        nbts,   // out: # backtracks
			double         row,    // start in this rectangle row
			double         col,    // start in this rectangle column
			RandomSource  rand);  // random gen, to choose among equal paths
	
	protected boolean backtrace(
			long       escore, // in: expected score
			boolean           fill,   // in: use mini-fill?
			boolean           usecp,  // in: use checkpoints?
			SwResult      res,    // out: store results (edits and scores) here
			double        off,    // out: store diagonal projection of origin
			double         row,    // start in this rectangle row
			double         col,    // start in this rectangle column
			double         maxiter,// max # extensions to try
			double        niter,  // # extensions tried
			RandomSource  rnd)    // random gen, to choose among equal paths
		{
			bter_.initBt(
				escore,              // in: alignment score
				row,                 // in: start in this row
				col,                 // in: start in this column
				fill,                // in: use mini-fill?
				usecp,               // in: use checkpoints?
				cperTri_,            // in: triangle-shaped mini-fills?
				rnd);                // in: random gen, to choose among equal paths
			double nrej = 0;
			if(bter_.emptySolution()) {
				return false;
			} else {
				return bter_.nextAlignment(maxiter, res, off, nrej, niter, rnd);
			}
		}
	
	public static double ALPHA_SIZE = 5;
	public enum AlignerState {
		STATE_UNINIT,  // init() hasn't been called yet
		STATE_INITED,  // init() has been called, but not align()
		STATE_ALIGNED, // align() has been called
	}
	public SwAligner(OutputStream dpLog, booleanean firstRead) {
		sseU8fw_ =	sseU8rc_=sseI16fw_=sseI16rc_=btnstack_=btcells_=btncand_=rfwbuf_=btncanddone_=6;
		state_ = AlignerState.STATE_UNINIT;
		initedRead_=readSse16_=initedRef_=false;
		btncanddoneSucc_=btncanddoneFail_=colstop_=lastsolcol_=cural_=0;
		dpLog_ = dpLog;
		firstRead_ = firstRead;
	}
	
	public void printResultsStacked(SwResult res, OutputStream os) {
		res.alres.printStacked(rd_, os);
	}
	
	public booleanean done() {
		return cural_ == btncand_.size();
	}
	
	public booleanean initedRef() {
		return initedRef_;
	}
	
	public booleanean initedRead() {
		return initedRead_;
	}
	
	public void reset() {
		initedRef_ = initedRead_ = false;
	}
	
	public double numAlignmentsReported() { return cural_; }
	
	public void merge(
			SSEMetrics sseU8ExtendMet,
			SSEMetrics sseU8MateMet,
			SSEMetrics sseI16ExtendMet,
			SSEMetrics sseI16MateMet,
			long   nbtfiltst,
			long   nbtfiltsc,
			long   nbtfiltdo)
			{
		sseU8ExtendMet.merge(sseU8ExtendMet_);
		sseU8MateMet.merge(sseU8MateMet_);
		sseI16ExtendMet.merge(sseI16ExtendMet_);
		sseI16MateMet.merge(sseI16MateMet_);
		nbtfiltst += nbtfiltst_;
		nbtfiltsc += nbtfiltsc_;
		nbtfiltdo += nbtfiltdo_;
	}
	
	public void resetCounters() {
		sseU8ExtendMet_.reset();
		sseU8MateMet_.reset();
		sseI16ExtendMet_.reset();
		sseI16MateMet_.reset();
		nbtfiltst_ = nbtfiltsc_ = nbtfiltdo_ = 0;
	}
	
	public double size() {
		return dpRows() * (rff_ - rfi_);
	}
	
	public void initRead(
			BTDnaString rdfw, // forward read sequence
			BTDnaString rdrc, // revcomp read sequence
			BTString qufw,    // forward read qualities
			BTString qurc,    // reverse read qualities
			double rdi,              // offset of first read char to align
			double rdf,              // offset of last read char to align
			Scoring sc)       // scoring scheme
	{
		int nceil = sc.nCeil.f<Integer>((double)rdfw.length());
		rdfw_    = rdfw;      // read sequence
		rdrc_    = rdrc;      // read sequence
		qufw_    = qufw;      // read qualities
		qurc_    = qurc;      // read qualities
		rdi_     = rdi;        // offset of first read char to align
		rdf_     = rdf;        // offset of last read char to align
		sc_      = sc;        // scoring scheme
		nceil_   = nceil;      // max # Ns allowed in ref portion of aln
		readSse16_ = false;    // true . sse16 from now on for this read
		initedRead_ = true;

		if(dpLog_ != null) {
			if(!firstRead_) {
				dpLog_.write('\n');
			}
			dpLog_.write(rdfw.toZBuf() + '\t' + qufw.toZBuf());
		}
		firstRead_ = false;
	}
	
	public void initRef(
			Boolean fw,               // whether to forward or revcomp read is aligning
			long refidx,         // id of reference aligned against
			DPRect rect,    // DP rectangle
			char rf,              // reference sequence
			double rfi,            // offset of first reference char to align to
			double rff,            // offset of last reference char to align to
			long reflen,        // length of reference sequence
			Scoring sc,     // scoring scheme
			long minsc,        // minimum score
			Boolean enable8,          // use 8-bit SSE if possible?
			double cminlen,        // minimum length for using checkpointing scheme
			double cpow2,          // interval b/t checkpointed diags; 1 << this
			Boolean doTri,            // triangular mini-fills?
			Boolean extend)           // is this a seed extension?
			{
		double readGaps = sc.maxReadGaps(minsc, rdfw_.length());
		double refGaps  = sc.maxRefGaps(minsc, rdfw_.length());
		rdgap_       = readGaps;  // max # gaps in read
		rfgap_       = refGaps;   // max # gaps in reference
		state_       = AlignerState.STATE_INITED;
		fw_          = fw;       // orientation
		rd_          = fw ? rdfw_ : rdrc_; // read sequence
		qu_          = fw ? qufw_ : qurc_; // quality sequence
		refidx_      = refidx;   // id of reference aligned against
		rf_          = rf;       // reference sequence
		rfi_         = rfi;      // offset of first reference char to align to
		rff_         = rff;      // offset of last reference char to align to
		reflen_      = reflen;   // length of entire reference sequence
		rect_        = rect;    // DP rectangle
		minsc_       = minsc;    // minimum score
		cural_       = 0;        // idx of next alignment to give out
		initedRef_   = true;     // indicate we've initialized the ref portion
		enable8_     = enable8;  // use 8-bit SSE if possible?
		extend_      = extend;   // true iff this is a seed extension
		cperMinlen_  = cminlen;  // reads shorter than this won't use checkpointer
		cperPerPow2_ = cpow2;    // interval b/t checkpointed diags; 1 << this
		cperEf_      = true;     // whether to checkpoint H, E, and F
		cperTri_     = doTri;    // triangular mini-fills?
		bter_.initRef(
			fw_ ? rdfw_.buf() : // in: read sequence
				  rdrc_.buf(), 
			fw_ ? qufw_.buf() : // in: quality sequence
				  qurc_.buf(),
			rd_.length(),       // in: read sequence length
			rf_ + rfi_,          // in: reference sequence
			rff_ - rfi_,         // in: in-rectangle reference sequence length
			reflen,              // in: total reference sequence length
			refidx_,             // in: reference id
			rfi_,                // in: reference offset
			fw_,                 // in: orientation
			rect_,               // in: DP rectangle
			cper_,              // in: checkpointer
			sc_,                // in: scoring scheme
			nceil_);             // in: N ceiling
		// Record the reference sequence in the log
		if(dpLog_ != null) {
			dpLog_.write( '\t');
			dpLog_.write( refidx_ + ',');
			dpLog_.write( reflen_ + ',');
			dpLog_.write( minsc_ + ',');
			dpLog_.write( (fw ? '+' : '-') + ',');
			rect_.write(dpLog_);
			dpLog_.write( ',');
			for(long i = rfi_; i < rff_; i++) {
				dpLog_.write( mask2dna[(int)rf[i]];
			}
		}
	}
	
	public void initRef(
			boolean fw,               // whether to forward or revcomp read is aligning
			long refidx,         // reference aligned against
			DPRect rect,    // DP rectangle
			BitPairReference refs, // Reference strings
			long reflen,        // length of reference sequence
			Scoring sc,     // scoring scheme
			long minsc,        // minimum score
			boolean enable8,          // use 8-bit SSE if possible?
			double cminlen,        // minimum length for using checkpointing scheme
			double cpow2,          // interval b/t checkpointed diags; 1 << this
			boolean doTri,            // triangular mini-fills?
			boolean extend,           // true iff this is a seed extension
			double  upto,          // count the number of Ns up to this offset
			double nsUpto)        // output: the number of Ns up to 'upto'
			{
		long rfi = rect.refl;
		long rff = rect.refr + 1;
		assert_gt(rff, rfi);
		// Capture an extra reference character outside the rectangle so that we
		// can check matches in the next column over to the right
		rff++;
		// rflen = full length of the reference substring to consider, including
		// overhang off the boundaries of the reference sequence
		double rflen = (double)(rff - rfi);
		// Figure the number of Ns we're going to add to either side
		double leftNs  =
			(rfi >= 0               ? 0 : (double)Math.abs((long) rfi);
		leftNs = min(leftNs, rflen);
		double rightNs =
			(rff <= (long)reflen ? 0 : (double)Math.abs((long)(rff - reflen));
		rightNs = min(rightNs, rflen);
		// rflenInner = length of just the portion that doesn't overhang ref ends
		double rflenInner = rflen - (leftNs + rightNs);

		boolean haveRfbuf2 = false;
		EList<char> rfbuf2(rflen);
		// This is really slow, so only do it some of the time
		if((rand() % 10) == 0) {
			long rfii = rfi;
			for(double i = 0; i < rflen; i++) {
				if(rfii < 0 || (long)rfii >= reflen) {
					rfbuf2.push_back(4);
				} else {
					rfbuf2.push_back(refs.getBase(refidx, (double)rfii));
				}
				rfii++;
			}
			haveRfbuf2 = true;
		}
		// rfbuf_ = int list large enough to accommodate both the reference
		// sequence and any Ns we might add to either side.
		rfwbuf_.resize((rflen + 16) / 4);
		int offset = refs.getStretch(
			rfwbuf_.ptr(),               // buffer to store words in
			refidx,                      // which reference
			(rfi < 0) ? 0 : (double)rfi, // starting offset (can't be < 0)
			rflenInner);// for BitPairReference::getStretch()
		rf_ = (String)rfwbuf_.ptr() + offset;
		// Shift ref chars away from 0 so we can stick Ns at the beginning
		if(leftNs > 0) {
			// Slide everyone down
			for(double i = rflenInner; i > 0; i--) {
				rf_[i+leftNs-1] = rf_[i-1];
			}
			// Add Ns
			for(double i = 0; i < leftNs; i++) {
				rf_[i] = 4;
			}
		}
		if(rightNs > 0) {
			// Add Ns to the end
			for(double i = 0; i < rightNs; i++) {
				rf_[i + leftNs + rflenInner] = 4;
			}
		}

		// Count Ns and convert reference characters into A/C/G/T masks.  Ambiguous
		// nucleotides (IUPAC codes) have more than one mask bit set.  If a
		// reference scanner was provided, use it to opportunistically resolve seed
		// hits.
		nsUpto = 0;
		for(double i = 0; i < rflen; i++) {
			// rf_[i] gets mask version of refence char, with N=16
			if(i < upto && rf_[i] > 3) {
				nsUpto++;
			}
			rf_[i] = (1 << rf_[i]);
		}
		// Correct for having captured an extra reference character
		rff--;
		initRef(
			fw,          // whether to forward or revcomp read is aligning
			refidx,      // id of reference aligned against
			rect,        // DP rectangle
			rf_,         // reference sequence, wrapped up in BTString object
			0,           // use the whole thing
			(double)(rff - rfi), // ditto
			reflen,      // reference length
			sc,          // scoring scheme
			minsc,       // minimum score
			enable8,     // use 8-bit SSE if possible?
			cminlen,     // minimum length for using checkpointing scheme
			cpow2,       // interval b/t checkpointed diags; 1 << this
			doTri,       // triangular mini-fills?
			extend);     // true iff this is a seed extension
	}
	
	public ungappedAlign(
			BTDnaString      rd,     // read sequence (could be RC)
			BTString         qu,     // qual sequence (could be rev)
			Coord            coord,  // coordinate aligned to
			BitPairReference refs,   // Reference strings
			double                  reflen, // length of reference sequence
			Scoring          sc,     // scoring scheme
			boolean                    ohang,  // allow overhang?
			long                minsc,  // minimum score
			SwResult               res)    // put alignment result here
	{
		double len = rd.length();
		int nceil = sc.nCeil.f<Integer>((double)len);
		int ns = 0;
		long rfi = coord.off();
		long rff = rfi + (long)len;
		long refidx = coord.ref();
		// Figure the number of Ns we're going to add to either side
		double leftNs = 0;
		if(rfi < 0) {
			if(ohang) {
				leftNs = (double)(-rfi);
			} else {
				return 0;
			}
		}
		double rightNs = 0;
		if(rff > (long)reflen) {
			if(ohang) {
				rightNs = (double)(rff - (long)reflen);
			} else {
				return 0;
			}
		}
		if((leftNs + rightNs) > (double)nceil) {
			return 0;
		}
		// rflenInner = length of just the portion that doesn't overhang ref ends
		double rflenInner = len - (leftNs + rightNs);

		boolean haveRfbuf2 = false;
		EList<char> rfbuf2(len);
		// This is really slow, so only do it some of the time
		if((rand() % 10) == 0) {
			long rfii = rfi;
			for(double i = 0; i < len; i++) {
				if(rfii < 0 || (double)rfii >= reflen) {
					rfbuf2.push_back(4);
				} else {
					rfbuf2.push_back(refs.getBase(refidx, (double)rfii));
				}
				rfii++;
			}
			haveRfbuf2 = true;
		}

		// rfbuf_ = int list large enough to accommodate both the reference
		// sequence and any Ns we might add to either side.
		rfwbuf_.resize((len + 16) / 4);
		int offset = refs.getStretch(
			rfwbuf_.ptr(),               // buffer to store words in
			refidx,                      // which reference
			(rfi < 0) ? 0 : (double)rfi, // starting offset (can't be < 0)
			rflenInner);// for BitPairReference::getStretch()

		rf_ = (String)rfwbuf_.ptr() + offset;
		// Shift ref chars away from 0 so we can stick Ns at the beginning
		if(leftNs > 0) {
			// Slide everyone down
			for(double i = rflenInner; i > 0; i--) {
				rf_[i+leftNs-1] = rf_[i-1];
			}
			// Add Ns
			for(double i = 0; i < leftNs; i++) {
				rf_[i] = 4;
			}
		}
		if(rightNs > 0) {
			// Add Ns to the end
			for(double i = 0; i < rightNs; i++) {
				rf_[i + leftNs + rflenInner] = 4;
			}
		}

		// Count Ns and convert reference characters into A/C/G/T masks.  Ambiguous
		// nucleotides (IUPAC codes) have more than one mask bit set.  If a
		// reference scanner was provided, use it to opportunistically resolve seed
		// hits.
		long score = 0;
		res.alres.reset();
		double rowi = 0;
		double rowf = len-1;
		if(sc.monotone) {
			for(double i = 0; i < len; i++) {
				// rf_[i] gets mask version of refence char, with N=16
				score += sc.score(rd[i], (int)(1 << rf_[i]), qu[i] - 33, ns);
				if(score < minsc || ns > nceil) {
					// Fell below threshold
					return 0;
				}
			}
			// Got a result!  Fill in the rest of the result object.
		} else {
			// Definitely ways to short-circuit this.  E.g. if diff between cur
			// score and minsc can't be met by matches.
			long floorsc = 0;
			long scoreMax = floorsc;
			double lastfloor = 0;
			rowi = MAX_SIZE_T;
			double sols = 0;
			for(double i = 0; i < len; i++) {
				score += sc.score(rd[i], (int)(1 << rf_[i]), qu[i] - 33, ns);
				if(score >= minsc && score >= scoreMax) {
					scoreMax = score;
					rowf = i;
					if(rowi != lastfloor) {
						rowi = lastfloor;
						sols++;
					}
				}
				if(score <= floorsc) {
					score = floorsc;
					lastfloor = i+1;
				}
			}
			if(ns > nceil || scoreMax < minsc) {
				// Too many Ns
				return 0;
			}
			if(sols > 1) {
				// >1 distinct solution in this diag; defer to DP aligner
				return -1;
			}
			score = scoreMax;
			// Got a result!  Fill in the rest of the result object.  
		}
		// Now fill in the edits

		EList<Edit> ned = res.alres.ned();
		double refns = 0;
		for(double i = rowi; i <= rowf; i++) {
			if(rf_[i] > 3 || rd[i] != rf_[i]) {
				// Add edit
				Edit e = new Edit((int)i,
				       mask2dna[1 << (int)rf_[i]],
				       "ACGTN"[(int)rd[i]],
				       EDIT_TYPE_MM);
				ned.push_back(e);
				if(rf_[i] > 3) {
					refns++;
				}
			}
		}
		res.alres.setScore(new AlnScore(score,
									(int)(rd.length() - ned.size()),
									(int)ned.size(), ns, 0));
		boolean fw = coord.fw();
		double trimEnd = (len-1) - rowf;
		res.alres.setShape(
			coord.ref(),  // ref id
			coord.off()+rowi, // 0-based ref offset
			reflen,       // length of reference sequence aligned to
			fw,           // aligned to Watson?
			len,          // read length
			true,         // pretrim soft?
			0,            // pretrim 5' end
			0,            // pretrim 3' end
			true,         // alignment trim soft?
			fw ? rowi : trimEnd,  // alignment trim 5' end
			fw ? trimEnd : rowi); // alignment trim 3' end
		res.alres.setRefNs(refns);

		BTDnaString editstr;
		Edit.toRef(rd, ned, editstr, true, rowi, trimEnd);
		if(refstr != editstr) {
			System.err.println( "Decoded nucleotides and edits don't match reference:" );
			System.err.println( "           score: " + res.alres.score().score());
			System.err.println( "           edits: ");
			Edit.print(cerr, ned);
			System.err.println();
			System.err.println( "    decoded nucs: " + rd + endl;
			System.err.println( "     edited nucs: " + editstr + endl;
			System.err.println( "  reference nucs: " + refstr + endl;

		}
		if(!fw) {
			// All edits are currently w/r/t upstream end; if read aligned to Crick
			// strand, invert them to be w/r/t 5' end instead.
			res.alres.invertEdits();
		}
		return 1;
	}
	
	public boolean align(long best) {
		state_ = AlignerState.STATE_ALIGNED;
		// Reset solutions lists
		btncand_.clear();
		btncanddone_.clear();
		btncanddoneSucc_ = btncanddoneFail_ = 0;
		best = Long.MIN_VALUE;
		sse8succ_ = sse16succ_ = false;
		int flag = 0;
		double rdlen = rdf_ - rdi_;
		boolean checkpointed = rdlen >= cperMinlen_;
		boolean gathered = false; // Did gathering happen along with alignment?
		if(sc_.monotone) {
			// End-to-end
			if(enable8_ && !readSse16_ && minsc_ >= -254) {
				// 8-bit end-to-end
				if(checkpointed) {
					best = alignGatherEE8(flag, false);
					if(flag == 0) {
						gathered = true;
					}
				} else {
					best = alignNucleotidesEnd2EndSseU8(flag, false);
					int flagtmp = 0;
					long besttmp = alignGatherEE8(flagtmp, true); // debug
				}
				sse8succ_ = (flag == 0);

				{
					int flag2 = 0;
					long best2 = alignNucleotidesEnd2EndSseI16(flag2, true);
					{
						int flagtmp = 0;
						long besttmp = alignGatherEE16(flagtmp, true);
					}
					sse16succ_ = (flag2 == 0);
				}

			} else {
				// 16-bit end-to-end
				if(checkpointed) {
					best = alignGatherEE16(flag, false);
					if(flag == 0) {
						gathered = true;
					}
				} else {
					best = alignNucleotidesEnd2EndSseI16(flag, false);
					int flagtmp = 0;
					long besttmp = alignGatherEE16(flagtmp, true);
				}
				sse16succ_ = (flag == 0);
			}
		} else {
			// Local
			flag = -2;
			if(enable8_ && !readSse16_) {
				// 8-bit local
				if(checkpointed) {
					best = alignGatherLoc8(flag, false);
					if(flag == 0) {
						gathered = true;
					}
				} else {
					best = alignNucleotidesLocalSseU8(flag, false);
					int flagtmp = 0;
					long besttmp = alignGatherLoc8(flagtmp, true);
				}
			}
			if(flag == -2) {
				// 16-bit local
				flag = 0;
				if(checkpointed) {
					best = alignNucleotidesLocalSseI16(flag, false);
					best = alignGatherLoc16(flag, false);
					if(flag == 0) {
						gathered = true;
					}
				} else {
					best = alignNucleotidesLocalSseI16(flag, false);
					int flagtmp = 0;
					long besttmp = alignGatherLoc16(flagtmp, true);
				}
				sse16succ_ = (flag == 0);
			} else {
				sse8succ_ = (flag == 0);
				int flag2 = 0;
				long best2 = alignNucleotidesLocalSseI16(flag2, true);
				{
					int flagtmp = 0;
					long besttmp = alignGatherLoc16(flagtmp, true);
				}
				sse16succ_ = (flag2 == 0);
			}
		}
		if(!checkpointed && ((Math.random() * Integer.MAX_VALUE) & 15) == 0 && sse8succ_ && sse16succ_) {
			SSEData d8  = fw_ ? sseU8fw_  : sseU8rc_;
			SSEData d16 = fw_ ? sseI16fw_ : sseI16rc_;
			for(double i = 0; i < d8.mat_.nrow(); i++) {
				for(double j = 0; j < colstop_; j++) {
					int h8  = d8.mat_.helt(i, j);
					int h16 = d16.mat_.helt(i, j);
					int e8  = d8.mat_.eelt(i, j);
					int e16 = d16.mat_.eelt(i, j);
					int f8  = d8.mat_.felt(i, j);
					int f16 = d16.mat_.felt(i, j);
					long h8s  =
						(sc_.monotone ? (h8  - 0xff  ) : h8);
					long h16s =
						(sc_.monotone ? (h16 - 0x7fff) : (h16 + 0x8000));
					long e8s  =
						(sc_.monotone ? (e8  - 0xff  ) : e8);
					long e16s =
						(sc_.monotone ? (e16 - 0x7fff) : (e16 + 0x8000));
					long f8s  =
						(sc_.monotone ? (f8  - 0xff  ) : f8);
					long f16s =
						(sc_.monotone ? (f16 - 0x7fff) : (f16 + 0x8000));
					if(h8s < minsc_) {
						h8s = minsc_ - 1;
					}
					if(h16s < minsc_) {
						h16s = minsc_ - 1;
					}
					if(e8s < minsc_) {
						e8s = minsc_ - 1;
					}
					if(e16s < minsc_) {
						e16s = minsc_ - 1;
					}
					if(f8s < minsc_) {
						f8s = minsc_ - 1;
					}
					if(f16s < minsc_) {
						f16s = minsc_ - 1;
					}
				}
			}
		}
		cural_ = 0;
		if(best == Long.MIN_VALUE || best < minsc_) {
			if(dpLog_ != null) {
				dpLog_.write( ",0,0");
			}
			return false;
		}
		if(!gathered) {
			// Look for solutions using SSE matrix
			if(sc_.monotone) {
				if(sse8succ_) {
					gatherCellsNucleotidesEnd2EndSseU8(best);

					if(sse16succ_) {
						cand_tmp_ = btncand_;
						gatherCellsNucleotidesEnd2EndSseI16(best);
						cand_tmp_.sort();
						btncand_.sort();
					}
				} else {
					gatherCellsNucleotidesEnd2EndSseI16(best);
				}
			} else {
				if(sse8succ_) {
					gatherCellsNucleotidesLocalSseU8(best);
					if(sse16succ_) {
						cand_tmp_ = btncand_;
						gatherCellsNucleotidesLocalSseI16(best);
						cand_tmp_.sort();
						btncand_.sort();
					}
				} else {
					gatherCellsNucleotidesLocalSseI16(best);
				}
			}
		}
		if(!btncand_.empty()) {
			btncand_.sort();
		}
		if(dpLog_ != null) {
			dpLog_.write( ",1," + best);
		}
		return !btncand_.empty();
	}
	
	public boolean nextAlignment(
			SwResult res,
			long minsc,
			RandomSource rnd){
		if(done()) {
			res.reset();
			return false;
		}
		double off = 0, nbts = 0;
		// For each candidate cell that we should try to backtrack from...
		double candsz = btncand_.size();
		double SQ = dpRows() >> 4;
		if(SQ == 0) SQ = 1;
		double rdlen = rdf_ - rdi_;
		boolean checkpointed = rdlen >= cperMinlen_;
		while(cural_ < candsz) {
			// Doing 'continue' anywhere in here simply causes us to move on to the
			// next candidate
			if(btncand_[cural_].score < minsc) {
				btncand_[cural_].fate = BT_CAND_FATE_FILT_SCORE;
				nbtfiltsc_++; cural_++; continue;
			}
			nbts = 0;
			double row = btncand_[cural_].row;
			double col = btncand_[cural_].col;
			if(sse16succ_) {
				SSEData d = fw_ ? sseI16fw_ : sseI16rc_;
				if(!checkpointed && d.mat_.reset_[row] && d.mat_.reportedThrough(row, col)) {
					// Skipping this candidate because a previous candidate already
					// moved through this cell
					btncand_[cural_].fate = BT_CAND_FATE_FILT_START;
					//System.err.println( "  skipped becuase starting cell was covered" << endl;
					nbtfiltst_++; cural_++; continue;
				}
			} else if(sse8succ_) {
				SSEData d = fw_ ? sseU8fw_ : sseU8rc_;
				if(!checkpointed && d.mat_.reset_[row] && d.mat_.reportedThrough(row, col)) {
					// Skipping this candidate because a previous candidate already
					// moved through this cell
					btncand_[cural_].fate = BT_CAND_FATE_FILT_START;
					//System.err.println( "  skipped becuase starting cell was covered" << endl;
					nbtfiltst_++; cural_++; continue;
				}
			}
			if(sc_.monotone) {
				boolean ret = false;
				if(sse8succ_) {
					int reseed = rnd.nextU32() + 1;
					rnd.init(reseed);
					res.reset();
					if(checkpointed) {
						double maxiter = MAX_SIZE_T;
						double niter = 0;
						ret = backtrace(
							btncand_[cural_].score, // in: expected score
							true,     // in: use mini-fill?
							true,     // in: use checkpoints?
							res,      // out: store results (edits and scores) here
							off,      // out: store diagonal projection of origin
							row,      // start in this rectangle row
							col,      // start in this rectangle column
							maxiter,  // max # extensions to try
							niter,    // # extensions tried
							rnd);     // random gen, to choose among equal paths
					} else {
						ret = backtraceNucleotidesEnd2EndSseU8(
							btncand_[cural_].score, // in: expected score
							res,    // out: store results (edits and scores) here
							off,    // out: store diagonal projection of origin
							nbts,   // out: # backtracks
							row,    // start in this rectangle row
							col,    // start in this rectangle column
							rnd);   // random gen, to choose among equal paths
					}
					// if(...) statement here should check not whether the primary
					// alignment was checkpointed, but whether a checkpointed
					// alignment was done at all.
					if(!checkpointed) {
						SwResult res2;
						double maxiter2 = MAX_SIZE_T;
						double niter2 = 0;
						boolean ret2 = backtrace(
							btncand_[cural_].score, // in: expected score
							true,     // in: use mini-fill?
							true,     // in: use checkpoints?
							res2,     // out: store results (edits and scores) here
							off,      // out: store diagonal projection of origin
							row,      // start in this rectangle row
							col,      // start in this rectangle column
							maxiter2, // max # extensions to try
							niter2,   // # extensions tried
							rnd);     // random gen, to choose among equal paths
						// After the first alignment, there's no guarantee we'll
						// get the same answer from both backtrackers because of
						// differences in how they handle marking cells as
						// reported-through.
					}
					if(sse16succ_ && !checkpointed) {
						SwResult res2;
						double off2, nbts2 = 0;
						rnd.init(reseed);
						boolean ret2 = backtraceNucleotidesEnd2EndSseI16(
							btncand_[cural_].score, // in: expected score
							res2,   // out: store results (edits and scores) here
							off2,   // out: store diagonal projection of origin
							nbts2,  // out: # backtracks
							row,    // start in this rectangle row
							col,    // start in this rectangle column
							rnd);   // random gen, to choose among equal paths

					}
					rnd.init(reseed+1); // debug/release pseudo-randoms in lock step
				} else if(sse16succ_) {
					int reseed = rnd.nextU32() + 1;
					res.reset();
					if(checkpointed) {
						double maxiter = MAX_SIZE_T;
						double niter = 0;
						ret = backtrace(
							btncand_[cural_].score, // in: expected score
							true,     // in: use mini-fill?
							true,     // in: use checkpoints?
							res,      // out: store results (edits and scores) here
							off,      // out: store diagonal projection of origin
							row,      // start in this rectangle row
							col,      // start in this rectangle column
							maxiter,  // max # extensions to try
							niter,    // # extensions tried
							rnd);     // random gen, to choose among equal paths
					} else {
						ret = backtraceNucleotidesEnd2EndSseI16(
							btncand_[cural_].score, // in: expected score
							res,    // out: store results (edits and scores) here
							off,    // out: store diagonal projection of origin
							nbts,   // out: # backtracks
							row,    // start in this rectangle row
							col,    // start in this rectangle column
							rnd);   // random gen, to choose among equal paths
					}

					// if(...) statement here should check not whether the primary
					// alignment was checkpointed, but whether a checkpointed
					// alignment was done at all.
					if(!checkpointed) {
						SwResult res2;
						double maxiter2 = MAX_SIZE_T;
						double niter2 = 0;
						boolean ret2 = backtrace(
							btncand_[cural_].score, // in: expected score
							true,     // in: use mini-fill?
							true,     // in: use checkpoints?
							res2,     // out: store results (edits and scores) here
							off,      // out: store diagonal projection of origin
							row,      // start in this rectangle row
							col,      // start in this rectangle column
							maxiter2, // max # extensions to try
							niter2,   // # extensions tried
							rnd);     // random gen, to choose among equal paths
						// After the first alignment, there's no guarantee we'll
						// get the same answer from both backtrackers because of
						// differences in how they handle marking cells as
						// reported-through.
					}
					rnd.init(reseed); // debug/release pseudo-randoms in lock step
				}
				if(ret) {
					btncand_[cural_].fate = BT_CAND_FATE_SUCCEEDED;
					break;
				} else {
					btncand_[cural_].fate = BT_CAND_FATE_FAILED;
				}
			} else {
				// Local alignment
				// Check if this solution is "dominated" by a prior one.
				// Domination is a heuristic designed to eliminate the vast
				// majority of valid-but-redundant candidates lying in the
				// "penumbra" of a high-scoring alignment.
				boolean dom = false;
				{
					double donesz = btncanddone_.size();
					double col = btncand_[cural_].col;
					double row = btncand_[cural_].row;
					for(double i = 0; i < donesz; i++) {
						double colhi = col, rowhi = row;
						double rowlo = btncanddone_[i].row;
						double collo = btncanddone_[i].col;
						if(colhi < collo) swap(colhi, collo);
						if(rowhi < rowlo) swap(rowhi, rowlo);
						if(colhi - collo <= SQ && rowhi - rowlo <= SQ) {
							// Skipping this candidate because it's "dominated" by
							// a previous candidate
							dom = true;
							break;
						}
					}
				}
				if(dom) {
					btncand_[cural_].fate = BT_CAND_FATE_FILT_DOMINATED;
					nbtfiltdo_++;
					cural_++;
					continue;
				}
				boolean ret = false;
				if(sse8succ_) {
					int reseed = rnd.nextU32() + 1;
					res.reset();
					rnd.init(reseed);
					if(checkpointed) {
						double maxiter = MAX_SIZE_T;
						double niter = 0;
						ret = backtrace(
							btncand_[cural_].score, // in: expected score
							true,     // in: use mini-fill?
							true,     // in: use checkpoints?
							res,      // out: store results (edits and scores) here
							off,      // out: store diagonal projection of origin
							row,      // start in this rectangle row
							col,      // start in this rectangle column
							maxiter,  // max # extensions to try
							niter,    // # extensions tried
							rnd);     // random gen, to choose among equal paths
					} else {
						ret = backtraceNucleotidesLocalSseU8(
							btncand_[cural_].score, // in: expected score
							res,    // out: store results (edits and scores) here
							off,    // out: store diagonal projection of origin
							nbts,   // out: # backtracks
							row,    // start in this rectangle row
							col,    // start in this rectangle column
							rnd);   // random gen, to choose among equal paths
					}
					// if(...) statement here should check not whether the primary
					// alignment was checkpointed, but whether a checkpointed
					// alignment was done at all.
					if(!checkpointed) {
						SwResult res2;
						double maxiter2 = MAX_SIZE_T;
						double niter2 = 0;
						boolean ret2 = backtrace(
							btncand_[cural_].score, // in: expected score
							true,     // in: use mini-fill?
							true,     // in: use checkpoints?
							res2,     // out: store results (edits and scores) here
							off,      // out: store diagonal projection of origin
							row,      // start in this rectangle row
							col,      // start in this rectangle column
							maxiter2, // max # extensions to try
							niter2,   // # extensions tried
							rnd);     // random gen, to choose among equal paths
						// After the first alignment, there's no guarantee we'll
						// get the same answer from both backtrackers because of
						// differences in how they handle marking cells as
						// reported-through.
						// TODO: I find that sometimes there is disagreement here
						// where the alignments are in the same place with
						// identical scores, but one is more soft-trimmed than the other
						//assert(cural_ > 0 || !ret || res.alres == res2.alres);
					}
					if(!checkpointed && sse16succ_) {
						SwResult res2;
						double off2, nbts2 = 0;
						rnd.init(reseed); // same b/t backtrace calls
						boolean ret2 = backtraceNucleotidesLocalSseI16(
							btncand_[cural_].score, // in: expected score
							res2,   // out: store results (edits and scores) here
							off2,   // out: store diagonal projection of origin
							nbts2,  // out: # backtracks
							row,    // start in this rectangle row
							col,    // start in this rectangle column
							rnd);   // random gen, to choose among equal paths
					}
					rnd.init(reseed+1); // debug/release pseudo-randoms in lock step
				} else if(sse16succ_) {
					int reseed = rnd.nextU32() + 1;
					res.reset();
					if(checkpointed) {
						double maxiter = MAX_SIZE_T;
						double niter = 0;
						ret = backtrace(
							btncand_[cural_].score, // in: expected score
							true,     // in: use mini-fill?
							true,     // in: use checkpoints?
							res,      // out: store results (edits and scores) here
							off,      // out: store diagonal projection of origin
							row,      // start in this rectangle row
							col,      // start in this rectangle column
							maxiter,  // max # extensions to try
							niter,    // # extensions tried
							rnd);     // random gen, to choose among equal paths
					} else {
						ret = backtraceNucleotidesLocalSseI16(
							btncand_[cural_].score, // in: expected score
							res,    // out: store results (edits and scores) here
							off,    // out: store diagonal projection of origin
							nbts,   // out: # backtracks
							row,    // start in this rectangle row
							col,    // start in this rectangle column
							rnd);   // random gen, to choose among equal paths
					}

					// if(...) statement here should check not whether the primary
					// alignment was checkpointed, but whether a checkpointed
					// alignment was done at all.
					if(!checkpointed) {
						SwResult res2;
						double maxiter2 = MAX_SIZE_T;
						double niter2 = 0;
						boolean ret2 = backtrace(
							btncand_[cural_].score, // in: expected score
							true,     // in: use mini-fill?
							true,     // in: use checkpoints?
							res2,     // out: store results (edits and scores) here
							off,      // out: store diagonal projection of origin
							row,      // start in this rectangle row
							col,      // start in this rectangle column
							maxiter2, // max # extensions to try
							niter2,   // # extensions tried
							rnd);     // random gen, to choose among equal paths
						// After the first alignment, there's no guarantee we'll
						// get the same answer from both backtrackers because of
						// differences in how they handle marking cells as
						// reported-through.
					}
					rnd.init(reseed); // same b/t backtrace calls
				}
				if(ret) {
					btncand_[cural_].fate = BT_CAND_FATE_SUCCEEDED;
					btncanddone_.push_back(btncand_[cural_]);
					btncanddoneSucc_++;
					break;
				} else {
					btncand_[cural_].fate = BT_CAND_FATE_FAILED;
					btncanddone_.push_back(btncand_[cural_]);
					btncanddoneFail_++;
				}
			}
			cural_++;
		} // while(cural_ < btncand_.size())
		if(cural_ == btncand_.size()) {
			return false;
		}
		if(!fw_) {
			// All edits are currently w/r/t upstream end; if read aligned
			// to Crick strand, we need to invert them so that they're
			// w/r/t the read's 5' end instead.
			res.alres.invertEdits();
		}
		cural_++;
		return true;
	}
	
	
}
