package com.uwb.bt2j.indexer;

import java.io.*;
import java.util.Scanner;

import com.uwb.bt2j.indexer.blockwise.KarkkainenBlockwiseSA;
import com.uwb.bt2j.indexer.filebuf.FileBuf;
import com.uwb.bt2j.indexer.types.EList;
import com.uwb.bt2j.indexer.types.RandomSource;
import com.uwb.bt2j.indexer.types.RefReadInParams;
import com.uwb.bt2j.indexer.types.RefRecord;
import com.uwb.bt2j.indexer.types.SideLocus;
import com.uwb.bt2j.indexer.util.EbwtParams;
import com.uwb.bt2j.indexer.util.IndexTypes;

public abstract class EBWT <TStr>{
	public static final String gEbwt_ext = "bt2";
	public String gLastIOErrMsg;
	public int    _overrideOffRate;
	public boolean       _verbose;
	public boolean       _passMemExc;
	public boolean       _sanity;
	public boolean       fw_;     // true iff this is a forward index
	public File       _in1;    // input fd for primary index file
	public File       _in2;    // input fd for secondary index file
	public String     _in1Str; // filename for primary index file
	public String     _in2Str; // filename for secondary index file
	public String     _inSaStr;  // filename for suffix-array file
	public String     _inBwtStr; // filename for BWT file
	public long  _zOff;
	public long  _zEbwtByteOff;
	public long   _zEbwtBpOff;
	public long  _nPat;  /// number of reference texts
	public int  _nFrag; /// number of fragments
	long[] _plen;
	int[] _rstarts;
	long _fchr;
	long _ftab;
	long _eftab;
	long[] _offs;
	byte _ebwt;
	public boolean       _useMm;        /// use memory-mapped files to hold the index
	public boolean       useShmem_;     /// use shared memory to hold large parts of the index
	public EList<String> _refnames; /// names of the reference sequences
	public String mmFile1_;
	public String mmFile2_;
	public EbwtParams _eh;
	public boolean packed_;
	public static final long default_bmax = IndexTypes.OFF_MASK;
	public static final long default_bmaxMultSqrt = IndexTypes.OFF_MASK;
	public static final long default_bmaxDivN = 4;
	public static final int      default_dcv = 1024;
	public static final boolean     default_noDc = false;
	public static final boolean     default_useBlockwise = true;
	public static final int default_seed = 0;
	public static final int      default_lineRate = 6;
	public static final int      default_offRate = 5;
	public static final int      default_offRatePlus = 0;
	public static final int      default_ftabChars = 10;
	public static final boolean     default_bigEndian = false;
	
	public enum EbwtFlags {
		EBWT_COLOR(2),
		EBWT_ENTIRE_REV(4);
		private int x;
		EbwtFlags(int y){x = y;}
	}
	
	public EBWT(String in,
				int color,
				int needEntireReverse,
				boolean fw,
				int overrideOffRate, // = -1,
				int offRatePlus, // = -1,
				boolean useMm, // = false,
				boolean useShmem, // = false,
				boolean mmSweep, // = false,
				boolean loadNames, // = false,
				boolean loadSASamp, // = true,
				boolean loadFtab, // = true,
				boolean loadRstarts, // = true,
				boolean verbose, // = false,
				boolean startVerbose, // = false,
				boolean passMemExc, // = false,
				boolean sanityCheck) {
		_overrideOffRate = overrideOffRate;
		_verbose = verbose;
		_passMemExc = passMemExc;
		_sanity = sanityCheck;
		fw_ = fw;
		_in1 = null;
		_in2 = null;
		_zOff = IndexTypes.OFF_MASK;
		_zEbwtByteOff = IndexTypes.OFF_MASK;
		_zEbwtBpOff = -1;
		_nPat = 0;
		_nFrag = 0;
		//_plen = 1;
		//_rstarts = 1;
		_fchr = 1;
		_ftab = 1;
		_eftab = 1;
		//_offs = 1;
		_ebwt = 1;
		_useMm = false;
		useShmem_ = false;
		_refnames = new EList(1);
		mmFile1_ = null;
		mmFile2_ = null;
		packed_ = false;
		_useMm = useMm;
		useShmem_ = useShmem;
		_in1Str = in + ".1." + gEbwt_ext;
		_in2Str = in + ".2." + gEbwt_ext;
		readIntoMemory(
				color,       // expect index to be colorspace?
				fw ? -1 : needEntireReverse, // need REF_READ_REVERSE
				loadSASamp,  // load the SA sample portion?
				loadFtab,    // load the ftab & eftab?
				loadRstarts, // load the rstarts array?
				true,        // stop after loading the header portion?
				_eh,        // params
				mmSweep,     // mmSweep
				loadNames,   // loadNames
				startVerbose); // startVerbose
		// If the offRate has been overridden, reflect that in the
		// _eh._offRate field
		if(offRatePlus > 0 && _overrideOffRate == -1) {
			_overrideOffRate = _eh._offRate + offRatePlus;
		}
		if(_overrideOffRate > _eh._offRate) {
			_eh.setOffRate(_overrideOffRate);
		}
	}
	
	public EBWT(TStr exampleStr,
				boolean packed,
				int color,
				int needEntireReverse,
				int lineRate,
				int offRate,
				int ftabChars,
				int nthreads,
				String file,   // base filename for EBWT files
				boolean fw,
				boolean useBlockwise,
				long bmax,
				long bmaxSqrtMult,
				long bmaxDivN,
				int dcv,
				EList<FileBuf> is,
				EList<RefRecord> szs,
				long sztot,
				RefReadInParams refparams,
				int seed,
				int overrideOffRate,
				boolean doSaFile,
				boolean doBwtFile,
				boolean verbose,
				boolean passMemExc,
				boolean sanityCheck) throws FileNotFoundException, IOException
	{
		_eh(
				joinedLen(szs),
				lineRate,
				offRate,
				ftabChars,
				color,
				refparams.reverse == ReadDir.REF_READ_REVERSE
		);
		_overrideOffRate = overrideOffRate;
		_verbose = verbose;
		_passMemExc = passMemExc;
		_sanity = sanityCheck;
		fw_ = fw;
		_in1 = null;
		_in2 = null;
		_zOff = IndexTypes.OFF_MASK;
		_zEbwtByteOff = IndexTypes.OFF_MASK;
		_zEbwtBpOff = -1;
		_nPat = 0;
		_nFrag = 0;
		//_plen = 1;
		//_rstarts = 1;
		_fchr = 1;
		_ftab = 1;
		_eftab = 1;
		//_offs = 1;
		_ebwt = 1;
		_useMm = false;
		useShmem_ = false;
		//_refnames = 1;
		mmFile1_ = null;
		mmFile2_ = null;
		_in1Str = file + ".1." + gEbwt_ext;
		_in2Str = file + ".2." + gEbwt_ext;
		packed_ = packed;
		// Open output files
		FileOutputStream fout1 = new FileOutputStream(_in1Str);
		if(!fout1.getChannel().isOpen()) {
			System.err.println("Could not open index file for writing: \"" + _in1Str + "\n" + "Please make sure the directory exists and that permissions allow writing by Bowtie.");
		}
		FileOutputStream fout2 = new FileOutputStream(_in2Str);
		if(!fout2.getChannel().isOpen()) {
			System.err.println("Could not open index file for writing: \"" + _in2Str + "\"" + "\n" +
					"Please make sure the directory exists and that permissions allow writing by Bowtie." + "\n");
		}
		_inSaStr = file + ".sa";
		_inBwtStr = file + ".bwt";
		FileInputStream saOut = null, bwtOut = null;
		if(doSaFile) {
			saOut = new FileInputStream(_inSaStr);
			if(saOut.available() <= 0) {
				System.err.println("Could not open suffix-array file for writing: \"" + _inSaStr + "Please make sure the directory exists and that permissions allow writing by Bowtie.");
			}
		}
		if(doBwtFile) {
			bwtOut = new FileInputStream(_inBwtStr);
			if(bwtOut.available() <= 0) {
				System.err.println("Could not open suffix-array file for writing: \"" + _inBwtStr + "\"" + "\n" +
						"Please make sure the directory exists and that permissions allow writing by" + "\n"
						+ "Bowtie.");
			}
		}
		// Build SA(T) and BWT(T) block by block
		initFromVector<TStr>(
			is,
			szs,
			sztot,
			refparams,
			fout1,
			fout2,
			file,
			saOut,
			bwtOut,
			nthreads,
			useBlockwise,
			bmax,
			bmaxSqrtMult,
			bmaxDivN,
			dcv,
			seed,
			verbose);
		// Close output files
		fout1.flush();
		/**
		 long tellpSz1 = (long)fout1.tellp();
		 VMSG_NL("Wrote " + fout1.tellp() + " bytes to primary EBWT file: " + _in1Str);
		 fout1.close();
		 boolean err = false;
		 if(tellpSz1 > fileSize(_in1Str)) {
		 err = true;
		 System.err.println("Index is corrupt: File size for " + _in1Str + " should have been " + tellpSz1
		 + " but is actually " + fileSize(_in1Str) + "." + "\n");;
		 }
		 fout2.flush();
		 
		 long tellpSz2 = (long)fout2.tellp();
		 VMSG_NL("Wrote " + fout2.tellp() + " bytes to secondary EBWT file: " + _in2Str);
		 fout2.close();
		 if(tellpSz2 > fileSize(_in2Str)) {
		 err = true;
		 System.err.println("Index is corrupt: File size for " + _in2Str + " should have been " + tellpSz2
		 + " but is actually " + fileSize(_in2Str) + "." + "\n");;
		 }
		 
		 if(saOut != null) {
		 // Check on suffix array output file size
		 long tellpSzSa = (long)saOut.tellp();
		 VMSG_NL("Wrote " + tellpSzSa + " bytes to suffix-array file: " + _inSaStr);
		 saOut.close();
		 if(tellpSzSa > fileSize(_inSaStr)) {
		 err = true;
		 System.err.println("Index is corrupt: File size for " + _inSaStr + " should have been " + tellpSzSa
		 + " but is actually " + fileSize(_inSaStr) + "." + "\n");;
		 }
		 }
		 
		 if(bwtOut != null) {
		 // Check on suffix array output file size
		 long tellpSzBwt = (long)bwtOut.tellp();
		 VMSG_NL("Wrote " + tellpSzBwt + " bytes to BWT file: " + _inBwtStr);
		 bwtOut.close();
		 if(tellpSzBwt > fileSize(_inBwtStr)) {
		 err = true;
		 System.err.println("Index is corrupt: File size for " + _inBwtStr + " should have been " + tellpSzBwt
		 + " but is actually " + fileSize(_inBwtStr) + "." + "\n");;
		 }
		 }
		 
		 if(err) {
		 System.err.println("Please check if there is a problem with the disk or if disk is full." + "\n");;
		 throw 1;
		 }
		 **/
		// Reopen as input streams
		//VMSG_NL("Re-opening _in1 and _in2 as input streams");
		if(_sanity) {
			//VMSG_NL("Sanity-checking Bt2");
			assert(!isInMemory());
			readIntoMemory(
					color,                       // colorspace?
					fw ? -1 : needEntireReverse, // 1 . need the reverse to be reverse-of-concat
					true,                        // load SA sample (_offs[])?
					true,                        // load ftab (_ftab[] & _eftab[])?
					true,                        // load r-starts (_rstarts[])?
					false,                       // just load header?
					null,                        // Params object to fill
					false,                       // mm sweep?
					true,                        // load names?
					false);                      // verbose startup?
			sanityCheckAll(refparams.reverse);
			evictFromMemory();
			assert(!isInMemory());
		}
		//VMSG_NL("Returning from Ebwt constructor");
	}
	
	
	public TStr join(EList<TStr> l, int seed) {
		RandomSource rand = new RandomSource(); // reproducible given same seed
		rand.init(seed);
		TStr ret;
		long guessLen = 0;
		for(int i = 0; i < l.size(); i++) {
			guessLen += l.size();
			//guessLen += length(l.get(i));
		}
		ret.resize(guessLen);
		long off = 0;
		for(int i = 0; i < l.size(); i++) {
			TStr s = l[i];
			for(int j = 0; j < s.size(); j++) {
				ret.set(s.get(j), off++);
			}
		}
		return ret;
	}
	public String join(
			EList<FileBuf> l,
			EList<RefRecord> szs,
			long sztot,
			RefReadInParams refparams,
			int seed) {
		RandomSource rand = new RandomSource(); // reproducible given same seed
		rand.init(seed);
		RefReadInParams rpcp = refparams;
		TStr ret;
		long guessLen = sztot;
		ret.resize(guessLen);
		long dstoff = 0;
		for(int i = 0; i < l.size(); i++) {
			// For each sequence we can pull out of istream l[i]...
			boolean first = true;
			while(!l.get(i).eof()) {
				RefRecord rec = fastaRefReadAppend(l.get(i), first, ret, dstoff, rpcp);
				first = false;
				long bases = rec.len;
				if(bases == 0) continue;
			}
		}
		return ret.toString();
	}
	
	public long joinedLen(EList<RefRecord> szs) {
		long ret = 0;
		for(int i = 0; i < szs.size();i++) {
			ret += szs.get(i).len;
		}
		return ret;
	}
	
	public void joinToDisk(
			EList<FileBuf> l,
			EList<RefRecord> szs,
			long sztot,
			RefReadInParams refparams,
			TStr ret,
			OutputStream out1,
			OutputStream out2
	) {
		RefReadInParams rpcp = refparams;
		// Not every fragment represents a distinct sequence - many
		// fragments may correspond to a single sequence.  Count the
		// number of sequences here by counting the number of "first"
		// fragments.
		this._nPat = 0;
		this._nFrag = 0;
		for(int i = 0; i < szs.size(); i++) {
			if(szs.get(i).len > 0) this._nFrag++;
			if(szs.get(i).first && szs.get(i).len > 0) this._nPat++;
		}
		//_rstarts.reset();
		writeU<long>(out1, this._nPat, this.toBe());
		// Allocate plen[]
		try {
			this._plen.init(new long[this._nPat], this._nPat);
		} catch(Exception e) {
			System.err.println("Out of memory allocating plen[] in Ebwt.join()"
					+ " at " + __FILE__ + ":" + __LINE__ + "\n");
			throw e;
		}
		// For each pattern, set plen
		long npat = -1;
		for(long i = 0; i < szs.size(); i++) {
			if(szs.get((int)i).first && szs.get((int)i).len > 0) {
				if(npat >= 0) {
					writeU<long>(out1, this.plen()[(int)npat], this.toBe());
				}
				this.plen()[(int)++npat] = (szs.get((int)i).len + szs.get((int)i).off);
			} else {
				// edge case, but we could get here with npat == -1
				// e.g. when building from a reference of all Ns
				if (npat < 0) npat = 0;
				this.plen()[(int)npat] += (szs.get((int)i).len + szs.get((int)i).off);
			}
		}
		writeU<long>(out1, this.plen()[(int)npat], this.toBe());
		// Write the number of fragments
		writeU<long>(out1, this._nFrag, this.toBe());
		long seqsRead = 0;
		long dstoff = 0;
		// For each filebuf
		for(int i = 0; i < l.size(); i++) {
			boolean first = true;
			long patoff = 0;
			// For each *fragment* (not necessary an entire sequence) we
			// can pull out of istream l[i]...
			while(!l.get(i).eof()) {
				String name;
				// Push a new name onto our vector
				_refnames.push_back("");
				RefRecord rec = fastaRefReadAppend(
						l.get(i), first, ret, dstoff, rpcp, _refnames.back());
				first = false;
				long bases = rec.len;
				if(rec.first && rec.len > 0) {
					if(_refnames.back().length() == 0) {
						// If name was empty, replace with an index
						try
						{
							FileOutputStream stm = new FileOutputStream("");
							stm.write((int)seqsRead);
							_refnames.back() = stm.str();
						}
						catch (IOException e)
						{
						
						}
						
					}
				} else {
					// This record didn't actually start a new sequence so
					// no need to add a name
					////assert_eq(0, _refnames.back().length());
					_refnames.pop_back();
				}
				// Increment seqsRead if this is the first fragment
				if(rec.first && rec.len > 0) seqsRead++;
				if(bases == 0) continue;
				//assert_leq(bases, this.plen()[seqsRead-1]);
				// Reset the patoff if this is the first fragment
				if(rec.first) patoff = 0;
				patoff += rec.off; // add fragment's offset from end of last frag.
				// Adjust rpcps
				//long seq = seqsRead-1;
				// This is where rstarts elements are written to the output stream
				//writeU32(out1, oldRetLen, this.toBe()); // offset from beginning of joined String
				//writeU32(out1, seq,       this.toBe()); // sequence id
				//writeU32(out1, patoff,    this.toBe()); // offset into sequence
				patoff += bases;
			}
			l.get(i).reset();
		}
	}
	
	public void buildToDisk(
			InorderBlockwiseSA<TStr> sa,
			TStr s,
			OutputStream out1,
			OutputStream out2,
			OutputStream saOut,
			OutputStream bwtOut
	) {
		EbwtParams eh = _eh;
		
		//assert(eh.repOk());
		//assert_eq(s.length()+1, sa.size());
		//assert_eq(s.length(), eh._len);
		//assert_gt(eh._lineRate, 3);
		//assert(sa.suffixItrIsReset());
		
		long len = eh._len;
		long ftabLen = eh._ftabLen;
		long sideSz = eh._sideSz;
		long ebwtTotSz = eh._ebwtTotSz;
		long fchr[] = {0, 0, 0, 0, 0};
		EList<long> ftab(1);
		long zOff = IndexTypes.OFF_MASK;
		
		// Save # of occurrences of each character as we walk along the bwt
		long occ[] = {0, 0, 0, 0};
		long occSave[] = {0, 0, 0, 0};
		
		// Record rows that should "absorb" adjacent rows in the ftab.
		// The absorbed rows represent suffixes shorter than the ftabChars
		// cutoff.
		byte absorbCnt = 0;
		EList<byte> absorbFtab(1);
		try {
			//VMSG_NL("Allocating ftab, absorbFtab");
			ftab.resize((int)ftabLen);
			ftab.fillZero();
			absorbFtab.resize(ftabLen);
			absorbFtab.fillZero();
		} catch(Exception e) {
			System.err.println("Out of memory allocating ftab[] or absorbFtab[] "
					+ "in Ebwt.buildToDisk() at " + __FILE__ + ":"
					+ __LINE__ + "\n");;
			throw e;
		}
		
		// Allocate the side buffer; holds a single side as its being
		// constructed and then written to disk.  Reused across all sides.
		/**
		 #ifdef SIXTY4_FORMAT
		 EList<long> ebwtSide(1);
		 #else
		 EList<byte> ebwtSide(1);
		 #endif
		 try {
		 #ifdef SIXTY4_FORMAT
		 ebwtSide.resize(sideSz >> 3);
		 #else
		 ebwtSide.resize(sideSz);
		 #endif
		 } catch(bad_alloc &e) {
		 System.err.println("Out of memory allocating ebwtSide[] in "
		 + "Ebwt.buildToDisk() at " + __FILE__ + ":"
		 + __LINE__ + "\n");;
		 throw e;
		 }
		 **/
		// Points to the base offset within ebwt for the side currently
		// being written
		long side = 0;
		
		// Whether we're assembling a forward or a reverse bucket
		boolean fw;
		long sideCur = 0;
		fw = true;
		
		// Have we skipped the '$' in the last column yet?
		//assert_ONLY(boolean dollarSkipped = false);
		
		long si = 0;   // String offset (chars)
		//assert_ONLY(long lastSufInt = 0);
		//assert_ONLY(boolean inSA = true); // true iff saI still points inside suffix
		// array (as opposed to the padding at the
		// end)
		// Iterate over packed bwt bytes
		//VMSG_NL("Entering Ebwt loop");
		//assert_ONLY(long beforeEbwtOff = (long)out1.tellp()); // @double-check - pos_type, std.streampos
		
		// First integer in the suffix-array output file is the length of the
		// array, including $
		if(saOut != null) {
			// Write length word
			writeU<long>(*saOut, len+1, this.toBe());
		}
		
		// First integer in the BWT output file is the length of BWT(T), including $
		if(bwtOut != null) {
			// Write length word
			writeU<long>(*bwtOut, len+1, this.toBe());
		}
		
		while(side < ebwtTotSz) {
			// Sanity-check our cursor into the side buffer
			//assert_geq(sideCur, 0);
			//assert_lt(sideCur, (int)eh._sideBwtSz);
			//assert_eq(0, side % sideSz); // 'side' must be on side boundary
			ebwtSide[sideCur] = 0; // clear
			//assert_lt(side + sideCur, ebwtTotSz);
			// Iterate over bit-pairs in the si'th character of the BWT
	#ifdef SIXTY4_FORMAT
			for(int bpi = 0; bpi < 32; bpi++, si++)
	#else
			for(int bpi = 0; bpi < 4; bpi++, si++)
	#endif
			{
				int bwtChar;
				boolean count = true;
				if(si <= len) {
					// Still in the SA; extract the bwtChar
					long saElt = sa.nextSuffix();
					// Write it to the optional suffix-array output file
					if(saOut != null) {
						writeU<long>(*saOut, saElt, this.toBe());
					}
					// TODO: what exactly to write to the BWT output file?  How to
					// represent $?  How to pack nucleotides into bytes/words?
					
					// (that might have triggered sa to calc next suf block)
					if(saElt == 0) {
						// Don't add the '$' in the last column to the BWT
						// transform; we can't encode a $ (only A C T or G)
						// and counting it as, say, an A, will mess up the
						// LR mapping
						bwtChar = 0; count = false;
						//assert_ONLY(dollarSkipped = true);
						zOff = si; // remember the SA row that
						// corresponds to the 0th suffix
					} else {
						bwtChar = (int)(s[saElt-1]);
						//assert_lt(bwtChar, 4);
						// Update the fchr
						fchr[bwtChar]++;
					}
					// Update ftab
					if((len-saElt) >= (long)eh._ftabChars) {
						// Turn the first ftabChars characters of the
						// suffix into an integer index into ftab.  The
						// leftmost (lowest index) character of the suffix
						// goes in the most significant bit pair if the
						// integer.
						long sufInt = 0;
						for(int i = 0; i < eh._ftabChars; i++) {
							sufInt += 2;
							//assert_lt((long)i, len-saElt);
							sufInt |= (unsigned char)(s[saElt+i]);
						}
						// Assert that this prefix-of-suffix is greater
						// than or equal to the last one (true b/c the
						// suffix array is sorted)
						#ifndef NDEBUG
						if(lastSufInt > 0) //assert_geq(sufInt, lastSufInt);
							lastSufInt = sufInt;
						#endif
						// Update ftab
						//assert_lt(sufInt+1, ftabLen);
						ftab[sufInt+1]++;
						if(absorbCnt > 0) {
							// Absorb all short suffixes since the last
							// transition into this transition
							absorbFtab[sufInt] = absorbCnt;
							absorbCnt = 0;
						}
					} else {
						// Otherwise if suffix is fewer than ftabChars
						// characters long, then add it to the 'absorbCnt';
						// it will be absorbed into the next transition
						//assert_lt(absorbCnt, 255);
						absorbCnt++;
					}
					// Suffix array offset boundary? - update offset array
					if((si & eh._offMask) == si) {
						//assert_lt((si >> eh._offRate), eh._offsLen);
						// Write offsets directly to the secondary output
						// stream, thereby avoiding keeping them in memory
						writeU<long>(out2, saElt, this.toBe());
					}
				} else {
					// Strayed off the end of the SA, now we're just
					// padding out a bucket
					#ifndef NDEBUG
					if(inSA) {
						// Assert that we wrote all the characters in the
						// String before now
						//assert_eq(si, len+1);
						inSA = false;
					}
					#endif
							// 'A' used for padding; important that padding be
							// counted in the occ[] array
							bwtChar = 0;
				}
				if(count) occ[bwtChar]++;
				// Append BWT char to bwt section of current side
				if(fw) {
					// Forward bucket: fill from least to most
	#ifdef SIXTY4_FORMAT
					ebwtSide[sideCur] |= ((long)bwtChar + (bpi + 1));
					if(bwtChar > 0) //assert_gt(ebwtSide[sideCur], 0);
	#else
					pack_2b_in_8b(bwtChar, ebwtSide[sideCur], bpi);
					//assert_eq((ebwtSide[sideCur] >> (bpi*2)) & 3, bwtChar);
	#endif
				} else {
					// Backward bucket: fill from most to least
	#ifdef SIXTY4_FORMAT
					ebwtSide[sideCur] |= ((long)bwtChar + ((31 - bpi) + 1));
					if(bwtChar > 0) //assert_gt(ebwtSide[sideCur], 0);
	#else
					pack_2b_in_8b(bwtChar, ebwtSide[sideCur], 3-bpi);
					//assert_eq((ebwtSide[sideCur] >> ((3-bpi)*2)) & 3, bwtChar);
	#endif
				}
			} // end loop over bit-pairs
			//assert_eq(dollarSkipped ? 3 : 0, (occ[0] + occ[1] + occ[2] + occ[3]) & 3);
	#ifdef SIXTY4_FORMAT
			//assert_eq(0, si & 31);
	#else
			//assert_eq(0, si & 3);
	#endif
			
			sideCur++;
			if(sideCur == (int)eh._sideBwtSz) {
				sideCur = 0;
				long *cpptr = reinterpret_cast<long*>(ebwtSide.ptr());
				// Write 'A', 'C', 'G' and 'T' tallies
				side += sideSz;
				//assert_leq(side, eh._ebwtTotSz);
	#ifdef BOWTIE_64BIT_INDEX
				cpptr[(sideSz >> 3)-4] = endianizeU<long>(occSave[0], this.toBe());
				cpptr[(sideSz >> 3)-3] = endianizeU<long>(occSave[1], this.toBe());
				cpptr[(sideSz >> 3)-2] = endianizeU<long>(occSave[2], this.toBe());
				cpptr[(sideSz >> 3)-1] = endianizeU<long>(occSave[3], this.toBe());
	#else
				cpptr[(sideSz >> 2)-4] = endianizeU<long>(occSave[0], this.toBe());
				cpptr[(sideSz >> 2)-3] = endianizeU<long>(occSave[1], this.toBe());
				cpptr[(sideSz >> 2)-2] = endianizeU<long>(occSave[2], this.toBe());
				cpptr[(sideSz >> 2)-1] = endianizeU<long>(occSave[3], this.toBe());
	#endif
				occSave[0] = occ[0];
				occSave[1] = occ[1];
				occSave[2] = occ[2];
				occSave[3] = occ[3];
				// Write backward side to primary file
				out1.write((char *)ebwtSide.ptr(), sideSz);
			}
		}
		//VMSG_NL("Exited Ebwt loop");
		//assert_neq(zOff, IndexTypes.OFF_MASK);
		if(absorbCnt > 0) {
			// Absorb any trailing, as-yet-unabsorbed short suffixes into
			// the last element of ftab
			absorbFtab[ftabLen-1] = absorbCnt;
		}
		// Assert that our loop counter got incremented right to the end
		//assert_eq(side, eh._ebwtTotSz);
		// Assert that we wrote the expected amount to out1
		//assert_eq(((long)out1.tellp() - beforeEbwtOff), eh._ebwtTotSz); // @double-check - pos_type
		// assert that the last thing we did was write a forward bucket
		
		//
		// Write zOff to primary stream
		//
		writeU<long>(out1, zOff, this.toBe());
		
		//
		// Finish building fchr
		//
		// Exclusive prefix sum on fchr
		for(int i = 1; i < 4; i++) {
			fchr[i] += fchr[i-1];
		}
		//assert_eq(fchr[3], len);
		// Shift everybody up by one
		for(int i = 4; i >= 1; i--) {
			fchr[i] = fchr[i-1];
		}
		fchr[0] = 0;
		if(_verbose) {
			for(int i = 0; i < 5; i++)
				System.out.println( "fchr[" + "ACGT$" + i + "]: " + fchr[i] + "\n");
		}
		// Write fchr to primary file
		for(int i = 0; i < 5; i++) {
			writeU<long>(out1, fchr[i], this.toBe());
		}
		
		//
		// Finish building ftab and build eftab
		//
		// Prefix sum on ftable
		long eftabLen = 0;
		//assert_eq(0, absorbFtab[0]);
		for(long i = 1; i < ftabLen; i++) {
			if(absorbFtab[i] > 0) eftabLen += 2;
		}
		//assert_leq(eftabLen, (long)eh._ftabChars*2);
		eftabLen = eh._ftabChars*2;
		EList<long> eftab(1);
		try {
			eftab.resize(eftabLen);
			eftab.fillZero();
		} catch(Exception e) {
			System.err.println("Out of memory allocating eftab[] "
					+ "in Ebwt.buildToDisk() at " + __FILE__ + ":"
					+ __LINE__ + "\n");;
			throw e;
		}
		long eftabCur = 0;
		for(long i = 1; i < ftabLen; i++) {
			long lo = ftab[i] + Ebwt.ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i-1);
			if(absorbFtab[i] > 0) {
				// Skip a number of short pattern indicated by absorbFtab[i]
				long hi = lo + absorbFtab[i];
				//assert_lt(eftabCur*2+1, eftabLen);
				eftab[eftabCur*2] = lo;
				eftab[eftabCur*2+1] = hi;
				ftab[i] = (eftabCur++) ^ IndexTypes.OFF_MASK; // insert pointer into eftab
				//assert_eq(lo, Ebwt.ftabLo(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i));
				//assert_eq(hi, Ebwt.ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, i));
			} else {
				ftab[i] = lo;
			}
		}
		////assert_eq(Ebwt.ftabHi(ftab.ptr(), eftab.ptr(), len, ftabLen, eftabLen, ftabLen-1), len+1);
		// Write ftab to primary file
		for(long i = 0; i < ftabLen; i++) {
			writeU<long>(out1, ftab[i], this.toBe());
		}
		// Write eftab to primary file
		for(long i = 0; i < eftabLen; i++) {
			writeU<long>(out1, eftab[i], this.toBe());
		}
		
		// Note: if you'd like to sanity-check the Ebwt, you'll have to
		// read it back into memory first!
		assert(!isInMemory());
		//VMSG_NL("Exiting Ebwt.buildToDisk()");
	}
	
	public void joinedToTextOff(
			long qlen,
			long off,
			int tidx,
			long textoff,
			long tlen,
			boolean rejectStraddle,
			boolean straddled){
		int top = 0;
		int bot = _nFrag; // 1 greater than largest addressable element
		int elt = (int) IndexTypes.OFF_MASK;
		// Begin binary search
		while(true) {
			elt = top + ((bot - top) >> 1);
			long lower = rstarts()[elt*3];
			long upper;
			if(elt == _nFrag-1) {
				upper = _eh._len;
			} else {
				upper = rstarts()[((elt+1)*3)];
			}
			long fraglen = upper - lower;
			if(lower <= off) {
				if(upper > off) { // not last element, but it's within
					// off is in this range; check if it falls off
					if(off + qlen > upper) {
						straddled = true;
						if(rejectStraddle) {
							// it falls off; signal no-go and return
							tidx = (int)IndexTypes.OFF_MASK;
							return;
						}
					}
					// This is the correct text idx whether the index is
					// forward or reverse
					tidx = rstarts()[(elt*3)+1];
					// it doesn't fall off; now calculate textoff.
					// Initially it's the number of characters that precede
					// the alignment in the fragment
					long fragoff = off - rstarts()[(elt*3)];
					if(!this.fw_) {
						fragoff = fraglen - fragoff - 1;
						fragoff -= (qlen-1);
					}
					// Add the alignment's offset into the fragment
					// ('fragoff') to the fragment's offset within the text
					textoff = fragoff + rstarts()[(elt*3)+2];
					break; // done with binary search
				} else {
					// 'off' belongs somewhere in the region between elt
					// and bot
					top = elt;
				}
			} else {
				// 'off' belongs somewhere in the region between top and
				// elt
				bot = elt;
			}
			// continue with binary search
		}
		tlen =plen()[tidx];
	}
	
	public long walkLeft(long row, long steps) {
		SideLocus l;
		if(steps > 0) l.initFromRow(row, _eh, ebwt());
		while(steps > 0) {
			if(row == _zOff) return IndexTypes.OFF_MASK;
			long newrow = mapLF(l);
			row = newrow;
			steps--;
			if(steps > 0) l.initFromRow(row, _eh, ebwt());
		}
		return row;
	}
	
	public long getOffset(int row) {
		if(row == _zOff) return 0;
		if((row & _eh._offMask) == row) return offs()[row >> _eh._offRate];
		long jumps = 0;
		SideLocus l;
		l.initFromRow(row, _eh, ebwt());
		while(true) {
			int newrow = (int)mapLF(l);
			jumps++;
			row = newrow;
			if(row == _zOff) {
				return jumps;
			} else if((row & _eh._offMask) == row) {
				return jumps + this.offs()[row >> _eh._offRate];
			}
			l.initFromRow(row, _eh, ebwt());
		}
	}
	
	public long getOffset(int elt, boolean fw, long hitlen) {
		long off = getOffset(elt);
		if(!fw) {
			off = _eh._len - off - 1;
			off -= (hitlen-1);
		}
		return off;
	}
	
	public boolean contains(String str, long[] otop, long[] obot) {
		SideLocus tloc, bloc;
		if(str == "") {
			if(otop != null && obot != null) {
				otop = null;
				obot = null;
			}
			return true;
		}
		int c = str.charAt(str.length()-1);
		long top = 0, bot = 0;
		if(c < 4) {
			top = fchr()[c];
			bot = fchr()[c+1];
		} else {
			boolean set = false;
			for(int i = 0; i < 4; i++) {
				if(fchr()[c] < fchr()[c+1]) {
					if(set) {
						return false;
					} else {
						set = true;
						top = fchr()[c];
						bot = fchr()[c+1];
					}
				}
			}
		}
		tloc.initFromRow(top, eh(), ebwt());
		bloc.initFromRow(bot, eh(), ebwt());
		for(int i = str.length()-2; i >= 0; i--) {
			c = str.charAt(i);
			if(c <= 3) {
				top = mapLF(tloc, c);
				bot = mapLF(bloc, c);
			} else {
				long sz = bot - top;
				int c1 = mapLF1(top, tloc);
				bot = mapLF(bloc, c1);
				if(bot - top < sz) {
					// Encountered an N and could not proceed through it because
					// there was more than one possible nucleotide we could replace
					// it with
					return false;
				}
			}
			if(i > 0) {
				tloc.initFromRow(top, eh(), ebwt());
				bloc.initFromRow(bot, eh(), ebwt());
			}
		}
		if(otop != null && obot != null) {
			otop = top; obot = bot;
		}
		return bot > top;
	}
	
	public String adjustEbwtBase(String cmdline, String ebwtFileBase, boolean verbose) {
		String str = ebwtFileBase;
		File in = new File((str + ".1." + gEbwt_ext));
		
		if(verbose) System.out.println( "Trying " + str);
		if(!in.exists())
			if(verbose) System.out.println( "  didn't work" );
		if(System.getenv("BOWTIE2_INDEXES") != null) {
			str = System.getenv("BOWTIE2_INDEXES") + "/" + ebwtFileBase;
			if(verbose) System.out.println( "Trying " + str);
			in=new File((str + ".1." + gEbwt_ext));
		}
		if(!in.exists()) {
			System.err.println("Could not locate a Bowtie index corresponding to basename \"" + ebwtFileBase + "\"" );
		}
		return str;
	}
	
	public long    zOff()         { return _zOff; }
	public long    zEbwtByteOff() { return _zEbwtByteOff; }
	public long    zEbwtBpOff()   { return _zEbwtBpOff; }
	public long    nPat()         { return _nPat; }
	public long    nFrag()        { return _nFrag; }
	public long   ftab()              { return _ftab; }
	public long   eftab()             { return _eftab; }
	public long[]   offs()              { return _offs; }
	public long[]   plen()              { return _plen; }
	public int[]   rstarts()           { return _rstarts; }
	public boolean        sanityCheck()  { return _sanity; }
	public EList<String> refnames()        { return _refnames; }
	public boolean        fw()           { return fw_; }
	
	public void readIntoMemory(
			int color,
			int needEntireRev,
			boolean loadSASamp,
			boolean loadFtab,
			boolean loadRstarts,
			boolean justHeader,
			EbwtParams params,
			boolean mmSweep,
			boolean loadNames,
			boolean startVerbose) {
		
	}
	
	public void writeFromMemory(boolean justHeader, OutputStream out1, OutputStream out2) {
		EbwtParams eh = this._eh;
		//assert(eh.repOk());
		long be = this.toBe();
		//assert(out1.good());
		//assert(out2.good());
		
		// When building an Ebwt, these header parameters are known
		// "up-front", i.e., they can be written to disk immediately,
		// before we join() or buildToDisk()
		writeI<long>(out1, 1, be); // endian hint for priamry stream
		writeI<long>(out2, 1, be); // endian hint for secondary stream
		writeU<long>(out1, eh._len,          be); // length of string (and bwt and suffix array)
		writeI<long>(out1, eh._lineRate,     be); // 2^lineRate = size in bytes of 1 line
		writeI<long>(out1, 2,                be); // not used
		writeI<long>(out1, eh._offRate,      be); // every 2^offRate chars is "marked"
		writeI<long>(out1, eh._ftabChars,    be); // number of 2-bit chars used to address ftab
		long flags = 1;
		if(eh._color) flags |= EBWT_COLOR;
		if(eh._entireReverse) flags |= EBWT_ENTIRE_REV;
		writeI<long>(out1, -flags, be); // BTL: chunkRate is now deprecated
		
		if(!justHeader) {
			//assert(rstarts() != null);
			//assert(offs() != null);
			//assert(ftab() != null);
			//assert(eftab() != null);
			//assert(isInMemory());
			// These Ebwt parameters are known after the inputs strings have
			// been joined() but before they have been built().  These can
			// written to the disk next and then discarded from memory.
			writeU<long>(out1, this._nPat,      be);
			for(long i = 0; i < this._nPat; i++)
				writeU<long>(out1, this.plen()[i], be);
			//assert_geq(this._nFrag, this._nPat);
			writeU<long>(out1, this._nFrag, be);
			for(long i = 0; i < this._nFrag*3; i++)
				writeU<long>(out1, this.rstarts()[i], be);
			
			// These Ebwt parameters are discovered only as the Ebwt is being
			// built (in buildToDisk()).  Of these, only 'offs' and 'ebwt' are
			// terribly large.  'ebwt' is written to the primary file and then
			// discarded from memory as it is built; 'offs' is similarly
			// written to the secondary file and discarded.
			out1.write((const char *)this.ebwt(), eh._ebwtTotLen);
			writeU<long>(out1, this.zOff(), be);
			long offsLen = eh._offsLen;
			for(long i = 0; i < offsLen; i++)
				writeU<long>(out2, this.offs()[i], be);
			
			// 'fchr', 'ftab' and 'eftab' are not fully determined until the
			// loop is finished, so they are written to the primary file after
			// all of 'ebwt' has already been written and only then discarded
			// from memory.
			for(int i = 0; i < 5; i++)
				writeU<long>(out1, this.fchr()[i], be);
			for(long i = 0; i < eh._ftabLen; i++)
				writeU<long>(out1, this.ftab()[i], be);
			for(long i = 0; i < eh._eftabLen; i++)
				writeU<long>(out1, this.eftab()[i], be);
		}
	}
	
	public int readFlags(String instr) {
		ifstream in;
		// Initialize our primary and secondary input-stream fields
		in.open((instr + ".1." + gEbwt_ext).c_str(), ios_base.in | ios.binary);
		if(!in.is_open()) {
			throw EbwtFileOpenException("Cannot open file " + instr);
		}
		assert(in.is_open());
		assert(in.good());
		Boolean switchEndian = false;
		long one = readU<long>(in, switchEndian); // 1st word of primary stream
		if(one != 1) {
			//assert_eq((1u<<24), one);
			//assert_eq(1, endianSwapU32(one));
			switchEndian = true;
		}
		readU<long>(in, switchEndian);
		readI<long>(in, switchEndian);
		readI<long>(in, switchEndian);
		readI<long>(in, switchEndian);
		readI<long>(in, switchEndian);
		long flags = readI<long>(in, switchEndian);
		return (int)flags;
	}
	
	public static int readFlags(String instr) {
		Boolean switchEndian; // dummy; caller doesn't care
		#ifdef BOWTIE_MM
			char *mmFile[] = { null, null };
		#endif
		if(_in1Str.length() > 0) {
			if(_verbose || startVerbose) {
				System.err.println("About to open input files: ");
				logTime(cerr);
			}
			// Initialize our primary and secondary input-stream fields
			if(_in1 != null) fclose(_in1);
			if(_verbose || startVerbose) System.err.println( "Opening \"" +  _in1Str.c_str() + "\"");
			if((_in1 = fopen(_in1Str.c_str(), "rb")) == null) {
				System.err.println( "Could not open index file " + _in1Str.c_str() );
			}
			if(loadSASamp) {
				if(_in2 != null) fclose(_in2);
				if(_verbose || startVerbose) System.err.println( "Opening \"" + _in2Str.c_str() + "\"");
				if((_in2 = fopen(_in2Str.c_str(), "rb")) == null) {
					System.err.println( "Could not open index file " + _in2Str.c_str());
				}
			}
			if(_verbose || startVerbose) {
				System.err.println( "  Finished opening input files: ";
				logTime(cerr);
			}
			
		#ifdef BOWTIE_MM
			if(_useMm /*&& !justHeader*/) {
					char *names[] = {_in1Str.c_str(), _in2Str.c_str()};
				int fds[] = { fileno(_in1), fileno(_in2) };
				for(int i = 0; i < (loadSASamp ? 2 : 1); i++) {
					if(_verbose || startVerbose) {
						System.err.println( "  Memory-mapping input file " + (i+1) + ": ");
						//logTime(cerr);
					}
					struct stat sbuf;
					if (stat(names[i], &sbuf) == -1) {
						//perror("stat");
						System.err.println( "Error: Could not stat index file " + names[i] + " prior to memory-mapping");
						//throw 1;
					}
					mmFile[i] = (char*)mmap((void *)0, (int)sbuf.st_size,
							PROT_READ, MAP_SHARED, fds[(int)i], 0);
					if(mmFile[i] == (void *)(-1)) {
						perror("mmap");
						System.err.println( "Error: Could not memory-map the index file " << names[i] << endl;
						//throw 1;
					}
					if(mmSweep) {
						int sum = 0;
						for(off_t j = 0; j < sbuf.st_size; j += 1024) {
							sum += (int) mmFile[i][j];
						}
						if(startVerbose) {
							System.err.println( "  Swept the memory-mapped ebwt index file 1; checksum: " + sum + ": ");
							//logTime(cerr);
						}
					}
				}
				mmFile1_ = mmFile[0];
				mmFile2_ = loadSASamp ? mmFile[1] : null;
			}
		#endif
		}
		#ifdef BOWTIE_MM
			else if(_useMm && !justHeader) {
			mmFile[0] = mmFile1_;
			mmFile[1] = mmFile2_;
		}
		if(_useMm && !justHeader) {
			assert(mmFile[0] == mmFile1_);
			assert(mmFile[1] == mmFile2_);
		}
		#endif
		
		if(_verbose || startVerbose) {
			System.err.println( "  Reading header: ");
			//logTime(cerr);
		}
		
		// Read endianness hints from both streams
		long bytesRead = 0;
		switchEndian = false;
		long one = readU<long>(_in1, switchEndian); // 1st word of primary stream
		bytesRead += 4;
		if(loadSASamp) {
		#ifndef NDEBUG
			assert_eq(one, readU<long>(_in2, switchEndian)); // should match!
		#else
			readU<long>(_in2, switchEndian);
		#endif
		}
		if(one != 1) {
			assert_eq((1u<<24), one);
			assert_eq(1, endianSwapU32(one));
			switchEndian = true;
		}
		
		// Can't switch endianness and use memory-mapped files; in order to
		// support this, someone has to modify the file to switch
		// endiannesses appropriately, and we can't do this inside Bowtie
		// or we might be setting up a race condition with other processes.
		if(switchEndian && _useMm) {
			System.err.println( "Error: Can't use memory-mapped files when the index is the opposite endianness");
			//throw 1;
		}
		
		// Reads header entries one by one from primary stream
		long len          = readU<long>(_in1, switchEndian);
		bytesRead += OFF_SIZE;
		long  lineRate     = readI<long>(_in1, switchEndian);
		bytesRead += 4;
			/*int32_t  linesPerSide =*/ readI<long>(_in1, switchEndian);
		bytesRead += 4;
		long  offRate      = readI<long>(_in1, switchEndian);
		bytesRead += 4;
		// TODO: add isaRate to the actual file format (right now, the
		// user has to tell us whether there's an ISA sample and what the
		// sampling rate is.
		long  ftabChars    = readI<long>(_in1, switchEndian);
		bytesRead += 4;
		// chunkRate was deprecated in an earlier version of Bowtie; now
		// we use it to hold flags.
		long flags = readI<long>(_in1, switchEndian);
		Boolean entireRev = false;
		if(flags < 0 && (((-flags) & EBWT_COLOR) != 0)) {
			if(color != -1 && !color) {
				System.err.println( "Error: -C was not specified when running bowtie, but index is in colorspace.  If"
						+ "your reads are in colorspace, please use the -C option.  If your reads are not"
						+ "in colorspace, please use a normal index (one built without specifying -C to"
						+ "bowtie-build).");
				//throw 1;
			}
			color = 1;
		} else if(flags < 0) {
			if(color != -1 && color) {
				System.err.println( "Error: -C was specified when running bowtie, but index is not in colorspace.  If"
						+ "your reads are in colorspace, please use a colorspace index (one built using"
						+ "bowtie-build -C).  If your reads are not in colorspace, don't specify -C when"
						+ "running bowtie.");
				throw 1;
			}
			color = 0;
		}
		if(flags < 0 && (((-flags) & EBWT_ENTIRE_REV) == 0)) {
			if(needEntireRev != -1 && needEntireRev != 0) {
				System.err.println( "Error: This index is compatible with 0.* versions of Bowtie, but not with 2.*"
						+ "versions.  Please build or download a version of the index that is compitble"
						+ "with Bowtie 2.* (i.e. built with bowtie-build 2.* or later)");
				//throw 1;
			}
		} else entireRev = true;
		bytesRead += 4;
		
		// Create a new EbwtParams from the entries read from primary stream
		EbwtParams *eh;
		Boolean deleteEh = false;
		if(params != null) {
			params.init(len, lineRate, offRate, ftabChars, color, entireRev);
			if(_verbose || startVerbose) params.print(cerr);
			eh = params;
		} else {
			eh = new EbwtParams(len, lineRate, offRate, ftabChars, color, entireRev);
			deleteEh = true;
		}
		
		// Set up overridden suffix-array-sample parameters
		long offsLen = eh._offsLen;
		long offsSz = eh._offsSz;
		long offRateDiff = 0;
		long offsLenSampled = offsLen;
		if(_overrideOffRate > offRate) {
			offRateDiff = _overrideOffRate - offRate;
		}
		if(offRateDiff > 0) {
			offsLenSampled >>= offRateDiff;
			if((offsLen & ~(OFF_MASK << offRateDiff)) != 0) {
				offsLenSampled++;
			}
		}
		
		// Can't override the offrate or isarate and use memory-mapped
		// files; ultimately, all processes need to copy the sparser sample
		// into their own memory spaces.
		if(_useMm && (offRateDiff)) {
			System.err.println( "Error: Can't use memory-mapped files when the offrate is overridden");
			//throw 1;
		}
		
		// Read nPat from primary stream
		this._nPat = readI<long>(_in1, switchEndian);
		bytesRead += OFF_SIZE;
		_plen.reset();
		// Read plen from primary stream
		if(_useMm) {
		#ifdef BOWTIE_MM
			_plen.init((long*)(mmFile[0] + bytesRead), _nPat, false);
			bytesRead += _nPat*OFF_SIZE;
			fseeko(_in1, _nPat*OFF_SIZE, SEEK_CUR);
		#endif
		} else {
			try {
				if(_verbose || startVerbose) {
					System.err.println( "Reading plen (" + this._nPat + "): ";
					logTime(cerr);
				}
				_plen.init(new long[_nPat], _nPat, true);
				if(switchEndian) {
					for(long i = 0; i < this._nPat; i++) {
						plen()[i] = readU<long>(_in1, switchEndian);
					}
				} else {
					int r = MM_READ(_in1, (void*)(plen()), _nPat*OFF_SIZE);
					if(r != (int)(_nPat*OFF_SIZE)) {
						System.err.println( "Error reading _plen[] array: " + r + ", " + _nPat*OFF_SIZE );
						//throw 1;
					}
				}
			} catch(bad_alloc& e) {
				System.err.println( "Out of memory allocating plen[] in Ebwt.read()"
						+ " at " + __FILE__ + ":" + __LINE__ );
				//throw e;
			}
		}
		
		Boolean shmemLeader;
		
		// TODO: I'm not consistent on what "header" means.  Here I'm using
		// "header" to mean everything that would exist in memory if we
		// started to build the Ebwt but stopped short of the build*() step
		// (i.e. everything up to and including join()).
		if(justHeader) goto done;
		
		this._nFrag = readU<long>(_in1, switchEndian);
		bytesRead += OFF_SIZE;
		if(_verbose || startVerbose) {
			System.err.println( "Reading rstarts (" + this._nFrag*3 + "): ");
			//logTime(cerr);
		}
		assert_geq(this._nFrag, this._nPat);
		_rstarts.reset();
		if(loadRstarts) {
			if(_useMm) {
		#ifdef BOWTIE_MM
				_rstarts.init((long*)(mmFile[0] + bytesRead), _nFrag*3, false);
				bytesRead += this._nFrag*OFF_SIZE*3;
				fseeko(_in1, this._nFrag*OFF_SIZE*3, SEEK_CUR);
		#endif
			} else {
				_rstarts.init(new long[_nFrag*3], _nFrag*3, true);
				if(switchEndian) {
					for(long i = 0; i < this._nFrag*3; i += 3) {
						// fragment starting position in joined reference
						// string, text id, and fragment offset within text
						this.rstarts()[i]   = readU<long>(_in1, switchEndian);
						this.rstarts()[i+1] = readU<long>(_in1, switchEndian);
						this.rstarts()[i+2] = readU<long>(_in1, switchEndian);
					}
				} else {
					int r = MM_READ(_in1, (void *)rstarts(), this._nFrag*OFF_SIZE*3);
					if(r != (int)(this._nFrag*OFF_SIZE*3)) {
						System.err.println( "Error reading _rstarts[] array: " + r + ", " + (this._nFrag*OFF_SIZE*3) );
						//throw 1;
					}
				}
			}
		} else {
			// Skip em
			//assert(rstarts() == null);
			bytesRead += this._nFrag*OFF_SIZE*3;
			fseeko(_in1, this._nFrag*OFF_SIZE*3, SEEK_CUR);
		}
		
		_ebwt.reset();
		if(_useMm) {
		#ifdef BOWTIE_MM
			_ebwt.init((uint8_t*)(mmFile[0] + bytesRead), eh._ebwtTotLen, false);
			bytesRead += eh._ebwtTotLen;
			fseek(_in1, eh._ebwtTotLen, SEEK_CUR);
		#endif
		} else {
			// Allocate ebwt (big allocation)
			if(_verbose || startVerbose) {
				System.err.println( "Reading ebwt (" + eh._ebwtTotLen + "): ");
				//logTime(cerr);
			}
			Boolean shmemLeader = true;
			if(useShmem_) {
				uint8_t *tmp = null;
				shmemLeader = ALLOC_SHARED_U8(
						(_in1Str + "[ebwt]"), eh._ebwtTotLen, &tmp,
						"ebwt[]", (_verbose || startVerbose));
				assert(tmp != null);
				_ebwt.init(tmp, eh._ebwtTotLen, false);
				if(_verbose || startVerbose) {
					System.err.println( "  shared-mem " + (shmemLeader ? "leader" : "follower") );
				}
			} else {
				try {
					_ebwt.init(new uint8_t[eh._ebwtTotLen], eh._ebwtTotLen, true);
				} catch(bad_alloc& e) {
					System.err.println( "Out of memory allocating the ebwt[] array for the Bowtie index.  Please try"
							+ "again on a computer with more memory." );
					//throw 1;
				}
			}
			if(shmemLeader) {
				// Read ebwt from primary stream
				long bytesLeft = eh._ebwtTotLen;
					char *pebwt = (char*)this.ebwt();
				
				while (bytesLeft>0){
					int r = MM_READ(_in1, (void *)pebwt, bytesLeft);
					if(MM_IS_IO_ERR(_in1,r,bytesLeft)) {
						System.err.println( "Error reading _ebwt[] array: " + r + ", "
								+ bytesLeft + gLastIOErrMsg );
						//throw 1;
					}
					pebwt += r;
					bytesLeft -= r;
				}
				if(switchEndian) {
					uint8_t *side = this.ebwt();
					for(int i = 0; i < eh._numSides; i++) {
							long *cums = reinterpret_cast<long*>(side + eh._sideSz - OFF_SIZE*2);
						cums[0] = endianSwapU(cums[0]);
						cums[1] = endianSwapU(cums[1]);
						side += this._eh._sideSz;
					}
				}
		#ifdef BOWTIE_SHARED_MEM
				if(useShmem_) NOTIFY_SHARED(ebwt(), eh._ebwtTotLen);
		#endif
			} else {
				// Seek past the data and wait until master is finished
				fseeko(_in1, eh._ebwtTotLen, SEEK_CUR);
		#ifdef BOWTIE_SHARED_MEM
				if(useShmem_) WAIT_SHARED(ebwt(), eh._ebwtTotLen);
		#endif
			}
		}
		
		// Read zOff from primary stream
		_zOff = readU<long>(_in1, switchEndian);
		bytesRead += OFF_SIZE;
		assert_lt(_zOff, len);
		
		try {
			// Read fchr from primary stream
			if(_verbose || startVerbose) System.err.println( "Reading fchr (5)" );
			_fchr.reset();
			if(_useMm) {
		#ifdef BOWTIE_MM
				_fchr.init((long*)(mmFile[0] + bytesRead), 5, false);
				bytesRead += 5*OFF_SIZE;
				fseek(_in1, 5*OFF_SIZE, SEEK_CUR);
		#endif
			} else {
				_fchr.init(new long[5], 5, true);
				for(int i = 0; i < 5; i++) {
					this.fchr()[i] = readU<long>(_in1, switchEndian);
					assert_leq(this.fchr()[i], len);
					//assert(i <= 0 || this.fchr()[i] >= this.fchr()[i-1]);
				}
			}
			assert_gt(this.fchr()[4], this.fchr()[0]);
			// Read ftab from primary stream
			if(_verbose || startVerbose) {
				if(loadFtab) {
					System.err.println( "Reading ftab (" + eh._ftabLen + "): ");
					//logTime(cerr);
				} else {
					System.err.println( "Skipping ftab (" + eh._ftabLen + "): ");
				}
			}
			_ftab.reset();
			if(loadFtab) {
				if(_useMm) {
		#ifdef BOWTIE_MM
					_ftab.init((long*)(mmFile[0] + bytesRead), eh._ftabLen, false);
					bytesRead += eh._ftabLen*OFF_SIZE;
					fseeko(_in1, eh._ftabLen*OFF_SIZE, SEEK_CUR);
		#endif
				} else {
					_ftab.init(new long[eh._ftabLen], eh._ftabLen, true);
					if(switchEndian) {
						for(long i = 0; i < eh._ftabLen; i++)
							this.ftab()[i] = readU<long>(_in1, switchEndian);
					} else {
						int r = MM_READ(_in1, (void *)ftab(), eh._ftabLen*OFF_SIZE);
						if(r != (int)(eh._ftabLen*OFF_SIZE)) {
							System.err.println( "Error reading _ftab[] array: " + r + ", " + (eh._ftabLen*OFF_SIZE) );
							//throw 1;
						}
					}
				}
				// Read etab from primary stream
				if(_verbose || startVerbose) {
					if(loadFtab) {
						System.err.println( "Reading eftab (" + eh._eftabLen + "): ");
						//logTime(cerr);
					} else {
						System.err.println( "Skipping eftab (" + eh._eftabLen + "): ");
					}
					
				}
				_eftab.reset();
				if(_useMm) {
		#ifdef BOWTIE_MM
					_eftab.init((long*)(mmFile[0] + bytesRead), eh._eftabLen, false);
					bytesRead += eh._eftabLen*OFF_SIZE;
					fseeko(_in1, eh._eftabLen*OFF_SIZE, SEEK_CUR);
		#endif
				} else {
					_eftab.init(new long[eh._eftabLen], eh._eftabLen, true);
					if(switchEndian) {
						for(long i = 0; i < eh._eftabLen; i++)
							this.eftab()[i] = readU<long>(_in1, switchEndian);
					} else {
						int r = MM_READ(_in1, (void *)this.eftab(), eh._eftabLen*OFF_SIZE);
						if(r != (int)(eh._eftabLen*OFF_SIZE)) {
							System.err.println( "Error reading _eftab[] array: " + r + ", " + (eh._eftabLen*OFF_SIZE) );
							//throw 1;
						}
					}
				}
				for(long i = 0; i < eh._eftabLen; i++) {
					if(i > 0 && this.eftab()[i] > 0) {
						//assert_geq(this.eftab()[i], this.eftab()[i-1]);
					} else if(i > 0 && this.eftab()[i-1] == 0) {
						//assert_eq(0, this.eftab()[i]);
					}
				}
			} else {
				//assert(ftab() == null);
				//assert(eftab() == null);
				// Skip ftab
				bytesRead += eh._ftabLen*OFF_SIZE;
				fseeko(_in1, eh._ftabLen*OFF_SIZE, SEEK_CUR);
				// Skip eftab
				bytesRead += eh._eftabLen*OFF_SIZE;
				fseeko(_in1, eh._eftabLen*OFF_SIZE, SEEK_CUR);
			}
		} catch(bad_alloc& e) {
			System.err.println( "Out of memory allocating fchr[], ftab[] or eftab[] arrays for the Bowtie index. "
					+ "Please try again on a computer with more memory.");
			//throw 1;
		}
		
		// Read reference sequence names from primary index file (or not,
		// if --refidx is specified)
		if(loadNames) {
			while(true) {
				char c = '\0';
				if(MM_READ(_in1, (void *)(&c), (int)1) != (int)1) break;
				bytesRead++;
				if(c == '\0') break;
				else if(c == '\n') {
					this._refnames.push_back("");
				} else {
					if(this._refnames.size() == 0) {
						this._refnames.push_back("");
					}
					this._refnames.back().push_back(c);
				}
			}
		}
		
		_offs.reset();
		if(loadSASamp) {
			bytesRead = 4; // reset for secondary index file (already read 1-sentinel)
			
			shmemLeader = true;
			if(_verbose || startVerbose) {
				System.err.println( "Reading offs (" + offsLenSampled + std.setw(2) + OFF_SIZE*8 + "-bit words): ");
				//logTime(cerr);
			}
			
			if(!_useMm) {
				if(!useShmem_) {
					// Allocate offs_
					try {
						_offs.init(new long[offsLenSampled], offsLenSampled, true);
					} catch(bad_alloc& e) {
						System.err.println( "Out of memory allocating the offs[] array  for the Bowtie index." + "\n"
								+ "Please try again on a computer with more memory." );
						//throw 1;
					}
				} else {
						long *tmp = null;
					shmemLeader = ALLOC_SHARED_U(
							(_in2Str + "[offs]"), offsLenSampled*OFF_SIZE, &tmp,
							"offs", (_verbose || startVerbose));
					_offs.init((long*)tmp, offsLenSampled, false);
				}
			}
			
			if(_overrideOffRate < 32) {
				if(shmemLeader) {
					// Allocate offs (big allocation)
					if(switchEndian || offRateDiff > 0) {
						//assert(!_useMm);
							long blockMaxSz = (2 * 1024 * 1024); // 2 MB block size
							long blockMaxSzU = (blockMaxSz >> (OFF_SIZE/4 + 1)); // # U32s per block
							char *buf;
						try {
							buf = new char[blockMaxSz];
						} catch(std.bad_alloc& e) {
							System.err.println( "Error: Out of memory allocating part of _offs array: '" + e.what() + "'" );
							//throw e;
						}
						for(long i = 0; i < offsLen; i += blockMaxSzU) {
							long block = min<long>(blockMaxSzU, offsLen - i);
							int r = MM_READ(_in2, (void *)buf, block << (OFF_SIZE/4 + 1));
							if(r != (int)(block << (OFF_SIZE/4 + 1))) {
								System.err.println( "Error reading block of _offs[] array: " + r + ", " + (block + (OFF_SIZE/4 + 1)) );
								//throw 1;
							}
							long idx = i >> offRateDiff;
							for(long j = 0; j < block; j += (1 << offRateDiff)) {
								//assert_lt(idx, offsLenSampled);
								this.offs()[idx] = ((long*)buf)[j];
								if(switchEndian) {
									this.offs()[idx] = endianSwapU(this.offs()[idx]);
								}
								idx++;
							}
						}
						delete[] buf;
					} else {
						if(_useMm) {
		#ifdef BOWTIE_MM
							_offs.init((long*)(mmFile[1] + bytesRead), offsLen, false);
							bytesRead += offsSz;
							fseeko(_in2, offsSz, SEEK_CUR);
		#endif
						} else {
							// Workaround for small-index mode where MM_READ may
							// not be able to handle read amounts greater than 2^32
							// bytes.
							long bytesLeft = offsSz;
								char *offs = (char *)this.offs();
							
							while(bytesLeft > 0) {
								int r = MM_READ(_in2, (void*)offs, bytesLeft);
								if(MM_IS_IO_ERR(_in2,r,bytesLeft)) {
									System.err.println( "Error reading block of _offs[] array: "
											+ r + ", " + bytesLeft + gLastIOErrMsg );
									//throw 1;
								}
								offs += r;
								bytesLeft -= r;
							}
						}
					}
		#ifdef BOWTIE_SHARED_MEM
					if(useShmem_) NOTIFY_SHARED(offs(), offsLenSampled*OFF_SIZE);
		#endif
				} else {
					// Not the shmem leader
					fseeko(_in2, offsLenSampled*OFF_SIZE, SEEK_CUR);
		#ifdef BOWTIE_SHARED_MEM
					if(useShmem_) WAIT_SHARED(offs(), offsLenSampled*OFF_SIZE);
		#endif
				}
			}
		}
		
		this.postReadInit(*eh); // Initialize fields of Ebwt not read from file
		if(_verbose || startVerbose) print(cerr, *eh);
		
		// The fact that _ebwt and friends actually point to something
		// (other than null) now signals to other member functions that the
		// Ebwt is loaded into memory.
		
		done: // Exit hatch for both justHeader and !justHeader
		
		// Be kind
		if(deleteEh) delete eh;
		if(_in1 != null) {
			rewind(_in1);
		}
		if(_in2 != null) {
			rewind(_in2);
		}
	}
	
	public void sanityCheckUpToSide(long upToSide) {
		long occ[] = {0, 0, 0, 0};
		long cur = 0; // byte pointer
		EbwtParams eh = _eh;
		boolean fw = false;
		while(cur < (long)(upToSide * eh._sideSz)) {
			for(int i = 0; i < eh._sideBwtSz; i++) {
				byte by = ebwt()[cur + (fw ? i : eh._sideBwtSz-i-1)];
				for(int j = 0; j < 4; j++) {
					// Unpack from lowest to highest bit pair
					int twoBit = unpack_2b_from_8b(by, fw ? j : 3-j);
					occ[twoBit]++;
				}
			}
			
			occ_save[0] = occ[0];
			occ_save[1] = occ[1];
			occ_save[2] = occ[2];
			occ_save[3] = occ[3];
			cur += eh._sideSz;
		}
	}
	
	public void sanityCheckAll(int reverse) {
		EbwtParams eh = this._eh;
		assert(isInMemory());
		// Check ftab
		for(long i = 1; i < eh._ftabLen; i++) {
			//assert_geq(this.ftabHi(i), this.ftabLo(i-1));
			//assert_geq(this.ftabLo(i), this.ftabHi(i-1));
			//assert_leq(this.ftabHi(i), eh._bwtLen+1);
		}
		//assert_eq(this.ftabHi(eh._ftabLen-1), eh._bwtLen);
		
		// Check offs
		TIndexOff seenLen = (eh._bwtLen + 31) >> ((long)5);
		TIndexOff *seen;
		try {
			seen = new TIndexOff[seenLen]; // bitvector marking seen offsets
		} catch(Exception e) {
			System.err.println( "Out of memory allocating seen[] at " + __FILE__ + ":" + __LINE__ );
			throw e;
		}
		memset(seen, 0, OFF_SIZE * seenLen);
		long offsLen = eh._offsLen;
		for(long i = 0; i < offsLen; i++) {
			//assert_lt(this.offs()[i], eh._bwtLen);
			TIndexOff w = this.offs()[i] >> 5;
			TIndexOff r = this.offs()[i] & 31;
			//assert_eq(0, (seen[w] >> r) & 1); // shouldn't have been seen before
			seen[w] |= (1 << r);
		}
		delete[] seen;
		
		// Check nPat
		//assert_gt(this._nPat, 0);
		
		// Check plen, flen
		for(long i = 0; i < this._nPat; i++) {
			//assert_geq(this.plen()[i], 0);
		}
		
		// Check rstarts
		if(this.rstarts() != null) {
			for(long i = 0; i < this._nFrag-1; i++) {
				//assert_gt(this.rstarts()[(i+1)*3], this.rstarts()[i*3]);
				if(reverse == REF_READ_REVERSE) {
					//assert(this.rstarts()[(i*3)+1] >= this.rstarts()[((i+1)*3)+1]);
				} else {
					//assert(this.rstarts()[(i*3)+1] <= this.rstarts()[((i+1)*3)+1]);
				}
			}
		}
		
		// Check ebwt
		sanityCheckUpToSide(eh._numSides);
		//VMSG_NL("Ebwt.sanityCheck passed");
	}
	
	public void checkOrigs(EList<String> os, boolean color, boolean mirror) {
		SString<char> rest;
		restore(rest);
		long restOff = 0;
		int i = 0, j = 0;
		if(mirror) {
			// TODO: FIXME
			return;
		}
		while(i < os.size()) {
			int olen = os[i].length();
			int lastorig = -1;
			for(; j < olen; j++) {
				int joff = j;
				if(mirror) joff = olen - j - 1;
				if((int)os[i][joff] == 4) {
					// Skip over Ns
					lastorig = -1;
					if(!mirror) {
						while(j < olen && (int)os[i][j] == 4) j++;
					} else {
						while(j < olen && (int)os[i][olen-j-1] == 4) j++;
					}
					j--;
					continue;
				}
				if(lastorig == -1 && color) {
					lastorig = os[i][joff];
					continue;
				}
				if(color) {
					//assert_neq(-1, lastorig);
					//assert_eq(dinuc2color[(int)os[i][joff]][lastorig], rest[restOff]);
				} else {
					//assert_eq(os[i][joff], rest[restOff]);
				}
				lastorig = (int)os[i][joff];
				restOff++;
			}
			if(j == os[i].length()) {
				// Moved to next sequence
				i++;
				j = 0;
			} else {
				// Just jumped over a gap
			}
		}
	}
	
	public void postReadInit(EbwtParams eh) {
		long sideNum     = _zOff / eh._sideBwtLen;
		long sideCharOff = _zOff % eh._sideBwtLen;
		long sideByteOff = sideNum * eh._sideSz;
		_zEbwtByteOff = sideCharOff >> 2;
		_zEbwtBpOff = sideCharOff & 3;
		_zEbwtByteOff += sideByteOff;
	}
	
	public void print(OutputStream out, EbwtParams eh) throws IOException {
		eh.print(out); // print params
		out.write(("Ebwt (" + (isInMemory()? "memory" : "disk") + "):" + "\n"
				+ "    zOff: "         + _zOff + "\n"
				+ "    zEbwtByteOff: " + _zEbwtByteOff + "\n"
				+ "    zEbwtBpOff: "   + _zEbwtBpOff + "\n"
				+ "    nPat: "  + _nPat + "\n"
				+ "    plen: ").getBytes());
		if(plen() == null) {
			out.write(( "null" + "\n").getBytes());
		} else {
			out.write(( "non-null, [0] = " + plen()[0] + "\n").getBytes());
		}
		out.write( ("    rstarts: ").getBytes());
		if(rstarts() == null) {
			out.write( ("null" + "\n").getBytes());
		} else {
			out.write(( "non-null, [0] = " + rstarts()[0] + "\n").getBytes());
		}
		out.write( ("    ebwt: ").getBytes());
		if(ebwt() == null) {
			out.write( "null" + "\n");
		} else {
			out.write( "non-null, [0] = " + ebwt()[0] + "\n");
		}
		out.write( "    fchr: ");
		if(fchr() == null) {
			out.write( "null" + "\n");
		} else {
			out.write( "non-null, [0] = " + fchr()[0] + "\n");
		}
		out.write( "    ftab: ");
		if(ftab() == null) {
			out.write( "null" + "\n");
		} else {
			out.write( "non-null, [0] = " + ftab()[0] + "\n");
		}
		out.write( "    eftab: ");
		if(eftab() == null) {
			out.write( "null" + "\n");
		} else {
			out.write( "non-null, [0] = " + eftab()[0] + "\n");
		}
		out.write( "    offs: ");
		if(offs() == null) {
			out.write( "null" + "\n");
		} else {
			out.write( "non-null, [0] = " + offs()[0] + "\n");
		}
	}
	
	public long countBt2Side(SideLocus l, int c) {
		byte[] side = l.side(this.ebwt());
		long cCnt = countUpTo(l, c);
		if(c == 0 && l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
			// Adjust for the fact that we represented $ with an 'A', but
			// shouldn't count it as an 'A' here
			if((l._sideByteOff + l._by > _zEbwtByteOff) ||
					(l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
			{
				cCnt--; // Adjust for '$' looking like an 'A'
			}
		}
		long ret;
		// Now factor in the occ[] count at the side break
		byte[] acgt8 = side + _eh._sideBwtSz;
		long[] acgt = reinterpret_cast<long>(acgt8);
		ret = acgt[c] + cCnt + this.fchr()[c];
		return ret;
	}
	
	public void countBt2SideRange(SideLocus l, long num, long[] cntsUpto, long[] cntsIn, EList<boolean> masks) {
		countUpToEx(l, cntsUpto);
		WITHIN_FCHR_DOLLARA(cntsUpto);
		WITHIN_BWT_LEN(cntsUpto);
		byte[] side = l.side(this.ebwt());
		if(l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
			// Adjust for the fact that we represented $ with an 'A', but
			// shouldn't count it as an 'A' here
			if((l._sideByteOff + l._by > _zEbwtByteOff) ||
					(l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
			{
				cntsUpto[0]--; // Adjust for '$' looking like an 'A'
			}
		}
		// Now factor in the occ[] count at the side break
		long acgt = reinterpret_cast<long>(side + _eh._sideBwtSz);
		cntsUpto[0] += (acgt[0] + this.fchr()[0]);
		cntsUpto[1] += (acgt[1] + this.fchr()[1]);
		cntsUpto[2] += (acgt[2] + this.fchr()[2]);
		cntsUpto[3] += (acgt[3] + this.fchr()[3]);
		masks[0].resize(num);
		masks[1].resize(num);
		masks[2].resize(num);
		masks[3].resize(num);
		WITHIN_FCHR_DOLLARA(cntsUpto);
		WITHIN_FCHR_DOLLARA(cntsIn);
		// 'cntsUpto' is complete now.
		// Walk forward until we've tallied the entire 'In' range
		long nm = 0;
		// Rest of this side
		nm += countBt2SideRange2(l, true, num - nm, cntsIn, masks, nm);
		SideLocus lcopy = l;
		while(nm < num) {
			// Subsequent sides, if necessary
			lcopy.nextSide(this._eh);
			nm += countBt2SideRange2(lcopy, false, num - nm, cntsIn, masks, nm);
			WITHIN_FCHR_DOLLARA(cntsIn);
		}
		WITHIN_FCHR_DOLLARA(cntsIn);
	}
	
	public long countBt2SideRange2(
			SideLocus l,
			boolean startAtLocus,
			long num,
			long[] arrs,
			EList<boolean> masks,
			long maskOff){
		long nm = 0; // number of nucleotides tallied so far
		int iby = 0;      // initial byte offset
		int ibp = 0;      // initial base-pair offset
		if(startAtLocus) {
			iby = l._by;
			ibp = l._bp;
		} else {
			// Start at beginning
		}
		int by = iby, bp = ibp;
		byte[] side = l.side(this.ebwt());
		while(nm < num) {
			int c = (side[by] >> (bp * 2)) & 3;
			masks[0][maskOff + nm] = masks[1][maskOff + nm] =
					masks[2][maskOff + nm] = masks[3][maskOff + nm] = false;
			// Note: we tally $ just like an A
			arrs[c]++; // tally it
			masks[c][maskOff + nm] = true; // not dead
			nm++;
			if(++bp == 4) {
				bp = 0;
				by++;
				if(by == (int)this._eh._sideBwtSz) {
					// Fell off the end of the side
					break;
				}
			}
		}
		WITHIN_FCHR_DOLLARA(arrs);
		return nm;
	}
	
	public int rowL(SideLocus l) {
		// Extract and return appropriate bit-pair
		return unpack_2b_from_8b(l.side(this.ebwt())[l._by], l._bp);
	}
	
	public int rowL(long i) {
		// Extract and return appropriate bit-pair
		SideLocus l;
		l.initFromRow(i, _eh, ebwt());
		return rowL(l);
	}
	
	
	
	public void countBt2SideEx(SideLocus l, long[] arrs) {
		countUpToEx(l, arrs);
		if(l._sideByteOff <= _zEbwtByteOff && l._sideByteOff + l._by >= _zEbwtByteOff) {
			// Adjust for the fact that we represented $ with an 'A', but
			// shouldn't count it as an 'A' here
			if((l._sideByteOff + l._by > _zEbwtByteOff) ||
					(l._sideByteOff + l._by == _zEbwtByteOff && l._bp > _zEbwtBpOff))
			{
				arrs[0]--; // Adjust for '$' looking like an 'A'
			}
		}
		WITHIN_FCHR(arrs);
		WITHIN_BWT_LEN(arrs);
		// Now factor in the occ[] count at the side break
		byte[] side = l.side(this.ebwt());
		byte[] acgt16 = side + this._eh._sideSz - OFF_SIZE*4;
		long[] acgt = reinterpret_cast<long>(acgt16);
		arrs[0] += (acgt[0] + this.fchr()[0]);
		arrs[1] += (acgt[1] + this.fchr()[1]);
		arrs[2] += (acgt[2] + this.fchr()[2]);
		arrs[3] += (acgt[3] + this.fchr()[3]);
		WITHIN_FCHR(arrs);
	}
	
	public long countUpTo(SideLocus l, int c) {
		// Count occurrences of c in each 64-bit (using bit trickery);
		// Someday countInU64() and pop() functions should be
		// vectorized/SSE-ized in case that helps.
		long cCnt = 0;
		byte[] side = l.side(this.ebwt());
		int i = 0;
		for(; i + 7 < l._by; i += 8) {
			cCnt += countInU64(c, (long)side[i]);
		}
		// Count occurences of c in the rest of the side (using LUT)
		for(; i < l._by; i++) {
			cCnt += cCntLUT_4[0][c][side[i]];
		}
		// Count occurences of c in the rest of the byte
		if(l._bp > 0) {
			cCnt += cCntLUT_4[(int)l._bp][c][side[i]];
		}
		return cCnt;
	}
	
	public void countUpToEx(SideLocus l, long[] arrs) {
		int i = 0;
		// Count occurrences of each nucleotide in each 64-bit word using
		// bit trickery; note: this seems does not seem to lend a
		// significant boost to performance in practice.  If you comment
		// out this whole loop (which won't affect correctness - it will
		// just cause the following loop to take up the slack) then runtime
		// does not change noticeably. Someday the countInU64() and pop()
		// functions should be vectorized/SSE-ized in case that helps.
		byte[] side = l.side(this.ebwt());
		
		for(; i+7 < l._by; i += 8) {
			countInU64Ex((long)side[i], arrs);
		}
		// Count occurences of nucleotides in the rest of the side (using LUT)
		// Many cache misses on following lines (~20K)
		for(; i < l._by; i++) {
			arrs[0] += cCntLUT_4[0][0][side[i]];
			arrs[1] += cCntLUT_4[0][1][side[i]];
			arrs[2] += cCntLUT_4[0][2][side[i]];
			arrs[3] += cCntLUT_4[0][3][side[i]];
		}
		// Count occurences of c in the rest of the byte
		if(l._bp > 0) {
			arrs[0] += cCntLUT_4[(int)l._bp][0][side[i]];
			arrs[1] += cCntLUT_4[(int)l._bp][1][side[i]];
			arrs[2] += cCntLUT_4[(int)l._bp][2][side[i]];
			arrs[3] += cCntLUT_4[(int)l._bp][3][side[i]];
		}
	}
	
	public boolean isInMemory() {
		if(ebwt() != null) {
			return true;
		}
		return false;
	}
	
	public boolean isEvicted() {
		return !isInMemory();
	}
	
	public void loadIntoMemory(
			int color,
			int needEntireReverse,
			boolean loadSASamp,
			boolean loadFtab,
			boolean loadRstarts,
			boolean loadNames,
			boolean verbose){
		readIntoMemory(
				color,       // expect index to be colorspace?
				needEntireReverse, // require reverse index to be concatenated reference reversed
				loadSASamp,  // load the SA sample portion?
				loadFtab,    // load the ftab (_ftab[] and _eftab[])?
				loadRstarts, // load the r-starts (_rstarts[])?
				false,       // stop after loading the header portion?
				null,        // params
				false,       // mmSweep
				loadNames,   // loadNames
				verbose);    // startVerbose
	}
	
	public void evictFromMemory() {
		_zEbwtByteOff = IndexTypes.OFF_MASK;
		_zEbwtBpOff = -1;
	}
	
	public static long fileSize(String name) {
		File f = new File(name);
		return f.length();
	}
	
	public static int pop32(long x) {
		// Lots of cache misses on following lines (>10K)
		x = x - ((x >> 1) & 0x55555555);
		x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
		x = (x + (x >> 4)) & 0x0F0F0F0F;
		x = x + (x >> 8);
		x = x + (x >> 16);
		x = x + (x >> 32);
		return (int)(x & 0x3F);
	}
	
	public static int countInU64(int c, long dw) {
		long c_table[] = {
				0xffffffff,
				0xaaaaaaaa,
				0x55555555,
				0x00000000
		};
		long c0 = c_table[c];
		long x0 = dw ^ c0;
		long x1 = (x0 >> 1);
		long x2 = x1 & (0x55555555);
		long x3 = x0 & x2;
		long tmp = pop32(x3);
		return (int) tmp;
	}
	
	public EbwtParams eh() {
		return _eh;
	}
	
	public long ftabSeqToInt(BTDnaString seq, int off, boolean rev) {
		int fc = _eh._ftabChars;
		int lo = off, hi = lo + fc;
		long ftabOff = 0;
		for(int i = 0; i < fc; i++) {
			boolean fwex = fw();
			if(rev) fwex = !fwex;
			// We add characters to the ftabOff in the order they would
			// have been consumed in a normal search.  For BWT, this
			// means right-to-left order; for BWT' it's left-to-right.
			int c = (fwex ? seq[lo + i] : seq[hi - i - 1]);
			if(c > 3) {
				return long.MAX_VALUE;
			}
			ftabOff += 2;
			ftabOff |= c;
		}
		return ftabOff;
	}
	
	public long ftabHi(BTDnaString seq, int off) {
		return ftabHi(ftabSeqToInt(seq, off, false));
	}
	
	public long ftabHi(long i) {
		return EBWT.ftabHi(
				ftab(),
				eftab(),
				_eh._len,
				_eh._ftabLen,
				_eh._eftabLen,
				i);
	}
	
	public static long ftabHi(
			long ftab[],
			long eftab[],
			long len,
			long ftabLen,
			long eftabLen,
			long i) {
		if(ftab[i] <= len) {
			return ftab[i];
		} else {
			long efIdx = ftab[i] ^ com.uwb.bt2j.indexer.util.IndexTypes.OFF_MASK;
			return eftab[efIdx*2+1];
		}
	}
	
	public static long ftabLo(
			long ftab[],
			long eftab[],
			long len,
			long ftabLen,
			long eftabLen,
			long i) {
		if(ftab[i] <= len) {
			return ftab[i];
		} else {
			long efIdx = ftab[i] ^ IndexTypes.OFF_MASK;
			return eftab[efIdx*2];
		}
	}
	
	public long ftabLo(long i) {
		return EBWT.ftabLo(
				ftab(),
				eftab(),
				_eh._len,
				_eh._ftabLen,
				_eh._eftabLen,
				i);
	}
	
	public long ftabLo(BTDnaString seq, int off) {
		return ftabLo(ftabSeqToInt(seq, off, false));
	}
	
	public boolean ftabLoHi(
			BTDnaString seq, // sequence to extract from
			int off,             // offset into seq to begin extracting
			boolean rev,               // reverse while extracting
			long top,
			long bot)
	{
		long fi = ftabSeqToInt(seq, off, rev);
		if(fi == long.MAX_VALUE) {
			return false;
		}
		top = ftabHi(fi);
		bot = ftabLo(fi+1);
		return true;
	}
	
	public static Pair<EBWT, EBWT> fromString(
			String str,
			boolean packed,
			int color,
			int reverse,
			boolean bigEndian,
			int lineRate,
			int offRate,
			int ftabChars,
			String file,
			boolean useBlockwise,
			long bmax,
			long bmaxSqrtMult,
			long bmaxDivN,
			int dcv,
			int seed,
			boolean verbose,
			boolean autoMem,
			boolean sanity){
		EList<String> strs = new EList(1);
		strs.push_back(str);
		return fromStrings<TStr>(
				strs,
		packed,
				color,
				reverse,
				bigEndian,
				lineRate,
				offRate,
				ftabChars,
				file,
				useBlockwise,
				bmax,
				bmaxSqrtMult,
				bmaxDivN,
				dcv,
				seed,
				verbose,
				autoMem,
				sanity);
	}
	public static Pair<EBWT, EBWT> fromStrings(
			String strs,
			boolean packed,
			int color,
			int reverse,
			boolean bigEndian,
			int lineRate,
			int offRate,
			int ftabChars,
			String file,
			boolean useBlockwise,
			long bmax,
			long bmaxSqrtMult,
			long bmaxDivN,
			int dcv,
			int seed,
			boolean verbose,
			boolean autoMem,
			boolean sanity){
		EList<FileBuf> is = new EList(1);
		RefReadInParams refparams = new RefReadInParams(color, REF_READ_FORWARD, false, false);
		String ss = "";
		for(long i = 0; i < strs.size(); i++) {
			ss += ">" + i + "\n" + strs[i] + "\n";
		}
		auto_ptr<FileBuf> fb = new FileBuf(ss.get());
		is.push_back(fb.get());
		// Vector for the ordered list of "records" comprising the input
		// sequences.  A record represents a stretch of unambiguous
		// characters in one of the input sequences.
		EList<RefRecord> szs = new EList(1);
		Pair<Long,Long> sztot;
		sztot = BitPairReference.szsFromFasta(is, file, bigEndian, refparams, szs, sanity);
		// Construct Ebwt from input strings and parameters
		Ebwt ebwtFw = new Ebwt(
				TStr,
				packed,
				refparams.color ? 1 : 0,
				-1,           // fw
				lineRate,
				offRate,      // suffix-array sampling rate
				ftabChars,    // number of chars in initial arrow-pair calc
				file,         // basename for .?.ebwt files
				true,         // fw?
				useBlockwise, // useBlockwise
				bmax,         // block size for blockwise SA builder
				bmaxSqrtMult, // block size as multiplier of sqrt(len)
				bmaxDivN,     // block size as divisor of len
				dcv,          // difference-cover period
				is,           // list of input streams
				szs,          // list of reference sizes
				sztot.first,  // total size of all unambiguous ref chars
				refparams,    // reference read-in parameters
				seed,         // pseudo-random number generator seed
				-1,           // override offRate
				verbose,      // be talkative
				autoMem,      // pass exceptions up to the toplevel so that we can adjust memory settings automatically
				sanity);      // verify results and internal consistency
		refparams.reverse = reverse;
		szs.clear();
		sztot = BitPairReference.szsFromFasta(is, file, bigEndian, refparams, szs, sanity);
		// Construct Ebwt from input strings and parameters
		Ebwt ebwtBw = new Ebwt(
				TStr(),
				packed,
				refparams.color ? 1 : 0,
				reverse == REF_READ_REVERSE,
				lineRate,
				offRate,      // suffix-array sampling rate
				ftabChars,    // number of chars in initial arrow-pair calc
				file + ".rev",// basename for .?.ebwt files
				false,        // fw?
				useBlockwise, // useBlockwise
				bmax,         // block size for blockwise SA builder
				bmaxSqrtMult, // block size as multiplier of sqrt(len)
				bmaxDivN,     // block size as divisor of len
				dcv,          // difference-cover period
				is,           // list of input streams
				szs,          // list of reference sizes
				sztot.first,  // total size of all unambiguous ref chars
				refparams,    // reference read-in parameters
				seed,         // pseudo-random number generator seed
				-1,           // override offRate
				verbose,      // be talkative
				autoMem,      // pass exceptions up to the toplevel so that we can adjust memory settings automatically
				sanity);      // verify results and internal consistency
		return new Pair(ebwtFw, ebwtBw);
	}
	
	public boolean isPacked1() {
		return packed_;
	}
	
	public void szsToDisk(EList<RefRecord> szs, OutputStream os, int reverse) {
		long seq = 0;
		long off = 0;
		long totlen = 0;
		for(int i = 0; i < szs.size(); i++) {
			if(szs[i].len == 0) continue;
			if(szs[i].first) off = 0;
			off += szs[i].off;
			if(szs[i].first && szs[i].len > 0) seq++;
			long seqm1 = seq-1;
			long fwoff = off;
			if(reverse == ReadDir.REF_READ_REVERSE) {
				// Invert pattern idxs
				seqm1 = _nPat - seqm1 - 1;
				// Invert pattern idxs
				fwoff = plen()[seqm1] - (off + szs[i].len);
			}
			writeU<long>(os, totlen, this.toBe()); // offset from beginning of joined string
			writeU<long>(os, seqm1,  this.toBe()); // sequence id
			writeU<long>(os, fwoff,  this.toBe()); // offset into sequence
			totlen += szs[i].len;
			off += szs[i].len;
		}
	}
	
	public void initFromVector(
			EList<FileBuf> is,
			EList<RefRecord> szs,
			long sztot,
			RefReadInParams refparams,
			FileOutputStream out1,
			FileOutputStream out2,
			String outfile,
			FileOutputStream saOut,
			FileOutputStream bwtOut,
			int nthreads,
			boolean useBlockwise,
			long bmax,
			long bmaxSqrtMult,
			long bmaxDivN,
			int dcv,
			int seed,
			boolean verbose){
		// Compose text strings into single String
		if(this.verbose()) {
			String tmp = "Calculating joined length";
			this.verbose(tmp);
		}
		TStr s; // holds the entire joined reference after call to joinToDisk
		long jlen;
		jlen = joinedLen(szs);
		////assert_geq(jlen, sztot);
		if(this.verbose()) {
			String tmp = "Writing header";
			this.verbose(tmp);
		}
		writeFromMemory(true, out1, out2);
		try {
			if(this.verbose()) {
				String tmp = "Reserving space for joined String";
				this.verbose(tmp);
			}
			s.resize(jlen);
			VMSG_NL("Joining reference sequences");
			if(refparams.reverse == REF_READ_REVERSE) {
				//Timer timer(cout, "  Time to join reference sequences: ", _verbose);
				joinToDisk(is, szs, sztot, refparams, s, out1, out2);
				//Timer timer(cout, "  Time to reverse reference sequence: ", _verbose);
				//EList<RefRecord> tmp(1);
				s.reverse();
				reverseRefRecords(szs, tmp, false, verbose);
				szsToDisk(tmp, out1, refparams.reverse);
			} else {
				//Timer timer(cout, "  Time to join reference sequences: ", _verbose);
				joinToDisk(is, szs, sztot, refparams, s, out1, out2);
				szsToDisk(szs, out1, refparams.reverse);
			}
			// Joined reference sequence now in 's'
		} catch(Exception e) {
			// If we throw an allocation exception in the try block,
			// that means that the joined version of the reference
			// String itself is too larger to fit in memory.  The only
			// alternatives are to tell the user to give us more memory
			// or to try again with a packed representation of the
			// reference (if we haven't tried that already).
			System.err.println("Could not allocate space for a joined String of " + jlen + " elements.");;
			if(!isPacked() && _passMemExc) {
				// Pass the exception up so that we can retry using a
				// packed String representation
				throw e;
			}
			// There's no point passing this exception on.  The fact
			// that we couldn't allocate the joined String means that
			// --bmax is irrelevant - the user should re-run with
			// ebwt-build-packed
			if(isPacked1()) {
				System.err.println("Please try running bowtie-build on a computer with more memory.");
			} else {
				System.err.println("Please try running bowtie-build in packed mode (-p/--packed) or in automatic" + "\n"
						+ "mode (-a/--auto), or try again on a computer with more memory.");
			}
			if(sizeof(null) == 4) {
				System.err.println("If this computer has more than 4 GB of memory, try using a 64-bit executable;" + "\n"
						+ "this executable is 32-bit." + "\n");
			}
			throw 1;
		}
		// Succesfully obtained joined reference String
		if(bmax != IndexTypes.OFF_MASK) {
			VMSG_NL("bmax according to bmax setting: " + bmax);
		}
		else if(bmaxSqrtMult != IndexTypes.OFF_MASK) {
			bmax *= bmaxSqrtMult;
			VMSG_NL("bmax according to bmaxSqrtMult setting: " + bmax);
		}
		else if(bmaxDivN != IndexTypes.OFF_MASK) {
			bmax = max<long>(jlen / bmaxDivN, 1);
			VMSG_NL("bmax according to bmaxDivN setting: " + bmax);
		}
		else {
			bmax = (long)sqrt(s.length());
			VMSG_NL("bmax defaulted to: " + bmax);
		}
		int iter = 0;
		boolean first = true;
		streampos out1pos = out1.tellp();
		streampos out2pos = out2.tellp();
		// Look for bmax/dcv parameters that work.
		while(true) {
			if(!first && bmax < 40 && _passMemExc) {
				System.err.println("Could not find approrpiate bmax/dcv settings for building this index.");
				if(!isPacked1()) {
					// Throw an exception exception so that we can
					// retry using a packed String representation
					//throw bad_alloc();
				} else {
					System.err.println("Already tried a packed String representation.");
				}
				System.err.println("Please try indexing this reference on a computer with more memory.");
				if(sizeof(void*) == 4) {
					System.err.println("If this computer has more than 4 GB of memory, try using a 64-bit executable;" + "\n"
							+ "this executable is 32-bit.");
				}
				throw 1;
			}
			if(!first) {
				out1.seekp(out1pos);
				out2.seekp(out2pos);
			}
			if(dcv > 4096) dcv = 4096;
			if((iter % 6) == 5 && dcv < 4096 && dcv != 0) {
				dcv += 1; // double difference-cover period
			} else {
				bmax -= (bmax >> 2); // reduce by 25%
			}
			//VMSG("Using parameters --bmax " + bmax);
			if(dcv == 0) {
				//VMSG_NL(" and *no difference cover*");
			} else {
				//VMSG_NL(" --dcv " + dcv);
			}
			iter++;
			try {
				{
					//VMSG_NL("  Doing ahead-of-time memory usage test");
					// Make a quick-and-dirty attempt to force a bad_alloc iff
					// we would have thrown one eventually as part of
					// constructing the DifferenceCoverSample
					dcv += 1;
					long sz = (long)DifferenceCoverSample<TStr>.simulateAllocs(s, dcv >> 1);
					if(nthreads > 1) sz *= (nthreads + 1);
					AutoArray<byte> tmp(sz, 1);
					dcv >>= 1;
					// Likewise with the KarkkainenBlockwiseSA
					sz = (long)KarkkainenBlockwiseSA<TStr>.simulateAllocs(s, bmax);
					AutoArray<byte> tmp2(sz, 1);
					// Now throw in the 'ftab' and 'isaSample' structures
					// that we'll eventually allocate in buildToDisk
					AutoArray<long> ftab(_eh._ftabLen * 2, 1);
					AutoArray<byte> side(_eh._sideSz, 1);
					// Grab another 20 MB out of caution
					AutoArray<int> extra(20*1024*1024, 1);
					// If we made it here without throwing bad_alloc, then we
					// passed the memory-usage stress test
					//VMSG("  Passed!  Constructing with these parameters: --bmax " + bmax + " --dcv " + dcv);
					if(isPacked()) {
						//VMSG(" --packed");
					}
					//VMSG_NL("");
				}
				//VMSG_NL("Constructing suffix-array element generator");
				KarkkainenBlockwiseSA<TStr> bsa(s, bmax, nthreads, dcv, seed, _sanity, _passMemExc, _verbose, outfile);
				assert(bsa.suffixItrIsReset());
				////assert_eq(bsa.size(), s.length()+1);
				//VMSG_NL("Converting suffix-array elements to index image");
				buildToDisk(bsa, s, out1, out2, saOut, bwtOut);
				out1.flush(); out2.flush();
				boolean failed = out1.fail() || out2.fail();
				if(saOut != null) {
					saOut.flush();
					failed = failed || saOut.fail();
				}
				if(bwtOut != null) {
					bwtOut.flush();
					failed = failed || bwtOut.fail();
				}
				if(failed) {
					System.err.println("An error occurred writing the index to disk.  Please check if the disk is full." + "\n");
					//throw 1;
				}
				break;
			} catch(Exception e) {
				if(_passMemExc) {
					//VMSG_NL("  Ran out of memory; automatically trying more memory-economical parameters.");
				} else {
					System.err.println("Out of memory while constructing suffix array.  Please try using a smaller" + "\n"
							+ "number of blocks by specifying a smaller --bmax or a larger --bmaxdivn");
					//throw 1;
				}
			}
			first = false;
		}
		//assert(repOk());
		// Now write reference sequence names on the end
		////assert_eq(this._refnames.size(), this._nPat);
		for(long i = 0; i < this._refnames.size(); i++) {
			out1 + this._refnames[i] + "\n");;
		}
		out1 + '\0';
		out1.flush(); out2.flush();
		if(out1.fail() || out2.fail()) {
			System.err.println("An error occurred writing the index to disk.  Please check if the disk is full." + "\n");;
			throw 1;
		}
		//VMSG_NL("Returning from initFromVector");
	}
	
	public long tryOffset(long elt) {
		if(elt == _zOff) return 0;
		if((elt & _eh._offMask) == elt) {
			long eltOff = elt >> _eh._offRate;
			long off = offs()[eltOff];
			return off;
		} else {
			// Try looking at zoff
			return com.uwb.bt2j.indexer.util.IndexTypes.OFF_MASK;
		}
	}
	
	public long tryOffset(long elt, boolean fw, long hitlen) {
		long off = tryOffset(elt);
		if(off != com.uwb.bt2j.indexer.util.IndexTypes.OFF_MASK && !fw) {
			off = _eh._len - off - 1;
			off -= (hitlen-1);
		}
		return off;
	}
	
	public boolean is_read_err(int fdesc, int ret, int count) {
		if (ret < 0) {
			gLastIOErrMsg = "ERRNO: " + errno + " ERR Msg:" + strerror(errno) + "\n";
			return true;
		}
		return false;
	}
	
	public boolean is_fread_err(File file_hd, int ret, int count) {
		if(file_hd.canRead()) {
			gLastIOErrMsg = "Error Reading File!";
			return true;
		}
		return false;
	}
	
	public long[] fchr() {
		return _fchr.get();
	}
	
	public byte ebwt() {
		return _ebwt.get();
	}
	
	public long mapLF(SideLocus l) {
		long ret;
		int c = rowL(l);
		ret = countBt2Side(l, c);
		return ret;
	}
	
	public long mapLF1(long row, SideLocus l, int c) {
		if(rowL(l) != c || row == _zOff) return com.uwb.bt2j.indexer.util.IndexTypes.OFF_MASK;
		long ret = countBt2Side(l, c);
		return ret;
	}
	
	private PrintStream log() {
		return System.out;
	}
	
	public boolean verbose() {
		return _verbose;
	}
	
	private void verbose(String s) {
		if(verbose()) {
			log().println(s);
			log().flush();
		}
	}
	
	public boolean readEbwtColor(String instr) {
		long flags = Ebwt.readFlags(instr);
		if(flags < 0 && (((-flags) & EBWT_COLOR) != 0)) {
			return true;
		} else {
			return false;
		}
	}
	
	public boolean readEntireReverse(String instr) {
		long flags = Ebwt.readFlags(instr);
		if(flags < 0 && (((-flags) & EBWT_ENTIRE_REV) != 0)) {
			return true;
		} else {
			return false;
		}
	}
	
	public void readEbwtRefNames(File fin, EList<String> refnames) {
		// Read endianness hints from both streams
		Boolean switchEndian = false;
		long one = readU<long>(fin, switchEndian); // 1st word of primary stream
		if(one != 1) {
			//assert_eq((1u<<24), one);
			switchEndian = true;
		}
		
		// Reads header entries one by one from primary stream
		long len          = readU<long>(fin, switchEndian);
		long  lineRate     = readI<long>(fin, switchEndian);
		/*int32_t  linesPerSide =*/ readI<long>(fin, switchEndian);
		long  offRate      = readI<long>(fin, switchEndian);
		long  ftabChars    = readI<long>(fin, switchEndian);
		// BTL: chunkRate is now deprecated
		long flags = readI<long>(fin, switchEndian);
		Boolean color = false;
		Boolean entireReverse = false;
		if(flags < 0) {
			color = (((-flags) & EBWT_COLOR) != 0);
			entireReverse = (((-flags) & EBWT_ENTIRE_REV) != 0);
		}
		
		// Create a new EbwtParams from the entries read from primary stream
		EbwtParams eh(len, lineRate, offRate, ftabChars, color, entireReverse);
		
		long nPat = readI<long>(fin, switchEndian); // nPat
		fseeko(fin, nPat*OFF_SIZE, SEEK_CUR);
		
		// Skip rstarts
		long nFrag = readU<long>(fin, switchEndian);
		fseeko(fin, nFrag*OFF_SIZE*3, SEEK_CUR);
		
		// Skip ebwt
		fseeko(fin, eh._ebwtTotLen, SEEK_CUR);
		
		// Skip zOff from primary stream
		readU<long>(fin, switchEndian);
		
		// Skip fchr
		fseeko(fin, 5 * OFF_SIZE, SEEK_CUR);
		
		// Skip ftab
		fseeko(fin, eh._ftabLen*OFF_SIZE, SEEK_CUR);
		
		// Skip eftab
		fseeko(fin, eh._eftabLen*OFF_SIZE, SEEK_CUR);
		
		// Read reference sequence names from primary index file
		while(true) {
			char c = '\0';
			int read_value = 0;
			read_value = fgetc(fin);
			if(read_value == EOF) break;
			c = (char)read_value;
			if(c == '\0') break;
			else if(c == '\n') {
				refnames.push_back("");
			} else {
				if(refnames.size() == 0) {
					refnames.push_back("");
				}
				refnames.back().push_back(c);
			}
		}
		if(refnames.back().empty()) {
			refnames.pop_back();
		}
		
		// Be kind
		fseeko(fin, 0, SEEK_SET);
	}
	
	public void readEbwtRefNames(String isntr, EList<String> refnames) {
		File fin = new File(instr + ".1." + gEbwt_ext);
		// Initialize our primary and secondary input-stream fields
		if(fin == null) {
			throw new Exception("Cannot open file " + instr);
		}
		readEbwtRefnames(fin, refnames);
	}
	
	public void restore(String s) {
		s.resize(_eh._len);
		long jumps = 0;
		long i = _eh._len; // should point to final SA elt (starting with '$')
		SideLocus l = new SideLocus(i, _eh, ebwt());
		while(i != _zOff) {
			//if(_verbose) cout << "restore: i: " << i << endl;
			// Not a marked row; go back a char in the original string
			long newi = mapLF(l);
			s[_eh._len - jumps - 1] = rowL(l);
			i = newi;
			l.initFromRow(i, _eh, ebwt());
			jumps++;
		}
	}
}
