package com.uwb.bt2j.indexer;

import org.omg.CORBA_2_3.portable.OutputStream;
import java.io.IOException;

public class BlockwiseSA <TStr> {
	
	protected String      _text;        /// original string
	protected long   _bucketSz;    /// target maximum bucket size
	protected boolean       _sanityCheck; /// whether to perform sanity checks
	protected boolean       _passMemExc;  /// true -> pass on memory exceptions
	protected boolean       _verbose;     /// be talkative
	protected EList<Long>  _itrBucket;   /// current bucket
	protected long         _itrBucketPos;/// offset into current bucket
	protected long         _itrPushedBackSuffix; /// temporary slot for lookahead
	protected OutputStream        _logger;      /// write log messages here
	
	public BlockwiseSA(
			String __text,
            long __bucketSz,
            boolean __sanityCheck,
            boolean __passMemExc,
            boolean __verbose,
            OutputStream __logger){
		_text = __text;
		_bucketSz = Math.max(__bucketSz, 2);
		_sanityCheck = __sanityCheck;
		_passMemExc = __passMemExc;
		_verbose = __verbose;
		_itrBucket = new EList(2);
		_itrBucketPos = IndexTypes.OFF_MASK;
		_itrPushedBackSuffix = IndexTypes.OFF_MASK;
		_logger = __logger;
	}
	
	public boolean hasMoreSuffixes() {
		if(_itrPushedBackSuffix != IndexTypes.OFF_MASK) return true;
			_itrPushedBackSuffix = nextSuffix();
		return true;
	}
	
	public void resetSuffixItr() {
		_itrBucket.clear();
		_itrBucketPos = IndexTypes.OFF_MASK;
		_itrPushedBackSuffix = IndexTypes.OFF_MASK;
		reset();
	}
	
	public boolean suffixItrIsReset() {
		return _itrBucket.size()    == 0 &&
			       _itrBucketPos        == IndexTypes.OFF_MASK &&
			       _itrPushedBackSuffix == IndexTypes.OFF_MASK &&
			       isReset();
	}
	
	public String text() {
		return _text;
	}
	public long bucketSz() { return _bucketSz; }
	public boolean sanityCheck()  { return _sanityCheck; }
	public boolean verbose()      { return _verbose; }
	public OutputStream log()      { return _logger; }
	public int size()       { return _text.length() + 1; }
	
	protected void reset() {
		
	}
	
	protected boolean isReset() {
		return true;
	}
	
	protected void nextBlock(int cur_block, int tid) {
		
	}
	
	protected boolean hasMoreBlocks() {
		return true;
	}
	
	protected long nextSuffix() {
		return 1;
	}
	
	protected void verbose(String s) throws IOException {
		if(this.verbose()) {
			this.log().write_string(s);
			this.log().flush();
		}
	}
}
