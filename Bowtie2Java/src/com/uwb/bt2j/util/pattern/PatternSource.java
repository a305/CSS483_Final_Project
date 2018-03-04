package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;
import com.uwb.bt2j.indexer.Formats;
import com.uwb.bt2j.indexer.Formats.FileFormat;

public abstract class PatternSource {
	protected PatternParams pp_;
	protected long readCnt_;
	
	public PatternSource(PatternParams p) {
		pp_ = p;
		readCnt_ = 0;
	}
	
	public abstract Pair<boolean, Integer> nextBatch(PerThreadReadBuf pt, boolean batch_a, boolean lock) {
		
	}
	
	public abstract boolean parse(Read ra, Read rb, long rdid);
	
	public abstract void reset() {
		readCnt_ = 0;
	}
	
	public static PatternSource patsrcFromStrings(PatternParams p, EList<String> qs) {
		switch(p.format) {
		case FASTA:       return new FastaPatternSource(qs, p);
		case FASTA_CONT:  return new FastaContinuousPatternSource(qs, p);
		case RAW:         return new RawPatternSource(qs, p);
		case FASTQ:       return new FastqPatternSource(qs, p);
		case INTERLEAVED: return new FastqPatternSource(qs, p, true /* interleaved */);
		case TAB_MATE5:   return new TabbedPatternSource(qs, p, false);
		case TAB_MATE6:   return new TabbedPatternSource(qs, p, true);
		case CMDLINE:     return new VectorPatternSource(qs, p);
		case QSEQ:        return new QseqPatternSource(qs, p);
		default: {
			System.err.println("Internal error; bad patsrc format: " + p.format);
		}
	}
	}
	
	public long readCount() {
		return readCnt_;
	}
}
