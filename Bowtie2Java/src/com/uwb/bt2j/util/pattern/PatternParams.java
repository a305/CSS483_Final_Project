package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.indexer.Formats;
import com.uwb.bt2j.indexer.Formats.FileFormat;

public class PatternParams {
	Formats.FileFormat format;			  // file format
	boolean fileParallel;	  // true -> wrap files with separate PatternComposers
	double seed;		  // pseudo-random seed
	PerThreadReadBuf max_buf;		  // number of reads to buffer in one read
	boolean solexa64;		  // true -> qualities are on solexa64 scale
	boolean phred64;		  // true -> qualities are on phred64 scale
	boolean intQuals;		  // true -> qualities are space-separated numbers
	int trim5;            // amt to hard clip from 5' end
	int trim3;            // amt to hard clip from 3' end
	int sampleLen;		  // length of sampled reads for FastaContinuous...
	int sampleFreq;		  // frequency of sampled reads for FastaContinuous...
	double skip;		  // skip the first 'skip' patterns
	int nthreads;		  // number of threads for locking
	boolean fixName;
	
	public PatternParams(
			FileFormat format_,
			boolean fileParallel_,
			double seed_,
			PerThreadReadBuf max_buf_,
			boolean solexa64_,
			boolean phred64_,
			boolean intQuals_,
			int trim5_,
			int trim3_,
			int sampleLen_,
			int sampleFreq_,
			double skip_,
			int nthreads_,
			boolean fixName_) {
		format = format_;
		fileParallel = fileParallel_;
		seed=seed_;
		max_buf=max_buf_;
		solexa64=solexa64_;
		phred64=phred64_;
		intQuals=intQuals_;
		trim5=trim5_;
		trim3=trim3_;
		sampleLen=sampleLen_;
		sampleFreq=sampleFreq_;
		skip=skip_;
		nthreads=nthreads_;
		fixName=fixName_;
	}
}
