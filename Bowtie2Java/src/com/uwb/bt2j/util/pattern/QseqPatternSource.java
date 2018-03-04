package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

import javafx.util.Pair;

public abstract class QseqPatternSource extends CFilePatternSource{
	protected EList<String> qualToks_;

	public QseqPatternSource(PatternParams p, EList<String> infiles) {
		super(p, infiles);
	}

	@Override
	protected abstract Pair<Boolean, Integer> nextBatchFromFile(PerThreadReadBuf pt, boolean batch_a, int read_idx);

	@Override
	protected void resetForNextFile() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void nextBatch(PerThreadReadBuf pt, boolean batch_a, boolean lock) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public abstract boolean parse(Read ra, Read rb, long rdid);

}
