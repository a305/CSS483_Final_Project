package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

public class VectorPatternSource extends PatternSource{
	private double cur_;
	private double skip_;
	private boolean paired_;
	private EList<String> tokbuf_;
	private String nametmp_;
	
	public VectorPatternSource(PatternParams p, EList<String> seqs) {
		super(p);
		cur_ = p.skip;
		skip_ = p.skip;
		paired_ = false;
		
		// Install sequences in buffers, ready for immediate copying in
		// nextBatch().  Formatting of the buffer is just like
		// TabbedPatternSource.
		int seqslen = seqs.size();
		for(int i = 0; i < seqslen; i++) {
			tokbuf_.clear();
			tokenize(seqs[i], ":", tokbuf_, 2);
			// Get another buffer ready
			bufs_.expand();
			bufs_.back().clear();
			// Install name
			nametmp_ = String.valueOf((long)i);
			bufs_.back().install(nametmp_);
			bufs_.back().append('\t');
			// Install sequence
			bufs_.back().append(tokbuf_[0].c_str());
			bufs_.back().append('\t');
			// Install qualities
			if(tokbuf_.size() > 1) {
				bufs_.back().append(tokbuf_[1].c_str());
			} else {
				double len = tokbuf_[0].length();
				for(int i = 0; i < len; i++) {
					bufs_.back().append('I');
				}
			}
		}
	}
	
	private Pair<Boolean, Integer> nextBatchImpl(PerThreadReadBuf pt, boolean batch_a) {
		pt.setReadId(cur_);
		EList<Read> readbuf = batch_a ? pt.bufa_ : pt.bufb_;
		double readi = 0;
		for(; readi < pt.max_buf_ && cur_ < bufs_.size(); readi++, cur_++) {
			readbuf[readi].readOrigBuf = bufs_[cur_];
		}
		readCnt_ += readi;
		return make_pair(cur_ == bufs_.size(), readi);
	}
	
	@Override
	public void nextBatch(PerThreadReadBuf pt, boolean batch_a, boolean lock) {
		
	}

	@Override
	public boolean parse(Read ra, Read rb, long rdid) {
		// Very similar to TabbedPatternSource

		// Light parser (nextBatchFromFile) puts unparsed data
		// into Read& r, even when the read is paired.
		int c = '\t';
		double cur = 0;
		double buflen = ra.readOrigBuf.length();
		
		// Loop over the two ends
		for(int endi = 0; endi < 2 && c == '\t'; endi++) {
			Read r = ((endi == 0) ? ra : rb);
			// Parse name if (a) this is the first end, or
			// (b) this is tab6
			if(endi < 1 || paired_) {
				// Parse read name
				c = ra.readOrigBuf[cur++];
				while(c != '\t' && cur < buflen) {
					r.name.append(c);
					c = ra.readOrigBuf[cur++];
				}
				if(cur >= buflen) {
					return false; // record ended prematurely
				}
			} else if(endi > 0) {
				// if this is the second end and we're parsing
				// tab5, copy name from first end
				rb.name = ra.name;
			}

			// Parse sequence
			c = ra.readOrigBuf[cur++];
			int nchar = 0;
			while(c != '\t' && cur < buflen) {
				if(isalpha(c)) {
					if(nchar++ >= pp_.trim5) {
						r.patFw.append(asc2dna[c]); // ascii to int
					}
				}
				c = ra.readOrigBuf[cur++];
			}
			if(cur >= buflen) {
				return false; // record ended prematurely
			}
			// record amt trimmed from 5' end due to --trim5
			r.trimmed5 = (int)(nchar - r.patFw.length());
			// record amt trimmed from 3' end due to --trim3
			r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
			
			// Parse qualities
			c = ra.readOrigBuf[cur++];
			int nqual = 0;
			while(c != '\t' && c != '\n' && c != '\r') {
				if(c == ' ') {
					wrongQualityFormat(r.name);
					return false;
				}
				char cadd = charToPhred33(c, false, false);
				if(++nqual > pp_.trim5) {
					r.qual.append(cadd);
				}
				if(cur >= buflen) break;
				c = ra.readOrigBuf[cur++];
			}
			if(nchar > nqual) {
				tooFewQualities(r.name);
				return false;
			} else if(nqual > nchar) {
				tooManyQualities(r.name);
				return false;
			}
			r.qual.trimEnd(pp_.trim3);
		}
		ra.parsed = true;
		if(!rb.parsed && !rb.readOrigBuf.empty()) {
			return parse(rb, ra, rdid);
		}
		return true;
	}

	@Override
	public void reset() {
		readCnt_ = 0;
		cur_ = skip_;
		paired_ = false;
	}
}
