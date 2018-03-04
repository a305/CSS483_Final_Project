package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

import javafx.util.Pair;

public abstract class PatternComposer {
	public PatternComposer(PatternParams p) {
		
	}
	
	public abstract void reset();
	public abstract Pair<Boolean, Integer> nextBatch(PerThreadReadBuf pt);
	public abstract boolean parse(Read ra, Read rb, long rdid);
	public static PatternComposer setupPatternComposer(
			EList<String> si,
			EList<String> m1,
			EList<String> m2,
			EList<String> m12,
			EList<String> q,
			EList<String> q1,
			EList<String> q2,
			PatternParams p,
			boolean verbose) {
		EList<PatternSource> a  = new EList<PatternSource>();
		EList<PatternSource> b  = new EList<PatternSource>();
		EList<PatternSource> ab = new EList<PatternSource>();
		// Create list of pattern sources for paired reads appearing
		// interleaved in a single file
		for(int i = 0; i < m12.size(); i++) {
			EList<String>* qs = m12;
			EList<String> tmp;
			if(p.fileParallel) {
				// Feed query files one to each PatternSource
				qs = &tmp;
				tmp.push_back(m12[i]);
			}
			ab.push_back(PatternSource.patsrcFromStrings(p, qs));
			if(!p.fileParallel) {
				break;
			}
		}

		// Create list of pattern sources for paired reads
		for(int i = 0; i < m1.size(); i++) {
			EList<String> qs = m1;
			EList<String> tmpSeq;
			EList<String> tmpQual;
			if(p.fileParallel) {
				// Feed query files one to each PatternSource
				qs = tmpSeq;
				tmpSeq.push_back(m1[i]);
			}
			a.push_back(PatternSource.patsrcFromStrings(p, qs));
			if(!p.fileParallel) {
				break;
			}
		}

		// Create list of pattern sources for paired reads
		for(int i = 0; i < m2.size(); i++) {
			EList<String> qs = m2;
			EList<String> tmpSeq;
			EList<String> tmpQual;
			if(p.fileParallel) {
				// Feed query files one to each PatternSource
				qs = &tmpSeq;
				tmpSeq.push_back(m2[i]);
			}
			b.push_back(PatternSource.patsrcFromStrings(p, qs));
			if(!p.fileParallel) {
				break;
			}
		}
		// All mates/mate files must be paired

		// Create list of pattern sources for the unpaired reads
		for(int i = 0; i < si.size(); i++) {
			EList<String> qs = si;
			PatternSource patsrc = null;
			EList<String> tmpSeq;
			EList<String> tmpQual;
			if(p.fileParallel) {
				// Feed query files one to each PatternSource
				qs = tmpSeq;
				tmpSeq.push_back(si[i]);
			}
			patsrc = PatternSource.patsrcFromStrings(p, qs);
			a.push_back(patsrc);
			b.push_back(null);
			if(!p.fileParallel) {
				break;
			}
		}

		PatternComposer patsrc = null;
		if(m12.size() > 0) {
			patsrc = new SoloPatternComposer(p, ab);
		} else {
			patsrc = new DuelPatternComposer(p, a, b);
		}
		return patsrc;
	}
	
	public static void freeEListPMembers(EList<PatternSource> e) {
		
	}
}
