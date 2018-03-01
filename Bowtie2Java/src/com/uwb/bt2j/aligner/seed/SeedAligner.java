package com.uwb.bt2j.aligner.seed;

import javafx.util.Pair;

public class SeedAligner {
	public enum SeedAlignerActions {
		SA_ACTION_TYPE_RESET(1),
				SA_ACTION_TYPE_SEARCH_SEED(2), // 2
				SA_ACTION_TYPE_FTAB(3),        // 3
				SA_ACTION_TYPE_FCHR(4),        // 4
				SA_ACTION_TYPE_MATCH(5),       // 5
				SA_ACTION_TYPE_EDIT(6);         // 6
				private int x;
		SeedAlignerActions(int y) {
			x = y;
		}
	}
	
	public void instantiateSeq(Read read, BTDnaString seq, BTString qual, int len, int depth, Boolean fw) {
		// Fill in 'seq' and 'qual'
		int seedlen = len;
		if((int)read.length() < seedlen) seedlen = (int)read.length();
		seq.resize(len);
		qual.resize(len);
		// If fw is false, we take characters starting at the 3' end of the
		// reverse complement of the read.
		for(int i = 0; i < len; i++) {
			seq.set(read.patFw.windowGetDna(i, fw, depth, len), i);
			qual.set(read.qual.windowGet(i, fw, depth, len), i);
		}
	}
	
	public Pair<Integer, Integer> instantiateSeeds(
			EList<Seed> seeds,
			double off,
			int per,
			Read read,
			Scoring pens,
			Boolean nofw,
			Boolean norc,
			AlignmentCacheIFace cache,
			SeedResults sr,
			SeedSearchMetrics met,
			Pair<Integer, Integer> instFw,
			Pair<Integer, Integer> instRc) {
		offIdx2off_.clear();
		int len = seeds[0].len; // assume they're all the same length
		// Calc # seeds within read interval
		int nseeds = 1;
		if((int)read.length() - (int)off > len) {
			nseeds += ((int)read.length() - (int)off - len) / per;
		}
		for(int i = 0; i < nseeds; i++) {
			offIdx2off_.push_back(per * i + (int)off);
		}
		Pair<Integer, Integer> ret;
		ret.first = 0;  // # seeds that require alignment
		ret.second = 0; // # seeds that hit in cache with non-empty results
		sr.reset(read, offIdx2off_, nseeds);
		// For each seed position
		for(int fwi = 0; fwi < 2; fwi++) {
			Boolean fw = (fwi == 0);
			if((fw && nofw) || (!fw && norc)) {
				// Skip this orientation b/c user specified --nofw or --norc
				continue;
			}
			// For each seed position
			for(int i = 0; i < nseeds; i++) {
				int depth = i * per + (int)off;
				int seedlen = seeds[0].len;
				// Extract the seed sequence at this offset
				// If fw == true, we extract the characters from i*per to
				// i*(per-1) (exclusive).  If fw == false, 
				instantiateSeq(
					read,
					sr.seqs(fw)[i],
					sr.quals(fw)[i],
					std::min<int>((int)seedlen, (int)read.length()),
					depth,
					fw);
				QKey qk(sr.seqs(fw)[i] ASSERT_ONLY(, tmpdnastr_));
				// For each search strategy
				EList<InstantiatedSeed>& iss = sr.instantiatedSeeds(fw, i);
				for(int j = 0; j < (int)seeds.size(); j++) {
					iss.expand();
					assert_eq(seedlen, seeds[j].len);
					InstantiatedSeed* is = &iss.back();
					if(seeds[j].instantiate(
						read,
						sr.seqs(fw)[i],
						sr.quals(fw)[i],
						pens,
						depth,
						i,
						j,
						fw,
						*is))
					{
						// Can we fill this seed hit in from the cache?
						ret.first++;
						if(fwi == 0) { instFw.first++; } else { instRc.first++; }
					} else {
						// Seed may fail to instantiate if there are Ns
						// that prevent it from matching
						met.filteredseed++;
						iss.pop_back();
					}
				}
			}
		}
		return ret;
	}
}
