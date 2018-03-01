package com.uwb.bt2j.aligner.seed;

public class SeedSearchMetrics {
	public long seedsearch;   // # times we executed strategy in InstantiatedSeed
	public long nrange;       // # ranges found
	public long nelt;         // # range elements found
	public long possearch;    // # offsets where aligner executed >= 1 strategy
	public long intrahit;     // # offsets where current-read cache gave answer
	public long interhit;     // # offsets where across-read cache gave answer
	public long filteredseed; // # seed instantiations skipped due to Ns
	public long ooms;         // out-of-memory errors
	public long bwops;        // Burrows-Wheeler operations
	public long bweds;        // Burrows-Wheeler edits
	public long bestmin0;     // # times the best min # edits was 0
	public long bestmin1;     // # times the best min # edits was 1
	public long bestmin2;     // # times the best min # edits was 2
	public SeedSearchMetrics() {
		reset();
	}
	
	public void merge(SeedSearchMetrics m, boolean getLock) {
		seedsearch   += m.seedsearch;
		nrange       += m.nrange;
		nelt         += m.nelt;
		possearch    += m.possearch;
		intrahit     += m.intrahit;
		interhit     += m.interhit;
		filteredseed += m.filteredseed;
		ooms         += m.ooms;
		bwops        += m.bwops;
		bweds        += m.bweds;
		bestmin0     += m.bestmin0;
		bestmin1     += m.bestmin1;
		bestmin2     += m.bestmin2;
	}
	
	public void reset() {
		seedsearch =
				nrange =
				nelt =
				possearch =
				intrahit =
				interhit =
				filteredseed =
				ooms =
				bwops =
				bweds =
				bestmin0 =
				bestmin1 =
				bestmin2 = 0;
	}
}
