package com.uwb.bt2j.aligner;

public class OuterLoopMetrics {
		long reads;   // total reads
		long bases;   // total bases
		long srreads; // same-read reads
		long srbases; // same-read bases
		long freads;  // filtered reads
		long fbases;  // filtered bases
		long ureads;  // unfiltered reads
		long ubases;  // unfiltered bases
		
		public OuterLoopMetrics() {
			reset();
		}
		
		public void reset() {
			reads = bases = srreads = srbases =
					freads = fbases = ureads = ubases = 0;
		}
		
		public void merge(OuterLoopMetrics m) {
			reads += m.reads;
			bases += m.bases;
			srreads += m.srreads;
			srbases += m.srbases;
			freads += m.freads;
			fbases += m.fbases;
			ureads += m.ureads;
			ubases += m.ubases;
		}
}
