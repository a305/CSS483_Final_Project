package com.uwb.bt2j.aligner;

public class SSEMetrics {
	public long dp;       // DPs tried
	public long dpsat;    // DPs saturated
	public long dpfail;   // DPs failed
	public long dpsucc;   // DPs succeeded
	public long col;      // DP columns
	public long cell;     // DP cells
	public long inner;    // DP inner loop iters
	public long fixup;    // DP fixup loop iters
	public long gathsol;  // DP gather solution cells found
	public long bt;       // DP backtraces
	public long btfail;   // DP backtraces failed
	public long btsucc;   // DP backtraces succeeded
	public long btcell;   // DP backtrace cells traversed
	public long corerej;  // DP backtrace core rejections
	public long nrej;     // DP backtrace N rejections
	
	public SSEMetrics() {
		reset();
	}
	
	public void clear() {
		reset();
	}
	
	public void reset() {
		dp = dpsat = dpfail = dpsucc = 
				col = cell = inner = fixup =
				gathsol = bt = btfail = btsucc = btcell =
				corerej = nrej = 0;
	}
	
	public void merge(SSEMetrics o) {
		dp       += o.dp;
		dpsat    += o.dpsat;
		dpfail   += o.dpfail;
		dpsucc   += o.dpsucc;
		col      += o.col;
		cell     += o.cell;
		inner    += o.inner;
		fixup    += o.fixup;
		gathsol  += o.gathsol;
		bt       += o.bt;
		btfail   += o.btfail;
		btsucc   += o.btsucc;
		btcell   += o.btcell;
		corerej  += o.corerej;
		nrej     += o.nrej;
	}
	
	
}
