package com.uwb.bt2j.aligner.sink;

public class ALNSinkSam {
	protected final SamConfig samc_;
	protected BTDnaString dseq_;
	protected BTString dqual_;
	
	public ALNSinkSam() {
		
	}
	
	public void append(
			BString o,
			StackedAln staln,
			double threadId,
			final Read rd1,
			final Read rd2,
			final long rdid,
			AlnRes rs1,
			AlnRes rs2,
			final AlnSetSumm summ,
			final SeedAlSumm ssm1,
			final SeedAlSumm ssm2,
			final AlnFlags flags1,
			final AlnFlags flags2,
			final PerReadMetrics prm,
			final Mapq mapq,
			final Scoring sc,
			Boolean report2) {
		if(rd1 != null) {
			appendMate(o, staln, rd1, rd2, rdid, rs1, rs2, summ, ssm1, ssm2,
			           flags1, prm, mapq, sc);
		}
		if(rd2 != null && report2) {
			appendMate(o, staln, rd2, rd1, rdid, rs2, rs1, summ, ssm2, ssm1,
			           flags2, prm, mapq, sc);
		}
	}
	
	protected void appendMate() {
		
	}
}
