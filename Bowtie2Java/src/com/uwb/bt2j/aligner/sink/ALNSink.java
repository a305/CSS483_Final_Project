package com.uwb.bt2j.aligner.sink;

import java.io.DataOutputStream;
import java.io.IOException;

public class ALNSink {
	protected OutputQueue oq_;
	protected final EList<String> refnames_;
	protected Boolean quiet_;
	protected ReportingMetrics met_;
	
	public ALNSink(){
		
	}
	
	public void append() {
		
	}
	
	public void reportHits(
			BString o,
			StackedAln staln,
			double threadId,
			final Read rd1,
			final Read rd2,
			final long rdid,
			final EList<double> select1,
			final EList<double> select2,
			EList<AlnRes> rs1,
			EList<AlnRes> rs2,
			Boolean maxed,
			final AlnSetSumm summ,
			final SeedAlSumm ssm1,
			final SeedAlSumm ssm2,
			final AlnFlags flags1,
			final AlnFlags flags2,
			final PerReadMetrics prm,
			final Mapq mapq,
			final Scoring sc,
			Boolean reportBoth,
			Boolean getLock) {
		AlnFlags flagscp1, flagscp2;
		if(flags1 != NULL) {
			flagscp1 = flags1;
			flags1 = flagscp1;
			flagscp1.setPrimary(true);
		}
		if(flags2 != NULL) {
			flagscp2 = flags2;
			flags2 = flagscp2;
			flagscp2.setPrimary(true);
		}
		if(select2 != NULL) {
			// Handle case 4
			AlnRes r1pri = ((rs1 != NULL) ? rs1.get(select1[0]) : NULL);
			AlnRes r2pri = ((rs2 != NULL) ? rs2.get((select2)[0]) : NULL);
			append(o, staln, threadId, rd1, rd2, rdid, r1pri, r2pri, summ,
			       ssm1, ssm2, flags1, flags2, prm, mapq, sc, false);
			flagscp1.setPrimary(false);
			flagscp2.setPrimary(false);
			for(double i = 1; i < select1.size(); i++) {
				AlnRes r1 = ((rs1 != NULL) ? rs1.get(select1[i]) : NULL);
				append(o, staln, threadId, rd1, rd2, rdid, r1, r2pri, summ,
				       ssm1, ssm2, flags1, flags2, prm, mapq, sc, false);
			}
			if(reportBoth) {
				for(double i = 1; i < select2->size(); i++) {
					AlnRes r2 = ((rs2 != NULL) ? rs2.get((select2)[i]) : NULL);
					append(o, staln, threadId, rd2, rd1, rdid, r2, r1pri, summ,
						   ssm2, ssm1, flags2, flags1, prm, mapq, sc, false);
				}
			}
		} else {
			// Handle cases 1-3 and 5
			for(double i = 0; i < select1.size(); i++) {
				AlnRes r1 = ((rs1 != NULL) ? rs1.get(select1[i]) : NULL);
				AlnRes r2 = ((rs2 != NULL) ? rs2.get(select1[i]) : NULL);
				append(o, staln, threadId, rd1, rd2, rdid, r1, r2, summ,
				       ssm1, ssm2, flags1, flags2, prm, mapq, sc, true);
				if(flags1 != NULL) {
					flagscp1.setPrimary(false);
				}
				if(flags2 != NULL) {
					flagscp2.setPrimary(false);
				}
			}
		}
	}
	
	public void reportUnaligned(
			BString o,
			StackedAln staln,
			double threadId,
			final Read rd1,
			final Read rd2,
			final long rdid,
			final AlnSetSumm summ,
			final SeedAlSumm ssm1,
			final SeedAlSumm ssm2,
			final AlnFlags flags1,
			final AlnFlags flags2,
			final PerReadMetrics prm,
			final Mapq mapq,
			final Scoring sc,
			Boolean report2,
			Boolean getLock) {
		append(o, staln, threadId, rd1, rd2, rdid, NULL, NULL, summ,
			       ssm1, ssm2, flags1, flags2, prm, mapq, sc, report2);
	}
	
	public void printAlSumm() {
		
	}
	
	public void reportSeedSummary() {
		
	}
	
	public void reportEmptySeedSummary() {
		
	}
	
	public void mergeMetrics(final ReportingMetrics met) {
		met_.merge(met);
	}
	
	public OutputQueue outq() {
		return oq_;
	}
	
	public void finish(double repThresh, Boolean discord, Boolean mixed, Boolean hadoopOut) {
		if(!quiet_) {
			printAlSumm(met_,repThread, discord, mixed, hadoopOut);
		}
	}
	
	public static DataOutputStream printPct(DataOutputStream os, long num, long denom) throws IOException {
		double pct= 0f;
		if(denom != 0)	{
			pct = 100.0 * (double)num / (double)denom;
		}
		
		os.writeDouble(pct);
		return os;
	}
}