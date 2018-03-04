package com.uwb.bt2j.aligner.sink;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import com.uwb.bt2j.util.OutputQueue;

public class ALNSink {
	protected OutputQueue oq_;
	protected final EList<String> refnames_;
	protected boolean quiet_;
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
			boolean maxed,
			final AlnSetSumm summ,
			final SeedAlSumm ssm1,
			final SeedAlSumm ssm2,
			final AlnFlags flags1,
			final AlnFlags flags2,
			final PerReadMetrics prm,
			final Mapq mapq,
			final Scoring sc,
			boolean reportBoth,
			boolean getLock) {
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
			boolean report2,
			boolean getLock) {
		append(o, staln, threadId, rd1, rd2, rdid, NULL, NULL, summ,
			       ssm1, ssm2, flags1, flags2, prm, mapq, sc, report2);
	}
	
	public static OutputStream printPct(OutputStream os, long num, long denom) throws IOException {
		double pct = 0.0f;
		if(denom != 0) { pct = 100.0 * (double)num / (double)denom; }
		os.write((pct + "%").getBytes());
		return os;
	}
	
	public void printAlSumm(
			ReportingMetrics met,
			int repThresh,   // threshold for uniqueness, or max if no thresh
			boolean discord,       // looked for discordant alignments
			boolean mixed,         // looked for unpaired alignments where paired failed?
			boolean hadoopOut)     // output Hadoop counters?
			{
		// NOTE: there's a filtering step at the very beginning, so everything
		// being reported here is post filtering

		boolean canRep = repThresh != Integer.MAX_VALUE;
		if(hadoopOut) {
			System.err.println( "reporter:counter:Bowtie,Reads processed," + met.nread);
		}
		long totread = met.nread;
		if(totread > 0) {
			System.err.println( "" + met.nread + " reads; of these:" );
		} else {
			System.err.println( "" + totread + " reads");
		}
		long totpair = met.npaired;
		if(totpair > 0) {
			// Paired output
			System.err.print( "  " + totpair + " (");
			printPct(System.err, totpair, totread);
			System.err.print( ") were paired; of these:");

			// Concordants
			System.err.print( "    " + met.nconcord_0 + " (");
			printPct(System.err, met.nconcord_0, met.npaired);
			System.err.println( ") aligned concordantly 0 times");
			if(canRep) {
				// Print the number that aligned concordantly exactly once
				System.err.print( "    " + met.nconcord_uni1 + " (");
				printPct(System.err, met.nconcord_uni1, met.npaired);
				System.err.println( ") aligned concordantly exactly 1 time");
				
				// Print the number that aligned concordantly more than once but
				// fewer times than the limit
				
				System.err.print( "    " + met.nconcord_uni2+met.nconcord_rep + " (");
				printPct(System.err, met.nconcord_uni2+met.nconcord_rep, met.npaired);
				System.err.println( ") aligned concordantly >1 times");
			} else {
				// Print the number that aligned concordantly exactly once
				System.err.print( "    " + met.nconcord_uni1 + " (");
				printPct(System.err, met.nconcord_uni1, met.npaired);
				System.err.println( ") aligned concordantly exactly 1 time");

				// Print the number that aligned concordantly more than once
				System.err.print( "    " + met.nconcord_uni2 + " (");
				printPct(System.err, met.nconcord_uni2, met.npaired);
				System.err.println( ") aligned concordantly >1 times" );
			}
			if(discord) {
				// TODO: what about discoardant and on separate chromosomes?
			
				// Bring out the unaligned pair total so we can subtract discordants
				System.err.println( "    ----" );
				System.err.print( "    " + met.nconcord_0
				     + " pairs aligned concordantly 0 times; of these:");
				// Discordants
				System.err.print( "      " + met.ndiscord + " (");
				printPct(System.err, met.ndiscord, met.nconcord_0);
				System.err.println( ") aligned discordantly 1 time");
			}
			long ncondiscord_0 = met.nconcord_0 - met.ndiscord;
			if(mixed) {
				// Bring out the unaligned pair total so we can subtract discordants
				System.err.println( "    ----" );
				System.err.println( "    " + ncondiscord_0
				     + " pairs aligned 0 times concordantly or discordantly; of these:");
				System.err.println( "      " + (ncondiscord_0 * 2) + " mates make up the pairs; of these:");
				System.err.print( "        " + met.nunp_0_0 + " " + "(");
				printPct(System.err, met.nunp_0_0, ncondiscord_0 * 2);
				System.err.println( ") aligned 0 times");
				if(canRep) {
					// Print the number that aligned exactly once
					System.err.print( "        " + met.nunp_0_uni1 + " (");
					printPct(System.err, met.nunp_0_uni1, ncondiscord_0 * 2);
					System.err.println( ") aligned exactly 1 time" );

					// Print the number that aligned more than once but fewer times
					// than the limit
					System.err.print( "        " + met.nunp_0_uni2+met.nunp_0_rep + " (");
					printPct(System.err, met.nunp_0_uni2+met.nunp_0_rep, ncondiscord_0 * 2);
					System.err.println( ") aligned >1 times" );
				} else {
					// Print the number that aligned exactly once
					System.err.print( "        " + met.nunp_0_uni1 + " (");
					printPct(System.err, met.nunp_0_uni1, ncondiscord_0 * 2);
					System.err.println( ") aligned exactly 1 time");

					// Print the number that aligned more than once but fewer times
					// than the limit
					System.err.print( "        " + met.nunp_0_uni2 + " (");
					printPct(System.err, met.nunp_0_uni2, ncondiscord_0 * 2);
					System.err.println( ") aligned >1 times");
				}
			}
		}
		long totunpair = met.nunpaired;
		if(totunpair > 0) {
			// Unpaired output
			System.err.print( "  " + totunpair + " (");
			printPct(System.err, totunpair, totread);
			System.err.println( ") were unpaired; of these:" );
			
			System.err.print( "    " + met.nunp_0 + " (");
			printPct(System.err, met.nunp_0, met.nunpaired);
			System.err.println( ") aligned 0 times");
			if(hadoopOut) {
				System.err.println( "reporter:counter:Bowtie 2,Unpaired reads with 0 alignments,"
				     + met.nunpaired);
			}
			
			if(canRep) {
				// Print the number that aligned exactly once
				System.err.print( "    " + met.nunp_uni1 + " (");
				printPct(System.err, met.nunp_uni1, met.nunpaired);
				System.err.println( ") aligned exactly 1 time");

				// Print the number that aligned more than once but fewer times
				// than the limit
				System.err.print( "    " + met.nunp_uni2+met.nunp_rep + " (");
				printPct(System.err, met.nunp_uni2+met.nunp_rep, met.nunpaired);
				System.err.println( ") aligned >1 times");
			} else {
				// Print the number that aligned exactly once
				System.err.print( "    " + met.nunp_uni1 + " (");
				printPct(System.err, met.nunp_uni1, met.nunpaired);
				System.err.println( ") aligned exactly 1 time");

				// Print the number that aligned more than once
				System.err.print( "    " + met.nunp_uni2 + " (");
				printPct(System.err, met.nunp_uni2, met.nunpaired);
				System.err.println( ") aligned >1 times");
			}
		}
		long tot_al_cand = totunpair + totpair*2;
		long tot_al =
			(met.nconcord_uni + met.nconcord_rep)*2 +
			(met.ndiscord)*2 +
			met.nunp_0_uni +
			met.nunp_0_rep + 
			met.nunp_uni +
			met.nunp_rep;
		printPct(System.err, tot_al, tot_al_cand);
		System.err.println( " overall alignment rate");
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
	
	public void finish(double repThresh, boolean discord, boolean mixed, boolean hadoopOut) {
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
