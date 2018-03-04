package com.uwb.bt2j.aligner;

import java.io.OutputStream;

import com.uwb.bt2j.indexer.OutFileBuf;

public class PerfMetrics {
	
	
	public PerfMetrics() {
		first = true;
		reset();
	}
	
	public void reset() {
		olm.reset();
		sdm.reset();
		wlm.reset();
		swmSeed.reset();
		swmMate.reset();
		rpm.reset();
		dpSse8Seed.reset();   // 8-bit SSE seed extensions
		dpSse8Mate.reset();   // 8-bit SSE mate finds
		dpSse16Seed.reset();  // 16-bit SSE seed extensions
		dpSse16Mate.reset();  // 16-bit SSE mate finds
		nbtfiltst = 0;
		nbtfiltsc = 0;
		nbtfiltdo = 0;
		
		olmu.reset();
		sdmu.reset();
		wlmu.reset();
		swmuSeed.reset();
		swmuMate.reset();
		rpmu.reset();
		dpSse8uSeed.reset();  // 8-bit SSE seed extensions
		dpSse8uMate.reset();  // 8-bit SSE mate finds
		dpSse16uSeed.reset(); // 16-bit SSE seed extensions
		dpSse16uMate.reset(); // 16-bit SSE mate finds
		nbtfiltst_u = 0;
		nbtfiltsc_u = 0;
		nbtfiltdo_u = 0;
	}
	
	public void merge(
			 OuterLoopMetrics ol,
			 SeedSearchMetrics sd,
			 WalkMetrics wl,
			 SwMetrics swSeed,
			 SwMetrics swMate,
			 ReportingMetrics rm,
			 SSEMetrics dpSse8Ex,
			 SSEMetrics dpSse8Ma,
			 SSEMetrics dpSse16Ex,
			 SSEMetrics dpSse16Ma,
			long nbtfiltst_,
			long nbtfiltsc_,
			long nbtfiltdo_) {
		if(ol != null) {
			olmu.merge(ol);
		}
		if(sd != null) {
			sdmu.merge(sd);
		}
		if(wl != null) {
			wlmu.merge(wl);
		}
		if(swSeed != null) {
			swmuSeed.merge(swSeed);
		}
		if(swMate != null) {
			swmuMate.merge(swMate);
		}
		if(rm != null) {
			rpmu.merge(rm);
		}
		if(dpSse8Ex != null) {
			dpSse8uSeed.merge(dpSse8Ex);
		}
		if(dpSse8Ma != null) {
			dpSse8uMate.merge(dpSse8Ma);
		}
		if(dpSse16Ex != null) {
			dpSse16uSeed.merge(dpSse16Ex);
		}
		if(dpSse16Ma != null) {
			dpSse16uMate.merge(dpSse16Ma);
		}
		nbtfiltst_u += nbtfiltst_;
		nbtfiltsc_u += nbtfiltsc_;
		nbtfiltdo_u += nbtfiltdo_;
	}
	
	public void reportInterval(
			OutFileBuf o,        // file to send output to
			boolean metricsStderr,   // additionally output to stderr?
			boolean total,           // true -> report total, otherwise incremental
			BTString name) // non-null name pointer if is per-read record
		{
			OutputStream stderrSs;
			time_t curtime = time(0);
			String buf);
			if(first) {
				String str =
					/*  1 */ "Time"           "\t"
					/*  2 */ "Read"           "\t"
					/*  3 */ "Base"           "\t"
					/*  4 */ "SameRead"       "\t"
					/*  5 */ "SameReadBase"   "\t"
					/*  6 */ "UnfilteredRead" "\t"
					/*  7 */ "UnfilteredBase" "\t"
					
					/*  8 */ "Paired"         "\t"
					/*  9 */ "Unpaired"       "\t"
					/* 10 */ "AlConUni"       "\t"
					/* 11 */ "AlConRep"       "\t"
					/* 12 */ "AlConFail"      "\t"
					/* 13 */ "AlDis"          "\t"
					/* 14 */ "AlConFailUni"   "\t"
					/* 15 */ "AlConFailRep"   "\t"
					/* 16 */ "AlConFailFail"  "\t"
					/* 17 */ "AlConRepUni"    "\t"
					/* 18 */ "AlConRepRep"    "\t"
					/* 19 */ "AlConRepFail"   "\t"
					/* 20 */ "AlUnpUni"       "\t"
					/* 21 */ "AlUnpRep"       "\t"
					/* 22 */ "AlUnpFail"      "\t"
					
					/* 23 */ "SeedSearch"     "\t"
					/* 24 */ "NRange"         "\t"
					/* 25 */ "NElt"           "\t"
					/* 26 */ "IntraSCacheHit" "\t"
					/* 27 */ "InterSCacheHit" "\t"
					/* 28 */ "OutOfMemory"    "\t"
					/* 29 */ "AlBWOp"         "\t"
					/* 30 */ "AlBWBranch"     "\t"
					/* 31 */ "ResBWOp"        "\t"
					/* 32 */ "ResBWBranch"    "\t"
					/* 33 */ "ResResolve"     "\t"
					/* 34 */ "ResReport"      "\t"
					/* 35 */ "RedundantSHit"  "\t"

					/* 36 */ "BestMinEdit0"   "\t"
					/* 37 */ "BestMinEdit1"   "\t"
					/* 38 */ "BestMinEdit2"   "\t"

					/* 39 */ "ExactAttempts"  "\t"
					/* 40 */ "ExactSucc"      "\t"
					/* 41 */ "ExactRanges"    "\t"
					/* 42 */ "ExactRows"      "\t"
					/* 43 */ "ExactOOMs"      "\t"

					/* 44 */ "1mmAttempts"    "\t"
					/* 45 */ "1mmSucc"        "\t"
					/* 46 */ "1mmRanges"      "\t"
					/* 47 */ "1mmRows"        "\t"
					/* 48 */ "1mmOOMs"        "\t"

					/* 49 */ "UngappedSucc"   "\t"
					/* 50 */ "UngappedFail"   "\t"
					/* 51 */ "UngappedNoDec"  "\t"

					/* 52 */ "DPExLt10Gaps"   "\t"
					/* 53 */ "DPExLt5Gaps"    "\t"
					/* 54 */ "DPExLt3Gaps"    "\t"

					/* 55 */ "DPMateLt10Gaps" "\t"
					/* 56 */ "DPMateLt5Gaps"  "\t"
					/* 57 */ "DPMateLt3Gaps"  "\t"

					/* 58 */ "DP16ExDps"      "\t"
					/* 59 */ "DP16ExDpSat"    "\t"
					/* 60 */ "DP16ExDpFail"   "\t"
					/* 61 */ "DP16ExDpSucc"   "\t"
					/* 62 */ "DP16ExCol"      "\t"
					/* 63 */ "DP16ExCell"     "\t"
					/* 64 */ "DP16ExInner"    "\t"
					/* 65 */ "DP16ExFixup"    "\t"
					/* 66 */ "DP16ExGathSol"  "\t"
					/* 67 */ "DP16ExBt"       "\t"
					/* 68 */ "DP16ExBtFail"   "\t"
					/* 69 */ "DP16ExBtSucc"   "\t"
					/* 70 */ "DP16ExBtCell"   "\t"
					/* 71 */ "DP16ExCoreRej"  "\t"
					/* 72 */ "DP16ExNRej"     "\t"

					/* 73 */ "DP8ExDps"       "\t"
					/* 74 */ "DP8ExDpSat"     "\t"
					/* 75 */ "DP8ExDpFail"    "\t"
					/* 76 */ "DP8ExDpSucc"    "\t"
					/* 77 */ "DP8ExCol"       "\t"
					/* 78 */ "DP8ExCell"      "\t"
					/* 79 */ "DP8ExInner"     "\t"
					/* 80 */ "DP8ExFixup"     "\t"
					/* 81 */ "DP8ExGathSol"   "\t"
					/* 82 */ "DP8ExBt"        "\t"
					/* 83 */ "DP8ExBtFail"    "\t"
					/* 84 */ "DP8ExBtSucc"    "\t"
					/* 85 */ "DP8ExBtCell"    "\t"
					/* 86 */ "DP8ExCoreRej"   "\t"
					/* 87 */ "DP8ExNRej"      "\t"

					/* 88 */ "DP16MateDps"     "\t"
					/* 89 */ "DP16MateDpSat"   "\t"
					/* 90 */ "DP16MateDpFail"  "\t"
					/* 91 */ "DP16MateDpSucc"  "\t"
					/* 92 */ "DP16MateCol"     "\t"
					/* 93 */ "DP16MateCell"    "\t"
					/* 94 */ "DP16MateInner"   "\t"
					/* 95 */ "DP16MateFixup"   "\t"
					/* 96 */ "DP16MateGathSol" "\t"
					/* 97 */ "DP16MateBt"      "\t"
					/* 98 */ "DP16MateBtFail"  "\t"
					/* 99 */ "DP16MateBtSucc"  "\t"
					/* 100 */ "DP16MateBtCell"  "\t"
					/* 101 */ "DP16MateCoreRej" "\t"
					/* 102 */ "DP16MateNRej"    "\t"

					/* 103 */ "DP8MateDps"     "\t"
					/* 104 */ "DP8MateDpSat"   "\t"
					/* 105 */ "DP8MateDpFail"  "\t"
					/* 106 */ "DP8MateDpSucc"  "\t"
					/* 107 */ "DP8MateCol"     "\t"
					/* 108 */ "DP8MateCell"    "\t"
					/* 109 */ "DP8MateInner"   "\t"
					/* 110 */ "DP8MateFixup"   "\t"
					/* 111 */ "DP8MateGathSol" "\t"
					/* 112 */ "DP8MateBt"      "\t"
					/* 113 */ "DP8MateBtFail"  "\t"
					/* 114 */ "DP8MateBtSucc"  "\t"
					/* 115 */ "DP8MateBtCell"  "\t"
					/* 116 */ "DP8MateCoreRej" "\t"
					/* 117 */ "DP8MateNRej"    "\t"

					/* 118 */ "DPBtFiltStart"  "\t"
					/* 119 */ "DPBtFiltScore"  "\t"
					/* 120 */ "DpBtFiltDom"    "\t"
				
				if(name != null) {
					if(o != null) o.writeString("Name\t");
					if(metricsStderr) stderrSs.write("Name\t".getBytes());
				}
				
				if(o != null) o.writeString(str);
				if(metricsStderr) stderrSs.write(str.getBytes());
				first = false;
			}
			
			if(total) mergeIncrementals();
			
			// 0. Read name, if needed
			if(name != null) {
				if(o != null) {
					o.writeString(name.toZBuf());
					o.write('\t');
				}
				if(metricsStderr) {
					stderrSs.write(name + '\t');
				}
			}
				
			// 1. Current time in secs
			buf = String.valueOf(curtime);
			if(metricsStderr) stderrSs.write(buf + '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			OuterLoopMetrics ol = total ? olm : olmu;
			
			// 2. Reads
			buf = String.valueOf(ol.reads);
			if(metricsStderr) stderrSs.write(buf + '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 3. Bases
			buf = String.valueOf(ol.bases);
			if(metricsStderr) stderrSs.write(buf + '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 4. Same-read reads
			buf = String.valueOf(ol.srreads);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 5. Same-read bases
			buf = String.valueOf(ol.srbases);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 6. Unfiltered reads
			buf = String.valueOf(ol.ureads);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 7. Unfiltered bases
			buf = String.valueOf(ol.ubases);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }

			const ReportingMetrics& rp = total ? rpm : rpmu;

			// 8. Paired reads
			buf = String.valueOf(rp.npaired);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 9. Unpaired reads
			buf = String.valueOf(rp.nunpaired);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 10. Pairs with unique concordant alignments
			buf = String.valueOf(rp.nconcord_uni);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 11. Pairs with repetitive concordant alignments
			buf = String.valueOf(rp.nconcord_rep);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 12. Pairs with 0 concordant alignments
			buf = String.valueOf(rp.nconcord_0);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 13. Pairs with 1 discordant alignment
			buf = String.valueOf(rp.ndiscord);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 14. Mates from unaligned pairs that align uniquely
			buf = String.valueOf(rp.nunp_0_uni);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 15. Mates from unaligned pairs that align repetitively
			buf = String.valueOf(rp.nunp_0_rep);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 16. Mates from unaligned pairs that fail to align
			buf = String.valueOf(rp.nunp_0_0);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 17. Mates from repetitive pairs that align uniquely
			buf = String.valueOf(rp.nunp_rep_uni);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 18. Mates from repetitive pairs that align repetitively
			buf = String.valueOf(rp.nunp_rep_rep);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 19. Mates from repetitive pairs that fail to align
			buf = String.valueOf(rp.nunp_rep_0);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 20. Unpaired reads that align uniquely
			buf = String.valueOf(rp.nunp_uni);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 21. Unpaired reads that align repetitively
			buf = String.valueOf(rp.nunp_rep);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 22. Unpaired reads that fail to align
			buf = String.valueOf(rp.nunp_0);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }

			SeedSearchMetrics sd = total ? sdm : sdmu;
			
			// 23. Seed searches
			buf = String.valueOf(sd.seedsearch);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 24. Seed ranges found
			buf = String.valueOf(sd.nrange);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 25. Seed elements found
			buf = String.valueOf(sd.nelt);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 26. Hits in 'current' cache
			buf = String.valueOf(sd.intrahit);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 27. Hits in 'local' cache
			buf = String.valueOf(sd.interhit);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 28. Out of memory
			buf = String.valueOf(sd.ooms);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 29. Burrows-Wheeler ops in aligner
			buf = String.valueOf(sd.bwops);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 30. Burrows-Wheeler branches (edits) in aligner
			buf = String.valueOf(sd.bweds);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			WalkMetrics wl = total ? wlm : wlmu;
			
			// 31. Burrows-Wheeler ops in resolver
			buf = String.valueOf(wl.bwops);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 32. Burrows-Wheeler branches in resolver
			buf = String.valueOf(wl.branches);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 33. Burrows-Wheeler offset resolutions
			buf = String.valueOf(wl.resolves);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 34. Offset reports
			buf = String.valueOf(wl.reports);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			// 35. Redundant seed hit
			buf = String.valueOf(total ? swmSeed.rshit : swmuSeed.rshit);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }

			// 36. # times the best (out of fw/rc) minimum # edits was 0
			buf = String.valueOf(total ? sdm.bestmin0 : sdmu.bestmin0);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 37. # times the best (out of fw/rc) minimum # edits was 1
			buf = String.valueOf(total ? sdm.bestmin1 : sdmu.bestmin1);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 38. # times the best (out of fw/rc) minimum # edits was 2
			buf = String.valueOf(total ? sdm.bestmin2 : sdmu.bestmin2);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			// 39. Exact aligner attempts
			buf = String.valueOf(total ? swmSeed.exatts : swmuSeed.exatts);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 40. Exact aligner successes
			buf = String.valueOf(total ? swmSeed.exsucc : swmuSeed.exsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 41. Exact aligner ranges
			buf = String.valueOf(total ? swmSeed.exranges : swmuSeed.exranges);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 42. Exact aligner rows
			buf = String.valueOf(total ? swmSeed.exrows : swmuSeed.exrows);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 43. Exact aligner OOMs
			buf = String.valueOf(total ? swmSeed.exooms : swmuSeed.exooms);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }

			// 44. 1mm aligner attempts
			buf = String.valueOf(total ? swmSeed.mm1atts : swmuSeed.mm1atts);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 45. 1mm aligner successes
			buf = String.valueOf(total ? swmSeed.mm1succ : swmuSeed.mm1succ);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 46. 1mm aligner ranges
			buf = String.valueOf(total ? swmSeed.mm1ranges : swmuSeed.mm1ranges);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 47. 1mm aligner rows
			buf = String.valueOf(total ? swmSeed.mm1rows : swmuSeed.mm1rows);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 48. 1mm aligner OOMs
			buf = String.valueOf(total ? swmSeed.mm1ooms : swmuSeed.mm1ooms);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }

			// 49 Ungapped aligner success
			buf = String.valueOf(total ? swmSeed.ungapsucc : swmuSeed.ungapsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 50. Ungapped aligner fail
			buf = String.valueOf(total ? swmSeed.ungapfail : swmuSeed.ungapfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 51. Ungapped aligner no decision
			buf = String.valueOf(total ? swmSeed.ungapnodec : swmuSeed.ungapnodec);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }

			// 52. # seed-extend DPs with < 10 gaps
			buf = String.valueOf(total ? swmSeed.sws10 : swmuSeed.sws10);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 53. # seed-extend DPs with < 5 gaps
			buf = String.valueOf(total ? swmSeed.sws5 : swmuSeed.sws5);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 54. # seed-extend DPs with < 3 gaps
			buf = String.valueOf(total ? swmSeed.sws3 : swmuSeed.sws3);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }

			// 55. # seed-extend DPs with < 10 gaps
			buf = String.valueOf(total ? swmMate.sws10 : swmuMate.sws10);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 56. # seed-extend DPs with < 5 gaps
			buf = String.valueOf(total ? swmMate.sws5 : swmuMate.sws5);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 57. # seed-extend DPs with < 3 gaps
			buf = String.valueOf(total ? swmMate.sws3 : swmuMate.sws3);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			SSEMetrics dpSse16s = total ? dpSse16Seed : dpSse16uSeed;
			
			// 58. 16-bit SSE seed-extend DPs tried
			buf = String.valueOf(dpSse16s.dp);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 59. 16-bit SSE seed-extend DPs saturated
			buf = String.valueOf(dpSse16s.dpsat);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 60. 16-bit SSE seed-extend DPs failed
			buf = String.valueOf(dpSse16s.dpfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 61. 16-bit SSE seed-extend DPs succeeded
			buf = String.valueOf(dpSse16s.dpsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 62. 16-bit SSE seed-extend DP columns completed
			buf = String.valueOf(dpSse16s.col);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 63. 16-bit SSE seed-extend DP cells completed
			buf = String.valueOf(dpSse16s.cell);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 64. 16-bit SSE seed-extend DP inner loop iters completed
			buf = String.valueOf(dpSse16s.inner);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 65. 16-bit SSE seed-extend DP fixup loop iters completed
			buf = String.valueOf(dpSse16s.fixup);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 66. 16-bit SSE seed-extend DP gather, cells with potential solutions
			buf = String.valueOf(dpSse16s.gathsol);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 67. 16-bit SSE seed-extend DP backtrace attempts
			buf = String.valueOf(dpSse16s.bt);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 68. 16-bit SSE seed-extend DP failed backtrace attempts
			buf = String.valueOf(dpSse16s.btfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 69. 16-bit SSE seed-extend DP succesful backtrace attempts
			buf = String.valueOf(dpSse16s.btsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 70. 16-bit SSE seed-extend DP backtrace cells
			buf = String.valueOf(dpSse16s.btcell);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 71. 16-bit SSE seed-extend DP core-diag rejections
			buf = String.valueOf(dpSse16s.corerej);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 72. 16-bit SSE seed-extend DP N rejections
			buf = String.valueOf(dpSse16s.nrej);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			SSEMetrics dpSse8s = total ? dpSse8Seed : dpSse8uSeed;
			
			// 73. 8-bit SSE seed-extend DPs tried
			buf = String.valueOf(dpSse8s.dp);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 74. 8-bit SSE seed-extend DPs saturated
			buf = String.valueOf(dpSse8s.dpsat);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 75. 8-bit SSE seed-extend DPs failed
			buf = String.valueOf(dpSse8s.dpfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 76. 8-bit SSE seed-extend DPs succeeded
			buf = String.valueOf(dpSse8s.dpsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 77. 8-bit SSE seed-extend DP columns completed
			buf = String.valueOf(dpSse8s.col);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 78. 8-bit SSE seed-extend DP cells completed
			buf = String.valueOf(dpSse8s.cell);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 79. 8-bit SSE seed-extend DP inner loop iters completed
			buf = String.valueOf(dpSse8s.inner);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 80. 8-bit SSE seed-extend DP fixup loop iters completed
			buf = String.valueOf(dpSse8s.fixup);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 81. 16-bit SSE seed-extend DP gather, cells with potential solutions
			buf = String.valueOf(dpSse8s.gathsol);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 82. 16-bit SSE seed-extend DP backtrace attempts
			buf = String.valueOf(dpSse8s.bt);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 83. 16-bit SSE seed-extend DP failed backtrace attempts
			buf = String.valueOf(dpSse8s.btfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 84. 16-bit SSE seed-extend DP succesful backtrace attempts
			buf = String.valueOf(dpSse8s.btsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 85. 16-bit SSE seed-extend DP backtrace cells
			buf = String.valueOf(dpSse8s.btcell);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 86. 16-bit SSE seed-extend DP core-diag rejections
			buf = String.valueOf(dpSse8s.corerej);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 87. 16-bit SSE seed-extend DP N rejections
			buf = String.valueOf(dpSse8s.nrej);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			SSEMetrics dpSse16m = total ? dpSse16Mate : dpSse16uMate;
			
			// 88. 16-bit SSE mate-finding DPs tried
			buf = String.valueOf(dpSse16m.dp);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 89. 16-bit SSE mate-finding DPs saturated
			buf = String.valueOf(dpSse16m.dpsat);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 90. 16-bit SSE mate-finding DPs failed
			buf = String.valueOf(dpSse16m.dpfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 91. 16-bit SSE mate-finding DPs succeeded
			buf = String.valueOf(dpSse16m.dpsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 92. 16-bit SSE mate-finding DP columns completed
			buf = String.valueOf(dpSse16m.col);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 93. 16-bit SSE mate-finding DP cells completed
			buf = String.valueOf(dpSse16m.cell);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 94. 16-bit SSE mate-finding DP inner loop iters completed
			buf = String.valueOf(dpSse16m.inner);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 95. 16-bit SSE mate-finding DP fixup loop iters completed
			buf = String.valueOf(dpSse16m.fixup);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 96. 16-bit SSE mate-finding DP gather, cells with potential solutions
			buf = String.valueOf(dpSse16m.gathsol);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 97. 16-bit SSE mate-finding DP backtrace attempts
			buf = String.valueOf(dpSse16m.bt);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 98. 16-bit SSE mate-finding DP failed backtrace attempts
			buf = String.valueOf(dpSse16m.btfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 99. 16-bit SSE mate-finding DP succesful backtrace attempts
			buf = String.valueOf(dpSse16m.btsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 100. 16-bit SSE mate-finding DP backtrace cells
			buf = String.valueOf(dpSse16m.btcell);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 101. 16-bit SSE mate-finding DP core-diag rejections
			buf = String.valueOf(dpSse16m.corerej);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 102. 16-bit SSE mate-finding DP N rejections
			buf = String.valueOf(dpSse16m.nrej);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			SSEMetrics dpSse8m = total ? dpSse8Mate : dpSse8uMate;
			
			// 103. 8-bit SSE mate-finding DPs tried
			buf = String.valueOf(dpSse8m.dp);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 104. 8-bit SSE mate-finding DPs saturated
			buf = String.valueOf(dpSse8m.dpsat);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 105. 8-bit SSE mate-finding DPs failed
			buf = String.valueOf(dpSse8m.dpfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 106. 8-bit SSE mate-finding DPs succeeded
			buf = String.valueOf(dpSse8m.dpsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 107. 8-bit SSE mate-finding DP columns completed
			buf = String.valueOf(dpSse8m.col);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 108. 8-bit SSE mate-finding DP cells completed
			buf = String.valueOf(dpSse8m.cell);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 109. 8-bit SSE mate-finding DP inner loop iters completed
			buf = String.valueOf(dpSse8m.inner);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 110. 8-bit SSE mate-finding DP fixup loop iters completed
			buf = String.valueOf(dpSse8m.fixup);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 111. 16-bit SSE mate-finding DP gather, cells with potential solutions
			buf = String.valueOf(dpSse8m.gathsol);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 112. 16-bit SSE mate-finding DP backtrace attempts
			buf = String.valueOf(dpSse8m.bt);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 113. 16-bit SSE mate-finding DP failed backtrace attempts
			buf = String.valueOf(dpSse8m.btfail);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 114. 16-bit SSE mate-finding DP succesful backtrace attempts
			buf = String.valueOf(dpSse8m.btsucc);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 115. 16-bit SSE mate-finding DP backtrace cells
			buf = String.valueOf(dpSse8m.btcell);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 116. 16-bit SSE mate-finding DP core rejections
			buf = String.valueOf(dpSse8m.corerej);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 117. 16-bit SSE mate-finding N rejections
			buf = String.valueOf(dpSse8m.nrej);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			
			// 118. Backtrace candidates filtered due to starting cell
			buf = String.valueOf(total ? nbtfiltst : nbtfiltst_u);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 119. Backtrace candidates filtered due to low score
			buf = String.valueOf(total ? nbtfiltsc : nbtfiltsc_u);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			// 120. Backtrace candidates filtered due to domination
			buf = String.valueOf(total ? nbtfiltdo : nbtfiltdo_u);
			if(metricsStderr) stderr.write( buf << '\t');
			if(o != null) { o.writeString(buf); o.write('\t'); }
			

			if(o != null) { o.write('\n'); }
			if(metricsStderr) System.err.println(stderrSs);
			if(!total) mergeIncrementals();
		}
}
