package com.uwb.bt2j.aligner.seed;

import com.uwb.bt2j.aligner.Edit;
import com.uwb.bt2j.aligner.PerReadMetrics;
import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.aligner.Scoring;
import com.uwb.bt2j.aligner.cache.QKey;
import com.uwb.bt2j.aligner.cache.QVal;
import com.uwb.bt2j.indexer.Ebwt;
import com.uwb.bt2j.indexer.SideLocus;
import com.uwb.bt2j.util.strings.BTDnaString;
import com.uwb.bt2j.util.strings.BTString;
import com.uwb.bt2j.util.types.EList;

import javafx.util.Pair;

public class SeedAligner {
	public Ebwt ebwtFw_;
	public Ebwt ebwtBw_;
	public Scoring sc_;
	public InstantiatedSeed s_;
	public Read read_;
	public BTDnaString seq_;
	public BTString qual_;
	public double off_;
	public boolean fw_;
	public EList<Edit> edits_;
	public AlignmentCacheIface ca_;
	public EList<Double> offIdx2off_;
	public long bwops_;
	public long bwedits_;
	public BTDnaString tmprfdnastr_;
	public BTDnaString tmpdnastr_;
	
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
		int len = seeds.get(0).len; // assume they're all the same length
		// Calc # seeds within read interval
		int nseeds = 1;
		if((int)read.length() - (int)off > len) {
			nseeds += ((int)read.length() - (int)off - len) / per;
		}
		for(int i = 0; i < nseeds; i++) {
			offIdx2off_.push_back((double) (per * i + (int)off));
		}
		Pair<Integer, Integer> ret = new Pair(0,0);
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
				int seedlen = seeds.get(0).len;
				// Extract the seed sequence at this offset
				// If fw == true, we extract the characters from i*per to
				// i*(per-1) (exclusive).  If fw == false, 
				instantiateSeq(
					read,
					sr.seqs(fw).get(i),
					sr.quals(fw).get(i),
					Integer.min((int)seedlen, (int)read.length()),
					depth,
					fw);
				QKey qk = new QKey(sr.seqs(fw).get(i));
				// For each search strategy
				EList<InstantiatedSeed> iss = sr.instantiatedSeeds(fw, i);
				for(int j = 0; j < (int)seeds.size(); j++) {
					iss.expand();
					InstantiatedSeed is = iss.back();
					if(seeds.get(j).instantiate(
						read,
						sr.seqs(fw).get(i),
						sr.quals(fw).get(i),
						pens,
						depth,
						i,
						j,
						fw,
						is))
					{
						// Can we fill this seed hit in from the cache?
						ret = new Pair(ret.getKey() + 1, ret.getValue());
						if(fwi == 0) { instFw = new Pair(instFw.getKey() + 1, instFw.getValue()); } 
						else { instRc = new Pair(instRc.getKey() + 1, instRc.getValue()); }
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
	
	public void searchAllSeeds(
			EList<Seed> seeds,
			Ebwt ebwtFw,
			Ebwt ebwtBw,
			Read read,
			Scoring pens,
			AlignmentCacheIface cache,
			SeedResults sr,
			SeedSearchMetrics met,
			PerReadMetrics prm) {
		ebwtFw_ = ebwtFw;
		ebwtBw_ = ebwtBw;
		sc_ = pens;
		read_ = read;
		ca_ = cache;
		bwops_ = bwedits_ = 0;
		long possearches = 0, seedsearches = 0, intrahits = 0, interhits = 0, ooms = 0;
		// For each instantiated seed
		for(int i = 0; i < (int)sr.numOffs(); i++) {
			double off = sr.idx2off(i);
			for(int fwi = 0; fwi < 2; fwi++) {
				Boolean fw = (fwi == 0);
				EList<InstantiatedSeed> iss = sr.instantiatedSeeds(fw, i);
				if(iss.empty()) {
					// Cache hit in an across-read cache
					continue;
				}
				QVal qv;
				seq_  = sr.seqs(fw).get(i);  // seed sequence
				qual_ = sr.quals(fw).get(i); // seed qualities
				off_  = off;              // seed offset (from 5')
				fw_   = fw;               // seed orientation
				// Tell the cache that we've started aligning, so the cache can
				// expect a series of on-the-fly updates
				int ret = cache.beginAlign(seq_, qual_, qv);
				if(ret == -1) {
					// Out of memory when we tried to add key to map
					ooms++;
					continue;
				}
				Boolean abort = false;
				if(ret == 0) {
					// Not already in cache
					possearches++;
					for(int j = 0; j < iss.size(); j++) {
						// Set seq_ and qual_ appropriately, using the seed sequences
						// and qualities already installed in SeedResults
						s_ = iss.get(j);
						// Do the search with respect to seq_, qual_ and s_.
						if(!searchSeedBi()) {
							// Memory exhausted during search
							ooms++;
							abort = true;
							break;
						}
						seedsearches++;
					}
					if(!abort) {
						qv = cache.finishAlign();
					}
				} else {
					// Already in cache
					intrahits++;
				}
				if(qv.valid()) {
					sr.add(
						qv,    // range of ranges in cache
						cache.current(), // cache
						i,     // seed index (from 5' end)
						fw);   // whether seed is from forward read
				}
			}
		}
		prm.nSeedRanges = sr.numRanges();
		prm.nSeedElts = sr.numElts();
		prm.nSeedRangesFw = sr.numRangesFw();
		prm.nSeedRangesRc = sr.numRangesRc();
		prm.nSeedEltsFw = sr.numEltsFw();
		prm.nSeedEltsRc = sr.numEltsRc();
		prm.seedMedian = (sr.medianHitsPerSeed() + 0.5);
		prm.seedMean = sr.averageHitsPerSeed();

		prm.nSdFmops += bwops_;
		met.seedsearch += seedsearches;
		met.nrange += sr.numRanges();
		met.nelt += sr.numElts();
		met.possearch += possearches;
		met.intrahit += intrahits;
		met.interhit += interhits;
		met.ooms += ooms;
		met.bwops += bwops_;
		met.bweds += bwedits_;
	}
	
	public boolean sanityPartial(
			Ebwt ebwtFw,
			Ebwt ebwtBw,
			BTDnaString seq,
			int dep,
			double len,
			boolean do1mm,
			long topFw,
			long botfw,
			long topbw,
			long botbw) {
		tmpdnastr_.clear();
		for(int i = dep; i < len; i++) {
			tmpdnastr_.append(seq.get(i));
		}
		long top_fw = 0, bot_fw = 0;
		ebwtFw.contains(tmpdnastr_, top_fw, bot_fw);
		if(do1mm && ebwtBw != null) {
			tmpdnastr_.reverse();
			long top_bw = 0, bot_bw = 0;
			ebwtBw.contains(tmpdnastr_, top_bw, bot_bw);
		}
		return true;
	}
	
	private void INIT_LOCS(long top, long bot, SideLocus tloc, SideLocus bloc, Ebwt e) {
		if(bot - top == 1) {
			tloc.initFromRow(top, e.eh(), e.ebwt());
			bloc.invalidate();
		} else {
			SideLocus.initFromTopBot(top, bot, e.eh(), e.ebwt(), tloc, bloc);
		}
	}
	
	public double exactSweep(
			Ebwt ebwt,
			Read read,
			Scoring sc,
			Boolean nofw,
			Boolean norc,
			double mineMax,
			double mineFw,
			double mineRc,
			Boolean repex,
			SeedResults hits,
			SeedSearchMetrics met) {
		long top = 0, bot = 0;
		SideLocus tloc, bloc;
		int len = read.length();
		double nelt = 0;
		for(int fwi = 0; fwi < 2; fwi++) {
			boolean fw = (fwi == 0);
			if( fw && nofw) continue;
			if(!fw && norc) continue;
			BTDnaString seq = fw ? read.patFw : read.patRc;
			int ftabLen = ebwt.eh().ftabChars();
			int dep = 0;
			double nedit = 0;
			Boolean done = false;
			while(dep < len && !done) {
				top = bot = 0;
				double left = len - dep;
				Boolean doFtab = ftabLen > 1 && left >= (double)ftabLen;
				if(doFtab) {
					// Does N interfere with use of Ftab?
					for(int i = 0; i < (double)ftabLen; i++) {
						String c = seq.get(len-dep-1-i);
						if(c.compareTo("3") > 0) {
							doFtab = false;
							break;
						}
					}
				}
				if(doFtab) {
					// Use ftab
					ebwt.ftabLoHi(seq, len - dep - ftabLen, false, top, bot);
					dep += (double)ftabLen;
				} else {
					// Use fchr
					String c = seq.get(len-dep-1);
					if(c.compareTo("4") < 0) {
						top = ebwt.fchr().get(c);
						bot = ebwt.fchr().get(c+1);
					}
					dep++;
				}
				if(bot <= top) {
					nedit++;
					if(nedit >= mineMax) {
						if(fw) { mineFw = nedit; } else { mineRc = nedit; }
						break;
					}
					continue;
				}
				INIT_LOCS(top, bot, tloc, bloc, ebwt);
				// Keep going
				while(dep < len) {
					String c = seq.get(len-dep-1);
					if(c.compareTo("3") > 0) {
						top = bot = 0;
					} else {
						if(bloc.valid()) {
							bwops_ += 2;
							top = ebwt.mapLF(tloc);
							bot = ebwt.mapLF(bloc);
						} else {
							bwops_++;
							top = ebwt.mapLF1(top, tloc, c);
							if(top == IndexTypes.OFF_MASK) {
								top = bot = 0;
							} else {
								bot = top+1;
							}
						}
					}
					if(bot <= top) {
						nedit++;
						if(nedit >= mineMax) {
							if(fw) { mineFw = nedit; } else { mineRc = nedit; }
							done = true;
						}
						break;
					}
					INIT_LOCS(top, bot, tloc, bloc, ebwt);
					dep++;
				}
				if(done) {
					break;
				}
				if(dep == len) {
					// Set the minimum # edits
					if(fw) { mineFw = nedit; } else { mineRc = nedit; }
					// Done
					if(nedit == 0 && bot > top) {
						if(repex) {
							// This is an exact hit
							long score = len * sc.match();
							if(fw) {
								hits.addExactEeFw(top, bot, null, null, fw, score);
							} else {
								hits.addExactEeRc(top, bot, null, null, fw, score);
							}
						}
						nelt += (bot - top);
					}
					break;
				}
				dep++;
			}
		}
		return nelt;
	}
	
	public Boolean oneMmSearch(
			Ebwt        ebwtFw, // BWT index
			Ebwt        ebwtBw, // BWT' index
			Read        read,   // read to align
			Scoring     sc,     // scoring
			long            minsc,  // minimum score
			Boolean               nofw,   // don't align forward read
			Boolean               norc,   // don't align revcomp read
			Boolean               local,  // 1mm hits must be legal local alignments
			Boolean               repex,  // report 0mm hits?
			Boolean               rep1mm, // report 1mm hits?
			SeedResults       hits,   // holds all the seed hits (and exact hit)
			SeedSearchMetrics met) {   // metrics{
		int len = read.length();
		int nceil = sc.nCeil.f<Integer>((double)len);
		double ns = read.ns();
		if(ns > 1) {
			// Can't align this with <= 1 mismatches
			return false;
		} else if(ns == 1 && !rep1mm) {
			// Can't align this with 0 mismatches
			return false;
		}
		double halfFw = len >> 1;
		double halfBw = len >> 1;
		if((len & 1) != 0) {
			halfBw++;
		}
		SideLocus tloc, bloc;
		long t[], b[];   // dest BW ranges for BWT
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		long tp[], bp[]; // dest BW ranges for BWT'
		tp[0] = tp[1] = tp[2] = tp[3] = 0;
		bp[0] = bp[1] = bp[2] = bp[3] = 0;
		long top = 0, bot = 0, topp = 0, botp = 0;
		// Align fw read / rc read
		Boolean results = false;
		for(int fwi = 0; fwi < 2; fwi++) {
			Boolean fw = (fwi == 0);
			if( fw && nofw) continue;
			if(!fw && norc) continue;
			// Align going right-to-left, left-to-right
			int lim = rep1mm ? 2 : 1;
			for(int ebwtfwi = 0; ebwtfwi < lim; ebwtfwi++) {
				Boolean ebwtfw = (ebwtfwi == 0);
				Ebwt ebwt  = (ebwtfw ? ebwtFw : ebwtBw);
				Ebwt ebwtp = (ebwtfw ? ebwtBw : ebwtFw);
				BTDnaString seq =
					(fw ? (ebwtfw ? read.patFw : read.patFwRev) :
					      (ebwtfw ? read.patRc : read.patRcRev));
				BTString qual =
					(fw ? (ebwtfw ? read.qual    : read.qualRev) :
					      (ebwtfw ? read.qualRev : read.qual));
				int ftabLen = ebwt.eh().ftabChars();
				double nea = ebwtfw ? halfFw : halfBw;
				// Check if there's an N in the near portion
				Boolean skip = false;
				for(int dep = 0; dep < nea; dep++) {
					if(seq.get(len-dep-1).compareTo("3") > 0) {
						skip = true;
						break;
					}
				}
				if(skip) {
					continue;
				}
				double dep = 0;
				// Align near half
				if(ftabLen > 1 && (double)ftabLen <= nea) {
					// Use ftab to jump partway into near half
					Boolean rev = !ebwtfw;
					ebwt.ftabLoHi(seq, len - ftabLen, rev, top, bot);
					if(rep1mm) {
						ebwtp.ftabLoHi(seq, len - ftabLen, rev, topp, botp);
					}
					if(bot - top == 0) {
						continue;
					}
					int c = seq.get(len - ftabLen);
					t[c] = top; b[c] = bot;
					tp[c] = topp; bp[c] = botp;
					dep = ftabLen;
					// initialize tloc, bloc??
				} else {
					// Use fchr to jump in by 1 pos
					int c = seq[len-1];
					top = topp = tp[c] = ebwt.fchr()[c];
					bot = botp = bp[c] = ebwt.fchr()[c+1];
					if(bot - top == 0) {
						continue;
					}
					dep = 1;
					// initialize tloc, bloc??
				}
				INIT_LOCS(top, bot, tloc, bloc, *ebwt);
				Boolean do_continue = false;
				for(; dep < nea; dep++) {
					int rdc = seq[len - dep - 1];
					tp[0] = tp[1] = tp[2] = tp[3] = topp;
					bp[0] = bp[1] = bp[2] = bp[3] = botp;
					if(bloc.valid()) {
						bwops_++;
						t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
						ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
						SANITY_CHECK_4TUP(t, b, tp, bp);
						top = t[rdc]; bot = b[rdc];
						if(bot <= top) {
							do_continue = true;
							break;
						}
						topp = tp[rdc]; botp = bp[rdc];
					} else {
						bwops_++;
						top = ebwt.mapLF1(top, tloc, rdc);
						if(top == IndexTypes.OFF_MASK) {
							do_continue = true;
							break;
						}
						bot = top + 1;
						t[rdc] = top; b[rdc] = bot;
						tp[rdc] = topp; bp[rdc] = botp;
						// topp/botp stay the same
					}
					INIT_LOCS(top, bot, tloc, bloc, *ebwt);
				}
				if(do_continue) {
					continue;
				}
				// Align far half
				for(; dep < len; dep++) {
					int rdc = seq[len-dep-1];
					int quc = qual[len-dep-1];
					if(rdc > 3 && nceil == 0) {
						break;
					}
					tp[0] = tp[1] = tp[2] = tp[3] = topp;
					bp[0] = bp[1] = bp[2] = bp[3] = botp;
					int clo = 0, chi = 3;
					Boolean match = true;
					if(bloc.valid()) {
						bwops_++;
						t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
						ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
						SANITY_CHECK_4TUP(t, b, tp, bp);
						match = rdc < 4;
						top = t[rdc]; bot = b[rdc];
						topp = tp[rdc]; botp = bp[rdc];
					} else {
						bwops_++;
						clo = ebwt.mapLF1(top, tloc);
						match = (clo == rdc);
						if(clo < 0) {
							break; // Hit the $
						} else {
							t[clo] = top;
							b[clo] = bot = top + 1;
						}
						bp[clo] = botp;
						tp[clo] = topp;
						chi = clo;
					}
					//assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
					if(rep1mm && (ns == 0 || rdc > 3)) {
						for(int j = clo; j <= chi; j++) {
							if(j == rdc || b[j] == t[j]) {
								// Either matches read or isn't a possibility
								continue;
							}
							// Potential mismatch - next, try
							double depm = dep + 1;
							long topm = t[j], botm = b[j];
							long topmp = tp[j], botmp = bp[j];
							long tm[4], bm[4];   // dest BW ranges for BWT
							tm[0] = t[0]; tm[1] = t[1];
							tm[2] = t[2]; tm[3] = t[3];
							bm[0] = b[0]; bm[1] = t[1];
							bm[2] = b[2]; bm[3] = t[3];
							long tmp[4], bmp[4]; // dest BW ranges for BWT'
							tmp[0] = tp[0]; tmp[1] = tp[1];
							tmp[2] = tp[2]; tmp[3] = tp[3];
							bmp[0] = bp[0]; bmp[1] = tp[1];
							bmp[2] = bp[2]; bmp[3] = tp[3];
							SideLocus tlocm, blocm;
							INIT_LOCS(topm, botm, tlocm, blocm, *ebwt);
							for(; depm < len; depm++) {
								int rdcm = seq[len - depm - 1];
								tmp[0] = tmp[1] = tmp[2] = tmp[3] = topmp;
								bmp[0] = bmp[1] = bmp[2] = bmp[3] = botmp;
								if(blocm.valid()) {
									bwops_++;
									tm[0] = tm[1] = tm[2] = tm[3] =
									bm[0] = bm[1] = bm[2] = bm[3] = 0;
									ebwt.mapBiLFEx(tlocm, blocm, tm, bm, tmp, bmp);
									SANITY_CHECK_4TUP(tm, bm, tmp, bmp);
									topm = tm[rdcm]; botm = bm[rdcm];
									topmp = tmp[rdcm]; botmp = bmp[rdcm];
									if(botm <= topm) {
										break;
									}
								} else {
									bwops_++;
									topm = ebwt.mapLF1(topm, tlocm, rdcm);
									if(topm == IndexTypes.OFF_MASK) {
										break;
									}
									botm = topm + 1;
									// topp/botp stay the same
								}
								INIT_LOCS(topm, botm, tlocm, blocm, *ebwt);
							}
							if(depm == len) {
								// Success; this is a 1MM hit
								double off5p = dep;  // offset from 5' end of read
								double offstr = dep; // offset into patFw/patRc
								if(fw == ebwtfw) {
									off5p = len - off5p - 1;
								}
								if(!ebwtfw) {
									offstr = len - offstr - 1;
								}
								Edit e((double)off5p, j, rdc, EDIT_TYPE_MM, false);
								results = true;
								long score = (len - 1) * sc.match();
								// In --local mode, need to double-check that
								// end-to-end alignment doesn't violate  local
								// alignment principles.  Specifically, it
								// shouldn't to or below 0 anywhere in the middle.
								int pen = sc.score(rdc, (int)(1 << j), quc - 33);
								score += pen;
								Boolean valid = true;
								if(local) {
									long locscore_fw = 0, locscore_bw = 0;
									for(double i = 0; i < len; i++) {
										if(i == dep) {
											if(locscore_fw + pen <= 0) {
												valid = false;
												break;
											}
											locscore_fw += pen;
										} else {
											locscore_fw += sc.match();
										}
										if(len-i-1 == dep) {
											if(locscore_bw + pen <= 0) {
												valid = false;
												break;
											}
											locscore_bw += pen;
										} else {
											locscore_bw += sc.match();
										}
									}
								}
								if(valid) {
									valid = score >= minsc;
								}
								if(valid) {
	#ifndef NDEBUG
									BTDnaString rf = tmprfdnastr_;
									rf.clear();
									edits_.clear();
									edits_.push_back(e);
									if(!fw) Edit.invertPoss(edits_, len, false);
									Edit::toRef(fw ? read.patFw : read.patRc, edits_, rf);
									if(!fw) Edit.invertPoss(edits_, len, false);
	#endif
									long toprep = ebwtfw ? topm : topmp;
									long botrep = ebwtfw ? botm : botmp;
									hits.add1mmEe(toprep, botrep, e, null, fw, score);
								}
							}
						}
					}
					if(bot > top && match) {
						if(dep == len-1) {
							// Success; this is an exact hit
							if(ebwtfw && repex) {
								if(fw) {
									results = true;
									long score = len * sc.match();
									hits.addExactEeFw(
										ebwtfw ? top : topp,
										ebwtfw ? bot : botp,
										null, nul, fw, score);
								} else {
									results = true;
									long score = len * sc.match();
									hits.addExactEeRc(
										ebwtfw ? top : topp,
										ebwtfw ? bot : botp,
										null, null, fw, score);
								}
							}
							break; // End of far loop
						} else {
							INIT_LOCS(top, bot, tloc, bloc, *ebwt);
						}
					} else {
						break; // End of far loop
					}
				} // for(; dep < len; dep++)
			} // for(int ebwtfw = 0; ebwtfw < 2; ebwtfw++)
		} // for(int fw = 0; fw < 2; fw++)
		return results;
	}
	
	public Boolean searchSeedBi() {
		return searchSeedBi(
				0, 0,
				0, 0, 0, 0,
				new SideLocus(), new SideLocus(),
				s_.cons[0], s_.cons[1], s_.cons[2], s_.overall,
				null);
	}
	
	public Boolean searchSeedBi(
			int step,             // depth into steps_[] array
			int depth,            // recursion depth
			long topf,        // top in BWT
			long botf,        // bot in BWT
			long topb,        // top in BWT'
			long botb,        // bot in BWT'
			SideLocus tloc,       // locus for top (perhaps unititialized)
			SideLocus bloc,       // locus for bot (perhaps unititialized)
			Constraint c0,        // constraints to enforce in seed zone 0
			Constraint c1,        // constraints to enforce in seed zone 1
			Constraint c2,        // constraints to enforce in seed zone 2
			Constraint overall,   // overall constraints to enforce
			DoublyLinkedList<Edit> prevEdit  // previous edit
			) {
		InstantiatedSeed s = s_;

		if(step == (int)s.steps.size()) {
			// Finished aligning seed
			if(!reportHit(topf, botf, topb, botb, seq_.length(), prevEdit)) {
				return false; // Memory exhausted
			}
			return true;
		}

		int off;
		long tp[4], bp[4]; // dest BW ranges for "prime" index
		if(step == 0) {
			// Just starting
			off = s.steps[0];
			Boolean ltr = off > 0;
			off = abs(off)-1;
			// Check whether/how far we can jump using ftab or fchr
			int ftabLen = ebwtFw_.eh().ftabChars();
			if(ftabLen > 1 && ftabLen <= s.maxjump) {
				if(!ltr) {
					off = off - ftabLen + 1;
				}
				ebwtFw_.ftabLoHi(seq_, off, false, topf, botf);

				if(ebwtBw_ != null) {
					ebwtBw_.ftabLoHi(seq_, off, false, topb, botb);
				}
				if(botf - topf == 0) return true;

				step += ftabLen;
			} else if(s.maxjump > 0) {
				// Use fchr
				int c = (seq_)[off];
				topf = topb = ebwtFw_.fchr()[c];
				botf = botb = ebwtFw_.fchr()[c+1];
				if(botf - topf == 0) return true;
				step++;
			} else {
				topf = topb = 0;
				botf = botb = ebwtFw_.fchr()[4];
			}
			if(step == (int)s.steps.size()) {
				// Finished aligning seed
				if(!reportHit(topf, botf, topb, botb, seq_.length(), prevEdit)) {
					return false; // Memory exhausted
				}
				return true;
			}
			nextLocsBi(tloc, bloc, topf, botf, topb, botb, step);

		}
		long t[4], b[4]; // dest BW ranges
		Constraint zones[3] = { c0, c1, c2 };
		for(int i = step; i < (int)s.steps.size(); i++) {
			off = s.steps[i];
			Boolean ltr = off > 0;
			Ebwt ebwt = ltr ? ebwtBw_ : ebwtFw_;
			if(ltr) {
				tp[0] = tp[1] = tp[2] = tp[3] = topf;
				bp[0] = bp[1] = bp[2] = bp[3] = botf;
			} else {
				tp[0] = tp[1] = tp[2] = tp[3] = topb;
				bp[0] = bp[1] = bp[2] = bp[3] = botb;
			}
			t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
			if(bloc.valid()) {
				// Range delimited by tloc/bloc has size >1.  If size == 1,
				// we use a simpler query (see if(!bloc.valid()) blocks below)
				bwops_++;
				ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
			}
			long tf = ltr ? tp : t, tb = ltr ? t : tp;
			long bf = ltr ? bp : b, bb = ltr ? b : bp;
			off = abs(off)-1;
			//
			Boolean leaveZone = s.zones[i].first < 0;
			//Boolean leaveZoneIns = zones_[i].second < 0;
			Constraint cons    = zones[abs(s.zones[i].first)];
			//Constraint& insCons = *zones[abs(s.zones[i].second)];
			int c = (seq_)[off];
			int q = (qual_)[off];
			// Is it legal for us to advance on characters other than 'c'?
			if(!(cons.mustMatch() && !overall.mustMatch()) || c == 4) {
				// There may be legal edits
				Boolean bail = false;
				if(!bloc.valid()) {
					// Range delimited by tloc/bloc has size 1
					long ntop = ltr ? topb : topf;
					bwops_++;
					int cc = ebwt.mapLF1(ntop, tloc);
					if(cc < 0) bail = true;
					else { t[cc] = ntop; b[cc] = ntop+1; }
				}
				if(!bail) {
					if((cons.canMismatch(q, sc_) && overall.canMismatch(q, sc_)) || c == 4) {
						Constraint oldCons = cons, oldOvCons = overall;
						SideLocus oldTloc = tloc, oldBloc = bloc;
						if(c != 4) {
							cons.chargeMismatch(q, sc_);
							overall.chargeMismatch(q, sc_);
						}
						// Can leave the zone as-is
						if(!leaveZone || (cons.acceptable() && overall.acceptable())) {
							for(int j = 0; j < 4; j++) {
								if(j == c || b[j] == t[j]) continue;
								// Potential mismatch
								nextLocsBi(tloc, bloc, tf[j], bf[j], tb[j], bb[j], i+1);
								int loff = off;
								if(!ltr) loff = (int)(s.steps.size() - loff - 1);
								Edit edit(off, j, c, EDIT_TYPE_MM, false);
								DoublyLinkedList<Edit> editl;
								editl.payload = edit;
								if(prevEdit != null) {
									prevEdit.next = editl;
									editl.prev = prevEdit;
								}
								bwedits_++;
								if(!searchSeedBi(
									i+1,     // depth into steps_[] array
									depth+1, // recursion depth
									tf[j],   // top in BWT
									bf[j],   // bot in BWT
									tb[j],   // top in BWT'
									bb[j],   // bot in BWT'
									tloc,    // locus for top (perhaps unititialized)
									bloc,    // locus for bot (perhaps unititialized)
									c0,      // constraints to enforce in seed zone 0
									c1,      // constraints to enforce in seed zone 1
									c2,      // constraints to enforce in seed zone 2
									overall, // overall constraints to enforce
									editl))  // latest edit
								{
									return false;
								}
								if(prevEdit != null) prevEdit.next = null;
							}
						} else {
							// Not enough edits to make this path
							// non-redundant with other seeds
						}
						cons = oldCons;
						overall = oldOvCons;
						tloc = oldTloc;
						bloc = oldBloc;
					}
				} // if(!bail)
			}
			if(c == 4) {
				return true; // couldn't handle the N
			}
			if(leaveZone && (!cons.acceptable() || !overall.acceptable())) {
				// Not enough edits to make this path non-redundant with
				// other seeds
				return true;
			}
			if(!bloc.valid()) {
				// Range delimited by tloc/bloc has size 1
				long top = ltr ? topb : topf;
				bwops_++;
				t[c] = ebwt.mapLF1(top, tloc, c);
				if(t[c] == IndexTypes.OFF_MASK) {
					return true;
				}
				b[c] = t[c]+1;
			}
			if(b[c] == t[c]) {
				return true;
			}
			topf = tf[c]; botf = bf[c];
			topb = tb[c]; botb = bb[c];
			if(i+1 == (int)s.steps.size()) {
				// Finished aligning seed
				if(!reportHit(topf, botf, topb, botb, seq_.length(), prevEdit)) {
					return false; // Memory exhausted
				}
				return true;
			}
			nextLocsBi(tloc, bloc, tf[c], bf[c], tb[c], bb[c], i+1);
		}
		return true;
	}
	
	public void nextLocsBi(
			SideLocus tloc,            // top locus
			SideLocus bloc,            // bot locus
			long topf,              // top in BWT
			long botf,              // bot in BWT
			long topb,              // top in BWT'
			long botb,              // bot in BWT'
			int step                    // step to get ready for
			) {
		if(step == (int)s_.steps.size()) return; // no more steps!
		// Which direction are we going in next?
		if(s_.steps[step] > 0) {
			// Left to right; use BWT'
			if(botb - topb == 1) {
				// Already down to 1 row; just init top locus
				tloc.initFromRow(topb, ebwtBw_.eh(), ebwtBw_.ebwt());
				bloc.invalidate();
			} else {
				SideLocus.initFromTopBot(
					topb, botb, ebwtBw_.eh(), ebwtBw_.ebwt(), tloc, bloc);
			}
		} else {
			// Right to left; use BWT
			if(botf - topf == 1) {
				// Already down to 1 row; just init top locus
				tloc.initFromRow(topf, ebwtFw_.eh(), ebwtFw_.ebwt());
				bloc.invalidate();
			} else {
				SideLocus.initFromTopBot(
					topf, botf, ebwtFw_.eh(), ebwtFw_.ebwt(), tloc, bloc);
			}
		}
	}
	
	public Boolean extendAndReportHit(
			long topf,                     // top in BWT
			long botf,                     // bot in BWT
			long topb,                     // top in BWT'
			long botb,                     // bot in BWT'
			int len,                      // length of hit
			DoublyLinkedList<Edit> prevEdit) {
		double nlex = 0, nrex = 0;
		long t[4], b[4];
		long tp[4], bp[4];
		SideLocus tloc, bloc;
		if(off_ > 0) {
			Ebwt *ebwt = ebwtFw_;
			// Extend left using forward index
			BTDnaString seq = fw_ ? read_.patFw : read_.patRc;
			// See what we get by extending 
			long top = topf, bot = botf;
			t[0] = t[1] = t[2] = t[3] = 0;
			b[0] = b[1] = b[2] = b[3] = 0;
			tp[0] = tp[1] = tp[2] = tp[3] = topb;
			bp[0] = bp[1] = bp[2] = bp[3] = botb;
			SideLocus tloc, bloc;
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
			for(double ii = off_; ii > 0; ii--) {
				double i = ii-1;
				// Get char from read
				int rdc = seq.get(i);
				// See what we get by extending 
				if(bloc.valid()) {
					bwops_++;
					t[0] = t[1] = t[2] = t[3] =
					b[0] = b[1] = b[2] = b[3] = 0;
					ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
					SANITY_CHECK_4TUP(t, b, tp, bp);
					int nonz = -1;
					Boolean abort = false;
					for(int j = 0; j < 4; j++) {
						if(b[i] > t[i]) {
							if(nonz >= 0) {
								abort = true;
								break;
							}
							nonz = j;
							top = t[i]; bot = b[i];
						}
					}
					if(abort || nonz != rdc) {
						break;
					}
				} else {
					bwops_++;
					int c = ebwt.mapLF1(top, tloc);
					if(c != rdc) {
						break;
					}
					bot = top + 1;
				}
				if(++nlex == 255) {
					break;
				}
				INIT_LOCS(top, bot, tloc, bloc, *ebwt);
			}
		}
		double rdlen = read_.length();
		double nright = rdlen - off_ - len;
		if(nright > 0 && ebwtBw_ != null) {
			Ebwt ebwt = ebwtBw_;
			// Extend right using backward index
			const BTDnaString& seq = fw_ ? read_.patFw : read_.patRc;
			// See what we get by extending 
			long top = topb, bot = botb;
			t[0] = t[1] = t[2] = t[3] = 0;
			b[0] = b[1] = b[2] = b[3] = 0;
			tp[0] = tp[1] = tp[2] = tp[3] = topb;
			bp[0] = bp[1] = bp[2] = bp[3] = botb;
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
			for(double i = off_ + len; i < rdlen; i++) {
				// Get char from read
				int rdc = seq.get(i);
				// See what we get by extending 
				if(bloc.valid()) {
					bwops_++;
					t[0] = t[1] = t[2] = t[3] =
					b[0] = b[1] = b[2] = b[3] = 0;
					ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
					SANITY_CHECK_4TUP(t, b, tp, bp);
					int nonz = -1;
					Boolean abort = false;
					for(int j = 0; j < 4; j++) {
						if(b[i] > t[i]) {
							if(nonz >= 0) {
								abort = true;
								break;
							}
							nonz = j;
							top = t[i]; bot = b[i];
						}
					}
					if(abort || nonz != rdc) {
						break;
					}
				} else {
					bwops_++;
					int c = ebwt.mapLF1(top, tloc);
					if(c != rdc) {
						break;
					}
					bot = top + 1;
				}
				if(++nrex == 255) {
					break;
				}
				INIT_LOCS(top, bot, tloc, bloc, *ebwt);
			}
		}
		return reportHit(topf, botf, topb, botb, len, prevEdit);
	}
	
	public Boolean reportHit(
			long topf,                     // top in BWT
			long botf,                     // bot in BWT
			long topb,                     // top in BWT'
			long botb,                     // bot in BWT'
			int len,                      // length of hit
			DoublyLinkedList<Edit> prevEdit)  // previous edit
	{
		// Add information about the seed hit to AlignmentCache.  This
		// information eventually makes its way back to the SeedResults
		// object when we call finishAlign(...).
		BTDnaString rf = tmprfdnastr_;
		rf.clear();
		edits_.clear();
		if(prevEdit != null) {
			prevEdit.toList(edits_);
			Edit.sort(edits_);
			Edit.toRef(seq_, edits_, rf);
		} else {
			rf = seq_;
		}
		// Sanity check: shouldn't add the same hit twice.  If this
		// happens, it may be because our zone Constraints are not set up
		// properly and erroneously return true from acceptable() when they
		// should return false in some cases.
		if(!ca_.addOnTheFly(rf, topf, botf, topb, botb)) {
			return false;
		}
		return true;
	}
}
