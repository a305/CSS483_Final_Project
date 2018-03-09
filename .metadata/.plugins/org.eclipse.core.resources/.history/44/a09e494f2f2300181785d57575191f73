package com.uwb.bt2j.indexer;

import com.uwb.bt2j.indexer.types.EList;
import com.uwb.bt2j.indexer.types.S2bDnaString;

public class MultikeyQSort <TStr, TPos, THost, T> {
    public void swap(TStr s, int slen, TPos a, TPos b) {
        swap(s[a],s[b]);
    }

    public void vecswap(TStr s, int slen, TPos i, TPos j, TPos n, TPos begin, TPos end) {
        while(n-- > 0) {
		TPos a = i+n;
		TPos b = j+n;
		swap(s, slen, a, b);
	    }
    }

    public void vecswap2(TStr s, int slen, TStr s2, TPos i, TPos j, TPos n, TPos begin, TPos end) {
        while(n-- > 0) {
		TPos a = i+n;
		TPos b = j+n;
		swap(s, slen, a, b);
        swap(s2, slen, a, b);
	    }
    }

    public boolean assertPartitionedSuf(THost host, long[] s, int hi, int pivot, int begin, int end, int depth) {
        int hlen = host.length();
	    int state = 0; // 0 . 1st = section, 1 . < section, 2 . > section, 3 . 2nd = section
	for(int i = begin; i < end; i++) {
		switch(state) {
		case 0:
			if       (CHAR_AT_SUF(i, depth) < pivot)  { state = 1; break; }
			else if  (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
            break;
		case 1:
			if       (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
			else if  (CHAR_AT_SUF(i, depth) == pivot) { state = 3; break; }
            break;
		case 2:
			if       (CHAR_AT_SUF(i, depth) == pivot) { state = 3; break; }
            break;
		}
	}
	return true;
    }

    public boolean assertPartitionedSuf2(THost host, long[] s, int hi, int pivot, int begin, int end, int depth) {
        int hlen = host.length();
	    int state = 0; // 0 . 1st = section, 1 . < section, 2 . > section, 3 . 2nd = section
	for(int i = begin; i < end; i++) {
		switch(state) {
		case 0:
            if       (CHAR_AT_SUF(i, depth) == pivot) { state = 1; break; }
			else if  (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
		case 1:
			if       (CHAR_AT_SUF(i, depth) > pivot)  { state = 2; break; }
		}
	}
	return true;
    }

    public void mkeyQSortSuf(T host, int heln, long[] s, int slen, int hi, int begin, int end, int depth, int upto) {
        // Helper for making the recursive call; sanity-checks arguments to
	// make sure that the problem actually got smaller.
	#define MQS_RECURSE_SUF(nbegin, nend, ndepth) { \
		assert(nbegin > begin || nend < end || ndepth > depth); \
		if(ndepth < upto) { /* don't exceed depth of 'upto' */ \
			mkeyQSortSuf(host, hlen, s, slen, hi, nbegin, nend, ndepth, upto); \
		} \
	}

	int a, b, c, d, /*e,*/ r;
	int n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF); // pick pivot, swap it into [begin]
	int v = CHAR_AT_SUF(begin, depth); // v <- randomly-selected pivot value
	a = b = begin;
	c = d = end-1;
	while(true) {
		// Invariant: everything before a is = pivot, everything
		// between a and b is <
		int bc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v >= (bc = CHAR_AT_SUF(b, depth))) {
			if(v == bc) {
				SWAP(s, a, b); a++;
			}
			b++;
		}
		// Invariant: everything after d is = pivot, everything
		// between c and d is >
		int cc = 0; // shouldn't have to init but gcc on Mac complains
		while(b <= c && v <= (cc = CHAR_AT_SUF(c, depth))) {
			if(v == cc) {
				SWAP(s, c, d); d--;
			}
			c--;
		}
		if(b > c) break;
		SWAP(s, b, c);
		b++;
		c--;
	}
    r = Math.min(a-begin, b-a);
    VECSWAP(s, begin, b-r,   r);  // swap left = to center
	r = Math.min(d-c, end-d-1);
    VECSWAP(s, b,     end-r, r);  // swap right = to center
	r = b-a; // r <- # of <'s
	if(r > 0) {
		MQS_RECURSE_SUF(begin, begin + r, depth); // recurse on <'s
	}
	// Do not recurse on ='s if the pivot was the off-the-end value;
	// they're already fully sorted
	if(v != hi) {
		MQS_RECURSE_SUF(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
	}
	r = d-c; // r <- # of >'s excluding those exhausted
	if(r > 0 && v < hi-1) {
		MQS_RECURSE_SUF(end-r, end, depth); // recurse on >'s
	}
    }

    public void mkeyQSortSuf(T host, long[] s, int slen, int hi, boolean verbose, boolean sanityCheck, int upto) {
        int hlen = host.length();
	if(sanityCheck) sanityCheckInputSufs(s, slen);
	mkeyQSortSuf(host, hlen, s, slen, hi, (int)0, slen, (int)0, upto);
	if(sanityCheck) sanityCheckOrderedSufs(host, hlen, s, slen, upto);
    }

    public void mkeyQSortSuf2(T host, int hlen, long[] s, int slen, long[] s2, int hi, int _begin, int _end, int _depth, int upto, EList<Integer> boundaries) {
        ELList<QSortRange, 3, 1024> block_list;
    while(true) {
        int begin = 0, end = 0, depth = 0;
        if(block_list.size() == 0) {
            begin = _begin;
            end = _end;
            depth = _depth;
        } else {
            if(block_list.back().size() > 0) {
                begin = block_list.back()[0].begin;
                end = block_list.back()[0].end;
                depth = block_list.back()[0].depth;
                block_list.back().erase(0);
            } else {
                block_list.resize(block_list.size() - 1);
                if(block_list.size() == 0) {
                    break;
                }
            }
        }
        if(depth == upto) {
            if(boundaries != null) {
                boundaries.push_back(end);
            }
            continue;
        }
        int a, b, c, d, /*e,*/ r;
        int n = end - begin;
        if(n <= 1) { // 1-element list already sorted
            if(n == 1 && boundaries != null) {
                boundaries.push_back(end);
            }
            continue;
        }
        CHOOSE_AND_SWAP_PIVOT(SWAP2, CHAR_AT_SUF); // pick pivot, swap it into [begin]
        int v = CHAR_AT_SUF(begin, depth); // v <- randomly-selected pivot value

        a = b = begin;
        c = d = /*e =*/ end-1;
        while(true) {
            // Invariant: everything before a is = pivot, everything
            // between a and b is <
            int bc = 0; // shouldn't have to init but gcc on Mac complains
            while(b <= c && v >= (bc = CHAR_AT_SUF(b, depth))) {
                if(v == bc) {
                    SWAP2(s, s2, a, b); a++;
                }
                b++;
            }
            // Invariant: everything after d is = pivot, everything
            // between c and d is >
            int cc = 0; // shouldn't have to init but gcc on Mac complains
            while(b <= c && v <= (cc = CHAR_AT_SUF(c, depth))) {
                if(v == cc) {
                    SWAP2(s, s2, c, d); d--; /*e--;*/
                }
                //else if(c == e && v == hi) e--;
                c--;
            }
            if(b > c) break;
            SWAP2(s, s2, b, c);
            b++;
            c--;
        }
        r = Math.min(a-begin, b-a);
        VECSWAP2(s, s2, begin, b-r,   r);  // swap left = to center
        r = Math.min(d-c, end-d-1);
        VECSWAP2(s, s2, b,     end-r, r);  // swap right = to center
        r = b-a; // r <- # of <'s
        block_list.expand();
        block_list.back().clear();
        if(r > 0) { // recurse on <'s
            block_list.back().expand();
            block_list.back().back().begin = begin;
            block_list.back().back().end = begin + r;
            block_list.back().back().depth = depth;
        }
        // Do not recurse on ='s if the pivot was the off-the-end value;
        // they're already fully sorted
        //if(v != hi) { // recurse on ='s
            block_list.back().expand();
            block_list.back().back().begin = begin + r;
            block_list.back().back().end = begin + r + (a-begin) + (end-d-1);
            block_list.back().back().depth = depth + 1;
        //}
        r = d-c;   // r <- # of >'s excluding those exhausted
        if(r > 0 /*&& v < hi-1*/) { // recurse on >'s
            block_list.back().expand();
            block_list.back().back().begin = end - r;
            block_list.back().back().end = end;
            block_list.back().back().depth = depth;
        }
    }
    }

    public void mkeyQSortSuf2(T host, long[] s, int slen, long[] s2, int hi, boolean verbose, boolean sanityCheck, int upto, EList<Integer> boundaries) {
        int hlen = host.length();
    if(sanityCheck) sanityCheckInputSufs(s, slen);
    long[] sOrig = null;
    if(sanityCheck) {
        sOrig = new TIndexOffU[slen];
        memcpy(sOrig, s, IndexTypes.OFF_SIZE * slen);
    }
    mkeyQSortSuf2(host, hlen, s, slen, s2, hi, (int)0, slen, (int)0, upto, boundaries);
    if(sanityCheck) {
        sanityCheckOrderedSufs(host, hlen, s, slen, upto);
    }
    }

    public boolean sufDcLt(T1 host, T2 s1, T2 s2, DifferenceCoverSample<T1> dc, boolean sanityCheck) {
        int diff = dc.tieBreakOff(s1, s2);
	boolean ret = dc.breakTie(s1+diff, s2+diff) < 0;
	return ret;
    }

    public void qsortSufDc(T host, int hlen, long[] s, int slen, DifferenceCoverSample<T> dc, int begin, int end, boolean sanityCheck) {
        int n = end - begin;
	if(n <= 1) return;                 // 1-element list already sorted
	// Note: rand() didn't really cut it here; it seemed to run out of
	// randomness and, after a time, returned the same thing over and
	// over again
	int a = (rand() % n) + begin; // choose pivot between begin and end
	SWAP(s, end-1, a); // move pivot to end
	int cur = 0;
	for(int i = begin; i < end-1; i++) {
		if(sufDcLt(host, s[i], s[end-1], dc, sanityCheck)) {
			SWAP(s, i, begin + cur);
			cur++;
		}
	}
	// Put pivot into place
	SWAP(s, end-1, begin+cur);
	if(begin+cur > begin) qsortSufDc(host, hlen, s, slen, dc, begin, begin+cur);
	if(end > begin+cur+1) qsortSufDc(host, hlen, s, slen, dc, begin+cur+1, end);
    }

    public void mkeyQSortSufDcU8(T1 host1, T2 host, int hlen, long[] s, int slen, DifferenceCoverSample<T1> dc, int hi, boolean verbose, boolean sanityCheck) {
        if(sanityCheck) sanityCheckInputSufs(s, slen);
	mkeyQSortSufDcU8(host1, host, hlen, s, slen, dc, hi, 0, slen, 0, sanityCheck);
	if(sanityCheck) sanityCheckOrderedSufs(host1, hlen, s, slen, IndexTypes.OFF_MASK);
    }
    
    public boolean sufDcLtU8(T1 host1, T2 host, int hlen, int s1, int s2, DifferenceCoverSample<T1> dc, boolean sanityCheck) {
    	hlen += 0;
    	size_t diff = dc.tieBreakOff((TIndexOffU)s1, (TIndexOffU)s2);
    	assert_lt(diff, dc.v());
    	assert_lt(diff, hlen-s1);
    	assert_lt(diff, hlen-s2);
    	if(sanityCheck) {
    		for(size_t i = 0; i < diff; i++) {
    			assert_eq(host[s1+i], host1[s2+i]);
    		}
    	}
    	bool ret = dc.breakTie((TIndexOffU)(s1+diff), (TIndexOffU)(s2+diff)) < 0;
    	// Sanity-check return value using dollarLt
    #ifndef NDEBUG
    	bool ret2 = sstr_suf_lt(host1, s1, hlen, host, s2, hlen, false);
    	assert(!sanityCheck || ret == ret2);
    #endif
    	return ret;
    }
    
    public void qsortSufDcU8(T1 host1, T2 host, int hlen, long[] s, int slen, DifferenceCoverSample<T1> dc, int begin, int end, boolean sanityCheck) {
    	size_t n = end - begin;
    	if(n <= 1) return;                 // 1-element list already sorted
    	// Note: rand() didn't really cut it here; it seemed to run out of
    	// randomness and, after a time, returned the same thing over and
    	// over again
    	size_t a = (rand() % n) + begin; // choose pivot between begin and end
    	assert_lt(a, end);
    	assert_geq(a, begin);
    	SWAP(s, end-1, a); // move pivot to end
    	size_t cur = 0;
    	for(size_t i = begin; i < end-1; i++) {
    		if(sufDcLtU8(host1, host, hlen, s[i], s[end-1], dc, sanityCheck)) {
    #ifndef NDEBUG
    			if(sanityCheck) {
    				assert(sstr_suf_lt(host1, s[i], hlen, host1, s[end-1], hlen, false));
    			}
    			assert_lt(begin + cur, end-1);
    #endif
    			SWAP(s, i, begin + cur);
    			cur++;
    		}
    	}
    	// Put pivot into place
    	assert_lt(cur, end-begin);
    	SWAP(s, end-1, begin+cur);
    	if(begin+cur > begin) qsortSufDcU8(host1, host, hlen, s, slen, dc, begin, begin+cur);
    	if(end > begin+cur+1) qsortSufDcU8(host1, host, hlen, s, slen, dc, begin+cur+1, end);
    }
    
    public short get_uint8(TStr t, int off) {
    	return t[off];
    }
    
    public short get_uint8(S2bDnaString t, int off) {
    	return (short) t[off];
    }
    
    public int char_at_suf_u8(TStr host, int hlen, lonh[] s, int si, int off, short hi) {
    	return ((off+s[si]) < hlen) ? get_uint8(host, off+s[si]) : (hi);
    }
    
    public void selectionSortSufDcU8(T1 host1, T2 host, int hlen, long[] s, int slen, DifferenceCoverSample<T1> dc, int hi, int begin, int end, int depth, boolean sanityCheck) {
    	#define ASSERT_SUF_LT(l, r) \
    	if(sanityCheck && \
    	   !sstr_suf_lt(host1, s[l], hlen, host1, s[r], hlen, false)) { \
    		assert(false); \
    	}

    	assert_gt(end, begin+1);
    	assert_leq(end-begin, SELECTION_SORT_CUTOFF);
    	assert_eq(hi, 4);
    	size_t v = dc.v();
    	if(end == begin+2) {
    		size_t off = dc.tieBreakOff(s[begin], s[begin+1]);
    		if(off + s[begin] >= hlen ||
    		   off + s[begin+1] >= hlen)
    		{
    			off = OFF_MASK;
    		}
    		if(off != OFF_MASK) {
    			if(off < depth) {
    				qsortSufDcU8<T1,T2>(host1, host, hlen, s, slen, dc,
    				                    begin, end, sanityCheck);
    				// It's helpful for debugging if we call this here
    				if(sanityCheck) {
    					sanityCheckOrderedSufs(host1, hlen, s, slen,
    					                       OFF_MASK, begin, end);
    				}
    				return;
    			}
    			v = off - depth + 1;
    		}
    	}
    	assert_leq(v, dc.v());
    	size_t lim = v;
    	assert_geq(lim, 0);
    	for(size_t i = begin; i < end-1; i++) {
    		size_t targ = i;
    		size_t targoff = depth + s[i];
    		for(size_t j = i+1; j < end; j++) {
    			assert_neq(j, targ);
    			size_t joff = depth + s[j];
    			size_t k;
    			for(k = 0; k <= lim; k++) {
    				assert_neq(j, targ);
    				uint8_t jc = (k + joff < hlen)    ? get_uint8(host, k + joff)    : hi;
    				uint8_t tc = (k + targoff < hlen) ? get_uint8(host, k + targoff) : hi;
    				assert(jc != hi || tc != hi);
    				if(jc > tc) {
    					// the jth suffix is greater than the current
    					// smallest suffix
    					ASSERT_SUF_LT(targ, j);
    					break;
    				} else if(jc < tc) {
    					// the jth suffix is less than the current smallest
    					// suffix, so update smallest to be j
    					ASSERT_SUF_LT(j, targ);
    					targ = j;
    					targoff = joff;
    					break;
    				} else if(k == lim) {
    					// Check whether either string ends immediately
    					// after this character
    					assert_leq(k + joff + 1, hlen);
    					assert_leq(k + targoff + 1, hlen);
    					if(k + joff + 1 == hlen) {
    						// targ < j
    						assert_neq(k + targoff + 1, hlen);
    						ASSERT_SUF_LT(targ, j);
    						break;
    					} else if(k + targoff + 1 == hlen) {
    						// j < targ
    						ASSERT_SUF_LT(j, targ);
    						targ = j;
    						targoff = joff;
    						break;
    					}
    				} else {
    					// They're equal so far, keep going
    				}
    			}
    			// The jth suffix was equal to the current smallest suffix
    			// up to the difference-cover period, so disambiguate with
    			// difference cover
    			if(k == lim+1) {
    				assert_neq(j, targ);
    				if(sufDcLtU8(host1, host, hlen, s[j], s[targ], dc, sanityCheck)) {
    					// j < targ
    					assert(!sufDcLtU8(host1, host, hlen, s[targ], s[j], dc, sanityCheck));
    					ASSERT_SUF_LT(j, targ);
    					targ = j;
    					targoff = joff;
    				} else {
    					assert(sufDcLtU8(host1, host, hlen, s[targ], s[j], dc, sanityCheck));
    					ASSERT_SUF_LT(targ, j); // !
    				}
    			}
    		}
    		if(i != targ) {
    			ASSERT_SUF_LT(targ, i);
    			// swap i and targ
    			TIndexOffU tmp = s[i];
    			s[i] = s[targ];
    			s[targ] = tmp;
    		}
    		for(size_t j = i+1; j < end; j++) {
    			ASSERT_SUF_LT(i, j);
    		}
    	}
    	if(sanityCheck) {
    		sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
    	}
    }
    
    public void bucketSortSufDcU8(T1 host1, T2 host, int hlen, long[] s, int slen, DifferenceCoverSample<T1> dc, int hi, int _begin, int _end, int _depth, boolean sanityCheck) {
    	// 5 64-element buckets for bucket-sorting A, C, G, T, $
        TIndexOffU* bkts[4];
        for(size_t i = 0; i < 4; i++) {
            bkts[i] = new TIndexOffU[4 * 1024 * 1024];
        }
        ELList<size_t, 5, 1024> block_list;
        bool first = true;
        while(true) {
            size_t begin = 0, end = 0;
            if(first) {
                begin = _begin;
                end = _end;
                first = false;
            } else {
                if(block_list.size() == 0) {
                    break;
                }
                if(block_list.back().size() > 1) {
                    end = block_list.back().back(); block_list.back().pop_back();
                    begin = block_list.back().back();
                } else {
                    block_list.resize(block_list.size() - 1);
                    if(block_list.size() == 0) {
                        break;
                    }
                }
            }
            size_t depth = block_list.size() + _depth;
            assert_leq(end-begin, BUCKET_SORT_CUTOFF);
            assert_eq(hi, 4);
            if(end <= begin + 1) { // 1-element list already sorted
                continue;
            }
            if(depth > dc.v()) {
                // Quicksort the remaining suffixes using difference cover
                // for constant-time comparisons; this is O(k*log(k)) where
                // k=(end-begin)
                qsortSufDcU8<T1,T2>(host1, host, hlen, s, slen, dc, begin, end, sanityCheck);
                continue;
            }
            if(end-begin <= SELECTION_SORT_CUTOFF) {
                // Bucket sort remaining items
                selectionSortSufDcU8(host1, host, hlen, s, slen, dc, hi,
                                     begin, end, depth, sanityCheck);
                if(sanityCheck) {
                    sanityCheckOrderedSufs(host1, hlen, s, slen,
                                           OFF_MASK, begin, end);
                }
                continue;
            }
            size_t cnts[] = { 0, 0, 0, 0, 0 };
            for(size_t i = begin; i < end; i++) {
                size_t off = depth + s[i];
                uint8_t c = (off < hlen) ? get_uint8(host, off) : hi;
                assert_leq(c, 4);
                if(c == 0) {
                    s[begin + cnts[0]++] = s[i];
                } else {
                    bkts[c-1][cnts[c]++] = s[i];
                }
            }
            assert_eq(cnts[0] + cnts[1] + cnts[2] + cnts[3] + cnts[4], end - begin);
            size_t cur = begin + cnts[0];
            if(cnts[1] > 0) { memcpy(&s[cur], bkts[0], cnts[1] << (OFF_SIZE/4 + 1)); cur += cnts[1]; }
            if(cnts[2] > 0) { memcpy(&s[cur], bkts[1], cnts[2] << (OFF_SIZE/4 + 1)); cur += cnts[2]; }
            if(cnts[3] > 0) { memcpy(&s[cur], bkts[2], cnts[3] << (OFF_SIZE/4 + 1)); cur += cnts[3]; }
            if(cnts[4] > 0) { memcpy(&s[cur], bkts[3], cnts[4] << (OFF_SIZE/4 + 1)); }
            // This frame is now totally finished with bkts[][], so recursive
            // callees can safely clobber it; we're not done with cnts[], but
            // that's local to the stack frame.
            block_list.expand();
            block_list.back().clear();
            block_list.back().push_back(begin);
            for(size_t i = 0; i < 4; i++) {
                if(cnts[i] > 0) {
                    block_list.back().push_back(block_list.back().back() + cnts[i]);
                }
            }
        }
        // Done
        
        for(size_t i = 0; i < 4; i++) {
            delete [] bkts[i];
        }
    }
    
    public void mkeyQSortSufDcU8(T1 host1, T2 host, int hlen, long[] s, int slen, DifferenceCoverSample<T1> dc, int hi, int _begin, int _end, int _depth, boolean sanityCheck) {
    	// Helper for making the recursive call; sanity-checks arguments to
    	// make sure that the problem actually got smaller.
    	#define MQS_RECURSE_SUF_DC_U8(nbegin, nend, ndepth) { \
    		assert(nbegin > begin || nend < end || ndepth > depth); \
    		mkeyQSortSufDcU8(host1, host, hlen, s, slen, dc, hi, nbegin, nend, ndepth, sanityCheck); \
    	}
    	assert_leq(begin, slen);
    	assert_leq(end, slen);
    	size_t n = end - begin;
    	if(n <= 1) return; // 1-element list already sorted
    	if(depth > dc.v()) {
    		// Quicksort the remaining suffixes using difference cover
    		// for constant-time comparisons; this is O(k*log(k)) where
    		// k=(end-begin)
    		qsortSufDcU8<T1,T2>(host1, host, hlen, s, slen, dc, begin, end, sanityCheck);
    		if(sanityCheck) {
    			sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
    		}
    		return;
    	}
    	if(n <= BUCKET_SORT_CUTOFF) {
    		// Bucket sort remaining items
    		bucketSortSufDcU8(host1, host, hlen, s, slen, dc,
    		                  (uint8_t)hi, begin, end, depth, sanityCheck);
    		if(sanityCheck) {
    			sanityCheckOrderedSufs(host1, hlen, s, slen, OFF_MASK, begin, end);
    		}
    		return;
    	}
    	size_t a, b, c, d, r;
    	CHOOSE_AND_SWAP_PIVOT(SWAP1, CHAR_AT_SUF_U8); // choose pivot, swap to begin
    	int v = CHAR_AT_SUF_U8(begin, depth); // v <- pivot value
    	#ifndef NDEBUG
    	{
    		bool stillInBounds = false;
    		for(size_t i = begin; i < end; i++) {
    			if(depth < (hlen-s[i])) {
    				stillInBounds = true;
    				break;
    			} else { /* already fell off this suffix */ }
    		}
    		assert(stillInBounds); // >=1 suffix must still be in bounds
    	}
    	#endif
    	a = b = begin;
    	c = d = end-1;
    	while(true) {
    		// Invariant: everything before a is = pivot, everything
    		// between a and b is <
    		int bc = 0; // shouldn't have to init but gcc on Mac complains
    		while(b <= c && v >= (bc = CHAR_AT_SUF_U8(b, depth))) {
    			if(v == bc) {
    				SWAP(s, a, b); a++;
    			}
    			b++;
    		}
    		// Invariant: everything after d is = pivot, everything
    		// between c and d is >
    		int cc = 0; // shouldn't have to init but gcc on Mac complains
    		//bool hiLatch = true;
    		while(b <= c && v <= (cc = CHAR_AT_SUF_U8(c, depth))) {
    			if(v == cc) {
    				SWAP(s, c, d); d--;
    			}
    			//else if(hiLatch && cc == hi) { }
    			c--;
    		}
    		if(b > c) break;
    		SWAP(s, b, c);
    		b++;
    		c--;
    	}
    	assert(a > begin || c < end-1);                      // there was at least one =s
    	assert_lt(d-c, n); // they can't all have been > pivot
    	assert_lt(b-a, n); // they can't all have been < pivot
    	r = min(a-begin, b-a); VECSWAP(s, begin, b-r,   r);  // swap left = to center
    	r = min(d-c, end-d-1); VECSWAP(s, b,     end-r, r);  // swap right = to center
    	r = b-a; // r <- # of <'s
    	if(r > 0) {
    		MQS_RECURSE_SUF_DC_U8(begin, begin + r, depth); // recurse on <'s
    	}
    	// Do not recurse on ='s if the pivot was the off-the-end value;
    	// they're already fully sorted
    	if(v != hi) {
    		MQS_RECURSE_SUF_DC_U8(begin + r, begin + r + (a-begin) + (end-d-1), depth+1); // recurse on ='s
    	}
    	r = d-c; // r <- # of >'s excluding those exhausted
    	if(r > 0 && v < hi-1) {
    		MQS_RECURSE_SUF_DC_U8(end-r, end, depth); // recurse on >'s
    	}
    }
}

