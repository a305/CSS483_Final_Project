package com.uwb.bt2j.indexer;

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
}

