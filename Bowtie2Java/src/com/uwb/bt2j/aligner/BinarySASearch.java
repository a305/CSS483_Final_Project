package com.uwb.bt2j.aligner;

import com.uwb.bt2j.indexer.EList;

public class BinarySASearch<TStr,TSufElt> {
	public long binarySASearch(
			TStr host,
			int qry,
			EList<TSufElt> sa){
		int lLcp = 0, rLcp = 0; // greatest observed LCPs on left and right
		int l = 0, r = (int)sa.size()+1; // binary-search window
		int hostLen = (int)host.length();
		while(true) {
			int m = (l+r) >> 1;
			if(m == l) {
				// Binary-search window has closed: we have an answer
				if(m > 0 && sa.get(m-1) == qry) {
					return Long.MAX_VALUE; // qry matches
				}
				return m; // Return index of right-hand suffix
			}
			int suf = sa.get(m-1);
			if(suf == qry) {
				return Long.MAX_VALUE; // query matches an elt of sa
			}
			int lcp = Math.min(lLcp, rLcp);

			// Keep advancing lcp, but stop when query mismatches host or
			// when the counter falls off either the query or the suffix
			while(suf+lcp < hostLen && qry+lcp < hostLen && host[suf+lcp] == host[qry+lcp]) {
				lcp++;
			}
			// Fell off the end of either the query or the sa elt?
			boolean fell = (suf+lcp == hostLen || qry+lcp == hostLen);
			if((fell && qry+lcp == hostLen) || (!fell && host[suf+lcp] < host[qry+lcp])) {
				// Query is greater than sa elt
				l = m;                 // update left bound
				lLcp = Math.max(lLcp, lcp); // update left lcp
			}
			else if((fell && suf+lcp == hostLen) || (!fell && host[suf+lcp] > host[qry+lcp])) {
				// Query is less than sa elt
				r = m;                 // update right bound
				rLcp = Math.max(rLcp, lcp); // update right lcp
			}
		}
		return Long.MAX_VALUE;
	}
}
