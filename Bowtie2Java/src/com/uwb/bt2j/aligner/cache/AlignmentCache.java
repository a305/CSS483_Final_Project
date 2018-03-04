package com.uwb.bt2j.aligner.cache;

import com.uwb.bt2j.indexer.EList;
import com.uwb.bt2j.indexer.IndexTypes;
import com.uwb.bt2j.indexer.PList;
import com.uwb.bt2j.indexer.PListSlice;
import com.uwb.bt2j.indexer.Pool;
import com.uwb.bt2j.indexer.RedBlack;
import com.uwb.bt2j.indexer.RedBlackNode;

class AlignmentCache {
  public static final int CACHE_PAGE_SZ = 16 * 1024;
  
  protected Pool pool_;
  protected RedBlack<QKey, QVal> qmap_;
  protected PList<QKey> qlist_;
  protected RedBlack<QKey, QVal> samap_;
  protected PList<Long> salist_;
  protected boolean shared_;
  protected double version_;
  
  public AlignmentCache(long bytes, boolean shared) {
	  pool_ = new Pool(bytes, CACHE_PAGE_SZ);
	  qmap_ = new RedBlack(CACHE_PAGE_SZ,3);
	  qlist_ = new PList(3);
	  samap_ = new RedBlack(CACHE_PAGE_SZ,3);
	  salist_ = new PList(3);
	  shared_ = shared;
	  version_ = 0;
  }
  
  public void queryQVal(
		  QVal qv,
			EList<SATuple> satups,
			int nrange,
			int nelt,
			boolean getLock)
  {
	  if(shared_ && getLock) {
			queryQvalImpl(qv, satups, nrange, nelt);
		} else {
			queryQvalImpl(qv, satups, nrange, nelt);
		}
  }
  
  public final Boolean empty() {
	  boolean ret = qmap_.empty();
		return ret;
  }
  
  public QVal add(
		  QKey qk,
			boolean added,
			boolean getLock) {
	  if(shared_ && getLock) {
			return addImpl(qk, added);
		} else {
			return addImpl(qk, added);
		}
  }
  
  public boolean addOnTheFly(
		  QVal qv,         // qval that points to the range of reference substrings
			QKey sak, // the key holding the reference substring
			int topf,    // top range elt in BWT index
			int botf,    // bottom range elt in BWT index
			int topb,    // top range elt in BWT' index
			int botb,    // bottom range elt in BWT' index
			boolean getLock){
	  if(shared_ && getLock) {
			return addOnTheFlyImpl(qv, sak, topf, botf, topb, botb);
		} else {
			return addOnTheFlyImpl(qv, sak, topf, botf, topb, botb);
		}
  }
  
  public void clear() {
	  pool_.clear();
		qmap_.clear();
		qlist_.clear();
		samap_.clear();
		salist_.clear();
		version_++;
  }
  
  public double qNumKeys() {
	  return qmap_.size();
  }
  
  public double saNumKeys() {
	  return samap_.size();
  }
  
  public double qSize() {
	  return qlist_.size();
  }
  
  public double saSize() {
	  return salist_.size();
  }
  
  public Pool pool() {
	  return pool_;
  }
  
  public Boolean shared() {
	  return shared_;
  }
  
  public double version() {
	  return version_;
  }
  
  private void queryQvalImpl(
		  QVal qv,
			EList<SATuple> satups,
			int nrange,
			int nelt){
	  int refi = qv.offset();
	  long reff = refi + qv.numRanges();
		// For each reference sequence sufficiently similar to the
		// query sequence in the QKey...
		for(int i = refi; i < reff; i++) {
			// Get corresponding SAKey, containing similar reference
			// sequence & length
			QKey sak = qlist_.get(i);
			// Shouldn't have identical keys in qlist_
			// Get corresponding SANode
			RedBlackNode<QKey,SAVal> n = samap_.lookup(sak);
			SAVal sav = n.payload;
			if(sav.len > 0) {
				nrange++;
				satups.expand();
				satups.back().init(sak, sav.topf, sav.topb, new PListSlice(salist_, sav.i, sav.len));
				nelt += sav.len;
			}
		}
  }
  
  private boolean addOnTheFlyImpl(
		  	QVal qv,         // qval that points to the range of reference substrings
			QKey sak, // the key holding the reference substring
			int topf,    // top range elt in BWT index
			int botf,    // bottom range elt in BWT index
			int topb,    // top range elt in BWT' index
			int botb)    // bottom range elt in BWT' index
  {
	  boolean added = true;
		// If this is the first reference sequence we're associating with
		// the query sequence, initialize the QVal.
		if(!qv.valid()) {
			qv.init(qlist_.size(), 0, 0);
		}
		qv.addRange(botf-topf); // update tally for # ranges and # elts
		if(!qlist_.add(pool(), sak)) {
			return false; // Exhausted pool memory
		}

		RedBlackNode<QKey,SAVal> s = samap_.add(pool(), sak, added);
		if(s == null) {
			return false; // Exhausted pool memory
		}
		if(added) {
			s.payload.i = salist_.size();
			s.payload.len = (botf - topf);
			s.payload.topf = topf;
			s.payload.topb = topb;
			for(int j = 0; j < (botf-topf); j++) {
				if(!salist_.add(pool(), IndexTypes.OFF_MASK)) {
					// Change the payload's len field
					s.payload.len = j;
					return false; // Exhausted pool memory
				}
			}
		}
		// Now that we know all allocations have succeeded, we can do a few final
		// updates
		
		return true; 
  }
  
  private QVal addImpl(QKey qk, boolean added) {
	  RedBlackNode<QKey,QVal> n = qmap_.add(pool(), qk, added);
		return (n != null ? n.payload : null);
  }
}
