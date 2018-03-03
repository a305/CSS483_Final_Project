package com.uwb.bt2j.aligner;

/* #include <iostream>
#include "ds.h"
#include "read.h"
#include "threading.h"
#include "mem_ids.h"
#include "simple_func.h"
#include "btypes.h" */

long CACHE_PAGE_SZ = (16 * 1024)
//typedef PLisPListSlice<TIndexOffU, CACHE_PAGE_SZ> TSlice;

class QKey() {
	public long seq; // sequence
	public int len; // length of sequence
	
	/**
	 * Initialize invalid QKey.
	 */
	public QKey() { reset(); }

	/**
	 * Initialize QKey from DNA string.
	 */
	public QKey(BTDnaString s, BTDnaString tmp) {
		init(s, tmp);
	}
	
	/**
	 * Initialize QKey from DNA string.  Rightmost character is placed in the
	 * least significant bitpair.
	 */
	public Boolean init(BTDnaString s, BTDnaString tmp)
	{
		seq = 0;
		len = (int)s.length();
		tmp.clear();
		if(len > 32) {
			len = 0xffffffff;
			return false; // wasn't cacheable
		} else {
			// Rightmost char of 's' goes in the least significant bitpair
			for(long i = 0; i < 32 && i < s.length(); i++) {
				int c = (int)s.get(i);
				//assert_range(0, 4, c);
				if(c == 4) {
					len = 0xffffffff;
					return false;
				}
				seq = (seq << 2) | s.get(i);
			}
			tmp.toString();
			assert(sstr_eq(tmp, s));
			//assert_leq(len, 32);
			return true; // was cacheable
		}
	}
	
	/**
	 * Return true if the read substring is cacheable.
	 */
	public Boolean cacheable() { return len != 0xffffffff; }
	
	/**
	 * Reset to uninitialized state.
	 */
	public void reset() { seq = 0; len = 0xffffffff; }

	/**
	 * True . my key is less than the given key.
	 */
	public String Output(QKey o) {
		return seq.toString();
	}

	/**
	 * True . my key is greater than the given key.
	 */
	public Boolean operator>(QKey o) {
		return !(*this < o || *this == o);
	}

	/**
	 * True . my key is equal to the given key.
	 */
	public Boolean Equals(QKey o) {
		return seq == o.seq && len == o.len;
	}


	/**
	 * True . my key is not equal to the given key.
	 */
	public Boolean NotEquals(QKey o) {
		return !(*this == o);
	}
	
////#ifndef NDEBUG
	/**
	 * Check that this is a valid, initialized QKey.
	 */
	Boolean repOk()  {
		return len != 0xffffffff;
	}
////#endif
}

/**
 * Payload for the query multimap: a range of elements in the reference
 * string list.
 */
class QVal {
	
	protected TIndexOffU i_;      // idx of first elt in qlist
	protected TIndexOffU rangen_; // # ranges (= # associated reference substrings)
	protected TIndexOffU eltn_;   // # elements (total)

	public QVal() { reset(); }

	/**
	 * Return the offset of the first reference substring in the qlist.
	 */
	public TIndexOffU offset()  { return i_; }

	/**
	 * Return the number of reference substrings associated with a read
	 * substring.
	 */
	public TIndexOffU numRanges()  {
		assert(valid());
		return rangen_;
	}

	/**
	 * Return the number of elements associated with all associated
	 * reference substrings.
	 */
	public TIndexOffU numElts()  {
		assert(valid());
		return eltn_;
	}
	
	/**
	 * Return true iff the read substring is not associated with any
	 * reference substrings.
	 */
	public Boolean empty()  {
		assert(valid());
		return numRanges() == 0;
	}

	/**
	 * Return true if the QVal is valid.
	 */
	public Boolean valid()  { return rangen_ != OFF_MASK; }
	
	/**
	 * Reset to invalid state.
	 */
	public void reset() {
		i_ = 0; rangen_ = eltn_ = OFF_MASK;
	}
	
	/**
	 * Initialize Qval.
	 */
	public void init(TIndexOffU i, TIndexOffU ranges, TIndexOffU elts) {
		i_ = i; rangen_ = ranges; eltn_ = elts;
	}
	
	/**
	 * Tally another range with given number of elements.
	 */
	public void addRange(TIndexOffU numElts) {
		rangen_++;
		eltn_ += numElts;
	}
	
////#ifndef NDEBUG
	/**
	 * Check that this QVal is internally consistent and consistent
	 * with the contents of the given cache.
	 */
	public Boolean repOk( AlignmentCache ac){
		if(rangen_ > 0) {
			assert(i_, ac.qSize());
			assert(i_ + rangen_, ac.qSize());
		}
		assert(eltn_, rangen_);
		return true;
	}
////#endif
}

/**
 * Key for the suffix array multimap: the reference substring and its
 * length.  Same as QKey so I typedef it.
 */
//typedef QKey SAKey;

/**
 * Payload for the suffix array multimap: (a) the top element of the
 * range in BWT, (b) the offset of the first elt in the salist, (c)
 * length of the range.
 */
class SAVal {
	
	public TIndexOffU topf;  // top in BWT
	public TIndexOffU topb;  // top in BWT'
	public TIndexOffU i;     // idx of first elt in salist
	public TIndexOffU len;   // length of range

	SAVal() : topf(), topb(), i(), len(OFF_MASK) { }

	/**
	 * Return true iff the SAVal is valid.
	 */
	public Boolean valid() { return len != OFF_MASK; }

////#ifndef NDEBUG
	/**
	 * Check that this SAVal is internally consistent and consistent
	 * with the contents of the given cache.
	 */
	public Boolean repOk( AlignmentCache ac){
		assert(len == 0 || i < ac.saSize());
		assert(i + len, ac.saSize());
		return true;
	}
////#endif
	
	/**
	 * Initialize the SAVal.
	 */
	public void init(
		TIndexOffU tf,
		TIndexOffU tb,
		TIndexOffU ii,
		TIndexOffU ln)
	{
		topf = tf;
		topb = tb;
		i = ii;
		len = ln;
	}
}

/**
 * One data structure that encapsulates all of the cached information
 * associated with a particular reference substring.  This is useful
 * for summarizing what info should be added to the cache for a partial
 * alignment.
 */
class SATuple {
	
	// bot/length of SA range equals offs.size()
	public QKey    key;  // sequence key
	public TIndexOffU topf;  // top in BWT index
	public TIndexOffU topb;  // top in BWT' index
	public PListSlice<TIndexOffU, CACHE_PAGE_SZ>   offs; // offsets

	public SATuple() { reset(); };

	public SATuple(QKey k, TIndexOffU tf, TIndexOffU tb, PListSlice<TIndexOffU, CACHE_PAGE_SZ> o) {
		init(k, tf, tb, o);
	}
	
	public void init(QKey k, TIndexOffU tf, TIndexOffU tb, PListSlice<TIndexOffU, CACHE_PAGE_SZ> o) {
		key = k; topf = tf; topb = tb; offs = o;
	}

	/**
	 * Initialize this SATuple from a subrange of the SATuple 'src'.
	 */
	public void init( SATuple src, long first, long last) {
		//assert_neq(OFF_MASK, src.topb);
		key = src.key;
		topf = (TIndexOffU)(src.topf + first);
		topb = OFF_MASK; // unknown!
		offs.init(src.offs, first, last);
	}
	
////#ifndef NDEBUG
	/**
	 * Check that this SATuple is internally consistent and that its
	 * PLisPListSlice<TIndexOffU, CACHE_PAGE_SZ> is consistent with its backing PList.
	 */
	public Boolean repOk()  {
		assert(offs.repOk());
		return true;
	}
////#endif

	/**
	 * Function for ordering SATuples.  This is used when prioritizing which to
	 * explore first when extending seed hits into full alignments.  Smaller
	 * ranges get higher priority and we use 'top' to break ties, though any
	 * way of breaking a tie would be fine.
	 */
	public Boolean GreaterThan( SATuple o)  {
		if(offs.size() < o.offs.size()) {
			return true;
		}
		if(offs.size() > o.offs.size()) {
			return false;
		}
		return topf < o.topf;
	}
	
	public Boolean LessThan( SATuple o)  {
		if(offs.size() < o.offs.size()) {
			return false;
		}
		if(offs.size() > o.offs.size()) {
			return true;
		}
		return topf > o.topf;
	}
	
	public Boolean EqualTo( SATuple o)  {
		return key == o.key && topf == o.topf && topb == o.topb && offs == o.offs;
	}

	public void reset() { topf = topb = OFF_MASK; offs.reset(); }
	
	/**
	 * Set the length to be at most the original length.
	 */
	public void setLength(long nlen) {
		//assert_leq(nlen, offs.size());
		offs.setLength(nlen);
	}
	
	/**
	 * Return the number of times this reference substring occurs in the
	 * reference, which is also the size of the 'offs' PListSlice<TIndexOffU, CACHE_PAGE_SZ>.
	 */
	public long size()  { return offs.size(); }
}

/**
 * Encapsulate the data structures and routines that itute a
 * particular cache, i.e., a particular stratum of the cache system,
 * which might comprise many strata.
 *
 * Each thread has a "current-read" AlignmentCache which is used to
 * build and store subproblem results as alignment is performed.  When
 * we're finished with a read, we might copy the cached results for
 * that read (and perhaps a bundle of other recently-aligned reads) to
 * a higher-level "across-read" cache.  Higher-level caches may or may
 * not be shared among threads.
 *
 * A cache consists chiefly of two multimaps, each implemented as a
 * Red-Black tree map backed by an EList.  A 'version' counter is
 * incremented every time the cache is cleared.
 */
class AlignmentCache<S> {
	
	protected Pool                   pool_;   // dispenses memory pages
	protected RedBlack<QKey, QVal>   qmap_;   // map from query substrings to reference substrings
	protected PList<QKey, CACHE_PAGE_SZ> qlist_;  // list of reference substrings
	protected RedBlack<QKey, SAVal> samap_;  // map from reference substrings to SA ranges
	protected PList<TIndexOffU, CACHE_PAGE_SZ> salist_; // list of SA ranges
	
	protected Boolean     shared_;  // true . this cache is global
	protected MUTEX_T mutex_m;    // mutex used for syncronization in case the the cache is shared.
	protected int version_; // cache version

	/* typedef RedBlackNode<QKey,  QVal>  QNode;
	typedef RedBlackNode<QKey, SAVal> SANode;

	typedef PList<QKey, CACHE_PAGE_SZ> TQList;
	typedef PList<TIndexOffU, CACHE_PAGE_SZ> TSAList; */

	public AlignmentCache(
		long bytes,
		Boolean shared) :
		pool_(bytes, CACHE_PAGE_SZ, CA_CAT),
		qmap_(CACHE_PAGE_SZ, CA_CAT),
		qlist_(CA_CAT),
		samap_(CACHE_PAGE_SZ, CA_CAT),
		salist_(CA_CAT),
		shared_(shared),
		mutex_m(),
		version_(0) { }

	/**
	 * Given a QVal, populate the given EList of SATuples with records
	 * describing all of the cached information about the QVal's
	 * reference substrings.
	 */
	//template <int S>
	public void queryQval(
		 QVal qv,
		EList<SATuple, S> satups,
		long nrange,
		long nelt,
		Boolean getLock = true)
	{
		if(shared_ && getLock) {
			ThreadSafe ts(mutex_m);
			queryQvalImpl(qv, satups, nrange, nelt);
		} else {
			queryQvalImpl(qv, satups, nrange, nelt);
		}
	}

	/**
	 * Return true iff the cache has no entries in it.
	 */
	public Boolean empty()  {
		Boolean ret = qmap_.empty();
		assert(!ret || qlist_.empty());
		assert(!ret || samap_.empty());
		assert(!ret || salist_.empty());
		return ret;
	}
	
	/**
	 * Add a new query key ('qk'), usually a 2-bit encoded substring of
	 * the read) as the key in a new Red-Black node in the qmap and
	 * return a pointer to the node's QVal.
	 *
	 * The expectation is that the caller is about to set about finding
	 * associated reference substrings, and that there will be future
	 * calls to addOnTheFly to add associations to reference substrings
	 * found.
	 */
	public QVal add(
		 QKey qk,
		Boolean added)
	{
		add(qk, added, true);
	}
	
	public QVal add(
		 QKey qk,
		Boolean added,
		Boolean getLock)
	{
		if(shared_ && getLock) {
			ThreadSafe ts(mutex_m);
			return addImpl(qk, added);
		} else {
			return addImpl(qk, added);
		}
	}
	
	/**
	 * Add a new association between a read sequnce ('seq') and a
	 * reference sequence ('')
	 */
	public Boolean addOnTheFly(
		QVal qv,         // qval that points to the range of reference substrings
		 QKey sak, // the key holding the reference substring
		TIndexOffU topf,    // top range elt in BWT index
		TIndexOffU botf,    // bottom range elt in BWT index
		TIndexOffU topb,    // top range elt in BWT' index
		TIndexOffU botb,    // bottom range elt in BWT' index
		Boolean getLock = true)
	{
		if(shared_ && getLock) {
			ThreadSafe ts(mutex_m);
			return addOnTheFlyImpl(qv, sak, topf, botf, topb, botb);
		} else {
			return addOnTheFlyImpl(qv, sak, topf, botf, topb, botb);
		}
	}

	/**
	 * Clear the cache, i.e. turn it over.  All HitGens referring to
	 * ranges in this cache will become invalid and the corresponding
	 * reads will have to be re-aligned.
	 */
	public void clear() {
		ThreadSafe ts(mutex_m);
		pool_.clear();
		qmap_.clear();
		qlist_.clear();
		samap_.clear();
		salist_.clear();
		version_++;
	}

	/**
	 * Return the number of keys in the query multimap.
	 */
	public long qNumKeys()  { return qmap_.size(); }

	/**
	 * Return the number of keys in the suffix array multimap.
	 */
	public long saNumKeys()  { return samap_.size(); }

	/**
	 * Return the number of elements in the reference substring list.
	 */
	public long qSize()  { return qlist_.size(); }

	/**
	 * Return the number of elements in the SA range list.
	 */
	public long saSize()  { return salist_.size(); }

	/**
	 * Return the pool.
	 */
	public Pool pool() { return pool_; }
	
	/**
	 * Return the lock object.
	 */
	public MUTEX_T lock() {
		return mutex_m;
	}

	/**
	 * Return true iff this cache is shared among threads.
	 */
	public Boolean shared()  { return shared_; }
	
	/**
	 * Return the current "version" of the cache, i.e. the total number
	 * of times it has turned over since its creation.
	 */
	public int version()  { return version_; }

	//template <int S>
	private void queryQvalImpl(
		 QVal qv,
		EList<SATuple, S> satups,
		long nrange,
		long nelt)
	{
		assert(qv.repOk(this));
		 long refi = qv.offset();
		 long reff = refi + qv.numRanges();
		// For each reference sequence sufficiently similar to the
		// query sequence in the QKey...
		for(long i = refi; i < reff; i++) {
			// Get corresponding QKey, containing similar reference
			// sequence & length
			QKey sak = qlist_.get(i);
			// Shouldn't have identical keys in qlist_
			assert(i == refi || qlist_.get(i) != qlist_.get(i-1));
			// Get corresponding SANode
			RedBlackNode<QKey, SAVal> n = samap_.lookup(sak);
			assert(n != null);
			 SAVal sav = n.payload;
			assert(sav.repOk(this));
			if(sav.len > 0) {
				nrange++;
				satups.expand();
				satups.back().init(sak, sav.topf, sav.topb, PListSlice<TIndexOffU, CACHE_PAGE_SZ>(salist_, sav.i, sav.len));
				nelt += sav.len;
////#ifndef NDEBUG
				// Shouldn't add consecutive identical entries too satups
				if(i > refi) {
					 SATuple b1 = satups.back();
					 SATuple b2 = satups[satups.size()-2];
					assert(b1.key != b2.key || b1.topf != b2.topf || b1.offs != b2.offs);
				}
////#endif
			}
		}
	}
	
	/**
	 * Add a new association between a read sequnce ('seq') and a
	 * reference sequence ('')
	 */
	private Boolean addOnTheFlyImpl(
		QVal qv,         // qval that points to the range of reference substrings
		 QKey sak, // the key holding the reference substring
		TIndexOffU topf,    // top range elt in BWT index
		TIndexOffU botf,    // bottom range elt in BWT index
		TIndexOffU topb,    // top range elt in BWT' index
		TIndexOffU botb)
	{
		bool added = true;
		// If this is the first reference sequence we're associating with
		// the query sequence, initialize the QVal.
		if(!qv.valid()) {
			qv.init((uint32_t)qlist_.size(), 0, 0);
		}
		qv.addRange(botf-topf); // update tally for # ranges and # elts
		if(!qlist_.add(pool(), sak)) {
			return false; // Exhausted pool memory
		}
		
		for(size_t i = qv.offset(); i < qlist_.size(); i++) {
			if(i > qv.offset()) {
				assert(qlist_.get(i) != qlist_.get(i-1));
			}
		}
	
		assert_eq(qv.offset() + qv.numRanges(), qlist_.size());
		RedBlackNode<QKey, SAVal> s = samap_.add(pool(), sak, &added);
		if(s == NULL) {
			return false; // Exhausted pool memory
		}
		assert(s->key.repOk());
		if(added) {
			s->payload.i = (TIndexOffU)salist_.size();
			s->payload.len = botf - topf;
			s->payload.topf = topf;
			s->payload.topb = topb;
			for(size_t j = 0; j < (botf-topf); j++) {
				if(!salist_.add(pool(), OFF_MASK)) {
					// Change the payload's len field
					s->payload.len = (TIndexOffU)j;
					return false; // Exhausted pool memory
				}
			}
			assert(s->payload.repOk(*this));
		}
		// Now that we know all allocations have succeeded, we can do a few final
		// updates
		
		return true; 
	}   // bottom range elt in BWT' index

	/**
	 * Add a new query key ('qk'), usually a 2-bit encoded substring of
	 * the read) as the key in a new Red-Black node in the qmap and
	 * return a pointer to the node's QVal.
	 *
	 * The expectation is that the caller is about to set about finding
	 * associated reference substrings, and that there will be future
	 * calls to addOnTheFly to add associations to reference substrings
	 * found.
	 */
	private QVal addImpl(
		 QKey qk,
		Boolean added)
	{
		assert(qk.cacheable());
		RedBlackNode<QKey,  QVal> n = qmap_.add(pool(), qk, added);
		return (n != null ? n.payload : null);
	}
}

/**
 * Interface used to query and update a pair of caches: one thread-
 * local and unsynchronized, another shared and synchronized.  One or
 * both can be null.
 */
class AlignmentCacheIface<S> {
	
	protected QKey qk_;  // key representation for current read substring
	protected QVal qv_; // pointer to value representation for current read substring
	protected QVal qvbuf_; // buffer for when key is uncacheable but we need a qv
	protected Boolean cacheable_; // true iff the read substring currently being aligned is cacheable
	
	protected long rangen_; // number of ranges since last alignment job began
	protected long eltsn_;  // number of elements since last alignment job began

	protected AlignmentCache current_; // cache dedicated to the current read
	protected AlignmentCache local_;   // local, unsynchronized cache
	protected AlignmentCache shared_;  // shared, synchronized cache
	public BTDnaString tmpdnastr_;

	public AlignmentCacheIface(
		AlignmentCache current,
		AlignmentCache local,
		AlignmentCache shared) :
		qk_(),
		qv_(null),
		cacheable_(false),
		rangen_(0),
		eltsn_(0),
		current_(current),
		local_(local),
		shared_(shared)
	{
		assert(current_ != null);
	}

//#if 0
	/**
	 * Query the relevant set of caches, looking for a QVal to go with
	 * the provided QKey.  If the QVal is found in a cache other than
	 * the current-read cache, it is copied into the current-read cache
	 * first and the QVal pointer for the current-read cache is
	 * returned.  This function never returns a pointer from any cache
	 * other than the current-read cache.  If the QVal could not be
	 * found in any cache OR if the QVal was found in a cache other
	 * than the current-read cache but could not be copied into the
	 * current-read cache, null is returned.
	 */
	 public QVal queryCopy( QKey qk) {
		 queryCopy(qk, true);
	 }
	public QVal queryCopy( QKey qk, Boolean getLock) {
		assert(qk.cacheable());
		AlignmentCache caches[3] = { current_, local_, shared_ };
		for(int i = 0; i < 3; i++) {
			if(caches[i] == null) continue;
			QVal qv = caches[i].query(qk, getLock);
			if(qv != null) {
				if(i == 0) return qv;
				if(!current_.copy(qk, qv, caches[i], getLock)) {
					// Exhausted memory in the current cache while
					// attempting to copy in the qk
					return null;
				}
				QVal curqv = current_.query(qk, getLock);
				assert(curqv != null);
				return curqv;
			}
		}
		return null;
	}

	/**
	 * Query the relevant set of caches, looking for a QVal to go with
	 * the provided QKey.  If a QVal is found and which is non-null,
	 * *which is set to 0 if the qval was found in the current-read
	 * cache, 1 if it was found in the local across-read cache, and 2
	 * if it was found in the shared across-read cache.
	 */
	// Was inline
	public QVal query(
		 QKey qk,
		AlignmentCache which)
	{
		query(qk, which, true);
	}
	
	// Was inline
	public QVal query(
		 QKey qk,
		AlignmentCache which,
		Boolean getLock)
	{
		assert(qk.cacheable());
		AlignmentCache caches[3] = { current_, local_, shared_ };
		for(int i = 0; i < 3; i++) {
			if(caches[i] == null) continue;
			QVal qv = caches[i].query(qk, getLock);
			if(qv != null) {
				if(which != null) which = caches[i];
				return qv;
			}
		}
		return null;
	}
//#endif

	/**
	 * This function is called whenever we start to align a new read or
	 * read substring.  We make key for it and store the key in qk_.
	 * If the sequence is uncacheable, we don't actually add it to the
	 * map but the corresponding reference substrings are still added
	 * to the qlist_.
	 *
	 * Returns:
	 *  -1 if out of memory
	 *  0 if key was found in cache
	 *  1 if key was not found in cache (and there's enough memory to
	 *    add a new key)
	 */
	public int beginAlign(
		 BTDnaString seq,
		 BTString qual,
		QVal qv,              // out: filled in if we find it in the cache
		Boolean getLock = true)
	{
		assert(repOk());
		qk_.init(seq, tmpdnastr_);
		//if(qk_.cacheable() && (qv_ = current_.query(qk_, getLock)) != null) {
		//	// qv_ holds the answer
		//	assert(qv_.valid());
		//	qv = *qv_;
		//	resetRead();
		//	return 1; // found in cache
		//} else
		if(qk_.cacheable()) {
			// Make a QNode for this key and possibly add the QNode to the
			// Red-Black map; but if 'seq' isn't cacheable, just create the
			// QNode (without adding it to the map).
			qv_ = current_.add(qk_, cacheable_, getLock);
		} else {
			qv_ = qvbuf_;
		}
		if(qv_ == null) {
			resetRead();
 			return -1; // Not in memory
		}
		qv_.reset();
		return 0; // Need to search for it
	}
	
	/**
	 * Called when is finished aligning a read (and so is finished
	 * adding associated reference strings).  Returns a copy of the
	 * final QVal object and resets the alignment state of the
	 * current-read cache.
	 *
	 * Also, if the alignment is cacheable, it commits it to the next
	 * cache up in the cache hierarchy.
	 */
	public QVal finishAlign(Boolean getLock = true) {
		if(!qv_.valid()) {
			qv_.init(0, 0, 0);
		}
		// Copy this pointer because we're about to reset the qv_ field
		// to null
		QVal qv = qv_;
		// Commit the contents of the current-read cache to the next
		// cache up in the hierarchy.
		// If qk is cacheable, then it must be in the cache
//#if 0
		if(qk_.cacheable()) {
			AlignmentCache caches[3] = { current_, local_, shared_ };
			AlignmentCache which;
			QVal qv2 = query(qk_, which, true);
			assert(qv2 == qv);
			assert(which == current_);
			for(int i = 1; i < 3; i++) {
				if(caches[i] != null) {
					// Copy this key/value pair to the to the higher
					// level cache and, if its memory is exhausted,
					// clear the cache and try again.
					caches[i].clearCopy(qk_, qv_, current_, getLock);
					break;
				}
			}
		}
//#endif
		// Reset the state in this iface in preparation for the next
		// alignment.
		resetRead();
		assert(repOk());
		return qv;
	}

	/**
	 * A call to this member indicates that the caller has finished
	 * with the last read (if any) and is ready to work on the next.
	 * This gives the cache a chance to reset some of its state if
	 * necessary.
	 */
	public void nextRead() {
		current_.clear();
		resetRead();
		assert(!aligning());
	}
	
	/**
	 * Return true iff we're in the middle of aligning a sequence.
	 */
	public Boolean aligning()  {
		return qv_ != null;
	}
	
	/**
	 * Clears both the local and shared caches.
	 */
	public void clear() {
		if(current_ != null) current_.clear();
		if(local_   != null) local_.clear();
		if(shared_  != null) shared_.clear();
	}
	
	/**
	 * Add an alignment to the running list of alignments being
	 * compiled for the current read in the local cache.
	 */
	 public Boolean addOnTheFly(
		 BTDnaString rfseq, // reference sequence close to read seq
		TIndexOffU topf,            // top in BWT index
		TIndexOffU botf,            // bot in BWT index
		TIndexOffU topb,            // top in BWT' index
		TIndexOffU botb            // bot in BWT' index
		)      // true . lock is not held by caller
	{
		addOnTheFly(rfseq, topf, botf, topb, botb, true);
	}
	
	public Boolean addOnTheFly(
		 BTDnaString rfseq, // reference sequence close to read seq
		TIndexOffU topf,            // top in BWT index
		TIndexOffU botf,            // bot in BWT index
		TIndexOffU topb,            // top in BWT' index
		TIndexOffU botb,            // bot in BWT' index
		Boolean getLock)      // true . lock is not held by caller
	{
		
		assert(aligning());
		assert(repOk());
		BTDnaString tmp;
		QKey sak(rfseq, tmp);
		//assert(sak.cacheable());
		if(current_.addOnTheFly(qv_, sak, topf, botf, topb, botb, getLock)) {
			rangen_++;
			eltsn_ += (botf-topf);
			return true;
		}
		return false;
	}

	/**
	 * Given a QVal, populate the given EList of SATuples with records
	 * describing all of the cached information about the QVal's
	 * reference substrings.
	 */
	//template<int S>
	public void queryQval(
		 QVal qv,
		EList<SATuple, S> satups,
		long nrange,
		long nelt,
		Boolean getLock = true)
	{
		current_.queryQval(qv, satups, nrange, nelt, getLock);
	}

	/**
	 * Return a pointer to the current-read cache object.
	 */
	public AlignmentCache currentCache()  { return current_; }
	
	public long curNumRanges()  { return rangen_; }
	public long curNumElts()    { return eltsn_;  }
	
//#ifndef NDEBUG
	/**
	 * Check that AlignmentCacheIface is internally consistent.
	 */
	public Boolean repOk()  {
		assert(current_ != null);
		//assert_geq(eltsn_, rangen_);
		if(qv_ == null) {
			//assert_eq(0, rangen_);
			//assert_eq(0, eltsn_);
		}
		return true;
	}
//#endif
	
	/**
	 * Return the alignment cache for the current read.
	 */
	 public AlignmentCache current() {
		return current_;
	}

	/**
	 * Reset fields encoding info about the in-process read.
	 */
	protected void resetRead() {
		cacheable_ = false;
		rangen_ = eltsn_ = 0;
		qv_ = null;
	}
}

//#endif /*ALIGNER_CACHE_H_*/