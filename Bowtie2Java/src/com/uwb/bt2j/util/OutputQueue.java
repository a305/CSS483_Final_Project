package com.uwb.bt2j.util;

public class OutputQueue {
	protected OutFileBuf obuf_;
	protected long cur_;
	protected long nstarted_;
	protected long nfinished_;
	protected long nflushed_;
	protected EList<BTString> lines_;
	protected EList<Boolean> started_;
	protected Boolean reorder_;
	protected Boolean threadSafe_;
	protected double nthreads_;
	protected BTString perThreadBuf;
	protected int perThreadCounter;
	protected int perThreadBufSize_;
	
	public OutputQueue(OutFileBuf& obuf,
			Boolean reorder,
			double nthreads,
			Boolean threadSafe,
			int perThreadBufSize,
			long rdid) {
		nstarted_=0;
		if(!reorder)
		{
			perThreadBuf = new BTString[nthreads_];
			perThreadCounter = new int[nthreads_];
			double i = 0;
			for(i=0;i<nthreads_;i++)
			{
				perThreadBuf[i] = new BTString[perThreadBufSize_];
				perThreadCounter[i] = 0;
			}
		}
	}
	
}
	private void beginReadImpl(long rdid, double threadId) {
		if(reorder_) {
			if(rdid - cur_ >= lines_.size()) {
				// Make sure there's enough room in lines_, started_ and finished_
				size_t oldsz = lines_.size();
				lines_.resize(rdid - cur_ + 1);
				started_.resize(rdid - cur_ + 1);
				finished_.resize(rdid - cur_ + 1);
				for(size_t i = oldsz; i < lines_.size(); i++) {
					started_[i] = finished_[i] = false;
				}
			}
			started_[rdid - cur_] = true;
			finished_[rdid - cur_] = false;
	}
		
	public void beginRead(long rdid, double threadId) {
		if(reorder_ && threadSafe_) {
			ThreadSafe ts(mutex_m);
			beginReadImpl(rdid, threadId);
		} else {
			beginReadImpl(rdid, threadId);
		}
	}
	
	private void finishReadImpl(BTString rec, long rdid, double threadId) {
		if(reorder_) {
			lines_[rdid - cur_] = rec;
			nfinished_++;
			finished_[rdid - cur_] = true;
			flush(false, false); // don't force; already have lock
		} else {
			// obuf_ is the OutFileBuf for the output file
			int i = 0;
			for(i=0; i < perThreadBufSize_; i++)
			{
				obuf_.writeString(perThreadBuf[threadId][i]);
				//TODO: turn these into atomics
				nfinished_++;
				nflushed_++;
			}
			perThreadCounter[threadId] = 0;
		}
	}
	
	public void finishRead(BTString rec, long rdid, double threadId) {
		if(reorder_ || perThreadCounter[threadId] >= perThreadBufSize_) {
			if(threadSafe_) {
				ThreadSafe ts(mutex_m);
				finishReadImpl(rec, rdid, threadId);
			} else {
				finishReadImpl(rec, rdid, threadId);
			}
		}
		if(!reorder_) {
			perThreadBuf[threadId][perThreadCounter[threadId]++] = rec;
		}
	}
	
	private void flushImpl(Boolean force) {
		if(!reorder_) {
			double i = 0;
			int j = 0;
			for(i=0;i<nthreads_;i++)
			{
				for(j=0;j<perThreadCounter[i];j++)
				{
					obuf_.writeString(perThreadBuf[i][j]);
					nfinished_++;
					nflushed_++;
				}
				perThreadCounter[i]=0;
			}
			return;
		}
		double nflush = 0;
		while(nflush < finished_.size() && finished_[nflush]) {
			nflush++;
		}
		// Waiting until we have several in a row to flush cuts down on copies
		// (but requires more buffering)
		if(force || nflush >= NFLUSH_THRESH) {
			for(size_t i = 0; i < nflush; i++) {
				obuf_.writeString(lines_[i]);
			}
			lines_.erase(0, nflush);
			started_.erase(0, nflush);
			finished_.erase(0, nflush);
			cur_ += nflush;
			nflushed_ += nflush;
		}
	}
	
	public void flush(Boolean force, Boolean getLock) {
		if(getLock && threadSafe_) {
			ThreadSafe ts(mutex_m);
			flushImpl(force);
		} else {
			flushImpl(force);
		}
	}
	
	public double size() {
		return lines_.size();
	}
	
	public long numFlushed() {
		return nflushed_;
	}
	
	public long numStarted() {
		return nstarted_;
	}
	
	public long numFinished() {
		return nfinished_;
	}
	
	class OutputQueueMark {
		protected OutputQueue q_;
		protected final BTString rec_;
		protected long rdid_;
		protected double threadId_;
		
		public OutputQueueMark(OutputQueue q, BTString rec, long rdid, double threadId) {
			q_ = q;
			rec_ = rec;
			rdid_ = rdid;
			threadId_ = threadId;
			q_.beginRead(rec_,  rdid_, threadId_);
		}
	}
}