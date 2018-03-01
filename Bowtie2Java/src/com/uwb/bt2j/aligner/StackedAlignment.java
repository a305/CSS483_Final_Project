package com.uwb.bt2j.aligner;

public class StackedAlignment {
	protected Boolean inited_;
	protected double trimLS_;
	protected double trimLH_;
	protected double trimRS_;
	protected double trimRH_;
	protected EList<char> stackRef_;
	protected EList<char> stackRel_;
	protected EList<char> stackRead_;
	protected Boolean cigDistMm_;
	protected Boolean cigCalc_;
	protected EList<char> cigOp_;
	protected EList<double> cigRun_;
	protected Boolean mdzCalc_;
	protected EList<char> mdzOp_;
	protected EList<char> mdzChr_;
	protected EList<double> mdzRun_;
	
  public StackedAlignment() {
  stackRef_ = stackRel_ = stackRead_ = cigOp_ = cigRun_ = mdzOp_ = mdzChr_ = mdzRun_ = 7;
    reset();
  }
  
  public void reset() {
  inited_ = false;
		trimLS_ = trimLH_ = trimRS_ = trimRH_ = 0;
		stackRef_.clear();
		stackRel_.clear();
		stackRead_.clear();
		cigDistMm_ = cigCalc_ = false;
		cigOp_.clear();
		cigRun_.clear();
		mdzCalc_ = false;
		mdzOp_.clear();
		mdzChr_.clear();
		mdzRun_.clear();
  }
  
  public Boolean inited() {
	  return inited_;
  }
  
  public void init(BTDnaString s, EList<Edit> ed, double trimLS, double trimLH, double trimRS, double trimRH) {
	  trimLS_ = trimLS;
		trimLH_ = trimLH;
		trimRS_ = trimRS;
		trimRH_ = trimRH;
		stackRef_.clear();
		stackRel_.clear();
		stackRead_.clear();
		double rdoff = trimLS;
		for(double i = 0; i < ed.size(); i++) {
			double pos = ed[i].pos + trimLS;
			while(rdoff < pos) {
				int c = s[rdoff++];
				stackRef_.push_back("ACGTN"[c]);
				stackRel_.push_back('=');
				stackRead_.push_back("ACGTN"[c]);
			}
			if(ed[i].isMismatch()) {
				int c = s[rdoff++];
				stackRef_.push_back(ed[i].chr);
				stackRel_.push_back('X');
				stackRead_.push_back("ACGTN"[c]);
			} else if(ed[i].isRefGap()) {
				int c = s[rdoff++];
				stackRef_.push_back('-');
				stackRel_.push_back('I');
				stackRead_.push_back("ACGTN"[c]);
			} else if(ed[i].isReadGap()) {
				stackRef_.push_back(ed[i].chr);
				stackRel_.push_back('D');
				stackRead_.push_back('-');
			}
		}
		while(rdoff < s.length() - trimRS) {
			int c = s[rdoff++];
			stackRef_.push_back("ACGTN"[c]);
			stackRel_.push_back('=');
			stackRead_.push_back("ACGTN"[c]);
		}
		inited_ = true;
  }
  
  public void leftAlign(Boolean pastMms) {
	  Boolean changed = false;
		double ln = stackRef_.size();
		// Scan left-to-right
		for(double i = 0; i < ln; i++) {
			int rel = stackRel_[i];
			if(rel != '=' && rel != 'X') {
				// Neither a match nor a mismatch - must be a gap
				double glen = 1;
				// Scan further right to measure length of gap
				for(double j = i+1; j < ln; j++) {
					if(rel != (int)stackRel_[j]) break;
					glen++;
				}
				// We've identified a gap of type 'rel' (D = deletion or read
				// gap, I = insertion or ref gap) with length 'glen'.  Now we
				// can try to slide it to the left repeatedly.
				double l = i - 1;
				double r = l + glen;
				EList<char> gp  = ((rel == 'I') ? stackRef_ : stackRead_);
				EList<char> ngp = ((rel == 'I') ? stackRead_ : stackRef_);
				while(l > 0 && ngp[l] == ngp[r]) {
					if(!pastMms && stackRel_[l] == 'X') {
						break;
					}
					swap(gp[l], gp[r]);
					swap(stackRel_[l], stackRel_[r]);
					assert_neq('-', gp[r]);
					assert_eq('-', gp[l]);
					l--; r--;
					changed = true;
				}
				i += (glen-1);
			}
		}
		if(changed) {
			cigCalc_ = mdzCalc_ = false;
		}
  }
  
  public Boolean buildCigar(Boolean xeq) {
		if(cigCalc_) {
			return false; // already done
		}
		cigOp_.clear();
		cigRun_.clear();
		if(trimLS_ > 0) {
			cigOp_.push_back('S');
			cigRun_.push_back(trimLS_);
		}
		double ln = stackRef_.size();
		for(double i = 0; i < ln; i++) {
			char op = stackRel_[i];
			if(!xeq && (op == 'X' || op == '=')) {
				op = 'M';
			}
			double run = 1;
			for(; i + run < ln; run++) {
				char op2 = stackRel_[i + run];
				if(!xeq && (op2 == 'X' || op2 == '=')) {
					op2 = 'M';
				}
				if(op2 != op) {
					break;
				}
			}
			i += (run-1);
			cigOp_.push_back(op);
			cigRun_.push_back(run);
		}
		if(trimRS_ > 0) {
			cigOp_.push_back('S');
			cigRun_.push_back(trimRS_);
		}
		cigCalc_ = true;
		return true;
  }
  
  public Boolean buildMdz() {
		if(mdzCalc_) {
			return false; // already done
		}
		mdzOp_.clear();
		mdzChr_.clear();
		mdzRun_.clear();
		double ln = stackRef_.size();
		for(double i = 0; i < ln; i++) {
			char op = stackRel_[i];
			if(op == '=') {
				double run = 1;
				double ninserts = 0;
				// Skip over matches and insertions (ref gaps)
				for(; i+run < ln; run++) {
					if(stackRel_[i + run] == '=') {
						// do nothing
					} else if(stackRel_[i + run] == 'I') {
						ninserts++;
					} else {
						break;
					}
				}
				i += (run - 1);
				mdzOp_.push_back('='); // = X or G
				mdzChr_.push_back('-');
				mdzRun_.push_back(run - ninserts);
			} else if(op == 'X') {
				mdzOp_.push_back('X'); // = X or G
				mdzChr_.push_back(stackRef_[i]);
				mdzRun_.push_back(1);
			} else if(op == 'D') {
				mdzOp_.push_back('G'); // = X or G
				mdzChr_.push_back(stackRef_[i]);
				mdzRun_.push_back(1);
			}
		}
		mdzCalc_ = true;
		return true;
  }
  
  public void writeCigar(BTString o, char oc) {
		EList<char> op = cigOp_;
		EList<double> run = cigRun_;
		if(o != null || occ != null) {
			char buf[128];
			for(double i = 0; i < op.size(); i++) {
				double r = run[i];
				if(r > 0) {
					itoa10<double>(r, buf);
					if(o != null) {
						o.append(buf);
						o.append(op[i]);
					}
					if(occ != null) {
						COPY_BUF();
						*occ = op[i];
						occ++;
					}
				}
			}
			if(occ != null) {
				*occ = '\0';
			}
		}
  }
  
  public void writeMdz(BTString o, char oc) {
	  char buf[128];
		Boolean mm_last = false;
		Boolean rdgap_last = false;
		Boolean first_print = true;
		EList<char> op = mdzOp_;
		EList<char> ch = mdzChr_;
		EList<double> run = mdzRun_;
		for(double i = 0; i < op.size(); i++) {
			double r = run[i];
			if(r > 0) {
				if(op[i] == '=') {
					// Write run length
					itoa10<double>(r, buf);
					if(o != null)  { o.append(buf); }
					if(occ != null) { COPY_BUF(); }
					first_print = false;
					mm_last = false;
					rdgap_last = false;
				} else if(op[i] == 'X') {
					if(o != null) {
						if(rdgap_last || mm_last || first_print) {
							o.append('0');
						}
						o.append(ch[i]);
					}
					if(occ != null) {
						if(rdgap_last || mm_last || first_print) {
							occ = '0';
							occ++;
						}
						occ = ch[i];
						occ++;
					}
					first_print = false;
					mm_last = true;
					rdgap_last = false;
				} else if(op[i] == 'G') {
					if(o != null) {
						if(mm_last || first_print) {
							o.append('0');
						}
						if(!rdgap_last) {
							o.append('^');
						}
						o.append(ch[i]);
					}
					if(occ != null) {
						if(mm_last || first_print) {
							occ = '0'; occ++;
						}
						if(!rdgap_last) {
							occ = '^'; occ++;
						}
						*occ = ch[i];
						occ++;
					}
					first_print = false;
					mm_last = false;
					rdgap_last = true;
				}
			} // if r > 0
		} // for loop over ops
		if(mm_last || rdgap_last) {
			if(o   != null) { o.append('0'); }
			if(occ != null) { occ = '0'; occ++; }
		}
		if(occ != null) { occ = '\0'; }
  }
}
