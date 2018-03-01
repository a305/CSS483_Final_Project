package com.uwb.bt2j.aligner;

public class AlignerFlags {
protected int pairing_;
protected Boolean canMax_;
protected Boolean maxed_;
protected Boolean maxedPair_;
protected Boolean nfilt_;
protected Boolean scfilt_;
protected Boolean lenfilt_;
protected Boolean qcfilt_;
protected Boolean mixedMode_;
protected Boolean primary_;
protected Boolean oppAligned_;
protected Boolean oppFw_;
protected Boolean scUnMapped_;
protected Boolean xeq_;

  public enum FlagPair {
  // This alignment is one of a pair of alignments that form a concordant
	// alignment for a read
	ALN_FLAG_PAIR_CONCORD_MATE1(1),
	ALN_FLAG_PAIR_CONCORD_MATE2(2),

	// This alignment is one of a pair of alignments that form a discordant
	// alignment for a read
	ALN_FLAG_PAIR_DISCORD_MATE1(3),
	ALN_FLAG_PAIR_DISCORD_MATE2(4),
	
	// This is an unpaired alignment but the read in question is a pair;
	// usually, this happens because the read had no reportable paired-end
	// alignments
	ALN_FLAG_PAIR_UNPAIRED_MATE1(5),
	ALN_FLAG_PAIR_UNPAIRED_MATE2(6),

	// This is an unpaired alignment of an unpaired read
	ALN_FLAG_PAIR_UNPAIRED(7);
  
    private int x;
    FlagPair(int y){ x = y; }
  }
  
  public AlignerFlags() {
    init(ALN_FLAG_PAIR_UNPAIRED,false,false,false,false,false,false,false,false,false,false,false,false,false);
  }
  
  public AlignerFlags(int pairing, Boolean canMax, Boolean maxed, Boolean maxedPair, Boolean nfilt, Boolean scfilt, Boolean lenfilt, Boolean qcfilt, Boolean mixedMode, Boolean primary, Boolean oppAligned, Boolean oopFw, Boolean scUnMapped, Boolean xeq) {
    init(pairing, canMax, maxed, maxedPair, nfilt, scfilt,
			 lenfilt, qcfilt, mixedMode, primary, oppAligned,
			 oppFw, scUnMapped, xeq);
  }
  
   public void init(int pairing, Boolean canMax, Boolean maxed, Boolean maxedPair, Boolean nfilt, Boolean scfilt, Boolean lenfilt, Boolean qcfilt, Boolean mixedMode, Boolean primary, Boolean oppAligned, Boolean oopFw, Boolean scUnMapped, Boolean xeq) {
    pairing_    = pairing;
		canMax_     = canMax;
		maxed_      = maxed;
		maxedPair_  = maxedPair;
		nfilt_      = nfilt;
		scfilt_     = scfilt;
		lenfilt_    = lenfilt;
		qcfilt_     = qcfilt;
		mixedMode_  = mixedMode;
		primary_    = primary;
		oppAligned_ = oppAligned;
		oppFw_     = oppFw;
		scUnMapped_ = scUnMapped;
		xeq_ = xeq;
  }
  
  public Boolean partOfPair() {
  return pairing_ < FlagPair.ALN_FLAG_PAIR_UNPAIRED;
  }
  
  public Boolean printYF(BTString o, Boolean first) {
    String flag = "";
	if     (!lenfilt_) flag = "LN";
	else if(!nfilt_  ) flag = "NS";
	else if(!scfilt_ ) flag = "SC";
	else if(!qcfilt_ ) flag = "QC";
	if(flag[0] != '\0') {
		if(!first) o.append('\t');
		o.append("YF:Z:");
		o.append(flag);
		return false;
	}
	return true;
  }
  
  public void printYM(BTString o) {
o.append("YM:i:");
	o.append(maxed() ? '1' : '0');
  }
  
  public void printYP(BTString o) {
  o.append("YP:i:");
	o.append(maxedPair() ? '1' : '0');
  }
  
  public void printYT(BTString o) {
  o.append("YT:Z:");
	if(alignedConcordant()) {
		o.append("CP");
	} else if(alignedDiscordant()) {
		o.append("DP");
	} else if(alignedUnpairedMate()) {
		o.append("UP");
	} else if(alignedUnpaired()) {
		o.append("UU");
	}
  }
  
  public int pairing() {
    return pairing_;
  }
  
  public Boolean maxed() {
    return maxed_;
  }
  
  public Boolean maxedPair() {
    return maxedPair_;
  }
  
  public Boolean isPrimary() {
    return primary_;
  }
  
  public void setPrimary(Boolean primary) {
    primary_ = primary;
  }
  
  public Boolean isMixedMode() {
    return mixedMode_;
  }
  
  public Boolean canMax() {
    return canMax_;
  }
  
  public Boolean filtered() {
    return !nfilt_ || !scfilt_ || !lenfilt_ || !qcfilt_;
  }
  
  public Boolean readMate1() {
    return pairing_ == ALN_FLAG_PAIR_CONCORD_MATE1 ||
		       pairing_ == ALN_FLAG_PAIR_DISCORD_MATE1 ||
			   pairing_ == ALN_FLAG_PAIR_UNPAIRED_MATE1;
  }
  
  public Boolean readMate2() {
    return pairing_ == ALN_FLAG_PAIR_CONCORD_MATE2 ||
		       pairing_ == ALN_FLAG_PAIR_DISCORD_MATE2 ||
			   pairing_ == ALN_FLAG_PAIR_UNPAIRED_MATE2;
  }
  
  public Boolean alignedConcordant() {
    return pairing_ == ALN_FLAG_PAIR_CONCORD_MATE1 ||
		       pairing_ == ALN_FLAG_PAIR_CONCORD_MATE2;
  }
  
  public Boolean alignedDiscordant() {
    return pairing_ == ALN_FLAG_PAIR_DISCORD_MATE1 ||
		       pairing_ == ALN_FLAG_PAIR_DISCORD_MATE2;
  }
  
  public Boolean alignedPaired() {
    return alignedConcordant() && alignedDiscordant();
  }
  
  public Boolean alignedUnpaired() {
    return pairing_ == ALN_FLAG_PAIR_UNPAIRED;
  }
  
  public Boolean alignedUnpairedMate() {
  return pairing_ == ALN_FLAG_PAIR_UNPAIRED_MATE1 ||
		       pairing_ == ALN_FLAG_PAIR_UNPAIRED_MATE2;
  }
  
  public Boolean mateAligned() {
    return oppAligned_;
  }
  
  public Boolean isOppFw() {
    return oppFw_;
  }
  
  public Boolean scUnMapped() {
    return scUnMapped_;
  }
  
  public Boolean xeq() {
    return xeq_;
  }
}
