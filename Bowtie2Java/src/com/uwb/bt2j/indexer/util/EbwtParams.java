package com.uwb.bt2j.indexer.util;

import java.io.IOException;
import java.io.OutputStream;

public class EbwtParams {
	public long _len;
	public long _bwtLen;
	public long _sz;
	public long _bwtSz;
	public int  _lineRate;
	public int  _origOffRate;
	public int  _offRate;
	public long _offMask;
	public int  _ftabChars;
	public int _eftabLen;
	public int _eftabSz;
	public long _ftabLen; 
	public long _ftabSz; 
	public long _offsLen;
	public long _offsSz;
	public int _lineSz; 
	public int _sideSz; 
	public int _sideBwtSz; 
	public int _sideBwtLen; 
	public long _numSides;
	public long _numLines;
	public long _ebwtTotLen;
	public long _ebwtTotSz;
	public boolean     _color;
	public boolean     _entireReverse;
	
	public EbwtParams(
		long len,
		int lineRate,
		int offRate,
		int ftabChars,
		boolean color,
		boolean entireReverse){
		init(len, lineRate, offRate, ftabChars, color, entireReverse);
	}
	
	public EbwtParams(EbwtParams eh) {
		init(eh._len, eh._lineRate, eh._offRate,
			     eh._ftabChars, eh._color, eh._entireReverse);
	}
	
	public void init(
			long len,
			int lineRate,
			int offRate,
			int ftabChars,
			boolean color,
			boolean entireReverse) {
		_color = color;
		_entireReverse = entireReverse;
		_len = len;
		_bwtLen = _len + 1;
		_sz = (len+3)/4;
		_bwtSz = (len/4 + 1);
		_lineRate = lineRate;
		_origOffRate = offRate;
		_offRate = offRate;
		_offMask = IndexTypes.OFF_MASK << _offRate;
		_ftabChars = ftabChars;
		_eftabLen = _ftabChars*2;
		_eftabSz = _eftabLen*IndexTypes.OFF_SIZE;
		_ftabLen = (1 << (_ftabChars*2))+1;
		_ftabSz = _ftabLen*IndexTypes.OFF_SIZE;
		_offsLen = (_bwtLen + (1 << _offRate) - 1) >> _offRate;
		_offsSz = (long)_offsLen*IndexTypes.OFF_SIZE;
		_lineSz = 1 << _lineRate;
		_sideSz = _lineSz * 1 /* lines per side */;
		_sideBwtSz = _sideSz - IndexTypes.OFF_SIZE*4;
		_sideBwtLen = _sideBwtSz*4;
		_numSides = (_bwtSz+(_sideBwtSz)-1)/(_sideBwtSz);
		_numLines = _numSides * 1 /* lines per side */;
		_ebwtTotLen = _numSides * _sideSz;
		_ebwtTotSz = _ebwtTotLen;
	}
	
	public long len()           	{ return _len; }
	public long lenNucs()   	    { return _len + (_color ? 1 : 0); }
	public long bwtLen()        	{ return _bwtLen; }
	public long sz()           		{ return _sz; }
	public long bwtSz()         	{ return _bwtSz; }
	public int  lineRate()      	{ return _lineRate; }
	public int  origOffRate()   	{ return _origOffRate; }
	public int  offRate()       	{ return _offRate; }
	public long offMask()       	{ return _offMask; }
	public int  ftabChars()     	{ return _ftabChars; }
	public int  eftabLen()      	{ return _eftabLen; } 
	public int  eftabSz()       	{ return _eftabSz; } 
	public long ftabLen()       	{ return _ftabLen; }
	public long ftabSz()        	{ return _ftabSz; }
	public long offsLen()       	{ return _offsLen; }
	public long offsSz()        	{ return _offsSz; }
	public int  lineSz()        	{ return _lineSz; } 
	public int  sideSz()        	{ return _sideSz; } 
	public int  sideBwtSz()     	{ return _sideBwtSz; } 
	public int  sideBwtLen()    	{ return _sideBwtLen; } 
	public long numSides()     		{ return _numSides; }
	public long numLines()      	{ return _numLines; }
	public long ebwtTotLen()    	{ return _ebwtTotLen; }
	public long ebwtTotSz()     	{ return _ebwtTotSz; }
	public boolean color()          { return _color; }
	public boolean entireReverse()  { return _entireReverse; }
	
	public void setOffRate(int __offRate) {
		_offRate = __offRate;
		_offMask = IndexTypes.OFF_MASK << _offRate;
		_offsLen = (_bwtLen + (1 << _offRate) - 1) >> _offRate;
		_offsSz = (long)_offsLen * IndexTypes.OFF_SIZE;
	}
	
	public void print(OutputStream out) throws IOException {
		out.write(("Headers:" + "\n"
	    + "    len: "          + _len + "\n"
	    + "    bwtLen: "       + _bwtLen + "\n"
	    + "    sz: "           + _sz + "\n"
	    + "    bwtSz: "        + _bwtSz + "\n"
	    + "    lineRate: "     + _lineRate + "\n"
	    + "    offRate: "      + _offRate + "\n"
	    + "    offMask: 0x"    + Long.toHexString(_offMask) + "\n"
	    + "    ftabChars: "    + _ftabChars + "\n"
	    + "    eftabLen: "     + _eftabLen + "\n"
	    + "    eftabSz: "      + _eftabSz + "\n"
	    + "    ftabLen: "      + _ftabLen + "\n"
	    + "    ftabSz: "       + _ftabSz + "\n"
	    + "    offsLen: "      + _offsLen + "\n"
	    + "    offsSz: "       + _offsSz + "\n"
	    + "    lineSz: "       + _lineSz + "\n"
	    + "    sideSz: "       + _sideSz + "\n"
	    + "    sideBwtSz: "    + _sideBwtSz + "\n"
	    + "    sideBwtLen: "   + _sideBwtLen + "\n"
	    + "    numSides: "     + _numSides + "\n"
	    + "    numLines: "     + _numLines + "\n"
	    + "    ebwtTotLen: "   + _ebwtTotLen + "\n"
	    + "    ebwtTotSz: "    + _ebwtTotSz + "\n"
	    + "    color: "        + _color + "\n"
	    + "    reverse: "      + _entireReverse + "\n").getBytes());
	}
}
