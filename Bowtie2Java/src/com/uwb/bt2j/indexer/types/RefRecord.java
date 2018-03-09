package com.uwb.bt2j.indexer.types;

import java.io.File;
import java.io.FileInputStream;
import java.io.OutputStream;

import com.uwb.bt2j.indexer.filebuf.BitpairOutFileBuf;
import com.uwb.bt2j.indexer.filebuf.FileBuf;
import com.uwb.bt2j.indexer.util.WordIO;

import javafx.util.Pair;

public class RefRecord {
	protected long off;
	protected long len;
	protected boolean first;
	
	public enum ReadDir {
		REF_READ_FORWARD, // don't reverse reference sequence
				REF_READ_REVERSE,     // reverse entire reference sequence
				REF_READ_REVERSE_EACH // reverse each unambiguous stretch of reference
	}
	
	public RefRecord() {}
	public RefRecord(long _off, long _len, boolean _first) {
		off = _off;
		len = _len;
		first = _first;
	}
	
	public RefRecord(File in, boolean swap) {
		if(!in.canRead()) {
			System.err.println("Error reading RefRecord offset from FILE");
		}
		if(!in.canRead()) {
			System.err.println("Error reading RefRecord offset from FILE");
		}
		first = in.canRead() ? true : false;
	}
	
	public void write(OutputStream out, boolean be) {
		WordIO.writeU<Long>(out, off, be);
		WordIO.writeU<Long>(out, len, be);
		out.write(first ? 1 : 0);
	}
	
	public RefRecord fastaRefReadSize(FileBuf in, RefReadInParams rparms, boolean first, BitpairOutFileBuf bpout) {
		int c;
		static int lastc = '>'; // last character seen

		// RefRecord params
		TIndexOffU len = 0; // 'len' counts toward total length
		// 'off' counts number of ambiguous characters before first
		// unambiguous character
		int off = 0;

		// Pick off the first carat and any preceding whitespace
		if(first) {
			assert(!in.eof());
			lastc = '>';
			c = in.getPastWhitespace();
			if(in.eof()) {
				// Got eof right away; emit warning
				cerr << "Warning: Empty input file" << endl;
				lastc = -1;
				return RefRecord(0, 0, true);
			}
			assert(c == '>');
		}

		first = true;
		// Skip to the end of the id line; if the next line is either
		// another id line or a comment line, keep skipping
		if(lastc == '>') {
			// Skip to the end of the name line
			do {
				if((c = in.getPastNewline()) == -1) {
					// No more input
					cerr << "Warning: Encountered empty reference sequence" << endl;
					lastc = -1;
					return RefRecord(0, 0, true);
				}
				if(c == '>') {
					cerr << "Warning: Encountered empty reference sequence" << endl;
				}
				// continue until a non-name, non-comment line
			} while (c == '>');
		} else {
			first = false; // not the first in a sequence
			off = 1; // The gap has already been consumed, so count it
			if((c = in.get()) == -1) {
				// Don't emit a warning, since this might legitimately be
				// a gap on the end of the final sequence in the file
				lastc = -1;
				return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
			}
		}

		// Now skip to the first DNA character, counting gap characters
		// as we go
		int lc = -1;
		while(true) {
			int cat = asc2dnacat[c];
			if(rparms.nsToAs && cat >= 2) c = 'A';
			if(cat == 1) {
				break; // to read-in loop
			} else if(cat >= 2) {
				if(lc != -1 && off == 0) off++;
				lc = -1;
				off++; // skip over gap character and increment
			} else if(c == '>') {
				if(off > 0 && lastc == '>') {
					cerr << "Warning: Encountered reference sequence with only gaps" << endl;
				} else if(lastc == '>') {
					cerr << "Warning: Encountered empty reference sequence" << endl;
				}
				lastc = '>';
				//return RefRecord(off, 0, false);
				return RefRecord((TIndexOffU)off, 0, first);
			}
			c = in.get();
			if(c == -1) {
				// End-of-file
				if(off > 0 && lastc == '>') {
					cerr << "Warning: Encountered reference sequence with only gaps" << endl;
				} else if(lastc == '>') {
					cerr << "Warning: Encountered empty reference sequence" << endl;
				}
				lastc = -1;
				//return RefRecord(off, 0, false);
				return RefRecord((TIndexOffU)off, 0, first);
			}
		}
		assert_eq(1, asc2dnacat[c]); // C must be unambiguous base

		// in now points just past the first character of a sequence
		// line, and c holds the first character
		while(c != -1 && c != '>') {
			if(rparms.nsToAs && asc2dnacat[c] >= 2) c = 'A';
			uint8_t cat = asc2dnacat[c];
			int cc = toupper(c);
			if(rparms.bisulfite && cc == 'C') c = cc = 'T';
			if(cat == 1) {
				// It's a DNA character
				assert(cc == 'A' || cc == 'C' || cc == 'G' || cc == 'T');
				// Check for overflow
				if((TIndexOffU)(len + 1) < len) {
					throw RefTooLongException();
				}
				// Consume it
				len++;
				// Output it
				if(bpout != NULL) {
					// output nucleotide
					bpout->write(asc2dna[c]);
				}
				lc = asc2dna[(int)c];
			} else if(cat >= 2) {
				// It's an N or a gap
				lastc = c;
				assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
				return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
			} else {
				// Not DNA and not a gap, ignore it
	#ifndef NDEBUG
				if(!isspace(c)) {
					cerr << "Unexpected character in sequence: ";
					if(isprint(c)) {
						cerr << ((char)c) << endl;
					} else {
						cerr << "(" << c << ")" << endl;
					}
				}
	#endif
			}
			c = in.get();
		}
		lastc = c;
		return RefRecord((TIndexOffU)off, (TIndexOffU)len, first);
	}
	
	public Pair<Integer,Integer> fastaRefReadSizes(EList<FileBuf> in, EList<RefRecord> recs, RefReadInParams rparms, BitpairOutFileBuf bpout, long numSeqs) {
		long unambigTot = 0;
		int bothTot = 0;
		// For each input istream
		for(int i = 0; i < in.size(); i++) {
			bool first = true;
			assert(!in[i]->eof());
			// For each pattern in this istream
			while(!in[i]->eof()) {
				RefRecord rec;
				try {
					rec = fastaRefReadSize(*in[i], rparms, first, bpout);
					if((unambigTot + rec.len) < unambigTot) {
						throw RefTooLongException();
					}
				}
				catch(RefTooLongException& e) {
					cerr << e.what() << endl;
					throw 1;
				}
				// Add the length of this record.
				if(rec.first) numSeqs++;
				unambigTot += rec.len;
				bothTot += rec.len;
				bothTot += rec.off;
				first = false;
				if(rec.len == 0 && rec.off == 0 && !rec.first) continue;
				recs.push_back(rec);
			}
			// Reset the input stream
			in[i].reset();
		}
		return new Pair(
			unambigTot, // total number of unambiguous DNA characters read
			bothTot); // total number of DNA characters read, incl. ambiguous ones
	}
	
	public void reverseRefRecords(EList<RefRecord> src, EList<RefRecord> dst, boolean recursive, boolean verbose) {
		dst.clear();
		{
			EList<RefRecord> cur;
			for(int i = (int)src.size()-1; i >= 0; i--) {
				boolean first = (i == (int)src.size()-1 || src.get(i+1).first);
				// Clause after the || on next line is to deal with empty FASTA
				// records at the end of the 'src' list, which would be wrongly
				// omitted otherwise.
				if(src.get(i).len || (first && src.get(i).off == 0)) {
					cur.push_back(new RefRecord(0, src.get(i).len, first));
					first = false;
				}
				if(src.get(i).off) cur.push_back(new RefRecord(src.get(i).off, 0, first));
			}
			for(int i = 0; i < (int)cur.size(); i++) {
				if(i < (int)cur.size()-1 && cur.get(i).off != 0 && !cur.get(i+1).first) {
					dst.push_back(new RefRecord(cur.get(i).off, cur.get(i+1).len, cur.get(i).first));
					i++;
				} else {
					dst.push_back(cur.get(i));
				}
			}
		}
	}
	
	public static RefRecord fastaRefReadAppend(FileBuf in, boolean first, TStr dst, long dstoff, RefReadInParams rparms, String name) {
		int c;
		char lastc = '>';
		if(first) {
			c = in.getPastWhitespace();
			if(c != '>') {
				cerr << "Reference file does not seem to be a FASTA file" << endl;
				throw 1;
			}
			lastc = c;
		}

		// RefRecord params
		int len = 0;
		int off = 0;
		first = true;

		int ilen = dstoff;

		// Chew up the id line; if the next line is either
		// another id line or a comment line, keep chewing
		int lc = -1; // last-DNA char variable for color conversion
		c = lastc;
		if(c == '>' || c == '#') {
			do {
				while (c == '#') {
					if((c = in.getPastNewline()) == -1) {
						lastc = -1;
						goto bail;
					}
				}
				assert_eq('>', c);
				while(true) {
					c = in.get();
					if(c == -1) {
						lastc = -1;
						goto bail;
					}
					if(c == '\n' || c == '\r') {
						while(c == '\r' || c == '\n') c = in.get();
						if(c == -1) {
							lastc = -1;
							goto bail;
						}
						break;
					}
					if (name) name->push_back(c);
				}
				// c holds the first character on the line after the name
				// line
				if(c == '>') {
					// If there's another name line immediately after this one,
					// discard the previous name and start fresh with the new one
					if (name) name->clear();
				}
			} while (c == '>' || c == '#');
		} else {
			ASSERT_ONLY(int cc = toupper(c));
			assert(cc != 'A' && cc != 'C' && cc != 'G' && cc != 'T');
			first = false;
		}

		// Skip over an initial stretch of gaps or ambiguous characters.
		// For colorspace we skip until we see two consecutive unambiguous
		// characters (i.e. the first unambiguous color).
		while(true) {
			int cat = asc2dnacat[c];
			if(rparms.nsToAs && cat >= 2) {
				c = 'A';
			}
			int cc = toupper(c);
			if(rparms.bisulfite && cc == 'C') c = cc = 'T';
			if(cat == 1) {
				// This is a DNA character
				if(rparms.color) {
					if(lc != -1) {
						// Got two consecutive unambiguous DNAs
						break; // to read-in loop
					}
					// Keep going; we need two consecutive unambiguous DNAs
					lc = asc2dna[(int)c];
					// The 'if(off > 0)' takes care of the case where
					// the reference is entirely unambiguous and we don't
					// want to incorrectly increment off.
					if(off > 0) off++;
				} else {
					break; // to read-in loop
				}
			} else if(cat >= 2) {
				if(lc != -1 && off == 0) {
					off++;
				}
				lc = -1;
				off++; // skip it
			} else if(c == '>') {
				lastc = '>';
				goto bail;
			}
			c = in.get();
			if(c == -1) {
				lastc = -1;
				goto bail;
			}
		}
		if(first && rparms.color && off > 0) {
			// Handle the case where the first record has ambiguous
			// characters but we're in color space; one of those counts is
			// spurious
			off--;
		}
		// in now points just past the first character of a sequence
		// line, and c holds the first character
		while(true) {
			// Note: can't have a comment in the middle of a sequence,
			// though a comment can end a sequence
			int cat = asc2dnacat[c];
			if(cat == 1) {
				// Consume it
				if(!rparms.color || lc != -1) len++;
				// Add it to referenece buffer
				if(rparms.color) {
					dst.set((char)dinuc2color[asc2dna[(int)c]][lc], dstoff++);
				} else if(!rparms.color) {
					dst.set(asc2dna[c], dstoff++);
				}
				lc = asc2dna[(int)c];
			}
			c = in.get();
			if(rparms.nsToAs && asc2dnacat[c] >= 2) c = 'A';
			if (c == -1 || c == '>' || c == '#' || asc2dnacat[c] >= 2) {
				lastc = c;
				break;
			}
			if(rparms.bisulfite && toupper(c) == 'C') c = 'T';
		}

	  bail:
		// Optionally reverse the portion that we just appended.
		// ilen = length of buffer before this last sequence was appended.
		if(rparms.reverse == REF_READ_REVERSE_EACH) {
			// Find limits of the portion we just appended
			int nlen = dstoff;
			dst.reverseWindow(ilen, nlen);
		}
		return new RefRecord((long)off, (long)len, first);
	}
}
