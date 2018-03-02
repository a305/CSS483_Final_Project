package com.uwb.bt2j.aligner.dp;

import com.uwb.bt2j.aligner.AlignmentScore;

public class DynamicProgrammingNucFrame {
	double   nedsz;    // size of the nucleotide edit list at branch (before
    // adding the branch edit)
	double   aedsz;    // size of ambiguous nucleotide edit list at branch
	double   celsz;    // size of cell-traversed list at branch
	double   row;      // row of cell where branch occurred
	double   col;      // column of cell where branch occurred
	double   gaps;     // number of gaps before branch occurred
	double   readGaps; // number of read gaps before branch occurred
	double   refGaps;  // number of ref gaps before branch occurred
	AlignmentScore score;    // score where branch occurred
	int      ct;       // table type (oall, rdgap or rfgap)
	
	public void init(
			double   nedsz_,
			double   aedsz_,
			double   celsz_,
			double   row_,
			double   col_,
			double   gaps_,
			double   readGaps_,
			double   refGaps_,
			AlignmentScore score_,
			int      ct_) {
		nedsz    = nedsz_;
		aedsz    = aedsz_;
		celsz    = celsz_;
		row      = row_;
		col      = col_;
		gaps     = gaps_;
		readGaps = readGaps_;
		refGaps  = refGaps_;
		score    = score_;
		ct       = ct_;
	}
}
