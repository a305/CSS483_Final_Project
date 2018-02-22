package com.uwb.bt2j.util;
class Formats {
  enum file_format {
    FASTA = 1,
    FASTA_CONT,
	  FASTQ,
	  INTERLEAVED,
	  TAB_MATE5,
	  TAB_MATE6,
	  RAW,
	  CMDLINE,
	  QSEQ
  }
  
  public static final String[] file_format_names = {
    "Invalid!",
	  "FASTA",
	  "FASTA sampling",
	  "FASTQ",
	  "Tabbed mated",
	  "Raw",
	  "Command line",
  	"Chain file",
	  "Random",
	  "Qseq"
  };
}
