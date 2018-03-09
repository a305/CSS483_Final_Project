package com.uwb.bt2j.indexer.util;
public class Formats {
  public enum FileFormat {
    FASTA(1),
    FASTA_CONT(2),
	  FASTQ(3),
	  INTERLEAVED(4),
	  TAB_MATE5(5),
	  TAB_MATE6(6),
	  RAW(7),
	  CMDLINE(8),
	  QSEQ(9);
	  
	  private int y;
	  FileFormat(int x){y = x;}
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
