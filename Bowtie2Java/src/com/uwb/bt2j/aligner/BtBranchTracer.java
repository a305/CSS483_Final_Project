package com.uwb.bt2j.aligner;
import com.uwb.bt2j.util.EList;
import com.uwb.bt2j.util.types.ELList;

class BtBranchTracer {
  protected BtBranchProblem prob_;
  protected EFactory<BtBranch> bs_;
  protected ELList<pair<double,double>> seenPaths_;
  protected ELSet<double> sawcell_;
  protected EList<pair<double,double>> unsorted_;
  protected EList<double> sorted1_;
  protected EList<double> sorted2_;
  protected EList<double> solutions_;
  protected Boolean sortedSel_;
  protected double cur_;
  protected double nmm_;
  protected double nnmm_;
  protected double nrdop_;
  protected double nrfop_;
  protected double nrdex_;
  protected double nrfex_;
  protected double nmmPrune_;
  protected double nnmmPrune_;
  protected double nrdopPrune_;
  protected double nrfopPrune_;
  protected double nrdexPrune_;
  protected double nrfexPrune_;
  protected double row_;
  protected double col_;
  protected Boolean doTri_;
  protected EList<CpQuad> sq_;
  protected EList<double> ndep_;
  protected ELList<CpQuad> tri_;

  public BtBranchTracer() {
    
  }
  
  public void add() {
    
  }
  
  public void addSolution() {
  
  }
  
  public void examineBranch() {
  
  }
  
  public void addOffshoots() {
    
  }
  
  public double best() {
  
  }
  
  public Boolean empty() {
    
  }
  
  public double size() {
  
  }
  
  public Boolean emptySolution() {
  
  }
  
  public double sizeSolution() {
  
  }
  
  public void flushUnsorted() {
    
  }
  
  public void initRef() {
    
  }
  
  public void initBt() {
    
  }
  
  public Boolean nextAlignment() {
    
  }
  
  public Boolean inited() {
  
  }
  
  public Boolean doTri() {
  
  }
  
  public void triangleFill() {
  
  }
  
  public void squareFill() {
  
  }
  
  protected Boolean nextAlignmentBacktrace() {
    
  }
  
  protected Boolean nextAlignmentFill() {
    
  }
  
  protected Boolean trySolutions() {
  
  }
  
  protected int trySolution() {
  
  }
}