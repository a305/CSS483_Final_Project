package com.uwb.bt2j.aligner;

public class RunningStat {
	private int m_n;
	private double m_tot;
	private double m_oldM, m_newM, m_oldS, m_newS;
	
	public RunningStat() {
		m_n = 0;
		m_tot = 0.0;
	}
	
	public void clear() {
		m_n = 0;
		m_tot = 0.0;
	}
	
	public void push(float x) {
		m_n++;
		m_tot += x;
		// See Knuth TAOCP vol 2, 3rd edition, page 232
		if (m_n == 1) {
			m_oldM = m_newM = x;
			m_oldS = 0.0;
		} else {
			m_newM = m_oldM + (x - m_oldM)/m_n;
			m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
			// set up for next iteration
			m_oldM = m_newM;
			m_oldS = m_newS;
		}
	}
	
	public int num() {
		return m_n;
	}

	public double tot() {
		return m_tot;
	}

	public double mean() {
		return (m_n > 0) ? m_newM : 0.0;
	}

	public double variance() {
		return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
	}

	public double stddev() {
		return Math.sqrt(variance());
	}
}
