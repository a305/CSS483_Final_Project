package com.uwb.bt2j.indexer;

public class LarssonSadakane{
	public int I, V[], r, h;
	
	public static int[] suffixArray_LS(int[] S) {
		int n = S.length;
		int[] SA = new int[n];
		int[] ISA = new int[n];
		radixSort_LS(S,SA,ISA);
		boolean iteration;
		int h = 1;
		do {
			iteration = false;
			for(int i = 0; i < n; i = ISA[SA[i]] + 1) {
				if(i < ISA[SA[i]]) {
					ternary_split_quicksort_LS(SA,i,ISA[SA[i]], ISA, h);
					iteration = true;
				}
			}
			h *= 2;
		} while(iteration);
		return SA;
	}
	
	private static void radixSort_LS(int[] S, int[] SA, int[] ISA) {
		int n = S.length;
		int m = 0;
		for(int i = 0; i < n; i++)
			if(m < S[i])
				m = S[i];
		m++;
		int[] count = new int[m];
		for(int i = 0; i < m; i++)
			count[i] = 0;
		for(int i = 0; i < n; i++)
			count[S[i]]++;
		for(int i = 1; i < m; i++)
			count[i] += count[i - 1];
		for(int i = n - 1; 0 <= i; i--)
			ISA[i] = count[S[i]] - 1;
		for(int i = n - 1; 0 <= i; i--){
			SA[count[S[i]] - 1] = i;
			count[S[i]]--;
		}
	}
	
	private static void ternary_split_quicksort_LS(int[] SA, int l, int r, int[] ISA, int h) {
		if(l > r)
			return;
		if(l == r){
			ISA[SA[l]] = l;
			return;
		}
		
		int v = ISA[SA[l + (int) Math.floor(Math.random() * (r - l + 1))] + h];
		int i = l;
		int mi = l;
		int j = r;
		int mj = r;
		int tmp;
		for(;;) {
			for(; i <= j && ISA[SA[i] + h] <= v; i++){
				if(ISA[SA[i] + h] == v) {
					tmp = SA[i];
					SA[i] = SA[mi];
					SA[mi] = tmp;
					mi++;
				}
			}
			for(; i <= j && v <= ISA[SA[j] + h]; j--){
				if(ISA[SA[j] + h] == v) {
					tmp = SA[j];
					SA[j] = SA[mj];
					SA[mj] = tmp;
					mj--;
				}
			}
			if(i > j)
				break;
			tmp = SA[i];
			SA[i] = SA[j];
			SA[j] = tmp;
			i++;
			j--;
		}
		for(mi--,i--; l <= mi; mi--, i--) {
			tmp = SA[i];
			SA[i] = SA[mi];
			SA[mi] = tmp;
		}
		for(mj++,j++; mj <= r; mj++, j++) {
			tmp = SA[j];
			SA[j] = SA[mj];
			SA[mj] = tmp;
		}
		ternary_split_quicksort_LS(SA, l, i, ISA, h);
		for(int k = i + 1; k < j; k++)
			ISA[SA[k]] = j - 1;
		ternary_split_quicksort_LS(SA, j, r, ISA, h);
	}
}
