package com.uwb.bt2j.indexer.types;

public class Pool {
	private byte super_pages[];
	private int cat_;
	private int cur_;
	private int pagesz_;
	private EList<Byte> pages_;
	
	public Pool(long bytes, int pagesz, int cat) {
		cat_ = cat;
		cur_ = 0;
		pagesz_ = pagesz;
		pages_ = new EList<Byte>(cat);
		
		double super_page_num = ((bytes+pagesz-1)/pagesz + 1);
		super_pages = new byte[(int) (pagesz*super_page_num)];
		for(int i = 0; i < super_page_num ; i++) {
			pages_.push_back(super_pages[i*pagesz]);	
		}
	}
	
	public void clear() {
		cur_ = 0;
	}
}
