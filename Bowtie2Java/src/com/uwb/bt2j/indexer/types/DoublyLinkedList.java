package com.uwb.bt2j.indexer.types;

public class DoublyLinkedList<T> {
	DoublyLinkedList<T> prev;
	DoublyLinkedList<T> next;
	T payload;
	
	public DoublyLinkedList() {	}
	
	public void toList(EList<T> l) {
		// Add this and all subsequent elements
				DoublyLinkedList<T> cur = this;
				while(cur != null) {
					l.push_back(cur.payload);
					cur = cur.next;
				}
				// Add all previous elements
				cur = prev;
				while(cur != null) {
					l.push_back(cur.payload);
					cur = cur.prev;
				}
	}
}
