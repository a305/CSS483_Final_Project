package com.uwb.bt2j.indexer;

public class RedBlackNode <K,P>{
	public RedBlackNode parent;
	public RedBlackNode left;
	public RedBlackNode right;
	public boolean red;
	public K key;
	public P payload;
	
	public RedBlackNode grandparent() {
		return parent != null ? parent.parent : null;
	}
	
	public RedBlackNode uncle() {
		if(parent == null) return null; // no parent
		if(parent.parent == null) return null; // parent has no siblings
		return (parent.parent.left == parent) ? parent.parent.right : parent.parent.left;
	}
	
	public boolean isLeftChild() {
		return parent.left == this; 
	}
	
	public boolean isRightChild() {
		return parent.right == this; 
	}
	
	public void replaceChild(RedBlackNode ol, RedBlackNode nw) {
		if(left == ol) {
			left = nw;
		} else {
			right = nw;
		}
	}
	
	public int numChildren() {
		return ((left != null) ? 1 : 0) + ((right != null) ? 1 : 0);
	}
}
