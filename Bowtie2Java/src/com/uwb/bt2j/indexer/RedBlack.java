package com.uwb.bt2j.indexer;

public class RedBlack <K, P>{
	public int keys_;
	public int cur_;
	public int curPage_;
	public int perPage_;
	public RedBlackNode root_;
	public EList<RedBlackNode> pages_;
	public int intenseRepOkCnt_;
	
	public RedBlack(int pageSz, int cat) {
		//perPage_ = pageSz/
		pages_ = new EList(cat);
		clear();
	}
	
	public RedBlackNode lookup(K key) {
		RedBlackNode cur = root_;
		while(cur != null) {
			if(cur == key) return cur;
			if(cur < key) {
				cur = cur.right;
			} else {
				cur = cur.left;
			}
		}
		return null;
	}
	
	public RedBlackNode add(Pool p, K key, Boolean added) {
		// Look for key; if it's not there, get its parent
				RedBlackNode cur = root_;
				RedBlackNode parent = null;
				boolean leftChild = true;
				while(cur != null) {
					if(cur == key) {
						// Found it; break out of loop with cur != null
						break;
					}
					parent = cur;
					if(cur < key) {
						if((cur = cur.right) == null) {
							// Fell off the bottom of the tree as the right
							// child of parent 'lastCur'
							leftChild = false;
						}
					} else {
						if((cur = cur.left) == null) {
							// Fell off the bottom of the tree as the left
							// child of parent 'lastCur'
							leftChild = true;
						}
					}
				}
				if(cur != null) {
					// Found an entry; assert if we weren't supposed to
					if(added != null) added = false;
				} else {
					if(!addNode(p, cur)) {
						// Exhausted memory
						return null;
					}
					// Initialize new node
					cur.key = key;
					cur.left = cur.right = null;
					cur.red = true; // red until proven black
					keys_++;
					if(added != null) added = true;
					// Put it where we know it should go
					addNode(cur, parent, leftChild);
				}
				return cur; // return the added or found node
	}
	
	public void clear() {
		cur_ = curPage_ = 0;
		root_ = null;
		keys_ = 0;
		intenseRepOkCnt_ = 0;
		pages_.clear();
	}
	
	public int size() {
		return keys_;
	}
	
	public boolean empty() {
		return keys_ == 0;
	}
	
	public boolean addNode(Pool p, RedBlackNode node) {
		// Allocation of the first page
				if(pages_.size() == 0) {
					if(addPage(p) == null) {
						node = null;
						return false;
					}
				}
				if(cur_ == perPage_) {
					if(curPage_ == pages_.size()-1 && addPage(p) == null) {
						return false;
					}
					cur_ = 0;
					curPage_++;
				}
				node = pages_[curPage_][cur_];
				cur_++;
				return true;
	}
	
	protected void leftRotate(RedBlackNode n) {
		RedBlackNode r = n.right;
		n.right = r.left;
		if(n.right != null) {
			n.right.parent = n;
		}
		r.parent = n.parent;
		n.parent = r;
		r.left = n;
		if(r.parent != null) {
			r.parent.replaceChild(n, r);
		}
		if(root_ == n) root_ = r;
	}
	
	protected void rightRotate(RedBlackNode n) {
		RedBlackNode r = n.left;
		n.left = r.right;
		if(n.left != null) {
			n.left.parent = n;
		}
		r.parent = n.parent;
		n.parent = r;
		r.right = n;
		if(r.parent != null) {
			r.parent.replaceChild(n, r);
		}
		if(root_ == n) root_ = r;
	}
	
	public void addNode(RedBlackNode n, RedBlackNode parent, boolean leftChild) {
		if(parent == null) {
			// Case 1: inserted at root
			root_ = n;
			root_.red = false; // root must be black
			n.parent = null;
		} else {
			// Add new node to tree
			if(leftChild) {
				parent.left = n;
			} else {
				parent.right = n;
			}
			n.parent = parent;
			int thru = 0;
			while(true) {
				thru++;
				parent = n.parent;
				if(parent == null && n.red) {
					n.red = false;
				}
				if(parent == null || !parent.red) {
					break;
				}
				RedBlackNode uncle = n.uncle();
				RedBlackNode gparent = n.grandparent();
				boolean uncleRed = (uncle != null ? uncle.red : false);
				if(uncleRed) {
					// Parent is red, uncle is red; recursive case
					parent.red = uncle.red = false;
					gparent.red = true;
					n = gparent;
					continue;
				} else {
					if(parent.isLeftChild()) {
						// Parent is red, uncle is black, parent is
						// left child
						if(!n.isLeftChild()) {
							n = parent;
							leftRotate(n);
						}
						n = n.parent;
						n.red = false;
						n.parent.red = true;
						rightRotate(n.parent);
					} else {
						// Parent is red, uncle is black, parent is
						// right child.
						if(!n.isRightChild()) {
							n = parent;
							rightRotate(n);
						}
						n = n.parent;
						n.red = false;
						n.parent.red = true;
						leftRotate(n.parent);
					}
				}
				break;
			}
		}
	}
	
	public RedBlackNode addPage(Pool p) {
		RedBlackNode n = (RedBlackNode)p.alloc();
		if(n != null) {
			pages_.push_back(n);
		}
		return n;
	}
}
