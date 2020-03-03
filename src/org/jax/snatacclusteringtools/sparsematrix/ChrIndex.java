package org.jax.snatacclusteringtools.sparsematrix;


public class ChrIndex {

	private String _chr;
	private int _index;
	
	public ChrIndex(String chr, int index) {
		_chr = chr;
		_index = index;
	}
	
	public String getChr() {
		return _chr;
	}
	
	public int getIndex() {
		return _index;
	}
}
