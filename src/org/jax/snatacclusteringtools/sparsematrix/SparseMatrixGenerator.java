package org.jax.snatacclusteringtools.sparsematrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.TreeMap;

import org.jax.snatacclusteringtools.util.Util;

public class SparseMatrixGenerator {
	
	public static void main(String[] args) {

		if(args.length != 4) {
			System.out.println("Usage: generatesparsematrix countfiles totalcountfiles topn outfile");
			System.exit(0);
		}
		String chrfilelist = args[0];
		String totalfilelist = args[1];
		int n = Integer.parseInt(args[2]);
		String output = args[3];
		SparseMatrixGenerator smg = new SparseMatrixGenerator();
		smg.generateSparseMatrix(chrfilelist, totalfilelist, n, output);
	}
	
	public void generateSparseMatrix(String chrfilelist, String totalfilelist, int n, String output) {
		try {
			Util u = new Util();
			String[] chrfiles = u.getFiles(chrfilelist);
			String[] totalfiles = u.getFiles(totalfilelist);
			TreeMap<String, Integer>  topchrindices = getTopChrIndices(totalfiles, n);
			writeSparseMatrix(chrfiles, topchrindices, output);

		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}

	private void writeSparseMatrix(String[] chrfiles, TreeMap<String, Integer> topregions, String outfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write("Cell ID");
		bw.write("\t");
		bw.write("Row Index");
		bw.write("\t");
		bw.write("Column Index");
		bw.write("\t");
		bw.write("Count");
		bw.write("\n");
		TreeMap<String, Integer> cellpositions = new TreeMap<String, Integer>(); //Rows
		int newrowid = 0;
		for(int i = 0; i < chrfiles.length; i++) {
			BufferedReader br = new BufferedReader(new FileReader(chrfiles[i]));
			br.readLine(); //Header
			while(br.ready()) {
				String[] split = br.readLine().split("\t");
				String chr = split[0];
				int index = Integer.parseInt(split[1]);
				int count = Integer.parseInt(split[2]);
				String cellid = split[3];
				
				if(!cellpositions.containsKey(cellid)) {
					cellpositions.put(cellid, newrowid++);
				}
				
				int currowid = cellpositions.get(cellid);
				
				String colkey = chr+":"+index;
				if(topregions.containsKey(colkey)) {
					int curcolid = topregions.get(colkey);
					if(count > 0) {
						bw.write(cellid);
						bw.write("\t");
						bw.write(Integer.toString(currowid));
						bw.write("\t");
						bw.write(Integer.toString(curcolid));
						bw.write("\t");
						bw.write(Integer.toString(count));
						bw.write("\n");
					}
				}
				
			}
			
			br.close();
		}
		
		bw.flush();
		bw.close();
		
	}
	
	private TreeMap<String, Integer> getTopChrIndices(String[] totalfiles, int n) throws IOException {
		
		int totaladded = 0;
		int smallestn = Integer.MAX_VALUE;
		TreeMap<Integer, LinkedList<ChrIndex>> tracker = new TreeMap<Integer, LinkedList<ChrIndex>>();
		
		for(int i = 0; i < totalfiles.length; i++) {
			
			BufferedReader br = new BufferedReader(new FileReader(totalfiles[i]));
			br.readLine(); //Header
			while(br.ready()) {
				String[] split = br.readLine().split("\t");
				String chr = split[0];
				int index = Integer.parseInt(split[1]);
				int totalcount = Integer.parseInt(split[2]);
				//int uniquecount = Integer.parseInt(split[3]);
				
				if(totaladded < n) {
					if(!tracker.containsKey(totalcount)) {
						tracker.put(totalcount, new LinkedList<ChrIndex>());
					}
					tracker.get(totalcount).add(new ChrIndex(chr, index));
					smallestn = Math.min(totalcount, smallestn);
					totaladded++;
				}
				else {
					//Only add if it beats the smallest value in the top N
					if(totalcount == smallestn) {
						tracker.get(totalcount).add(new ChrIndex(chr, index)); //Might go over the top due to ties
						totaladded++;
					}
					else if(totalcount > smallestn) {
						int totalremoved = tracker.get(smallestn).size();
						//Only remove if we are closer to n by doing so
						if((smallestn - totalremoved) + 1 >= n) {
							tracker.remove(smallestn);
							totaladded -= totalremoved;
							
						}
						
						if(!tracker.containsKey(totalcount)) {
							tracker.put(totalcount, new LinkedList<ChrIndex>());
						}
						tracker.get(totalcount).add(new ChrIndex(chr, index));
						totaladded++;
						
						smallestn = tracker.firstKey();
					}
				}
			}
			
			br.close();
			
		}
		
		TreeMap<String, Integer> rv = new TreeMap<String, Integer>();
		int returncount = totaladded;
		for(Iterator<LinkedList<ChrIndex>> it = tracker.values().iterator(); it.hasNext();) {
			LinkedList<ChrIndex> next = it.next();
			for(Iterator<ChrIndex> it2 = next.iterator(); it2.hasNext();) {
				ChrIndex curidx = it2.next();
				if(returncount <= n) {
					rv.put(curidx.getChr()+":"+curidx.getIndex(), returncount);
				}
				returncount--;
			}
		}
		
		
		return rv;
	}

}
