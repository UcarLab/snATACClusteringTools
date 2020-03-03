package org.jax.snatacclusteringtools.peakprocessing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;

import org.jax.snatacclusteringtools.util.Util;

import java.util.TreeMap;
import java.util.TreeSet;

public class SummitPeakMerger {
	
	//TODO Implement an alternative method that uses the union of all summits instead of  "merging" summits
	public static void main(String[] args) {
		if(args.length != 5) {
			System.out.println("SummitPeakMerger.jar fileofpeakfiles chromosomesizelist extendamt outdir outfileprefix");
			System.exit(0);
		}
		SummitPeakMerger spm = new SummitPeakMerger();
		int extend = Integer.parseInt(args[2]);
		spm.processPeaks(args[0], args[1], extend, args[3], args[4]);
	}
	
	public void processPeaks(String in, String chromsizelist, int extend, String outdir, String prefix) {
		try {
			Util u = new Util();
			String[] files = u.getFiles(in);
			TreeMap<String, Integer> chromsizes = u.readChromSizes(chromsizelist);
			mergeSummits(files, chromsizes, extend, outdir, prefix);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void writeSummits(TreeMap<String, TreeMap<Integer, Double>> mergedpeaks, int extend, String output) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		for(Iterator<Entry<String, TreeMap<Integer, Double>>> chrit = mergedpeaks.entrySet().iterator(); chrit.hasNext();) {
			Entry<String, TreeMap<Integer, Double>> next = chrit.next();
			String chr = next.getKey();
			TreeMap<Integer, Double> positions = next.getValue();
			int index = 1;
			for(Iterator<Entry<Integer, Double>> posit = positions.entrySet().iterator(); posit.hasNext();) {
				Entry<Integer, Double> posnext = posit.next();
				int summit = posnext.getKey();
				double score = posnext.getValue();
				bw.write(chr);
				bw.write("\t");
				bw.write(Integer.toString(summit));
				bw.write("\t");
				bw.write(Integer.toString(summit+1));
				bw.write("\t");
				bw.write(Integer.toString(index++));
				bw.write("\t");
				bw.write(Double.toString(1/score)); //Inverse for summit scores to maintain higher is better to reuse this method
				bw.write("\n");
			}
			
		}
		
		bw.flush();
		bw.close();
	}
	
	private void writePeaks(TreeMap<String, TreeMap<Integer, Double>> mergedpeaks, int extend, String output) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		for(Iterator<Entry<String, TreeMap<Integer, Double>>> chrit = mergedpeaks.entrySet().iterator(); chrit.hasNext();) {
			Entry<String, TreeMap<Integer, Double>> next = chrit.next();
			String chr = next.getKey();
			TreeMap<Integer, Double> positions = next.getValue();
			int index = 1;
			for(Iterator<Entry<Integer, Double>> posit = positions.entrySet().iterator(); posit.hasNext();) {
				Entry<Integer, Double> posnext = posit.next();
				int summit = posnext.getKey();
				double score = posnext.getValue();
				bw.write(chr);
				bw.write("\t");
				bw.write(Integer.toString(summit-extend));
				bw.write("\t");
				bw.write(Integer.toString(summit+extend));
				bw.write("\t");
				bw.write(Integer.toString(index++));
				bw.write("\t");
				bw.write(Double.toString(score));
				bw.write("\n");
			}
			
		}
		
		bw.flush();
		bw.close();
	}

	
	public void mergeSummits(String[] files, TreeMap<String, Integer> chromsizes, int extend, String outdir, String prefix) throws IOException {
		@SuppressWarnings("unchecked")
		TreeMap<String, TreeMap<Integer, Double>>[] allmergedsumits = new TreeMap[files.length];
		
		for(int i = 0; i < files.length; i++) {
			TreeMap<String, TreeMap<Integer, Double>> summits = readSummits(files[i], chromsizes);
			TreeMap<String, TreeMap<Integer, Double>> mergedsummits = mergeSummits(summits, extend);
			TreeMap<String, TreeMap<Integer, Double>> quantileadjusted = getQuantileAdjustedSummits(mergedsummits);
			allmergedsumits[i] = quantileadjusted;
			//writeSummits(quantileadjusted, extend, files[i]+"_merged.bed");
			//writePeaks(quantileadjusted, extend, files[i]+"_merged_peaks.bed");
		}
		
		if(files.length > 0) {
			TreeMap<String, TreeMap<Integer, Double>> combinedsummits = combineQuantileSummits(allmergedsumits);
			TreeMap<String, TreeMap<Integer, Double>> mergedcombinedsummits = mergeSummits(combinedsummits, extend);
			writeSummits(mergedcombinedsummits, extend, outdir+"/"+prefix+"_summits.bed");
			writePeaks(mergedcombinedsummits, extend, outdir+"/"+prefix+"_peaks.bed");
		}
	}
	

	private TreeMap<String, TreeMap<Integer, Double>> combineQuantileSummits(TreeMap<String, TreeMap<Integer, Double>>[] allmergedsummits) {
		TreeMap<String, TreeMap<Integer, Double>> rv = new TreeMap<String, TreeMap<Integer, Double>>();
		TreeSet<String> chromosomes = getChromosomes(allmergedsummits);
		for(Iterator<String> chrit = chromosomes.iterator(); chrit.hasNext();) {
			String curchr = chrit.next();
			TreeMap<Integer, Double> curtree = new TreeMap<Integer, Double>();
			rv.put(curchr, curtree);
			
			for(int i = 0; i < allmergedsummits.length; i++) {
				TreeMap<Integer, Double> curmergedtree = allmergedsummits[i].get(curchr);
				if(curmergedtree != null) {
					for(Iterator<Entry<Integer, Double>> eit = curmergedtree.entrySet().iterator(); eit.hasNext();) {
						Entry<Integer, Double> enext = eit.next();
						int position = enext.getKey();
						double quantilescore = enext.getValue();
						if(curtree.containsKey(position)) {
							curtree.put(position, Math.max(curtree.get(position), quantilescore));
						}
						else {
							curtree.put(position, quantilescore);
						}
					}
				}
			}
		}
		return rv;
	}
	
	private TreeSet<String> getChromosomes(TreeMap<String, TreeMap<Integer, Double>>[] allmergedsummits) {
		TreeSet<String> rv = new TreeSet<String>();
		for(int i = 0; i < allmergedsummits.length; i++) {
			rv.addAll(allmergedsummits[i].keySet());
		}
		return rv;
	}
	
	
	private TreeMap<String, TreeMap<Integer, Double>> getQuantileAdjustedSummits(TreeMap<String, TreeMap<Integer, Double>> mergedpeaks) {
		TreeMap<Double, Double> adjustedmapping = getAdjustedMapping(mergedpeaks);
		TreeMap<String, TreeMap<Integer, Double>> rv = new TreeMap<String, TreeMap<Integer, Double>>();
		for(Iterator<Entry<String, TreeMap<Integer, Double>>> chrentryit = mergedpeaks.entrySet().iterator(); chrentryit.hasNext();) {
			Entry<String, TreeMap<Integer, Double>> chrentry = chrentryit.next();
			String chr = chrentry.getKey();
			TreeMap<Integer, Double> chrvalues = new TreeMap<Integer, Double>();
			rv.put(chr, chrvalues);
			for(Iterator<Entry<Integer, Double>> eit = chrentry.getValue().entrySet().iterator(); eit.hasNext();) {
				Entry<Integer, Double> nextentry = eit.next();
				chrvalues.put(nextentry.getKey(), adjustedmapping.get(nextentry.getValue()));
			}
		}
		return rv;
	}
	
	//Scores assume that least score is better
	private	TreeMap<Double, Double> getAdjustedMapping(TreeMap<String, TreeMap<Integer, Double>> mergedpeaks){
		LinkedList<Double> allscores = new LinkedList<Double>();
		for(Iterator<TreeMap<Integer, Double>> it = mergedpeaks.values().iterator(); it.hasNext();) {
			TreeMap<Integer, Double> next = it.next();
			for(Iterator<Double> scoreit = next.values().iterator(); scoreit.hasNext();) {
				allscores.add(scoreit.next());
			}
		}
		allscores.sort(new Comparator<Double>() {
			@Override
			public int compare(Double o1, Double o2) {
				if(o1 < o2) {
					return -1;
				}
				else if(o1 > o2){
					return 1;
				}
				return 0;
			}});
		
		TreeMap<Double, Double> rv = new TreeMap<Double, Double>();
		int currank = 1;
		int totalscores = allscores.size();
		for(Iterator<Double> it = allscores.iterator(); it.hasNext();) {
			Double next = it.next();
			if(rv.containsKey(next)) {
				currank++;
			}
			else {
				double newscore = (currank++)/(double)totalscores;
				rv.put(next, newscore);
			}
		}
		return rv;
	}
	

	
	private TreeMap<String, TreeMap<Integer, Double>> readSummits(String file, TreeMap<String, Integer> chromosomes) throws IOException{
		TreeMap<String, TreeMap<Integer, Double>> rv = new TreeMap<String, TreeMap<Integer, Double>>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		while(br.ready()) {
			String[] split = br.readLine().split("\t");
			String chr = split[0];
			int position = Integer.parseInt(split[1]);
			double score = Double.parseDouble(split[4]);
			
			if(!chromosomes.containsKey(chr)) {
				continue;
			}
			
			if(!rv.containsKey(chr)) {
				rv.put(chr, new TreeMap<Integer, Double>());
			}
			TreeMap<Integer, Double> chrposition = rv.get(chr);
			
			if(chrposition.containsKey(position)) {
				double prevscore = chrposition.get(position);
				if(prevscore < score) {
					chrposition.put(position, score);
				}
			}
			else {
				chrposition.put(position, score);
			}
		}
		br.close();
		return rv;
	}
	
	
	private TreeMap<String, TreeMap<Integer, Double>> mergeSummits(TreeMap<String, TreeMap<Integer, Double>> summits, int extend) {
		TreeMap<String, TreeMap<Integer, Double>> rv = new TreeMap<String, TreeMap<Integer, Double>>();
		for(Iterator<String> chrit = summits.keySet().iterator(); chrit.hasNext();) {
			String curchr = chrit.next();
			if(summits.containsKey(curchr)) {
				rv.put(curchr, mergeChrSummits(summits.get(curchr), extend));
			}
		}
		return rv;
	}
	
	private TreeMap<Integer, Double> mergeChrSummits(TreeMap<Integer, Double> summits, int extend) {
		int[] sortedsummits = getSortedSummits(summits);
		TreeMap<Integer, Integer> sortedindex = getSortedSummitIndex(sortedsummits);
		//Higher score is better
		TreeMap<Double, LinkedList<Integer>> summitscores = getSummitsByScore(summits, true);
		
		TreeSet<Integer> completed = new TreeSet<Integer>();
		TreeMap<Integer, Double> finalsummits = new TreeMap<Integer, Double>();

		
		for(Iterator<Entry<Double, LinkedList<Integer>>> it = summitscores.entrySet().iterator(); it.hasNext();) {
			Entry<Double, LinkedList<Integer>> next = it.next();
			LinkedList<Integer> positions = next.getValue();
			double score = -next.getKey(); //Convert the negative score back to positive
			for(Iterator<Integer> posit = positions.iterator(); posit.hasNext();) {
				int nextpos = posit.next();
				if(!completed.contains(nextpos)) {
					
					//1. Scan forward and backwards and mark those within the extend amount as completed
					int sp = sortedindex.get(nextpos);//start position
					
					//Forward scan
					int fcutoff = nextpos + extend;
					for(int i = sp; i < sortedsummits.length; i++) {
						int fsp = sortedsummits[i];
						if(fsp > fcutoff) {
							break;
						}
						completed.add(fsp);
					}
					
					//Backward scan
					int bcutoff = nextpos - extend;
					for(int i = sp; i > -1; i--) {
						int bsp = sortedsummits[i];
						if(bsp < bcutoff) {
							break;
						}
						completed.add(bsp);
					}
					
					//2. add current position to final summits list and mark as complete
					finalsummits.put(nextpos, score);
					completed.add(nextpos);
				}
			}
			
		}
		
		return finalsummits;
	}
	
	
	//Note, makes the scores negative so that minimum score is the best score in the tree
	private TreeMap<Double, LinkedList<Integer>> getSummitsByScore(TreeMap<Integer, Double> summits, boolean higherisbetter){
		TreeMap<Double, LinkedList<Integer>> rv = new TreeMap<Double, LinkedList<Integer>>();
		for(Iterator<Entry<Integer, Double>> it = summits.entrySet().iterator(); it.hasNext();) {
			Entry<Integer, Double> next = it.next();
			int summit = next.getKey();
			double qval = next.getValue();
			if(higherisbetter) {
				qval = -qval;
			}
			if(!rv.containsKey(qval)) {
				rv.put(qval, new LinkedList<Integer>());
			}
			rv.get(qval).add(summit);
		}
		return rv;
	}
	
	private TreeMap<Integer, Integer> getSortedSummitIndex(int[] sortedsummits){
		TreeMap<Integer, Integer> rv = new TreeMap<Integer, Integer>();
		for(int i = 0; i < sortedsummits.length; i++) {
			rv.put(sortedsummits[i], i);
		}
		return rv;
	}
	
	private int[] getSortedSummits(TreeMap<Integer, Double> summits) {
		int[] rv = new int[summits.size()];
		int index = 0;
		for(Iterator<Integer> it = summits.keySet().iterator(); it.hasNext();) {
			rv[index++] = it.next();
		}
		return rv;
	}
	

}
