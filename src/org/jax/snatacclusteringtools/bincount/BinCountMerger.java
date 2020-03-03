package org.jax.snatacclusteringtools.bincount;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;

import org.jax.snatacclusteringtools.util.Util;

import java.util.TreeMap;

public class BinCountMerger {
	
	public static void main(String[] args) {
		if(args.length != 4) {
			System.out.println("Usage: mergebincounts directorylist chromosomesizes binsize outdir");
			System.exit(0);
		}
		String directorylist = args[0];
		String chromsizes = args[1];
		int binsize = Integer.parseInt(args[2]);
		String outdir = args[3];
		BinCountMerger pm = new BinCountMerger();
		pm.mergePartitions(directorylist, chromsizes, binsize, outdir);
	}
	
	public void mergePartitions(String directorylist, String chromlist, int binsize, String outdir) {
		Util u = new Util();
		try {
			String[] directories = u.getFiles(directorylist);
			String[] samplenames = u.inferSampleNames(directories);
			TreeMap<String, Integer> chromosomes = u.readChromSizes(chromlist);
			for(Iterator<String> it = chromosomes.keySet().iterator(); it.hasNext();) {
				String curchr = it.next();
				mergeBinCounts(directories, samplenames, curchr, binsize, outdir);
			}
			
			u.writeCompletionFile(outdir, "CompletionSummary", "Complete!");

			
		} catch (IOException e) {
			try {
				u.writeErrorFile(e, outdir);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			e.printStackTrace();
		}
	}
	
	
	private void mergeBinCounts(String[] directories, String[] samplenames, String chr, int binsize, String outdir) throws IOException {
		LinkedList<String> totalchrfiles = new LinkedList<String>();
		File mergefile = new File(outdir+"/MergedBinCounts_"+chr+"_"+binsize+".txt");
		BufferedWriter bw = new BufferedWriter(new FileWriter(mergefile));
		bw.write("Chromosome");
		bw.write("\t");
		bw.write("Bin Index ("+binsize+")");
		bw.write("\t");
		bw.write("# Inserts Overlapping Bin");
		bw.write("\t");
		bw.write("Sample_Cell Id");
		bw.write("\n");
		for(int i = 0; i < directories.length; i++) {
			String curfile = directories[i]+"/BinCounts_"+chr+"_"+binsize+".txt";
			if((new File(curfile)).exists()){
				String sampleid = samplenames[i]+"_";
				appendBinCountChrData(bw, sampleid, curfile);
				totalchrfiles.add(directories[i]+"/TotalBinCounts_"+chr+"_"+binsize+".txt");
			}
			else {
				System.out.println("File does not exist. Skipping. "+curfile);
			}
		}
		bw.flush();
		bw.close();
		if(totalchrfiles.isEmpty()) {
			mergefile.delete();
		}
		else {
			writeTotalBinCounts(totalchrfiles.toArray(new String[0]), chr, binsize, outdir+"/TotalMergedBinCounts_"+chr+"_"+binsize+".txt");
		}
	}
	
	private void writeTotalBinCounts(String[] files, String chr, int binsize, String outfile) throws IOException {
		if(files.length == 0) {
			return;
		}
		TreeMap<Integer, int[]> totalcounts = getTotalCounts(files);
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		
		bw.write("Chromosome");
		bw.write("\t");
		bw.write("Bin Index ("+binsize+")");
		bw.write("\t");
		bw.write("# Inserts Overlapping Bin");
		bw.write("\t");
		bw.write("# of Unique Cells");
		bw.write("\n");
		for(Iterator<Entry<Integer, int[]>> it = totalcounts.entrySet().iterator(); it.hasNext();) {
			Entry<Integer, int[]> next = it.next();
			int bin = next.getKey();
			int[] counts = next.getValue();
			bw.write(chr);
			bw.write("\t");
			bw.write(Integer.toString(bin));
			bw.write("\t");
			bw.write(Integer.toString(counts[0]));
			bw.write("\t");
			bw.write(Integer.toString(counts[1]));
			bw.write("\n");
		}
		
		bw.flush();
		bw.close();
	}
	
	private TreeMap<Integer, int[]> getTotalCounts(String[] files) throws IOException {
		TreeMap<Integer, int[]> totalcounts = new TreeMap<Integer, int[]>();
		for(int i = 0; i < files.length; i++) {
			BufferedReader br = new BufferedReader(new FileReader(files[i]));
			br.readLine();
			while(br.ready()) {
				String[] split = br.readLine().split("\t");
				int bin = Integer.parseInt(split[1]);
				int numinserts = Integer.parseInt(split[2]);
				int numcells = Integer.parseInt(split[3]);
				if(!totalcounts.containsKey(bin)) {
					totalcounts.put(bin, new int[2]);
				}
				int[] a = totalcounts.get(bin);
				a[0] += numinserts;
				a[1] += numcells;
			}
			br.close();
		}
		return totalcounts;
	}
	
	private void appendBinCountChrData(BufferedWriter bw, String sampleid, String file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		br.readLine(); //Skip header
		while(br.ready()) {
			String[] split = br.readLine().split("\t");
			String chr = split[0];
			String bin = split[1];
			String numinserts = split[2];
			String cellid = split[3];
			
			bw.write(chr);
			bw.write("\t");
			bw.write(bin);
			bw.write("\t");
			bw.write(numinserts);
			bw.write("\t");
			bw.write(sampleid+cellid);
			bw.write("\n");
		}
		br.close();
	}
	
}
