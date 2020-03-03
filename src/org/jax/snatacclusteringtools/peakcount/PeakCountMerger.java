package org.jax.snatacclusteringtools.peakcount;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.jax.snatacclusteringtools.util.Util;

public class PeakCountMerger {
	
	public static void main(String[] args) {
		if(args.length != 4) {
			System.out.println("Usage: mergepeakcounts filelist totalfilelist peakextendsize outdir");
			System.exit(0);
		}
		String filelist = args[0];
		String totalfilelist = args[1];

		int binsize = Integer.parseInt(args[2]);
		String outdir = args[3];
		PeakCountMerger pm = new PeakCountMerger();
		pm.mergePeakCounts(filelist, totalfilelist, binsize, outdir);
	}
	
	public void mergePeakCounts(String filelist, String totalfilelist, int binsize, String outdir) {
		try {
			Util u = new Util();
			String[] files = u.getFiles(filelist);
			String[] totalfiles = u.getFiles(totalfilelist);
			String[] samplenames = u.inferSampleNames(totalfiles);

			mergeCounts(files, totalfiles, samplenames, binsize, outdir);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	private void mergeCounts(String[] files, String[] totalfiles, String[] samplenames, int extend, String outdir) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outdir+"/MergedPeakCounts_"+extend+".txt"));
		bw.write("Chromosome");
		bw.write("\t");
		bw.write("Peak Index ("+extend+")");
		bw.write("\t");
		bw.write("# Inserts Overlapping Bin");
		bw.write("\t");
		bw.write("Sample_Cell Id");
		bw.write("\n");
		for(int i = 0; i < files.length; i++) {
			String sampleid = samplenames[i]+"_";
			String curfile = files[i];
			appendPartitionData(bw, sampleid, curfile);
		}
		bw.flush();
		bw.close();
		writeTotalParitionFile(totalfiles, extend, outdir+"/TotalMergedPeakCounts_"+extend+".txt");
	}
	
	private void writeTotalParitionFile(String[] files,int extend, String outfile) throws IOException {
		TreeMap<String, int[]> totalcounts = getTotalCounts(files);
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		
		bw.write("Chromosome");
		bw.write("\t");
		bw.write("Peak Index ("+extend+")");
		bw.write("\t");
		bw.write("# Inserts Overlapping Bin");
		bw.write("\t");
		bw.write("# of Unique Cells");
		bw.write("\n");
		for(Iterator<Entry<String, int[]>> it = totalcounts.entrySet().iterator(); it.hasNext();) {
			Entry<String, int[]> next = it.next();
			String key = next.getKey();
			int[] counts = next.getValue();
			String[] split = key.split(":");
			bw.write(split[0]);
			bw.write("\t");
			bw.write(split[1]);
			bw.write("\t");
			bw.write(Integer.toString(counts[0]));
			bw.write("\t");
			bw.write(Integer.toString(counts[1]));
			bw.write("\n");
		}
		
		bw.flush();
		bw.close();
	}
	
	private TreeMap<String, int[]> getTotalCounts(String[] files) throws IOException {
		TreeMap<String, int[]> totalcounts = new TreeMap<String, int[]>();
		for(int i = 0; i < files.length; i++) {
			BufferedReader br = new BufferedReader(new FileReader(files[i]));
			br.readLine();
			while(br.ready()) {
				String[] split = br.readLine().split("\t");
				String chr = split[0];
				int bin = Integer.parseInt(split[1]);
				int numinserts = Integer.parseInt(split[2]);
				int numcells = Integer.parseInt(split[3]);
				String key = chr+":"+bin;
				if(!totalcounts.containsKey(key)) {
					totalcounts.put(key, new int[2]);
				}
				int[] a = totalcounts.get(key);
				a[0] += numinserts;
				a[1] += numcells;
			}
			br.close();
		}
		return totalcounts;
	}
	
	private void appendPartitionData(BufferedWriter bw, String sampleid, String file) throws IOException {
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
