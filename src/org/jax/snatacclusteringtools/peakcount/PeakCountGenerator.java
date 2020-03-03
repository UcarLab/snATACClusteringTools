package org.jax.snatacclusteringtools.peakcount;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map.Entry;

import org.jax.snatacclusteringtools.util.Util;

import java.util.TreeMap;
import java.util.TreeSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class PeakCountGenerator {
	
	public static void main(String[] args) {
		if(args.length != 5) {
			System.out.println("Usage generatepeakcounts bamfile barcodes summits extend outdir");
			System.exit(0);
		}
		String bamfile = args[0];
		String barcodes = args[1];
		String peaks = args[2];
		int extend = Integer.parseInt(args[3]);
		String outdir = args[4];
		PeakCountGenerator rip = new PeakCountGenerator();
		try {
			rip.generatePeakCounts(bamfile, barcodes, peaks, extend, outdir);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}	

	@SuppressWarnings("unchecked")
	public void generatePeakCounts(String bamfile, String cellbarcodes, String peaks, int extend, String outdir) throws IOException {
		
		SamReaderFactory factory = SamReaderFactory.makeDefault()
	              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	              .validationStringency(ValidationStringency.SILENT);
		
		final SamReader reader = factory.open(new File(bamfile));

		if(!(reader.getFileHeader().getSortOrder().getComparatorInstance() instanceof SAMRecordCoordinateComparator)) {
			System.out.println("The input BAM file must be coordinate sorted.");
			reader.close();
			System.exit(0);
		}
		
		SAMRecordIterator it = reader.iterator();
		Util u = new Util();
		TreeMap<String, String> cellbarcodemap = u.readCellBarcodes(cellbarcodes);
		TreeMap<String, Integer> cellidtoint = new TreeMap<String, Integer>();
		
		TreeMap<String, int[]> peaksummits = readPeakSummits(peaks);
		
		TreeMap<String, int[]> peakcounts = new TreeMap<String, int[]>();
		TreeMap<String, TreeSet<Integer>[]> uniquecellcounts = new TreeMap<String, TreeSet<Integer>[]>();
		TreeMap<String, TreeMap<String, int[]>> peakcountsbycell = new TreeMap<String, TreeMap<String, int[]>>();
		
		int uniquecellcount = 0;
		System.out.println("Reading BAM file.");
		while(it.hasNext()){
			SAMRecord next = it.next();
			boolean isnegative = next.getReadNegativeStrandFlag();
			int insertsize = next.getInferredInsertSize();

			if(!isnegative && u.isValidRead(next) && insertsize > 0 && insertsize < 900) {
				
				String chr = next.getReferenceName();
				if(!peaksummits.containsKey(chr)) {
					continue;
				}
					
				String barcode = u.getBarcode(next);
				if(barcode == null) {
					continue;
				}
					
				String cellid = cellbarcodemap.get(barcode);
				if(cellid == null) {
					continue;
				}
				
				
				
				if(!peakcountsbycell.containsKey(cellid)) {
					peakcountsbycell.put(cellid, new TreeMap<String, int[]>());
					cellidtoint.put(cellid, uniquecellcount++);
				}
				
				int intcellid = cellidtoint.get(cellid);
				
				TreeMap<String, int[]> cellpeakcounts = peakcountsbycell.get(cellid);
					
				int[] sortedsummits = peaksummits.get(chr);
				if(!peakcounts.containsKey(chr)) {
					peakcounts.put(chr, new int[sortedsummits.length]);
					uniquecellcounts.put(chr, new TreeSet[sortedsummits.length]);
				}
				
				if(!cellpeakcounts.containsKey(chr)) {
					cellpeakcounts.put(chr, new int[sortedsummits.length]);
				}
				
				int start = next.getAlignmentStart();
				int end = start+insertsize-1;
				
				int[] range = getOverlappingPeakRange(start, end, sortedsummits, extend);
				int[] curcounts = peakcounts.get(chr);
				int[] curcellcounts = peakcountsbycell.get(cellid).get(chr);
				
				TreeSet<Integer>[] uniquetrees = uniquecellcounts.get(chr);
				
				for(int i = range[0]; i <= range[1]; i++) {
					curcounts[i]++;
					curcellcounts[i]++;
					
					TreeSet<Integer> curtree = uniquetrees[i];
					if(curtree == null) {
						curtree = new TreeSet<Integer>();
						uniquetrees[i] = curtree;
					}
					curtree.add(intcellid);
				}
			}
			
		}
		
		it.close();
		reader.close();
		
		System.out.println("Writing peak overlaps");
		writeCellInPeakCounts(peakcountsbycell, extend, outdir);
		writeTotalCounts(peakcounts, uniquecellcounts, extend, outdir);
		
		
		StringBuilder output = new StringBuilder();

		output.append("Unique Cell Count: ");
		output.append(Integer.toString(uniquecellcount));
		output.append("\n");

		u.writeCompletionFile(outdir, "CompletionSummary", output.toString());
		
		
		System.out.println("Complete.");

	}
	
	private void writeCellInPeakCounts(TreeMap<String, TreeMap<String, int[]>> cellcounts, int extend, String outdir) throws IOException {
		BufferedWriter curbw = new BufferedWriter(new FileWriter(outdir+"/PeakCounts_"+extend+".txt"));
		
		curbw.write("Chromosome");
		curbw.write("\t");
		curbw.write("Peak Index ("+Integer.toString(extend)+")");
		curbw.write("\t");
		curbw.write("# Inserts Overlapping Peak");
		curbw.write("\t");
		curbw.write("Cell Id");
		curbw.write("\n");
		
		for(Iterator<Entry<String, TreeMap<String, int[]>>> it = cellcounts.entrySet().iterator(); it.hasNext();) {
			Entry<String, TreeMap<String, int[]>> next = it.next();
			String cellid = next.getKey();
			for(Iterator<Entry<String, int[]>> countit = next.getValue().entrySet().iterator(); countit.hasNext();) {
				Entry<String, int[]> curcounts = countit.next();
				String curchr = curcounts.getKey();
				int[] curpeakcounts = curcounts.getValue();
				for(int i = 0; i < curpeakcounts.length; i++) {
					if(curpeakcounts[i] > 0) {
						curbw.write(curchr);
						curbw.write("\t");
						curbw.write(Integer.toString(i));
						curbw.write("\t");
						curbw.write(Integer.toString(curpeakcounts[i]));
						curbw.write("\t");
						curbw.write(cellid);
						curbw.write("\n");
					}
				}

			}
		}
		
		curbw.flush();
		curbw.close();
	}
	
	private void writeTotalCounts(TreeMap<String, int[]> totalcounts, TreeMap<String, TreeSet<Integer>[]> uniquecounts, int extend, String outdir) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outdir+"/TotalPeakCounts_"+extend+".txt"));
		
		bw.write("Chromosome");
		bw.write("\t");
		bw.write("Peak Index ("+Integer.toString(extend)+")");
		bw.write("\t");
		bw.write("# Inserts Overlapping Peak");
		bw.write("\t");
		bw.write("# of Unique Cells");
		bw.write("\n");
		
		for(Iterator<Entry<String, int[]>> it = totalcounts.entrySet().iterator(); it.hasNext();) {
			Entry<String, int[]> next = it.next();
			String curchr = next.getKey();
			int[] curtotalcounts = next.getValue();
			TreeSet<Integer>[] curu = uniquecounts.get(curchr);
			
			for(int i = 0;i < curtotalcounts.length; i++) {
				int numunique = 0;
				if(curu[i] != null) {
					numunique = curu[i].size();
				}
				if(curtotalcounts[i] > 0 && numunique > 0) {
					bw.write(curchr);
					bw.write("\t");
					bw.write(Integer.toString(i));
					bw.write("\t");
					bw.write(Integer.toString(curtotalcounts[i]));
					bw.write("\t");
					bw.write(Integer.toString(numunique));
					bw.write("\n");
				}
			}
		}
		
		bw.flush();
		bw.close();
	}
	
	
	private int[] getOverlappingPeakRange(int fposition, int rposition, int[] sortedsummits, int extend) {
		int start = fposition-extend;
		int end = rposition+extend;
		int overlapstart = getClosestIndex(start, sortedsummits)[0];
		int overlapend = getClosestIndex(end, sortedsummits)[1];

		for(int i = overlapstart; i < overlapend; i++) {
			if(start >= sortedsummits[i] && start <= sortedsummits[i+1]) {
				if(start == sortedsummits[i]) {
					overlapstart = i;
				}
				else {
					overlapstart = i+1;
				}
				break;
			}
		}
		
		for(int i = overlapend; i > overlapstart; i--) {
			if(end <= sortedsummits[i] && end >= sortedsummits[i-1]) {
				if(end == sortedsummits[i]) {
					overlapend = i;
				}
				else {
					overlapend = i-1;
				}
				break;
			}
		}
		
		if(overlapstart == overlapend) {
			if(end < sortedsummits[overlapstart] || start > sortedsummits[overlapstart]) {
				overlapstart = 1;
				overlapend = 0;
			}
		}
		
		return new int[] { overlapstart, overlapend };
	}
	
	private int[] getClosestIndex(int position, int[] sortedpositions) {
		int mini = 0;
		int maxi = sortedpositions.length-1;
		
		if(position < sortedpositions[mini]) {
			return new int[] { mini, mini };
		}
		
		if(position > sortedpositions[maxi]) {
			return new int[] { maxi, maxi };
		}
		
		
		while(maxi-mini > 1) {
			int curi = mini+(maxi-mini)/2;
			int midpos = sortedpositions[curi];
			if(position < midpos) {
				maxi = curi;
			}
			else if(position > midpos) {
				mini = curi;
			}
			else {
				mini = curi;
				maxi = curi;
				break;
			}
		}
		
		if(maxi-mini == 1) {
			if(position == sortedpositions[mini]) {
				maxi = mini;
			}
			else if(position == sortedpositions[maxi]) {
				mini = maxi;
			}
			
		}
		
		return new int[] {mini, maxi};
	}
	
	private TreeMap<String, int[]> readPeakSummits(String file) throws IOException{
		TreeMap<String, TreeSet<Integer>> tree = new TreeMap<String, TreeSet<Integer>>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		while(br.ready()) {
			String[] split = br.readLine().split("\t");
			String chr = split[0];
			int position = Integer.parseInt(split[1]);
			
			if(!tree.containsKey(chr)) {
				tree.put(chr, new TreeSet<Integer>());
			}
			tree.get(chr).add(position);
		}
		br.close();
		
		TreeMap<String, int[]> rv = new TreeMap<String, int[]>();
		
		for(Iterator<Entry<String, TreeSet<Integer>>> it = tree.entrySet().iterator(); it.hasNext();) {
			Entry<String, TreeSet<Integer>> next = it.next();
			rv.put(next.getKey(), getSortedArrayFromTree(next.getValue()));
		}
		
		return rv;
	}
	
	private int[] getSortedArrayFromTree(TreeSet<Integer> treevals) {
		int[] rv = new int[treevals.size()];
		int index = 0;
		for(Iterator<Integer> it = treevals.iterator(); it.hasNext();) {
			rv[index++] = it.next();
		}
		return rv;
	}
	

	
}
