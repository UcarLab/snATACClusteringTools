package org.jax.snatacclusteringtools.bincount;

import java.io.BufferedWriter;
import java.io.File;
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

public class BinCountGenerator {
	
	public static void main(String[] args) {
		if(args.length != 5) {
			System.out.println("Usage: generatebincounts bamfile cellbarcodes chromsizes binsize outdir");
			System.exit(0);
		}
		String bamfile = args[0];
		String cellbarcodes = args[1];
		String chromsizes = args[2];
		int binsize = Integer.parseInt(args[3]);
		String outdir = args[4];
		
		BinCountGenerator pc = new BinCountGenerator();
		try {
			 pc.generateBinCounts(bamfile, cellbarcodes, chromsizes, binsize, outdir);
		} catch (IOException e) {
			e.printStackTrace();
			try {
				Util u = new Util();
				u.writeErrorFile(e, outdir);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			e.printStackTrace();
		}
	}
	
	public void generateBinCounts(String bamfile, String cellbarcodes, String chromsizes, int binsize, String outdir) throws IOException{

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
		TreeMap<String, Integer> chromsizesmap = u.readChromSizes(chromsizes);
		TreeMap<String, String> cellbarcodemap = u.readCellBarcodes(cellbarcodes);
		
		String curchr = null;
		TreeMap<String, TreeMap<Integer, Integer>> cellindex = null;
		TreeMap<Integer, Integer> totalcellcounts = null;
		
		TreeMap<String, Integer> cellidint = new TreeMap<String, Integer>();
		TreeMap<Integer, TreeSet<Integer>> uniquecellcounts = new TreeMap<Integer, TreeSet<Integer>>();
		
		System.out.println("Reading BAM file.");
		int totalcount = 0;
		int totalreads = 0;
		int barcodeindex = 0;
		while(it.hasNext()){
			if(totalreads % 10000000 == 0) {
				System.out.println(Integer.toString(totalreads));
			}
			
			SAMRecord next = it.next();
			totalreads++;
			boolean isnegative = next.getReadNegativeStrandFlag();
			int insertsize = next.getInferredInsertSize();

			if(!isnegative && u.isValidRead(next) && insertsize > 0) {
				
				String chr = next.getReferenceName();
				if(chromsizesmap.containsKey(chr)) {
					
					//Write previous chromosome file if exploring a new chromosome
					if(!chr.equals(curchr)) {
						if(curchr != null) {
							writeBinCounts(curchr, cellindex, binsize, outdir);
							writeTotalCounts(curchr, totalcellcounts, uniquecellcounts, binsize, outdir);
						}
						cellindex = new TreeMap<String, TreeMap<Integer, Integer>>();
						totalcellcounts = new TreeMap<Integer, Integer>();
						uniquecellcounts = new TreeMap<Integer, TreeSet<Integer>>();
						curchr = chr;
					}
					
					String barcode = u.getBarcode(next);

					if(barcode == null) {
						continue;
					}
					
					String cellid = cellbarcodemap.get(barcode);
					
					if(cellid == null) {
						continue;
					}
					
					int start = next.getAlignmentStart();
					int end = start+insertsize-1;
					
					int si = start/binsize;
					int ei = end/binsize;
					
					if(end-start > 900) {
						continue;
					}
					
					if(!cellidint.containsKey(cellid)) {
						cellidint.put(cellid, barcodeindex++);
					}
					
					if(!cellindex.containsKey(cellid)) {
						cellindex.put(cellid, new TreeMap<Integer, Integer>());
					}
					int curcbi = cellidint.get(cellid);


					TreeMap<Integer, Integer> sparsecount = cellindex.get(cellid);
					for(int i = si; i <= ei; i++) {
						if(sparsecount.containsKey(i)) {
							sparsecount.put(i, sparsecount.get(i)+1);
						}
						else {
							sparsecount.put(i, 1);
						}
						
						if(totalcellcounts.containsKey(i)) {
							totalcellcounts.put(i,  totalcellcounts.get(i)+1);
						}
						else {
							totalcellcounts.put(i, 1);
						}
						
						if(!uniquecellcounts.containsKey(i)) {
							uniquecellcounts.put(i, new TreeSet<Integer>());
						}
						uniquecellcounts.get(i).add(curcbi);
						
					}
					totalcount++;
				}
				
			}
			
		}
		
		System.out.println("Reading BAMs complete.");
		
		StringBuilder output = new StringBuilder();
		output.append("Total Reads: ");
		output.append(Integer.toString(totalreads));
		output.append("\n");

		output.append("Reads Counted: ");
		output.append(Integer.toString(totalcount));
		output.append("\n");

		output.append("Cells Counted: ");
		output.append(cellidint.size());
		output.append("\n");

		u.writeCompletionFile(outdir, "CompletionSummary", output.toString());
		
		it.close();
		reader.close();

		System.out.println("Complete.");
	}
	
	/*private void writeBarcodes(TreeMap<String, Integer> barcodes, String outdir) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outdir+"/BinCounts_barcodes.txt"));
		
		for(Iterator<String> it = barcodes.keySet().iterator(); it.hasNext();) {
			bw.write(it.next()+"\n");
		}
		bw.flush();
		bw.close();
	}*/
	
	private void writeBinCounts(String chr, TreeMap<String, TreeMap<Integer, Integer>> cellcounts, int binsize, String outdir) throws IOException {
		
		BufferedWriter curbw = new BufferedWriter(new FileWriter(outdir+"/BinCounts_"+chr+"_"+binsize+".txt"));
		
		curbw.write("Chromosome");
		curbw.write("\t");
		curbw.write("Bin Index ("+Integer.toString(binsize)+")");
		curbw.write("\t");
		curbw.write("# Inserts Overlapping Bin");
		curbw.write("\t");
		curbw.write("Cell Id");
		curbw.write("\n");
		
		for(Iterator<Entry<String, TreeMap<Integer, Integer>>> it = cellcounts.entrySet().iterator(); it.hasNext();) {
			Entry<String, TreeMap<Integer, Integer>> next = it.next();
			String barcode = next.getKey();
			for(Iterator<Entry<Integer, Integer>> countit = next.getValue().entrySet().iterator(); countit.hasNext();) {
				Entry<Integer, Integer> curcounts = countit.next();
				curbw.write(chr);
				curbw.write("\t");
				curbw.write(Integer.toString(curcounts.getKey()));
				curbw.write("\t");
				curbw.write(Integer.toString(curcounts.getValue()));
				curbw.write("\t");
				curbw.write(barcode);
				curbw.write("\n");
			}
		}
		
		
		curbw.flush();
		curbw.close();
	}
	
	private void writeTotalCounts(String chr, TreeMap<Integer, Integer> totalcounts, TreeMap<Integer, TreeSet<Integer>> uniquecounts, int binsize, String outdir) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outdir+"/TotalBinCounts_"+chr+"_"+binsize+".txt"));
		
		bw.write("Chromosome");
		bw.write("\t");
		bw.write("Bin Index ("+Integer.toString(binsize)+")");
		bw.write("\t");
		bw.write("# Inserts Overlapping Bin");
		bw.write("\t");
		bw.write("# of Unique Cells");
		bw.write("\n");
		
		for(Iterator<Entry<Integer, Integer>> countit = totalcounts.entrySet().iterator(); countit.hasNext();) {
			Entry<Integer, Integer> curcounts = countit.next();
			bw.write(chr);
			bw.write("\t");
			bw.write(Integer.toString(curcounts.getKey()));
			bw.write("\t");
			bw.write(Integer.toString(curcounts.getValue()));
			bw.write("\t");
			bw.write(Integer.toString(uniquecounts.get(curcounts.getKey()).size()));
			bw.write("\n");
		}
		
		bw.flush();
		bw.close();
	}
	
	

	
	
}
