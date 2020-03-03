package org.jax.snatacclusteringtools.splitbam;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import org.jax.snatacclusteringtools.util.Util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BAMSplitter {
	
	public static void main(String[] args) {
		if(args.length != 4) {
			System.out.println("Usage: splitbams bamfile cellbarcodes cellclusters outdir");
			System.exit(0);
		}
		String bamfile = args[0];
		String cellbarcodes = args[1];
		String cellcluster = args[2];
		String outdir = args[3];
		
		BAMSplitter bcs = new BAMSplitter();
		try {
			bcs.splitByCluster(bamfile, cellbarcodes, cellcluster, outdir);
		} catch (IOException e) {
			try {
				Util u = new Util();
				u.writeErrorFile(e, outdir);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			e.printStackTrace();
		}
	}
	
	public void splitByCluster(String bamfile, String cellbarcodes, String cellclusterfile, String outdir) throws IOException{

		SamReaderFactory factory = SamReaderFactory.makeDefault()
	              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	              .validationStringency(ValidationStringency.SILENT);
		
		final SamReader reader = factory.open(new File(bamfile));

		SAMFileHeader header = reader.getFileHeader();

		
		if(!(header.getSortOrder().getComparatorInstance() instanceof SAMRecordCoordinateComparator)) {
			System.out.println("The input BAM file must be coordinate sorted.");
			reader.close();
			System.exit(0);
		}
		
		
		SAMRecordIterator it = reader.iterator();
		Util u = new Util();
		TreeMap<String, String> cellbarcodemap = u.readCellBarcodes(cellbarcodes);
		ClusterInfo ci = readCellClusters(cellclusterfile);
		TreeMap<String, Integer> cellclusters = ci.cellclustermap;
		TreeSet<Integer> clusters = ci.clusters;
		
		
		TreeMap<Integer, SAMFileWriter> writers = new TreeMap<Integer, SAMFileWriter>();
		SAMFileWriterFactory wfactory = new SAMFileWriterFactory();
		for(Iterator<Integer> clusterit = clusters.iterator(); clusterit.hasNext();) {
			int curcluster = clusterit.next();
			File f = new File(outdir+"/clusterbam_"+curcluster+".bam");
			if(f.exists()) {
				System.out.println("Error:"+f.getAbsolutePath()+" exists.  Exiting.");
				System.exit(1);
			}
			SAMFileWriter curwriter = wfactory.makeBAMWriter(header, true, new File(outdir+"/clusterbam_"+curcluster+".bam"));
			writers.put(curcluster, curwriter);
		}

		System.out.println("Splitting BAM files into clusters.");

		while(it.hasNext()){
			SAMRecord next = it.next();
			
			String barcode = u.getBarcode(next);

			if(barcode == null) {
				continue;
			}
				
			String cellid = cellbarcodemap.get(barcode);
				
			if(cellid == null) {
				continue;
			}
			
			Integer cluster = cellclusters.get(cellid);
			
			if(cluster == null) {
				System.out.println("Missing cluster for cell "+cellid);
			}
			
			SAMFileWriter curwriter = writers.get(cluster);
			curwriter.addAlignment(next);
			
		}

		it.close();
		reader.close();
		
		for(Iterator<SAMFileWriter> wit = writers.values().iterator(); wit.hasNext();) {
			SAMFileWriter next = wit.next();
			next.close();
		}
		
		writeCompletionFile(outdir, clusters, cellclusters);

		System.out.println("Complete.");
	}
		
	
	private void writeCompletionFile(String outdir, TreeSet<Integer> clusters, TreeMap<String, Integer> cellclusters) throws IOException {
		Util u = new Util();
		
		
		int[] clustercounts = new int[clusters.pollLast()+1];
		
		for(Iterator<Integer> it = cellclusters.values().iterator(); it.hasNext();){
			clustercounts[it.next()]++;
		}
		
		StringBuilder sb = new StringBuilder();
		sb.append("Cluster ID\tNumber of Cells\n");

		for(int i = 0; i < clustercounts.length; i++) {
			sb.append("Cluster "+Integer.toString(i)+"\t"+Integer.toString(clustercounts[i])+"\n");
		}
		
		u.writeCompletionFile(outdir, "CompletionSummary", sb.toString());

	}
	
	private ClusterInfo readCellClusters(String cellclusterfile) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(cellclusterfile));
		TreeMap<String, Integer> clustermap = new TreeMap<String, Integer>();
		TreeSet<Integer> clusters = new TreeSet<Integer>();
		while(br.ready()) {
			String line = br.readLine();
			String[] split = line.split("\t|,");
			String cellid = split[0];
			int clusterid = Integer.parseInt(split[1]);
			clustermap.put(cellid, clusterid);
			clusters.add(clusterid);
		}
		br.close();
		ClusterInfo rv = new ClusterInfo();
		rv.cellclustermap = clustermap;
		rv.clusters = clusters;
		return rv;
	}
	
	
	private class ClusterInfo {
		
		TreeMap<String, Integer> cellclustermap;
		TreeSet<Integer> clusters;
		
	}
	
	
}
