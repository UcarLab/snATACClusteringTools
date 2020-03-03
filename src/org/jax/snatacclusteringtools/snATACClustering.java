package org.jax.snatacclusteringtools;

import org.jax.snatacclusteringtools.bincount.BinCountGenerator;
import org.jax.snatacclusteringtools.bincount.BinCountMerger;
import org.jax.snatacclusteringtools.peakcount.PeakCountGenerator;
import org.jax.snatacclusteringtools.peakcount.PeakCountMerger;
import org.jax.snatacclusteringtools.peakprocessing.SummitPeakMerger;
import org.jax.snatacclusteringtools.sparsematrix.SparseMatrixGenerator;
import org.jax.snatacclusteringtools.splitbam.BAMSplitter;

public class snATACClustering {
	
	public static void main(String[] args) {
		if (args.length < 1) {
			System.out.println("Usage: snATACClustering.jar command arguments");
			System.out.println("Commands:");
			System.out.println("generatebincounts");
			System.out.println("mergebincounts");
			System.out.println("generatesparsematrix");
			System.out.println("splitbams");
			System.out.println("mergepeaksummits");
			System.out.println("generatepeakcounts");
			System.out.println("mergepeakcounts");
			System.exit(0);
		}
		
		new snATACClustering().runCommand(args);
	}

	private void runCommand(String[] args) {
		String command = args[0];
		String[] newargs = getCommandArguments(args);
		
		switch (command) {
		case "generatebincounts":
			BinCountGenerator.main(newargs);
			break;
		case "mergebincounts":
			BinCountMerger.main(newargs);
			break;
		case "generatesparsematrix":
			SparseMatrixGenerator.main(newargs);
			break;
		case "splitbams":
			BAMSplitter.main(newargs);
			break;
		case "processpeaks":
			SummitPeakMerger.main(newargs);
			break;
		case "generatepeakcounts":
			PeakCountGenerator.main(newargs);
			break;
		case "mergepeakcounts":
			PeakCountMerger.main(newargs);
			break;
		default:
			System.out.println("Command not found.");
			break;
		}
		
	}

	private String[] getCommandArguments(String[] origargs) {
		String[] rv = new String[origargs.length - 1];
		for (int i = 0; i < rv.length; i++) {
			rv[i] = origargs[i + 1].trim();
		}
		return rv;
	}
}
