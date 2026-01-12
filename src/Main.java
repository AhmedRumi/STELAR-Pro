import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import preprocessing.GeneTrees;
import utils.Config;
import core.InferenceDP;
import core.SpeciesTreeScorer;
import tree.RangeBipartition;
import tree.MixedBipartition;
import tree.Tree;

/**
 * Main entry point for phylogeny project with GeneTrees processing.
 * 
 * This implementation uses the actual GeneTrees class to parse and process
 * gene trees from Newick format input files.
 */
public class Main {

    /**
     * Main method that handles command line arguments and orchestrates the analysis.
     */
    public static void main(String[] args) throws IOException {

        String inputFilePath = null;
        String outputFilePath = null;
        String computationMode = null;
        String expansionMethod = null;
        String distanceMethod = null;
        boolean verboseExpansion = false;
        boolean disableExpansion = false;
        String branchSupport = null;
        double lambda = 0.5;
        boolean useMixedBipartitions = true;  // Cross-tree recombination flag
        String speciesTreePath = null;  // For score-only mode

        // Parse command line arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i") && i + 1 < args.length) {
                inputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            } else if (args[i].equals("-o") && i + 1 < args.length) {
                outputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            } else if ((args[i].equals("-c") || args[i].equals("--score")) && i + 1 < args.length) {
                speciesTreePath = args[i + 1];
                i++; // Skip next argument as it's the species tree path
            } else if (args[i].equals("-m") && i + 1 < args.length) {
                computationMode = args[i + 1];
                i++; // Skip next argument as it's the mode
            } else if (args[i].equals("-e") && i + 1 < args.length) {
                expansionMethod = args[i + 1];
                i++; // Skip next argument as it's the expansion method
            } else if (args[i].equals("-d") && i + 1 < args.length) {
                distanceMethod = args[i + 1];
                i++; // Skip next argument as it's the distance method
            } else if (args[i].equals("-v")) {
                verboseExpansion = true;
            } else if (args[i].equals("--no-expansion")) {
                disableExpansion = true;
            } else if (args[i].equals("-s") && i + 1 < args.length) {
                branchSupport = args[i + 1];
                i++; // Skip next argument as it's the support type
            } else if (args[i].equals("--lambda") && i + 1 < args.length) {
                try {
                    lambda = Double.parseDouble(args[i + 1]);
                    i++; // Skip next argument as it's the lambda value
                } catch (NumberFormatException e) {
                    System.err.println("Error: Invalid lambda value '" + args[i + 1] + "'");
                    System.exit(-1);
                }
            } else if (args[i].equals("--use-mixed") || args[i].equals("--extend-candidates")) {
                useMixedBipartitions = true;
            } else if (args[i].equals("--no-mixed")) {
                useMixedBipartitions = false;
            }
        }

        // Validate required arguments
        // Score-only mode: -i gene_trees.tre -c species_tree.tre (no -o needed)
        // Inference mode: -i gene_trees.tre -o output.tre
        if (inputFilePath == null || (outputFilePath == null && speciesTreePath == null)) {
            System.out.println("Usage:");
            System.out.println("  Inference mode: java Main -i <gene_trees> -o <output_file> [options]");
            System.out.println("  Score mode:     java Main -i <gene_trees> -c <species_tree> [options]");
            System.out.println("");
            System.out.println("Options:");
            System.out.println("  -c <tree>     Calculate triplet score between gene trees and given species tree");
            System.out.println("  --score <tree> Same as -c");
            System.out.println("  -m <mode>     Computation mode: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
            System.out.println("  -e <method>   Expansion method: NONE, DISTANCE_ONLY, CONSENSUS_ONLY, DISTANCE_CONSENSUS, FULL");
            System.out.println("  -d <method>   Distance method: UPGMA, NEIGHBOR_JOINING, BOTH");
            System.out.println("  -s <support>  Branch support: NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL");
            System.out.println("  --lambda <val> Lambda parameter for branch support (default: 0.5)");
            System.out.println("  -v            Verbose expansion output");
            System.out.println("  --no-expansion Disable bipartition expansion");
            System.out.println("  --use-mixed   Enable cross-tree recombination (default: ON)");
            System.out.println("  --no-mixed    Disable cross-tree recombination");
            System.exit(-1);
        }

        // Validate input file exists
        File inputFile = new File(inputFilePath);
        if (!inputFile.exists()) {
            System.err.println("Error: Input file '" + inputFilePath + "' does not exist.");
            System.exit(-1);
        }

        // Set computation mode if specified
        if (computationMode != null) {
            try {
                Config.COMPUTATION_MODE = Config.ComputationMode.valueOf(computationMode);
            } catch (IllegalArgumentException e) {
                System.err.println("Error: Invalid computation mode '" + computationMode + "'");
                System.err.println("Valid modes: CPU_SINGLE, CPU_PARALLEL, GPU_PARALLEL");
                System.exit(-1);
            }
        }
        
        // Configure bipartition expansion
        if (disableExpansion) {
            utils.BipartitionExpansionConfig.EXPANSION_METHOD = utils.BipartitionExpansionConfig.ExpansionMethod.NONE;
        } else if (expansionMethod != null) {
            try {
                utils.BipartitionExpansionConfig.EXPANSION_METHOD = 
                    utils.BipartitionExpansionConfig.ExpansionMethod.valueOf(expansionMethod);
            } catch (IllegalArgumentException e) {
                System.err.println("Error: Invalid expansion method '" + expansionMethod + "'");
                System.err.println("Valid methods: NONE, DISTANCE_ONLY, CONSENSUS_ONLY, DISTANCE_CONSENSUS, FULL");
                System.exit(-1);
            }
        }
        
        // Set distance method if specified
        if (distanceMethod != null) {
            try {
                utils.BipartitionExpansionConfig.DISTANCE_METHOD = 
                    utils.BipartitionExpansionConfig.DistanceMethod.valueOf(distanceMethod);
            } catch (IllegalArgumentException e) {
                System.err.println("Error: Invalid distance method '" + distanceMethod + "'");
                System.err.println("Valid methods: UPGMA, NEIGHBOR_JOINING, BOTH");
                System.exit(-1);
            }
        }
        
        // Set verbose expansion if specified
        if (verboseExpansion) {
            utils.BipartitionExpansionConfig.VERBOSE_EXPANSION = true;
            System.out.println("Verbose expansion output enabled.");
        }

        // Determine mode
        boolean scoreMode = (speciesTreePath != null);
        
        System.out.println("Input file: " + inputFilePath);
        if (scoreMode) {
            System.out.println("Mode: SCORE (calculate triplet score for given species tree)");
            System.out.println("Species tree: " + speciesTreePath);
        } else {
            System.out.println("Mode: INFERENCE (find optimal species tree)");
            System.out.println("Output file: " + outputFilePath);
        }
        System.out.println("Computation mode: " + Config.COMPUTATION_MODE);
        if (!scoreMode) {
            System.out.println("Expansion method: " + utils.BipartitionExpansionConfig.EXPANSION_METHOD);
            if (utils.BipartitionExpansionConfig.isDistanceExpansionEnabled()) {
                System.out.println("Distance method: " + utils.BipartitionExpansionConfig.DISTANCE_METHOD);
            }
            System.out.println("Cross-tree recombination (--use-mixed): " + (useMixedBipartitions ? "ENABLED" : "disabled"));
            if (branchSupport != null) {
                System.out.println("Branch support: " + branchSupport);
                System.out.println("Lambda parameter: " + lambda);
            }
        }

        long startTime = System.nanoTime();

        // Process gene trees
        GeneTrees geneTrees = new GeneTrees(inputFilePath);
        geneTrees.readTaxaNames(); // Ensure taxaMap is initialized
        geneTrees.readGeneTrees(null);

        // ================================================================
        // SCORE MODE: Calculate triplet score for given species tree
        // ================================================================
        if (scoreMode) {
            // Validate species tree file exists
            File speciesFile = new File(speciesTreePath);
            if (!speciesFile.exists()) {
                System.err.println("Error: Species tree file '" + speciesTreePath + "' does not exist.");
                System.exit(-1);
            }
            
            // Read species tree newick
            String speciesNewick = null;
            try (BufferedReader reader = new BufferedReader(new FileReader(speciesTreePath))) {
                speciesNewick = reader.readLine();
                if (speciesNewick != null) {
                    speciesNewick = speciesNewick.trim();
                }
            }
            
            if (speciesNewick == null || speciesNewick.isEmpty()) {
                System.err.println("Error: Species tree file is empty.");
                System.exit(-1);
            }
            
            System.out.println("\nParsing species tree...");
            Tree speciesTree = new Tree(speciesNewick, geneTrees.taxaMap);
            System.out.println("Species tree parsed: " + speciesTree.leavesCount + " leaves");
            
            // Calculate triplet score
            SpeciesTreeScorer scorer = new SpeciesTreeScorer(geneTrees);
            double score = scorer.calculateScore(speciesTree);
            
            long endTime = System.nanoTime();
            double duration = (endTime - startTime) / 1_000_000_000.0;
            
            System.out.println("\n========================================");
            System.out.println("TRIPLET_SCORE: " + score);
            System.out.println("========================================");
            System.out.println("Time taken: " + duration + " seconds");
            System.out.println("Score calculation completed successfully!");
            
            return;
        }

        // ================================================================
        // INFERENCE MODE: Find optimal species tree via DP
        // ================================================================
        
        // Generate candidate bipartitions with cross-tree recombination extension
        System.out.println("Generating candidate bipartitions...");
        
        // Generate mixed bipartitions via cross-tree recombination
        // These will only be used in DP if --use-mixed flag is set
        List<RangeBipartition> candidates = geneTrees.generateExtendedCandidateBipartitions(useMixedBipartitions);
        System.out.println("Total candidate bipartitions (gene tree): " + candidates.size());
        
        // Report mixed bipartitions generated by cross-tree recombination
        List<MixedBipartition> mixedBips = geneTrees.getMixedBipartitions();
        if (mixedBips != null && !mixedBips.isEmpty()) {
            System.out.println("Mixed bipartitions from cross-tree recombination: " + mixedBips.size());
            long crossTree = mixedBips.stream().filter(MixedBipartition::isCrossTree).count();
            System.out.println("  - Cross-tree (truly new): " + crossTree);
            System.out.println("  - Same-tree: " + (mixedBips.size() - crossTree));
        }

        // Run inference
        InferenceDP inference = new InferenceDP(geneTrees, candidates);
        
        // Enable mixed bipartitions in DP if flag is set
        if (useMixedBipartitions && mixedBips != null && !mixedBips.isEmpty()) {
            System.out.println("\nEnabling mixed bipartitions in DP inference...");
            inference.enableMixedBipartitions(mixedBips);
        }
        
        double score = inference.solve();
        Tree resultTree = inference.reconstructTree();

        // Calculate branch support if requested
        if (branchSupport != null && !branchSupport.equals("NONE")) {
            System.out.println("\nCalculating branch support...");
            
            core.BranchSupportCalculator.BranchAnnotationType annotationType;
            try {
                switch (branchSupport.toUpperCase()) {
                    case "POSTERIOR":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.POSTERIOR_ONLY;
                        break;
                    case "DETAILED":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.DETAILED;
                        break;
                    case "LENGTH":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.BRANCH_LENGTH_ONLY;
                        break;
                    case "BOTH":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.POSTERIOR_AND_LENGTH;
                        break;
                    case "PVALUE":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.PVALUE_ONLY;
                        break;
                    case "ALL":
                        annotationType = core.BranchSupportCalculator.BranchAnnotationType.ALL;
                        break;
                    default:
                        System.err.println("Error: Invalid branch support type '" + branchSupport + "'");
                        System.err.println("Valid types: NONE, POSTERIOR, DETAILED, LENGTH, BOTH, PVALUE, ALL");
                        System.exit(-1);
                        return;
                }
            } catch (Exception e) {
                System.err.println("Error parsing branch support type: " + e.getMessage());
                System.exit(-1);
                return;
            }
            
            core.BranchSupportCalculator supportCalculator = 
                new core.BranchSupportCalculator(geneTrees, resultTree, lambda, annotationType);
            
            // Validate quartet frequencies for debugging (optional)
            if (verboseExpansion) {
                supportCalculator.validateQuartetFrequencies();
            }
            
            // Annotate branches
            supportCalculator.annotateBranches();
            
            // Print statistics
            core.BranchSupportCalculator.BranchSupportStatistics stats = 
                supportCalculator.calculateStatistics();
            System.out.println("\n" + stats.toString());
        }

        // Write output
        try (FileWriter writer = new FileWriter(outputFilePath)) {
            writer.write(resultTree.getNewickFormat());
        }

        long endTime = System.nanoTime();
        double duration = (endTime - startTime) / 1_000_000_000.0; // Convert to seconds

        System.out.println("\n========================================");
        System.out.println("OPTIMAL_TRIPLET_SCORE: " + score);
        System.out.println("========================================");
        System.out.println("Time taken: " + duration + " seconds");
        System.out.println("Program completed successfully!");
        System.out.println("Output written to: " + outputFilePath);
    }

    /**
     * Processes gene trees using the GeneTrees class and returns analysis results.
     * 
     * @param inputFilePath Path to the input file containing gene trees in Newick format
     * @return Formatted string with analysis results
     * @throws FileNotFoundException if the input file cannot be read
     */
    private static String processGeneTrees(String inputFilePath) throws FileNotFoundException {
        System.out.println("Initializing GeneTrees...");
        
        // Create GeneTrees object and read taxa names
        GeneTrees geneTrees = new GeneTrees(inputFilePath);
        var taxaMap = geneTrees.readTaxaNames();

        
        System.out.println("Reading and parsing gene trees...");
        
        // Read and process all gene trees
        geneTrees.readGeneTrees(null); // No distance matrix needed for basic analysis
        
        // Debug output
        // debugOutput(geneTrees);

        System.out.println(geneTrees.geneTrees.get(0).isRooted);
        
        // Test InferenceDP algorithm
        System.out.println("Testing InferenceDP algorithm...");
        List<RangeBipartition> candidates = new ArrayList<>(geneTrees.rangeBipartitions.keySet());
        
        if (!candidates.isEmpty()) {
            InferenceDP dp = new InferenceDP(geneTrees, candidates);
            double maxScore = dp.solve();
            
            System.out.println("DP Algorithm completed with maximum score: " + maxScore);
            
            Tree reconstructedTree = dp.reconstructTree();
            if (reconstructedTree != null && reconstructedTree.root != null) {
                System.out.println("Tree reconstruction successful");
                return reconstructedTree.getNewickFormat();
            }
        }
        
        return "";
    }
    
    /**
     * Debug function that prints detailed analysis information to console
     */
    private static void debugOutput(GeneTrees geneTrees) {
        // Gather analysis information
        int geneTreeCount = geneTrees.geneTrees.size();
        int taxaCount = geneTrees.realTaxaCount;
        // int uniquePartitions = geneTrees.triPartitions.size();
        int uniqueRangeBipartitions = geneTrees.rangeBipartitions.size();
        
        System.out.println("Processing complete:");
        System.out.println("  - Gene trees processed: " + geneTreeCount);
        System.out.println("  - Taxa found: " + taxaCount);
        // System.out.println("  - Unique tripartitions: " + uniquePartitions);
        System.out.println("  - Unique RangeBipartitions: " + uniqueRangeBipartitions);
        
        // Print taxa names
        System.out.print("Taxa names: ");
        for (int i = 0; i < geneTrees.taxonIdToLabel.length; i++) {
            if (i > 0) System.out.print(", ");
            System.out.print(geneTrees.taxonIdToLabel[i]);
        }
        System.out.println();
        
        // Print RangeBipartitions with counts
        System.out.println("RangeBipartitions:");
        for (var entry : geneTrees.rangeBipartitions.entrySet()) {
            System.out.println("  " + entry.getKey().toString() + " : " + entry.getValue());
        }
    }

    /**
     * Writes analysis results to the specified output file.
     * 
     * @param outputFilePath Path to the output file
     * @param content Content to write to the file
     * @throws IOException if there's an error writing to the file
     */
    private static void writeResults(String outputFilePath, String content) throws IOException {
        FileWriter writer = new FileWriter(outputFilePath);
        writer.write(content);
        writer.close();
    }
}
