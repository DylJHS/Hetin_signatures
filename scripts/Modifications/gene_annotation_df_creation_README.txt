
The gene_annotation_df_creation.RMD script generates a gene annotation table that includes chromosomal arm labels (p or q) for all annotated protein-coding genes in the human genome (GRCh38). This annotation is essential for analyses where the genomic position of genes relative to chromosomal instability (CIN) features is important—such as identifying and filtering cis-acting genes (i.e., those located on the same chromosome arm as the target CIN feature) during model interpretation.

The output is used to support interpretability filtering in feature importance analyses of a machine learning model trained to predict CIN features (e.g., arm-level aneuploidy or pericentromeric CNVs) from gene expression data.

Purpose
In the context of this project, certain CIN features are defined per chromosome arm (e.g., "5p" or "17q"). When interpreting feature importance from machine learning predictions, it is often desirable to:

- Exclude genes located on the same arm as the target feature to focus on trans-effects.
- Optionally retain or highlight cis-acting genes separately, particularly if they represent mechanistically relevant candidates.

To do this reliably, each gene must be annotated not just by chromosome, but by chromosomal arm.

Input Files
Gene Annotation GTF File
- Filename: Homo_sapiens.GRCh38.114.chr.gtf.gz
- Source: Ensembl
- Description: Contains gene-level annotations including gene names, genomic start and end positions, and chromosome identifiers, for the GRCh38 genome build.
- Usage: Extracts gene name, chromosome, and start position to determine chromosomal location and arm assignment.

Centromere Position File
- Filename: centromeres.txt.gz
- Source: UCSC Genome Browser – hg38 database
- Description: Provides the centromere midpoints for each chromosome, which are used as boundaries to separate the p and q arms.
- Usage: Used to define the transition point between the p and q arms for each chromosome.

Processing Steps
Filter Gene Entries:

The GTF file is filtered to retain only rows where the feature type is "gene".

Extract Metadata:

- Gene names and Ensembl IDs are parsed from the attributes column.

- Chromosome names and gene start positions are retained.

Standardise Chromosomes:

- Only genes located on standard chromosomes (1–22, X, Y) are retained.

- Any genes on scaffolds or alternative contigs are excluded.

Merge with Centromere Data:

- The UCSC centromere midpoint file is used to assign each gene to either the p arm (if the start coordinate is less than or equal to the centromere midpoint) or the q arm (if the start coordinate is greater than the midpoint).

- This logic is applied across all chromosomes.

Deduplicate and Finalise:

- For genes with multiple entries, the script retains a unique row per gene name.

The final DataFrame includes columns for gene name, chromosome, start and end coordinates, and arm label.

Outputs
The script produces two annotation files:

gene_annotation_df:
A full gene annotation table with the following fields:
- gene_name
- chr (e.g. "chr1")
- start
- end
- arm (either "p" or "q")

SOI_gene_annotation_df:
- A filtered version of the full table, containing only the genes in the Set of Interest (SOI)—a predefined list of heterochromatin-associated genes used as input features in the CIN prediction model.

Use Case
The outputs from this script are used for the model feature importance interpretation, specifically to:

- Filter out cis-located genes from feature importance results.
- Create separate heatmaps or plots showing either all genes or only trans-acting genes per CIN feature.
- Annotate feature importance tables with genomic context for biological interpretation.
