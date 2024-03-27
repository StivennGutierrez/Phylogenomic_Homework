# 🧬🔍 ABS 18S-rDNA _Trypanosoma cruzi_ isolates 📟
Understanding the epidemiology and impact of Tripanosomatid species is crucial for disease control and prevention. Current diagnostic methods face challenges due to the high genetic diversity of these species. Conventional PCR and high-resolution melting genotyping have limitations in resolving diversity. Amplicon-based sequencing (ABS) emerges as a promising solution, utilizing universal primers for simultaneous DNA amplification of multiple species. ABS, using Illumina sequencing, has shown advantages in diversity detection. However, the potential of long-read Nanopore sequencing remains unexplored. The 18S rDNA gene proves vital for ABS, enabling accurate identification and genotyping of Tripanosomatids.

#### :hospital:The data utilized in this pipeline consists of DNA extracted from samples of 23 patients who tested positive for Chagas disease through qPCR, referred from the Fundación Cardio Infantil (FCI) hospital.

The resulting library was loaded onto an R9.4 flow cell (FLO-MIN106) and sequenced on the MinIonMk1C instrument with MinKNOW software V22.10.7. 

## Basecalling 📲
Basecalling was conducted using Guppy v6.3.9, and reads below a minimum quality score of seven were excluded from further analysis.

## Nanostat for _T. cruzi_ Sequencing Quality Assessment 📈
**Nanostat** is a versatile tool for analyzing FASTQ and SAM/BAM files. It offers a comprehensive set of metrics to assess the quality of sequencing data, including base calling accuracy, read length distribution, and GC content. When applied to FASTQ archives of ***T. cruzi*** sequencing obtained through Oxford Nanopore, Nanostat can provide valuable insights into the quality of the sequencing data and aid in identifying potential issues that may affect downstream analyses.

1) **Installation:**

To utilize Nanostat for quality assessment, start by installing it with a simple command. Execute the following in your terminal:

`curl -sSL https://raw.githubusercontent.com/sanger-pathogens/nanostat/master/install.sh | sh`

2) **Running Nanostat on _T. cruzi_ FASTQ Files:**

Once Nanostat is installed, you can assess the quality of your _T. cruzi_ sequencing data with the following command:

`nanostat input.fastq`

3) **Interpreting Nanostat Output:**

🚨Nanostat generates a comprehensive report summarizing quality metrics specific to _T. cruzi_ sequencing. Key metrics include:

**Read N50:** This signifies the median read length, crucial for evaluating sequencing quality. Higher values suggest superior sequencing quality.

**Base Calling Accuracy:** Indicates the proportion of accurately called bases during sequencing. High accuracy is vital for reliable downstream analyses.

**GC Content:** Measures the proportion of guanine (G) and cytosine (C) bases. _T. cruzi_ typically has a GC content of around 45%.

>Analyzing these metrics allows an assessment of overall sequencing quality and identification of potential areas for improvement. A low read N50 may indicate sample preparation or sequencing condition issues, while low base calling accuracy could point to problems with the sequencing device or data analysis pipeline.

## Filtlong for Improved _T. cruzi_ Sequencing Data Quality ✂️
Filtlong is a specialized long-read filtering tool crafted for Oxford Nanopore sequencing data, with a focus on enhancing the quality of _T. cruzi_ sequencing. Tailored to address challenges posed by _T. cruzi_ genomes, known for repetitive regions and homopolymer stretches, Filtlong employs a sophisticated algorithm to eliminate low-quality and chimeric reads from FASTQ archives, consequently refining overall data quality and file size.

1) **Installation:**

Initiate the installation of Filtlong with the following command:

`pip install filtlong`

2) **Download _T. cruzi_ Reference Genome:**

To facilitate the identification of chimeric reads, Filtlong requires a reference genome. Download the _T. cruzi_ reference genome from reputable sources like NCBI or [TriTrypDB](https://tritrypdb.org/tritrypdb/app/record/genomic-sequence/TcX10_chr1).

3) **Filtering T. cruzi FASTQ Files:**
Once Filtlong is installed and the reference genome is obtained, filter T. cruzi FASTQ files with the following command:

`filtlong -i input.fastq -o filtered.fastq -r reference.fasta`

🚨Advantages of Filtlong for _T. cruzi_ Sequencing:

**Improved Quality:** Filtlong excels in the removal of low-quality reads, mitigating the risk of errors in subsequent analyses.

**Reduced Chimeric Reads:** By effectively eliminating chimeric reads—formed from fragments of different DNA molecules—Filtlong streamlines assembly and downstream analyses.

**Increased Read Size:** Through the removal of low-quality and chimeric reads, Filtlong enhances the average read length, contributing to improved assembly contiguity and accuracy.

## Centrifuge for Efficient Taxonomic Classification 📚
Centrifuge, a swift taxonomic classifier, employs a compressed k-mer database for efficient analysis of Oxford Nanopore sequencing data. Notably effective for _T. cruzi_ sequencing, it precisely assigns reads to the species level, aiding in accurate identification. The system's innovative indexing scheme, rooted in the Burrows-Wheeler transform and Ferragina-Manzini index, optimizes metagenomic classification for both speed and memory efficiency. 

1) **Installation:**
Begin by installing Centrifuge through the following commands:

`git clone https://github.com/vgteam/centrifuge.git`

`cd centrifuge`

`./build.sh`

2) **Downloading Centrifuge Reference Database:**
  
Centrifuge requires a reference database for taxonomic classification. Download the latest Centrifuge reference database from the official Centrifuge website.

3) **Running Centrifuge on _T. cruzi_ FASTQ Files:**

After installing Centrifuge and obtaining the reference database, classify your T. cruzi FASTQ files with the following command:

- Activate module on cluster:

`module load centrifuge/1.0.3`

- Preparation and file permissions:

`cat *fast.gz > Z4-rojo.fastq.gz`

`gunzip Z4-rojo.fastq.gz`

`chmod 777 Z4-rojo.fastq`

- Additional filtering for reads:

`NanoPlot --fastq Z4-rojo.fastq`  ❗Review the distribution of reads

`filtlong --min_length 800 --keep_percent 90 --min_mean_q 10 Z4-rojo.fastq.gz > Z4-rojo_filter.fastq.gz`  ❗Filter reads by size and quality

`NanoFilt -z -q 10 -l 600 Z4-rojo.fastq.gz > Z4-rojo_filter.fastq.gz`  ❗Filter reads by size and quality

`NanoPlot --fastq Z4-rojo_filter.fastq`  ❗Check distribution of reads after filtering

1. Database indexing:

`centrifuge-build -p 4 --conversion-table seqid2taxid.txt \`  ❗Specify only the organism name and its corresponding taxonomic number in the GFASTA file

`--taxonomy-tree  /datacnmat01/biologia/gimur/dbtaxonomizr/nodes.dmp --name-table /datacnmat01/biologia/gimur/dbtaxonomizr/names.dmp \
18s_base.fa 18s_base`

2. Running 143:

`centrifuge -x /datacnmat01/biologia/gimur/dbtaxonomizr/cruzi/18s_base /datagimur/GIMUR2/5_MinION_data/Corridas_2023/R143-ABS/20230203_1909_MC-110797_FAS85573_5d708e48/fastq_pass/barcode86/Z4-rojo.fastq \`

`--report-file Z4-rojo.txt -S  Z4-rojo.csv -k 1 --min-totallen 100 --qc-filter`

3. Output:

After Centrifuge provides the assignment in a .csv file, open the file in Excel to apply the corresponding filter based on scores in Seq_Length and Hit_Length within the following range:

► Seq-Length: 800-1100

► Hit-Length: =>450

🚨Additionally, categories that were not assigned to _T. cruzi_ (e.g., species) are disregarded.

🚨Finally, assignments of _T. cruzi_ strains other than TcI were validated using Blastn on NCBI and TriTrypDB, in accordance with the respective sequence reference code.

## 🌄 _T. cruzi_'s Chromatic Spectrum: Relative Abundances Illuminated 📊
First and foremost, it is crucial to establish the appropriate database for visualizing relative abundances. For this purpose, four columns are considered (corresponding barcode, strain name, associated DTU, and the count of each sequence presence).

![image](https://github.com/StivennGutierrez/BIO-analysis/assets/128840301/9763c0b3-16a7-4127-8cff-b07b1fb3d0b2)

**ggplot2** in R plays a pivotal role in visualizing _T. cruzi_ strains from Chagas disease patient samples obtained through Oxford Nanopore sequencing. ggplot2's versatility in creating bar charts, line charts, and scatterplots, tailored for _T. cruzi_ sequencing data, stands out. It serves as a potent tool for generating insightful visualizations, portraying the relative abundances of strains. This visual exploration contributes to a nuanced comprehension of the distribution of _T. cruzi_ strains, enhancing our understanding of Chagas disease dynamics.

To do this we run the following commands within the Rstudio environment:

Required libraries:

`library(dplyr)`

`library(ggplot2)`

Load database and organize by factors:

`DTUs_fci$Asignación <- factor(DTUs_fci$Asignación, levels = c("Tcruzi_Brazil_A4_TcI", "Tcruzi_DA_TcI", "Tcruzi_Dm28c_TcI", "Tcruzi_Esmeraldo_cl3_TcII","Tcruzi_G_TcI","Tcruzi_strain_231_TcI","Tcruzi_Tc1","Tcruzi_Y_TcII","TcruziSilvio_X10_cl1_TcI","Trypanosoma rangeli Choachi"))`

Filter data by excluding rows with NA values in assignment:

`DTUs_fci_filtered <- DTUs_fci[!is.na(DTUs_fci$Asignación), ]`

Create chart with filtered data:

`ggplot(DTUs_fci_filtered, aes(x = Barcode, y = as.numeric(Abundancia), fill = Asignación)) +
  geom_bar(position = "fill", stat = "identity", col = NA) +
  scale_fill_viridis_d(option = "inferno") +  # Utilizar la escala de colores viridis
  theme_bw() +
  ylab("Abundancia (%)") +
  xlab("Barcodes positivos") +
  labs(fill = "Cepa-DTU") +
  ggtitle("Pacientes ECHA - FCI") +
  theme(axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 12,face = "italic"))`

Following this command line you will get a graph as follows:

![Abundancia_Barcodes](https://github.com/StivennGutierrez/BIO-analysis/assets/128840301/554dcb8f-0746-4994-b49b-a3e5fb2a49a5)

## Creation of the multifasta 🔠

From the consensus FastQ, extract the verified and positive sequences for _T. cruzi_ [18S_seqs](https://github.com/StivennGutierrez/BIO-analysis/blob/main/all_seqs.fasta), considering their respective identification code.

## Harmonizing Sequences: MAFFT Magic 🎹 

MAFFT is a versatile multiple sequence alignment (MSA) tool that is particularly well-suited for aligning long-read sequencing data, such as that generated by Oxford Nanopore sequencers. For T. cruzi sequencing, MAFFT can effectively align long reads from different T. cruzi strains or closely related organisms, enabling the identification of conserved and variable regions within the T. cruzi genome.

`module load mafft/7.407` !We load the version installed in our terminal, in my case is v7.407!

`mafft --thread 15 R143_multifasta.fasta > R143_output_mafft_align.fasta` !**--thread 15**: This setting indicates the number of threads or processors to be used to perform the alignment. In this case, it is specifying that 15 threads be used to speed up the alignment process. You can adjust this number according to the amount of processing resources available on your system¡.

## SNP Safari: Unveiling _T. cruzi_'s Genomic Mysteries 🎯
We are now embarking on a thorough SNP analysis utilizing two potent tools: SNP-sites and SNP-dists :deciduous_tree:. Tailored for the intricate exploration of _T. cruzi_'s genetic makeup, [SNP-sites](http://sanger-pathogens.github.io/snp-sites/) excels in identifying Single Nucleotide Polymorphisms (SNPs) within aligned sequences specific to this parasite. Complementing this, [SNP-dists](https://github.com/tseemann/snp-dists) takes it a step further by calculating genetic distances based on SNP variations among _T. cruzi_ sequences. This dynamic duo enables researchers to extract and scrutinize SNP information within a specialized phylogenetic framework, providing invaluable insights into the evolutionary relationships and genetic diversity of _T. cruzi_.

First let's detect the SNPs in the debugged multifasta file with the following command: 

`module load snp-sites/2.4.1` !We load the tool in our terminal environment¡

`snp-sites R143_output_mafft_align.fasta > R143_SNPsites.fasta` ❗To obtain a fasta file

`snp-sites -v -o R143_SNPsites.vcf R143_output_mafft_align.fasta` ❗To obtain a vcf file

🚨It will generate a **.fasta** file report for our multifasta alignment with SNPs sites. !You can obtain files in alternative formats including **.aln** (alignment), **.vcf** (which provides the position of each SNP in the reference sequence and its occurrence in other samples), and **Relaxed Phylip format** (which presents all SNP sites in a format suitable for RAxML and other tree building applications)¡.

Then, to calculate the distances for the detected SNPs we will use **snp-dists** with the following command:

`module load snp-dists 0.8.2`

`snp-dists R143_SNPsites.fasta > R143_SNPdists.tsv`

🚨It will generate a **.tsv** file report for our multifasta alignment. !The .tsv format matrix produced by snp-dists contains pairwise genetic distances between sequences, calculated based on SNP variations. Each sequence is represented by a row and column in the matrix, and the value at the intersection of a row and column represents the genetic distance between those two sequences. This matrix serves as a quantitative measure of the genetic dissimilarities among sequences, enabling the comparison and analysis of their evolutionary relationships¡.

>As a tip you could visualize the distance matrix with [Heatmapper](http://www.heatmapper.ca/), this web tool is very useful for visualizing research data, since it allows you to generate heatmaps by clustering methods such as hierarchical clustering and k-means🌌, for example:

![heatmap_143](https://github.com/StivennGutierrez/BIO-analysis/assets/128840301/614e1c2f-3577-438f-a845-d573ec4b5db2)

## Phylogenetic Odyssey: Tracing _T. cruzi_'s Genetic Roots 🎨

Now, we embark on the construction of a phylogenetic tree, a powerful tool shedding light on the genetic diversity, epidemiology, and strain distribution of _T. cruzi_. Employing the widely adopted **maximum likelihood algorithm**, we delve into assessing the likelihood that the proposed model and hypothesized evolutionary history align with the observed sequencing data. This approach allows us to unravel the intricate web of evolutionary relationships within _T. cruzi_, providing valuable insights into its genomic landscape and epidemiological dynamics.

There are several tools that provide the option of generating phylogenetic trees such as Mega or Ugene, however for this occasion we are going to use two software tools.

First, we must generate a matrix that concatenates the previously obtained alignment, for this we will use the [FASconCAT](https://github.com/PatrickKueck/FASconCAT) tool.
>**FASconCAT** is a software tool that is widely utilized in phylogenetic analyses to merge multiple gene alignments into a unified supermatrix alignment. Its primary purpose is to combine individual gene alignments into a comprehensive alignment, which can be further analyzed using methods like maximum likelihood or Bayesian inference. By consolidating the alignments, FASconCAT facilitates the examination of evolutionary relationships and the extraction of meaningful insights from the data.

:warning:However, prior to executing the command, it is essential to ensure that our .fasta file containing the SNP locations and the script for the tool are isolated in a dedicated folder.

Now that everything is organized, we can proceed by executing the following command to open the FASconCAT interface:

`./FASconCAT-G_v1.05.1.pl`

In this way different options are displayed, this time we are going to focus on configuring the command with the following settings: 

1) Nexus: Blocks `n`
2) Phylip: Relaxed `p`
and `s` to run.

🚨The output file typically contains the concatenated sequences from all the input alignments in several formats. !Phylip (PHYLogenetic Inference Package) format is a widely used file format for representing aligned sequences in phylogenetic analyses.!

Thus, with our matrix in .phy format we proceed to generate the tree with the [IQ-TREE](http://www.iqtree.org/) tool. 

>**IQ-TREE** is a powerful software tool for phylogenetic analysis, providing fast and reliable maximum likelihood (ML) tree estimation. It employs the IQ-TREE algorithm, a fast likelihood optimization method, to explore tree space and find the best-fit tree. IQ-TREE supports large datasets and conducts ultrafast bootstrap analyses to assess statistical support. It also offers model selection based on criteria like AIC and BIC, allowing users to choose the optimal substitution model. With its advanced features and efficient analysis of phylogenetic relationships, IQ-TREE is widely used in evolutionary biology.

After we have the IQ-TREE installed and attached to a path command, we proceed to run the command for our matrix as follows:

`module load iqtree/1.6.12`

`iqtree -s FcC_supermatrix.phy -m TEST -bb 1000 -pre R143_MLtree`

🚨IQ-TREE provides several key results for phylogenetic analysis. It generates a phylogenetic tree representing evolutionary relationships among the input sequences, which can be visualized in Newick tree format. The tool performs model selection to identify the best-fit substitution model, crucial for accurate inference. Branch support values indicate statistical reliability, estimated through bootstrap analysis. A log file contains detailed information about the analysis, including optimization progress and model parameters. IQ-TREE can generate a consensus tree summarizing topology and support values from bootstrap replicates. Finally, it reports estimated model parameters like substitution rates and equilibrium frequencies. !In this case we will focus on using the .treefile format, which can be used in several phylogenetic tree visualization tools¡

>Another tip is that to visualize the bootstrap values and root our tree, a very useful and quite customizable tool is [ITOL](https://itol.embl.de/), for example we can generate clade groupings according to our criteria and inferred phylogenetic relationships: 

![rect_tree_154](https://github.com/StivennGutierrez/BIO-analysis/assets/128840301/25ca0d00-b5f6-47c9-92a1-6db151053306)

## SNP Scoreboard: Tracking Genetic Changes in _T. cruzi_ 🌆
Counting SNPs per sequence in _T. cruzi_ is a powerful tool for characterizing genetic diversity, understanding transmission dynamics, and supporting epidemiological surveillance efforts. This information is invaluable for public health interventions, outbreak investigations, and the development of targeted control strategies.

In our terminal it is first important to load the VCFtools module:

`module load vcftools`

then if we evaluate the SNP statistics from the previously obtained **.vcf file**:

`vcf-stats R143_SNPsites.vcf > R143_SNPsites.txt`

Finally, from the .txt file, extract the information on the total number of SNPs per sequence along with their respective identification names. Subsequently, organize this data into an Excel spreadsheet for visualizing the SNP count per sequence.

could be plotted as follows:

![barchart](https://github.com/StivennGutierrez/BIO-analysis/assets/128840301/998b9a39-0bc3-4a5b-8c05-542d00cee4fc)
