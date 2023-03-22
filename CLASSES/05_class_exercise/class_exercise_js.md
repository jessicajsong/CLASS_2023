class\_exercise
================
Jessica Song
3/17/2023

# Load the libraries you need

# Load functions you need “my\_class\_functions”

# load in your peak files for each replicate of each protein

## import peaks

# Here I am starting to analyze my data for my proteins of interest:

ARID3A ATF3 BHLHE40 CEPBP H3K27me3 HDAC2 SIN3A SIN3B SUZ12 TBL1XR1

# First I will read in each replicate file

``` r
# import files
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/jeso3380"
peak_path <- "CLASS_2023/CLASSES/03_Nextflow/01_my_chipseq/results/bwa/mergedLibrary/macs/broadPeak"
broadpeakfilepath <- file.path(basepath, peak_path)
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

# print out a table of the number of peaks in each file:
peak_numbers <- sapply(peak_list, length) %>% as.data.frame(row.names = T)
```

    ## Warning in as.data.frame.integer(., row.names = T): 'row.names' is not a
    ## character vector of length 30 -- omitting it. Will be an error!

``` r
names(peak_numbers) <- c("number_of_peaks")
peak_numbers <- peak_numbers %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")
write_csv(peak_numbers, "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/number_peaks_df.csv")
```

# Now I am going to create consensus peaks for each protein

``` r
dbps <- unique(sapply(names(peak_list), function(x) {
   unlist(strsplit(x, "_"))[1]
}))

consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)
names(consensus_list) <- dbps

sapply(consensus_list, length)
```

    ##   ARID3A     ATF3  BHLHE40    CEBPB H3K27me3    HDAC2    SIN3A    SIN3B 
    ##    23669     2876     1234    12519      209    26326     9458     5882 
    ##    SUZ12  TBL1XR1 
    ##     1235    12960

``` r
number_consensus_peaks <- sapply(consensus_list, length) %>% 
  as.data.frame() %>%
  rownames_to_column( var = "dbp") %>%
  dplyr::rename(number_consensus_peaks = ".")

peak_num <- left_join(peak_numbers, number_consensus_peaks)
```

    ## Joining with `by = join_by(dbp)`

``` r
write_csv(peak_numbers, "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/number_peaks_df.csv")

# export consensus peaks to results folder
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/jeso3380"
consensus_path <- "/CLASS_2023/CLASSES/05_class_exercise/results/"
exportpath <- file.path(basepath, consensus_path)

for(i in 1:length(consensus_list)) {
rtracklayer::export(consensus_list[[i]], paste0(exportpath, names(consensus_list)[i], "_consensus_peaks.bed") )}
```

# Now I am going to make my consensus peaks compatible with UCSC genome browser

``` r
consensus_file_list <- list.files("/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results", full.names = T, pattern = ".bed")

peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))

names(peaks) <- dbps

canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")

peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))

new_filenames <- paste0("/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/00_consensus_peaks/", names(peaks), "_consensus.bed")

for(i in 1:length(peaks)) {
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}

headers <- paste0("track type=bed name=", names(peaks))
headers
```

    ##  [1] "track type=bed name=ARID3A"   "track type=bed name=ATF3"    
    ##  [3] "track type=bed name=BHLHE40"  "track type=bed name=CEBPB"   
    ##  [5] "track type=bed name=H3K27me3" "track type=bed name=HDAC2"   
    ##  [7] "track type=bed name=SIN3A"    "track type=bed name=SIN3B"   
    ##  [9] "track type=bed name=SUZ12"    "track type=bed name=TBL1XR1"

``` r
new_filenames <- paste0("/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/", names(peaks), ".bed")
new_filenames
```

    ##  [1] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/ARID3A.bed"  
    ##  [2] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/ATF3.bed"    
    ##  [3] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/BHLHE40.bed" 
    ##  [4] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/CEBPB.bed"   
    ##  [5] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/H3K27me3.bed"
    ##  [6] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/HDAC2.bed"   
    ##  [7] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/SIN3A.bed"   
    ##  [8] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/SIN3B.bed"   
    ##  [9] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/SUZ12.bed"   
    ## [10] "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/results/UCSC_consensus_peaks/TBL1XR1.bed"

``` r
# print out consensus peak files in a results/UCSC directory
for(i in 1:length(peaks)) {
  writeLines(headers[[i]], new_filenames[[i]])
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}

for(i in 1:length(peaks)) {
  writeLines(headers[[i]], new_filenames[[i]])
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}
```

# Now I want to compare a protein with a previous analysis

``` r
# sorry John I moved things around from this point on

# I loaded in my SUZ12 consensus peak file and the consensus peaks match up exactly with the SUZ12 consensus peaks from 2021!
```

# Now I am going to determine how my peaks for each protein overlap annotations of the genome

# First I will find the overlaps between my consensus peaks with promoters of lncRNA and mRNA promoters

``` r
# reset file paths
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/jeso3380"
peak_path <- "CLASS_2023/CLASSES/05_class_exercise/00_consensus_peaks"
consensusPeakPath <- file.path(basepath, peak_path)
setwd("/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features")

# get the file list
  
consensus_peaks_files <- list.files(consensusPeakPath, 
                                             pattern = "*.bed",
                                             full.names = TRUE)

# make Granges, add DBP names, clean up the file names
consensus_peaks <- lapply(consensus_peaks_files, rtracklayer::import)
names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/00_consensus_peaks/|_consensus.bed","", consensus_peaks_files)

#load in gencode, get list of mRNA and lncRNA genes
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
table(gencode_gr$type)
```

    ## 
    ##           gene     transcript           exon            CDS    start_codon 
    ##          60609         227462        1372308         761508          87662 
    ##     stop_codon            UTR Selenocysteine 
    ##          79913         310193            119

``` r
rtracklayer::export(gencode_genes, "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/gene_annotations/gencode_genes.gtf")

mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"] 
rtracklayer::export(mrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/gene_annotations/mrna_genes.gtf")
table(gencode_genes$gene_type)
```

    ## 
    ##                          IG_C_gene                    IG_C_pseudogene 
    ##                                 14                                  9 
    ##                          IG_D_gene                          IG_J_gene 
    ##                                 37                                 18 
    ##                    IG_J_pseudogene                      IG_pseudogene 
    ##                                  3                                  1 
    ##                          IG_V_gene                    IG_V_pseudogene 
    ##                                144                                188 
    ##                             lncRNA                              miRNA 
    ##                              16849                               1881 
    ##                           misc_RNA                            Mt_rRNA 
    ##                               2212                                  2 
    ##                            Mt_tRNA             polymorphic_pseudogene 
    ##                                 22                                 42 
    ##               processed_pseudogene                     protein_coding 
    ##                              10171                              19965 
    ##                         pseudogene                           ribozyme 
    ##                                 18                                  8 
    ##                               rRNA                    rRNA_pseudogene 
    ##                                 52                                500 
    ##                             scaRNA                              scRNA 
    ##                                 49                                  1 
    ##                             snoRNA                              snRNA 
    ##                                942                               1901 
    ##                               sRNA                                TEC 
    ##                                  5                               1061 
    ##                          TR_C_gene                          TR_D_gene 
    ##                                  6                                  4 
    ##                          TR_J_gene                    TR_J_pseudogene 
    ##                                 79                                  4 
    ##                          TR_V_gene                    TR_V_pseudogene 
    ##                                106                                 33 
    ##   transcribed_processed_pseudogene     transcribed_unitary_pseudogene 
    ##                                495                                130 
    ## transcribed_unprocessed_pseudogene    translated_processed_pseudogene 
    ##                                923                                  2 
    ##  translated_unprocessed_pseudogene                 unitary_pseudogene 
    ##                                  2                                 98 
    ##             unprocessed_pseudogene                           vaultRNA 
    ##                               2631                                  1

``` r
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 
rtracklayer::export(lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/gene_annotations/lncrna_genes.gtf")

mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
rtracklayer::export(mrna_lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/gene_annotations/mrna_lncrna_genes.gtf")

# annotation file for mRNA and lncRNA + set it up as a data frame
lncrna_mrna_genes <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/gene_annotations/mrna_lncrna_genes.gtf")
lncrna_mrna_genes_df <- lncrna_mrna_genes %>% as.data.frame()

# annotation file for promoters of mRNA and lncRNA genes
lncrna_mrna_promoters <- promoters(lncrna_mrna_genes, upstream = 1000, downstream = 1000)

rtracklayer::export(lncrna_mrna_promoters, "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/gene_annotations/lncrna_mrna_promoters.gtf")

lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
table(mrna_lncrna_genes$gene_type)
```

    ## 
    ##         lncRNA protein_coding 
    ##          16849          19965

``` r
mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]

# Looking at overlaps

# Load in peaks
number_peaks_df <- data.frame("dbp" = names(consensus_peaks),
                           "num_peaks" = sapply(consensus_peaks, length))

# How much of the genome is covered by all peaks for a DBP
number_peaks_df$total_peak_length <- sapply(consensus_peaks, function(x) sum(width(x)))

# Use count_peaks_per_feature function to find number of overlaps at each promoter (cols=promoters, rows=DPS)
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, type = "counts")

#row sum for each DBP
number_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

# peaks overlapping promoters for each protein
number_peaks_df
```

    ##               dbp num_peaks total_peak_length peaks_overlapping_promoters
    ## ARID3A     ARID3A    355035         139753095                       69015
    ## ATF3         ATF3     43140          45746805                       39885
    ## BHLHE40   BHLHE40     18510           7227585                        9420
    ## CEBPB       CEBPB    187785          61641225                       26130
    ## H3K27me3 H3K27me3      3762          13958532                        2772
    ## HDAC2       HDAC2    473868         598874652                      200214
    ## SIN3A       SIN3A    170244         126087462                      156780
    ## SIN3B       SIN3B    105876          65013840                      120258
    ## SUZ12       SUZ12     22230           7770042                       24966
    ## TBL1XR1   TBL1XR1    233280         105279300                      136854

## results:

\#1) What can you determine from these overlaps? \# ATF3, SIN3A, SIN3B,
and SUZ12 bind promoters pretty heavily compared to the other proteins.

# Now I want to compare the overlaps with lncRNA and mRNA promoters seperately

``` r
# Break up promoters into two groups - lncRNA and mRNA by indexing
number_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])
number_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

write_csv(number_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/results/number_peaks_df.csv")

# peaks overlapping mRNA vs lncRNA promoters for each protein
number_peaks_df
```

    ##               dbp num_peaks total_peak_length peaks_overlapping_promoters
    ## ARID3A     ARID3A    355035         139753095                       69015
    ## ATF3         ATF3     43140          45746805                       39885
    ## BHLHE40   BHLHE40     18510           7227585                        9420
    ## CEBPB       CEBPB    187785          61641225                       26130
    ## H3K27me3 H3K27me3      3762          13958532                        2772
    ## HDAC2       HDAC2    473868         598874652                      200214
    ## SIN3A       SIN3A    170244         126087462                      156780
    ## SIN3B       SIN3B    105876          65013840                      120258
    ## SUZ12       SUZ12     22230           7770042                       24966
    ## TBL1XR1   TBL1XR1    233280         105279300                      136854
    ##          peaks_overlapping_lncrna_promoters peaks_overlapping_mrna_promoters
    ## ARID3A                                18885                            50130
    ## ATF3                                   7980                            31905
    ## BHLHE40                                2385                             7035
    ## CEBPB                                  7935                            18195
    ## H3K27me3                                810                             1962
    ## HDAC2                                 49158                           151056
    ## SIN3A                                 31716                           125064
    ## SIN3B                                 22284                            97974
    ## SUZ12                                  3942                            21024
    ## TBL1XR1                               30078                           106776

## results:

# 1) What is the difference in overlaps between mRNA and lncRNA promoters

The overlaps between mRNA and lncRNA promoters differs by a factor of 2,
which makes sense since there are two times more mRNA promoters compared
to lncRNA promoters. This trend holds true for the proteins I looked at,
except for ATF3, SIN3A, SIN3B, SUZ12, TBL1XR1 which seem to bind mRNA
promoters more than lncRNA promoters.

# I am curious if my proteins are transcription factors so I will use the annotations

# in a cell paper I found and see

``` r
url <- "https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx"
destination_for_url <- "/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/results/TF_annotations.xlsx"
download.file(url, destination_for_url)

human_tfs <- readxl::read_excel("/scratch/Shares/rinnclass/CLASS_2023/jeso3380/CLASS_2023/CLASSES/05_class_exercise/01_peak_features/results/TF_annotations.xlsx",
                                sheet = 2, skip = 1)
```

    ## Warning: Expecting logical in M1006 / R1006C13: got 'Contains a SANT and
    ## multiple DNA-binding C2H2 domains. Motif is 99% AA ID from mouse (Transfac).'

    ## Warning: Expecting logical in M1021 / R1021C13: got 'Close ortholog (PP1RA)
    ## binds to mRNA; single-stranded DNA (ssDNA); poly(A) and poly(G) homopolymers
    ## (Uniprot)'

    ## Warning: Expecting logical in M1542 / R1542C13: got 'Contains 1 SANT domain'

    ## Warning: Expecting logical in M1543 / R1543C13: got 'Contains 2 Myb DBDs.
    ## Sources of Hocomoco/Transfac motifs are unclear. However these sequences look
    ## similar to in vitro sites selected by SELEX (PMID:11082045)'

    ## Warning: Expecting logical in M1544 / R1544C13: got 'Although CHD2 has weak
    ## similarity to a Myb domain (PMID:9326634), it's more closely related to the
    ## non-DNA-binding SANT domain based on our alignment analysis. The data showing
    ## that show that CHD2 binding histone H3.3 (PMID:22569126) further support the
    ## conclusion that the Myb domain is probably a SANT domain facilitating the
    ## histone interaction'

    ## Warning: Expecting logical in M1545 / R1545C13: got 'Contains a single SANT
    ## domain, no evidence for sequence-specific DNA binding'

    ## Warning: Expecting logical in M1546 / R1546C13: got 'Contains 2 Myb DBDs'

    ## Warning: Expecting logical in M1547 / R1547C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1548 / R1548C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1549 / R1549C13: got 'Contains a single SANT
    ## domain, no evidence for sequence-specific DNA binding'

    ## Warning: Expecting logical in M1550 / R1550C13: got 'Domain is truncated, and
    ## there is nothing known about this gene'

    ## Warning: Expecting logical in M1551 / R1551C13: got 'Contains a single SANT
    ## domain, no evidence for sequence-specific DNA binding'

    ## Warning: Expecting logical in M1552 / R1552C13: got 'MIER2's Myb domain is more
    ## similar to the non-DNA-binding SANT domain'

    ## Warning: Expecting logical in M1553 / R1553C13: got 'MIER3's Myb domain is more
    ## similar to the non-DNA-binding SANT domain'

    ## Warning: Expecting logical in M1554 / R1554C13: got 'Contains 1 SANT domain,
    ## and a SANTA domain'

    ## Warning: Expecting logical in M1555 / R1555C13: got 'Contains a single Myb-like
    ## domain with an insertion in the middle. It is ambiguous whether Myb-like
    ## domains are DNA or protein binding. Since it has a single domain it's likely
    ## non-specific, but future experiments should be performed to assay it's
    ## specificity'

    ## Warning: Expecting logical in M1556 / R1556C13: got 'Contains 3 Myb DBDs'

    ## Warning: Expecting logical in M1557 / R1557C13: got 'Contains 3 Myb DBDs'

    ## Warning: Expecting logical in M1558 / R1558C13: got 'Contains 3 Myb DBDs'

    ## Warning: Expecting logical in M1559 / R1559C13: got 'Contains a single Myb-like
    ## domain. Mouse ortholog has motif'

    ## Warning: Expecting logical in M1560 / R1560C13: got 'MYSM1 has been shown to
    ## bind DNA ? interaction with DNA requires the MYSM1 Myb but not the SWIRM domain
    ## (PMID:17428495). Domain sequence alignment places it near DNA-binding Myb
    ## domains but scores slightly higher as a SANT rather than Myb domain based on
    ## Prosite patterns. Given that most Myb proteins that bind DNA sequence
    ## specifically have multiple Myb domains in an array this protein could bind DNA
    ## sequence non-specifically with it?s single Myb domain. Future experiments
    ## should assay MYSM1?s specificity'

    ## Warning: Expecting logical in M1561 / R1561C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1562 / R1562C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1564 / R1564C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1565 / R1565C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains'

    ## Warning: Expecting logical in M1566 / R1566C13: got 'Contains 2 SANT domains,
    ## and no other putative DNA-binding domains. RCOR3 SANT domains are known to
    ## facilitate PPIs'

    ## Warning: Expecting logical in M1567 / R1567C13: got 'SMARCA1 contains a
    ## truncated Myb-like and SANT domain. Given the presence of the Myb-like domain,
    ## and other domains known to associated with DNA (DEAD box helicase) it likely
    ## associates with DNA non-sequence-specifically'

    ## Warning: Expecting logical in M1568 / R1568C13: got 'Contains a SANT, and
    ## Myb-like domain'

    ## Warning: Expecting logical in M1569 / R1569C13: got 'Contains 1 SANT domain,
    ## and no other putative DNA-binding domains. Motif logos look like bZIP dimeric
    ## binding sites, and are thus likely specificifities of SMARCC1 interactors'

    ## Warning: Expecting logical in M1570 / R1570C13: got 'Contains 1 SANT domain,
    ## and no other putative DNA-binding domains. Motif logos ares likely
    ## specificifities of SMARCC2 interactors'

    ## Warning: Expecting logical in M1571 / R1571C13: got 'Contains only Myb DBDs'

    ## Warning: Expecting logical in M1572 / R1572C13: got 'Contains 1 SANT domain'

    ## Warning: Expecting logical in M1573 / R1573C13: got 'TADA2B contains a single
    ## SANT domain and is thus unlikely to bind DNA'

    ## Warning: Expecting logical in M1574 / R1574C13: got 'Contains a single Myb
    ## domain (with slightly less simialrity to a SANT domain.) This domain has been
    ## shown to be involved in PPIs but this may not be mutually exclusive with
    ## DNA-binding. The sequence-specificity of CCDC79 should be investigated in the
    ## future'

    ## Warning: Expecting logical in M1575 / R1575C13: got 'Contains 1 Myb domain, and
    ## has structural evidence of DNA-binding'

    ## Warning: Expecting logical in M1576 / R1576C13: got 'Motif is inferred from
    ## mouse (92% DBD AA ID)'

    ## Warning: Expecting logical in M1577 / R1577C13: got 'TERF2IP contains a single
    ## Myb-like domain. While it's unclear if TERF2IP (Human Rap1) contacts DNA
    ## directly it has been shown to affect the DNA binding activity of TRF2'

    ## Warning: Expecting logical in M1578 / R1578C13: got 'This protein contains Myb,
    ## and Myb-like domains and is annotated as a Pol1 terminator. TTF1 DNA-binding
    ## has been demonstrated in vitro (PMID: 7597036), but it's specificity has not
    ## been determined'

    ## Warning: Expecting logical in M1579 / R1579C13: got 'Contains 1 Myb DBD'

    ## Warning: Expecting logical in M1580 / R1580C13: got 'Contains a GATA and SANT
    ## domain. Unclear whether the GATA domain is a bona fide DBD as the MTA/RERE
    ## family domains are atypical to human GATA domains (see alignment). In CIS-BP
    ## there is one protein from C.elegans that shares domain homology and binds a
    ## GATA motif (elg-27, ChIP-seq). The GATA ZnF domain of MTA1 is required for it's
    ## interaction with RBBP4 and RBBP7 (PMID:18067919). Full-length protein has been
    ## tried in HT-SELEX and did not yield a motif'

    ## Warning: Expecting logical in M1581 / R1581C13: got 'Contains a GATA and SANT
    ## domain. Unclear whether the GATA domain is a bona fide DBD as the MTA/RERE
    ## family domains are atypical to human GATA domains (see alignment). In CIS-BP
    ## there is one protein from C.elegans that shares domain homology and binds a
    ## GATA motif (elg-27, ChIP-seq). Full-length protein has been tried in HT-SELEX,
    ## and DBD has been tried on PBM - neither yielded motifs'

    ## Warning: Expecting logical in M1582 / R1582C13: got 'Contains a GATA and SANT
    ## domain. Unclear whether the GATA domain is a bona fide DBD as the MTA/RERE
    ## family domains are atypical to human GATA domains (see alignment). In CIS-BP
    ## there is one protein from C.elegans that shares domain homology and binds a
    ## GATA motif (elg-27, ChIP-seq). Hasn't been tried in any in vitro assays'

    ## Warning: Expecting logical in M1583 / R1583C13: got 'Contains a GATA and SANT
    ## domain. Unclear whether the GATA domain is a bona fide DBD as the MTA/RERE
    ## family domains are atypical to human GATA domains (see alignment). In CIS-BP
    ## there is one protein from C.elegans that shares domain homology and binds a
    ## GATA motif (elg-27, ChIP-seq). Has been tried as a DBD in HT-SELEX but did not
    ## yield a motif'

    ## Warning: Expecting logical in M1791 / R1791C13: got 'CNOT3 is a part of the
    ## CCR4-NOT complex involved in mRNA decay'

    ## Warning: Expecting logical in M1932 / R1932C13: got '"Prosite identifies a
    ## low-confidence Myb-like domain (e.g. can?t decide between Myb and SANT) so it?s
    ## probably not a TF"'

    ## New names:
    ## • `` -> `...4`

``` r
# rename 4th column of human_tfs
names(human_tfs)[4] <- "is_tf"

# intersect between my dbp names and if it's a TF
length(which(tolower(number_peaks_df$dbp) %in% tolower(human_tfs$Name)))
```

    ## [1] 8

``` r
# merge the TF annotation file to the number_peaks_df

human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(number_peaks_df$dbp), 1:4]

names(human_tfs) <- c("ensembl_id",
                      "dbp",
                      "dbd",
                      "tf")

number_peaks_df <- merge(number_peaks_df, human_tfs, all.x = T)
dim(number_peaks_df[is.na(number_peaks_df$tf),])
```

    ## [1] 2 9

``` r
number_peaks_df <- number_peaks_df[1:9]
write_csv(number_peaks_df, "results/number_peaks_df.csv")
number_peaks_df
```

    ##         dbp num_peaks total_peak_length peaks_overlapping_promoters
    ## 1    ARID3A    355035         139753095                       69015
    ## 2      ATF3     43140          45746805                       39885
    ## 3   BHLHE40     18510           7227585                        9420
    ## 4     CEBPB    187785          61641225                       26130
    ## 5  H3K27me3      3762          13958532                        2772
    ## 6     HDAC2    473868         598874652                      200214
    ## 7     SIN3A    170244         126087462                      156780
    ## 8     SIN3B    105876          65013840                      120258
    ## 9     SUZ12     22230           7770042                       24966
    ## 10  TBL1XR1    233280         105279300                      136854
    ##    peaks_overlapping_lncrna_promoters peaks_overlapping_mrna_promoters
    ## 1                               18885                            50130
    ## 2                                7980                            31905
    ## 3                                2385                             7035
    ## 4                                7935                            18195
    ## 5                                 810                             1962
    ## 6                               49158                           151056
    ## 7                               31716                           125064
    ## 8                               22284                            97974
    ## 9                                3942                            21024
    ## 10                              30078                           106776
    ##         ensembl_id         dbd   tf
    ## 1  ENSG00000116017 ARID/BRIGHT  Yes
    ## 2  ENSG00000162772        bZIP  Yes
    ## 3  ENSG00000134107        bHLH  Yes
    ## 4  ENSG00000172216        bZIP  Yes
    ## 5             <NA>        <NA> <NA>
    ## 6  ENSG00000196591     Unknown   No
    ## 7  ENSG00000169375     Unknown   No
    ## 8  ENSG00000127511     Unknown   No
    ## 9  ENSG00000178691     Unknown   No
    ## 10            <NA>        <NA> <NA>

# Now I am going to test if there is more binding over gene bodies than promoters

# I will seperate lncRNA and mRNA gene bodies to find the overlaps

``` r
# Find overlaps with gene_bodies
genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                consensus_peaks, 
                                                type = "counts")
number_peaks_df$peaks_overlapping_genebody <- 
  rowSums(genebody_peak_counts)

# lncRNA gene bodies 
number_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA gene bodies
number_peaks_df$peaks_overlapping_mrna_genebody <- 
  rowSums(genebody_peak_counts[,mrna_gene_ids])

write_csv(number_peaks_df, "results/number_peaks_df.csv")
number_peaks_df
```

    ##         dbp num_peaks total_peak_length peaks_overlapping_promoters
    ## 1    ARID3A    355035         139753095                       69015
    ## 2      ATF3     43140          45746805                       39885
    ## 3   BHLHE40     18510           7227585                        9420
    ## 4     CEBPB    187785          61641225                       26130
    ## 5  H3K27me3      3762          13958532                        2772
    ## 6     HDAC2    473868         598874652                      200214
    ## 7     SIN3A    170244         126087462                      156780
    ## 8     SIN3B    105876          65013840                      120258
    ## 9     SUZ12     22230           7770042                       24966
    ## 10  TBL1XR1    233280         105279300                      136854
    ##    peaks_overlapping_lncrna_promoters peaks_overlapping_mrna_promoters
    ## 1                               18885                            50130
    ## 2                                7980                            31905
    ## 3                                2385                             7035
    ## 4                                7935                            18195
    ## 5                                 810                             1962
    ## 6                               49158                           151056
    ## 7                               31716                           125064
    ## 8                               22284                            97974
    ## 9                                3942                            21024
    ## 10                              30078                           106776
    ##         ensembl_id         dbd   tf peaks_overlapping_genebody
    ## 1  ENSG00000116017 ARID/BRIGHT  Yes                     297705
    ## 2  ENSG00000162772        bZIP  Yes                      55785
    ## 3  ENSG00000134107        bHLH  Yes                      17430
    ## 4  ENSG00000172216        bZIP  Yes                     150975
    ## 5             <NA>        <NA> <NA>                       3978
    ## 6  ENSG00000196591     Unknown   No                     483570
    ## 7  ENSG00000169375     Unknown   No                     219150
    ## 8  ENSG00000127511     Unknown   No                     142344
    ## 9  ENSG00000178691     Unknown   No                      28494
    ## 10            <NA>        <NA> <NA>                     242244
    ##    peaks_overlapping_lncrna_genebody peaks_overlapping_mrna_genebody
    ## 1                              67905                          229800
    ## 2                              10830                           44955
    ## 3                               3495                           13935
    ## 4                              34365                          116610
    ## 5                               1314                            2664
    ## 6                             104184                          379386
    ## 7                              37998                          181152
    ## 8                              23220                          119124
    ## 9                               4716                           23778
    ## 10                             50148                          192096

## results:

# 1) Do my proteins have more overlaps with promoters or genebodies?

ARID3A, CEBPB, HDAC2, and TBL1XR1 seem to overlap with gene bodies more
than promoters, whereas ATF3, H3K27me3, SIN3A, SIN3B and SUZ12 bind to
promoters and gene bodies about the same.

# It is nice and all to find overlaps, but I am interested in how many proteins

# bind a specific promoter. I will use my handy “occurence” parameter in

# " count peaks per feature"

``` r
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, 
                                               type = "occurrence")
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))
write.table(promoter_peak_occurence, "results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))
write_csv(peak_occurence_df, "results/peak_occurence_dataframe.csv")
```

## results: I find the max number of proteins on a promoter to be 8

# Now I want to start plotting my results

# First I will see if there is a realtionship between peak number and total DNA covered

``` r
ggplot(number_peaks_df, aes(x = num_peaks, 
                         y = total_peak_length)) +
  geom_point() 
```

![](class_exercise_js_files/figure-gfm/total%20peak%20length%20vs%20number%20peaks-1.png)<!-- -->

# Now I want to color my plot by whether the protein is a TF or not.

``` r
ggplot(number_peaks_df, aes(x = num_peaks, 
                         y = total_peak_length,
                        color = tf )) +
  geom_point()
```

![](class_exercise_js_files/figure-gfm/total%20peak%20length%20vs%20number%20peaks,%20+%20tf-1.png)<!-- -->

# I want to make a histogram of the number of peaks for each of my proteins

``` r
ggplot(number_peaks_df, aes(x = num_peaks)) +
  geom_histogram(bins = 80)
```

![](class_exercise_js_files/figure-gfm/number%20of%20peaks%20for%20each%20protein%20-%20histogram-1.png)<!-- -->

# Now I want to facet this by the type of DNA binding domain my protein has.

``` r
ggplot(number_peaks_df, aes(x = num_peaks, fill = dbd)) +
  geom_histogram(bins = 80)
```

![](class_exercise_js_files/figure-gfm/number_peaks%20w%20dbd%20-%20histagram-1.png)<!-- -->

# Cool now I am ready to send my result to my collaborator as a

# Knitted document
