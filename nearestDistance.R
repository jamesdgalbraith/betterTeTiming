# Check if MAFFT and BLAST are installed
mafft_path_check <- Sys.which("mafft")
if (mafft_path_check != "") {
  mafft_path <- unname(mafft_path_check) 
} else {
  stop("MAFFT not found in PATH, please ensure MAFFT is installed and in your PATH")
}
blast_path_check <- Sys.which("blastn")
if (blast_path_check != "") {
  blast_path <- unname(blast_path_check) 
} else {
  stop("BLAST+ not found in PATH, please ensure NCBI BLAST+ is installed and in your PATH")
}

# Read in options
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-g", "--genome"), default=NA, type = "character",
              help="Path to genome FASTA file"),
  make_option(c("-l", "--library"), default=NA, type = "character",
              help="Path to EarlGrey repeat FASTA file"),
  make_option(c("-a", "--annotation_gff"), default=NA, type = "character",
              help="Path to EarlGrey gff"),
  make_option(c("-o", "--outdir"), default=NA, type = "character",
              help="Output directory"),
  make_option(c("-t", "--threads"), default=1, type = "integer",
              help="Prefix to add to gene names"),
  make_option(c("-L", "--min_len"), default=250, type = "integer",
              help="Minimum length of repeats for a superfamily to be considered"),
  make_option(c("-N", "--min_n"), default=250, type = "integer",
              help="Minimum number or repeats for superfamily to be considered"),
  make_option(c("-C", "--min_cov"), default=80, type = "integer",
              help="Minimum coverage a repeat must have from another repeat to be included")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check required variables provided
if (!is.na(opt$genome)) {
  genome <- opt$genome
} else {
  stop("Path to genome assembly FASTA must be provided. See script usage (--help)")
}

if (!is.na(opt$library)) {
  library <- opt$library
} else {
  stop("Path to repeat library FASTA must be provided. See script usage (--help)")
}

if (!is.na(opt$annotation_gff)) {
  annotation_gff <- opt$annotation_gff
} else {
  stop("Path to Earl Grey GFF must be provided. See script usage (--help)")
}

if (!is.na(opt$outdir)) {
  outdir <- opt$outdir
  if (!dir.exists(outdir)){
    dir.create(outdir, recursive = T)
  }
} else {
  stop("Path to out directory must be provided. See script usage (--help)")
}

# Import packages
suppressMessages(library(tidyverse))
suppressMessages(library(purrr))
suppressMessages(library(parallel))
suppressMessages(library(plyranges))
suppressMessages(library(Biostrings))
suppressMessages(library(ape))
suppressMessages(library(ips))

# Function for filtering repeats and extracting sequences from genomes
gff2seq <- function(genome_path, library_path, gff_path, min_len, min_n){
  
  # read in genome
  message("Reading genome")
  genome_seq <- readDNAStringSet(genome_path)
  names(genome_seq) <- sub(" .*", "", names(genome_seq))
  # Make genome index
  genome_idx <- tibble(seqnames = names(genome_seq), start = 1, end = width(genome_seq)) %>%
    as_granges()
  message("Reading gff")
  # Read in TEs, fix Penelopes, select necessary columns, 
  eg_gff <- plyranges::read_gff(gff_path) %>%
    plyranges::filter(seqnames %in% seqnames(genome_idx)) %>%
    dplyr::select(type, ID) %>%
    dplyr::mutate(type = ifelse(type == "LINE/Penelope", "PLE/Penelope", as.character(type))) %>%
    dplyr::mutate(subclass = sub("/.*", "", type),
                  superfamily = sub("-.*", "", sub(".*/", "", type))) %>%
    dplyr::filter(subclass != "Other") %>%
    dplyr::filter(subclass %in% c("LINE", "LTR", "DNA", "PLE"),  # just for testing, alter when complete
                  width >=min_len) %>%
    dplyr::mutate(ID = tolower(ID))
  message("Reading library")
  # read in Earl Grey families, remove any not in GFF
  eg_seq <- Biostrings::readDNAStringSet(library_path)
  names(eg_seq) <- tolower(sub("#.*", "", names(eg_seq)))
  eg_seq <- eg_seq[names(eg_seq) %in% eg_gff$ID,]
  message("Filtering gff")
  # filter out TEs over 150% size of consensus sequences, may alter for customisablity
  filtered_gff <- eg_gff %>%
    dplyr::as_tibble() %>%
    dplyr::inner_join(tibble(ID = names(eg_seq), con_width = width(eg_seq)), by = "ID") %>%
    dplyr::filter(width <= 1.5* con_width)
  
  # identify supefamilies present and numerous for calculation
  sf_n <- as_tibble(as.data.frame(table(eg_gff$type))) %>%
    dplyr::rename(type = Var1) %>%
    dplyr::mutate(type = as.character(type)) %>%
    arrange(-Freq) %>%
    filter(Freq >= min_n)
  message("Getting repeat seq")
  all_seq <- BSgenome::getSeq(genome_seq, eg_gff)
  names(all_seq) <- paste0(seqnames(eg_gff), ":", ranges(eg_gff))
  
  filtered_gff <- filtered_gff %>%
    as_tibble() %>%
    dplyr::mutate(names = paste0(seqnames, ":", start, "-", end))
  
  return(list(filtered_gff, all_seq, sf_n))
}

# Function for performing all vs all blast for each superfamily and filtering results
getSeqBlastn <- function(eg_gff, eg_seq, subclass_n, outdir, n_threads, min_cov){
  
  comp_blast <- dplyr::tibble()
  comp_seq <- Biostrings::DNAStringSet()
  if(!dir.exists(paste0(outdir, "/db/"))){
    dir.create(paste0(outdir, "/db/"), recursive = T)
  }
  
  
  for (i in seq_along(subclass_n$type)) {
    
    message(paste0(subclass_n$type[i], " ", i, " of ", nrow(subclass_n)))
    # per superfamily, testing on CR1s
    # get sequence
    message("\t getting seq ")
    subset_eg_seq <- eg_seq[names(eg_seq) %in% eg_gff[eg_gff$type == subclass_n$type[i],]$names]
    comp_seq <- c(comp_seq, subset_eg_seq)
    
    # assign db path and check dir exists
    db_path <- paste0(outdir, "/db/", gsub("\\/", "_", subclass_n$type[i]))
    # write sequence to file
    writeXStringSet(x = subset_eg_seq,
                    filepath =  paste0(db_path, ".fasta"))
    system(paste0("makeblastdb -in ", db_path, ".fasta -dbtype nucl"), intern = TRUE)
    # perform all vs all blast
    message("\t BLASTing ")
    system(paste0("blastn -query ", db_path, ".fasta",
                  " -db ", db_path, ".fasta",
                  " -outfmt \"6 qseqid sseqid qstart qend sstart send nident mismatch qcovs length qlen slen\"",
                  " -max_target_seqs 6 -out ", db_path, ".blastn -task dc-megablast -num_threads ", n_threads))
    # blast, take top hits, remove self, ensure coverage >0.5, calculate dist
    message("\t processing blast out")
    blast_out <- read_tsv(paste0(db_path, ".blastn"),
                          col_names = c("qseqid", "sseqid", "qstart", "qend",
                                        "sstart", "send", "match", "mismatch",
                                        "qcovs", "length", "qlen", "slen"), show_col_types = F) %>%
      filter(qseqid != sseqid) %>%
      dplyr::filter(qcovs >= min_cov) %>%
      dplyr::group_by(qseqid) %>%
      dplyr::slice(1)
    
    # compile blast data
    comp_blast <- rbind(blast_out, comp_blast)
    
    # remove db and fasta
    message("\t cleaning up ")
    delfiles <- dir(path=paste0(outdir, "/db/") ,pattern=paste0(gsub("\\/", "_", subclass_n$type[i]), ".*"))
    file.remove(file.path(paste0(outdir, "/db/"), delfiles))
    
    
  }
  
  return(list(comp_blast, comp_seq))
  
}

# Function for calculating distance
alnDist <- function(seq_in){
  
  # get unaligned sequences and convert to ape
  unaln_ape <- ape::as.DNAbin(seq_in)
  # align with mafft
  aln_ape <- ips::mafft(unaln_ape, method = "localpair", thread = 1, exec = mafft_path, options = "--adjustdirection")
  # calculate distances
  dist_df <- base::data.frame(names = names(seq_in)[1],
                              kdist = ape::dist.dna(aln_ape, model="K80", pairwise.deletion = TRUE)[1],
                              jcdist = ape::dist.dna(aln_ape, model="JC69", pairwise.deletion = TRUE)[1],
                              rawdist = ape::dist.dna(aln_ape, model="raw", pairwise.deletion = TRUE)[1])
  return(dist_df)
}

# Get repeat sequence and classes
repeat_data <- gff2seq(opt$genome, opt$library, opt$annotation_gff, opt$min_len, opt$min_n)

# Run all vs all blast for each superfamily and return compiled top hits and seqs
blastOut <- getSeqBlastn(repeat_data[[1]], repeat_data[[2]], repeat_data[[3]], opt$outdir, opt$threads, opt$min_cov)
comp_blast <- blastOut[[1]]
comp_seq <- blastOut[[2]]

# write sequences and best hit to file
readr::write_tsv(comp_blast, paste0(opt$outdir, "/blast_best_hits.tsv"))
Biostrings::writeXStringSet(comp_seq, paste0(opt$outdir, "/genomic_te_sequences.fasta") )

# Add sequence data to tibble
to_align_tbl <- inner_join(comp_blast, tibble(qseqid = names(comp_seq), qseq = as.character(comp_seq)), by = "qseqid") %>%
  inner_join(tibble(sseqid = names(comp_seq), sseq = as.character(comp_seq)), by = "sseqid")

# make list of pairs to align
message("Compiling sequence pairs")
to_align <- purrr::map(seq_along(to_align_tbl$qseq), ~ {
  DNAStringSet(c(to_align_tbl$qseq[.x], to_align_tbl$sseq[.x]))
})

# FURRR LOOP/FUNCTION TO CALCULATE AND RETURN DISTANCE
message("Aligning and calculating distance")
kdist_list <- mclapply(to_align, alnDist, mc.cores = opt$threads)
kdist_tbl <- purrr::list_rbind(kdist_list) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(kdist = base::round(kdist, 2),
                jcdist = base::round(jcdist, 2),
                rawdist = base::round(rawdist, 2))
readr::write_tsv(kdist_tbl, paste0(opt$outdir, "/final_kdist.tsv"))

# Join with gff and write to file
dplyr::inner_join(repeat_data[[1]], kdist_tbl, by = "names") %>%
  plyranges::select(-names) %>%
  readr::write_tsv(paste0(opt$outdir, "/kdist_", sub(".gff", ".tsv", sub(".*\\/", "", opt$annotation_gff))))