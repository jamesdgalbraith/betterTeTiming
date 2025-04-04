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
    dplyr::as_tibble() %>%
    dplyr::select_if(names(.) %in% c("seqnames", "start", "end", "strand", "type", "ID", "KIMURA80")) %>%
    plyranges::as_granges() %>%
    dplyr::mutate(type = ifelse(type == "LINE/Penelope", "PLE/Penelope", as.character(type))) %>%
    dplyr::mutate(subclass = sub("/.*", "", type),
                  superfamily = sub("-.*", "", sub(".*/", "", type))) %>%
    dplyr::filter(subclass != "Other") %>%
    dplyr::filter(subclass %in% c("LINE", "LTR", "DNA", "PLE", "RC", "SINE"),  # just for testing, alter when complete
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

# Get repeat sequence and classes
repeat_data <- gff2seq(opt$genome, opt$library, opt$annotation_gff, opt$min_len, opt$min_n)

# Run all vs all blast for each superfamily and return compiled top hits and seqs
blastOut <- getSeqBlastn(repeat_data[[1]], repeat_data[[2]], repeat_data[[3]], opt$outdir, opt$threads, opt$min_cov)
comp_blast <- blastOut[[1]]
comp_seq <- blastOut[[2]]
filtered_gff <- repeat_data[[1]]
base::remove(blastOut)
base::remove(repeat_data)
gc()
uniq_blast <- comp_blast %>%
  dplyr::select(qseqid, sseqid) %>%
  base::unique()

# write sequences and best hit to file
readr::write_tsv(comp_blast, paste0(opt$outdir, "/blast_best_hits.tsv"))
Biostrings::writeXStringSet(comp_seq, paste0(opt$outdir, "/genomic_te_sequences.fasta") )

# ID hits in both directions
comp_blast$strand <- ifelse(comp_blast$sstart < comp_blast$send, "+", "-")
fwd <- comp_blast[comp_blast$strand == "+",]
both_dir <- base::unique(fwd[fwd$qseqid %in% comp_blast[comp_blast$strand == "-",]$qseqid,]$qseqid)
base::remove(fwd)

# Add sequence data to tibble
to_align_tbl <- inner_join(uniq_blast, tibble(qseqid = names(comp_seq), qseq = as.character(comp_seq)), by = "qseqid") %>%
  inner_join(tibble(sseqid = names(comp_seq), sseq = as.character(comp_seq)), by = "sseqid") %>%
  as.data.frame()
to_align_v <- base::unique(to_align_tbl$qseqid)

# Distance calc function
message('Calculating genetic distances')
alnDist <- function(seqid){
  
  to_assess <- comp_blast[comp_blast$qseqid == seqid,]
  # Determine if multistranded
  if (seqid %in% both_dir) {
    # Determine strands
    to_assess <- to_assess %>%
      dplyr::mutate(seqnames = qseqid,
                    start = qstart,
                    end = qend)
    # calculate coverage in direction
    fwd <- to_assess[to_assess$strand == "+",] %>% 
      plyranges::as_granges() %>%
      plyranges::reduce_ranges_directed()
    rev <- to_assess[to_assess$strand == "-",] %>% 
      plyranges::as_granges() %>%
      plyranges::reduce_ranges_directed()
    # select longest
    if(sum(width(fwd)) > sum(width(rev))){
      longest_hit <- GenomicRanges::GRanges(seqnames = fwd$qseqid[1],
                                            ranges = IRanges::IRanges(min(fwd$qstart), max(fwd$qend)))
    } else {
      longest_hit <- GenomicRanges::GRanges(seqnames = re$qseqid[1],
                                            ranges = IRanges::IRanges(min(rev$qstart), max(rev$qend)))
    }
    
  } else {
    # If in one piece seqlect blast hit region
    longest_hit <- GenomicRanges::GRanges(seqnames = to_assess$qseqid[1],
                                          ranges = IRanges::IRanges(min(to_assess$qstart), max(to_assess$qend)))
  }
  
  # Combine longest extracxt and best hit from repeat data
  paired_seq <- c(BSgenome::getSeq(comp_seq, longest_hit),
                  repeat_data[[2]][names(repeat_data[[2]]) == to_assess$sseqid[1]])  
  # get unaligned sequences and convert to ape
  unaln_ape <- ape::as.DNAbin(paired_seq)
  # align with mafft
  aln_ape <- ips::mafft(unaln_ape, method = "localpair", thread = 1, exec = mafft_path, options = "--adjustdirection")
  # calculate distances
  dist_df <- as.data.frame(longest_hit)
  dist_df$kdist <- ape::dist.dna(aln_ape, model="K80", pairwise.deletion = TRUE)[1]
  dist_df$jcdist <- ape::dist.dna(aln_ape, model="JC69", pairwise.deletion = TRUE)[1]
  dist_df$rawdist <- ape::dist.dna(aln_ape, model="raw", pairwise.deletion = TRUE)[1]
  return(dist_df)
  
}

# Run distance function
kdist_list <- mclapply(to_align_v, alnDist, mc.cores = opt$threads)
message("Compiling kdist data")

# Compile distance info
message("Compiling list into dataframe")
kdist_tbl <- purrr::list_rbind(kdist_list) 

# Removing distance list and garbage collection
base::remove(kdist_list)
gc()

message("Converting dataframe to tibble")
kdist_tbl <- dplyr::as_tibble(kdist_tbl)

message("Renaming columns")
kdist_tbl <- kdist_tbl %>%
  dplyr::rename(names = seqnames,
                inner_start = start,
                inner_end = end)

message("Rounding distances")
kdist_tbl <- kdist_tbl %>%
  dplyr::mutate(names = as.character(names),
                kdist = base::round(kdist, 4),
                jcdist = base::round(jcdist, 4),
                rawdist = base::round(rawdist, 4))

message("Wrirint distance tsv to file")
readr::write_tsv(kdist_tbl, paste0(opt$outdir, "/final_kdist.tsv"))

# Join with gff and write to file
message("Joining gff and removing unnecessary columns")
filtered_gff <- dplyr::inner_join(filtered_gff, kdist_tbl, by = "names") %>%
  dplyr::select(-names, -con_width, -subclass, -superfamily, -strand.y, -width.x, -width.y)

message("Renaming strand")
filtered_gff <- filtered_gff %>% dplyr::rename(strand = strand.x)

message("Adjusting coordinates and removing inner coordinates")
filtered_gff <- filtered_gff %>% 
  dplyr::mutate(eg_start = start, eg_end = end,
                start = start + inner_start - 1,
                end = end - width.x + inner_end,
                source = "BetterTiming", score = ".") %>% # Add source and score
  dplyr::select(-inner_start, -inner_end)

message("Cobverteding to Granges object")
filtered_gff <- filtered_gff %>% 
  plyranges::as_granges()

message("Writing gff to file")
plyranges::write_gff3(filtered_gff, paste0(opt$outdir, "/kdist_", sub(".*\\/", "", opt$annotation_gff)))
