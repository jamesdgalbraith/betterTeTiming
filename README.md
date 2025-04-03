# Better TE Timing

A Rscript aiming to better approximate the recency of TE insertions by calculating genetic distance from the most similar TE

Current dependancies
- MAFFT
- BLAST+
- R (with optparse purrr biostrings bsgenome ape ips tidyverse plyranges installed)


For easy set up with mamba: `mamba create -n betterTiming conda-forge::r-optparse conda-forge::r-purrr bioconda::bioconductor-biostrings bioconda::bioconductor-bsgenome conda-forge::r-ape conda-forge::r-ips conda-forge::r-tidyverse bioconda::bioconductor-plyranges bioconda::mafft bioconda::blast`


To run from terminal:
```
conda activate betterTiming

Rscript nearestDistance.R \
    -g </path/to/genome_assembly.fasta>  \
    -l </path/to/earlgrey_library.fasta> \
    -a </path/to/earlgrey_annotation.gff> \
    -o </path/to/out/directory> \
    -t <number of threads>
```


The gff should have the subclass and superfamily of TEs formatted in the RepeatMasker/Repbase style in the type column (e.g. `LINE/CR1`), and the family/subfamily of the TE as described in the TE library's FASTA headers as ID in the attribute column (e.g. `ID=rnd_1-family_1`)
