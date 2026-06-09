# Better TE Timing

## Quick start/how to run
A project aiming to better approximate the recency of TE insertions by calculating genetic distance from the most similar TE

Current dependancies
- Python
    - Biopython
    - pyfaidx
- BLAST+ or minimap2

For easy set up with mamba/conda:
`mamba create -n betterTiming conda-forge::biopython bioconda::pyfaidx bioconda::minimap2 bioconda::blast`

To run the minimap2 version from terminal:
```
conda activate betterTiming
python /home/jgalbrai/ross/scales/comparative/GenomeDynamicsOfScaleInsects/scripts/better_TE_timing_minimap2.py \
        -t <number of threads> \
        -i </path/to/earlgrey_annotation.gff> \
        -o </path/to/betterTEtiming_annotation.gff> \
        -g </path/to/genome_assembly.fasta>
```

There are more options allowing you to tweak with minimap2 and BLASTN options, more can be added in furture if desired.

I've tested this successfully on a macOS Tahoe 26.5 (M2 processor) and on a Ubuntu 22.04.5 LTS server (Intel(R) Xeon(R) Gold 5416S processor).

The gff should have the subclass and superfamily of TEs formatted in the RepeatMasker/Repbase style in the type column (e.g. `LINE/CR1`), and the family/subfamily of the TE as described in the TE library's FASTA headers as NAME in the attribute column (e.g. `NAME=rnd_1-family_1`)

Try to make sure you don't have any nested TEs: if a nested repeat and a parent repeat belong to the same superfamily, they'll align to each other (and hence have genetic distance of 0).

## Rationale/reasoning behind the idea
At present, RepeatMasker and EarlGrey calculate the Kimura distance between a repeat and the consensus repeat they identified it with using BLAST, treating the consensus of a hypothetical ancestor. This hypothetical ancestor is the consensus of numerous repeats from across the genome which exhibit enough similarity to be grouped together by RECON or RepeatScout. These clusters could consist of a repeat expansion at a single time point, in which case the current approach is appropriate. Assuming they are not under selection and are evolving neutrally, all the repeats RepeatMasker/EarlGrey finds will have accumulated SNPs at approximately the same rate, and hence old expansions will appear old and recent expansions recent. However, after many years of manual curation of multiple sequences alignments, I've found these groups are more likley to consist of multiple repeat expansions with clear "subfamilies" within the larger family or of ongoing expansions, with repeats with incrementaly decreasing similarity from two or three very similar repeats (likely the newest insertions).

As such, while the RepeatModeler's approach has long used this as an approximation for the age/timing of the insertions I feel this is highly flawed. A newly inserted repeat will have 0 mutations compared to its "parent", yet the genetic distance reported by RepeatModeler is unlikely to be zero, as the consensus repeat is formed from numerous repeats. If manually curated using the commonly used suggestions of [Wicker et al. 2007](https://doi.org/10.1038/nrg2165), these families of repeats will have up to 80% divergence, a pairwise distance of 0.2! The approach Better TE timing (BTT) instead calculates the genetic distance between a repeat insertion anf the most similar repeat insertion from the same repeat superfamily. This should provide a better estimate of how recently the repeat was inserted, however I recognise there are a lot of assumptions this method relies on and (like every model) there are flaws. The most similar repeat is identified using either BLASTN or minimap2, and the genetic disatance between the two sequences calculated using the number of matches vs mismatches reported by BLASTN/minimap2. At present Jukes-Cantor distance and the default genetic distance metric used by the tool are reported. 

The figure below hopefully explains this rationale graphically:
![Image explaining RepeatModeler's approach to calculating genetic distance compared to BTT's](https://github.com/jamesdgalbraith/betterTeTiming/blob/main/BTT_Explained.png)

p.s. If someone can provide a decent rationale as to why Kimura 2-parameter distance should be used instead other than "we've done it forever" please let me know in the issues and I'll try and incorporate it as an option is a future release. Until then I'll use the "simplest is best" argument and stick to Jukes-Cantor.

p.p.s. `better_TE_timing_blastn_repeat_library.py` replicates the RepeatMasker approach using BLASTN for testing. To run this you'll additionally need the repeat library used to annotate the genome.