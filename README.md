# PPRmatcher
Julia (https://julialang.org) code to match pentatricopeptide repeat (PPR) protein sequences to likely RNA targets.
The software requires the ArgParse.jl and Statistics.jl packages to be installed.<br>

The software uses scoring table to relate the the 5th/last amino acids of each PPR motif with the aligned nucleotide.
The original idea is from Barkan et al. Plos Genetics 8:e1002910. Several scoring tables are available to choose from.<br>
Yan.tsv uses data from Yan et al. Nucleic Acids Res. 47:3728-3738; this is experimental data on binding of synthetic P-type PPR proteins to oligonucleotides in vitro.<br>
Kobayashi.tsv used data from Plant Cell Physiol. 60:862-874; this is from observed frequencies of different 5th/last combinations in alignments of natural PLS-type proteins to their target RNAs.<br>
The Millman .tsv files are similarly derived from observed frequencies of different 5th/last combinations in alignments of natural PLS-type proteins to their target RNAs, but separate scoring tables have been derived for each type of PPR motif.<br>

To run the software, install julia (instructions here: https://julialang.org/downloads/) and the ArgParse and Statistics packages (instructions here: https://datatofish.com/install-package-julia/)<br>

PPRmatcher.jl requires three arguments:
1. One or more PPR sequences in the output format from PPRfinder (https://github.com/ian-small/PPRfinder)
2. One or more target seqeunces in fasta (or multifasta) format
3. One or more scoring tables from the scoring_table directory

usage: PPRmatcher.jl [-r REFERENCE] [-e EDIT_SITE] [-s STRAND]
                     [-o OUTFILE] [-h] PPR(s) target(s)
                     scoring_table(s)...

positional arguments:
  PPR(s)                path to file of PPRfinder motifs to search
                        with
  target(s)             path to file of FASTA format target
                        sequence(s) to search in
  scoring_table(s)      path to scoring table or directory of scoring
                        table files

optional arguments:
  -r, --reference REFERENCE
                        path to reference sequence for calculating
                        score distribution
  -e, --edit_site EDIT_SITE
                        edit site position within the sequence(s)
                        (type: Int64, default: 0)
  -s, --strand STRAND   which target strand(s) to match to (default:
                        "F")
  -o, --outfile OUTFILE
                        path to outfile for saving results
  -h, --help            show this help message and exit


Examples of use cases:

#1 finding the probable target of mRPF2 in the nad6 transcript (the command below assumes that the pwd is PPRmatcher/test)<br>
julia ../src/PPRmatcher.jl mRPF2.motifs.txt nad6.fasta ../scoring_tables/Yan.tsv

#2 finding the best matching editing site for the synthetic PPR dsn3PLS (the command below assumes that the pwd is PPRmatcher/test)<br>
julia ../src/PPRmatcher.jl -e 41 dsn3PLS.motifs.txt cp_editing_sites.fasta ../scoring_tables/Kobayashi.tsv
