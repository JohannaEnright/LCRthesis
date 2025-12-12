This file contains all the code which is referenced in my final honour's thesis project. Much of this project was published in the journal of Molecular Biology and Evolution and can be found at the following address: https://doi.org/10.1093/molbev/msad084 (Enright, Dickson, & Golding, 2023). 
This project contains additional code from that project used to gather results for further exploration and analysis of low complexity regions (LCRs),
which are compositionally biased, highly repetitive regions in proteins which can impart critical cellular functions, or can be implicated in diseases, depending on the protein and its cellular role.
The project looks at LCRs in the entire proteome and genome of five model organisms: Homo sapiens, Saccharomyces cerevisiae, Drosophila melanogaster, Caenorhabditis elegans, and Arabidopsis thaliana,
and compares the entropy of protein LCRs to their encoding DNA sequences and vice versa, to help gain evolutionary insights, and better understand the role of these repetitive sequences.

The main script components include:

1. calculate_entropy_biological_lcrs.py
   This script was used to determine the entropies of LCRs from real biological sequences and requires the following inputs:
      i) a genbank file for a particular organism, which can be multiple genbank files for different chromosomes, concatenated together
      ii) the direction of comparison, which can be "pro_dna", "pro_cod", "dna_pro", "dna_cod", "cod_dna", or "cod_pro".
          For example, pro_dna scans for LCRs at the protein level and compares the sequence entropy to that of its underlying DNA sequence,
          while dna_cod scans for LCRs at the DNA level and compares the entropy to that of the codon entropy it encodes.
      iii) segA parameter W ,which is the window length as an integer
      iv) segA parameter K1, which is the "trigger complexity" as a float.
          Trigger complexity is the maximum entropy threshold within the sliding window that triggers segA to stop and expand the window in search of the full LCR.
      v) segA parameter K2, which is the "extension complexity" as a float.
         Upon finding an LCR of the given window length below the trigger complexity threshold, segA will scan outwards from this window until it determines the full LCR with an entropy below the extension complexity threshold.

   This script outputs a plain text files parsed on whitespace, name as inputfilename_Output.
   It contains the information for protein ID, genbank accession number, the DNA strand the LCR was found on, if there were introns found in the gene, if the LCR spanned an exon juncion, the start position of the gene, the end position of the gene, the starting amino acid position of the protin LCR, the end amino acid position of the protein LCR, the entropy of the LCR from the sequence type you were searching for, the starting LCR position on the chromosome, the end LCR position on the chromosome, and the entropy of the corresponding sequence type which was being compared. The key information which was used for further analysis was entropy y ~ entropy x.
   
   Other plain text files are also outputed and include more information on the LCRs, information on proteins which did not contain LCRs (for counting/checking purposes), and information if the LCRs/genes were more atypical. For example, if the LCR was found in a gene which was located on both strands of the DNA, or if the LCR spanned an exon junction.
   
2. GiveCombSysCallArgs.sh
   This is a bash file which is used to feed CombinedUsingSysCall_Codon3.py a combination of set parameters for segA. This was important to ensure that we were observing the same trends given slightly different LCR identifying parameters, considering the somewhat arbitrary nature of these parameters. In total, it will feed CombinedUsingSysCall_Codon3.py 27 unique parameter combinations as shown in the supplementary section of the publication (Enright, Dickson, & Golding, 2023).
   
3. average_protein_length.py
   This script requires a genbank file as input, and will return the average protein length for that file.
   Again, it also accepts multiple concatenated genbank files, hence it was used to calculate the average protein length for each organism.
   This information was used in the generaterandomproteome script.

4. change_aa_props_degeneracy.py
    Some amino acids can be encoded by multiple codons. For example, serine can be encoded by six unique codons and thus has a degeneracy of six,
    while methionine for example can only be encoded by AUG.
    This script allows you to alter the proportions of amino acids of specified degeneracies. 
    This script takes the following inputs:
       i) filename of a file outputed by nt_cod_aa_bias.py which contains whole proteome amino acid proportions.
       ii) a string of integers separated by a space corresponding to amino acid degeneracies whose proportions you wish to alter.
          For example "2" or "2 3 4".
       iii) the percentage you wish to alter the selected amino acid proportions by. For example 120 will increase the proportions by 20%.
     This output file of altered proportions can then be used as input for generate_biologically_biased_proteome.py to simulate protein and DNA sequences based on these proportions.
     This allows for the experimentation of how preferences for amino acids of particular degeneracies can influence the correlation of LCR entropies.

 5. change_codon_props.py
     This script works similarly to change_top_props.py in that it changes certain codon proportions by a specified amount so that you can later simulate sequences based on these biased proportions.
     The input requirements are:
        i) a filename of letter proportions outputed by the nt_cod_aa_bias.py script.
        ii) the codon type whose proportion you wish to alter - "1", "2", or "3", corresponding to mononucleic, dinucleic, or trinucleic codons.
        iii) n - the number of codons, from the most commmon to least common, that you wish to alter
        iv) the percent you wish to alter them by. For example '50' would half the selected codon proportions
     This output file can then be used as input for generate_biologically_biased_proteome.py to experiment with how codon biases of specific codon types
     in an organism can influence the correlation of LCR entropies.


 6. change_top_props.py
     This script is meant to take the most common letters in an alphabet (either nucleotide, codon, or amino acid) and alter them by a given percentage.
     The inputs are: start integer (exp. 1), end integer (exp. 5), and percent (exp. 200). This will increase the top 5 most common letters by 200%.
     And will write the changed proportions back to a file. This output file can then be used on generate_biologically_biased_proteome.py to experiment with how letter biases in an organism can influence the number and nature of LCRs in a proteome.

 7. count_pattern_types.sh
     This is a bash script which takes the name of an output file from periodicity_and_entropy.py.
     It counts the number of times were mono-, di-, and tri-, (up to nine) repeats etc. occuring within an LCR in a proteome/genome.
     This gives insight into if certain repeat types are more common or less common in an organism.
     
8. codon_class_categorize_lcrs.py
     This script requires the _Seq file outputed by CombinedUsingSysCall_Codon3.py which contains the sequences of all of the detected LCRs for a given organism, and a file output location.
     It categorizes each LCR into one of three categories:
       i) containing three unique nucleotides
       ii) containing two unique nucleotides
       iii) containing one unique nucleotide and outputs this information as counts of LCRs in each category.
     This information is used for the calculation of the preference coefficient, see Enright, Dickson, & Golding (2023).

 9. count_codon_types_all_cds.py
     This script takes a genbank (or multiple) file as input. It outputs a text file with information for the total codons of each codon class summed for all sequences in a genome.
     This information is used for the calculation of the preference coefficient, see Enright, Dickson, & Golding (2023).

 10. generate_biologically_biased_proteome.py
     This script creates a genbank file with 100 000 protein sequences and encoding DNA sequences to see if LCRs can be created simply from either an organisms bias
     towards/against certain nucleotides, codons, and/or amino acids. It takes a file as outputed from nt_pro_cod_genome_bias.py which contains either the nucleotide, codon, or amino acid
     proportions and a species name (just for naming purposes) and simulates a proteome. This can then be used to count the number of LCRs that occur and compare their nature to the actual biollogical LCRs.
     
 11. generate_random_proteome.py
     This script is to simulate 100 000 sequences outputed in fasta file format. 
     The input values are the simulation type ('random', 'biased' or 'mutate'), species ('sc', 'dm', 'hs', 'ce', 'at'), fasta file of sequences to mutate (only if you choose the 'mutate' option),
     write file (name of file to write to).
     The species input is necessary to generate sequences which have the average length, the correct codon class proportions, and the same proportion of LCRs as the corresponding species.
     The 'random' simulation will randomly generate sequences selecting from codons with weightings equal to the biological proportions, and is designed to show that LCRs are not a result of random chance.
     The 'biased' simulation will generate sequences selecting from codons with weightings equal to the biological proportions, where the codon that was previously selected has a higher weighting which was chosen so to create the same proportion of LCRs as found in the respective species.
     The 'mutate' simulation will add 1000 mutations to the sequences generated from the 'biased' simulation.
     Only the 'random' simulation was used in the final publication, the other simulations were later improved upon in the improved_codon_class_lcr_simulation.py and the inc_slip_sim_species_specific.py scripts.
     
 12. improved_codon_class_lcr_simulation.py
    It will simulate 100 000 LCR sequences using an increasing DNA polyermase slippage model that also takes into account the codonclass proportions (class 1,2 or 3).
    Codon class proportion biases are based off of precalculated preference coefficients for each species, see Enright, Dickson, & Golding (2023).
    This model is reffered to as the 'Slip+CC' model.
    The input requirements are a file of secies specific codon bias proportions Outputed by the nt_cod_aa_bias.py script and species type ('hs', 'sc', 'dm', 'at', 'ce').
    It generates sequences using the average protein length, codon proportion, and preference coefficient for biological relevance.
     
 13. inc_slip_sim_species_specific.py
     This simulation generates 100 000 sequences protein and corresponding DNA sequences, outputed as a genbank file. It creates LCRs simulating the DNA polyermase slippage mechanism,
     where the probability of resampling the same codon increases with increasing number of repeats of that codon.
     Sequences are generated such that they include the same proportion of LCRs, the same protein length, and same codon bias as are biologically found in the organism.
     The script requires the inputs:
      i) simulation type ('random', 'biased' or 'mutate') where 'random' randomly generates sequences not using an increasing slippage model,
         'biased' generates sequences using an increasing slippage model (referred to as the 'Slip model' in the publication),
          and 'mutate' adds 1000 random mutations to the third position of codons from the generated sequences (reffered to as 'Slip + Syn model' in the publication),
          using biologially relevant weights (5:1) for transition/tranversion mutations.
      ii) species ('hs', 'sc', 'dm', 'at', 'ce')
      iii) file of codon proportions as generated from nt_pro_cod_genome_bias.py

 14. make_periodicity_entropy_plots.R
     This is an R script which takes creates scatter plots with a linear regression and calculated correlation coefficient and confidence interval
     from entropy results of periodic repeats of biological organisms.
     The inputs are: i) the direction ("pro_dna", "pro_cod", "dna_pro", "dna_cod", "cod_pro", "cod_dna") as "x_y"
                     ii) the entropy outputs from periodicity_and_entropy.py as a filename.
     
 15. making_scatter_correlation_plots.R
     This is an R script which takes creates scatter plots with a linear regression and calculated correlation coefficient and confidence interval from entropy results of biological organisms.
     The inputs are: i) the direction ("pro_dna", "pro_cod", "dna_pro", "dna_cod", "cod_pro", "cod_dna") as "x_y"
                     ii) the entropy outputs from calculate_entropy_biological_lcrs.py as a filename.
     
 16. nt_pro_cod_genome_bias.py
     This script takes a genbank file (or multiple concatenated) as input and outputs 4 files: each containing the nucleotide frequency,
     the nucleotide in coding sequence frequency, the codon frequency, and the amino acid frequency.
     The codon frequency output was used for seqeunce simulations such as improved_codon_class_lcr_simulation.py and inc_slip_sim_species_specific.py.
     
 17. periodicity_and_entropy.py
     This script requires a _Seq file outputed by CombinedUsingSysCall_Codon3.py containing the sequences of the detected biological LCRs and the sequence type, either 'protein' or 'dna'.
     Optional inputs include:
       i) wholeseq="True" or "False, which outputs either the entropies of only the repeat segments ("False") or the entropies of the entire LCR sequence ("True").
       ii) The lengths that each repeat has to be to be considered a repeat. For example, for mono-amino acid repeats, the releat size needed to be 4 (rep1p="4"),
           and the length of tri-nucleotide repeats needed to be 4 (rep3d="4").
    It outputs 2 files:
       i) one contains the DNA, codon, and protein entropy information for LCR sequences which did not contain periodic repeats.
       ii) the second contains this information for sequences which do contain periodic reapeats, the index positions where the repeat begins and ends,
           as well as the repeat type, for example mono-, di-, or tri- amino acid repeats.
    This information was used to compare the correlation of entropies between LCRs which contained periodic repeats and those which did not.

