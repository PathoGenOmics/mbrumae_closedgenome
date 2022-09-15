# Genomic analysis of _Mycobacterium brumae_ sustains its nonpathogenic and immunogenic phenotype  
__Chantal Renau-Mínguez<sup>1</sup>,__ 
__Paula Herrero-Abadía<sup>2</sup>,__ 
__Vicente Sentandreu<sup>3</sup>,__ 
__Paula Ruiz-Rodriguez<sup>1</sup>,__ 
__Eduard Torrents<sup>4,5</sup>,__ 
__Álvaro Chiner-Oms<sup>6</sup>,__ 
__Manuela Torres-Puente<sup>6</sup>,__ 
__Iñaki Comas<sup>6</sup>,__ 
__Esther Julián<sup>2*</sup>__
__and Mireia Coscolla<sup>1*</sup>__
<br>
<sub>

1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain  

2. Genetics and Microbiology Department, Faculty of Biosciences, Autonomous University of Barcelona, 08193, Bellaterra, Barcelona, Spain 

3. Genomics Unit, Central Service for Experimental Research (SCSIE), University of Valencia, Spain  

4. Bacterial Infections and Antimicrobial Therapies Group, Institute for Bioengineering of Catalonia (IBEC), Baldiri Reixac 15-21, 08028 Barcelona, Spain  

5. Microbiology Section, Department of Genetics, Microbiology and Statistics, Biology Faculty, Universitat de Barcelona, 08028 Barcelona, Spain  

6. Instituto de Biomedicina de Valencia (IBV), CSIC, 46010, Valencia, Spain  </sub>

<sub> * Correspondence:  <sub>

<sub> mireia.coscolla@uv.es (Mireia Coscolla); Esther.Julian@uab.cat (Esther Julián) <sub>

## Aim of  this repository
The main purpose of this repository is to display the scripts made for that release, in order to achieve reproducibility.

## Scripts
This folder contains multiple subfolders with scripts made for certain purpose.
### blast_analysis  
Scripts for the analysis to get protein identity of *Mycobacterium tuberculosis* H37Rv in the analyzed genomes of interest. 
### duplicated_genes  
Scripts for the analysis to get genes with less than 300bp repeated in order to exclude this genes in Illumina genomic analysis. 
- clean_genes.py
- multifasta.py
- process_mummer_output.py
- command_mummer.sh
- duplicated_genescoord.tsv

## Closed genome
This folder contains multiple files related to the analyzed closed genome: *Mycobacterium brumae* ATCC 51384<sup>T</sup>  

__brumae.fasta:__ fasta file with the genomic sequence of *Mycobacterium brumae* ATCC 51384<sup>T</sup>  

__brumae.gff:__ annotation file used in this study.  

__brumae.sqn:__ file for submission to NCBI.  