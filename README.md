![alt text](https://github.com/gchchung/TeloSearchLR/blob/main/logo_copy_v3.svg)
TeloSearchLR (**telo**mere **search** using **l**ong sequencing **r**eads) is a Python script for aiding the identification of telomeric repeat motifs.

## Contents
[Installation](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#installation)

[Usage](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#usage)

1. [Telomeric repeat motif discovery using reads from NCBI SRA](https://github.com/gchchung/TeloSearchLR?tab=readme-ov-file#1-telomeric-repeat-motif-discovery-using-reads-from-ncbi-sequence-read-archive)
2. [Telomeric repeat motif discovery in your own library](https://github.com/gchchung/TeloSearchLR?tab=readme-ov-file#2-telomeric-repeat-motif-discovery-in-your-own-library)
3. [Unusually long telomeric repeat motifs](https://github.com/gchchung/TeloSearchLR?tab=readme-ov-file#3-unusually-long-telomeric-repeat-motifs)
4. [Telomeric repeat motif discovery, but sorting the motifs by repeat period first](https://github.com/gchchung/TeloSearchLR?tab=readme-ov-file#4-telomeric-repeat-motif-discovery-but-sorting-the-motifs-by-repeat-period-first)
5. [Testing to see if a specific motif is repeated at read ends](https://github.com/gchchung/TeloSearchLR?tab=readme-ov-file#5-testing-to-see-if-a-specfic-motif-is-repeated-at-read-ends)

[Commands and options](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#commands-and-options)

[How it works](https://github.com/gchchung/TeloSearchLR/blob/main/README.md#how-it-works)

- [Tested cases](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#tested-cases-2024-11)

[Contributing](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#contributing)

[License](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#license)

[Please cite](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#please-cite)

[Contact](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#contact)

## Installation

### Conda
Install in a new environment named "telosearchlr-env".
```bash
conda create -n telosearchlr-env -c bioconda -y telosearchlr    # Create a new conda environment with TeloSearchLR
conda activate telosearchlr-env                                 # Activate this environment
conda install -c conda-forge -y python-kaleido                  # In the next release we will fix this: we didn't include "python-kaleido"
TeloSearchLR.py -h                                              # Test the installation by calling the help function
conda deactivate                                                # De-activate the conda environment before using other environments
```

### Docker
```bash
docker pull gchchung/telosearchlr:telosearchlr_v1.0.0   # Pull the TeloSearchLR image from Docker repository ```gchchung/telosearchlr```.
docker run telosearchlr_v1.0.0 TeloSearchLR.py -h       # Test by asking for the help message.
```

### From source
```bash
git clone https://github.com/gchchung/TeloSearchLR.git  # clone repo
cd TeloSearchLR

# after you have installed dependencies, you can run
python TeloSearch.py -h
```


## Usage
### 1. Telomeric repeat motif discovery using reads from NCBI Sequence Read Archive
|required       | description                                                         |
|---------------|---------------------------------------------------------------------|
|-f             | FASTA file of the reads (STR)                                       |
|-k             | shortest period to search for (INT)                                 |
|-K             | longest period to search for (<=500 nt if -t is unspecified) (INT)  |
|-m             | numerical rank of the most frequent motif to plot (INT)             |
|-M             | numerical rank of the least frequent motif to plot (INT)            |
|-n             | number of nucleotides to plot the repeat occupancy (INT)            |

Download a PacBio genomic sequencing library from *Caenorhabditis elegans* by [Yoshimura & al (2019)](https://pubmed.ncbi.nlm.nih.gov/31123080/) using ```fasterq-dump``` from [sra-tools](https://github.com/ncbi/sra-tools). Then, get TeloSearchLR to find and rank repeat motifs between 4 and 20-nt long (-k 4 -K 20) from the terminal 1000 nts (-t 1000) of reads 2000 nts or over. Plot the occupancy patterns of the top 100 motifs (-m 1 -M 100) in the first and last 6000 nt of all reads (-n 6000) longer than 12,000 nts.
```bash
# get the reads
fasterq-dump --fasta SRR7594465

# run TeloSearchLR
python TeloSearchLR.py -f SRR7594465.fasta -k 4 -K 20 -t 1000 -m 1 -M 100 -n 6000
```
The [output of this](https://github.com/gchchung/TeloSearchLR/blob/main/repeatPattern.m1.M100.png) can be found in the *.results folder. As expected, the known telomeric repeat of *C. elegans*, TTAGGC, is almost exclusively found at the 3' ends of reads and as the reverse complement at the 5' ends of reads. It is also the second most frequent terminal repeat motif (second row of plots). Typically, telomeric repeat motifs are among the most frequent repeat motifs, but may not be when the telomeres are short.   

### 2. Telomeric repeat motif discovery in your own library
The required flags are the same as above, but you must convert your reads into a FASTA file, usually from a FASTQ. A possible method is by using the Unix sed.
```bash
# convert YOUR_LIBRARY.fastq to YOUR_LIBRARY.fasta 
sed -n '1~4s/^@/>/p;2~4p' YOUR_LIBRARY.fastq > YOUR_LIBRARY.fasta

# run TeloSearchLR
python3 TeloSearchLR.py -f YOUR_LIBRARY.fasta -k 4 -K 20 -m 1 -M 100 -n 6000
```

### 3. Unusually long telomeric repeat motifs
Typically, telomeres maintained by telomerase have short motifs (<=30 bps), but telomeres maintained by ALT can have very long motifs. To find repeat motifs longer than 500 bps, change the -t paramenter to 2 * [the maximum period you'd like to consider]. This is because when *t* = 1000, the longest tandem repeat that can fit in 1000 nt is a 500-nt repeat (2*500 = 1000).
|option    | description                                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------------------------|
|-t            | the terminal region (in bps) to rank repeat motifs. It must be at least 2 times the *K*-value (INT, default 1000)       |

In our manuscript, we described how telomeric motifs were found in two *Strongyloides stercoralis* libraries. *Strongyloides* nematode species are thought to have evolutionarily lost telomerase ([Mota & al 2024](https://pubmed.ncbi.nlm.nih.gov/38316773/), Fig. 1), and the low chromosome count of *S. stercoralis* ([Kounosu & al 2023](https://pubmed.ncbi.nlm.nih.gov/38008120/)) makes this species ideal to search for novel telomeres.

```bash
# download the S. stercoralis libraries
fasterq-dump --fasta SRR25177361
fasterq-dump --fasta SRR25177362

# concatenate the two libraries
cat SRR25177361.fasta SRR25177362.fasta > Sstercoralis.fasta

# run TeloSearchLR, search for 21 to 1000-mers of repeat motifs and graph the terminal 8000 nts
# because we are testing up to 1000-mers, t-value has to be at least 2000
python TeloSearchLR.py -f Sstercoralis.fasta -k 21 -K 1000 -m 1 -M 100 -t 2000 -n 8000
```


### 4. Telomeric repeat motif discovery, but sorting the motifs by repeat period first
With the -e flag, each TeloSearchLR run can also sort the discovered motifs first by repeat period, then by occupancy. This mode may be useful if you already suspect the telomeric motif to be between *k* and *K* nucleotides long. Required options are the same as above plus the -e.
|required       | description                               |
|---------------|-------------------------------------------|
|-e             | exhaustive mode: sort motifs by period first then by occupancy. |

Using the same *C. elegans* library as in example 1, we run TeloSearchLR in "exhaustive" mode. 
```bash
# run TeloSearchLR in exhaustive mode - sorting the motifs first by period then by occupancy
python TeloSearchLR.py -f SRR7594465.fasta -k 4 -K 20 -t 1000 -m 1 -M 100 -n 6000 -e
```

### 5. Testing to see if one specfic motif is repeated at read ends
Suppose you sequence a novel organism (library: NOVEL_ORG.fasta) whose closest known relative has a telomeric motif of TTAGGG. You would like to quickly check if this novel organism also has TTAGGG telomeres. TeloSearchLR offers a "single-motif" mode that allows you to check only one motif to see if it is telomeric or not. 

|required       | description                                                         |
|---------------|---------------------------------------------------------------------|
|-f             | FASTA file of the reads (STR)                                       |
|-n             | number of nucleotides to plot the repeat occupancy (INT)            |
|-s             | the repeat motif to test, eg. ```-s TTAGGG``` (STR)                 |
|-T             | a TideHunter search result file in tabular format, eg. a *TideHunterTable.txt file from an earlier TeloSearchLR run (STR)  |

We require a [TideHunter search output in tabular format](https://github.com/yangao07/TideHunter?tab=readme-ov-file#to-generate-consensus-sequences-in-tabular-format) in this mode. This can come from an independent run of TideHunter, or an earlier run of TeloSearchLR, as long as motif period was included in the TideHunter search. Here we demonstrate how to check if a 6-mer TTAGGG motif is a telomeric motif by setting -p, -P, and -m all to 6 in TideHunter.
```bash
# run TideHunter. TideHunter search result in NOVEL_ORG_6mer_repeats.TideHunterTable.txt
# the -p 6 -P 6 -m 6 means TideHunter will only search for and display repeats of period 6 (nt)
TideHunter -f 2 -p 6 -P 6 -m 6 NOVEL_ORG.fasta > NOVEL_ORG_6mer_repeats.TideHunterTable.txt

# run TeloSearchLR
python TeloSearchLR.py -f NOVEL_ORG.fasta -n 4000 -s TTAGGG -T NOVEL_ORG_6mer_repeats.TideHunterTable.txt
```

## Commands and options
```text
Usage:   python TeloSearchLR.py [options]
                
Options:
  Run modes:
    (default)                      "occupancy mode", repeat motifs ranked by occupancy
    -e --exhaustive                enable "exhaustive mode", motifs ranked by period AND occupancy
    -s --single_motif       STR    enable "single-motif mode", specify the motif whose occupancy is to be plotted
  Required for all modes:
    -f --fasta              STR    long-read sequencing library file name
    -n --num_of_nucleotides INT    number of nucleotides to plot motif occupancy
  Required for occupancy and exhaustive modes:
    -k --k_value            INT    shortest repeat period (>0) to consider by TideHunter
    -K --K_value            INT    longest repeat period (≥k) to consider by TideHunter
    -m --mth_pattern        INT    most frequent motif to plot (>0)
    -M --Mth_pattern        INT    least frequent motif to plot (≥m)
  Required for single-motif mode:
    -T --TideHunter         STR    a TideHunter (≥v1.5.4) tabular output
  Other options:
    -t --terminal           INT    terminal number of nucleotides to consider for ranking motif occupancy [1000]
    -c --cores              INT    number of threads to use for TideHunter [4]
    -p --path               STR    path of TideHunter (if not already in $PATH)
    -v --version                   display the version number and quit
    -h --help                      display this help message and quit
```

## How it works
We take advantage of the fact that telomeric sequences only appear at the ends of long genomic reads prepared using standard procedures. Thus, we can focus our search for telomeric motifs on the terminal regions (*t* bp) of sequencing reads. Intact telomeres will always be read in its canonical (G-rich strand) orientation at the 3' ends of reads, and as the reverse complement on the 5' ends of reads.
![image](https://github.com/user-attachments/assets/d0043e21-668d-4d07-89fe-ab8704fdc519)

TeloSearch uses [TideHunter](https://github.com/yangao07/TideHunter) to search for repeats in the first and last *t* bp of all reads ≥ 2*t* bp. You need to specify the shortest (-k) and the longest (-K) repeat periods to search.

![image](https://github.com/user-attachments/assets/fb1650a3-1274-4414-b881-cde8ca0aaf75)

After this initial search, TeloSearchLR ranks the discovered motifs by their occupancy: the number of bases they occupy in the first and last *t* bps. The motifs are counted together if they are reverse complements (eg. ATTG, CAAT), or if they are cyclical permutations (eg. ATGC, TGCA, GCAT, CATG).

![image](https://github.com/user-attachments/assets/616b915b-34ca-4598-8e4f-2f5baa38e6be)

TeloSearchLR calls TideHunter again to search for repeat motifs across entire reads. You specify which motif rankings in the first and last *t* bps TeloSearchLR should graph (-m and -M).

![image](https://github.com/user-attachments/assets/d18b9ac7-facc-4146-9f1e-c0b1ed7a2f10)

If entire telomeres were captured on long sequencing reads, you should be able to visually identify candidate telomeric motifs by their **terminal stranded occupancy**.

![image](https://github.com/user-attachments/assets/587d06fd-0001-4e9e-b484-a1e16e6aab73)

Occupancies are represented in three plots: the occupancy by position at the 5' *n* nucleotides, the occupancy by position at the 3' *n* nucleotides, and a bar plots summarizing the occupancies at the terminal *n* nucleotides and the middle of the reads.

![image](https://github.com/user-attachments/assets/9a6fcf8c-6899-446f-8217-903f89d16c2f)

We have seen instances with very long telomeres that the true telomeric motif does not have an obvious terminal stranded occupancy ([Chung & al 2024](https://doi.org/10.1101/2024.10.29.617943), Fig. 4). This is because if the telomeres are *longer* than the reads, the telomeric motifs are no longer only found on one or the other end. **Thus, this algorithm works best if sequencing reads are much longer than the average telomere in your sample (and telomere lengths are species- and strain-specific).**

If you can identify the terminal stranded occupancy, you can make direct inferences on the library's telomere length. This is because the shape of the occupancy curve reflects the distribution of telomeres captured by the long sequencing reads. Thus, this may be used as a quick way to estimate telomere lengths without specialized protocols (such as in [Karimian & al 2024](https://pubmed.ncbi.nlm.nih.gov/38603523/) or in [Schmidt & al 2024](https://pubmed.ncbi.nlm.nih.gov/38890299/)).

![image](https://github.com/user-attachments/assets/fb2a3dad-8c2d-4ac4-b328-b54561356b7b)

For all the repeat motifs plotted, a corresponding text file with the plotted values is provided so you may use your favourite software (eg. R ggplot2, Python matplotlib) to do the plotting. 

### Tested cases (2024-11)
We have tested TeloSearchLR on libraries from metazoans, plants and fungi. In several beetle and wasp species, we successfully identified new telomeric motifs where short-read searches and other methods did not work (motifs in red below).

![image](https://github.com/user-attachments/assets/09e5508f-9e1c-47a5-9f77-23bb34c48727)
(from Table 2, [Chung & al 2024](https://doi.org/10.1101/2024.10.29.617943))


## To do

- Improving the SVG output. PNG output may be deprecated in a future release as SVG is immediately publication-ready
- Possibly move to matplotlib for graphing
- A fork for short paired-end reads?

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.


## License

[NYU License](https://github.com/gchchung/TeloSearchLR/blob/main/LICENSE)

## Please cite
Chung G, Piano F, Gunsalus KC. 2025. TeloSearchLR: an algorithm to detect novel telomere repeat motifs using long sequencing reads. *G3 (Bethesda)* jkaf062.
 [link](https://academic.oup.com/g3journal/article/15/6/jkaf062/8102967).

## Contact
George Chung (gc95@nyu.edu), Research Scientist, New York University
