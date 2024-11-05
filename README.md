![alt text](https://github.com/gchchung/TeloSearchLR/blob/main/logo_copy_v3.svg)
TeloSearchLR (**telo**mere **search** using **l**ong sequencing **r**eads) is a Python script for aiding the identification of telomeric repeat motifs.

## Contents
[Installation](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#installation)

[Usage](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#usage)

[Commands and options](https://github.com/gchchung/TeloSearchLR/tree/main?tab=readme-ov-file#commands-and-options)

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
conda install -c conda-forge -y python-kaleido                  # In the next release we will fix this, we did not include "python-kaleido"
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
### How it works
![Figures_01 (2024-10-11) panel A medium](https://github.com/user-attachments/assets/579c422d-8b8a-4dbb-add0-6f097551c3e6)



(to do)

### Telomeric repeat motif discovery using reads from NCBI Sequence Read Archive
For this, TeloSearchLR requires these parameters.
|paramenter    | description                                                         |
|--------------|---------------------------------------------------------------------|
|-f            | FASTA file of the reads (STR)                                       |
|-k            | shortest period to search for (INT)                                 |
|-K            | longest period to search for (<=500 nt if -t is unspecified) (INT)  |
|-m            | numerical rank of the most frequent motif to plot (INT)             |
|-M            | numerical rank of the least frequent motif to plot (INT)            |
|-n            | number of nucleotides to plot the repeat occupancy (INT)            |

Let's try downloading a PacBio genomic sequencing library from *Caenorhabditis elegans* by [Yoshimura & al (2019)](https://pubmed.ncbi.nlm.nih.gov/31123080/) using ```fasterq-dump``` from [sra-tools](https://github.com/ncbi/sra-tools). Then, get TeloSearchLR to find and rank repeat motifs between 4 and 20-nt long (-k 4 -K 20) from the terminal 1000 nts (-t 1000). Plot the occupancy patterns of the top 100 motifs (-m 1 -M 100) in the first and last 6000 nt of all reads (-n 6000) longer than 12,000 nts.
```bash
fasterq-dump --fasta SRR7594465
python TeloSearchLR.py -f SRR7594465.fasta -k 4 -K 20 -t 1000 -m 1 -M 100 -n 6000
```
The [output of this](https://github.com/gchchung/TeloSearchLR/blob/main/repeatPattern.m1.M100.png) can be found in the *.results folder. The known telomeric repeat of *C. elegans*, TTAGGC, is the second most frequent repeat motif (second row of plots) and has the typical stranded occupancy pattern of telomeric repeat motifs. 

### Telomeric repeat motif discovery in your own library
Sequencing library reads are generally in the FASTQ format, which you must first convert into a FASTA file. A possible method is by using the Unix sed.
```bash
sed -n '1~4s/^@/>/p;2~4p' YOUR_LIBRARY.fastq > YOUR_LIBRARY.fasta
```
Then run TeloSearch.py on YOUR_LIBRARY.fasta.
```bash
python3 TeloSearch.py -f YOUR_LIBRARY.fasta -k 4 -K 20 -m 1 -M 100 -n 6000
```

### Unusually long telomeric repeat motifs
As run above, the longest telomeric motif that can be detected is 500 bps. To find repeat motifs longer than this, change the -t paramenter.
|paramenter    | description                                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------------------------|
|-t            | the terminal region (in bps) to rank repeat motifs. The *K*-value can be at most 1/2 this *t*-value (INT, default 1000) |

### Telomeric repeat motif discovery, with a known repeat period length
(to do)


### Testing to see if a repeat motif shows stranded occupancy at the read ends
(to do)


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

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.


## License

[NYU License](https://github.com/gchchung/TeloSearchLR/blob/main/LICENSE)

## Please cite
Chung & al. (2024) TeloSearchLR: an algorithm to detect novel telomere repeat motifs using long sequencing reads. bioRxiv [DOI: 10.1101/2024.10.29.617943](https://doi.org/10.1101/2024.10.29.617943).

## Contact
George Chung (gc95@nyu.edu)
