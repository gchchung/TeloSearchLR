![alt text](https://github.com/gchchung/TeloSearchLR/blob/main/logo_copy_v3.svg)
TeloSearchLR (**telo**mere **search** using **l**ong sequencing **r**eads) is a Python script for aiding the identification of telomeric repeat motifs.

## Installation

### Conda
Install in a new environment named "telosearchlr-env".
```bash
# Create a new conda environment, install TeloSearchLR and dependencies, then activate the environment
conda create -n telosearchlr-env -c bioconda -y telosearchlr
conda activate telosearchlr-env

# In the next release we will fix this, but our original conda recipe didn't include "python-kaleido"
conda install -c conda-forge -y python-kaleido

# Test the installation by calling the help function
TeloSearchLR.py -h

# De-activate the conda environment before using other environments
conda deactivate
```
Note that if you install via conda, you will not need to call the python command.


### Docker
```bash
# Pull the TeloSearchLR image from Docker repository ```gchchung/telosearchlr```. Test by asking for the help message.
docker pull gchchung/telosearchlr:telosearchlr_v1.0.0

# Then display the help message.
docker run telosearchlr_v1.0.0 TeloSearchLR.py -h
```

### From source
```bash
# clone repo
git clone https://github.com/gchchung/TeloSearchLR.git
cd TeloSearchLR

# install dependencies (not shown), then run script
python TeloSearch.py -h
```


## Usage
### Using a library from NCBI Sequence Read Archive
Download a test dataset - a genomic sequencing library from *Caenorhabditis elegans* generated using PacBio by [Yoshimura & al (2019)](https://pubmed.ncbi.nlm.nih.gov/31123080/). You will need to have ```fasterq-dump``` from [sra-tools](https://github.com/ncbi/sra-tools) installed.
```bash
fasterq-dump --fasta SRR7594465
```
Wait for the download to complete, then run TeloSearch.py on the downloaded library.
```bash
python TeloSearchLR.py -f SRR7594465.fasta -k 4 -K 20 -t 1000 -m 1 -M 100 -n 6000
```
The algorithm will look for tandem repeats of period 4-20 bp (```-k 4 -K 20```) in the terminal 1000 bps (specified through ```-t 1000```) on reads ≥ 2000bp, and rank each tandem repeat motif based on its occupancy in those terminal 1000 bps.  The algorithm then plots the occupancy of the top 100 patterns ranked this way (```-m 1 -M 100```), at the first and last 6000 bps of all reads 12,000 bps or longer (```-n 6000```). The [output of this](https://github.com/gchchung/TeloSearchLR/blob/main/repeatPattern.m1.M100.png) can be found in the *.results folder. The known telomeric repeat of *C. elegans*, TTAGGC, is the second most frequent repeat motif (second row of plots) and has the typical stranded occupancy pattern of telomeric repeat motifs. 

### Using your own sequencing library
Sequencing library reads are generally in the FASTQ format, which you must first convert into a FASTA file. A possible method is by using the Unix sed.
```bash
sed -n '1~4s/^@/>/p;2~4p' YOUR_LIBRARY.fastq > YOUR_LIBRARY.fasta
```
Then run TeloSearch.py on YOUR_LIBRARY.fasta.
```bash
python3 TeloSearch.py -f YOUR_LIBRARY.fasta -k 4 -K 20 -m 1 -M 100 -n 6000
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
## Sample use cases
### Telomeric repeat motif discovery
In this use case, TeloSearchLR requires these parameters.
|paramenter    | description                                                         |
|--------------|---------------------------------------------------------------------|
|-f            | FASTA file of the reads (STR)                                       |
|-k            | shortest period to search for (INT)                                 |
|-K            | longest period to search for (<=500 nt if -t is unspecified) (INT)  |
|-m            | numerical rank of the most frequent motif to plot (INT)             |
|-M            | numerical rank of the least frequent motif to plot (INT)            |
|-n            | number of nucleotides to plot the repeat occupancy (INT)            |                      

### Unusually long telomeric repeat motifs
As run above, the longest telomeric motif that can be detected is 500 bps. To find repeat motifs longer than this, change the -t paramenter.
|paramenter    | description                                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------------------------|
|-t            | the terminal region (in bps) to rank repeat motifs. The *K*-value can be at most 1/2 this *t*-value (INT, default 1000) |


### Telomeric repeat motif discovery, with a known repeat period length

### Testing to see if a repeat motif shows stranded occupancy at the read ends

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.


## License

[NYU License](https://github.com/gchchung/TeloSearchLR/blob/main/LICENSE)

## Please cite
Chung & al. (2024) TeloSearchLR: an algorithm to detect novel telomere repeat motifs using long sequencing reads. bioRxiv [DOI: 10.1101/2024.10.29.617943](https://doi.org/10.1101/2024.10.29.617943).

## Contact
George Chung (gc95@nyu.edu)
