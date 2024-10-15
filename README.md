# TeloSearchLR

TeloSearchLR (**telo**mere **search** using **l**ong sequencing **r**eads) is a Python script for aiding the identificaiton of telomeric repeat motifs.

## Installation

Using conda,

```bash
conda install _________________
```

## Usage

To get started, download a test dataset - a genomic sequencing library from *Caenorhabditis elegans* generated using PacBio (using ```fasterq-dump``` from [sra-tools](https://github.com/ncbi/sra-tools)).
```bash
# Download a test Caenorhabditis elegans PacBio long-read dataset from SRA using sra-tools
fasterq-dump --fasta SRR7594465
```
Wait a bit until the download is complete, then run TeloSearch.py.
```bash
# Run TeloSearchLR
python3 TeloSearchLR.py -f SRR7594465.fasta -k 4 -K 20 -m 1 -M 100 -n 4000
```
The algorithm will look for tandem repeats of period 4-20 bp (```-k 4 -K 20```)in reads ≥ 8000bp (2×4000, specified through ```-n 4000```) and rank each tandem repeat motif based on its occupancy in first and last 1000 bp.  The algorithm then plots the occupancy of the top 100 patterns ranked this way (```-m 1 -M 100```).

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

Please make sure to update tests as appropriate.

## License

NYU License

## References
Please cite this repository. A preprint describing this algorithm will be posted shortly.

## Contact
George Chung (gc95@nyu.edu)
