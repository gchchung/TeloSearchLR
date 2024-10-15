# TeloSearchLR

TeloSearchLR (**TELO**mere **SEARCH** using **L**ong sequencing **R**eads a Python script for aiding the identificaiton of telomeric repeat motifs.

## Installation

Using conda,

```bash
conda install _________________
```

## Usage

```bash
# Download a test Caenorhabditis elegans PacBio long-read dataset from SRA using sra-tools
fasterq-dump --fasta SRR7594465
```
```python
python3 TeloSearchLR.py -f SRR7594465.fasta -k 4 -K 20 -m 1 -M 100 -n 4000

```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

NYU License
