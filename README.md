# HaploPhase-dev
Haplotype phasing algorithm for diploid genome(under development)

# Description
Our phasing algorithm can adapt the assembly result from various TGS assemblers and generate phased assemblies. Currently we have evaluated against [Verkko](https://github.com/marbl/verkko) and [Hifiasm](https://github.com/chhylp123/hifiasm).

# Dependencies
## Dependencies for main algorithm
Manually install dependencies:
- [Kmer2SNP](https://github.com/yanboANU/Kmer2SNP)
- [HiC-Pro](https://github.com/nservant/HiC-Pro)

And python modules:
- [graph-tool](https://graph-tool.skewed.de)
- [numpy](https://numpy.org)
- [gfapy](https://github.com/ggonnella/gfapy)
- [matplotlib](https://matplotlib.org)
- [pysam](https://github.com/pysam-developers/pysam)


## Dependencies for evaluation against reference genome
Manually install dependencies:
- [MUMMER](https://anaconda.org/bioconda/mummer)
- [GenomeViz](https://github.com/moshi4/pyGenomeViz)
- [matplotlib](https://matplotlib.org)
- [seaborn](https://seaborn.pydata.org/installing.html)
- [pandas](https://pandas.pydata.org/docs/getting_started/install.html)
- [numpy](https://numpy.org)

# Usage

```sh
src/main.py <gfa_file> <asm_file> <bed> <matrix> <snp_file> <snp_ext_file> <output_dir>
```

`gfa_file` and `asm_file` correspond to the assembly graph and assembled contigs from assembler, respectively. `bed` and `matrix` correspond to the alignment result from HiC-Pro, `snp_file` and `snp_ext_file` correspond to the called snps from Kmer2SNP, and `output_dir` corresponds to the final output directory.

## Inputs
we expect the LA reads (such as Oxford Nanopore Duplex or PacBio Hifi reads), UL reads (such as Oxford Nanopore Simplex reads), and HiC reads as inputs. Homozygous coverage of LA reads is required by Kmer2SNP prior phasing, which can be estimated by comparing against similar species or perform a pre-assembly.

## Outputs
The program stores all output files in `<output_dir>`, which is set by the user.

`<output_dir>/bin0.fasta` and `<output_dir>/bin1.fasta` store the phased contigs belong to haplotype 0 and haplotype 1, respectively. `<output_dir>/binN.fasta` stores the unphased contigs.

If the assembly graph is provided, `<output_dir>/tree_bin0.fasta`, `<output_dir>/tree_bin1.fasta`, and `<output_dir>/tree_binN.fasta` include the phasing result after running Steiner tree module.

## Support Verkko
When incorporating Verkko assembler, you can obtain the uncompressed assembly graph via [Verkko issue](https://github.com/marbl/verkko/issues/140).

```sh
verkko -d output --hifi duplex.fastq --nano simplex.fastq
```

## Support Hifiasm
When incorporating Hifiasm assembler, please use `prefix.p_utg.gfa` as input. Assembled contigs in fasta format can be obtain by following the instruction described at the [Hifiasm](https://github.com/chhylp123/hifiasm) main page.

```sh
hifiasm -o <output> -t 16 --ul simplex.fastq duplex.fastq
```
# Evaluation
Various evaluation scripts to generate haplotype phasing scatter plot and chromosome-level phasing result are provided under `plot_scripts`, with incorporation of [QUAST](https://github.com/ablab/quast). All the scripts can be executed under [Jupyter](https://jupyter.org).

## Contact
Currently, the program is under active development, please raise a bug by directly contacting john.luo@anu.edu.au or upload via `Issues`.