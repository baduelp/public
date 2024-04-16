# SPLITREADER 

This public repository contains the scripts for running the SPLITREADER pipeline initially developped by [Quadrana et al. eLife 2016](https://doi.org/10.7554/eLife.15716) and further improved by [Baduel et al. MMB 2021](https://doi.org/10.1007/978-1-0716-1134-0_15). 

If you have any question or comment please reach me through the channels listed [there](http://pbaduel.com/about). 

## [/SPLITREADER](/SPLITREADER) 
This folder contains the scripts to detect non-reference transposable element insertions from short-read sequencing data as described in [Baduel et al. MMB 2021](https://doi.org/10.1007/978-1-0716-1134-0_15). <br/>
Accessory files and wrapper scripts are also available in [/thaliana](/SPLITREADER/thaliana) to run the SPLITREADER pipeline on the _A. thaliana_ genome. <br/>

## [/SPLITREADER_V1](/SPLITREADER_V1)
This folder contains a snakemake version of the SPLITREADER pipeline produced by [@aurelpetit](https://github.com/aurelpetit) for enhanced portability and ease of use. 

![workflow](SPLITREADER_v1/dag_workflow.svg)
### Clone the SPLITREADER branch

git clone https://github.com/baduelp/public/SPLITREADER_V1.git

### Installation Conda

- **Windows:** [Installation Conda on Windows](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)

- **Linux:** [Installation Conda on Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

### Run SPLITREADER 

#### First Step: 

Fill the file `config.yml` with full paths.

#### Second step: 

You need the following files to run the SPLITREADER pipeline :

**in $workspace_dir/$ref_dir**

- $genome.fasta file with the reference genome sequence (test: *A.thaliana* => **[link to the genome file](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_chromosome_files))**

**in $workspace_dir/TE_sequence**

- superfamily_TSD.txt : *Tab-delimited file with superfamilies in the 1st column and their respective Target Site Duplication (TSD) lengths in the 2nd column.*
- TE-list.txt : *file with a list of the names of the TEs annotated in the reference genome.*
- $TE_lib.fasta : *file with the sequence of all the TEs annotated in the reference genome.*
- TEfamily-superfamily.txt : *file with TE names in the 1st column and their respective superfamily in the 2nd column.*
- $TE_annotation.gff *file of the TEs annotated.*

#### Third step:

1. **Create the conda environment :**

    ```bash
    conda env create -f SPLITREADER.yml
    conda activate SPLITREADER
    ```

2. **Run the snakefile :**

    ```bash
    snakemake --snakefile SPLITREADER.snk --cores X
    ```

## References

**Baduel P**, Quadrana L, Colot V. Efficient detection of transposable element insertion polymorphisms between genomes using short-read sequencing data. Plant Transposable Elements, Methods in Molecular Biology, [10.1007/978-1-0716-1134-0_15](https://doi.org/10.1007/978-1-0716-1134-0_15), 04/2021.


---

