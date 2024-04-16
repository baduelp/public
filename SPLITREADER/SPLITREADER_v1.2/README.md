## SPLITREADER_V1.2

### Clone the SPLITREADER branch

git clone https://github.com/baduelp/public/SPLITREADER.git

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
