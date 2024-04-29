## SPLITREADER_v1.2

### Install Conda

- **Windows:** [Installation Conda on Windows](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)

- **Linux:** [Installation Conda on Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

### Clone the SPLITREADER branch

git clone https://github.com/baduelp/public/SPLITREADER.git


### Run SPLITREADER 

#### First Step: 

Create a $workspace_dir from where you will run the SPLITREADER pipeline.

In this $workspace_dir, you will need the following setup files:

**in $workspace_dir/Reference**

- $genome.fasta file with the reference genome sequence (for *A. thaliana* the TAIR10 reference genome fasta can be downloaded **[here](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_chromosome_files))**

**in $workspace_dir/TE_sequence**

- superfamily_TSD.txt : *Tab-delimited file with superfamilies in the 1st column and their respective Target Site Duplication (TSD) lengths in the 2nd column.*
- TE-list.txt : *file with a list of the names of the TEs annotated in the reference genome.*
- $TE_lib.fasta : *file with the sequence of all the TEs annotated in the reference genome.*
- TEfamily-superfamily.txt : *file with TE names in the 1st column and their respective superfamily in the 2nd column.*
- $TE_annotation.gff *file of the TEs annotated.*
- 
#### Second step: 

Fill in the `config.yml` file with the full paths corresponding to the location of your setup and input files ($workspace_dir, $ref_dir, etc) and define the names of the filename variables ($genome, $TE_lib, etc) to match those above.

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
