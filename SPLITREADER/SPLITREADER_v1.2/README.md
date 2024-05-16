## SPLITREADER_v1.2

### Installing Conda

- **Windows:** [Installation Conda on Windows](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)

- **Linux:** [Installation Conda on Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

### Getting the SPLITREADER

`git clone https://github.com/baduelp/public.git`

### Setting and running the SPLITREADER 

#### First Step: organizing the workspace with the needed meta files

Create a `$WORKSPACE_DIR` from where you will run the SPLITREADER pipeline.

In this `$WORKSPACE_DIR`, you will need the following meta files:

**in `WORKSPACE_DIR/Reference`**

- `$GENOME.fasta` file with the reference genome sequence (for *A. thaliana* the TAIR10 reference genome fasta can be downloaded **[here](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_chromosome_files))**

**in `WORKSPACE_DIR/TE_sequence`**

- `superfamily_TSD.txt` : *Tab-delimited file with superfamilies in the 1st column and their respective Target Site Duplication (TSD) lengths in the 2nd column.*
- `TE-list.txt` : *file with a list of the names of the TEs annotated in the reference genome.*
- `$TE_LIB.fasta` : *file with the sequence of all the TEs annotated in the reference genome.*
- `TEfamily-superfamily.txt` : *file with TE names in the 1st column and their respective superfamily in the 2nd column.*
- `$TE_ANNOTATION.gff`: *file of the TEs annotated.*

#### Second step: setting up the variables

Fill in the `config.yml` file with the full paths corresponding to the location of your input and meta files ( `$WORKSPACE_DIR`, `$REFERENCE_DIR`, etc) and define the names of the filename variables ( `$GENOME`, `$TE_LIB`, etc) to match those above.

You also have the choice, in `$CHOICE`, of using the Negative coverage filter. If `$NEGATIVE_COVERAGE` is set to `False`, you can put any value for `$DPrefmin` as it will not be used.

#### Third step: running the pipeline

1. **Creating the conda environment :**

    ```bash
    conda env create -f SPLITREADER.yml
    conda activate SPLITREADER
    ```

2. **Running the snakefile :**

    ```bash
    snakemake --snakefile SPLITREADER.snk --cores X
    ```
`--cores X` : X being the number of cores you allocate for the run. 

--------------

#### About the multi-threading :

In `config.yml` you can give a value for `$SNAKEMAKE_THREADS`, that specify the number of threads that some rules can use. For instance, if `$SNAKEMAKE_THREADS: 10` each multi_threadable rule will use 10 threads.

> [!IMPORTANT]
> Be careful with the number of cores you give via `--cores X` :
> - If `--cores 10` and `$SNAKEMAKE_THREADS: 10`, then each multi-threadable rule will use 10 threads, but 10 threads are allocated here. So only one sample of your data will be treated at once.
> - If `--cores 10` and `$SNAKEMAKE_THREADS: 2`, then each multi-threadable rule will use 2 threads, and 10 threads are allocated, so 5 samples will be treated at once. 
