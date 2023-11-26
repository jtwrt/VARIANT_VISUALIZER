# VARIANT_VISUALIZER

> README is work in progress!

Visualizing genomic regions, transcripts and proteins, their known functional regions, regulatory elements and the mutations affecting them.

## Setup

### Prepare config.yml

At a minimum, the following values need to be defined in the `config.yml` (general section):
- `init_bedtools`: Bash command that enables use of the `bedtools` command in the command line.
- `ucsc_liftover`: Path of the UCSC-Liftover executable.

### Prepare conda enviroment

Create conda enviroment:
```bash
conda env create -f conda_environment.yml
```
Activate the enviroment:
```bash
conda activate variant_visualizer
```

### Download dependencies

Run `setup_dependencies.py`. The script automatically downloads all dependencies not mentioned in the previous step and prepares them for further use.

```bash
python setup_dependencies.py
```

### Generate clusters

Generate genomic clusters by running `setup_clusters.py` passing the number of parallel processes as argument.
Generating all clusters is currently very resource intensive. For each parallel process, 20 GB of memory are advised. Use as many parallel processes as you can afford.

```bash
python setup_clusters.py 1
```

### Set up the ClusterIndex

Run `setup_index.py` to generate a index of the previously generated clusters.

```bash
python setup_index.py
```

## Usage

Read `vignette.ipynb` for exampels and detailed explaination on how this package can be used. 