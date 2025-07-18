# PlastidPipeline

A software pipeline for the automated, reproducible and well-documented production of assembled, annotated and standardized plastid genomes, developed by Siddharth Annaldasula and Katja Reichel.

Full documentation: [PlastidPipeline Google Docs](https://docs.google.com/document/d/1kSNYbWYWll7QQv1zzU0TqP_HkZW93I5C3cAuKrxtv6g/edit?usp=sharing)

To run the pipeline on the HPC (remember to adapt the conda prefix & profile paths):
snakemake --use-conda --use-envmodules --resources chrome=1 --conda-prefix /path/to/envs --configfile config/config.yaml --profile ~/plastidpipeline/resources

For installation:
1. Install conda + mamba
2. Install snakemake into its own "snakemake" environment
3. Git clone this repository & move into it
4. Run once (remember to adapt the conda prefix path):
`snakemake --use-conda --conda-prefix /path/to/envs --conda-create-envs-only --configfile config/config.yaml -c1`

No warranty – always visually check the output.
