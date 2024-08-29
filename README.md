# PlastidPipeline

Current name: Plastid Pipeline

The current (in development) documentation: [PlastidPipeline Google Docs](https://docs.google.com/document/d/1kSNYbWYWll7QQv1zzU0TqP_HkZW93I5C3cAuKrxtv6g/edit?usp=sharing)


snakemake --profile resources/simple/ --configfile config/config.yaml --resources chrome=1 --use-envmodules --use-conda --conda-prefix /scratch/sjannalda/projects/PlastidTutorial/envs --restart-times 3

To run the pipeline on the HPC: `snakemake --profile resources/simple/ --configfile config/config.yaml --resources chrome=1 --use-envmodules --use-conda`
