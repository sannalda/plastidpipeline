#!/bin/bash

#SBATCH --mail-user=sjannalda@zedat.fu-berlin.de
#SBATCH --job-name=PlastidPipeline
#SBATCH --mail-type=end
#SBATCH --mem=100
#SBATCH --time=0-2:00:00
#SBATCH --qos=standard
#SBATCH --output=PlastidPipeline.out
#SBATCH --error=PlastidPipeline.err

snakemake --profile resources/simple/ --configfile config/config.yaml --use-envmodules --use-conda 


snakemake --profile resources/simple/ --configfile config/config.yaml --use-envmodules --use-conda

$ /home/sjannalda/opt/sratoolkit.3.0.0-centos_linux64/bin/prefetch ERR5554746

 /home/sjannalda/opt/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --split-files ERR5554746