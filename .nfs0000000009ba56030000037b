#!/bin/bash

#SBATCH --job-name=PlastidPipeline
#SBATCH --mail-type=end
#SBATCH --mem=100
#SBATCH --time=1-0:00:00
#SBATCH --qos=standard
#SBATCH --output=PlastidPipeline.out
#SBATCH --error=PlastidPipeline.err

snakemake --profile resources/simple/ --configfile config/config.yaml --use-envmodules --use-conda

###SBATCH --mail-user=sjannalda@zedat.fu-berlin.de
