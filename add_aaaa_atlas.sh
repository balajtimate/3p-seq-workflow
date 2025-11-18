#!/bin/bash

#SBATCH --job-name=add_aaaa_atlas
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --output=add_aaaa_atlas.o%j
#SBATCH --error=add_aaaa_atlas.e%j

source /scicore/home/zavolan/balajt0000/.bashrc
conda activate paqr3
python /scicore/home/zavolan/balajt0000/Work/PTA_immune/3p-seq-workflow/modules/atlas_motif_AAAA_analysis/add_atlas_AAAA.py \
      -i /scicore/home/zavolan/balajt0000/Work/PTA_immune/3p-seq-workflow/modules/atlas_motif_AAAA_analysis/atlas.clusters.3.0.GRCh38.GENCODE_42.bed \
      -f /scicore/home/zavolan/balajt0000/Work/PTA_immune/3p-seq-workflow/modules/atlas_motif_AAAA_analysis/GRCh38.primary_assembly.genome.fa \
      -o /scicore/home/zavolan/balajt0000/Work/PTA_immune/3p-seq-workflow/modules/atlas_motif_AAAA_analysis/atlas.clusters.3.0.GRCh38.GENCODE_42.with_AAAA.bed
