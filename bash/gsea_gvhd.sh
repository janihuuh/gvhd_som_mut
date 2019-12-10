#!/bin/bash
me=$(whoami)
RNKPATH=/Users/$me/Dropbox/gvhd_scrnaseq/
GSEA=/Users/$me/Dropbox/applications/gsea-3.0.jar ## gsea
# GSEA=/Users/janihuuh/Dropbox/gvhd_scrnaseq/applications/gsea/GSEA_4.0.2/gsea-cli.sh

GMT=/Users/$me/Dropbox/applications/gsea/h.all.v6.2.symbols.gmt # path to gmt file to count the enrichment against
OUTDIR=$RNKPATH/results/gsea/ #out directory




### Total
RNKFILE=$RNKPATH/results/mutated_clonotype/cluster_markers_clonotype_for_gsea.rnk ## file with genes ranked by fold-change
LABEL=clonotype # prefix-label for subfolder of $OUTDIR

java  -cp $GSEA \
      -Xmx5g xtools.gsea.GseaPreranked \
      -rpt_label $LABEL \
      -rnk $RNKFILE \
      -gmx $GMT \
      -out $OUTDIR \
      -plot_top_x 250 \
      -collapse false \
      -mode Max_probe \
      -norm meandiv \
      -scoring_scheme weighted \
      -include_only_symbols true \
      -make_sets true \
      -rnd_seed 149 \
      -zip_report false \
      -gui false \
      -nperm 1000 \
      -set_min 5 \
      -set_max 500



### Cluster 0 vs 1

RNKFILE=$RNKPATH/results/mutated_clonotype/cluster_0v1_markers_clonotype_for_gsea.rnk ## file with genes ranked by fold-change
LABEL=cluster0v1 # prefix-label for subfolder of $OUTD
java  -cp $GSEA \
      -Xmx5g xtools.gsea.GseaPreranked \
      -rpt_label $LABEL \
      -rnk $RNKFILE \
      -gmx $GMT \
      -out $OUTDIR \
      -plot_top_x 250 \
      -collapse false \
      -mode Max_probe \
      -norm meandiv \
      -scoring_scheme weighted \
      -include_only_symbols true \
      -make_sets true \
      -rnd_seed 149 \
      -zip_report false \
      -gui false \
      -nperm 1000 \
      -set_min 5 \
      -set_max 500


### Cluster 0
RNKFILE=$RNKPATH/results/mutated_clonotype/cluster0_markers_clonotype_for_gsea.rnk ## file with genes ranked by fold-change
LABEL=cluster0 # prefix-label for subfolder of $OUTDIR

java  -cp $GSEA \
      -Xmx5g xtools.gsea.GseaPreranked \
      -rpt_label $LABEL \
      -rnk $RNKFILE \
      -gmx $GMT \
      -out $OUTDIR \
      -plot_top_x 250 \
      -collapse false \
      -mode Max_probe \
      -norm meandiv \
      -scoring_scheme weighted \
      -include_only_symbols true \
      -make_sets true \
      -rnd_seed 149 \
      -zip_report false \
      -gui false \
      -nperm 1000 \
      -set_min 5 \
      -set_max 500
