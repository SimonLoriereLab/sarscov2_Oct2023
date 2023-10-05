# Contextual sequence IDs taken from nextstrain build global GISAID from ~18/09/2023
xz -d -c metadata_tsv_2023_09_16.tar.xz| tar -xOvf - metadata.tsv | cut -f 1,5,6,7,14 > meta_shrink.tsv
perl ../bin/simplifyPango.pl meta_shrink.tsv > meta_shrink2.tsv
grep -f <(gotree labels -i ../results/align.treefile) meta_shrink2.tsv > meta_shrink2_tree.tsv
bin/simplifyPango.pl data/meta_shrink2_tree.tsv > data/meta_shrink2_tree_short.tsv

wget https://github.com/cov-lineages/pango-designation/raw/e6a71b06d7b2e0649f222c52af4bb38bd7af3f3c/lineage_notes.txt
perl bin/simplifyPango.pl data/meta_shrink.tsv  data/lineage_notes.txt > data/meta_shrink3.tsv

grep -f <(gotree labels -i results/align_fbp.treefile) data/meta_shrink3.tsv > data/meta_shrink3_tree.tsv


grep -f <(gotree labels -i align_fbp.treefile) gisaid_pangolin.tsv | cut -f 1,2 -d "," > gisaid_pangolin_subsample.tsv
shell --bind /pasteur ../singularity/evolbioinfo-table2itol-fa4b43c.img
table2itol.R -i ID -s "," gisaid_pangolin_subsample.tsv

../bin/simplifyPango2.pl  gisaid_pangolin_subsample.tsv  ../data/lineage_notes.txt > gisaid_pangolin_subsample_simple.tsv
shell --bind /pasteur ../singularity/evolbioinfo-table2itol-fa4b43c.img
table2itol.R -i ID -s ","  gisaid_pangolin_subsample_simple.tsv

perl ../bin/domain2branch.pl iTOL_domains-short_pango.txt > iTOL_domains-short_pango_tips.txt
