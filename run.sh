# Infering the tree from contextual sequences + closest sequences
nextflow run -c gofasta.config gofasta.nf
# Adding JN.1 sequence
nextflow run -c gofasta.config add_jn1.nf

# Contextual sequence IDs taken from nextstrain build global GISAID from ~18/09/2023
xz -d -c data/metadata_tsv_2023_09_16.tar.xz| tar -xOvf - metadata.tsv | cut -f 1,5,6,7,14 > data/meta_shrink.tsv
wget https://github.com/cov-lineages/pango-designation/raw/e6a71b06d7b2e0649f222c52af4bb38bd7af3f3c/lineage_notes.txt
perl bin/simplifyPango.pl data/meta_shrink.tsv  lineage_notes.txt > data/meta_shrink3.tsv
grep -f <(gotree labels -i results/align_fbp.treefile) data/meta_shrink3.tsv > data/meta_shrink3_tree.tsv

# Producing annotation files for itol
grep -f <(gotree labels -i results/align_fbp.treefile) results/gisaid_pangolin.tsv | cut -f 1,2 -d "," > results/gisaid_pangolin_subsample.tsv
singularity pull docker://evolbioinfo/table2itol:fa4b43c
singularity exec --bind /pasteur table2itol-fa4b43c.img table2itol.R -i ID -s "," results/gisaid_pangolin_subsample.tsv
bin/simplifyPango2.pl  results/gisaid_pangolin_subsample.tsv  lineage_notes.txt > gisaid_pangolin_subsample_simple.tsv
singularity exec  --bind /pasteur table2itol-fa4b43c.img table2itol.R -i ID -s ","  results/gisaid_pangolin_subsample_simple.tsv
perl bin/domain2branch.pl results/iTOL_domains-short_pango.txt > results/iTOL_domains-short_pango_tips.txt

# Statistics up to October
wget -O data/lineage_notes_oct.txt https://raw.githubusercontent.com/cov-lineages/pango-designation/1bf41233556859320bad75e269f595690294098a/lineage_notes.txt
xz -d -c data/metadata_tsv_2023_10_28.tar.xz| tar -xOvf - metadata.tsv | cut -f 1,5,6,7,14 > data/meta_shrink_oct.tsv
perl bin/simplifyPango.pl data/meta_shrink_oct.tsv data/lineage_notes_oct.txt > data/meta_shrink3_oct.tsv

# Statistics up to November
wget -O data/lineage_notes_nov.txt https://raw.githubusercontent.com/cov-lineages/pango-designation/e4da6b33910e7694ccdbd91892e53f26bf85ccca/lineage_notes.txt
xz -d -c data/metadata_tsv_2023_11_27.tar.xz| tar -xOvf - metadata.tsv | cut -f 1,5,6,7,14 > data/meta_shrink_nov.tsv
perl bin/simplifyPango.pl data/meta_shrink_nov.tsv data/lineage_notes_nov.txt > data/meta_shrink3_nov.tsv

# Statistics up to almost December
wget -O data/lineage_notes_dec.txt https://raw.githubusercontent.com/cov-lineages/pango-designation/a0837ab84b0d2066a910afef1642c6ccfd10da12/lineage_notes.txt
xz -d -c data/metadata_tsv_2023_12_11.tar.xz| tar -xOvf - metadata.tsv | cut -f 1,5,6,7,14 > data/meta_shrink_dec.tsv
perl bin/simplifyPango.pl data/meta_shrink_dec.tsv data/lineage_notes_dec.txt > data/meta_shrink3_dec.tsv

# Statistics up to January
wget -O data/lineage_notes_jan.txt https://raw.githubusercontent.com/cov-lineages/pango-designation/a0507852280782373a70a21e272d14228352e9d8/lineage_notes.txt
xz -d -c data/metadata_tsv_2024_01_15.tar.xz| tar -xOvf - metadata.tsv | cut -f 1,5,6,7,14 > data/meta_shrink_jan.tsv
perl bin/simplifyPango.pl data/meta_shrink_jan.tsv data/lineage_notes_jan.txt > data/meta_shrink3_jan.tsv
