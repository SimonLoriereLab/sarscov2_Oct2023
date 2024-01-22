nextflow.enable.dsl=2

// Download nextclade reference datasets
process getRefData {

label 'getref'

output:
file "ncref"

script:
"""
nextclade dataset get --name 'sars-cov-2' --output-dir 'ncref'
"""
}

// Align jn1 sequences using nextalign and
// reference nextclade dataset
process alignJN1 {

label 'align'

input:
file ncref
file jn1

output:
file "jn1.aligned.fasta.gz"

script:
"""
nextalign run  -j ${task.cpus} --input-ref ${ncref}/reference.fasta  --input-gene-map ${ncref}/genemap.gff --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --output-all . --output-basename jn1 ${jn1}
gzip jn1.aligned.fasta
"""
}


// Mask specific positions of the alignment
process mask {
   publishDir "results_jn1", mode: 'copy'

   label 'goalignbigmem'

   input:
   file fasta

   output:
   file "*_masked.fasta.gz"

   script:
   """
   gunzip -c $fasta | goalign mask -s 0 -l 55 --ignore-identical 1 \
     | goalign mask -s 29803 -l 101 \
     | goalign mask --pos 186,1058,2093,3036,3129,6989,8021,10322,10740,11073,13407,14785,19683,20147,21136,24033,24377,25562,26143,26460,26680,28076,28825,28853,29699,4049,13401,11082,15323,21574,21845,21986,22992,23011,23062,23603,24409 \
     | gzip -c - > ${fasta.baseName}_masked.fasta.gz
   """
}

process AddJN1ToTree {
    publishDir "results_jn1", mode: "copy"

    label 'gotree'

    input:
    path tree
    path toadd

    output:
    path "starttree_jn1.treefile"

    script:
    """
    gotree repopulate -i $tree -g $toadd -o starttree_jn1.treefile
    """
}




// Build reference tree
process Phylogeny {
    publishDir "results_jn1", mode: 'link'

    label 'phylogeny'

    input:
    path tree
    path align

    output:
    path "*.treefile"
    path "align"
    
    script:
    """
    gunzip -c $align > align
    iqtree -ntmax ${task.cpus} -s align -m GTR -ninit 2 -n 2 -me 0.05 -nt AUTO -redo -ninit 10 -n 4 -t $tree --seed 123456789
    """
}

workflow {

fasta = Channel.fromPath("results/*_masked*")
jn1 = Channel.fromPath("data/JN.1.fasta")
starttree = Channel.fromPath("results/align.treefile")
addJN1 = Channel.fromPath("data/JN.1_add.txt")

starttreejn1 = AddJN1ToTree(starttree,addJN1)
ncref=getRefData()
jn1al=mask(alignJN1(ncref, jn1))

inalign = jn1al.mix(fasta).collect()

finaltree = Phylogeny(starttreejn1, inalign)
}