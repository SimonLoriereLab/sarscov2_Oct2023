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

// Remove sequences that are <29000nt or have >1% N
// And splits the dataset in files containing 100000 sequences
process filterGISAID {

label 'filter'

input:
file gisaid
file ba286names
file contextnames

output:
file "*_filter*.fasta.gz"
file "ba286.fasta.gz"
file "context.fasta.gz"

script:
"""
xz -d -c ${gisaid} | tar -xOvf - sequences.fasta| filterSequences.py -s - -S 100000 -o ${gisaid.baseName}_filter -B ${ba286names} -b ba286.fasta.gz -C ${contextnames} -c context.fasta.gz -t ${task.cpus-2}
"""
}

// Annotate input sequences using pangolin
process Pangolin {
    label 'pangolin'
    
    input:
    file seq

    output:
    path "lineage_report.csv"

    script:
    """
    PATH=/opt/conda/bin/:\$PATH
    pangolin --usher '${seq}' -t ${task.cpus} --outfile lineage_report.csv
    """
}

// Align input sequences using nextalign and
// reference nextclade dataset
process alignGISAID {

label 'align'

input:
file ncref
file seq

output:
file "*.aligned.fasta.gz"

script:
"""
nextalign run  -j ${task.cpus} --input-ref ${ncref}/reference.fasta  --input-gene-map ${ncref}/genemap.gff --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --output-all . --output-basename ${seq.baseName}.fasta ${seq}
gzip ${seq.baseName}*.aligned.fasta
"""
}

// Align contextual sequences using nextalign and
// reference nextclade dataset
process alignContext {

label 'align'

input:
file ncref
file context

output:
file "context.aligned.fasta.gz"

script:
"""
nextalign run  -j ${task.cpus} --input-ref ${ncref}/reference.fasta  --input-gene-map ${ncref}/genemap.gff --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --output-all . --output-basename context ${context}
gzip context.aligned.fasta
"""
}

// Align ba286 sequences using nextalign and
// reference nextclade dataset
process alignBA286 {

label 'align'

input:
file ncref
file ba286

output:
file "ba286.aligned.fasta.gz"

script:
"""
nextalign run  -j ${task.cpus} --input-ref ${ncref}/reference.fasta  --input-gene-map ${ncref}/genemap.gff --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --output-all . --output-basename ba286 ${ba286}
gzip ba286.aligned.fasta
"""
}

// Create the database containing contextual + gisaid filtered sequences
process ConcatDB {
    
    input:
    path gisaid
    path context
    
    output:
    path "db_all.fasta.gz"
    
    script:
    """
    cat $gisaid $context > db_all.fasta.gz
    """
}

// Search the closest sequences from the query in the database
// using gofasta
process GoFasta {
    
    label 'gofasta'
    
    publishDir "results", mode: 'link'
    
    input:
    file db
    file query
    
    output:
    file "*closest.txt"
    
    script:
    """
    gunzip -c $query > query
    gunzip -c $db | gofasta closest -t ${task.cpus} --query query -n 200 --target stdin -o ${query.baseName}_closest.txt
    """
}

// Process the results of gofasta, and deduplicate output sequence names
// because several query sequences can have intersecting sets of closest sequences
process UniqueClosest {

    label 'uniq'
    
    publishDir "results", mode : 'link'
    
    input:
    file closest
    
    output:
    file "*_unique.txt"
    
    script:
    """
    cut -f 2 -d ',' ${closest} | sed 's/;/\\n/g' | grep -v "closest" | sort -u > ${closest.baseName}_unique.txt
    """
}

// Given the name of he input sequences, retrieve their sequence in the input fasta file
process ExtractSequences {
    
    label 'extract'
    
    publishDir "results", mode : 'link'
    
    input:
    file fasta
    file names
    
    output:
    file "*_sequences.fasta.gz"
    
    script:
    """
    extractSequences.py -i ${fasta}  -n ${names} -o ${names.baseName}_sequences.fasta.gz
    """
}

// Mask specific positions of the alignment
process mask {
   publishDir "results/", mode: 'copy'

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

// Generate bootstrap alignments
process BootAligns {
    publishDir "results/bootaligns/", mode: 'copy'

    label 'goalignbigmemcpus'

    input:
    path fasta

    output:
    path "boot_*"

    script:
    """
    gunzip -c $fasta | goalign build seqboot -n 100 -o boot_ -S --gz --seed 123456789 -t ${task.cpus}
    """
}

// Build reference tree
process Phylogeny {
    publishDir "results/", mode: 'copy'

    label 'phylogeny'

    publishDir "results", mode: 'link'

    input:
    path align

    output:
    path "*.treefile"

    script:
    """
    gunzip -c $align > align
    #iqtree -nt ${task.cpus} --safe -s align -m GTR+G4+FO 
    iqtree -ntmax ${task.cpus} -s align -m GTR -ninit 2 -n 2 -me 0.05 -nt AUTO -redo -ninit 10 -n 4  --seed 123456789
    rm align
    """
}

// Build bootstrap trees
process BootTrees {
    publishDir "results/boottrees", mode: 'copy'

    label 'phylogeny'

    input:
    path align

    output:
    path "*.treefile"

    script:
    """
    gunzip -c $align > boot.fasta
    iqtree -nt ${task.cpus} --safe -s boot.fasta -m GTR+G4+FO  --seed 123456789
    rm boot.fasta
    """
}

// Compute bootstrap supports, given reference + bootstrap trees
process ComputeSupports {
    publishDir "results/", mode: 'copy'

    label 'gotree'

    input:
    path reftree
    path boottrees

    output:
    path "${reftree.baseName}_fbp.treefile"

    script:
    """
    gotree compute support fbp -i $reftree -b $boottrees -o ${reftree.baseName}_fbp.treefile -t ${task.cpus}
    """
}

// Full workflow
workflow {
   // This file should be downloaded from GISAID
   gisaid = file("data/sequences_fasta_2023_09_16.tar.xz")
   ba286names = file("data/BA.2.86_all.txt")
   contextnames = file("data/context.txt")
 
   // Download nextclade reference data
   ncref = getRefData()
   // Filter bad quality sequences + split in batches of 100000
   splitgisaid = filterGISAID(gisaid, ba286names, contextnames)
   // Align the sequences against reference using nextaliagn
   galign = alignGISAID(ncref, splitgisaid[0].flatten())
   // Align the BA286 sequences against the reference using nextalign
   ba286align = alignBA286(ncref, splitgisaid[1])
   // Align the Contextual sequences against the reference using nextalign
   contextalign = alignContext(ncref, splitgisaid[2])
   // Annotation of all the sequences with pangolin
   Pangolin(splitgisaid[0].flatten().mix(splitgisaid[1]).mix(splitgisaid[2])).collectFile(name: "gisaid_pangolin.tsv").subscribe{
   it -> it.copyTo("results/")
   }
   // concatenaate GISAID filtered + context 
   db = ConcatDB(galign.collect(), contextalign)
   // Get the closest sequences from ba286 in the database
   closest = GoFasta(db, ba286align)
   // We take the unique closest sequence names
   unique = UniqueClosest(closest)
   // And extract their sequences from the full dataset
   uniquesequences=ExtractSequences(db,unique)
   // We mask some positions of the alignment
   align = mask(uniquesequences.mix(ba286align).mix(contextalign)).collect()
   // We build bootstrap alignments
   bootaligns = BootAligns(align)
   // We infer bootstrap trees
   boottrees  = BootTrees(bootaligns.flatten())
   boottreecoll = boottrees.collectFile(name: 'boottrees.treefile')
   boottreecoll.subscribe{f -> f.copyTo('results/boottrees/')}
   // We infer the reference phylogeny
   reftree = Phylogeny(align)
   // We compute bootstrap supports
   ComputeSupports(reftree, boottreecoll)
}
