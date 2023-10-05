nextflow.enable.dsl=2

process getRefData {

label 'getref'

output:
file "ncref"

script:
"""
nextclade dataset get --name 'sars-cov-2' --output-dir 'ncref'
"""
}

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

workflow {
   query = file("data/gisaid_hcov-19_2023_09_07_07.fasta")
   gisaid = file("data/sequences_fasta_2023_09_16.tar.xz")
   ba286names = file("data/BA.2.86_all.txt")
   contextnames = file("data/context.txt")
   
   ncref = getRefData()
   splitgisaid = filterGISAID(gisaid, ba286names, contextnames)
   galign = alignGISAID(ncref, splitgisaid[0].flatten())
   ba286align = alignBA286(ncref, splitgisaid[1])
   contextalign = alignContext(ncref, splitgisaid[2])

   Pangolin(splitgisaid[0].flatten().mix(splitgisaid[1]).mix(splitgisaid[2])).collectFile(name: "gisaid_pangolin.tsv").subscribe{
   it -> it.copyTo("results/")
   }

   db = ConcatDB(galign.collect(), contextalign)
   closest = GoFasta(db, ba286align)
   unique = UniqueClosest(closest)
   
   uniquesequences=ExtractSequences(db,unique)

   align = mask(uniquesequences.mix(ba286align).mix(contextalign)).collect()
   bootaligns = BootAligns(align)
   boottrees  = BootTrees(bootaligns.flatten())
   boottreecoll = boottrees.collectFile(name: 'boottrees.treefile')
   boottreecoll.subscribe{f -> f.copyTo('results/boottrees/')}
   reftree = Phylogeny(align)

   ComputeSupports(reftree, boottreecoll)
}
