#!/usr/bin/perl -w

use strict;

print "TREE_COLORS\n";

print "SEPARATOR SPACE\n";

print "DATASET_LABEL short_pango_tips\n";
print "DATA\n";

my $start = 0;
while(<>){

    if(/^DATA$/){
	$start = 1;
    }elsif($start eq 1){
	my @cols = split(/\t/);
	my $name = $cols[0];
	my $annot = $cols[2];
	my @annots = split(/\|/,$annot);
	my $color = $annots[3];
	print "$name branch $color normal 1\n";	
    }
}
#915|777 branch #00ff00 dashed 5

#NODE_ID TYPE COLOR LABEL_OR_STYLE SIZE_FACTOR

#Examples
#internal node with solid branches colored blue and twice the standard width
#9031|9606 clade #0000ff normal 2
#internal node with dashed branches colored red and one half the standard width
#601|340 clade #ff0000 dashed 0.5
#a single internal branch colored green, dashed and 5 times the normal width

#colored range covering all leaves of an internal node,  colored red and with label 'Eukaryota'
#184922|9606 range #ff0000 Eukaryota
#examples of colored ranges from iTOL's Tree of Life
#2190|2287 range #aaffaa Archaea
#623|1502 range #aaaaff Bacteria

#leaf label for node 9606 will be displayed in green, bold and twice the regular font size
#9606 label #00ff00 bold 2

#leaf label for node 9031 will be displayed in yellow, bold italic and half the regular font size
#9031 label #ffff00 bold-italic 0.5

#leaf label for node 8015 will be displayed in blue
#8015 label #0000ff

#leaf label for node 9606 will have a semi-transparent red background
#9606 label_background rgba(255,0,0,0.5)
