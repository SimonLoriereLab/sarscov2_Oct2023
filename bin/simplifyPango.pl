#!/usr/bin/perl -w
use strict;
use Fatal;


my %annots;
my %pango;
my %counts;

my $meta = $ARGV[0]; # meta data file 
my $lineage = $ARGV[1]; # lineage file from https://github.com/cov-lineages/pango-designation/blob/e6a71b06d7b2e0649f222c52af4bb38bd7af3f3c/lineage_notes.txt
my $head = "";
my %alias; # lineage aliases
my %reversealias;
my %simplify;
my %full;

sub shrink{
    my ($lin, $alias) = @_;

    my @aliascols = split(/\./,$alias);
    my $last = 4;
    if($#aliascols<4){
	$last = $#aliascols;
    }
    
    if($lin =~ /^BQ\.1$/ || $lin =~ /^BQ\.1\./){
	$last = 10;
    }
    if($lin =~ /^BQ\.1\.1$/ || $lin =~ /^BQ\.1\.1\./){
	$last = 11;
    }
    if($lin =~ /^FD\.1\.1$/ || $lin =~ /^FD\.1\.1\./){
	$last = 5;
    }
    if($lin =~ /^BA\.2\.86$/ || $lin =~ /^BA\.2\.86\./ || $alias =~ /^B\.1\.1\.529\.2\.86$/ || $alias =~ /^B\.1\.1\.529\.2\.86\./){
	print STDERR  $lin." ".$alias."\n";
	$last = 5;
    }
    if($lin =~ /^JN\.1$/ || $lin =~ /^JN\.1\./ || $alias =~ /^B\.1\.1\.529\.2\.86\.1\.1$/ || $alias =~ /^B\.1\.1\.529\.2\.86\.1\.1\./){
	print STDERR  $lin." ".$alias."\n";
	$last = 7;
    }
    if($lin =~ /^B\.1\.617\.2$/ || $lin =~ /^B\.1\.617\.2\./ || $lin=~/^AY$/ || $lin=~/^AY\./){
	$last = 3;
    }
    if($lin =~ /^XBB\.[^.]+\./ || $alias =~ /^XBB\.[^.]+\./){
	$last = 2;
    }
    if($lin =~ /^XBB\.1\.22\.1$/ || $lin =~ /^XBB\.1\.22\.1\./ || $lin=~/^FY/){
	$last = 3;
    }
    if($lin =~ /^FD\.1\.1$/ || $lin =~ /^FD\.1\.1\./){
	$last = 5;
    }
    if($lin =~ /^FE\.1$/ || $lin =~ /^FE\.1\./){
	$last = 4;
    }
    if($lin =~ /^FL\.1$/ || $lin =~ /^FL\.1\./){
	$last = 4;
    }
    if($lin =~ /^FK\.1$/ || $lin =~ /^FK\.1\./){
	$last = 13;
    }
    if($lin =~ /^EG\.5$/ || $lin =~/^EG\.5\./){
	$last = 4;
    }
    if($lin =~ /^EG\.5\.1$/ || $lin =~ /^EG\.5\.1\./){
	$last = 5;
    }
    if($lin =~ /^BN\.1$/ || $lin =~ /^BN\.1\./){
	$last = 7;
    }
    if($lin =~ /^CH\.1$/ || $lin =~ /^CH\.1\./){
	$last = 10;
    }
    if($lin =~ /^BA\.2\.75$/ || $lin =~ /^BA\.2\.75\./){
	$last = 5;
    }
    if($lin =~ /^BA\.2\.12\.1$/ || $lin =~ /^BA\.2\.12\.1\./ || $lin =~/^BG\./){
	$last = 6;
    }

    if($last>$#aliascols){
	$last = $#aliascols;
    }

    my $simpl = join(".",@aliascols[0..$last]);
    
    return $simpl;
}


open(LIN,$lineage);
while(<LIN>){
    chomp;
    my @cols = split(/\t/);
    my $lin = $cols[0];
    $lin=~s/^\*//;
    if($cols[1] =~ /Alias of ([^,\s]+)/i){
	my $al = $1;
	$alias{$lin} = $al;
	$reversealias{$al}=$lin;
	my $simpl = shrink($lin,$al);
	$simplify{$al}=$simpl;
	$simplify{$lin} = $simpl;
    }elsif($cols[1] =~ /^(XBB[^,\s]+)/i){
	my $al = $1;
	$alias{$lin} = $al;
	$reversealias{$al}=$lin;
	my $simpl = shrink($lin,$al);
	$simplify{$al}=$simpl;
	$simplify{$lin} = $simpl;	
    }else{
	$alias{$lin} = $lin;
	$reversealias{$lin}=$lin;
	my $last = 4;

	my $simpl = shrink($lin,$lin);
	$simplify{$lin} = $simpl;
    }
}
close(LIN);

open(IN,$meta);
$_=<IN>;
chomp;
$head = $_;
while(<IN>){
    chomp;
    my @cols=split(/\t/);
    my $name = $cols[0];

    if($#cols > 3){
	my $pango = $cols[4];
	
	$pango=~s/^B\.1\.88\.1$/B\.1\.88/;

	if(length($pango)>0 && $pango ne "Unassigned"){
	    my $annot = $simplify{$pango};
	    my $fullannot = $annot;
	    if(defined $reversealias{$annot}){
		$annot = $reversealias{$annot};
	    }
	    $pango{$name} = $annot;
	    $full{$name} = $fullannot;
	    $counts{$annot}+=1;
	    $annots{$name}=join("\t",@cols);
	}
    }
}
close(IN);

print $head."\tshort_pango\tfull_pango\n";
for my $name (keys(%pango)){
    my $p = $pango{$name};
    my $c = $counts{$p};
    my $a = $annots{$name};
    my $f = $full{$name};

    print $a."\t".$p."\t".$f."\n";
}
