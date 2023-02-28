#!/usr/bin/perl
use warnings;
use strict;

#Inputs
#0 = list of fusion genes in format ####TODO
#1 = list of overlapping fusion genes (to be excluded):  gene1\tgene2\n
#2 = manta-produced VCF file
#3 = bam file
#4 = output file
#5 = background output file
#    for fusions connecting any possible combination of genes in the fusion list

open(my $FUSIONS,$ARGV[0]) or die("couldn't read fusions file: " . $ARGV[0] . "\n");
open(my $OVERLAPS,$ARGV[1]) or die("couldn't read overlap file: " . $ARGV[1] . "\n");
open(my $INPUTVCF, "gunzip -c $ARGV[2] |") or die("couldn't read gzipped file: ". $ARGV[2] . "\n");
my $bam = $ARGV[3];
open(my $OUTFILE,">$ARGV[4]") or die("couldn't write to outfile: " . $ARGV[4] . "\n");
open(my $BCKFILE,">$ARGV[5]") or die("couldn't write to outfile: " . $ARGV[5] . "\n");

my %out_lines; #hash outputs to prevent dup lines

#read in all the fusions, store their coords
my @fusion;
my %genes;
my $fcount=0;
while(<$FUSIONS>)
{
  chomp($_);
  my @line = split(/\t/,$_);
  #reading in chr1    3069210 3438621 chr1    2228694 2310119 PRDM16_SKI      .       +       +
  for(my $i=0; $i<scalar(@line); $i++)
  {
    $fusion[$fcount][$i] = $line[$i];
  }

  #store a hash of every gene with coords and fusion partners
  my @gname = split("_",$line[6]);
  #left side
  if(exists($genes{$gname[0]})){
      $genes{$gname[0]}{"partners"} = $genes{$gname[0]}{"partners"} . "," . $gname[1];
  } else {
      $genes{$gname[0]}{"chr"} = $line[0];
      $genes{$gname[0]}{"start"} = $line[1];
      $genes{$gname[0]}{"stop"} = $line[2];
      #include self in partner list for easy filtering below
      $genes{$gname[0]}{"partners"} = $gname[0] . "," . $gname[1];
  }
  
  #right side
  if(exists($genes{$gname[1]})){
      $genes{$gname[1]}{"partners"} = $genes{$gname[1]}{"partners"} . "," . $gname[0];
  } else {
      $genes{$gname[1]}{"chr"} = $line[3];
      $genes{$gname[1]}{"start"} = $line[4];
      $genes{$gname[1]}{"stop"} = $line[5];
      #include self in partner list for easy filtering below
      $genes{$gname[1]}{"partners"} = $gname[1] . "," . $gname[0];
  }
  $fcount++;
}

#read in all the genes that overlap, store their pairings
my %overlaps;
while(<$OVERLAPS>)
{
  chomp($_);
  my @F = split(/\t/,$_);
  $overlaps{$F[0] . "_" . $F[1]} = 1;
  $overlaps{$F[1] . "_" . $F[0]} = 1;
}

#do these coordinates intersect?
sub intersects{
    my ($chrom1, $start1, $stop1, $chrom2, $start2, $stop2) = @_;
    if($chrom1 eq $chrom2){
        if(($start1 <= $stop2 && $start1 >= $start2) || ($start2 <= $stop1 && $start2 >= $start1)){
            return 1;
        }
    }
    return 0;
}

#extract surrounding reads from the bam file, use them to determine whether a site is suspect
sub quality_stats{
    my ($chr1,$pos1,$bam) = @_;
    #use 250 bp window around the site
    my $pos = $chr1 . ":" . ($pos1-250) . "-" . ($pos1+250);

    unless (-e $bam) { #if no bam, can't check qual
        return join("\t",("NA","NA","NA","NA"));
    }

    open(my $BAM, "samtools view $bam $pos | ") or die("couldn't read bam file: ". $ARGV[3] . "\n");
    my %sa_chrs;
    my %mate_chrs;
    while(<$BAM>){
        my $line = $_;
        chomp($line);
        my @F = split("\t",$line);
        #if surrounding reads have supplementary mappings to at least three different chromosomes, toss it
        if($line =~ /SA\:Z\:(.+)/){
            my @sastrings = split(";",$1); #should only have one SA because bam was pre-filtered
            my @safields = split(",",$sastrings[0]);
            $sa_chrs{$safields[0]} = 1;
        }
        #if surrounding reads have mate pairs on at least two other chromosomes, toss it
        if($F[6] ne "=" && $F[6] ne $F[1]){
            $mate_chrs{$F[6]} = 1;
        }
    }
    close($BAM);
    # print STDERR "sa_chrs: " . keys(%sa_chrs) . "  mate_chrs: " . keys(%mate_chrs) . "\n";
    # if(keys(%sa_chrs) > 2 || keys(%mate_chrs) > 1){
    #     return 0;
    # };
    # return 1;
    return join("\t",scalar(keys(%sa_chrs)),scalar(keys(%mate_chrs)));
}

#reorder fusion such that lexically smaller chromosome/numerically smaller position is always first, for deduping purposes
sub order_print{
    my ($OUTPUT_HANDLE, $chr1, $sp1, $chr2, $sp2, $name1, $name2, @rest) = @_;
    my $flip = 0;
    if(($chr1 cmp $chr2) == -1){ #no rearanging needed, just print
        $flip = 0;
    } elsif(($chr1 cmp $chr2) == 1){ #flip
        $flip = 1;
    } else { #same chr, so tiebreak with pos
        if($sp1 > $sp2){
            $flip = 1;
        }
    }

    my $out_line;
    if($flip){
        $out_line = join("\t", ($chr2, $sp2, $chr1, $sp1, $name2, $name1, @rest));
    } else {
        $out_line = join("\t", ($chr1, $sp1, $chr2, $sp2, $name1, $name2, @rest));
    }
    unless(exists($out_lines{$out_line})){ #remove duplicates
        print $OUTPUT_HANDLE $out_line . "\n";
        $out_lines{$out_line} = 1;
    }
}


#read in the vcf, go line by line and see if it matches either our fusion list or
#is a "background" hit that matches two genes in the list in an unexpected way
while(<$INPUTVCF>)
{
    chomp($_);
    next if $_ =~ /^#/; #skip header

    my @line = split(/\s+/,$_);
    my $chr1=$line[0];
    my $st1=$line[1]-1;
    my $sp1=$line[1];

    my $chr2;
    my $st2;
    my $sp2;

    #different sv types require different handling
    if($line[2] =~ /^MantaBND/){
        my @line2 = split(/[][]/,$line[4]);
        my @line3 = split(/\:/,$line2[1]);
        $chr2=$line3[0];
        $st2=$line3[1]-1;
        $sp2=$line3[1];
    } elsif ($line[2] =~ /^MantaDEL/ || $line[2] =~ /^MantaDUP/){
        my @line2 = split(/\;/,$line[7]);
        $chr2=$line[0];
        if($line2[0] =~ /(?:^|;)END=(\d+)/){ #matches start of line or after ;
            $sp2=$1;
        } else {
            die("END coords weren't found in the INFO field for line:" . join("\t",@line) . "\n");
        }
        $st2=$sp2-1;
    #skipping MantaINS, as it's hard to see how a short insertion would
    #create a fusion in a way that's interpretable here
    } elsif ($line[2] =~ /^MantaINS/) {
        next;
    } else {
        die($line[2] . "is not a recognized alteration type\n");
    }

    my $paircount = 0;
    if($line[7] =~ /(?:^|;)PAIR_COUNT=(\d+)/){
        $paircount = $1;
    }

    #loop through each fusion gene, see if it matches
    #yeah, it's inefficient, but it's not bad if the lists are short
    my $count = 0;
    for(my $j=0; $j<@fusion; $j++)
    {
        my $match1 = 0;
        my $match2 = 0;
        my $matchrec1 = 0;
        my $matchrec2 = 0;
        my @fusion_names = split("_",$fusion[$j][6]);

        #do ends match?
        if(intersects($chr1, $st1, $sp1, $fusion[$j][0], $fusion[$j][1], $fusion[$j][2])){
            $match1=1;
        }
        if(intersects($chr1, $st1, $sp1, $fusion[$j][3], $fusion[$j][4], $fusion[$j][5])){
            $matchrec1=1;
        }
        if(intersects($chr2, $st2, $sp2, $fusion[$j][3], $fusion[$j][4], $fusion[$j][5])){
            $match2=1;
        }
        if(intersects($chr2, $st2, $sp2, $fusion[$j][0], $fusion[$j][1], $fusion[$j][2])){
            $matchrec2=1;
        }

        #print STDERR join("\t",($match1,$matchrec1,$match2,$matchrec2,@fusion_names)) . "\n";
        #no match, skip ahead
        if(($match1 + $matchrec1 + $match2 + $matchrec2) < 1){
            next;
        }

        #if genes overlap naturally, remove it
        if(exists($overlaps{join("_",@fusion_names)})){
            next;
        }
        # #check for edge case where fusion genes overlap and at least one breakpoint is in both genes
        # #this is now covered by overlap case above
        # if(($match1 && $matchrec1) || ($match2 && $matchrec2)){
        #     next;
        # }

        #clean match
        if(($match1 && $match2)){
            if($count > 1){
                print STDERR "WARNING: multiple matches for line:\n" . join("\t",@line) . "\n";
            }
            #verify both ends are "clean regions"
            order_print($OUTFILE, $chr1, $sp1, $chr2, $sp2, $fusion_names[0], $fusion_names[1], $paircount, 
                        @line[2..4], $line[7], quality_stats($chr1, $sp1, $bam), quality_stats($chr2, $sp2, $bam));
            $count++;
        } elsif ($matchrec1 && $matchrec2){
            if($count > 1){
                print STDERR "WARNING: multiple matches for line:\n" . join("\t",@line) . "\n";
            }
            #verify both ends are "clean" regions
            order_print($OUTFILE, $chr1, $sp1, $chr2, $sp2, $fusion_names[1], $fusion_names[0], $paircount, 
                        @line[2..4], $line[7], quality_stats($chr1, $sp1, $bam), quality_stats($chr2, $sp2, $bam));
            $count++;
        }
        
        ### now do background matching for this region
        sub checkGenes {
            my ($chrom, $start, $stop, $othergene) = @_;
            my @matching_genes;
            #loop through every gene
            foreach my $genename (keys %genes){
                #does this gene match the hanging end?
                if(intersects($chrom, $start, $stop, $genes{$genename}{"chr"}, $genes{$genename}{"start"}, $genes{$genename}{"stop"})){
                    # #check for case where the genes overlap, remove these
                    if(exists($overlaps{$othergene . "_" . $genename})){
                        next;
                    }
                    my @partners = split(",",$genes{$genename}{"partners"});
                    #don't add it if we already know about this fusion pairing (it's in our list)
                    #print STDERR "found a match to " . $othergene . " - " . $genename . "\n";
                    #print STDERR join("\t",@partners) . "\n";
                    unless($othergene ~~ @partners){
                        push(@matching_genes, $genename);
                    }
                }
            }
            #print STDERR "matching genes: " . join("\t",@matching_genes) . "\n";
            
            if(@matching_genes > 0){
                return(join(",",@matching_genes))
            }
            return 0;
        }

        
        if($match1){ #find match for second end
            my $matches = checkGenes($chr2, $st2, $sp2, $fusion_names[0]);
            if($matches){
                #print OUTFILE "match1\n";
                order_print($BCKFILE, $chr1, $sp1, $chr2, $sp2, $fusion_names[0], $matches, $paircount, 
                            @line[2..4], $line[7], quality_stats($chr1, $sp1, $bam), quality_stats($chr2, $sp2, $bam));
                #print OUTFILE join("\t",(@line[0..7],$fusion_names[0],$matches,$paircount)) . "\n";
                next;
            }
        }
        if($match2){ #find match for first end
            my $matches = checkGenes($chr1, $st1, $sp1, $fusion_names[1]);
            if($matches){
                #print OUTFILE "match2\n";
                order_print($BCKFILE,$chr1, $sp1, $chr2, $sp2, $matches, $fusion_names[1], $paircount, 
                            @line[2..4], $line[7], quality_stats($chr1, $sp1, $bam), quality_stats($chr2, $sp2, $bam));
                #print OUTFILE join("\t",(@line[0..7],$matches,$fusion_names[1],$paircount)) . "\n";
                next;
            }
        }
        
        if($matchrec1){ #find rec match for second end
            my $matches = checkGenes($chr2, $st2, $sp2, $fusion_names[1]);
            if($matches){
                #print OUTFILE "match3\n";
                order_print($BCKFILE, $chr1, $sp1, $chr2, $sp2, $fusion_names[1], $matches, $paircount, 
                            @line[2..4], $line[7],quality_stats($chr1, $sp1, $bam), quality_stats($chr2, $sp2, $bam));
                #print OUTFILE join("\t",(@line[0..7],$matches,$fusion_names[1],$paircount)) . "\n";
                next;
            }
        }
        if($matchrec2){ #find rec match for first end
            my $matches = checkGenes($chr1, $st1, $sp1, $fusion_names[0]);
            if($matches){
                #print OUTFILE "match4\n";
                order_print($BCKFILE, $chr1, $sp1, $chr2, $sp2, $matches, $fusion_names[0], $paircount, 
                            @line[2..4], $line[7],quality_stats($chr1, $sp1, $bam), quality_stats($chr2, $sp2, $bam));
                #print OUTFILE join("\t",(@line[0..7],$fusion_names[0],$matches,$paircount)) . "\n";
                next;
            }
        }
    }
}
