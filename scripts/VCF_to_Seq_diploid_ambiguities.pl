#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#This script produces unformatted sequence files for each chromosome arm in the reference genome. The shifted all-sites vcf files and INDEL files must be in the same directory as this file

#INPUT

my @Chromosomes = ('Yhet' , 'mtDNA' , 'Chr2L' , 'ChrX' , 'Chr3L' , 'Chr4', 'Chr2R' , 'Chr3R');
# , 'Uextra' , '2RHet' , '2LHet' , '3LHet' , '3RHet' , 'U' , 'XHet' , 'Wolbachia');
my @Samples = ();
my @Files = ();
# JCF my @AllFiles = ();
my $cmd = '';

my $k = 0;

########################################################################
# JCF edits: remove windows line breaks (show as ^M char in vim) from original file
#
# %s/^M//g
#
# Where ^M is entered by holding down Ctrl and typing v followed by m, and then releasing Ctrl. This is sometimes abbreviated as ^V^M, but note that you must enter it as described in the previous sentence, rather than typing it out literally.
#
# Don't grab files from dir, take from command line
#open directory
#opendir DIR, "." or die "couldn't open directory\n";
#@AllFiles = readdir(DIR);
#closedir DIR;

my @AllFiles = (@ARGV);

########################################################################


#find correct input files
for ($a = 0; $a < @AllFiles; $a++){
  if ($AllFiles[$a] =~ m'_shifted.vcf'){
    push (@Files, (split '_shifted', $AllFiles[$a])[0]);
  }
}

#@Files=('SD88N');

#die "@Files\n";

for ($k = 0; $k < @Files; $k++) {

  my $RefChrom = 1;
  my $Contigs = 8;   #Set this to the number of chromosome arms in your reference genome (i.e., if you have euchromatin from 2R, 2L, 3R, 3L, and X then this would be 5

  my $Name = $Files[$k];
    my $VCFin = $Files[$k] . '_shifted.vcf';
    #my $VCFin = $Files[$k] . '_sites.vcf';
  my $INDELSin = $Files[$k] . '_INDELS.vcf';
    $cmd = "gunzip " . $VCFin . '.gz';
    system($cmd);

  for ($RefChrom = 1; $RefChrom <= $Contigs; $RefChrom++) {
    my $outfile = $Name . '_' . $Chromosomes[$RefChrom - 1] . '_diploid.fas';
    #$cmd = 'purge';
    #system($cmd);

    my @line = ();
    my $ChromX = 0;
    my $PositionX = 0;
    my $RefX = 0;
    my $AltX = 0;
    my $QualX = 0;
    my $IndividualX = 0;
    my $Genotype = 0;
    my @IndelChromX = ();
    my @IndelPositionX = ();
    my @IndelRefX = ();
    my @IndelAltX = ();
    my $Chrom2 = 0;
    my $Position2 = 0;
    my $Ref2 = 0;
    my $Alt2 = 0;
    my $Qual2 = 0;
    my $Individual2 = 0;
    my $AlleleAlt_Depth2 = 0;
    my $AlleleRef_Depth2 = 0;
    my $Genotype2 = 0;
    my $gap = 0;
    my $offset = 0;

    open C, ">$outfile";

    my $m = 0;
    my $a = 0;
    my $prevposX = 0;

    #open X, "</Volumes/4_TB_RAID_Set_1/dgrp/round1/indel_vcfs/$INDELSin";
    open X, "<$INDELSin";
    while (<X>){
      chomp;
      last if m/^$/;
      @line = split;
      next if ($line[0] =~ m/#/);
      last if ($line[0] > $RefChrom);
      if ($line[0] == $RefChrom){
        push @IndelChromX, $line[0];
        push @IndelPositionX, $line[1];
        push @IndelRefX, length$line[3];
        push @IndelAltX, length$line[4];
      }
    }
    close X;
	
    undef @line;

    open Y, "<$VCFin";
    while (<Y>){
      chomp;
      last if m/^$/;
      @line = split;
      next if ($line[0] =~ m/#/);
      last if ($line[0] > $RefChrom);
      if ($line[0] == $RefChrom){
      if ($line[4] eq '.') {
        $ChromX = $line[0];
        $PositionX = $line[1];
        $RefX = $line[3];
        $AltX = $line[4];
        $QualX = $line[5];
        $IndividualX = "0,0";
        $Genotype = '0/0';
      }
      else {
        $ChromX = $line[0];
        $PositionX = $line[1];
        $RefX = $line[3];
        $AltX = $line[4];
        $QualX = $line[5];
        $IndividualX = ((split /:/, $line[9])[1]);
        $Genotype = ((split /:/, $line[9])[0]);
      }

      if ($a < @IndelChromX) {
        if ($PositionX < ($IndelPositionX[$a] - 2)) {
        $Chrom2 = $ChromX;
        $Position2 = $PositionX;
        $Ref2 = $RefX;
        $Alt2 = $AltX;
        $Qual2 = $QualX;
        $AlleleRef_Depth2 = ((split /,/, $IndividualX)[0]);
        $AlleleAlt_Depth2 = ((split /,/, $IndividualX)[1]);
        $Genotype2 = $Genotype;
      }
      elsif (($PositionX >= ($IndelPositionX[$a] - 2)) && ($IndelRefX[$a] == 1)) {
        if ($PositionX <= ($IndelPositionX[$a] + 3)) {
        }
        else {
          $Chrom2 = $ChromX;
          $Position2 = $PositionX;
          $Ref2 = $RefX;
          $Alt2 = $AltX;
          $Qual2 = $QualX;
          $AlleleRef_Depth2 = ((split /,/, $IndividualX)[0]);
          $AlleleAlt_Depth2 = ((split /,/, $IndividualX)[1]);
          $Genotype2 = $Genotype;
          $a++;
        }
      }
      elsif (($PositionX >= ($IndelPositionX[$a] - 2)) && ($IndelAltX[$a] == 1)) {
        $offset = $IndelRefX[$a];
        if ($PositionX <= ($offset + 2 + $IndelPositionX[$a])) {
        }
        else {
          $Chrom2 = $ChromX;
          $Position2 = $PositionX;
          $Ref2 = $RefX;
          $Alt2 = $AltX;
          $Qual2 = $QualX;
          $AlleleRef_Depth2 = ((split /,/, $IndividualX)[0]);
          $AlleleAlt_Depth2 = ((split /,/, $IndividualX)[1]);
          $Genotype2 = $Genotype;
          $a++;
        }
      }
    }
    else {
      $Chrom2 = $ChromX;
      $Position2 = $PositionX;
      $Ref2 = $RefX;
      $Alt2 = $AltX;
      $Qual2 = $QualX;
      $AlleleRef_Depth2 = ((split /,/, $IndividualX)[0]);
      $AlleleAlt_Depth2 = ((split /,/, $IndividualX)[1]);
      $Genotype2 = $Genotype;
    }
    #Now, lets write sequences
    if ($m == 0) {
      if (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 eq '.')){
        print C "$Ref2";
      }
      elsif (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 ne '.') && ($Genotype2 eq'1/1')){
        print C "$Alt2";
      }
      elsif (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 ne '.') && ($Genotype2 eq'0/1')){
        if ($Ref2 eq 'C') {
          if ($Alt2 eq 'A') {
            print C "M";
          }
          elsif ($Alt2 eq 'G') {
            print C "S";
          }
          elsif ($Alt2 eq 'T') {
            print C "Y";
          }
        }
        elsif ($Ref2 eq 'G') {
          if ($Alt2 eq 'A') {
            print C "R";
          }
          elsif ($Alt2 eq 'C') {
            print C "S";
          }
          elsif ($Alt2 eq 'T') {
            print C "K";
          }
        }
        elsif ($Ref2 eq 'T') {
          if ($Alt2 eq 'A') {
            print C "W";
          }
          elsif ($Alt2 eq 'G') {
            print C "K";
          }
          elsif ($Alt2 eq 'C') {
            print C "Y";
          }
        }
        elsif ($Ref2 eq 'A') {
          if ($Alt2 eq 'C') {
            print C "M";
          }
          elsif ($Alt2 eq 'G') {
            print C "R";
          }
          elsif ($Alt2 eq 'T') {
            print C "W";
          }
        }
      }
      else {
        print C "N";
      }
    }
    else {
      if ($Position2 == $prevposX) {
    }
    elsif (($Position2 - $prevposX) > 1) {
      my $gap = 1;
      while ($gap < ($Position2 - $prevposX)) {
        print C "N";
        $gap++;
      }
      if (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 eq '.')){
        print C "$Ref2";
      }
      elsif (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 ne '.') && ($Genotype2 eq'1/1')){
        print C "$Alt2";
      }
      elsif (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 ne '.') && ($Genotype2 eq'0/1')){
        if ($Ref2 eq 'C') {
          if ($Alt2 eq 'A') {
            print C "M";
          }
          elsif ($Alt2 eq 'G') {
            print C "S";
          }
          elsif ($Alt2 eq 'T') {
            print C "Y";
          }
        }
        elsif ($Ref2 eq 'G') {
          if ($Alt2 eq 'A') {
            print C "R";
          }
          elsif ($Alt2 eq 'C') {
            print C "S";
          }
          elsif ($Alt2 eq 'T') {
            print C "K";
          }
        }
        elsif ($Ref2 eq 'T') {
          if ($Alt2 eq 'A') {
            print C "W";
          }
          elsif ($Alt2 eq 'G') {
            print C "K";
          }
          elsif ($Alt2 eq 'C') {
            print C "Y";
          }
        }
        elsif ($Ref2 eq 'A') {
          if ($Alt2 eq 'C') {
            print C "M";
          }
          elsif ($Alt2 eq 'G') {
            print C "R";
          }
          elsif ($Alt2 eq 'T') {
            print C "W";
          }
        }
      }
      else {
        print C "N";
      }
    }
    else {
       if (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 eq '.')){
         print C "$Ref2";
       }
       elsif (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 ne '.') && ($Genotype2 eq'1/1')){
         print C "$Alt2";
       }
       elsif (($Qual2 ne '.') && ($Qual2 >= 32) && ($Alt2 ne '.') && ($Genotype2 eq'0/1')){
         if ($Ref2 eq 'C') {
           if ($Alt2 eq 'A') {
             print C "M";
           }
           elsif ($Alt2 eq 'G') {
             print C "S";
           }
           elsif ($Alt2 eq 'T') {
             print C "Y";
           }
        }
        elsif ($Ref2 eq 'G') {
          if ($Alt2 eq 'A') {
            print C "R";
          }
          elsif ($Alt2 eq 'C') {
            print C "S";
          }
          elsif ($Alt2 eq 'T') {
            print C "K";
          }
        }
        elsif ($Ref2 eq 'T') {
          if ($Alt2 eq 'A') {
            print C "W";
          }
          elsif ($Alt2 eq 'G') {
            print C "K";
          }
          elsif ($Alt2 eq 'C') {
            print C "Y";
          }
        }
        elsif ($Ref2 eq 'A') {
          if ($Alt2 eq 'C') {
            print C "M";
          }
          elsif ($Alt2 eq 'G') {
            print C "R";
          }
          elsif ($Alt2 eq 'T') {
            print C "W";
          }
        }
      }
      else {
        print C "N";
      }
    }
  }
  $prevposX = $Position2;
  $m++;
}
}
  }
    $cmd = "gzip " . $VCFin;
    system($cmd);
}
