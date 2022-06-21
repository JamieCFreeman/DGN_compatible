#!/usr/bin/perl -w


# Purpose:
#  Take in fasta files and format to Pool lab 1000 bp wrapped .fas1k files.
#
# Edited version of a lab script that globbed all .fas files in a dir to take one file at a time
# 	as a command line argument, which allows for parallelization, and to work with full paths.
#



if(@ARGV != 1){
	die "This script requires one argument in this order: name of fasta file to wrap\n";
}
# Takes one argument: a fasta file
my $fasta_file = $ARGV[0];

$sequence='';
open F, '<', $fasta_file or die "cannot find $fasta_file\n";
while(<F>){
	chomp;
	$sequence=$_;
}
@sequence=split //, $sequence;
$seq_length=@sequence;
$num_rows=int($seq_length/1000);
$count=0;
@output=split '.f', $fasta_file;
$output_file=$output[0].'.fas1k';
open O, ">$output_file";
for($j=0;$j<$num_rows;$j++){
	print O @sequence[$count .. $count+999];
	print O "\n";
	$count+=1000;
}
$get=$seq_length%1000;
$start=$seq_length-$get;
$end=$seq_length-1;
print O @sequence[$start .. $end];
close O;
