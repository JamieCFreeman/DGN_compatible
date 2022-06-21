@fasta_files=glob("*.fas");

for($i=0;$i<@fasta_files;$i++){
	$sequence='';
	open F, "<$fasta_files[$i]" or die "cannot find $fasta_files[$i]\n";
	while(<F>){
		chomp;
		$sequence=$_;
	}
	@sequence=split //, $sequence;
	$seq_length=@sequence;
	$num_rows=int($seq_length/1000);
	$count=0;
	@output=split '.f', $fasta_files[$i];
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
}