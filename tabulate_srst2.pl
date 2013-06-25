#!/usr/bin/perl


# e.g. shigella-samples.csv
open(TRUE_ST,$ARGV[0]);

while(<TRUE_ST>)
{
	chomp;
	
	my@foo=split(/\,/);
	$samplename=$foo[0];
	$this_st=$foo[1];
	if($this_st eq "" || $this_st eq "?" || $this_st eq "unknown" ||
	$this_st =~ /SLV/ || $this_st =~ /DLV/)
	{
		$this_st="NA";
	}
	$st_true{$samplename}=$this_st;
}

close(TRUE_ST);

# e.g. ../data/pub_mlst/Shigella-sonnei/ecoli_ST_table.txt
open(ST_FILE,$ARGV[1]);
$count=0;
while(<ST_FILE>)
{
	chomp;
	
	if($count==0)
	{
		@gene_order=split(/\t/,uc($_));
		shift(@gene_order);
	}
	else
	{
		my@foo=split(/\t/);
		$st_map{$foo[1]}{$foo[2]}{$foo[3]}{$foo[4]}{$foo[5]}{$foo[6]}{$foo[7]}=$foo[0];
		$st_map2{$foo[0]}="$foo[1]\t$foo[2]\t$foo[3]\t$foo[4]\t$foo[5]\t$foo[6]\t$foo[7]";
	}
	
	$count++;
}

close(ST_FILE);



open(LIST,$ARGV[2]);

while(<LIST>)
{
	chomp;

	$file=$_;

	open(FILE,$file);

	if($file =~ /\//)
	{
		my@goo=split(/\//,$file);
		$file=$goo[scalar(@goo)-1];
	}

	my@foo=split(/\./,$file);
	$sample=$foo[0];
	$allele=uc($foo[1]);

	$hash_depth{$sample}{$allele}=$avg_depth;
	$count=0;

	while(<FILE>)
	{
		chomp;
		if($count>0)
		{
			$score=$avg_depth=$depth_a=$depth_z="NA";
			($allelenum,$score,$avg_depth,$depth_a,$depth_z)=split(/\t/);

			$depth_a="NA" if ($depth_a eq "");
			$depth_z="NA" if ($depth_z eq "");

			$hash_depth{$sample}{$allele}{$allelenum}=$avg_depth;
			$hash_depth_edge{$sample}{$allele}{$allelenum}="$depth_a\t$depth_z";
			$hash{$sample}{$allele}{$allelenum}=$score;
		}
		$count++;
	}

	close(FILE);i
	
}

close(LIST);

#header
print "Sample\tGene\tAllele_SRST2\tAllele_True\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\tConcordant_alleles\n";

TWO: for $sample (keys %hash)
{
	%this_true_hash=();
	$this_true_st=$st_true{$sample};
	@this_true_alleles=split(/\t/,$st_map2{$this_true_st});

	for($i=0;$i<scalar(@gene_order);$i++)
	{
		if(defined $st_true{$sample})
		{
			$this_true_hash{$gene_order[$i]}=$this_true_alleles[$i];
		}
		else
		{
			print STDERR "No 'true' ST for $sample!\n";
			$this_true_hash{$gene_order[$i]}="NA";
		}
	}

	for $allele (keys %{ $hash{$sample} })
	{
		$counter=0;

		ONE: for $allelenum ( sort {$hash{$sample}{$allele}{$a} <=> $hash{$sample}{$allele}{$b}} keys %{ $hash{$sample}{$allele} } )
		{
			if($counter==0)
			{
				$best_score=$hash{$sample}{$allele}{$allelenum};
			}
		
			$score=$hash{$sample}{$allele}{$allelenum};
			$avg_depth=$hash_depth{$sample}{$allele}{$allelenum};
			$edge_depth=$hash_depth_edge{$sample}{$allele}{$allelenum};

			if($counter>0 && $best_score<$score)
			{
				next ONE;
			}

			if(defined $this_true_hash{$allele})
			{
				($a,$b)=split(/\-/,$allelenum);
				if($b eq $this_true_hash{$allele})
				{
					$concord=1;
				}
				elsif($this_true_hash{$allele} eq "NA")
				{
					$concord="NA";
				}
				else
				{
					$concord=0;
				}

				print "$sample\t$allele\t$allelenum\t$this_true_hash{$allele}\t$score\t$avg_depth\t$edge_depth\t$concord\n";
			}
			else
			{
				print "$sample\t$allele\t$allelenum\tNA\t$score\t$avg_depth\t$edge_depth\tNA\n";
			}

		$counter++;

		}
	}
}


