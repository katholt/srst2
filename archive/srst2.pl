#!/usr/bin/perl


#
# SRST2 - Short Read Sequence Typer (v2)
# 
# Author - Michael Inouye (minouye@unimelb.edu.au)
#
# Dependencies:
#    bowtie2       http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
#    SAMtools      http://samtools.sourceforge.net
#    R             http://www.r-project.org 
#



# Path for the bowtie2 suite
# $path_bowtie2="/software/bowtie2-2.1.0";
$path_bowtie2="";

# Path for SAMtools
#$path_samtools="/software/samtools-0.1.18";
$path_samtools="";

# Defaults for options
$in_file="-1";
$out_file="-1";
$index_file="-1";
$index_file_list="-1";
$verbose=0;
$se=1;
$in_file_type="-q";
$qual_map=1;   
$qual_base=20;
$edge_a=$edge_z=2;

$prob_err=0.01;

# Parse the input arguments, use '-h' for help
for($i=0;$i<scalar(@ARGV);$i++)
{
	if($ARGV[$i] eq "-i" ||
	$ARGV[$i] eq "--input")
	{
		shift(@ARGV);
		$in_file=$ARGV[$i];
	}
	elsif($ARGV[$i] eq "-o" ||
	$ARGV[$i] eq "--output")
	{
		shift(@ARGV);
		$out_file=$ARGV[$i];
	}
	elsif($ARGV[$i] eq "-b" ||
        $ARGV[$i] eq "--build")
        {
                shift(@ARGV);
                $index_file=$ARGV[$i];
        }
	elsif($ARGV[$i] eq "-bl" ||
        $ARGV[$i] eq "--build-list")
        {
                shift(@ARGV);
                $index_file_list=$ARGV[$i];
        }
	elsif($ARGV[$i] eq "-pe" ||
        $ARGV[$i] eq "--paired-end")
        {
		# Input is paired-end reads, not single-end
                $se=0;
        }
	elsif($ARGV[$i] eq "-v" ||
        $ARGV[$i] eq "--verbose")
        {
                $verbose=1;
        }
	elsif($ARGV[$i] eq "-t" ||
        $ARGV[$i] eq "--input-type")
        {
                shift(@ARGV);
                $in_file_type=$ARGV[$i];
        }
	elsif($ARGV[$i] eq "-z" ||
        $ARGV[$i] eq "--other")
        {
                shift(@ARGV);
                $other_bowtie2_options=$ARGV[$i];
        }
	elsif($ARGV[$i] eq "-mapq")
        {
                shift(@ARGV);
                $qual_map=$ARGV[$i];
        }
	elsif($ARGV[$i] eq "-baseq")
        {
                shift(@ARGV);
                $qual_base=$ARGV[$i];
        }
	elsif($ARGV[$i] eq "-h" ||
	$ARGV[$i] eq "--help")
	{
		print STDERR "\n";
		
		print STDERR "Usage:\n";
		print STDERR "\tsrst2.pl [other-arguments] -i [input-file] -o [output-file]\n";
		print STDERR "\n";
		print STDERR "Options:\n";
		print STDERR "\t-i\tInput file ([file1],[file2] if paired-end)\n";
		print STDERR "\t-t\tInput file type (-q or -qseq or -f)\n";
		print STDERR "\t-b\tFasta file to build bowtie2 index\n";
		print STDERR "\t-bl\tList of fasta files for building bowtie2 indicies\n";
		print STDERR "\t-pe\tReads are paired-end (default: single-end)\n";
		print STDERR "\t-z\tOther options for bowtie2\n";
		print STDERR "\t-baseq\tFilter bases less than this PHRED quality (default: 20)\n";
		print STDERR "\t-o\tOutput file\n";
		print STDERR "\n";

		exit(1);
	}
}

print STDERR "\n";

# Check arguments
if($in_file ne "-1")
{
	if($se==1 && -e $in_file)
	{
		print STDERR "\tInput file\t\t$in_file\n";
	}
	elsif($se==0)
	{
		($a,$b)=split(/\,/,$in_file);
		
		if(-e $a && -e $b)
		{
			print STDERR "\tInput file1\t\t$a\n";
			print STDERR "\tInput file2\t\t$b\n";
		}
		else
		{
			print STDERR "Input file has not been correctly defined.\n";
                	exit(1);
		}
	}
	else
	{
		print STDERR "Input file has not been correctly defined.\n";
        	exit(1);
	}
}
else
{
	print STDERR "Input file has not been correctly defined.\n";
	exit(1);
}

if($out_file ne "-1")
{
	print STDERR "\tOutput file prefix\t$out_file\n";
}
else
{
	print STDERR "Output file prefix (-o) must be defined!";
	exit(1);
}

if($in_file_type eq "-q" ||
   $in_file_type eq "-qseq" ||
   $in_file_type eq "-f")
{
        print STDERR "\tInput file type\t\t$in_file_type\n";
}
else
{
	print STDERR "Type of input file has not been correctly defined.\n";
        exit(1);
}

if($index_file ne "-1")
{
        print STDERR "\tIndex file\t\t$index_file\n";
}
elsif($index_file_list ne "-1")
{
	print STDERR "\tIndex list\t\t$index_file_list\n";
}
else
{
	print STDERR "Index file or list has not been correctly defined.\n";
        exit(1);
}

print STDERR "\n";



### Build a bowtie2 index from the given input

if($index_file ne "-1")
{
	$index_file[0]=$index_file;

	# Start building the command
	#$command="$path_bowtie2/bowtie2-build";
	$command="bowtie2-build";
	
	$command .= " $index_file";
	$command .= " $index_file";

	$built_index=$index_file . ".1.bt2";
	if(-e $built_index)
	{
		print STDERR "Index for $index_file is already built...\n";
	}
	else
	{
		print STDERR "Building bowtie2 index for $index_file...\n";
		`$command`;    # Run the command
	}
}
else
{
	# User has provide a list of indices to build, parse that list
	open(LIST_INDEX,"$index_file_list");	
	$counter=0;

	while(<LIST_INDEX>)
	{
		chomp;
	
		$index_file[$counter]=$_;
		
		$command="$path_bowtie2/bowtie2-build";
        
	        $command .= " $index_file[$counter]";
        	$command .= " $index_file[$counter]";

		$built_index=$index_file[$counter] . ".1.bt2";
        	
		if(-e $built_index)
        	{
                	print STDERR "Index for $index_file[$counter] is already built...\n";
        	}
        	else
        	{
                	print STDERR "Building bowtie2 index for $index_file[$counter]...\n";
                	`$command`;   # Run the command
        	}

		$counter++;
	}

	close(LIST_INDEX);

}


#############################################################

### Run bowtie2 on the newly built index/indices

#############################################################

# Create output prefix
$out_file_orig=$out_file;

for($counter=0;$counter<scalar(@index_file);$counter++)
{

if($index_file[$counter] =~ /\//)
{
	my@goo=split(/\//,$index_file[$counter]);
	$suffix=$goo[scalar(@goo)-1];
}
else
{
	$suffix=$index_file[$counter];
}

$out_file=$out_file_orig . ".$suffix.srst2";

print STDERR "Output prefix set to $out_file\n";


# Start a command and add bowtie2 parameters
$command="$path_bowtie2/bowtie2";

if($se==1)
{
	$command .= " -U $in_file";
}
else
{
	($a,$b)=split(/\,/,$in_file);
	$command .= " -1 $a -2 $b";
}

if($out_file ne "-1")
{
	$command .= " -S $out_file";
} 
else
{
	print STDERR "Output file prefix needs to be defined!\n";
	exit(1);
}

$command .= " $in_file_type";
$command .= " --very-sensitive-local";
$command .= " --no-unal";

# Search for and report all alignments
$command .= " -a";

# The index to be aligned to
$command .= " -x $index_file[$counter]";

$command .= " $other_bowtie2_options";



print STDERR "Aligning reads to index $index_file[$counter] using bowtie2...\n";
`$command`;    # Run bowtie2 alignment



### Modify Bowtie's SAM formatted output so that we get secondary alignments in downstream pileup
open(SAM,"$out_file");
open(SAM2,">$out_file.mod");

while(<SAM>)
{
	chomp;

	if($_ !~ /^\@/)
	{
		my@foo=split(/\t/,$_,3);
		#my@foo=split(/\t/,$_,7);

		$flag=$foo[1];
		#$cigar=$foo[5];
		#$start=$foo[3];

		# Bitwise edit flag for secondary alignments
		#$flag = $flag ^ 256;
		$flag = $flag - 256 if ($flag & 256);

		#$cigar =~ s/S/M/g;

		print SAM2 "$foo[0]\t$flag\t$foo[2]\n";
		#print SAM2 "$foo[0]\t$flag\t$foo[2]\t$foo[3]\t$foo[4]\t$cigar\t$foo[6]\n";
	}
	else
	{
		print SAM2 "$_\n";
	}
}

close(SAM);
close(SAM2);


### Analyse output with SAMtools
print STDERR "Processing Bowtie2 output with SAMtools...\n";
print STDERR "Generate and sort BAM file...\n";
$out_file_sam1=$out_file . ".bam";
#$command="$path_samtools/samtools view -b -o $out_file_sam1 -q $qual_map -S $out_file.mod";
$command="samtools view -b -o $out_file_sam1 -q $qual_map -S $out_file.mod";
`$command`;

$out_file_sam2=$out_file . ".sorted";
#$command="$path_samtools/samtools sort $out_file_sam1 $out_file_sam2";
$command="samtools sort $out_file_sam1 $out_file_sam2";
`$command`;

#$command="$path_samtools/samtools faidx $index_file[$counter]";
$command="samtools faidx $index_file[$counter]";
`$command`;


print STDERR "Generate pileup...\n";
$out_file_sam3=$out_file . ".pileup";
$faidx_file="$index_file[$counter]";
#$command="$path_samtools/samtools mpileup -L 1000 -f $faidx_file -Q $qual_base -q $qual_map $out_file_sam2.bam > $out_file_sam3";
$command="samtools mpileup -L 1000 -f $faidx_file -Q $qual_base -q $qual_map $out_file_sam2.bam > $out_file_sam3";
`$command`;



### Process SAMtools output
print STDERR "Processing SAMtools output...\n";

# Get sequence lengths for reference alleles - important for scoring
open(INDEX_FAI,"$index_file[$counter].fai");
$count=0;
while(<INDEX_FAI>)
{
	chomp;

	my@foo=split(/\t/);
	$size{$foo[0]}=$foo[1];

	$count++;	
}
close(INDEX_FAI);

# Use pileup for binomial-based scoring
open(PILEUP,$out_file_sam3);
$count=1;
$prob_success=1-$prob_err;    # Set by user, default is $prob_err == 0.01
%hash_alignment=();
$total_depth=0;
$depth_a=0;   # Gene 5' depth
$depth_z=0;   # Gene 3' depth
$num_indel=0;
$this_allele="-1";
$max_depth=1;

while(<PILEUP>)
{
	chomp;

	my@foo=split(/\s+/);
	$allele=$foo[0];
	$nuc_num=$foo[1];
	$nuc=$foo[2];
	$nuc_depth=$foo[3];
	$allele_size=$size{$foo[0]};
	@aligned_bases=split(//,$foo[4]);
	
	if($count>1 && $allele ne $this_allele)
	{
		# Not the first line and allele alignment has changed in pileup
		$avg_depth=int($total_depth/$count);
		$avg_a=int($depth_a/$edge_a);    # Avg depth at 5' end, num basepairs determined by $edge_a
		$avg_z=int($depth_z/$edge_z);    # 3'

		$hash_max_depth{$this_allele}=$max_depth;

		$hash_edge_depth{$this_allele}="$avg_a\t$avg_z";
		
		### Penalize insertions/deletions and truncations 
		$num_missing=abs($this_allele_size-($this_nuc_num-1))+$num_indel;
		
		for($j=1;$j<=$num_missing;$j++)
		{	
			# Maintain a penalty of at least 5 mismatches if there's a gap
			if($avg_depth>=5)
			{
				$min_penalty=$avg_depth;
			}
			else
			{
				$min_penalty=5;
			}
			
			# Save in hash for later processing in R
			$hash_alignment{$this_allele}[$this_nuc_num+$j]="0\t$min_penalty\t$prob_success";
		}

		if(defined $avg_depth_allele{$allele})
		{
			$avg_depth_allele{$allele} = $avg_depth_allele{$allele} + $avg_depth;
		}
		else
		{
			$avg_depth_allele{$allele}=$avg_depth;
		}
		
		# Reset counters and indicators
		$count=1;
		$total_depth=0;
		$depth_a=$depth_z=0;
		$num_indel=0;
		$this_allele=$allele;
		$this_nuc_num=$nuc_num;
		$max_depth=0;
	}
	
        if($count==1)
        {
		# First line of file or of new allele
                $this_allele=$allele;
                $this_nuc_num=$nuc_num;
                $this_allele_size=$allele_size;
		$max_depth=$nuc_depth;
        }

        if($allele eq $this_allele && $nuc_num != $count)
        {
		# Same allele but alignment skips basepairs
                $num_indel=$num_indel + abs($count-$nuc_num);
                $count=$nuc_num;
                $this_nuc_num=$nuc_num;
        }

	if($allele eq $this_allele)
	{
		# Same allele, calculate depths and update counters
		if($count>0 && $count<($edge_a+1))
		{
			$depth_a=$depth_a+$nuc_depth;
		}
		if(abs($count-$allele_size)<$edge_z)
		{
			$depth_z=$depth_z+$nuc_depth;
		}

		if($nuc_depth>$max_depth)
		{
			$hash_max_depth{$allele}=$nuc_depth;
			$max_depth=$nuc_depth;
		}

		$total_depth=$total_depth+$nuc_depth;
		$this_nuc_num++;
		$count++;
	}


	$num_match=0;

	ONE: for($i=0;$i<scalar(@aligned_bases);$i++)
	{
		if ($aligned_bases[$i] eq "^")
		{
			# Signifies start of a read, next char is mapping quality (skip it)
			$i=$i+2;
		}

		if($aligned_bases[$i] eq "." ||
		$aligned_bases[$i] eq ",")
		{
			$num_match++;
		}
	}
	
	$num_mismatch=$nuc_depth - $num_match;

	# Hash for later processing in R
	$hash_alignment{$allele}[$nuc_num]="$num_match\t$num_mismatch\t$prob_success";
}

close(PILEUP);



# Table for R processing
open(TABLE,">$out_file_sam3.table");

for $allele ( keys %hash_alignment )
{
	for $nuc_num ( 1 .. $#{ $hash_alignment{$allele} })
	{
		if(defined $hash_alignment{$allele}[$nuc_num])
		{
			$this=$hash_alignment{$allele}[$nuc_num];
			($a,$b,$c,$d)=split(/\t/,$this);
			
			if($a > 0 || $b > 0)
			{
				print TABLE "$allele\t$nuc_num\t$this\t$hash_max_depth{$allele}\n";
			}
		}
		else
		{
			#print TABLE "$allele\t$nuc_num\tNA\tNA\tNA\tNA\n";
		}
	}
}

close(TABLE);




### Generate R script for processing table file

print STDERR "Scoring alleles...\n";
open(FILER,">$out_file_sam3.table.R");

###############################################
print FILER <<EOM;

a<-read.table("$out_file_sam3.table")
totalbp<-dim(a)[1]
loop<-1:totalbp
j<-1
pvals<-NULL
score_table<-NULL

for(i in loop)
{
pvals[j] <- binom.test(c(a\$V3[i],a\$V4[i]),p=a\$V5[i],alternative="less")\$p.value

# Weight pvalue by (depth/max_depth)
weight <- (a\$V3[i]+a\$V4[i]) / a\$V6[i]
pvals[j] <- pvals[j] * weight

pvals[j] <- -log10(pvals[j])

if(i==totalbp || a\$V1[i] != a\$V1[i+1])
{
# Fit linear model to observed Pval distribution vs expected Pval distribution (QQ plot)
obs_pvals <- sort(pvals,decreasing=TRUE)
exp_pvals <- c(1:length(obs_pvals))
exp_pvals2 <- -(log10(exp_pvals / (length(exp_pvals)+1)))
qq_lm<-lm(obs_pvals ~ exp_pvals2)
slope<-qq_lm\$coef[2]   # Slope is score
allele<-toString(a\$V1[i])
score_table<-rbind(cbind(allele,slope),score_table)
pvals<-NULL
j<-0
}
j<-j+1
}

write.table(score_table,file="$out_file_sam3.table.scores")

EOM
###############################################

close(FILER);

# Process the data in R
`R CMD BATCH $out_file_sam3.table.R $out_file_sam3.table.Rout`;


# Process the R output and add depths etc
print STDERR "Processing allele scores...\n";
open(SCORES,"$out_file_sam3.table.scores");
open(SCORES2,">$out_file_sam3.table.scores2");
$count=0;
while(<SCORES>)
{
	chomp;
	
	if($count == 0)
	{
		print SCORES2 "Allele\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\n";
	}
	else
	{
		$_ =~ s/"//g;
		$_ =~ s/ /\t/g;
		($a,$b,$c)=split(/\t/);

		if(defined $hash_edge_depth{$b})
		{
			$edge_depth=$hash_edge_depth{$b};
		}
		else
		{
			$edge_depth="NA\tNA";
		}

		if(defined $avg_depth_allele{$allele})
		{
			$this_depth=$avg_depth_allele{$allele}
		}
		else
		{
			$this_depth="NA";
		}

		print SCORES2 "$b\t$c\t$this_depth\t$edge_depth\n";
	}

	$count++;
}

close(SCORES);
close(SCORES2);
`mv $out_file_sam3.table.scores2 $out_file_sam3.table.scores`;


print STDERR "Finished processing for index file $index_file[$counter]...\n";
print STDERR "\n";

### Go back and do next index file (if there is one)
}

print STDERR "SRST2 has finished.\n";
