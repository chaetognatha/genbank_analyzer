# /usr/bin/perl

use strict;
use warnings;

if (!open(GENBANK, "genbank.txt")) {
	 die "Check input file 1\n";
}
my @genbank = <GENBANK>;
close (GENBANK);

my $f = join("", @genbank);
my $c = 0;
my $trembl_counter=0;
my @starts;
my @stops;
my $cds_length=0;
my $interspace=0;

#calculate amount of CDS
while ($f =~ /CDS\D+(\d+)\.\.(\d+)/g){
  $c++;
  push (@starts, $1);
  push (@stops, $2);
}

#put starts and stops in multiarray
my @table = (\@starts, \@stops);
my $cds_count = scalar(@starts);
for (my $i=0;$i<$cds_count;$i++){
			$cds_length = $cds_length + (${$table[1]}[$i] - ${$table[0]}[$i]);
			if (defined ${$table[1]}[$i+1]){#when the value is not defined, exit loop
				$interspace = $interspace + (${$table[1]}[$i+1] - ${$table[0]}[$i]);
			} else {last;}
}
#extract TrEMBL
if (!open(TREM, ">sequencesTrembl.txt")) {
	 die "Check input file 2\n";
}

my @translation = split(/UniProtKB\//, $f); #cutting so interesting stuff is at beginning of elements
foreach (@translation){
	if ($_ =~ /TrEMBL:(\w+)/){#extract accession number
		print TREM ">$1\n"; #print first line, fasta format
		$trembl_counter++; #count the trembl CDSs
	}
	if ($_ =~ /translation="(\w+[\n\s+\w+]+)/){ #get the aa sequence
		my $x = $1;
		$x =~ s/\s+//g; #get rid of whitespace
	  $x =~ s/([^\n]{80})/$1\n/g; #insert newline every 80 characters
		print TREM "$x\n"; #put sequence into file and end in a newline
	}
	if ($_ =~ /ORIGIN/){ #breaks out of loop at ORIGIN which is last in file
		last;
	}
}
my $percent_cds = ($cds_length/($cds_length+$interspace))*100;
$percent_cds = int $percent_cds;
my $empty_cds = $cds_count-$trembl_counter;
print "$cds_count CDS entries found.\n
Total summed length of CDS: $cds_length nts.\n
Total summed length of intergene regions: $interspace nts.\n
Percentage of CDS:$percent_cds%\n
$trembl_counter of CDS have TrEMBL entries and $empty_cds do not.\n";
