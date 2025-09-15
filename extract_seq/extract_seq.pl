 #!usr/bin/perl -w

use experimental 'smartmatch';
#use warnings;

open FILE1, $ARGV[0] or die; # reference fasta file

open FILE2, $ARGV[1] or die; # Need to be a file with one id per line. If it is a fasta need to be grepped (grep '>') before.

#my $output = $ARGV[2];
#open OUT, ">$output" or die;

my %sequences;
my $id;
my $id_2;
my @list;

while (<FILE1>) {
	chomp;

	if (/^>/) {
		$id = $_;
		$id =~ s/>//;
#		print "ID $_\n";
		#$sequences{$id} = '';

	} else {

		$sequences{$id} .= $_;
#		print "SEQ $_\n";

	}

}

while (<FILE2>) {
        chomp;
        push @list, $_;

		}

foreach $id (sort keys %sequences) {
        my @id = split(/\s/,$id);
        $id_2 = $id[0];
        if ($id_2 ~~ @list) {

                print ">$id\n$sequences{$id}\n";
#               print "$name\n$taxa{$name}\n";
    }

 }

#close OUT;
close FILE1;
close FILE2;