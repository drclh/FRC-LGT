#!c:\perl -w

use strict;
use warnings;
use diagnostics;
use Statistics::ChisqIndep;
use Text::NSP::Measures::2D::Fisher::right;
use Getopt::Long;

my $infile;
my $outfile;
my $sample_included = 'yes';

GetOptions (
	"infile=s" => \$infile,
	"outfile=s" => \$outfile,
	"sample_included:s" => \$sample_included
) or die "Error in command line arguments\n";

my %samp_count;
my %bact_count;
my $stats_ref;
my @orgs;
my %cog_names;

open (my $samp_fh, "<", $infile);

while (<$samp_fh>) { # parse sample abundance data into $samp_count
	chomp (my ($func_id, $func_name, @counts) = split /\t/);
	if ($func_id eq 'Func_id') {
		@orgs = @counts;
	} else {
		$cog_names{$func_id} = $func_name unless $cog_names{$func_id};
		foreach my $i (0..$#orgs) {
			$samp_count{$func_id}{'counts'}{$orgs[$i]} += $counts[$i];
			$samp_count{$func_id}{'counts'}{'all_ref'} += $counts[$i];
		}
	}
}
close $samp_fh;

opendir(my $dir_dh, ".");
my @files = grep {/^abund\d+\.txt/} readdir($dir_dh);
close $dir_dh;

foreach (@files) { # parse bacteria abundance data into $bact_count
	chomp;
	open (my $fh, "<", $_);
	while (<$fh>) {
		chomp (my ($func_id, $func_name, @counts) = split /\t/);
		if ($func_id eq 'Func_id') {
			@orgs = @counts;
		} else {
			$cog_names{$func_id} = $func_name unless $cog_names{$func_id};
			foreach my $i (0..$#orgs) {
				$bact_count{$func_id}{'counts'}{$orgs[$i]} += $counts[$i];
				$bact_count{$func_id}{'counts'}{'all_bact'} += $counts[$i];
			}
		}
	}
	close $fh;
}

$stats_ref = &calculateOdds({%samp_count}, {%bact_count});

open (my $out_fh, ">", $outfile);
print $out_fh "Reference Organism\tCOG\tCOG Function\t#COGXXX in reference\t#COGs in reference\t#COGXXXX in bact\t#COGs in bact\tOdds Ratio\tChi Squared Value\tp value\tFisher's Exact Test (1-tailed right)\n";

grep {
	my $cog = $_;
	my $name = $cog_names{$cog};
	grep {
		my $org = $_;
		my $odds_ratio = sprintf( "%.3f", $stats_ref->{$cog}{$org}{'odds ratio'} );
		my $chi_value = sprintf( "%.3f", $stats_ref->{$cog}{$org}{'chi'}{'chisq_statistic'} );
		my $chi_p_value = sprintf( "%.3f", $stats_ref->{$cog}{$org}{'chi'}{'p_value'} );
		my $fisher_right = sprintf( "%.3f", $stats_ref->{$cog}{$org}{'fisher_right'} );
		my $ref_COG = $stats_ref->{$cog}{$org}{'counts'}{'ref_COG'};
		my $ref_all = $stats_ref->{$cog}{$org}{'counts'}{'ref_all'};
		my $bact_COG = $stats_ref->{$cog}{$org}{'counts'}{'bact_COG'};
		my $bact_all = $stats_ref->{$cog}{$org}{'counts'}{'bact_all'};
		print $out_fh "$org\t$cog\t$name\t$ref_COG\t$ref_all\t$bact_COG\t$bact_all\t$odds_ratio\t$chi_value\t$chi_p_value\t$fisher_right\n";
	} sort keys %{ $stats_ref->{$cog} };
} sort keys %{ $stats_ref };

my $temp;

sub calculateOdds {
	my ($references, $all) = @_;
	my (%a, %b, %c, $d);
	my %stats;
  
	grep { # parse references
		my $cog = $_;
		grep {
			my $org = $_;
			my $name = $references->{$cog}{'name'};
			$a{$cog}{$org} = $references->{$cog}{'counts'}{$org};
			$b{$org} += $references->{$cog}{'counts'}{$org};
		} sort keys %{ $references->{$cog}{'counts'} };
	} sort keys %{ $references };

	grep { # parse all bacteria
		my $cog = $_;
		grep {
			my $org = $_;
			$c{$cog} += $all->{$cog}{'counts'}{$org};
			$d += $all->{$cog}{'counts'}{$org};
		} sort keys %{ $all->{$cog}{'counts'} };
	} sort keys %{ $all };

	grep { # calculate odds ratio, chi squared
		my $cog = $_;
		grep {
			my $org = $_;
			
			my ($local_a, $local_b, $local_c) = ($a{$cog}{$org}, $b{$org}, $c{$cog});
			
			if ($sample_included eq 'yes') {
				$local_c -= $local_a;
				$d -= $local_b;
			}
			
			$local_c = 1 if $local_c == 0;
			$stats{$cog}{$org}{'odds ratio'} = sprintf ("%.2f", ($local_a/$local_b) / ($local_c/$d));
			my $chisq = new Statistics::ChisqIndep;
			my @obs = ([$local_a, $local_b], [$local_c, $d]);
			$chisq->load_data(\@obs);

			my $n1p = $local_a + $local_b;
			my $np1 = $local_c;
			my $npp = $n1p + $d;
			my $right_value = calculateStatistic( 
				n11=>$local_a,
		                n1p=>$n1p,
                                np1=>$np1,
                                npp=>$npp
                        );
			
			$stats{$cog}{$org}{'counts'} = {
				'ref_COG' => $local_a,
				'ref_all' => $local_b,
				'bact_COG' => $local_c,
				'bact_all' => $d
			};
			$stats{$cog}{$org}{'chi'} = $chisq;
			$stats{$cog}{$org}{'fisher_right'} = $right_value;

		} sort keys %{ $a{$cog} };
	} sort keys %a;
	
	return { %stats };
}

#partial pod of Text::NSP::Measures::2D::Fisher::right that I included for my own reference

=head1 DESCRIPTION

Assume that the frequency count data associated with a bigram
<word1><word2> is stored in a 2x2 contingency table:

          word2   ~word2
  word1    n11      n12 | n1p
 ~word1    n21      n22 | n2p
           --------------
           np1      np2   npp

where n11 is the number of times <word1><word2> occur together, and
n12 is the number of times <word1> occurs with some word other than
word2, and n1p is the number of times in total that word1 occurs as
the first word in a bigram.

The odds ratio computes the ratio of the number of times that
the words in a bigram occur together (or not at all) to the
number of times the words occur individually. It is the cross
product of the diagonal and the off-diagonal.

Thus, ODDS RATIO = n11*n22/n21*n12

if n21 and/or n12 is 0, then each zero value is "smoothed" to one to
avoid a zero in the denominator.

=over

=cut
