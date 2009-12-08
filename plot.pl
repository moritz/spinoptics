#!/usr/bin/perl
use strict;
use warnings;
use lib 'lib';
use ReadMatrix qw(real_matrix);
use Data::Dumper;
use Getopt::Long;
use 5.010;

GetOptions("upto=s" => \my $upto)
    or warn "some of these options look wrong\n";

my @x;
my @d;
my $isfirst = 1;
for my $d (@ARGV) {
    my $number;
    if ($d =~ m/(\d+)/) {
        $number = $1;
    } else {
        die "Don't know what to do with $d (doesn't contain a number)";
    }

    my @files = glob "data/$number/*.dat" or die "No files in `data/$number/*'";

    my $i = 0;
    for my $fn (@files) {
        $fn =~ /(energy|alpha|phi)([\d.]+)\.dat/ or die "can't work with file name `$fn'";
        state $what //= $1;
        my $datapoint = $2;
        if ($what ne $1) {
            die "Can't mix different variables on the X axis ($what vs. $1)";
        }
        if (defined($upto) && ($upto < $datapoint)) {
            last;
        }
        my $m = tpq_for_file($fn);
        if ($isfirst) {
            push @x, $datapoint;
        } else {
            if (abs($x[$i]-$datapoint) > 1e5) {
                die "Can't merge datasets with different X axes\n"
                    . "  for the ${i}th set I expected $x[$i] and got $datapoint";
            }
        }
        my $agg = 0;
        if ($d =~ /uu/) {
            $agg += $m->[0][2];
        } elsif ($d =~ /dd/) {
            $agg += $m->[1][3];
        } elsif ($d =~ /ud/) {
            $agg += $m->[1][2];
        } elsif ($d =~ /du/ ) {
            $agg += $m->[0][3];
        } elsif ($d =~ /diff/) {
            $agg += $m->[0][2] - $m->[1][3];
        } elsif ($d =~ /rel/) {
            $agg += ($m->[0][2] - $m->[1][3]) / (($m->[0][2] + $m->[1][3]) || 1);
        }
        push @{$d[$i]}, $agg;
    } continue {
        $i++;
    }
    $isfirst = 0;
}
use Data::Dumper;
#print Dumper \@x;
#print Dumper \@d;

for my $i (0..$#x) {
    say $x[$i], ' ', join(' ', @{$d[$i]});

}


sub tpq_for_file {
    my $fn = shift;
    open my $handle, '<', $fn or die "Can't open file `$fn' for reading: $!";
    my $m;
    while (<$handle>) {
        if (/^final tpq/) {
            $m = real_matrix($_);
            last;
        }
    }
    close $handle or warn $!;
    return $m;
}

# vim: sw=4 ts=4 expandtab
