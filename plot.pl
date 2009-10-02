#!/usr/bin/perl
use strict;
use warnings;
use lib 'lib';
use ReadMatrix qw(real_matrix);
use Data::Dumper;

my $d = shift(@ARGV) // die "Usage: $0 <datadir>";

my @files = glob "data/$d/*.dat" or die "No files in `data/$d/*'";

my @d;
for my $fn (@files) {
    $fn =~ /phi(\d\d)\.dat/ or die "can't work with file name `$fn'";
    my $phi = $1;
    my $m = tpq_for_file($fn);
    print $phi, " ", $m->[0][7], $/;
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
