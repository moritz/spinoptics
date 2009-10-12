#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use Data::Dumper;

my @hosts = glob "wvbh07{2,3,3,6,7,8}";
my $parallel_jobs = @hosts;

my $dir = 'data/' . int(rand() * 10_000);
die "bad luck" if -e $dir;

mkdir $dir or die "Can't mkdir `$dir': $!";

print "Writing data to `$dir'\n";

my $b = 0;

my $pm = Parallel::ForkManager->new($parallel_jobs);

my $count = -1;
for my $phi (0..90) {
    $count++;
    my $pid = $pm->start and next;
    my $host = $hosts[$count % @hosts];
    print "($host) Phi = ", $phi, " degrees\n";
        my $fn = sprintf "%s/bz%+.2f,phi%02d.dat", $dir, $b, $phi;
        my @args = (
            -b => $b,
            -e => 0.2,
            -o => $fn,
            -r => 0.05,
            -p => $phi,
            '-q',
            -n => 19,
        );
        my $ret = system('ssh', '-x', $host, './run.sh', @args);
        if ($ret != 0) {
            warn "can't run ssh $host run.sh: $?, $!\n";
            warn "re-running it locally...\n";
            system './cppspin', @args;
        }
    $pm->finish;
}

print "finished run `$dir'\n";

# vim: ft=perl sw=4 ts=4 expandtab
