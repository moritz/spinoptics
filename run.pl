#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use Parallel::ForkManager;
use Data::Dumper;


my @hosts = glob "wvbh07{0,1,2,3,3,5,8} wvbh06{6,9} wthp009 wthp010 wthp10{4,4,5,5,6,6}";
my $parallel_jobs = @hosts;
my $revoke;
$revoke = 1 if $ARGV[0] && $ARGV[0] eq 'revoke';

my $what =  shift(@ARGV) // 'phi';


my $dir = 'data/' . int(rand() * 10_000);
if (@ARGV) {
    $dir = "data/$ARGV[0]";
    print "Writing data to `$dir'\n";

} elsif (!$revoke) {
    die "bad luck" if -e $dir;
    mkdir $dir or die "Can't mkdir `$dir': $!";
    print "Writing data to `$dir'\n";
}

my %defaults = (
    -b => 0,
    -e => 2.0,
    -r => 0.02,
    -p => 40,
    -n => 21,
);

my %vars = (
    alpha => {
        from    => 0,
        to      => 1.0,
        step    => 0.01,
	option  => '-r',
    },
    phi => {
        name    => 'phi',
        from    => 0,
        to      => 90,
	step	=> 0.1,
	format  => 'phi%04.1f',
	option  => '-p',
    },
);


my %v = %{$vars{$what} // {}};
my $pm = Parallel::ForkManager->new($parallel_jobs);

my $count = -1;
if ($revoke) {
    for my $h (@hosts) {
        my $pid = $pm->start and next;
        system 'ssh', '-x', $h, 'killall', 'cppspin';
        $pm->finish;
    }
} else {
    for (my $var = $v{from}; $var <= $v{to}; $var += $v{step}) {
        $count++;
        my $pid = $pm->start and next;
        my $host = $hosts[$count % @hosts];
        my $fn = sprintf "%s/$v{format}.dat", $dir, $var;
        if (-e $fn) {
            say "$fn exists already, skipping...";
            $pm->finish;
            next;
        }
        my %args = %defaults;
        $args{'-o'} = $fn;
        $args{$v{option}} = $var;
        my @args = ( %args, '-q');
        printf "START: (%s) $v{format}\n", $host, $var;
        my $ret = system('ssh', '-x', $host, './run.sh', @args);
        printf "END:   (%s) $v{format}\n", $host, $var;
        if ($ret != 0) {
            warn "can't run ssh $host run.sh: $?, $!\n";
            warn "re-running it locally...\n";
            system './cppspin', @args;
        }
        $pm->finish;
    }
}
$pm->wait_all_children;

print "finished run `$dir'\n";

# vim: ft=perl sw=4 ts=4 expandtab

