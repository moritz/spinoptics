#! /usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(max sum);

# a smaller value than 1e-8 does not make sense because of the data file precision
# my $epsilon = 5e-7;
my $epsilon = 1e-5;

my @tpq_sets;

my $max_index = 0;
die "Usage: $0 tpq+b-file tpq-b-file\n" unless @ARGV == 2;
for my $fn (@ARGV) {
    push @tpq_sets, read_file($fn);
}

sub read_file {
    my $fn = shift;
    open my $h, '<', $fn or die "Can't open '$fn' for reading: $!@";
    my @items;
    while (<$h>) {
        s/\[.*?\]//g;
        push @items, $1 while m/(\b\d([\d.eEdD-]*\d\b)?)/g;
    }
    my $count = scalar @items;
    my $size = int sqrt($count) + 0.5;
    die "Read $count numbers from $fn, and that's not a square\n"
        unless $size * $size == $count;
    my @matrix = ([(undef) x ($size + 1)]);

    for (0..($size-1)) {
        my $s = $_ * $size;
        push @matrix, [undef, @items[$s .. $s+$size-1]];
    }
    $max_index = $size;
    return \@matrix;
}

my $problems = 0;

# my @a = @{$tpq_sets[0]};
# print "$tpq_sets[0]\n";
# print "$a[1][1]\n";
# exit(1);


my %checks;
sub check1() {
    my ($set, $x, $y) = @_;
    return ($set, $y, $x);
}
$checks{"time rev"} = \&check1;


# for spin_hall:
# my @leads = (1,2,5,7);
# my %other_spin = (1, 3, 3, 1, 2, 4, 4, 2, 5, 6, 6, 5, 7, 8, 8, 7);
# for spin_optics: (!!! check again)
my @leads = (1,3,5,7);
my %other_spin = (1, 2, 2, 1, 3, 4, 4, 3, 5, 6, 6, 5, 7, 8, 8, 7);

sub check2() {
    my ($set, $x, $y) = @_;

    my $invspin_x = $other_spin{$x};
    my $invspin_y = $other_spin{$y};

#    print "$x, $y -> $invspin_y, $invspin_x \n";
    return ($set, $invspin_y, $invspin_x);
}
$checks{"time + spin rev"} = \&check2; 


sub check3() {
    my ($set, $x, $y) = @_;

    my $other_field = ($set+1)%2;

    # no spinor inversion:
    my $invspin_x = $x;
    my $invspin_y = $y;
    # spinor inversion:
    # my $invspin_x = $other_spin{$x};
    # my $invspin_y = $other_spin{$y};

#    print "$x, $y -> $invspin_y, $invspin_x \n";
    return ($other_field, $invspin_y, $invspin_x);
}
$checks{"time + B-field rev"} = \&check3;

sub check4() {
    my ($set, $x, $y) = @_;

    my $other_field = ($set+1)%2;

    my $invspin_x = $other_spin{$x};
    my $invspin_y = $other_spin{$y};
#    my $invspin_x = $x;    # test: no spinor inversion
#    my $invspin_y = $y;

#    print "$x, $y -> $invspin_y, $invspin_x \n";
    return ($other_field, $invspin_y, $invspin_x);
}
$checks{"time + spin + B-field rev"} = \&check4;

my $summary = "";
foreach my $check (keys(%checks)) {
    print "performing symmetry check: $check\n";
    my $fpointer = $checks{$check};
    my $failcount = 0;
    my $count = 0;

    for my $x (1 .. $max_index) {
        for my $y (1 .. $max_index) {
            my ($set, $tx, $ty) = &$fpointer(0, $x, $y);
            my @tpq1 = @{$tpq_sets[0]};
            my @tpq2 = @{$tpq_sets[$set]};
            # print "set = $set, tpq[1,1] = $tpq1[1][1]\n";

            my $diff = abs($tpq1[$x][$y] - $tpq2[$tx][$ty]);
            # print "diff = $diff\n";

            if ($diff >= $epsilon) {
                print " comparision tpq1[$x, $y] with tpq2[$tx, $ty]:  diff = $diff\n";
                print "tpq1[$x,$y] = $tpq1[$x][$y]\n";
                print "tpq2[$tx,$ty] = $tpq2[$tx][$ty]\n";
                $failcount++;
            }
            $count++;
        }
    }
    my $msg = "check $check: $failcount failures of $count\n";
    $summary .= $msg;
    print "end of $msg\n"; 
}



#! sum rules: this needs correction for spin dependency
my @tpq = @{$tpq_sets[0]};

# $tpq[0] = [0];  # the index [0] hasn't been used yet
# my @tx = map { ${$_}[0] = 0; sum @$_ } @tpq;
# print @tx;


for my $x (@leads) {
    my ($sumq11, $sumq12, $sumq21, $sumq22) = (0,0,0,0);
    my ($sump11, $sump12, $sump21, $sump22) = (0,0,0,0);
    for my $y (@leads) {
        # print "$y $other_spin{$y}\n";
        $sumq11 += $tpq[$x][$y];
        $sumq12 += $tpq[$x][$other_spin{$y}];
        $sumq21 += $tpq[$other_spin{$x}][$y];
        $sumq22 += $tpq[$other_spin{$x}][$other_spin{$y}];
        $sump11 += $tpq[$y][$x];
        $sump12 += $tpq[$y][$other_spin{$x}];
        $sump21 += $tpq[$other_spin{$y}][$x];
        $sump22 += $tpq[$other_spin{$y}][$other_spin{$x}];
    }
    print "\np = $x:\n";
    print "sum_q  with spin(up->up) = $sumq11\n";
    print "sum_q  with spin(up->down) = $sumq12\n";
    print "sum_q  with spin(down->up) = $sumq21\n";
    print "sum_q  with spin(down->down) = $sumq22\n";
    my $sumq_all = $sumq11 + $sumq12 + $sumq21 + $sumq22;
    print "sum_q sum_spin(i,j) = $sumq_all\n";
    # this value appears in tsum.dat:
    # print $sumq11 + $sumq12 . "\n";

    print "\nq = $x:\n";
    print "sum_p  with spin(up->up) = $sump11\n";
    print "sum_p  with spin(up->down) = $sump12\n";
    print "sum_p  with spin(down->up) = $sump21\n";
    print "sum_p  with spin(down->down) = $sump22\n";
    my $sump_all = $sump11 + $sump12 + $sump21 + $sump22;
    print "sum_p sum_spin(i,j) = $sump_all\n";
    # print $sump22 + $sump12 . "\n";
    
    my $d1 = $sump_all - $sumq_all;
    print "j=$x: sum_i_sig_sig' t[i_sig, j_sig'] - sum_j_sig_sig' t[i_sig, j_sig'] = " 
       . $d1 . "\n";
    my $msg;
    if (abs($d1) > $epsilon) { 
        $msg = "spin independent sum rule:  not met for index $x\n";
        print "$msg\n";
        $problems++;
    } else {
        $msg = "spin independent sum rule: ok for index $x\n";
    }
    $summary .= $msg;
}




# if ($problems){
#    print "Some symmetries were broken\n";
#    exit 1;
# } else {
#     print "Everything is fine.\n";
# }
print "\nSummary:\n";
print $summary;






# vim: sw=4 ts=4 expandtab

