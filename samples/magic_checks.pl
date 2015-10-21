#!/bin/sh
#
# Auto-test function: adapted (and simplified) from the pc_auto-test function
# of the pencil-code
#
# Run the right perl version:
if [ -x /usr/local/bin/perl ]; then
  perl=/usr/local/bin/perl
elif [ -x /usr/bin/perl ]; then
  perl=/usr/bin/perl
else
  perl=`which perl| sed 's/.*aliased to *//'`
fi

exec $perl -x -S $0 "$@"
#! /Good_Path/perl -w

use strict;
use Getopt::Long;
use Cwd;
use POSIX qw(floor);
use vars qw{ %failed $test_status };


my $help  = 0;
my $default_level = 1;
my $reference_out = 'reference.out';
my (%opts);                     # variables written by GetOptions
my @testdirs; 

#
# Here is defined the runs you want to compute by priority levels
# Default level is zero
#
my %tests = (
             0 => [qw(
                         dynamo_benchmark
                         varProps
                         boussBenchSat
                    )],
             1 => [qw(
                         testRestart
                         testTruncations
                         testMapping
                         testOutputs
                         testRadialOutputs
                    )],
             2 => [qw(
                         hydro_bench_anel
                    )],
             3 => [qw(
                         fluxPerturbation
                         isothermal_nrho3
                         dynamo_benchmark_condICrotIC
                         varCond
                    )],
            );

my $itest = 1;
my $ntests      = 0;                 # total number of tests run
my $failed      = 0;
my $indi_fmt = '%-17s'; 


my $usage=
"Usage:  magic_checks.pl [options]
  Test compilation and results on a set of sample directories (or on the
  list of directories given in the command line) to validate MAGIC.
  Uses Perl because we need to do pattern matching to extract the relevant
  output lines.

  Adapted from the auto-test subroutines of the pencil-code developed by
  W. Dobler

Options:
  -h,  --help              \tShow usage overview
  -c,  --clean             \tClean the directories when it is finished
  -a,  --all               \tAll auto-tests are computed
       --level=LEV         \tRun only tests from level LEV
       --max-level=LEV     \tRun all tests below with level <= LEV (default: 0)
       --hybrid            \tRun the hybrid version
       --use-mkl           \tUse the MKL for FFTs and Lapack calls
       --use-cmake         \tUse Cmake instead of make

Example:
  magic_checks.pl --clean
  magic_checks.pl --level=1
  magic_checks.pl --clean --all
  magic_checks.pl --clean --all --hybrid --use-cmake --use-mkl
";


eval {
    Getopt::Long::config("bundling"); # make single-letter opts. case-sensitive
};
GetOptions(\%opts,
           qw( -h   --help
               -c   --clean
               -a   --all
                    --level=s
                    --max-level=s
                    --hybrid
		    --use-mkl
		    --use-cmake
                    )
          ) or $help=1, die "Aborting.\n";

if ($opts{'h'} || $opts{'help'}) { $help=1; die "$usage\n"; }
my $clean=($opts{'C'} || $opts{'clean'} || 0 );
my $all=($opts{'a'} || $opts{'all'} || 0 );
my $level=($opts{'level'});
my $max_level=($opts{'max-level'});
my $hybrid=($opts{'hybrid'} || 0 );
my $cmake=($opts{'use-cmake'} || 0 );
my $mkl=($opts{'use-mkl'} || 0 );
my $MAGIC_HOME = "";
my $OMP_NUM_THREADS = 1;
my $execdir = "$MAGIC_HOME/tmp";

if (!exists($ENV{MAGIC_HOME})) {
    die "MAGIC_HOME is not set!\n\tRun one of the sourceme scripts in the MAGIC base directory.!";
} else {
    $MAGIC_HOME="$ENV{MAGIC_HOME}";
}

if ((!exists($ENV{OMP_NUM_THREADS})) && $hybrid) {
    die "You must set OMP_NUM_THREADS for the hybrid version first.";
} else {
    $OMP_NUM_THREADS = "$ENV{OMP_NUM_THREADS}";
}

# Make sure we are in the top directory and have the right PATH
die "Need to set environment variable MAGIC_HOME\n"
    unless (defined($MAGIC_HOME));
my $topdir = "$MAGIC_HOME";
chdir $topdir;
$ENV{PATH} .= ":$MAGIC_HOME/bin";

# Timing 
my %t_global = (
                'compile'   => 0,
                'run' => 0,
               );
my $t0_global = get_time();         # remember start time
my $t1;



# Define sample dirs
my @sampdirs;
if ($all) {
    $max_level=3;
    @sampdirs = get_all_sampdirs($max_level);
}
else {
    if (defined($level)) {
        @sampdirs = get_sampdirs($level);
    } else {
        $max_level = $default_level unless defined($max_level);
        @sampdirs = get_all_sampdirs($max_level);
    }
}
@testdirs =
  ( map { "$topdir/samples/$_" } @sampdirs,
  );

$ntests = @testdirs;

# Compile
print "Compilation.. ";
print "\n";
if ( $cmake ) {
    $execdir = "$MAGIC_HOME/tmp";
    mkdir "$topdir/tmp";
    chdir "$topdir/tmp";
    if ($hybrid) {
        if ($mkl){
            `cmake .. -DUSE_FFTLIB=MKL -DUSE_MKL=yes -DOPENMP=yes`;
        }
        else {
            `cmake .. -DUSE_FFTLIB=JW -DUSE_MKL=no -DOPENMP=yes`;
        }
    }
    else {
        if ($mkl){
            `cmake .. -DUSE_FFTLIB=MKL -DUSE_MKL=yes -DOPENMP=no`;
        }
        else {
            `cmake .. -DUSE_FFTLIB=JW -DUSE_MKL=no -DOPENMP=no`;
        }
    }
    `make -j`;
    $t1 = get_time();
    $t_global{'compile'} += ($t1-$t0_global);
    chdir $topdir;
}
else {
    $execdir = "$MAGIC_HOME/src";
    chdir "$topdir/src";
    `make clean`; # clean directory
    `make -j`;
    $t1 = get_time();
    $t_global{'compile'} += ($t1-$t0_global);
}
 
# Run tests consecutively
for my $d (@testdirs) {
    test_rundir("$d");
}

# Cleaning
print "Clean.. ";
if ( $cmake ) {
    `rm -fr $topdir/tmp`;
}
else {
    chdir "$topdir/src";
    `make clean`;
}




# ---------------------------------------------------------------------- #
sub get_all_sampdirs {
# Return a flat list of all sample dirs that belong to a level <= $level .
    my $level = shift;

    my @dirs;
    my @levels = grep { $_ <= $level } keys %tests;
    foreach my $lev (sort {$a <=> $b} @levels) {
        push @dirs, @{$tests{$lev}};
    }

    return @dirs;
}
# ---------------------------------------------------------------------- #
sub get_sampdirs {
# Return list of sample dirs that belong to $level .
    my $level = shift;

    my @dirs;
    my @levels = grep { $_ == $level } keys %tests;
    foreach my $lev (sort {$a <=> $b} @levels) {
        push @dirs, @{$tests{$lev}};
    }

    return @dirs;
}
# ---------------------------------------------------------------------- #
sub test_rundir {
    my $dir = shift;
    my $t0  = get_time();
    my $t1;

    $test_status = 0;           #  so far, everything is OK

    # Go to directory and identify it
    if (! -d $dir) {
        print STDERR "No such directory: $dir\n";
        return;
    }
    chdir $dir;
    my $cwd = cwd();
    print "\n$cwd:";
    print " ($itest/$ntests)" if ($ntests>1);
    print "\n";
    $itest++;

    #0. Make sure nothing is here yet
    `magic_clean`;
    `ln -s $execdir/magic.exe .`;

    `cp $MAGIC_HOME/bin/run_magic.sh .`;

    #2. Run
    my $t2 = get_time();
    if ( -e 'runMe.sh' ) {
        if ( $hybrid ) {
            `bash runMe.sh hybrid $OMP_NUM_THREADS &> /dev/null`;
        } else {
            `bash runMe.sh &> /dev/null`;
        }
    } else { 
        if ( $hybrid ) {
            `./run_magic.sh hybrid $OMP_NUM_THREADS`; 
        } else {
            `./run_magic.sh &> /dev/null`; 
        }
    }
    my $exitcode=$?;
    print "exitcode = $exitcode";
    test_results($dir);

    my $t3 = get_time();
    $t_global{'run'} += ($t3-$t2);

    # Summarize timings in human-readable form
    my $t4 = get_time();
    print "    Time used: ",
    s_to_hms(time_diff($t0,$t3), 44),
    "\n";

    if ($clean && $exitcode==0) {
        `magic_clean`;
    };
}
# ---------------------------------------------------------------------- #
sub get_time {
# Wrapper around time() that uses Time::HiRes if possible
    my $time = 0;
    eval {
        require Time::HiRes;
        $time = Time::HiRes::time();
    };
    return $time || time();
}
# ---------------------------------------------------------------------- #
sub time_diff {
# Return difference of times if both args are defined, undef otherwise
    my $t1 = shift;
    my $t2 = shift;

    if (defined($t1) && defined($t2)) {
        return $t2-$t1;
    } else {
        return undef;
    }
}
# ---------------------------------------------------------------------- #
sub s_to_hms {
# Convert no. of seconds to [ddd][hh:]mm:ss string
    my $secs  = shift;
    my $width = (shift || 0);

    my $string;

    # Not much to do if arg is undef:
    if (! defined($secs)) {
        $string = 'undef';
    } else {
        my $ss = $secs % 60;
        my $mm = floor($secs/60) % 60;
        my $hh = floor($secs/3600) % 24;
        my $dd = floor($secs/86400);

        $string = sprintf("%02d:%02d", $mm,$ss);
        if ($hh) { $string = sprintf("%02d:", $hh) . $string };
        if ($dd) { $string = sprintf("%dd", $dd) . $string };
    }

    if (length($string) < $width) {
        $string = (" " x ($width-length($string))) . $string;
    };

    return $string;
}
# ---------------------------------------------------------------------- #
sub test_results {
# Analyze results from code
    my $dir = shift;

    my ($rdmsg,$diagn);
    print "    Validating results.. ";
    printf "$indi_fmt ", '';

    my @output;

    @output = read_lines('e_kin.test', $rdmsg);
    if ($rdmsg) {
        my_ok(0,1,"[$rdmsg]",$dir,"results");
        return;
    }

    my @ref_output = read_lines($reference_out,$rdmsg);

    if ($#output != $#ref_output) {
        printf("\n    Warning: this case did not run correctly... \n");
        printf("    Check your stacksize limits before trying it again\n");
        $failed{$dir} = "did not run";
        return;
    }

    if ($rdmsg) {
        print " [$rdmsg]\n";
    } else {
        my $comp = compare_results(\@ref_output,\@output,$diagn);
        my_ok($comp,1,$diagn,$dir,"results");
    }
}
# ---------------------------------------------------------------------- #
sub read_lines {
# Read file an return hash of non-empty lines
    my $file = shift;
    my @lines = ();
    my $msg = "";

    {
        local $/ = undef;       # read in whole file
        if (open (REF, "< $file")) {
            @lines = grep { /\S/ } split(/[\n\r]+/,<REF>);
            # Remove leading comment sign from header line:
            $lines[0] =~ s/^(\s*)#/$1 /;
        } else {
            $msg = "Couldn't open $file";
        }
    }
    $_[0] = $msg;
    @lines;
}
# ---------------------------------------------------------------------- #
sub compare_results {
# Compare two arrays of lines linewise; if all lines are the same, return
# 1; if not, return 0 and write report to third argument
    my $arr1 = shift;
    my $arr2 = shift;

    my @pref = ("< ", "> ", "  "); # prefixes to mark differing lines
    my $n1 = $#$arr1;
    my $n2 = $#$arr2;
    my $N = ($n2>$n1 ? $n2 : $n1);
    my $diagn = "";
    my $equal = 1;
    for my $i (0..$N) {
        my $line1 = ($$arr1[$i] || ""); chomp $line1;
        my $line2 = ($$arr2[$i] || ""); chomp $line2;
        unless (compare_lines_fuzzily($line1,$line2)) {
            $diagn .= $pref[0] . $line1 . "\n" if ($line1);
            $diagn .= $pref[1] . $line2 . "\n" if ($line2);
            $equal = 0;
        } else {
            if ($i == 0) {      # Keep header line for easier reading
                $diagn .= $pref[2] . $line1 . "\n";
            }
        }
    }
    $_[0] = $diagn;             # the remaining argument (was the third one)
    $equal;
}
# ---------------------------------------------------------------------- #
sub compare_lines_fuzzily {
# Compare the numerical values in two lines fuzzily. Return 1 if numbers
# are approximately equal (differ by 1 or less in the last decimal), zero
# otherwise.
    my $xline = shift;
    my $yline = shift;

    my $equal = 1;

    my @x = split(/\s+/,$xline);
    my @y = split(/\s+/,$yline);
    @x = grep(/([0-9]|NaN|Inf)/i, @x);      # weed out empty entries
    @y = grep(/([0-9]|NaN|Inf)/i, @y);

    return 0 unless ($#x == $#y); # lengths must match
    foreach my $i (0..$#x) {
        my ($x,$y) = ($x[$i], $y[$i]);
        $equal = 0 unless (compare_numbers_fuzzily($x,$y));
    }

    return $equal;
}
# ---------------------------------------------------------------------- #
sub compare_numbers_fuzzily {
# Compare two numbers and return true if they differ only by one in the
# last digit, else false.
    my $x = shift;
    my $y = shift;

    # Short-circuit NaN/Inf -- they will be taken care of as they will make
    # any comparison fail
    return 1 if (   ($x =~ /^([+-]?(NaN|Inf))$/i)
                 || ($y =~ /^([+-]?(NaN|Inf))$/i) );

    # Regexp for general float, groups into pre-exponent and
    # post-exponent part
    my $numeric = '((?:[+-]?)(?=\d|\.\d)\d*(?:\.\d*)?)((?:[EeDd](?:[+-]?\d+))?)';
    my ($x1,$x2) = ($x =~ /^\s*$numeric\s*$/);
    $x2 = '' unless defined($x2);
    $x  =~ s/[dD]/e/;
    #$x1 =~ s/[dD]/e/;
    $x2 =~ s/[dD]/e/;

    my ($y1,$y2) = ($y =~ /^\s*$numeric\s*$/);
    $y2 = '' unless defined($y2);
    $y =~ s/[dD]/e/;
    #$y1 =~ s/[dD]/e/;
    $y2 =~ s/[dD]/e/;

    # If too small (1e-20), don't test:
    return 1 if (abs("$x1$x2") < 1e-18 and "$x1$x2" != 0.);
    # If diff too small (1e-15), don't test:
    return 1 if (abs("$x1$x2"-"$y1$y2") < 1e-14);

    # Are $x, $y really numeric?
    unless(defined($x1) && defined($y1)) {
        warn "Not a numerical value: <$x> or <$y>\n";
        return 0;
    }

    # Short-circuit identical numbers (so Perl's rounding errors to double
    # precision can't screw them up)
    return 1 if ("$x1#$x2" eq "$y1#$y2");

    # Choose the longer string for decrementing/incrementing, so we
    # correctly compare stuff like [1.3599, 1.36] (succeeds) and [1.3598,
    # 1.36] (fails) in both directions.
    if (length($x1) < length($y1)) {
        ($x,$x1,$x2, $y,$y1,$y2) = ($y,$y1,$y2, $x,$x1,$x2);
    }

    # Do comparison
    $x1 =~ /([+-]?[0-9]*)\.?([0-9]*)/;
    my ($xint,$xfrac) = ($1, $2);

    die "Cannot split <$x1> into int/fract parts (should not happen)"
      unless defined($xint);
    my $delta;
    if (length($xfrac) > 0) {
        $delta = '0.' . '0'x(length($xfrac)-1) . '1';
    } else {
        $delta = 1;
    }
    my $x1_p = $x1 + $delta;
    my $x1_m = $x1 - $delta;
    ($x1_m, $x1_p) = sort  { $a <=> $b } $x1_m, $x1_p;
    my ($x_m, $x_p) = ("${x1_m}${x2}", "${x1_p}${x2}");

    # Are the constructed numbers $x_m, $x_p really numerical?
    eval { my $dummy = 1 + $x_m + $x_p };
    if ($@) {
        warn "$@";
        return 0;
    }

    if (($y >= $x_m) && ($y <= $x_p)) {
        return 1;
    } else {
        return 0;
    }
}
# ---------------------------------------------------------------------- #
sub my_ok {
# Similar to Test's and Test::Harness' ok() function: consider success if
# first two argumets are equal, otherwise report a problem and print the
# third argument (which should normally contain the output from the shell
# calls that we are normally testing here).
# Args #4 and #5 are the current run directory and the phase (compilation,
# starting, running) we are in.

    my $arg1 = shift;
    my $arg2 = shift;
    my $mesg = (shift || "<No diagnostics>");
    chomp($mesg);
    my $dir = shift;
    my $phase = shift;
    my $dt    = shift;

    my $timestr;
    $timestr = '';

    # Allow for calls like `ok(0)' or `ok(1)':
    if (!defined($arg2)) {
        $arg2 = $arg1 ? $arg1 : 1;
    }
    if ($arg1 eq $arg2) {
        print " ok      $timestr\n";
    } else {
        print " not ok: $timestr\n$mesg\n";
        $test_status = 1;
        # report only first failure:
        $failed{$dir} = $phase unless defined($failed{$dir});
    }

}

## Summarize results
END {
    unless ($help) {
        print "\n" . "-" x 70 . "\n";

        # Print failure header that can be identified by pencil-test
        if ($failed || %failed) {
            print "### auto-test failed ###\n";
        }

        # Failed outside individual tests (e.g. other auto test is running)
        if ($failed) {
            print "Crash...";
        }

        # Failed during some of the individual tests
        if (%failed) {
            print "Failed ", scalar(keys %failed),
              " test(s) out of $ntests:\n";
            while (my ($dir,$phase) = each %failed) {
                print "  $dir ($phase)\n";
            }
        } else {
            if ($ntests == 1) {
                print "Test succeeded.\n";
            } elsif ($ntests > 1) {
                print "All $ntests tests succeeded.\n";
            }
            elsif (($ntests < 1) && ! $failed) {
                print "There was no test to run???\n";
            }
        }

        # Print timing numbers
        my @t = times();
        print "\nCPU time (including compilation): ",
              s_to_hms($t[2], 7) . "u ",
              s_to_hms($t[3]   ) . "s\n";

        my $t1_global = time(); # end time
        print "Total wall-clock time:            ",
                 s_to_hms(time_diff($t0_global,$t1_global), 7),
          " = ", s_to_hms($t_global{'compile'}),
          " + ", s_to_hms($t_global{'run'}),
          "\n";

    }

}
