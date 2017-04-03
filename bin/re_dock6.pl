#!/usr/bin/env perl
# 
################################################################################
#
# MIT License
# 
# Copyright (c) 2017 Michael Roach
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
################################################################################
#
# Description: A Rigid Exhaustive extension for UCSF Dock6. See documentation
#              for requirements, installation, and usage. 
#
################################################################################

use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use File::Temp qw/tempfile tempdir/;
use Getopt::Long;
use Time::Piece;
use threads;
use threads::shared;
use Thread::Semaphore;

my $re_dock6_version = "RE_Dock6 version 1.0\n";

my $usage = "
re_dock6.pl  -r  <rec.mol2>  -b  <box.pdb>  -l  <ligands.mol2>  [options]

REQUIRED PARAMETERS:
-rec / -r               Ready-to-dock receptor structure in mol2 format.
-box / -b               Box file in pdb format.
-lig / -l               Molecule(s) to dock in either mol2 or gzipped mol2 (.mol2.gz) format.

OPTIONAL PARAMETERS:
-proc / -p  <INT>       Number of processors/threads to use. Default = 1.
-dir / -d               ABSOLUTE PATH to the directory for output, default is the current 
                        working directory.
-outfile / -o           Specify a file name for the scored conforations, default = scored.mol2.
-tmpdir / -t            ABSOLUTE PATH to the directory in which to create the temp folder  
                        (default = make the temp folder in current directory).
-force / -f             Force overwrite of output file (if it exists). Default = re_dock6 will
                        resume docking and skip already docked ligands.
-version / -v           Print version and exit.
-help / -h              Print this message and exit.

PRESET EFFORT SETTINGS:
-e_sparse / -e_s        Sparse generation of ligand poses. Use this for large, open binding
                        sites.
-e_quick / -e_q         Quick docking. A bit less dense and less energy minimisation than default.
-e_thorough / -e_t      Very dense generation of ligand poses, if the default setting is not
                        quite enough (e.g. very tight binding pocket). Use with caution.

ADVANCED USAGE: DOCK-RELATED PARAMETERS:
-n_emin / -n_e  <INT>   Number of steps of energy minimization. default = 50, sensible range 
                        = 25-200.
-n_confs / -n_c  <INT>  Number of conformers for Dock to output for EACH ligand (default = 1).

ADVANCED USAGE: CONFORMER GENERATION PARAMETERS:
-c_rots / -c_r  <INT>   Number of x,y,z rotations per step to perform (more = more conformations,
                        default = 3). Total rotations per step = n_rot^3.
-c_step / -c_s  <FLOAT> 
                        Step distance across grid box (less = more conformations, default = 1.0,
                        sensible range = 0.5-2.5).
-c_tether / -c_t <FLOAT>
                        Tether conformation generation to within a specified radius of the box
                        center. Default = 3.0, sensible range = 2.0-5.0 for normal docking, or
                        some arbitrary large number to carpetbomb everything in the box.
";

#---PARAMETERS---
my $curdir = getcwd;
my $outfile = "scored";
my $tmp_dir;

# Required parameters
my $recmol2;
my $box;
my $lib;

# optional parameters
my $outdir;
my $tmpdir;
my $force;
my $sparse;
my $quick;
my $thorough;
my $confs = 1;
my $e_min = 50;
my $n_rot = 3;
my $d_step = 1;
my $tether = 3;
my $threads = 1;

# misc
my $version;
my $help;
my @LIG;
my $cur_lig;
my $lig_number = 0;
my $die=0;
my %skip_ligands;

my $command = "$0 @ARGV";

GetOptions (
            "rec=s" => \$recmol2,
            "box=s" => \$box,
            "lig=s" => \$lib,
            "proc=i" => \$threads,
            "dir=s" => \$outdir,
            "tmpdir=s" => \$tmpdir,
            "force" => \$force,
            "outfile=s" => \$outfile,
            "n_confs=i" => \$confs,
            "n_emin=i" => \$e_min,
            "e_sparse" => \$sparse,
            "e_quick" => \$quick,
            "e_thorough" => \$thorough,
            "c_rots=i" => \$n_rot,
            "c_step=f" => \$d_step,
            "c_tether=f" => \$tether,
            "version" => \$version,
            "help" => \$help
			)
or die($usage);

if($help){
    die $usage;
}

if($version){
    die $re_dock6_version;
}

if ( !defined($recmol2) || !defined($box) || !defined($lib) ) {
    die $usage;
}

my $presets = 0;
foreach ($sparse, $quick, $thorough){
    $presets++ if ($_);
}

if ($presets > 1){
    die "$usage\nPlease choose no more than one of -e_sparse, -e_quick, or -e_thorough\n";
}

if ($sparse){
    $d_step = 2;
    $tether = 2.5;
} elsif ($quick){
    $e_min = 25;
    $d_step = 1.5;
    $tether = 2;
} elsif ($thorough){
    $n_rot = 4;
    $d_step = 0.75;
    $tether = 5;
}

# set up log file
my $logfile = "RE_DOCK6.log";
my $t = localtime;
my $faillog = "RE_DOCK6_FAILURES_" . $t->dmy . "_" . $t->hms . ".log";
open(my $LOG, ">", $logfile) or err("Failed to open $logfile for writing");
msg("INFO: Starting RE_DOCK6. Command = $command");
msg("INFO: Log file = $curdir/$logfile");

# TODO, more checks for inputs

# quick and dirty check that files exist
if (!(-s $recmol2)){
    msg("ERROR: no such file '$recmol2'");
    $die=1;
} else {
    msg("INFO: Receptor file = $recmol2");
}
if (!(-s $box)){
    msg("ERROR: no such file '$box'");
    $die=1;
} else {
    msg("INFO: Box file = $box");
}
if (!(-s $lib)){
    msg("ERROR: no such file '$lib'");
    $die=1;
} else {
    msg("INFO: Ligand(s) file = $lib");
}

# setup output directory
if ($outdir){
    if (!(-d $outdir)){
        msg("ERROR: no such directory '$outdir'");
        $die=1;
    }
} else {
    $outdir = "$curdir";
}
msg("INFO: Output directory = $outdir");

# set up temp directory
my $template = "tmp_re_dock6_XXXXXX";
if ($tmpdir){
    if (!(-d $tmpdir)){
        msg("ERROR: no such directory '$tmpdir'");
        $die=1;
    } else {
        $tmp_dir = tempdir ( $template, DIR => $tmpdir, CLEANUP => 1 );
    }
} else { 
    $tmp_dir = tempdir ( $template, DIR => $curdir, CLEANUP => 1 );
}
msg("INFO: Temp directory = $tmp_dir");


# set up for multithreading
my $available_threads = Thread::Semaphore->new($threads);
my $writing_to_out = Thread::Semaphore->new(1);

msg("INFO: Using $threads threads");

# check that dock is installed
my $dockexe = `dock6 -v 2>&1`;
if (!($dockexe)){
    msg("ERROR: dock6 doesn't appear to be installed or is not in the system PATH (using the command: \`dock6 -v 2>&1\`)");
    $die=1;
} else {
    if ($dockexe =~ /(DOCK v[\d\.]+)/){
        $dockexe = $1;
    } else {
        $dockexe = "Unsure of Dock version";
    }
    msg("INFO: Dock is installed, version = $dockexe");
}


# check that grid is installed
my $gridexe = `grid 2>&1`;
if (!($gridexe)){
    msg("ERROR: grid (part of UCSF dock6) doesn't appear to be installed or is not in the system PATH (using the command: \`grid 2>&1\`)");
    $die=1;
}
msg("INFO: Grid is installed");

# make sure we can find the parameters folder
if (!(-d "$Bin/../parameters")) {
    msg("ERROR: paramters folder not found using search '-d \$bin/../parameters'
    Add the RE_Dock6 bin folder to the system PATH (don't symlink), parameters folder should be ~/re_dock6/parameters");
    $die=1;
}
my $params = "$Bin/../parameters";

# all checks are finished, die if any issues
if ($die){
    err("ERROR: One or more errors encountered, exiting...\n");
}

#---IO SETUP---
$outfile = "$outfile.mol2" if ($outfile !~ /\.mol2$/);
if (-f "$outdir/$outfile"){
    if ($force){
        msg("INFO: Overwriting previous ouput file $outdir/$outfile");
        unlink "$outdir/$outfile";
    } else {
        msg("INFO: Out file $outdir/$outfile already exists, resuming previous docking run");
        read_previous_output();
    }
}

# settings

msg("PARAM: Ligand x, y, z rotations:             $n_rot");
msg("PARAM: Ligand placement step distance:       $d_step");
msg("PARAM: Tether distance from box center:      $tether");
msg("PARAM: Energy minimisation steps:            $e_min");
msg("PARAM: Scored conformers output per ligand:  $confs");


# prep for grid
my $prefix = $recmol2;
$prefix =~ s/\.mol2$//;

# run grid if the grid files are missing
if (!(-f "$prefix.nrg") || !(-f "$prefix.bmp")) {
    make_grid_files();
}

# open the library for reading
if ($lib =~ /\.gz$/) {
    open(LIGMOL2, "gunzip -c $lib |") or err("can't open pipe to $lib");
} else {
    open(LIGMOL2, "<", $lib) or err("could not open $lib for reading");
}

#---MAIN DOCKING LOOP---
while(<LIGMOL2>){
    if ($_ =~ /@<TRIPOS>MOLECULE/){
        my $ligname = <LIGMOL2>;
        $ligname =~ s/\s//g;
        if (!($cur_lig)){
            $cur_lig = $ligname;
        }
        if ($ligname ne $cur_lig){
            if ($skip_ligands{$cur_lig}){
                msg("INFO: Skipping docked ligand $cur_lig");
            } else {
                # will hold within prep_and_dock untill a thread frees up
                prep_and_dock();
            }
            $cur_lig = $ligname;
            $lig_number++;
            @LIG = ();
            push @LIG, "@<TRIPOS>MOLECULE\n$ligname\n";
        } else {
            push @LIG, "@<TRIPOS>MOLECULE\n$ligname\n";
        }
    } else {
        push @LIG, $_ if ($cur_lig);
    }
}
err("ERROR: Failed to get ligands from file $recmol2, check file format") if (!($cur_lig));
# dock the final ligand
if ($skip_ligands{$cur_lig}){
    msg("INFO: Skipping docked ligand $cur_lig");
} else {
    prep_and_dock();
}

# wait on remaining jobs
while(){
    if (threads->list()){
        sleep 0.5;
    } else {
        last;
    }
}

if (-s $faillog){
    msg("WARNING: Not all compounds docked successfully. This is common for large screens.
        Dock output for failed compunds is in $faillog");
}
msg("FINISHED: Docking pipeline has finished");

exit(0);

#---SUBROUTINES---

sub read_previous_output {
    msg("INFO: Reading previous re_dock6 out file and gathering ligand IDs to skip");
    open(my $PREV, "$outdir/$outfile") or err("failed to open $outdir/$outfile for reading");
    # TODO, run a check to make sure this is indeed a dock6 scored .mol2 format file
    while(<$PREV>){
        if ($_ =~ /@<TRIPOS>MOLECULE/){
            my $ln = <$PREV>;
            $ln =~ s/\s//g;
            $skip_ligands{$ln} = 1;
        }
    }
}

sub prep_and_dock {
    $available_threads->down(1);
    msg("INFO: Docking $cur_lig");
    write_lig($lig_number);
    prep_dock($lig_number);
    threads->create(\&dock_lig, $lig_number, $cur_lig);
}

sub make_grid_files {
    msg("INFO: Grid files not found for receptor structure, generating grid files");
    open(GRDT, "<", "$params/grid.in") or err("couldn't open $params/grid.in for reading");
    open(GRD, ">", "$tmp_dir/grid.in") or err("couldn't open $tmp_dir/grid.in for writing");
    while(<GRDT>){
        $_ =~ s/##RECEPTORMOL2##/$recmol2/;
        $_ =~ s/##BOXPDB##/$box/;
        $_ =~ s/##VDWDEFN##/$params\/vdw.defn/;
        $_ =~ s/##GRID##/$prefix/;
        print GRD $_;
    }
    close(GRDT);
    close(GRD);
    runcmd("grid -i $tmp_dir/grid.in -o $tmp_dir/grid.log");
}

sub write_lig{
    open(LIG, ">", "$tmp_dir/$_[0].tmp.mol2") or err("couldn't open $tmp_dir/$_[0].tmp.mol2 for writing");
    print LIG @LIG;
    close(LIG);
}

sub dock_lig {
    my $ln = $_[0];
    my $lname = $_[1];
    # run conf ensemble
    qrun("$Bin/mol2_conf_ensemble $recmol2 $box $tmp_dir/$ln.tmp.mol2 $tmp_dir/$ln.ensmbl.mol2 $n_rot $d_step $tether 1> $tmp_dir/$ln.log 2>&1");

    # run dock6 to score orientations
    qrun("dock6 -i $tmp_dir/dock.$ln.in 1> $tmp_dir/$ln.log 2>&1");

    # write output
    $writing_to_out->down(1);

    # check output
    if (!(-s "$tmp_dir/$ln\_ranked.mol2")){
        msg("WARNING: Docking of $lname appears to have failed, check $faillog for hints");
        qrun("cat $tmp_dir/$ln.log >> $faillog");
    }
    qrun("cat $tmp_dir/$ln\_ranked.mol2 >> $outdir/$outfile");
    $writing_to_out->up(1);

    # clean up temp files
    foreach my $file ("$tmp_dir/$ln.tmp.mol2", "$tmp_dir/$ln.ensmbl.mol2", 
                      "$tmp_dir/dock.$ln.in", "$tmp_dir/$ln\_ranked.mol2", "$tmp_dir/ln.log"){
        unlink $file;
    }
    
    msg("INFO: Finished docking $lname");
    $available_threads->up(1);
    threads->detach();
}

sub prep_dock {
    my $ln = $_[0];
    open(DOCKT, "<", "$params/dock.in") or err("ERROR: couldn't open $params/dock.in for reading");
    open(DOCK, ">", "$tmp_dir/dock.$ln.in") or err("ERROR: couldn't open $tmp_dir/dock.$ln.in for writing");
    while(<DOCKT>){
        $_ =~ s/##LIGANDSMOL2##/$tmp_dir\/$ln.ensmbl.mol2/;
        $_ =~ s/##GRID##/$prefix/;
        $_ =~ s/##EMIN##/$e_min/;
        $_ =~ s/##VDWDEFN##/$params\/vdw.defn/;
        $_ =~ s/##FLEXDEFN##/$params\/flex.defn/;
        $_ =~ s/##FLEXDRIVETBL##/$params\/flex_drive.tbl/;
        $_ =~ s/##LIGANDOUT##/$tmp_dir\/$ln/;
        $_ =~ s/##CONFS##/$confs/;
        print DOCK $_;
    }
    close(DOCKT);
    close(DOCK);
}


sub msg {
    my $t = localtime;
    my $line = "[ " . $t->dmy . " " . $t->hms . " ] @_\n";
    print STDERR $line;
    print $LOG $line;
}


sub err {
    msg(@_);
    msg("FAILURE: Docking did not finish successfully");
    exit(1);
}


sub runcmd {
    msg("RUNNING:", @_);
    system(@_)==0 or err("ERROR: Failed to run command:", @_);
}


sub qrun {
    system(@_)==0 or err("ERROR: Failed to run command:", @_, "\n    Check the log files in $tmp_dir");
}



