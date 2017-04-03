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
# Description: Script for generating a box.pdb file, intended for use with 
#              UCSF Dock 6, re_dock6 etc.
# 
################################################################################

use strict;
use warnings;
use Getopt::Long;
use List::Util qw( min max sum );

my $usage = "
Usage:
perl  lig2box.pl  -l  ligand.mol2  -p  5  -o  box.pdb [-e]

REQUIRED ARGUMENTS:
-lig / -l       This is either the known ligand in it's bound (crystal structure) 
                conformation, or an ensemble of high-scoring docked ligands, or 
                a .mol2 file containing only the active site residues.

-padding / -p   Extra padding in angstroms for the box for all 6 directions 
                (+/- x y z). Recomend 2-3 when using active site residues and 
                3-6 when using a ligand.

-out / -o       Specify the out file. 

OPTIONAL ARGUMENTS:
-ensemble / -e  Use this flag to tell the script to use all available entries in 
                the mol2 file, not just the first. Needed  when using an ensemble 
                of high-scoring docked poses.

";

my $lig;
my $box;
my $padding;
my $ensemble;

GetOptions (
            "lig=s" => \$lig,
            "out=s" => \$box,
            "padding=i" => \$padding,
            "ensemble" => \$ensemble
			)
or die("\n###\nError in command line arguments\n###\n", $usage);

if ( !($lig) || !($box) || !($padding) ){
    die "$usage\nError, missing command line argument(s)\n";
}

open(my $IN, $lig) or die "$usage\nError opening $lig for reading\n";
open(my $OUTPDB, ">$box") or die "$usage\nError opening $box for writing\n";

my $start = 0;
my @x_coords;
my @y_coords;
my @z_coords;

# read in all the atom coords in the mol2 file
ITR: while(<$IN>) {
    if ($_ =~ /@<TRIPOS>BOND/) {
        if ($ensemble) {
            $start = 0;
            next ITR;
        } else {
            last ITR;
        }
    }
    if ($start == 1) {
        $_ =~ s/^\s+//;
        my @line = split(/\s+/, $_);
        push @x_coords, $line[2];
        push @y_coords, $line[3];
        push @z_coords, $line[4];
        next ITR;
    }
    if ($_ =~ /@<TRIPOS>ATOM/) {
        $start = 1;
        next ITR;
    }
}
close($IN);

# get center of mass and dimensions with padding
my $x_cent = sprintf ("%.2f", sum(@x_coords)/@x_coords);
my $y_cent = sprintf ("%.2f", sum(@y_coords)/@y_coords);
my $z_cent = sprintf ("%.2f", sum(@z_coords)/@z_coords);
my $x_dim;
my $y_dim;
my $z_dim;

if ( abs($x_cent - min(@x_coords)) > abs($x_cent - max(@x_coords))) {
    $x_dim = (2 * abs($x_cent - min(@x_coords))) + (2 * $padding);
} else {
    $x_dim = (2 * abs(max(@x_coords)-$x_cent)) + (2 * $padding);
}
if ( abs($y_cent - min(@y_coords)) > abs($y_cent - max(@y_coords))) {
    $y_dim = (2 * abs($y_cent - min(@y_coords))) + (2 * $padding);
} else {
    $y_dim = (2 * abs(max(@y_coords)-$y_cent)) + (2 * $padding);
}
if ( abs($z_cent - min(@z_coords)) > abs($z_cent - max(@z_coords))) {
    $z_dim = (2 * abs($z_cent - min(@z_coords))) + (2 * $padding);
} else {
    $z_dim = (2 * abs(max(@z_coords)-$z_cent)) + (2 * $padding);
}

my $xmin = sprintf ("%.2f", ($x_cent - (0.5 * $x_dim)));
my $xmax = sprintf ("%.2f", ($x_cent + (0.5 * $x_dim)));
my $ymin = sprintf ("%.2f", ($y_cent - (0.5 * $y_dim)));
my $ymax = sprintf ("%.2f", ($y_cent + (0.5 * $y_dim)));
my $zmin = sprintf ("%.2f", ($z_cent - (0.5 * $z_dim)));
my $zmax = sprintf ("%.2f", ($z_cent + (0.5 * $z_dim)));

$x_dim = sprintf ("%.2f", $x_dim);
$y_dim = sprintf ("%.2f", $y_dim);
$z_dim = sprintf ("%.2f", $z_dim);

# write to file
print $OUTPDB "HEADER    CORNERS OF BOX
REMARK    CENTER (X Y Z)   $x_cent  $y_cent  $z_cent
REMARK    DIMENSIONS (X Y Z)   $x_dim  $y_dim  $z_dim\n";

print $OUTPDB (pack  ("A7A4A2A4A6A9A8A8A8A1",("ATOM", 1, " ", "DUA", "BOX", "1", $xmin, $ymin, $zmin, "1"))), "\n";
print $OUTPDB (pack  ("A7A4A2A4A6A9A8A8A8A1",("ATOM", 2, " ", "DUB", "BOX", "1", $xmax, $ymin, $zmin, "1"))), "\n";
print $OUTPDB (pack  ("A7A4A2A4A6A9A8A8A8A1",("ATOM", 3, " ", "DUC", "BOX", "1", $xmax, $ymin, $zmax, "1"))), "\n";
print $OUTPDB (pack  ("A7A4A2A4A6A9A8A8A8A1",("ATOM", 4, " ", "DUD", "BOX", "1", $xmin, $ymin, $zmax, "1"))), "\n";
print $OUTPDB (pack  ("A7A4A2A4A6A9A8A8A8A1",("ATOM", 5, " ", "DUE", "BOX", "1", $xmin, $ymax, $zmin, "1"))), "\n";
print $OUTPDB (pack  ("A7A4A2A4A6A9A8A8A8A1",("ATOM", 6, " ", "DUF", "BOX", "1", $xmax, $ymax, $zmin, "1"))), "\n";
print $OUTPDB (pack  ("A7A4A2A4A6A9A8A8A8A1",("ATOM", 7, " ", "DUG", "BOX", "1", $xmax, $ymax, $zmax, "1"))), "\n";
print $OUTPDB (pack  ("A7A4A2A4A6A9A8A8A8A1",("ATOM", 8, " ", "DUH", "BOX", "1", $xmin, $ymax, $zmax, "1"))), "\n";

print $OUTPDB "CONECT    1    2    4    5
CONECT    2    1    3    6
CONECT    3    2    4    7
CONECT    4    1    3    8
CONECT    5    1    6    8
CONECT    6    2    5    7
CONECT    7    3    6    8
CONECT    8    4    5    7";

close($OUTPDB);
