# RE_DOCK6: A Rigid Exhaustive extension for UCSF Dock 6 #

RE_DOCK6 is an extension for UCSF Dock 6 that allows it to dock in a rigid exhaustive manner. While Dock already has a method for rigid ligand docking, this method is much faster and doesn't require the generation of target spheres.

RE_DOCK6 is intended to be easy to use; it takes in as little as a prepared receptor structure, prepared ligand file to dock, and a box file. The best scoring pose(s) are output in standard Dock .mol2 format. There is multi-threading support and a number of configurable parameters for conformation ensemble generation.

```
re_dock6.pl  -r  <rec.mol2>  -b  <box.pdb>  -l  <ligand(s).mol2>  [other options]
```

Check the WIKI for detailed usage.

### IMPORTANT NOTE ###

UCSF Dock 6 is NOT included in this git. The MIT licence provided with this git covers the source code in this git only and does NOT cover any UCSF Dock6 files. Before you can use RE_DOCK6 you will need to obtain a UCSF Dock6 license, and download and install Dock6.

[ UCSF Dock 6 URL ](http://dock.compbio.ucsf.edu/DOCK_6/index.html)

### HOW TO CITE ###
If RE_DOCK6 has helped in your work then you can cite it like so, but please don't forget to cite UCSF Dock6 (and Chimera if you use it).

**RE_DOCK6:**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.439792.svg)](https://doi.org/10.5281/zenodo.439792)

**UCSF Dock6:**

W. J. Allen, T. E. Balius, S. Mukherjee, S. R. Brozell, D. T. Moustakas, P. T. Lang, D. A. Case, I. D. Kuntz, R. C. Rizzo. J. Comput. Chem. 2015, 36, 1132–1156. [ DOI: 10.1002/jcc.23905 ](http://dx.doi.org/10.1002/jcc.23905)

**UCSF Chimera:**

[ Citing UCSF Chimera ](https://www.cgl.ucsf.edu/chimera/docs/credits.html)


## INSTALLATION ##

### Installing UCSF Dock 6 (Ubuntu/Debian) ###

* Install dependencies (if missing)

```bash
# dependencies you probably already have #
sudo apt-get install make g++

# less common dependencies #
sudo apt-get install flex gfortran bison
```

* Download the source code then extract, configure and compile

```bash
tar xzvf dock.6.7_source.tar.gz
cd dock6/install
./configure gnu
make all
```

* Add the Dock6/bin directory to your PATH, for instance by adding a line to .bashrc

`PATH=$PATH:/path/to/dock6/bin`

### Installing RE\_DOCK6 (Ubuntu/Debian) ###

* Clone or download and extract this git.

* Compile the conformation ensemble generator using your favorite Fortran compiler (I've currently only tested gfortran on Ubuntu, there is a simple makefile in re_dock6/src).

```bash
cd RE_DOCK6/src
# use the makefile (edit for different Fortran compiler)
make
make install
make clean

# or you can do it manually
# compile
gfortran -O3 -o ../bin/mol2_conf_ensemble boxio.f95 mol2io_min.f95 mol2_conf_ensemble.f95
# cleanup
rm -f *.o *.mod *.MOD
```

* Add the bin directory to your PATH, again you can modify your .bashrc

`PATH=$PATH:/path/to/RE_DOCK6/bin`

That's it! The scripts don't use any uncommon Perl modules but if your Linux distro is missing them just install them via cpan or cpanminus. `cpan <modulename>` / `cpanm <modulename>`. You may need to open a new terminal before dock6, grid, lig2box.pl, and re_dock6.pl will be in your PATH.


## Other Required Software ###
* [ UCSF Chimera ](https://www.cgl.ucsf.edu/chimera/)

Chimera is a protein structure viewer/editor, used here for preparing the receptor protein structure and ligand files for use with Dock 6.
To install, just download and run the installer, it's very easy and is my personal favourite protein viewer.

* [ openbabel ](https://github.com/openbabel/openbabel.git)

Openbabel is used for generating multiconformers of ligands for use with screening. If you have another program you would prefer to use for this step then openbabel would not be needed, but you may have compatibility issues. RE_DOCK6 uses the molecule names in the .mol2 ligand file for grouping multi-conformers; multi-conformers will need to have the same name or else RE_DOCK6 will treat them as different ligands.

### Making a local install of openbabel (Ubuntu/Debian) ###

```bash
sudo apt-get install libeigen3-dev

git clone https://github.com/openbabel/openbabel.git
mv openbabel openbabel-build
mkdir openbabel
cd openbabel-build
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../openbabel
make -j2
make install
cd ..

# clean up #
rm -rf openbabel-build

# add to PATH via .bashrc #
cd openbabel/bin
echo "PATH=\$PATH:$PWD" >> $HOME/.bashrc
```
