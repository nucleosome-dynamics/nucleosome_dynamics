# Install

For running the latest 'Nucleosome Dynamics CLI' version, simply clone the master branch of the present repository. It contains the R scripts in `/bin` and `/statistics` folders.

```sh
git clone http://mmb.irbbarcelona.org/gitlab/NuclDynamics/nucleosome_dynamics_CLI 
``` 

Some functionalities of these scripts depend on third-party software that needs an independent installation:

- 1. nucleR R package
- 2. NucDyn R package
- 3. UCSC wig utils


<a name="nucleR"></a>
#### 1. nucleR R package
Nucleosome positioning predictions in 'Nucleosome Dynamics CLI' relies on the specific methods implemented in 'nucleR' package. There are several options for installing nucleR. Check them all [here](http://mmb.pcb.ub.es/gitlab/NuclDynamics/nucleR).

<a name="NucDyn"></a>
#### 2. NucDyn R package
Nucleosome architecture comparison in 'Nucleosome Dynamics CLI' depend on the specific methods implemented in 'NucDyn' package. There are also several options for installing NuclDyn. Check them all [here](http://mmb.pcb.ub.es/gitlab/NuclDynamics/NucDyn).

<a name="UCSC"></a>
#### 3. UCSC wig utils
WIG format conversions (wig, bigwig) in 'Nucleosome Dynamics CLI' relies on the stand-alone binaries maintainced by [UCSC](https://genome.ucsc.edu/goldenPath/help/bigWig.html). For installing them, download the executables of interest (`wigToBigWig`, `bigWigToWig` and `fetchChromSizes`) for your operating system from the official binary repository. For a Linux 64-bit machine do:

```sh
mkdir wig_utils;
cd wig_utils;
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig &&  \ 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig  && \ 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
```

Make sure you save them in a directory included in your **$PATH**, or simply add your installation path directory to it by editing in your `.bashrc` home directory the following line:

```sh
export PATH="$PATH:/absolute/path/to/wig_utils"
```
