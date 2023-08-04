# Precompiled installation packages
_note this list is not comprehensive, is not meant to validate nor endorse a particular third-party build over others and may not always lead to the newest version_

## Multiple operating systems
Binaries designed to be compatible with [Julia](https://github.com/JuliaLang/julia) are regularly provided in
https://github.com/staticfloat/OpenBLASBuilder/releases . Note that their 64bit builds with INTERFACE64=1 have
`_64` appended to the function symbol names, so your code needs to refer to e.g. gemm_64_ rather than gemm_

the [Conda-Forge](https://github.com/conda-forge) project maintains packages for the conda package manager at
https://github.com/conda-forge/openblas-feedstock  
## FreeBSD
 > pkg install openblas

 see https://www.freebsd.org/ports/index.html
## Linux
### Debian/Ubuntu/Mint/Kali
 OpenBLAS package is available in default repositories and can act as default BLAS in system

Example installation commands:
```
  $ sudo apt update
  $ apt search openblas
  $ sudo apt install libopenblas-dev
  $ sudo update-alternatives --config libblas.so.3
```
 Alternatively, if distributor's package proves unsatisfactory, you may try latest version of OpenBLAS, [Following guide in OpenBLAS FAQ](https://github.com/xianyi/OpenBLAS/wiki/faq#debianlts)
 
### openSuSE/SLE: 
 Recent OpenSUSE versions include OpenBLAS in default repositories and also permit OpenBLAS to act as replacement of system-wide BLAS.

 Example installation commands:
```
  $ sudo zypper ref
  $ zypper se openblas
  $ sudo zypper in openblas-devel
  $ sudo update-alternatives --config libblas.so.3
``` 
Should you be using older OpenSUSE or SLE that provides no OpenBLAS, you can attach optional or experimental openSUSE repository as a new package source to acquire recent build of OpenBLAS following [instructions on openSUSE software site](https://software.opensuse.org/package/openblas)

### Fedora/CentOS/RHEL
Fedora provides OpenBLAS in default installation repositories.

To install it try following:
```
  $ dnf search openblas
  $ dnf install openblas-devel
```
For CentOS/RHEL/Scientific Linux packages are provided via [Fedora EPEL repository](https://fedoraproject.org/wiki/EPEL)

After adding repository and repository keys installation is pretty straightforward:
```
  $ yum search openblas
  $ yum install openblas-devel
```
No alternatives mechanism is provided for BLAS, and packages in system repositories are linked against NetLib BLAS or ATLAS BLAS libraries. You may wish to re-package RPMs to use OpenBLAS instead [as described here](https://fedoraproject.org/wiki/How_to_create_an_RPM_package)

### Mageia
Mageia offers ATLAS and NetLIB LAPACK in base repositories.
You can build your own OpenBLAS replacement, and once installed in /opt
TODO: populate /usr/lib64 /usr/include accurately to replicate netlib with update-alternatives

### Arch/Manjaro/Antergos
```
  $ sudo pacman -S openblas
```

## Solaris (via pkgsrc)


## OSX
https://www.macports.org/ports.php?by=name&substr=openblas

`brew install openblas`

or using the conda package manager from
https://github.com/conda-forge/miniforge#download
(which also has packages for the new M1 cpu)

 `conda install openblas`

## Windows
http://sourceforge.net/projects/openblas/files

https://www.nuget.org/packages?q=openblas
