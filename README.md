# Welcome to the LDAK GitHub Pages

Please note that this page focuses on how to download LDAK; if you are instead looking for advice on how to run LDAK, please visit either www.ldak-kvik.com (for documentation related to LDAK-KVIK) or www.dougspeed.com (for documentation on all other features).

If you find this page confusing, you may prefer to instead follow the instructions on www.dougspeed.com/two-simple-analyses.

# How to obtain LDAK

There are three ways to obtain LDAK; you can download the LDAK executable, you can install LDAK via conda, or you can compile your own version of LDAK from the source code. We generally recommend the first option (downloading the LDAK executable), while the third option (compiling from source code) is technically challenging for those without a computer science background.

Please note that the Linux LDAK executable is most regularly updated (by contrast, the Mac executable, conda version and source code tend to be a few months out of date).


# 1A - Download the Linux executable:
(note that if using a Mac, you should instead follow the instructions in 1B)

If you plan to run LDAK on a Linux computer, please download the file ldak6.1.linux; you can either click on the name of the file at the top of this page, then find the download button, or you can use the following command

```
wget https://github.com/dougspeed/LDAK/raw/main/ldak6.1.linux
```
# 1B - Download the Mac executable:
(note that if using Linux, you should instead follow the instructions in 1A)

If you plan to run LDAK on a Mac, please download the file ldak6.1.mac; you can either click on the name of the file at the top of this page, then find the download button, or you can use the following command

```
wget https://github.com/dougspeed/LDAK/raw/main/ldak6.1.linux
```
You may have noticed there is also a beta Linux version of LDAK, but please only use this version if asked.

# 2 - Install LDAK via conda (Linux systems only):

```
# create a new environment and install ldak6
conda create -n ldak_env -c genomedk ldak6

# load the new environment
conda activate ldak_env

# check ldak runs if you type the following
ldak6
```
# 3A - Compile a Linux version of LDAK from source code:
(note that if using a Mac, you should instead follow the instructions in 3B)

Please download and extract the file Source_Code.zip (avialable at the top of this page). If you open a terminal window, and navigate to the Source_Code folder, then you can hopefully compile LDAK using a command such as

```
gcc -O3 -o ldak6.1 ldak_nomkl.c libqsopt.linux.a -lblas -llapack -lm -lz
```
Sometimes it is necessary to add the option "-no-pie", while adding "-Wformat-overflow=0" reduces the warnings. 
If succcessful, you should have created a file called ldak6.1, which you can then run by typing
```
chmod a+x ldak6.1
./ldak6.1
```
Please note that this version of LDAK will likely be slower than the pre-compiled Linux version of LDAK, because it does not utilize the Intel MKL libraries (see 3C for more details).

# 3B - Compile a Mac version of LDAK from source code:
(note that if using Linux, you should instead follow the instructions in 3A)

Please download and extract the file Source_Code.zip (avialable at the top of this page). If you open a terminal window, and navigate to the Source_Code folder, then you can hopefully compile LDAK using a command such as

```
gcc -O3 -o ldak6.1 ldak_nomkl.c libqsopt.mac.a -lblas -llapack -lm -lz
```
If succcessful, you should have created a file called ldak6.1, which you can then run by typing
```
chmod a+x ldak6.1
./ldak6.1
```
Note that you may first have to allow your computer to execute non-Apple software. You can do this by opening "System Settings", then clicking on "Privacy & Security", followed by "Developer Tools", and ticking the box that allows apps downloaded from "Anywhere".

# 3C - Compile a Linux version of LDAK that utilizes the Intel MKL libraries:

The pre-compiled Linux version of LDAK uses the Intel MKL libraries, which provide highly-efficient algebraic routines (e.g., for multiplying and decomposing matrices). I currently use the libraries within the 2024 version of oneAPI; the latest version of oneAPI can be downloaded free from  
https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html

I then compile LDAK using the commands

```
source intel/oneapi/setvars.sh
gcc --static -static-libgcc -O3 -o ldak6.1.linux ldak/ldak_mkl.c ldak/libqsopt.linux.a -m64 -Wl,--start-group ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_gnu_thread.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -lz -I${MKLROOT}/include -fopenmp
```
Note that the first command depends on where you installed oneAPI. Further, I generated the second command based on help from the Intel MKL Link Advisor (https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html).

# LDAK documentation

For advice on running LDAK-KVIK, please visit www.ldak-kvik.com

For documentation on all other features of LDAK, please see www.dougspeed.com

# Asking for help

If you have problems running LDAK, we prefer if you open an issue at https://github.com/dougspeed/LDAK/issues. However, you can instead email doug \<at\> qgg \<dot\> au \<dot\> .dk.

# Contributors

These pages and the LDAK software are maintained and developed by Doug Speed and Jasper Hof.
