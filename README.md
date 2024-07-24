# Welcome to the GitHub pages for LDAK

LDAK is our software for analyzing data from genome-wide association studies (GWAS). It includes state-of-the-art tools for association analysis, constructing prediction models and heritability analysis.

# How to obtain LDAK

You can download the LDAK executable by clicking on the file links above (this will open a new page, then you can click the download button on the right). If you plan to use LDAK on a Linux computer, please download the file ldak6.linux; if you plan to use LDAK on a MAC, please download the file ldak5.2.mac (we will soon release a MAC version of LDAK6). There is also a beta Linux version of LDAK, however, please only use this file if asked by me or Isabella.

Alternatively, you can install LDAK via conda, using the following command

```
#command to come here
```

The final option is to compile LDAK yourself. You should first download the source code folder (either use git clone, or click on the green button that says "Code", then download zip). You can then compile LDAK using a command such as

```
gcc -O3 -o ldak6 ldak.c libqsopt.linux.a -lblas -llapack -lm -lz -fopenmp
chmod a+x ldak6
```

# LDAK documentation

For advice on running LDAK-KVIK, please visit www.ldak-kvik.com

For documentation on all other features of LDAK, please see www.dougspeed.com

# Asking for help

If you have problems running LDAK, we prefer if you open an issue at https://github.com/dougspeed/LDAK/issues. However, you can instead email doug <at> qgg <dot> au <dot> .dk.
