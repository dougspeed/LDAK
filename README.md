# Welcome to the GitHub pages for LDAK

LDAK is our software for analyzing data from genome-wide association studies (GWAS). It includes state-of-the-art tools for association analysis, constructing prediction models and heritability analysis. 

# How to obtain LDAK

1 - The easiest way to obtain LDAK is by downloading the executable from the repository at the top of this page.

If you plan to run LDAK on a Linux computer, please download the file ldak6.linux; you can either click on the link, then find the download button, or you can use the following command

```
wget https://github.com/dougspeed/LDAK/blob/main/ldak6.linux
```

If you plan to run LDAK on a MAC, please download the file ldak6.mac; you can either click on the link, then find the download button, or you can use the following command

```
wget https://github.com/dougspeed/LDAK/blob/main/ldak6.mac
```
Note that there is also a beta Linux version of LDAK, but please only use this version if asked.

2 - Alternatively, you can install LDAK via conda, using the following command

```
#command to come here
```

3 - The final option is to compile LDAK yourself. You should first download the source code from the repository at the top of this page (either use git clone, or press the green button that says "Code", then select "Download zip"). You can then compile LDAK using a command such as

```
gcc -O3 -o ldak6 ldak.c libqsopt.linux.a -lblas -llapack -lm -lz -fopenmp
chmod a+x ldak6
```

# LDAK documentation

For advice on running LDAK-KVIK, please visit www.ldak-kvik.com

For documentation on all other features of LDAK, please see www.dougspeed.com

# Asking for help

If you have problems running LDAK, we prefer if you open an issue at https://github.com/dougspeed/LDAK/issues. However, you can instead email doug \<at\> qgg \<dot\> au \<dot\> .dk.

# Contributors

These pages and the LDAK software are maintained and developed by Doug Speed and Jasper Hof.
