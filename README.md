# Welcome to the LDAK GitHub Pages

Please note that this page focuses on how to download LDAK; if you are instead looking for advice on how to run LDAK, please visit either www.ldak-kvik.com (for documentation related to LDAK-KVIK) or www.dougspeed.com (for documentation on all other features of LDAK).

If you find this page very confusing, you may prefer to instead follow the instructions on www.dougspeed.com/two-simple-analyses.

# How to obtain LDAK

There are three ways to obtain LDAK; you can download the LDAK executable, you can install LDAK via conda, or you can compile your own version of LDAK from the source code. We generally recommend the first way (i.e., downloading the LDAK executable).

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

3 - The final option is to compile LDAK yourself. You should first download the source code from the repository at the top of this page (either use git clone, or press the green button that says "Code", then select "Download zip"). You can then compile LDAK using a command such as

```
gcc -O3 -o ldak6.1 ldak.c libqsopt.linux.a -lblas -llapack -lm -lz -fopenmp
chmod a+x ldak6.1
```

# LDAK documentation

For advice on running LDAK-KVIK, please visit www.ldak-kvik.com

For documentation on all other features of LDAK, please see www.dougspeed.com

# Asking for help

If you have problems running LDAK, we prefer if you open an issue at https://github.com/dougspeed/LDAK/issues. However, you can instead email doug \<at\> qgg \<dot\> au \<dot\> .dk.

# Contributors

These pages and the LDAK software are maintained and developed by Doug Speed and Jasper Hof.
