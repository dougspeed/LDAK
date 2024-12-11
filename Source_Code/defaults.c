/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

//////////////////////////

//Set random seeds, and a few parameters, such as those used by multiple modes and/or required for parsefiles.c

///////////////////////////

if(seed==-9999)	//seed random number generators
{
//set random seed
srand((unsigned)clock()+(unsigned)time(NULL)+(unsigned)getpid());

//set seed for normal generator
zigset_safe(rand());

//set seed for uniform generator
init_genrand(rand());
}
else	//use provided seed
{
srand(seed);
zigset_safe(seed);
init_genrand(seed);
}

#if MKL==1
mkl_set_num_threads(maxthreads);

omp_set_dynamic(0);     //stops number of threads varying
omp_set_num_threads(maxthreads);
#endif

#if MKL==2
//openblas_set_num_threads(maxthreads);
bli_thread_set_num_threads(maxthreads);

omp_set_dynamic(0);     //stops number of threads varying
omp_set_num_threads(maxthreads);
#endif

///////////////////////////

//fill bedzeros, bedones and bedtwos
for(i=0;i<256;i++)
{
bedzeros[i]=0;
bedones[i]=0;
bedtwos[i]=0;

for(j=0;j<4;j++)
{
switch((i >> 2*j) & 3)
{
case 3: bedzeros[i]++;break;
case 2: bedones[i]++;break;
case 0: bedtwos[i]++;
}
}
}

///////////////////////////

if(nonsnp==-9999)
{
nonsnp=0;
if(dtype==3||dtype==4||(dtype==5&&genprobs==1))
{printf("Warning, predictors are assumed to be biallelic SNP allele counts (i.e., take values within [0,2]); if this is not the case, you should use \"--SNP-data NO\"\n\n");}
}

if(num_subs==-9999){num_subs=0;}

if(hwestand==-9999)
{
if(nonsnp==0){hwestand=1;}
else{hwestand=0;}
}

if(encoding==-9999){encoding=1;}

if(kindetails==-9999){kindetails=1;}

if(num_regs==-9999){num_regs=0;}

if(kingz==-9999){kingz=0;}
if(kinraw==-9999){kinraw=0;}

if(diagonal==-9999){diagonal=0;}
if(discenv==-9999){discenv=0;}

if(allone==-9999){allone=1;}

if(strcmp(indhers,"blank")!=0){ignoreweights=1;power=-1;}

if(strcmp(probsfile,"blank")!=0){ignoreweights=1;power=-1;}

if(checkroot==-9999){checkroot=1;}

if(manypreds==-9999){manypreds=0;}
if(manysamples==-9999){manysamples=0;}

///////////////////////////

