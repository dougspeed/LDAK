//To-do list:
//Can not do spa with dosage data

//Hmm, must be something!

//Some notes:

//we assume resp, covars, X, are not large enough to require size_t when indexing
//normal (two-sided) pvalue for T is erfc(|T| root(1/2)), chisq pvalue for S is erfc(S^.5 root(1/2))
//cdf(x)=1-.5*erfc(x root(1/2))=.5*erfc(-x root(1/2))

//reml integrates out fixed effs, then solve vars, then given vars, finds maximising fixed + random effs
//when vars fixed, integrating out fixed effects same (up to constant) as setting fixed effects to mle
//therefore random effects are same whether obtain with fixed effects set to mle or integrated out
//and vice versa, fixed effects same whether obtain with random effects set to mle or integrated out

/*
Copyright 2024 Doug Speed.
 
    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.
*/

///////////////////////////

//Main file - see mode.c for details of the different modes

///////////////////////////

//The pre-compiled Linux version of LDAK uses Intel MKL libraries; I compile it using 
//source intel/oneapi/setvars.sh
//gcc -O3 -static -o ldak6.linux source/ldak.c source/libqsopt.linux.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -lz -m64 -I${MKLROOT}/include -fopenmp -L/home/doug/opt/lib -I/home/doug/opt/include -L/home/doug/opt/glibc/lib -I/home/doug/opt/glibc/include

////////

//If you wish to compile yourself, but do not have Intel Libraries, you must change the variable MKL to zero (look a few lines below)

//On my personal (LINUX) laptop I can compile a dynamic version using
//gcc source/ldak.c source/libqsopt.linux.a -o ldak.out -lblas -llapack -lm -lz

//for a static version, consider something like
//gcc source/ldak.c source/libqsopt.linux.a -o ldak.out -lc libraries/liblapack.a libraries/libblas.a libraries/libc.a libraries/libgfortran.a libraries/libm.a libraries/libz.a -static

//On a MAC, compile using 
//gcc source/ldak.c source/libqsopt.mac.a -o ldak.out -lblas -llapack -lm -lz
//you may have to add --framework accelerate to this this command, and/or first run xcode-select --install

///////////////////////////

#define MET 0	//0 for qsopt, 1 for glpk (for glpk, edit compilation line to include source/glpk.h)
#define MKL 1	//1 to compile with mkl, 0 to compile without mkl, 2 to compile with AOCL
//(changing MKL changes one line below and a few in defaults.c)

//library includes

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <zlib.h>
#include <limits.h>

#if MET==0
#include"qsopt.h"
#else
#include"glpk.h"
#endif

#if MKL==1
#include<mkl.h>
#include <omp.h>
#include <sys/sysinfo.h>
#endif

#if MKL==2
#include <omp.h>
#include <blis.h>
#endif

//will export functions from fortran

extern void dgemm_();
extern void sgemm_();
extern void dgemv_();
extern void sgemv_();
extern void dsymm_();
extern void ssymm_();
extern void dsymv_();
extern void ssymv_();
extern void dpotrf_();
extern void spotrf_();
extern void dpotri_();
extern void spotri_();
extern void dpotrs_();
extern void spotrs_();
extern void dsytrf_();
extern void ssytrf_();
extern void dsytri_();
extern void ssytri_();
extern void dsytrs_();
extern void ssytrs_();
extern void dsyev_();
extern void ssyev_();
extern void dsyevx_();
extern void dgesvd_();
extern void dgesdd_();
extern void dtrmm_();
extern void strmm_();
extern void dtrmv_();
extern void strmv_();

extern void sspmv_();
extern void spptrf_();
extern void spptri_();
extern void spptrs_();

struct sorting_double{double value;int index;};
struct sorting_string{char *ptr;int index;};

//these defines are used in nnls.c
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define ABS(x) ((x) >= 0 ? (x) : -(x))

#include "mt19937.c"
#include "norm.c"
#include "gamma.c"
#include "nm_gamma.c"
#include "nnls.c"
#include "toms462.c"

#include "permute.c"
#include "oddsnends.c"
#include "sort.c"
#include "decomp.c"
#include "regfuns.c"
#include "cgd.c"
#include "spa.c"

#include "filedata.c"
#include "dataops.c"
#include "filemain.c"
#include "filekins.c"

#include "weightfuns.c"
#include "kinfuns.c"
#include "genefuns.c"

#include "remlhe.c"
#include "remlmulti.c"
#include "remlfast.c"
#include "blup.c"
#include "he.c"
#include "family.c"

#include "remllinear.c"
#include "remlgene.c"
#include "remladv.c"

#include "sumfuns.c"


//make break point for when using gbd
int breakfun()
{
return(0);
}

///////////////////////////

int main (int argc, const char * argv[])
{
//this line makes the buffer flush
setlinebuf(stdout);

if(MKL==0){printf("Please note that this MAC version is not compiled with MKL, so will be slower than the LINUX version\n\n");}

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\nLDAK - Software for obtaining Linkage Disequilibrium Adjusted Kinships and Loads More\nVersion 6 - Help pages at http://www.ldak.org\n-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n\n");

if(sizeof(unsigned char)!=1||sizeof(unsigned short)!=2||sizeof(unsigned)!=4||sizeof(float)!=4||sizeof(double)!=8||sizeof(size_t)!=8)
{printf("Error, this compilation uses %zd / %zd / %zd / %zd / %zd / %zd bytes to save an unsigned char / unsigned short / unsigned integer / float / size_t (LDAK requires 1 / 2 / 4 / 4 / 8 / 8)\n", sizeof(unsigned char), sizeof(unsigned short), sizeof(unsigned), sizeof(float), sizeof(double), sizeof(size_t));exit(1);}

//declare variables
#include "declare.c"

//deal with command line arguments (come in pairs)
#include "readargs.c"
//this sets mode - see modes.c for a list of modes
//it also takes care of mode=134 (decides whether it should become mode 131 or 132)

//append filenames and check exist
#include "append.c"

if(mode==134)	//work out whether switching to mode 131 or 132
{
#include "switcha.c"
}

if(mode==140)	//check not using multiple phenotypes
{
#include "switchb.c"
}

//check whether arguments are consistent
#include "consistent.c"

//read section, partition or gene details - will set a few more parameters
if(mode==102||mode==103||mode==104||mode==112||mode==113||mode==137||mode==138||mode==139||mode==192||mode==193||mode==194)
{
#include "readdetails.c"
}

//set random seeds and a few parameters - generally those used by multiple modes
#include "defaults.c"
      
//check data/regions/kinships etc
#include "parsefiles.c"

//mode-specific checks
#include "required.c"

//tell the user what will happen
#include "param.c"

///////////////////////////

if(use_data!=5)	//get num_samples_use, and get size of data
{
#include "getnums.c"

if(use_data==1)	//using data (normally)
{
printf("Data contain %d samples and %d predictors", num_samples, num_preds);
if(num_samples_use<num_samples||num_preds_use<num_preds)
{printf("; will be using %d and %d", num_samples_use, num_preds_use);}
printf("\n\n");
}

if(mode==131)	//deal with families, trios and duos
{
if(families==1)	//find families
{
#include "findfamilies.c"
}
if(trios==1)	//find trios (and also families)
{
#include "findtrios.c"
}
if(duos!=0)	//find duos (and also families)
{
#include "findduos.c"
}
}
}
//will deal with use_data=5 (merging) later

if(strcmp(idsfile,"blank")==0)	//without ids, cant use respfile, covarfile or envfile
{
flag=0;
if(strcmp(respfile,"blank")!=0)
{printf("Warning, phenotypes were provided but will not be used\n");strcpy(respfile,"blank");flag=1;}
if(strcmp(covarfile,"blank")!=0)
{printf("Warning, covariates were provided but will not be used\n");strcpy(covarfile,"blank");flag=1;}
if(strcmp(envfile,"blank")!=0)
{printf("Warning, environmentals were provided but will not be used\n");strcpy(envfile,"blank");flag=1;}
if(flag==1){printf("\n");}
}

///////////////////////////

if(strcmp(topfile,"blank")!=0)	//using tops - get annotations - already have tkeeppreds
{
tchr=malloc(sizeof(double)*num_tops);
tbp=malloc(sizeof(double)*num_tops);
tpreds=malloc(sizeof(char *)*num_tops);
tal1=malloc(sizeof(char)*num_tops);
tal2=malloc(sizeof(char)*num_tops);
for(j=0;j<num_tops;j++)
{
tchr[j]=allchr[tkeeppreds[j]];
tbp[j]=allbp[tkeeppreds[j]];
copy_string(tpreds,j,allpreds[tkeeppreds[j]]);
tal1[j]=allal1[tkeeppreds[j]];
tal2[j]=allal2[tkeeppreds[j]];
}
}	//end of using tops

if(num_regs>0)	//sort regions - must have modes 121/123/124
{
rkeeppreds=malloc(sizeof(int)*num_preds);
regindex=malloc(sizeof(int*)*num_regs);
rnum_preds_use=read_regions(regindex, rkeeppreds, num_regs, regpref, num_preds, allpreds, predorder, num_preds_use, keeppreds);
rpreds=malloc(sizeof(char*)*rnum_preds_use);
ral1=malloc(sizeof(char)*rnum_preds_use);
ral2=malloc(sizeof(char)*rnum_preds_use);
for(j=0;j<rnum_preds_use;j++)
{
copy_string(rpreds,j,allpreds[rkeeppreds[j]]);
ral1[j]=allal1[rkeeppreds[j]];
ral2[j]=allal2[rkeeppreds[j]];
}
}	//end of num_regs>0

///////////////////////////

if(use_data==1||use_data==6)	//reduce predictors to data_length and keeppreds_use
{
#include "setdl.c"
}

////////

if(use_data==1||use_data==6)	//(normal case) sort scalings and pvalues, perhaps reduce, then sort summaries
{
//some of these allocations are  unnecessary (especially for use_data=6)
centres=malloc(sizeof(double)*data_length);
mults=malloc(sizeof(double)*data_length);
sqdevs=malloc(sizeof(double)*data_length);
weights=malloc(sizeof(double)*data_length);
pvalues=malloc(sizeof(double)*data_length);

for(j=0;j<data_length;j++){centres[j]=-9999;mults[j]=-9999;sqdevs[j]=-9999;weights[j]=1;pvalues[j]=2;}

if(strcmp(centresfile,"blank")!=0)
{
printf("Reading centres from %s\n", centresfile);
read_centres(centresfile, centres, data_length, preds, al1, al2, bimfile);
printf("First few centres are: %s %.4e", preds[0], centres[0]);
for(j=1;j<data_length;j++){if(j<5){printf(" | %s %.4e", preds[j], centres[j]);}}
printf("\n\n");
}

if(strcmp(weightsfile,"blank")!=0)	//read weights from weightsfile
{
printf("Reading weights from %s\n", weightsfile);
read_weightsfile(weightsfile, weights, data_length, preds, keeppreds_use, num_preds, allpreds, bimfile);
printf("First few weights are: %s %.4f", preds[0], weights[0]);
for(j=1;j<data_length;j++){if(j<5){printf(" | %s %.4f", preds[j], weights[j]);}}
printf("\n\n");
}

if(strcmp(pvafile,"blank")!=0)
{
printf("Reading p-values from %s\n", pvafile);
read_pvafile(pvafile, pvalues, data_length, preds);
printf("First few p-values are: %s %.4e", preds[0], pvalues[0]);
for(j=1;j<data_length;j++){if(j<5){printf(" | %s %.4e", preds[j], pvalues[j]);}}
printf("\n\n");
}

if(strcmp(impfile,"blank")!=0)
{
printf("Reading importances from %s\n", impfile);
read_impfile(impfile, pvalues, data_length, preds);
printf("First few importances are: %s %.4e", preds[0], -pvalues[0]);
for(j=1;j<data_length;j++){if(j<5){printf(" | %s %.4e", preds[j], -pvalues[j]);}}
printf("\n\n");
}

if((mode==105&&reduce==1)||mode==112||mode==114||(mode==127&&reduce==1)||(mode==128&&reduce==1)||mode==136||mode==140||(mode==141&&reduce==1)||mode==151||mode==152||mode==153||mode==154||mode==159)	//can reduce
{
count=reduce_values(weights,data_length,keeppreds_use,chr,preds,cm,bp,cmbp,al1,al2,centres,pvalues,0,NULL);
if(count==0)
{printf("Error, there are no predictors with non-zero weight\n\n");exit(1);}
if(count<data_length)
{printf("By ignoring those with weight zero, number of predictors reduced to %d", count);
if(mode==105||mode==141){printf(" (to avoid this, use \"--reduce NO\")");}
printf("\n\n");
}
data_length=count;
}

if(strcmp(indhers,"blank")!=0)	//put relative per-snp heritabilities into weights and set her
{
read_indhers(indhers, weights, data_length, preds);
count=reduce_values(weights,data_length,keeppreds_use,chr,preds,cm,bp,cmbp,al1,al2,centres,pvalues,0,NULL);
if(count==0)
{printf("Error, there are no predictors with non-zero heritability\n\n");exit(1);}
if(count<data_length)
{printf("By ignoring those with non-positive heritability, number of predictors reduced to %d\n", count);}
data_length=count;

her=0;for(j=0;j<data_length;j++){her+=weights[j];}
for(j=0;j<data_length;j++){weights[j]=weights[j]/her;}

if(herscale!=-9999){her*=herscale;}

printf("Sum of per-predictor heritabilities: %.6f\n", her);
if(mode==151||mode==152||mode==153||mode==154)	//must control high values (ok for megaprs)
{
if(her>=1)
{printf("Warning, this must be less than one; will rescale pre-predictor heritabilities so that they sum to 0.9\n");her=0.9;}
if(her>.9)
{printf("Warning, very high values can affect convergence; will rescale pre-predictor heritabilities so that they sum to 0.9\n");her=0.9;}
}
printf("\n");
}

if(strcmp(sumsfile,"blank")!=0)	//must be mode 131, 132, 138, 140, 158, 159 or 172
{
nss=malloc(sizeof(double)*data_length);
chis=malloc(sizeof(double)*data_length);
rhos=malloc(sizeof(double)*data_length);

printf("Reading summary statistics from %s\n", sumsfile);
read_sumsfile(sumsfile, nss, chis, rhos, data_length, preds, al1, al2, bimfile, amb, fixn, scaling, 0, -9999);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis[0], nss[0]);
for(j=1;j<data_length;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis[j], nss[j]);}}
printf("\n\n");
}
}	//end of use_data=1

///////////////////////////

if(num_resps_use>0)	//get resp (possibly a fake response)
{
resp=malloc(sizeof(double)*num_samples_use*num_resps_use);
respcounts=malloc(sizeof(int)*num_resps_use);
respinds=malloc(sizeof(int)*num_samples_use*num_resps_use);

if(strcmp(respfile,"blank")!=0)	//phenotype provided (if dichot=1, have already checked responses are binary)
{read_respfile(respfile, resp, num_samples_use, ids3, num_resps_use, keepresps, num_resps, respcounts, respinds, missingvalue, (mode==130||mode==132||prev!=-9999), pad, (mode==229||mode==230));}
else	//make a fake phenotype (must have num_resps_use=1, mode 115, 121, 122, 125, 138. 140 or 172)
{
if(mode==121||mode==138||mode==140)	//will make a realistic phenotype (not entirely sure necessary)
{
if(prev!=-9999)	//binary - will have ascer
{
value=round(ascer*num_samples_use);
for(i=0;i<num_samples_use;i++){resp[i]=(i<value)-value/num_samples_use;}
}
else	//quantitative - use normal quantiles
{
for(i=0;i<num_samples_use;i++){value=(double)(i+1)/(num_samples_use+1);resp[i]=normal_inv(value);}
}

respcounts[0]=num_samples_use;
}
else	//make blank one
{
for(i=0;i<num_samples_use;i++){resp[i]=missingvalue;}
respcounts[0]=0;
}
}	//end of no phenotypes provided
}	//end of using resp

if(mode==229||mode==230)	//read second phenotype
{
keepresps2=malloc(sizeof(int)*1);
resp2=malloc(sizeof(double)*num_samples_use);
respcounts2=malloc(sizeof(int)*1);
respinds2=malloc(sizeof(int)*num_samples_use);

keepresps2[0]=mpheno2-1;
read_respfile(respfile, resp2, num_samples_use, ids3, 1, keepresps2, num_resps, respcounts2, respinds2, missingvalue, (mode==230), pad, 2);
}

////////

if(num_fixed>0)	//using covar - modes 121/122/123/124/125/126/127/128/130/131/132/133/138/140/145/151/152/153/154/169/170/172/173/175/194
{
covar=malloc(sizeof(double)*num_samples_use*num_fixed);

if(num_covars>0)
{
//first covariate is intercept
for(i=0;i<num_samples_use;i++){covar[i]=1;}

if(strcmp(covarfile,"blank")!=0)	//normal covariates
{read_covarfile(covarfile, covar+num_samples_use, num_samples_use, ids3, num_covars-1, missingvalue);}
}

if(strcmp(envfile,"blank")!=0)	//environmental variables
{read_envfile(envfile, covar+num_covars*num_samples_use, num_samples_use, ids3, num_envs, discenv, ids1, ids2);}

if(strcmp(topfile,"blank")!=0)	//have tops
{
tcentres=malloc(sizeof(double)*num_tops);
tvars=malloc(sizeof(double)*num_tops);
read_tops(datafile, covar+(num_covars+num_envs)*num_samples_use, num_samples_use, ids3, num_tops, tkeeppreds, tcentres, tvars, famfile, famhead, dtype, binary, genskip, genheaders, genprobs, bedbytes, num_preds, allpreds, missingvalue, nonsnp, sumsfile);

if(strcmp(sumsfile,"blank")!=0)	//read in summaries for tops - currently redudant, as not allowed
{
tnss=malloc(sizeof(double)*num_tops);
tchis=malloc(sizeof(double)*num_tops);
trhos=malloc(sizeof(double)*num_tops);

printf("Reading summary statistics for top predictors from %s\n", sumsfile);
read_sumsfile(sumsfile, tnss, tchis, trhos, num_tops, tpreds, tal1, tal2, bimfile, amb, fixn, scaling, 0, -9999);
printf("First few stats and ns are: %s %.3f %.1f", tpreds[0], tchis[0], tnss[0]);
for(j=1;j<num_tops;j++){if(j<3){printf(" | %s %.3f %.1f", tpreds[j], tchis[j], tnss[j]);}}
printf("\n\n");
}
}

if(strcmp(cofile,"blank")!=0)	//mode must be 122, 125 or 172
{
printf("Reading coefficients for %d covariates from %s\n\n", num_fixed, cofile);
thetas=malloc(sizeof(double)*num_fixed);
if(countcols(cofile)==1){read_values(cofile,thetas,num_fixed,NULL,1,0,0);}
else{read_values(cofile,thetas,num_fixed,NULL,2,1,0);}
}
}

if(mode==132)	//deal with offsets
{
offsets=malloc(sizeof(double)*num_samples_use);

if(strcmp(offsetfile,"blank")!=0)	//read from file
{
read_offsetfile(offsetfile, offsets, num_samples_use, ids3, ids1, ids2);
}
else	//set to zero
{
for(i=0;i<num_samples_use;i++){offsets[i]=0;}
}
}

///////////////////////////

if(num_regs>0)	//deal with regions - sort scalings and summaries, then read data, set rdata_length
{
rcentres=malloc(sizeof(double)*rnum_preds_use);
rmults=malloc(sizeof(double)*rnum_preds_use);
rsqdevs=malloc(sizeof(double)*rnum_preds_use);
rweights=malloc(sizeof(double)*rnum_preds_use);

for(j=0;j<rnum_preds_use;j++){rcentres[j]=-9999;rmults[j]=-9999;rweights[j]=1;}

if(strcmp(weightsfile,"blank")!=0)	//read weights from weightsfile and reduce
{
printf("Reading weights for region predictors from %s\n", weightsfile);
read_weightsfile(weightsfile, rweights, rnum_preds_use, rpreds, rkeeppreds, num_preds, allpreds, bimfile);
printf("First few weights are: %s %.6f", rpreds[0], rweights[0]);
for(j=1;j<rnum_preds_use;j++){if(j<5){printf(" | %s %.6f", rpreds[j], rweights[j]);}}
printf("\n\n");

count=reduce_values(rweights,rnum_preds_use,rkeeppreds,NULL,rpreds,NULL,NULL,NULL,ral1,ral2, NULL, NULL, num_regs, regindex);
//reduce_values will report an error if count==0
if(count<rnum_preds_use)
{printf("By ignoring those with weight zero, the number of (unique) region predictors is reduced to %d\n\n", count);}
rnum_preds_use=count;
}

if(strcmp(sumsfile,"blank")!=0)	//read in summaries for regions (must all be present)
{
rnss=malloc(sizeof(double)*rnum_preds_use);
rchis=malloc(sizeof(double)*rnum_preds_use);
rrhos=malloc(sizeof(double)*rnum_preds_use);

printf("Reading summary statistics for %d region predictors from %s\n", rnum_preds_use, sumsfile);
read_sumsfile(sumsfile, rnss, rchis, rrhos, rnum_preds_use, rpreds, ral1, ral2, bimfile, amb, fixn, scaling, 0, -9999);

printf("First few stats and ns are: %s %.3f %.1f", rpreds[0], rchis[0], rnss[0]);
for(j=1;j<rnum_preds_use;j++){if(j<3){printf(" | %s %.3f %.1f", rpreds[j], rchis[j], rnss[j]);}}
printf("\n\n");
}

//read in regions data and standardize
rdata_warn(rnum_preds_use,num_samples_use);
rdata=malloc(sizeof(double)*num_samples_use*rnum_preds_use);

if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
(void)read_data_fly(datafile, dtype, rdata, NULL, num_samples_use, keepsamps, 0, rnum_preds_use, rkeeppreds, datainputgz, 0, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
if(binary==0){gzclose(datainputgz);}

stand_data(rdata, rcentres, rmults, rsqdevs, num_samples_use, rnum_preds_use, missingvalue, power, 0, hwestand, rweights, 1, rpreds);

//adjust regindex in order to remove trivial predictors, and possibly prune 
rdata_length=prune_regions(num_regs, regindex, rmults, rprune, rdata, num_samples_use);

if(mode==121&&shortcut==1)	//warn if too many region predictors
{
if(rdata_length>num_samples_use)
{printf("Warning, there are more region predictors (%d) than samples (%d); usually it is better to include them as kinships instead\n\n", rdata_length, num_samples_use);}
if(rdata_length>2000&&rdata_length<num_samples_use)
{printf("Warning, the number of region predictors is quite high (%d); consider using \"--region-prune\" to thin these or including them as kinships instead\n\n", rdata_length);}
}
}	//end of num_regs>0

///////////////////////////

//do some memory warnings now (because reading kins might take a while)
if(mode==121&&shortcut==0){anal_warn(2*num_samples_use, num_samples_use);}
if(mode==126){anal_warn(num_samples_use/2, num_samples_use);}
if(mode==133){anal_warn((2+memsave)*num_samples_use, num_samples_use);}

if(num_kins>0&&(mode==120||mode==121||mode==123||mode==124||mode==126||mode==131||mode==133||mode==138))	//deal with kinships - will store in Mkins or Mkins_single
{
if(memsave==0)	//for calc-sim-grms, he, pcgc and fast-reml don't need double precision
{
if(mode==120||mode==123||mode==124||mode==126){kin_warn(num_kins, num_samples_use, 0, 1);}
else{kin_warn(num_kins, num_samples_use, 0, 0);}
}

if(mode==120||mode==123||mode==124||mode==126){Mkins_single=malloc(sizeof(float *)*num_kins);}
else{Mkins=malloc(sizeof(double *)*num_kins);}

kintraces=malloc(sizeof(double)*num_kins);

if(memsave==0)	//read kinships (and get traces direct)
{
for(k=0;k<num_kins;k++)
{
if(mode==120||mode==123||mode==124||mode==126)
{
Mkins_single[k]=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kins(kinstems[k], NULL, Mkins_single[k], 1.0, num_samples_use, ids3, 0, maxthreads);
sum=0;for(i=0;i<num_samples_use;i++){sum+=Mkins_single[k][(size_t)i*num_samples_use+i];}
kintraces[k]=sum/num_samples_use;
}
else
{
Mkins[k]=malloc(sizeof(double*)*num_samples_use*num_samples_use);
read_kins(kinstems[k], Mkins[k], NULL, 1.0, num_samples_use, ids3, 0, maxthreads);
sum=0;for(i=0;i<num_samples_use;i++){sum+=Mkins[k][(size_t)i*num_samples_use+i];}
kintraces[k]=sum/num_samples_use;
}
}
}
else	//just get traces
{
for(k=0;k<num_kins;k++)
{
kintraces[k]=read_kin_trace(kinstems[k], num_samples_use, ids3, NULL, 0, diagonal, maxthreads);
printf("\n");
}
}

count=0;
for(k=0;k<num_kins;k++)
{
if(fabs(kintraces[k]-1)>.01)
{printf("Warning, Kinship %d has average diagonal value %.6f (usually this is 1)\n", k+1, kintraces[k]);count++;}
}
if(count>0){printf("\n");}
}

////////

if(num_kins==1&&((mode==121&&shortcut==1)||mode==131||mode==138))	//will be using a decomposition
{
eigen_warn(num_samples_use);
U=malloc(sizeof(double)*num_samples_use*num_samples_use);
E=malloc(sizeof(double)*num_samples_use);

if(strcmp(eigenfile,"blank")==0)	//decompose
{
if(memsave==0)	//are there cases when we decompose and must keep kinship? - I think not
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){U[(size_t)i2*num_samples_use+i]=Mkins[0][(size_t)i2*num_samples_use+i];}
}
}
else
{read_kins(kinstems[0], U, NULL, 1.0, num_samples_use, ids3, 0, maxthreads);}

printf("Performing an eigen-decomposition; this can take a while\nIf preferred, you can perform the analysis in two steps; first replace the main argument with \"--decompose\" (this will save the eigen-decomposition) then restart adding \"--eigen\"\n\n");

lwork=-1;
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

decomp_warn(lwork);
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug info %d\n\n", info);exit(1);}
free(work);
}
else	//read it from file
{read_eigens(eigenfile, U, E, num_samples_use, ids3, -9999, maxthreads);}
}

///////////////////////////
//start of work
///////////////////////////

if(mode==101)	//cut-weights
{
if(nothin==0)	//reduce data_length to thinned SNPs 
{
//usedpreds starts as -1 (can be considered) or -2 (filtered out)
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=-2;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds_use[j]]=-1;}

if(window_kb>25||window_length>25)	//do in two steps
{
if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("Will prune in two passes; first with windows of size 0.01cM, then %.4fcM\n\n", window_cm);}
else
{printf("Will prune in two passes; first with windows of size 10kb, then %.2fkb\n\n", window_kb);}
}
else
{printf("Will prune in two passes; first with windows of length 10, then %d\n\n", window_length);}
bitsize2=bitsize;window_kb2=window_kb;window_length2=window_length;
if(window_kb!=-9999){window_kb=10;}
if(window_length!=-9999){window_length=10;}
flag=1;
#include "prune.c"
bitsize=bitsize2;window_kb=window_kb2;window_length=window_length2;flag=2;
flag=2;
#include "prune.c"
}
else	//do in one
{
flag=0;
#include "prune.c"
}

free(usedpreds);
}
if(nothin==1)	//probably pointless, but print "duplicates"
{
sprintf(filename,"%sthin.dups",folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor Replacement\n"); 
for(j=0;j<num_preds_use;j++){fprintf(output, "%s %s\n", allpreds[keeppreds[j]], allpreds[keeppreds[j]]);}
fclose(output);
}
//for nothin=2, nothing to do

cut_sections(data_length, chr, cmbp, section_kb, section_length, buffer_kb, buffer_length, num_samples_use, folder, window_kb, window_length, window_cm, datafile, bimfile, extract, nothin);
}	//end of mode=101

if(mode==102||mode==104)	//calc-weights
{
#include "weights.c"
}

if(mode==103||mode==104)	//join-weights
{
join_sections(folder, num_sections, sstarts, sends, sstarts2, sends2, allpreds, num_preds, keeppreds, num_preds_use, keeppreds_use, spread);
}

if(mode==105)	//adjust-weights
{
#include "tagging.c"
}

////////

if(mode==106||mode==107)	//thin or thin-tops
{
//usedpreds starts as -1 (can be considered) or -2 (filtered out)
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=-2;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds_use[j]]=-1;}

if(window_kb>25||window_length>25)	//do in two steps
{
if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("Will prune in two passes; first with windows of size 0.01cM, then %.4fcM\n\n", window_cm);}
else
{printf("Will prune in two passes; first with windows of size 10kb, then %.2fkb\n\n", window_kb);}
}
else
{printf("Will prune in two passes; first with windows of length 10, then %d\n\n", window_length);}
bitsize2=bitsize;window_kb2=window_kb;window_length2=window_length;
if(window_kb!=-9999){window_kb=10;}
if(window_length!=-9999){window_length=10;}
flag=1;
#include "prune.c"
bitsize=bitsize2;window_kb=window_kb2;window_length=window_length2;flag=2;
flag=2;
#include "prune.c"
}
else
{
flag=0;
#include "prune.c"
}

free(usedpreds);
}	//end of mode=106/107

if(mode==108||mode==109)	//find-tags or remove-tags - each target SNP becomes a gene
{
#include "finding.c"
}

///////////////////////////

if(mode==111)	//cut-kins - make kinships partitions
{
cut_partitions(data_length, chr, preds, part_length, bychr, num_parts, partpref, checkpart, folder, datafile, bimfile, extract);
}

if(mode==112||mode==114)	//calculate kinships
{
if(strcmp(invsfile,"blank")==0)
{
#include "kinsa.c"
}
else
{
#include "kinsb.c"
}
}

if(mode==113)	//join-kins - have set kinstems (in consistent.c)
{
sprintf(outfile,"%skinships.all",folder);
manip_kins(outfile, num_kins, kinstems, kinsums, ids1, ids2, ids3, num_samples_use, kingz, kinraw, 0, NULL, maxthreads);
}

////////

if(mode==115)	//filter
{
#include "filter.c"
}

if(mode==116)	//add-grm
{
manip_kins(outfile, num_kins, kinstems, kinsums, ids1, ids2, ids3, num_samples_use, kingz, kinraw, 1, NULL, maxthreads);
}

if(mode==117)	//sub-grm
{
if(extract==1)
{
#include "subkin.c"
}	//end of subtracting predictors
else
{
manip_kins(outfile, num_kins, kinstems, kinsums, ids1, ids2, ids3, num_samples_use, kingz, kinraw, 2, NULL, maxthreads);
}
}

if(mode==118)	//convert-gz
{
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kinsgz(kinstems[0], kins_single, num_samples_use, ids3, partial);
write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 4);
free(kins_single);
}

if(mode==119)	//convert-raw
{
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kinsraw(kinstems[0], kins_single, num_samples_use, ids3);
write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 4);
free(kins_single);
}

////////

if(mode==120)	//computing correlations between kinships - must have single kins and memsave=0
{
cors=malloc(sizeof(double)*num_kins*num_kins);

for(k=0;k<num_kins;k++)
{
cors[k+k*num_kins]=1;
for(k2=k+1;k2<num_kins;k2++)
{
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;count=0;
for(i2=1;i2<num_samples_use;i2++)
{
for(i=0;i<i2;i++)
{
sum+=Mkins_single[k][(size_t)i2*num_samples_use+i];
sum2+=Mkins_single[k2][(size_t)i2*num_samples_use+i];
sumsq+=pow(Mkins_single[k][(size_t)i2*num_samples_use+i],2);
sumsq2+=pow(Mkins_single[k2][(size_t)i2*num_samples_use+i],2);
sumsq3+=Mkins_single[k][(size_t)i2*num_samples_use+i]*Mkins_single[k2][(size_t)i2*num_samples_use+i];
count++;
}
}
mean=sum/count;mean2=sum2/count;
value=(sumsq3-count*mean*mean2)/pow(sumsq-count*mean*mean,.5)/pow(sumsq2-count*mean2*mean2,.5);
cors[k+k2*num_kins]=value;cors[k2+k*num_kins]=value;
printf("Correlation between Matrices %d and %d is %.6f\n", k+1, k2+1, value);
}	//end of k2 loop
}
printf("\n");

sprintf(filename,"%s.cors.grms",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(k=0;k<num_kins;k++)
{
for(k2=0;k2<num_kins;k2++){fprintf(output,"%.10f ", cors[k+k2*num_kins]);}
fprintf(output,"\n");
}
fclose(output);

printf("Correlations saved in %s\n\n", filename);

free(cors);
}	//end of mode=120

///////////////////////////

if(mode==121)	//reml
{
//allocate - begin by guessing memory used by allocations within remlmulti.c
//if(shortcut==0){anal_warn(2*num_samples_use, num_samples_use);}

sweights=malloc(sizeof(double)*num_samples_use);
order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

if(num_regs>0)	//will be using Xs
{
X=malloc(sizeof(double)*num_samples_use*rdata_length);
Xstarts=malloc(sizeof(int)*num_regs);
Xends=malloc(sizeof(int)*num_regs);
Xsums=malloc(sizeof(double)*num_regs);
if(strcmp(sumsfile,"blank")!=0)
{
Xnss=malloc(sizeof(double)*rdata_length);
Xrhos=malloc(sizeof(double)*rdata_length);
}
}

if(strcmp(sampwfile,"blank")!=0)	//get sample weights
{read_sampwfile(sampwfile, sweights, num_samples_use, ids3, ids1, ids2);}
else
{
for(i=0;i<num_samples_use;i++){sweights[i]=1;}
}

for(m=0;m<num_resps_use;m++)	//loop through phenotypes
{
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Z (maybe permuted) and maybe X (not permuted) - will fill Y (maybe permuted) below
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}
if(strcmp(sumsfile,"blank")==0)
{Xtotal=fill_X(X, Xstarts, Xends, Xsums, NULL, NULL, num_samples_use, num_regs, regindex, rdata, NULL, NULL);}
else
{Xtotal=fill_X(X, Xstarts, Xends, Xsums, Xnss, Xrhos, num_samples_use, num_regs, regindex, rdata, rnss, rrhos);}

if(mpheno!=-1)	//single phenotype
{
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]];}
strcpy(filename,outfile);
}
else	//analysing phenotypes in turn
{
printf("Analysing Response %d of %d\n", m+1, num_resps_use);
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]+m*num_samples_use];}
sprintf(filename, "%s.%d", outfile, m+1);
}

if(num_kins==0&&num_regs>0&&shortcut==1)	//do quick version - might have summary stats
{
if(strcmp(sumsfile,"blank")==0)
{fast_reml(num_samples_use, num_covars, num_tops, num_envs, num_regs, Y, Z, X, Xtotal, Xstarts, Xends, Xsums, NULL, NULL, prev, respcounts[m], constrain, shrink, strip, tol, maxiter, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, rkeeppreds, regindex, rcentres, rmults, rweights, tkeeppreds, tcentres);}
else
{fast_reml(num_samples_use, num_covars, 0, 0, num_regs, Y, Z, X, Xtotal, Xstarts, Xends, Xsums, Xnss, Xrhos, prev, -9999, constrain, shrink, strip, tol, maxiter, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, rkeeppreds, regindex, rcentres, rmults, rweights, tkeeppreds, tcentres);}
}
else	//do full version
{multi_reml(num_samples_use, num_covars, num_envs, num_tops, num_kins, num_regs, Y, Z, Mkins, NULL, kintraces, kinsums, sweights, X, Xtotal, Xstarts, Xends, Xsums, prev, respcounts[m], hersfile, hestart, shortcut, U, E, 0, -9999, discenv, oversfile, ssums, constrain, tol, maxiter, memsave, maxthreads, kinstems, ids3, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, rkeeppreds, regindex, rcentres, rmults, rweights, tkeeppreds, tcentres, NULL);}
}	//end of m loop

free(sweights);free(order);free(Y);free(Z);
if(num_regs>0){free(X);free(Xstarts);free(Xends);free(Xsums);}
if(num_regs>0&&strcmp(sumsfile,"blank")!=0){free(Xnss);free(Xrhos);}
}

////////

if(mode==122)	//calc-blups
{
#include "timesa.c"
}

////////

if(mode==123||mode==124)	//he and pcgc regression
{
//allocate
order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

if(num_regs>0)	//will be using X
{
X=malloc(sizeof(double)*num_samples_use*rdata_length);
Xstarts=malloc(sizeof(int)*num_regs);
Xends=malloc(sizeof(int)*num_regs);
Xsums=malloc(sizeof(double)*num_regs);
}

if(num_blocks==-1||num_blocks>num_samples_use){num_blocks=num_samples_use;}

for(m=0;m<num_resps_use;m++)	//loop through phenotypes
{
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Z (maybe permuted) and maybe X (not permuted) - will fill Y (maybe permuted) below
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}
Xtotal=fill_X(X, Xstarts, Xends, Xsums, NULL, NULL, num_samples_use, num_regs, regindex, rdata, NULL, NULL);

if(mpheno!=-1)	//single phenotype
{
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]];}
strcpy(filename,outfile);
}
else	//analysing phenotypes in turn
{
printf("Analysing Response %d of %d\n", m+1, num_resps_use);
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]+m*num_samples_use];}
sprintf(filename, "%s.%d", outfile, m+1);
}

he_reg(num_samples_use, num_subs, num_covars, num_envs, num_tops, num_kins, num_regs, subindex, Y, Z, Mkins_single, kintraces, kinsums, X, Xtotal, Xstarts, Xends, Xsums, prev, respcounts[m], discenv, oversfile, ssums, num_blocks, memsave, kinstems, ids3, filename, (mode==124), trun);
}	//end of m loop

free(order);free(Y);free(Z);
if(num_regs>0){free(X);free(Xstarts);free(Xends);free(Xsums);}
}

////////

if(mode==125)	//reml-pred
{
#include "remlpred.c"
}

////////

if(mode==126)	//fast-reml - not using region nor summaries
{
//allocate - begin by guessing memory used by allocations within remlmulti.c
//anal_warn(num_samples_use/2, num_samples_use);

sweights=malloc(sizeof(double)*num_samples_use);
order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

if(strcmp(sampwfile,"blank")!=0)	//get sample weights
{read_sampwfile(sampwfile, sweights, num_samples_use, ids3, ids1, ids2);}
else
{
for(i=0;i<num_samples_use;i++){sweights[i]=1;}
}

for(m=0;m<num_resps_use;m++)	//loop through phenotypes
{
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Z (maybe permuted) - will fill Y (maybe permuted) below
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

if(mpheno!=-1)	//single phenotype
{
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]];}
strcpy(filename,outfile);
}
else	//analysing phenotypes in turn
{
printf("Analysing Response %d of %d\n", m+1, num_resps_use);
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]+m*num_samples_use];}
sprintf(filename, "%s.%d", outfile, m+1);
}

multi_reml(num_samples_use, num_covars, num_envs, num_tops, num_kins, num_regs, Y, Z, NULL, Mkins_single, kintraces, kinsums, sweights, NULL, 0, NULL, NULL, NULL, prev, respcounts[m], hersfile, hestart, 2, NULL, NULL, num_vects, ldlt, 0, oversfile, ssums, constrain, tol, maxiter, memsave, maxthreads, kinstems, ids3, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, NULL, NULL, NULL, NULL, NULL, tkeeppreds, tcentres, NULL);
}	//end of m loop

free(sweights);free(order);free(Y);free(Z);
if(num_regs>0){free(X);free(Xstarts);free(Xends);free(Xsums);}
}

////////

if(mode==127||mode==128)	//fast he and pcgc regression
{
#include "fast.c"
}

////////

if(mode==129||mode==130)	//estimation of family heritability
{
//allocate and fill response and covariates

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);
ids4=malloc(sizeof(char*)*num_samples_use);

for(m=0;m<num_resps_use;m++)	//loop through phenotypes
{
//index non-missing phenotypes
count=0;
for(i=0;i<num_samples_use;i++)
{
if(respinds[i+m*num_samples_use]==1){order[count]=i;copy_string(ids4,count,ids3[i]);count++;}
}
if(count==0){printf("Error, Phenotype %d has no non-missing phenotypes\n\n", m+1);exit(1);}
if(count<3){printf("Error, Phenotype %d has only %d non-missing phenotypes; unable to continue\n\n", m+1, count);exit(1);}

//fill Z (maybe permuted) - will fill Y (maybe permuted) below
if(permute==1){permute_int(order,count);}
for(i=0;i<count;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*count]=covar[order[i]+j*num_samples_use];}
}

if(mpheno!=-1)	//single phenotype
{
for(i=0;i<count;i++){Y[i]=resp[order[i]];}
strcpy(filename,outfile);
}
else	//analysing phenotypes in turn
{
printf("Analysing Response %d of %d\n", m+1, num_resps_use);
for(i=0;i<count;i++){Y[i]=resp[order[i]+m*num_samples_use];}
sprintf(filename, "%s.%d", outfile, m+1);
}

mle_fam(count, num_covars, Y, Z, relfile, prev, -9999, tol, maxiter, ids4, filename, constrain, (mode==130), cordups, num_blocks, 0);

for(i=0;i<count;i++){free(ids4[i]);}
}	//end of m loop

free(order);free(Y);free(Z);
free(ids4);
}

if(mode==229||mode==230)	//estimation of family coheritability
{
//allocate and fill response and covariates

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use*2);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);
ids4=malloc(sizeof(char*)*num_samples_use);

//index non-missing phenotypes (for both responses)
count=0;
for(i=0;i<num_samples_use;i++)
{
if(respinds[i]==1&&respinds2[i]==1){order[count]=i;copy_string(ids4,count,ids3[i]);count++;}
}
if(count==0){printf("Error, there are no samples with values for both phenotypes\n\n");exit(1);}
if(count<3){printf("Error, there are only %d samples with values for both phenotypes; unable to continue\n\n", count);exit(1);}

//fill Y and Z (maybe permuted)
if(permute==1){permute_int(order,count);}
for(i=0;i<count;i++)
{
Y[i]=resp[order[i]];
Y[i+count]=resp2[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*count]=covar[order[i]+j*num_samples_use];}
}

mle_fam(count, num_covars, Y, Z, relfile, prev, prev2, tol, maxiter, ids4, outfile, constrain, (mode==230), cordups, num_blocks, 1);

free(order);free(Y);free(Z);
for(i=0;i<count;i++){free(ids4[i]);}free(ids4);
}

///////////////////////////

if(mode==131)	//linear (num_kins will be 0 or 1)
{
if(strcmp(envfile,"blank")==0&&strcmp(locofile,"blank")==0&&families==0&&trios==0&&duos==0)	//basic case
{
#include "singlea.c"
}
if(strcmp(envfile,"blank")!=0)	//one environmental variable
{
#include "singleb.c"
}
if(strcmp(locofile,"blank")!=0)	//LOCO analysis
{
#include "singlec.c"
}
if(families==1)	//families analysis
{
#include "singlee.c"
}
if(trios==1)	//trio analysis
{
#include "singlef.c"
}
if(duos!=0)	//duos analysis
{
#include "singleg.c"
}
}

if(mode==132)	//logistic
{
if(strcmp(locofile,"blank")==0)	//basic case
{
#include "singled.c"
}
else	//LOCO analysis
{
#include "singlec.c"
}
}

////////

if(mode==133)	//solve-null - will always have shortcut=0 and num_resps_use=1
{
//anal_warn((2+memsave)*num_samples_use, num_samples_use);

sweights=malloc(sizeof(double)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

vstarts=malloc(sizeof(double)*num_kins);

//set sweights to one
for(i=0;i<num_samples_use;i++){sweights[i]=1;}

//fill Y and Z
for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[i];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

//estimate the variances
printf("Estimating the genetic variances\n");
multi_reml(num_samples_use, num_covars, 0, num_tops, num_kins, 0, Y, Z, Mkins, NULL, kintraces, kinsums, sweights, NULL, 0, NULL, NULL, NULL, prev, respcounts[0], hersfile, hestart, 0, NULL, NULL, 0, -9999, 0, oversfile, NULL, constrain, tol, maxiter, memsave, maxthreads, kinstems, ids3, filename, ids1, ids2, num_preds, allpreds, allal1, allal2, NULL, NULL, NULL, NULL, NULL, tkeeppreds, tcentres, vstarts);

//add up kinships and save (for simplicity, will read kinships from file)
printf("Creating the merged kinship matrix\n\n");
manip_kins(outfile, num_kins, kinstems, kinsums, ids1, ids2, ids3, num_samples_use, kingz, kinraw, 3, vstarts, maxthreads);
printf("Merged kinship matrix saved with stem %s; complete the mixed-model linear regression by using \"--linear\" with \"--grm %s\" (with the same covariates, top predictors and sample filterings used now)\n\n", outfile, outfile);

free(sweights);free(Y);free(Z);
free(vstarts);
}	//end of solve-null

////////

if(mode==136)	//cut-genes
{
//get an upper limit on number of genes/chunks
if(strcmp(genefile,"blank")!=0)	//using genefile - easy
{num_genes=countrows(genefile)-check_head(genefile,"Gene","Name",1);}

if(chunks!=-9999)	//using weights
{
sum=weights[0];count=1;
for(j=1;j<data_length;j++)
{
sum+=weights[j];
if(chr[j]!=chr[j-1]){count++;}
}
num_genes=sum/chunks*(1+overlap)+2*count;
}

if(chunksbp==1)	//each basepair observed becomes a gene - easy
{num_genes=data_length;}

if(chunksbp>1)
{
sum=0;count=1;
for(j=1;j<data_length;j++)
{
if(chr[j]!=chr[j-1]){sum+=bp[j-1];count++;}
}
sum+=bp[data_length-1];

if(sum/chunksbp*(1+overlap)+2*count>pow(2,31))
{printf("Error, there are too many chunks; to continue, you must increase chunksbp (currently %d)\n\n", chunksbp);exit(1);}
num_genes=sum/chunksbp*(1+overlap)+2*count;
}

gchr=malloc(sizeof(int)*num_genes);
gbp1=malloc(sizeof(double)*num_genes);
gbp2=malloc(sizeof(double)*num_genes);
gnames=malloc(sizeof(char*)*num_genes);
gstrand=malloc(sizeof(int)*num_genes);
gstarts=malloc(sizeof(int)*num_genes);
gends=malloc(sizeof(int)*num_genes);
for(g=0;g<num_genes;g++)
{
if(chunks!=-9999||chunksbp!=-9999){gnames[g]=malloc(sizeof(char)*100);}
gstarts[g]=-9999;gends[g]=-9999;
}

cut_genes(gchr, gbp1, gbp2, gnames, gstrand, gstarts, gends, data_length, chr, preds, bp, weights, genefile, chunks, chunksbp, up_buffer, down_buffer, minweight, overlap, pvafile, pvalues, part_length, bychr, 0, folder, datafile, bimfile, extract, -9999);

for(g=0;g<num_genes;g++){free(gnames[g]);}free(gnames);
free(gchr);free(gbp1);free(gbp2);free(gstrand);free(gstarts);free(gends);
}

////////

if(mode==137||mode==138)	//calc-genes-kins and calc-genes-reml
{
#include "genes.c"
}

if(mode==139)	//join-genes-reml
{
join_genes_reml(folder, num_genes, gnames, gstarts, gends, gparts, gchr, gbp1, gbp2, gamp, gam1, gam2, cut1, cut2, data_length, preds, keeppreds_use, 0, magma);
}

////////

if(mode==140)	//kvik-step3
{
//cut-genes

//get an upper limit on number of genes/chunks
if(strcmp(genefile,"blank")!=0)	//using genefile - easy
{num_genes=countrows(genefile)-check_head(genefile,"Gene","Name",1);}

if(chunks!=-9999)	//using weights
{
sum=weights[0];count=1;
for(j=1;j<data_length;j++)
{
sum+=weights[j];
if(chr[j]!=chr[j-1]){count++;}
}
num_genes=sum/chunks*(1+overlap)+2*count;
}

if(chunksbp==1)	//each basepair observed becomes a gene - easy
{num_genes=data_length;}

if(chunksbp>1)
{
sum=0;count=1;
for(j=1;j<data_length;j++)
{
if(chr[j]!=chr[j-1]){sum+=bp[j-1];count++;}
}
sum+=bp[data_length-1];

if(sum/chunksbp*(1+overlap)+2*count>pow(2,31))
{printf("Error, there are too many chunks; to continue, you must increase chunksbp (currently %d)\n\n", chunksbp);exit(1);}
num_genes=sum/chunksbp*(1+overlap)+2*count;
}

gchr=malloc(sizeof(int)*num_genes);
gbp1=malloc(sizeof(double)*num_genes);
gbp2=malloc(sizeof(double)*num_genes);
gnames=malloc(sizeof(char*)*num_genes);
gstrand=malloc(sizeof(int)*num_genes);
gstarts=malloc(sizeof(int)*num_genes);
gends=malloc(sizeof(int)*num_genes);
for(g=0;g<num_genes;g++)
{
if(chunks!=-9999||chunksbp!=-9999){gnames[g]=malloc(sizeof(char)*100);}
gstarts[g]=-9999;gends[g]=-9999;
}

cut_genes(gchr, gbp1, gbp2, gnames, gstrand, gstarts, gends, data_length, chr, preds, bp, weights, genefile, chunks, chunksbp, up_buffer, down_buffer, minweight, overlap, pvafile, pvalues, part_length, bychr, 3, outfile, datafile, bimfile, extract, -9999);

for(g=0;g<num_genes;g++){free(gnames[g]);}free(gnames);
free(gchr);free(gbp1);free(gbp2);free(gstrand);free(gstarts);free(gends);

////////

//calc-genes and join genes

//allocate and get gene details
#include "readdetails.c"

//perform analysis
#include "genes.c"

join_genes_reml(outfile, num_genes, gnames, gstarts, gends, gparts, gchr, gbp1, gbp2, gamp, gam1, gam2, cut1, cut2, data_length, preds, keeppreds_use, 1, magma);
}
///////////////////////////

if(mode==141)	//calc-tagging
{
#include "tagging.c"
}

if(mode==142)	//join-tagging
{
#include "join.c"
}

if(mode==143)	//merge-tagging
{
#include "combine.c"
}

if(mode==144)	//reduce-tagging
{
#include "reduce.c"
}

if(mode==145)	//calc-overlaps
{
#include "tagging.c"
}

////////

if(mode==146)	//sum-hers
{
#include "sumsa.c"

if(strcmp(matfile,"blank")!=0)	//also want per-predictor heritabilities
{
#include "exps.c"
}
}

if(mode==147)	//sum-cors
{
#include "sumsb.c"
}

if(mode==149)	//calc-exps
{
#include "exps.c"
}

if(mode==150)	//calc-posts
{
#include "posts.c"
}

///////////////////////////

if(mode==151||mode==152||mode==153||mode==154)	//blup, bolt, bayesr and elastic
{
#include "predmaster.c"
}

if(mode==156)	//calc-cors
{
#include "cors.c"
}

if(mode==157)	//join-cors
{
#include "concat.c"
}

if(mode==158)	//pseudo-summaries
{
#include "timesd.c"
}

if(mode==159)	//mega-prs
{
#include "litesa.c"
}

if(mode==160)	//validate
{
#include "timese.c"
}

///////////////////////////

if(mode==161)	//pca
{
if(strcmp(eigenfile,"blank")==0)
{
kin_warn(1, num_samples_use, 0, 0);
kins=malloc(sizeof(double*)*num_samples_use*num_samples_use);
}

U=malloc(sizeof(double)*num_samples_use*axes);
E=malloc(sizeof(double)*axes);

if(strcmp(eigenfile,"blank")==0)	//must compute pcs
{
read_kins(kinstems[0], kins, NULL, 1.0, num_samples_use, ids3, 0, maxthreads);

printf("Performing the decomposition; this can take a while\n\n");

vl=-1;vu=-1;	//vl and vu required but not used
count=num_samples_use-axes+1;
iwork=malloc(sizeof(int)*5*num_samples_use);
ifail=malloc(sizeof(int)*num_samples_use);
lwork=-1;
dsyevx_("V", "I", "U", &num_samples_use, kins, &num_samples_use, &vl, &vu, &count, &num_samples_use, &tol, &found, E, U, &num_samples_use, &wkopt, &lwork, iwork, ifail, &info);
if(info!=0){printf("PCA error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

decomp_warn(lwork);
work=malloc(sizeof(double)*lwork);
dsyevx_("V", "I", "U", &num_samples_use, kins, &num_samples_use, &vl, &vu, &count, &num_samples_use, &tol, &found, E, U, &num_samples_use, work, &lwork, iwork, ifail, &info);
if(info!=0){printf("PCA error 2; please tell Doug\n\n");exit(1);}

free(kins);
free(iwork);free(ifail);free(work);
}
else	//can get from eigen-decomposition
{read_eigens(eigenfile, U, E, num_samples_use, ids3, axes, maxthreads);}

sprintf(filename, "%s.vect", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
for(i=0;i<num_samples_use;i++)
{
fprintf(output, "%s %s ", ids1[i], ids2[i]);
for(k=0;k<axes;k++){fprintf(output, "%.6f ", U[(size_t)(axes-1-k)*num_samples_use+i]);}
fprintf(output, "\n");
}
fclose(output);

sprintf(filename2, "%s.values", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
for(k=0;k<axes;k++){fprintf(output2, "%.6f\n", E[axes-1-k]);}
fclose(output2);

sprintf(filename3, "%s.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
fprintf(output3, "Kinship %s\n", kinstems[0]);
fclose(output3);

printf("First %d eigenvectors saved in %s, with eigenvalues saved in %s\n\n", axes, filename, filename2);

free(U);free(E);
}	//end of mode=161

////////

if(mode==162)	//calc-pca-loads
{
#include "timesb.c"
}

if(mode==163)	//decompose and save
{
kin_warn(1, num_samples_use, 0, 0);

U=malloc(sizeof(double)*num_samples_use*num_samples_use);
E=malloc(sizeof(double)*num_samples_use);

read_kins(kinstems[0], U, NULL, 1.0, num_samples_use, ids3, 0, maxthreads);

printf("Performing the decomposition; this can take a while\n\n");

lwork=-1;
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

decomp_warn(lwork);
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_samples_use, U, &num_samples_use, E, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug\n\n");exit(1);}

write_eigen(outfile, U, E, num_samples_use, ids1, ids2, kinstems[0], eigenraw);

free(U);free(E);free(work);
}	//end of mode=163

if(mode==164)	//adjust-grm
{
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);

Z=malloc(sizeof(double)*num_samples_use*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);

Z_single=malloc(sizeof(float)*num_samples_use*num_fixed);
ZTZ_single=malloc(sizeof(float)*num_fixed*num_fixed);
kinZ_single=malloc(sizeof(float)*num_samples_use*num_fixed);
ZTkinZ_single=malloc(sizeof(float)*num_fixed*num_fixed);
W_single=malloc(sizeof(float)*num_samples_use*num_fixed);
WZTkinZ_single=malloc(sizeof(float)*num_samples_use*num_fixed);

read_kins(kinstems[0], NULL, kins_single, 1.0, num_samples_use, ids3, 0, maxthreads);

//fill Z
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

//get inverseZTZ
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &num_samples_use, &alpha, Z, &num_samples_use, Z, &num_samples_use, &beta, ZTZ, &num_fixed);
(void)eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

//now copy Z and ZTZ into single versions
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z_single[i+j*num_samples_use]=Z[i+j*num_samples_use];}
}
for(j=0;j<num_fixed;j++)
{
for(k=0;k<num_fixed;k++){ZTZ_single[j+k*num_fixed]=ZTZ[j+k*num_fixed];}
}

//then kinZ, ZTkinZ, W=ZinvZTZ and WZTkinZ
alpha_single=1.0;beta_single=0.0;
sgemm_("N", "N", &num_samples_use, &num_fixed, &num_samples_use, &alpha_single, kins_single, &num_samples_use, Z_single, &num_samples_use, &beta_single, kinZ_single, &num_samples_use);
sgemm_("T", "N", &num_fixed, &num_fixed, &num_samples_use, &alpha_single, Z_single, &num_samples_use, kinZ_single, &num_samples_use, &beta_single, ZTkinZ_single, &num_fixed);
sgemm_("N", "N", &num_samples_use, &num_fixed, &num_fixed, &alpha_single, Z_single, &num_samples_use, ZTZ_single, &num_fixed, &beta_single, W_single, &num_samples_use);
sgemm_("N", "N", &num_samples_use, &num_fixed, &num_fixed, &alpha_single, W_single, &num_samples_use, ZTkinZ_single, &num_fixed, &beta_single, WZTkinZ_single, &num_samples_use);

//ready to adjust kins
printf("Regressing the kinship matrix on the fixed effects\n\n");

//first subtract WZTkin and kinZWT
alpha_single=-1.0;beta_single=1.0;
sgemm_("N", "T", &num_samples_use, &num_samples_use, &num_fixed, &alpha_single, W_single, &num_samples_use, kinZ_single, &num_samples_use, &beta_single, kins_single, &num_samples_use);
sgemm_("N", "T", &num_samples_use, &num_samples_use, &num_fixed, &alpha_single, kinZ_single, &num_samples_use, W_single, &num_samples_use, &beta_single, kins_single, &num_samples_use);

//then add WZTkinZ WT
alpha_single=1.0;beta_single=1.0;
sgemm_("N", "T", &num_samples_use, &num_samples_use, &num_fixed, &alpha_single, WZTkinZ_single, &num_samples_use, W_single, &num_samples_use, &beta_single, kins_single, &num_samples_use);

write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 1);

if(kindetails==1)	//copy over details and adjust
{
sprintf(cmd, "cp %s.grm.details %s.grm.details", kinstems[0], outfile);
system(cmd);
sprintf(cmd, "cp %s.grm.adjust %s.grm.adjust", kinstems[0], outfile);
system(cmd);
}

//make root
sprintf(filename,"%s.grm.root",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
fprintf(output,"Kinship %s\n", kinstems[0]);
if(strcmp(covarfile,"blank")!=0){fprintf(output,"Covariates %s\n", covarfile);}
else{fprintf(output,"Covariates none\n");}
if(strcmp(topfile,"blank")!=0){fprintf(output,"Top_Predictors %s\n", topfile);}
else{fprintf(output,"Top_Predictors none\n");}
if(strcmp(envfile,"blank")!=0){fprintf(output,"Environments %s\n", envfile);}
else{fprintf(output,"Environments none\n");}
fclose(output);

free(kins_single);
free(Z);free(ZTZ);free(ZTZ2);free(ZTZ3);
free(Z_single);free(ZTZ_single);free(kinZ_single);free(ZTkinZ_single);free(W_single);free(WZTkinZ_single);
}	//end of mode=164

////////

if(mode==166||mode==167||mode==168||mode==169||mode==170)	//truncate-grm, pca-grm, square-grm and gxemms
{
#include "others.c"
}

///////////////////////////

if(mode==171)	//calc-stats
{
#include "stats.c"
}

if(mode==172)	//calc-blups
{
#include "timesc.c"
}

if(mode==173)	//make-phenos
{
#include "phens.c"
}

if(mode==174)	//make-snps
{
if(trios==1){gens=malloc(sizeof(int)*num_inds);}
if(quads==1){gens=malloc(sizeof(int)*8);}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%s.bed", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
gen=108;fwrite(&gen, sizeof(unsigned char), 1, output2);
gen=27;fwrite(&gen, sizeof(unsigned char), 1, output2);
gen=1;fwrite(&gen, sizeof(unsigned char), 1, output2);

//how many ind in each population?
count=(num_inds-1)/pops+1;
if(count<num_inds){printf("There will be %d individuals in each population\n\n", count);}

bittotal=(num_snps-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>num_snps){bitend=num_snps;}
bitlength=bitend-bitstart;

printf("Making SNPs for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Making SNPs for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

if(trios==0&&quads==0)	//regular case
{
for(j=0;j<bitlength;j++)
{
gen=0;
for(i=0;i<num_inds;i++)
{
if(i%count==0)	//start a new population
{
unifrand=genrand_real1();
maf=unifrand*(maf2-maf1)+maf1;
}

//flag indicates if we need to make a new genotype
if(i%dups==0)	//first individual in a family - always make a new genotype
{flag=1;}
else	//subsequent individual in a family
{
if(i%dups==1)	//decide whether it matches previous (if i%dups>1, will have set flag already)
{
unifrand=genrand_real1();
if(unifrand<closeness){flag=0;}
else{flag=1;}
}
}

if(flag==1)	//generate a new genotype
{
unifrand=genrand_real1();

gen2=0;					//genotype 2
if(unifrand<1-pow(maf,2)){gen2=2;}	//genotype 1
if(unifrand<pow(1-maf,2)){gen2=3;}	//genotype 0
}
gen+=(gen2<<(2*(i%4)));
if(i%4==3||i==num_inds-1){fwrite(&gen, sizeof(unsigned char), 1, output2);gen=0;}
}	//end of i loop
}
}	//end of regular case

if(trios==1)
{
for(j=0;j<bitlength;j++)
{
unifrand=genrand_real1();
maf=unifrand*(maf2-maf1)+maf1;

//first make all values for snp
for(i=0;i<num_inds;i+=3)
{
//make first allele of father
unifrand=genrand_real1();
gen=(unifrand<maf);

//now finish and save genotype
unifrand=genrand_real1();
gens[i]=gen+(unifrand<maf);

//make first allele of mother
unifrand=genrand_real1();
gen2=(unifrand<maf);

//now finish and save genotype
unifrand=genrand_real1();
gens[i+1]=gen2+(unifrand<maf);

//now save child
gens[i+2]=gen+gen2;
}

//can now save all values for snp
gen=0;
for(i=0;i<num_inds;i++)
{
if(gens[i]==0){gen2=3;}
if(gens[i]==1){gen2=2;}
if(gens[i]==2){gen2=0;}

gen+=(gen2<<(2*(i%4)));
if(i%4==3||i==num_inds-1){fwrite(&gen, sizeof(unsigned char), 1, output2);gen=0;}
}
}
}	//end of trios=1

if(quads==1)
{
for(j=0;j<bitlength;j++)
{
unifrand=genrand_real1();
maf=unifrand*(maf2-maf1)+maf1;

for(i=0;i<num_inds;i+=4)
{
//make haplotypes for father
unifrand=genrand_real1();
gens[0]=(unifrand<maf);
unifrand=genrand_real1();
gens[1]=(unifrand<maf);

//make haplotypes for mother
unifrand=genrand_real1();
gens[2]=(unifrand<maf);
unifrand=genrand_real1();
gens[3]=(unifrand<maf);

//make haplotypes for son
unifrand=genrand_real1();
gens[4]=gens[(unifrand<maf)];
unifrand=genrand_real1();
gens[5]=gens[2+(unifrand<maf)];

//make haplotypes for daughter
unifrand=genrand_real1();
gens[6]=gens[(unifrand<maf)];
unifrand=genrand_real1();
gens[7]=gens[2+(unifrand<maf)];

//now save
gen=0;

//father
if(gens[0]+gens[1]==0){gen2=3;}
if(gens[0]+gens[1]==1){gen2=2;}
if(gens[0]+gens[1]==2){gen2=0;}
gen+=gen2;

//mother
if(gens[2]+gens[3]==0){gen2=3;}
if(gens[2]+gens[3]==1){gen2=2;}
if(gens[2]+gens[3]==2){gen2=0;}
gen+=(gen2<<2);

//son
if(gens[4]+gens[5]==0){gen2=3;}
if(gens[4]+gens[5]==1){gen2=2;}
if(gens[4]+gens[5]==2){gen2=0;}
gen+=(gen2<<4);

//daughter
if(gens[6]+gens[7]==0){gen2=3;}
if(gens[6]+gens[7]==1){gen2=2;}
if(gens[6]+gens[7]==2){gen2=0;}
gen+=(gen2<<6);

fwrite(&gen, sizeof(unsigned char), 1, output2);
}
}
}	//end of quads=1
}	//end of bit loop
printf("\n");
fclose(output2);

sprintf(filename3, "%s.bim", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
mark=0;mark2=1;
for(j=0;j<num_snps;j++)
{
if(1+j*nchrom/num_snps>mark2)	//new chromosome
{mark=j;mark2++;}
fprintf(output3,"%d SNP%d %.2f %.0f A C\n", mark2, j+1, 0.01*(j-mark+1), 10000.0*(j-mark+1));
}
fclose(output3);

sprintf(filename4, "%s.fam", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename4);exit(1);}
if(trios==0&&quads==0)
{
for(i=0;i<num_inds;i++){fprintf(output4,"FAMILY%d IND%d 0 0 %d 0\n", i/dups+1, i%dups+1, i%2+1);}
}
if(trios==1)
{
for(i=0;i<num_inds;i+=3)
{
fprintf(output4,"FAMILY%d IND1 0 0 2 0\n", i/3+1);
fprintf(output4,"FAMILY%d IND2 0 0 1 0\n", i/3+1);
fprintf(output4,"FAMILY%d IND3 IND1 IND2 %d 0\n", i/3+1, i%2+1);
}
}
if(quads==1)
{
for(i=0;i<num_inds;i+=4)
{
fprintf(output4,"FAMILY%d IND1 0 0 2 0\n", i/4+1);
fprintf(output4,"FAMILY%d IND2 0 0 1 0\n", i/4+1);
fprintf(output4,"FAMILY%d IND3 IND1 IND2 2 0\n", i/4+1);
fprintf(output4,"FAMILY%d IND4 IND1 IND2 1 0\n", i/4+1);
}
}
fclose(output4);

sprintf(filename5, "%s.covar", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename5);exit(1);}
fprintf(output5,"FID IID Sex Age\n");
if(trios==0&&quads==0)
{
for(i=0;i<num_inds;i++){fprintf(output5,"FAMILY%d IND%d %d %d\n", i/dups+1, i%dups+1, i%2+1, 20+i%10);}
}
if(trios==1)
{
for(i=0;i<num_inds;i+=3)
{
fprintf(output5,"FAMILY%d IND1 2 %d\n", i/3+1, 40+i%10);
fprintf(output5,"FAMILY%d IND2 1 %d\n", i/3+1, 40+i%10);
fprintf(output5,"FAMILY%d IND3 %d %d\n", i/3+1, i%2+1, 20+i%10);
}
}
if(quads==1)
{
for(i=0;i<num_inds;i+=4)
{
fprintf(output5,"FAMILY%d IND1 2 %d\n", i/4+1, 40+i%10);
fprintf(output5,"FAMILY%d IND2 1 %d\n", i/4+1, 40+i%10);
fprintf(output5,"FAMILY%d IND3 2 %d\n", i/4+1, 20+i%10);
fprintf(output5,"FAMILY%d IND4 1 %d\n", i/4+1, 20+i%10);
}
}
fclose(output5);

sprintf(filename6,"%s.pairs", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}

if(trios==1)
{
for(i=0;i<num_inds;i+=3)
{
fprintf(output6,"FAMILY%d IND1 FAMILY%d IND3 0.5\n", i/3+1, i/3+1);
fprintf(output6,"FAMILY%d IND2 FAMILY%d IND3 0.5\n", i/3+1, i/3+1);
}
}
if(quads==1)
{
for(i=0;i<num_inds;i+=4)
{
fprintf(output6,"FAMILY%d IND1 FAMILY%d IND3 0.5\n", i/3+1, i/3+1);
fprintf(output6,"FAMILY%d IND2 FAMILY%d IND3 0.5\n", i/3+1, i/3+1);
fprintf(output6,"FAMILY%d IND1 FAMILY%d IND4 0.5\n", i/3+1, i/3+1);
fprintf(output6,"FAMILY%d IND2 FAMILY%d IND4 0.5\n", i/3+1, i/3+1);
fprintf(output6,"FAMILY%d IND3 FAMILY%d IND4 0.5\n", i/3+1, i/3+1);
}
}
if(dups>1)
{
for(i=0;i<num_inds;i+=dups)
{
for(i2=0;i2<dups-1;i2++)
{
for(i3=i2+1;i3<dups;i3++)
{
if(i+i2<num_inds&&i+i3<num_inds){fprintf(output6,"FAMILY%d IND%d FAMILY%d IND%d %.6f\n", i/dups+1, i2+1, i/dups+1, i3+1, closeness);}
}
}
}
}
if(trios==0&&quads==0&&dups==1){fprintf(output6,"No Related Pairs\n");}
fclose(output6);

printf("Simulated data saved in %s, %s and %s, with pretend covariates in %s", filename2, filename3, filename4, filename5);
if(trios==1||quads==1||dups>1){printf(", and a list of related pairs in %s", filename6);}
printf("\n\n");

if(trios==1||quads==1){free(gens);}
}	//end of mode=174

if(mode==175)	//calc-inflation
{
#include "inflate.c"
}

////////

if(mode==176)	//jackknife
{
#include "jackknife.c"
}

if(mode==177)	//cut-folds
{
usedids=malloc(sizeof(int)*num_samples_use);
order=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){order[i]=i;}
permute_int(order,num_samples_use);

count=(num_samples_use-1)/num_folds+1;
printf("Dividing the %d samples into %d folds with size (approximately) %d\n\n", num_samples_use, num_folds, count);

for(k=0;k<num_folds;k++)
{
for(i=0;i<num_samples_use;i++){usedids[i]=0;}
for(i=0;i<num_samples_use;i++)
{
if(i>=k*count&&i<(k+1)*count){usedids[order[i]]=1;}
}

sprintf(filename,"%s.train%d", outfile, k+1);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
sprintf(filename2,"%s.test%d", outfile, k+1);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==1){fprintf(output2,"%s %s\n", ids1[i], ids2[i]);}
else{fprintf(output,"%s %s\n", ids1[i], ids2[i]);}
}
fclose(output);
fclose(output2);
}	//end of k loop

printf("Training folds saved in %s.train1, ..., %s.train%d;  testing folds in %s.test1, .., %s.test%d\n\n", outfile, outfile, num_folds, outfile, outfile, num_folds);

free(usedids);free(order);
}	//end of mode=177

if(mode==178)	//find-gaussian
{
#include "grid.c"
}

if(mode==179)	//winners-curse
{
#include "winners.c"
}

///////////////////////////

if(mode==181||mode==182||mode==183||mode==184||mode==185)	//make-data
{
#include "merge.c"
}

if(mode==186||mode==187||mode==188||mode==189)	//condense data
{
#include "condense.c"
}	//end of condense data

if(mode==190)	//calc-sim-data
{
#include "merge.c"
}

///////////////////////////

if(mode==191)	//cut-gre
{
#include "gre.c"
}

if(mode==192)	//cut-gre
{
#include "gre.c"
}

if(mode==193)	//join-gre
{
#include "gre.c"
}

if(mode==194)	//solve-gre
{
#include "gre.c"
}

///////////////////////////

if(mode==201)	//speed-tests
{
count=4;
if(num_samples_use>=20000){count=2;}
if(num_samples_use>=40000){count=1;}

//read in kinship
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kins(kinstems[0], NULL, kins_single, 1.0, num_samples_use, ids3, 0, maxthreads);

//save in kins
kin_warn(1, num_samples_use, 0, 0);
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){kins[(size_t)i2*num_samples_use+i]=kins_single[(size_t)i2*num_samples_use+i];}
}

Rvsing=malloc(sizeof(float)*num_samples_use*num_vects);
E=malloc(sizeof(double)*num_samples_use);
ipiv=malloc(sizeof(int)*num_samples_use);

//first, want to do single cholesky and solve
printf("First %d cholesky solves for %d individuals\n\n", count, num_samples_use);

sprintf(cmd,"date > %s.time.schol.start", outfile);
system(cmd);

for(p=0;p<count;p++)
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){kins_single[(size_t)i2*num_samples_use+i]=kins[(size_t)i2*num_samples_use+i];}
}

for(g=0;g<num_vects;g++)
{
for(i=0;i<num_samples_use;i++){Rvsing[i+g*num_samples_use]=rnorm_safe();}
}

//do cholesky
spotrf_("U", &num_samples_use, kins_single, &num_samples_use, &info);
if(info!=0){printf("Error, Cholesky failed (info %d, length %d); please tell Doug\n", info, num_samples_use);exit(1);}

//solve
spotrs_("U", &num_samples_use, &num_vects, kins_single, &num_samples_use, Rvsing, &num_samples_use, &info);
if(info!=0){printf("Error first Cholesky solve failed, please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}
}

sprintf(cmd,"date > %s.time.schol.end", outfile);
system(cmd);

for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++)
{kins_single[(size_t)i2*num_samples_use+i]=kins[(size_t)i2*num_samples_use+i];}
}

//now do count double LDLT and inverse
printf("Now %d LDLT inverses for %d individuals\n\n", count, num_samples_use);

sprintf(cmd,"date > %s.time.dldlt.start", outfile);
system(cmd);

for(p=0;p<count;p++)
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++)
{kins[(size_t)i2*num_samples_use+i]=kins_single[(size_t)i2*num_samples_use+i];}
}

//do ldlt
lwork=-1;
dsytrf_("U", &num_samples_use, kins, &num_samples_use, ipiv, &wkopt, &lwork, &info);
if(info!=0)
{printf("Error, LDLT priming failed; please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsytrf_("U", &num_samples_use, kins, &num_samples_use, ipiv, work, &lwork, &info);
if(info!=0)
{printf("Error, LDLT decomposition failed, please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}

dsytri_("U", &num_samples_use, kins, &num_samples_use, ipiv, work, &info);
if(info!=0)
{printf("Error, LDLT inverting failed, please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}
free(work);
}

sprintf(cmd,"date > %s.time.dldlt.end", outfile);
system(cmd);

//now do count eigens and inverse
printf("Now %d eigendecomps for %d individuals\n\n", count, num_samples_use);

sprintf(cmd,"date > %s.time.eigen.start", outfile);
system(cmd);

for(p=0;p<count;p++)
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++)
{kins[(size_t)i2*num_samples_use+i]=kins_single[(size_t)i2*num_samples_use+i];}
}

lwork=-1;
dsyev_("V", "U", &num_samples_use, kins, &num_samples_use, E, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_samples_use, kins, &num_samples_use, E, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug info %d\n\n", info);exit(1);}
free(work);
}

sprintf(cmd,"date > %s.time.eigen.end", outfile);
system(cmd);

free(kins_single);
free(kins);
free(Rvsing);
free(ipiv);
free(E);
}	//end of speed-tests

////////

if(mode==202)	//speed-tests
{
//read in kinship
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kins(kinstems[0], NULL, kins_single, 1.0, num_samples_use, ids3, 0, maxthreads);

//save in kins
kin_warn(1, num_samples_use, 0, 0);
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);

E=malloc(sizeof(double)*num_samples_use);

//allocate and make random vectors
Rv=malloc(sizeof(double)*num_samples_use*num_vects);
Rv2=malloc(sizeof(double)*num_samples_use*num_vects);
for(g=0;g<num_vects;g++)
{
for(i=0;i<num_samples_use;i++){Rv[i+g*num_samples_use]=rnorm_safe();}
}

//first, want to do 10 cholesky and solve
printf("First 10 cholesky solves for %d individuals\n\n", num_samples_use);

sprintf(cmd,"date > %s.time.csolve.start", outfile);
system(cmd);

for(p=0;p<10;p++)
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){kins[(size_t)i2*num_samples_use+i]=kins_single[(size_t)i2*num_samples_use+i];}
}

//do cholesky
dpotrf_("U", &num_samples_use, kins, &num_samples_use, &info);
if(info!=0){printf("Error, Cholesky failed (info %d, length %d); please tell Doug\n", info, num_samples_use);exit(1);}

//solve
dpotrs_("U", &num_samples_use, &num_vects, kins, &num_samples_use, Rv2, &num_samples_use, &info);
if(info!=0){printf("Error first Cholesky solve failed, please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}
}

sprintf(cmd,"date > %s.time.csolve.end", outfile);
system(cmd);

//now do 10 cholesky and inverse
printf("Now 10 cholesky inverses for %d individuals\n\n", num_samples_use);

sprintf(cmd,"date > %s.time.cinv.start", outfile);
system(cmd);

for(p=0;p<10;p++)
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){kins[(size_t)i2*num_samples_use+i]=kins_single[(size_t)i2*num_samples_use+i];}
}

//do cholesky
dpotrf_("U", &num_samples_use, kins, &num_samples_use, &info);
if(info!=0){printf("Error, Cholesky failed (info %d, length %d); please tell Doug\n", info, num_samples_use);exit(1);}

//solve
dpotri_("U", &num_samples_use, kins, &num_samples_use, &info);
if(info!=0){printf("Error second Cholesky solve failed, please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}
}

sprintf(cmd,"date > %s.time.cinv.end", outfile);
system(cmd);

//now do 2 eigens and inverse
printf("Now 2 eigendecomps for %d individuals\n\n", num_samples_use);

sprintf(cmd,"date > %s.time.eigen.start", outfile);
system(cmd);

for(p=0;p<2;p++)
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){kins[(size_t)i2*num_samples_use+i]=kins_single[(size_t)i2*num_samples_use+i];}
}

lwork=-1;
dsyev_("V", "U", &num_samples_use, kins, &num_samples_use, E, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_samples_use, kins, &num_samples_use, E, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug info %d\n\n", info);exit(1);}
free(work);
}

sprintf(cmd,"date > %s.time.eigen.end", outfile);
system(cmd);

//finally, want to do 10 cholesky and solve
printf("Final 10 cholesky solves for %d individuals\n\n", num_samples_use);

sprintf(cmd,"date > %s.time.ccheck.start", outfile);
system(cmd);

for(p=0;p<10;p++)
{
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){kins[(size_t)i2*num_samples_use+i]=kins_single[(size_t)i2*num_samples_use+i];}
}

//do cholesky
dpotrf_("U", &num_samples_use, kins, &num_samples_use, &info);
if(info!=0){printf("Error, Cholesky failed (info %d, length %d); please tell Doug\n", info, num_samples_use);exit(1);}

//solve
dpotrs_("U", &num_samples_use, &num_vects, kins, &num_samples_use, Rv2, &num_samples_use, &info);
if(info!=0){printf("Error first Cholesky solve failed, please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}
}

sprintf(cmd,"date > %s.time.ccheck.end", outfile);
system(cmd);

free(kins_single);
free(kins);
}	//end of speed-tests2

////////

if(mode==203)	//speed-tests3
{
//read in kinship
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kins(kinstems[0], NULL, kins_single, 1.0, num_samples_use, ids3, 0, maxthreads);

//save in binary file
sprintf(filename,"%s.kin.temp", outfile);
if((output=fopen(filename,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(i=0;i<num_samples_use;i++){fwrite(kins_single+(size_t)i*num_samples_use, sizeof(float), i+1, output);}
fclose(output);

free(kins_single);

scount=(size_t)num_samples_use*(num_samples_use+1)/2;
kins_packed=malloc(sizeof(float)*scount);

if(num_vects!=-9999)	//do linear system
{
//allocate and make random vectors
Rvsing=malloc(sizeof(float)*num_samples_use*num_vects);
Rvsing2=malloc(sizeof(float)*num_samples_use*num_vects);
Rvsing3=malloc(sizeof(float)*num_samples_use*num_vects);
for(g=0;g<num_vects;g++)
{
for(i=0;i<num_samples_use;i++){Rvsing[i+g*num_samples_use]=rnorm_safe();}
}

//write a start time file
sprintf(cmd,"date > %s.dime.solve.start", outfile);
system(cmd);

for(p=0;p<10;p++)
{
printf("Loop %d solve\n", p+1);

//read back in kinship
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
fseeko(input, 0, SEEK_SET);

scount2=0;
for(i=0;i<num_samples_use;i++)
{
if(fread(kins_packed+scount2, sizeof(float), i+1, input)!=i+1)
{printf("Error reading row %d out of %d\n\n", i+1, num_samples_use);exit(1);}
scount2+=i+1;
}
fclose(input);

if(p==1)	//confirm last inverse worked
{
alpha_single=1.0;beta_single=0.0;
sspmv_("U", &num_samples_use, &alpha_single, kins_packed, Rvsing2, &one, &beta_single, Rvsing3, &one);

for(i=0;i<5;i++){printf("%d %f %f\n", i+1, Rvsing[i], Rvsing3[i]);}
for(i=num_samples_use-5;i<num_samples_use;i++){printf("%d %f %f\n", i+1, Rvsing[i], Rvsing3[i]);}
}

//read in randoms
for(g=0;g<num_vects;g++)
{
for(i=0;i<num_samples_use;i++){Rvsing2[i+g*num_samples_use]=Rvsing[i+g*num_samples_use];}
}

//do cholesky
spptrf_("U", &num_samples_use, kins_packed, &info);
if(info!=0){printf("Error, Cholesky failed (info %d, length %d); please tell Doug\n", info, num_samples_use);exit(1);}

//solve
spptrs_("U", &num_samples_use, &num_vects, kins_packed, Rvsing2, &num_samples_use, &info);
if(info!=0){printf("Error first Cholesky solve failed, please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}
}

//write an end time file
sprintf(cmd,"date > %s.dime.solve.end", outfile);
system(cmd);

free(Rvsing);free(Rvsing2);free(Rvsing3);
}	//end of solving

sprintf(cmd,"rm %s", filename);
system(cmd);

free(kins_single);
}	//end of speed-tests3

///////////////////////////
//free stuff
///////////////////////////

//allocations from readdetails.c

if(mode==102||mode==103||mode==104){free(sstarts);free(sends);free(sstarts2);free(sends2);}
if(mode==112){free(pstarts);free(pends);}
if(mode==137||mode==138||mode==139||mode==140)
{
for(g=0;g<num_genes;g++){free(gnames[g]);}free(gnames);
free(gstarts);free(gends);free(gparts);free(gchr);free(gbp1);free(gbp2);free(gpvas);
}
if(mode==192||mode==193||mode==194){free(pstarts);free(pends);}

//datafile, kinships and resp allocations (from parsefiles.c and top of ldak.c)

if(use_data==1||use_data==2||use_data==3||use_data==4||use_data==5)
{
for(k=0;k<num_files;k++){free(datastems[k]);free(bimstems[k]);free(famstems[k]);}
free(datastems);free(bimstems);free(famstems);
}

if(num_kins>0)
{
for(k=0;k<num_kins;k++){free(kinstems[k]);}free(kinstems);
free(kinnums);free(kinsums);

if(mode==120||mode==121||mode==123||mode==124||mode==126||mode==131||mode==133||mode==138)
{
if(mode==120||mode==123||mode==124||mode==126){if(memsave==0){for(k=0;k<num_kins;k++){free(Mkins_single[k]);}}free(Mkins_single);}
else{if(memsave==0){for(k=0;k<num_kins;k++){free(Mkins[k]);}}free(Mkins);}
free(kintraces);
}
}

if(num_resps_use>0){free(keepresps);free(resp);free(respcounts);free(respinds);}
if(mode==229||mode==230){free(keepresps2);free(resp2);free(respcounts2);free(respinds2);}

//more allocations from parsefiles.c

if(strcmp(oversfile,"blank")!=0)
{
for(k=0;k<num_kins;k++){free(ssums[k]);}
free(ssums);
}

if(strcmp(taglist,"blank")!=0)
{
for(k=0;k<num_tags;k++){free(tagstems[k]);}
free(tagstems);
}

if(strcmp(matlist,"blank")!=0)
{
for(k=0;k<num_tags;k++){free(matstems[k]);}
free(matstems);
}

if(strcmp(catfile,"blank")!=0){free(keepparts);free(keepparts2);}

if(strcmp(corslist,"blank")!=0)
{
for(k=0;k<num_cors;k++){free(corstems[k]);}
free(corstems);
}

//allocations from getnums.c 

if(strcmp(idsfile,"blank")!=0)
{
for(i=0;i<num_samples_use;i++){free(ids1[i]);free(ids2[i]);free(ids3[i]);}
free(ids1);free(ids2);free(ids3);

if(num_subs>0)
{
for(s=0;s<num_subs;s++){free(subindex[s]);}
free(subindex);
}

if(use_data==1||use_data==2||use_data==4){free(keepsamps);}
}

if(use_data==1||use_data==2||use_data==3||use_data==6)
{
free(allchr);free(allcm);free(allbp);free(allal1);free(allal2);
for(j=0;j<num_preds;j++){free(allpreds[j]);}free(allpreds);
free(predorder);free(keeppreds);
}

//global allocations from setdl.c - will have done mode-specific ones above

if(use_data==1||use_data==6)
{
free(keeppreds_use);
free(chr);free(cm);free(bp);free(cmbp);free(al1);free(al2);
for(j=0;j<data_length;j++){free(preds[j]);}free(preds);
}

//tops and regions (from top of ldak.c)

if(strcmp(topfile,"blank")!=0)
{
free(tkeeppreds);
free(tchr);free(tbp);free(tal1);free(tal2);
for(j=0;j<num_tops;j++){free(tpreds[j]);}free(tpreds);
free(tcentres);free(tvars);
if(strcmp(sumsfile,"blank")!=0){free(tnss);free(tchis);free(trhos);}
}

if(num_regs>0)
{
free(rkeeppreds);
for(r=0;r<num_regs;r++){free(regindex[r]);}free(regindex);
free(ral1);free(ral2);
for(j=0;j<rnum_preds_use;j++){free(rpreds[j]);}free(rpreds);
free(rcentres);free(rmults);free(rsqdevs);free(rweights);
if(strcmp(sumsfile,"blank")!=0){free(rnss);free(rchis);free(rrhos);}
free(rdata);
}

//more allocations from top of ldak.c

if(mode==131&&(families==1||trios==1||duos!=0)){free(famindex);free(famcounts);}

if(use_data==1||use_data==6)
{
free(centres);free(mults);free(sqdevs);free(weights);free(pvalues);
if(strcmp(sumsfile,"blank")!=0){free(nss);free(chis);free(rhos);}
}

if(num_fixed>0)
{
free(covar);
if(strcmp(cofile,"blank")!=0){free(thetas);}
}

if(mode==132){free(offsets);}

if(num_kins==1&&((mode==121&&shortcut==1)||mode==131||mode==138)){free(U);free(E);}

///////////////////////////

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\nMission completed. All your basepair are belong to us :)\n-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n");

return(0);
}	//end of main

///////////////////////////

