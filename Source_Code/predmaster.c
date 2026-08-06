/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Wrapper for making prediction models

///////////////////////////

//deal with progress file
if(substep<0){sprintf(filename,"%s.progress",outfile);}
else{sprintf(filename,"%s.pheno%d.progress",outfile, substep);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);} 
fclose(output);

//will do all analyses in bits
bittotal=(data_length-1)/bitsize+1;

if(substep==-1)	//save top of multivariate.bin
{
sprintf(filename2,"%s.multivariate.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fwrite(&num_samples_use, sizeof(int), 1, output2);
fwrite(&data_length, sizeof(int), 1, output2);
fwrite(&num_resps_use, sizeof(int), 1, output2);
fwrite(&bittotal, sizeof(int), 1, output2);

fclose(output2);
}

if(substep>0)	//check size and top of cutting.bin
{
sprintf(filename2,"%s.multivariate.bin", outfile);
if(just_check(filename2)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step CUT\"\n\n", filename2);exit(1);}

if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}

//read first four elements
fseeko(input2, 0, SEEK_END);
if(ftello(input2)<(off_t)sizeof(int)*4)
{printf("Error readint %s; should have size at least %jd, but instead has size %jd\n\n", filename2, (off_t)sizeof(int)*4, ftello(input2));exit(1);}

fseeko(input2, 0, SEEK_SET);

if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading number of samples from %s\n\n", filename2);exit(1);}
if(readint!=num_samples_use)
{printf("Error, number of samples in %s (%d) does not match number being used now (%d)\n\n", filename2, readint, num_samples_use);exit(1);}

if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading number of predictors from %s\n\n", filename2);exit(1);}
if(readint!=data_length)
{printf("Error, number of predictors in %s (%d) does not match number being used now (%d)\n\n", filename2, readint, data_length);exit(1);}

if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading number of responses from %s\n\n", filename2);exit(1);}
if(readint!=num_resps_use)
{printf("Error, number of chunk in %s (%d) does not match number being used now (%d)\n\n", filename2, readint, num_resps_use);exit(1);}

if(fread(&readint, sizeof(int), 1, input2)!=1)
{printf("Error reading number of chunks from %s\n\n", filename2);exit(1);}
if(readint!=bittotal)
{printf("Error, number of chunk in %s (%d) does not match number being used now (%d)\n\n", filename2, readint, bittotal);exit(1);}

//check full size
fseeko(input2, 0, SEEK_END);
if(ftello(input2)!=(off_t)sizeof(int)*(6+num_samples_use+num_resps_use)+sizeof(double)*data_length*(6+num_resps_use+bitsize)+sizeof(double)*num_resps_use*(2+2*num_resps_use+2*num_samples_use))
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename2, (off_t)sizeof(int)*(6+num_samples_use+num_resps_use)+sizeof(double)*data_length*(6+num_resps_use+bitsize)+sizeof(double)*num_resps_use*(2+2*num_resps_use+2*num_samples_use), ftello(input2));exit(1);}

fclose(input2);
}

if(skipcv==0)	//get num_train, num_test and cv samples (substep=-1 adds these to multivariate.bin; substep>0 reads these from multivariate.bin)
{
keeptrain=malloc(sizeof(int)*num_samples_use);
keeptest=malloc(sizeof(int)*num_samples_use);

#include "setcv.c"
}

//allocate powers and load up possible values
if(power!=-9999){num_pows=1;}
else
{
if(strcmp(powfile,"blank")==0){num_pows=5;}
else{num_pows=countrows(powfile);}
}

powers=malloc(sizeof(double)*num_pows);

if(power!=-9999){powers[0]=power;}
else
{
if(strcmp(powfile,"blank")==0)
{powers[0]=-1.0;powers[1]=-0.75;powers[2]=-0.5;powers[3]=-0.25;powers[4]=0;}
else
{
printf("Reading %d power values from %s\n\n", num_pows, powfile);
read_values(powfile,powers,num_pows,NULL,1,0,0);
}
}

if(checkped==1)	//check number of pedigree predictors
{
if(nped>data_length)
{
nped=data_length;
printf("Warning, the number of pedigree predictors has been reduced to %d (the total number of predictors)\n\n", nped);
}

if(nped>1024){printf("Warning, there are many pedigree predictors (%d); this can substantially increase memory usage and runtmime\n\n", nped);}
}

if(ncal>0)	//check number of calibration predictors
{
if(ncal>data_length)
{
ncal=data_length;
printf("Warning, the number of calibration predictors has been reduced to %d (the total number of predictors)\n\n", data_length);
}

if(ncal>32){printf("Warning, there are many calibration predictors (%d); this can substantially increase memory usage and runtime\n\n", ncal);}
}

if(ncomp>0)	//check number of comparison predictors
{
if(ncomp>data_length)
{
ncomp=data_length;
printf("Warning, the number of comparison predictors has been reduced to %d (the total number of predictors)\n\n", ncomp);
}
}

//set max number of relatives
if(strcmp(relfile,"blank")==0)
{
maxpairs=num_samples_use*1000;
if(maxpairs>1e9){maxpairs=1e9;}
}
else
{maxpairs=countrows(relfile);}

//set number of pedigree heritabilities
num_hers1=20;

//set number of vb heritabilities
num_hers2=3;

////////

//do allocations required for testing for structure

sweights=malloc(sizeof(double)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use*num_resps_use);
Yadj=malloc(sizeof(double)*num_samples_use*num_resps_use);
Yadj2=malloc(sizeof(double)*num_samples_use*num_resps_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed*num_resps_use);
thetasds=malloc(sizeof(double)*num_fixed*num_resps_use);
thetapvas=malloc(sizeof(double)*num_fixed*num_resps_use);

if(dichot==1)
{
nullprobs=malloc(sizeof(double)*num_samples_use*num_resps_use);
nullweights=malloc(sizeof(double)*num_samples_use*num_resps_use);
}

gaussian=malloc(sizeof(double)*num_samples_use*nmcmc);

if(ncal>0)	//might not need, but fairly small
{
cindex=malloc(sizeof(int)*data_length);
cindex2=malloc(sizeof(int)*ncal);
cdata=malloc(sizeof(double)*num_samples_use*ncal*num_resps_use);
ccentres=malloc(sizeof(double)*ncal);
cmults=malloc(sizeof(double)*ncal);
csqdevs=malloc(sizeof(double)*ncal);
crates=malloc(sizeof(double)*ncal);
cinfos=malloc(sizeof(double)*ncal);
}

if(ncomp>0)	//again, might not need, but very small
{
eindex=malloc(sizeof(int)*data_length);
eindex2=malloc(sizeof(int)*ncomp);
ecentres=malloc(sizeof(double)*ncomp);
emults=malloc(sizeof(double)*ncomp);
esqdevs=malloc(sizeof(double)*ncomp);
erates=malloc(sizeof(double)*ncomp);
einfos=malloc(sizeof(double)*ncomp);
}

//do initializations required for testing for structure

if(strcmp(sampwfile,"blank")!=0)	//get square roots of sample weights (W^1/2) //not currently used
{
readdoubles=malloc(sizeof(double)*num_samples_use);
read_sampwfile(sampwfile, readdoubles, num_samples_use, ids3, ids1, ids2);

//scale roots so the mean inverse is one (might help stability and helps interpretation)
sum=0;for(i=0;i<num_samples_use;i++){sum+=pow(readdoubles[i],-1);}
mean=sum/num_samples_use;
for(i=0;i<num_samples_use;i++){readdoubles[i]*=mean;}

for(i=0;i<num_samples_use;i++){sweights[i]=pow(readdoubles[i],0.5);}
free(readdoubles);
}
else
{
for(i=0;i<num_samples_use;i++){sweights[i]=1;}
}

//fill Y (maybe multiplying by W^1/2)
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Y[i+m*num_samples_use]=resp[i+m*num_samples_use]*sweights[i];}
}

//fill Z (maybe multiplying by W^1/2)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<num_samples_use;i++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use]*sweights[i];}
}

for(m=0;m<num_resps_use;m++)	//solve null models (and pad any missing values)
{
if(dichot==0)	//get thetas, adjusted response and pad missing
{
reg_covar_lin_missing(Y+m*num_samples_use, Z, num_samples_use, num_fixed, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, Yadj+m*num_samples_use, missingvalue);
}
else		//get thetas, nullprobs, nullweights and adjusted response (Y-mu)/W, padding missing
{
reg_covar_log_missing(Y+m*num_samples_use, Z, num_samples_use, num_fixed, NULL, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, nullprobs+m*num_samples_use, nullweights+m*num_samples_use, Yadj+m*num_samples_use, 0.001, 100, missingvalue, 1);
}

if(verbose==1&&substep<0)	//save coefficients
{
if(mpheno!=-1){sprintf(filename3,"%s.coeff", outfile);}
else{sprintf(filename3,"%s.pheno%d.coeff", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(dichot==0){fprintf(output3, "Component Effect SE P\n");}
else{fprintf(output3, "Component LogOR SE P\n");}
fprintf(output3, "Intercept %.4e %.4e %.4e\n", thetas[0+m*num_fixed], thetasds[0+m*num_fixed], thetapvas[0+m*num_fixed]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.4e %.4e %.4e\n",j, thetas[j+m*num_fixed], thetasds[j+m*num_fixed], thetapvas[j+m*num_fixed]);}
fclose(output3);
}
}	//end of m loop

//fill noise
for(g=0;g<nmcmc;g++)
{
for(i=0;i<num_samples_use;i++){gaussian[i+g*num_samples_use]=rnorm_safe();}
}

if(ncal>0)	//deal with calibration predictors (might not use)
{
//pick predictors at random (some could be trivial)
for(j=0;j<data_length;j++){cindex[j]=j;}
permute_int(cindex,data_length);
qsort(cindex,ncal,sizeof(int), compare_int);

if(chr[cindex[ncal-1]]==chr[cindex[0]]){printf("Warning, all calibration predictors are on Chromosome %d (this indicates that very few predictors are not on this chromosome)\n\n", chr[cindex[0]]);}

for(j=0;j<ncal;j++){cindex2[j]=keeppreds_use[cindex[j]];}

(void)read_data_fly(datafile, dtype, cdata, NULL, num_samples_use, keepsamps, 0, ncal, cindex2, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(cdata, ccentres, cmults, csqdevs, crates, cinfos, num_samples_use, ncal, missingvalue, -1, 0, 0, NULL, 1);

count=0;for(j=0;j<ncal;j++){count+=(cmults[j]==-9999);}
if(count==ncal){printf("Error, all of the %d calibration predictors are trivial\n\n", ncal);}
if(count>0){printf("Warning, %d of the %d calibration predictors are trivial\n\n", count, ncal);}

if(strcmp(sampwfile,"blank")!=0)	//multiply by W^1/2 - not currently used
{copy_matrix(num_samples_use, ncal, cdata, cdata, 1, sweights);}

for(m=1;m<num_resps_use;m++)	//make copies of cdata
{copy_matrix(num_samples_use, ncal, cdata, cdata+m*ncal*num_samples_use, 0, NULL);}

if(num_fixed>1)	//adjust predictors for covariates
{
for(m=0;m<num_resps_use;m++)
{
if(dichot==0){reg_covar_matrix(cdata+m*ncal*num_samples_use, Z, num_samples_use, ncal, num_fixed);}
else{reg_covar_weighted(cdata+m*ncal*num_samples_use, Z, num_samples_use, ncal, num_fixed, nullweights+m*num_samples_use);}
}
}
}

//set starting values for useped and highstruct
useped=fastgwa;
highstruct=0;

if(checkped==1)	//test for structure, and maybe make pedigree prs / perform fastGWA (only here if performing LOCO)
{
if(substep<0)
{
#include "findrels.c"
}
else	//can read from file
{
sprintf(filename,"%s.structure",outfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step CUT\"\n\n", filename);exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input,"Average_Sq_Correlation %lf ", &value)!=1)
{printf("Error reading Row 1 of %s, suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step CUT\"\n\n", filename);exit(1);}
if(fscanf(input,"ChiSq_Statistic %lf ", &value2)!=1)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step CUT\"\n\n", filename);exit(1);}
if(fscanf(input,"MAI %lf ", &value4)!=1)
{printf("Error reading Row 3 of %s, suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--sub-step CUT\"\n\n", filename);exit(1);}
fclose(input);

value3=erfc(pow(value2,.5)*M_SQRT1_2);
printf("Average cross-chromosome squared correlation is %.4e (chi-squared test statistic %.4f; p-value %.4e); estimated maximum average inflation of test statistics is %.4f\n\n", value, value2, value3, value4);

if(value3>=maisig||value4<=maithresh)	//low structure
{printf("This is relatively low structure, so will set the scaling factor to one\n\n");}
else	//high structure (must have dichot=0)
{
if(mode==151){printf("This is relatively high structure, so will use the GRAMMAR-GAMMA test statistic scaling factor\n\n");}
else{printf("This is relatively high structure, so will estimate the test statistic scaling factor (instead of setting it to one)\n\n");}
highstruct=1;
}
}
}
else
{
if(verbose==1)
{
sprintf(filename,"%s.structure",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Not Tested\n");
fclose(output);
}
}

///////////////////////////

if(useped==0)	//continuing
{
//set and save parameters (will set num_pows and num_try, and allocate powers tryps and tryf2s)
#include "setparamsa.c"

if(loco==1)	//set usecal, useridge and usecomp
{
//decide if using calibration predictors (normally only when high structure, but also when not checking pedigree)
if(highstruct==1||checkped==0){usecal=1;}
else{usecal=0;}

//decide if using ridge models (required if using calibration, unless restricting to ridge models)
if(usecal==1&&mode!=151){useridge=1;}
else{useridge=0;}

//decide if using comparison predictors (whenever using loco, unless usecal=1 and mode=151)
if(usecal==1&&mode==151){usecomp=0;}
else{usecomp=1;}
}
else{usecal=0;useridge=0;usecomp=0;}

if(revher==-9999)	//decide whether to revise her using MCMC REML
{
if(her==-9999&&(usecal==1||num_samples_use<40000||strcmp(sampwfile,"blank")!=0)){revher=1;}
else{revher=0;}
}

if(skipcv==0)	//set num_small and num_full (can only have usecal=1 and/or usecomp=1 with loco=1)
{
num_small=num_try+revher*num_hers2*(1+nmcmc);
num_full=1+loco*num_chr+usecal*ncal+useridge*num_chr;
}
else	//only set these so allocations work
{
num_small=num_try;
num_full=num_try;
}

//allocate most variables (will do rhe-specific and larger vb-specific variables later)

total=num_small*num_resps_use;
if(num_full*num_resps_use>total){total=num_full*num_resps_use;}

//this allows for allocations of data, data2 and cors
if(dichot==0){data_warn(bitsize,2*num_samples_use+data_length);}
else{data_warn(bitsize,2*num_samples_use+data_length*num_resps_use);}

//estimate memory for rhe (ignores allocation of data and data2)
if(her!=-9999){value=0;}
else{value=(double)num_samples_use/1024/1024/1024*8*(nmcmc+num_resps_use)*ndivs*num_pows;}

//estimate memory for vb (ignores allocation of data)
value2=(double)data_length/1024/1024/1024*8*(2*total+2*num_resps_use)+num_samples_use/1024/1024/1024*8*(total+2*num_resps_use);

max=value;
if(value2>max){max=value2;}
if(max>1){printf("Warning, to perform the analysis requires approximately %.1f Gb; sorry, this can not be reduced\n\n", max);}

data=malloc(sizeof(double)*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);
if(dichot==0)
{
datasqs=malloc(sizeof(double)*data_length);
cors=malloc(sizeof(double)*bitsize*data_length);
}
else
{
datasqs=malloc(sizeof(double)*data_length*num_resps_use);
cors=malloc(sizeof(double)*bitsize*data_length*num_resps_use);
}

if(skipcv==0)
{
chrindex=malloc(sizeof(int)*num_full);
chrprops=malloc(sizeof(double)*num_full*num_resps_use);
}

YTdata=malloc(sizeof(double)*bitsize*total);

hers=malloc(sizeof(double)*num_resps_use);
hers2=malloc(sizeof(double)*num_resps_use);
hersold=malloc(sizeof(double)*num_resps_use);
exps=malloc(sizeof(double)*data_length*num_resps_use);

if(revher==1)
{
tryhers=malloc(sizeof(double)*num_hers2*num_resps_use);
polates=malloc(sizeof(double)*num_hers2*num_resps_use);
}

lambdas=malloc(sizeof(double)*num_try);
lambdas2=malloc(sizeof(double)*num_try);
lambdas3=malloc(sizeof(double)*num_try);
lambdas4=malloc(sizeof(double)*num_try);

pens=malloc(sizeof(double)*total);
likes=malloc(sizeof(double)*total);
likesold=malloc(sizeof(double)*total);

cgammas=malloc(sizeof(double)*num_resps_use);
csds=malloc(sizeof(double)*num_resps_use);
ceffs=malloc(sizeof(double)*num_resps_use);

bitrun=malloc(sizeof(int)*bittotal*num_resps_use);
bitdo=malloc(sizeof(int)*bittotal);
bitactive=malloc(sizeof(int)*total);
bitdet1=malloc(sizeof(int)*bittotal);
bitdet2=malloc(sizeof(int)*bittotal);
bitdiffs=malloc(sizeof(double)*bittotal*total);
bitpens=malloc(sizeof(double)*bittotal*total);

Mtops=malloc(sizeof(int)*num_resps_use);
Mbests=malloc(sizeof(int)*num_resps_use);
Mincs=malloc(sizeof(int)*num_resps_use);
Mscales=malloc(sizeof(double)*num_resps_use);
Mmses=malloc(sizeof(double)*num_resps_use);
Mmses2=malloc(sizeof(double)*num_resps_use);
Mneffs=malloc(sizeof(double)*num_resps_use);

if(multi==1)
{
Umat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Umat2=malloc(sizeof(double)*num_resps_use*num_resps_use);
}

if(skipcv==0)	//sort out chrindexes - the chromosome corresponding to each phenotype/calibration residual
{
//start with genome-wide model
chrindex[0]=-1;

if(loco==1)	//add loco models
{
chrindex[1]=chr[0];
count=2;
for(j=1;j<data_length;j++)
{
if(chr[j]!=chrindex[count-1]){chrindex[count]=chr[j];count++;}
}

if(usecal==1)	//now add calibration models
{
for(j=0;j<ncal;j++){chrindex[1+num_chr+j]=chr[cindex[j]];}
}

if(useridge==1)	//and ridge loco models
{
chrindex[1+num_chr+ncal]=chr[0];
count=1+num_chr+ncal+1;
for(j=1;j<data_length;j++)
{
if(chr[j]!=chrindex[count-1]){chrindex[count]=chr[j];count++;}
}
}
}
}

////////

if(dichot==0)	//will use scaled phenotypes
{
for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

Mscales[m]=pow(var,-0.5);
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=(Yadj[i+m*num_samples_use]-mean)*Mscales[m];}
}
}
else	//use scaled phenotypes and scaled weights
{
for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

Mscales[m]=pow(var,-0.5);
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=(Yadj[i+m*num_samples_use]-mean)*Mscales[m];}

sum=0;for(i=0;i<num_samples_use;i++){sum+=pow(nullweights[i+m*num_samples_use],-1);}
mean=sum/num_samples_use;
for(i=0;i<num_samples_use;i++){nullweights[i+m*num_samples_use]*=mean;}
}
}

//save a copy of phenotypes (if RHE finds big predictors, these are regressed out of Yadj2)
for(m=0;m<num_resps_use;m++) 
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use];}
}

if(her!=-9999)	//have her, power and maybe indhers - get exps and correlations (must have multi=0)
{
for(m=0;m<num_resps_use;m++){Mtops[m]=0;hers[m]=her;}

#include "getindhersa.c"
}
else	//do not have her (and probably not power)
{
if(multi==0)	//basic analysis - estimate her and power, and get exps and correlations
{
#include "getindhersb.c"
}
else	//multivariate analysis - substep is either -9999, -1 or a positive integer
{
if(substep<0)	//estimates her and power, multivariate terms, and get exps and correlations
{
#include "getindhersb.c"

if(substep==-1)	//save values
{
sprintf(filename2,"%s.multivariate.bin", outfile);
if((output2=fopen(filename2,"ab"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

fwrite(centres, sizeof(double), data_length, output2);
fwrite(mults, sizeof(double), data_length, output2);
fwrite(sqdevs, sizeof(double), data_length, output2);
fwrite(rates, sizeof(double), data_length, output2);
fwrite(infos, sizeof(double), data_length, output2);

fwrite(Mtops, sizeof(int), num_resps_use, output2);
fwrite(hers, sizeof(double), num_resps_use, output2);
fwrite(hers2, sizeof(double), num_resps_use, output2);
fwrite(exps, sizeof(double), data_length*num_resps_use, output2);

fwrite(Umat, sizeof(double), num_resps_use*num_resps_use, output2);
fwrite(Umat2, sizeof(double), num_resps_use*num_resps_use, output2);
fwrite(Yadj, sizeof(double), num_samples_use*num_resps_use, output2);
fwrite(Yadj2, sizeof(double), num_samples_use*num_resps_use, output2);

fwrite(datasqs, sizeof(double), data_length, output2);
fwrite(cors, sizeof(double), bitsize*data_length, output2);

fclose(output2);
}
}
else	//read values from file
{
sprintf(filename2,"%s.multivariate.bin", outfile);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

fseeko(input2, (off_t)sizeof(int)*(6+num_samples_use), SEEK_SET);

if(fread(centres, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading predictor means from %s\n\n", filename2);exit(1);}
if(fread(mults, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading predictor scalings from %s\n\n", filename2);exit(1);}
if(fread(sqdevs, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading predictor variances from %s\n\n", filename2);exit(1);}
if(fread(rates, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading predictor call rates from %s\n\n", filename2);exit(1);}
if(fread(infos, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading predictor info scores from %s\n\n", filename2);exit(1);}

if(fread(Mtops, sizeof(int), num_resps_use, input2)!=num_resps_use)
{printf("Error reading best powers from %s\n\n", filename2);exit(1);}
if(fread(hers, sizeof(double), num_resps_use, input2)!=num_resps_use)
{printf("Error reading phenotype heritabilities from %s\n\n", filename2);exit(1);}
if(fread(hers2, sizeof(double), num_resps_use, input2)!=num_resps_use)
{printf("Error reading phenotype heritabilities from %s\n\n", filename2);exit(1);}
if(fread(exps, sizeof(double), data_length*num_resps_use, input2)!=data_length*num_resps_use)
{printf("Error reading per-predictor heritabilities from %s\n\n", filename2);exit(1);}

if(fread(Umat, sizeof(double), num_resps_use*num_resps_use, input2)!=num_resps_use*num_resps_use)
{printf("Error reading transformation matrix from %s\n\n", filename2);exit(1);}
if(fread(Umat2, sizeof(double), num_resps_use*num_resps_use, input2)!=num_resps_use*num_resps_use)
{printf("Error reading transformation inverse from %s\n\n", filename2);exit(1);}
if(fread(Yadj, sizeof(double), num_samples_use*num_resps_use, input2)!=num_samples_use*num_resps_use)
{printf("Error reading first transformed phenotypes from %s\n\n", filename2);exit(1);}
if(fread(Yadj2, sizeof(double), num_samples_use*num_resps_use, input2)!=num_samples_use*num_resps_use)
{printf("Error reading second transformed phenotypes from %s\n\n", filename2);exit(1);}

if(fread(datasqs, sizeof(double), data_length, input2)!=data_length)
{printf("Error reading predictor variances from %s\n\n", filename2);exit(1);}
if(fread(cors, sizeof(double), bitsize*data_length, input2)!=bitsize*data_length)
{printf("Error reading predictor-predictor correlations from %s\n\n", filename2);exit(1);}

fclose(input2);

printf("Have read details from \"--sub-step CUT\"\n");
printf("Phenotype %d has estimated heritability %.4f\n\n", substep, hers[substep-1]);
}
}
}

time(&midtime);
printf("Time check: have so far spent %.2f hours\n\n", (double)(midtime-starttime)/60/60);

if(substep<0){sprintf(filename,"%s.progress",outfile);}
else{sprintf(filename,"%s.pheno%d.progress",outfile, substep);}
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Time check: have so far spent %.2f hours\n", (double)(midtime-starttime)/60/60);
fclose(output);

if(skipcv==0)	//sort out chrprops
{
//chrprops is sum of exps not on target chromosome for test models
for(m=0;m<num_resps_use;m++)
{
for(p=0;p<num_full;p++)
{
sum=0;for(j=0;j<data_length;j++){sum+=exps[j+m*data_length]*(chr[j]!=chrindex[p]);}
chrprops[p+m*num_full]=sum;
}
}
}

if(revher==1)	//work out what heritabilities to test when performing MC REML
{
for(m=0;m<num_resps_use;m++)
{
if(hers[m]>0.01)
{tryhers[m*num_hers2]=0.5*hers[m];tryhers[1+m*num_hers2]=0.75*hers[m];tryhers[2+m*num_hers2]=hers[m];}
else
{tryhers[m*num_hers2]=0.01;tryhers[1+m*num_hers2]=0.03;tryhers[2+m*num_hers2]=0.05;}
}
}

if(dougvar==99){exit(1);}

////////

//set lambdas for effect sizes (values corresponding to exp(betaj^2)=1 - will later scale)
for(p=0;p<num_try;p++)
{
//start at zero, then change when required
lambdas[p]=0;lambdas2[p]=0;lambdas3[p]=0;lambdas4[p]=0;

if(mode==151)	//ridge - exp(betaj^2) = lam
{lambdas[p]=1.0;}

if(mode==152)	//bolt - exp(betaj^2) = plam + p2lam2, with p2lam2/plam = f2/(1-f2) - p and p2 within (0,1)
{
lambdas[p]=(1-tryf2s[p])/tryps[p];
lambdas2[p]=tryf2s[p]/tryp2s[p];
}

if(mode==153)	//bayesr - exp(betaj^2) = plam + p2lam2 + p3lam3 + p4lam4
{
if(pointmass==1)	//sparse version
{
value=tryp2s[p]/100+tryp3s[p]/10+tryp4s[p];
lambdas[p]=0.0;
lambdas2[p]=0.01/value;
lambdas3[p]=0.1/value;
lambdas4[p]=1.0/value;
}
else
{
value=tryps[p]/1000+tryp2s[p]/100+tryp3s[p]/10+tryp4s[p];
lambdas[p]=0.001/value;
lambdas2[p]=0.01/value;
lambdas3[p]=0.1/value;
lambdas4[p]=1.0/value;
}
}

if(mode==154)	//elastic - exp(betaj^2) = 2p/lam^2 + 2p2/lam2^2 + p3lam3 (alt 2(p+p2)/lam^2 + p3lam3)
{
if(tryp3s[p]==0)	//lasso (will have p=p2=1/2)
{lambdas[p]=1.0;lambdas2[p]=1.0;}
if(tryp3s[p]==1)	//ridge
{lambdas3[p]=1.0;}
if(tryp3s[p]>0||tryp3s[p]<1)	//elastic
{
lambdas[p]=pow(4*tryps[p]/(1-tryf2s[p]),.5);
lambdas2[p]=pow(4*tryp2s[p]/(1-tryf2s[p]),.5);
lambdas3[p]=tryf2s[p]/tryp3s[p];
}
}
}	//end of p loop

////////

/*
if(dougvar==98)	//compute invV Y by conjugate gradient and vb
{
#include "contest.c"
exit(1);
}
*/

//more allocations (must refresh total because likely changed in getindhers)

total=num_small*num_resps_use;
if(num_full*num_resps_use>total){total=num_full*num_resps_use;}

effs=malloc(sizeof(double)*data_length*total);
effs2=malloc(sizeof(double)*data_length*num_resps_use);
effs3=malloc(sizeof(double)*data_length*num_resps_use);
probs=malloc(sizeof(double)*data_length*total);
residuals=malloc(sizeof(double)*num_samples_use*total);
residuals2=malloc(sizeof(double)*num_samples_use*total);
residuals3=malloc(sizeof(double)*num_samples_use*num_resps_use);
changes=malloc(sizeof(double)*bitsize*total);

if(substep<0&&loco==1)	//save root file
{
sprintf(filename3,"%s.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Analysis KVIK\n");
if(dichot==0){fprintf(output3,"Regression Linear\n");}
else{fprintf(output3,"Regression Logistic\n");}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Phenotypes %s\n", respfile);
if(strcmp(covarfile,"blank")!=0){fprintf(output3,"Covariates %s\n", covarfile);}
else{fprintf(output3,"Covariates none\n");}
if(strcmp(factorfile,"blank")!=0){fprintf(output3,"Factors %s\n", factorfile);}
else{fprintf(output3,"Factors none\n");}
if(mpheno!=-1){fprintf(output3,"Phenotype_Number %d\n", mpheno);}
else{fprintf(output3,"Num_Phenotypes %d\n", num_resps_use);}
fprintf(output3,"Num_Samples_Used %d\n", num_samples_use);
fprintf(output3,"Num_Predictors_Used %d\n", data_length);
fprintf(output3,"Num_Chromosomes %d\n", num_chr);
if(multi==0){fprintf(output3,"Multivariate NO\n");}
else{fprintf(output3,"Multivariate YES\n");}
fclose(output3);
}

if(substep==-9999)
{
#include "quickpred.c"
}
else	//either have substep=-1 (cutting) or >0 (analyzing phenotype substep)
{
if(substep==-1)
{
printf("Analysis completed; please rerun this command once for each phenotype (i.e., first using \"--sub-step 1\", then \"--sub-step 2\", etc)\n\n");
}
else
{
//convenient to set hers to zero for phenotypes we are ignoring
for(m=0;m<num_resps_use;m++)
{
if(substep!=m+1){hers[m]=0;}
}

#include "quickpred.c"

printf("Analysis completed for Phenotype %d; once all phenotypes are finished, please rerun this command using \"--sub-step JOIN\"\n\n", substep);
}
}

free(tryps);free(tryp2s);free(tryp3s);free(tryp4s);free(tryf2s);
free(data);free(data2);free(datasqs);free(cors);
if(skipcv==0){free(chrindex);free(chrprops);}
free(YTdata);
free(hers);free(hers2);free(hersold);free(exps);
if(revher==1){free(tryhers);free(polates);}
free(lambdas);free(lambdas2);free(lambdas3);free(lambdas4);
free(pens);free(likes);free(likesold);
free(cgammas);free(csds);free(ceffs);
free(bitrun);free(bitdo);free(bitactive);free(bitdet1);free(bitdet2);free(bitdiffs);free(bitpens);
free(Mtops);free(Mbests);free(Mincs);free(Mscales);free(Mmses);free(Mmses2);free(Mneffs);
if(multi==1){free(Umat);free(Umat2);}

free(effs);free(effs2);free(effs3);free(probs);free(residuals);free(residuals2);free(residuals3);free(changes);
}	//end of continuing

////////

if(skipcv==0){free(keeptrain);free(keeptest);}
free(powers);
free(sweights);free(Y);free(Yadj);free(Yadj2);free(Z);
free(thetas);free(thetasds);free(thetapvas);
if(dichot==1){free(nullprobs);free(nullweights);}
free(gaussian);
if(ncal>0){free(cindex);free(cindex2);free(cdata);free(ccentres);free(cmults);free(csqdevs);free(crates);free(cinfos);}
if(ncomp>0){free(eindex);free(eindex2);free(ecentres);free(emults);free(esqdevs);free(erates);free(einfos);}

///////////////////////////

