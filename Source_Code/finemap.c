/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Fine-mapping
//start out assuming ind-level data, then later modify
//have summary statistics

///////////////////////////

if(cvar==-9999)	//set prior variance based on sample size
{cvar=29.72/num_samples_use;}
printf("The standardized effect sizes of causal predictors are assumed to come from a Gaussian distribution with mean zero and variance %.4e; change the latter using \"--causal-variance\"\n\n", cvar);

//work out how many bits and bitmax

have set window_cm or something

//will later read in prs if using (will only have one chromosome) - hope ok to not pad

prs=malloc(sizeof(double)*num_samples_use);
for(i=0;i<num_samples_use;i++){prs[i]=0;}

//allocate variables

data_warn(bitmax,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitmax);

Y=malloc(sizeof(double)*num_samples_use);
Yadj=malloc(sizeof(double)*num_samples_use);
Yadj2=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

effs=malloc(sizeof(double)*data_length*2);
effs2=malloc(sizeof(double)*data_length);
probs=malloc(sizeof(double)*data_length);
residuals=malloc(sizeof(double)*num_samples_use);
residuals2=malloc(sizeof(double)*num_samples_use);
changes=malloc(sizeof(double)*bitsize*2);

//fill Y (subtracting prs)
for(i=0;i<num_samples_use;i++){Y[i+m*num_samples_use]=resp[i];}

//fill Z
for(j=0;j<num_fixed;j++)
{
for(i=0;i<num_samples_use;i++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

//solve null model (and pad any missing values)
if(dichot==0)	//get thetas, adjusted response and pad missing
{
reg_covar_lin_missing(Y, Z, num_samples_use, num_fixed, thetas, thetasds, thetapvas, Yadj, missingvalue);
}
else		//get thetas, nullprobs, nullweights and adjusted response, padding missing
{
reg_covar_log_missing(Y, Z, num_samples_use, num_fixed, NULL, thetas, thetasds, thetapvas, nullprobs, nullweights, Yadj, 0.001, 100, missingvalue, 1);
}

//subtract prs

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);






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
hers2=malloc(sizeof(double)*num_resps_use);
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

//save a copy of phenotypes (Yadj2 will always match Yadj, unless RHE finds big predictors, in which case these are regressed out of Yadj2)
for(m=0;m<num_resps_use;m++) 
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use];}
}

if(her!=-9999)	//have her and power - set exps based on indhers or by scaling data (also compute correlations)
{
for(m=0;m<num_resps_use;m++){Mtops[m]=0;hers[m]=her;}

#include "getindhersa.c"
}
else	//estimate hers, exps and probably power, and if multi=1 get multivariate terms (also compute statistics and correlations)
{
if(substep<0)	//compute from scratch
{
#include "getindhersb.c"

if(substep==-1)	//save to file
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
else	//read from file
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
{printf("Error reading predictor correlations from %s\n\n", filename2);exit(1);}

fclose(input2);

printf("Have read details from \"--sub-step CUT\"\n");
printf("Phenotype %d has estimated heritability %.4f\n\n", substep, hers[substep-1]);
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

if(mode==154)	//elastic - exp(betaj^2) = 2p/lam^2 + 2p2/lam2^2 + p3lam3
{
if(tryp3s[p]==0)	//lasso
{lambdas[p]=1.0;lambdas2[p]=1.0;}
if(tryp3s[p]==1)	//ridge
{lambdas3[p]=1.0;}
if(tryp3s[p]>0||tryp3s[p]<1)	//elastic
{
lambdas[p]=pow(2*(1-tryp3s[p])/(1-tryf2s[p]),.5);
lambdas2[p]=pow(2*(1-tryp3s[p])/(1-tryf2s[p]),.5);
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
if(strcmp(sampwfile,"blank")==0){fprintf(output3,"Sample_Weights NO\n");}
else{fprintf(output3,"Sample_Weights YES\n");}
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
free(hers);free(hersold);free(exps);
if(revher==1){free(tryhers);free(polates);}
free(lambdas);free(lambdas2);free(lambdas3);free(lambdas4);
free(pens);free(likes);free(likesold);
free(cgammas);free(csds);free(ceffs);
free(bitrun);free(bitdo);free(bitactive);free(bitdet1);free(bitdet2);free(bitdiffs);free(bitpens);
free(Mtops);free(Mbests);free(Mincs);free(Mscales);free(Mmses);free(Mmses2);free(Mneffs);
if(multi==1){free(hers2);free(Umat);free(Umat2);}

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

