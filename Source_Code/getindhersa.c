/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Get per-predictor heritabilities based on indhers or by scaling data - also get cors and datasqs and load up eindexes

///////////////////////////

//a few allocations

R=malloc(sizeof(double)*num_samples_use*num_resps_use);
RTdata=malloc(sizeof(double)*bitsize*num_resps_use);

//fill start of R with adjusted phenotypes
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){R[(size_t)m*num_samples_use+i]=Yadj[i+m*num_samples_use];}
}

if(strcmp(indhers,"blank")==0){printf("Computing per-predictor heritabilities\n\n");}
else{printf("Computing predictor means and variances\n\n");}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Computing per-predictor heritabilities\n");
fclose(output);

//ready for bit loop
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%200==0)
{
printf("Processing predictors in Chunk %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Processing predictors in Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read data, compute statistics, standardize and set missing to zero
if(dtype==1)	//fast way
{(void)read_bed_wrapper(datafile, data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, keepsamps, bitlength, keeppreds_use+bitstart, num_samples, num_preds, missingvalue, bedzeros, bedones, bedtwos, 2, maxthreads);}
else	//slow way
{
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1);
}

if(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999)	//perform qc
{
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
maf=centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]);
value=sqdevs[bitstart+j]/centres[bitstart+j]/(1-centres[bitstart+j]/2);}
if(minmaf!=-9999&&maf<minmaf){mults[bitstart+j]=-9999;}
if(maxmaf!=-9999&&maf>maxmaf){mults[bitstart+j]=-9999;}
if(minvar!=-9999&&sqdevs[bitstart+j]<minvar){mults[bitstart+j]=-9999;}
if(minobs!=-9999&&rates[bitstart+j]<minobs){mults[bitstart+j]=-9999;}
if(mininfo!=-9999&&value<mininfo){mults[bitstart+j]=-9999;}
}
}

//get cors and datasqs
if(dichot==0)	//using non-weighted regression, so only require one set
{
if(num_fixed>1)	//adjust for covariates
{reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);

for(j=bitstart;j<bitend;j++){datasqs[j]=cors[(size_t)j*bitsize+j-bitstart];}
}
else	//using weighted regression, so require one set per phenotype
{
for(m=0;m<num_resps_use;m++)
{
if(num_fixed>1)	//adjust for covariates - its ok if have already adjusted data for covariates
{reg_covar_weighted(data, Z, num_samples_use, bitlength, num_fixed, nullweights+m*num_samples_use);}

//put weighted version of data into data2
copy_matrix(num_samples_use, bitlength, data, data2, 1, nullweights+m*num_samples_use);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data2, &num_samples_use, &beta, cors+(size_t)(bitstart+m*data_length)*bitsize, &bitsize);

for(j=bitstart;j<bitend;j++){datasqs[j+m*data_length]=cors[(size_t)(j+m*data_length)*bitsize+j-bitstart];}
}
}
}	//end of bit loop
printf("\n");

count=0;for(j=0;j<num_preds_use;j++){count+=(mults[j]==-9999);}
if(count==num_preds_use)
{
if(minmaf==-9999&&maxmaf==-9999&&minvar==-9999&&minobs==-9999&&mininfo==-9999)
{printf("Error, all predictors are trivial (showed no variation)\n\n");}
else
{printf("Error, all predictors failed quality control\n\n");}
exit(1);
}
if(count>0)
{
if(minmaf==-9999&&maxmaf==-9999&&minvar==-9999&&minobs==-9999&&mininfo==-9999)
{printf("Warning, %d predictors were excluded because they were trivial (showed no variation)\n\n", count);}
else
{printf("Warning, %d predictors were excluded because they failed quality control\n\n", count);}
}

//get exps for first phenotype based on weights and power
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}
else{exps[j]=0;}
}

//make sure exps sum to one
sum=0;for(j=0;j<data_length;j++){sum+=exps[j];}
for(j=0;j<data_length;j++){exps[j]=exps[j]/sum;}

for(m=1;m<num_resps_use;m++)	//copy exps into other phenotypes
{
for(j=0;j<data_length;j++){exps[j+m*data_length]=exps[j];}
}

////////

free(R);free(RTdata);

///////////////////////////

