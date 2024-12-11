/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Computing kinships - normal case

///////////////////////////

if(dosage!=-9999)	//get list of males
{
keepsamps_males=malloc(sizeof(int)*num_samples_use);

count=countrows(malesfile);
printf("Reading list of %d samples from %s\n", count, malesfile);
wantids=malloc(sizeof(char*)*count);
read_ids(malesfile, NULL, NULL, wantids, count, NULL, 0, 0);

count2=find_strings(wantids, count, ids3, num_samples_use, NULL, keepsamps_males, NULL, NULL, NULL, NULL, 1);
if(count2==0){printf("Error, none of the %d samples are in %s\n\n", num_samples_use, malesfile);exit(1);}
if(count2==num_samples_use){printf("Error, all of the %d samples are in  %s\n\n", num_samples_use, malesfile);exit(1);}
if(count2==count){printf("All males are in the data\n\n");}
else{printf("%d of the males are in the data\n\n", count2);}

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
}

//allocate variables (and set kins to zero)

if(single==0)
{
data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
}
else
{
data_warn2(bitsize,num_samples_use*3/2);
data=malloc(sizeof(double)*num_samples_use*bitsize);
data_single=malloc(sizeof(float)*num_samples_use*bitsize);
}

if(onlydets==0)
{
if(single==0)
{
anal_warn(num_samples_use, num_samples_use);
kins=calloc((size_t)num_samples_use*num_samples_use,sizeof(double));
}
else
{
anal_warn(num_samples_use, num_samples_use/2);
kins_single=calloc((size_t)num_samples_use*num_samples_use,sizeof(single));
}
}

exps=malloc(sizeof(double)*data_length);

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
if(mode==112){sprintf(filename,"%s/progress.%d",folder,partition);}
else{sprintf(filename,"%s.progress",outfile);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(onlydets==0){printf("Calculating kinships for Chunk %d of %d\n", bit+1, bittotal);}
else{printf("Calculating details for Chunk %d of %d\n", bit+1, bittotal);}
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(onlydets==0){fprintf(output,"Calculating kinships for Chunk %d of %d\n", bit+1, bittotal);}
else{fprintf(output,"Calculating details for Chunk %d of %d\n", bit+1, bittotal);}
fclose(output);

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);

if(dosage!=-9999)	//change for males
{
wcount=0;
for(j=0;j<bitlength;j++)
{
for(i=0;i<count2;i++)
{
i2=keepsamps_males[i];
if(data[(size_t)j*num_samples_use+i2]==1){wcount++;}
data[(size_t)j*num_samples_use+i2]*=dosage/2;
}
}
if(wcount>0){printf("Warning, there were %d ocassions when a predictor had value 1 for a male\n\n", wcount);}
}

if(david==0)
{
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, power, strcmp(centresfile,"blank")!=0, hwestand, weights+bitstart, 1);
}
else
{
//first subtract one and set missing to zero
#pragma omp parallel for private(j,i) schedule (static)
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue){data[(size_t)j*num_samples_use+i]--;}
else{data[(size_t)j*num_samples_use+i]=0;}
}
}

//now compute Ms for all predictors, and save in mults
#pragma omp parallel for private(j,sum,sumsq,i) schedule (static)
for(j=0;j<bitlength;j++)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++)
{sum+=data[(size_t)j*num_samples_use+i];sumsq+=pow(data[(size_t)j*num_samples_use+i],2);}
mults[bitstart+j]=0.5+.5*(pow(sum,2)-sumsq)/num_samples_use/(num_samples_use-1);
}

//finally scale predictors by (2(1-Ms))^-.5
#pragma omp parallel for private(j,i,value) schedule (static)
for(j=0;j<bitlength;j++)
{
value=pow(2*(1-mults[bitstart+j]),-.5);
if(value!=value){printf("Error snp %d and %f\n", j+1, mults[bitstart+j]);exit(1);}
for(i=0;i<num_samples_use;i++){data[(size_t)j*num_samples_use+i]*=value;}
}
}

if(onlydets==0)
{
if(single==0)
{
alpha=1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &bitlength, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, kins, &num_samples_use);
}
else
{
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use;i++){data_single[(size_t)j*num_samples_use+i]=(float)data[(size_t)j*num_samples_use+i];}
}

alpha_single=1.0;beta_single=1.0;
sgemm_("N", "T", &num_samples_use, &num_samples_use, &bitlength, &alpha_single, data_single, &num_samples_use, data_single, &num_samples_use, &beta_single, kins_single, &num_samples_use);
}
}
}	//end of bit loop
printf("\n");

count=0;for(j=0;j<data_length;j++){count+=(mults[j]!=-9999);}
if(count==0){printf("Error, all %d predictors trivial\n\n", data_length);exit(1);}

if(david==1)	//have (X-1)^TD(X-1) - must have single=0
{
//compute sum (.5-Ms)/(1-Ms)
value=0;
for(j=0;j<data_length;j++){value+=(.5-mults[j])/(1-mults[j]);}
printf("value is %f\n", value);

//add value to kinshps, then divide by number of predictors
#pragma omp parallel for private(i2,i) schedule (static)
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++)
{kins[(size_t)i2*num_samples_use+i]=(kins[(size_t)i2*num_samples_use+i]+value)/data_length;}
}
}

//compute expected heritabilities (will scale so they equal one when saving)
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}
else{exps[j]=0;}
}

if(onlydets==0)	//save kins
{
if(mode==112){sprintf(outfile,"%skinships.%d", folder, partition);}
if(single==0){write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 1, preds, keeppreds_use, centres, mults, weights, al1, al2, exps, data_length, datafile, power, kingz, kinraw, 1);}
else{write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 1, preds, keeppreds_use, centres, mults, weights, al1, al2, exps, data_length, datafile, power, kingz, kinraw, 1);}
}
else	//only save details
{
//get sum of exps
value=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){value+=exps[j];}
}

sprintf(filename2,"%s.grm.details", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor Index Centre Scaling Weight A1 A2 Exp_Heritability Share\n");
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
fprintf(output2, "%s %d %.6f %.6f %.6f %c %c %.6f %.6e\n", preds[j], keeppreds_use[j]+1, centres[j], mults[j], weights[j], al1[j], al2[j], exps[j], exps[j]/value);
}
}
fclose(output2);
}

if(dosage!=-9999){free(keepsamps_males);}
free(data);if(single==1){free(data_single);}
if(single==0){free(kins);}else{free(kins_single);}
free(exps);
if(binary==0){gzclose(datainputgz);}

///////////////////////////

