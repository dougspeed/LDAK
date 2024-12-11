/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//compute correlations between two sets of predictors

///////////////////////////

if(bitsize>numa){bitsize=numa;}

//allocate variables

data_warn3(bitsize+numb,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*numb);

Z=malloc(sizeof(double)*num_samples_use*num_fixed);

anal_warn2(bitsize, numb);
cors=malloc(sizeof(double)*numb*bitsize);

acentres=malloc(sizeof(double)*numa);
amults=malloc(sizeof(double)*numa);
asqdevs=malloc(sizeof(double)*numa);
arates=malloc(sizeof(double)*numa);
ainfos=malloc(sizeof(double)*numa);
bcentres=malloc(sizeof(double)*numb);
bmults=malloc(sizeof(double)*numb);
bsqdevs=malloc(sizeof(double)*numb);
brates=malloc(sizeof(double)*numb);
binfos=malloc(sizeof(double)*numb);

apreds=malloc(sizeof(char*)*numa);
bpreds=malloc(sizeof(char*)*numb);

//fill covariates and predictor names

for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}
for(j=0;j<numa;j++){apreds[j]=allpreds[keepa[j]];}
for(j=0;j<numb;j++){bpreds[j]=allpreds[keepb[j]];}

//deal with progress and on-the-fly files
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%s.rjk2.average", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor\tAverage_Squared_Correlation\tNumber_Pairs\n");
fclose(output2);

if(savepairs==1)
{
sprintf(filename3,"%s.rjk.full", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fclose(output3);
}

//read and standardize data for listb (maybe regressing out covariates)
printf("Reading data for %d listb predictors\n", numb);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Reading data for %d listb predictors\n", numb);
fclose(output);

if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

current=read_data_fly(datafile, dtype, data2, NULL, num_samples_use, keepsamps, 0, numb, keepb, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);

stand_data(data2, bcentres, bmults, bsqdevs, brates, binfos, num_samples_use, numb, missingvalue, -1, 0, 0, NULL, 1);

if(num_covars>1)	//regress out covariates - must then standardize again
{
reg_covar_matrix(data2, Z, num_samples_use, numb, num_covars);
stand_matrix_nomiss(data2, num_samples_use, num_samples_use, numb);
}

if(binary==0){gzclose(datainputgz);}
printf("\n");

////////

//prepare for reading data (again)
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//will save sum of all squared correlations
sum2=0;
scount=0;

//ready for bit loop
bittotal=(numa-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>numa){bitend=numa;}
bitlength=bitend-bitstart;

printf("Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

//read and standardize bit of data for lista (maybe regressing out covariates)
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keepa, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);

stand_data(data, acentres+bitstart, amults+bitstart, asqdevs+bitstart, arates+bitstart, ainfos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1);

if(num_covars>1)	//regress out covariates - must then standardize again
{
reg_covar_matrix(data, Z, num_samples_use, bitlength, num_covars);
stand_matrix_nomiss(data, num_samples_use, num_samples_use, bitlength);
}

//compute correlations
alpha=1.0/num_samples_use;beta=0.0;
dgemm_("T", "N", &numb, &bitlength, &num_samples_use, &alpha, data2, &num_samples_use, data, &num_samples_use, &beta, cors, &numb);

//save average squares (restricting to different chromosomes)
count2=0;
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
for(j=0;j<bitlength;j++)
{
if(amults[bitstart+j]!=-9999)
{
sum=0;count=0;
for(j2=0;j2<numb;j2++)
{
if(bmults[j2]!=-9999&&allchr[keepa[bitstart+j]]!=allchr[keepb[j2]])
{sum+=pow(cors[(size_t)j*numb+j2],2);count++;}
}
if(count>0){fprintf(output2,"%s\t%.6e\t%d\n", allpreds[keepa[bitstart+j]], sum/count, count);}
else{fprintf(output2,"%s\tNA\t0\n", allpreds[keepa[bitstart+j]]);count2++;}
sum2+=sum;
scount+=count;
}
else
{fprintf(output2,"%s\tNA\tNA\n", allpreds[keepa[bitstart+j]]);count2++;}
}
fclose(output2);

if(count2>0){printf("Error, there were no different-chromosome listb predictors for %d of the lista predictors\n\n", count2);exit(1);}

if(savepairs==1)
{
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(j=0;j<bitlength;j++)
{
sum=0;count=0;
for(j2=0;j2<numb;j2++){fprintf(output3,"%.6e ",cors[(size_t)j*numb+j2]);}
fprintf(output3,"\n");
}
fclose(output3);
}
}	//end of bit loop
printf("\n");

if(binary==0){gzclose(datainputgz);}

sprintf(filename4,"%s.lista.used", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(j=0;j<numa;j++){fprintf(output4,"%s\n", apreds[j]);}
fclose(output4);

sprintf(filename5,"%s.listb.used", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(j=0;j<numb;j++){fprintf(output5,"%s\n", bpreds[j]);}
fclose(output5);

printf("The average squared correlations are saved in %s", filename2);
if(savepairs==1){printf(", with all squared correlations saved in %s", filename3);}
printf("\n\n");

//work out overall average, then compare to mean under the null distribution
value=sum2/scount;
mean=1.0/(num_samples_use-1);
var=2*(num_samples_use-2)*pow(num_samples_use-1,-2)/(num_samples_use+1)/scount;
value2=pow(value-mean,2)/var;
value3=erfc(pow(value2,.5)*M_SQRT1_2);

printf("Average (different-chromosome) squared correlation is %4e; chi-squared test statistic %4f; p-value %4e; estimated maximum average inflation of test statistics is %f\n\n", value, value2, value3, (value-mean)*num_samples_use);

////////

//free allocations from setdl.c
free(keepa);
free(keepb);

//frees from above
free(data);free(data2);
free(Z);
free(cors);
free(acentres);free(amults);free(asqdevs);free(arates);free(ainfos);
free(bcentres);free(bmults);free(bsqdevs);free(brates);free(binfos);
free(apreds);free(bpreds);

///////////////////////////

