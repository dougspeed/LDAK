/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Do multiplications for calculating blups

///////////////////////////

//sort out step
step=(50000/bitsize);if(step<20){step=20;}
if(step>20){step=10*(step/10);}if(step>50){step=20*(step/20);}
if(step>100){step=50*(step/50);}if(step>300){step=100*(step/100);}
if(step>1000){step=500*(step/500);}

total=num_kins+num_scores+(num_tops>0);

//allocate variables
data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

if(num_kins>0){bluprands=malloc(sizeof(double)*num_samples_use*num_kins);}

effects=malloc(sizeof(double*)*(total+3));	//last three are sum, centres and predictor used?
guesses=malloc(sizeof(double*)*(total+3));	//last three are sum of genetics, covs and sum of both
for(k=0;k<total+3;k++)
{effects[k]=malloc(sizeof(double)*data_length);guesses[k]=malloc(sizeof(double)*num_samples_use);}

if(num_kins>0)	//read random effects
{read_rand(blupfile, bluprands, num_kins, num_samples_use, ids3);}

//set effects and guesses to zero
for(k=0;k<total+3;k++)
{
for(j=0;j<data_length;j++){effects[k][j]=0;}
for(i=0;i<num_samples_use;i++){guesses[k][i]=0;}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
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

if(bit%step==0)
{
printf("Calculating blups for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating blups for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);

for(j=0;j<bitlength;j++)	//get raw effect size, then contributions to projections
{
//for kinships, effect size is multsj sum (Xij bluprandsi) / kinsums, where Xij = (Sij-centrej)*multsj
for(k=0;k<num_kins;k++)
{
if(blupfactors[k][bitstart+j]!=0)
{
value=0;
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{value+=(data[(size_t)j*num_samples_use+i]-blupcentres[k][bitstart+j])*bluprands[i+k*num_samples_use];}
}
effects[k][bitstart+j]=value*pow(blupfactors[k][bitstart+j],2)/kinsums[k];
}}

//for regions and tops, already have raw effect sizes
for(k=num_kins;k<total;k++){effects[k][bitstart+j]=blupfactors[k][bitstart+j];}

//now add on contributions
for(k=0;k<total;k++)
{
if(effects[k][bitstart+j]!=0)
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{guesses[k][i]+=(data[(size_t)j*num_samples_use+i]-blupcentres[k][bitstart+j])*effects[k][bitstart+j];}
}
}}
}	//end of j loop
}	//end of bit loop
printf("\n");

////////

//get sums of effects, centres and whether used
for(j=0;j<data_length;j++)
{
for(k=0;k<total;k++)
{
if(effects[k][j]!=0)
{effects[total][j]+=effects[k][j];effects[total+1][j]=blupcentres[k][j];effects[total+2][j]=1;}
}
}

//get sum of genetics
for(i=0;i<num_samples_use;i++)
{
for(k=0;k<total;k++){guesses[total][i]+=guesses[k][i];}
}

//get contribution of covariates - will always have
alpha=1.0;beta=0.0;
dgemv_("N", &num_samples_use, &num_fixed, &alpha, covar, &num_samples_use, thetas, &one, &beta, guesses[total+1], &one);

//add sum of genetic to covariates
for(i=0;i<num_samples_use;i++){guesses[total+2][i]=guesses[total][i]+guesses[total+1][i];}

//save
write_blups(effects, data_length, guesses, num_samples_use, num_kins, num_scores, num_tops, preds, al1, al2, ids1, ids2, outfile, resp, missingvalue);

//free allocations from setdl.c
for(k=0;k<total;k++){free(blupcentres[k]);free(blupfactors[k]);}free(blupcentres);free(blupfactors);
free(bluprands);

//frees from above
free(data);
for(k=0;k<total+3;k++){free(effects[k]);free(guesses[k]);}free(effects);free(guesses);
if(binary==0){gzclose(datainputgz);}

///////////////////////////

