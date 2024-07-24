/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Get per-predictor heritabilities based on indhers or by scaling data - also get savethetas

///////////////////////////

if(strcmp(indhers,"blank")!=0)	//already have exps - must have only one phenotype
{
for(j=0;j<data_length;j++){exps[j]=weights[j];}
}
else
{
printf("Computing per-predictor heritabilities\n\n");

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

//read data and standardize
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1, preds+bitstart);

if(num_fixed>1)	//get savethetas
{
if(dichot==0)	//only one set of thetas
{reg_covar_thetas(savethetas+bitstart*num_fixed, data, Z, num_samples_use, bitlength, num_fixed, NULL, 0);}
else	//can be multiple sets
{
for(m=0;m<num_resps_use;m++)
{reg_covar_thetas(savethetas+(bitstart+m*data_length)*num_fixed, data, Z, num_samples_use, bitlength, num_fixed, nullweights+m*num_samples_use, 0);}
}
}
}	//end of bit loop
printf("\n");
}

//get exps for first phenotype based on weights and power
count=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
count++;
}
else{exps[j]=0;}
}
if(count==0){printf("Error, all predictors are trivial\n\n");exit(1);}

//make sure exps sum to one
sum=0;for(j=0;j<data_length;j++){sum+=exps[j];}
for(j=0;j<data_length;j++){exps[j]=exps[j]/sum;}

for(m=1;m<num_resps_use;m++)	//copy exps into other phenotypes
{
for(j=0;j<data_length;j++){exps[j+m*data_length]=exps[j];}
}

///////////////////////////

