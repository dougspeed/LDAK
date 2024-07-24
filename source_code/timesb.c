/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Do multiplications for calculating pca loadings

///////////////////////////

//sort out step
step=(50000/bitsize);if(step<20){step=20;}
if(step>20){step=10*(step/10);}if(step>50){step=20*(step/20);}
if(step>100){step=50*(step/50);}if(step>300){step=100*(step/100);}
if(step>1000){step=500*(step/500);}

//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

blupmids=malloc(sizeof(double)*bitsize*axes);
blupprojs=malloc(sizeof(double)*num_samples_use*axes);
effects=malloc(sizeof(double*)*axes);
guesses=malloc(sizeof(double*)*axes);
for(k=0;k<axes;k++)
{effects[k]=malloc(sizeof(double)*data_length);guesses[k]=malloc(sizeof(double)*num_samples_use);}

//set effects and guesses to zero (redundant for some columns of effects)
for(k=0;k<axes;k++)
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
printf("Calculating loadings for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating loadings for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);

//standardize based on blupcentres[0] and blupfactors[0]
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{data[(size_t)j*num_samples_use+i]=(data[(size_t)j*num_samples_use+i]-blupcentres[0][bitstart+j])*blupfactors[0][bitstart+j];}
else
{data[(size_t)j*num_samples_use+i]=0;}
}
}

//want scaled effect sizes XT bluprands / kinsums / eigens; also want X times these

//get XT bluprands
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &axes, &num_samples_use, &alpha, data, &num_samples_use, bluprands, &num_samples_use, &beta, blupmids, &bitlength);

//divide these by kinsums and eigens
for(k=0;k<axes;k++)
{
for(j=0;j<bitlength;j++){blupmids[j+k*bitlength]=blupmids[j+k*bitlength]/kinsums[0]/blupvalues[k];}
}

//load in raw effect sizes
for(k=0;k<axes;k++)
{
for(j=0;j<bitlength;j++){effects[k][bitstart+j]=blupmids[j+k*bitlength]*blupfactors[0][bitstart+j];}
}

//get contribution to projections
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_samples_use, &axes, &bitlength, &alpha, data, &num_samples_use, blupmids, &bitlength, &beta, blupprojs, &num_samples_use);

//add these to guesses
for(k=0;k<axes;k++)
{
for(i=0;i<num_samples_use;i++){guesses[k][i]+=blupprojs[i+k*num_samples_use];}
}
}	//end of bit loop
printf("\n");

////////

//save
write_loads(effects, data_length, guesses, num_samples_use, axes, blupcentres, preds, al1, al2, ids1, ids2, outfile);

//free allocations from setdl.c
free(blupcentres[0]);free(blupcentres);
free(blupfactors[0]);free(blupfactors);
free(bluprands);free(blupvalues);

//frees from above
free(data);
free(blupmids);free(blupprojs);
for(k=0;k<axes;k++){free(effects[k]);free(guesses[k]);}free(effects);free(guesses);
if(binary==0){gzclose(datainputgz);}

///////////////////////////

