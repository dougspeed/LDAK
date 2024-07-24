/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Subtract predictors from a kinship

///////////////////////////

//allocate variables

data_warn2(bitsize,num_samples_use);

data=malloc(sizeof(double)*num_samples_use*bitsize);
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

anal_warn(num_samples_use, num_samples_use);
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);

//read details from kinstem[0]
sprintf(filename,"%s.grm.details",kinstems[0]);
count=countrows(filename)-1;
kpreds=malloc(sizeof(char*)*count);
kindex=malloc(sizeof(int)*count);
kcentres=malloc(sizeof(double)*count);
kmults=malloc(sizeof(double)*count);
kweights=malloc(sizeof(double)*count);
kal1=malloc(sizeof(char)*count);
kal2=malloc(sizeof(char)*count);
kexps=malloc(sizeof(double)*count);
read_details(filename, kpreds, kindex, kcentres, kmults, kweights, kal1, kal2, kexps, count);

//get data indexes, fill centres, mults and weights, blanking kmults when found
indexer=malloc(sizeof(int)*data_length);
count2=find_strings(preds, data_length, kpreds, count,  NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count2<data_length)
{printf("Error, please tell Doug; Code X27 %d and %d\n\n", count2, data_length);exit(1);}

for(j=0;j<data_length;j++)
{
j2=indexer[j];
if(al1[j]!=kal1[j2]&&al1[j]!=kal2[j2])
{printf("Error, alleles for %s (%c and %c) not consistent with those in %s (%c and %c)\n\n", preds[j], al1[j], al2[j], filename, kal1[j2], kal2[j2]);exit(1);}
if(al2[j]!=kal1[j2]&&al2[j]!=kal2[j2])
{printf("Error, alleles for %s (%c and %c) not consistent with those in %s (%c and %c)\n\n", preds[j], al1[j], al2[j], filename, kal1[j2], kal2[j2]);exit(1);}

if(al1[j]==kal1[j2]){centres[j]=kcentres[j2];}
else{centres[j]=2-kcentres[j2];}
mults[j]=kmults[j2];
weights[j]=kweights[j2];
kmults[j2]=-9999;
}	//end of j loop
free(indexer);

//ready to subtract - kins initialized to first kinship matrix
read_kins(kinstems[0], kins, NULL, kinsums[0], num_samples_use, ids3, 0, maxthreads);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Subtracting predictors for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Subtracting predictors for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, -9999, 0, -9999, weights+bitstart, 2, preds+bitstart);

alpha=-1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &bitlength, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, kins, &num_samples_use);
}	//end of bit loop
printf("\n");

////////

//get datafile and power
sprintf(filename,"%s.grm.adjust",kinstems[0]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")==0){flag=0;}
else{flag=1;}
if(fscanf(input, "%s %s ", readstring, readstring3)!=2)
{printf("Error reading Row 2 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Power")==0){power=atof(readstring3);}
else{power=-9999;}
fclose(input);

if(flag==0)
{write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 1, kpreds, kindex, kcentres, kmults, kweights, kal1, kal2, kexps, count, readstring2, power, kingz, kinraw, 1);}
else
{write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 1, kpreds, kindex, kcentres, kmults, kweights, kal1, kal2, kexps, count, NULL, power, kingz, kinraw, 1);}

free(data);
if(binary==0){gzclose(datainputgz);}
free(kins);
for(j=0;j<count;j++){free(kpreds[j]);}free(kpreds);
free(kindex);free(kcentres);free(kmults);free(kweights);free(kal1);free(kal2);free(kexps);

///////////////////////////

