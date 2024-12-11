/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Computing kinships - have a predictor-predictor inverse

///////////////////////////

//work out how much memory we use - try to be clever with allocations

data_warn(data_length,num_samples_use);
if(data_length>num_samples_use){anal_warn(data_length, data_length/2);}
else{anal_warn(num_samples_use, num_samples_use/2);}

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//read the data and fill data_single
printf("Reading the data\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Reading the data\n");
fclose(output);

data=malloc(sizeof(double)*num_samples_use*data_length);
data_single=malloc(sizeof(float)*num_samples_use*data_length);

if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}

(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, 0, data_length, keeppreds_use, datainputgz, 0, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres, mults, sqdevs, rates, infos, num_samples_use, data_length, missingvalue, power, strcmp(centresfile,"blank")!=0, hwestand, weights, 1);

count=0;for(j=0;j<data_length;j++){count+=(mults[j]!=-9999);}
if(count==0){printf("Error, all %d predictors trivial\n\n", data_length);exit(1);}

//load data into data_single
for(j=0;j<data_length;j++)
{
for(i=0;i<num_samples_use;i++)
{data_single[(size_t)j*num_samples_use+i]=(float)data[(size_t)j*num_samples_use+i];
}
}

/*
test=fopen("data.txt","w");
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<data_length;j++){fprintf(test, "%f ", data_single[j*num_samples+i]);}
fprintf(test, "\n");
}
fclose(test);
*/

//can now free data
free(data);
if(binary==0){gzclose(datainputgz);}

//read the inverse - have checked order and size is correct
printf("Reading the inverse\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Reading the inverse\n");
fclose(output);

inv_single=malloc(sizeof(float)*data_length*data_length);
for(j=0;j<data_length;j++)
{
for(j2=0;j2<data_length;j2++){inv_single[j2*data_length+j]=0;}
}

sprintf(filename2,"%s.inverse.bin",invsfile);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);
for(j2=0;j2<data_length;j2++)
{
if(fread(inv_single+(size_t)j2*data_length, sizeof(float), j2+1, input2)!=j2+1)
{printf("Error reading Column %d of %s\n\n", j2+1, filename2);exit(1);}
}
fclose(input2);

/*
FILE *test;
test=fopen("inv.txt","w");
for(j=0;j<data_length;j++)
{
for(j2=0;j2<data_length;j2++){fprintf(test, "%f ", inv_single[j2*data_length+j]);}
fprintf(test,"\n");
}
fclose(test);
*/

//get data_single * inv_single
printf("Multiplying the data by the inverse\n\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Multiplying the data by the inverse\n");
fclose(output);

data_single2=malloc(sizeof(float)*num_samples_use*data_length);

alpha_single=1.0;beta_single=0.0;
ssymm_("R", "U", &num_samples_use, &data_length, &alpha_single, inv_single, &data_length, data_single, &num_samples_use, &beta_single, data_single2, &num_samples_use);

//can now free inv_single
free(inv_single);

/*
test=fopen("datainv.txt","w");
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<data_length;j++){fprintf(test, "%f ", data_single2[j*num_samples+i]);}
fprintf(test, "\n");
}
fclose(test);
*/

//get data_single2 x data_single
printf("Computing the kinships\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Computing the kinships\n");
fclose(output);

kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);

alpha_single=1.0;beta_single=0.0;
sgemm_("N", "T", &num_samples_use, &num_samples_use, &data_length, &alpha_single, data_single2, &num_samples_use, data_single, &num_samples_use, &beta_single, kins_single, &num_samples_use);

//save kins
printf("Saving the kinships\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Saving the kinships\n");
fclose(output);

write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, data_length, datafile, power, kingz, kinraw, 1);

free(data_single);free(data_single2);
free(kins_single);

///////////////////////////

