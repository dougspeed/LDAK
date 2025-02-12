/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Perform GRE - modes 191, 192, 193

///////////////////////////

if(mode==191)	//cut into partitions
{
if(chr[data_length-1]>chr[0])
{printf("Warning, the predictors span multiple chromosomes; usually GRE is run separately for each chromosome\n\n");}

if(part_length==-9999)	//set so that calc step (with bitsize=part_length) require about 15G memory
{
value=(double)data_length/1024*data_length/1024*8/1024;
if(value<16){value=16;}
part_length=15*1024/8*1024/num_samples_use*1024/2;
if(part_length>data_length){part_length=data_length;}
value=(double)part_length/1024*num_samples_use/1024*8/1024*2;
printf("The partition length is set to %d, so that the memory required when using \"--calc-gre\" will be approximately %.1f Gb; to change this value, use \"--partition-length\"\n\n", part_length, value);
}

//work out how many partitions and default size
num_parts=(data_length-1)/part_length+1;
count=(data_length-1)/num_parts+1;

//open file which stores partition details and write header lines, then partitions
sprintf(filename,"%spartition.details",folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Datafiles %s %s\n", datafile, bimfile);
if(extract==1){fprintf(output, "Using Filtered Predictors\n");}
else{fprintf(output, "Using All Predictors\n");}
if(strcmp(bsampfile,"blank")!=0||strcmp(csampfile,"blank")!=0){fprintf(output, "Using Filtered Samples\n");}
else{fprintf(output, "Using All Samples\n");}
fprintf(output,"Partition Start_Predictor End_Predictor\n");

for(q=0;q<num_parts;q++)
{
start=q*count;
end=(q+1)*count;
if(end>data_length){end=data_length;}
fprintf(output, "%d %d %d\n", count+1, start+1, end);
}
fclose(output);

//also save number of partitions
sprintf(filename2,"%spartition.number",folder);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check write permission granted and a folder of this name does not exist\n\n",filename2);exit(1);}
fprintf(output2,"%d\n", num_parts);
fclose(output2);

printf("The %d predictors have been split into %d partitions with (approx) length %d; details saved in %spartition.details\n\n", data_length, num_parts, part_length, folder);

}	//end of mode=191

////////

if(mode==192)	//compute correlations for specified partition
{
//bitsize is always the partition length
bitsize=pends[partition-1]-pstarts[partition-1];

//allocate variables
data_warn(2*bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);

anal_warn(bitsize, bitsize);
cors=malloc(sizeof(double)*bitsize*bitsize);

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress and on-the-fly files
sprintf(filename,"%sprogress.%d",folder,partition);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%scorrelations.%d",folder, partition);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

//ready for bit loop
bittotal=(data_length-pstarts[partition-1]-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=pstarts[partition-1]+bit*bitsize;
bitend=pstarts[partition-1]+(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
fclose(output2);
if((output2=fopen(filename2,"ab"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//read data into data2
current=read_data_fly(datafile, dtype, data2, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data2, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1);

if(bit==0)	//copy data2 into data - bitlength will equal bitsize
{
for(j=0;j<bitsize;j++)
{
for(i=0;i<num_samples_use;i++){data[(size_t)j*num_samples_use+i]=data2[(size_t)j*num_samples_use+i];}
}
}

//get correlations
alpha=1.0/num_samples_use;beta=0.0;
dgemm_("T", "N", &bitsize, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data2, &num_samples_use, &beta, cors, &bitsize);

//save them
for(j=0;j<bitlength;j++){fwrite(cors+(size_t)j*bitsize, sizeof(double), bitsize, output2);}
}	//end of bit loop
printf("\n");

fclose(output2);

free(data);free(data2);
free(cors);
}	//end of mode=192

////////

if(mode==193)	//join correlations
{
//see which partitions not present (or wrong size)
sprintf(filename,"%spartition.missings", folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

count3=0;
for(q=0;q<num_parts;q++)
{
count=pends[q]-pstarts[q];
count2=data_length-pstarts[q];
sprintf(filename2,"%scorrelations.%d", folder, q+1);
if(just_check(filename2)!=0)
{
if(count3<5){printf("Error, %s does not exist indicating that \"--calc-gre\" did not finish for Partition %d\n", filename2, q+1);}
fprintf(output, "%d\n", q+1);
count3++;
}
else
{
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}

fseeko(input2, 0, SEEK_END);
if(ftello(input2)!=(off_t)sizeof(double)*count*count2)
{printf("Error reading %s; should have size %jd (%d x %d predictors), but instead has size %jd\n\n", filename2, (off_t)sizeof(double)*count*count2, count, count2, ftello(input2));exit(1);}

fclose(input2);
}
}
if(count3==0){fprintf(output,"All partitions complete\n");}
fclose(output);

if(count3==1){printf("\nIn total, 1 partition did not complete; this is listed in %s\n\n", filename);exit(1);}
if(count3>1){printf("\nIn total, %d partitions did not complete; these are listed in %s\n\n", count3, filename);exit(1);}

////////

//allocate variables
anal_warn(data_length, data_length);
cors=malloc(sizeof(double)*data_length*data_length);
ipiv=malloc(sizeof(int)*data_length);

losts=malloc(sizeof(int)*data_length);

//read in correlations
for(q=0;q<num_parts;q++)
{
count=pends[q]-pstarts[q];
count2=data_length-pstarts[q];
printf("Reading correlations from Partition %d of %d\n", q+1, num_parts);

sprintf(filename2,"%scorrelations.%d", folder, q+1);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
fseeko(input2, 0, SEEK_SET);

for(j=0;j<count2;j++)
{
if(fread(cors+(size_t)(pstarts[q]+j)*data_length+pstarts[q], sizeof(double), count, input2)!=count)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-gre\"\n\n", j+1, filename2);exit(1);}
}
fclose(input2);
}
printf("\n");

//work out which predictors ignored (these will get self correlation 1, but correlation 0 with all others)
for(j=0;j<data_length;j++){losts[j]=0;}

//predictors with zero self correlation are trivial (will have correlation zero for all other elements)
count=0;
for(j=0;j<data_length;j++)
{
if(cors[(size_t)j*data_length+j]==0)
{
cors[(size_t)j*data_length+j]=1;
losts[j]=1;
count++;
}
}
if(count>0){printf("There are %d trivial predictors; these will be ignored\n\n", count);}

//check for near perfect correlations - simply move left-right, top-bottom (remember, we only use upper triangle)
printf("Searching for pairs of predictors that have squared correlation greater than %.6f\n", maxcor);

count2=0;
for(j2=0;j2<data_length;j2++)
{
for(j=0;j<j2;j++)
{
if(pow(cors[(size_t)j2*data_length+j],2)>maxcor)	//will remove j2
{
losts[j2]=2;
count2++;
//blank out upper triangle entries for j2 (their self correlation will remain 1)
for(j3=0;j3<j2;j3++){cors[(size_t)j2*data_length+j3]=0;}
for(j3=j2+1;j3<data_length;j3++){cors[(size_t)j3*data_length+j2]=0;}
break;
}
}
}

if(count2==0){printf("No predictors were removed\n\n");}
if(count2==1){printf("One predictor was removed\n\n");}
if(count2>1){printf("%d predictors were removed\n\n",count2);}

//perform decomposition
printf("Performing LDLT Decomposition\n");

lwork=-1;
dsytrf_("U", &data_length, cors, &data_length, ipiv, &wkopt, &lwork, &info);
if(info!=0){printf("Error, LDLT priming failed; please tell Doug (info %d, length %d)\n\n", info, data_length);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsytrf_("U", &data_length, cors, &data_length, ipiv, work, &lwork, &info);
if(info!=0){printf("Error, LDLT decomposition failed; please tell Doug (info %d, length %d)\n\n", info, data_length);exit(1);}
free(work);

//save decomposition and indexes of remaining predictors
printf("Saving the decomposition\n\n");

sprintf(filename3,"%sLDLT.bin",folder);
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(j=0;j<data_length;j++){fwrite(cors+(size_t)j*data_length, sizeof(double), j+1, output3);}
fwrite(ipiv, sizeof(int), data_length, output3);
fclose(output3);

sprintf(filename4,"%sLDLT.predictors",folder);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(j=0;j<data_length;j++)
{
if(losts[j]==0){fprintf(output4,"%s 0 Retained\n", preds[j]);}
if(losts[j]==1){fprintf(output4,"%s 1 Trivial\n", preds[j]);}
if(losts[j]==2){fprintf(output4,"%s 2 Duplicate\n", preds[j]);}
}
fclose(output4);

if(sinv==1)	//get and save inverse (ignoring lost predictors)
{
printf("Computing the inverse\n");

work=malloc(sizeof(double)*data_length);
dsytri_("U", &data_length, cors, &data_length, ipiv, work, &info);
if(info!=0){printf("Error, inverse failed; please tell Doug (info %d, length %d)\n\n", info, data_length);exit(1);}
free(work);

//save the upper triangle
printf("Saving the inverse\n");

sprintf(filename5,"%scors.inverse.bin",folder);
if((output5=fopen(filename5,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}

writefloats=malloc(sizeof(float)*data_length);
for(j2=0;j2<data_length;j2++)
{
if(losts[j2]==0)
{
count=0;
for(j=0;j<=j2;j++)
{
if(losts[j]==0){writefloats[count]=cors[(size_t)j2*data_length+j];count++;}
}
fwrite(writefloats, sizeof(float), count, output5);
}
}
fclose(output5);

sprintf(filename6,"%scors.inverse.predictors",folder);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
for(j=0;j<data_length;j++)
{
if(losts[j]==0){fprintf(output6, "%s\n", preds[j]);}
}
fclose(output6);

printf("The inverse is saved in files with prefix %scors\n\n", folder);
}

free(cors);free(ipiv);
free(losts);
}	//end of mode=193

////////

if(mode==194)	//solve
{
//check the response
count=0;for(i=0;i<num_samples_use;i++){count+=(resp[i]==missingvalue);}
if(count>0){printf("Error, phenotypes are missing for %d of the %d samples; you should either restart the analysis filtering out these samples, or use \"--dentist YES\" to pad missing values\n\n", count, num_samples_use);exit(1);}

//check the decomposition
sprintf(filename,"%sLDLT.bin",folder);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(double)*data_length*(data_length+1)/2+sizeof(int)*data_length)
{printf("Error reading %s; should have size %jd (%d x %d predictors), but instead has size %jd\n\n", filename, (off_t)sizeof(double)*data_length*(data_length+1)/2+sizeof(int)*data_length, data_length, data_length, ftello(input));exit(1);}
fclose(input);

if(bitsize==-9999)	//set bitsize
{
value=(double)data_length/1024*data_length/1024*8/1024;
if(value<1){value=1;}
bitsize=value*1024/8*1024/num_samples_use*1024;
if(bitsize>data_length){bitsize=data_length;}
printf("The bit-size will be set to %d (you can change this using \"--bit-size\")\n\n", bitsize);
}

//work out the total memory - the max of reading data and storing decomposition
value=(double)bitsize/1024*num_samples_use/1024*8/1024;
value2=(double)data_length/1024*data_length/1024*8/1024;
value3=value;
if(value2>value){value3=value2;}

if(value3>1)
{
if(value>value2+1)
{printf("Warning, to perform the analysis requires approximate %.1f Gb; you can reduce this to %.1f Gb using \"--bit-size\" (currently %d)\n\n", value, value2, bitsize);}
else
{printf("Warning, to perform the analysis requires approximately %.1f Gb\n\n", value2);}
}

////////

//allocate some variables (save data and decomposition for later)

losts=malloc(sizeof(int)*data_length);

Y=malloc(sizeof(double)*num_samples_use);
Yadj=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

YTdata=malloc(sizeof(double)*data_length);
YTdata2=malloc(sizeof(double)*data_length);

//read which predictors we will use
sprintf(filename,"%sLDLT.predictors",folder);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}

for(j=0;j<data_length;j++)
{
if(fscanf(input,"%s %d %s ", readstring, losts+j, readstring)!=3)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--join-gre\"\n\n", j+1, filename);exit(1);}
}
fclose(input);

//fill Y and Z, get thetas and standardized residuals (could be missing values)
for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[i];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}
reg_covar_lin_missing(Y, Z, num_samples_use, num_fixed, thetas, thetasds, thetapvas, Yadj, missingvalue);
stand_matrix_nomiss(Y, num_samples_use, num_samples_use, 1);

//save
sprintf(filename,"%s.coeff", greout);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Effect SE P\n");
fprintf(output, "Intercept %.4e %.4e %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d %.4e %.4e %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output);

////////

//now allocate data

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

data=malloc(sizeof(double)*num_samples_use*bitsize);

//deal with progress and on-the-fly files
sprintf(filename,"%s.progress",greout);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

fclose(output);
printf("Computing XTY for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Computing XTY for Chunk %d of %d\n", bit+1, bittotal);

//read data for chunk, and standardize
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1);

//compute YTdata for chunk
alpha=1.0/num_samples_use;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, Yadj, &one, &beta, YTdata+bitstart, &one);
}	//end of bit loop
printf("\n");

fclose(output);

//count number of predictors being used, and set YTdata to zero for remainder
count=0;
for(j=0;j<data_length;j++)
{
if(losts[j]==0){count++;}
else{YTdata[j]=0;}
}

//are now finished with data
free(data);
if(binary==0){gzclose(datainputgz);}

////////

//allocate and read decomposition
cors=malloc(sizeof(double)*data_length*data_length);
ipiv=malloc(sizeof(int)*data_length);

printf("Reading the LDLT decomposition\n\n");
sprintf(filename,"%sLDLT.bin",folder);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}

fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(double)*data_length*(data_length+1)/2+sizeof(int)*data_length)
{printf("Error reading %s; should have size %jd (%d x %d predictors), but instead has size %jd\n\n", filename, (off_t)sizeof(double)*data_length*(data_length+1)/2+sizeof(int)*data_length, data_length, data_length, ftello(input));exit(1);}

fseeko(input, 0, SEEK_SET);
for(j=0;j<data_length;j++)
{
if(fread(cors+(size_t)j*data_length, sizeof(double), j+1, input)!=j+1)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--join-gre\"\n\n", j+1, filename);exit(1);}
;}
if(fread(ipiv, sizeof(int), data_length, input)!=data_length)
{printf("Error reading the final row of %s, suggesting the file has been changed since creation with \"--join-gre\"\n\n", filename);exit(1);}
fclose(input);

//now compute inv cors YTdata
printf("Estimating joint effect sizes\n\n");

for(j=0;j<data_length;j++){YTdata2[j]=YTdata[j];}
dsytrs_("U", &data_length, &one, cors, &data_length, ipiv, YTdata2, &num_samples_use, &info);
if(info!=0){printf("Error, LDLT solve failed, please tell Doug (info %d, length %d)\n\n", info, num_samples_use);exit(1);}

//compute estimate = (n marginal beta * joint beta - m)/(n-m)
sum=0;for(j=0;j<data_length;j++){sum+=YTdata[j]*YTdata2[j];}
value=(num_samples_use*sum-count)/(num_samples_use-count);

//compute variance = complicated
value2=pow(num_samples_use-count,-2)*(2*count*(1-value)+4*value*num_samples_use)*(1-value);

sprintf(filename,"%s.estimate",greout);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Heritability %.6f\nSE %.6f\n", value, pow(value2,.5));
fclose(output);

printf("The estimated heritability is %.4f (SE %.4f), saved in %s\n\n", value, pow(value2,.5), filename);

//save scorefile
sprintf(filename2,"%s.score",greout);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Center Marginal Joint\n");
for(j=0;j<data_length;j++)
{
if(losts[j]==0){fprintf(output2, "%s %c %c %.6f %.6f %.6f\n", preds[j], al1[j], al2[j], centres[j], YTdata[j]*mults[j], YTdata2[j]*mults[j]);}
}
fclose(output2);

free(cors);free(ipiv);

free(losts);
free(Y);free(Yadj);free(Z);
free(thetas);free(thetasds);free(thetapvas);
free(YTdata);free(YTdata2);
}

///////////////////////////

