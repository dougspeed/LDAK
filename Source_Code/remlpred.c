/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//BLUP predictions using only kinship matrices

///////////////////////////

//get numbers of blup and predictor individuals
numa=countrows(blupfile);
numb=num_samples_use-numa;
printf("Computing predictions for %d samples\n\n", numb);

indexer=malloc(sizeof(int)*numa);
indexer2=malloc(sizeof(int)*numb);

//get details of blup individuals
wantids=malloc(sizeof(char*)*numa);
read_ids(blupfile, NULL, NULL, wantids, numa, NULL, 0, 0);

//get details of predictor individuals
kinids=malloc(sizeof(char*)*numb);
count=find_strings(ids3, num_samples_use, wantids, numa, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count!=numa){printf("Error finding indexes AA, please tell Doug\n\n");exit(1);}

usedids=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){usedids[i]=0;}
for(i=0;i<numa;i++){usedids[indexer[i]]=1;}

count=0;
for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==0){copy_string(kinids,count,ids3[i]);count++;}
}
free(usedids);

if(bitsize>numb){bitsize=numb;}

//allocate variables
anal_warn(numa, numb/2);
kins_single=malloc(sizeof(float)*numa*numb);

bluprands=malloc(sizeof(double)*numa*num_kins);
guesses=malloc(sizeof(double*)*(num_kins+3));	//last three are sum of genetics, covs and sum of both
for(k=0;k<num_kins+3;k++){guesses[k]=malloc(sizeof(double)*numb);}
Z=malloc(sizeof(double)*numb*num_covars);

bluprand_single=malloc(sizeof(float)*numa);
guess_single=malloc(sizeof(float)*numb);

Mcurrent=malloc(sizeof(int)*maxthreads);
Minput=malloc(sizeof(FILE *)*maxthreads);
Mdatatemp=malloc(sizeof(float*)*maxthreads);

//read random effects
read_rand(blupfile, bluprands, num_kins, numa, wantids);

//set guesses to zero
for(k=0;k<num_kins+3;k++)
{
for(i=0;i<numb;i++){guesses[k][i]=0;}
}

//load covariates (for predictor samples)
count=find_strings(ids3, num_samples_use, kinids, numb, indexer2, NULL, NULL, NULL, NULL, NULL, 3);
if(count!=numb){printf("Error finding indexes BB, please tell Doug\n\n");exit(1);}
for(i=0;i<numb;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*numb]=covar[indexer2[i]+j*num_samples_use];}
}

////////

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

for(k=0;k<num_kins;k++)
{
printf("Calculating predictions for Kinship Matrix %d of %d\n", k+1, num_kins);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating predictions for Kinship Matrix %d of %d\n", k+1, num_kins);
fclose(output);

//get indexes of individuals we want
sprintf(filename2, "%s.grm.id", kinstems[k]);
count=countrows(filename2);
wantids2=malloc(sizeof(char*)*count);
read_ids(filename2, NULL, NULL, wantids2, count, NULL, 0, 0);

count2=find_strings(wantids, numa, wantids2, count, NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count2!=numa){printf("Error finding indexes CC, please tell Doug\n\n");exit(1);}
count2=find_strings(kinids, numb, wantids2, count, NULL, indexer2, NULL, NULL, NULL, NULL, 3);
if(count2!=numb){printf("Error finding indexes DD, please tell Doug\n\n");exit(1);}

//usedids indicates whether not reading (0), part of numa (>0) or part of numb (<0)
usedids=malloc(sizeof(int)*count);
for(i=0;i<count;i++){usedids[i]=0;}
for(i=0;i<numa;i++){usedids[indexer[i]]=1+i;}
for(i=0;i<numb;i++){usedids[indexer2[i]]=-1-i;}

//open kinship and check size
sprintf(filename2, "%s.grm.bin", kinstems[k]);
if((Minput[0]=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n", filename2);exit(1);}
fseeko(Minput[0], 0, SEEK_END);
if(ftello(Minput[0])!=(off_t)sizeof(float)*count*(count+1)/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename2, (off_t)sizeof(float)*count*(count+1)/2, ftello(Minput[0]));exit(1);}

#pragma omp parallel for private(thread,threadstart,threadend,i,i2) schedule (static,1)
for(thread=0;thread<maxthreads;thread++)
{
threadstart=pow((double)thread/maxthreads,.5)*count;
threadend=pow((double)(thread+1)/maxthreads,.5)*count;

if((Minput[thread]=fopen(filename2,"rb"))==NULL)
{printf("Error opening %s\n\n", filename2);exit(1);}

fseeko(Minput[thread], 0, SEEK_SET);
Mcurrent[thread]=0;

//read one "row" at a time
Mdatatemp[thread]=malloc(sizeof(float)*count);
for(i=threadstart;i<threadend;i++)
{
if(count>40000&&thread==0&&i%20000==0&&i*maxthreads<count)
{printf("Reading kinships for Sample %d out of %d\n", i*maxthreads+1, count);}

if(usedids[i]>0)	//reading and part of numa
{
if(i!=Mcurrent[thread])	//get to start of Row i
{fseeko(Minput[thread], (off_t)sizeof(float)*i*(i+1)/2, SEEK_SET);}
if(fread(Mdatatemp[thread], sizeof(float), i+1, Minput[thread])!=i+1)
{printf("Error reading Row %d of %s\n\n", i+1, filename2);exit(1);};
Mcurrent[thread]=i+1;

for(i2=0;i2<i;i2++)
{
if(usedids[i2]<0)	//part of numb
{kins_single[(size_t)(-usedids[i2]-1)*numa+(usedids[i]-1)]=Mdatatemp[thread][i2];}
}
}

if(usedids[i]<0)	//reading and part of numb
{
if(i!=Mcurrent[thread])	//get to start of Row i
{fseeko(Minput[thread], (off_t)sizeof(float)*i*(i+1)/2, SEEK_SET);}
if(fread(Mdatatemp[thread], sizeof(float), i+1, Minput[thread])!=i+1)
{printf("Error reading Row %d of %s\n\n", i+1, filename2);exit(1);};
Mcurrent[thread]=i+1;

for(i2=0;i2<i;i2++)
{
if(usedids[i2]>0)	//part of numa
{kins_single[(size_t)(-usedids[i]-1)*numa+(usedids[i2]-1)]=Mdatatemp[thread][i2];}
}
}
}	//end of i loop

free(Mdatatemp[thread]);
}	//end of thread loop
printf("\n");

for(i=0;i<count;i++){free(wantids2[i]);}free(wantids2);
free(usedids);

//now do multiplication
alpha_single=1.0;beta_single=0.0;
for(i=0;i<numa;i++){bluprand_single[i]=bluprands[i+k*numa];}
sgemv_("T", &numa, &numb, &alpha_single, kins_single, &numa, bluprand_single, &one, &beta_single, guess_single, &one);
for(i=0;i<numb;i++){guesses[k][i]=guess_single[i];}
}	//end of k loop

//get sum of genetics
for(i=0;i<numb;i++)
{
for(k=0;k<num_kins;k++){guesses[num_kins][i]+=guesses[k][i];}
}

//get contribution of covariates
alpha=1.0;beta=0.0;
dgemv_("N", &numb, &num_fixed, &alpha, Z, &numb, thetas, &one, &beta, guesses[num_kins+1], &one);

//add sum of genetic to covariates
for(i=0;i<numb;i++){guesses[num_kins+2][i]=guesses[num_kins][i]+guesses[num_kins+1][i];}

//save
count=find_strings(kinids, numb, ids3, num_samples_use, NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count!=numb){printf("Error finding indexes EE, please tell Doug\n\n");exit(1);}

sprintf(filename3,"%s.pred",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"ID1\tID2\tPhenotype\tCovariates\tGenetics\tTotal\n");
for(i=0;i<numb;i++)
{
fprintf(output3, "%s\t%s\t", ids1[indexer[i]], ids2[indexer[i]]);
if(resp[indexer[i]]!=missingvalue){fprintf(output3, "%.6f\t", resp[indexer[i]]);}
else{fprintf(output3, "NA\t");}
fprintf(output3, "%.6f\t%.6f\t%.6f\n", guesses[num_kins+1][i], guesses[num_kins][i], guesses[num_kins+2][i]);
}
fclose(output3);

sprintf(filename4,"%s.pred.full",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"ID1\tID2");
for(k=0;k<num_kins;k++){fprintf(output4,"\tKinship_%d", k+1);}
fprintf(output4,"\n");
for(i=0;i<numb;i++)
{
fprintf(output4, "%s\t%s\t", ids1[indexer[i]], ids2[indexer[i]]);
for(k=0;k<num_kins;k++){fprintf(output4, "%.6f\t", guesses[k][i]);}
fprintf(output4, "\n");
}
fclose(output4);

printf("Predictions in %s (and %s)\n\n", filename3, filename4);

free(indexer);free(indexer2);
for(i=0;i<numa;i++){free(wantids[i]);}free(wantids);
for(i=0;i<numb;i++){free(kinids[i]);}free(kinids);
free(kins);
free(bluprands);
for(k=0;k<num_kins+3;k++){free(guesses[k]);}free(guesses);
free(Z);
free(bluprand_single);free(guess_single);
free(Mcurrent);free(Minput);free(Mdatatemp);

///////////////////////////

