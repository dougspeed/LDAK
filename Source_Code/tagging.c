/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Calculating tagging (and checks and coverage) and adjusting weights
//if using pathways (genpreds!=NULL), will have reduce=0, savemat=0 and coverage=0

///////////////////////////

//get labels
catlabels=malloc(sizeof(char *)*(num_parts+1));
for(q=0;q<num_parts+1;q++){catlabels[q]=malloc(sizeof(char)*500);}

if(num_parts>0)
{
if(strcmp(labfile,"blank")!=0)	//read from file
{
if((input=fopen(labfile,"r"))==NULL)
{printf("Error opening %s\n\n",labfile);exit(1);}
for(q=0;q<num_parts;q++)
{
catlabels[q]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", catlabels[q])!=1)
{printf("Error reading Row %d of %s\n\n", q+1, labfile);exit(1);}

if(strcmp(catlabels[q],"Base")==0||strcmp(catlabels[q],"Background")==0)
{printf("Error reading Row %d of %s; labels can not be either \"Base\" or \"Background\"\n\n", q+1, labfile);exit(1);}

for(q2=0;q2<q;q2++)	//check different to previous
{
if(strcmp(catlabels[q],catlabels[q2])==0)
{printf("Error reading %s; the label %s appears twice\n\n", labfile, catlabels[q]);exit(1);}
}
}
fclose(input);
}
else	//give generic names
{
if(parttype==0)
{
if(strcmp(genpreds,"blank")==0)
{
for(q=0;q<num_parts;q++){sprintf(catlabels[q],"Annotation_%d",q+1);}
}
else
{
for(q=0;q<num_parts;q++){sprintf(catlabels[q],"Pathway_%d",q+1);}
}
}
else
{
for(q=0;q<num_parts;q++){sprintf(catlabels[q],"Partition_%d",q+1);}
}
}
}

if(parttype==0)	//always have a base
{strcpy(catlabels[num_parts],"Base");}
else	//add background, even if not used
{strcpy(catlabels[num_parts],"Background");}

//addgenic indicates whether using genic snps
if(strcmp(genpreds,"blank")==0){addgenic=0;}
else{addgenic=1;}

//get the heritability model
model_warn(data_length*3/2, num_parts+1);
pindexes=malloc(sizeof(int *)*(num_parts+2));
pweights=malloc(sizeof(double *)*(num_parts+2));
usedpreds=malloc(sizeof(int)*data_length);

addpart=get_her_model(num_parts, partpref, pindexes, pweights, data_length, keeppreds_use, num_preds, allpreds, predorder, parttype, backpart, allone);

//get_her_model only allocates the first num_parts+1 columns of pindexes and pweights
pindexes[num_parts+1]=malloc(sizeof(int)*data_length);
pweights[num_parts+1]=malloc(sizeof(double)*data_length);

if(addpart==2)	//there are some predictors not in a partition
{
if(mode==145){printf("Warning, it is unusual for there to be predictors that are not in a partition\n\n");}

if(reduce==1||mode==145)	//squeeze down
{
count=0;
for(j=0;j<data_length;j++)
{
flag=0;for(q=0;q<num_parts;q++){flag+=pindexes[q][j];}
if(flag>0)
{
if(count!=j)
{
chr[count]=chr[j];cmbp[count]=cmbp[j];al1[count]=al1[j];al2[count]=al2[j];
free(preds[count]);copy_string(preds,count,preds[j]);
free(along1[count]);copy_string(along1,count,along1[j]);free(along2[count]);copy_string(along2,count,along2[j]);
keeppreds_use[count]=keeppreds_use[j];
centres[count]=centres[j];mults[count]=mults[j];sqdevs[count]=sqdevs[j];rates[count]=rates[j];infos[count]=infos[j];
weights[count]=weights[j];pvalues[count]=pvalues[j];
for(q=0;q<num_parts;q++){pindexes[q][count]=pindexes[q][j];pweights[q][count]=pweights[q][j];}
}
count++;
}}
for(j=count;j<data_length;j++){free(preds[j]);free(along1[j]);free(along2[j]);}
data_length=count;
}

//reset addpart
addpart=0;
}

if(strcmp(genpreds,"blank")!=0)	//get genic predictors
{
count=countrows(genpreds);
printf("Reading %d genic predictors from %s\n", count, genpreds);
wantpreds=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*count);

read_strings(genpreds, wantpreds, count, NULL, 1, 0);
count2=find_strings(preds, data_length, wantpreds, count, indexer, NULL, NULL, genpreds, NULL, NULL, 3);
if(count2==0){printf("Error, none of these are in the data\n\n");exit(1);}
if(count2<count){printf("Warning, only %d of these are in the data\n\n", count2);}
else{printf("All of these are in the data\n\n");}

for(j=0;j<data_length;j++){pindexes[num_parts+addpart][j]=0;pweights[num_parts+addpart][j]=0;}
for(j=0;j<count2;j++){pindexes[num_parts+addpart][indexer[j]]=1;pweights[num_parts+addpart][indexer[j]]=1;}

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);
free(indexer);

//see if predictors consistent with annotations
count2=0;
for(q=0;q<num_parts;q++)
{
count=0;
for(j=0;j<data_length;j++){count+=(pindexes[q][j]==1&&pindexes[num_parts+addpart][j]==0);}
if(count>0)
{
if(count2<10)
{printf("Warning, %d of the predictors in Annotation %d are not in the list of genic predictors\n", count, q+1);}
count2++;
}
}
if(count2>0){printf("\n");}
}

//record how many categories each predictor in (do not include genic)
for(j=0;j<data_length;j++)
{
usedpreds[j]=0;for(q=0;q<num_parts+addpart;q++){usedpreds[j]+=pindexes[q][j];}
}

if(mode==141)	//get lists of predictors to print and to compute ssums over
{
windex=malloc(sizeof(int)*data_length);
vindex=malloc(sizeof(int)*data_length);

if(strcmp(printfile,"blank")!=0)	//print those in printfile
{
count=countrows(printfile);
printf("Reading list of %d regression predictors from %s\n", count, printfile);
wantpreds=malloc(sizeof(char*)*count);
read_strings(printfile, wantpreds, count, NULL, 1, 0);
indexer=malloc(sizeof(int)*count);
count2=find_strings(preds, data_length, wantpreds, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, none of these are in the heritability model\n\n");exit(1);}
if(count2<count){printf("Warning, only %d of these are in the heritability model\n", count2);}
printf("\n");

for(j=0;j<data_length;j++){windex[j]=0;}
for(j=0;j<count2;j++){windex[indexer[j]]=1;}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}
else	//will print all
{
for(j=0;j<data_length;j++){windex[j]=1;}
}

if(strcmp(herfile,"blank")!=0)	//when computing heritabilities, use those in herfile
{
count=countrows(herfile);
printf("Reading list of %d heritability predictors from %s\n", count, herfile);
wantpreds=malloc(sizeof(char*)*count);
read_strings(herfile, wantpreds, count, NULL, 1, 0);
indexer=malloc(sizeof(int)*count);
count2=find_strings(preds, data_length, wantpreds, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, none of these are in the heritability model\n\n");exit(1);}
if(count2<count){printf("Warning, only %d of these are in the heritability model\n", count2);}
printf("\n");

for(j=0;j<data_length;j++){vindex[j]=0;}
for(j=0;j<count2;j++){vindex[indexer[j]]=1;}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}
else	//will use all
{
for(j=0;j<data_length;j++){vindex[j]=1;}
}
}

////////

if(mode==105)	//print out new weights and partitions
{
for(q=0;q<num_parts+addpart;q++)
{
count=0;for(j=0;j<data_length;j++){count+=(usedpreds[j]==q+1);}
if(count>0)
{
if(q==0){printf("%d predictors are in exactly 1 partition\n", count);}
else{printf("%d predictors are in exactly %d partitions\n", count, q+1);}
}
}
printf("\n");

if(allone==1)
{
for(j=0;j<data_length;j++){usedpreds[j]=(usedpreds[j]>0);}
}

sprintf(filename,"%s.weights",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(j=0;j<data_length;j++){fprintf(output,"%s %.6f\n", preds[j], weights[j]/usedpreds[j]);}
fclose(output);

for(q=0;q<num_parts+addpart;q++)
{
sprintf(filename2,"%s%d",outfile, q+1);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
for(j=0;j<data_length;j++)
{
if(pweights[q][j]>0){fprintf(output2,"%s\n", preds[j]);}
}
fclose(output2);
}

printf("Adjusted weights saved in %s, with partition predictors in %s1", filename, outfile);
if(num_parts+addpart==2){printf(" and %s2", outfile);}
if(num_parts+addpart>2){printf(" ... %s%d", outfile, num_parts+addpart);}
printf("\n\n");
}	//end of mode=105

////////

if(mode==141)	//get tagging
{
if(bitsize==-9999)	//optimal bitsize is average number of predictors in each window (always using window_kb)
{
scount=0;
k=1;
for(j=0;j<data_length;j++)
{
while(1)	//move along until out of range, tallying each data_length predictors encountered
{
if(k==data_length){break;}
if(chr[k]!=chr[j]||cmbp[k]-cmbp[j]>1000*window_kb){break;}
k++;
}
scount+=k-j;
}
bitsize=(int)((scount-1)/data_length)+1;

if(bitsize<20){bitsize=20;}
if(bitsize>8000){bitsize=8000;}
printf("The bit-size will be set to %d (you can change this using \"--bit-size\")\n\n", bitsize);
}

step=(50000/bitsize);if(step<20){step=20;}
if(step>20){step=10*(step/10);}if(step>50){step=20*(step/20);}
if(step>100){step=50*(step/50);}if(step>300){step=100*(step/100);}
if(step>1000){step=500*(step/500);}

//work out bitmax
bitmax=bitsize;
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
if(bitend2-bitstart>bitmax){bitmax=bitend2-bitstart;}
}

//allocate variables
data_warn3(bitmax,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitmax);
exps=malloc(sizeof(double)*data_length);

anal_warn(bitmax, bitmax);
cors=malloc(sizeof(double)*bitmax*bitmax);

//tallies store number of tagging predictors, total tagging, relative heritability tagged (and coverage)
tally1=malloc(sizeof(double)*data_length);
tally2=malloc(sizeof(double)*data_length);
tally3=malloc(sizeof(double)*data_length*(num_parts+addpart+addgenic));
if(cover==1){tally4=malloc(sizeof(double)*data_length*(num_parts+addpart+addgenic));}
for(j=0;j<data_length;j++)
{
tally1[j]=0;
tally2[j]=0;
for(q=0;q<num_parts+addpart+addgenic;q++){tally3[j+q*data_length]=0;}
if(cover==1)
{
for(q=0;q<num_parts+addpart+addgenic;q++){tally4[j+q*data_length]=0;}
}
}

if(strcmp(genpreds,"blank")==0)
{
ssums=malloc(sizeof(double*)*(num_parts+addpart));
for(q=0;q<num_parts+addpart;q++){ssums[q]=malloc(sizeof(double)*(num_parts+addpart+1));}
}
else	//will use a reduced version of ssums
{
ssums=malloc(sizeof(double*)*3);
ssums[0]=malloc(sizeof(double)*4);
ssums[1]=malloc(sizeof(double)*4);
ssums[2]=malloc(sizeof(double)*4);
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

//deal with progress and on-the-fly files
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

if(strcmp(genpreds,"blank")==0){sprintf(filename2,"%s.tagging", outfile);}
else{sprintf(filename2,"%s.base.tagging", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1 A2 Neighbours Tagging ");
if(hwestand==1){fprintf(output2,"Weight MAF Categories Exp_Heritability");}
else{fprintf(output2,"Weight Variance Categories Exp_Heritability");}

if(strcmp(genpreds,"blank")==0)
{
for(q=0;q<num_parts+addpart;q++){fprintf(output2," %s", catlabels[q]);}
fprintf(output2,"\n");
}
else{fprintf(output2," Base Genic\n");}

if(savemat==1)
{
sprintf(filename3,"%s.matrix", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor");
for(q=0;q<num_parts+addpart;q++){fprintf(output3," %s", catlabels[q]);}
fprintf(output3,"\n");
}

if(cover==1)
{
sprintf(filename4,"%s.coverage", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Predictor");
for(q=0;q<num_parts;q++){fprintf(output4," Partition_%d", q+1);}
if(addpart==1){fprintf(output4," Background");}
fprintf(output4,"\n");
}

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
bitlength=bitend2-bitstart;

if(bit%step==0)
{
printf("Calculating tagging for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Calculating tagging for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
if(savemat==1)
{
fclose(output3);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
}
if(cover==1)
{
fclose(output4);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
}
}

shuffle=0;
for(j=0;j<end-bitstart;j++)	//using values already in data, so shuffle back
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(bitstart-start+j)*num_samples_use+i];}
shuffle++;
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, bitstart+shuffle, bitend2, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data+(size_t)shuffle*num_samples_use, centres+bitstart+shuffle, mults+bitstart+shuffle, sqdevs+bitstart+shuffle, rates+bitstart+shuffle, infos+bitstart+shuffle, num_samples_use, bitlength-shuffle, missingvalue, -1, 0, 0, NULL, 1);

for(j=bitstart+shuffle;j<bitend2;j++)	//want exps to contain expected heritability
{
if(mults[j]!=-9999)	//non trivial
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}
else{exps[j]=0;}
}

//get correlation
alpha=1.0/num_samples_use;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors, &bitlength);

//now loop through predictors in chunk, getting required sums and printing
for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)	//non trivial
{
//diagonals
tally1[j]++;
if(windex[j]==1){tally2[j]++;}
for(q=0;q<num_parts+addpart+addgenic;q++){tally3[j+q*data_length]+=pweights[q][j]*exps[j];}

if(cover==1)
{
for(q=0;q<num_parts+addpart;q++){tally4[j+q*data_length]=pindexes[q][j];}
}

//off-diagonals
for(k=j+1;k<bitend2;k++)
{
if(chr[k]!=chr[j]){break;}
if(cmbp[k]-cmbp[j]>1000*window_kb){break;}

if(mults[k]!=-9999)	//non-trivial and within range
{
value=pow(cors[(size_t)(j-bitstart)*bitlength+(k-bitstart)],2);
if(unbias==1){value2=value-(1.0-value)/(num_samples_use-2);}
else{value2=value;}

tally1[j]++;
tally1[k]++;
if(windex[j]==1&&windex[k]==1)	//both are regression SNPs
{tally2[j]+=value2;tally2[k]+=value2;}
if(value>=mincor)
{
for(q=0;q<num_parts+addpart+addgenic;q++)
{tally3[j+q*data_length]+=pweights[q][k]*value2*exps[k];
tally3[k+q*data_length]+=pweights[q][j]*value2*exps[j];}
}

if(cover==1)
{
for(q=0;q<num_parts+addpart;q++)
{
if(pindexes[q][k]==1&&value>tally4[j+q*data_length]){tally4[j+q*data_length]=value;}
}
}
}}	//end of using k and k loop

if(windex[j]==1)	//regression snp, so print tagging
{
fprintf(output2, "%s %c %c %d %.3f ", preds[j], al1[j], al2[j], (int)tally1[j], tally2[j]);
if(strcmp(weightsfile,"blank")!=0){fprintf(output2, "%.4f ", weights[j]);}
else{fprintf(output2, "1 ");}
if(hwestand==1)
{fprintf(output2, "%.6f %d %.4f", centres[j]/2+(centres[j]>1)*(1-centres[j]), usedpreds[j], exps[j]);}
else
{fprintf(output2, "%.6f %d %.4f", pow(sqdevs[j],2), usedpreds[j], exps[j]);}

if(strcmp(genpreds,"blank")==0)
{
for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %.4f", tally3[j+q*data_length]);}
fprintf(output2, "\n");
}
else{fprintf(output2, " %.4f %.4f\n", tally3[j+num_parts*data_length], tally3[j+(num_parts+1)*data_length]);}
}

if(savemat==1&&vindex[j]==1)
{
fprintf(output3, "%s ", preds[j]);
for(q=0;q<num_parts+addpart;q++){fprintf(output3, " %.4f", pweights[q][j]*exps[j]);}
fprintf(output3, "\n");
}
if(cover==1)
{
fprintf(output4, "%s ", preds[j]);
for(q=0;q<num_parts+addpart;q++){fprintf(output4, " %.4f", tally4[j+q*data_length]);}
fprintf(output4, "\n");
}
}}	//end of using j and j loop

start=bitstart;if(bitend2>end){end=bitend2;}
}	//end of bit loop
printf("\n");

if(strcmp(genpreds,"blank")!=0)	//print pathway taggings as a separate files, then print names
{
sprintf(filename5,"%s.pathway.tagging", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}

for(q=0;q<num_parts;q++)
{
fprintf(output5,"%s", catlabels[q]);
for(j=0;j<data_length;j++)
{
if(windex[j]==1)
{
if(fourdp==0){fprintf(output5," %.1f", tally3[j+q*data_length]);}
else{fprintf(output5," %.4f", tally3[j+q*data_length]);}
}
}
fprintf(output5,"\n");
}
fclose(output5);

//printf("Compressing the file %s\n\n", filename5);
//sprintf(cmd, "gzip -f %s", filename5);
//system(cmd);

/*
for(q=0;q<num_parts;q++)
{
sprintf(filename5,"%s.pathway.%d", outfile, q+1);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(j=0;j<data_length;j++)
{
if(windex[j]==1){fprintf(output5,"%.4f\n", tally3[j+q*data_length]);}
}
fclose(output5);
sprintf(cmd, "gzip -f %s", filename5);
system(cmd);

sprintf(filename5,"%s.bathway.%d", outfile, q+1);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(j=0;j<data_length;j++)
{
if(windex[j]==1){fprintf(output5,"%.1f\n", tally3[j+q*data_length]);}
}
fclose(output5);
sprintf(cmd, "gzip -f %s", filename5);
system(cmd);
}

sprintf(filename5,"%s.pathway.names", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(q=0;q<num_parts;q++){fprintf(output5," %s\n", catlabels[q]);}
fclose(output5);
*/
}

if(strcmp(genpreds,"blank")==0)	//compute full ssums (using only heritability snps) and add to bottom of tagging file
{
for(q=0;q<num_parts+addpart;q++)
{
//ssums[q][q2]indicates how much q2 contributes to q
for(q2=0;q2<num_parts+addpart;q2++){ssums[q][q2]=0;}
//ssums[q][num_parts+addpart] is the number of snps in category q (previously was sum of its snps heritabilities)
ssums[q][num_parts+addpart]=0;

for(j=0;j<data_length;j++)	//will include trivial snps (but these will have exp zero)
{
if(pindexes[q][j]==1&&vindex[j]==1)	//predictor is part of q and a heritability predictor
{
//add on contribution from q2
for(q2=0;q2<num_parts+addpart;q2++){ssums[q][q2]+=pweights[q2][j]*exps[j];}
//increase tally (previously added on expected heritability)
ssums[q][num_parts+addpart]++;
}
}
}

for(q2=0;q2<num_parts+addpart;q2++)
{
if(parttype==0)
{
if(q2==num_parts){fprintf(output2, "The relative contribution of the Base to each category");}
else{fprintf(output2, "The relative contribution of Annotation %d to each category", q2+1);}
}
else
{
if(q2==num_parts){fprintf(output2, "The relative contribution of the Background to each category");}
else{fprintf(output2, "The relative contribution of Partition %d to each category", q2+1);}
}

for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %.4f", ssums[q][q2]);}
fprintf(output2,"\n");
}

//get numbers of non-trivial heritability and regression snps
count=0;for(j=0;j<data_length;j++){count+=(mults[j]!=-9999&&windex[j]==1);}
count2=0;for(j=0;j<data_length;j++){count2+=(mults[j]!=-9999&&vindex[j]==1);}

fprintf(output2, "There are %d reference %d regression %d heritability predictors", data_length, count, count2);
for(q=0;q<num_parts+addpart;q++){fprintf(output2, " %d", (int)ssums[q][num_parts+addpart]);}
fprintf(output2, "\n");
}
else	//compute reduced version of ssums and print separately
{
sprintf(filename5,"%s.pathway.sums", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}

for(q=0;q<num_parts;q++)
{
//work out contributions for base, genic and annotation q
ssums[0][0]=0;ssums[0][1]=0;ssums[0][2]=0;ssums[0][3]=0;
ssums[1][0]=0;ssums[1][1]=0;ssums[1][2]=0;ssums[1][3]=0;
ssums[2][0]=0;ssums[2][1]=0;ssums[2][2]=0;ssums[2][3]=0;

for(j=0;j<data_length;j++)
{
if(pindexes[num_parts][j]==1&&vindex[j]==1)	//predictor is part of base and a heritability predictor
{
//add on contributions to base, from the base, genic and annotation, plus number of snps
ssums[0][0]+=pweights[num_parts][j]*exps[j];
ssums[0][1]+=pweights[num_parts+1][j]*exps[j];
ssums[0][2]+=pweights[q][j]*exps[j];
ssums[0][3]++;
}

if(pindexes[num_parts+1][j]==1&&vindex[j]==1)	//predictor is part of genic and a heritability predictor
{
//add on contributions to base, from the base, genic and annotation, plus number of snps
ssums[1][0]+=pweights[num_parts][j]*exps[j];
ssums[1][1]+=pweights[num_parts+1][j]*exps[j];
ssums[1][2]+=pweights[q][j]*exps[j];
ssums[1][3]++;
}

if(pindexes[q][j]==1&&vindex[j]==1)	//predictor is part of q and a heritability predictor
{
//add on contributions to annotation q, from the base, genic and annotation, plus number of snps
ssums[2][0]+=pweights[num_parts][j]*exps[j];
ssums[2][1]+=pweights[num_parts+1][j]*exps[j];
ssums[2][2]+=pweights[q][j]*exps[j];
ssums[2][3]++;
}
}

fprintf(output5,"%s", catlabels[q]);
fprintf(output5," %.4f %.4f %.4f %d", ssums[0][0], ssums[0][1], ssums[0][2], (int)ssums[0][3]);
fprintf(output5," %.4f %.4f %.4f %d", ssums[1][0], ssums[1][1], ssums[1][2], (int)ssums[1][3]);
fprintf(output5," %.4f %.4f %.4f %d\n", ssums[2][0], ssums[2][1], ssums[2][2], (int)ssums[2][3]);
}
fclose(output5);
}

fclose(output2);
if(savemat==1){fclose(output3);}
if(cover==1){fclose(output4);}

if(strcmp(genpreds,"blank")==0)
{
printf("Taggings saved in %s", filename2);
if(savemat==1&&cover==0){printf(", with heritability matrix in %s", filename3);}
if(savemat==0&&cover==1){printf(", with coverages in %s", filename4);}
if(savemat==1&&cover==1){printf(", with heritability matrix in %s and coverages in %s", filename3, filename4);}
printf("\n\n");
}
else{printf("Base taggings saved in %s, with pathway taggings in %s\n\n", filename2, filename5);}

free(windex);free(vindex);
free(data);
free(exps);
free(cors);
free(tally1);free(tally2);free(tally3);if(cover==1){free(tally4);}
if(strcmp(genpreds,"blank")==0){for(q=0;q<num_parts+addpart+addgenic;q++){free(ssums[q]);}free(ssums);}
else{free(ssums[0]);free(ssums[1]);free(ssums[2]);free(ssums);}
if(binary==0){gzclose(datainputgz);}
}	//end of mode=141

////////

if(mode==145)	//calc-overlap
{
//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
exps=malloc(sizeof(double)*data_length);

ssums=malloc(sizeof(double*)*(num_parts+addpart));
for(q=0;q<num_parts+addpart;q++){ssums[q]=malloc(sizeof(double)*(num_parts+addpart+2));}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress and on-the-fly files
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

printf("Reading data for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Reading data for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 0);
}	//end of bit loop
printf("\n");

//compute ssums - slightly different to the ssums used with tagging files
//ssums[q][q2] indicates fraction of category q2 that contributes to category q predictors
//ssums[q][num_parts+addpart] is number of predictors in category
//ssums[q][num_parts+addpart+1] is expected share of heritability from category q predictors

//first get exps based on just weights and power (ignoring pweights)
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}
else{exps[j]=0;}
}

//now multiply pweights by exps, then make sure scale to one
for(q=0;q<num_parts+addpart;q++)
{
for(j=0;j<data_length;j++){pweights[q][j]=pweights[q][j]*exps[j];}
sum=0;for(j=0;j<data_length;j++){sum+=pweights[q][j];}
for(j=0;j<data_length;j++){pweights[q][j]=pweights[q][j]/sum;}
}

//can now set ssums
for(q=0;q<num_parts+addpart;q++)
{
for(q2=0;q2<num_parts+addpart;q2++){ssums[q][q2]=0;}
ssums[q][num_parts+addpart]=0;

for(j=0;j<data_length;j++)	//will include trivial snps (but these will have exp zero)
{
if(pindexes[q][j]==1)	//predictor is part of q
{
for(q2=0;q2<num_parts+addpart;q2++){ssums[q][q2]+=pweights[q2][j];}
ssums[q][num_parts+addpart]++;
}
}

ssums[q][num_parts+addpart+1]=ssums[q][num_parts+addpart]/data_length;
}

sprintf(filename2,"%s.overs",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

for(q=0;q<num_parts+addpart;q++){fprintf(output2,"%s ", catlabels[q]);}
fprintf(output2,"\n");
for(q2=0;q2<num_parts+addpart+2;q2++)
{
for(q=0;q<num_parts+addpart;q++){fprintf(output2,"%.6f ", ssums[q][q2]);}
fprintf(output2,"\n");
}
fclose(output2);

printf("Overlap proportions saved in %s\n\n", filename2);

free(data);
free(exps);
for(q=0;q<num_parts+addpart;q++){free(ssums[q]);}free(ssums);
if(binary==0){gzclose(datainputgz);}
}	//end of mode=145

////////

for(q=0;q<num_parts+1;q++){free(catlabels[q]);}free(catlabels);
for(q=0;q<num_parts+2;q++){free(pindexes[q]);free(pweights[q]);}free(pindexes);free(pweights);
free(usedpreds);

///////////////////////////

