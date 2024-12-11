/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Merge datasets based on samples or predictors - can not be here if dtype=2

///////////////////////////

//deal with samples

Xnall=malloc(sizeof(int)*num_files);
Xids1=malloc(sizeof(char**)*num_files);
Xids2=malloc(sizeof(char**)*num_files);
Xids3=malloc(sizeof(char**)*num_files);
Xnuse=malloc(sizeof(int)*num_files);
Xks=malloc(sizeof(int*)*num_files);
Xks2=malloc(sizeof(int*)*num_files);

//first read in all ids and check unique
for(k=0;k<num_files;k++)
{
Xnall[k]=countrows(famstems[k])-famhead;
if(Xnall[k]==0){printf("Error, %s contains no samples\n\n", famstems[k]);exit(1);}
printf("Reading IDs for %d samples from %s\n", Xnall[k], famstems[k]);
Xids1[k]=malloc(sizeof(char*)*Xnall[k]);
Xids2[k]=malloc(sizeof(char*)*Xnall[k]);
Xids3[k]=malloc(sizeof(char*)*Xnall[k]);
read_ids(famstems[k], Xids1[k], Xids2[k], Xids3[k], Xnall[k], NULL, famhead, 0);
idsorder=malloc(sizeof(int)*Xnall[k]);
check_dups(Xids3[k],Xnall[k],famstems[k],idsorder,1);
free(idsorder);
}

//set csamps (0 if merging predictors, 1 if merging samples) and get num_samples and allids3

if(num_files==1)	//always have csamps 1
{csamps=1;}
else	//see whether allowed csamps 0
{
//get list of unique samples
num_samples=0;for(k=0;k<num_files;k++){num_samples+=Xnall[k];}
allids3=malloc(sizeof(char*)*num_samples);
count=0;
for(k=0;k<num_files;k++)
{
for(i=0;i<Xnall[k];i++){copy_string(allids3,count,Xids3[k][i]);count++;}
}

//set csamps assuming no duplicates, then change if there are
csamps=0;
wantids=malloc(sizeof(char*)*num_samples);
for(i=0;i<num_samples;i++){copy_string(wantids,i,allids3[i]);}
qsort(wantids,num_samples,sizeof(char *), compare_string);
for(i=1;i<num_samples;i++)
{
if(strcmp(wantids[i-1],wantids[i])==0){csamps=1;break;}
}

for(i=0;i<num_samples;i++){free(wantids[i]);}free(wantids);
if(csamps==1)
{
for(i=0;i<num_samples;i++){free(allids3[i]);}free(allids3);
}
}

if(mode==190&&csamps==0){printf("Error, there are no samples common to all %d datasets\n\n", k+1);exit(1);}

if(csamps==1)	//get list of overlapping samples (or give error if none)
{
//set num_samples and allids3 based on first dataset
num_samples=Xnall[0];
allids3=malloc(sizeof(char*)*num_samples);
for(i=0;i<num_samples;i++){copy_string(allids3,i,Xids3[0][i]);}

indexer=malloc(sizeof(int)*num_samples);
for(k=1;k<num_files;k++)	//update num_samples and allids3 based on which samples also in dataset k+1
{
count=num_samples;
num_samples=find_strings(allids3, count, Xids3[k], Xnall[k], indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(num_samples==0)
{
if(mode!=190){printf("Error, it is only possible to merge datasets that contain distinct samples, or that have samples common to all datasets\n\n");}
else{printf("Error, there are no samples common to all %d datasets\n\n", k+1);}
exit(1);
}

//squeeze down samples
for(i=0;i<num_samples;i++)
{
if(indexer[i]!=i){free(allids3[i]);copy_string(allids3,i, allids3[indexer[i]]);}
}
for(i=num_samples;i<count;i++){free(allids3[i]);}
}
free(indexer);
}

if(num_files>1)
{
if(csamps==1){printf("\nThere are %d samples common to all %d datasets\n", num_samples, num_files);}
else{printf("\nIn total, the %d datasets contain %d samples\n", num_files, num_samples);}
}

////////

//consider filterings
usedids=malloc(sizeof(int)*num_samples);
for(i=0;i<num_samples;i++){usedids[i]=1;}

if(strcmp(bsampfile,"blank")!=0)	//keep samples
{
count=countrows(bsampfile);
printf("Reading list of %d samples to keep from %s\n", count, bsampfile);
wantids=malloc(sizeof(char*)*count);
read_ids(bsampfile, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, can not find any of these\n\n");exit(1);}
if(count2<count){printf("Warning, can only find %d of these\n", count2);}
printf("\n");
for(i=0;i<count2;i++){usedids[indexer[i]]++;}
for(i=0;i<num_samples;i++){usedids[i]=(usedids[i]==2);}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);
}

if(strcmp(csampfile,"blank")!=0)	//remove samples
{
count=countrows(csampfile);
printf("Reading list of %d samples to remove from %s", count, csampfile);
if(strcmp(bsampfile,"blank")!=0){printf(" (takes priority over \"--keep\")");}
printf("\n");
wantids=malloc(sizeof(char*)*count);
read_ids(csampfile, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Warning, can not find any of these\n");}
if(count2<count){printf("Warning, can only find %d of these\n", count2);}
printf("\n");
for(i=0;i<count2;i++){usedids[indexer[i]]=0;}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);
}

//reduce
num_samples_use=0;
for(i=0;i<num_samples;i++){num_samples_use+=(usedids[i]==1);}
if(num_samples_use==0){printf("Error, after filtering samples, none remain\n\n");exit(1);}
if(num_samples_use<num_samples){printf("Will be using %d of the %d samples\n", num_samples_use, num_samples);}
printf("\n");

//decide if we can use quickmerge (can use for bed files if keeping all samples and fam files are the same)
qmerge=0;
if(mode==181&&csamps==1&&dtype==1&&num_samples_use==num_samples)
{
qmerge=1;
for(k=0;k<num_files;k++)
{
if(Xnall[k]!=num_samples_use){qmerge=0;break;}
for(i=0;i<Xnall[k];i++)
{
if(strcmp(Xids3[k][i],allids3[i])!=0){qmerge=0;break;}
}
}
}

//fill up ids
ids1=malloc(sizeof(char*)*num_samples_use);
ids2=malloc(sizeof(char*)*num_samples_use);
ids3=malloc(sizeof(char*)*num_samples_use);

//ids3 are easy
count=0;
for(i=0;i<num_samples;i++)
{
if(usedids[i]==1){copy_string(ids3,count,allids3[i]);count++;}
}
free(usedids);

//match up ids1/2 which also gets indexes
for(k=0;k<num_files;k++)
{
Xks[k]=malloc(sizeof(int)*Xnall[k]);
Xks2[k]=malloc(sizeof(int)*Xnall[k]);
Xnuse[k]=find_strings(Xids3[k], Xnall[k], ids3, num_samples_use, Xks[k], Xks2[k], NULL, NULL, NULL, NULL, 3);
if(Xnuse[k]==0){printf("Error, after filtering samples, none remain in Dataset %d\n\n", k+1);exit(1);}

for(i=0;i<Xnuse[k];i++)
{
i2=Xks[k][i];i3=Xks2[k][i];
if(csamps==0||k==0)	//load up details for this sample
{
copy_string(ids1,i3,Xids1[k][i2]);
copy_string(ids2,i3,Xids2[k][i2]);
}
else	//in theory, there can be mismatches
{
if(strcmp(Xids1[k][i2],ids1[i3])!=0||strcmp(Xids2[k][i2],ids2[i3])!=0)
{printf("Error matching sample IDs (%s %s and %s %s); this has never happened before, so please tell Doug\n\n", Xids1[k][i2], ids1[i3], Xids2[k][i2], ids2[i3]);exit(1);}
}
}
}

//now get additional details (unless mode=190)

pid=malloc(sizeof(char*)*num_samples_use);
mid=malloc(sizeof(char*)*num_samples_use);
schar=malloc(sizeof(char*)*num_samples_use);
pchar=malloc(sizeof(char*)*num_samples_use);
for(i=0;i<num_samples_use;i++)
{
pid[i]=malloc(sizeof(char)*2);strcpy(pid[i],"0");
mid[i]=malloc(sizeof(char)*2);strcpy(mid[i],"0");
schar[i]=malloc(sizeof(char)*2);strcpy(schar[i],"0");
pchar[i]=malloc(sizeof(char)*3);strcpy(pchar[i],"NA");
}

if(mode!=190&&famhead==0)	//read in values (if using fam files with at least 6 columns)
{
flag=0;
for(k=0;k<num_files;k++)
{
if(countcols(famstems[k])<6)
{
printf("Warning, %s has only %d columns; famfiles should generally have (at least) 6 columns\n\n", famstems[k], countcols(famstems[k]));
flag=1;break;
}
}

if(flag==0)
{
erra=0;errb=0;errc=0;errd=0;
for(k=0;k<num_files;k++)
{
wantids=malloc(sizeof(char*)*Xnuse[k]);

read_strings(famstems[k], wantids, Xnuse[k], Xks[k], 3, 0);
for(i=0;i<Xnuse[k];i++)
{
i2=Xks2[k][i];
if(erra==0){erra=check_copy(pid, i2, wantids[i], "0", ids1[i2], ids2[i2], famstems[k], 1);}
free(wantids[i]);
}
read_strings(famstems[k], wantids, Xnuse[k], Xks[k], 4, 0);
for(i=0;i<Xnuse[k];i++)
{
i2=Xks2[k][i];
if(errb==0){errb=check_copy(mid, i2, wantids[i], "0", ids1[i2], ids2[i2], famstems[k], 2);}
free(wantids[i]);
}
read_strings(famstems[k], wantids, Xnuse[k], Xks[k], 5, 0);
for(i=0;i<Xnuse[k];i++)
{
i2=Xks2[k][i];
if(errc==0){errc=check_copy(schar, i2, wantids[i], "0", ids1[i2], ids2[i2], famstems[k], 3);}
free(wantids[i]);
}
read_strings(famstems[k], wantids, Xnuse[k], Xks[k], 6, 0);
for(i=0;i<Xnuse[k];i++)
{
i2=Xks2[k][i];
if(errd==0){errd=check_copy(pchar, i2, wantids[i], "NA", ids1[i2], ids2[i2], famstems[k], 4);}
free(wantids[i]);
}

free(wantids);
}
}

//reset for those that failed
for(i=0;i<num_samples_use;i++)
{
if(erra==1){strcpy(pid[i],"0");}
if(errb==1){strcpy(mid[i],"0");}
if(errc==1){strcpy(schar[i],"0");}
if(errd==1){strcpy(pchar[i],"NA");}
}
}	//end of reading extra details

//can now free a few things
for(k=0;k<num_files;k++)
{
for(i=0;i<Xnall[k];i++){free(Xids1[k][i]);free(Xids2[k][i]);free(Xids3[k][i]);}
free(Xids1[k]);free(Xids2[k]);free(Xids3[k]);
}
free(Xids1);free(Xids2);free(Xids3);
for(i=0;i<num_samples;i++){free(allids3[i]);}free(allids3);

///////////////////////////

//now predictors - if filtering, will do so on-the-fly (must be careful not to use count and count2)

XNall=malloc(sizeof(int*)*num_files);
Xchr=malloc(sizeof(int*)*num_files);
Xcm=malloc(sizeof(double*)*num_files);
Xbp=malloc(sizeof(double*)*num_files);
Xpreds=malloc(sizeof(char**)*num_files);
Xalong1=malloc(sizeof(char**)*num_files);
Xalong2=malloc(sizeof(char**)*num_files);
Xal1=malloc(sizeof(char*)*num_files);
Xal2=malloc(sizeof(char*)*num_files);
XNuse=malloc(sizeof(int)*num_files);
Xkp=malloc(sizeof(int*)*num_files);
Xkp2=malloc(sizeof(int*)*num_files);

//can not have both extract and snp
count=0;
if(strcmp(bpredfile,"blank")!=0)
{
count=countrows(bpredfile);
printf("Reading list of %d predictors to extract from %s\n\n", count, bpredfile);
wantpreds=malloc(sizeof(char*)*count);
read_strings(bpredfile, wantpreds, count, NULL, 1, 0);
}
if(strcmp(onesnp,"blank")!=0)
{
count=1;
wantpreds=malloc(sizeof(char*));
wantpreds[0]=malloc(sizeof(char)*(strlen(onesnp)+1));
strcpy(wantpreds[0],onesnp);
}

count2=0;
if(strcmp(cpredfile,"blank")!=0)	//exclude predictors
{
count2=countrows(cpredfile);
printf("Reading list of %d predictors to exclude from %s (takes priority over any other predictor filterings)\n", count2, cpredfile);
wantpreds2=malloc(sizeof(char*)*count2);
read_strings(cpredfile, wantpreds2, count2, NULL, 1, 0);
}

////////

for(k=0;k<num_files;k++)
{
if(strcmp(bimstems[k],datastems[k])!=0)	//have a separate bimfile
{
if(countcols(bimstems[k])!=6){printf("Error, %s should have six columns (not %d); columns should be chr, name, genetic distance, bp, A1, A2\n\n", bimstems[k], countcols(bimstems[k]));exit(1);}

XNall[k]=countrows_plus(bimstems[k], 6);
printf("Reading details for %d predictors from %s\n", XNall[k], bimstems[k]);

Xchr[k]=malloc(sizeof(int)*XNall[k]);
Xcm[k]=malloc(sizeof(double)*XNall[k]);
Xbp[k]=malloc(sizeof(double)*XNall[k]);
Xpreds[k]=malloc(sizeof(char*)*XNall[k]);
Xalong1[k]=malloc(sizeof(char*)*XNall[k]);
Xalong2[k]=malloc(sizeof(char*)*XNall[k]);
Xal1[k]=malloc(sizeof(char)*XNall[k]);
Xal2[k]=malloc(sizeof(char)*XNall[k]);

//predictors only need to be ordered if data not binary
read_bimfile(bimstems[k], Xchr[k], Xpreds[k], Xcm[k], Xbp[k], Xalong1[k], Xalong2[k], Xal1[k], Xal2[k], XNall[k], 1, -9999, binary, exlong);
}
else	//no bimfile, must be dtype=5
{
//assume at most maxpreds SNPs
printf("Reading predictor details from %s\n", bimstems[k]);
Xchr[k]=malloc(sizeof(int)*maxpreds);
Xcm[k]=malloc(sizeof(double)*maxpreds);
Xbp[k]=malloc(sizeof(double)*maxpreds);
Xpreds[k]=malloc(sizeof(char*)*maxpreds);
Xalong1[k]=malloc(sizeof(char*)*maxpreds);
Xalong2[k]=malloc(sizeof(char*)*maxpreds);
Xal1[k]=malloc(sizeof(char)*maxpreds);
Xal2[k]=malloc(sizeof(char)*maxpreds);
for(j=0;j<maxpreds;j++){Xchr[k][j]=oxchr;Xcm[k][j]=0;}

//predictors need to be ordered
XNall[k]=read_genfile(datastems[k], Xpreds[k], Xbp[k], Xalong1[k], Xalong2[k], Xal1[k], Xal2[k], Xnall[k], genskip, genheaders, genprobs, maxpreds, 1, 0, exlong);

if(XNall[k]>maxpreds)	//then file larger than allocated, will have to read again
{
for(j=0;j<maxpreds;j++){free(Xpreds[k][j]);free(Xalong1[k][j]);free(Xalong2[k][j]);}free(Xpreds[k]);free(Xalong1[k]);free(Xalong2[k]);
free(Xchr[k]);free(Xcm[k]);free(Xbp[k]);free(Xal1[k]);free(Xal2[k]);
Xchr[k]=malloc(sizeof(int)*XNall[k]);
Xcm[k]=malloc(sizeof(double)*XNall[k]);
Xbp[k]=malloc(sizeof(double)*XNall[k]);
Xpreds[k]=malloc(sizeof(char*)*XNall[k]);
Xalong1[k]=malloc(sizeof(char*)*XNall[k]);
Xalong2[k]=malloc(sizeof(char*)*XNall[k]);
Xal1[k]=malloc(sizeof(char)*XNall[k]);
Xal2[k]=malloc(sizeof(char)*XNall[k]);
for(j=0;j<XNall[k];j++){Xchr[k][j]=oxchr;Xcm[k][j]=0;}

(void)read_genfile(datastems[k], Xpreds[k], Xbp[k], Xalong1[k], Xalong2[k], Xal1[k], Xal2[k], Xnall[k], genskip, genheaders, genprobs, XNall[k], 1, 0, exlong);
}

if(XNall[k]>=20000){printf("%s contains genotypes for %d SNPs\n", datastems[k], XNall[k]);}
}

if(exlong==0)	//see if any long alleles
{
count3=0;for(j=0;j<XNall[k];j++){count3+=(Xbp[k][j]>0&&(strlen(Xalong1[k][j])+strlen(Xalong2[k][j])>2));}
if(count3>0){printf("Warning, %d predictors have multi-character alleles (you can remove these using \"--exclude-long-alleles YES\")\n", count3);}
}

////////

//check there are some predictors with positive bp
count3=0;for(j=0;j<XNall[k];j++){count3+=(Xbp[k][j]>0);}
if(count3==0)
{
if(exlong==0){printf("Error, there are no predictors with non-positive basepairs in %s\n\n", bimstems[k]);}
else{printf("Error, there are no predictors with non-positive basepairs and single-character alleles in %s\n\n", bimstems[k]);}
exit(1);
}

//deal with onechr - exclude by giving basepair -1
if(onechr>0)
{
count3=0;
for(j=0;j<XNall[k];j++)
{
if(Xchr[k][j]==onechr){count3++;}
else{Xbp[k][j]=-1;}
}
if(count3==0)
{printf("Error, there are no Chromosome %d predictors in %s\n", onechr, bimstems[k]);exit(1);}
printf("There are %d Chromosome %d predictors in %s\n", count3, onechr, bimstems[k]);
}
if(onechr==-1)
{
count3=0;
for(j=0;j<XNall[k];j++)
{
if(Xchr[k][j]>0&&Xchr[k][j]<23){count3++;}
else{Xbp[k][j]=-1;}
}
if(count3==0)
{printf("Error, there are no autosomal predictors in %s\n", bimstems[k]);exit(1);}
printf("There are %d autosomal predictors in %s\n", count3, bimstems[k]);
}
if(onechr==-3)
{
count3=0;
for(j=0;j<XNall[k];j++)
{
if(Xchr[k][j]>0&&Xchr[k][j]%2==1){count3++;}
else{Xbp[k][j]=-1;}
}
if(count3==0)
{printf("Error, there are no odd-chromosome predictors in %s\n", bimstems[k]);exit(1);}
printf("There are %d odd-chromosome predictors in %s\n", count3, bimstems[k]);
}
if(onechr==-2)
{
count3=0;
for(j=0;j<XNall[k];j++)
{
if(Xchr[k][j]>0&&Xchr[k][j]%2==0){count3++;}
else{Xbp[k][j]=-1;}
}
if(count3==0)
{printf("Error, there are no even-chromosome predictors in %s\n", bimstems[k]);exit(1);}
printf("There are %d even-chromosome predictors in %s\n", count3, bimstems[k]);
}

//get predictor order and check for duplicates names
predorder=malloc(sizeof(int)*XNall[k]);
check_dups(Xpreds[k],XNall[k],NULL,predorder, 0);

//get to first retained predictor (must be at least one)
found=0;while(Xbp[k][predorder[found]]<=0){found++;}

//check if names are unique
wcount=0;
for(j=found+1;j<XNall[k];j++)
{
if(Xbp[k][predorder[j]]>0)
{
if(strcmp(Xpreds[k][predorder[j]],Xpreds[k][predorder[found]])==0)	//name matches
{
if(wcount<5)
{
if(exsame==0){printf("Warning, at least two predictors are called %s (you can remove these using \"--exclude-same-names YES\")\n", allpreds[predorder[j]]);}
else
{
printf("Warning, at least two predictors are called %s; these will be ignored\n", allpreds[predorder[j]]);
Xbp[k][predorder[found]]=-1;
Xbp[k][predorder[j]]=-1;
}
}
wcount++;

//if three or more with same name, don't want to keep warning
while(j<XNall[k]-1)
{
j++;
if(Xbp[k][predorder[j]]>0)
{
if(strcmp(Xpreds[k][predorder[j]],Xpreds[k][predorder[found]])!=0){break;}
if(exsame==1){Xbp[k][predorder[j]]=-1;}
}
}
}
found=j;
}}
if(wcount>5){printf("In total, %d predictor names appear more than once\n", wcount);}

//check there are still some predictors with positive bp
count3=0;for(j=0;j<XNall[k];j++){count3+=(Xbp[k][j]>0);}
if(count3==0){printf("Error, there are no predictors with unique names\n\n");exit(1);}

if(count>0||count2>0)	//take care of extract, onesnp and exclude - exclude by giving basepair -1
{
if(count>0)	//extract
{
indexer=malloc(sizeof(int)*XNall[k]);
usedpreds=malloc(sizeof(int)*XNall[k]);

count3=find_strings(Xpreds[k], XNall[k], wantpreds, count, indexer, NULL, NULL, NULL, predorder, NULL, 3);
if(count3==0)
{
if(strcmp(bpredfile,"blank")!=0){printf("Error, none of the %d predictors in %s are in %s\n\n", count, bpredfile, bimstems[k]);exit(1);}
else{printf("Error, Predictor %s is not in %s\n\n", wantpreds[0], bimstems[k]);exit(1);}
}
if(count3<count){printf("Warning, only %d of the %d predictors in %s are in %s\n", count3, count, bpredfile, bimstems[k]);}

for(j=0;j<XNall[k];j++){usedpreds[j]=0;}
for(j=0;j<count3;j++){usedpreds[indexer[j]]=1;}
for(j=0;j<XNall[k];j++)
{
if(usedpreds[j]==0){Xbp[k][j]=-1;}
}
free(indexer);free(usedpreds);
}

if(count2>0)	//exclude
{
indexer=malloc(sizeof(int)*XNall[k]);
count3=find_strings(Xpreds[k], XNall[k], wantpreds2, count2, indexer, NULL, NULL, NULL, predorder, NULL, 3);
if(count3==0){printf("Warning, none of the %d predictors in %s are in %s\n\n", count2, cpredfile, bimstems[k]);}
if(count3>0&&count3<count2){printf("Warning, only %d of the %d predictors in %s are in %s\n\n", count3, count2, cpredfile, bimstems[k]);}
for(j=0;j<count3;j++){Xbp[k][indexer[j]]=-1;}
free(indexer);
}

count3=0;for(j=0;j<XNall[k];j++){count3+=(Xbp[k][j]>0);}
if(count3==0){printf("Error, after filtering predictors, none remain in %s\n\n", bimstems[k]);exit(1);}
}

free(predorder);
printf("\n");
}	//end of k loop

if(count>0)
{
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);
}
if(count2>0)
{
for(j=0;j<count2;j++){free(wantpreds2[j]);}free(wantpreds2);
}

////////

//get num_preds_use and predictor details
if(mode!=190&&csamps==1)	//not allowed predictors in multiple datasets (after allowing for filtering)
{
//get details for all predictors
num_preds_use=0;
for(k=0;k<num_files;k++)
{
for(j=0;j<XNall[k];j++){num_preds_use+=(Xbp[k][j]>0);}
}

allchr=malloc(sizeof(int)*num_preds_use);
allcm=malloc(sizeof(double)*num_preds_use);
allbp=malloc(sizeof(double)*num_preds_use);
allpreds=malloc(sizeof(char*)*num_preds_use);
allalong1=malloc(sizeof(char*)*num_preds_use);
allalong2=malloc(sizeof(char*)*num_preds_use);
allal1=malloc(sizeof(char)*num_preds_use);
allal2=malloc(sizeof(char)*num_preds_use);

count=0;
for(k=0;k<num_files;k++)
{
for(j=0;j<XNall[k];j++)
{
if(Xbp[k][j]>0)
{
allchr[count]=Xchr[k][j];
allcm[count]=Xcm[k][j];
allbp[count]=Xbp[k][j];
copy_string(allpreds,count,Xpreds[k][j]);
copy_string(allalong1,count,Xalong1[k][j]);
copy_string(allalong2,count,Xalong2[k][j]);
allal1[count]=Xal1[k][j];
allal2[count]=Xal2[k][j];
count++;
}
}
}

//check for duplicates (use wantpreds so that preserve order of samples in allpreds)
wantpreds=malloc(sizeof(char*)*num_preds_use);
for(j=0;j<num_preds_use;j++){copy_string(wantpreds,j,allpreds[j]);}
qsort(wantpreds,num_preds_use,sizeof(char *), compare_string);
for(j=1;j<num_preds_use;j++)
{
if(strcmp(wantpreds[j-1],wantpreds[j])==0){printf("Error, Predictor %s appears in multiple datasets (this is only allowed when merging datasets containing distinct samples)\n\n", wantpreds[j]);exit(1);}
}
for(j=0;j<num_preds_use;j++){free(wantpreds[j]);}free(wantpreds);
}
else	//must be predictors common to all datasets (after allowing for filtering)
{
//set predictor details based on first dataset
num_preds_use=0;for(j=0;j<XNall[0];j++){num_preds_use+=(Xbp[0][j]>0);}
allchr=malloc(sizeof(int)*num_preds_use);
allcm=malloc(sizeof(double)*num_preds_use);
allbp=malloc(sizeof(double)*num_preds_use);
allpreds=malloc(sizeof(char*)*num_preds_use);
allalong1=malloc(sizeof(char*)*num_preds_use);
allalong2=malloc(sizeof(char*)*num_preds_use);
allal1=malloc(sizeof(char)*num_preds_use);
allal2=malloc(sizeof(char)*num_preds_use);

count=0;
for(j=0;j<XNall[0];j++)
{
if(Xbp[0][j]>0)
{
allchr[count]=Xchr[0][j];
allcm[count]=Xcm[0][j];
allbp[count]=Xbp[0][j];
copy_string(allpreds,count,Xpreds[0][j]);
copy_string(allalong1,count,Xalong1[0][j]);
copy_string(allalong2,count,Xalong2[0][j]);
allal1[count]=Xal1[0][j];
allal2[count]=Xal2[0][j];
count++;
}
}

for(k=1;k<num_files;k++)	//update based on which predictors also in dataset k+1 (allowing for filtering)
{
wantpreds=malloc(sizeof(char*)*XNall[k]);
count2=0;
for(j=0;j<XNall[k];j++)
{
if(Xbp[k][j]>0){copy_string(wantpreds,count2,Xpreds[k][j]);count2++;}
}
indexer=malloc(sizeof(int)*num_preds_use);

count=num_preds_use;
num_preds_use=find_strings(allpreds, count, wantpreds, count2, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(num_preds_use==0)
{
if(mode==190){printf("Error, there are no predictors common to the first %d datasets\n\n", k+1);}
else{printf("Error, there are no predictors common to the first %d datasets (this is only allowed when merging datasets containing overlapping samples)\n\n", k+1);}
exit(1);
}

//squeeze down predictors
for(j=0;j<num_preds_use;j++)
{
if(indexer[j]!=j)
{
allchr[j]=allchr[indexer[j]];
allcm[j]=allcm[indexer[j]];
allbp[j]=allbp[indexer[j]];
free(allpreds[j]);copy_string(allpreds,j, allpreds[indexer[j]]);
free(allalong1[j]);copy_string(allalong1,j, allalong1[indexer[j]]);
free(allalong2[j]);copy_string(allalong2,j, allalong2[indexer[j]]);
allal1[j]=allal1[indexer[j]];
allal2[j]=allal2[indexer[j]];
}
}
for(j=num_preds_use;j<count;j++){free(allpreds[j]);free(allalong1[j]);free(allalong2[j]);}

for(j=0;j<count2;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}
}

if(num_files==1)
{
if(num_samples_use<Xnall[0]||num_preds_use<XNall[0]){printf("Will be using %d samples and %d predictors\n\n", num_samples_use, num_preds_use);}
}
else
{
if(mode!=190&&csamps==1){printf("In total, the %d datasets contain %d predictors\n\n", num_files, num_preds_use);}
else{printf("There are %d predictors common to all %d datasets\n\n", num_preds_use, num_files);}
}

//preds must contain allpreds, sorted by position

chr=malloc(sizeof(int)*num_preds_use);
cm=malloc(sizeof(double)*num_preds_use);
bp=malloc(sizeof(double)*num_preds_use);
preds=malloc(sizeof(char*)*num_preds_use);
along1=malloc(sizeof(char*)*num_preds_use);
along2=malloc(sizeof(char*)*num_preds_use);
al1=malloc(sizeof(char)*num_preds_use);
al2=malloc(sizeof(char)*num_preds_use);

//see if already sorted
flag=0;
for(j=1;j<num_preds_use;j++)
{
if(allchr[j]<allchr[j-1]||(allchr[j]==allchr[j-1]&&allbp[j]<allbp[j-1])){flag=1;break;}
}

if(flag==0)	//already sorted
{
for(j=0;j<num_preds_use;j++)
{
chr[j]=allchr[j];cm[j]=allcm[j];bp[j]=allbp[j];
copy_string(preds,j,allpreds[j]);
copy_string(along1,j,allalong1[j]);
copy_string(along2,j,allalong2[j]);
al1[j]=allal1[j];al2[j]=allal2[j];
}
}
else	//necessary to sort
{
//get largest chromosome
count=0;
for(j=0;j<num_preds_use;j++)
{
if(allchr[j]>count){count=allchr[j];}
}

dptrs=malloc(sizeof(struct sorting_double)*num_preds_use);
count2=0;
for(k=0;k<=count;k++)	//find and sort predictors on chr k
{
count3=0;
for(j=0;j<num_preds_use;j++)
{
if(allchr[j]==k){dptrs[count3].value=allbp[j];dptrs[count3].index=j;count3++;}
}

if(count3>0)
{
qsort(dptrs, count3, sizeof(struct sorting_double), compare_sorting_double);
for(j=0;j<count3;j++)
{
j2=dptrs[j].index;
chr[count2+j]=allchr[j2];cm[count2+j]=allcm[j2];bp[count2+j]=allbp[j2];
copy_string(preds,count2+j,allpreds[j2]);
copy_string(along1,count2+j,allalong1[j2]);
copy_string(along2,count2+j,allalong2[j2]);
al1[count2+j]=allal1[j2];al2[count2+j]=allal2[j2];
}
count2+=count3;
}
}

free(dptrs);
}

free(allchr);free(allcm);free(allbp);free(allal1);free(allal2);
for(j=0;j<num_preds_use;j++){free(allpreds[j]);free(allalong1[j]);free(allalong2[j]);}free(allpreds);free(allalong1);free(allalong2);

//warn if there are duplicate positions
found=0;
wcount=0;
for(j=1;j<num_preds_use;j++)
{
if(chr[j]==chr[found]&&bp[j]==bp[found])	//a duplicate
{
if(wcount<5){printf("Warning, there are multiple predictors at Chromosome %d, Basepair %.1f)\n", chr[found], bp[found]);}
wcount++;
}
else{found=j;}
}
if(wcount>5){printf("In total, there are multiple predictors at %d positions\n", wcount);}
if(wcount>0){printf("\n");}

////////

//get indexes and check details consistent (if binary=0, indexes must be ordered)
predorder=malloc(sizeof(int)*num_preds_use);
check_dups(preds,num_preds_use,NULL,predorder,0);

if(mode!=190&&csamps==1)	//each predictor in only one dataset
{
for(k=0;k<num_files;k++)
{
Xkp[k]=malloc(sizeof(int)*XNall[k]);
Xkp2[k]=malloc(sizeof(int)*XNall[k]);

wantpreds=malloc(sizeof(char*)*XNall[k]);
indexer=malloc(sizeof(int)*XNall[k]);
indexer2=malloc(sizeof(int)*XNall[k]);
count=0;
for(j=0;j<XNall[k];j++)
{
if(Xbp[k][j]>0){copy_string(wantpreds,count,Xpreds[k][j]);indexer[count]=j;count++;}
}

XNuse[k]=find_strings(wantpreds, count, preds, num_preds_use, indexer2, Xkp2[k], NULL, NULL, NULL, predorder, 3);
if(XNuse[k]==0){printf("Doug error 482NA - only found %d of %d in %d\n", XNuse[k], num_preds_use, k+1);exit(1);}
for(j=0;j<XNuse[k];j++){Xkp[k][j]=indexer[indexer2[j]];}

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);free(indexer2);
}
}
else	//each predictor is in all datasets
{
for(k=0;k<num_files;k++)
{
Xkp[k]=malloc(sizeof(int)*XNall[k]);
Xkp2[k]=malloc(sizeof(int)*XNall[k]);

wantpreds=malloc(sizeof(char*)*XNall[k]);
indexer=malloc(sizeof(int)*XNall[k]);
indexer2=malloc(sizeof(int)*XNall[k]);
count=0;
for(j=0;j<XNall[k];j++)
{
if(Xbp[k][j]>0){copy_string(wantpreds,count,Xpreds[k][j]);indexer[count]=j;count++;}
}

XNuse[k]=find_strings(wantpreds, count, preds, num_preds_use, indexer2, Xkp2[k], NULL, NULL, NULL, predorder, 3);
if(XNuse[k]!=num_preds_use){printf("Doug error 482NCB - only found %d of %d in %d\n", XNuse[k], num_preds_use, k+1);exit(1);}
for(j=0;j<XNuse[k];j++){Xkp[k][j]=indexer[indexer2[j]];}

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);free(indexer2);

//check predictor details match or consistent with those in dataset 1 (probably redundant for dataset 1)
for(j=0;j<XNuse[k];j++)
{
j2=Xkp[k][j];j3=Xkp2[k][j];
if(Xchr[k][j2]!=chr[j3]){printf("Error, the chromosome of Predictor %s differs between datasets (%d and %d)\n\n", preds[j3], Xchr[k][j2], chr[j3]);exit(1);}

//no need to check genetic distances for calc-sim-data
if(mode!=190&&Xcm[k][j2]!=cm[j3]){printf("Error, the genetic distance of Predictor %s differs between datasets (%.4f and %.4f)\n\n", preds[j3], Xcm[k][j2], cm[j3]);exit(1);}

if(Xbp[k][j2]!=bp[j3]){printf("Error, the basepair of Predictor %s differs between datasets (%.1f and %.1f)\n\n", preds[j3], Xbp[k][j2], bp[j3]);exit(1);}
if((strcmp(Xalong1[k][j2],along1[j3])!=0||strcmp(Xalong2[k][j2],along2[j3])!=0)&&(strcmp(Xalong1[k][j2],along2[j3])!=0||strcmp(Xalong2[k][j2],along1[j3])!=0))
{printf("Error, the alleles of Predictor %s are not consistent across datasets (%s & %s and %s & %s)\n\n", preds[j3], Xalong1[k][j2], Xalong2[k][j2], along1[j3], along2[j3]);exit(1);}
}
}
}

free(predorder);

if(binary==0)	//make sure Xkp are in order
{
for(k=0;k<num_files;k++)
{
for(j=1;j<XNuse[k];j++)
{
if(Xkp[k][j]<=Xkp[k][j-1]){printf("Doug Error 385JXX - %d predictors in Dataset %d not ordered (%d %d)\n\n", XNuse[k], k+1, Xkp[k][j], Xkp[k][j-1]);exit(1);}
}
}
}

//free some things we no longer need
for(k=0;k<num_files;k++)
{
for(j=0;j<XNall[k];j++){free(Xpreds[k][j]);}free(Xpreds[k]);
free(Xchr[k]);free(Xcm[k]);free(Xbp[k]);free(Xal1[k]);free(Xal2[k]);
}
free(Xchr);free(Xcm);free(Xbp);free(Xpreds);free(Xal1);free(Xal2);

//must later free alongs and Xkps

///////////////////////////
//ready to look at data
///////////////////////////

if(qmerge==1)	//do separately
{
#include "rapid.c"
}
else
{
if(bitsize>num_preds_use){bitsize=num_preds_use;}

if(genprobs<2||mode==190){data_warn2(bitsize,2*num_samples_use);}
else{data_warn2(bitsize,4*num_samples_use);}

data=malloc(sizeof(double)*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);

Xdatainputgz=malloc(sizeof(gzFile)*num_files);
Xcurrent=malloc(sizeof(int)*num_files);
for(k=0;k<num_files;k++)
{
if(binary==0){open_datagz(Xdatainputgz+k, datastems[k], Xnall[k], genskip, genheaders, genprobs);}
Xcurrent[k]=0;
}

if(genprobs>=2&&mode!=190)
{
ps=malloc(sizeof(float*)*2);
ps[0]=malloc(sizeof(float)*num_samples_use*bitsize);
ps[1]=malloc(sizeof(float)*num_samples_use*bitsize);
ps2=malloc(sizeof(float*)*2);
ps2[0]=malloc(sizeof(float)*num_samples_use*bitsize);
ps2[1]=malloc(sizeof(float)*num_samples_use*bitsize);
}
else{ps=NULL;ps2=NULL;}

centres=malloc(sizeof(double)*num_preds_use);
mults=malloc(sizeof(double)*num_preds_use);
sqdevs=malloc(sizeof(double)*num_preds_use);
rates=malloc(sizeof(double)*num_preds_use);
infos=malloc(sizeof(double)*num_preds_use);

nums=malloc(sizeof(int*)*6);
for(k=0;k<6;k++)
{
nums[k]=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){nums[k][j]=0;}
}

//set nums to zero
for(k=0;k<6;k++)
{
for(j=0;j<num_preds_use;j++){nums[k][j]=0;}
}

////////

//open output files

sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

if(mode!=190)
{
sprintf(filename2,"%s.stats",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor\tA1\tA2\tA1_Mean\tMAF\tCall_Rate\tInfo\n");
}
else	//calculating correlations
{
sprintf(filename2,"%s.pred.cors", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Predictor A1\tA2\tNon_Missing\tA1_Mean_D1\tMAF_D1\tA1_Mean_D2\tMAF_D2\tCorrelation\n");
}

if(mode==181)	//bed
{
sprintf(filename3,"%s.bed", outfile);
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
onechar=108;fwrite(&onechar, sizeof(unsigned char), 1, output3);
onechar=27;fwrite(&onechar, sizeof(unsigned char), 1, output3);
onechar=1;fwrite(&onechar, sizeof(unsigned char), 1, output3);
}
if(mode==182)	//sp
{
sprintf(filename3,"%s.sp", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}
if(mode==183)	//sped
{
sprintf(filename3,"%s.sped", outfile);
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}
if(mode==184)	//speed
{
sprintf(filename3,"%s.speed", outfile);
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}
if(mode==185)	//gen
{
sprintf(filename3,"%s.gen", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}

////////

//ready for bit loop
bittotal=(num_preds_use-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>num_preds_use){bitend=num_preds_use;}
bitlength=bitend-bitstart;

if(bit%10==0)
{
if(mode!=190){printf("Making data for Chunk %d of %d\n", bit+1, bittotal);}
else{printf("Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);}
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mode!=190){fprintf(output, "Making data for Chunk %d of %d\n", bit+1, bittotal);}
else{fprintf(output, "Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);}
fclose(output);

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
}

for(k=0;k<num_files;k++)
{
//bitstart2 and bitend2 indicate which of Xkp2[k] are >= bitstart and < bitend
bitstart2=0;
while(bitstart2<XNuse[k])
{
if(Xkp2[k][bitstart2]>=bitstart){break;}
bitstart2++;
}
bitend2=bitstart2;
while(bitend2<XNuse[k])
{
if(Xkp2[k][bitend2]>=bitend){break;}
bitend2++;
}
bitlength2=bitend2-bitstart2;

if(bitlength2>0)	//will be reading from dataset k
{
if(genprobs>=2&&mode!=190)	//read data and probabilities - this is rarely used, and must be dtype=5
{Xcurrent[k]=read_data_fly(datastems[k], 5, data2, ps2, Xnuse[k], Xks[k], bitstart2, bitend2, Xkp[k], Xdatainputgz[k], Xcurrent[k], Xnall[k], XNall[k], genskip, genheaders, genprobs, NULL, missingvalue, threshold, minprob, nonsnp, -9999);}
else	//read just data
{
Xcurrent[k]=read_data_fly(datastems[k], dtype, data2, NULL, Xnuse[k], Xks[k], bitstart2, bitend2, Xkp[k], Xdatainputgz[k], Xcurrent[k], Xnall[k], XNall[k], genskip, genheaders, genprobs, NULL, missingvalue, threshold, minprob, nonsnp, maxthreads);}

//see if any predictors need turning (should only be necessary if predictors in multiple datasets)
#pragma omp parallel for private(j2, j, i2) schedule(static)
for(j2=0;j2<bitlength2;j2++)
{
j=Xkp2[k][bitstart2+j2]-bitstart;
if(strcmp(Xalong1[k][Xkp[k][bitstart2+j2]],along1[bitstart+j])!=0)	//turn
{
for(i2=0;i2<Xnuse[k];i2++)
{
if(data2[(size_t)j2*Xnuse[k]+i2]!=missingvalue)
{
data2[(size_t)j2*Xnuse[k]+i2]=2-data2[(size_t)j2*Xnuse[k]+i2];
if(genprobs>=2&&mode!=190)
{ps2[0][(size_t)j2*Xnuse[k]+i2]=1-ps2[0][(size_t)j2*Xnuse[k]+i2]-ps2[1][(size_t)j2*Xnuse[k]+i2];}
}}
}
}

if(mode!=190||k==0)	//load values into data
{
#pragma omp parallel for private(j2, j, i2) schedule(static)
for(j2=0;j2<bitlength2;j2++)
{
j=Xkp2[k][bitstart2+j2]-bitstart;
for(i2=0;i2<Xnuse[k];i2++){data[(size_t)j*num_samples_use+Xks2[k][i2]]=data2[(size_t)j2*Xnuse[k]+i2];}
if(genprobs>=2&&mode!=190)
{
for(i2=0;i2<Xnuse[k];i2++)
{ps[0][(size_t)j*num_samples_use+Xks2[k][i2]]=ps2[0][(size_t)j2*Xnuse[k]+i2];
ps[1][(size_t)j*num_samples_use+Xks2[k][i2]]=ps2[1][(size_t)j2*Xnuse[k]+i2];}
}
}	//end of j2 loop
}	//end of storing	
else	//must be 190 and k=1 so compare and save - will have Xnuse[k]=num_samples_use
{
for(j2=0;j2<bitlength2;j2++)
{
j=Xkp2[k][bitstart2+j2]-bitstart;

sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;indcount=0;
for(i2=0;i2<num_samples_use;i2++)
{
i=Xks2[k][i2];
if(data[(size_t)j*num_samples_use+i]!=missingvalue&&data2[(size_t)j2*Xnuse[k]+i2]!=missingvalue)
{
sum+=data[(size_t)j*num_samples_use+i];sum2+=data2[(size_t)j2*Xnuse[k]+i2];
sumsq+=pow(data[(size_t)j*num_samples_use+i],2);sumsq2+=pow(data2[(size_t)j2*Xnuse[k]+i2],2);
sumsq3+=data[(size_t)j*num_samples_use+i]*data2[(size_t)j2*Xnuse[k]+i2];
indcount++;
}}

fprintf(output2, "%s\t%s\t%s\t%d\t", preds[bitstart+j], along1[bitstart+j], along2[bitstart+j], indcount);
if(indcount>0)
{
mean=sum/indcount;mean2=sum2/indcount;
if(nonsnp==0)
{fprintf(output2, "%.6f\t%.6f\t%.6f\t%.6f\t", mean, mean/2+(mean>1)*(1-mean), mean2, mean2/2+(mean2>1)*(1-mean2));}
else
{fprintf(output2,"%.6f\t-1\t%.6f\t-1\t", mean, mean2);}
}
if(indcount>1)
{
if(sumsq>indcount*mean*mean&&sumsq2>indcount*mean2*mean2)
{
value=(sumsq3/indcount-mean*mean2)*pow(sumsq/indcount-mean*mean,-.5)*pow(sumsq2/indcount-mean2*mean2,-.5);
if(value>.99995||value<-0.99995){fprintf(output2, "%d\n",(value>0)-(value<0));}
else{fprintf(output2, "%.4f\n", value);}
}
else{fprintf(output2, "0\n");}
}
else{fprintf(output2, "NA\n");}
}	//end of j2 loop
}	//end of compare
}}	//end of bitlength2>0 and k loop

if(encoding!=1)	//change coding - and maybe also alleles
{change_coding(data, al1+bitstart, al2+bitstart, num_samples_use, bitlength, encoding, missingvalue);}

if(mode!=190)	//get stats and save
{
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 0);

if(minprob==0)	//info scores are not valid
{
for(j=0;j<bitlength;j++){infos[bitstart+j]=1;}
}

if(passqc==1)	//perform QC
{
for(j=0;j<bitlength;j++)
{
maf=centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]);
if(minmaf!=-9999&&maf<minmaf){nums[1][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(maxmaf!=-9999&&maf>maxmaf){nums[2][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(minvar!=-9999&&sqdevs[bitstart+j]<minvar){nums[3][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(minobs!=-9999&&rates[bitstart+j]<minobs){nums[4][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(mininfo!=-9999&&infos[bitstart+j]<mininfo){nums[5][bitstart+j]=1;nums[0][bitstart+j]=1;}
}
}

//save predictor stats (allowing for QC)
for(j=0;j<bitlength;j++)
{
if(nums[0][bitstart+j]==0)
{
fprintf(output2, "%s\t%s\t%s\t", preds[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
if(rates[bitstart+j]>0)
{
mean=centres[bitstart+j];
maf=mean/2+(mean>1)*(1-mean);
fprintf(output2,"%.6f\t",mean);
if(nonsnp==0){fprintf(output2, "%.6f\t", maf);}
else{fprintf(output2, "NA\t");}
}
else{fprintf(output2, "NA\tNA\t");}
if(rates[bitstart+j]==0||rates[bitstart+j]==1){fprintf(output2, "%d\t", (int)rates[bitstart+j]);}
else{fprintf(output2, "%.6f\t", rates[bitstart+j]);}
if(genprobs<2){fprintf(output2, "NA\n");}
else{fprintf(output2, "\t%.4f\n", infos[bitstart+j]);}
}
}

//save data (allowing for QC)
for(j=0;j<bitlength;j++)
{
if(nums[0][bitstart+j]==0)
{
#include "savedata.c"
}
}
}
}	//end of bit loop
printf("\n");

fclose(output2);
if(mode!=190){fclose(output3);}
for(k=0;k<num_files;k++)
{
if(binary==0){gzclose(Xdatainputgz[k]);}
}

////////

if(mode!=190)	//write bimfile, samplefile/famfile and .qc
{
count=0;for(j=0;j<num_preds_use;j++){count+=(nums[0][j]==0);}
if(count==0)
{
sprintf(cmd,"rm %s", filename3);
system(cmd);
printf("Error, no predictors passed quality control, so no data were saved\n\n");

if((output=fopen(filename,"w"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fprintf(output, "Error, all %d chunks were read, but no predictors passed quality control (so no files were written)\n\n", bittotal);
fclose(output);
exit(1);
}

sprintf(filename4,"%s.bim", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(j=0;j<num_preds_use;j++)
{
if(nums[0][j]==0)
{
fprintf(output4, "%d %s ", chr[j], preds[j]);
if(cm[j]==0){fprintf(output4, "0 ");}
else{fprintf(output4, "%.6f ", cm[j]);}
fprintf(output4, "%ld %s %s\n", (long int)bp[j], along1[j], along2[j]);
}
}
fclose(output4);

if(mode==185)	//sample
{
sprintf(filename5,"%s.sample", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "ID_1 ID_2 missing sex\n");
for(i=0;i<num_samples_use;i++)
{fprintf(output5, "%s %s 0 %s\n", ids1[i], ids2[i], schar[i]);}
fclose(output5);
}
else
{
sprintf(filename5,"%s.fam", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(i=0;i<num_samples_use;i++)
{fprintf(output5, "%s %s %s %s %s %s\n", ids1[i], ids2[i], pid[i], mid[i], schar[i], pchar[i]);}
fclose(output5);
}

if(passqc==1)	//save qc info
{
sprintf(filename6,"%s.qc", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6, "Predictor\tFail_Any\tFail_Min_MAF\tFail_Max_MAF\tFail_Min_Var\tFail_Min_Obs\tFail_Min_Info\n");
for(j=0;j<num_preds_use;j++)
{
fprintf(output6, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", preds[j], nums[0][j], nums[1][j], nums[2][j], nums[3][j], nums[4][j], nums[5][j]);}
fclose(output6);
}
}

if(mode!=190)
{
printf("Data for %d samples and %d predictors saved in %s, %s and %s, with statistics in %s", num_samples_use, count, filename3, filename4, filename5, filename2);
if(passqc==1){printf("; %s indicates which predictors failed quality control",filename6);}
printf("\n\n");
}
else{printf("Correlations saved in %s\n\n", filename2);}

if(mode==185)
{
sprintf(cmd, "gzip -f %s", filename3);
printf("Trying to gzip %s; if this fails to complete, run the command\n%s\n\n", filename3, cmd);
system(cmd);
}
}	//end of qmerge=0

//now free loads of stuff (have already freed quite a few details)

for(k=0;k<num_files;k++){free(Xks[k]);free(Xks2[k]);}
free(Xnall);free(Xnuse);free(Xks);free(Xks2);
for(k=0;k<num_files;k++)
{
for(j=0;j<XNall[k];j++){free(Xalong1[k][j]);free(Xalong2[k][j]);}
free(Xalong1[k]);free(Xalong2[k]);free(Xkp[k]);free(Xkp2[k]);
}
free(XNall);free(Xalong1);free(Xalong2);free(XNuse);free(Xkp);free(Xkp2);

for(i=0;i<num_samples_use;i++){free(ids1[i]);free(ids2[i]);free(ids3[i]);}
free(ids1);free(ids2);free(ids3);
for(i=0;i<num_samples_use;i++){free(pid[i]);free(mid[i]);free(schar[i]);free(pchar[i]);}
free(pid);free(mid);free(schar);free(pchar);
for(j=0;j<num_preds_use;j++){free(preds[j]);free(along1[j]);free(along2[j]);}free(preds);free(along1);free(along2);
free(chr);free(cm);free(bp);free(al1);free(al2);

if(qmerge==0)
{
free(data);free(data2);
free(Xdatainputgz);free(Xcurrent);
if(genprobs>=2&&mode!=190){free(ps[0]);free(ps[1]);free(ps);free(ps2[0]);free(ps2[1]);free(ps2);}
free(centres);free(mults);free(sqdevs);free(rates);free(infos);
for(k=0;k<6;k++){free(nums[k]);}free(nums);
}

///////////////////////////

