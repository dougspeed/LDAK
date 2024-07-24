/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Merge (or compute correlations between) datasets

///////////////////////////

//deal with samples
Xnall=malloc(sizeof(int)*num_files);
Xids1=malloc(sizeof(char**)*num_files);
Xids2=malloc(sizeof(char**)*num_files);
Xids3=malloc(sizeof(char**)*num_files);
Xnuse=malloc(sizeof(int)*num_files);
Xks=malloc(sizeof(int*)*num_files);
Xks2=malloc(sizeof(int*)*num_files);

//first read in all ids
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

if(quickmerge==1)
{
//check that all sample ids are the same size and in the same order
flag=1;
for(k=1;k<num_files;k++)
{
if(Xnall[k]!=Xnall[0]){flag=0;printf("here\n");break;}
for(i=0;i<Xnall[k];i++)
{
if(strcmp(Xids1[k][i],Xids1[0][i])!=0||strcmp(Xids2[k][i],Xids2[0][i])!=0){flag=0;printf("%d\n", i+1);break;}
}
}

if(flag==0)
{printf("Error, you can only use \"--quick-merge YES\" when all datasets contain the same samples in the same order\n\n");exit(1);}

//set ids3 (keeping original order)
num_samples=Xnall[0];
allids3=malloc(sizeof(char*)*num_samples);
for(i=0;i<num_samples;i++){copy_string(allids3,i,Xids3[0][i]);}
}
else	//must set num_samples and allids3 to intersect / union
{
num_samples=Xnall[0];
allids3=malloc(sizeof(char*)*num_samples);
for(i=0;i<num_samples;i++){copy_string(allids3,i,Xids3[0][i]);}
qsort(allids3,num_samples,sizeof(char *), compare_string);

for(k=1;k<num_files;k++)
{
count=num_samples;
wantids=malloc(sizeof(char*)*count);
for(i=0;i<count;i++){copy_string(wantids,i,allids3[i]);free(allids3[i]);}free(allids3);

if(comsamps==1)	//get intersect
{
allids3=malloc(sizeof(char*)*count);
num_samples=inter_ids(wantids, count, Xids3[k], Xnall[k], allids3, 2);
if(num_samples==0)
{printf("Error, there are no samples common to the first %d datasets\n\n", k+1);exit(1);}
}
else	//get union
{
allids3=malloc(sizeof(char*)*(count+Xnall[k]));
num_samples=uni_ids(wantids, count, Xids3[k], Xnall[k], allids3, 2+(k==1));
}

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
}	//end of k loop

if(num_files>1)
{
if(comsamps==1){printf("\nThere are %d samples common to the %d datasets\n", num_samples, num_files);}
else{printf("\nThere are %d samples across the %d datasets\n", num_samples, num_files);}
}
}	//end of quickmerge=0

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
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, NULL, NULL, 2);
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
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, NULL, NULL, 2);
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
usedids=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){usedids[i]=0;}

for(k=0;k<num_files;k++)
{
Xks[k]=malloc(sizeof(int)*Xnall[k]);
Xks2[k]=malloc(sizeof(int)*Xnall[k]);
if(quickmerge==1)
{
Xnuse[k]=num_samples_use;
for(i=0;i<num_samples_use;i++){Xks[k][i]=i;Xks2[k][i]=i;}
}
else
{
Xnuse[k]=find_strings(Xids3[k], Xnall[k], ids3, num_samples_use, Xks[k], Xks2[k], famstems[k], NULL, NULL, NULL, 1);
if(Xnuse[k]==0){printf("Error, after filtering samples, none remain in Dataset %d\n\n", k+1);exit(1);}
}
for(i=0;i<Xnuse[k];i++)
{
i2=Xks2[k][i];
if(usedids[i2]==0)	//have not already got details for this sample
{
copy_string(ids1,i2,Xids1[k][Xks[k][i]]);
copy_string(ids2,i2,Xids2[k][Xks[k][i]]);
usedids[i2]=1;
}
else	//in theory, there can be mismatches
{
if(strcmp(Xids1[k][Xks[k][i]],ids1[i2])!=0||strcmp(Xids2[k][Xks[k][i]],ids2[i2])!=0)
{printf("Error matching sample IDs (%s %s and %s %s), please tell Doug\n\n", Xids1[k][Xks[k][i]], ids1[i2], Xids2[k][Xks[k][i]], ids2[i2]);exit(1);}
}
}
}

count=0;for(i=0;i<num_samples_use;i++){count+=usedids[i];}
if(count!=num_samples_use){printf("Doug error 43566FH %d and %d\n\n", count, num_samples_use);exit(1);}

free(usedids);

//now get additional details
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

if(famhead==0)	//read in values (if using fam files with at least 6 columns)
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

//reset for those which failed
for(i=0;i<num_samples_use;i++)
{
if(erra==1){strcpy(pid[i],"0");}
if(errb==1){strcpy(mid[i],"0");}
if(errc==1){strcpy(schar[i],"0");}
if(errd==1){strcpy(pchar[i],"NA");}
}
}	//end of k loop

for(k=0;k<num_files;k++)
{
for(i=0;i<Xnall[k];i++){free(Xids1[k][i]);free(Xids2[k][i]);free(Xids3[k][i]);}
free(Xids1[k]);free(Xids2[k]);free(Xids3[k]);
}
free(Xids1);free(Xids2);free(Xids3);
for(i=0;i<num_samples;i++){free(allids3[i]);}free(allids3);

///////////////////////////

//now predictors - if filtering, will do so on-the-fly
XNall=malloc(sizeof(int*)*num_files);
Xchr=malloc(sizeof(int*)*num_files);
Xcm=malloc(sizeof(double*)*num_files);
Xbp=malloc(sizeof(double*)*num_files);
Xpreds=malloc(sizeof(char**)*num_files);
Xal1=malloc(sizeof(char*)*num_files);
Xal2=malloc(sizeof(char*)*num_files);
XNuse=malloc(sizeof(int)*num_files);
Xkp=malloc(sizeof(int*)*num_files);
Xkp2=malloc(sizeof(int*)*num_files);

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
XNall[k]=countrows(bimstems[k]);
printf("Reading details for %d predictors from %s\n", XNall[k], bimstems[k]);
Xchr[k]=malloc(sizeof(int)*XNall[k]);
Xcm[k]=malloc(sizeof(double)*XNall[k]);
Xbp[k]=malloc(sizeof(double)*XNall[k]);
Xpreds[k]=malloc(sizeof(char*)*XNall[k]);
Xal1[k]=malloc(sizeof(char)*XNall[k]);
Xal2[k]=malloc(sizeof(char)*XNall[k]);
read_bimfile(bimstems[k], Xchr[k], Xpreds[k], Xcm[k], Xbp[k], Xal1[k], Xal2[k], XNall[k], NULL, 1, -9999, (binary==1&&usenames==1), multi);
}
else	//will be reading from gen file
{
//assume at most maxpreds SNPs
printf("Reading predictor details from %s\n", bimstems[k]);
Xchr[k]=malloc(sizeof(int)*maxpreds);
Xcm[k]=malloc(sizeof(double)*maxpreds);
Xbp[k]=malloc(sizeof(double)*maxpreds);
Xpreds[k]=malloc(sizeof(char*)*maxpreds);
Xal1[k]=malloc(sizeof(char)*maxpreds);
Xal2[k]=malloc(sizeof(char)*maxpreds);
for(j=0;j<maxpreds;j++){Xchr[k][j]=oxchr;Xcm[k][j]=0;}

XNall[k]=read_genfile(bimstems[k], Xpreds[k], Xbp[k], Xal1[k], Xal2[k], Xnall[k], genskip, genheaders, genprobs, maxpreds, 1);
if(XNall[k]>maxpreds)	//then file larger than allocated, will have to read again
{
for(j=0;j<maxpreds;j++){free(Xpreds[k][j]);}free(Xpreds[k]);
free(Xchr[k]);free(Xbp[k]);free(Xal1[k]);free(Xal2[k]);
Xchr[k]=malloc(sizeof(int)*XNall[k]);
Xcm[k]=malloc(sizeof(double)*XNall[k]);
Xbp[k]=malloc(sizeof(double)*XNall[k]);
Xpreds[k]=malloc(sizeof(char*)*XNall[k]);
Xal1[k]=malloc(sizeof(char)*XNall[k]);
Xal2[k]=malloc(sizeof(char)*XNall[k]);
for(j=0;j<XNall[k];j++){Xchr[k][j]=oxchr;Xcm[k][j]=0;}
XNall[k]=read_genfile(bimstems[k], Xpreds[k], Xbp[k], Xal1[k], Xal2[k], Xnall[k], genskip, genheaders, genprobs, XNall[k], 1);
}
if(XNall[k]>=20000){printf("%s contains genotypes for %d SNPs\n", bimstems[k], XNall[k]);}
}

////////

//check for non-positives
count3=0;for(j=0;j<XNall[k];j++){count3+=(Xbp[k][j]<=0);}
if(count3==XNall[k]){printf("Error, no predictors remain in  %s\n\n", bimstems[k]);exit(1);}

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


//check for duplicates
predorder=malloc(sizeof(int)*XNall[k]);
check_dups(Xpreds[k],XNall[k],NULL,predorder,2);

found=0;
while(Xbp[k][predorder[found]]<=0){found++;}

wcount=0;
for(j=found+1;j<XNall[k];j++)
{
if(Xbp[k][predorder[j]]>0)
{
if(strcmp(Xpreds[k][predorder[found]],Xpreds[k][predorder[j]])==0)	//matches
{
if(exsame==0)
{
printf("Error, the predictor name %s appears twice in %s", Xpreds[k][predorder[j]], bimstems[k]);
if(usenames==0){printf("; to ignore these use \"--exclude-same YES\"");}
printf("\n\n");
exit(1);
}
if(wcount<5){printf("Warning, at least two predictors in %s have the same name (%s); these will be ignored\n", bimstems[k], Xpreds[k][predorder[j]]);}
Xbp[k][predorder[found]]=-1;
Xbp[k][predorder[j]]=-1;
wcount++;

//if three or more with same name, don't want to keep warning
while(j<XNall[k]-1)
{
j++;
if(Xbp[k][predorder[j]]>0)
{
if(strcmp(Xpreds[k][predorder[found]],Xpreds[k][predorder[j]])!=0){break;}
Xbp[k][predorder[j]]=-1;
}
}
}
found=j;
}}
if(wcount>5){printf("In total, %d predictor names appeared twice\n", wcount);}
if(wcount>0){printf("\n");}

if(count>0||count2>0)	//take care of extract, onesnp and exclude - exclude by giving basepair -1
{
if(count>0)	//extract
{
indexer=malloc(sizeof(int)*XNall[k]);
usedpreds=malloc(sizeof(int)*XNall[k]);

count3=find_strings(Xpreds[k], XNall[k], wantpreds, count, indexer, NULL, NULL, NULL, predorder, NULL, 3);
if(count3==0)
{
if(strcmp(bpredfile,"blank")!=0)
{printf("Error, none of the predictors in %s are in %s\n", bpredfile, bimstems[k]);exit(1);}
else
{printf("Error, Predictor %s is not in %s\n", wantpreds[0], bimstems[k]);exit(1);}
}
if(count3<count)
{printf("Warning, only %d of the predictors in %s are in %s\n", count3, bpredfile, bimstems[k]);}

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
if(count3==0)
{printf("Warning, none of the predictors in %s are in %s\n", cpredfile, bimstems[k]);}
if(count3>0&&count3<count2)
{printf("Warning, only %d of the predictors in %s are in %s\n", count3, cpredfile, bimstems[k]);}
for(j=0;j<count3;j++){Xbp[k][indexer[j]]=-1;}
free(indexer);
}

count3=0;for(j=0;j<XNall[k];j++){count3+=(Xbp[k][j]>0);}
if(count3==0)
{printf("Error, after filtering predictors, no predictors remain in %s\n\n", bimstems[k]);exit(1);}
}

free(predorder);
}	//end of k loop
printf("\n");

if(count>0)
{
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);
}
if(count2>0)
{
for(j=0;j<count2;j++){free(wantpreds2[j]);}free(wantpreds2);
}

////////

count2=0;
for(k=0;k<num_files;k++)
{
found=0;
while(Xbp[k][found]<=0){found++;}

count=0;flag=0;
for(j=found+1;j<XNall[k];j++)
{
if(Xbp[k][j]>0)
{
if(Xchr[k][j]==Xchr[k][found]&&Xbp[k][j]==Xbp[k][found])	//matches
{
if(exdups==1)	 //remove j, found stays same
{Xbp[k][j]=-1;}
count++;flag=1;
}
else	//different - will need to update found
{
if(flag==1)	//found was a duplicate
{
if(exdups==1){Xbp[k][found]=-1;}
count++;flag=0;
}
found=j;
}
}
}	//end of j loop

if(flag==1)	//last found was a duplicate
{
if(exdups==1){Xbp[k][found]=-1;}
count++;
}

if(count>0){printf("There are %d predictors with the same position in %s\n", count, bimstems[k]);}
count2+=count;

count=0;for(j=0;j<XNall[k];j++){count+=(Xbp[k][j]>0);}
if(count==0){printf("Error, after removing these, no predictors remain in %s\n\n", bimstems[k]);exit(1);}
}	//end of k loop
if(count2>0){printf("\n");}

if(usenames==0&&num_files>1&&count2>0&&exdups==0)
{
printf("Error, it is not possible to merge based on position when multiple predictors have the same position\n");
printf("You should either use \"--extract\" and/or \"--exclude\" to (selectively) filter predictors, use \"--exclude-dups YES\" to remove all predictors with the same position, or use \"--use-names YES\" to merge based on predictor names\n\n");
exit(1);
}

////////

//now set num_preds_use and chr/preds/bp/al1/al2 and get indexes - will sort cm later
if(usenames==0)	//sort based on position
{
//start by setting details based on first file
num_preds_use=0;for(j=0;j<XNall[0];j++){num_preds_use+=(Xbp[0][j]>0);}
chr=malloc(sizeof(int)*num_preds_use);
bp=malloc(sizeof(double)*num_preds_use);
preds=malloc(sizeof(char*)*num_preds_use);
al1=malloc(sizeof(char)*num_preds_use);
al2=malloc(sizeof(char)*num_preds_use);

count=0;
for(j=0;j<XNall[0];j++)
{
if(Xbp[0][j]>0)
{
chr[count]=Xchr[0][j];bp[count]=Xbp[0][j];
copy_string(preds,count,Xpreds[0][j]);
al1[count]=Xal1[0][j];al2[count]=Xal2[0][j];
count++;
}
}

for(k=1;k<num_files;k++)	//search remaining datasets for intersect / union
{
//save predictor details
count=num_preds_use;
allchr=malloc(sizeof(int)*count);
allbp=malloc(sizeof(double)*count);
allpreds=malloc(sizeof(char*)*count);
allal1=malloc(sizeof(char)*count);
allal2=malloc(sizeof(char)*count);

for(j=0;j<count;j++)
{
allchr[j]=chr[j];allbp[j]=bp[j];
copy_string(allpreds,j,preds[j]);
allal1[j]=al1[j];allal2[j]=al2[j];
}
for(j=0;j<num_preds_use;j++){free(preds[j]);}free(preds);
free(chr);free(bp);free(al1);free(al2);

if(compreds==1)	//get intersect
{
chr=malloc(sizeof(int)*count);
bp=malloc(sizeof(double)*count);
preds=malloc(sizeof(char*)*count);
al1=malloc(sizeof(char)*count);
al2=malloc(sizeof(char)*count);
num_preds_use=inter_preds(allchr, allbp, allpreds, allal1, allal2, count, Xchr[k], Xbp[k], Xpreds[k], Xal1[k], Xal2[k], XNall[k], chr, bp, preds, al1, al2, k);
if(num_preds_use==0)
{printf("Error, there are no predictors common to the first %d datasets\n\n", k+1);exit(1);}
}
else	//get union
{
chr=malloc(sizeof(int)*(count+XNall[k]));
bp=malloc(sizeof(double)*(count+XNall[k]));
preds=malloc(sizeof(char*)*(count+XNall[k]));
al1=malloc(sizeof(char)*(count+XNall[k]));
al2=malloc(sizeof(char)*(count+XNall[k]));
num_preds_use=uni_preds(allchr, allbp, allpreds, allal1, allal2, count, Xchr[k], Xbp[k], Xpreds[k], Xal1[k], Xal2[k], XNall[k], chr, bp, preds, al1, al2, k);
}

for(j=0;j<count;j++){free(allpreds[j]);}free(allpreds);
free(allchr);free(allbp);free(allal1);free(allal2);
}	//end of k loop

if(compreds==1){printf("There are %d predictors common to the %d datasets\n\n", num_preds_use, num_files);}
else{printf("There are %d predictors across the %d datasets\n\n", num_preds_use, num_files);}

////////

//get indexes
for(k=0;k<num_files;k++)
{
Xkp[k]=malloc(sizeof(int)*XNall[k]);
Xkp2[k]=malloc(sizeof(int)*XNall[k]);

count=0;found=0;
for(j=0;j<XNall[k];j++)
{
if(Xbp[k][j]>0)	//see if j being used
{
while(chr[count]<Xchr[k][j]||(chr[count]==Xchr[k][j]&&bp[count]<Xbp[k][j]))
{
if(count==num_preds_use-1){break;}
count++;
}
if(chr[count]==Xchr[k][j]&&bp[count]==Xbp[k][j])
{Xkp[k][found]=j;Xkp2[k][found]=count;found++;count++;}
if(count==num_preds_use){break;}
}}
XNuse[k]=found;
if(XNuse[k]==0){printf("Error, will not be using any predictors from Dataset %d\n\n", k+1);exit(1);}
}	//end of k loop
}	//end of usenames=0
else	//sort based on name
{
//first work out which predictors to use

if(compreds==1)	//want intersect
{
//set names based on first file, then sort them
num_preds_use=0;for(j=0;j<XNall[0];j++){num_preds_use+=(Xbp[0][j]>0);}
preds=malloc(sizeof(char*)*num_preds_use);
count=0;
for(j=0;j<XNall[0];j++)
{
if(Xbp[0][j]>0){copy_string(preds,count,Xpreds[0][j]);count++;}
}
qsort(preds,num_preds_use,sizeof(char *), compare_string);

for(k=1;k<num_files;k++)	//search remaining datasets
{
indexer=malloc(sizeof(int)*num_preds_use);
indexer2=malloc(sizeof(int)*num_preds_use);
count=find_strings(preds, num_preds_use, Xpreds[k], XNall[k], indexer, indexer2, NULL, NULL, NULL, NULL, 2);
if(count==0)
{printf("Error, there are no predictors common to the first %d datasets\n\n", k+1);exit(1);}

//remove those not found or with negative bps (only need to check bps in Dataset k)
count2=0;
for(j=0;j<count;j++)
{
j2=indexer[j];
j3=indexer2[j];
if(Xbp[k][j3]>0)
{
if(count2!=j2)	//squeeze down
{free(preds[count2]);copy_string(preds,count2,preds[j2]);}
count2++;
}
}

if(count2==0)
{printf("Error, after filtering predictors, there are none common to the first %d datasets\n\n", k+1);exit(1);}

for(j=count2;j<num_preds_use;j++){free(preds[j]);}
num_preds_use=count2;
free(indexer);free(indexer2);
}

printf("There are %d predictors common to the %d datasets\n\n", num_preds_use, num_files);
}
else	//want union
{
//read in all names
count=0;
for(k=0;k<num_files;k++)
{
for(j=0;j<XNall[k];j++){count+=(Xbp[k][j]>0);}
}

allpreds=malloc(sizeof(char*)*count);
count=0;
for(k=0;k<num_files;k++)
{
for(j=0;j<XNall[k];j++)
{
if(Xbp[k][j]>0){copy_string(allpreds,count,Xpreds[k][j]);count++;}
}
}

//find unique names
printf("Sorting names of %d predictors; this can take a while\n", count);
qsort(allpreds,count,sizeof(char *), compare_string);

num_preds_use=1;
for(j=1;j<count;j++){num_preds_use+=(strcmp(allpreds[j],allpreds[j-1])!=0);}

preds=malloc(sizeof(char*)*num_preds_use);
copy_string(preds,0,allpreds[0]);
num_preds_use=1;
for(j=1;j<count;j++)
{
if(strcmp(allpreds[j],allpreds[j-1])!=0)	//it's new
{copy_string(preds,num_preds_use,allpreds[j]);num_preds_use++;}
}

for(j=0;j<count;j++){free(allpreds[j]);}free(allpreds);

printf("In total, there are %d unique predictor names\n\n", num_preds_use);
}

////////

//get details for unique names, and also indexes
chr=malloc(sizeof(int)*num_preds_use);
bp=malloc(sizeof(double)*num_preds_use);
al1=malloc(sizeof(char)*num_preds_use);
al2=malloc(sizeof(char)*num_preds_use);
for(k=0;k<num_files;k++)
{
Xkp[k]=malloc(sizeof(int)*XNall[k]);
Xkp2[k]=malloc(sizeof(int)*XNall[k]);
}

usedpreds=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){usedpreds[j]=0;}

for(k=0;k<num_files;k++)
{
XNuse[k]=find_strings(Xpreds[k], XNall[k], preds, num_preds_use, Xkp[k], Xkp2[k], NULL, NULL, NULL, NULL, 1);

for(j=0;j<XNuse[k];j++)
{
j2=Xkp[k][j];
j3=Xkp2[k][j];
if(usedpreds[j3]==0)	//it's new
{
chr[j3]=Xchr[k][j2];bp[j3]=Xbp[k][j2];
al1[j3]=Xal1[k][j2];al2[j3]=Xal2[k][j2];
usedpreds[j3]=k+1;
}
else	//check it matches
{
if(chr[j3]!=Xchr[k][j2]||bp[j3]!=Xbp[k][j2])
{printf("Error, the position of Predictor %s in Dataset %d (%d:%.2f) does not match that in Dataset %d (%d:%.2f)\n\n", preds[j3], usedpreds[j3], chr[j3], bp[j3], k+1, Xchr[k][j2], Xbp[k][j2]);exit(1);}

if((al1[j3]!=Xal1[k][j2]&&al1[j3]!=Xal2[k][j2])||(al2[j3]!=Xal1[k][j2]&&al2[j3]!=Xal2[k][j2]))
{printf("Error, the alleles of Predictor %s in Dataset %d (%c and %c) do not match those in Dataset %d (%c and %c)\n\n", preds[j3], usedpreds[j3], al1[j3], al2[j3], k+1, Xal1[k][j2], Xal2[k][j2]);exit(1);}
}
}
}

free(usedpreds);

//copy over predictor details
allchr=malloc(sizeof(int)*num_preds_use);
allbp=malloc(sizeof(double)*num_preds_use);
allpreds=malloc(sizeof(char*)*num_preds_use);
allal1=malloc(sizeof(char)*num_preds_use);
allal2=malloc(sizeof(char)*num_preds_use);

for(j=0;j<num_preds_use;j++)
{
allchr[j]=chr[j];allbp[j]=bp[j];
copy_string(allpreds,j,preds[j]);
allal1[j]=al1[j];allal2[j]=al2[j];
}
for(j=0;j<num_preds_use;j++){free(preds[j]);}

//now refill, sorted by bp
dptrs=malloc(sizeof(struct sorting_double)*num_preds_use);
usedpreds=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){usedpreds[j]=-1;}

//find max chromosome number
count=1;
for(j=0;j<num_preds_use;j++)
{
if(allchr[j]>=count){count=allchr[j]+1;}
}

//loop through chromosomes
count3=0;
for(k=0;k<count;k++)
{
count2=0;
for(j=0;j<num_preds_use;j++)
{
if(allchr[j]==k){dptrs[count2].value=allbp[j];dptrs[count2].index=j;count2++;}
}

if(count2>0)	//there are predictors on this chromosome, so sort and load up
{
qsort(dptrs, count2, sizeof(struct sorting_double), compare_sorting_double);

for(j=0;j<count2;j++)
{
j2=dptrs[j].index;
chr[count3]=allchr[j2];bp[count3]=allbp[j2];
copy_string(preds,count3,allpreds[j2]);
al1[count3]=allal1[j2];al2[count3]=allal2[j2];
usedpreds[j2]=count3;
count3++;
}
}
}

//update Xkp2
for(k=0;k<num_files;k++)
{
for(j=0;j<XNuse[k];j++){Xkp2[k][j]=usedpreds[Xkp2[k][j]];}
}

free(allchr);free(allbp);free(allal1);free(allal2);
for(j=0;j<num_preds_use;j++){free(allpreds[j]);}free(allpreds);
free(dptrs);free(usedpreds);

//finally, must make sure Xkp2 is in order
indexer=malloc(sizeof(int)*num_preds_use);
indexer2=malloc(sizeof(int)*num_preds_use);
usedpreds=malloc(sizeof(int)*num_preds_use);

for(k=0;k<num_files;k++)
{
for(j=0;j<num_preds_use;j++){usedpreds[j]=-1;}
for(j=0;j<XNuse[k];j++){indexer[j]=Xkp[k][j];indexer2[j]=Xkp2[k][j];usedpreds[Xkp2[k][j]]=j;}
count=0;
for(j=0;j<num_preds_use;j++)
{
if(usedpreds[j]!=-1){Xkp[k][count]=indexer[usedpreds[j]];Xkp2[k][count]=indexer2[usedpreds[j]];count++;}
}
}

free(indexer);free(indexer2);free(usedpreds);
}	//end of usenames=1

////////

//get cms
cm=malloc(sizeof(double)*num_preds_use);
for(j=0;j<num_preds_use;j++){cm[j]=-9999;}

flag=0;
for(k=0;k<num_files;k++)
{
if(flag==0)
{
for(count=0;count<XNuse[k];count++)
{
j=Xkp[k][count];j2=Xkp2[k][count];
if(cm[j2]==-9999){cm[j2]=Xcm[k][j];}
if(cm[j2]!=Xcm[k][j])
{
printf("Warning, %s has at least two genetic distances (%f and %f), so all genetic distances will be set to 0\n\n", preds[j], cm[j2], Xcm[k][j]);
for(j=0;j<num_preds_use;j++){cm[j]=0;}
flag=1;break;
}
}}
}

if(quickmerge==1)	//check each predictor in exactly one dataset
{
usedpreds=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){usedpreds[j]=0;}

for(k=0;k<num_files;k++)
{
for(j=0;j<XNuse[k];j++)
{
j2=Xkp2[k][j];
if(usedpreds[j2]!=0){printf("Error, Predictor %s is in multiple datasets; you can only use \"--quick-merge YES\" when each predictor is in only one dataset\n\n", preds[j2]);exit(1);}
usedpreds[j2]=1;
}
}

free(usedpreds);
}

for(k=0;k<num_files;k++)
{
for(j=0;j<XNall[k];j++){free(Xpreds[k][j]);}free(Xpreds[k]);
free(Xchr[k]);free(Xcm[k]);free(Xbp[k]);
}
free(Xchr);free(Xcm);free(Xbp);free(Xpreds);

///////////////////////////
//ready to look at data
///////////////////////////

if(quickmerge==1)	//do separately
{
#include "rapid.c"
}
else
{
if(bitsize>num_preds_use){bitsize=num_preds_use;}

data_warn2(bitsize,2*num_samples_use);

data=malloc(sizeof(double)*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);

Xdatainputgz=malloc(sizeof(gzFile)*num_files);
Xcurrent=malloc(sizeof(int)*num_files);
for(k=0;k<num_files;k++)
{
if(binary==0){open_datagz(Xdatainputgz+k, datastems[k], Xnall[k], genskip, genheaders, genprobs);}
Xcurrent[k]=0;
}

if(genprobs>1&&mode!=190)
{
data_warn2(bitsize,num_samples_use);
ps=malloc(sizeof(float*)*2);
ps[0]=malloc(sizeof(float)*num_samples_use*bitsize);
ps[1]=malloc(sizeof(float)*num_samples_use*bitsize);
ps2=malloc(sizeof(float*)*2);
ps2[0]=malloc(sizeof(float)*num_samples_use*bitsize);
ps2[1]=malloc(sizeof(float)*num_samples_use*bitsize);
}

//will record which are lost to qc
if(passall==0)	//doing one lot of qc
{
nums=malloc(sizeof(int*)*6);
for(k=0;k<6;k++)
{
nums[k]=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){nums[k][j]=0;}
}
}
if(passall==1)	//doing qc for each dataset
{
nums=malloc(sizeof(int*)*(1+num_files));
for(k=0;k<1+num_files;k++)
{
nums[k]=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){nums[k][j]=0;}
}
}
if(passall==-9999)	//not using qc (so all pass)
{
nums=malloc(sizeof(int*));
nums[0]=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){nums[0][j]=0;}
}

//will record whether predictor encountered before
usedpreds=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){usedpreds[j]=0;}

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

bittotal=(num_preds_use-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>num_preds_use){bitend=num_preds_use;}
bitlength=bitend-bitstart;

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

//start everything at missing
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use;i++)
{
data[(size_t)j*num_samples_use+i]=missingvalue;
if(genprobs>1&&mode!=190){ps[0][(size_t)j*num_samples_use+i]=0;ps[1][(size_t)j*num_samples_use+i]=0;}
}
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
if(genprobs>1&&mode!=190)	//read data and probabilities - this is rarely used, and must be dtype=5
{Xcurrent[k]=read_data_fly(datastems[k], dtype, data2, ps2, Xnuse[k], Xks[k], bitstart2, bitend2, Xkp[k], Xdatainputgz[k], Xcurrent[k], Xnall[k], XNall[k], genskip, genheaders, genprobs, bedbytes, missingvalue, threshold, minprob, nonsnp, 1);}
else	//read just data
{Xcurrent[k]=read_data_fly(datastems[k], dtype, data2, NULL, Xnuse[k], Xks[k], bitstart2, bitend2, Xkp[k], Xdatainputgz[k], Xcurrent[k], Xnall[k], XNall[k], genskip, genheaders, genprobs, bedbytes, missingvalue, threshold, minprob, nonsnp, maxthreads);}

if(passall==1)	//see which predictors pass qc
{
for(j2=0;j2<bitlength2;j2++)
{
j=Xkp2[k][bitstart2+j2]-bitstart;
sum=0;sumsq=0;indcount=0;
for(i2=0;i2<Xnuse[k];i2++)
{
if(data2[(size_t)j2*Xnuse[k]+i2]!=missingvalue)
{sum+=data2[(size_t)j2*Xnuse[k]+i2];sumsq+=pow(data2[(size_t)j2*Xnuse[k]+i2],2);indcount++;}
}
if(indcount>0){mean=sum/indcount;var=sumsq/indcount-pow(mean,2);}
else{mean=0;var=0;}
maf=mean/2+(mean>1)*(1-mean);	//will not use when indcount=0 (or if nonsnp=1)

if(minmaf!=-9999&&(indcount==0||maf<minmaf)){nums[1+k][bitstart+j]=1;}
if(maxmaf!=-9999&&(indcount==0||maf>maxmaf)){nums[1+k][bitstart+j]=1;}
if(minvar!=-9999&&(indcount<2||var<minvar)){nums[1+k][bitstart+j]=1;}
if(minobs!=-9999&&indcount<minobs*Xnuse[k]){nums[1+k][bitstart+j]=1;}

if(mininfo!=-9999)
{
if(compute_info(data2+(size_t)j2*Xnuse[k], ps2[0]+j2*Xnuse[k], ps2[1]+j2*Xnuse[k], NULL, Xnuse[k], missingvalue)<=mininfo)
{nums[1+k][bitstart+j]=1;}
}
nums[0][bitstart+j]+=nums[1+k][bitstart+j];
}
}

//see if any predictors need turning - cany only be case if predictors in multiple datasets
for(j2=0;j2<bitlength2;j2++)
{
j=Xkp2[k][bitstart2+j2]-bitstart;
if(Xal1[k][Xkp[k][bitstart2+j2]]!=al1[bitstart+j])	//turn
{
for(i2=0;i2<Xnuse[k];i2++)
{
if(data2[(size_t)j2*Xnuse[k]+i2]!=missingvalue)
{
data2[(size_t)j2*Xnuse[k]+i2]=2-data2[(size_t)j2*Xnuse[k]+i2];
if(genprobs>1&&mode!=190)
{ps2[0][(size_t)j2*Xnuse[k]+i2]=1-ps2[0][(size_t)j2*Xnuse[k]+i2]-ps2[1][(size_t)j2*Xnuse[k]+i2];}
}}
}
}

if(mode!=190||k==0)	//store values into data
{
for(j2=0;j2<bitlength2;j2++)
{
j=Xkp2[k][bitstart2+j2]-bitstart;

if(usedpreds[bitstart+j]==0)	//not seen in a previous dataset, so can always copy new values in
{
for(i2=0;i2<Xnuse[k];i2++){data[(size_t)j*num_samples_use+Xks2[k][i2]]=data2[(size_t)j2*Xnuse[k]+i2];}
if(genprobs>1&&mode!=190)
{
for(i2=0;i2<Xnuse[k];i2++)
{ps[0][(size_t)j*num_samples_use+Xks2[k][i2]]=ps2[0][(size_t)j2*Xnuse[k]+i2];
ps[1][(size_t)j*num_samples_use+Xks2[k][i2]]=ps2[1][(size_t)j2*Xnuse[k]+i2];}
}
usedpreds[bitstart+j]=1;
}
else	//predictor seen before, so must check new values consistent
{
for(i2=0;i2<Xnuse[k];i2++)
{
i=Xks2[k][i2];
if(data2[(size_t)j2*Xnuse[k]+i2]!=missingvalue)	//have a new value
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)	//have a value already - check new one matches
{
if(fabs(data2[(size_t)j2*Xnuse[k]+i2]-data[(size_t)j*num_samples_use+i])>.01)
{printf("Error, Sample %s %s %d has conflicting values (%.6f and %.6f) for %s %c %c %c %c %d %d\n\n", ids1[i], ids2[i], i2, data[(size_t)j*num_samples_use+i], data2[(size_t)j2*Xnuse[k]+i2], preds[bitstart+j], al1[bitstart+j], al2[bitstart+j], Xal1[k][bitstart2+j2], Xal2[k][bitstart2+j2], bitstart+j, bitstart2+j2);exit(1);}
}
else	//no value yet - copy in
{
data[(size_t)j*num_samples_use+i]=data2[(size_t)j2*Xnuse[k]+i2];
if(genprobs>1&&mode!=190)
{
ps[0][(size_t)j*num_samples_use+i]=ps2[0][(size_t)j2*Xnuse[k]+i2];
ps[1][(size_t)j*num_samples_use+i]=ps2[1][(size_t)j2*Xnuse[k]+i2];
}
}
}}
}
}	//end of j2 loop
}	//end of storing	
else	//must be 190 and k=1 so compare and save - will have Xnuse[k]=num_samples_use
{
for(j2=0;j2<bitlength2;j2++)
{
j=Xkp2[k][bitstart2+j2]-bitstart;

sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;indcount=0;
for(i2=0;i2<Xnuse[k];i2++)
{
i=Xks2[k][i2];
if(data[(size_t)j*num_samples_use+i]!=missingvalue&&data2[(size_t)j2*Xnuse[k]+i2]!=missingvalue)
{
sum+=data[(size_t)j*num_samples_use+i];sum2+=data2[(size_t)j2*Xnuse[k]+i2];
sumsq+=pow(data[(size_t)j*num_samples_use+i],2);sumsq2+=pow(data2[(size_t)j2*Xnuse[k]+i2],2);
sumsq3+=data[(size_t)j*num_samples_use+i]*data2[(size_t)j2*Xnuse[k]+i2];
indcount++;
}}

fprintf(output2, "%s\t%c\t%c\t%d\t", preds[bitstart+j], al1[bitstart+j], al2[bitstart+j], indcount);
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

if(encoding!=1)	//changing coding - and maybe change alleles
{
for(j=0;j<bitlength;j++)
{
if(encoding==2)	//dominant - switch 1s for 2s
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]==1){data[(size_t)j*num_samples_use+i]=2;}
}
}
if(encoding==3)	//recessive - switch 1s for 0s
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]==1){data[(size_t)j*num_samples_use+i]=0;}
}
}
if(encoding==4)	//heterogeneous - 0/1/2 go to 0/2/0
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{data[(size_t)j*num_samples_use+i]=2*(data[(size_t)j*num_samples_use+i]==1);}
}
al1[bitstart+j]='M';al2[bitstart+j]='O';
}
if(encoding==5)	//ensure minor allele is a1 allele
{
sum=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{sum+=data[(size_t)j*num_samples_use+i];indcount++;}
}
if(sum/indcount>1)	//switch
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{data[(size_t)j*num_samples_use+i]=2-data[(size_t)j*num_samples_use+i];}
}
readchar=al1[bitstart+j];
al1[bitstart+j]=al2[bitstart+j];
al2[bitstart+j]=readchar;
}
}
if(encoding==6)	//missing or not missing
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)j*num_samples_use+i]=(data[(size_t)j*num_samples_use+i]==missingvalue);}
al1[bitstart+j]='M';al2[bitstart+j]='P';
}
}
}	//end of encoding

if(mode!=190)	//have stored, so now get stats and save
{
for(j=0;j<bitlength;j++)
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{sum+=data[(size_t)j*num_samples_use+i];sumsq+=pow(data[(size_t)j*num_samples_use+i],2);indcount++;}
}

if(indcount>0){mean=sum/indcount;var=sumsq/indcount-pow(mean,2);}
else{mean=0;var=0;}
maf=mean/2+(mean>1)*(1-mean);	//will not use when indcount=0 (or if nonsnp=1)
if(genprobs>1){value=compute_info(data+(size_t)j*num_samples_use, ps[0]+j*num_samples_use, ps[1]+j*num_samples_use, NULL, num_samples_use, missingvalue);}

if(passall==0)	//will be filtering
{
if(minmaf!=-9999&&(indcount==0||maf<minmaf)){nums[1][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(maxmaf!=-9999&&(indcount==0||maf>maxmaf)){nums[2][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(minvar!=-9999&&(indcount<2||var<minvar)){nums[3][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(minobs!=-9999&&indcount<minobs*num_samples_use){nums[4][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(mininfo!=-9999&&value<=mininfo){nums[5][bitstart+j]=1;nums[0][bitstart+j]=1;}
}

if(nums[0][bitstart+j]==0)	//will be keeping this predictor
{
fprintf(output2, "%s\t%c\t%c\t%.6f", preds[bitstart+j], al1[bitstart+j], al2[bitstart+j], mean);
if(nonsnp==0){fprintf(output2, "\t%.6f", mean/2+(mean>1)*(1-mean));}
else{fprintf(output2, "\t-1");}
if(indcount==0||indcount==num_samples_use){fprintf(output2, "\t%d", (indcount==num_samples_use));}
else{fprintf(output2, "\t%.6f", (double)indcount/num_samples_use);}
if(genprobs>1){fprintf(output2, "\t%.4f\n", value);}
else{fprintf(output2, "\t-1\n");}

if(mode==184)	//work out min and max and print header columns
{
if(nonsnp==0){minfloat=0;maxfloat=2;}
else
{
minfloat=missingvalue;maxfloat=missingvalue;
for(i=0;i<num_samples_use;i++)
{
value=data[(size_t)j*num_samples_use+i];
if(value!=missingvalue)
{
if(minfloat==missingvalue){minfloat=value;maxfloat=value;}
if(value<minfloat){minfloat=value;}
if(value>maxfloat){maxfloat=value;}
}
}
if(minfloat==missingvalue){minfloat=0;maxfloat=2;}
if(minfloat==maxfloat){maxfloat++;}
}
if(speedlong==0){writefloat=1;}
else{writefloat=2;}
fwrite(&writefloat, sizeof(float), 1, output3);
fwrite(&minfloat, sizeof(float), 1, output3);
fwrite(&maxfloat, sizeof(float), 1, output3);
writefloat=0;
for(k=0;k<13;k++){fwrite(&writefloat, sizeof(float), 1, output3);}
}

if(mode==185)	//print headers
{fprintf(output3, "%d %s %ld %c %c ", chr[bitstart+j], preds[bitstart+j], (long int)bp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);}

for(i=0;i<num_samples_use;i++)
{
value=data[(size_t)j*num_samples_use+i];

if(mode==181)
{
if(i%4==0){onechar=0;}
if(value==0){onechar+=(3<<(2*(i%4)));}
if(value==1){onechar+=(2<<(2*(i%4)));}
//if(value==2){onechar+=(0<<(2*(i%4)));}
if(value==missingvalue){onechar+=(1<<(2*(i%4)));}
if(i%4==3||i==num_samples_use-1){fwrite(&onechar, sizeof(unsigned char), 1, output3);}
}
if(mode==182)
{
if(value!=missingvalue)
{
if(fabs(value-round(value))<0.00005){fprintf(output3,"%d ", (int)round(value));}
else{fprintf(output3,"%.4f ", value);}
}
else
{fprintf(output3,"NA ");}
}
if(mode==183)
{
writefloat=(float)value;
fwrite(&writefloat, sizeof(float), 1, output3);
}
if(mode==184)
{
if(speedlong==0)	//use 256
{
if(value==missingvalue){onechar=255;}
else{onechar=round((value-minfloat)/(maxfloat-minfloat)*254);}
fwrite(&onechar, sizeof(unsigned char), 1, output3);
}
else	//use 65536
{
if(value==missingvalue){oneshort=65535;}
else{oneshort=round((value-minfloat)/(maxfloat-minfloat)*65534);}
fwrite(&oneshort, sizeof(unsigned short), 1, output3);
}
}
if(mode==185)
{
value=ps[0][(size_t)j*num_samples_use+i];
if(value==0||value==1){fprintf(output3,"%d ", (int)value);}
else{fprintf(output3,"%.3f ", value);}
value=ps[1][(size_t)j*num_samples_use+i];
if(value==0||value==1){fprintf(output3,"%d ", (int)value);}
else{fprintf(output3,"%.3f ", value);}
value=1-ps[0][(size_t)j*num_samples_use+i]-ps[1][(size_t)j*num_samples_use+i];
if(value==0||value==1){fprintf(output3,"%d ", (int)value);}
else{fprintf(output3,"%.3f ", value);}
}
}	//end of i loop
if(mode==182||mode==185){fprintf(output3,"\n");}
}	//end of using this predictor
}}	//end of j loop and mode!=190
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
printf("Error, no predictors passed quality control; no files were written\n\n");
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
fprintf(output4, "%ld %c %c\n", (long int)bp[j], al1[j], al2[j]);
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

if(passall==0)
{
sprintf(filename,"%s.qc", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tFail_Any\tFail_Min_MAF\tFail_Max_MAF\tFail_Min_Var\tFail_Min_Obs\tFail_Min_Info\n");
for(j=0;j<num_preds_use;j++)
{
fprintf(output, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", preds[j], nums[0][j], nums[1][j], nums[2][j], nums[3][j], nums[4][j], nums[5][j]);}
fclose(output);
}
if(passall==1)
{
sprintf(filename,"%s.qc", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tFail_Overall\t");
for(k=0;k<num_files;k++){fprintf(output, "Fail_Dataset%d\t", k+1);}
fprintf(output, "\n");
for(j=0;j<num_preds_use;j++)
{
fprintf(output, "%s\t%d\t", preds[j], nums[0][j]);
for(k=0;k<num_files;k++){fprintf(output, "%d\t", nums[1+k][j]);}
fprintf(output, "\n");
}
fclose(output);
}
}

if(mode!=190)
{
printf("Data for %d samples and %d predictors saved in %s, %s and %s, with statistics in %s", num_samples_use, count, filename3, filename4, filename5, filename2);
if(passall==0)
{printf("; %s indicates how predictors failed quality control",filename);}
if(passall==1)
{printf("; %s indicates in which datasets predictors failed quality control",filename);}
printf("\n\n");
}
else{printf("Correlations saved in %s\n\n", filename2);}

if(mode==185)
{
sprintf(cmd, "gzip -f %s", filename3);
printf("Trying to gzip %s; if this fails to complete, run the command\n%s\n\n", filename3, cmd);
system(cmd);
}
}	//end of quickmerge=0

//now free loads of stuff (have already freed quite a few details)

for(k=0;k<num_files;k++){free(Xks[k]);free(Xks2[k]);}
free(Xnall);free(Xnuse);free(Xks);free(Xks2);
for(k=0;k<num_files;k++){free(Xal1[k]);free(Xal2[k]);free(Xkp[k]);free(Xkp2[k]);}
free(XNall);free(Xal1);free(Xal2);free(Xkp);free(Xkp2);free(XNuse);

for(i=0;i<num_samples_use;i++){free(ids1[i]);free(ids2[i]);free(ids3[i]);}
free(ids1);free(ids2);free(ids3);
for(i=0;i<num_samples_use;i++){free(pid[i]);free(mid[i]);free(schar[i]);free(pchar[i]);}
free(pid);free(mid);free(schar);free(pchar);
for(j=0;j<num_preds_use;j++){free(preds[j]);}free(preds);
free(chr);free(cm);free(bp);free(al1);free(al2);

if(quickmerge==0)
{
free(data);free(data2);
free(Xdatainputgz);free(Xcurrent);
if(genprobs>1&&mode!=190){free(ps[0]);free(ps[1]);free(ps);free(ps2[0]);free(ps2[1]);free(ps2);}
if(passall==0){for(k=0;k<6;k++){free(nums[k]);}free(nums);}
if(passall==1){for(k=0;k<1+num_files;k++){free(nums[k]);}free(nums);}
if(passall==-9999){free(nums[0]);free(nums);}
free(usedpreds);
}

///////////////////////////

