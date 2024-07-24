/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Get names of samples / predictors and apply filterings - will not be here if merging

///////////////////////////

if(use_data==1||use_data==2||use_data==3||use_data==4||num_kins>0||mode==121||mode==123||mode==124||mode==129||mode==130||mode==229||mode==230||mode==177)	//will be using sample ids (and will set idsfile)
{
if(num_kins>0)	//get (sorted) intersect of kinship matrices
{
sprintf(filename,"%s.grm.id",kinstems[0]);
count=countrows(filename);
kinids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, kinids, count, NULL, 0, 0);

if(mode==113||(mode==121&&memsave==1&&hestart==1)||(mode==123&&memsave==1)||(mode==124&&memsave==1)||(mode==126&&memsave==1&&hestart==1)||(mode==133&&memsave==1&&hestart==1))	//all id files must match
{
for(k=1;k<num_kins;k++)
{
sprintf(filename2,"%s.grm.id",kinstems[k]);
count2=countrows(filename2);
if(count2!=count)
{
if(mode==113){printf("Error, %s does not match %s; each kinship matrix should consider the same samples (in the same order)\nConsider using \"--add-grm\" instead\n\n", filename, filename2);exit(1);}
if(mode==121||mode==126||mode==133){printf("Error, %s does not match %s; when using \"--memory-save YES\", each kinship matrix should correspond to the same samples (in the same order), otherwise you must add \"--he-starts NO\"\n\n", filename, filename2);exit(1);}
if(mode==123||mode==124){printf("Error, %s does not match %s; when using \"--memory-save YES\", each kinship matrix should consider the same samples (in the same order)\n\n", filename, filename2);exit(1);}
}

kinids2=malloc(sizeof(char*)*count2);
read_ids(filename2, NULL, NULL, kinids2, count2, NULL, 0, 0);

for(i=0;i<count2;i++)
{
if(strcmp(kinids[i],kinids2[i])!=0)
{
if(mode==113){printf("Error, %s does not match %s; each kinship matrix should consider the same samples (in the same order)\nConsider using \"--add-grm\" instead\n\n", filename, filename2);exit(1);}
if(mode==121||mode==126||mode==133){printf("Error, %s does not match %s; when using \"--memory-save YES\", each kinship matrix should ideally consider the same samples (in the same order), otherwise you must add \"--he-starts NO\"\n\n", filename, filename2);exit(1);}
if(mode==123||mode==124){printf("Error, %s does not match %s; when using \"--memory-save YES\", each kinship matrix should consider the same samples (in the same order)\n\n", filename, filename2);exit(1);}
}}
for(i=0;i<count2;i++){free(kinids2[i]);}free(kinids2);
}
}

qsort(kinids, count, sizeof(char *), compare_string);
total=count;

//now find overlap - will be redundant when id files match
for(k=1;k<num_kins;k++)
{
sprintf(filename2,"%s.grm.id",kinstems[k]);
count2=countrows(filename2);
kinids2=malloc(sizeof(char*)*count2);
read_ids(filename2, NULL, NULL, kinids2, count2, NULL, 0, 0);

count3=total;
kinids3=malloc(sizeof(char*)*count3);
for(i=0;i<count3;i++){kinids3[i]=kinids[i];}
total=inter_ids(kinids3, count3, kinids2, count2, kinids, 2);

for(i=0;i<count2;i++){free(kinids2[i]);}free(kinids2);
for(i=0;i<count3;i++){free(kinids3[i]);}free(kinids3);
}

if(total==0){printf("Error, there are no samples common to the %d kinship matrices\n\n", num_kins);exit(1);}
}	//end of num_kins>0

////////

//set idsfile and get num_samples

if(num_kins>0&&mode!=122&&mode!=131&&mode!=138)	//start with intersect of kinships
{
sprintf(idsfile,"%s.grm.id",kinstems[0]);
num_samples=total;
if(num_kins==1){printf("Reading IDs for %d samples from %s\n\n", num_samples, idsfile);}
else{printf("There are %d samples common to the %d kinship matrices\n\n", num_samples, num_kins);}

count=countrows(idsfile);
wantids=malloc(sizeof(char*)*count);
read_ids(idsfile, NULL, NULL, wantids, count, NULL, 0, 0);
indexer=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, kinids, num_samples, indexer, NULL, NULL, NULL, NULL, NULL, 1);
if(count2<num_samples){printf("Error 78C; please tell Doug %d %d\n\n", count2, num_samples);exit(1);}
allids1=malloc(sizeof(char*)*num_samples);
allids2=malloc(sizeof(char*)*num_samples);
allids3=malloc(sizeof(char*)*num_samples);
read_ids(filename, allids1, allids2, allids3, num_samples, indexer, 0, 0);
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);
}
else
{
if(use_data==1||use_data==2||use_data==3||use_data==4)	//start with famfile
{
strcpy(idsfile,famfile);
idshead=famhead;
num_samples=countrows(idsfile)-idshead;
printf("Reading IDs for %d samples from %s\n\n", num_samples, idsfile);
allids1=malloc(sizeof(char*)*num_samples);
allids2=malloc(sizeof(char*)*num_samples);
allids3=malloc(sizeof(char*)*num_samples);
read_ids(idsfile, allids1, allids2, allids3, num_samples, NULL, idshead, 0);
}
else	//start with response - must be null reml/he/pcgc, or no data/kins cut-folds, or family
{
strcpy(idsfile,respfile);
idshead=check_head_ids(idsfile,0);
num_samples=countrows(idsfile)-idshead;
printf("Reading IDs for %d samples from %s\n\n", num_samples, idsfile);
allids1=malloc(sizeof(char*)*num_samples);
allids2=malloc(sizeof(char*)*num_samples);
allids3=malloc(sizeof(char*)*num_samples);
read_ids(idsfile, allids1, allids2, allids3, num_samples, NULL, idshead, 0);
}
}

//check no duplicates and get idsorder
idsorder=malloc(sizeof(int)*num_samples);
check_dups(allids3,num_samples,idsfile,idsorder,1);

if(mode==122||mode==125)	//check have all individuals in blupfile
{
count=countrows(blupfile);
wantids=malloc(sizeof(char*)*count);
read_ids(blupfile, NULL, NULL, wantids, count, NULL, 0, 0);

count2=find_strings(allids3, num_samples, wantids, count, NULL, NULL, NULL, NULL, idsorder, NULL, 3);
if(count2==0){printf("Error, none of the %d samples in %s are in %s\n\n", count, blupfile, idsfile);exit(1);}
if(count2<count){printf("Error, only %d of the %d samples in %s are in %s\n\n", count2, count, blupfile, idsfile);exit(1);}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
}

////////

//usedids indicates which samples will be used
usedids=malloc(sizeof(int)*num_samples);

//start by assuming all num_samples are used
for(i=0;i<num_samples;i++){usedids[i]=1;}

if(strcmp(bsampfile,"blank")!=0)	//keep samples
{
count=countrows(bsampfile);
printf("Reading list of %d samples to keep from %s\n", count, bsampfile);
wantids=malloc(sizeof(char*)*count);
read_ids(bsampfile, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, idsorder, NULL, 3);
if(count2==0){printf("Error, none of these are in %s\n\n", idsfile);exit(1);}
if(count2<count){printf("Warning, only %d of these are in %s\n", count2, idsfile);}
printf("\n");
for(i=0;i<count2;i++){usedids[indexer[i]]++;}
for(i=0;i<num_samples;i++){usedids[i]=(usedids[i]==2);}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);
}

if(mode==162)	//calc-pca-loads - check have all individuals in pcastem, and reduce to these
{
sprintf(filename,"%s.vect",pcastem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, idsorder, NULL, 3);
if(count2==0){printf("Error, none of the %d samples in %s are in %s\n\n", count, filename, idsfile);exit(1);}
if(count2<count){printf("Error, only %d of the %d samples in %s are in %s\n\n", count2, count, filename, idsfile);exit(1);}
for(i=0;i<num_samples;i++){usedids[i]=0;}
for(i=0;i<count2;i++){usedids[indexer[i]]++;}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);
}

if(num_subs>0&&mode!=171)	//weightings, he or pcgc - reduce to samples in at least one subset
{
count=0;
for(s=0;s<num_subs;s++)
{
sprintf(filename,"%s%d",subpref,s+1);
count+=countrows(filename);
}
wantids=malloc(sizeof(char*)*count);
count=0;
for(s=0;s<num_subs;s++)
{
sprintf(filename,"%s%d",subpref,s+1);
count2=countrows(filename);
read_ids(filename, NULL, NULL, wantids+count, count2, NULL, 0, 0);
count+=count2;
}

indexer=malloc(sizeof(int)*count);
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, idsorder, NULL, 3);
if(count2==0){printf("Error, none of the samples in the %d subsets are in %s\n\n", num_subs, idsfile);exit(1);}
for(i=0;i<count2;i++){usedids[indexer[i]]++;}
for(i=0;i<num_samples;i++){usedids[i]=(usedids[i]==2);}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);
}

if(strcmp(csampfile,"blank")!=0)	//remove samples
{
count=countrows(csampfile);
printf("Reading list of %d samples to remove from %s", count, csampfile);
if(strcmp(bsampfile,"blank")!=0){printf(" (note that \"--remove\" takes priority over \"--keep\")");}
printf("\n");
wantids=malloc(sizeof(char*)*count);
read_ids(csampfile, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(allids3, num_samples, wantids, count, indexer, NULL, NULL, NULL, idsorder, NULL, 3);
if(count2<count)
{
if(count2==0){printf("Warning, none of these are in %s\n", idsfile);}
else{printf("Warning, only %d of these are in %s\n", count2, idsfile);}
}
printf("\n");
for(i=0;i<count2;i++){usedids[indexer[i]]=0;}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);
}

count=0;for(i=0;i<num_samples;i++){count+=(usedids[i]==1);}
if(count==0){printf("Error, after filtering, no samples remain\n\n");exit(1);}

if((mode==106||mode==107||mode==108||mode==109||mode==112||mode==114||mode==116||mode==117||mode==118||mode==119||mode==120||mode==121||mode==123||mode==124||mode==126||mode==127||mode==128||(mode==131&&trios==0&&duos==0)||mode==132||mode==133||mode==137||mode==138||mode==140||mode==141||mode==145||mode==151||mode==152||mode==153||mode==154||mode==156||mode==160||mode==161||mode==163||mode==164||mode==166||mode==167||mode==168||mode==169||mode==170||mode==171||mode==177||mode==181||mode==182||mode==183||mode==184||mode==185||mode==186||mode==187||mode==188||mode==189||mode==190||mode==191||mode==192||mode==193||mode==194)&&strcmp(respfile,"blank")!=0&&pad==0)
//reduce individuals to those with phenotypes (not for filter, blups, reml-pred, family or scores, nor some linear)
{check_respfile(respfile, usedids, num_samples, allids3, num_resps_use, keepresps, num_resps, missingvalue);}

////////

//see which samples remain

num_samples_use=0;
for(i=0;i<num_samples;i++){num_samples_use+=(usedids[i]==1);}

if(num_samples_use<3)
{
if(mode!=172)
{printf("Error, unable to continue with fewer than three samples; come on, you can do better ;)\n\n");exit(1);}
}

if((mode==102||mode==104)&&num_samples_use>10000&&manysamples==0)
{printf("Error, it is normally not necessary to use very many samples when calculating weightings; instead, we recommend using \"--keep rand.5000\", where the file \"rand.5000\" contains the IDs of 5000 samples, picked at random (e.g., you could make this file using the UNIX command \"shuf -n 5000 %s > rand.5000\")\nNote that you can override this error by adding \"--allow-many-samples YES\" (but the resulting analysis may be very slow)\n\n", famfile);exit(1);}

if((mode==106||mode==107)&&num_samples_use>10000&&manysamples==0)
{printf("Error, it is normally not necessary to use very many samples when thinning; instead, we recommend using \"--keep rand.5000\", where the file \"rand.5000\" contains the IDs of 5000 samples, picked at random (e.g., you could make this file using the UNIX command \"shuf -n 5000 %s > rand.5000\")\nNote that you can override this error by adding \"--allow-many-samples YES\" (but the resulting analysis may be very slow)\n\n", famfile);exit(1);}

if(mode==126&&num_samples_use<10000)
{printf("Warning, \"--fast-reml\" is designed for problems with very many samples (say, over 20000); otherwise it is normally better to use \"--reml\"\n\n");}

if(mode==138&&strcmp(sumsfile,"blank")!=0&&num_samples_use>10000&&manysamples==0)
{printf("Error, it is normally not necessary to use very many samples when performing gene/chunk-based REML with summary statistics; instead, we recommend using \"--keep rand.5000\", where the file \"rand.5000\" contains the IDs of 5000 samples, picked at random (e.g., you could make this file using the UNIX command \"shuf -n 5000 %s > rand.5000\")\nNote that you can override this error by adding \"--allow-many-samples YES\" (but the resulting analysis may be very slow)\n\n", famfile);exit(1);}

if(mode==140&&num_samples_use>10000&&manysamples==0)	//reduce to 5000 samples
{
printf("Warning, this analysis does not require very many samples, so LDAK will reduce to 5000, picked at random (you can force LDAK to use all samples by adding \"--allow-many-samples YES\", but the resulting analysis may be very slow)\n\n");

order=malloc(sizeof(int)*num_samples_use);
count=0;
for(i=0;i<num_samples;i++)
{
if(usedids[i]==1){order[count]=i;count++;}
}
for(i=0;i<num_samples;i++){usedids[i]=0;}

num_samples_use=5000;
permute_int(order,count);
qsort(order,num_samples_use,sizeof(int), compare_int);
for(i=0;i<num_samples_use;i++){usedids[order[i]]=1;}

free(order);
}

if(mode==141&&num_samples_use>10000&&manysamples==0)
{printf("Error, it is normally not necessary to use very many samples when calculating taggings; instead, we recommend using \"--keep rand.5000\", where the file \"rand.5000\" contains the IDs of 5000 samples, picked at random (e.g., you could make this file using the UNIX command \"shuf -n 5000 %s > rand.5000\")\nNote that you can override this error by adding \"--allow-many-samples YES\" (but the resulting analysis may be very slow)\n\n", famfile);exit(1);}

if(mode==156&&num_samples_use>10000&&manysamples==0)
{printf("Error, it is normally not necessary to use very many samples when calculating correlations; instead, we recommend using \"--keep rand.5000\", where the file \"rand.5000\" contains the IDs of 5000 samples, picked at random (e.g., you could make this file using the UNIX command \"shuf -n 5000 %s > rand.5000\")\nNote that you can override this error by adding \"--allow-many-samples YES\" (but the resulting analysis may be very slow)\n\n", famfile);exit(1);}

if(mode==160&&num_samples_use<500)
{printf("Warning, it is difficult to reliably measure the accuracy of prediction models with so few samples (ideally, there should be at least 500 samples)\n\n");}

if((mode==126||mode==127||mode==128)&&num_vects==-9999)	//can set num_vects for fastreml, fasthe and fastpcgc
{
if(num_samples_use>20000){num_vects=20;}
else{num_vects=100;}

printf("Will use %d random vectors (to change this, use \"--repetitions\")\n\n", num_vects);
}

if(mode==138||mode==140)	//deal with nclump and limit
{
if(saveall==1){nclump=1;}
else
{
if(num_samples_use<5000){nclump=1;}
else{nclump=num_samples_use/5000;}
}

if(nclump>1){printf("Will only save estimated genetic effects for %d samples, which is sufficient for clumping (you can force LDAK to save effects for all samples by adding \"--save-all YES\")\n\n", (num_samples_use-1)/nclump+1);}

if(limit==-9999)
{
if(strcmp(sumsfile,"blank")!=0&&num_samples_use<2000){limit=2;}
else{limit=0;}
}

if(limit>0){printf("Will exclude genes whose estimated heritability is more than %.4f times an estimate from least-squares regression (to change this threshold use \"--her-limit\")\n\n", limit);}
}

if((mode==151||mode==152||mode==153||mode==154)&&nmcmc==-9999)	//can set nmcmc
{
if(num_samples_use>40000){nmcmc=3;}
else{nmcmc=10;}
}

ids1=malloc(sizeof(char*)*num_samples_use);
ids2=malloc(sizeof(char*)*num_samples_use);
ids3=malloc(sizeof(char*)*num_samples_use);

count=0;
for(i=0;i<num_samples;i++)
{
if(usedids[i]==1)
{
copy_string(ids1,count,allids1[i]);
copy_string(ids2,count,allids2[i]);
copy_string(ids3,count,allids3[i]);
count++;
}
}
for(i=0;i<num_samples;i++){free(allids1[i]);free(allids2[i]);free(allids3[i]);}
free(allids1);free(allids2);free(allids3);free(idsorder);free(usedids);

////////

//some mode specific checks

if(mode==117&&extract==1)	//check have data for all samples
{
count=countrows(famfile)-famhead;
wantids=malloc(sizeof(char*)*count);
read_ids(famfile, NULL, NULL, wantids, count, NULL, 0, 0);

count2=find_strings(ids3, num_samples_use, wantids, count, NULL, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Error, can not find data for any of these\n\n");exit(1);}
if(count2<num_samples_use)
{printf("Error, data are only available for %d of these\n\n", count2);exit(1);}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
}

if(mode==122||mode==125||mode==162)	//if filtering, check still have individuals in blupfile / pcastems
{
if(strcmp(bsampfile,"blank")!=0||strcmp(csampfile,"blank")!=0)
{
if(mode==122||mode==125){strcpy(filename,blupfile);}
else{sprintf(filename,"%s.vect",pcastem);}
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

count2=find_strings(wantids, count, ids3, num_samples_use, NULL, NULL, NULL, NULL, NULL, NULL, 3);
if(count2<count)
{printf("Error, after filtering samples, only %d of the %d samples in %s remain\n\n", count2, count, filename);exit(1);}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
}
}

if(mode==125)	//check some extra individuals (those we predict for)
{
count=countrows(blupfile);
if(count==num_samples_use)
{printf("Error, all %d samples are in %s; \"--reml-pred\" can only predict phenotypes for samples who were not used when running \"--reml\" (or \"--fast-reml\")\n\n", num_samples_use, blupfile);exit(1);}
}

if((mode==131||mode==138)&&num_kins>0)	//check have kinships for all
{
count=find_strings(ids3, num_samples_use, kinids, total, NULL, NULL, NULL, NULL, NULL, NULL, 1);
if(count==0){printf("Error, can not find kinships for any of these\n\n");exit(1);}
if(count<num_samples_use)
{printf("Error, kinships are only available for %d of the samples (if this is intentional, use \"--keep\" or \"--remove\" to exclude the samples without kinships from the analysis)\n\n", count);exit(1);}
}

if(mode==161||mode==167)	//check number of axes
{
if(axes>num_samples_use){printf("Error, the number of axes (%d) can not be larger than the number of samples (%d)\n\n", axes, num_samples_use);exit(1);}
}

////////

if(num_subs>0)	//fill subsets - modes 102/104/123/124/171 - for modes 123/124, some subsets can be empty
{
subindex=malloc(sizeof(int*)*num_subs);
count2=num_subs;
for(s=0;s<num_subs;s++)
{
sprintf(filename, "%s%d", subpref, s+1);
count=countrows(filename);
printf("Reading subset of %d samples from %s\n", count, filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

subindex[s]=malloc(sizeof(int)*(1+count));
subindex[s][0]=find_strings(ids3, num_samples_use, wantids, count, subindex[s]+1, NULL, NULL, filename, NULL, NULL, 3);
if(subindex[s][0]==0)
{
if(mode==123||mode==124)
{printf("Warning, can not find any of these samples in %s, so this subset will be excluded\n", idsfile);count2--;}
else
{printf("Error, can not find any of these samples in %s\n\n", idsfile);exit(1);}
}
if(subindex[s][0]>0&&subindex[s][0]<count)
{printf("Warning, will be using only %d of these\n", subindex[s][0]);}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
}

if((mode==123||mode==124)&&count2==1){printf("Error, it is not possible to continue as only one subset remains\n\n");exit(1);}
if((mode==123||mode==124)&&count2<num_subs){printf("%d subsets remain\n", count2);}
printf("\n");

if(mode==123||mode==124)	//must check no overlap
{
usedids=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){usedids[i]=0;}
for(s=0;s<num_subs;s++)
{
for(i=0;i<subindex[s][0];i++){usedids[subindex[s][1+i]]++;}
}
count=0;for(i=0;i<num_samples_use;i++){count+=(usedids[i]>1);}
if(count>0)
{printf("Error, %d samples feature in more than one subset; when using \"--he\" or \"--pcgc\", subsets must not overlap\n\n", count);exit(1);}
free(usedids);
}
}

if(num_subs==0&&(mode==102||mode==104||mode==123||mode==124))	//convenient to make one subset
{
num_subs=1;
subindex=malloc(sizeof(int*));
subindex[0]=malloc(sizeof(int)*(1+num_samples_use));
subindex[0][0]=num_samples_use;
for(i=0;i<num_samples_use;i++){subindex[0][1+i]=i;}
}

////////

if(mode==191)	//save samples
{
sprintf(filename,"%scorrelations.samples",folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(i=0;i<num_samples_use;i++){fprintf(output,"%s %s\n", ids1[i], ids2[i]);}
fclose(output);
}

if(mode==192||mode==193||mode==194)	//check match those saved and sensible
{
sprintf(filename,"%scorrelations.samples",folder);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when using \"--cut-gre\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=num_samples_use)
{printf("\nError, the number of samples remaining (%d) does not match the number in %s (%d), indicating different sample filterings were used with \"--cut-gre\"\n\n", num_samples_use, filename, count);exit(1);}

printf("Reading list of %d samples used with \"--cut-gre\" from %s\n\n", count, filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

for(i=0;i<num_samples_use;i++)
{
if(strcmp(ids3[i],wantids[i])!=0)
{printf("Error, the samples remaining do not match those in %s, indicating different sample filterings were used with \"--cut-gre\"\n\n", filename);exit(1);}
}

for(i=0;i<num_samples_use;i++){free(wantids[i]);}free(wantids);
}

////////

if(use_data==1||use_data==2||use_data==4)	//set keepsamps, the indexes corresponding to ids3
{
num_samples=countrows(famfile)-famhead;
wantids=malloc(sizeof(char*)*num_samples);
read_ids(famfile, NULL, NULL, wantids, num_samples, NULL, famhead, 0);

keepsamps=malloc(sizeof(int)*num_samples_use);
count=find_strings(wantids, num_samples, ids3, num_samples_use, keepsamps, NULL, NULL, NULL, NULL, NULL, 3);
if(count<num_samples_use){printf("Error 3DE; please tell Doug %d %d\n\n", count, num_samples_use);exit(1);}
for(i=0;i<num_samples;i++){free(wantids[i]);}free(wantids);
}

if(use_data==3)	//using data only through toppreds - check some samples in common
{
num_samples=countrows(famfile)-famhead;
wantids=malloc(sizeof(char*)*num_samples);
read_ids(famfile, NULL, NULL, wantids, num_samples, NULL, famhead, 0);

count=find_strings(wantids, num_samples, ids3, num_samples_use, NULL, NULL, NULL, NULL, NULL, NULL, 3);
if(count==0){printf("Error, none of the samples are in %s\n\n", famfile);exit(1);}
if(count<num_samples_use){printf("Warning, only %d of the samples are in %s\n\n", count, famfile);}
for(i=0;i<num_samples;i++){free(wantids[i]);}free(wantids);
}

if(num_kins>0)	//free allocations from above
{
for(i=0;i<total;i++){free(kinids[i]);}free(kinids);
}
}	//end of using sample ids

///////////////////////////

if(use_data==1||use_data==2||use_data==3||use_data==6)	//will be using predictor ids, so set num_preds, num_preds_use and get predictor details
{
//read in all predictors
if(strcmp(bimfile,datafile)!=0)	//have a separate bimfile
{
num_preds=countrows(bimfile);
printf("Reading details for %d predictors from %s\n", num_preds, bimfile);
allchr=malloc(sizeof(int)*num_preds);
allcm=malloc(sizeof(double)*num_preds);
allbp=malloc(sizeof(double)*num_preds);
allpreds=malloc(sizeof(char*)*num_preds);
allal1=malloc(sizeof(char)*num_preds);
allal2=malloc(sizeof(char)*num_preds);
read_bimfile(bimfile, allchr, allpreds, allcm, allbp, allal1, allal2, num_preds, NULL, 1, window_cm, 0, multi);
}
else	//must be dtype=5, no bimfile, read from genfile (so can not have genetic distances)
{
//assume at most maxpreds SNPs
printf("Reading predictor details from %s\n", bimfile);
allchr=malloc(sizeof(int)*maxpreds);
allcm=malloc(sizeof(double)*maxpreds);
allbp=malloc(sizeof(double)*maxpreds);
allpreds=malloc(sizeof(char*)*maxpreds);
allal1=malloc(sizeof(char)*maxpreds);
allal2=malloc(sizeof(char)*maxpreds);
for(j=0;j<maxpreds;j++){allchr[j]=oxchr;allcm[j]=0;}

num_preds=read_genfile(bimfile, allpreds, allbp, allal1, allal2, num_samples, genskip, genheaders, genprobs, maxpreds, 1);
if(num_preds>maxpreds)	//then file larger than allocated, will have to read again
{
for(j=0;j<maxpreds;j++){free(allpreds[j]);}free(allpreds);
free(allchr);free(allbp);free(allal1);free(allal2);
allchr=malloc(sizeof(int)*num_preds);
allbp=malloc(sizeof(double)*num_preds);
allpreds=malloc(sizeof(char*)*num_preds);
allal1=malloc(sizeof(char)*num_preds);
allal2=malloc(sizeof(char)*num_preds);
for(j=0;j<num_preds;j++){allchr[j]=oxchr;allcm[j]=0;}
read_genfile(bimfile, allpreds, allbp, allal1, allal2, num_samples, genskip, genheaders, genprobs, num_preds, 1);
}
if(num_preds>=20000){printf("%s contains genotypes for %d SNPs\n", bimfile, num_preds);}
}
printf("\n");

//check for duplicates and get predorder
predorder=malloc(sizeof(int)*num_preds);
check_dups(allpreds,num_preds,bimfile,predorder,1);

////////

if(strcmp(topfile,"blank")!=0)	//check tops in the data (and no duplicates), and set tkeeppreds
{
wantpreds=malloc(sizeof(char*)*num_tops);
read_strings(topfile, wantpreds, num_tops, NULL, 1, 0);
tkeeppreds=malloc(sizeof(int)*num_tops);
count=find_strings(allpreds, num_preds, wantpreds, num_tops, tkeeppreds, NULL, NULL, topfile, predorder, NULL, 3);
if(count==0){printf("Error, none of the %d top predictors are in the data\n\n", num_tops);exit(1);}
if(count<num_tops)
{printf("Error, only %d of the %d top predictors are in the data\n\n", count, num_tops);exit(1);}
for(j=0;j<num_tops;j++){free(wantpreds[j]);}free(wantpreds);
}

if(mode==173&&num_causals>num_preds)
{printf("Error, there are only %d predictors in the data, so it will not be possible to generate phenotypes with %d causal predictors\n\n", num_preds, num_causals);exit(1);}

////////

//now set usedpreds

if(mode==117&&extract==1)	//subbing predictors from a kinship matrix
{
//first get list of predictors used in kinship matrix
sprintf(filename, "%s.grm.details", kinstems[0]);
count=countrows(filename)-1;
kinpreds=malloc(sizeof(char*)*count);
read_strings(filename, kinpreds, count, NULL, 1, 1);

//usedpreds first indicates which kinship predictors we will subtract
usedpreds=malloc(sizeof(int)*count);

//flag indicates if no kinship predictors to subtract
flag=0;

if(strcmp(bpredfile,"blank")!=0)	//extracting predictors - will subtract those not in bpredfile
{
count2=countrows(bpredfile);
printf("Reading list of %d predictors to retain from %s\n", count2, bpredfile);
wantpreds=malloc(sizeof(char*)*count2);
read_strings(bpredfile, wantpreds, count2, NULL, 1, 0);

indexer=malloc(sizeof(int)*count2);
count3=find_strings(kinpreds, count, wantpreds, count2, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count3==0){printf("Error, none of these contribute towards the kinship matrix\n\n");exit(1);}
if(count3==count)
{printf("Warning, this list contains all %d kinship predictors, so no predictors will be subtracted\n\n", count);flag=1;}
else
{
if(count3<count2){printf("Warning, only %d of these contribute towards the kinship matrix\n\n", count3);}
else{printf("All of these contribute towards the kinship matrix\n\n");}
}

for(j=0;j<count;j++){usedpreds[j]=1;}
for(j=0;j<count3;j++){usedpreds[indexer[j]]=0;}
for(j=0;j<count2;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}
else	//excluding predictors - will subtract those in cpredfile
{
count2=countrows(cpredfile);
printf("Reading list of %d predictors to subtract from %s\n", count2, cpredfile);
wantpreds=malloc(sizeof(char*)*count2);
read_strings(cpredfile, wantpreds, count2, NULL, 1, 0);

indexer=malloc(sizeof(int)*count2);
count3=find_strings(kinpreds, count, wantpreds, count2, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count3==count)
{printf("Error, this list contains all %d kinship predictors; the resulting kinship matrix would be null\n\n", count);exit(1);}
if(count3==0){printf("Warning, none of these contribute towards the kinship matrix, so no predictors will be subtracted\n\n");flag=1;}
else
{
if(count3<count2){printf("Warning, only %d of these contribute towards the kinship matrix\n", count3);}
else{printf("All of these contribute towards the kinship matrix\n");}
}

for(j=0;j<count;j++){usedpreds[j]=0;}
for(j=0;j<count3;j++){usedpreds[indexer[j]]=1;}
for(j=0;j<count2;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}

if(flag==1)	//no predictors to subtract, so simply copy over
{
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);
read_kins(kinstems[0], kins, NULL, 1, num_samples_use, ids3, 0, maxthreads);
write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 1);
free(kins);
sprintf(cmd, "cp %s.grm.details %s.grm.details", kinstems[0], outfile);
system(cmd);
sprintf(cmd, "cp %s.grm.adjust %s.grm.adjust", kinstems[0], outfile);
system(cmd);
exit(1);
}

//finally, reset usedpreds so that it aligns with allpreds
count2=0;for(j=0;j<count;j++){count2+=usedpreds[j];}
wantpreds=malloc(sizeof(char*)*count2);
count2=0;
for(j=0;j<count;j++)
{
if(usedpreds[j]==1){copy_string(wantpreds,count2,kinpreds[j]);count2++;}
}
for(j=0;j<count;j++){free(kinpreds[j]);}free(kinpreds);

indexer=malloc(sizeof(int)*count2);
count3=find_strings(allpreds, num_preds, wantpreds, count2, indexer, NULL, NULL, NULL, predorder, NULL, 3);
if(count3==0)
{printf("Error, none of the %d predictors to subtract are in the data, indicating different data files were used when calculating kinships\n\n", count2);exit(1);}
if(count3<count2)
{printf("Error, only %d of the %d predictors to subtract are in the data, indicating different data files were used when calculating kinships\n\n", count3, count2);exit(1);}
printf("Will be subtracting %d of the %d kinship predictors\n\n", count2, count);

free(usedpreds);
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(j=0;j<count2;j++){usedpreds[indexer[j]]=1;}
for(j=0;j<count2;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}	//end of subtracting predictors

////////

if(mode!=117||extract!=1)	//normal case
{
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=1;}

//blank out predictors with non-positive basepairs (prob due to multi-character alleles)
count=0;
for(j=0;j<num_preds;j++)
{
if(allbp[j]<=0){usedpreds[j]=0;count++;}
if(count==num_preds){printf("Error, no predictors have non-positive basepairs or single-character alleles\n\n");exit(1);}
}

//do predictor filtering
(void)extraction(usedpreds, num_preds, allpreds, allchr, predorder, bpredfile, cpredfile, onechr, onesnp, bimfile);

if(mode==107)	//reduce to predictors with P<cutoff
{
head=check_head(pvafile,"Predictor","SNP",0);
count=countrows(pvafile)-head;
wantpreds=malloc(sizeof(char*)*count);
pvalues=malloc(sizeof(double)*count);
read_strings(pvafile, wantpreds, count, NULL, 1, head);
read_values(pvafile, pvalues, count, NULL, 2, head, 1);

count2=0;
for(j=0;j<count;j++)
{
if(pvalues[j]!=-9999&&pvalues[j]<=cutoff)
{
if(count2!=j){free(wantpreds[count2]);copy_string(wantpreds,count2,wantpreds[j]);}
count2++;
}
}
if(count2==0){printf("Error, none of the %d predictors in %s have P<= %.2e\n\n", count, pvafile, cutoff);exit(1);}
printf("Of the %d predictors in %s, %d have P <= %.2e\n", count, pvafile, count2, cutoff);

indexer=malloc(sizeof(int)*count2);
count3=find_strings(allpreds, num_preds, wantpreds, count2, indexer, NULL, NULL, pvafile, predorder, NULL, 3);
if(count3==0){printf("Error, none of these are in the data\n\n");exit(1);}

count4=0;for(j=0;j<count3;j++){count4+=usedpreds[indexer[j]];}
if(count4==0){printf("Error, after filtering predictors, none remain\n\n");exit(1);}
if(count4<count3){printf("Warning, only %d of these are in the data\n\n", count4);}
else{printf("All of these are in the data\n\n");}

for(j=0;j<count3;j++){usedpreds[indexer[j]]++;}
for(j=0;j<num_preds;j++){usedpreds[j]=(usedpreds[j]==2);}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(pvalues);free(indexer);
}

if(strcmp(invsfile,"blank")!=0)	//reduce to predictors in inverse file (must be mode 114)
{
sprintf(filename,"%s.inverse.predictors",invsfile);
count=countrows(filename);
printf("Reading list of %d predictors from %s\n", count, filename);
wantpreds=malloc(sizeof(char*)*count);
read_strings(filename, wantpreds, count, NULL, 1, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(allpreds, num_preds, wantpreds, count, indexer, NULL, NULL, filename, predorder, NULL, 3);
if(count2==0){printf("Error, none of these are in the data\n\n");exit(1);}
if(count2<count){printf("Error, only %d of these are in the data\n\n", count2);exit(1);}

count2=0;for(j=0;j<count;j++){count2+=usedpreds[indexer[j]];}
if(count2==0){printf("Error, after filtering predictors, none remain\n\n");exit(1);}
if(count2<count){printf("Error, after filtering predictors, only %d of these remain\n\n", count2);exit(1);}

for(j=1;j<count;j++)
{
if(indexer[j]<indexer[j-1]){printf("Error, the order of predictors in %s does not match the order in %s (e.g., Predictors %s and %s)\n\n", filename, bimfile, allpreds[indexer[j-1]], allpreds[indexer[j]]);exit(1);}
}

for(j=0;j<count;j++){usedpreds[indexer[j]]++;}
for(j=0;j<num_preds;j++){usedpreds[j]=(usedpreds[j]==2);}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(pvalues);free(indexer);
}

if((mode==151||mode==152||mode==153||mode==154)&&strcmp(topfile,"blank")!=0)	//remove top predictors
{
for(j=0;j<num_tops;j++){usedpreds[tkeeppreds[j]]=0;}
}

if(mode==173&&strcmp(weightsfile,"blank")!=0)	//reduce to predictors with non-zero weights
{
head=check_head(weightsfile,"Predictor","SNP",0);
count=countrows(weightsfile)-head;
printf("Identifying predictors with non-zero weight from %s\n", weightsfile);
wantpreds=malloc(sizeof(char*)*count);
weights=malloc(sizeof(double)*count);
read_strings(weightsfile, wantpreds, count, NULL, 1, head);
read_values(weightsfile, weights, count, NULL, 2, head, 1);

count2=0;
for(j=0;j<count;j++)
{
if(weights[j]>0)
{
if(count2!=j){free(wantpreds[count2]);copy_string(wantpreds,count2,wantpreds[j]);}
count2++;
}
}
if(count2==0){printf("Error, none of the %d predictors have non-zero weight\n\n", count);exit(1);}
printf("Of the %d predictors, %d have non-zero weight\n", count, count2);

indexer=malloc(sizeof(int)*count2);
count3=find_strings(allpreds, num_preds, wantpreds, count2, indexer, NULL, NULL, weightsfile, predorder, NULL, 3);
if(count3==0){printf("Error, none of these are in the data\n\n");exit(1);}

count4=0;for(j=0;j<count3;j++){count4+=usedpreds[indexer[j]];}
if(count4==0){printf("Error, after filtering predictors, none remain\n\n");exit(1);}
if(count4<count3){printf("Warning, only %d of these are in the data\n\n", count4);}
else{printf("All of these are in the data\n\n");}

for(j=0;j<count3;j++){usedpreds[indexer[j]]++;}
for(j=0;j<num_preds;j++){usedpreds[j]=(usedpreds[j]==2);}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(weights);free(indexer);
}	//end of mode=173 and using weights
}	//end of mode!=117|extract!=1

////////

//set keeppreds
keeppreds=malloc(sizeof(int)*num_preds);
num_preds_use=0;
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]==1){keeppreds[num_preds_use]=j;num_preds_use++;}
}

if(num_preds_use==0){printf("Error, after filtering predictors, none remain\n\n");exit(1);}

//work out num_chr
num_chr=1;
for(j=1;j<num_preds_use;j++)	
{
if(allchr[keeppreds[j]]!=allchr[keeppreds[j-1]]){num_chr++;}
}

//do some checks and set some values

if(window_length==-1){window_length=num_preds_use;}

if(loco==1&&num_chr==1)
{printf("Error, there is only one chromosome, so it is not possible to use \"--LOCO YES\"\n\n");exit(1);}

if((mode==151||mode==152||mode==153||mode==154)&&num_preds_use>=1000000&&manypreds==0)
{printf("Error, it is normally not necessary to use more than about 500,000 predictors when computing the LOCO PRS; instead, we recommend using \"--extract highQC.predictors\", where the file \"highQC.predictors\" contains the IDs of common, high-quality predictors\nNote that you can override this error by adding \"--allow-many-predictors YES\" (but the resulting analysis may be very slow)\n\n");exit(1);}

if(mode==173)
{
if(num_causals>num_preds_use)
{
printf("Error, it will not be possible to generate phenotypes with %d causal predictors, because only %d predictors remain after filtering", num_causals, num_preds_use);
if(strcmp(weightsfile,"blank")!=0){printf(" and removing those with weight zero");}
printf("\n\n");
exit(1);
}

if(num_causals==-1){num_causals=num_preds_use;}
}

if(mode==191&&num_samples_use<=num_preds_use)
{printf("Error, the number of samples (%d) should be higher than the number of predictors (%d)\n\n", num_samples_use, num_preds_use);exit(1);}

////////

/*
if(strcmp(locofile,"blank")!=0)	//check sample size from root file
{
//open root file and skip 7 rows
sprintf(filename,"%s.root", locofile);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
for(j=0;j<7;j++)
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

//read and check num_samples_used
if(fscanf(input, "%s %d ", readstring, &readint)!=2)
{printf("Error reading Row 8 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Num_Samples_Used")!=0)
{printf("Error, reading Row 8 of %s (should begin \"Num_Samples_Used\"), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);exit(1);}
if(num_samples_use!=readint)
{printf("Error, %d samples were used when making the LOCO PRS, but %d are used now, indicating different sample filterings are being used\n\n", readint, num_samples_use);exit(1);}

fclose(input);
}
*/

///////////////////////////

if(mode==101||mode==111||mode==136||mode==191)	//save predictors
{
if(mode==101){sprintf(filename,"%sweights.predictors",folder);}
if(mode==111){sprintf(filename,"%skinships.predictors",folder);}
if(mode==136){sprintf(filename,"%sgenes.predictors",folder);}
if(mode==191){sprintf(filename,"%scorrelations.predictors",folder);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(j=0;j<num_preds_use;j++){fprintf(output,"%s\n", allpreds[keeppreds[j]]);}
fclose(output);
}

if(mode==102||mode==103||mode==104||mode==112||mode==113||mode==137||mode==138||mode==139||mode==192||mode==193||mode==194)	//check match those saved and sensible
{
if(mode==102||mode==103||mode==104)
{sprintf(filename,"%sweights.predictors",folder);strcpy(writestring,"\"--cut-weights\"");}
if(mode==112||mode==113)
{sprintf(filename,"%skinships.predictors",folder);strcpy(writestring,"\"--cut-kins\"");}
if(mode==137||mode==138||mode==139)
{sprintf(filename,"%sgenes.predictors",folder);strcpy(writestring,"\"--cut-genes\"");}
if(access(filename, R_OK)!=0)
{printf("Error, unable to read %s, this should have been created when using %s\n\n", filename, writestring);exit(1);}
if(mode==192||mode==193||mode==194)
{sprintf(filename,"%scorrelations.predictors",folder);strcpy(writestring,"\"--cut-gre\"");}

if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when using %s\n\n", filename, writestring);exit(1);}

count=countrows(filename);
if(count!=num_preds_use)
{printf("Error, the number of predictors remaining (%d) does not match the number in %s (%d), indicating different predictor filterings were used with %s\n\n", num_preds_use, filename, count, writestring);exit(1);}

printf("Reading list of %d predictors used with %s from %s\n", count, writestring, filename);
wantpreds=malloc(sizeof(char*)*count);
read_strings(filename, wantpreds, count, NULL, 1, 0);

for(j=0;j<num_preds_use;j++)
{
if(strcmp(allpreds[keeppreds[j]],wantpreds[j])!=0)
{printf("Error, the predictors remaining do not match those in %s, indicating different predictor filterings were used with %s\n\n", filename, writestring);exit(1);}
}

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);
printf("\n");
}

free(usedpreds);
}	//end of use_data=1/2/3

///////////////////////////

