/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

//////////////////////////

//Find trios - will update ids, num_samples_use and keepsamps

///////////////////////////

printf("Searching for offspring-father-mother trios (using parental IDs from the 3rd and 4th columns of %s\n\n", famfile);

pid=malloc(sizeof(char*)*num_samples_use);
mid=malloc(sizeof(char*)*num_samples_use);
retain=malloc(sizeof(int)*num_samples_use);
usedids=malloc(sizeof(int)*num_samples_use);
indexer=malloc(sizeof(int)*3*num_samples_use);
indexer2=malloc(sizeof(int)*3*num_samples_use);

//usedids first indicates which have both parents listed
read_parents(famfile, pid, mid, num_samples_use, ids3, famhead);
for(i=0;i<num_samples_use;i++){usedids[i]=(strcmp(pid[i],"0")!=0&&strcmp(mid[i],"0")!=0);}

count=0;for(i=0;i<num_samples_use;i++){count+=usedids[i];}

if(count==0){printf("Error, none of the %d samples have both parental IDs\n\n", num_samples_use);exit(1);}
printf("There are %d samples with both parental IDs\n\n", count);

//usedids now updated to reflect which also have phenotypes present
check_respfile(respfile, usedids, num_samples_use, ids3, num_resps_use, keepresps, num_resps, missingvalue);

//squeeze down pid and mid, saving index
count=0;
for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==1)
{
if(count!=i)	//squeeze down
{
free(pid[count]);copy_string(pid,count,pid[i]);
free(mid[count]);copy_string(mid,count,mid[i]);
}
retain[count]=i;
count++;
}
}

if(count==0){printf("Error, none of the %d individuals have phenotypes and parental IDs\n\n", num_samples_use);exit(1);}

//count how many parents in data
for(i=0;i<count;i++){usedids[i]=0;}

count2=find_strings(pid, count, ids3, num_samples_use, indexer, NULL, NULL, NULL, NULL, NULL, 3);
for(i=0;i<count2;i++){usedids[indexer[i]]++;}

count2=find_strings(mid, count, ids3, num_samples_use, indexer, NULL, NULL, NULL, NULL, NULL, 3);
for(i=0;i<count2;i++){usedids[indexer[i]]++;}

//reduce pid and mid again, saving index
count2=0;
for(i=0;i<count;i++)
{
if(usedids[i]==2)
{
if(count2!=i)	//squeeze down
{
free(pid[count2]);copy_string(pid,count2,pid[i]);
free(mid[count2]);copy_string(mid,count2,mid[i]);
}
indexer[count2]=retain[i];
count2++;
}
}

if(count2==0){printf("Error, none of the %d individuals have phenotypes and their parents in the data\n\n", num_samples_use);exit(1);}

//add indexes of fathers
count3=find_strings(pid, count2, ids3, num_samples_use, NULL, indexer+count2, NULL, NULL, NULL, NULL, 3);
if(count3!=count2){printf("Error 7DEF; please tell Doug %d %d\n\n", count2, count3);exit(1);}

//and indexes of mothers
count3=find_strings(mid, count2, ids3, num_samples_use, NULL, indexer+2*count2, NULL, NULL, NULL, NULL, 3);
if(count3!=count2){printf("Error 7DEG; please tell Doug %d %d\n\n", count2, count3);exit(1);}

//update IDs (this works because first count2 elements of indexer are in order)
for(i=0;i<count2;i++)
{
if(indexer[i]!=i)
{
free(ids1[i]);copy_string(ids1,i,ids1[indexer[i]]);
free(ids2[i]);copy_string(ids2,i,ids2[indexer[i]]);
free(ids3[i]);copy_string(ids3,i,ids3[indexer[i]]);
}
}

for(i=count2;i<num_samples_use;i++){free(ids1[i]);free(ids2[i]);free(ids3[i]);}

//should free pid and mid before we update num_samples_use
for(i=0;i<num_samples_use;i++){free(pid[i]);free(mid[i]);}free(pid);free(mid);

//update num_samples_use and keepsamps
num_samples_use=count2;
for(i=0;i<3*count2;i++){indexer2[i]=keepsamps[indexer[i]];}
free(keepsamps);
keepsamps=malloc(sizeof(int)*3*num_samples_use);
for(i=0;i<3*count2;i++){keepsamps[i]=indexer2[i];}

//get numbers of families (code copied from findfamilies.c)
sptrs=malloc(sizeof(struct sorting_string)*num_samples_use);
for(i=0;i<num_samples_use;i++){sptrs[i].ptr=ids1[i];sptrs[i].index=i;}
qsort(sptrs, num_samples_use, sizeof(struct sorting_string), compare_sorting_string);

//load famindex (indicates which family each individual belongs to))
famindex=malloc(sizeof(int)*num_samples_use);
famcounts=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){famcounts[i]=0;}

num_fams=1;
famindex[sptrs[0].index]=0;
famcounts[0]=1;
for(i=1;i<num_samples_use;i++)
{
if(strcmp(sptrs[i].ptr,sptrs[i-1].ptr)!=0){num_fams++;}
famindex[sptrs[i].index]=num_fams-1;
famcounts[num_fams-1]++;
}

printf("There are %d offspring with phenotypes and both parents in the data\n\n", num_samples_use);

if(num_fams==1){printf("Error, all offspring belongs to the same family (check the first column of %s contains family IDs)\n\n", famfile);exit(1);}

printf("The offspring span %d families\n\n", num_fams);

free(retain);free(usedids);free(indexer);free(indexer2);
free(sptrs);

///////////////////////////

