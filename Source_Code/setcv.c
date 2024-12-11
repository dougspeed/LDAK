/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Deal with CV samples

///////////////////////////

if(cvprop!=-9999)	//pick cvprop at random to be test
{
num_test=cvprop*num_samples_use;
num_train=num_samples_use-num_test;

//pick samples at random
for(i=0;i<num_samples_use;i++){keeptrain[i]=i;}
permute_int(keeptrain,num_samples_use);
for(i=0;i<num_test;i++){keeptest[i]=keeptrain[num_train+i];}
qsort(keeptrain,num_train,sizeof(int), compare_int);
qsort(keeptest,num_test,sizeof(int), compare_int);
}
else	//using ids from bvsfile
{
count=countrows(bvsfile);
printf("Reading list of %d cross-validation samples from %s\n", count, bvsfile);
wantids=malloc(sizeof(char*)*count);
read_ids(bvsfile, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
num_test=find_strings(ids3, num_samples_use, wantids, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(num_test==0){printf("Error, can not find any of these samples in the data\n\n");exit(1);}
if(num_test<count){printf("Warning, only %d of these are in the data\n", num_test);}
num_train=num_samples_use-num_test;

usedids=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){usedids[i]=0;}
for(i=0;i<num_test;i++){usedids[indexer[i]]=1;}
count2=0;count3=0;
for(i=0;i<num_samples_use;i++)
{
if(usedids[i]==0){keeptrain[count2]=i;count2++;}
else{keeptest[count3]=i;count3++;}
}
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);free(usedids);
}

//check sufficient non-missing phenotypes
for(m=0;m<num_resps_use;m++)
{
count=0;
for(i=0;i<num_train;i++){count+=(resp[keeptrain[i]+m*num_samples_use]!=missingvalue);}
count2=0;
for(i=0;i<num_test;i++){count2+=(resp[keeptest[i]+m*num_samples_use]!=missingvalue);}

if(mpheno!=-1)
{
printf("When performing cross-validation, will use %d samples to train models and %d to test their accuracy\n", count, count2);
if(count<3){printf("Error, unable to continue with fewer than three training samples\n\n");exit(1);}
if(count2<3){printf("Error, unable to continue with fewer than three test samples\n\n");exit(1);}
}
else
{
printf("For phenotype %d, will use %d samples to train models and %d to test their accuracy\n", m+1, count, count2);
if(count<3){printf("Error, unable to continue with fewer than three training samples\n\n");exit(1);}
if(count2<3){printf("Error, unable to continue with fewer than three test samples\n\n");exit(1);}
}
}
printf("\n");

///////////////////////////

