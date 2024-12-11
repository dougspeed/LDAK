/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

//////////////////////////

//Find families

///////////////////////////

printf("Searching for families (using the family IDs in the first column of %s)\n\n", famfile);

//sort based on ids1
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

if(num_fams==num_samples_use){printf("Error, each of the %d samples belongs to a different family (check the first column of %s contains family IDs)\n\n" ,num_samples_use, famfile);exit(1);}
if(num_fams==1){printf("Error, all %d samples belongs to the same family (check the first column of %s contains family IDs)\n\n" ,num_samples_use, famfile);exit(1);}

printf("In total, the %d samples belong to %d families\n\n", num_samples_use, num_fams);

free(sptrs);

///////////////////////////

