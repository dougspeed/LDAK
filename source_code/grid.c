/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Grid search to find gaussian distribution(s)

///////////////////////////

count=countcols(likefile);

if(count==2)	//performing 1-d search
{
mvals=malloc(sizeof(double)*num_means);
for(j=0;j<num_means;j++){mvals[j]=minmean+j*(maxmean-minmean)/(num_means-1);}

svals=malloc(sizeof(double)*num_sds);
for(k=0;k<num_sds;k++){svals[k]=(k+1)*maxsd/num_sds;}

num_pows=countrows(likefile);
powers=malloc(sizeof(double)*num_pows);
likes=malloc(sizeof(double)*num_pows);
read_values(likefile,powers,num_pows,NULL,1,0,0);
read_values(likefile,likes,num_pows,NULL,2,0,0);

perfs=malloc(sizeof(double*)*num_means);
for(j=0;j<num_means;j++){perfs[j]=malloc(sizeof(double)*num_sds);}

//find max likes, then subtract max and get exponent
value=likes[0];
for(m=1;m<num_pows;m++)
{
if(likes[m]>value){value=likes[m];}
}
for(m=0;m<num_pows;m++){likes[m]=exp(likes[m]-value);}

//loop through all pairs of means and sds
for(j=0;j<num_means;j++)
{
//need to subtract something to avoid massive values
value=fabs(powers[0]-mvals[j]);
for(m=1;m<num_pows;m++)
{
if(fabs(powers[m]-mvals[j])<value){value=fabs(powers[m]-mvals[j]);}
}

for(k=0;k<num_sds;k++)
{
//get sums and sumsqs for all data points
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;
for(m=0;m<num_pows;m++)
{
value2=exp(-.5*(pow(powers[m]-mvals[j],2)-pow(value,2))*pow(svals[k],-2));
sum+=likes[m];sum2+=value2;sumsq+=pow(likes[m],2);sumsq2+=pow(value2,2);sumsq3+=likes[m]*value2;
}

if(omitone==0)	//get correlation across all values
{perfs[j][k]=(num_pows*sumsq3-sum*sum2)*pow(num_pows*sumsq-sum*sum,-.5)*pow(num_pows*sumsq2-sum2*sum2,-.5);}
else	//get max correlation with one value removed
{
perfs[j][k]=-2;
for(m=0;m<num_pows;m++)
{
value2=exp(-.5*(pow(powers[m]-mvals[j],2)-pow(value,2))*pow(svals[k],-2));
sum-=likes[m];sum2-=value2;sumsq-=pow(likes[m],2);sumsq2-=pow(value2,2);sumsq3-=likes[m]*value2;
value3=((num_pows-1)*sumsq3-sum*sum2)*pow((num_pows-1)*sumsq-sum*sum,-.5)*pow((num_pows-1)*sumsq2-sum2*sum2,-.5);
if(isinf(value3)!=1&&value3>perfs[j][k]){perfs[j][k]=value3;}
sum+=likes[m];sum2+=value2;sumsq+=pow(likes[m],2);sumsq2+=pow(value2,2);sumsq3+=likes[m]*value2;
}
}
}}	//end of j and k loops

//get max
topm=0;tops=0;value=perfs[0][0];
for(j=0;j<num_means;j++)
{
for(k=0;k<num_sds;k++)
{
if(perfs[j][k]>value){topm=j;tops=k;value=perfs[j][k];}
}}

//get mode and save
found=0;
for(m=1;m<num_pows;m++)
{
if(likes[m]>likes[found]){found=m;}
}

sprintf(filename,"%s.power",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Component Mode Mean SD Correlation\nPower %.4f %.4f %.4f %.4f\n", powers[found], mvals[topm], svals[tops], perfs[topm][tops]);
fclose(output);

printf("Best power is %.4f (SD %.4f), saved in %s\n\n", mvals[topm], svals[tops], filename);

free(powers);free(likes);
free(mvals);free(svals);
for(j=0;j<num_means;j++){free(perfs[j]);}free(perfs);
}	//end of 1-d search

if(count==3)	//performing 2-d search
{
mvals=malloc(sizeof(double)*num_means);
for(j=0;j<num_means;j++){mvals[j]=minmean+j*(maxmean-minmean)/(num_means-1);}

svals=malloc(sizeof(double)*num_sds);
for(k=0;k<num_sds;k++){svals[k]=(k+1)*maxsd/num_sds;}

num_pows=countrows(likefile);
powers=malloc(sizeof(double)*num_pows);
powers2=malloc(sizeof(double)*num_pows);
likes=malloc(sizeof(double)*num_pows);
read_values(likefile,powers,num_pows,NULL,1,0,0);
read_values(likefile,powers2,num_pows,NULL,2,0,0);
read_values(likefile,likes,num_pows,NULL,3,0,0);

perfs=malloc(sizeof(double*)*num_means*num_means);
for(j=0;j<num_means*num_means;j++){perfs[j]=malloc(sizeof(double)*num_sds*num_sds);}

//find max likes, then subtract max and get exponent
value=likes[0];
for(m=1;m<num_pows;m++)
{
if(likes[m]>value){value=likes[m];}
}
for(m=0;m<num_pows;m++){likes[m]=exp(likes[m]-value);}

//loop through all pairs of mean+mean and sd+sd
for(j=0;j<num_means;j++)
{
printf("Testing values %d out of %d\n", j+1, num_means);

//need to subtract something to avoid massive values
value=fabs(powers[0]-mvals[j]);
for(m=1;m<num_pows;m++)
{
if(fabs(powers[m]-mvals[j])<value){value=fabs(powers[m]-mvals[j]);}
}

for(j2=0;j2<num_means;j2++)
{
//need to subtract a second something to avoid massive values
value4=fabs(powers2[0]-mvals[j2]);
for(m=1;m<num_pows;m++)
{
if(fabs(powers2[m]-mvals[j2])<value4){value4=fabs(powers2[m]-mvals[j2]);}
}

for(k=0;k<num_sds;k++)
{
for(k2=0;k2<num_sds;k2++)
{
//get sums and sumsqs for all data points
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;
for(m=0;m<num_pows;m++)
{
value2=exp(-.5*(pow(powers[m]-mvals[j],2)-pow(value,2))*pow(svals[k],-2)-.5*(pow(powers2[m]-mvals[j2],2)-pow(value4,2))*pow(svals[k2],-2));
sum+=likes[m];sum2+=value2;sumsq+=pow(likes[m],2);sumsq2+=pow(value2,2);sumsq3+=likes[m]*value2;
}

//now get max correlation with one value removed
perfs[j+j2*num_means][k+k2*num_sds]=-2;
for(m=0;m<num_pows;m++)
{
value2=exp(-.5*(pow(powers[m]-mvals[j],2)-pow(value,2))*pow(svals[k],-2)-.5*(pow(powers2[m]-mvals[j2],2)-pow(value4,2))*pow(svals[k2],-2));
sum-=likes[m];sum2-=value2;sumsq-=pow(likes[m],2);sumsq2-=pow(value2,2);sumsq3-=likes[m]*value2;
value3=((num_pows-1)*sumsq3-sum*sum2)*pow((num_pows-1)*sumsq-sum*sum,-.5)*pow((num_pows-1)*sumsq2-sum2*sum2,-.5);
if(isinf(value3)!=1&&value3>perfs[j+j2*num_means][k+k2*num_sds]){perfs[j+j2*num_means][k+k2*num_sds]=value3;}
sum+=likes[m];sum2+=value2;sumsq+=pow(likes[m],2);sumsq2+=pow(value2,2);sumsq3+=likes[m]*value2;
}
}}	//end of k and k2 loops
}}	//end of j and j2 loops
printf("\n");

//get max
topm=0;tops=0;value=perfs[0][0];
for(j=0;j<num_means*num_means;j++)
{
for(k=0;k<num_sds*num_sds;k++)
{
if(perfs[j][k]>value){topm=j;tops=k;value=perfs[j][k];}
}}

//get mode and save
found=0;
for(m=1;m<num_pows;m++)
{
if(likes[m]>likes[found]){found=m;}
}

sprintf(filename,"%s.power",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Component Mode Mean SD Correlation\nPower1 %.4f %.4f %.4f %.4f\nPower2 %.4f %.4f %.4f %.4f\n", powers[found], mvals[topm%num_means], svals[tops%num_sds], perfs[topm][tops], powers2[found], mvals[topm/num_means], svals[tops/num_sds], perfs[topm][tops]);
fclose(output);

printf("Best powers are %.4f (SD %.4f) and %.4f (SD %.4f), saved in %s\n\n", mvals[topm%num_means], svals[tops%num_sds], mvals[topm/num_means], svals[tops/num_sds], filename);

free(powers);free(likes);
free(mvals);free(svals);
for(j=0;j<num_means*num_means;j++){free(perfs[j]);}free(perfs);
}	//end of 2-d search

///////////////////////////

