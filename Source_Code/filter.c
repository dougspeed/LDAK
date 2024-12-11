/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Filter based on relatedness

///////////////////////////

//read kinships
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
read_kins(kinstems[0], NULL, kins_single, 1.0, num_samples_use, ids3, 0, maxthreads);

//standardize
sum=0;sum2=0;
for(i=0;i<num_samples_use;i++)
{
sum+=kins_single[(size_t)i*num_samples_use+i];
for(i2=0;i2<num_samples_use;i2++){sum2+=kins_single[(size_t)i*num_samples_use+i2];}
}
mean=sum/num_samples_use;
mean2=sum2/num_samples_use/num_samples_use;

printf("The average diagonal kinship is %.4f, the average off-diagonal kinship is %.4f\n", mean, (sum2-sum)/num_samples_use/(num_samples_use-1));
if(kinstand==1)
{
if(fabs(mean-1)>.00005||fabs((sum2-sum)/num_samples_use/(num_samples_use-1))>.00005)
{printf("Values will be standardized so these are 1 and 0; to avoid this, use \"--kin-stand NO\"\n");}
value=pow(num_samples_use,-1);
for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<num_samples_use;i2++)
{kins_single[(size_t)i*num_samples_use+i2]=(kins_single[(size_t)i*num_samples_use+i2]-mean2)/(mean-mean2)+value;}
}
}
printf("\n");

//get off-diagonal min
min=kins_single[0];
for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<i;i2++)
{
if(kins_single[(size_t)i*num_samples_use+i2]<min){min=kins_single[(size_t)i*num_samples_use+i2];}
}
}

//for convenience, blank diagonals
for(i=0;i<num_samples_use;i++){kins_single[(size_t)i*num_samples_use+i]=min-1;}

if(minrel!=-9999)	//find related samples
{
usedids=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){usedids[i]=0;}

sprintf(filename2,"%s.pairs", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

count2=0;
for(i=0;i<num_samples_use;i++)
{
for(i2=0;i2<i;i2++)
{
if(kins_single[(size_t)i*num_samples_use+i2]>minrel)
{
usedids[i]=1;usedids[i2]=1;

fprintf(output2,"%s %s %s %s %.6f\n", ids1[i], ids2[i], ids1[i2], ids2[i2], kins_single[(size_t)i*num_samples_use+i2]);
count2++;
}
}}
fclose(output2);

count=0;for(i=0;i<num_samples_use;i++){count+=(usedids[i]);}

sprintf(filename,"%s.related", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(i=0;i<num_samples_use;i++)
{
if(usedids[i]){fprintf(output,"%s %s\n", ids1[i], ids2[i]);}
}
fclose(output);

printf("List of %d related samples saved in %s, and %d related pairs in %s\n\n", count, filename, count2, filename2);
free(usedids);
}
else	//filtering
{
if(maxrel==-9999)	//set to min kinship
{
if(min>=0)
{printf("Error, there are no kinships below zero, so you must use \"--max-rel\" to specify a cutoff\n\n");exit(1);}
printf("The smallest observed kinship is %.6f\n\n", min);
maxrel=-min;

sprintf(filename,"%s.maxrel", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"%f\n", maxrel);
fclose(output);
}

//get (off-diagonal) maxes and highs
highs=malloc(sizeof(int)*num_samples_use);
maxes=malloc(sizeof(double)*num_samples_use);
for(i=0;i<num_samples_use;i++)
{
maxes[i]=min-1;
for(i2=0;i2<num_samples_use;i2++)
{
if(kins_single[(size_t)i*num_samples_use+i2]>maxes[i])
{highs[i]=i2;maxes[i]=kins_single[(size_t)i*num_samples_use+i2];}
}
}

//now start removing
losts=malloc(sizeof(int)*num_samples_use);
for(i=0;i<num_samples_use;i++){losts[i]=0;}
value=maxrel+1;
count=0;
while(count<num_samples_use-1)
{
//which of maxes is highest
value=min-1;
for(i=0;i<num_samples_use;i++)
{
if(losts[i]==0&&maxes[i]>value){j=i;j2=highs[i];value=maxes[i];}
}
if(value<=maxrel){break;}

//will lose j and keep j2 - with no other info,  will pick based on iteration
i=j;i2=j2;
if(count%2==0){j=i;j2=i2;}
else{j=i2;j2=i;}
if(resp[i]!=missingvalue&&resp[i2]==missingvalue){j=i2;j2=i;}
if(resp[i]==missingvalue&&resp[i2]!=missingvalue){j=i;j2=i2;}

if(strcmp(respfile,"blank")!=0)
{
printf("Remove %s %s (has kinship %.6f with %s %s); Phenotypes ", ids1[j], ids2[j], value, ids1[j2], ids2[j2]);
if(resp[j]==missingvalue){printf("NA and ");}
else{printf("%.3f and ", resp[j]);}
if(resp[j2]==missingvalue){printf("NA\n");}
else{printf("%.3f\n", resp[j2]);}
}
else
{printf("Remove %s %s (has kinship %.6f with %s %s)\n", ids1[j], ids2[j], value, ids1[j2], ids2[j2]);}
losts[j]=1;
count++;

//blank kins for j, then sort highs and maxes
for(i=0;i<num_samples_use;i++)
{kins_single[(size_t)j*num_samples_use+i]=min-1;kins_single[(size_t)i*num_samples_use+j]=min-1;}
for(i=0;i<num_samples_use;i++)
{
if(losts[i]==0&&highs[i]==j)	//recompute
{
maxes[i]=min-1;
for(i2=0;i2<num_samples_use;i2++)
{
if(kins_single[(size_t)i*num_samples_use+i2]>maxes[i])
{highs[i]=i2;maxes[i]=kins_single[(size_t)i*num_samples_use+i2];}
}
}}
}	//end of while loop
if(count>0){printf("\n");}

sprintf(filename,"%s.keep", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
sprintf(filename2,"%s.lose", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
for(i=0;i<num_samples_use;i++)
{
if(losts[i]==0){fprintf(output,"%s %s\n", ids1[i], ids2[i]);}
else{fprintf(output2,"%s %s\n", ids1[i], ids2[i]);}
}
fclose(output);fclose(output2);

printf("Filtering complete: %d samples kept (%s), %d samples lost (%s)\n\n", num_samples_use-count, filename, count, filename2);

free(highs);free(losts);free(maxes);
}

free(kins_single);

///////////////////////////

