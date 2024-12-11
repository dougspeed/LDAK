/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Functions for cutting and pruning genes and joining results from calc-gene-reml

///////////////////////////

int cut_genes(int *gchr, double *gbp1, double *gbp2, char **gnames, int *gstrand, int *gstarts, int *gends, int length, int *chr, char **preds, double *bp, double *weights, char *genefile, double chunks, int chunksbp, double up_buffer, double down_buffer, double minweight, int overlap, char *pvafile, double *pvalues, int part_length, int bychr, int type, char *folder, char *datafile, char *bimfile, int extract, int num_seek)
{
//type=0 - standard cutting, type=1 - condensing, type=2 - finding tags, type=3 - same as type=0 but folder is a filename
//for types 2 gbps might be gcms (but no consequences)
int g, j, count, count2, count3, count4, misscount, found;
int mark, mark2, *usedpreds, num_parts, *mins;
double value, left_buffer, right_buffer, *cumsum, weightsum, minpvalue, minpvalue2, *dists, *dists1, *dists2;

char filename[500], filename2[500], filename3[500];
FILE *output, *output2, *output3;


//get max number of geneschunks and set chr and start and end bp (cumsum) of each gene/chunk
//with genefile or chunksbp, include both start and end bp; with chunks exclude start and include end

if(type==0||type==1||type==3)
{
if(strcmp(genefile,"blank")!=0)	//using genefile
{count=countrows(genefile)-check_head(genefile,"Gene","Name",1);
read_genefile(genefile, gnames, gchr, gbp1, gbp2, gstrand);}

if(chunks!=-9999)	//using weights
{
//get cumsum (gbp1 and gbp2 will refer to this, rather than bp)
cumsum=malloc(sizeof(double)*length);
cumsum[0]=weights[0];
for(j=1;j<length;j++){cumsum[j]=cumsum[j-1]+weights[j];}

gchr[0]=chr[0];gbp1[0]=0;gbp2[0]=chunks/(1+overlap);count=1;
if(overlap==1){gchr[1]=chr[0];gbp1[1]=0;gbp2[1]=chunks;count=2;}
for(j=0;j<length;j++)
{
if(chr[j]==gchr[count-1])	//same chromosome
{
while(cumsum[j]>=gbp1[count-1]+chunks/(1+overlap))	//add a new chunk
{gchr[count]=chr[j];gbp1[count]=gbp1[count-1]+chunks/(1+overlap);gbp2[count]=gbp1[count]+chunks;count++;}
}
else	//new chr - make new chunk (or two)
{
gchr[count]=chr[j];gbp1[count]=cumsum[j-1];gbp2[count]=gbp1[count]+chunks/(1+overlap);count++;
if(overlap==1){gchr[count]=chr[j];gbp1[count]=cumsum[j-1];gbp2[count]=gbp1[count]+chunks;count++;}

//go back to ensure snp j included (only possible if its weight larger than chunk)
j--;
}
}
}	//end of chunks!=-9999

if(chunksbp==1)	//each basepair observed becomes a gene (any duplicates will be sorted later)
{
for(j=0;j<length;j++){gchr[j]=chr[j];gbp1[j]=bp[j];gbp2[j]=bp[j];}
count=length;
}

if(chunksbp>1)
{
gchr[0]=chr[0];gbp1[0]=1;gbp2[0]=chunksbp/(1+overlap);count=1;
if(overlap==1){gchr[1]=chr[0];gbp1[1]=1;gbp2[1]=chunksbp;count=2;}
for(j=0;j<length;j++)
{
if(chr[j]==gchr[count-1])	//same chromosome
{
while(bp[j]>gbp1[count-1]+chunksbp/(1+overlap))	//add a new chunk
{gchr[count]=chr[j];gbp1[count]=gbp1[count-1]+chunksbp/(1+overlap);gbp2[count]=gbp1[count]+chunksbp-1;count++;}
}
else	//new chr - make new chunk (or two)
{
gchr[count]=chr[j];gbp1[count]=1;gbp2[count]=chunksbp/(1+overlap);count++;
if(overlap==1){gchr[count]=chr[j];gbp1[count]=1;gbp2[count]=chunksbp;count++;}

//go back to ensure snp j included (possible if its basepair is larger than chunkbp)
j--;
}
}
}	//end of chunksbp>1

}	//end of type=0, 1 or 2

if(type==2)	//have already set details from findfile
{count=num_seek;}

////////

//now see which snps inside each gene/chunk (remember, for weights, exclude start, include end)
//gstarts and gends have been set to -9999
if(strcmp(genefile,"blank")!=0)
{printf("Finding which predictor are in each gene (this can take a while if using many genes and wide buffers)\n\n");}

mark=0;
misscount=0;count2=0;
for(g=0;g<count;g++)
{
if(strcmp(genefile,"blank")!=0&&g%100000==0){printf("Processing Gene %d out of %d\n", g+1, count);}

left_buffer=0;right_buffer=0;
if(strcmp(genefile,"blank")!=0)	//update buffers
{
if(gstrand[g]==1){left_buffer=up_buffer;right_buffer=down_buffer;}
else{left_buffer=down_buffer;right_buffer=up_buffer;}
}

//first get to chromosome at least as large
while(chr[mark]<gchr[g]&&mark<length-1){mark++;}

if(chr[mark]==gchr[g])	//have found the correct chr
{
while(mark>0)	//see if left predictor included
{
if(chr[mark-1]<gchr[g]){break;}
if(chunks!=-9999)	//using weights - stop if mark-1 outside start
{
if(cumsum[mark-1]<=gbp1[g]){break;}
}
else	//using bp - stop if mark-1 outside start
{
if(bp[mark-1]<gbp1[g]-left_buffer){break;}
}
//have not broken, so move back one (and try moving left again)
mark--;
}	//end of backtracking

while(mark<length-1)	//see if need to move right
{
if(chr[mark+1]>gchr[g]){break;}
if(chunks!=-9999)	//using weights - stop if already inside start
{
if(cumsum[mark]>gbp1[g]){break;}
}
else	//using bp - stop if already inside start
{
if(bp[mark]>=gbp1[g]-left_buffer){break;}
}
mark++;
}	//end of moving right

//see if we are within - then find mark2 and set start and end for gene
if(chunks!=-9999)	//using weights
{
if(chr[mark]==gchr[g]&&cumsum[mark]>gbp1[g]&&cumsum[mark]<=gbp2[g])
{
mark2=mark+1;
while(mark2<length)	//stop if mark2 outside end
{
if(chr[mark2]>gchr[g]||cumsum[mark2]>gbp2[g]){break;}
mark2++;
}
gstarts[g]=mark;gends[g]=mark2;
}
}
else	//using bp (or maybe cm)
{
if(chr[mark]==gchr[g]&&bp[mark]>=gbp1[g]-left_buffer&&bp[mark]<=gbp2[g]+right_buffer)
{
mark2=mark+1;
while(mark2<length)	//stop if mark2 outside end
{
if(chr[mark2]>gchr[g]||bp[mark2]>gbp2[g]+right_buffer){break;}
mark2++;
}
gstarts[g]=mark;gends[g]=mark2;
}
}

if(strcmp(genefile,"blank")!=0)	
{
if(gstarts[g]==-9999)	//on chr but didn't find gene
{
if(misscount<5)
{printf("Warning, could not find any predictors within Gene %s (Chr%d:%.0f-%.0f)\n", gnames[g], gchr[g], gbp1[g], gbp2[g]);}
misscount++;
}
else{count2++;}
}
}	//end of at correct chr
}	//end of loop through genes

if(strcmp(genefile,"blank")!=0)
{
if(count2==0){printf("Error, unable to find any of the %d genes\n\n", count);exit(1);}
if(misscount>5){printf("For the chromosomes represented, in total %d genes contained no predictors\n", misscount);}
printf("\n");
}

////////

if(strcmp(genefile,"blank")!=0)	//find distance of each predictor from genes
{
printf("Calculating how far each predictor is from each gene (this can also take a while if using many genes and wide buffers)\n");
mins=malloc(sizeof(int)*length);
dists=malloc(sizeof(double)*length);
dists1=malloc(sizeof(double)*length);	//distance up
dists2=malloc(sizeof(double)*length);	//distance down
for(j=0;j<length;j++)
{mins[j]=-1;dists[j]=up_buffer+down_buffer+1;dists1[j]=down_buffer+1;dists2[j]=up_buffer+1;}

for(g=0;g<count;g++)
{
if(g%100000==0){printf("Processing Gene %d out of %d\n", g+1, count);}
if(gstarts[g]!=-9999)
{
for(j=gstarts[g];j<gends[g];j++)
{
value=fabs(bp[j]-.5*gbp1[g]-.5*gbp2[g])-.5*(gbp2[g]-gbp1[g]);
if(value<dists[j]){dists[j]=value;mins[j]=g;}
if(value<=0){dists1[j]=0;dists2[j]=0;}

if(gstrand[g]==1)	//upstream is left, downstream is right
{
if(bp[j]<gbp1[g]&&gbp1[g]-bp[j]<dists1[j]){dists1[j]=gbp1[g]-bp[j];}
if(bp[j]>gbp2[g]&&bp[j]-gbp2[g]<dists2[j]){dists2[j]=bp[j]-gbp2[g];}
}
else	//left is downstream, right upstream
{
if(bp[j]<gbp1[g]&&gbp1[g]-bp[j]<dists2[j]){dists2[j]=gbp1[g]-bp[j];}
if(bp[j]>gbp2[g]&&bp[j]-gbp2[g]<dists1[j]){dists1[j]=bp[j]-gbp2[g];}
}
}}
}	//end of g loop
printf("\n");

if(overlap==0)	//deal with overlap between genes
{
for(g=0;g<count;g++)
{
if(gstarts[g]!=-9999)
{
mark=gstarts[g];
while(mark<gends[g])
{
if(mins[mark]==g){break;}
mark++;
}
if(mark<gends[g])	//then at least one predictor assigned to this gene, see if any more
{
mark2=mark+1;
while(mark2<gends[g])
{
if(mins[mark2]!=g){break;}
mark2++;
}
gstarts[g]=mark;gends[g]=mark2;
}
else	//no predictors assigned to gene
{gstarts[g]=-9999;gends[g]=-9999;}
}}
}
}

////////

if(chunks!=-9999)	//replace cumsums with bp
{
for(g=0;g<count;g++)
{
if(gstarts[g]!=-9999){gbp1[g]=bp[gstarts[g]];gbp2[g]=bp[gends[g]-1];}
}
}

if(chunks!=-9999||chunksbp!=-9999)	//remove any identical chunks
{
g=0;
while(g<count)
{
if(gstarts[g]!=-9999){mark=g;break;}
g++;
}

for(g=mark+1;g<count;g++)
{
if(gstarts[g]!=-9999)
{
if(gstarts[g]==gstarts[mark]&&gends[g]==gends[mark])	//same predictors, so swallow up g in mark
{gbp2[mark]=gbp2[g];gstarts[g]=-9999;gends[g]=-9999;}
else{mark=g;}
}
}
}

////////

usedpreds=malloc(sizeof(int)*length);
for(j=0;j<length;j++){usedpreds[j]=0;}

if(type==0||type==3)	//print out
{
if(type==0){sprintf(filename,"%sgenes.details",folder);}
else{sprintf(filename,"%s.genes.details",folder);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Datafiles %s %s\n", datafile, bimfile);
if(extract==0){fprintf(output, "Using All Predictors\n");}
else{fprintf(output, "Using Filtered Predictors\n");}
fprintf(output,"Gene_Name Start_Predictor End_Predictor Partition Weight Gene_Chr Gene_Start Gene_End Start_Predictor_BP End_Predictor_BP Min_Pvalue Min_Pvalue_Bonferroni\n");

found=0;	//number of genes found
count2=0;	//number of predictors in largest gene
num_parts=1;	//number of partitions
count3=0;	//number of predictors in partition
for(g=0;g<count;g++)
{
if(gstarts[g]!=-9999)
{
//get weightsum and min pvalue (even if not using pvalues)
weightsum=0;minpvalue=1;count4=0;
for(j=gstarts[g];j<gends[g];j++)
{
weightsum+=weights[j];
if(pvalues[j]!=2){count4++;}
if(pvalues[j]<minpvalue){minpvalue=pvalues[j];}
}
minpvalue2=1-pow(1-minpvalue,count4);
if(minpvalue<1e-10){minpvalue2=minpvalue*count4;}

if(weightsum>=minweight)
{
for(j=gstarts[g];j<gends[g];j++){usedpreds[j]++;}

//see if should update num_parts, and count number in this partition
if(part_length!=-9999)
{
if(count3>part_length){num_parts++;count3=0;}
count3+=gends[g]-gstarts[g];
}
if(bychr==1)
{
if(found==0){mark=chr[gstarts[g]];}
if(chr[gstarts[g]]>mark){num_parts++;count3=0;}
count3+=gends[g]-gstarts[g];
mark=chr[gstarts[g]];
}

if(chunks!=-9999||chunksbp!=-9999)	//set gene name
{sprintf(gnames[g],"Chunk_%d", found+1);}

fprintf(output, "%s %d %d %d %.4f %d %.0f %.0f %.0f %.0f ", gnames[g], gstarts[g]+1, gends[g], num_parts, weightsum, gchr[g], gbp1[g], gbp2[g], bp[gstarts[g]], bp[gends[g]-1]);
if(strcmp(pvafile,"blank")==0){fprintf(output, "NA NA\n");}
else{fprintf(output, "%.2e %.2e\n", minpvalue, minpvalue2);}

//squeeze down gnames, so that later can check for duplicates
if(g!=found)
{
free(gnames[found]);
copy_string(gnames,found,gnames[g]);
}

if(gends[g]-gstarts[g]>count2){count2=gends[g]-gstarts[g];}
found++;
}	//end of big enough
}}	//end of using g and g loop

fclose(output);

if(strcmp(genefile,"blank")!=0)	//check whether gnames unique
{
qsort(gnames,found,sizeof(char *), compare_string);

for(g=1;g<found;g++)
{
if(strcmp(gnames[g],gnames[g-1])==0)
{printf("Warning, the gene names are not unique (e.g., %s appears at least twice), so it will not be possible to clump the results of the gene-based association analysis\n\n", gnames[g]);break;}
}
}

//store the genic predictors
if(type==0){sprintf(filename2,"%sgenes.predictors.used",folder);}
else{sprintf(filename2,"%s.genes.predictors.used",folder);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

count3=0;
for(j=0;j<length;j++)
{
if(usedpreds[j]>0){fprintf(output2,"%s\n", preds[j]);count3++;}
}
fclose(output2);

if(strcmp(genefile,"blank")!=0)	//and also the distances from genes 
{
if(type==0){sprintf(filename3,"%sgenes.distances",folder);}
else{sprintf(filename3,"%s.genes.distances",folder);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor Distance Dist_Up Dist_Down\n");

for(j=0;j<length;j++)
{
if(usedpreds[j]>0){fprintf(output3,"%s %.1f %d %d\n", preds[j], dists[j], (int)dists1[j], (int)dists2[j]);}
}
fclose(output3);
}

if(found==0)
{
if(strcmp(genefile,"blank")!=0){printf("Error, none of the genes have non-zero weight\n\n");exit(1);}
else{printf("Error, none of the chunks have non-zero weight\n\n");exit(1);}
}

if(strcmp(genefile,"blank")!=0){printf("%d of the %d genes were found, ", found, count);}
else{printf("%d chunks were found, ", found);}
printf("spanning %d unique predictors ", count3);
if(num_parts==1){printf("and divided into 1 partition; ");}
else{printf("and divided into %d partitions; ", num_parts);}
if(strcmp(genefile,"blank")!=0){printf("the longest gene contains %d predictors\n\n", count2);}
else{printf("the longest chunk contains %d predictors\n\n", count2);}

printf("Details saved in %s and %s", filename, filename2);
if(strcmp(genefile,"blank")!=0){printf(", with distances of predictors from genes in %s", filename3);}
printf("\n\n");
}	//end of type=0 or 3

free(usedpreds);
if(chunks!=-9999){free(cumsum);}
if(strcmp(genefile,"blank")!=0){free(mins);free(dists);free(dists1);free(dists2);}

return(0);
}	//end of cut_genes

///////////////////////////

int prune_gene(int *retain, double gprune, double *data, int ns, int length, double *mults, double *weights, double *data2, double *cors, double *rhos, double *chis, double *var2)
//save space by using pointers for data2 and cors
//just below is code for checking test stat consistency
{
int i, j, j2, count;
double sum, sumsq, mean, var, value, alpha, beta;


//first set retain to exclude trivial, then count how many remain
for(j=0;j<length;j++){retain[j]=(mults[j]!=-9999&&weights[j]>0);}
count=0;for(j=0;j<length;j++){count+=retain[j];}

if(gprune<=1||rhos!=NULL)	//get correlations
{
//get correlations

for(j=0;j<length;j++)
{
if(retain[j]==1)
{
sum=0;sumsq=0;
for(i=0;i<ns;i++){sum+=data[(size_t)j*ns+i];sumsq+=pow(data[(size_t)j*ns+i],2);}
mean=sum/ns;
var=sumsq/ns-pow(mean,2);
value=pow(var,-.5);
for(i=0;i<ns;i++){data2[(size_t)j*ns+i]=(data[(size_t)j*ns+i]-mean)*value;}
}
else
{
for(i=0;i<ns;i++){data2[(size_t)j*ns+i]=0;}
}
}

alpha=1.0/ns;beta=0.0;
dgemm_("T", "N", &length, &length, &ns, &alpha, data2, &ns, data2, &ns, &beta, cors, &length);
}

if(gprune<=1)	//prune and update count
{
for(j=0;j<length;j++)
{
if(retain[j]==1)
{
for(j2=j+1;j2<length;j2++)
{
if(retain[j2]==1)
{
if(pow(cors[(size_t)j*length+j2],2)>gprune)	//always remove second in each pair
{retain[j2]=0;}
}}
}}

count=0;for(j=0;j<length;j++){count+=retain[j];}
}

if(rhos!=NULL)	//get weighted variance explained
{
*var2=0;
for(j=0;j<length;j++)
{
if(retain[j]==1)
{
sum=0;
for(j2=0;j2<length;j2++)
{
if(retain[j2]==1){sum+=pow(cors[(size_t)j*length+j2],2);}
}
*var2+=pow(rhos[j],2)/sum;
}}
}
else{*var2=1;}

return(count);
}	//end of prune_gene

/*
//some codes for test consistency

if(pclump==3)	//test each predictor for consistency
{
part1=malloc(sizeof(double)*length);
part1b=malloc(sizeof(double)*length);
part2=malloc(sizeof(double)*length*length);
part2b=malloc(sizeof(double)*length);
part3=malloc(sizeof(double)*length*length);
dptrs=malloc(sizeof(struct sorting_double)*length);
stats=malloc(sizeof(double)*length);

for(j=0;j<length;j++)
{
if(retain[j]==1)
{
//find (up to) ten most correlated snps
count=0;
for(j2=0;j2<length;j2++)
{
if(retain[j2]==1&&j2!=j&&pow(cors[(size_t)j*length+j2],2)>.1)
{dptrs[count].value=-pow(cors[(size_t)j*length+j2],2);dptrs[count].index=j2;count++;}
}

//printf("for %d there are %d\n", j+1, count);
if(count>10){count=10;}
if(count>0)	//can continue
{
qsort(dptrs, count, sizeof(struct sorting_double), compare_sorting_double);
for(j2=0;j2<count;j2++)
{
part1[j2]=cors[(size_t)j*length+dptrs[j2].index];
part1b[j2]=cors[(size_t)j*length+dptrs[j2].index];
for(j3=0;j3<count;j3++){part2[j2*count+j3]=cors[(size_t)dptrs[j2].index*length+dptrs[j3].index];}
part3[j2]=pow(chis[dptrs[j2].index],.5);
if(rhos[dptrs[j2].index]<0){part3[j2]=-part3[j2];}
}

//compute T=inv(part2)*part1
(void)eigen_invert(part2, count, part2b, 1, part1, 1);

//estimate is part3 * T
mean=0;for(j2=0;j2<count;j2++){mean+=part3[j2]*part1[j2];}

//var is 1 - part1 * T
var=1;for(j2=0;j2<count;j2++){mean+=part1b[j2]*part1[j2];}

//test is (orig-mean)/var
value=pow(chis[j],.5);
if(rhos[j]<0){value=-value;}
stats[j]=pow(value-mean,2)/var;
}
}}
}

{
dptrs=malloc(sizeof(struct sorting_double)*length);

for(j=0;j<length;j++)
{
if(retain[j]==1)
{
//find (up to) ten most correlated snps
count=0;
for(j2=0;j2<length;j2++)
{
if(retain[j2]==1&&j2!=j&&pow(cors[(size_t)j*length+j2],2)>.1)
{dptrs[count].value=-pow(cors[(size_t)j*length+j2],2);dptrs[count].index=j2;count++;}
}

printf("for %d there are %d\n", j+1, count);
if(count>10){count=10;}
if(count>0)	//can continue
{
qsort(dptrs, count, sizeof(struct sorting_double), compare_sorting_double);
for(j2=0;j2<count;j2++)
{
part1[j2]=cors[(size_t)j*length+dptrs[j2].index];
part1b[j2]=cors[(size_t)j*length+dptrs[j2].index];
for(j3=0;j3<count;j3++){part2[j2*count+j3]=cors[(size_t)dptrs[j2].index*length+dptrs[j3].index];}
part3[j2]=pow(chis[dptrs[j2].index],.5);
if(rhos[dptrs[j2].index]<0){part3[j2]=-part3[j2];}
}

//compute T=inv(part2)*part1
(void)eigen_invert(part2, count, part2b, 1, part1, 1);

//estimate is part3 * T
mean=0;for(j2=0;j2<count;j2++){mean+=part3[j2]*part1[j2];}

//var is 1 - part1 * T
var=1;for(j2=0;j2<count;j2++){mean+=part1b[j2]*part1[j2];}

//test is (orig-mean)/var
value=pow(chis[j],.5);
if(rhos[j]<0){value=-value;}
value2=pow(value-mean,2)/var;

if(value2>10.827)	//corresponds to p<0.001
{printf("lost %d stat %f\n", j+1, value2);retain[j]=0;}
}
}}

free(dptrs);
}
*/

////////

double fill_gene(double *X, double *Xnss, double *Xrhos, int ns, int length, int *order, int *retain, double *data, double *nss, double *rhos)
{
int i, i2, j, count;
double sumsq;


count=0;
sumsq=0;
for(j=0;j<length;j++)
{
if(retain[j]==1)
{
for(i=0;i<ns;i++)
{
i2=i;
if(order!=NULL){i2=order[i];}
X[i+count*ns]=data[(size_t)j*ns+i2];
sumsq+=pow(X[i+count*ns],2);
}
if(Xnss!=NULL){Xnss[count]=nss[j];Xrhos[count]=rhos[j];}
count++;
}
}

return(sumsq/ns);
}

///////////////////////////

int join_genes_reml(char *folder, int num_genes, char **gnames, int *gstarts, int *gends, int *gparts, int *gchr, double *gbp1, double *gbp2, double gamp, double gam1, double gam2, double cut1, double cut2, int num_preds_use, char **preds, int *keeppreds_use, int type, int magma)
//type=0 - normal, type=1 - same, except folder is a filename
{
int g, g2, i, j, p, q, count, count2, count3, info;
int num_parts, gene_perms, start, end, new, new2, *usedpreds;
double frac, alpha, beta, startbp, endbp;
double *gdetails, *pvalues, *permstore, *maxes1, *maxes2, *doubletemp;

double sum, sumsq, value, mean, var, *data, *cors;

int readint;
char readchar, readstring[500];
char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], cmd[500];
FILE *input, *output2, *output3, *output4, *output5;


//set num_parts
num_parts=gparts[num_genes-1];

if(type==0)	//check all files present
{
count=0;
for(q=0;q<num_parts;q++)
{
sprintf(filename, "%sremls.%d", folder, q+1);
if(just_check(filename)!=0)
{printf("Error, %s does not exist\n", filename);count++;exit(1);}

sprintf(filename, "%sprs.%d.sp", folder, q+1);
if(just_check(filename)!=0)
{printf("Error, %s does not exist\n", filename);count++;exit(1);}

sprintf(filename, "%sprs.%d.bim", folder, q+1);
if(just_check(filename)!=0)
{printf("Error, %s does not exist\n", filename);count++;exit(1);}

sprintf(filename, "%sprs.%d.fam", folder, q+1);
if(just_check(filename)!=0)
{printf("Error, %s does not exist\n", filename);count++;exit(1);}
}

if(count>0)
{
if(num_parts==1)
{printf("Please first analyse using \"--calc-genes-reml\" then try again\n\n");exit(1);}
if(count==num_parts)
{printf("Please analyse each partition using \"--calc-genes-reml\" then try again\n\n");exit(1);}
if(count==1)
{printf("Please analyse the missing partition using \"--calc-genes-reml\" then try again\n\n");exit(1);}
printf("Please analyse the missing partitions using \"--calc-genes-reml\" then try again\n\n");exit(1);
}
}

//get gene_perms from first partition
if(type==0){sprintf(filename, "%sremls.1", folder);}
else{sprintf(filename, "%s.remls.1", folder);}
if(countcols(filename)<9)
{printf("Error, %s should have at least nine columns (not %d)\n\n", filename, countcols(filename));exit(1);}
gene_perms=countcols(filename)-9;

//check dimensions of all partitions
for(q=0;q<num_parts;q++)
{
if(type==0){sprintf(filename, "%sremls.%d", folder, q+1);}
else{sprintf(filename, "%s.remls.%d", folder, q+1);}
if(countrows_min(filename,2)<2)
{printf("Error, %s should have at least two rows\n\n", filename);exit(1);}
if(countcols(filename)<9)
{printf("Error, %s should have at least nine columns (not %d)\n\n", filename, countcols(filename));exit(1);}

if(countcols(filename)!=9+gene_perms)
{
if(gene_perms==1){printf("Error, 1 permutation was used for Partition 1, ");}
else{printf("Error, %d permutations were used for Partition 1, ", gene_perms);}
printf("but %d for Partition %d\n\n", countcols(filename)-9, q+1);exit(1);
}

count=0;
for(g=0;g<num_genes;g++){count+=(gparts[g]==q+1);}
if(countrows(filename)!=count+1)
{printf("Error, %s should have %d rows (not %d)\n\n", filename, count+1, countrows(filename));exit(1);}
}

////////

//read in for each partition in turn

gdetails=malloc(sizeof(double)*num_genes*7);
pvalues=malloc(sizeof(double)*num_genes);
if(gene_perms>0){permstore=malloc(sizeof(double)*num_genes*gene_perms);}
maxes1=malloc(sizeof(double)*(gene_perms+1));
maxes2=malloc(sizeof(double)*(gene_perms+1));

g=0;
for(q=0;q<num_parts;q++)
{
if(type==0){sprintf(filename, "%sremls.%d", folder, q+1);}
else{sprintf(filename, "%s.remls.%d", folder, q+1);}
count=countrows(filename)-1;

//open and skip header row
if((input=fopen(filename,"r"))==NULL)
{printf("Error open %s\n\n", filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

for(j=0;j<count;j++)
{
//read first two columns
if(fscanf(input, "%d %s ", &readint, readstring)!=2)
{printf("Error reading start of Row %d of %s\n\n", j+2, filename);exit(1);}
if(readint!=g+1){printf("Error reading %s; Row %d should correspond to Gene/Chunk %d (not %d)\n\n", filename, j+1, g+1, readint);exit(1);}
if(strcmp(readstring,gnames[g])!=0)
{printf("Error reading %s; Row %d should correspond to %s (not %s)\n\n", filename, j+1, gnames[g], readstring);exit(1);}

//read next seven columns - note that values can be NA
for(p=0;p<7;p++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s, suggesting the file has been changed since creation with \"--calc-genes-reml\"\n\n", p+3, j+1, filename);exit(1);}
if(strcmp(readstring,"NA")==0){gdetails[g+p*num_genes]=-9999;}
else{gdetails[g+p*num_genes]=atof(readstring);}
}

for(p=0;p<gene_perms;p++)
{
/*
if(fscanf(input, "%[0-9eE.+-]%c", readstring, &readchar)==2){permstore[g+p*num_genes]=atof(readstring);}
else
{
if(fscanf(input, "%s%c", readstring, &readchar)!=2)
{printf("Error reading Permutation %d from Row %d of %s, suggesting the file has been changed since creation with \"--calc-genes-reml\"\n\n", p+1, j+1, filename);exit(1);}
if(strcmp(readstring,"NA")!=0)
{printf("Error reading %s; Permutation %d of Row %d is unrecognisable (%s), suggesting the file has been changed since creation with \"--calc-genes-reml\"\n\n\n\n", filename, p+1, j+1, readstring);exit(1);}
permstore[g+p*num_genes]=-9999;
}
*/

if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Permutation %d from Row %d of %s, suggesting the file has been changed since creation with \"--calc-genes-reml\"\n\n", p+1, j+1, filename);exit(1);}
if(strcmp(readstring,"NA")==0)
{permstore[g+p*num_genes]=-9999;}
else
{
if(sscanf(readstring, "%lf ", permstore+g+p*num_genes)!=1)
{printf("Error reading Permutation %d from Row %d of %s, suggesting the file has been changed since creation with \"--calc-genes-reml\"\n\n", p+1, j+1, filename);exit(1);}
}
}

g++;
}
fclose(input);
}	//end of partition loop

////////

if(gamp!=-9999)	//values provided for frac, alpha and beta
{frac=gamp;alpha=gam1;beta=gam2;}
else	//will compute from results
{
if(gene_perms==0)	//use alpha=beta=0.5 then set frac based on observed fraction
{
alpha=0.5;beta=0.5;
count=0;count2=0;
for(g=0;g<num_genes;g++)
{
if(gdetails[g]>0)	//has length>0
{
if(gdetails[g+5*num_genes]>0){count2++;}
count++;
}}
if(count>0){frac=(double)count2/count;}
else{frac=0.5;}
}
else	//instead set frac, alpha and beta based on permutations
{
doubletemp=malloc(sizeof(double)*num_genes*gene_perms);
count=0;count2=0;
for(g=0;g<num_genes;g++)
{
if(gdetails[g]>0)	//has length>0
{
for(p=0;p<gene_perms;p++)
{
if(permstore[g+p*num_genes]>0){doubletemp[count2]=permstore[g+p*num_genes];count2++;}
count++;
}
}}
if(count>0){frac=(double)count2/count;}
else{frac=0.5;}
if(count2>0)
{
qsort(doubletemp, count2, sizeof(double), compare_double_rev);
count3=count/10;
if(count3>count2){count3=count2;}
info=nm_optim(&alpha, &beta, doubletemp, count2, count3);
if(info!=0){printf("Error %d running nm_optim, please tell Doug\n\n", info);exit(1);}
}
else{alpha=0.5;beta=0.5;}
free(doubletemp);
}	//end of gene_perms>0
}

//want to store max LRT and min p-value
for(p=0;p<gene_perms+1;p++){maxes1[p]=0;maxes2[p]=2;}

//each new pvalue is frac*(1-Igamma(stat*beta, alpha))
for(g=0;g<num_genes;g++)	
{
if(gdetails[g]>0)	//has length>0
{
pvalues[g]=frac*(1-gamain(gdetails[g+5*num_genes]*beta, alpha, &info));
if(pvalues[g]==0||info>2)	//typically means lrt too large, so use pchisq(2*beta*stat,1)
{pvalues[g]=erfc(pow(2*beta*gdetails[g+5*num_genes],.5)*M_SQRT1_2);}
if(gdetails[g+5*num_genes]<=0){pvalues[g]=frac+.5*(1-frac);}	//non-positive LRT

if(gdetails[g+5*num_genes]>maxes1[gene_perms])
{maxes1[gene_perms]=gdetails[g+5*num_genes];maxes2[gene_perms]=pvalues[g];}

for(p=0;p<gene_perms;p++)
{
if(permstore[g+p*num_genes]>maxes1[p])	//highest seen so far
{
maxes1[p]=permstore[g+p*num_genes];
maxes2[p]=frac*(1-gamain(maxes1[p]*beta, alpha, &info));
if(maxes2[p]==0||info>2)	//typically means lrt too large, so use pchisq(2*beta*stat,1)
{maxes2[p]=erfc(pow(2*beta*maxes1[p],.5)*M_SQRT1_2);}
}
}
}
else	//so have length=0
{pvalues[g]=2;}
}

////////

//now save

if(type==0){sprintf(filename2,"%sremls.all", folder);}
else{sprintf(filename2,"%s.remls.all", folder);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Gene_Name Gene_Chr Gene_Start Gene_End Length Heritability SE Null_Likelihood Alt_Likelihood LRT_Stat LRT_P_Raw LRT_P_Perm\n");

for(g=0;g<num_genes;g++)
{
fprintf(output2, "%s %d %.0f %.0f %d ", gnames[g], gchr[g], gbp1[g], gbp2[g], (int)gdetails[g]);
if(gdetails[g]>0)
{
if(gdetails[g+2*num_genes]!=-9999){fprintf(output2, "%.6f %.6f ", gdetails[g+num_genes], gdetails[g+2*num_genes]);}
else{fprintf(output2, "%.6f NA ", gdetails[g+num_genes]);}
fprintf(output2, "%.4f %.4f %.4f %.4e %.4e\n", gdetails[g+3*num_genes], gdetails[g+4*num_genes], gdetails[g+5*num_genes], gdetails[g+6*num_genes], pvalues[g]);
}
else{fprintf(output2, "NA NA NA NA NA NA NA\n");}
}
fclose(output2);

if(type==0){sprintf(filename3,"%smaximums",folder);}
else{sprintf(filename3,"%s.maximums",folder);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
fprintf(output3, "Analysis Maximum_LRT LRT_P_Perm\nReal %.6f %.4e\n", maxes1[gene_perms], maxes2[gene_perms]);
for(p=0;p<gene_perms;p++){fprintf(output3,"Permutation_%d %.4f %.4e\n", p+1, maxes1[p], maxes2[p]);}
fclose(output3);

if(type==0){sprintf(filename4,"%sgammas",folder);}
else{sprintf(filename4,"%s.gammas",folder);}
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename4);exit(1);}
fprintf(output4, "Fraction %.6f\nShape %.6f\nRate %.6f\n", frac, alpha, beta);
fclose(output4);

//save p-values
if(type==0){sprintf(filename5,"%sprs.all.pvalues",folder);}
else{sprintf(filename5,"%s.prs.all.pvalues",folder);}
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename5);exit(1);}
for(g=0;g<num_genes;g++){fprintf(output5,"%s %.4e\n", gnames[g], pvalues[g]);}
fclose(output5);

if(type==0)	//copy and concatenate prs files
{
sprintf(cmd, "cp %sprs.1.sp %sprs.all.sp", folder, folder);
system(cmd);
sprintf(cmd, "cp %sprs.1.bim %sprs.all.bim", folder, folder);
system(cmd);
sprintf(cmd, "cp %sprs.1.fam %sprs.all.fam", folder, folder);
system(cmd);

for(q=1;q<num_parts;q++)
{
sprintf(cmd, "cat %sprs.%d.sp >> %sprs.all.sp", folder, q+1, folder);
system(cmd);
sprintf(cmd, "cat %sprs.%d.bim >> %sprs.all.bim", folder, q+1, folder);
system(cmd);
}
}

printf("Joined results saved in %s", filename2);
if(gene_perms==0){printf(" (Column 12 contains genomic-corrected p-values), ");}
if(gene_perms==1){printf(" (Column 12 contains p-values based on 1 permutation), ");}
if(gene_perms>1){printf(" (Column 12 contains p-values based on %d permutations), ", gene_perms);}
printf("with details of the assumed null distribution saved in %s, and maximums in %s\n\n", filename4, filename3);

if(type==0){printf("To clump results, use \"--thin-tops\" with \"--sp %sprs.all\", \"--pvalues %s\" and \"--SNP-data NO\" (you should also use \"--cutoff\", \"--window-prune\" and \"--window-kb\" to specify clumping parameters)\n\n", folder, filename5);}
else{printf("To clump results, use \"--thin-tops\" with \"--sp %s.prs.all\", \"--pvalues %s\" and \"--SNP-data NO\" (you should also use \"--cutoff\", \"--window-prune\" and \"--window-kb\" to specify clumping parameters)\n\n", folder, filename5);}

////////

if(magma==1)	//get correlations
{
if(type==0){sprintf(filename,"%sprs.all.fam",folder);}
else{sprintf(filename,"%s.prs.all.fam",folder);}
count=countrows(filename);

printf("Computing correlations between estimated genetic contributions\n\n");
anal_warn(count+num_genes,num_genes);

data=malloc(sizeof(double)*count*num_genes);
cors=malloc(sizeof(double)*num_genes*num_genes);

if(type==0){sprintf(filename,"%sprs.all.sp",folder);}
else{sprintf(filename,"%s.prs.all.sp",folder);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error open %s\n\n", filename);exit(1);}
for(g=0;g<num_genes;g++)
{
for(i=0;i<count;i++)
{
if(fscanf(input, "%lf ", data+(size_t)g*count+i)!=1)
{printf("Error reading Row %d of %s\n\n", g+1, filename);exit(1);}
}
}
fclose(input);

for(g=0;g<num_genes;g++)
{
sum=0;sumsq=0;
for(i=0;i<count;i++){sum+=data[(size_t)g*count+i];sumsq+=pow(data[(size_t)g*count+i],2);}
mean=sum/count;
var=sumsq/count-pow(mean,2);

if(var>0){value=pow(var,-.5);}
else{value=0;}
for(i=0;i<count;i++){data[(size_t)g*count+i]=(data[(size_t)g*count+i]-mean)*value;}
}

//compute correlations
alpha=1.0/count;beta=0.0;
dgemm_("T", "N", &num_genes, &num_genes, &count, &alpha, data, &count, data, &count, &beta, cors, &num_genes);

if(type==0){sprintf(filename2, "%sprs.all.magma", folder);}
else{sprintf(filename2, "%s.prs.all.magma", folder);}
if((output2=fopen(filename2, "w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}

for(g=0;g<num_genes;g++)
{
value=normal_inv(1-pvalues[g]);
fprintf(output2, "%s %d %d %d %d %d %d %d %.6f", gnames[g], gchr[g], (int)gbp1[g]+1, (int)gbp2[g], gends[g]-gstarts[g], gends[g]-gstarts[g], 1000, 1000, value);

for(g2=0;g2<g;g2++)
{
if(gchr[g]==gchr[g2]){fprintf(output2, " %.6f", cors[(size_t)g*num_genes+g2]);}
}
fprintf(output2, "\n");
}

printf("Within chromosome correlations saved in %s\n\n", filename2);

free(data);free(cors);
}

////////

if(cut1!=-9999)	//prepare regions for MultiBLUP
{
printf("Identifying genes/chunks with P < %f and merging with those with P < %f\n", cut1, cut2);

if(type==0){sprintf(filename2, "%sregion.details", folder);}
else{sprintf(filename2, "%s.region.details", folder);}
if((output2=fopen(filename2, "w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
fprintf(output2, "Number Chromosome Start_BP End_BP Start_Predictor End_Predictor\n");

if(gene_perms==0)	//replace gc corrected pvalues with raw
{
for(g=0;g<num_genes;g++){pvalues[g]=gdetails[g+6*num_genes];}
}

usedpreds=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){usedpreds[j]=0;}

count=0;count2=0;
for(g=0;g<num_genes;g++)
{
if(pvalues[g]!=-9999&&pvalues[g]<=cut1)	//gene is significant
{
start=gstarts[g];startbp=gbp1[g];
end=gends[g];endbp=gbp2[g];
pvalues[g]=-9999;

//move left - consider any gene which overlaps or touches start
new=g-1;new2=g;
while(new>=0)
{
if(gchr[new]!=gchr[g]){break;}	//(definitely) out of range
if(gends[new]>=start&&pvalues[new]!=-9999)	//touches, so consider including
{
if(gdetails[new]==0)	//include automatically - hopefully few of these
{
if(gstarts[new]<start){start=gstarts[new];startbp=gbp1[new];}
new2=new;
}
else	//include if significant
{
if(pvalues[new]<=cut2)
{
if(gstarts[new]<start){start=gstarts[new];startbp=gbp1[new];}
new2=new;
}}
}
new--;
}
for(j=new2;j<g;j++){pvalues[j]=-9999;}

//move right - consider any gene which overlaps or touches end
new=g+1;new2=g;
while(new<num_genes)
{
if(gchr[new]!=gchr[g]){break;}	//(definitely) out of range
if(gstarts[new]<=end&&pvalues[new]!=-9999)	//touches, so consider including
{
if(gdetails[new]==0)	//include automatically - hopefully few of these
{
if(gends[new]>end){end=gends[new];endbp=gbp2[new];}
new2=new;
}
else	//include if extends and significant
{
if(pvalues[new]<=cut2)
{
if(gends[new]>end){end=gends[new];endbp=gbp2[new];}
new2=new;
}}
}
new++;
}
for(j=g+1;j<=new2;j++){pvalues[j]=-9999;}

fprintf(output2, "%d %d %.0f %.0f %d %d\n",  count+1, gchr[g], startbp, endbp, start+1, end);

if(type==0){sprintf(filename3, "%sregion.%d", folder, count+1);}
else{sprintf(filename3, "%s.region.%d", folder, count+1);}
if((output3=fopen(filename3, "w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
for(j=start;j<end;j++){fprintf(output3,"%s\n",preds[j]);}
fclose(output3);

for(j=start;j<end;j++){count2++;usedpreds[j]++;}
count++;
}	//end of gene g significant
}	//end of g loop
fclose(output2);

count3=0;for(j=0;j<num_preds_use;j++){count3+=(usedpreds[j]>0);}

if(type==0){sprintf(filename3, "%sregion.0", folder);}
else{sprintf(filename3, "%s.region.0", folder);}
if((output3=fopen(filename3, "w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
for(j=0;j<num_preds_use;j++)
{
if(usedpreds[j]==0){fprintf(output3,"%s\n",preds[j]);}
}
fclose(output3);

if(type==0){sprintf(filename4, "%sregion.number", folder);}
else{sprintf(filename4, "%s.region.number", folder);}
if((output4=fopen(filename4, "w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename4);exit(1);}
fprintf(output4, "%d\n", count);
fclose(output4);

if(count==0){printf("There are no significant regions\n\n");}
if(count==1)
{
if(type==0){printf("There is one region, spanning %d predictors; details are stored with prefix %sregion.\n\n", count2, folder);}
else{printf("There is one region, spanning %d predictors; details are stored with prefix %s.region.\n\n", count2, folder);}
}
if(count>1)
{
printf("There are %d regions, spanning %d predictors", count, count2);
if(count3<count2){printf(" (%d unique)", count3);}
if(type==0){printf("; details are stored with prefix %sregion.\n\n", folder);}
else{printf("; details are stored with prefix %s.region.\n\n", folder);}
}

free(usedpreds);
}	//end of cut1!=-9999

free(gdetails);free(pvalues);
if(gene_perms>0){free(permstore);}
free(maxes1);free(maxes2);

return(0);
}	//end of join_genes_reml

///////////////////////////

