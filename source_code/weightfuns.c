/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Weighting functions

///////////////////////////

void cut_sections(int length, int *chr, double *cmbp, double section_kb, int section_length, double buffer_kb, int buffer_length, int ns, char *folder, double window_kb, int window_length, double window_cm, char *datafile, char *bimfile, int extract, int nothin)
{
int j, k, count, count2, count3, total;
int start, start2, end, end2, num_sections;
int numchr, *chrstarts, *chrlengths;
double value;

char filename[500], filename2[500];
FILE *output, *output2;


//get chromosome details
numchr=chr[length-1]+1;
chrstarts=malloc(sizeof(int)*numchr);
chrlengths=malloc(sizeof(int)*numchr);
for(j=0;j<numchr;j++){chrstarts[j]=-9999;chrlengths[j]=0;}
for(j=length-1;j>=0;j--){chrstarts[chr[j]]=j;}
for(j=0;j<length;j++){chrlengths[chr[j]]++;}

//open file which stores section details and write header lines
sprintf(filename,"%ssection.details",folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check write permission granted and a folder of this name does not exist\n\n",filename);exit(1);}

fprintf(output, "Datafiles %s %s\n", datafile, bimfile);
if(nothin==1)
{
if(extract==1){fprintf(output, "Using Filtered Predictors\n");}
else{fprintf(output, "Using All Predictors\n");}
}
else
{
if(extract==1){fprintf(output, "Filtering and Thinning\n");}
else{fprintf(output, "Using Thinned Predictors\n");}
}

if(window_kb!=-9999)
{
if(window_cm!=-9999){fprintf(output,"Window cM %.2f\n", window_cm);}
else{fprintf(output,"Window kb %.2f\n", window_kb);}
}
else
{fprintf(output,"Window Length %d\n",window_length);}

//now build sections
fprintf(output,"Section Analysed_Predictors Retained_Predictors\n");
num_sections=0;count2=0;count3=0;total=0;
for(j=0;j<numchr;j++)
{
if(chrlengths[j]>0)
{
//get number of sections on this chromosome
if(section_kb!=-9999){count=(cmbp[chrstarts[j]+chrlengths[j]-1]-cmbp[chrstarts[j]]-1)/(1000*section_kb)+1;}
else{count=(chrlengths[j]-1)/section_length+1;}

//see whether can get away with one section for the chr - set start2 and end2
if(buffer_kb!=-9999)
{
start2=0;
while(start2<chrlengths[j]-1)
{
if(cmbp[chrstarts[j]+start2+1]-cmbp[chrstarts[j]]>1000*buffer_kb){break;}
start2++;
}
end2=chrlengths[j];
while(end2>1)
{
if(cmbp[chrstarts[j]+chrlengths[j]-1]-cmbp[chrstarts[j]+end2-2]>1000*buffer_kb){break;}
end2--;
}
}
else{start2=buffer_length;end2=chrlengths[j]-buffer_length;}

//so is it small enough?
if(section_kb!=-9999)
{
if(cmbp[chrstarts[j]+end2-1]-cmbp[chrstarts[j]+start2]<=section_kb){count=1;}
}
else
{
if(end2-start2<=section_length){count=1;}
}

for(k=0;k<count;k++)
{
//set start2 and end2 - not guaranteed to be any inside
if(section_kb!=-9999)
{
value=(cmbp[chrstarts[j]+chrlengths[j]-1]-cmbp[chrstarts[j]]+1)/count;
for(start2=0;start2<chrlengths[j];start2++)
{
if(cmbp[chrstarts[j]+start2]-cmbp[chrstarts[j]]>=k*value){break;}
}
for(end2=1;end2<chrlengths[j];end2++)
{
if(cmbp[chrstarts[j]+end2]-cmbp[chrstarts[j]]>=(k+1)*value){break;}
}
}
else
{
start2=k*((chrlengths[j]-1)/count+1);
end2=(k+1)*((chrlengths[j]-1)/count+1);
if(end2>chrlengths[j]){end2=chrlengths[j];}
}

if(end2>start2)	//there are some predictors inside, so set start and end
{
if(buffer_kb!=-9999)
{
start=start2;
while(start>0)
{
if(cmbp[chrstarts[j]+start2]-cmbp[chrstarts[j]+start-1]>1000*buffer_kb){break;}
start--;
}
end=end2;
while(end<chrlengths[j])
{
if(cmbp[chrstarts[j]+end]-cmbp[chrstarts[j]+end2-1]>1000*buffer_kb){break;}
end++;
}
}
else	//using buffer_length
{
start=start2-buffer_length;
if(start<0){start=0;}
end=end2+buffer_length;
if(end>chrlengths[j]){end=chrlengths[j];}
}

fprintf(output,"%d %d-%d %d-%d\n", num_sections+1, chrstarts[j]+start+1, chrstarts[j]+end, start2-start+1, end2-start);
total+=end-start;
if(end-start>count2){count2=end-start;}
if(end-start>ns){count3++;}
num_sections++;
}	//end of some predictors inside
}	//end of k loop
}}	//end of chrlengths[j]>0 and j loop
fclose(output);

sprintf(filename2,"%ssection.number",folder);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check write permission granted and a folder of this name does not exist\n\n",filename2);exit(1);}
fprintf(output2,"%d\n", num_sections);
fclose(output2);

printf("The %d predictors have been split into %d sections (average length %d predictors, longest %d); details saved in %s\n", length, num_sections, total/num_sections, count2, filename);
if(count3==0)
{printf("Ideally each section would have length less than the number of samples (%d), and this is the case\n\n", ns);}
else
{printf("Ideally each section would have length less than the number of samples (%d); instead %d exceed this length\nConsider using \"--section-kb\", \"--section-length\" or \"--section-cm\" to reduce section sizes and perhaps \"--buffer-kb\", \"--buffer-length\" or \"--buffer-cm\" to reduce buffer sizes\n\n", ns, count3);}

printf("Section details saved in %ssection.details and %ssection.number\n\n", folder, folder);

free(chrlengths);free(chrstarts);
}	//end of cut_sections

///////////////////////////

void calc_correlations(double *cors, double *tally1, double *tally2, double *data, double missingvalue, int ns, int length, int num_subs, int **subindex, int *chr, double *cmbp, double window_kb, int window_length, double mincor, int lddecay, double decay, char *filename)
//tally1 is (effective) no of preds considered, tally2 is (effective) tagging
{
int i, i2, j, k, s, indcount;
int *usedpreds;
double sum, sumsq, mean, sd, ld, value, alpha, beta;
double *data2, *cors2;

FILE *output;


printf("Calculating correlations across %d predictors\n", length);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Calculating correlations across %d predictors\n", length);
fclose(output);

data2=malloc(sizeof(double)*length*ns);
usedpreds=malloc(sizeof(int)*length);
cors2=malloc(sizeof(double)*length*length);

//initialize things to zero
for(j=0;j<length;j++){tally1[j]=0;tally2[j]=0;usedpreds[j]=0;}
for(j=0;j<length;j++)
{
for(k=0;k<length;k++){cors[(size_t)k*length+j]=0;}
}

for(s=0;s<num_subs;s++)	//will be looping over num_subs
{
//fill data2 with standardized predictors
for(j=0;j<length;j++)
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<subindex[s][0];i++)
{
i2=subindex[s][1+i];
if(data[(size_t)j*ns+i2]!=missingvalue)
{sum+=data[(size_t)j*ns+i2];sumsq+=pow(data[(size_t)j*ns+i2],2);indcount++;}
}

if(indcount>1)	//using this predictor
{
mean=sum/indcount;
sd=pow((sumsq-indcount*mean*mean)/subindex[s][0],.5);
for(i=0;i<subindex[s][0];i++)
{
i2=subindex[s][1+i];
if(data[(size_t)j*ns+i2]==missingvalue){data2[(size_t)j*subindex[s][0]+i]=0;}
else{data2[(size_t)j*subindex[s][0]+i]=(data[(size_t)j*ns+i2]-mean)/sd;}
}
usedpreds[j]++;
}
else	//blank out
{
for(i=0;i<subindex[s][0];i++){data2[(size_t)j*subindex[s][0]+i]=0;}
}
}	//end of j loop

//get correlation
alpha=1.0/subindex[s][0];beta=0.0;
dgemm_("T", "N", &length, &length, subindex[s], &alpha, data2, subindex[s], data2, subindex[s], &beta, cors2, &length);

//put into cors if larger than current value, ld weighting if necessary
for(j=0;j<length;j++)
{
if(usedpreds[j]>0)	//non trivial so will have cor=1 (>mincor)
{
//diagonals (ld=1)
value=pow(cors2[(size_t)j*length+j],2);
tally1[j]+=1.0/num_subs;
cors[(size_t)j*length+j]=1;

//off diagonals
for(k=j+1;k<length;k++)
{
if(chr[k]!=chr[j]){printf("Chromosome error, please tell Doug\n");break;}
if(window_kb!=-9999&&cmbp[k]-cmbp[j]>1000*window_kb){break;}
if(window_length!=-9999&&k-j>window_length){break;}

if(usedpreds[k]>0)	//non-trivial and within range
{
if(lddecay==0){ld=1;}
else{ld=exp(-(cmbp[k]-cmbp[j])/1000*decay);}

value=pow(cors2[(size_t)k*length+j],2);
tally1[j]+=ld/num_subs;tally1[k]+=ld/num_subs;
if(value>=mincor)
{
if(value*ld>cors[(size_t)k*length+j]){cors[(size_t)k*length+j]=ld*value;cors[(size_t)j*length+k]=ld*value;}
}
}}	//end of using k and k loop
}}	//end of using j and j loop
}	//end of sub loop
free(data2);free(usedpreds);free(cors2);

//used to have a check whether any were very high - should not be necessary if thinning

//compute tally2, which contains the sum of squared correlations for each predictor
for(j=0;j<length;j++)
{
tally2[j]=0;
for(k=0;k<length;k++){tally2[j]+=cors[(size_t)k*length+j];}
}
}	//end of calc_correlations

///////////////////////////

#if MET==1

int solve_simplex(double *tally3, double *cors, double *infos, double *tally2, int length, int maxiter, double maxtime, char *filename, int fudge)
{
int j, k, count, total, flag, one=1;
int *row_ind, *col_ind;
double alpha, beta;
double *elements;

glp_prob *lp;
glp_smcp *parm;

FILE *output;


if(fudge==1)	//quick solve
{
for(j=0;j<length;j++)
{
if(tally2[j]>0){tally3[j]=infos[j]/tally2[j];}
else{tally3[j]=0;}
}
printf("Approximate weights computed\n\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Approximate weights computed\n");
fclose(output);

return(1);
}

////////

//solve simplex
//want to get cors |x| = infos
//have parameters y = (x,a), max -sum a
//so want A y < (info -info)
//where A = [ C -I] (<  info)
//	    [-C -I] (< -info)

printf("Attempting to solve the linear programming problem\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Attempting to solve the linear programming problem\n");
fclose(output);

//first initialise the ld workspace
lp=glp_create_prob();
glp_set_obj_dir(lp,GLP_MAX);
glp_add_rows(lp,2*length);
glp_add_cols(lp,2*length);

//set the right hand side limits (<info or <-info);
for(j=0;j<length;j++)
{glp_set_row_bnds(lp,1+j,GLP_UP, 0.0,infos[j]);}
for(j=length;j<2*length;j++)
{glp_set_row_bnds(lp,1+j,GLP_UP, 0.0,-infos[j]);}

//now set the limits on x (between 0 and 1) and a (positive) and the objective function coefficients
for(j=0;j<length;j++)
{
//glp_set_col_bnds(lp,1+j,GLP_LO, 0.0,0.0);
glp_set_col_bnds(lp,1+j,GLP_DB, 0.0,1.0);
glp_set_obj_coef(lp,1+j,0);
}
for(j=length;j<2*length;j++)
{glp_set_col_bnds(lp,1+j,GLP_LO, 0.0,0.0);
glp_set_obj_coef(lp,1+j,-1.0);}

//work out how many elements and allocate
total=2*length;
for(j=0;j<length;j++)
{
for(k=0;k<length;k++){total+=2*(cors[(size_t)k*length+j]>0);}
}
row_ind=malloc(sizeof(int)*(1+total));
col_ind=malloc(sizeof(int)*(1+total));
elements=malloc(sizeof(double)*(1+total));

//now fill up elements
count=0;

//for top left
for(j=0;j<length;j++){
for(k=0;k<length;k++){
if(cors[(size_t)k*length+j]>0)
{
row_ind[1+count]=j+1;
col_ind[1+count]=k+1;
elements[1+count]=cors[(size_t)k*length+j];
count++;
}
}}

//for bottom left
for(j=0;j<length;j++){
for(k=0;k<length;k++){
if(cors[(size_t)k*length+j]>0)
{
row_ind[1+count]=length+j+1;
col_ind[1+count]=k+1;
elements[1+count]=-cors[(size_t)k*length+j];
count++;
}
}}

//for top right - only need to fill diagonals (ones)
for(j=0;j<length;j++)
{
row_ind[1+count]=j+1;
col_ind[1+count]=length+j+1;
elements[1+count]=-1;
count++;
}

//for bottom right - only need to fill diagonals (ones)
for(j=0;j<length;j++)
{
row_ind[1+count]=length+j+1;
col_ind[1+count]=length+j+1;
elements[1+count]=-1;
count++;
}

//load these in
glp_load_matrix(lp, count, row_ind, col_ind, elements);

//initialise the parameters
parm=malloc(sizeof(glp_smcp));
glp_init_smcp(parm);
parm->it_lim=maxiter;
parm->tm_lim=maxtime*60*1000;
//parm->tol_bnd=tol;
//parm->tol_dj=tol;
//parm->tol_piv=tol;
parm->out_frq=10000;

//all loaded, so now solve the problem
flag=glp_simplex(lp,parm);
//glp_get_obj_val(lp);	//dont think we need this?

if(flag==0)	//success, so load in the weights
{
for(j=0;j<length;j++)
{
tally3[j]=glp_get_col_prim(lp,1+j);
if(tally3[j]<0){weights[j]=0;}
if(tally3[j]>1){weights[j]=1;}
}
}

glp_delete_prob(lp);
free(row_ind);free(col_ind);free(elements);free(parm);

if(flag==0)
{
printf("Linear programming problem solved successfully\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Linear programming problem solved successfully\n");
fclose(output);

return(0);
}

////////

//so if still here, must have been an error
if(flag==GLP_EITLIM)	//iteration limit reached
{
printf("\nWarning, the simplex method exceeded the iteration limit (%d), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but to (try to) avoid this approximation, rerun using \"--maxiter\" to raise the iteration limit\n\n", maxiter);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "Warning, the simplex method exceeded the iteration limit (%d), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but to (try to) avoid this approximation, rerun using \"--maxiter\" to raise the iteration limit\n", maxiter);
fclose(output);
}

if(flag==GLP_ETMLIM)	//time limit reached
{
printf("\nWarning, the simplex method exceeded the time limit (%.2f minutes), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but to (try to) avoid this approximation, rerun using \"--maxtime\" to raise the time limit\n\n", maxtime);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "Warning, the simplex method exceeded the time limit (%.2f minutes), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but to (try to) avoid this approximation, rerun using \"--maxtime\" to raise the time limit\n", maxtime);
fclose(output);
}

if(flag!=GLP_EITLIM&&flag!=GLP_ETMLIM)	//other error
{
printf("\nWarning, the simplex method failed (but not due to exceeding the iteration or time limits), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nPlease let Doug know (as this has never happened before!)\n\n");

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "Warning, the simplex method failed (but not due to exceeding the iteration or time limits), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nPlease let Doug know (as this has never happened before!)\n");
fclose(output);
}

//do quick solve
for(j=0;j<length;j++)
{
if(tally2[j]>0){tally3[j]=infos[j]/tally2[j];}
else{tally3[j]=0;}
}

return(2);
}	//end of solve_simplex

#else

int solve_simplex_new(double *tally3, double *cors, double *infos, double *tally2, int length, int maxiter, double maxtime, char *filename, int fudge)
{
int j, k, count;
int status, rval, *cmatcnt, *cmatbeg, *cmatind;
double *cmatval, *obj, *rhs, *lower, *upper, *x;
char *sense, **rownames, **colnames;

QSprob lp;

FILE *output;


if(fudge==1)	//quick solve
{
for(j=0;j<length;j++)
{
if(tally2[j]>0){tally3[j]=infos[j]/tally2[j];}
else{tally3[j]=0;}
}
printf("Approximate weights computed\n\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Approximate weights computed\n");
fclose(output);

return(1);
}

////////

//solve simplex
//want to get cors |x| = infos
//Cx < info+a and Cx >info-a with a +tive, then minimize sum(a)
//where A = [C -I] (< info)
//	    [C +I] (> info)

printf("Attempting to solve the linear programming problem\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Attempting to solve the linear programming problem\n");
fclose(output);

//obj, rhs, sense, lower, upper, colnames and rownames are fixed
obj=malloc(sizeof(double)*2*length);
rhs=malloc(sizeof(double)*2*length);
sense=malloc(sizeof(char)*2*length);
lower=malloc(sizeof(double)*2*length);
upper=malloc(sizeof(double)*2*length);

for(j=0;j<length;j++){obj[j]=0;obj[length+j]=1;}
for(j=0;j<length;j++){rhs[j]=infos[j];rhs[length+j]=infos[j];}
for(j=0;j<length;j++){sense[j]='|';sense[length+j]='G';}
for(j=0;j<2*length;j++){lower[j]=0;}
for(j=0;j<length;j++){upper[j]=1.0;upper[length+j]=QS_MAXDOUBLE;}

//work out how many nonzero elements
count=0;
for(j=0;j<length;j++)
{
for(k=0;k<length;k++)
{
if(cors[(size_t)k*length+j]>0){count++;}
}
}

//now fill working along columns, starting with C and C

cmatcnt=malloc(sizeof(int)*2*length);
cmatbeg=malloc(sizeof(int)*2*length);
cmatind=malloc(sizeof(int)*2*(count+length));
cmatval=malloc(sizeof(double)*2*(count+length));
rownames=malloc(sizeof(char*)*2*length);
for(j=0;j<2*length;j++){rownames[j]=malloc(sizeof(char)*100);sprintf(rownames[j],"X%d", j+1);}
colnames=malloc(sizeof(char*)*2*length);
for(k=0;k<2*length;k++){colnames[k]=malloc(sizeof(char)*100);sprintf(colnames[k],"Y%d", k+1);}

count=0;
for(k=0;k<length;k++)
{
cmatbeg[k]=count;
for(j=0;j<length;j++)
{
if(cors[(size_t)k*length+j]>0)
{cmatind[count]=j;cmatval[count]=cors[(size_t)k*length+j];count++;}
}
for(j=0;j<length;j++)
{
if(cors[(size_t)k*length+j]>0)
{cmatind[count]=length+j;cmatval[count]=cors[(size_t)k*length+j];count++;}
}
cmatcnt[k]=count-cmatbeg[k];
}
//now -I and I
for(k=0;k<length;k++)
{
cmatbeg[length+k]=count;
cmatind[count]=k;cmatval[count]=-1;count++;
cmatind[count]=length+k;cmatval[count]=1;count++;
cmatcnt[length+k]=count-cmatbeg[length+k];
}

lp = QSload_prob (NULL, 2*length, 2*length, cmatcnt, cmatbeg, cmatind, cmatval,
QS_MIN, obj, rhs, sense, lower, upper, (const char**) colnames, (const char**) rownames);
if (lp==(QSprob)NULL)
{
printf("Sorry, unable to load the LP problem, please tell Doug\n\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Sorry, unable to load the LP problem, please tell Doug\n");
fclose(output);
exit(1);
}

QSset_param(lp, QS_PARAM_SIMPLEX_MAX_ITERATIONS, maxiter);
QSset_param_double(lp, QS_PARAM_SIMPLEX_MAX_TIME, maxtime*60);

//all loaded, so now solve the problem
rval=QSopt_dual (lp, &status);

if(rval==0&&status==QS_LP_OPTIMAL)	//success, so load in the weights
{
x=malloc(sizeof(double)*2*length);
QSget_solution(lp, NULL, x, NULL, NULL, NULL);
for(j=0;j<length;j++)
{
tally3[j]=x[j];
if(tally3[j]<0){tally3[j]=0;}
if(tally3[j]>1){tally3[j]=1;}
}
free(x);
}

QSfree_prob(lp);
free(obj);free(rhs);free(sense);free(lower);free(upper);
free(cmatcnt);free(cmatbeg);free(cmatind);free(cmatval);
for(j=0;j<2*length;j++){free(rownames[j]);}free(rownames);
for(k=0;k<2*length;k++){free(colnames[k]);}free(colnames);

if(rval==0&&status==QS_LP_OPTIMAL)
{
printf("Linear programming problem solved successfully\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Linear programming problem solved successfully\n");
fclose(output);

return(0);
}

////////

//so if still here, must have been an error
if(status==QS_LP_ITER_LIMIT)	//iteration limit reached
{
printf("\nWarning, the simplex method exceeded the iteration limit (%d), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but to (try to) avoid this approximation, rerun using \"--maxiter\" to raise the iteration limit\n\n", maxiter);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "Warning, the simplex method exceeded the iteration limit (%d), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but to (try to) avoid this approximation, rerun using \"--maxiter\" to raise the iteration limit\n", maxiter);
fclose(output);
}

if(status==QS_LP_TIME_LIMIT)	//time limit reached
{
printf("\nWarning, the simplex method exceeded the time limit (%.2f minutes), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but to (try to) avoid this approximation, rerun using \"--maxtime\" to raise the time limit\n\n", maxtime);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "Warning, the simplex method exceeded the time limit (%.2f minutes), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but to (try to) avoid this approximation, rerun using \"--maxtime\" to raise the time limit\n", maxtime);
fclose(output);
}

if(status!=QS_LP_ITER_LIMIT&&status!=QS_LP_TIME_LIMIT)	//other error
{
printf("\nWarning, the simplex method failed (but not due to exceeding the iteration or time limits), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nPlease let Doug know (as this has never happened before!)\n\n");

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "Warning, the simplex method failed (but not due to exceeding the iteration or time limits), so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nPlease let Doug know (as this has never happened before!)\n");
fclose(output);
}

//do quick solve
for(j=0;j<length;j++)
{
tally3[j]=0;
if(tally2[j]>0){tally3[j]=infos[j]/tally2[j];}
}

return(2);
}	//end of solve_simplex_new

#endif

///////////////////////////

void join_sections(char *folder, int num_sections, int *sstarts, int *sends, int *sstarts2, int *sends2, char **allpreds, int num_preds, int *keeppreds, int num_preds_use, int *keeppreds_use, int spread)
{
int j, k, count, count2, count3;
int *replace, *usedpreds;
double *alltally1, *alltally2, *alltally3, *allinfos, *allcheck;

char readchar, readstring[500], readstring2[500];
char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500];
FILE *output, *input2, *output3, *output4, *input5, *output6;


//see which sections not present
sprintf(filename,"%ssection.missings", folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

count=0;
for(j=0;j<num_sections;j++)
{
sprintf(filename2,"%sweights.%d", folder, j+1);
if(just_check(filename2)!=0)
{
if(count<5){printf("Error, %s does not exist indicating that \"--calc-weights\" did not finish for Section %d\n", filename2, j+1);}
fprintf(output, "%d\n", j+1);
count++;
}
else
{
if(countrows(filename2)!=1+sends[j]-sstarts[j])
{printf("Error, %s should have %d rows (not %d)\n\n", filename2, 1+sends[j]-sstarts[j], countrows(filename2));exit(1);}
}
}
if(count==0){fprintf(output,"All sections complete\n");}
fclose(output);

if(count==1){printf("\nIn total, 1 section did not complete; this is listed in %s\nPlease allow it more time to complete, or if impatient restart using \"--quick-weights YES\"\n\n", filename);exit(1);}
if(count>1){printf("\nIn total, %d sections did not complete; these are listed in %s\nPlease allow them more time to complete, or if impatient restart using \"--quick-weights YES\"\n\n", count, filename);exit(1);}

////////

//read weights, keeping track of which approximated
if(num_sections==1){printf("Reading weights from 1 section\n");}
else{printf("Reading weights from %d sections\n", num_sections);}

sprintf(filename,"%ssection.approximations", folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

alltally1=malloc(sizeof(double)*num_preds);
alltally2=malloc(sizeof(double)*num_preds);
alltally3=malloc(sizeof(double)*num_preds);
allinfos=malloc(sizeof(double)*num_preds);
allcheck=malloc(sizeof(double)*num_preds);

//will be -1 if filtered out, -2 if thinned out
for(j=0;j<num_preds;j++)
{alltally1[j]=-1;alltally2[j]=-1;alltally3[j]=-1;allinfos[j]=-1;allcheck[j]=-1;}
for(j=0;j<num_preds_use;j++)
{alltally1[keeppreds[j]]=-2;alltally2[keeppreds[j]]=-2;alltally3[keeppreds[j]]=-2;
allinfos[keeppreds[j]]=-2;allcheck[keeppreds[j]]=-2;}

count=0;count2=0;count3=0;
for(j=0;j<num_sections;j++)
{
sprintf(filename2,"%sweights.%d", folder, j+1);
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n", filename2);exit(1);}

//deal with header row
if(fscanf(input2, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename2);exit(1);}
if(strcmp(readstring,"Weight")!=0&&strcmp(readstring,"Quick_Weight")!=0&&strcmp(readstring,"Approx_Weight")!=0)
{printf("Error, %s should begin \"Weight\", \"Quick_Weight\" or \"Approx_Weight\" (not %s), suggesting it was created using an old version of LDAK\n\n", filename2, readstring);exit(1);}

if(strcmp(readstring,"Approx_Weight")==0)
{
if(count<5){printf("Warning, approximate weights were computed for Section %d\n", j+1);}
fprintf(output, "%d\n", j+1);
count++;
}

if(strcmp(readstring2,"Eff_Neighbours")==0){count3++;}

//skip rest of header row, then ignore next start2
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input2, "%c", &readchar);}
for(k=0;k<sstarts2[j];k++)
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input2, "%c", &readchar);}
}

//now read the ones to keep
for(k=sstarts2[j];k<sends2[j];k++)
{
if(fscanf(input2, "%lf %lf %lf %lf %lf ", alltally3+keeppreds_use[count2], alltally1+keeppreds_use[count2], alltally2+keeppreds_use[count2], allinfos+keeppreds_use[count2], allcheck+keeppreds_use[count2])!=5)
{printf("Error reading Row %d of %s\n\n", 2+k, filename2);exit(1);}
count2++;
}

fclose(input2);
}	//end of j loop
if(count==0){fprintf(output,"All sections exact\n");}
fclose(output);

if(count==1)
{printf("\nApproximate weights were computed for 1 section, listed in %s; ", filename);}
if(count>1)
{printf("\nIn total, approximate weights were computed for %d sections, listed in %s; ", count, filename);}
if(count>0)
{printf("this is not normally a problem, but if desired, when calculating weights you can use \"--maxiter\" and \"--maxtime\" to raise the iteration and time limits\n\n");}

if(count3!=0&&count3!=num_sections)
{printf("Error, when calculating weights, \"--decay YES\" was used for %d sections, but not for %d; the same options should used for all sections\n\n", count3, num_sections-count3);exit(1);} 

////////

//print weights
sprintf(filename3,"%sweights.all", folder);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(count3==0){fprintf(output3,"Predictor Weight Neighbours Tagging Info Check\n");}
else{fprintf(output3,"Predictor Weight Eff_Neighbours Eff_Tagging Info Check\n");}

for(j=0;j<num_preds;j++)
{
if(alltally3[j]==-1){fprintf(output3, "%s NA NA NA NA NA\n", allpreds[j]);}
if(alltally3[j]==-2){fprintf(output3, "%s 0 0 0 0 0\n", allpreds[j]);}
if(alltally3[j]>=0)
{
if(count3==0){fprintf(output3, "%s %.6f %d %.3f %.3f %.3f\n", allpreds[j], alltally3[j], (int)alltally1[j], alltally2[j], allinfos[j], allcheck[j]);}
else{fprintf(output3, "%s %.6f %.3f %.3f %.3f %.3f\n", allpreds[j], alltally3[j], alltally1[j], alltally2[j], allinfos[j], allcheck[j]);}
}
}
fclose(output3);

//and weights for non-zero predictors
sprintf(filename4,"%sweights.short", folder);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(j=0;j<num_preds;j++)
{
if(alltally3[j]>0){fprintf(output4,"%s %.6f\n", allpreds[j], alltally3[j]);}
}
fclose(output4);

if(spread==1)	//get spread version
{
replace=malloc(sizeof(int)*num_preds);
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}

//open, skip header and read indexes
sprintf(filename5,"%sthin.dups", folder);
if((input5=fopen(filename5,"r"))==NULL)
{printf("Error opening %s\n\n", filename5);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input5, "%c", &readchar);}
for(j=0;j<num_preds;j++)
{
if(fscanf(input5, "%d ", replace+j)!=1){printf("Error reading Row %d of %s\n\n", j+2, filename5);exit(1);}
if(replace[j]!=-1)	//so not thinned out or trivial
{replace[j]--;usedpreds[replace[j]]++;}
}
fclose(input5);

//print out
sprintf(filename6,"%sweights.spread", folder);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
for(j=0;j<num_preds;j++)
{
if(replace[j]!=-1)
{
if(alltally3[replace[j]]>0)
{
fprintf(output6,"%s %.6f\n", allpreds[j], alltally3[replace[j]]/usedpreds[replace[j]]);}
}
}
fclose(output6);

free(replace);free(usedpreds);
}

printf("Merged weights saved in %s, with a condensed version in %s", filename3, filename4);
if(spread==1){printf(" and a spread version in %s", filename6);}
printf("\n\n");

free(alltally1);free(alltally2);free(alltally3);free(allinfos);free(allcheck);
}	//end of join_sections

///////////////////////////

int solve_quad(double *tally3, double *cors, double *infos, double *tally2, int length, char *filename, int fudge)
{
int j, k, flag;
double *cors2, *infos2;

FILE *output;


if(fudge==1)	//quick solve
{
for(j=0;j<length;j++)
{
if(tally2[j]>0){tally3[j]=infos[j]/tally2[j];}
else{tally3[j]=0;}
}
printf("Approximate weights computed\n\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Approximate weights computed\n");
fclose(output);

return(1);
}

printf("Attempting to solve the quadratic programming problem\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Attempting to solve the quadratic programming problem\n");
fclose(output);

cors2=malloc(sizeof(double)*length*length);
infos2=malloc(sizeof(double)*length);
for(j=0;j<length;j++)
{
for(k=0;k<length;k++){cors2[(size_t)k*length+j]=cors[(size_t)k*length+j];}
infos2[j]=infos[j];
}

flag=nnls_algorithm(cors2, length, length, infos2, tally3, NULL);

if(flag==0)
{
printf("Quadratic programming problem solved successfully\n");
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Quadratic programming problem solved successfully\n");
fclose(output);

free(cors2);free(infos2);

return(0);
}

////////

//so if still here, must have been an error
if(flag==1)	//iteration limit reached
{
printf("\nWarning, the quadratic solver exceeded the iteration limit, so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but if concerned you can rerun using \"--simplex YES\" to attempt to solve using the simplex method\n\n");

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "\nWarning, the quadratic solver exceeded the iteration limit, so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis is not normally a problem, but if concerned you can rerun using \"--simplex YES\" to attempt to solve using the simplex method\n\n");
fclose(output);
}

if(flag==2)	//memory allocation error
{
printf("\nWarning, the quadratic solver could not allocate sufficient memory, so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis error should not occur, so please let Doug know\n\n");

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "\nWarning, the quadratic solver could not allocate sufficient memory, so approximate solutions were found instead, equivalent to using \"--quick-weights YES\"\nThis error should not occur, so please let Doug know\n\n");
fclose(output);
}

//do quick solve
for(j=0;j<length;j++)
{
if(tally2[j]>0){tally3[j]=infos[j]/tally2[j];}
else{tally3[j]=0;}
}

return(2);
}	//end of solving quad

///////////////////////////

