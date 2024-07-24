/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Gene-based analyses (for num_kins=0, now start by regressing Y and X on Z)

///////////////////////////

//count number of gene in this partition and set genemax to size of longest gene
count4=0;genemax=0;
for(g=0;g<num_genes;g++)
{
if(gparts[g]==partition)
{
count4++;
if(gends[g]-gstarts[g]>genemax){genemax=gends[g]-gstarts[g];}
}
}

if(mode==137||mode==138){printf("The longest gene/chunk contains %d predictors\n\n", genemax);}

//allocate variables
data_warn(num_samples_use, genemax);
data=malloc(sizeof(double)*num_samples_use*genemax);

if(mode==137)
{
anal_warn(num_samples_use, num_samples_use);
kins=malloc(sizeof(double)*num_samples_use*num_samples_use);

exps=malloc(sizeof(double)*data_length);
}
else
{
if(num_kins==0){anal_warn(num_samples_use+genemax, genemax);}
else{anal_warn(5*num_samples_use+3*genemax, genemax);}

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);
X=malloc(sizeof(double)*num_samples_use*genemax);
XTCX=malloc(sizeof(double)*genemax*genemax);
Xbeta=malloc(sizeof(double)*num_samples_use);

if(strcmp(sumsfile,"blank")!=0)
{
Xnss=malloc(sizeof(double)*genemax);
Xrhos=malloc(sizeof(double)*genemax);
Xsqs=malloc(sizeof(double)*genemax);
}

retain=malloc(sizeof(int)*genemax);
stats=malloc(sizeof(double)*6);
if(gene_perms>0)
{
order2=malloc(sizeof(int)*num_samples_use);
permlike=malloc(sizeof(double)*gene_perms);
}

if(num_kins==0)
{
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);
thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*num_samples_use);
E2=malloc(sizeof(double)*genemax);
}
else	//most allocations done within remladv.c
{vstarts=malloc(sizeof(double)*3);}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

if(mode==138||mode==140)	//fill Y and Z, solve null models etc
{
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

if(num_kins==0)	//compute detZTZ, adjust response and get YTCY
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &num_samples_use, &alpha, Z, &num_samples_use, Z, &num_samples_use, &beta, ZTZ, &num_fixed);
detZTZ=eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

reg_covar_lin(Y, Z, num_samples_use, num_fixed, 0, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);

YTCY=0;for(i=0;i<num_samples_use;i++){YTCY+=pow(Yadj[i],2);}
}
else	//get starting heritabilities, varnull and likenull (saved in vstarts)
{
printf("Solving Null Model\n");
(void)adv_reml(num_samples_use, num_fixed, 0, Y, Z, U, E, kintraces, NULL, -9999, -9999, NULL, NULL, vstarts, tol, maxiter);
printf("\n");
}

if(gene_perms>0)	//fill order2
{
for(i=0;i<num_samples_use;i++){order2[i]=i;}
}
}

//deal with progress and on-the-fly files
if(mode==137||mode==138){sprintf(filename,"%sprogress.%d", folder, partition);}
else{sprintf(filename,"%s.progress", outfile);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

if(mode==138||mode==140)
{
if(mode==138){sprintf(filename2,"%sremls.%d", folder, partition);}
else{sprintf(filename2,"%s.remls.1", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Gene_Number Gene_Name Num_Predictors Heritability SD Null_Likelihood Alt_Likelihood LRT_Stat LRT_P_Raw");
for(p=0;p<gene_perms;p++){fprintf(output2, " Perm_%d", p+1);}
fprintf(output2, "\n");

if(mode==138){sprintf(filename3,"%sprs.%d.fam", folder, partition);}
else{sprintf(filename3,"%s.prs.all.fam", outfile);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(i=0;i<num_samples_use;i++)
{
if(i%nclump==0){fprintf(output3,"%s %s 0 0 0 0\n", ids1[i], ids2[i]);}
}
fclose(output3);

if(mode==138){sprintf(filename3,"%sprs.%d.bim", folder, partition);}
else{sprintf(filename3,"%s.prs.all.bim", outfile);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(g=0;g<num_genes;g++)
{
if(gparts[g]==partition){fprintf(output3,"%d %s 0 %.0f A B\n", gchr[g], gnames[g], gbp1[g]);}
}
fclose(output3);

if(mode==138){sprintf(filename3,"%sprs.%d.sp", folder, partition);}
else{sprintf(filename3,"%s.prs.all.sp", outfile);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}

////////

count2=0;
for(g=0;g<num_genes;g++)
{
total=gends[g]-gstarts[g];

if(gparts[g]==partition)
{
if(mode==137)
{
if(count2%200==0)
{
printf("Calculating kinships for Gene/Chunk %d of %d in Partition %d\n", count2+1, count4, partition);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Calculating kinships for Gene/Chunk %d of %d in Partition %d\n", count2+1, count4, partition);
}
}
else
{
if((num_kins==0&&count2%200==0)||(num_kins==1&&count2%20==0))
{
printf("Performing REML for Gene/Chunk %d of %d in Partition %d\n", count2+1, count4, partition);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Performing REML for Gene/Chunk %d of %d in Partition %d\n", count2+1, count4, partition);
fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
fclose(output3);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
}
}

shuffle=0;	//this is how many predictors we already have data for (can be more than total)
for(j=gstarts[g];j<end;j++)
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(j-start)*num_samples_use+i];}
shuffle++;
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, gstarts[g]+shuffle, gends[g], keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data+(size_t)shuffle*num_samples_use, centres+gstarts[g]+shuffle, mults+gstarts[g]+shuffle, sqdevs+gstarts[g]+shuffle, num_samples_use, total-shuffle, missingvalue, power, 0, hwestand, weights+gstarts[g]+shuffle, 1, preds+gstarts[g]+shuffle);

if((mode==138||mode==140)&&num_kins==0&&num_fixed>1)	//regress out covariates
{reg_covar_matrix(data+(size_t)shuffle*num_samples_use, Z, num_samples_use, total-shuffle, num_fixed);}

//ready to do work
if(mode==137)	//calc-kins-genes
{
//find number of non-trivial predictors
count3=0;for(j=gstarts[g];j<gends[g];j++){count3+=(mults[j]!=-9999&&weights[j]>0);}

if(count3==0)	//gene trivial
{
printf("Warning, all %d predictors in Gene/Chunk %s are trivial", total, gnames[g]);
if(strcmp(weightsfile,"blank")!=0){printf(" or have weight zero");}
printf("\n");
}
else	//non-trivial
{
//compute kinships
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &total, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, kins, &num_samples_use);

//compute expected heritabilities, scaled so they sum to one
for(j=gstarts[g];j<gends[g];j++)
{
if(mults[j]!=-9999)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}
else{exps[j]=0;}
}

sprintf(outfile,"%sgeneships.%d", folder, g+1);
write_kins(outfile, kins, NULL, num_samples_use, ids1, ids2, 1, preds+gstarts[g], keeppreds_use+gstarts[g], centres+gstarts[g], mults+gstarts[g], weights+gstarts[g], al1+gstarts[g], al2+gstarts[g], exps+gstarts[g], total, datafile, power, kingz, kinraw, 0);
}
}
else	//calc-genes-reml - note that will have summary statistics for all predictors
{
//find number of non-trivial predictors (indexed by retain)
//if using summary statistics, will also get a weighted variance explained
if(strcmp(sumsfile,"blank")==0){count3=prune_gene(retain, gprune, data, num_samples_use, total, mults+gstarts[g], weights+gstarts[g], X, XTCX, NULL, NULL, &var2);}
else{count3=prune_gene(retain, gprune, data, num_samples_use, total, mults+gstarts[g], weights+gstarts[g], X, XTCX, rhos+gstarts[g], chis+gstarts[g], &var2);}

if(count3==0)	//gene trivial
{
printf("Warning, all %d predictors in Gene/Chunk %s are trivial", total, gnames[g]);
if(strcmp(weightsfile,"blank")!=0){printf(" or have weight zero");}
printf("\n");

fprintf(output2, "%d %s 0 NA NA NA NA NA NA", g+1, gnames[g]);
for(p=0;p<gene_perms;p++){fprintf(output2, " NA");}
fprintf(output2, "\n");

for(i=0;i<num_samples_use;i++)
{
if(i%nclump==0){fprintf(output3,"0 ");}
}
fprintf(output3,"\n");
}
else	//non-trivial
{
//fill X and, if using, summary statistics
if(strcmp(sumsfile,"blank")==0)
{sumsq=fill_gene(X, NULL, NULL, num_samples_use, total, NULL, retain, data, NULL, NULL);}
else
{sumsq=fill_gene(X, Xnss, Xrhos, num_samples_use, total, NULL, retain, data, nss+gstarts[g], rhos+gstarts[g]);}

if(num_kins==0)	//compute and decompose XTCX
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &count3, &count3, &num_samples_use, &alpha, X, &num_samples_use, X, &num_samples_use, &beta, XTCX, &count3);

if(strcmp(sumsfile,"blank")!=0)	//extract diagonals of XTCX (will not have covariates, so XTCX=XTX)
{
for(j=0;j<count3;j++){Xsqs[j]=XTCX[j+j*count3];}
}

if(shrink<1)	//deflate off diagonal terms
{
for(j=0;j<count3;j++)
{
for(j2=j+1;j2<count3;j2++){XTCX[j+j2*count3]*=shrink;XTCX[j2+j*count3]*=shrink;}
}
}
if(shrink>1)	//inflate diagonal terms
{
for(j=0;j<count3;j++){XTCX[j+j*count3]*=shrink;}
}

//decompose XTCX
lwork=-1;
dsyev_("V", "U", &count3, XTCX, &count3, E2, &wkopt, &lwork, &info);
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &count3, XTCX, &count3, E2, work, &lwork, &info);
free(work);
}

//analyse real data
if(strcmp(sumsfile,"blank")==0)	//have phenotypes
{
if(num_kins==0)
{(void)gene_reml(num_samples_use, num_fixed, Yadj, YTCY, detZTZ, X, count3, XTCX, E2, sumsq, Xbeta, NULL, NULL, NULL, stats, limit, var2, tol, maxiter);}
else
{(void)adv_reml(num_samples_use, num_fixed, 1, Y, Z, U, E, kintraces, X, count3, sumsq, Xbeta, stats, vstarts, tol, maxiter);}
}
else	//using summary statistics
{(void)gene_reml(num_samples_use, -9999, NULL, YTCY, detZTZ, X, count3, XTCX, E2, sumsq, Xbeta, Xnss, Xrhos, Xsqs, stats, limit, var2, tol, maxiter);}

if(stats[4]==-9999&&stats[0]>0.01){printf("Warning, Gene/Chunk %s is excluded because its estimated heritability (%.4f) appears to be unrealistically high; this suggests an imperfectly-matched reference panel\n", gnames[g], stats[0]);}

for(p=0;p<gene_perms;p++)
{
//when permuting, always use individual-level data
permute_int(order2,num_samples_use);
sumsq=fill_gene(X, NULL, NULL, num_samples_use, total, order2, retain, data, NULL, NULL);

if(num_kins==0)
{permlike[p]=gene_reml(num_samples_use, num_fixed, Yadj, YTCY, detZTZ, X, count3, XTCX, E2, sumsq, NULL, NULL, NULL, NULL, NULL, 0, var2, tol, maxiter);}
else
{permlike[p]=adv_reml(num_samples_use, num_fixed, 1, Y, Z, U, E, kintraces, X, count3, sumsq, NULL, NULL, vstarts, tol, maxiter);}
}	//end of p loop

//print results
fprintf(output2, "%d %s %d ", g+1, gnames[g], count3);
if(stats[1]!=-9999){fprintf(output2, "%.6f %.6f ", stats[0], stats[1]);}
else{fprintf(output2, "%.6f NA ", stats[0]);}
fprintf(output2, "%.4f %.4f %.4f %.4e", stats[2], stats[3],  stats[4], stats[5]);
for(p=0;p<gene_perms;p++){fprintf(output2, " %.4f", permlike[p]);}
fprintf(output2, "\n");

for(i=0;i<num_samples_use;i++)
{
if(i%nclump==0){fprintf(output3,"%.4f ", Xbeta[i]);}
}
fprintf(output3,"\n");
}	//end of not trivial
}	//end of mode=138 or 139

count2++;
start=gstarts[g];if(gends[g]>end){end=gends[g];}
}}	//end of in right partition and g loop
printf("\n");

fclose(output);
if(mode==138||mode==140)
{
fclose(output2);
fclose(output3);
}

if(mode==137){printf("Kinship matrices saved with stem %sgeneships\n\n", folder);}
if(mode==138){printf("REML estimates saved in %sremls.%d\n\n", folder, partition);}
if(mode==140){printf("REML estimates saved in %s.remls.1\n\n", outfile);}

free(data);
if(mode==137){free(kins);free(exps);}
else
{
free(order);free(Y);free(Z);free(X);free(XTCX);free(Xbeta);
if(strcmp(sumsfile,"blank")!=0){free(Xnss);free(Xrhos);free(Xsqs);}
free(retain);free(stats);
if(gene_perms>0){free(order2);free(permlike);}
if(num_kins==0)
{free(ZTZ);free(ZTZ2);free(ZTZ3);free(thetas);free(thetasds);free(thetapvas);free(Yadj);free(E2);}
if(num_kins==1){free(vstarts);}
}
if(binary==0){gzclose(datainputgz);}

///////////////////////////

