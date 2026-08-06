/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//summary statistic fine-mapping

///////////////////////////

//set significance threshold (again)
value4=pow(normal_inv(cutoff/2),2);

if(num_sums3>1)  //replace focal summary statistics with mtag version
{
for(j=0;j<data_length;j++){nss[j]=nss4[j];chis[j]=chis4[j];rhos[j]=rhos4[j];}
}

//get XTY for each predictor
for(j=0;j<data_length;j++){YTdata[j]=rhos[j]*nss[j];}

//some allocations
dptrs=malloc(sizeof(struct sorting_double)*bitmax);
retain=malloc(sizeof(int)*bitmax);

printf("Fine-mapping the significant predictors\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Fine-mapping the significant predictors\n");
fclose(output);

//get number of windows with significant predictors
count4=0;
for(bit=0;bit<bittotal;bit++)
{
for(j=blockstarts[bit];j<blockends[bit];j++)
{
if(nss[j]>0&&chis[j]>value4){count4++;break;}
}
}
if(count4==0){printf("Error, there are no windows containing predictors with p-value below %.4e\n\n", cutoff);exit(1);}
else{printf("There are %d windows that contain predictors with p-value below %.4e\n\n", count4, cutoff);}

//set expected heritability of non-zero predictors
cvar=her/count4;
printf("her is %f, so causals have expected her %f\n", her, cvar);

//open results files
sprintf(filename3,"%s.clump",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor Statistic Window\n");

sprintf(filename4,"%s.reclump",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Predictor Probability Statistic Window\n");

sprintf(filename5,"%s.rereclump",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5,"Predictor Probability Statistic Window\n");

////////

count3=0;
ecount=0;
wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

count5=0;for(j=blockstarts[bit];j<blockends[bit];j++){count5+=(nss[j]>0&&chis[j]>value4);}
if(count5>0)
{
if(count3%10==0)
{
printf("Fine-mapping for Window %d of %d\n", count3+1, count4);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Fine-mapping for Window %d of %d\n", count3+1, count4);
fclose(output);

fclose(output3);
sprintf(filename3,"%s.clump",outfile);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}

fclose(output4);
sprintf(filename4,"%s.reclump",outfile);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}

fclose(output5);
sprintf(filename5,"%s.rereclump",outfile);
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}
}

//read correlations
if(num_sums3==1&&metasum==0) //easy case - only using focal correlations
{
s=sumpops[0];

//re-open cors
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//read lower triangle of focal correlations
fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

 //now extract correlations for predictors we are using (remember shrink, but no need for signs)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
cors[(size_t)k*bitlength+k]=1.0;
for(j=k+1;j<Dsizes[s][bit][1];j++)
{cors[(size_t)k*bitlength+j]=(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*shrink;}
}

fclose(input2);
}
else    //hard case - might use multiple correlations
{
//set cors to identity
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
cors[(size_t)k*bitlength+k]=1;
for(j=k+1;j<bitlength;j++){cors[(size_t)k*bitlength+j]=0;}
}

for(s=0;s<num_cors;s++)
{
if(cohers[s+bit*num_cors]>0&&Dsizes[s][bit][1]>0)	//this population contributes
{
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

fclose(input2);

//now add on off-diagonal correlations for predictors we are using (remember cohers, signs and shrink)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<Dsizes[s][bit][1];k++)
{
for(j=k+1;j<Dsizes[s][bit][1];j++)
{cors[(size_t)Duse2[s][bit][k]*bitlength+Duse2[s][bit][j]]+=cohers[s+bit*num_cors]*(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*Dsigns[s][bit][j]*Dsigns[s][bit][k]*shrink;}
}
}}
}

//order significant predictors based on test statistic
count=0;
for(j=0;j<bitlength;j++)
{
if(chis[bitstart+j]>value4){dptrs[count].value=chis[bitstart+j];dptrs[count].index=j;retain[count]=1;count++;}
}
qsort(dptrs, count, sizeof(struct sorting_double), compare_sorting_double_rev);

//now clump
for(j=0;j<count-1;j++)
{
if(retain[j]==1)    //blank correlated subsequent predictors
{
for(j2=j+1;j2<count;j2++)
{
if(dptrs[j2].index>dptrs[j].index)
{
if(pow(cors[(size_t)dptrs[j].index*bitlength+dptrs[j2].index],2)>wprune){retain[j2]=0;}
}
else
{
if(pow(cors[(size_t)dptrs[j2].index*bitlength+dptrs[j].index],2)>wprune){retain[j2]=0;}
}
}
}
}

//count and print out those remaining
count6=0;
for(j=0;j<count;j++)
{
if(retain[j]==1){fprintf(output3,"%s %.4f %d\n", preds[bitstart+dptrs[j].index], chis[bitstart+dptrs[j].index], bit+1);count6++;}
}

//prior probabilities for bit are exps[j]*escale, and should sum to count6
sum=0;for(j=bitstart;j<bitend;j++){sum+=exps[j];}
escale=pow(sum,-1)*count6;

//get starting effects for bit
#pragma omp parallel for private(j, value, value2) schedule(dynamic)
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//prior is exps[j]*escale N(0,cvar) + (1-exps[j]*escale)*N(0,0)
value=exps[j]*escale;
if(value>1){value=1.0;}
value2=1.0-value;
effs[j]=get_postmean(YTdata[j], cvar, 0, -9999, -9999, nss[j], resvar, value, value2, -9999, -9999, NULL, 4, NULL, NULL)/rjksums[j];
}
else	//set to zero
{effs[j]=0;}
}

//get average sample size
sum=0;for(j=bitstart;j<bitend;j++){sum+=nss[j];}
neff=sum/bitlength;

//will only compute rhosinvrhos when required
rhosinvrhos=-9999;

//flag=0 will indicate whether we need to do vb
flag=mcmcfinal;

if(mcmcfinal==1)	//first try fixed number of iterations, then save means
{
//use effs2 and probs2 to store sums of effect sizes and probs for each chain
for(p=0;p<num_chains;p++)
{
for(j=0;j<bitlength;j++){effs2[j+p*bitlength]=0;probs2[j+p*bitlength]=0;}
}

//put starting effects into effs3
for(p=0;p<num_chains;p++)
{
for(j=0;j<bitlength;j++){effs3[j+p*bitlength]=effs[j];}
}

//set variances to resvar
for(p=0;p<num_chains;p++){variances[p]=resvar;}

count2=0;
for(count=0;count<num_burn+num_mcmc;count++)
{
for(bitstart2=bitstart;bitstart2<bitend;bitstart2+=bitsize)
{
bitend2=bitstart2+bitsize;
if(bitend2>bitend){bitend2=bitend;}
bitlength2=bitend2-bitstart2;

//get XjXTbeta for all snps in (small) block
token=bitstart2-bitstart;
token2=bitend-bitend2;

//middle square - from bitstart2 to bitend2
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength2, &num_chains, &alpha, cors+(size_t)token*bitlength+token, &bitlength, effs3+token, &bitlength, &beta, residuals, &bitsize);

alpha=1.0;beta=1.0;
if(token>0)	//left rectangle - from zero to bitstart2
{dgemm_("N", "N", &bitlength2, &num_chains, &token, &alpha, cors+token, &bitlength, effs3, &bitlength, &beta, residuals, &bitsize);}
if(token2>0)	//right rectangle - from bitend2 to bitend
{dgemm_("T", "N", &bitlength2, &num_chains, &token2, &alpha, cors+(size_t)token*bitlength+(bitend2-bitstart), &bitlength, effs3+bitend2-bitstart, &bitlength, &beta, residuals, &bitsize);}

//sample new unifrands
for(j=0;j<128*num_chains;j++){manyrands[j]=genrand_real1();}

#pragma omp parallel for private(p,start,j,j2,sum,value,value2,postmean,postprob,postsamp,value3) schedule(dynamic)
for(p=0;p<num_chains;p++)
{
start=0;
for(j=bitstart2;j<bitend2;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j]-nss[j]*(residuals[j-bitstart2+p*bitsize]-effs3[j-bitstart+p*bitlength]);

//prior is exps[j]*escale N(0,cvar) + (1-exps[j]*escale)*N(0,0)
value=exps[j]*escale;
if(value>1){value=1.0;}
value2=1.0-value;
postmean=get_postsamp_pragma(sum, cvar, 0, -9999, -9999, nss[j], variances[p], value, value2, -9999, -9999, 4, &postprob, &postsamp, manyrands+p*128, &start);

//update residuals for remainder of block
value3=postsamp-effs3[j-bitstart+p*bitlength];
for(j2=j+1;j2<bitend2;j2++){residuals[j2-bitstart2+p*bitsize]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}

if(count>=num_burn){effs2[j-bitstart+p*bitlength]+=postmean;probs2[j-bitstart+p*bitlength]+=postprob;}
effs3[j-bitstart+p*bitlength]=postsamp;
}}	//end of j loop
}	//end of p loop
}	//end of bitstart2 loop

if(resfix==0)   //get new variances
{
//get Cbeta
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &num_chains, &alpha, cors, &bitlength, effs3, &bitlength, &beta, cors2, &bitlength);

for(p=0;p<num_chains;p++)
{
//compute t(beta) beta and t(beta) R beta
value=ddot_(&bitlength, effs3+p*bitlength, &one, effs3+p*bitlength, &one);
value2=ddot_(&bitlength, effs3+p*bitlength, &one, cors2+p*bitlength, &one);

if(value>1e-8&&value2>1e-8&&value>1.1*value2)    //resample variance from inverse gamma (provided not below 0.7)
{
if(rhosinvrhos==-9999)  //will need rhosinvrhos
{rhosinvrhos=cg_solve_lower_norm(cors, bitlength, rhos+bitstart, 1e-6);}

//need t(rhos) inv(cors) rhos - 2 t(rhos) beta + t(beta) cors beta
sum=rhosinvrhos-2*ddot_(&bitlength, effs3+p*bitlength, &one, rhos+bitstart, &one)+value2;
variances[p]=(0.5*sum*neff)/rgamma(.5*bitlength);
if(variances[p]<.7){variances[p]=1;}
}
else    //reset
{variances[p]=1;}
}
}
}	//end of count loop

//get mean effects and probs - store in first columns of effs2 and probs2 (must be careful if use pragma)
for(j=0;j<bitlength;j++)
{
sum=0;sum2=0;
for(p=0;p<num_chains;p++){sum+=effs2[j+p*bitlength];sum2+=probs2[j+p*bitlength];}
effs2[j]=sum/num_mcmc/num_chains;
probs2[j]=sum2/num_mcmc/num_chains;
}

//calculate ess for mean model - equal to 2YTXbeta - t(beta) XTXbeta (for multiple traits, this is ess of first trait only)
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &one, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);
//dsymv_("L", &bitlength, &alpha, cors, &bitlength, effs2, &one, &beta, cors2, &one);
//????? check here

ess[0]=2*ddot_(&bitlength, effs2, &one, rhos+bitstart, &one)-ddot_(&bitlength, effs2, &one, cors2, &one);

if(ess[0]<=1)	//model not terrible, so update effect sizes and probs
{
for(j=bitstart;j<bitend;j++){effs[j]=effs2[j-bitstart];probs[j]=probs2[j-bitstart];}
}
else	//model is suspect - leave effect sizes at starting estimates (probs will be zero)
{
printf("Warning, MCMC failed for Window %d, will revert to variational Bayes\n", bit+1);
flag=0;
ecount++;
wcount+=bitlength;
}
}

if(flag==0)	//either did not try mcmc, or it failed, so use vb
{
//put starting effects into effs2, and set vars to zero
for(j=0;j<bitlength;j++){effs2[j]=effs[bitstart+j];pars[j]=0;}

//set variance to resvar
variances[0]=resvar;

//calculate ess for model - equal to (2YTXbeta - t(beta) XTXbeta)/n
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &one, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);
//dsymv_("L", &bitlength, &alpha, cors, &bitlength, effs2, &one, &beta, cors2, &one);

ess[0]=2*ddot_(&bitlength, effs2, &one, rhos+bitstart, &one)-ddot_(&bitlength, effs2, &one, cors2, &one);

for(count=0;count<maxiter;count++)
{
for(bitstart2=bitstart;bitstart2<bitend;bitstart2+=bitsize)
{
bitend2=bitstart2+bitsize;
if(bitend2>bitend){bitend2=bitend;}
bitlength2=bitend2-bitstart2;

//get XjXTbeta for all snps in (small) block
token=bitstart2-bitstart;
token2=bitend-bitend2;

//middle square - from bitstart2 to bitend2
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength2, &one, &alpha, cors+(size_t)token*bitlength+token, &bitlength, effs2+token, &bitlength, &beta, residuals, &bitsize);

//dsymv_("L", &bitlength2, &alpha, cors+(size_t)token*bitlength+token, &bitlength, effs2+token, &one, &beta, residuals, &one);

alpha=1.0;beta=1.0;
if(token>0)	//left rectangle - from zero to bitstart2
{dgemm_("N", "N", &bitlength2, &one, &token, &alpha, cors+token, &bitlength, effs2, &bitlength, &beta, residuals, &bitsize);}
if(token2>0)	//right rectangle - from bitend2 to bitend
{dgemm_("T", "N", &bitlength2, &one, &token2, &alpha, cors+(size_t)token*bitlength+(bitend2-bitstart), &bitlength, effs2+bitend2-bitstart, &bitlength, &beta, residuals, &bitsize);}

for(j=bitstart2;j<bitend2;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j]-nss[j]*(residuals[j-bitstart2]-effs2[j-bitstart]);

//prior is exps[j]*escale N(0,cvar) + (1-exps[j]*escale)*N(0,0)
value=exps[j]*escale;
if(value>1){value=1.0;}
value2=1.0-value;
postmean=get_postmean(sum, cvar, 0, -9999, -9999, nss[j], variances[0], value, value2, -9999, -9999, NULL, 4, &postprob, &postvar);

//update residuals for remainder of block
value3=postmean-effs2[j-bitstart];
for(j2=j+1;j2<bitend2;j2++){residuals[j2-bitstart2]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}

effs2[j-bitstart]=postmean;
probs2[j-bitstart]=postprob;
pars[j-bitstart]=postvar;
}}	//end of j loop
}	//end of bitstart2 loop

//get Cbeta
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &one, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);

//dsymv_("L", &bitlength, &alpha, cors, &bitlength, effs2, &one, &beta, cors2, &one);
//????? check here

if(resfix==0)   //get new variance (have just got Cbeta)
{
//compute t(beta) beta and t(beta) R beta
value=ddot_(&bitlength, effs2, &one, effs2, &one);
value2=ddot_(&bitlength, effs2, &one, cors2, &one);

if(value>1e-8&&value2>1e-8&&value>1.1*value2)    //resample variance from inverse gamma (provided not below 0.7)
{
if(rhosinvrhos==-9999)  //will need rhosinvrhos
{rhosinvrhos=cg_solve_lower_norm(cors, bitlength, rhos+bitstart, 1e-6);}

//need t(rhos) inv(cors) rhos - 2 t(rhos) beta + t(beta) cors beta + sum var(beta)
sum=rhosinvrhos-2*ddot_(&bitlength, effs2, &one, rhos+bitstart, &one)+value2;
for(j=0;j<bitlength;j++){sum+=pars[j];}
variances[0]=sum*neff/bitlength;
if(variances[0]<.7){variances[0]=1;}
}
else    //reset
{variances[0]=1;}
}

//save old ess
ess2[0]=ess[0];

//compute new ess (have already computed Cbeta)
ess[0]=2*ddot_(&bitlength, effs2, &one, rhos+bitstart, &one)-ddot_(&bitlength, effs2, &one, cors2, &one);

//see whether converged
if(fabs(ess[0]-ess2[0])<tol){break;}
}	//end of count loop

if(count==maxiter){printf("Warning, Window %d did not converge within %d iterations\n", bit+1, maxiter);}

if(ess[0]<=1)	//model not terrible, so update effect sizes and probs
{
for(j=bitstart;j<bitend;j++){effs[j]=effs2[j-bitstart];probs[j]=probs2[j-bitstart];}
}
else	//model is suspect - leave effect sizes at starting estimates (probs will be zero)
{
printf("Warning, Variational Bayes failed for Window %d, will revert to starting estimates\n", bit+1);
if(mcmcfinal==0)
{
ecount++;
wcount+=bitlength;
}
}
}	//end of vb

//print out probs
for(j=bitstart;j<bitend;j++){fprintf(output4,"%s %.4e %.4f %d\n", preds[j], probs[j], chis[j], bit+1);}

//order predictors based on probabilties
for(j=0;j<bitlength;j++){dptrs[j].value=probs[bitstart+j];dptrs[j].index=j;retain[j]=1;}
qsort(dptrs, bitlength, sizeof(struct sorting_double), compare_sorting_double_rev);

//now clump, stopping when we have count6
count=0;
for(j=0;j<bitlength;j++)
{
if(retain[j]==1)    //blank correlated subsequent predictors, then print
{
for(j2=j+1;j2<bitlength;j2++)
{
if(dptrs[j2].index>dptrs[j].index)
{
if(pow(cors[(size_t)dptrs[j].index*bitlength+dptrs[j2].index],2)>wprune){retain[j2]=0;}
}
else
{
if(pow(cors[(size_t)dptrs[j2].index*bitlength+dptrs[j].index],2)>wprune){retain[j2]=0;}
}
}

fprintf(output5,"%s %.4e %.4f %d\n", preds[bitstart+dptrs[j].index], probs[bitstart+dptrs[j].index], chis[bitstart+dptrs[j].index], bit+1);
count++;
}
if(count==count6){break;}
}

count3++;
}   //end of count5>0
else    //set effects and probs to zero
{
for(j=bitstart;j<bitend;j++){effs[j]=0;probs[j]=0;}
}
}   //end of bit loop

printf("\n");
if(ecount>0){printf("%d windows failed (in total, %d predictors)\n\n", ecount, wcount);}

fclose(output3);
fclose(output4);
fclose(output5);

//some saves
sprintf(filename2,"%s.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{fprintf(output2, "%s %c %c %.6f %.4e\n", preds[j], al1[j], al2[j], centres[j], effs[j]*mults[j]);}
}
fclose(output2);

sprintf(filename3,"%s.probs",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

fprintf(output3,"Predictor Probability\n");
for(j=0;j<data_length;j++)
{fprintf(output3, "%s %.4e\n", preds[j], probs[j]);}
fclose(output3);

printf("Best model saved in %s, with posterior probabilities in %s\n\n", filename2, filename3);

free(dptrs);free(retain);

///////////////////////////

