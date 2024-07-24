/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Approximate he and pcgc

///////////////////////////

//sort out step
step=(50000/bitsize);if(step<20){step=20;}
if(step>20){step=10*(step/10);}if(step>50){step=20*(step/20);}
if(step>100){step=50*(step/50);}if(step>300){step=100*(step/100);}
if(step>1000){step=500*(step/500);}

//get labels
catlabels=malloc(sizeof(char *)*(num_parts+1));
for(q=0;q<num_parts+1;q++){catlabels[q]=malloc(sizeof(char)*500);}

if(num_parts>0)
{
if(strcmp(labfile,"blank")!=0)	//read from file
{
if((input=fopen(labfile,"r"))==NULL)
{printf("Error opening %s\n\n",labfile);exit(1);}
for(q=0;q<num_parts;q++)
{
catlabels[q]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", catlabels[q])!=1)
{printf("Error reading Row %d of %s\n\n", q+1, labfile);exit(1);}

if(strcmp(catlabels[q],"Base")==0||strcmp(catlabels[q],"Background")==0)
{printf("Error reading Row %d of %s; labels can not be either \"Base\" or \"Background\"\n\n", q+1, labfile);exit(1);}

for(q2=0;q2<q;q2++)	//check different to previous
{
if(strcmp(catlabels[q],catlabels[q2])==0)
{printf("Error reading %s; the label %s appears twice\n\n", labfile, catlabels[q]);exit(1);}
}
}
fclose(input);
}
else	//give generic names
{
if(parttype==0)
{
for(q=0;q<num_parts;q++){sprintf(catlabels[q],"Annotation_%d",q+1);}
}
else
{
for(q=0;q<num_parts;q++){sprintf(catlabels[q],"Partition_%d",q+1);}
}
}
}

if(parttype==0)	//always have a base
{strcpy(catlabels[num_parts],"Base");}
else	//add background, even if not used
{strcpy(catlabels[num_parts],"Background");}

//get the heritability model
model_warn(data_length*3/2, num_parts+1);
pindexes=malloc(sizeof(int *)*(num_parts+1));
pweights=malloc(sizeof(double *)*(num_parts+1));

addpart=get_her_model(num_parts, partpref, pindexes, pweights, data_length, keeppreds_use, num_preds, allpreds, predorder, parttype, backpart, allone);

if(addpart==2)	//were some redundant predictors, so squeeze down, and reset addpart
{
count=0;
for(j=0;j<data_length;j++)
{
flag=0;for(q=0;q<num_parts;q++){flag+=pindexes[q][j];}
if(flag>0)
{
if(count!=j)
{
chr[count]=chr[j];cmbp[count]=cmbp[j];al1[count]=al1[j];al2[count]=al2[j];
free(preds[count]);copy_string(preds,count,preds[j]);
keeppreds_use[count]=keeppreds_use[j];
centres[count]=centres[j];weights[count]=weights[j];pvalues[count]=pvalues[j];
for(q=0;q<num_parts;q++){pindexes[q][count]=pindexes[q][j];pweights[q][count]=pweights[q][j];}
}
count++;
}}
for(j=count;j<data_length;j++){free(preds[j]);}
data_length=count;

if(data_length<3)
{printf("Error, unable to continue with fewer than three predictors; come on, you can do better ;)\n\n");exit(1);}

addpart=0;
}

//find start and end of jackknife blocks
if(num_blocks>data_length){num_blocks=data_length;}
bstarts=malloc(sizeof(int)*num_blocks);
bends=malloc(sizeof(int)*num_blocks);
for(p=0;p<num_blocks;p++)
{bstarts[p]=(double)p/num_blocks*data_length;bends[p]=(double)(p+1)/num_blocks*data_length;}

//check each kinship matrix not entirely contained in a block
for(q=0;q<num_parts+addpart;q++)
{
count=0;
for(j=0;j<data_length;j++){count+=(pweights[q][j]!=0);}
for(p=0;p<num_blocks;p++)
{
count2=0;
for(j=bstarts[p];j<bends[p];j++){count2+=(pweights[q][j]!=0);}
if(count2==count){printf("Error, all %d predictors that contribute to Kinship Matrix %d are contained within Jackknife Block %d; most likely, you should remove this kinship matrix from the model\n\n", count, q+1, p+1);exit(1);}
}
}

////////

//set total, total2, total3, and bitsize
total=num_parts+addpart;
total2=num_vects+num_resps_use;
total3=total2*total;

//allocate variables
data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

value=(double)num_samples_use/1024/1024/1024*8*(2*total+total2+2*total3)+(double)bitsize/1024/1024/1024*8*(total2+total3);
if(value>1){printf("Warning, to perform the analysis requires approximately %.1f Gb\n\n", value);}

kinsums=malloc(sizeof(double)*total);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*num_samples_use);

R=malloc(sizeof(double)*num_samples_use*total2);
RTdata=malloc(sizeof(double)*total2*bitsize);
RTdata2=malloc(sizeof(double)*total3*bitsize);

MkinsD=malloc(sizeof(double)*num_samples_use*total);
MkinsD2=malloc(sizeof(double)*num_samples_use*total);
MkinsR=malloc(sizeof(double)*num_samples_use*total3);
MkinsR2=malloc(sizeof(double)*num_samples_use*total3);

KKtraces=malloc(sizeof(double)*total*total);
KKtraces2=malloc(sizeof(double)*total);
KYtraces=malloc(sizeof(double)*total);
KYtraces2=malloc(sizeof(double)*total);

exps=malloc(sizeof(double)*data_length);
ssums=malloc(sizeof(double*)*total);
for(q=0;q<total;q++){ssums[q]=malloc(sizeof(double)*(total+2));}

hers=malloc(sizeof(double)*(total+1));
hersds=malloc(sizeof(double)*(total+1));
shares=malloc(sizeof(double)*total);
sharesds=malloc(sizeof(double)*total);
cohers=malloc(sizeof(double)*total*total);
cohers2=malloc(sizeof(double)*total*total);
jacks=malloc(sizeof(double)*(total+1+total)*num_blocks);

////////

//fill Y and Z
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

//solve covariates (get covher and topher) and fill Yadj
if(mode==127)	//linear model - Yadj contains standardized residuals
{reg_covar_lin(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Yadj, 1, &covher, &topher);}
else	//logistic model - Yadj contains pcgc residuals
{reg_covar_log(Y, Z, num_samples_use, num_covars, num_tops, NULL, thetas, thetasds, thetapvas, Yadj, 0, &covher, &topher, prev, 0.001, 100);}

if(num_covars>1){printf("Proportion of variance explained by the %d covariates: %.4f\n", num_covars, covher);}
if(num_tops==1){printf("Proportion of variance explained by the top predictor: %.4f\n", topher);}
if(num_tops>1){printf("Proportion of variance explained by the %d top predictor: %.4f\n", num_tops, topher);}
if(num_fixed>1){printf("\n");}

//fill start of R with random values, then end with adjusted phenotypes
for(g=0;g<num_vects;g++)
{
for(i=0;i<num_samples_use;i++){R[i+g*num_samples_use]=rnorm_safe();}
}
for(i=0;i<num_samples_use;i++){R[i+num_vects*num_samples_use]=Yadj[i];}

//set MkinsD, MkinsD2, MkinsR and MkinsR2 to zero
for(q=0;q<total;q++)
{
for(i=0;i<num_samples_use;i++){MkinsD[(size_t)q*num_samples_use+i]=0;MkinsD2[(size_t)q*num_samples_use+i]=0;}
}
for(k=0;k<total3;k++)
{
for(i=0;i<num_samples_use;i++){MkinsR[(size_t)k*num_samples_use+i]=0;MkinsR2[(size_t)k*num_samples_use+i]=0;}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//open temp file
sprintf(filename2,"%s.temp", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

////////

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%step==0)
{
printf("Calculating traces for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating traces for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read data and standardize
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, power, 0, hwestand, weights+bitstart, 1, preds+bitstart);

if(num_fixed>1)	//adjust data for covariates
{reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}

//compute t(R) data
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &bitlength, &num_samples_use, &alpha, R, &num_samples_use, data, &num_samples_use, &beta, RTdata, &total2);

//load RTdata2, which has RTdata Wk for each vector then response, for kinship 1, then kinship 2, etc
#pragma omp parallel for private(j,q,g) schedule (static)
for(j=0;j<bitlength;j++)
{
for(q=0;q<total;q++)
{
for(g=0;g<total2;g++){RTdata2[q*total2+g+j*total3]=RTdata[g+j*total2]*pweights[q][bitstart+j];}
}
}

for(p=0;p<num_blocks;p++)	//loop through blocks, to find one(s) being used
{
if(bstarts[p]>=bitend){break;}

if(bends[p]>bitstart)	//some predictors within this block
{
start=bstarts[p];
if(start<bitstart){start=bitstart;}
end=bends[p];
if(end>bitend){end=bitend;}

//get contribution to diagonal terms
for(j=start;j<end;j++)
{
for(q=0;q<total;q++)
{
if(pweights[q][j]!=0)
{
for(i=0;i<num_samples_use;i++)
{MkinsD2[(size_t)q*num_samples_use+i]+=pow(data[(size_t)(j-bitstart)*num_samples_use+i],2)*pweights[q][j];}
}}
}

//get contribution to traces - do all kinships at once (wasteful if have non-overlapping partitions)
token=end-start;
alpha=1.0;beta=1.0;
dgemm_("N", "T", &num_samples_use, &total3, &token, &alpha, data+(size_t)(start-bitstart)*num_samples_use, &num_samples_use, RTdata2+(start-bitstart)*total3, &total3, &beta, MkinsR2, &num_samples_use);

if(bends[p]<=bitend)	//finished with this block
{
//save MkinsD2 and MkinsR2
for(q=0;q<total;q++){fwrite(MkinsD2+(size_t)q*num_samples_use, sizeof(double), num_samples_use, output2);}
for(k=0;k<total3;k++){fwrite(MkinsR2+(size_t)k*num_samples_use, sizeof(double), num_samples_use, output2);}

//add contribution of MkinsD2 to MkinsD then set to zero
for(q=0;q<total;q++)
{
for(i=0;i<num_samples_use;i++)
{MkinsD[(size_t)q*num_samples_use+i]+=MkinsD2[(size_t)q*num_samples_use+i];MkinsD2[(size_t)q*num_samples_use+i]=0;}
}

//add contribution of MkinsR2 to MkinsR then set to zero
for(k=0;k<total3;k++)
{
for(i=0;i<num_samples_use;i++)
{MkinsR[(size_t)k*num_samples_use+i]+=MkinsR2[(size_t)k*num_samples_use+i];MkinsR2[(size_t)k*num_samples_use+i]=0;}
}
}
}	//end of working on block p
}	//end of p loop
}	//end of bit loop
printf("\n");

fclose(output2);

//get the average diagonals of the kinship matrices
for(q=0;q<total;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)q*num_samples_use+i];}
kinsums[q]=sum/num_samples_use;
if(kinsums[q]==0){printf("Error, Kinship Matrix %d has trace zero, so it is not possible to continue\n\n", q+1);exit(1);}
}

//compute ssums - slightly different to the ssums used with tagging files
//ssums[q][q2] indicates fraction of category q2 that contributes to category q predictors
//ssums[q][total] is number of predictors in category
//ssums[q][total+1] is expected share of heritability from category q predictors

//first get exps based on just weights and power (ignoring pweights)
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}
else{exps[j]=0;}
}

//now multiply pweights by exps, then make sure scale to one
for(q=0;q<total;q++)
{
for(j=0;j<data_length;j++){pweights[q][j]=pweights[q][j]*exps[j];}
sum=0;for(j=0;j<data_length;j++){sum+=pweights[q][j];}
for(j=0;j<data_length;j++){pweights[q][j]=pweights[q][j]/sum;}
}

//can now set ssums
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){ssums[q][q2]=0;}
ssums[q][total]=0;

for(j=0;j<data_length;j++)	//will include trivial snps (but these will have exp zero)
{
if(pindexes[q][j]==1)	//predictor is part of q
{
for(q2=0;q2<total;q2++){ssums[q][q2]+=pweights[q2][j];}
ssums[q][total]++;
}
}

ssums[q][total+1]=ssums[q][total]/data_length;
}

//old code considered only fractions across non-trivial predictors
//sum=0;for(j=0;j<data_length;j++){sum+=(mults[j]!=-9999);}
//for(q=0;q<total;q++){ssums[q][total+1]=ssums[q][total]/sum;}

////////

//ready to estimate heritabilities, first for all predictors, then excluding a block

//open temp file for reading
if((output2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//want to save KKtraces across repetitions
sprintf(filename3,"%s.repetitions",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

for(g=0;g<num_vects;g++)
{
token=num_samples_use*total2;
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, MkinsR+(size_t)g*num_samples_use, &token, MkinsR+(size_t)g*num_samples_use, &token, &beta, KKtraces, &total);

for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){fprintf(output3,"%.6f ", KKtraces[q+q2*total]);}
fprintf(output3,"\n");
}
}
fclose(output3);

//set KKtraces to minus the cross-product of diagonals
alpha=-1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, MkinsD, &num_samples_use, MkinsD, &num_samples_use, &beta, KKtraces, &total);

//add average contribution from vectors
for(g=0;g<num_vects;g++)
{
token=num_samples_use*total2;
alpha=1.0/num_vects;beta=1.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, MkinsR+(size_t)g*num_samples_use, &token, MkinsR+(size_t)g*num_samples_use, &token, &beta, KKtraces, &total);
}

//KYtraces is phenotype terms from MkinsR x Y, minus MkinsD x Y^2
token=num_samples_use*total2;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &total, &alpha, MkinsR+(size_t)num_vects*num_samples_use, &token, Yadj, &one, &beta, KYtraces, &one);
for(q=0;q<total;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)q*num_samples_use+i]*pow(Yadj[i],2);}
KYtraces[q]-=sum;
}

//divide by traces and by two, and save in KYtraces2
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){KKtraces[q+q2*total]=KKtraces[q+q2*total]/kinsums[q]/kinsums[q2]/2;}
KYtraces[q]=KYtraces[q]/kinsums[q]/2;
KYtraces2[q]=KYtraces[q];
}

//get estimates
(void)eigen_invert(KKtraces, total, KKtraces2, 1, KYtraces, 1);
sum=0;for(q=0;q<total;q++){sum+=KYtraces[q];}
for(q=0;q<total;q++){hers[q]=KYtraces[q];}
hers[total]=sum;
for(q=0;q<total;q++){shares[q]=KYtraces[q]/sum;}

//get likelihoods and lrts
sum=0;sum2=0;for(i=0;i<num_samples_use;i++){sum+=pow(Yadj[i],2);sum2+=pow(Yadj[i],4);}
sumsq=(pow(sum,2)-sum2)/2;
smax=(size_t)num_samples_use*(num_samples_use-1)/2;
likenull=-.5*smax*(1+log(2*M_PI*sumsq/smax));
for(q=0;q<total;q++){sumsq-=KYtraces2[q]*KYtraces[q];}
like=-.5*smax*(1+log(2*M_PI*sumsq/smax));
lrtstat=2*(like-likenull);
lrtpva=erfc(pow(lrtstat,.5)*M_SQRT1_2);

//get heritability of each predictor
for(j=0;j<data_length;j++)
{
exps[j]=0;
for(q=0;q<total;q++){exps[j]+=pweights[q][j]*hers[q];}
}

////////

//get estimates excluding block p

for(p=0;p<num_blocks;p++)
{
if(p%100==0)
{
printf("Processing traces for Jackknife Block %d of %d\n", p+1, num_blocks);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Processing traces for Jackknife Block %d of %d\n", p+1, num_blocks);
fclose(output);
}

//read MkinsD2 and MkinsR2 - will already be at correct point in file
for(q=0;q<total;q++)
{
if(fread(MkinsD2+(size_t)q*num_samples_use, sizeof(double), num_samples_use, output2)!=num_samples_use)
{printf("Error reading values from %s\n\n", filename2);exit(1);}
}
for(k=0;k<total3;k++)
{
if(fread(MkinsR2+(size_t)k*num_samples_use, sizeof(double), num_samples_use, output2)!=num_samples_use)
{printf("Error reading values from %s\n\n", filename2);exit(1);}
}

//subtract MkinsD2 and MkinsR2 from MkinsD and MkinsR
for(q=0;q<total;q++)
{
for(i=0;i<num_samples_use;i++){MkinsD[(size_t)q*num_samples_use+i]-=MkinsD2[(size_t)q*num_samples_use+i];}
}
for(k=0;k<total3;k++)
{
for(i=0;i<num_samples_use;i++){MkinsR[(size_t)k*num_samples_use+i]-=MkinsR2[(size_t)k*num_samples_use+i];}
}

alpha=-1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, MkinsD, &num_samples_use, MkinsD, &num_samples_use, &beta, KKtraces, &total);

for(g=0;g<num_vects;g++)
{
token=num_samples_use*total2;
alpha=1.0/num_vects;beta=1.0;
dgemm_("T", "N", &total, &total, &num_samples_use, &alpha, MkinsR+(size_t)g*num_samples_use, &token, MkinsR+(size_t)g*num_samples_use, &token, &beta, KKtraces, &total);
}

token=num_samples_use*total2;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &total, &alpha, MkinsR+(size_t)num_vects*num_samples_use, &token, Yadj, &one, &beta, KYtraces, &one);
for(q=0;q<total;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)q*num_samples_use+i]*pow(Yadj[i],2);}
KYtraces[q]-=sum;
}

for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){KKtraces[q+q2*total]=KKtraces[q+q2*total]/kinsums[q]/kinsums[q2]/2;}
KYtraces[q]=KYtraces[q]/kinsums[q]/2;
}

(void)eigen_invert(KKtraces, total, KKtraces2, 1, KYtraces, 1);
sum=0;for(q=0;q<total;q++){sum+=KYtraces[q];}
for(q=0;q<total;q++){jacks[q+p*(total+1+total)]=KYtraces[q];}
jacks[total+p*(total+1+total)]=sum;
for(q=0;q<total;q++){jacks[total+1+q+p*(total+1+total)]=KYtraces[q]/sum;}

//add MkinsD2 and MkinsR2 back onto MkinsD and MkinsR
for(q=0;q<total;q++)
{
for(i=0;i<num_samples_use;i++){MkinsD[(size_t)q*num_samples_use+i]+=MkinsD2[(size_t)q*num_samples_use+i];}
}
for(k=0;k<total3;k++)
{
for(i=0;i<num_samples_use;i++){MkinsR[(size_t)k*num_samples_use+i]+=MkinsR2[(size_t)k*num_samples_use+i];}
}
}	//end of p loop
printf("\n");

fclose(output2);

//can delete temp file
sprintf(cmd,"rm -rf %s",filename2);
system(cmd);

sprintf(filename3,"%s.jackests", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(p=0;p<num_blocks;p++)
{
for(q=0;q<total;q++){fprintf(output3,"%f ", jacks[q+p*(total+1+total)]);}
fprintf(output3,"%f %d %d %d %d\n", jacks[total+p*(total+1+total)], bstarts[p], bends[p], chr[bstarts[p]], chr[bstarts[p]]);
}
fclose(output3);

for(q=0;q<total+1+total;q++)	//get sds
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=jacks[q+p*(total+1+total)];sumsq+=pow(jacks[q+p*(total+1+total)],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
if(q<total+1){hersds[q]=pow(var,.5);}
else{sharesds[q-total-1]=pow(var,.5);}
}

//get covariances for hers
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
sum+=jacks[q+p*(total+1+total)];
sum2+=jacks[q2+p*(total+1+total)];
sumsq+=jacks[q+p*(total+1+total)]*jacks[q2+p*(total+1+total)];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
cohers[q+q2*total]=var;
}}

//now for shares
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
sum+=jacks[total+1+q+p*(total+1+total)];
sum2+=jacks[total+1+q2+p*(total+1+total)];
sumsq+=jacks[total+1+q+p*(total+1+total)]*jacks[total+1+q2+p*(total+1+total)];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
cohers2[q+q2*total]=var;
}}

printf("Total heritability (SD): %.4f (%.4f)\n\n", hers[total], hersds[total]);

////////

//adjust for tops
for(q=0;q<total+1;q++){hers[q]*=(1-topher);hersds[q]*=(1-topher);}

//save stuff

sum=0;for(q=0;q<total;q++){sum+=kinsums[q];}

if(mode==127){sprintf(filename2,"%s.fasthe", outfile);}
else{sprintf(filename2,"%s.fastpcgc", outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Num_Kinships %d\nNum_Top_Predictors %d\nNum_Covariates %d\n", total, num_tops, num_covars);
fprintf(output2, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output2, "Total_Samples %d\nWith_Phenotypes %d\n", num_samples_use, respcounts[0]);
fprintf(output2, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output2, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output2, "Component Heritability SD Size Mega_Intensity SD\n");
for(q=0;q<total;q++){fprintf(output2, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", q+1, hers[q], hersds[q], kinsums[q], hers[q]/kinsums[q]*1000000, hersds[q]/kinsums[q]*1000000);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output2, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]+topher, hersds[total], sum, hers[total]/sum*1000000, hersds[total]/sum*1000000);
fclose(output2);

sprintf(filename3,"%s.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(mode==127){fprintf(output3, "Component Effect SD P\n");}
else{fprintf(output3, "Component Log_Odds SD P\n");}
fprintf(output3, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output3);

sprintf(filename4,"%s.share", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Share SD Expected Enrichment SD\n");
for(q=0;q<total;q++){fprintf(output4, "Share_K%d %.6f %.6f %.6f %.6f %.6f\n", q+1, shares[q], sharesds[q], kinsums[q]/sum, shares[q]/kinsums[q]*sum, sharesds[q]/kinsums[q]*sum);}
fclose(output4);

if(mode==127&&prev!=-9999)	//binary he
{
factor=get_factor(Y, num_samples_use, prev, -9999, outfile);

sprintf(filename5,"%s.fasthe.liab", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Top_Predictors %d\nNum_Covariates %d\n", total, num_tops, num_covars);
fprintf(output5, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", num_samples_use, respcounts[0]);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability SD Size Mega_Intensity SD\n");
for(q=0;q<total;q++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", q+1, hers[q]*factor, hersds[q]*factor, kinsums[q], hers[q]/kinsums[q]*1000000*factor, hersds[q]/kinsums[q]*1000000*factor);}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*factor+topher*factor, hersds[total]*factor, sum, hers[total]/sum*1000000*factor, hersds[total]/sum*1000000*factor);
fclose(output5);
}

if(mode==128)	//save marginals for pcgc
{
sprintf(filename5,"%s.fastpcgc.marginal", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Num_Kinships %d\nNum_Top_Predictors %d\nNum_Covariates %d\n", total, num_tops, num_covars);
fprintf(output5, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output5, "Total_Samples %d\nWith_Phenotypes %d\n", num_samples_use, respcounts[0]);
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
if(total==1){fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output5, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output5, "Component Heritability SD Size Mega_Intensity SD\n");
for(q=0;q<total;q++){fprintf(output5, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", q+1, hers[q]*(1-covher), hersds[q]*(1-covher), kinsums[q], hers[q]/kinsums[q]*1000000*(1-covher), hersds[q]/kinsums[q]*1000000*(1-covher));}
fprintf(output5, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output5, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*(1-covher)+topher*(1-covher), hersds[total]*(1-covher), sum, hers[total]/sum*1000000*(1-covher), hersds[total]/sum*1000000*(1-covher));
fclose(output5);
}

sprintf(filename6,"%s.cross", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
for(q=0;q<total;q++){fprintf(output6, "Her_K%d\t", q+1);}
fprintf(output6, "\n");
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++)
{fprintf(output6, "%.6f\t", cohers[q+q2*total]);}
fprintf(output6, "\n");
}
fclose(output6);

//now cats and enrichments - need to get linear combinations for each
sprintf(filename7,"%s.cats", outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
fprintf(output7, "Component Heritability SD\n");
for(q=0;q<total;q++)
{
sum=0;sum2=0;
for(q2=0;q2<total;q2++)
{
sum+=hers[q2]*ssums[q][q2];
for(q3=0;q3<total;q3++){sum2+=cohers[q2+q3*total]*ssums[q][q2]*ssums[q][q3];}
}
fprintf(output7, "Cat_K%d %.6f %.6f\n", q+1, sum, pow(sum2,.5));
}
fclose(output7);

sprintf(filename8,"%s.enrich", outfile);
if((output8=fopen(filename8,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename8);exit(1);}
fprintf(output8, "Component Share SD Expected Enrichment SD\n");
for(q=0;q<total;q++)
{
sum=0;sum2=0;
for(q2=0;q2<total;q2++)
{
sum+=shares[q2]*ssums[q][q2];
for(q3=0;q3<total;q3++){sum2+=cohers2[q2+q3*total]*ssums[q][q2]*ssums[q][q3];}
}
fprintf(output8, "Enrich_K%d %.6f %.6f %.6f %.6f %.6f\n", q+1, sum, pow(sum2,.5), ssums[q][total+1], sum/ssums[q][total+1], pow(sum2,.5)/ssums[q][total+1]);
}
fclose(output8);

sprintf(filename9,"%s.ind.hers", outfile);
if((output9=fopen(filename9,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename9);exit(1);}
fprintf(output9, "Predictor Heritability\n");
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){fprintf(output9, "%s %.6e\n", preds[j], exps[j]);}
}
fclose(output9);

sprintf(filename10,"%s.labels", outfile);
if((output10=fopen(filename10,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename10);exit(1);}
fprintf(output10, "Category Label\n");
for(q=0;q<total;q++){fprintf(output10, "K%d %s\n", q+1, catlabels[q]);}
fclose(output10);


printf("Main results saved in %s", filename2);
if(mode==127&&prev!=-9999){printf(", with a liability version saved in %s.liab", filename5);}
printf("\n\n");

////////

for(q=0;q<num_parts+1;q++){free(catlabels[q]);}free(catlabels);
for(q=0;q<num_parts+1;q++){free(pindexes[q]);free(pweights[q]);}free(pindexes);free(pweights);
free(bstarts);free(bends);
free(data);
free(kinsums);
free(order);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);free(Yadj);
free(R);free(RTdata);free(RTdata2);
free(MkinsD);free(MkinsD2);free(MkinsR);free(MkinsR2);
free(KKtraces);free(KKtraces2);free(KYtraces);free(KYtraces2);
free(exps);
for(q=0;q<total;q++){free(ssums[q]);}free(ssums);
free(hers);free(hersds);free(shares);free(sharesds);free(cohers);free(cohers2);free(jacks);
if(binary==0){gzclose(datainputgz);}

///////////////////////////

