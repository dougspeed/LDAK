/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Some kinship functions (mainly requested by others)

///////////////////////////

if(mode==166)	//truncate-grm
{
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);

read_kins(kinstems[0], NULL, kins_single, 1.0, num_samples_use, ids3, 0, maxthreads);

scount=0;
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<i2;i++)
{
if(kins_single[(size_t)i2*num_samples_use+i]<cutoff)
{kins_single[(size_t)i2*num_samples_use+i]=0;kins_single[(size_t)i*num_samples_use+i2]=0;scount++;}
}
}
if(scount==0){printf("Error, there are no off-diagonal allelic correlations less than %.6f\n\n", cutoff);exit(1);}
printf("There are %zu off-diagonal allelic correlations less than %.6f\n\n", scount, cutoff);

sprintf(filename,"%s.trun",outfile);
write_kins(filename, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 3);

printf("Truncated kinship matrix saved with stem %s\n\n", filename);

scount=0;
for(i2=0;i2<num_samples_use;i2++)
{
kins_single[(size_t)i2*num_samples_use+i2]=1;
for(i=0;i<i2;i++)
{
if(kins_single[(size_t)i2*num_samples_use+i]>0)
{kins_single[(size_t)i2*num_samples_use+i]=1;kins_single[(size_t)i*num_samples_use+i2]=1;scount++;}
}
}
if(scount==0){printf("Error, there are no off-diagonal allelic correlations at least %.6f\n\n", cutoff);exit(1);}
printf("There are %zu off-diagonal allelic correlations at least %.6f\n\n", scount, cutoff);

sprintf(filename2,"%s.cenv",outfile);
write_kins(filename2, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 3);

printf("Common environment kinship matrix saved with stem %s\n\n", filename2);

free(kins_single);
}

////////

if(mode==167)	//pca-grm
{
kin_warn(1, num_samples_use, 0, 0);
kins=malloc(sizeof(double*)*num_samples_use*num_samples_use);

anal_warn(num_samples_use, num_samples_use);
kins2=malloc(sizeof(double*)*num_samples_use*num_samples_use);
U=malloc(sizeof(double)*num_samples_use*axes);
U2=malloc(sizeof(double)*num_samples_use*axes);
E=malloc(sizeof(double)*axes);

//read kins, putting a copy in kins2
read_kins(kinstems[0], kins, NULL, 1.0, num_samples_use, ids3, 0, maxthreads);

for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){kins2[(size_t)i2*num_samples_use+i]=kins[(size_t)i2*num_samples_use+i];}
}

printf("Performing the decomposition; this can take a while\n\n");

vl=-1;vu=-1;	//vl and vu required but not used
count=num_samples_use-axes+1;
iwork=malloc(sizeof(int)*5*num_samples_use);
ifail=malloc(sizeof(int)*num_samples_use);
lwork=-1;
dsyevx_("V", "I", "U", &num_samples_use, kins, &num_samples_use, &vl, &vu, &count, &num_samples_use, &tol, &found, E, U, &num_samples_use, &wkopt, &lwork, iwork, ifail, &info);
if(info!=0){printf("PCA error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

decomp_warn(lwork);
work=malloc(sizeof(double)*lwork);
dsyevx_("V", "I", "U", &num_samples_use, kins, &num_samples_use, &vl, &vu, &count, &num_samples_use, &tol, &found, E, U, &num_samples_use, work, &lwork, iwork, ifail, &info);
if(info!=0){printf("PCA error 2; please tell Doug\n\n");exit(1);}

free(iwork);free(ifail);free(work);

printf("Will first calculate the PCA kinships, then subtract these from the original kinships\n\n");

//U2 contains eigenvectors multiplied by eigenvalues
for(j=0;j<axes;j++)
{
for(i=0;i<num_samples_use;i++){U2[(size_t)j*num_samples+i]=U[(size_t)j*num_samples+i]*E[j];}
}

//compute and save the PC kinships
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_samples_use, &axes, &alpha, U2, &num_samples_use, U, &num_samples_use, &beta, kins, &num_samples_use);
sprintf(filename,"%s.top",outfile);
write_kins(filename, kins, NULL, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 3);

//compute and save the remainder kinships
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){kins2[(size_t)i2*num_samples_use+i]-=kins[(size_t)i2*num_samples_use+i];}
}
sprintf(filename2,"%s.bottom",outfile);
write_kins(filename2, kins2, NULL, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 3);

printf("Kinship matrices saved with stems %s and %s\n\n", filename, filename2);

free(kins);free(kins2);free(U);free(U2);free(E);
}

////////

if(mode==168)	//square-grm
{
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);

read_kins(kinstems[0], NULL, kins_single, 1.0, num_samples_use, ids3, 0, maxthreads);

for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++)
{kins_single[(size_t)i2*num_samples_use+i]=pow(kins_single[(size_t)i2*num_samples_use+i],2);}
}

write_kins(outfile, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 3);

printf("Kinship matrix saved with stem %s\n\n", outfile);

free(kins_single);
}

////////

if(mode==169||mode==170)	//gxemm
{
if((mode==169&&noise==0)||mode==170)
{
kin_warn(1, num_samples_use, 0, 1);
kins_single=malloc(sizeof(float)*num_samples_use*num_samples_use);
}
else
{kins_single=malloc(sizeof(float)*num_samples_use);}

Z=malloc(sizeof(double)*num_samples_use*num_envs);

if(mode==170)
{
anal_warn(num_samples_use/2, num_samples_use);
kins_single2=malloc(sizeof(float)*num_samples_use*num_samples_use);
}

if((mode==169&&noise==0)||mode==170)
{read_kins(kinstems[0], NULL, kins_single, 1.0, num_samples_use, ids3, 0, maxthreads);}

for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_envs;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

if(mode==169)	//compute one new kinship (either genetic or noise)
{
printf("Constructing the kinship matrix\n\n");

if(noise==0)	//compute K ZZT
{
#pragma omp parallel for private(i2,i,value) schedule (static)
for(i2=0;i2<num_samples_use;i2++)
{
if(num_samples_use>40000&&i2<num_samples_use/maxthreads&&i2%20000==0&&i2*maxthreads<num_samples_use)
{printf("Processing Sample %d out of %d\n", i2*maxthreads+1, num_samples_use);}
for(i=0;i<num_samples_use;i++)
{
value=0;
for(j=0;j<num_envs;j++){value+=Z[i+j*num_samples_use]*Z[i2+j*num_samples_use];}
if(noise==0){kins_single[(size_t)i2*num_samples_use+i]*=value;}
else{kins_single[(size_t)i2*num_samples_use+i]=value;}
}
}

//normalize kinships
sum=0;sum2=0;
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){sum+=kins_single[(size_t)i2*num_samples_use+i];}
sum2+=kins_single[(size_t)i2*num_samples_use+i2];
}
mean=sum/num_samples_use/num_samples_use;
mean2=sum2/num_samples_use-mean;

for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++)
{kins_single[(size_t)i2*num_samples_use+i]=(kins_single[(size_t)i2*num_samples_use+i]-mean)/mean2;}
}

sprintf(filename,"%s.iid",outfile);
write_kins(filename, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 3);

sprintf(filename2,"%s.list",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "%s\n%s\n", kinstems[0], filename);
fclose(output2);

printf("IID kinship matrix saved with stem %s\n\nThe next step is to use \"--reml\", \"--he\" or \"--pcgc\" with \"--mgrm %s; remember to use \"--enviro\", to include the environmental variables you used in this step as fixed effects (you will also need to use \"--kinship-details NO\")\n\n", filename, filename2);
}
else	//compute I ZZT - save diagonal in start of kins_single
{
for(i=0;i<num_samples_use;i++)
{
value=0;for(j=0;j<num_envs;j++){value+=pow(Z[i+j*num_samples_use],2);}
kins_single[(size_t)i]=value;
}

//scale diagonal
sum=0;
for(i=0;i<num_samples_use;i++){sum+=kins_single[i];}
mean=sum/num_samples_use;
for(i=0;i<num_samples_use;i++){kins_single[i]=kins_single[i]/mean;}

sprintf(filename,"%s.noise",outfile);
write_kins(filename, NULL, kins_single, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 5);

printf("Noise kinship matrix saved with stem %s\n\n", filename);
}
}
else	//need pairs of kinships
{
for(j=0;j<num_envs;j++)
{
printf("Constructing kinship matrix pair %d out of %d\n", j+1, num_envs);

//compute K ZjZjT
#pragma omp parallel for private(i2,i) schedule (static)
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++)
{kins_single2[(size_t)i2*num_samples_use+i]=kins_single[(size_t)i2*num_samples_use+i]*Z[i+j*num_samples_use]*Z[i2+j*num_samples_use];}
}

//normalize kinships
sum=0;sum2=0;
for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++){sum+=kins_single2[(size_t)i2*num_samples_use+i];}
sum2+=kins_single2[(size_t)i2*num_samples_use+i2];
}
mean=sum/num_samples_use/num_samples_use;
mean2=sum2/num_samples_use-mean;

for(i2=0;i2<num_samples_use;i2++)
{
for(i=0;i<num_samples_use;i++)
{kins_single2[(size_t)i2*num_samples_use+i]=(kins_single2[(size_t)i2*num_samples_use+i]-mean)/mean2;}
}

sprintf(filename,"%s.genetic%d",outfile, j+1);
write_kins(filename, NULL, kins_single2, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 3);

//compute I ZjZjT - save diagonal in start of kins_single2
for(i=0;i<num_samples_use;i++){kins_single2[i]=pow(Z[i+j*num_samples_use],2);}

//scale diagonal
sum=0;
for(i=0;i<num_samples_use;i++){sum+=kins_single2[i];}
mean=sum/num_samples_use;
for(i=0;i<num_samples_use;i++){kins_single2[i]=kins_single2[i]/mean;}

sprintf(filename,"%s.noise%d",outfile, j+1);
write_kins(filename, NULL, kins_single2, num_samples_use, ids1, ids2, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, -9999, NULL, -9999, kingz, kinraw, 5);
}
printf("\n");

if(discenv==0)
{
sprintf(filename2,"%s.list",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "%s\n", kinstems[0]);
for(j=0;j<num_envs;j++){fprintf(output2, "%s.genetic%d\n", outfile, j+1);}
for(j=0;j<num_envs;j++){fprintf(output2, "%s.noise%d\n", outfile, j+1);}
fclose(output2);

printf("Kinship matrices saved with stems %s.genetic and %s.noise\n\nThe next step is to use \"--reml\", \"--he\" or \"--pcgc\" with \"--mgrm %s, to perform REML, Haseman Elston Regression or PCGC Regression; remember to use \"--enviro\", to include the environmental variables you used in this step as fixed effects (you will also need to use \"--kinship-details NO\")\n\n", outfile, outfile, filename2);
}
else
{
sprintf(filename2,"%s.list",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "%s\n", kinstems[0]);
for(j=0;j<num_envs;j++){fprintf(output2, "%s.genetic%d\n", outfile, j+1);}
for(j=0;j<num_envs-1;j++){fprintf(output2, "%s.noise%d\n", outfile, j+1);}
fclose(output2);

printf("Kinship matrices saved with stems %s.genetic and %s.noise\n\nThe next step is to use \"--reml\", \"--he\" or \"--pcgc\" with \"--mgrm %s, to perform REML, Haseman Elston Regression or PCGC Regression; remember to use both \"--enviro\", to include the environmental variables you used in this step as fixed effects, and \"--subgroups YES\" (you will also need to use \"--kinship-details NO\")\n\n", outfile, outfile, filename2);
}
}

free(kins_single);
free(Z);
if(mode==170){free(kins_single2);}
}

//////////////////////////

