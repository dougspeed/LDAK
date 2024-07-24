/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Solve weights

///////////////////////////

//set bitmax to size of longest section
if(mode==102){bitmax=sends[section-1]-sstarts[section-1];}
else
{
bitmax=0;
for(bit=0;bit<num_sections;bit++)
{
if(sends[bit]-sstarts[bit]>bitmax){bitmax=sends[bit]-sstarts[bit];}
}
}

//can now allocate variables
data_warn(bitmax,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitmax);

anal_warn(bitmax,2*bitmax);
cors=malloc(sizeof(double*)*bitmax*bitmax);
cors2=malloc(sizeof(double)*bitmax);
tally1=malloc(sizeof(double)*bitmax);
tally2=malloc(sizeof(double)*bitmax);
tally3=malloc(sizeof(double)*bitmax);
infos=malloc(sizeof(double)*data_length);

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

if(strcmp(infosfile,"blank")!=0)	//read infos from file
{
printf("Reading infos from %s\n", infosfile);
read_infosfile(infosfile, infos, data_length, preds);
printf("First few infos are: %s %.4f", preds[0], infos[0]);
for(j=1;j<data_length;j++){if(j<5){printf(" | %s %.4f", preds[j], infos[j]);}}
printf("\n\n");

//too fiddly to exclude predictors with zero info (as would have to change section boundaries)
count=0;for(j=0;j<data_length;j++){count+=(info>0);}
if(count==0){printf("Error, there are no predictors with non-zero info\n\n");exit(1);}
if(count<data_length){printf("Warning, %d predictors have info zero\n", data_length-count);}
printf("\n");
}
else	//set infos to one
{
for(j=0;j<data_length;j++){infos[j]=1;}
}

//deal with progress file
if(mode==102){sprintf(filename,"%sprogress.%d", folder, section);}
else{sprintf(filename,"%sprogress.all", folder);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}
fclose(output);

//ready for bit loop
for(bit=0;bit<num_sections;bit++)
{
bitstart=sstarts[bit];
bitend=sends[bit];
bitlength=bitend-bitstart;

if(bit>=section_start-1)	//calculate weights for this section
{
if(mode==104)
{
printf("Calculating weights for Section %d of %d\n", bit+1, num_sections);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating weights for Section %d of %d\n", bit+1, num_sections);
fclose(output);
}

shuffle=0;
for(j=0;j<end-bitstart;j++)	//using values already in data, so shuffle back
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(bitstart-start+j)*num_samples_use+i];}
shuffle++;
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, bitstart+shuffle, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);

calc_correlations(cors, tally1, tally2, data, missingvalue, num_samples_use, bitlength, num_subs, subindex, chr+bitstart, cmbp+bitstart, window_kb, window_length, mincor, lddecay, decay, filename);

if(simplex==0)
{
flag=solve_quad(tally3, cors, infos, tally2, bitlength, filename, fudge);
}
else
{
//calculate the weights (using fudge if necessary)
#if MET==0
flag=solve_simplex_new(tally3, cors, infos+bitstart, tally2, bitlength, maxiter, maxtime, filename, fudge);
#else
flag=solve_simplex(tally3, cors, infos+bitstart, tally2, bitlength, maxiter, maxtime, filename, fudge);
#endif
}

//see how well the weights do
alpha=1.0;beta=0.0;
dgemv_("N", &bitlength, &bitlength, &alpha, cors, &bitlength, tally3, &one, &beta, cors2, &one);

//print the results
sprintf(filename2, "%sweights.%d", folder, bit+1);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
if(flag==0){fprintf(output2, "Weight ");}
if(flag==1){fprintf(output2, "Quick_Weight ");}
if(flag==2){fprintf(output2, "Approx_Weight ");}
if(lddecay==0){fprintf(output2,"Neighbours Tagging Info Check\n");}
else{fprintf(output2,"Eff_Neighbours Eff_Tagging Info Check\n");}

for(j=0;j<bitlength;j++)
{fprintf(output2, "%.6f %.2f %.3f %.3f %.3f\n", tally3[j], tally1[j], tally2[j], infos[bitstart+j], cors2[j]);}
fclose(output2);
printf("Weightings for section %d saved in %s\n\n", bit+1, filename2);

start=bitstart;if(bitend>end){end=bitend;}

if(mode==102){break;}	//only doing this section
}}	//end of doing this section and bit loop

free(data);
free(infos);free(cors);free(cors2);free(tally1);free(tally2);free(tally3);
if(binary==0){gzclose(datainputgz);}

///////////////////////////

