/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//tell the user what's about to happen and offer advice

///////////////////////////

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n\n");

///////////////////////////

if(mode==101)
{
if(section_kb!=-9999)
{
if(section_cm!=-9999)
{printf("Cutting predictors into sections of (approximate) length %.4fcM; ", section_cm);}
else
{printf("Cutting predictors into sections of (approximate) length %.2fkb; ", section_kb);}
}
else
{printf("Cutting predictors into sections of (approximate) length %d predictors; ", section_length);}
if(buffer_kb!=-9999)
{
if(buffer_cm!=-9999)
{printf("either side, will be a buffer of %.4fcM\n", buffer_cm);}
else
{printf("either side, will be a buffer of %.2fkb\n", buffer_kb);}
}
else{printf("either side, will be a buffer of %d predictors\n", buffer_length);}
printf("You can change these settings using \"--section-kb\", \"--section-length\" or \"--section-cm\" and \"--buffer-kb\", \"buffer-length\" or \"--buffer-cm\"\n\n");

if(nothin==0)
{
if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("Will first thin so that there are no pairs of predictors within %.4fcM with correlation squared greater than %.4f", window_cm, wprune);}
else
{printf("Will first thin so that there are no pairs of predictors within %.2fkb with correlation squared greater than %.4f", window_kb, wprune);}
}
else
{
if(window_length!=-1){printf("Will first thin so that there are no pairs of predictors within %d with correlation squared greater than %.4f", window_length, wprune);}
else{printf("Will first thin so that there are no pairs of predictors on the same chromosome with correlation squared greater than %.4f", wprune);}
}
printf("; you can change these settings using \"--window-prune\" and \"--window-kb\", \"--window-length\" or \"--window-cm\"\n");
printf("Note that this thinning may take a few hours, but tends to make calculating weights both faster and more stable and is especially advised when analsying dense/imputed data; if you would prefer to skip, add \"--no-thin YES\", or if you have have previously thinned, add \"--no-thin DONE\" (LDAK will then expect to find the file %sthin.in)\n\n", folder);
}
if(nothin==1)
{printf("You have chosen not to first thin predictors; if this is because you have already thinned predictors, you should instead use \"--no-thin DONE\" (LDAK will then expect to find the file %sthin.in created by a previous run of \"--cut-weights\"), or use \"--extract\" to provide a list of thinned predictors\n\n", folder);}
if(nothin==2)
{printf("Will read a list of thinned predictors from %sthin.in\n\n", folder);}
}

if(mode==102)
{
printf("Calculating weights for Section %d (out of %d)\n\n", section, num_sections);
if(num_sections>1)
{printf("Note that it is normally easier to use \"--calc-weights-all\", which calculates weights for each section in turn\n\n");}

if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("When calculating how well each predictor is tagged, will consider correlations with predictors within %.2fcM", window_cm);}
else
{printf("When calculating how well each predictor is tagged, will consider correlations with predictors within %.2fkb", window_kb);}
}
else
{
if(window_length!=-1){printf("When calculating how well each predictor is tagged, will consider correlations with predictors within %d", window_length);}
else{printf("When calculating how well each predictor is tagged, will consider correlations with predictors on the same chromosome");}
}
if(mincor>0){printf(" (correlation squared values less than %.4f are considered noise and set to zero)", mincor);}
printf("; you can change these settings using \"--window-kb\", \"--window-length\" or \"--window-cm\"\n\n");

if(num_subs>1)
{
printf("When calculating correlations, will use the highest observed across %d subsets of samples, ", num_subs);
if(num_subs==2){printf("specified by %s1 and %s2\n\n", subpref, subpref);}
if(num_subs==3){printf("specified by %s1, %s2 and %s3\n\n", subpref, subpref, subpref);}
if(num_subs>3){printf("specified by %s1, ..., %s%d\n\n", subpref, subpref, num_subs);}
}

if(lddecay==1)
{printf("Will be modelling the decay of LD with distance, assuming a half-life of %.2fkb (this is not necessary if samples are \"unrelated\" and population homogeneous)\n", halflife);}

if(fudge==0)
{
if(simplex==0){printf("By default, LDAK now uses a quadratic solver; to revert to a linear solver use \"--simplex YES\"\n");}
printf("The solver will halt after %d iterations (option \"--max-iter\")", maxiter);
if(MET==0){printf(" or after %.2f minutes (option \"--max-time\")", maxtime);}
printf(", after which point it will compute approximate weights (equivalent to using \"--quick-weights YES\")\n\n");
}
else
{printf("Will compute approximate weights; ideally you should use this only if the default solver fails to complete\n\n");}
}

if(mode==103)
{
if(num_sections==1){printf("Joining weights across one section\n\n");}
else{printf("Joining weights across %d sections\n\n", num_sections);}
}

if(mode==104)
{
printf("Calculating weights for each of the %d sections in turn", num_sections);
if(section_start>0){printf(", starting from Section %d", section_start);}
printf("\n\n");

if(num_sections>1)
{printf("If your job is killed before all sections complete, you can continue from where it finished by adding \"--start-section\" (if progress is very slow, consider instead using \"--calc-weights\" to analyse sections individually)\n\n");}

if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("When calculating how well each predictor is tagged, will consider correlations with predictors within %.2fcM", window_cm);}
else
{printf("When calculating how well each predictor is tagged, will consider correlations with predictors within %.2fkb", window_kb);}
}
else
{
if(window_length!=-1){printf("When calculating how well each predictor is tagged, will consider correlations with predictors within %d", window_length);}
else{printf("When calculating how well each predictor is tagged, will consider correlations with predictors on the same chromosome");}
}
if(mincor>0){printf(" (correlation squared values less than %.4f are considered noise and set to zero)", mincor);}
printf("; you can change these settings using \"--window-kb\", \"--window-length\" or \"--window-cm\"\n\n");

if(num_subs>1)
{
printf("When calculating correlations, will use the highest observed across %d subsets of samples, ", num_subs);
if(num_subs==2){printf("specified by %s1 and %s2\n\n", subpref, subpref);}
if(num_subs==3){printf("specified by %s1, %s2 and %s3\n\n", subpref, subpref, subpref);}
if(num_subs>3){printf("specified by %s1, ..., %s%d\n\n", subpref, subpref, num_subs);}
}

if(lddecay==1)
{printf("Will be modelling the decay of LD with distance, assuming a half-life of %.2fkb (this is not necessary if samples are \"unrelated\" and population homogeneous)\n", halflife);}

if(fudge==0)
{
if(simplex==0){printf("By default, LDAK now uses a quadratic solver; to revert to a linear solver use \"--simplex YES\"\n");}
printf("The solver will halt after %d iterations (option \"--max-iter\")", maxiter);
if(MET==0){printf(" or after %.2f minutes (option \"--max-time\")", maxtime);}
printf(", after which point it will compute approximate weights (equivalent to using \"--quick-weights YES\")\n\n");
}
else
{printf("Will compute approximate weights; ideally use this only if the default solver fails to complete\n\n");}
}

if(mode==105)
{
printf("Dividing weights across %d partitions\n\n", num_parts);
}

////////

if(mode==106)
{
printf("Thinning predictors ");
if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("so that no pair within %.4fcM remains with correlation squared greater than %.4f\n\n", window_cm, wprune);}
else
{printf("so that no pair within %.2fkb remains with correlation squared greater than %.4f\n\n", window_kb, wprune);}
}
else
{
if(window_length!=-1){printf("so that no pair less than %d apart remains with correlation squared greater than %.4f\n\n", window_length, wprune);}
else{printf("so that no pair on the same chromosome remains with correlation squared greater than %.4f\n\n", wprune);}
}

if(strcmp(pvafile,"blank")!=0)
{printf("Predictors with higher p-values will be excluded first\n\n");}
else
{printf("If you add \"--pvalues\", predictors with higher p-values will be excluded first\n\n");}
}

if(mode==107)
{
printf("Will identify predictors with p-values less than %.2e, then prune ", cutoff);
if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("so that no pair within %.4fcM remains with correlation squared greater than %.4f ", window_cm, wprune);}
else
{printf("so that no pair within %.2fkb remains with correlation squared greater than %.4f ", window_kb, wprune);}
}
else
{
if(window_length!=-1){printf("so that no pair less than %d apart remains with correlation squared greater than %.4f\n\n", window_length, wprune);}
else{printf("so that no pair on the same chromosome remains with correlation squared greater than %.4f\n\n", wprune);}
}
printf("(predictors with higher p-values will be excluded first)\n\n");
}

if(mode==108)
{
printf("Finding tags for the predictors in %s; ", scorefile);
if(window_cm!=-9999){printf("tags must be within %.4fcM", window_cm);}
else{printf("tags must be within %.2fkb", window_kb);}
if(mincor>0){printf(" and have correlation squared at least %.4f\n\n", mincor);}
else{printf(", but can be of any strength (use \"--min-cor\" to require a minimum correlation squared)\n\n");}
}

if(mode==109)
{
printf("Removing tags of predictors specified in %s; ", targetfile);
if(window_cm!=-9999){printf("tags must be within %.4fcM", window_cm);}
else{printf("tags must be within %.2fkb", window_kb);}
printf(" and have correlation squared at least %.4f\n\n", mincor);
}

if(mode==110)
{
if(nonsnp==0){printf("Thinning predictors with MAF >= %.4f, ", cutoff);}
else{printf("Thinning predictors with variance >= %.4f, ", cutoff);}
if(window_kb!=-9999)
{
if(window_cm!=-9999)
{printf("so that no pair within %.4fcM remains with correlation squared greater than %.4f\n\n", window_cm, wprune);}
else
{printf("so that no pair within %.2fkb remains with correlation squared greater than %.4f\n\n", window_kb, wprune);}
}
else
{
if(window_length!=-1){printf("so that no pair less than %d apart remains with correlation squared greater than %.4f\n\n", window_length, wprune);}
else{printf("so that no pair on the same chromosome remains with correlation squared greater than %.4f\n\n", wprune);}
}
}

///////////////////////////

if(mode==111)
{
if(part_length!=-9999)
{printf("Cutting predictors into partitions with (approximate) length %d\n\n", part_length);}
if(bychr==1)
{printf("Cutting predictors into partitions according to chromosome\n\n");}
if(strcmp(partpref,"blank")!=0)
{
if(num_parts==1){printf("Cutting predictors into one partition according to %s1\n\n", partpref);}
if(num_parts==2){printf("Cutting predictors into two partitions according to %s1 and %s2\n\n", partpref, partpref);}
if(num_parts>2){printf("Cutting predictors into %d partitions according to %s1 ... %s%d\n\n", num_parts, partpref, partpref, num_parts);}
}

if(checkpart==0)
{printf("Will check that at least one predictor in each predictor list is in the data (if sure this is the case, you can skip this using \"--check-partitions NO\")\n\n");}
}

if(mode==112)
{
printf("Calculating kinships for Partition %d (out of %d)\n\n", partition, num_parts);

print_scaling(power,hwestand);

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
if(single==0){printf("If memory is an issue, add \"--single YES\" (then calculations will performed using single precistion)\n\n");}
}

if(mode==113)
{
if(num_parts==1){printf("Joining kinships across one partition\n\n");}
else{printf("Joining kinships across %d partitions\n\n", num_parts);}

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
}

if(mode==114)
{
printf("Calculating kinships directly\n\n");

print_scaling(power,hwestand);

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
if(single==0){printf("If memory is an issue, add \"--single YES\" (then calculations will performed using single precistion)\n\n");}
}

////////

if(mode==115)
{
if(minrel!=-9999){printf("Identifying samples with kinship at least %.4f with one other\n\n", minrel);}
else
{
if(maxrel==-9999)
{printf("Filtering samples so that no pair remains with kinship greater than (the negative of) the smallest kinship observed\n\n");}
else
{printf("Filtering samples so that no pair remains with kinship greater than %.4f\n\n", maxrel);}

if(strcmp(respfile,"blank")==0)
{printf("If you add \"--pheno\", samples with non-missing phenotypes will take priority\n\n");}
printf("If you prefer to identify samples that are related (instead of unrelated), use \"--min-rel\"\n\n");
}
}

if(mode==116)
{
if(num_kins==1)
{printf("Resaving a kinship matrix (note that is normal to provide multiple kinship matrices)\n\n");}
else
{printf("Adding the %d kinship matrices specified in %s\n\n", num_kins, kinlist);}

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
}

if(mode==117)
{
if(extract==1)
{
if(strcmp(cpredfile,"blank")!=0)
{printf("Subtracting from the kinship matrix the predictors listed in %s\n\n", cpredfile);}
else
{printf("Subtracting from the kinship matrix the predictors NOT listed in %s\n\n", bpredfile);}
}
else
{
if(num_kins==1)
{printf("Resaving a kinship matrix (note that it is normal either to provide multiple kinship matrices, or to use \"--extract\" or \"--exclude\" to subtract predictors from a kinship matrix)\n\n");}
else
{
if(num_kins==2){printf("Subtracting the second kinship matrix in %s from the first\n\n", kinlist);}
else{printf("Subtracting the last %d kinship matrices in %s from the first\n\n", num_kins-1, kinlist);}
}
}

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
}

if(mode==118)
{
printf("Converting a kinship matrix from gzip to binary format\n\n");
}

if(mode==119)
{
printf("Converting a kinship matrix from raw text to binary format\n\n");
}

if(mode==120)
{
printf("Calculating pairwise correlations between off-diagonal terms of %d kinship matrices\n\n", num_kins);
}

///////////////////////////

if(mode==121)
{
printf("Performing generalized REML ");
if(num_kins+num_regs==0){printf("(although no kinship matrices or regions provided!)\n\n");}
else
{
if(num_kins==1){printf("with one kinship matrix ");}
else{printf("with %d kinship matrices ", num_kins);}
if(num_regs==1){printf("and one region\n\n");}
else{printf("and %d regions\n\n", num_regs);}

if(num_regs>0){print_scaling(power,hwestand);}

if(constrain==0){printf("Variance components can be negative (to prevent this, use \"--constrain YES\")\n\n");}
else{printf("Variance components can not be negative (to allow this, use \"--constrain NO\")\n\n");}
}

if(mpheno==-1){printf("Will analyse each of the %d phenotypes in turn\n\n", num_resps);}

if(num_kins==0&&num_regs>0)
{
if(shrink<1){printf("Will scale off-diagonal values of the regional correlation matrix by %.4f (change this value using \"--shrink\")\n\n", shrink);}
if(shrink>1){printf("Will scale diagonal values of the regional correlation matrix by %.4f (change this value using \"--shrink\")\n\n", shrink);}
if(strip!=0){printf("Will remove a proportion %.2f of the variation of the regional correlation matrix (change this value using \"--strip\")\n\n", strip);}
}

if(strcmp(sumsfile,"blank")==0)
{
if(strcmp(covarfile,"blank")==0&&strcmp(factorfile,"blank")==0){printf("Consider using \"--covar\" and/or \"--factors\" to provide quantitative or categorical covariates\n\n");}
if(strcmp(topfile,"blank")==0){printf("You can use \"--top-preds\" to include (strongly-associated) predictors; these will be treated the same as fixed effects, except that the variance they explain will count towards total heritability\n\n");}
}

if(num_kins>1)
{
if(memsave==0){printf("If memory is an issue, use \"--memory-save YES\" and kinship matrices will be read on-the-fly\n\n");}
else{printf("Kinship matrices will be read on-the-fly\n\n");}
}

if(strcmp(envfile,"blank")!=0)
{
if(discenv==1){printf("Will also compute heritabilities for subgroups of samples, specified by the environmental variables\n\n");}
else{printf("If the environmental variables specify subgroups of samples, use \"--subgroups YES\"\n\n");}
}

if(hestart==1){printf("Will use Haseman-Elston Regression to set starting heritabilities (to instead set agnostically, use \"--he-starts NO\")\n\n");}
else{printf("Will set starting heritabilities agnostically (to instead use Haseman-Elston Regression use \"--he-starts YES\")\n\n");}
}

if(mode==122)
{
printf("Calculating BLUP effects based on %d kinship matrices, %d regions and %d top predictors\n\n", num_kins, num_scores, num_tops);
}

if(mode==123)
{
printf("Performing generalized Haseman Elston Regression ");
if(num_kins+num_regs==0){printf("(although no kinship matrices or regions provided!)\n\n");}
else
{
if(num_kins==1){printf("with one kinship matrix ");}
else{printf("with %d kinship matrices ", num_kins);}
if(num_regs==1){printf("and one region\n\n");}
else{printf("and %d regions\n\n", num_regs);}
}

if(num_regs>0){print_scaling(power,hwestand);}

if(mpheno==-1){printf("Will analyse each of the %d phenotypes in turn\n\n", num_resps);}

if(num_blocks==-1){printf("Will estimate SEs using leave-one-sample-out jackknifing (change this using \"--num-blocks\")\n\n");}
else{printf("Will estimate SEs using jackknifing with %d blocks (change this using \"--num-blocks\")\n\n", num_blocks);}

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" to provide quantitative covariates\n\n");}
if(strcmp(topfile,"blank")==0){printf("You can use \"--top-preds\" to include (strongly-associated) predictors; these will be treated the same as fixed effects, except that the variance they explain will count towards total heritability\n\n");}

if(num_kins>0)
{
if(memsave==0){printf("If memory is an issue, use \"--memory-save YES\" and kinship matrices will be read on-the-fly\n\n");}
else{printf("Kinship matrices will be read on-the-fly\n\n");}
}

if(strcmp(envfile,"blank")!=0)
{
if(discenv==1){printf("Will also compute heritabilities for subgroups of samples, specified by the environmental variables\n\n");}
else{printf("If the environmental variables specify subgroups of samples, use \"--subgroups YES\"\n\n");}
}
}

if(mode==124)
{
printf("Performing generalized PCGC regression ");
if(num_kins+num_regs==0){printf("(although no kinship matrices or regions provided!)\n\n");}
else
{
if(num_kins==1){printf("with one kinship matrix ");}
else{printf("with %d kinship matrices ", num_kins);}
if(num_regs==1){printf("and one region\n\n");}
else{printf("and %d regions\n\n", num_regs);}
}

if(num_regs>0){print_scaling(power,hwestand);}

if(mpheno==-1){printf("Will analyse each of the %d phenotypes in turn\n\n", num_resps);}

if(num_blocks==-1){printf("Will estimate SEs using leave-one-sample-out jackknifing (change this using \"--num-blocks\")\n\n");}
else{printf("Will estimate SEs using jackknifing with %d blocks (change this using \"--num-blocks\")\n\n", num_blocks);}

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" to provide quantitative covariates\n\n");}
if(strcmp(topfile,"blank")==0){printf("You can use \"--top-preds\" to include (strongly-associated) predictors; these will be treated the same as fixed effects, except that the variance they explain will count towards total heritability\n\n");}

if(num_kins>0)
{
if(memsave==0){printf("If memory is an issue, use \"--memory-save YES\" and kinship matrices will be read on-the-fly\n\n");}
else{printf("Kinship matrices will be read on-the-fly\n\n");}
}

if(strcmp(envfile,"blank")!=0)
{
if(discenv==1){printf("Will also compute heritabilities for subgroups of samples, specified by the environmental variables\n\n");}
else{printf("If the environmental variables specify subgroups of samples, use \"--subgroups YES\"\n\n");}
}
}

if(mode==125)
{
if(num_kins==1){printf("Calculating BLUP predictions based on one kinship matrix\n\n");}
else{printf("Calculating BLUP predictions based on %d kinship matrices\n\n", num_kins);}
}

if(mode==126)
{
printf("Performing fast generalized REML ");
if(num_kins==1){printf("with one kinship matrix ");}
else{printf("with %d kinship matrices ", num_kins);}
if(num_regs==1){printf("and one region\n\n");}
else{printf("and %d regions\n\n", num_regs);}

if(num_regs>0){print_scaling(power,hwestand);}

if(constrain==0){printf("Variance components can be negative (to prevent this, use \"--constrain YES\")\n\n");}
else{printf("Variance components can not be negative (to allow this, use \"--constrain NO\")\n\n");}

if(mpheno==-1){printf("Will analyse each of the %d phenotypes in turn\n\n", num_resps);}

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" to provide quantitative covariates\n\n");}
if(strcmp(topfile,"blank")==0){printf("You can use \"--top-preds\" to include (strongly-associated) predictors; these will be treated the same as fixed effects, except that the variance they explain will count towards total heritability\n\n");}

if(num_kins>1)
{
if(memsave==0){printf("If memory is an issue, use \"--memory-save YES\" and kinship matrices will be read on-the-fly\n\n");}
else{printf("Kinship matrices will be read on-the-fly\n\n");}
}

if(hestart==1){printf("Will use Haseman-Elston Regression to set starting heritabilities (to instead set agnostically, use \"--he-starts NO\")\n\n");}
else{printf("Will set starting heritabilities agnostically (to instead use Haseman-Elston Regression use \"--he-starts YES\")\n\n");}
}

////////

if(mode==127)
{
printf("Performing Randomized Haseman Elston Regression ");
if(num_parts>0)
{
if(parttype==0)
{
if(num_parts==1){printf(" using a base category and one annotation\n");}
else{printf(" using a base category and %d annotations\n", num_parts);}
}
else
{
if(num_parts==1){printf(" using one partition\n");}
else{printf(" using %d partitions\n", num_parts);}
}
}
else{printf("\n");}

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" to provide quantitative covariates\n\n");}
if(strcmp(topfile,"blank")==0){printf("You can use \"--top-preds\" to include (strongly-associated) predictors; these will be treated the same as fixed effects, except that the variance they explain will count towards total heritability\n\n");}

printf("Will estimate SEs using jackknifing with %d blocks (change this using \"--num-blocks\")\n\n", num_blocks);

if(num_parts>0&&strcmp(labfile,"blank")==0)
{
if(parttype==0){printf("Consider using \"--labels\" to provide a file containing the annotation names\n\n");}
else{printf("Consider using \"--labels\" to provide a file containing the partition names\n\n");}
}
}

if(mode==128)
{
printf("Performing Randomized PCGC Regression ");
if(num_parts>0)
{
if(parttype==0)
{
if(num_parts==1){printf(" using a base category and one annotation\n");}
else{printf(" using a base category and %d annotations\n", num_parts);}
}
else
{
if(num_parts==1){printf(" using one partition\n");}
else{printf(" using %d partitions\n", num_parts);}
}
}
else{printf("\n");}

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" to provide quantitative covariates\n\n");}
if(strcmp(topfile,"blank")==0){printf("You can use \"--top-preds\" to include (strongly-associated) predictors; these will be treated the same as fixed effects, except that the variance they explain will count towards total heritability\n\n");}

printf("Will estimate SEs using jackknifing with %d blocks (change this using \"--num-blocks\")\n\n", num_blocks);

if(num_parts>0&&strcmp(labfile,"blank")==0)
{
if(parttype==0){printf("Consider using \"--labels\" to provide a file containing the annotation names\n\n");}
else{printf("Consider using \"--labels\" to provide a file containing the partition names\n\n");}
}
}

////////

if(mode==129)
{
printf("Estimating heritability using relatives for a quantitative phenotype\n\n");

if(mpheno==-1){printf("Will analyse each of the %d phenotypes in turn\n\n", num_resps);}

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" to provide quantitative covariates\n\n");}
}

if(mode==130)
{
printf("Estimating heritability using relatives for a binary phenotype\n\n");

if(mpheno==-1){printf("Will analyse each of the %d phenotypes in turn\n\n", num_resps);}

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" to provide quantitative covariates\n\n");}
}

///////////////////////////

if(mode==131)
{
if(num_kins==0&&strcmp(envfile,"blank")==0&&strcmp(prsfile,"blank")==0&&families==0&&trios==0&&duos==0)	//basic case
{
if(num_resps_use==1){printf("Performing linear regression for one phenotype\n\n");}
else{printf("Performing linear regression for %d phenotypes\n\n", num_resps_use);}

if(strcmp(covarfile,"blank")!=0||strcmp(topfile,"blank")!=0||strcmp(factorfile,"blank")!=0)
{
if(adjpreds==0){printf("Will not adjust predictors for covariates (use \"--adjust-predictors YES\" to adjust all predictors, or use \"--adjust-predictors PARTIAL\" to adjust only the most significant predictors)\n\n");}
if(adjpreds==1){printf("Will only adjust the most significant predictors for covariates (use \"--adjust-predictors YES\" to adjust all predictors, or use \"--adjust-predictors NO\" to not adjust predictors)\n\n");}
if(adjpreds==2){printf("Will adjust all predictors for covariates (use \"--adjust-predictors NO\" to not adjust predictors, or use \"--adjust-predictors PARTIAL\" to adjust only the most significant predictors)\n\n");}
}

if(strcmp(sampwfile,"blank")==0){printf("To perform weighted linear regression, use \"--sample-weights\"\n\n");}
else
{
if(sandwich==0){printf("Will use the standard estimator of effect size variance; use \"--sandwich YES\" to switch to a sandwich estimator\n\n");}
else{printf("Will use the sandwich estimator of effect size variance; to switch to the standard estimator use \"--sandwich NO\"\n\n");}
}

if(spatest==0){printf("Will compute standard test statistics; use \"--spa-test YES\" to switch to a saddlepoint approximation, \n\n");}
else{printf("Will use a saddlepoint approximation for predictors with (non SPA) p-value below %.2e (use \"--spa-threshold\" to change the threshold, or use \"--spa-test NO\" to turn off the SPA)\n\n", spathresh);}
}

if(num_kins==1)
{
printf("Performing mixed-model linear regression\n\n");

if(exact==0){printf("Variance components will be fixed under the null model; to instead re-estimate for each predictor add \"--exact YES\"\n\n");}
else{printf("Variance components will be re-estimated for each predictor\n\n");}
}

if(strcmp(prsfile,"blank")!=0)
{
if(num_resps_use==1){printf("Performing linear regression including LOCO PRS as an offset for one phenotype\n\n");}
else{printf("Performing linear regression including LOCO PRS as an offset for %d phenotypes\n\n", num_resps_use);}

if(strcmp(covarfile,"blank")!=0||strcmp(topfile,"blank")!=0||strcmp(factorfile,"blank")!=0)
{
if(adjpreds==0){printf("Will not adjust predictors for covariates (use \"--adjust-predictors YES\" to adjust all predictors, or use \"--adjust-predictors PARTIAL\" to adjust only the most significant predictors)\n\n");}
if(adjpreds==1){printf("Will only adjust the most significant predictors for covariates (use \"--adjust-predictors YES\" to adjust all predictors, or use \"--adjust-predictors NO\" to not adjust predictors)\n\n");}
if(adjpreds==2){printf("Will adjust all predictors for covariates (use \"--adjust-predictors NO\" to not adjust predictors, or use \"--adjust-predictors PARTIAL\" to adjust only the most significant predictors)\n\n");}
}

if(spatest==0){printf("Will compute standard test statistics; use \"--spa-test YES\" to switch to a saddlepoint approximation, \n\n");}
else{printf("Will use a saddlepoint approximation for predictors with (non SPA) p-value below %.2e (use \"--spa-threshold\" to change the threshold, or use \"--spa-test NO\" to turn off the SPA)\n\n", spathresh);}
}

if(strcmp(envfile,"blank")!=0){printf("Performing linear regression with an environmental interaction\n\n");}

if(families==1){printf("Performing family-based linear regression\n\n");}

if(trios==1){printf("Performing trio-based linear regression\n\n");}

if(duos==1){printf("Performing parent-offspring linear regression\n\n");}
if(duos==2){printf("Performing father-offspring linear regression\n\n");}
if(duos==3){printf("Performing mother-offspring linear regression\n\n");}

if(strcmp(prsfile,"blank")==0)
{
flag=0;
if(strcmp(covarfile,"blank")==0&&strcmp(factorfile,"blank")==0){printf("Consider using \"--covar\" and/or \"--factors\" to provide quantitative or categorical covariates\n");flag=1;}
if(strcmp(topfile,"blank")==0&&strcmp(envfile,"blank")==0&&families==0&&trios==0&&duos==0)
{printf("You can use \"--top-preds\" to include (strongly-associated) predictors as extra covariates\n");flag=1;}
if(flag==1){printf("\n");}
}

print_qc(minmaf, maxmaf, minvar, minobs, mininfo, genprobs, 1);
}

if(mode==132)
{
if(strcmp(prsfile,"blank")==0)
{
if(num_resps_use==1){printf("Performing logistic regression for one phenotype\n\n");}
else{printf("Performing logistic regression for %d phenotypes\n\n", num_resps_use);}

if(scoretest==0){printf("Will use a Wald Test; use \"--score-test YES\" to switch to a Score Test (which is a bit less accurate, but much faster)\n\n");}
else{printf("Will use a Score Test; use \"--score-test NO\" to switch to a Wald Test (which is a bit more accurate, but much slower)\n\n");}

flag=0;
if(strcmp(covarfile,"blank")==0&&strcmp(factorfile,"blank")==0){printf("Consider using \"--covar\" and/or \"--factors\" to provide quantitative or categorical covariates\n");flag=1;}
if(strcmp(topfile,"blank")==0&&strcmp(envfile,"blank")==0&&families==0&&trios==0&&duos==0)
{printf("You can use \"--top-preds\" to include (strongly-associated) predictors as extra covariates\n");flag=1;}
if(flag==1){printf("\n");}
}
else
{
if(num_resps_use==1){printf("Performing quasi-logistic regression including LOCO PRS as an offset for one phenotype\n\n");}
else{printf("Performing quasi-logistic regression including LOCO PRS as an offset for %d phenotypes\n\n", num_resps_use);}
}

if(strcmp(covarfile,"blank")!=0||strcmp(topfile,"blank")!=0||strcmp(factorfile,"blank")!=0)
{
if(adjpreds==0){printf("Will not adjust predictors for covariates (use \"--adjust-predictors YES\" to adjust all predictors, or use \"--adjust-predictors PARTIAL\" to adjust only the most significant predictors)\n\n");}
if(adjpreds==1){printf("Will only adjust the most significant predictors for covariates (use \"--adjust-predictors YES\" to adjust all predictors, or use \"--adjust-predictors NO\" to not adjust predictors)\n\n");}
if(adjpreds==2){printf("Will adjust all predictors for covariates (use \"--adjust-predictors NO\" to not adjust predictors, or use \"--adjust-predictors PARTIAL\" to adjust only the most significant predictors)\n\n");}
}

if(spatest==0){printf("Will compute standard test statistics; use \"--spa-test YES\" to switch to a saddlepoint approximation, \n\n");}
else{printf("Will use a saddlepoint approximation for predictors with (non SPA) p-value below %.2e (use \"--spa-threshold\" to change the threshold, or use \"--spa-test NO\" to turn off the SPA)\n\n", spathresh);}

print_qc(minmaf, maxmaf, minvar, minobs, mininfo, genprobs, 1);
}

if(mode==133)
{
printf("Performing Step 1 of complex mixed-model linear regression ");
if(num_kins==1){printf("with one kinship matrix\n\n");}
else{printf("with %d kinship matrices\n\n", num_kins);}

if(constrain==0){printf("Variance components can be negative (to prevent this, use \"--constrain YES\")\n\n");}
else{printf("Variance components can not be negative (to allow this, use \"--constrain NO\")\n\n");}

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" or \"--factors\" to provide quantitative covariates\n\n");}
if(strcmp(topfile,"blank")==0){printf("You can use \"--top-preds\" to include (strongly-associated) predictors as extra covariates\n\n");}

if(num_kins>0)
{
if(memsave==0){printf("If memory is an issue, use \"--memory-save YES\" and kinship matrices will be read on-the-fly\n\n");}
else{printf("Kinship matrices will be read on-the-fly\n\n");}
}

if(hestart==1){printf("Will use Haseman-Elston Regression to set starting heritabilities (to instead set agnostically, use \"--he-starts NO\")\n\n");}
else{printf("Will set starting heritabilities agnostically (to instead use Haseman-Elston Regression use \"--he-starts Y\")\n\n");}
}

////////

if(mode==136)
{
if(strcmp(genefile,"blank")!=0)
{
printf("Cutting predictors into genes based on the annotations provided in %s\n\n", genefile);

if(up_buffer>0||down_buffer>0||gene_buffer>0)
{
if(gene_buffer>0){printf("Will include predictors within %d basepairs of each gene\n\n", gene_buffer);}
else{printf("Will include predictors within %d (up) and %d (down) basepairs of each gene\n\n", up_buffer, down_buffer);}
}
else{printf("Use \"--gene-buffer\" (or \"--up-buffer\" and \"--down-buffer\") to include predictors neighbouring genes\n\n");}

if(minweight!=1e-10)
{
if(strcmp(weightsfile,"blank")==0){printf("Will only consider genes containing at least %.2f predictors\n\n", minweight);}
else{printf("Will only consider genes with weighting at least %.2f\n\n", minweight);}
}
else{printf("Will include all genes (use \"--min-weight\" to exclude small genes)\n\n");}

if(overlap==1){printf("Will allow overlap between adjacent genes (to prevent this, use \"--overlap NO\")\n\n");}
else{printf("Will ensure predictors are assigned to at most one gene\n\n");}

if(part_length!=-9999){printf("Genes will be divided into partitions of (approximate) length %d predictors\n\n", part_length);}
if(bychr==1){printf("Genes will be divided into partitions according to chromosome\n\n");}
if(part_length==-9999&&bychr==0){printf("If you plan to use \"--calc-genes-kins\" or \"--calc-genes-reml\", consider adding \"--partition-length\" or \"--by-chr YES\" to divide genes into partitions (which will allow the next step to be parallelized)\n\n");}

if(strcmp(pvafile,"blank")==0){printf("If you add \"--pvalues\", LDAK will work out the lowest p-value for each gene\n\n");}
}
else	//so using chunks or chunks bp
{
if(chunks!=-9999)
{
if(strcmp(weightsfile,"blank")!=0)
{
if(overlap==1){printf("Cutting predictors into overlapping chunks with weighting %.2f (to avoid overlap, add \"--overlap NO\")\n\n", chunks);}
else{printf("Cutting predictors into non-overlapping chunks with weighting %.2f\n\n", chunks);}
}
else
{
if(overlap==1){printf("Cutting predictors into overlapping chunks of length %.2f (to avoid overlap, add \"--overlap NO\")\n\n", chunks);}
else{printf("Cutting predictors into non-overlapping chunks of length %.2f\n\n", chunks);}
}
}
else
{
if(overlap==1){printf("Cutting predictors into overlapping chunks of length %d basepairs (to avoid overlap, add \"--overlap NO\")\n\n", chunksbp);}
else{printf("Cutting predictors into non-overlapping chunks of length %d basepairs\n\n", chunksbp);}
}

if(minweight!=1e-10)
{
if(strcmp(weightsfile,"blank")==0){printf("Will only consider chunks containing at least %.2f predictors\n\n", minweight);}
else{printf("Will only consider chunks with weighting at least %.2f\n\n", minweight);}
}

if(part_length!=-9999){printf("Chunks will be divided into partitions of (approximate) length %d predictors\n\n", part_length);}
if(bychr==1){printf("Chunks will be divided into partitions according to chromosome\n\n");}
if(part_length==-9999&&bychr==0){printf("If you plan to use \"--calc-genes-kins\" or \"--calc-genes-reml\", consider adding \"--partition-length\" or \"--by-chr YES\" to divide chunks into partitions (which will allow the next step to be parallelized)\n\n");}


if(strcmp(pvafile,"blank")==0){printf("If you add \"--pvalues\", LDAK will work out the lowest p-value for each chunk\n\n");}
}
}

if(mode==137)
{
printf("Calculating kinships for genes/chunks in Partition %d (out of %d)\n\n", partition, num_parts);

print_scaling(power,hwestand);

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrices\n\n");}
}

if(mode==138)
{
printf("Performing gene/chunk-based REML");
if(num_kins==1){printf(" with a kinship matrix,");}
printf(" for Partition %d (out of %d)\n\n", partition, num_parts);

print_scaling(power,hwestand);

if(gene_perms>0)
{
if(gene_perms==1){printf("Will analyse the observed phenotypes and one permutation ");}
else{printf("Will analyse the observed phenotype and %d permutations", gene_perms);}
printf(" (change this number using \"--gene-permutations\")\n\n");
}

if(gprune<=1)
{printf("Will thin each gene/chunk so that there are no pairs of predictors with correlation squared greater than %.4f (change this threshold using \"--gene-prune\")\n\n", gprune);}

if(shrink<1)
{printf("Will scale off-diagonal values of each gene/chunk correlation matrix by %.4f (change this value using \"--shrink\")\n\n", shrink);}
if(shrink>1)
{printf("Will scale diagonal values of each gene/chunk correlation matrix by %.4f (change this value using \"--shrink\")\n\n", shrink);}
//if(strcmp(sumsfile,"blank")!=0&&shrink==1){printf("You can use \"--shrink\" to scale values of each gene/chunk correlation matrix\n\n");}

if(strcmp(sumsfile,"blank")==0)
{
if(strcmp(covarfile,"blank")==0&&strcmp(factorfile,"blank")==0){printf("Consider using \"--covar\" and/or \"--factors\" to provide quantitative or categorical covariates\n\n");}
if(strcmp(topfile,"blank")==0){printf("Consider using \"--top-preds\" to include (strongly-associated) predictors as extra covariates\n\n");}
}
}

if(mode==139)
{
if(num_parts==1){printf("Joining results from gene/chunk-based REML across one partition\n\n");}
else{printf("Joining results from gene/chunk-based REML across %d partitions\n\n", num_parts);}
}

////////

if(mode==140)
{
if(strcmp(genefile,"blank")!=0)
{
printf("Will first cut predictors into genes based on the annotations provided in %s\n\n", genefile);

if(up_buffer>0||down_buffer>0||gene_buffer>0)
{
if(gene_buffer>0){printf("Will include predictors within %d basepairs of each gene\n\n", gene_buffer);}
else{printf("Will include predictors within %d (up) and %d (down) basepairs of each gene\n\n", up_buffer, down_buffer);}
}
else{printf("Use \"--gene-buffer\" (or \"--up-buffer\" and \"--down-buffer\") to include predictors neighbouring genes\n\n");}

if(minweight!=1e-10)
{
if(strcmp(weightsfile,"blank")==0){printf("Will only consider genes containing at least %.2f predictors\n\n", minweight);}
else{printf("Will only consider genes with weighting at least %.2f\n\n", minweight);}
}
else{printf("Will include all genes (use \"--min-weight\" to exclude small genes)\n\n");}

if(overlap==1){printf("Will allow overlap between adjacent genes (to prevent this, use \"--overlap NO\")\n\n");}
else{printf("Will ensure predictors are assigned to at most one gene\n\n");}
}
else	//so using chunks or chunks bp
{
if(chunks!=-9999)
{
if(strcmp(weightsfile,"blank")!=0)
{
if(overlap==1){printf("Will first cut predictors into overlapping chunks with weighting %.2f (to avoid overlap, add \"--overlap NO\")\n\n", chunks);}
else{printf("Will first cut predictors into non-overlapping chunks with weighting %.2f\n\n", chunks);}
}
else
{
if(overlap==1){printf("Cutting predictors into overlapping chunks of length %.2f (to avoid overlap, add \"--overlap NO\")\n\n", chunks);}
else{printf("Cutting predictors into non-overlapping chunks of length %.2f\n\n", chunks);}
}
}
else
{
if(overlap==1){printf("Cutting predictors into overlapping chunks of length %d basepairs (to avoid overlap, add \"--overlap NO\")\n\n", chunksbp);}
else{printf("Cutting predictors into non-overlapping chunks of length %d basepairs\n\n", chunksbp);}
}

if(minweight!=1e-10)
{
if(strcmp(weightsfile,"blank")==0){printf("Will only consider chunks containing at least %.2f predictors\n\n", minweight);}
else{printf("Will only consider chunks with weighting at least %.2f\n\n", minweight);}
}
}

printf("Will then perform gene/chunk-based REML\n\n");

print_scaling(power,hwestand);

if(gene_perms>0)
{
if(gene_perms==1){printf("Will analyse the observed phenotypes and one permutation ");}
else{printf("Will analyse the observed phenotype and %d permutations", gene_perms);}
printf(" (change this number using \"--gene-permutations\")\n\n");
}

if(gprune<=1)
{printf("Will thin each gene/chunk so that there are no pairs of predictors with correlation squared greater than %.4f (change this threshold using \"--gene-prune\")\n\n", gprune);}

if(shrink<1)
{printf("Will scale off-diagonal values of each gene/chunk correlation matrix by %.4f (change this value using \"--shrink\")\n\n", shrink);}
if(shrink>1)
{printf("Will scale diagonal values of each gene/chunk correlation matrix by %.4f (change this value using \"--shrink\")\n\n", shrink);}
}

///////////////////////////

if(mode==141)
{
if(strcmp(genpreds,"blank")==0)
{
printf("Calculating a tagging file");
if(num_parts>0)
{
if(parttype==0)
{
if(num_parts==1){printf(" using a base category and one annotation");}
else{printf(" using a base category and %d annotations", num_parts);}
}
else
{
if(num_parts==1){printf(" using one partition");}
else{printf(" using %d partitions", num_parts);}
}
}
printf("\n\n");
}
else{printf("Calculating tagging files for %d pathways\n\n", num_parts);}

if(window_cm!=-9999)
{printf("Will calculate correlations between pairs of predictors within %.4fcM (change the window size using \"--window-cm\" or \"--window-kb\")", window_cm);}
else
{printf("Will calculate correlations between pairs of predictors within %.2fkb (change the window size using \"--window-cm\" or \"--window-kb\")", window_kb);}
if(mincor>0){printf("; correlation squared values less than %.4f are considered noise and set to zero (change this threshold using \"--min-cor\")", mincor);}
printf("\n\n");

if(strcmp(printfile,"blank")==0)
{printf("If you already know the subset of predictors you plan to use with \"--sum-hers\" or \"--sum-cors\", you can use \"--regression-predictors\" (this will reduce the size of the tagging file and save you having to compute the regression taggings separately)\n\n");}

if(strcmp(herfile,"blank")==0)
{printf("When using \"--sum-hers\" or \"--sum-cors\", heritability estimates will be computed across all predictors; to instead restrict to a subset of predictors use \"--heritability-predictors\"\n\n");}

if(strcmp(genpreds,"blank")==0&&savemat==0){printf("Add \"--save-matrix YES\" to also compute a heritability matrix (this is required to estimate the per-predictor heritabilities used for prediction)\n\n");}

if(num_parts>0&&strcmp(labfile,"blank")==0)
{
if(parttype==0){printf("Consider using \"--labels\" to provide a file containing the annotation names\n\n");}
else{printf("Consider using \"--labels\" to provide a file containing the partition names\n\n");}
}

if(strcmp(genpreds,"blank")!=0)
{
if(fourdp==0){printf("Will save taggings to one decimal place; to increase to four, use \"--full-accuracy YES\"\n\n");}
else{printf("Will save taggings to four decimal places; to reduce to one, use \"--full-accuracy NO\"\n\n");}
}
}

if(mode==142)
{
if(strcmp(taglist,"blank")!=0)
{
if(strcmp(matlist,"blank")==0)
{printf("Joining %d tagging files\n\nIf you added \"--save-matrix YES\" when calculating the tagging files then you should use \"--matlist\" to provide the corresponding list of heritability matrices\n\n", num_tags);}
else{printf("Joining %d tagging files, and the corresponding heritability matrices\n\n", num_tags);}
}
else
{printf("Joining %d sets of pathway tagging files\n\n", num_tags);}
}

if(mode==143)
{
if(strcmp(taglist,"blank")!=0)
{
printf("Merging %d tagging files\n\n", num_tags);

printf("Warning, this feature is mainly for model testing (e.g., to see whether combining two models significantly improves model fit or changes the overall estimate of heritability); for most other uses, including if you wish to estimate enrichments, you should construct the tagging file from scratch\n\n");
}
else
{printf("Merging %d sets of pathway tagging files\n\n", num_tags);}
}

if(mode==144)
{
if(strcmp(matfile,"blank")==0){printf("Reducing a tagging file; ");}
else{printf("Reducing a tagging file and a heritability matrix; ");}
if(parttype==0)
{
if(num_reds==1){printf("will be retaining only the base category\n\n");}
else{printf("will be retaining %d out of %d annotations, plus the base\n\n", num_reds-1, num_parts-1);}
}
else
{printf("will be retaining %d out of %d partitions\n\n", num_reds, num_parts);}

if(strcmp(matfile,"blank")==0)
{printf("Consider using \"--matrix\" to provide the corresponding heritability matrix\n\n");}
}

////////

if(mode==146)
{
if(strcmp(pathfile,"blank")==0)
{
if(num_parts==1){printf("Estimating heritability from summary statistics\n");}
else
{
if(parttype==0)
{
if(num_parts==2){printf("Estimating heritabilities from summary statistics with the base category and one annotation\n");}
else{printf("Estimating heritabilities from summary statistics with the base category and %d annotations\n", num_parts-1);}
}
else
{
if(divide==0){printf("Estimating heritabilities from summary statistics with %d partitions\n", num_parts);}
else{printf("Estimating heritability from summary statistics using %d sets of %d partitions\n", num_parts/divide, divide);}
}
}

if(strcmp(matfile,"blank")==0)
{printf("If you use \"--matrix\" to provide a heritability matrix, LDAK will additionally calculate the heritability contributed by each predictor\n\n");}
}
else{printf("Estimating heritability of %d pathways from summary statistics\n", num_parts);}

if(gcon==0&&cept==0){printf("Assuming no inflation of test statistics due to confounding");}
if(gcon==1&&cept==0){printf("Assuming multiplicative inflation of test statistics due to confounding");}
if(gcon==0&&cept==1){printf("Assuming additive inflation of test statistics due to confounding");}
if(gcon==1&&cept==1){printf("Assuming both additive and multiplicative inflation of test statistics due to confounding");}
printf(" (change this by adding \"--genomic-control YES\" or \"--intercept YES\")\n\n");

if(strcmp(altfile,"blank")==0){printf("The regression weights will be the inverse of the taggings in the tagging file; to provide different taggings, use \"--alternative-tags\"\n\n");}

if(strcmp(pathfile,"blank")==0&&strcmp(cvsfile,"blank")==0){printf("To perform cross-validation, use \"--cv-predictors\" (LDAK will then predict the test statistics of these predictors based on the test statistics of the remaining predictors)\n\n");}

if(chisol==1){printf("Will use a maximum likelihood solver; to revert to the (less-efficient) weighted least squared solver, use \"--chisq-solver NO\"\n\n");}

if(cutoff!=-9999){printf("Will be excluding predictors which individually explain at least %.4f of phenotypic variance\n\n", cutoff);}
else{printf("Warning, estimates of heritability can be biased by the presence of strong-effect loci, so either identify these yourself then use \"--exclude\", or use \"--cutoff\" to specify the maximum variance explained\n\n");}

if(num_blocks!=-9999){printf("Will estimate SEs using jackknifing with %d blocks (change this using \"--num-blocks\")\n\n", num_blocks);}
}

if(mode==147)
{
if(num_parts==1){printf("Estimating the genetic correlation from summary statistics\n");}
else
{
if(parttype==0)
{
if(num_parts==2){printf("Estimating genetic correlations from summary statistics with the base category and one annotation\n");}
else{printf("Estimating genetic correlations from summary statistics with the base category and %d annotations\n", num_parts-1);}
}
else
{printf("Estimating genetic correlations from summary statistics with %d partitions\n", num_parts);}
}

if(gcon==0&&cept==0){printf("Assuming no inflation of test statistics due to confounding");}
if(gcon==1&&cept==0){printf("Assuming multiplicative inflation of test statistics due to confounding");}
if(gcon==0&&cept==1){printf("Assuming additive inflation of test statistics due to confounding");}
if(gcon==1&&cept==1){printf("Assuming both additive and multiplicative inflation of test statistics due to confounding");}
printf(" (change this by adding \"--genomic-control YES\" or \"--intercept YES\")\n\n");

if(cutoff!=-9999){printf("Will be excluding predictors which individually explain at least %.6f of phenotypic variance\n\n", cutoff);}
else{printf("Warning, estimates of heritability can be biased by the presence of strong-effect loci, so either identify these yourself then use \"--exclude\", or use \"--cutoff\" to specify the maximum variance explained\n\n");}

if(num_blocks!=-9999){printf("Will estimate SEs using jackknifing with %d blocks (change this using \"--num-blocks\")\n\n", num_blocks);}
}

if(mode==149)
{
printf("Calculating estimates of the heritability tagged by each predictor\n\n");
}

if(mode==150)
{
printf("Calculating posterior effects of each predictor\n\n");
}

///////////////////////////

if(mode==151)
{
if(loco==0){printf("Constructing a ridge regression PRS\n\n");}
else{printf("Constructing LOCO Ridge Regression PRS\n\n");}

if(dichot==1){printf("Will analyze using a quasi-logistic model (use \"--binary NO\" to switch to a linear model)\n\n");} 

if(gctastep==-9999&&faststep==-9999)
{
if(power==-9999&&strcmp(powfile,"blank")==0){printf("Will consider five values for the predictor scaling (alpha = -1, -0.75, -0.5, -0.25 and 0); to instead specify the value, use \"--power\" (or use \"--powerfile\" to provide a range of values)\n\n");}
if(strcmp(powfile,"blank")!=0){printf("Will consider the predictor scalings in %s\n\n", powfile);}

if(fprs==0){printf("Will exclude the polygenic contribution if its estimated mean-squared error is above 0.995 (to always include, use \"--force-PRS YES\")\n\n");}
else{printf("Will always include the polygenic contribution, regardless of its estimated accuracy\n\n");}

if(fast==1){printf("When constructing PRS, will scan the data at most %d times (change this using \"--num-scans\")\n\n", nscan);}

if(her==-9999){printf("All heritability estimates must be between 0.01 and %.4f (change the upper bound using \"--max-her\")\n\n", maxher);}

if(nmcmc==-9999)
{printf("Will use either three or ten random vectors for Monte Carlo operations (decided based on the number of samples); change this number using \"--num-random-vectors\"\n\n");}
}

if(strcmp(covarfile,"blank")==0&&strcmp(factorfile,"blank")==0){printf("Consider using \"--covar\" and/or \"--factors\" to provide quantitative or categorical covariates\n\n");}

print_qc(minmaf, maxmaf, minvar, minobs, mininfo, genprobs, 1);

if(verbose==0){printf("LDAK will only save selected files (to save all files, use \"--verbose YES\")\n\n");}
}

if(mode==152)
{
if(loco==0){printf("Constructing a Bolt-LMM PRS\n\n");}
else{printf("Constructing LOCO Bolt-LMM PRS\n\n");}

if(dichot==1){printf("Will analyze using a quasi-logistic model (use \"--binary NO\" to switch to a linear model)\n\n");}

if(power==-9999&&strcmp(powfile,"blank")==0){printf("Will consider five values for the predictor scaling (alpha = -1, -0.75, -0.5, -0.25 and 0); to instead specify the value, use \"--power\" (or use \"--powerfile\" to provide a range of values)\n\n");}
if(strcmp(powfile,"blank")!=0){printf("Will consider the predictor scalings in %s\n\n", powfile);}

if(strcmp(fracfile,"blank")==0)
{printf("Will use the default prior parameter choices (saved in the file %s.parameters); to instead specify your own, use \"--parameters\"\n\n", outfile);}
else
{printf("Will use the %d sets of prior parameter choices provided in %s\n\n", countrows(fracfile), fracfile);}

if(skipcv==0)
{printf("Will select the best prior parameters via cross-validation, using %.2f randomly-picked test samples (use \"--cv-proportion\" to change this proportion, \"--cv-samples\" to explicitly specify the test samples, or \"--cv-skip\" to turn off cross-validation)\n\n", cvprop);}

if(fprs==0){printf("Will exclude the polygenic contribution if its estimated mean-squared error is above 0.995 (to always include, use \"--force-PRS YES\")\n\n");}
else{printf("Will always include the polygenic contribution, regardless of its estimated accuracy\n\n");}

if(fast==1){printf("When constructing PRS, will scan the data at most %d times (change this using \"--num-scans\")\n\n", nscan);}

if(her==-9999){printf("All heritability estimates must be between 0.01 and %.4f (change the upper bound using \"--max-her\")\n\n", maxher);}

if(nmcmc==-9999)
{printf("Will use either three or ten random vectors for Monte Carlo operations (decided based on the number of samples); change this number using \"--num-random-vectors\"\n\n");}

if(strcmp(covarfile,"blank")==0&&strcmp(factorfile,"blank")==0){printf("Consider using \"--covar\" and/or \"--factors\" to provide quantitative or categorical covariates\n\n");}

print_qc(minmaf, maxmaf, minvar, minobs, mininfo, genprobs, 1);

if(verbose==0){printf("LDAK will only save selected files (to save all files, use \"--verbose YES\")\n\n");}
}

if(mode==153)
{
if(loco==0){printf("Constructing a BayesR PRS\n\n");}
else{printf("Constructing LOCO BayesR PRS\n\n");}

if(dichot==1){printf("Will analyze using a quasi-logistic model (use \"--binary NO\" to switch to a linear model)\n\n");}

if(power==-9999&&strcmp(powfile,"blank")==0){printf("Will consider five values for the predictor scaling (alpha = -1, -0.75, -0.5, -0.25 and 0); to instead specify the value, use \"--power\" (or use \"--powerfile\" to provide a range of values)\n\n");}
if(strcmp(powfile,"blank")!=0){printf("Will consider the predictor scalings in %s\n\n", powfile);}

if(strcmp(fracfile,"blank")==0)
{printf("Will use the default prior parameter choices (saved in the file %s.parameters); to instead specify your own, use \"--parameters\"\n\n", outfile);}
else
{printf("Will use the %d sets of prior parameter choices provided in %s\n\n", countrows(fracfile), fracfile);}

if(skipcv==0)
{printf("Will select the best prior parameters via cross-validation, using %.2f randomly-picked test samples (use \"--cv-proportion\" to change this proportion, \"--cv-samples\" to explicitly specify the test samples, or \"--cv-skip\" to turn off cross-validation)\n\n", cvprop);}

if(fprs==0){printf("Will exclude the polygenic contribution if its estimated mean-squared error is above 0.995 (to always include, use \"--force-PRS YES\")\n\n");}
else{printf("Will always include the polygenic contribution, regardless of its estimated accuracy\n\n");}

if(fast==1){printf("When constructing PRS, will scan the data at most %d times (change this using \"--num-scans\")\n\n", nscan);}

if(her==-9999){printf("All heritability estimates must be between 0.01 and %.4f (change the upper bound using \"--max-her\")\n\n", maxher);}

if(nmcmc==-9999)
{printf("Will use either three or ten random vectors for Monte Carlo operations (decided based on the number of samples); change this number using \"--num-random-vectors\"\n\n");}

if(strcmp(covarfile,"blank")==0&&strcmp(factorfile,"blank")==0){printf("Consider using \"--covar\" and/or \"--factors\" to provide quantitative or categorical covariates\n\n");}

print_qc(minmaf, maxmaf, minvar, minobs, mininfo, genprobs, 1);

if(verbose==0){printf("LDAK will only save selected files (to save all files, use \"--verbose YES\")\n\n");}
}

if(mode==154)
{
if(loco==0){printf("Constructing an elastic net PRS\n\n");}
else{printf("Constructing LOCO Elastic Net PRS\n\n");}

if(dichot==1){printf("Will analyze using a quasi-logistic model (use \"--binary NO\" to switch to a linear model)\n\n");}

if(kvikstep==-9999)
{
if(power==-9999&&strcmp(powfile,"blank")==0){printf("Will consider five values for the predictor scaling (alpha = -1, -0.75, -0.5, -0.25 and 0); to instead specify the value, use \"--power\" (or use \"--powerfile\" to provide a range of values)\n\n");}
if(strcmp(powfile,"blank")!=0){printf("Will consider the predictor scalings in %s\n\n", powfile);}

if(strcmp(fracfile,"blank")==0)
{printf("Will use the default prior parameter choices (saved in the file %s.parameters); to instead specify your own, use \"--parameters\"\n\n", outfile);}
else
{printf("Will use the %d sets of prior parameter choices provided in %s\n\n", countrows(fracfile), fracfile);}

if(skipcv==0)
{printf("Will select the best prior parameters via cross-validation, using %.2f randomly-picked test samples (use \"--cv-proportion\" to change this proportion, \"--cv-samples\" to explicitly specify the test samples, or \"--cv-skip\" to turn off cross-validation)\n\n", cvprop);}

if(fprs==0){printf("Will exclude the polygenic contribution if its estimated mean-squared error is above 0.995 (to always include, use \"--force-PRS YES\")\n\n");}
else{printf("Will always include the polygenic contribution, regardless of its estimated accuracy\n\n");}

if(fast==1){printf("When constructing PRS, will scan the data at most %d times (change this using \"--num-scans\")\n\n", nscan);}

if(her==-9999){printf("All heritability estimates must be between 0.01 and %.4f (change the upper bound using \"--max-her\")\n\n", maxher);}

if(nmcmc==-9999)
{printf("Will use either three or ten random vectors for Monte Carlo operations (decided based on the number of samples); change this number using \"--num-random-vectors\"\n\n");}
}

if(strcmp(covarfile,"blank")==0&&strcmp(factorfile,"blank")==0){printf("Consider using \"--covar\" and/or \"--factors\" to provide quantitative or categorical covariates\n\n");}

print_qc(minmaf, maxmaf, minvar, minobs, mininfo, genprobs, 1);

if(verbose==0){printf("LDAK will only save selected files (to save all files, use \"--verbose YES\")\n\n");}
}

////////

if(mode==156)
{
if(strcmp(blockfile,"blank")==0)
{
if(window_cm!=-9999){printf("Calculating correlations between pairs of predictors within windows of size %.4fcM (change this size using \"--window-cm\" or \"--window-kb\")\n\n", window_cm);}
else{printf("Calculating correlations between pairs of predictors within %.2fkb (change this size using \"--window-cm\" or \"--window-kb\")\n\n", window_kb);}
}
else
{printf("Calculating correlations between pairs of predictors using the windows defined in %s (plus will automatically start a new window when the chromosome changes)\n\n", blockfile);}

if(minmaf>0){printf("Will exclude predictors with minor allele frequency below %.6f (change this threshold using \"--min-maf\")\n\n", minmaf);}

if(mincor>0){printf("Will set to zero correlations whose square is below %.2e (change this threshold using \"--min-cor\")\n\n", mincor);}

if(strip>0){printf("Will remove a proproption %.2f of the variation of each window correlation matrix (change this value using \"--strip\")\n\n", strip);}
}

if(mode==157)
{
printf("Joining %d correlations\n\n", num_cors);
}

if(mode==158)
{
printf("Generating partial summary statistics, corresponding to %.2f and %.2f of the total samples\n\n", subprop, 1-subprop);
}

if(mode==159)
{
printf("Constructing a MegaPRS PRS\n\n");

if(ptype==0){printf("will consider lasso, ridge regression, Bolt-LMM, BayesR and elastic net models; ");}
if(ptype==1){printf("will consider lasso-sparse models (use posterior modes instead of posterior means); ");}
if(ptype==2){printf("will consider lasso modelsc");}
if(ptype==3){printf("will consider ridge regression models;");}
if(ptype==4){printf("will consider Bolt models;");}
if(ptype==5){printf("will consider BayesR models;");}
if(ptype==6){printf("will consider BayesR-shrink models (replace point mass in prior with a Gaussian distribution);");}
if(ptype==7){printf("will consider elastic net models;");}
printf(" change this choice using \"--model\", followed by \"lasso-sparse\", \"lasso\", \"ridge\", \"bolt\", \"bayesr\", \"bayesr-shrink\", \"--elastic\" or \"--mega\" (we recommend \"--bayesr\")\n\n");

if(prsvar==1){printf("Will also construct %d jackknife models, each time reducing the sample size by %.2f\n\n", num_blocks, jackprop);}

if(strcmp(fracfile,"blank")==0)
{printf("Will use the default parameter choices (saved in the file %s.parameters); to instead specify your own, use \"--parameters\"\n\n", outfile);}
else{printf("Will use the %d sets of prior parameter choices provided in %s\n\n", countrows(fracfile), fracfile);}

if(skipcv==0)
{printf("Will select the best parameters via cross-validation, using %.2f randomly-picked test samples (use \"--cv-proportion\" to change this proportion or \"--cv-skip\" to turn off cross-validation)\n\n", cvprop);}

if(her==-9999){printf("The heritability estimates must be between 0.01 and %.4f (change the upper bound using \"--max-her\")\n\n", maxher);}

if(strcmp(ldfile,"blank")!=0)
{printf("Please note that you are no longer required to use \"--high-LD\", because LDAK can instead identify high LD regions from the correlations\n\n");}

if(checkfreq==1)
{printf("Will ignore predictors whose allele frequency differs significantly (P<1e-6) between the summary statistics and the correlations (to prevent this, use \"--check-frequencies NO\")\n\n");}

if(strcmp(sums2file,"blank")!=0)
{printf("We recommend that \"--summary\" provides full summary statistics and \"--summary2\" provides training summary statistics (the latter should be computed either using a subset of samples or with \"--pseudo-summaries\"). LDAK will then create training and full PRS; you can subsequently use \"--calc-scores\" to identify the most-accurate training model and the corresponding full model (the latter becomes the final model). For more details, see http://www.ldak.org/prediction\n\n");}

if(shrink<1){printf("Will scale predictor-predictor correlations by %.4f (change this value using \"--shrink\")\n\n", shrink);}
}

if(mode==160)
{
if(num_scores==1){printf("Measuring the accuracy of one PRS\n\n");}
else{printf("Measuring the accuracy of %d PRS\n\n", num_scores);}
}

///////////////////////////

if(mode==161)
{
if(strcmp(eigenfile,"blank")==0)
{printf("Calculating the top %d principal components of a kinship matrix\n\nIf you have previously decomposed the kinship matrix, use \"--eigen\" to extract the PCs from the eigen-decomposition\n\n", axes);}
else
{printf("Extracting the top %d principal components of a kinship matrix\n\n", axes);}
}

if(mode==162)
{
printf("Calculating predictor loadings for %d principal component axes\n\n", axes);
}

if(mode==163)
{
printf("Decomposing a kinship matrix\n\n");

if(strcmp(respfile,"blank")==0)
{printf("If you add \"--pheno\", samples with missing phenotypes will be excluded\n\n");}
}

if(mode==164)
{
printf("Adjusting a kinship matrix for %d covariates,  %d environmental variables and %d top predictors\n\n", num_quants, num_envs, num_tops);

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
}

if(mode==166)
{
printf("Truncating a kinship matrix; values below %.6f will be set to zero\n\n", cutoff);

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
}

if(mode==167)
{
printf("Dividing a kinship matrix based on its top %d principal components\n\n", axes);

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
}

if(mode==168)
{
printf("Squaring elements of a kinship matrix; good luck :)\n\n");

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
}

if(mode==169)
{
printf("Calculating a GxEMM-IID kinship matrix using %d environmental variables\n\n", num_envs);

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrix\n\n");}
}

if(mode==170)
{
printf("Calculating pairs of GxEMM-FREE kinship matrices for each of %d environmental variables\n\n", num_envs);

if(discenv==0)
{printf("If the environmental variables specify subgroups of samples, use \"--subgroups YES\"\n\n");}

if(kingz==0&&kinraw==0)
{printf("Add \"--kinship-gz YES\" and/or \"--kinship-raw YES\" to save gzipped and/or text versions of the kinship matrices\n\n");}
}

///////////////////////////

if(mode==171)
{
printf("Calculating predictor and sample statistics\n\n");
}

if(mode==172)
{
if(prsvar==0)
{
if(num_scores==1){printf("Calculating scores for one profile\n\n");}
else{printf("Calculating scores for %d profiles\n\n", num_scores);}
}
else{printf("Calculating a PRS and its corresponding SEs (the latter are based on %d jackknife PRS)\n\n", num_scores-1);}

if(power==0){printf("Please note that %s is assumed to contains raw effect sizes (e.g., those generated by \"--calc-blups\", \"--linear\", \"--ridge\", \"--bolt\", \"--bayesr\", \"--elastic\", \"--mega-prs\" or \"--calc-pca-loads\"); if it instead contains standardized effect sizes (e.g., those from \"--calc-posts\"), you should use \"--power -1\"\n\n", scorefile);}
else{print_scaling(power,hwestand);}

if(strcmp(respfile,"blank")==0&&strcmp(sumsfile,"blank")==0)
{printf("If you add \"--pheno\" (or \"--summary\"), LDAK will compute the correlation between scores and the phenotype\n\n");}

if(prsvar==0&&savecounts==0){printf("If you require counts (how many predictors contribute to each profile) add \"--save-counts\"\n\n");}
}

if(mode==173)
{
if(num_phenos==1){printf("Making one phenotype, with heritability %.4f and ", her);}
else{printf("Making %d phenotypes, each with heritability %.4f and ", num_phenos, her);}
if(num_causals==-1){printf("all predictors causal\n\n");}
else{printf("%d causal predictors\n\n", num_causals);}

print_scaling(power,hwestand);

if(strcmp(causalsfile,"blank")==0)
{printf("Causal predictors will be picked at random; if you would prefer to specify them, use \"--causals\"\n\n");}
else
{printf("Causal predictors are specified in %s\n\n", causalsfile);}

if(strcmp(effectsfile,"blank")==0)
{printf("Effect sizes (for scaled predictors) will be drawn from a Gaussian distribution; if you would prefer to specify them, use \"--effects\"\n\n");}
else
{printf("Effect sizes (for scaled predictors) are specified in %s\n\n", effectsfile);}

if(prev==-9999){printf("To generate a binary phenotype, use \"--prev\" to provide the prevalence of cases\n\n");}

if(bivar==-9999){printf("To generate pairs of correlated phenotypes, use \"--bivar\"\n\n");}
else{printf("Pairs of phenotypes will have correlation %.4f\n\n", bivar);}
}

if(mode==174)
{
printf("Making genotypes for %d individuals and %d SNPs, each with %.4f <= MAF <= %.4f and spanning %d chromosomes (change these values using \"--maf-low\", \"--maf-high\" and \"num-chr\")\n\n", num_inds, num_snps, maf1, maf2, nchrom);

if(famsize==1){printf("If you would like some individuals to be related, use \"--family-size\" and \"--relatedness\" (e.g., \"--family-size 2\" and \"--relatedness 0.5\" will generate pairs of full-siblings)\n\n");}
else{printf("There will be families of size %d with relatedness %.4f\n\n", famsize, closeness);}

if(pops==1){printf("If you would like to simulate multiple populations, use \"--populations\" (e.g., \"--populations 2\" will generate two independent subgroups)\n\n");}
else{printf("Will generate %d populations\n\n", pops);}
}

if(mode==175)
{
printf("Calculating correlations between predictors in %s and %s\n\n", predlista, predlistb);

if(savepairs==0){printf("LDAK will only report the average pairwise correlation for each predictor in %s; to report all pairwise correlations, add \"--save-pairs YES\"\n\n", predlista);}
}

if(mode==176)
{
printf("Estimating the standard deviation of measures of accuracy between predicted and observed values, using jackknifing with %d blocks\n\n", num_blocks);

if(auc==0){printf("Will compute the correlation, correlation squared, mean squared error and mean absolute error; if the phenotype is binary, add \"--AUC YES\" to also compute the area under curve\n\n");}
else{printf("Will compute the correlation, correlation squared, mean squared error, mean absolute error and area under curve\n\n");}
}

if(mode==177)
{
printf("Cutting samples into %d folds\n\n", num_folds);
}

if(mode==178)
{
if(countcols(likefile)==2){printf("Using a grid search to find the Gaussian distribution that best fits a set of datapoints and the corresponding likelihoods\n\n");}
else{printf("Using a grid search to find the two Gaussian distributions that best fit pairs of datapoints and the corresponding likelihoodsn\n\n");}

printf("Will consider %d values for the mean (ranging from %.4f to %.4f) and %d values for the sd (ranging from zero to %.4f); change these settings using \"--num-means\", \"--num-sds\", \"--min-mean\", \"--max-mean\", \"--max-sd\"\n\n", num_means, minmean, maxmean, num_sds, maxsd);

if(omitone==1){printf("Will compute correlations multiple times with one observation excluded (this is useful if one of the likelihoods is incorrect); to instead use all observations, use \"--omit-one NO\"\n\n");}
}

///////////////////////////

if(mode==181||mode==182||mode==183||mode==184||mode==185)
{
if(num_files==1){printf("Making a dataset in ");}
else{printf("Merging %d datasets and saving in ", num_files);}
if(mode==181){printf("binary PLINK format\n\n");}
if(mode==182){printf("SP format\n\n");}
if(mode==183){printf("sped (old speed) format\n\n");}
if(mode==184){printf("speed format\n\n");}
if(mode==185){printf("gzipped gen (chiamo) format\n\n");}

if(threshold!=-9999)
{
if(genprobs>1)
{
if(threshold!=0.5){printf("First, genotype probabilities will be converted to dosages (expected allele count), then dosages will be converted to hard genotypes as follows:\nDosage <= %f -> Genotype 0\n%f <= Dosage <= %f -> Genotype 1\nDosage >= %f -> Genotype 2\n\n", 1-threshold, threshold,2-threshold, 1+threshold);}
else{printf("First, genotype probabilities will be converted to dosages (expected allele count), then dosages will be converted to hard genotypes as follows:\nDosage < 0.5 -> Genotype 0\n0.5 <= Dosage <= 1.5 -> Genotype 1\nDosage > 1.5 -> Genotype 2\n\n");}
}
else
{
if(threshold!=0.5){printf("Values are assumed to be dosages, and will be converted to hard genotypes as follows:\nDosage <= %f -> Genotype 0\n%f <= Dosage <= %f -> Genotype 1\nDosage >= %f -> Genotype 2\n\n", 1-threshold, threshold,2-threshold, 1+threshold);}
else{printf("Values are assumed to be dosages, and will be converted to dosages (expected allele count), then dosages will be converted to hard genotypes as follows:\nDosage < 0.5 -> Genotype 0\n0.5 <= Dosage <= 1.5 -> Genotype 1\nDosage > 1.5 -> Genotype 2\n\n");}
}
}

if(minprob>=0.5)
{printf("Each genotype will be set to its most likely value, provided the certainty is at least %.2f\n\n", minprob);}
if(minprob==0)
{printf("Genotype values will be sampled at random based on state probabilities\n\n");}

if(encoding==2){printf("Will use a dominant coding of genotypes (0/1/2 becomes 0/2/2)\n\n");}
if(encoding==3){printf("Will use a recessive coding of genotypes (0/1/2 becomes 0/0/2)\n\n");}
if(encoding==4){printf("Will use a heterogeneous coding of genotypes (0/1/2 becomes 0/2/0)\n\n");}
if(encoding==5){printf("Will ensure the A1 allele has lower frequency than the A2 allele\n\n");}
if(encoding==6){printf("Predictor values will indicate whether the original genotype was missing or not\n\n");}

print_qc(minmaf, maxmaf, minvar, minobs, mininfo, genprobs, num_files);

if(mode==184)
{
if(speedlong==0)
{printf("Each predictor will be stored as one byte (accurate to 1/254); this is normally sufficient, but to increase this to 2 bytes (accurate to 1/65534) use \"--speed-long YES\"\n\n");}
else
{printf("Each predictor will be stored as two bytes (accurate to 1/65534); note that you can not use this format with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");}
}
}

////////

if(mode==186||mode==187||mode==188||mode==189)
{
if(mode==186){strcpy(writestring,"binary PLINK format\n");}
if(mode==187){strcpy(writestring,"SP format\n");}
if(mode==188){strcpy(writestring,"sped (old speed) format\n");}
if(mode==189){strcpy(writestring,"speed format\n");}

if(strcmp(genefile,"blank")!=0)
{
printf("Condensing predictors into genes based on the annotations provided in %s and saving in %s\n", genefile, writestring);

if(up_buffer>0||down_buffer>0||gene_buffer>0)
{
if(gene_buffer>0){printf("Will include predictors within %d basepairs of each gene\n\n", gene_buffer);}
else{printf("Will include predictors within %d (up) and %d (down) basepairs of each gene\n\n", up_buffer, down_buffer);}
}
else{printf("Use \"--gene-buffer\" (or \"--up-buffer\" and \"--down-buffer\") to include predictors neighbouring genes\n\n");}

if(minweight!=1e-10)
{
if(strcmp(weightsfile,"blank")==0){printf("Will only consider genes containing at least %.2f predictors\n\n", minweight);}
else{printf("Will only consider genes with weighting at least %.2f\n\n", minweight);}
}
else{printf("Will include all genes (use \"--min-weight\" to exclude small genes)\n\n");}

if(overlap==1){printf("Will allow overlap between adjacent genes (to prevent this, use \"--overlap NO\")\n\n");}
else{printf("Will ensure predictors are assigned to at most one gene\n\n");}

if(strcmp(pvafile,"blank")==0){printf("If you add \"--pvalues\", LDAK will work out the lowest p-value for each gene\n\n");}
}
else	//so using chunks or chunks bp
{
if(chunks!=-9999)
{
if(strcmp(weightsfile,"blank")!=0)
{
if(overlap==1){printf("Condensing predictors into overlapping chunks with weighting %.2f (to avoid overlap, add \"--overlap NO\")\n\n", chunks);}
else{printf("Condensing predictors into non-overlapping chunks with weighting %.2f\n\n", chunks);}
}
else
{
if(overlap==1){printf("Condensing predictors into overlapping chunks of length %.2f (to avoid overlap, add \"--overlap NO\")\n\n", chunks);}
else{printf("Condensing predictors into non-overlapping chunks of length %.2f\n\n", chunks);}
}
}
else
{
if(overlap==1){printf("Condensing predictors into overlapping chunks of length %d basepairs (to avoid overlap, add \"--overlap NO\")\n\n", chunksbp);}
else{printf("Condensing predictors into non-overlapping chunks of length %d basepairs\n\n", chunksbp);}
}

if(minweight!=1e-10)
{
if(strcmp(weightsfile,"blank")==0){printf("Will only consider chunks containing at least %.2f predictors\n\n", minweight);}
else{printf("Will only consider chunks with weighting at least %.2f\n\n", minweight);}
}

if(strcmp(pvafile,"blank")==0){printf("If you add \"--pvalues\", LDAK will work out the lowest p-value for each chunk\n\n");}
}

if(useminor==1)
{printf("SNPs will be coded based on the count of the minor allele (to instead always use the A1 allele use \"--count-minor NO\")\n\n");}

if(mode==186)
{printf("Warning, binary PLINK format can only store values 0, 1 or 2, so non-integer values will be rounded and those above 2 will be reduced\nTo avoid this, you should instead use \"--condense-sp\", \"--condense-sped\" or \"--condense-speed\"\n\n");}

if(mode==189)
{
if(speedlong==0)
{printf("Each predictor will be stored as one byte (accurate to 1/254); this is normally sufficient, but to increase this to 2 bytes (accurate to 1/65534) use \"--speed-long YES\"\n\n");}
else
{printf("Each predictor will be stored as two bytes (accurate to 1/65534); note that you can not use this format with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n");}
}
}

if(mode==190)
{
printf("Calculating correlations between two datasets\n\n");
}

///////////////////////////

if(mode==191)
{
{printf("Cutting predictors into partitions\n\n");}

if(strcmp(respfile,"blank")==0)
{printf("If you add \"--pheno\", samples with missing phenotypes will be excluded\n\n");}
}

if(mode==192)
{
printf("Calculating correlations for Partition %d (out of %d)\n\n", partition, num_parts);
}

if(mode==193)
{
printf("Joining up and decomposing the correlations\n\n");

printf("Will filter predictors so that no pair remains with correlation squared greater than %.4f (change this value using \"--max-cor\")\n\n", maxcor);
}

if(mode==194)
{
printf("Estimating heritability\n\n");

if(strcmp(covarfile,"blank")==0){printf("Consider using \"--covar\" to provide quantitative covariates\n\n");}
}

///////////////////////////

if(strcmp(workdir2,"blank")!=0){printf("The working directory is set to %s\n\n", workdir);}

if(seed!=-9999){printf("The random number generator seed is set to %d\n\n", seed);}

#if MKL==1
if(maxthreads==1){printf("To run the parallel version of LDAK, use \"--max-threads\" (this will only reduce runtime for some commands)\n\n");}
#endif

printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n\n");

///////////////////////////

