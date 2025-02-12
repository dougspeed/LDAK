/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Declare variables

///////////////////////////

//parameters set directly by user

int dougvar=0;
double dougvar2=-9999;

//mode indicates which main argument is called (see mode.c for list)
int mode=-9999;

//random seed
int seed=-9999;

//data input

char folder2[500]="blank", folder[500]="";
char outfile2[500]="blank", outfile[500]="";

char workdir2[500]="blank", workdir[500]="";

char udatafile2[500]="blank", udatafile[500]="";
char ubimfile2[500]="blank", ubimfile[500]="";
char ufamfile2[500]="blank", ufamfile[500]="";
char unamefile2[500]="blank", unamefile[500]="";
char datalist2[500]="blank", datalist[500]="";

int oxchr=-9999;
int genskip=-9999;
int genheaders=-9999;
int genprobs=-9999;
int nonsnp=-9999;

int bedzeros[256], bedones[256], bedtwos[256];
double missingvalue=-9999;

//data filtering

char bsampfile2[500]="blank", bsampfile[500]="";
char csampfile2[500]="blank", csampfile[500]="";
int num_subs=-9999;
char subpref2[500]="blank", subpref[500];

char bpredfile2[500]="blank", bpredfile[500]="";
char cpredfile2[500]="blank", cpredfile[500]="";
int onechr=-9999;
char onesnp[500]="blank";

double minmaf=-9999;
double maxmaf=-9999;
double minvar=-9999;
double minobs=-9999;
double mininfo=-9999;

//data scaling (and pvalues) and coding

char centresfile2[500]="blank", centresfile[500]="";
char weightsfile2[500]="blank", weightsfile[500]="";
int ignoreweights=-9999;
double power=-9999;
int hwestand=-9999;
char pvafile2[500]="blank", pvafile[500]="";
char impfile2[500]="blank", impfile[500]="";

int encoding=-9999;
double threshold=-9999;
double minprob=-9999;

//kinships, regions, responses, summaries and fixed

char kinname2[500]="blank", kinname[500]="";
char kinlist2[500]="blank", kinlist[500]="";
int kindetails=-9999;

int num_regs=-9999;
char regpref2[500]="blank", regpref[500];
double rprune=-9999;

char respfile2[500]="blank", respfile[500]="";
int mpheno=-9999;
int mpheno2=-9999;
int pad=-9999;

char sumsfile2[500]="blank", sumsfile[500]="";
char sums2file2[500]="blank", sums2file[500]="";

int fixn=-9999, fixn2=-9999;
int amb=-9999;
double scaling=-9999, scaling2=-9999;

double prev=-9999;
double prev2=-9999;
double ascer=-9999;

char covarfile2[500]="blank", covarfile[500]="";
char covarnums2[2000]="blank", covarnums[2000]="blank";
char covarnames2[2000]="blank", covarnames[2000]="blank";
char envfile2[500]="blank", envfile[500]="";
char topfile2[500]="blank", topfile[500]="";
char factorfile2[500]="blank", factorfile[500]="";
char povarfile2[500]="blank", povarfile[500]="";

char offsetfile2[500]="blank", offsetfile[500]="";

//calculating weights, thinning and finding/removing tags

int nothin=-9999;
double wprune=-9999;
double window_kb=-9999, window_kb2;
int window_length=-9999, window_length2;
double window_cm=-9999;

double section_kb=-9999;
int section_length=-9999;
double section_cm=-9999;
double buffer_kb=-9999;
int buffer_length=-9999;
double buffer_cm=-9999;

int section=-9999;
int section_start=-9999;

int lddecay=-9999;
double halflife=-9999;

int fudge=-9999;
int simplex=-9999;
double maxtime=-9999;
int spread=-9999;

char targetfile2[500]="blank", targetfile[500];

//calculating and manipulating kinships (partitions also used for gene/chunk-based reml, sumher, gre)

int part_length=-9999;
int bychr=-9999;
int num_parts=-9999;
char partpref2[500]="blank", partpref[500];
int checkpart=-9999;

int partition=-9999;
int kingz=-9999;
int kinraw=-9999;
int single=-9999;

double dosage=-9999;
char malesfile2[500]="blank", malesfile[500];
int onlydets=-9999;
char invsfile2[500]="blank", invsfile[500];
int david=-9999;

double maxrel=-9999, minrel=-9999;
int kinstand=-9999;

int partial=-9999;

//reml, blup and he/pcgc (shortcut used indirectly by other functions)

int diagonal=-9999;
char hersfile2[500]="blank", hersfile[500]="";
int hestart=-9999;
int shortcut=-9999;
int discenv=-9999;
char oversfile2[500]="blank", oversfile[500]="";

char remlfile2[500]="blank", remlfile[500]="";

int adjusted=-9999;
double trun=-9999;

int num_vects=-9999;
int ldlt=-9999;

char relfile[500], relfile2[500]="blank";
int cordups=-9999;

//association analysis (genes also used for condensing data)

int kvikstep=-9999;
int faststep=-9999;
int gctastep=-9999;
char prsfile2[500]="blank", prsfile[500]="";
int verbose=-9999;

int spatest=-9999;
int num_knots=-9999;
int num_bins=-9999;
int spaside=1;
double spathresh=-9999;
double spamax=-9999;

int families=-9999;
int trios=-9999;
int duos=-9999;

int adjpreds=-9999;
char sampwfile2[500]="blank", sampwfile[500]="";
int sandwich=-9999;
int exact=-9999;
int scoretest=-9999;
char transfile2[500]="blank", transfile[500]="";

char genefile2[500]="blank", genefile[500]="";
double chunks=-9999;
int chunksbp=-9999;
int gene_buffer=-9999;
int up_buffer=-9999, down_buffer=-9999;
double minweight=-9999;
int overlap=-9999;

int gene_perms=-9999;
int saveall=-9999;
double gprune=-9999;
double limit=-9999;

int magma=-9999;
double cut1=-9999, cut2=-9999;
double gamp=-9999, gam1=-9999, gam2=-9999;

//sumher (annotations also used by fast he and pcgc

int num_anns=-9999;
char annpref2[500]="blank", annpref[500];
char genpreds2[500]="blank", genpreds[500];
char labfile2[500]="blank", labfile[500];
int backpart=-9999;
int allone=-9999;
int reduce=-9999;

char printfile2[500]="blank", printfile[500]="";
char herfile2[500]="blank", herfile[500]="";
int unbias=-9999;
int savemat=-9999;
int cover=-9999;
int fourdp=-9999;

char taglist2[500]="blank", taglist[500]="";
char matlist2[500]="blank", matlist[500]="";
char pathlist2[500]="blank", pathlist[500]="";
int checkdups=-9999;

char tagfile2[500]="blank", tagfile[500]="";
char pathfile2[500]="blank", pathfile[500]="";
char altfile2[500]="blank", altfile[500]="";
char cvsfile2[500]="blank", cvsfile[500]="";
char catfile2[500]="blank", catfile[500]="";
char taufile2[500]="blank", taufile[500]="";

int checksums=-9999;
int gcon=-9999, cept=-9999;

int ldsc=-9999;
int chisol=-9999;
int tagone=-9999;

int divide=-9999;
int uptaus=-9999;
char powfile2[500]="blank", powfile[500]="";

int plet=-9999;

char matfile2[500]="blank", matfile[500]="";

char expfile2[500]="blank", expfile[500]="";
double cvar=-9999;

//individual-level data prediction, then megaprs

char indhers2[500]="blank", indhers[500]="";
double herscale=-9999;

int loco=-9999;
int dichot=-9999;
int multi=-9999;
int fast=-9999;
int fprs=-9999;
int fastgwa=-9999;
char fastfile2[500]="blank", fastfile[500]="";

int skipcv=-9999;
double cvprop=-9999;
char bvsfile2[500]="blank", bvsfile[500]="";

int ldpred=-9999;
int pointmass=-9999;

int ndivs=-9999;
int nmcmc=-9999;
double maxher=-9999;

int checkped=-9999;
int nped=-9999;
double maithresh=-9999;
double maisig=-9999;
int ncal=-9999;
int ncomp=-9999;
double cthresh=-9999;
int nscan=-9999;
int revher=-9999;

char blockfile2[500]="blank", blockfile[500]="";

char corslist2[500]="blank", corslist[500]="";

double subprop=-9999;

char corname2[500]="blank", corname[500]="";
char pseudostem2[500]="blank", pseudostem[500]="";

int ptype=-9999;
char bestfile2[500]="blank", bestfile[500]="";

int checkld=-9999;
char ldfile2[500]="blank", ldfile[500]="";

char fracfile2[500]="blank", fracfile[500]="";

int checkfreq=-9999;
int prsvar=-9999;
double jackprop=-9999;

//pca, decompose, adjust-grm and others

int axes=-9999;
char pcastem2[500]="blank", pcastem[500]="";

int eigenraw=-9999;
char eigenfile2[500]="blank", eigenfile[500]="";

int noise=-9999;

//stats, scores, making phenotypes, snps, inflation, jackknifing, folds, find gaussian

char scorefile2[500]="blank", scorefile[500]="";
char cofile2[500]="blank", cofile[500]="";
int savecounts=-9999;
char finalfile2[500]="blank", finalfile[500]="";

int num_phenos=-9999;
int num_causals=-9999;
double her=-9999;
double cher=-9999;
double bivar=-9999;
double bivar2=-9999;
double bivar3=-9999;
char probsfile2[500]="blank", probsfile[500];
char causalsfile2[500]="blank", causalsfile[500];
char effectsfile2[500]="blank", effectsfile[500];

int num_inds=-9999;
int num_snps=-9999;
double maf1=-9999, maf2=-9999;
int nchrom=-9999;
int famsize=-9999;
double closeness=-9999;
int quads=-9999;
int pops=-9999;
double missrate=-9999;

char predlista[500], predlista2[500]="blank";
char predlistb[500], predlistb2[500]="blank";
int savepairs=-9999;

char jackfile[500], jackfile2[500]="blank";
char proffile[500], proffile2[500]="blank";
int auc=-9999;

int num_folds=-9999;

char likefile[500], likefile2[500]="blank";
int num_means=-9999;
int num_sds=-9999;
double minmean=-9999, maxmean=-9999;
double maxsd=-9999;
int omitone=-9999;

//making and condensing data

int exsame=-9999;
int exlong=-9999;
int speedlong=-9999;

int useminor=-9999;

//gre

int sinv=-9999;
char greout[500], greout2[500]="blank";

//common options

int checkroot=-9999;
double mincor=-9999;
double maxcor=-9999;
double cutoff=-9999;
int constrain=-9999;
int num_blocks=-9999;
int permute=-9999;
double shrink=-9999;
double strip=-9999;

int bitsize=-9999, bitsize2;
double tol=-9999;
int maxiter=-9999;
int memsave=-9999;

int manysamples=-9999;
int manypreds=-9999;

//threading option

int maxthreads=1;

///////////////////////////

//variables which are used for a specific purpose (and set fairly early on)

//storing data

int dtype=-9999, binary=-9999, famhead=-9999, extract=0, use_data, num_files=0;
char **datastems, **bimstems, **famstems;
char datafile[500], bimfile[500], famfile[500];

int num_samples=-9999, num_samples_use=-9999, *keepsamps, *keepsamps_males, **subindex;
int num_preds=-9999, num_preds_use=-9999, *keeppreds, data_length=-9999, *keeppreds_use;

int bgen_preds=-9999, bgen_samples=-9999, bgen_comp=-9999, bgen_layout=-9999, bgen_ids=-9999;
size_t *bgen_indexes;
int maxpreds=500000;	//the apriori max number of predictors in gen files

char **ids1, **ids2, **ids3, **ids4, **allids1, **allids2, **allids3;
int *chr, *allchr, num_chr, num_chr2, num_chr3, *chrindex, *chrindex2, *chrindex3;
double *cm, *bp, *cmbp, *allcm, *allbp, *chrprops;
char **preds, **allpreds;
char **along1, **along2, **allalong1, **allalong2, *al1, *al2, *allal1, *allal2;

unsigned char **data_char;
float *data_single, *data_single2;
double *data, *data2, *data3;
double **bytelookup;

float *speedstarts, *speedscales, **ps, **ps2;
double *centres, *mults, *sqdevs, *rates, *infos, *weights, *pvalues;

int num_pows;
double *powers, *powers2;

//kinships, regions, responses, summaries and fixed

int num_kins, *kinnums;
float *kins_single, *kins_single2, **Mkins_single, *inv_single;
double *kins, *kins2, **Mkins, *kintraces, *kinsums;
char **kinstems;

int *kindex;
double *kcentres, *kmults, *kweights, *kexps;
char **kpreds, *kal1, *kal2;

int rnum_preds_use, *rkeeppreds, **regindex, rdata_length;
double *rdata, *rcentres, *rmults, *rsqdevs,  *rrates, *rinfos, *rweights;
char **rpreds, *ral1, *ral2;

int num_resps, num_resps_use, *keepresps, *keepresps2, *respcounts, *respcounts2;
double *resp, *resp2;

int sformat, sformat2, gotfreq;
double *nss, *nss2, *nss3, *rnss, *tnss, *snss, *snss2;
double *chis, *chis2, *chis3, *rchis, *tchis, *schis, *schis2;
double *rhos, *rhos2, *rhos3, *rrhos, *trhos, *srhos, *srhos2;
double *a1freq;

int num_quants, num_cats, num_covars, *keepcovars, num_envs, num_fixed, num_prs;
double *covar, *thetas, *thetasds, *thetapvas, covher, topher;

int num_tops, *tkeeppreds, *tchr;
double *tbp, *tcentres, *tvars;
char **tpreds, *tal1, *tal2;

double *offsets;

//calculating weights, thinnings and finding/removing tags

int num_sections;
double decay;
int *sstarts, *sends, *sstarts2, *sends2;
float *cors_single;
double *cors, *cors2, *datasqs, *tally1, *tally2, *tally3, *tally4;
int *replace;

int num_seek, *sindex, *stypes;

//calculating and manipulating kinships

int *pstarts, *pends;
int *highs, *losts;
double *maxes;

//reml, blup and he/pcgc

char blupfile2[500]="blank", blupfile[500]="";
char regfile2[500]="blank", regfile[500]="";

int **nums;
double **blupcentres, **blupfactors, *bluprands, *blupvalues, *blupmids, *blupprojs;
double **guesses, ***guesses2;

int *bstarts, *bends;
double *R, *RTdata, *RTdata2, *MkinsD, *MkinsD2, *MkinsR, *MkinsR2;
double *KKtraces, *KKtraces2, *KKtraces3, *KYtraces, *KYtraces2;

float *bluprand_single, *guess_single;

//association analysis

int num_fams, *famindex, *famcounts;

double *Pincs, *Pscales, *Pscales2, *Pscales3, *Pvarphens, *Plambdas, *Proots, *Pspamax;

int *tindex, *spastatus;
double *knots, *bins, **CGF0, **CGF1, **CGF2, **CGF3;
double cuts[6]={.1,.01,.001,.0001,.00001,5e-8};
double *sweights;
double DTD[9], DTD2[3], DTD3[4], DTS[3], DTS2[3], invDTD[4];
double *nullprobs, *nullprobs2, *nullweights, *nullweights2;

int kvikparity;
double *prs, *covprs;

int num_genes, genemax;
int *gchr, *gstrand, *gstarts, *gends, *gparts;
double *gbp1, *gbp2, *gpvas;
char **gnames;
int nclump;
double *permlike;

//tagging, sumher and bayes factors

int **pindexes, addpart, addpart2, addgenic, *windex, *vindex;
double **pweights;

int num_tags, num_anals, num_reds, parttype=-9999;
char **tagstems, **matstems;
int *pcounts, *keepparts, *keepparts2;

double *exps, *exps2;
char **catlabels;

int num_parts_use;
int ncv, *cvindex;
double *stags, **svars, **svars2, **ssums, **ssums2, **ssums_use, *cvexps;
char **spreds,  *sal1, *sal2;

int topm, tops;
double *influs, *powlikes, *mvals, *svals, **perfs;

double *pmeans, *pvars, *bfs, *podds;

//individual-level data prediction, then megaprs

int num_train, num_test, *keeptrain, *keeptest;
int num_try, num_hers1, num_hers2, num_small, num_full;
int useped, highstruct, usecal, usecomp;

int *pedtops;
double *pedscales, *pedhers, *pedprs, *pedgammas, *pedsds, *pedmse;

int *eindex, *eindex2;
double *ecentres, *emults, *esqdevs, *erates, *einfos;

int *cindex, *cindex2;
double *cdata, *ccentres, *cmults, *csqdevs, *crates, *cinfos;

int *dindex, *dindex2;
double *ddata, *dcentres, *dmults, *dsqdevs, *drates, *dinfos;

int maxpairs, num_rels, *firsts, *seconds;
double *kinships, *cX, *cR;

double *lambdas, *lambdas2, *lambdas3, *lambdas4;
double *cgammas, *csds, *ceffs;
double *effs, *effs2, *effs3, *probs, *probs2, *residuals, *residuals2, *residuals3, *changes, *ess, *ess2;


int *bitrun, *bitdo, *bitactive, *bitdet1, *bitdet2;
double *bitdiffs, *bitpens;

int *Mtops, *Mbests, *Mincs;
double *Mscales, *Mmses, *Mmses2, *Mneffs;
double *Gmat, *Emat, *Umat, *Umat2, *Wmat, *Wmat2, *Wmat3, *Wmat4, *Wmat5, *Vmat, *Vmat2, *Ymat;

//double *cY, *cX, *cR, *cP, *cHP, *cKP, *cVP, *ctemps, *calphas, *cbetas;
//double *results1a, *results1b, *results2a, *results2b;

size_t *blockindexes;
int *blockstarts, *blockends, *blocklengths, **blockuse;
double *rjksums, *rjkaves, *rjktemp;
double *randnorms, *randnorms2, *datarands, *datarands2, *subrhos, *restrhos;

int num_cors;
char **corstems;

int *trytypes;
double *tryhers, *trylams, *tryscales, *tryps, *tryp2s, *tryp3s, *tryp4s, *tryf2s, *tryvars;

int *highlds;
double *loads, *loads2;

//pca, decompose, adjust-grm and others

float *Z_single, *ZTZ_single, *kinZ_single, *ZTkinZ_single, *W_single, *WZTkinZ_single;

//stats, scores, making phenotypes, snps, inflation, folds, jackknifing and cors

int *presents, *hets;
int num_scores=0;
double *indsds, *predmeans, *predvars, *predcors, *predsds;
double **unders, **phens, *vsds, *veffs, *vadds, **liabs;

int numa, *keepa, numb, *keepb;
double *acentres, *amults, *bcentres, *bmults, *asqdevs, *bsqdevs, *arates, *brates, *ainfos, *binfos;
char **apreds, **bpreds;

int *sZ1, *sZ2;
double *sX, *sXTX, *sXTX2, *sW, *sT;

//making and condensing data

int *Xcurrent;
int *Xnall, *Xnuse, **Xks, **Xks2;
int *XNall, *XNuse, **Xkp, **Xkp2;
char ***Xids1, ***Xids2, ***Xids3;
int **Xchr;
double **Xcm, **Xbp;
char ***Xpreds, ***Xalong1, ***Xalong2, **Xal1, **Xal2;
FILE **Xinput;
gzFile *Xdatainputgz;

int passqc;
int csamps, qmerge;
int erra, errb, errc, errd;
char **pid, **mid, **bid, **schar, **pchar;
unsigned char onechar;
unsigned short oneshort;
float minfloat, maxfloat;

int rowlength, c0, c1, c2;
unsigned char startchars[3], *rowchars, *rowchars2;

//speed tests

float *Rvsing, *Rvsing2, *Rvsing3, *kins_packed;
double *Rv, *Rv2;
size_t scount2;

//common options

int bit, bittotal, bittotal2, bitstart, bitstart2, bitend, bitend2, bitlength, bitlength2, bitmax, step;

double *Y, *Yadj, *Yadj2, *YTdata, *YTdata2, YTCY;
double *Z, *Z2, *Z3, *ZTZ, *ZTZ2, *ZTZ3, detZTZ, *ZTdata;

double *U, *U2, *E, *E2, *UTY, *UTZ, *UTdata;

int Xtotal, *Xstarts, *Xends;
double *X, *XTCX, *XTCX2, *XTCX3, *XTCX4, *Xbeta, *Xsums, *Xnss, *Xrhos, *Xsqs;

//threading options

int thread, threadstart, threadend, *Mcurrent;
float **Mdatatemp;
FILE **Minput;

///////////////////////////

//generic working variables

size_t scount, smax;
int i, i2, i3, j, j2, j3, k, k2, g, m, m2, p, p2, q, q2, q3, r, s, count, count2, count3, count4;
int current, head, found, total, total2, total3, token, indcount, ecount, wcount, xcount, *ycounts;
int shuffle, start, end, best, worst, mark, mark2, mark3, mark4, gen, gen2, *gens, flag, cflag, eflag, hflag, *order, *order2, cols[6];
double sum, sum2, sum3, sumsq, sumsq2, sumsq3, sumsq4, mean, mean2, mean3, var, var2, var3, value, value2, value3, value4, value5;
double last, min, max, med, maf, varphen, gif, postmean, unifrand;
double *hers, *hersold, *hersds, *shares, *sharesds, *cohers, *cohers2, mat[4], mat2[2];
double likenull, like, lrtstat, lrtpva, *likes, *likesold;
double neff, neff2, neff3, *pens, factor;
double weightsum, minpvalue, minpvalue2, *jacks, *gaussian, *polates, *polates2;
double *vstarts, *stats, *stats2, *stats3, *stats4, **effects, *varexp;

int one=1, two=2, three=3, info, lwork, *iwork, *ipiv, *ifail;
float alpha_single, beta_single;
double alpha, beta, wkopt, *work, vl, vu;

int idshead, *idsorder, *predorder, *usedids, *usedpreds, *indexer, *indexer2, *retain;
char idsfile[500]="blank", **kinids, **kinids2, **kinids3, **wantids, **wantids2;
char **wantpreds, **wantpreds2, **kinpreds;

int readint, readint2, readint3, *readints, writeint;
float *readfloats, writefloat, *writefloats;
double readdouble, readdouble2, *readdoubles;
char readchar, readchar2, *rs;
char **readstrings, readstring[500], readstring2[500], readstring3[500], readstring4[500], readstring5[500], readstring6[500];
char **writestrings, writestring[500];
char cmd[2000], cmd2[2000];

char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500], filename7[500], filename8[500], filename9[500], filename10[500];
FILE *input, *input2, *input3;
FILE *output, *output2, *output3, *output4, *output5, *output6, *output7, *output8, *output9, *output10;

DIR *dir;
gzFile datainputgz;
struct stat statstruct;

struct sorting_double *dptrs;
struct sorting_string *sptrs;

///////////////////////////

