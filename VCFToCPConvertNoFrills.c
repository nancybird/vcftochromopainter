#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include "zlib.h"

#define PI 3.141593

/*****************************************************************

// to compile:     gcc -Wall -o VCFToCPConvertNoFrills VCFToCPConvertNoFrills.c -lm -lz

// usage:  ./VCFToCPConvertNoFrills -g vcf.infile.gz -r recom.infile -o out.file

// example:  ./VCFToCPConvertNoFrills -g /lustre8/home/ghellenthal-pg/VCFinput/YGDP-freeze1.phased.ancientsmerged_ALLPOPS_quality05_chr22_filtered.vcf.gz -r /home/nancybird-pg/maps/chr22.SHAPEIT.b38.gmap -o /lustre8/home/ghellenthal-pg/CPinput/YGDP-freeze1.phased.ancientsmerged_ALLPOPS_quality05_chr22.cversion.cp

// example:  ./VCFToCPConvertNoFrills -g /lustre8/home/ghellenthal-pg/temp/YGDP-freeze1.phased.ancientsmerged.chr5.vcf.gz -r /home/nancybird-pg/maps/chr5.SHAPEIT.b38.gmap -f /lustre8/home/ghellenthal-pg/datainfo/JapanDataAscertain.popfile.txt -t /lustre8/home/ghellenthal-pg/datainfo/JapanDataALLPOPSQuality05.idfile.txt -s /lustre8/home/ghellenthal-pg/NancyNeanderthalLiftover/CombinedArchaicSNPList_chr5.txt -o /lustre8/home/ghellenthal-pg/CPinput/YGDP-freeze1.phased.ancientsmerged_ALLPOPS_quality05_chr5.cversion.cp
// example:  ./VCFToCPConvertNoFrills -g /lustre8/home/ghellenthal-pg/temp/YGDP-freeze1.phased.ancientsmerged.chr21.vcf.gz -r /home/nancybird-pg/maps/chr21.SHAPEIT.b38.gmap -f /lustre8/home/ghellenthal-pg/datainfo/JapanDataAscertain.popfile.txt -t /lustre8/home/ghellenthal-pg/datainfo/JapanDataALLPOPSQuality05.idfile.txt -s /lustre8/home/ghellenthal-pg/NancyNeanderthalLiftover/CombinedArchaicSNPList_chr21.txt -o /lustre8/home/ghellenthal-pg/CPinput/YGDP-freeze1.phased.ancientsmerged_ALLPOPS_quality05_chr21.cversion.cp

*******************************************************************/

int reading(st, format, res)
    char **st, *format;
    void *res;
{
    int i;
    char *rs;
    rs = *st;
    for(i = 0; isspace(rs[i]); i++) ; 
    if (!rs[i]) return 0; 
    for(; !isspace(rs[i]); i++) ;
    if (rs[i]) rs[i++] = 0;  
    if (!sscanf(*st, format, res)) return 0; 
    *st += i;
    return 1;
}

void usage()
{
  printf("to run: use './VCFToCPConvertNoFrills' with following options:\n");
  printf ("       -g <vcf.filein>  (VCF.gz file;REQUIRED; no default)\n");
  printf ("       -r <recommap.filein>  (REQUIRED; no default)\n");
  printf ("       -t <labels.filein>  file listing id labels to keep\n");
  printf ("       -s <snplist.filein>  file listing SNPs to keep\n");
  printf ("       -n <int>  read in this many inds at a time, to save RAM (default=1000)\n");
  printf ("       -o <outfile-prefix>  (default = 'vcf.filein')\n");
  printf ("       --help  print this menu\n");
}

int main(int argc, char *argv[])
{
  int i,j,a;
  int geno_find, recom_find, idfile_find, snplistfile_find, donorlist_find, outfile_find, numindsmaxfind, num_found;
  int nsitesFULL, nsites, ninds, num_splits, line_max, line_check, nsitesTOKEEP, num_tri, num_unordered, num_odd, to_keep, startIND, endIND, snp_count, allele1_int, allele2_int, start, endspot, nind_ID, nsites_SNPLIST, nindsTOKEEP, num_monomorph, snp_found, snplist_start, num_not_in_snplist, ind_count, ndonors, keep_indASCERTAIN, nindsASCERTAIN, snplistfile_ordered;
  double prevPOS, recom_current, allelefreq;
  char *step;
  char line[2047];
  char * line2=malloc(5000000*sizeof(char));
  char waste[2047];
  char waste2[2047];
  FILE *fd, *fd2, *fd3, *fd4, *fd5, *fout, *fout2;
  char * filename = malloc(1000 * sizeof(char));
  char * filenameGEN = malloc(1000 * sizeof(char));
  char * filenameID = malloc(1000 * sizeof(char));
  char * filenameSNPLIST = malloc(1000 * sizeof(char));
  char * filenameDONORLIST = malloc(1000 * sizeof(char));
  char * filenameOUT = malloc(1000 * sizeof(char));
  srand((unsigned)time(NULL));

  /***********************************************************/
  // DEFAULT VALUES:

  int numindsmaxsize=1000;     // ONLY READ IN THIS NUMBER OF VCF INDS AT A TIME (TO MINIMIZE RAM REQUIREMENTS -- MAKES PROGRAM RUN A FACTOR OF (N/numinds.maxsize) SLOWER THAN IF YOU SET numinds.maxsize=N, where N is the total number of individuals in your ".vcf" file, but also uses a factor of (N/numinds.maxsize) less RAM                   chr22 -- numinds.maxsize=1000 seems to use at least 5Gb

  double small_recom_val=0.000000000000001;    // lower limit for small genmap rates

  char allele_0 = '0';
  char allele_1 = '1';

  /************************************************************/

  double Mb = 1000000;


  geno_find=0;
  recom_find=0;
  idfile_find=0;
  donorlist_find=0;
  snplistfile_find=0;
  outfile_find=0;
  numindsmaxfind=0;
  num_found=0;
  for (i=1; i < argc; i++)
    {
      if ((strcmp(argv[i],"-help")==0) || (strcmp(argv[i],"--help")==0) || (strcmp(argv[i],"-h")==0))
	{
	  usage();
	  exit(1);
	}
      if (strcmp(argv[i],"-g")==0)
	{
	  geno_find=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-r")==0)
	{
	  recom_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-t")==0)
	{
	  idfile_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-f")==0)
	{
	  donorlist_find=1;
	  num_found=num_found+1;
	}
     if (strcmp(argv[i],"-s")==0)
	{
	  snplistfile_find=1;
	  num_found=num_found+1;
	}
      if (strcmp(argv[i],"-n")==0)
	{
	  numindsmaxfind=1;
	  num_found=num_found+1;
	}
       if (strcmp(argv[i],"-o")==0)
	{
	  outfile_find=1;
	  num_found=num_found+1;
	}
    }
  if (argc != (num_found*2+1))
    {
      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
      usage();
      exit(1);
    }
  if ((geno_find==0) || (recom_find==0)) { printf("Error with command line -- each of -g and -r MUST be specified. Exiting...\n"); exit(1);}
  if (idfile_find==0) strcpy(filenameID,"NULL");
  if (snplistfile_find==0) strcpy(filenameSNPLIST,"NULL");
  if (donorlist_find==0) strcpy(filenameDONORLIST,"NULL");

  for (i=1; i < argc; i++)
    {
      if (strcmp(argv[i],"-g")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
 	  strcpy(filenameGEN,argv[(i+1)]);
	  if (outfile_find==0)
	    {
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout = gzopen(strcat(filenameOUT,".haps.gz"), "w");
	      strcpy(filenameOUT,argv[(i+1)]);
	      fout2 = fopen(strcat(filenameOUT,".recomrates"), "w");
	    }
	}
       if (strcmp(argv[i],"-r")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  strcpy(filename,argv[(i+1)]);
	}
      if (strcmp(argv[i],"-t")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  strcpy(filenameID,argv[(i+1)]);
	}
     if (strcmp(argv[i],"-f")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  strcpy(filenameDONORLIST,argv[(i+1)]);
	}
      if (strcmp(argv[i],"-s")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  strcpy(filenameSNPLIST,argv[(i+1)]);
	}
       if (strcmp(argv[i],"-n")==0)
	{
	  if (argv[(i+1)][0] == '-')
	    {
	      printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	      usage();
	      exit(1);
	    }
	  numindsmaxsize = atoi(argv[(i+1)]);
	  if (numindsmaxsize < 0)
	    {
	      printf("Number of inds to consider at a time, i.e. '-n' switch, must be >0. Exiting...\n");
	      exit(1);
	    }
	}
        if (strcmp(argv[i],"-o")==0)
	 {
	   if (argv[(i+1)][0] == '-')
	     {
	       printf("Something wrong with input command line (missing arguments?). Exiting....\n");
	       usage();
	       exit(1);
	     }
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout = gzopen(strcat(filenameOUT,".haps.gz"), "w");
	   strcpy(filenameOUT,argv[(i+1)]);
	   fout2 = fopen(strcat(filenameOUT,".recomrates"), "w");
	 }
    }
 if (fout == NULL) {printf("error opening closing file1\n"); exit(1);}
 if (fout2 == NULL) {printf("error opening closing file2\n"); exit(1);}

                  // (I) get #SNPs, positions and recombination rates:
 fd2 = fopen(filename,"r");
 if (fd2 == NULL) { printf("error opening recom map input file: %s\n",filename); exit(1);}
 fgets(line,2047,fd2);   // header
 nsitesFULL=0;
 while(!feof(fd2))
   {
     if (fgets(line,2047,fd2)!=NULL) nsitesFULL=nsitesFULL+1;
   }
 fclose(fd2);
 double * positionsFULL = malloc(nsitesFULL * sizeof(double));
 double * recom_mapFULL = malloc((nsitesFULL - 1) * sizeof(double));
 fd2 = fopen(filename,"r");
 if (fd2 == NULL) { printf("error opening recom map input file: %s\n",filename); exit(1);}
 fgets(line,2047,fd2);   // header
 for (j=0; j < (nsitesFULL-1); j++)
   {
     fgets(line,2047,fd2);
     step=line;
     reading(&step,"%lf",&positionsFULL[j]);    // basepair position
     reading(&step,"%lf",&recom_mapFULL[j]);
     if (recom_mapFULL[j] >= 0 && recom_mapFULL[j] <= small_recom_val)
       {
	 printf("recom rate very low at basepair %lf (%lf). Assuming recomb rate between this snp and next one is %lf....\n",positionsFULL[j],recom_mapFULL[j],small_recom_val);
	 recom_mapFULL[j]=small_recom_val;
       }
     if (recom_mapFULL[j]<0)
       {
	 printf("recom rate < 0 at basepair %lf. Assuming transition probability of 1 between this snp and next one....\n",positionsFULL[j]);
       }
   }
 fgets(line,2047,fd2);
 step=line;
 reading(&step,"%lf",&positionsFULL[j]);    // basepair position
 fclose(fd2);
           // check ordering of SNPs (only allowed to be less than previous position if recom_map<0 at position -- i.e. suggesting new chromosome):
  for (i=0; i < nsitesFULL; i++)
    {
      if (i > 0)
	{
	  if (positionsFULL[i] <= positionsFULL[(i-1)] && (recom_mapFULL[(i-1)]>=0))
	    {
	      printf("positions in %s are not listed in increasing order at basepairs %lf and %lf. Exiting....\n",filenameGEN,positionsFULL[(i-1)],positionsFULL[i]);
	      exit(1);
	    }
	}
    }
  printf("\n Number of SNPs in %s: %d\n",filename,nsitesFULL);
  
                      // (II) get #SNPs and INDS in VCF FILE:
  fd = gzopen(filenameGEN,"r");
  if (fd == NULL) { printf("error opening %s\n",filenameGEN); exit(1);}
  while(!gzeof(fd))
    {
      line_check=0; if (gzgets(fd,line2,5000000)!=NULL) line_check=1;
      //line_check=0; if (gzgets(fd,line,2047)!=NULL) line_check=1;
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameGEN); exit(1);}
      step=line2;
      reading(&step,"%s",waste);
      if (strcmp("#CHROM",waste)==0) break;
   }
  for (i=0; i < 8; i++) reading(&step,"%s",waste);
  for (i=0; i < 10000000; i+=1)
    {
      reading(&step,"%s",waste);
      if (strcmp(waste,waste2)==0) break;
      strcpy(waste2,waste);
    }
  gzclose(fd);
  ninds=i;
  printf("ninds=%d\n",ninds);
  fd = gzopen(filenameGEN,"r");
  if (fd == NULL) { printf("error opening %s\n",filenameGEN); exit(1);}
  while(!gzeof(fd))
    {
      line_check=0; if (gzgets(fd,line2,5000000)!=NULL) line_check=1;
      //line_check=0; if (gzgets(fd,line,2047)!=NULL) line_check=1;
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameGEN); exit(1);}
      step=line2;
      reading(&step,"%s",waste);
      if (strcmp("#CHROM",waste)==0) break;
   }
  for (i=0; i<8; i++) reading(&step,"%s",waste);
  char ** id_vec=malloc(ninds*sizeof(char *));
  for (i=0; i < ninds; i++) id_vec[i]=malloc(1000*sizeof(char));
  for (i=0; i < ninds; i++)
    {
      reading(&step,"%s",waste);
      strcpy(id_vec[i],waste);
    }
  line_max=ninds*8*20;
  char * bigline = malloc(line_max * sizeof(char));
  nsites=0;
  while(!gzeof(fd))
   {
     if (gzgets(fd,bigline,line_max)!=NULL) nsites=nsites+1;
   }
  gzclose(fd);
  printf(" Number of SNPs in %s: %d\n",filenameGEN,nsites);
  printf(" Number of inds: %d\n",ninds);
  
                  // (III) IF NECESSARY, GET POPS TO ASCERTAIN ON:
  ndonors=1;
  if (donorlist_find==1)
    {
      fd5 = fopen(filenameDONORLIST,"r");
      if (fd5 == NULL) { printf("error opening %s\n",filenameDONORLIST); exit(1);}
      ndonors=0;
      while(!feof(fd5))
	{
	  if (fgets(line,2047,fd5)!=NULL)
	    ndonors=ndonors+1;
	}
      fclose(fd5);
    }
  char ** donor_pop_vec=malloc(ndonors * sizeof(char *));
  for (i=0; i < ndonors; i++)
    donor_pop_vec[i]=malloc(1000*sizeof(char));
  if (donorlist_find==1)
    {
      fd5 = fopen(filenameDONORLIST,"r");
      if (fd5 == NULL) { printf("error opening %s\n",filenameDONORLIST); exit(1);}
      for (i=0; i < ndonors; i++)
	{
	  fgets(line,2047,fd5);
	  step=line;
	  reading(&step,"%s",waste);
	  strcpy(donor_pop_vec[i],waste);
	}
      fclose(fd5);
    }
  
                  // (IV) IF NECESSARY, GET INDS TO KEEP: 
  nind_ID=1;
  if (idfile_find==1)
    {
                         // open id file (to get total number of inds)
     fd4 = fopen(filenameID,"r");
     if (fd4 == NULL) { printf("error opening %s\n",filenameID); exit(1);}
     nind_ID=0;
     while(!feof(fd4))
       {
	 if (fgets(line,2047,fd4)!=NULL)
	   nind_ID=nind_ID+1;
       }
     fclose(fd4);
     if (nind_ID != ninds) {printf("number of inds in %s (%d) does not match number of inds in %s (%d)\n. Exiting....\n",filenameID,nind_ID,filenameGEN,ninds); exit(1);}  
    }
  char ** ind_label_vec=malloc(nind_ID * sizeof(char *));
  char ** pop_label_vec=malloc(nind_ID * sizeof(char *));
  for (i=0; i < nind_ID; i++)
    {
      ind_label_vec[i]=malloc(1000*sizeof(char));
      pop_label_vec[i]=malloc(1000*sizeof(char));
    }
  int * include_ind_vec=malloc(ninds*sizeof(int));
  int * include_ind_vecASCERTAIN=malloc(ninds*sizeof(int));
  for (i=0; i < ninds; i++)
    {
      include_ind_vec[i]=1;
      include_ind_vecASCERTAIN[i]=1;
    }
  if (idfile_find==1)
    {
      fd4 = fopen(filenameID,"r");
      if (fd4 == NULL) { printf("error opening %s\n",filenameID); exit(1);}
      for (i=0; i < nind_ID; i++)
	{
	  fgets(line,2047,fd4);
	  step=line;
	  reading(&step,"%s",waste);
	  strcpy(ind_label_vec[i],waste);
	  reading(&step,"%s",waste);
	  strcpy(pop_label_vec[i],waste);
	  reading(&step,"%s",waste);
	  if (strcmp(waste,"0")==0)
	    {
	      include_ind_vec[i]=0;
	      include_ind_vecASCERTAIN[i]=0;
	    }
	  if (include_ind_vecASCERTAIN[i]==1 && donorlist_find==1)
	    {
	      keep_indASCERTAIN=0;
	      for (j=0; j < ndonors; j++)
		{
		  if (strcmp(donor_pop_vec[j],pop_label_vec[i])==0)
		    {
		      keep_indASCERTAIN=1;
		      break;
		    }
		}
	      include_ind_vecASCERTAIN[i]=keep_indASCERTAIN;
	    }
	}
      fclose(fd4);
      for (i=0; i < ninds; i++)
	{
	  if (strcmp(ind_label_vec[i],id_vec[i])!=0)
	    {
	      printf("ID labels for ind %d do not match between %s and %s (%s %s)!! Exiting...\n",i+1,filenameGEN,filenameID,id_vec[i],ind_label_vec[i]);
	      exit(1);
	    }
	}
    }
  nindsTOKEEP=0;
  nindsASCERTAIN=0;
  for (i=0; i < ninds; i++)
    {
      if (include_ind_vec[i]==1) nindsTOKEEP=nindsTOKEEP+1;
      if (include_ind_vecASCERTAIN[i]==1) nindsASCERTAIN=nindsASCERTAIN+1;
    }
  num_splits=nindsTOKEEP/numindsmaxsize; 
  if ((num_splits*numindsmaxsize)<nindsTOKEEP) num_splits=num_splits+1;
  if (donorlist_find==1)
    {
      printf(" Ascertaining on SNPs with MAF>0 among %d inds from the following populations:\n",nindsASCERTAIN);
      for (i=0; i < ndonors; i++) printf("        %s\n",donor_pop_vec[i]);
    }
  
                  // (V) IF NECESSARY, GET SNPs TO KEEP:
  nsites_SNPLIST=1;
  if (snplistfile_find==1)
    {
      fd3 = fopen(filenameSNPLIST,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filenameSNPLIST); exit(1);}
      nsites_SNPLIST=0;
      while(!feof(fd3))
	{
	  if (fgets(line,2047,fd3)!=NULL)
	    nsites_SNPLIST=nsites_SNPLIST+1;
	}
      fclose(fd3);
    }
  double * pos_vecSNPLIST = malloc(nsites_SNPLIST*sizeof(double));
  if (snplistfile_find==1)
    {
      fd3 = fopen(filenameSNPLIST,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filenameSNPLIST); exit(1);}
      for (i=0; i < nsites_SNPLIST; i++)
	{
 	  fgets(line,2047,fd3);
	  step=line;
	  reading(&step,"%s",waste);
	  reading(&step,"%lf",&pos_vecSNPLIST[i]);
	}
      fclose(fd3);
    }
             // check if SNPS are in order (which will drastically speed up program):
  snplistfile_ordered=1;
  for (i=0; i < (nsites_SNPLIST-1); i++)
    {
      if (pos_vecSNPLIST[i]>pos_vecSNPLIST[i+1])
	{
	  snplistfile_ordered=0;
	  break;
	}
    }

               // (VI) FIND IF ANY TRIALLELIC, MONOMORPHIC (IN KEPT OR ASCERTAINED INDS), OUT-OF-ORDER SNPS, OR SNPS WHERE DATA IS NO GOOD:
 double * pos_vecALL=malloc(nsites*sizeof(double));
 int * to_keepSNP=malloc(nsites*sizeof(int));
 fd = gzopen(filenameGEN,"r");
 if (fd == NULL) { printf("error opening %s\n",filenameGEN); exit(1);}
 while(!gzeof(fd))
   {
     line_check=0; if (gzgets(fd,bigline,line_max)!=NULL) line_check=1;
     if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameGEN); exit(1);}
     step=bigline;
     reading(&step,"%s",waste);
     if (strcmp("#CHROM",waste)==0) break;
   }
 prevPOS=-1.0;
 num_tri=0;
 num_unordered=0;
 num_odd=0;
 num_monomorph=0;
 num_not_in_snplist=0;
 nsitesTOKEEP=0;
 snplist_start=0;
 for (j=0; j< nsites; j++)
   {
     line_check=0; if (gzgets(fd,bigline,line_max)!=NULL) line_check=1;
     if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameGEN); exit(1);}
     step=bigline;
     reading(&step,"%s",waste);
     reading(&step,"%lf",&pos_vecALL[j]);
     to_keep=1;
     if (pos_vecALL[j]<=prevPOS)
       {
	 to_keep=0;
	 num_unordered=num_unordered+1;
       }
     if (to_keep==1 && snplistfile_find==1)
       {
	 snp_found=0;
	 for (i=snplist_start; i < nsites_SNPLIST; i++)
	   {
	     if (pos_vecALL[j]==pos_vecSNPLIST[i])
	       {
		 snplist_start=i;
		 snp_found=1;
		 break;
	       }
	     if (snplistfile_ordered==1 && pos_vecALL[j]<pos_vecSNPLIST[i]) break;
	   }
	 for (i=snplist_start; i >=0; i--)
	   {
	     if (pos_vecALL[j]==pos_vecSNPLIST[i])
	       {
		 snplist_start=i;
		 snp_found=1;
		 break;
	       }
	     if (snplistfile_ordered==1 && pos_vecALL[j]>pos_vecSNPLIST[i]) break;
	   }
	 if (snp_found==0)
	   {
	     num_not_in_snplist=num_not_in_snplist+1;
	     to_keep=0;
	   }
       }
     prevPOS=pos_vecALL[j];
     reading(&step,"%s",waste);
     reading(&step,"%s",waste);
     reading(&step,"%s",waste);
     if (strcmp(waste,"A")!=0 && strcmp(waste,"C")!=0 && strcmp(waste,"T")!=0 && strcmp(waste,"G")!=0)
       {
	 printf("TRI-SNP? %s %d %d %lf\n",waste,j,nsites,pos_vecALL[j]);
	 to_keep=0;
	 num_tri=num_tri+1;
       }
     if (to_keep==1)
       {
	 reading(&step,"%s",waste);
	 reading(&step,"%s",waste);
	 reading(&step,"%s",waste);
	 reading(&step,"%s",waste);
	 allelefreq=0.0;
	 for (i=0; i < ninds; i++)
	   {
	     if (include_ind_vec[i]==1)
	       {
		 reading(&step,"%s",waste);
		 char allele1=waste[0];
		 char allele2=waste[2];
		 //if ((strcmp(waste[0],"1")!=0 && strcmp(waste[0],"0")!=0) || (strcmp(waste[2],"1")!=0 && strcmp(waste[2],"0")!=0))
		 //if ((strcmp(allele1,"1")!=0 && strcmp(allele1,"0")!=0) || (strcmp(allele2,"1")!=0 && strcmp(allele2,"0")!=0))
		 //if (allele1==allele_1) printf("yes\n");
		 //if (allele2==allele_0) printf("yes2\n");
		 if ((allele1!=allele_0 && allele1!=allele_1) || (allele2!=allele_0 && allele2!=allele_1))
		   {
		     printf("Alleles for ind %d at SNP %d (position %.1lf) in %s are not 0 or 1! (%c %c). Removing this SNP...\n",i+1,j+1,pos_vecALL[j],filenameGEN,waste[0],waste[2]);
		     to_keep=0;
		     num_odd=num_odd+1;
		     break;
		   }
		 if (include_ind_vecASCERTAIN[i]==1)
		   {
		     allele1_int=(int)(allele1)-48;
		     allele2_int=(int)(allele2)-48;
		     allelefreq=allelefreq+allele1_int+allele2_int;
		   }
	       }
	   }
	 allelefreq=allelefreq/(2.0*nindsASCERTAIN);
	 if (allelefreq==0.0 || allelefreq==1.0)
	   {
	     num_monomorph=num_monomorph+1;
	     to_keep=0;
	   }
       }
     to_keepSNP[j]=to_keep;
     nsitesTOKEEP=nsitesTOKEEP+to_keepSNP[j];
   }
 gzclose(fd);
 printf(" Number of inds to keep = %d\n",nindsTOKEEP);
 printf(" number of (at least) triallelic SNPs = %d\n",num_tri);
 if (snplistfile_find==1) printf(" number of SNPs not in %s = %d\n",filenameSNPLIST,num_not_in_snplist);
 printf(" number of SNPs out-of-order = %d\n",num_unordered);
 printf(" number of SNPs with non-0/1 data for any kept ind = %d\n",num_odd);
 printf(" number of monomorphic SNPs among ascertained inds = %d\n",num_monomorph);
 printf(" total SNPs to keep = %d\n",nsitesTOKEEP);

            // (VII) READ IN VCF HAPLOTYPES AND MAKE CHROMOPAINTER HAPLOTYPE INPUT FILE:
 double * pos_vec = malloc(nsitesTOKEEP*sizeof(double));
 int ** geno_mat = malloc((2*numindsmaxsize)*sizeof(int *));
 for (i=0; i < (2*numindsmaxsize); i++) geno_mat[i]=malloc(nsitesTOKEEP*sizeof(int));
 gzprintf(fout,"%d\n",2*nindsTOKEEP);
 gzprintf(fout,"%d\n",nsitesTOKEEP);
 gzprintf(fout,"P");
 snp_count=0;
 for (j=0; j < nsites; j++)
   {
     if (to_keepSNP[j]==1)
       {
	 gzprintf(fout," %.0lf",pos_vecALL[j]);
	 pos_vec[snp_count]=pos_vecALL[j];
	 snp_count=snp_count+1;
       }
   }
 gzprintf(fout,"\n");
 for (a=0; a < num_splits; a++)
   {
     startIND=a*numindsmaxsize;
     endIND=(a+1)*numindsmaxsize;
     if (endIND > ninds) endIND=nindsTOKEEP;
     printf("Converting run %d of %d (inds %d to %d, of %d)\n",a+1,num_splits,startIND+1,endIND,nindsTOKEEP);

                                           // READ IN:
     fd = gzopen(filenameGEN,"r");
     if (fd == NULL) { printf("error opening %s\n",filenameGEN); exit(1);}
     while(!gzeof(fd))
       {
	 line_check=0; if (gzgets(fd,bigline,line_max)!=NULL) line_check=1;
	 if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameGEN); exit(1);}
	 step=bigline;
	 reading(&step,"%s",waste);
	 if (strcmp("#CHROM",waste)==0) break;
       }
     snp_count=0;
     for (j=0; j < nsites; j++)
       {
	 line_check=0; if (gzgets(fd,bigline,line_max)!=NULL) line_check=1;
	 if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameGEN); exit(1);}
	 if (to_keepSNP[j]==1)
	   {
	     step=bigline;
	     for (i=0; i<9; i++) reading(&step,"%s",waste);
	     ind_count=0;
	     for (i=0; i < ninds; i++)
	       {
		 reading(&step,"%s",waste);
		 if (include_ind_vec[i]==1)
		   {
		     if (ind_count>=startIND)
		       {
			 char allele1=waste[0];
			 char allele2=waste[2];
			 allele1_int=(int)(allele1)-48;
			 allele2_int=(int)(allele2)-48;
			 //printf("%c %c %d %d\n",allele1,allele2,allele1_int,allele2_int);
			 geno_mat[2*(ind_count-startIND)][snp_count]=allele1_int;
			 geno_mat[2*(ind_count-startIND)+1][snp_count]=allele2_int;
		       }
		     ind_count=ind_count+1;
		     if (ind_count>=endIND) break;
		   }
	       }
	     snp_count=snp_count+1;
	   }
       }
      gzclose(fd);

                                           // PRINT OUT:
      for (i=0; i < (2*(endIND-startIND)); i++)
	{
	  for (j=0; j < nsitesTOKEEP; j++) gzprintf(fout,"%d",geno_mat[i][j]);
	  gzprintf(fout,"\n");
	}
   }
 gzclose(fout);
 
            // (VIII) MAKE CHROMOPAINTER RECOM-MAP INPUT FILE:
 double * recomrate = malloc((nsitesTOKEEP-1)*sizeof(double));
 start=0;
 for (i=0; i < nsitesTOKEEP; i++)
   {
     if (pos_vec[i]>=positionsFULL[(nsitesFULL-1)]) recomrate[i]=recom_mapFULL[(nsitesFULL-2)]/Mb;
     if (pos_vec[i]< positionsFULL[(nsitesFULL-1)] && pos_vec[i] >= positionsFULL[0])
       {
	 for (j=start; j < (nsitesFULL-1); j++)
	   {
	     if ((pos_vec[i] >= positionsFULL[j]) && (pos_vec[i] < positionsFULL[(j+1)]) && (pos_vec[(i+1)] <= positionsFULL[(j+1)]))
	       {
		 recomrate[i] = recom_mapFULL[j]/Mb; 
		 start = j;
		 break;
	       }
	     if ((pos_vec[i] >= positionsFULL[j]) && (pos_vec[i] < positionsFULL[(j+1)]) && (pos_vec[(i+1)] > positionsFULL[(j+1)]))
	       {
		 recom_current = recom_mapFULL[j]*(positionsFULL[(j+1)]-pos_vec[i]);
		 endspot = j+1;
		 if (endspot == nsitesFULL)
		   {
		     recom_current = recom_current + recom_mapFULL[j]*(pos_vec[(i+1)]-positionsFULL[(j+1)]);
		     recomrate[i] = (recom_current/(pos_vec[(i+1)]-pos_vec[i]))/Mb;
		     break;
		   }
		 while(pos_vec[(i+1)] > positionsFULL[(endspot+1)])
		   {
		     recom_current = recom_current + recom_mapFULL[endspot]*(positionsFULL[(endspot+1)]-positionsFULL[endspot]);
		     endspot = endspot+1;
		     if (endspot == nsitesFULL) break;
		   }
		 recom_current = recom_current + recom_mapFULL[endspot]*(pos_vec[(i+1)]-positionsFULL[endspot]);
		 recomrate[i] = (recom_current/(pos_vec[(i+1)]-pos_vec[i]))/Mb;
		 start = j;
		 break;
	       }
	   }
       }
     if (pos_vec[i] < positionsFULL[0]) printf("WARNING - basepair %lf of file is less than first basepair of genetic map! Will 'impute' using the rate from the first included SNP(s) in the map!",pos_vec[i]);
     if((pos_vec[i] < positionsFULL[0]) && (pos_vec[(i+1)] <= positionsFULL[1])) recomrate[i] = recom_mapFULL[0]/Mb;
     if((pos_vec[i] < positionsFULL[0]) && (pos_vec[(i+1)] > positionsFULL[1]))
       {
	 recom_current = recom_mapFULL[0]*(positionsFULL[1]-pos_vec[i]);
	 endspot = 1;
	 while(pos_vec[(i+1)] > positionsFULL[(endspot+1)])
	   {
	     recom_current = recom_current + recom_mapFULL[endspot]*(positionsFULL[(endspot+1)]-positionsFULL[endspot]);
	     endspot = endspot+1;
	     if (endspot == (nsitesFULL-1)) break;
	   }
	 recom_current = recom_current + recom_mapFULL[endspot]*(pos_vec[(i+1)]-positionsFULL[endspot]);
	 recomrate[i] = (recom_current/(pos_vec[(i+1)]-pos_vec[i]))/Mb;
       }
   }
 fprintf(fout2,"start.pos recom.rate.perbp\n");
 for (i=0; i < nsitesTOKEEP; i++)
   {
     if (i<(nsitesTOKEEP-1)) fprintf(fout2,"%.0lf %.15lf\n",pos_vec[i],recomrate[i]/100);
     if (i==(nsitesTOKEEP-1)) fprintf(fout2,"%.0lf 0\n",pos_vec[i]);
   }
 fclose(fout2);
 
 free(filenameGEN);
 free(filenameID);
 free(filenameSNPLIST);
 free(filenameDONORLIST);
 free(filename);
 free(filenameOUT);
 free(positionsFULL);
 free(recom_mapFULL);
 free(include_ind_vec);
 free(include_ind_vecASCERTAIN);
 for (i=0; i < ninds; i++) free(id_vec[i]);
 free(id_vec);
 for (i=0; i < nind_ID; i++)
   {
     free(ind_label_vec[i]);
     free(pop_label_vec[i]);
   }
 free(ind_label_vec);
 free(pop_label_vec);
 for (i=0; i < ndonors; i++)
   free(donor_pop_vec[i]);
 free(donor_pop_vec);
 free(bigline);
 free(pos_vecALL);
 for (i=0; i < (2*numindsmaxsize); i++) free(geno_mat[i]);
 free(geno_mat);
 free(pos_vec);
 free(recomrate);
 free(to_keepSNP);
 free(pos_vecSNPLIST);
 free(line2);
 
 printf("Finished!\n");
 
 return 0;
}

