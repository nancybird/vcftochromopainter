## vcftochromopainter
*please note this is a development version and may have lots of bugs*

author: Garrett Hellenthal

script for converting from vcf to chromopainterv2 format

// to compile:     gcc -Wall -o VCFToCPConvert VCFToCPConvert.c -lm -lz

// usage:  ./VCFToCPConvert -g vcf.infile.gz -r recom.infile -o out.file


to run: use './VCFToCPConvertNoFrills' with following options:
       
       -g <vcf.filein>  (VCF.gz file;REQUIRED; no default)
       
       -r <recommap.filein>  (REQUIRED; no default)
       
       -t <labels.filein>  file listing id labels to keep
       
       -s <snplist.filein>  file listing SNPs to keep
       
       -n <int>  read in this many inds at a time, to save RAM (default=1000)
       
       -o <outfile-prefix>  (default = 'vcf.filein')
       
       --help  print this menu

Map must be in format e.g.: 

Position(bp)	Rate(cM/Mb)	Map(cM)

55550	0	0

82571	2.98182894785537	0.080572

88169	2.08235083958557	0.092229

