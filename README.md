# vcftochromopainter
Perl script for converting from vcf to chromopainterv2 format

// to compile:     gcc -Wall -o VCFToCPConvertNoFrills VCFToCPConvertNoFrills.c -lm -lz
// usage:  ./VCFToCPConvertNoFrills -g vcf.infile.gz -r recom.infile -o out.file

Map must be in format e.g.: 

Position(bp)	Rate(cM/Mb)	Map(cM)
55550	0	0
82571	2.98182894785537	0.080572
88169	2.08235083958557	0.092229
285245	1.76189388865209	0.439456
629218	3.01969049896358	1.478148
629241	2.86956521738693	1.478214
630053	2.88669950738932	1.480558
