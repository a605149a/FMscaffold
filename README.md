#Compile

Before installing this version of StriDe, TBB library is needed. ( https://www.threadingbuildingblocks.org )
	
	(root)sudo yum install tbb
	1. ./autogen.sh 
	2. ./configure
	3. Add "-ltbb" into LIBS on "FMscaffold/StriDe/Makefile".
	4. make


#Parameter


Usage: FMscaffold/StriDe/stride scaffold {option} MatePairRead.fa

	option:
	-g, --asqg=NAME                  ASQGFILE.(skip mapping)
	-i, --insertSize=N               The insert size of matepair. (default: 3000)
	-n  --npair-threshold=N          The threshold of npair. (default: 1)
	-p, --prefix=NAME                prefix of FM-index of contigs. (bwt, rbwt, sai, rsai)
	-t, --threads=NUM                use NUM threads. (default: 1)
	-o, --out-prefix=NAME            use NAME as the prefix of the output files.
	--help                           display this help and exit.


#Example

	FMscaffold/StriDe/stride -t 30 index contig.fa
	FMscaffold/StriDe/stride scaffold -t 30 -i 3000 -p contig MatePairRead.fa

	
	(skip mapping)
	FMscaffold/StriDe/stride scaffold -g ASQGFILE MatePairRead.fa
