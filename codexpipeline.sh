#!/bin/bash
#Get options/arguments
function usage() { echo -e "Usage: $0 [-r|--R1 <R1 fastq file, required>] [-f|--R2 <R2 fastq file, required>] [-t|--tag <Barcode file, tabular, required>] [-p|-prj <project name directory>]" 1>&2; }

if [ $# -eq 0 ]
then
    usage
    exit 0
fi

function showHelp() {
cat << EOF

-h, --help              This help
-f, --R1 <file.fastq>   Fastq file unzipped with R1 reads; required
-r, --R2 <file.fastq>   Fastq file unzipped with R2 reads; required
-t, --tag <file.txt>    Barcode file; one line per well with the well name, the fragment 1 forward tag, then the fragment 1 reverse tag then the fragment 2 forward tag, then the fragm
ent 2 reverse tag separated by a tabulation; required
-p, --prj <directory>   Project name; Creates a directory were computing will be done; default BTC1
-m, --min <number>      Minimum number of reads for a fragment; default 10
-a, --minfr1 <number>   Minimum length for fragment 1; default 331
-b, --maxfr1 <number>   Maximum length for fragment 1; default 451
-c, --minfr2 <number>   Minimum length for fragment 2; default 313
-d, --maxfr2 <number>   Maximum length for fragment 2; default 434
-n, --minreads <number> Minimum reads matching a contig fragment; default 10
-e, --mincov <number>   Minimum contig fragment converage; default 5

First step - meerging R1 and R2 - need program PEAR under licence: https://www.h-its.org/downloads/pear-academic/
        Install Pear locally:
        Once downloaded, unzip the file pear-0.9.11-linux-x86_64.tar.gz by running:
                tar xzf pear-0.9.11-linux-x86_64.tar.gz
        then export it to your PATH:
                export PATH=\$PATH:\$PWD/pear-0.9.11-linux-x86_64/bin/

Second step - demultiplexing per well - no external program

Third step - demultiplexing per fragment type - need Blast - remove reads lower than 200 nucleotides

Four step - First step assembly per fragment - need trinity 

Fifth step - Second step assembly (assembly of both fragments) - need Blast - Contig filtering (see options)

EOF
}

OPTS=$(getopt -l "help,R1:,R2:,prj:,tag:,min:,minfr1:,maxfr1:,minfr2:,maxfr2:,minreads:,mincov:,ncbintdb:" -o "hf:r:p:t:m:a:b:c:d:n:e:i:" -a -- "$@")
eval set -- "$OPTS"

while true ; do
    case "$1" in
        -f|--R1)
            r1=$2
            if [ ! -e "${r1}" ]; then
                echo "requested file ${r1} doesn't exist" >&2
                exit 1
            fi
            shift;;
        -r|--R2)
            r2=$2
            if [ ! -e "${r2}" ]; then
                echo "requested file ${r2} doesn't exist" >&2
                exit 1
            fi
            shift;;
        -t|--tag)
            tag=$2
            if [ ! -e "${tag}" ]; then
                echo "requested file ${tag} doesn't exist" >&2
                exit 1
            fi
            shift;;
        -p|--prj)
            nameDir=$2
            shift;;
        -m|--min)
            ${minreads}=$2
            shift;;
        -a|--minfr1)
            ${minfr1}=$2
            shift;;
        -b|--maxfr1)
            ${maxfr1}=$2
            shift;;
        -c|--minfr2)
            ${minfr2}=$2
            shift;;
        -d|--maxfr2)
            ${maxfr2}=$2
            shift;;
        -n|--minreads)
            ${minreadsContig}=$2
            shift;;
        -e|--mincov)
            ${mincov}=$2
            shift;;
        -i|--ncbintdb)
            ${ncbintdb}=$2
            shift;;
        -h|--help)
            usage
            showHelp
            exit 0;;
        --)
            shift
            break;;
    esac
shift
done

if [ -z "${r1}" ] || [ -z "${r2}" ] || [ -z "${tag}" ]; then
    echo "Arguments are required. R1, R2 files and tag file."
    usage
    exit 0
fi

#variables needed
Libname=`echo $(basename $r1) | cut -d"_" -f1-2`

if [ -z $nameDir ]; then
    nameDir="BTC1"
fi

if [ -z ${minreads} ]; then
    minreads=10
fi

if [ -z ${minfr1} ]; then
    minfr1=$((331))
fi

if [ -z ${maxfr1} ]; then
    maxfr1=$((451))
fi

if [ -z ${minfr2} ]; then
    minfr2=$((313))
fi

if [ -z ${maxfr2} ]; then
    maxfr2=$((434))
fi

if [ -z ${minreadsContig} ]; then
    minreadsContig=10
fi

if [ -z ${mincov} ]; then
    mincov=5
fi

if [ -z ${ncbintdb} ]; then
    ncbintdb=/mnt/beegfs/shared_resources/NCBI_NT_DB/nt
fi

mkdir -p $nameDir
####First step: assembly of R1 and R2####
pearDir=$nameDir/First_step/
mkdir -p ${pearDir}
if [ ! -e  ${pearDir}/${Libname}.assembled.fastq ]; then
  #Run merging with default parameters
  pear -f $r1 -r $r2 -o ${pearDir}/$Libname
  #Summary
  echo "###########################################" >> ${pearDir}/${Libname}_step_1.txt
  echo $Libname >> ${pearDir}/${Libname}_step_1.txt
  echo "###########################################" >> ${pearDir}/${Libname}_step_1.txt
  echo "Total reads in R1 = "`zcat $r1 | wc -l | awk '{print $1/4}'` >> ${pearDir}/${Libname}_step_1.txt
  echo "Total assembled sequences = "`wc -l ${pearDir}/${Libname}.assembled.fastq | awk '{print $1/4}'` >> ${pearDir}/${Libname}_step_1.txt
  echo "Total reads unassembled = "`wc -l ${pearDir}/${Libname}.unassembled.forward.fastq | awk '{print $1/4}'` >> ${pearDir}/${Libname}_step_1.txt
fi

####Second step: demultiplexing, from Jawad python script "02a_codexSingle_v1.0.py"####
demDir=$nameDir/Second_step/
mkdir -p $demDir

if [ ! -e ${demDir}/${Libname}_step_2.txt ]; then
    cat ${pearDir}/${Libname}.assembled.fastq | paste - - - - | awk -F "\t" 'function revcomp(arg) {o = ""; for(i = length(arg); i > 0; i--){o = o c[substr(arg, i, 1)]} return(o)}
BEGIN{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\t"; while(getline < "'${tag}'" > 0){fwd[substr($2,1,5)]=1; fwd[substr($4,1,5)]; 
rev[substr($3,1,5)]=1; rev[substr($5,1,5)]=1; well[substr($2,1,5)"_"substr($3,1,5)]=$1; well[substr($3,1,5)"_"substr($2,1,5)]=$1; well[substr($4,1,5)"_"substr($5,1,5)]=$1; well[substr($5,1,5)"_"substr($4,1,5)]=$1}} \
{seqF=substr($2,1,5); seqR=revcomp(substr($2,length($2)-5+1,5));\
  if(fwd[seqF] && rev[seqR]){ \
    print $0 >> "'${demDir}'/'${Libname}'_"well[seqF"_"seqR]".fq";\
    FR += 1}\
  else if(fwd[seqR] && rev[seqF]){\
    print $0 >> "'${demDir}'/'${Libname}'_"well[seqR"_"seqF]".fq";\
    RF += 1}\
  else if(fwd[seqR] && ! rev[seqF]){\
    FnoR += 1;
    print $1 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $2 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $3 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $4 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
  }else if(! fwd[seqR] && rev[seqF]){\
    print $1 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $2 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $3 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $4 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    RnoF+= 1}\
  else if(! fwd[seqF] && rev[seqR]){\
    print $1 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $2 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $3 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $4 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    noFR += 1}\
  else if(! fwd[seqR] && rev[seqF]){\
    print $1 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $2 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $3 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $4 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    noRF += 1}\
  else{noFnoRnoRnoF += 1;\
    print $1 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $2 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $3 >> "'${demDir}'/'${Libname}'_nobarcode.fastq"; \
    print $4 >> "'${demDir}'/'${Libname}'_nobarcode.fastq";}
   }END{
    print "###########################################" >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "FR = \t\t\t "FR >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "RF = \t\t\t "RF >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "FnoR = \t\t\t "FnoR >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "RnoF = \t\t\t "RnoF >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "noFR = \t\t\t "noFR >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "noRF = \t\t\t "noRF >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "noFnoRnoRnoF = \t\t "noFnoRnoRnoF >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "totalOutcome = \t\t "FR+RF+FnoR+RnoF+noFR+noRF+noFnoRnoRnoF >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "###########################################" >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "total paires found = \t "FR+RF >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "total no paire found = \t "FnoR+RnoF+noFR+noRF+noFnoRnoRnoF >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "total seq = \t\t "NR >> "'$demDir'/'${Libname}'_step_2.txt";\
    print "###########################################" >> "'$demDir'/'${Libname}'_step_2.txt";}'
fi

####Third step: blast fragments on reads to segregate####
dem2Dir=$nameDir/Third_step/
mkdir -p ${dem2Dir}
#generate seq tags
#!/bin/bash

# Lookup table of degenerate IUPAC nucleotide codes.
declare -A deg2nuc=(
  ["R"]="A G"
  ["Y"]="C T"
  ["S"]="G C"
  ["W"]="A T"
  ["K"]="G T"
  ["M"]="A C"
  ["B"]="C G T"
  ["D"]="A G T"
  ["H"]="A C T"
  ["V"]="A C G"
  ["N"]="A C G T"
)

# Recursive function that replaces degenerate nucleotides with all combinations.
function generate {
    if [[ $1 =~ (.*)([RYSWiKMBDHVN])(.*) ]]; then
        local head=${BASH_REMATCH[1]}
        local tail=${BASH_REMATCH[3]}
        local -a seqs=()
        for nuc in ${deg2nuc[${BASH_REMATCH[2]}]}; do
            seqs+=($(generate "${head}${nuc}${tail}"))
        done
        echo "${seqs[@]}"
    else
        echo "$1"
    fi
}

# Generate the sequences and output in FASTA format
if [ ! -e  ${dem2Dir}/amorces.fa ]; then
    i=1
    for gen in `awk '{print $2}' ${tag} | sort -u`; do
        id="FR1_Fwd"
        sequence=$gen
        for seq in $(generate "$sequence"); do
            echo ">${id}_${i}" >> ${dem2Dir}/amorces.fa
            echo "$seq" >> ${dem2Dir}/amorces.fa
            ((i++))
        done
    done
    i=1
    for gen in `awk '{print $3}' ${tag} | sort -u`; do
        id="FR1_Rev"
        sequence=$gen
        for seq in $(generate "$sequence"); do
            echo ">${id}_${i}" >> ${dem2Dir}/amorces.fa
            echo "$seq" >> ${dem2Dir}/amorces.fa
            ((i++))
        done
    done
    i=1
    for gen in `awk '{print $4}' ${tag} | sort -u`; do
        id="FR2_Fwd"
        sequence=$gen
        for seq in $(generate "$sequence"); do
            echo ">${id}_${i}" >> ${dem2Dir}/amorces.fa
            echo "$seq" >> ${dem2Dir}/amorces.fa
            ((i++))
        done
    done
    i=1
    for gen in `awk '{print $5}' ${tag} | sort -u`; do
        id="FR2_Rev"
        sequence=$gen
        for seq in $(generate "$sequence"); do
            echo ">${id}_${i}" >> ${dem2Dir}/amorces.fa
            echo "$seq" >> ${dem2Dir}/amorces.fa
            ((i++))
        done
    done
fi
# makeblastdb
if [ ! -e  ${dem2Dir}/amorces.fa.nsq ]; then
    makeblastdb -dbtype nucl -in ${dem2Dir}/amorces.fa
fi
#Blastn tags on reads to segregate fragments
if [ ! -e ${dem2Dir}/${Libname}_step_3.txt ]; then
    for reads in `ls ${demDir}/${Libname}*.fq`;
    do
        well=$(basename `echo $reads | cut -d"_" -f4` .fq)
        if [ ! -e ${dem2Dir}/$(basename $reads .fq).fa ] && [ ! -e ${dem2Dir}/blast_$(basename $reads .fq) ]; then
            awk -F"\t" '{if($2>=200){print ">"$1; print $2;}}' $reads > ${dem2Dir}/$(basename $reads .fq).fa
            blastn -db ${dem2Dir}/amorces.fa -query ${dem2Dir}/$(basename $reads .fq).fa -outfmt 6 -word_size 7 -evalue 0.1 -perc_identity 95 | awk -F"\t" '
{OFS="\t"; split($2,strand,"_"); if(id==$1){if(strand[2]=="Fwd" && fwd<1){fr=fr","$2;pos=pos","$7"_"$8;fwd=fwd+1}\
if(strand[2]=="Rev" && rev<1){fr=fr","$2;pos=pos","$7"_"$8;rev=rev+1}}else{if(NR!=1) print id, fr, pos; \
if(strand[2]=="Fwd"){fwd=1;rev=0}else{rev=1;fwd=0}id=$1;fr=$2;pos=$7"_"$8;}}END{print id, fr, pos;}' | grep "," > ${dem2Dir}/blast_$(basename $reads .fq)
        fi
        awk -F "\t" 'function revcomp(arg) {o = ""; for(i = length(arg); i > 0; i--){o = o c[substr(arg, i, 1)]} return(o)}
        function revqual(qual){x=""; for(i=length(qual);i!=0;i--){x=x substr(qual,i,1)} return(x)}
        BEGIN{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; frag1=0; frag2=0; nofrag=0; OFS="\t"; 
        while(getline < "'${dem2Dir}/blast_$(basename $reads .fq)'" > 0){
            nb=split($2,fr,","); fragA[$1]=fr[1]; fragB[$1]=fr[2]; 
            split($3,position,","); posA[$1]=position[1]; posB[$1]=position[2]}}
        {split($1,s," "); seq=s[1]; split(fragA[seq],frA,"_"); split(fragB[seq],frB,"_"); split(posA[seq],psA,"_"); split(posB[seq],psB,"_");\
	if(frA[1] && frB[1]){
	    if(frA[1]==frB[1]){ \
                if(frA[1]=="FR1"){
                    if(psA[1]<psB[1]){
                        print $1 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR1.fastq"; \
                        print $2 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR1.fastq"; \
                        print $3 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR1.fastq"; \
                        print $4 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR1.fastq"; \
                    }else{
                        print $1 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR1.fastq"; \
                        print revcomp($2) >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR1.fastq"; \
                        print $3 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR1.fastq"; \
                        print revqual($4) >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR1.fastq"; \
                    }
                    frag1=frag1+1\
                }else{
                    if(psA[1]<psB[1]){
                        print $1 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR2.fastq"; \
                        print $2 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR2.fastq"; \
                        print $3 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR2.fastq"; \
                        print $4 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR2.fastq"; \
                    }else{
                        print $1 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR2.fastq"; \
                        print revcomp($2) >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR2.fastq"; \
                        print $3 >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR2.fastq"; \
                        print revqual($4) >> "'${dem2Dir}'/'$(basename $reads .fq)'_FR2.fastq"; \
                    }
                    frag2=frag2+1}
            }else{\
                nofrag=nofrag+1;}\
        }else{
             nofrag=nofrag+1;}\
        }END{\
        print "'${Libname}' '$well'" >> "'$dem2Dir'/'${Libname}'_step_3.txt";\
        print "###########################################" >> "'$dem2Dir'/'${Libname}'_step_3.txt";\
        print "Assigned fragments 1 = \t "frag1 >> "'$dem2Dir'/'${Libname}'_step_3.txt";\
        print "Assigned fragments 2 = \t "frag2 >> "'$dem2Dir'/'${Libname}'_step_3.txt";\
        print "Not Assigned fragments = \t "nofrag >> "'$dem2Dir'/'${Libname}'_step_3.txt";\
        print "###########################################" >> "'$dem2Dir'/'${Libname}'_step_3.txt";}' $reads
    done
fi
####Forth step : first step assembly####
trynDir=$nameDir/Forth_step/
mkdir -p ${trynDir}
if [ ! -e ${trynDir}/${Libname}_step_4.log ]; then
        for well in `awk '{print $1}' $tag`;
        do
        fragment1=${dem2Dir}/${Libname}_${well}_FR1.fastq
        fragment2=${dem2Dir}/${Libname}_${well}_FR2.fastq
        if [ -e $fragment1 ] && [ -e $fragment2 ]; then
                nbreadsFR1=`wc -l ${fragment1} | awk '{print $1/4}'`
                nbreadsFR2=`wc -l ${fragment2} | awk '{print $1/4}'`
                if [ $nbreadsFR1 -lt ${minreads} ]; then
                if [ $nbreadsFR2 -lt ${minreads} ]; then
                        rm ${fragment1} ${fragment2}
                        echo "Less than 10 reads found for each fragment. The well $well was removed from analysis" >> ${trynDir}/${Libname}_step_4.log
                else
                        rm ${fragment1}
                        echo "Less than 10 reads found for fragment 1; well $well. The assembly will be done only for fragment 2" >> ${trynDir}/${Libname}_step_4.log
                        cat ${fragment2} | paste - - - - | awk -F"\t" '{print $1; print $2}' | sed 's/@/>/' > ${trynDir}/${Libname}_${well}_FR2.fa
                        Trinity --seqType fa --max_memory 30G --single ${trynDir}/${Libname}_${well}_FR2.fa --CPU 2 --output ${trynDir}/trinity_${Libname}_${well}_FR2 --full_cleanup
                        minimap2 -ax sr ${trynDir}/trinity_${Libname}_${well}_FR2.Trinity.fasta ${trynDir}/${Libname}_${well}_FR2.fa | grep -v "^@" | awk '$2!=4 {print $3}' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2,($1/'$nbreadsFR2')*100}' > ${trynDir}/trinity_${Libname}_${well}_FR2.nbreads
                        blastn -db ${ncbintdb} -query ${trynDir}/trinity_${Libname}_${well}_FR2.Trinity.fasta -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -out ${trynDir}/trinity_${Libname}_${well}_FR2.nt -num_threads 2 -num_alignments 1
                fi
                else
                if [ $nbreadsFR2 -lt ${minreads} ]; then
                        rm $fragment2
                        echo "Less than 10 reads found for fragment 2; well $well. The assembly will be done only for fragment 1" >> ${trynDir}/${Libname}_step_4.log
                        cat ${fragment1} | paste - - - - | awk -F"\t" '{print $1; print $2}' | sed 's/@/>/' > ${trynDir}/${Libname}_${well}_FR1.fa
                        Trinity --seqType fa --max_memory 30G --single ${trynDir}/${Libname}_${well}_FR1.fa --CPU 2 --output ${trynDir}/trinity_${Libname}_${well}_FR1 --full_cleanup
                        minimap2 -ax sr ${trynDir}/trinity_${Libname}_${well}_FR1.Trinity.fasta ${trynDir}/${Libname}_${well}_FR1.fa | grep -v "^@" | awk '$2!=4 {print $3}' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2,$1,($1/'$nbreadsFR1')*100}' > ${trynDir}/trinity_${Libname}_${well}_FR1.nbreads
                        blastn -db ${ncbintdb} -query ${trynDir}/trinity_${Libname}_${well}_FR1.Trinity.fasta -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -out ${trynDir}/trinity_${Libname}_${well}_FR1.nt -num_threads 2 -num_alignments 1
                else
                        cat ${fragment1} | paste - - - - | awk -F"\t" '{print $1; print $2}' | sed 's/@/>/' > ${trynDir}/${Libname}_${well}_FR1.fa
                        cat ${fragment2} | paste - - - - | awk -F"\t" '{print $1; print $2}' | sed 's/@/>/' > ${trynDir}/${Libname}_${well}_FR2.fa
                        Trinity --seqType fa --max_memory 30G --single ${trynDir}/${Libname}_${well}_FR1.fa --CPU 2 --output ${trynDir}/trinity_${Libname}_${well}_FR1 --full_cleanup
                        minimap2 -ax sr ${trynDir}/trinity_${Libname}_${well}_FR1.Trinity.fasta ${trynDir}/${Libname}_${well}_FR1.fa  | grep -v "^@" | awk '$2!=4 {print $3}' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2,$1,($1/'$nbreadsFR1')*100}' > ${trynDir}/trinity_${Libname}_${well}_FR1.nbreads
                        blastn -db ${ncbintdb} -query ${trynDir}/trinity_${Libname}_${well}_FR1.Trinity.fasta -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -out ${trynDir}/trinity_${Libname}_${well}_FR1.nt -num_threads 2 -num_alignments 1
                        Trinity --seqType fa --max_memory 30G --single ${trynDir}/${Libname}_${well}_FR2.fa --CPU 2 --output ${trynDir}/trinity_${Libname}_${well}_FR2 --full_cleanup
                        minimap2 -ax sr ${trynDir}/trinity_${Libname}_${well}_FR2.Trinity.fasta ${trynDir}/${Libname}_${well}_FR2.fa  | grep -v "^@" | awk '$2!=4 {print $3}' | sort | uniq -c | sort -k1,1nr | awk '{OFS="\t"; print $2,$1,($1/'$nbreadsFR2')*100}' > ${trynDir}/trinity_${Libname}_${well}_FR2.nbreads
                        blastn -db ${ncbintdb} -query ${trynDir}/trinity_${Libname}_${well}_FR2.Trinity.fasta -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -out ${trynDir}/trinity_${Libname}_${well}_FR2.nt -num_threads 2 -num_alignments 1
                fi
                fi
        elif [ ! -e $fragment1 ] && [ ! -e $fragment2 ]; then
                echo "No reads found in well $well." >> ${trynDir}/${Libname}_step_4.log
        elif [ ! -e $fragment1 ]; then
                echo "No reads found in well $well fragment 1." >> ${trynDir}/${Libname}_step_4.log
        elif [ ! -e $fragment2 ]; then
                echo "No reads found in well $well fragment 2." >> ${trynDir}/${Libname}_step_4.log
        fi
    done
fi
####Fifth step : remove tags####
trynDir2=$nameDir/Fifth_step/
mkdir -p ${trynDir2}
for well in `awk '{print $1}' $tag`;
do
    #Fragment1: remove tag and generate stat tsv file
    frag1=${trynDir}/trinity_${Libname}_${well}_FR1.Trinity.fasta
    if [ -e ${frag1} ]; then
        #Remove tags off for FR1
        if [ ! -e ${trynDir2}/${Libname}_${well}_FR1.blastout ]; then
            blastn -db ${dem2Dir}/amorces.fa -query ${frag1} -outfmt 6 -word_size 7 -evalue 0.1 -perc_identity 95 | \
            awk '{if(id==$1){if(deb!=1){if($8<50){deb=1; print $0}} \
            if(fin!=1){if($8>50){fin=1; print $0}}}else{id=$1; deb=0; fin=0; if(id==$1){if(deb!=1){\
            if($8<50){deb=1; print $0}} if(fin!=1){if($8>50){fin=1; print $0}}}}}' | awk '{if(id==$1){cpt=cpt+1; l2=$0}else{if(cpt==2){print l1; print l2;} cpt=1; l1=$0; id=$1}}END{if(cpt==2){print l1; print l2;}}' \
            > ${trynDir2}/${Libname}_${well}_FR1.blastout
        fi
        if [ ! -e ${trynDir2}/${Libname}_${well}_FR1.notag.fa ]; then
            cat ${frag1} | paste - - | sed 's/>//' | awk -F "\t" 'function revcomp(arg) {o = ""; for(i = length(arg); i > 0; i--){o = o c[substr(arg, i, 1)]} return(o)}
            BEGIN{OFS="\t";c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ;\
            while(getline < "'${trynDir2}'/'${Libname}'_'${well}'_FR1.blastout" > 0){\
            if($8<50){deb[$1]=$8}else{fin[$1]=$7} split($2,t,"_"); if(($10>$9 && t[2]=="Rev")||($9>$10 && t[2]=="Fwd")){rev[$1]="yes"}else{rev[$1]="no"}}}
            {split($1,name," ");\
              if(deb[name[1]] && fin[name[1]]){
		print ">FR1_"$1;
                if(rev[$1]=="no"){
                    print substr($2,deb[name[1]]+1,length($2)-(deb[name[1]]+(length($2)-fin[name[1]]+1)))\
                }else{
                    print revcomp(substr($2,deb[name[1]]+1,length($2)-(deb[name[1]]+(length($2)-fin[name[1]]+1))))\
              }
            }}' >> ${trynDir2}/${Libname}_${well}_FR1.notag.fa
        fi
    fi
    #Generate stat tsv file
    if [ -e ${trynDir2}/${Libname}_${well}_FR1.notag.fa ]; then
        cat ${trynDir2}/${Libname}_${well}_FR1.notag.fa | paste - - | awk -F"\t" 'BEGIN{
        while(getline < "'${trynDir}'/trinity_'${Libname}'_'${well}'_FR1.nbreads" > 0){nbreads[">FR1_"$1]=$2; 
        perreads[">FR1_"$1]=$3;}
        while(getline < "'${trynDir}'/trinity_'${Libname}'_'${well}'_FR1.nt"  > 0){nt[">FR1_"$1]=$10} cpt=0} \
        {OFS="\t"; split($1,t," "); split(t[2],p,"="); l=p[2]; nreads=nbreads[t[1]]; 
        preads=perreads[t[1]]; ntmatch=nt[t[1]]; seq=$2; \
        print "'${Libname}'","'${well}'","FR1",l,nreads,preads,$1,seq >> "'${trynDir2}'/'${Libname}'_length_and_nbreads.tsv";\
        if(l>'$minfr1' && l<'$maxfr1' && nreads>'${minreadsContig}' && preads>'${mincov}'){\
        print "'${Libname}'","'${well}'","FR1",l,nreads,preads,$1,seq >> "'${trynDir2}'/'${Libname}'_length_and_nbreads.Filtered.tsv"; \
        print $1; print seq}}' >> ${trynDir2}/${Libname}_${well}_filtered.fa
        nblines=`wc -l ${trynDir2}/${Libname}_${well}_filtered.fa | awk '{print $1}'`
        if [ $nblines -eq 0 ]; then
            rm ${trynDir2}/${Libname}_${well}_filtered.fa
        fi
    fi
    #Fragment2: remove tag and generate stat tsv file
    frag2=${trynDir}/trinity_${Libname}_${well}_FR2.Trinity.fasta
    if [ -e ${frag2} ]; then
        #Remove tags off for FR2
        if [ ! -e ${trynDir2}/${Libname}_${well}_FR2.blastout ]; then
            blastn -db ${dem2Dir}/amorces.fa -query ${frag2} -outfmt 6 -word_size 7 -evalue 0.1 -perc_identity 95 | \
            awk '{if(id==$1){if(deb!=1){if($8<50){deb=1; print $0}} \
            if(fin!=1){if($8>50){fin=1; print $0}}}else{id=$1; deb=0; fin=0; if(id==$1){if(deb!=1){\
            if($8<50){deb=1; print $0}} if(fin!=1){if($8>50){fin=1; print $0}}}}}' | awk '{if(id==$1){cpt=cpt+1; l2=$0}else{if(cpt==2){print l1; print l2;} cpt=1; l1=$0; id=$1}}END{if(cpt==2){print l1; print l2;}}' \
            > ${trynDir2}/${Libname}_${well}_FR2.blastout
        fi
        if [ ! -e ${trynDir2}/${Libname}_${well}_FR2.notag.fa ]; then
            cat ${frag2} | paste - - | sed 's/>//' | awk -F "\t" 'function revcomp(arg) {o = ""; for(i = length(arg); i > 0; i--){o = o c[substr(arg, i, 1)]} return(o)}
            BEGIN{OFS="\t";c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ;\
            while(getline < "'${trynDir2}'/'${Libname}'_'${well}'_FR2.blastout" > 0){\
            if($8<50){deb[$1]=$8}else{fin[$1]=$7} split($2,t,"_"); if(($10>$9 && t[2]=="Rev")||($9>$10 && t[2]=="Fwd")){rev[$1]="yes"}else{rev[$1]="no"}}}
            {split($1,name," ");\
              if(deb[name[1]] && fin[name[1]]){
                print ">FR2_"$1;
                if(rev[$1]=="no"){
                    print substr($2,deb[name[1]]+1,length($2)-(deb[name[1]]+(length($2)-fin[name[1]]+1)))\
                }else{
                    print revcomp(substr($2,deb[name[1]]+1,length($2)-(deb[name[1]]+(length($2)-fin[name[1]]+1))))\
              }
            }}' >> ${trynDir2}/${Libname}_${well}_FR2.notag.fa
        fi
    fi
    #Generate stat tsv file
    if [ -e ${trynDir2}/${Libname}_${well}_FR2.notag.fa ]; then
        cat ${trynDir2}/${Libname}_${well}_FR2.notag.fa | paste - - | awk -F"\t" 'BEGIN{
        while(getline < "'${trynDir}'/trinity_'${Libname}'_'${well}'_FR2.nbreads" > 0){nbreads[">FR2_"$1]=$2; 
        perreads[">FR2_"$1]=$3;}
        while(getline < "'${trynDir}'/trinity_'${Libname}'_'${well}'_FR2.nt"  > 0){nt[">FR1_"$1]=$10} cpt=0} \
        {OFS="\t"; split($1,t," "); split(t[2],p,"="); l=p[2]; nreads=nbreads[t[1]]; 
        preads=perreads[t[1]]; seq=$2; ntmatch=nt[t[1]]; seq=$2;\
        print "'${Libname}'","'${well}'","FR2",l,nreads,preads,$1,seq >> "'${trynDir2}'/'${Libname}'_length_and_nbreads.tsv";\
        if(l>'$minfr2' && l<'$maxfr2' && nreads>'${minreadsContig}' && preads>'${mincov}'){\
        print "'${Libname}'","'${well}'","FR2",l,nreads,preads,$1,seq >> "'${trynDir2}'/'${Libname}'_length_and_nbreads.Filtered.tsv";\
        print $1; print seq}}' >> ${trynDir2}/${Libname}_${well}_filtered.fa
        nblines=`wc -l ${trynDir2}/${Libname}_${well}_filtered.fa | awk '{print $1}'`
        if [ $nblines -eq 0 ]; then
            rm ${trynDir2}/${Libname}_${well}_filtered.fa
        fi
    fi
done

cp ${trynDir2}/${Libname}_length_and_nbreads.tsv ${nameDir}/
echo -e "Library\tWell\tFragment\tLength\tNb of reads matching contig\tNb of reads matching contig / Nb of total read used for assembly\tMatch on NT\tsequence ID\tSequence" > ${nameDir}/${nameDir}_${Libname}_ALL_length_and_nbreads.Filtered.tsv
cat ${trynDir2}/${Libname}_length_and_nbreads.Filtered.tsv >> ${nameDir}/${nameDir}_${Libname}_ALL_length_and_nbreads.Filtered.tsv

####Final step: merging fragments####
#Get easy cases
awk -F"\t" '{if(id!=$1"_"$2){if(NR!=1){if(cpt==1){print ">'$nameDir'_"id"_"fr >> "'${nameDir}'/'${Libname}'_final.sequences.fa"; print $8 >> "'${nameDir}'/'${Libname}'_final.sequences.fa"
                          }else{if(cpt==2){print l1 >> "'${trynDir2}'/'${Libname}'case2contigs.out"; print l2 >> "'${trynDir2}'/'${Libname}'case2contigs.out"; 
                          }else{print id >> "'${trynDir2}'/'${Libname}'caseUp3contigs.out";}}} id=$1"_"$2; cpt=1; l1=$0; fr=$3;
            }else{cpt=cpt+1; l2=$0}}
                          END{if(cpt==1){print ">'$nameDir'_"id"_"fr >> "'${nameDir}'/'${Libname}'_final.sequences.fa"; print $8 >> "'${nameDir}'/'${Libname}'_final.sequences.fa"
                          }else{if(cpt==2){print l1 >> "'${trynDir2}'/'${Libname}'case2contigs.out"; print $0 >> "'${trynDir2}'/'${Libname}'case2contigs.out"; 
                          }else{print id >> "'${trynDir2}'/'${Libname}'caseUp3contigs.out";}}}' ${trynDir2}/${Libname}_length_and_nbreads.Filtered.tsv
#Case nb contigs = 2
for case2 in `awk '{print $1","$2}' ${trynDir2}/${Libname}case2contigs.out | sort -u`; do
    well=`echo $case2 | awk -F"," '{print $2}'`
    fasta=${trynDir2}/${Libname}_${well}_filtered.fa
    makeblastdb -dbtype nucl -in ${fasta}
    blastn -db ${fasta} -query ${fasta} -word_size 7 -evalue 0.1 -outfmt 6 -out ${trynDir2}/blast_${Libname}_${well}.out
        awk '$3==100 && $4==46' ${trynDir2}/blast_${Libname}_${well}.out > ${trynDir2}/blast_${Libname}_${well}.Filtered.out
    awk '$2=="'${well}'"' ${trynDir2}/${Libname}_length_and_nbreads.Filtered.tsv | paste - - | sed 's/>//g' | awk -F "\t" '
        function revcomp(arg) {o = ""; for(i = length(arg); i > 0; i--){o = o c[substr(arg, i, 1)]} return(o)}
        BEGIN{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\t"; 
        while(getline < "'${trynDir2}'/blast_'${Libname}'_'${well}'.Filtered.out" > 0){ couple[$1]=$2; posA[$1]=$7; posB[$1]=$8; posA[$2]=$9;  posB[$2]=$10;}
        while(getline < "'${fasta}'" > 0){if($1~/^>/){id=$0;}else{seq[id]=$0}}}
        {split($7,contig1," "); split($15,contig2," ");
        if($3==$11){
            print ">'$nameDir'_'${Libname}'_'${well}'_"$3"_1" >> "'$nameDir'/'${Libname}'_final.sequences.fa";
            print seq[">"$7] >> "'$nameDir'/'${Libname}'_final.sequences.fa";
            print ">'$nameDir'_'${Libname}'_'${well}'_"$11"_2" >> "'$nameDir'/'${Libname}'_final.sequences.fa";
            print seq[">"$15] >> "'$nameDir'/'${Libname}'_final.sequences.fa";
        }else{
            if(couple[contig1[1]]==contig2[1]){
                if(posB[contig1[1]]<50 && posB[contig2[1]]<50){
                    final_seq=revcomp(seq[">"$7])""substr(seq[">"$15],47,length(seq[">"$15]));
                }else if(posB[contig1[1]]<50 && posB[contig2[1]]>50){
                    final_seq=revcomp(seq[">"$7])""substr(revcomp(seq[">"$15]),47,length(seq[">"$15])); 
                }else if(posB[contig1[1]]>50 && posB[contig2[1]]<50){
                    final_seq=seq[">"$7]""substr(seq[">"$15],47,length(seq[">"$15]));
                }else{
                    final_seq=seq[">"$7]""substr(revcomp(seq[">"$15]),47,length(seq[">"$15]));
                }
                print ">'$nameDir'_'${Libname}'_'${well}'" >> "'$nameDir'/'${Libname}'_final.sequences.fa";
                print final_seq >> "'$nameDir'/'${Libname}'_final.sequences.fa";
            }
            if(!couple[contig1[1]]){
		print ">'$nameDir'_'${Libname}'_'${well}'_"$3 >> "'$nameDir'/'${Libname}'_final.sequences.fa";
                print seq[">"$7] >> "'$nameDir'/'${Libname}'_final.sequences.fa";
                print ">'$nameDir'_'${Libname}'_'${well}'_"$11 >> "'$nameDir'/'${Libname}'_final.sequences.fa";
                print seq[">"$15] >> "'$nameDir'/'${Libname}'_final.sequences.fa";
	    }
        }}'
done
#Case nb contigs > 3
for case3 in `cat ${trynDir2}/${Libname}caseUp3contigs.out`; do
    well=`echo $case3 | awk -F"_" '{print $3}'`
    fasta=${trynDir2}/${case3}_filtered.fa
    makeblastdb -dbtype nucl -in ${fasta}
    blastn -db ${fasta} -query ${fasta} -word_size 7 -evalue 0.1 -outfmt 6 -out ${trynDir2}/blast_${Libname}_${well}.out
    awk '$3==100 && $4==46 {if(!($1==id2 && $2==id1)){split($1,fr1,"_"); split($2,fr2,"_"); if(fr1[1]!=fr2[1]) print $0; id1=$1; id2=$2}}' ${trynDir2}/blast_${Libname}_${well}.out > ${trynDir2}/blast_${Libname}_${well}.Filtered.out
    nbmatch=`wc -l ${trynDir2}/blast_${Libname}_${well}.Filtered.out | awk '{print $1}'`
    if [ $nbmatch -eq 0 ]; then
        awk '{if($1~/^>/){split($1,s,">"); print ">'$nameDir'_'${Libname}'_'${well}'_"s[2]}else{print $0}}' ${trynDir2}/${Libname}_${well}_filtered.fa >> $nameDir/${Libname}_final.sequences.fa
    else
      cat ${trynDir2}/blast_${Libname}_${well}.Filtered.out | awk -F "\t" '
        function revcomp(arg) {o = ""; for(i = length(arg); i > 0; i--){o = o c[substr(arg, i, 1)]} return(o)}
        BEGIN{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\t"; 
        while(getline < "'${trynDir2}'/'${Libname}'_'${well}'_filtered.fa" > 0){split($1,t," "); if($1~/^>/){id=t[1];}else{seq[id]=$0}}}
        {if($8<50 && $10<50){
            final_seq=revcomp(seq[">"$1])""substr(seq[">"$2],47,length(seq[">"$2]));
        }else if($8<50 && $10>50){
            final_seq=revcomp(seq[">"$1])""substr(revcomp(seq[">"$2]),47,length(seq[">"$2])); 
        }else if($8>50 && $10<50){
            final_seq=seq[">"$1]""substr(seq[">"$2],47,length(seq[">"$2]));
        }else{
            final_seq=seq[">"$1]""substr(revcomp(seq[">"$2]),47,length(seq[">"$2]));
        }
        if('$nbmatch'==1){
            print ">'$nameDir'_'${Libname}'_'${well}'" >> "'$nameDir'/'${Libname}'_final.sequences.fa";
            print final_seq >> "'$nameDir'/'${Libname}'_final.sequences.fa";
        }else{
            cpt=cpt+1;
            print ">'$nameDir'_'${Libname}'_'${well}'_"cpt >> "'$nameDir'/'${Libname}'_final.sequences.fa";
            print final_seq >> "'$nameDir'/'${Libname}'_final.sequences.fa";
        }}'
    fi
done
