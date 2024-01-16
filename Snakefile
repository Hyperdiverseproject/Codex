#Snakefile
configfile: "config.yaml"
import pandas as pd
import glob
import os

# Function to pair R1 and R2 files
def get_read_pairs(directory):
    r1_files = glob.glob(os.path.join(directory, "*_R1.fastq.gz"))
    read_pairs = {}
    for r1 in r1_files:
        sample = os.path.basename(r1).split("_R1.fastq.gz")[0]
        r2 = os.path.join(directory, sample + "_R2.fastq.gz")
        read_pairs[sample] = {"R1": r1, "R2": r2}
    return read_pairs
read_pairs = get_read_pairs(config['reads_directory'])

# Read the well names from the barcode file
well= pd.read_table(config["barcodes"], header=None)[0].tolist()

#Run all
rule all:
    input:
        finalseq=expand("{outdir}/Final/{sample}_final.sequences.fa", outdir=config["output_directory"], sample=read_pairs.keys()),
    
# Rule for trimming R1 and R2 using fastp and merging aka STEP1
rule trim_merge_reads:
    input:
        r1=lambda wildcards: read_pairs[wildcards.sample]["R1"],
        r2=lambda wildcards: read_pairs[wildcards.sample]["R2"],
    output:
        r1trim="{outdir}/First_step/{sample}_R1_trim.fastq.gz",
        r2trim="{outdir}/First_step/{sample}_R2_trim.fastq.gz",
        read_merged="{outdir}/First_step/{sample}_merged.fastq.gz",
        read_merged_filtered="{outdir}/First_step/{sample}_filtered.fastq",
    threads: 28
    shell:
        """
        fastp -i {input.r1} -o {output.r1trim} -I {input.r2} -O {output.r2trim} -m --merged_out {output.read_merged} -w {threads};
        zcat {output.read_merged} | paste - - - - | awk -F"\\t" 'length($2)>=300 {{print $1; print $2; print $3; print $4;}}' > {output.read_merged_filtered};
        """

rule demultiplex_well:
    input:
        #assembled="{output_directory}/First_step/{sample}_filtered.fastq",
        assembled=config['output_directory'] + "/First_step/{sample}_filtered.fastq",
        tag=config['barcodes'],
    output:
        #reads="{output_directory}/Second_step/{sample}_{well}.fq",
        reads=[config['output_directory'] + "/Second_step/{sample}_" + w + ".fq" for w in well],
        #step2="{outdir}/Second_step/{sample}_step_2.txt",
    params:
        outdir=config['output_directory']
    shell:
        """
        demDir=`echo {input.assembled} | sed 's/_filtered.fastq//; s/First_step/Second_step/'`
        cat {input.assembled} | paste - - - - | awk -v demDir=$demDir -f demultiplex_script.awk {input.tag}
        for w in {well}; do
            touch "{params.outdir}/Second_step/{wildcards.sample}_${{w}}.fq"
        done
        """

# Rule for preparing tags database using BLAST aka STEP3.1
rule prepare_tag_database:
    input:
        tag=config['barcodes'],
    output:
        amorces_seq=expand("{outdir}/Third_step/amorces.fa", outdir=config["output_directory"]),
        db=expand("{outdir}/Third_step/amorces.fa.ndb", outdir=config["output_directory"]),
    threads: 1
    shell:
        """
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
        function generate {{
            if [[ $1 =~ (.*)([RYSWiKMBDHVN])(.*) ]]; then
                local head=${{BASH_REMATCH[1]}}
                local tail=${{BASH_REMATCH[3]}}
                local -a seqs=()
                for nuc in ${{deg2nuc[${{BASH_REMATCH[2]}}]}}; do
                    seqs+=($(generate "${{head}}${{nuc}}${{tail}}"))
                done
                echo "${{seqs[@]}}"
            else
                echo "$1"
            fi
        }}
        # Generate the sequences and output in FASTA format
        i=1
        for gen in `awk '{{print $2}}' {input.tag} | sort -u`; do
            id="FR1_Fwd"
            sequence=$gen
            for seq in $(generate "$sequence"); do
                echo ">${{id}}_${{i}}" >> {output.amorces_seq}
                echo "$seq" >> {output.amorces_seq}
                ((i++))
            done
        done
        i=1
        for gen in `awk '{{print $3}}' {input.tag} | sort -u`; do
            id="FR1_Rev"
            sequence=$gen
            for seq in $(generate "$sequence"); do
                echo ">${{id}}_${{i}}" >> {output.amorces_seq}
                echo "$seq" >> {output.amorces_seq}
                ((i++))
            done
        done
        i=1
        for gen in `awk '{{print $4}}' {input.tag} | sort -u`; do
            id="FR2_Fwd"
            sequence=$gen
            for seq in $(generate "$sequence"); do
                echo ">${{id}}_${{i}}" >> {output.amorces_seq}
                echo "$seq" >> {output.amorces_seq}
                ((i++))
            done
        done
        i=1
        for gen in `awk '{{print $5}}' {input.tag} | sort -u`; do
            id="FR2_Rev"
            sequence=$gen
            for seq in $(generate "$sequence"); do
                echo ">${{id}}_${{i}}" >> {output.amorces_seq}
                echo "$seq" >> {output.amorces_seq}
                ((i++))
            done
        done
        makeblastdb -dbtype nucl -in {output.amorces_seq}
        """

# Rule for generating fasta files aka STEP 3.2
rule generate_fasta:
    input:
        fastq=config['output_directory'] + "/Second_step/{sample}_{well}.fq",
    output:
        fasta=config['output_directory'] + "/Third_step/{sample}_{well}.fa",
    threads: 28
    shell:
        """
        if [ -s {input.fastq} ]; then
            awk -F"\\t" '{{if($2>=200){{print ">"$1; print $2;}}}}' {input.fastq} > {output.fasta}
        else
            touch {output.fasta}
        fi
        """

# Rule for generating blast files aka STEP 3.3
rule generate_blast:
    input:
        fasta=rules.generate_fasta.output.fasta,
        amorces=rules.prepare_tag_database.output.amorces_seq,
    output:
        blast="{outdir}/Third_step/blast_{sample}_{well}",
    threads: 28
    shell:
        """
        if [ -s "{input.fasta}" ]; then
            blastn -db {input.amorces} -query {input.fasta} -outfmt 6 -word_size 7 -evalue 0.1 -perc_identity 95 -num_threads {threads} | awk -F"\\t" '{{OFS="\\t"; split($2,strand,"_"); \
            if(id==$1){{if(strand[2]=="Fwd" && fwd<1){{fr=fr","$2;pos=pos","$7"_"$8;fwd=fwd+1}} if(strand[2]=="Rev" && rev<1){{fr=fr","$2;pos=pos","$7"_"$8;\
            rev=rev+1}}}}else{{if(NR!=1) print id, fr, pos; if(strand[2]=="Fwd"){{fwd=1;rev=0}}else{{rev=1;fwd=0}}id=$1;fr=$2;pos=$7"_"$8;}}}}\
            END{{print id, fr, pos;}}' | grep "," > {output.blast} || true
        else
            touch {output.blast}
        fi
        """

# Rule for segregating fragments aka STEP 3.4
rule segregate_fragments :
    input:
        blast="{outdir}/Third_step/blast_{sample}_{well}",
        reads="{outdir}/Second_step/{sample}_{well}.fq",
    params:
        step3="{outdir}/Third_step/{sample}_step_3.txt",
    output:
        reads_fragments1="{outdir}/Fourth_step/{sample}_{well}_FR1.fastq",
        reads_fragments2="{outdir}/Fourth_step/{sample}_{well}_FR2.fastq",
    threads: 28
    shell:
        """
        lib_well=$(basename {input.reads} .fq)
        touch {output.reads_fragments1}
        touch {output.reads_fragments2}
        awk -F "\\t" 'function revcomp(arg) {{o = ""; for(i = length(arg); i > 0; i--){{o = o c[substr(arg, i, 1)]}} return(o)}}
            function revqual(qual){{x=""; for(i=length(qual);i!=0;i--){{x=x substr(qual,i,1)}} return(x)}}
            BEGIN{{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; frag1=0; frag2=0; nofrag=0; OFS="\\t";
            while(getline < "{input.blast}" > 0){{
            nb=split($2,fr,","); fragA[$1]=fr[1]; fragB[$1]=fr[2];
            split($3,position,","); posA[$1]=position[1]; posB[$1]=position[2]}}}}
            {{split($1,s," "); seq=s[1]; split(fragA[seq],frA,"_"); split(fragB[seq],frB,"_"); split(posA[seq],psA,"_"); split(posB[seq],psB,"_");
            if(frA[1] && frB[1]){{
                if(frA[1]==frB[1]){{
                    if(frA[1]=="FR1"){{
                        if(psA[1]<psB[1]){{
                            print $1 >> "{output.reads_fragments1}";
                            print $2 >> "{output.reads_fragments1}";
                            print $3 >> "{output.reads_fragments1}";
                            print $4 >> "{output.reads_fragments1}";
                        }}else{{
                            print $1 >> "{output.reads_fragments1}";
                            print revcomp($2) >> "{output.reads_fragments1}";
                            print $3 >> "{output.reads_fragments1}";
                            print revqual($4) >> "{output.reads_fragments1}";
                        }}
                        frag1=frag1+1
                    }}else{{
                        if(psA[1]<psB[1]){{
                            print $1 >> "{output.reads_fragments2}";
                            print $2 >> "{output.reads_fragments2}";
                            print $3 >> "{output.reads_fragments2}";
                            print $4 >> "{output.reads_fragments2}";
                        }}else{{
                            print $1 >> "{output.reads_fragments2}";
                            print revcomp($2) >> "{output.reads_fragments2}";
                            print $3 >> "{output.reads_fragments2}";
                            print revqual($4) >> "{output.reads_fragments2}";
                        }}
                        frag2=frag2+1}}
                }}else{{
                    nofrag=nofrag+1;}}
            }}else{{
                nofrag=nofrag+1;}}
        }}END{{
            print "'${{lib_well}}'" >> "{params.step3}";
            print "###########################################" >> "{params.step3}";
            print "Assigned fragments 1 = \\t "frag1 >> "{params.step3}";
            print "Assigned fragments 2 = \\t "frag2 >> "{params.step3}";
            print "Not Assigned fragments = \\t "nofrag >> "{params.step3}";
            print "###########################################" >> "{params.step3}";}}' {input.reads}
        """
        
# Rule for assembling fragments aka STEP 4.1
rule trinity_run:
    input:
        fragmentsFR1="{outdir}/Fourth_step/{library}_{well}_FR1.fastq",
        fragmentsFR2="{outdir}/Fourth_step/{library}_{well}_FR2.fastq",
    output:
        trinity_FR1="{outdir}/Fourth_step/trinity_{library}_{well}_FR1.Trinity.fasta",
        trinity_FR2="{outdir}/Fourth_step/trinity_{library}_{well}_FR2.Trinity.fasta",
    threads: 28
    shell:
        """
        touch {output.trinity_FR1}
        touch {output.trinity_FR2}
        triout1="$(dirname {output.trinity_FR1})/$(basename {output.trinity_FR1} .Trinity.fasta)"
        triout2="$(dirname {output.trinity_FR2})/$(basename {output.trinity_FR2} .Trinity.fasta)"
        echo $triout1 $triout2
        if [ -s {input.fragmentsFR1} ]; then
            Trinity --seqType fq --max_memory 30G --CPU {threads} --single {input.fragmentsFR1} --output $triout1 --full_cleanup || true
        fi
        if [ -s {input.fragmentsFR2} ]; then
            Trinity --seqType fq --max_memory 30G --CPU {threads} --single {input.fragmentsFR2} --output $triout2 --full_cleanup || true
        fi
        """



# Rule for read depth aka STEP 4.2
rule minimap:
    input:
        fragmentsFR1="{outdir}/Fourth_step/{library}_{well}_FR1.fastq",
        fragmentsFR2="{outdir}/Fourth_step/{library}_{well}_FR2.fastq",
        trinity_FR1="{outdir}/Fourth_step/trinity_{library}_{well}_FR1.Trinity.fasta",
        trinity_FR2="{outdir}/Fourth_step/trinity_{library}_{well}_FR2.Trinity.fasta",
    output:
        minimap_FR1="{outdir}/Fourth_step/trinity_{library}_{well}_FR1.nbreads",
        minimap_FR2="{outdir}/Fourth_step/trinity_{library}_{well}_FR2.nbreads",
    threads: 28
    shell:
        """
        temp="$(dirname {input.fragmentsFR1})/$(basename {input.fragmentsFR1} .fastq).fa"
        cat {input.fragmentsFR1} | paste - - - - | awk -F"\\t" '{{print $1; print $2}}' | sed 's/@/>/' > $temp
        nbreads=`wc -l {input.fragmentsFR1} | awk '{{print $1/4}}'`
        if [ -s {input.fragmentsFR1} ] && [ -s $temp ]; then
            minimap2 -t {threads} -ax sr {input.trinity_FR1} $temp | grep -v "^@" | awk '$2!=4 {{print $3}}' | sort | uniq -c | sort -k1,1nr | awk '{{OFS="\\t"; print $2,$1,($1/'$nbreads')*100}}' > {output.minimap_FR1}
        else
            touch {output.minimap_FR1}
        fi
        temp="$(dirname {input.fragmentsFR2})/$(basename {input.fragmentsFR2} .fastq).fa"
        cat {input.fragmentsFR2} | paste - - - - | awk -F"\\t" '{{print $1; print $2}}' | sed 's/@/>/' > $temp
        nbreads=`wc -l {input.fragmentsFR2} | awk '{{print $1/4}}'`
        if [ -s {input.fragmentsFR2} ] && [ -s $temp ]; then
            minimap2 -t {threads} -ax sr {input.trinity_FR2} $temp | grep -v "^@" | awk '$2!=4 {{print $3}}' | sort | uniq -c | sort -k1,1nr | awk '{{OFS="\\t"; print $2,$1,($1/'$nbreads')*100}}' > {output.minimap_FR2}
        else
            touch {output.minimap_FR2}
        fi
        """

# Rule for match on NT aka STEP 4.3
rule blastOnDB:
    input:
        trinity_FR1="{outdir}/Fourth_step/trinity_{library}_{well}_FR1.Trinity.fasta",
        trinity_FR2="{outdir}/Fourth_step/trinity_{library}_{well}_FR2.Trinity.fasta",
    params:
        ncbintdb=config["ncbidb"]
    output:
        blast_FR1="{outdir}/Fourth_step/trinity_{library}_{well}_FR1.nt",
        blast_FR2="{outdir}/Fourth_step/trinity_{library}_{well}_FR2.nt",
    threads: 28
    run:
        if params.ncbintdb:
            shell(
                """
                if [ -s {input.trinity_FR1} ]; then
                    blastn -db {params.ncbintdb} -query {input.trinity_FR1} -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -out {output.blast_FR1} -num_threads {threads} -num_alignments 1
                else
                    touch {output.blast_FR1}
                fi
                if [ -s {input.trinity_FR2} ]; then
                    blastn -db {params.ncbintdb} -query {input.trinity_FR2} -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -out {output.blast_FR2} -num_threads {threads} -num_alignments 1
                else
                    touch {output.blast_FR2}
                fi
                """
            )
        else:
            shell(
                """
                echo 'ncbintdb not provided. Skipping blast on the nt databe.'
                touch {output.blast_FR1} {output.blast_FR2}
                """
            )

# Rule for removing tags aka STEP 5.1
rule remove_tag:
    input:
        frag1="{outdir}/Fourth_step/trinity_{library}_{well}_FR1.Trinity.fasta",
        frag2="{outdir}/Fourth_step/trinity_{library}_{well}_FR2.Trinity.fasta",
    params:
        amorce_seq=rules.prepare_tag_database.output.amorces_seq,
    output:
        blast_FR1="{outdir}/Fifth_step/trinity_{library}_{well}_FR1.blastout",
        blast_FR2="{outdir}/Fifth_step/trinity_{library}_{well}_FR2.blastout",
        notag_FR1="{outdir}/Fifth_step/trinity_{library}_{well}_FR1.notag.fa",
        notag_FR2="{outdir}/Fifth_step/trinity_{library}_{well}_FR2.notag.fa",
    threads: 28
    shell:
        """
        if [ -s {input.frag1} ]; then
            blastn -db {params.amorce_seq} -query {input.frag1} -outfmt 6 -word_size 7 -evalue 0.1 -perc_identity 95 | awk '{{if(id==$1){{if(deb!=1){{if($8<50){{deb=1; print $0}}}} \
            if(fin!=1){{if($8>50){{fin=1; print $0}}}}}}else{{id=$1; deb=0; fin=0; if(id==$1){{if(deb!=1){{\
            if($8<50){{deb=1; print $0}}}} if(fin!=1){{if($8>50){{fin=1; print $0}}}}}}}}}}' | awk '{{if(id==$1){{cpt=cpt+1; l2=$0}}else{{if(cpt==2){{print l1; print l2;}} cpt=1; l1=$0; id=$1}}}}END{{if(cpt==2){{print l1; print l2;}}}}' \
            > {output.blast_FR1}
            cat {input.frag1} | awk '{{OFS="\\t"; if($0~/^>/){{print id, seq; id=$1; seq="";}}else{{seq=seq""$0}}}}END{{print id,seq}}' | sed '1d' | sed 's/>//' | awk -F "\\t" 'function revcomp(arg) {{o = ""; for(i = length(arg); i > 0; i--){{o = o c[substr(arg, i, 1)]}} return(o)}}
            BEGIN{{OFS="\\t";c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ;\
            while(getline < "{output.blast_FR1}" > 0){{\
            if($8<50){{deb[$1]=$8}}else{{fin[$1]=$7}} split($2,t,"_"); if(($10>$9 && t[2]=="Rev")||($9>$10 && t[2]=="Fwd")){{rev[$1]="yes"}}else{{rev[$1]="no"}}}}}}
            {{split($1,name," ");\
            if(deb[name[1]] && fin[name[1]]){{\
                print ">FR1_"$1;\
            if(rev[$1]=="no"){{\
                print substr($2,deb[name[1]]+1,length($2)-(deb[name[1]]+(length($2)-fin[name[1]]+1)))\
            }}else{{\
                print revcomp(substr($2,deb[name[1]]+1,length($2)-(deb[name[1]]+(length($2)-fin[name[1]]+1))))\
            }}
            }}}}' >> {output.notag_FR1}
        else
            touch {output.blast_FR1} {output.notag_FR1}
        fi
        if [ -s {input.frag2} ]; then
            blastn -db {params.amorce_seq} -query {input.frag2} -outfmt 6 -word_size 7 -evalue 0.1 -perc_identity 95 | \
            awk '{{if(id==$1){{if(deb!=1){{if($8<50){{deb=1; print $0}}}} \
            if(fin!=1){{if($8>50){{fin=1; print $0}}}}}}else{{id=$1; deb=0; fin=0; if(id==$1){{if(deb!=1){{\
            if($8<50){{deb=1; print $0}}}} if(fin!=1){{if($8>50){{fin=1; print $0}}}}}}}}}}' | awk '{{if(id==$1){{cpt=cpt+1; l2=$0}}else{{if(cpt==2){{print l1; print l2;}} cpt=1; l1=$0; id=$1}}}}END{{if(cpt==2){{print l1; print l2;}}}}' \
            > {output.blast_FR2}
            cat {input.frag2} | awk '{{OFS="\\t"; if($0~/^>/){{print id, seq; id=$1; seq="";}}else{{seq=seq""$0}}}}END{{print id,seq}}' | sed '1d' | sed 's/>//' | awk -F "\\t" 'function revcomp(arg) {{o = ""; for(i = length(arg); i > 0; i--){{o = o c[substr(arg, i, 1)]}} return(o)}}
            BEGIN{{OFS="\\t";c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ;\
            while(getline < "{output.blast_FR2}" > 0){{\
            if($8<50){{deb[$1]=$8}}else{{fin[$1]=$7}} split($2,t,"_"); if(($10>$9 && t[2]=="Rev")||($9>$10 && t[2]=="Fwd")){{rev[$1]="yes"}}else{{rev[$1]="no"}}}}}}
            {{split($1,name," ");\
            if(deb[name[1]] && fin[name[1]]){{\
                print ">FR2_"$1;\
            if(rev[$1]=="no"){{\
                print substr($2,deb[name[1]]+1,length($2)-(deb[name[1]]+(length($2)-fin[name[1]]+1)))\
            }}else{{\
                print revcomp(substr($2,deb[name[1]]+1,length($2)-(deb[name[1]]+(length($2)-fin[name[1]]+1))))\
            }}
            }}}}' >> {output.notag_FR2}
        else
            touch {output.blast_FR2} {output.notag_FR2}
        fi
        """

rule get_tsv_files:
    input:
        #notag_FR1=config['output_directory'] + "/Fifth_step/trinity_{library}_{well}_FR1.notag.fa",
        notag_FR1=[config['output_directory'] + "/Fifth_step/trinity_{library}_" + w + "_FR1.notag.fa" for w in well],
        #notag_FR2=config['output_directory'] + "Fifth_step/trinity_{library}_{well}_FR2.notag.fa",
        notag_FR2=[config['output_directory'] + "/Fifth_step/trinity_{library}_" + w + "_FR2.notag.fa" for w in well],
        #nt_FR1=config['output_directory'] + "Fourth_step/trinity_{library}_{well}_FR1.nt",
        nt_FR1=[config['output_directory'] + "/Fourth_step/trinity_{library}_" + w + "_FR1.nt" for w in well],
        #nt_FR2=config['output_directory'] + "Fourth_step/trinity_{library}_{well}_FR2.nt",
        nt_FR2=[config['output_directory'] + "/Fourth_step/trinity_{library}_" + w + "_FR2.nt" for w in well],
        minimap_FR1=[config['output_directory'] + "/Fourth_step/trinity_{library}_" + w + "_FR1.nbreads" for w in well],
        minimap_FR2=[config['output_directory'] + "/Fourth_step/trinity_{library}_" + w + "_FR2.nbreads" for w in well],
    params:
        minfr1=config["minfr1"],
        maxfr1=config["maxfr1"],
        minfr2=config["minfr2"],
        maxfr2=config["maxfr2"],
        minreadsContig=config["minreadsContig"],
        mincov=config["mincov"],
        dir=config["output_directory"],
    output:
        #fasta_filtered=config['output_directory'] + "Fifth_step/{library}_{well}_filtered.fa",
        tsv=config['output_directory'] + "/Fifth_step/{library}_length_and_nbreads.tsv",
        tsv_filtered=config['output_directory'] + "/Fifth_step/{library}_length_and_nbreads.Filtered.tsv",
    threads: 1
    shell:
        """
        touch {output.tsv_filtered} {output.tsv}
        for notagFR1 in {input.notag_FR1}; 
        do well=`echo $notagFR1 | awk -F"_" '{{print $5}}'`;
            lib=`echo $notagFR1 | awk -F"_" '{{print $3"_"$4}}'`;
            minimapfr1="{params.dir}/Fourth_step/"$(basename $notagFR1 .notag.fa)".nbreads"
            fastafile="{params.dir}/Fifth_step/${{lib}}_${{well}}_filtered.fa"
            touch $fastafile
            if [ -s $notagFR1 ]; then
                cat $notagFR1 | paste - - | awk -F"\\t" 'BEGIN{{\
                while(getline < "'$minimapfr1'" > 0){{nbreads[">FR1_"$1]=$2; perreads[">FR1_"$1]=$3;}}}}\
                {{OFS="\\t"; split($1,t," "); seq=$2; l=length(seq); nreads=nbreads[t[1]]; preads=perreads[t[1]]; ntmatch=nt[t[1]];\
                print "'$lib'","'$well'","FR1",l,nreads,preads,ntmatch,$1,seq >> "{output.tsv}";\
                if(l>{params.minfr1} && l<{params.maxfr1} && nreads>{params.minreadsContig} && preads>{params.mincov}){{\
                print "'$lib'","'$well'","FR1",l,nreads,preads,ntmatch,$1,seq >> "{output.tsv_filtered}"; \
                print $1; print seq}}}}' >> $fastafile
            fi
        done
        for notagFR2 in {input.notag_FR2};
        do well=`echo $notagFR2 | awk -F"_" '{{print $5}}'`;
            lib=`echo $notagFR2 | awk -F"_" '{{print $3"_"$4}}'`;
            minimapfr2="{params.dir}/Fourth_step/"$(basename $notagFR2 .notag.fa)".nbreads"
            fastafile="{params.dir}/Fifth_step/${{lib}}_${{well}}_filtered.fa"
            touch $fastafile
            if [ -s $notagFR2 ]; then
                cat $notagFR2 | paste - - | awk -F"\\t" 'BEGIN{{\
                while(getline < "'$minimapfr2'" > 0){{nbreads[">FR2_"$1]=$2; perreads[">FR2_"$1]=$3;}}}}\
                {{OFS="\\t"; split($1,t," "); seq=$2; l=length(seq); nreads=nbreads[t[1]]; preads=perreads[t[1]]; ntmatch=nt[t[1]];\
                print "'$lib'","'$well'","FR2",l,nreads,preads,ntmatch,$1,seq >> "{output.tsv}";\
                if(l>{params.minfr2} && l<{params.maxfr2} && nreads>{params.minreadsContig} && preads>{params.mincov}){{\
                print "'$lib'","'$well'","FR2",l,nreads,preads,ntmatch,$1,seq >> "{output.tsv_filtered}"; \
                print $1; print seq}}}}' >> $fastafile
            fi
        done
        """

rule merging_fragments:
    input:
        tsv_filtered="{outdir}/Fifth_step/{library}_length_and_nbreads.Filtered.tsv",
    params:
        up3contigs="{outdir}/Fifth_step/{library}caseUp3contigs.out",
        eq2contigs="{outdir}/Fifth_step/{library}case2contigs.out",
        dir=config["output_directory"],
    output:
        finalseq="{outdir}/Final/{library}_final.sequences.fa",
    shell:
        """
        mkdir -p {params.dir}/Final/
        lib=`echo $(basename {input.tsv_filtered}) | awk -F"_" '{{print $1"_"$2}}'`;
        sort -k1,1 -k2,2 {input.tsv_filtered} > {params.dir}/Final/${{lib}}_length_and_nbreads.Filtered.tsv
        cat {input.tsv_filtered} | sort -k1,1 -k2,2 | awk -F"\\t" '{{\
            if(id!=$1"_"$2){{\
                if(NR!=1){{\
                    if(cpt==1){{\
                        print ">{params.dir}_"id"_"fr >> "{output.finalseq}"; \
                        print $9 >> "{output.finalseq}"\
                    }}else{{\
                        if(cpt==2){{\
                            print l1 >> "{params.eq2contigs}"; \
                            print l2 >> "{params.eq2contigs}"; \
                        }}else{{\
                            print id >> "{params.up3contigs}";\
                        }}\
                    }}\
                }}\
                id=$1"_"$2; cpt=1; l1=$0; fr=$3;\
            }}else{{\
                cpt=cpt+1; l2=$0\
            }}
        }}END{{\
            if(cpt==1){{\
                print ">{params.dir}_"id"_"fr >> "{output.finalseq}"; \
                print $9 >> "{output.finalseq}"\
            }}else{{\
                if(cpt==2){{\
                    print l1 >> "{params.eq2contigs}"; \
                    print $0 >> "{params.eq2contigs}"; \
                }}else{{\
                    print id >> "{params.up3contigs}";\
                }}\
            }}\
        }}'
        #Case nb contigs = 2
        for case2 in `awk '{{print $1","$2}}' {params.eq2contigs} | sort -u`; do
            well=`echo $case2 | awk -F"," '{{print $2}}'`
            fasta="{params.dir}/Fifth_step/${{lib}}_${{well}}_filtered.fa"
            blast="{params.dir}/Fifth_step/blast_${{lib}}_${{well}}.out"
            blast_filtered="{params.dir}/Fifth_step/blast_${{lib}}_${{well}}.filtered.out"
            makeblastdb -dbtype nucl -in ${{fasta}}
            blastn -db ${{fasta}} -query ${{fasta}} -word_size 7 -evalue 0.1 -outfmt 6 -out $blast
            awk '$3==100 && $4==46' $blast > $blast_filtered
            awk '$2=="'${{well}}'"' {input.tsv_filtered} | paste - - | sed 's/>//g' | awk -F "\\t" '\
            function revcomp(arg) {{o = ""; for(i = length(arg); i > 0; i--){{o = o c[substr(arg, i, 1)]}} return(o)}}\
            BEGIN{{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\\t"; \
            while(getline < "'$blast_filtered'" > 0){{couple[$1]=$2; posA[$1]=$7; posB[$1]=$8; posA[$2]=$9;  posB[$2]=$10;}}\
            while(getline < "'${{fasta}}'" > 0){{if($1~/^>/){{id=$0;}}else{{seq[id]=$0}}}}}}\
            {{split($7,contig1," "); split($15,contig2," ");\
            if($3==$11){{\
                print ">{params.dir}_'${{lib}}'_'${{well}}'_"$3"_1" >> "{output.finalseq}";\
                print seq[">"$7] >> "{output.finalseq}";\
                print ">{params.dir}_'${{lib}}'_'${{well}}'_"$11"_2" >> "{output.finalseq}";\
                print seq[">"$15] >> "{output.finalseq}";\
            }}else{{\
                if(couple[contig1[1]]==contig2[1]){{\
                    if(posB[contig1[1]]<50 && posB[contig2[1]]<50){{\
                        final_seq=revcomp(seq[">"$7])""substr(seq[">"$15],47,length(seq[">"$15]));\
                    }}\
                }}else if(posB[contig1[1]]<50 && posB[contig2[1]]>50){{\
                    final_seq=revcomp(seq[">"$7])""substr(revcomp(seq[">"$15]),47,length(seq[">"$15])); \
                }}else if(posB[contig1[1]]>50 && posB[contig2[1]]<50){{\
                    final_seq=seq[">"$8]""substr(seq[">"$17],47,length(seq[">"$15]));\
                }}else{{\
                    final_seq=seq[">"$8]""substr(revcomp(seq[">"$17]),47,length(seq[">"$15]));\
                }}\
                print ">{params.dir}_'${{lib}}'_'${{well}}'" >> "{output.finalseq}";\
                print final_seq >> "{output.finalseq}";\
            }}if(!couple[contig1[1]]){{\
                print ">{params.dir}_'${{lib}}'_'${{well}}'_"$3 >> "{output.finalseq}";\
                print seq[">"$8] >> "{output.finalseq}";\
                print ">{params.dir}_'${{lib}}'_'${{well}}'_"$12 >> "{output.finalseq}";\
                print seq[">"$17] >> "{output.finalseq}";\
            }}}}'
        done
        #Case nb contigs > 3
        for case3 in `cat {params.up3contigs}`; do
            well=`echo $case3 | awk -F"_" '{{print $3}}'`
            fasta="{params.dir}/Fifth_step/${{lib}}_${{well}}_filtered.fa"
            blast="{params.dir}/Fifth_step/blast_${{lib}}_${{well}}.out"
            blast_filtered="{params.dir}/Fifth_step/blast_${{lib}}_${{well}}.filtered.out"
            makeblastdb -dbtype nucl -in $fasta
            blastn -db $fasta -query $fasta -word_size 7 -evalue 0.1 -outfmt 6 -out $blast
            awk '$3==100 && $4==46 {{if(!($1==id2 && $2==id1)){{split($1,fr1,"_"); split($2,fr2,"_"); if(fr1[1]!=fr2[1]) print $0; id1=$1; id2=$2}}}}' $blast > $blast_filtered
            nbmatch=`wc -l $blast_filtered | awk '{{print $1}}'`
            if [ $nbmatch -eq 0 ]; then
                awk '{{if($1~/^>/){{split($1,s,">"); print ">{params.dir}_'${{lib}}'_'${{well}}'_"s[2]}}else{{print $0}}}}' $fasta >> {output.finalseq}
            else
                cat $blast_filtered | awk -F "\\t" '\
                function revcomp(arg) {{o = ""; for(i = length(arg); i > 0; i--){{o = o c[substr(arg, i, 1)]}} return(o)}}\
                BEGIN{{c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" ; OFS="\\t"; \
                while(getline < "'$fasta'" > 0){{split($1,t," "); if($1~/^>/){{id=t[1];}}else{{seq[id]=$0}}}}}}\
                {{if($8<50 && $10<50){{\
                    final_seq=revcomp(seq[">"$1])""substr(seq[">"$2],47,length(seq[">"$2]));\
                }}else if($8<50 && $10>50){{\
                    final_seq=revcomp(seq[">"$1])""substr(revcomp(seq[">"$2]),47,length(seq[">"$2])); \
                }}else if($8>50 && $10<50){{\
                    final_seq=seq[">"$1]""substr(seq[">"$2],47,length(seq[">"$2]));\
                }}else{{\
                    final_seq=seq[">"$1]""substr(revcomp(seq[">"$2]),47,length(seq[">"$2]));\
                }}\
                if('$nbmatch'==1){{\
                    print ">{params.dir}_'${{lib}}'_'${{well}}'" >> "{output.finalseq}";\
                    print final_seq >> "{output.finalseq}";\
                }}else{{\
                    cpt=cpt+1;\
                    print ">{params.dir}_'${{lib}}'_'${{well}}'_"cpt >> "{output.finalseq}";\
                    print final_seq >> "{output.finalseq}";\
                }}}}'
            fi
        done
        touch {output.finalseq}
        """
