# demultiplex_script.awk

BEGIN {
    FS="\t";
    OFS="\t";
    # Initialize arrays
    split("", fwd);
    split("", rev);
    split("", well);
    
    # Load barcodes and well assignments
    while (getline < ARGV[1] > 0) {
        fwd[substr($2,1,5)] = 1;
        fwd[substr($4,1,5)] = 1;
        rev[substr($3,1,5)] = 1;
        rev[substr($5,1,5)] = 1;
        well[substr($2,1,5)"_"substr($3,1,5)] = $1;
        well[substr($3,1,5)"_"substr($2,1,5)] = $1;
        well[substr($4,1,5)"_"substr($5,1,5)] = $1;
        well[substr($5,1,5)"_"substr($4,1,5)] = $1;
    }
    # remove the barcodes file from ARGV so it's not processed as input
    delete ARGV[1];
}

{
    # Process each sequence
    seqF = substr($2, 1, 5);
    seqR = revcomp(substr($2, length($2) - 5 + 1, 5));
    if (fwd[seqF] && rev[seqR]) {
        print $0 >> (demDir "_" well[seqF "_" seqR] ".fq");
        FR += 1;
    } else if (fwd[seqR] && rev[seqF]) {
        print $0 >> (demDir "_" well[seqR "_" seqF] ".fq");
        RF += 1;
    } else if ((fwd[seqF] && ! rev[seqR]) || (fwd[seqR] && ! rev[seqF])) {
        print $1 >> (demDir "_nobarcode.fastq");
        print $2 >> (demDir "_nobarcode.fastq");
        print $3 >> (demDir "_nobarcode.fastq");
        print $4 >> (demDir "_nobarcode.fastq");
        FnoR += 1;
    } else if ((! fwd[seqF] && rev[seqR]) || (! fwd[seqR] && rev[seqF])) {
        print $1 >> (demDir "_nobarcode.fastq");
        print $2 >> (demDir "_nobarcode.fastq");
        print $3 >> (demDir "_nobarcode.fastq");
        print $4 >> (demDir "_nobarcode.fastq");
        RnoF += 1;
    } else if (! fwd[seqF] && ! rev[seqR]) {
        print $1 >> (demDir "_nobarcode.fastq");
        print $2 >> (demDir "_nobarcode.fastq");
        print $3 >> (demDir "_nobarcode.fastq");
        print $4 >> (demDir "_nobarcode.fastq");
        noFR += 1;
    } else if (! fwd[seqR] && ! rev[seqF]) {
        print $1 >> (demDir "_nobarcode.fastq");
        print $2 >> (demDir "_nobarcode.fastq");
        print $3 >> (demDir "_nobarcode.fastq");
        print $4 >> (demDir "_nobarcode.fastq");
        noRF += 1;
    } else {
        print $1 >> (demDir "_nobarcode.fastq");
        print $2 >> (demDir "_nobarcode.fastq");
        print $3 >> (demDir "_nobarcode.fastq");
        print $4 >> (demDir "_nobarcode.fastq");
        noFnoRnoRnoF += 1;
    }
}

function revcomp(seq) {
    # Function to reverse complement a DNA sequence
    c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A";
    revseq = ""; 
    for(i = length(seq); i > 0; i--){
        revseq = revseq c[substr(seq, i, 1)]
    } 
    return(revseq)
}

END {
    # Final statistics or cleanup actions
    print "###########################################" >> (demDir "_step_2.txt");
    print "FR = \t\t\t "FR >> (demDir "_step_2.txt");
    print "RF = \t\t\t "RF >> (demDir "_step_2.txt");
    print "FnoR = \t\t\t "FnoR >> (demDir "_step_2.txt");
    print "RnoF = \t\t\t "RnoF >> (demDir "_step_2.txt");
    print "noFR = \t\t\t "noFR >> (demDir "_step_2.txt");
    print "noRF = \t\t\t "noRF >> (demDir "_step_2.txt");
    print "noFnoRnoRnoF = \t\t "noFnoRnoRnoF >> (demDir "_step_2.txt");
    print "totalOutcome = \t\t "FR+RF+FnoR+RnoF+noFR+noRF+noFnoRnoRnoF >> (demDir "_step_2.txt");
    print "###########################################" >> (demDir "_step_2.txt");
    print "total paires found = \t "FR+RF >> (demDir "_step_2.txt");
    print "total no paire found = \t "FnoR+RnoF+noFR+noRF+noFnoRnoRnoF >> (demDir "_step_2.txt");
    print "total seq = \t\t "NR >> (demDir "_step_2.txt");
    print "###########################################" >> (demDir "_step_2.txt");
}
