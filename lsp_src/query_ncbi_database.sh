$DATABASE=$1
$GENOME_NAME=$2
$FASTA=$3

blastp -db $DATABASE -query $FASTA -outfmt "6 std qcovs" -out sp_${GENOME_NAME}_vs_all.tab -num_threads 4
