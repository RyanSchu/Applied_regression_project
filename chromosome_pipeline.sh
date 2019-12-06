chromosomenumber=$1
genenames="${@:2}"

cd chr$chromosomenumber
# This function takes in a chromosome number and a gene to produce a smaller file of target values the gene

function creategenetargetexpressionfiles() {
chromosomenumber=$1
genename=$2
echo "Creating expression file for $genename"
head -1 chr$chromosomenumber.gene_expression.txt > ${genename}target.txt
grep $genename chr$chromosomenumber.gene_expression.txt >> ${genename}target.txt
echo "Expression target for genes created for chromosome $chromosomenumber for the gene $genename"
}

# Create a file that contains all of the SNP predictors(genotypes) that are less than the FDR threshold
function create_snp_predictor_file_given_threshold_for_gene(){
chromosomenumber=$1
shift
threshold=$1
genename=$2
head -1 chr$chromosomenumber.genoype.txt > ${genename}_snps_${threshold}.txt
echo "Getting SNPS for $genename with a threshold setting of $threshold"
awk '{if ($5 < 0.05) print $1}' chr$chromosomenumber.$genename.snps_MEQTL_PRUNED_$threshold | xargs -I{} grep {} chr$chromosomenumber.genoype.txt >> ${genename}_snps_${threshold}.txt
}


# Chromosome Pipeline is here

for genename in $genenames
do
    creategenetargetexpressionfiles $chromosomenumber $genename
    for threshold in 0.3 0.5 0.8
    do
        create_snp_predictor_file_given_threshold_for_gene $chromosomenumber $threshold $genename
        python3 ../create_data_for_r.py "${genename}target.txt" "${genename}_snps_${threshold}.txt" "${genename}_for_r_${threshold}.txt"
    done
done
cd ..
