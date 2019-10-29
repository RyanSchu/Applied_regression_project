For performing the anova on our data we will be using the MatrixEQTL package, the webpage of which can be found [here](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/).


# Pulling gene data for input
Using the gencode gene annotations found at https://www.gencodegenes.org/human/release_32lift37.html
```{bash}
zcat gencode.v32lift37.annotation.gtf.gz | tail -n +6 | awk '{gsub("\"","",$10);gsub(";","",$10); if ($12 == "\"protein_coding\";" && $3 == "gene") print $10 "\t" $1 "\t" $4 "\t" $5}' > protein_coding_list.txt
```
Next subset the data to contain only those genes present in our expression data

```{r}
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
observed_list<-read.table("Z:/data/GEU/expression/obsEXP_GEU_from_DGNdb.txt.gz",header = T)
theoretical_list<-fread("Z:/data/GEU/expression/protein_coding_list.txt",header = F)
#Sanitize gene names
observed_list$gene<-gsub("\\.[0-9]+","",observed_list$gene)
theoretical_list$V1<-gsub("\\.[0-9]+$","",theoretical_list$V1)
theoretical_list$V1<-gsub("\\.[0-9]+_[0-9]+$","",theoretical_list$V1)
observed_locations<-theoretical_list %>% filter(V1 %in% observed_list$gene)
colnames(observed_locations)<-c('gene','chr','start','stop')
for (i in c(1:22)){
  chrom<-"chr" %&% i
  tmp<-filter(observed_locations, chr == chrom)
  fwrite(tmp,"Z:/data/GEU/MEQTL/gene_locations/" %&% chrom %&% ".protein_coding_list")
  tmp2<-filter(observed_list, gene %in% observed_locations$gene)
  fwrite(tmp2,"Z:/data/GEU/MEQTL/ALL/" %&% chrom %&% ".gene_expression.txt")
}
```



# Pulling SNP locations
At this point we've already parsed our snps into dosage format and pulling these are fairly trivial
```{bash}
for i in {1..22};
do 
  zcat chr${i}.1KGP_MAF0.01.dosages.txt.gz | awk '{print $2 "\t" $1 "\t" $3}' > /home/rschubert1/data/GEU/MEQTL/snp_locations/chr.${i}.snp.loc.txt; 
done
```

To split this data into genotypes as well as also fairly trivial 
```{bash}
#Note transpose funciton taken from stackoverflow
transpose() {
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' "$1"
}

for i in {1..22}; 
do  
  zcat chr${i}.1KGP_MAF0.01.dosages.txt.gz | awk '{print $2}' | transpose > /home/rschubert1/data/GEU/MEQTL/ALL/chr${i}.genotype.txt
  cut -f 7- -d$'\t' >> /home/rschubert1/data/GEU/MEQTL/ALL/chr${i}.genotype.txt
done 
```

