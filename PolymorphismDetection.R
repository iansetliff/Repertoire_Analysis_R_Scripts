library(tigger)
library(dplyr)

# Load example Rep-Seq data and example germline database
data(sample_db, germline_ighv)

# Detect novel alleles
novel_df = findNovelAlleles(sample_db, germline_ighv)
# Extract and view the rows that contain successful novel allele calls
novel = selectNovel(novel_df)
glimpse(novel)

# Plot evidence of the first (and only) novel allele from the example data
plotNovel(sample_db, novel[1,])

# Infer the individual's genotype, using only unmutated sequences and checking
# for the use of the novel alleles inferred in the earlier step.
geno = inferGenotype(sample_db, find_unmutated = TRUE,
                     germline_db = germline_ighv, novel_df = novel_df)

# Save the genotype sequences to a vector
genotype_seqs = genotypeFasta(geno, germline_ighv, novel_df)

# Visualize the genotype and sequence counts
print(geno)

# Make a colorful visualization. Bars indicate presence, not proportion.
plotGenotype(geno, text_size = 10)

# Use the personlized genotype to determine corrected allele assignments
V_CALL_GENOTYPED = reassignAlleles(sample_db, genotype_seqs)

# Append the corrected calls to the original data.frame
sample_db = bind_cols(sample_db, V_CALL_GENOTYPED)

# Find the set of alleles in the original calls that were not in the genotype
not_in_genotype = sample_db$V_CALL %>%
  strsplit(",") %>%
  unlist() %>%
  unique() %>%
  setdiff(names(genotype_seqs))

# Determine the fraction of calls that were ambigious before/after correction
# and the fraction that contained original calls to non-genotype alleles. Note
# that by design, only genotype alleles are allowed in "after" calls.
data.frame(
  Ambiguous = c(mean(grepl(",",sample_db$V_CALL)),
                mean(grepl(",",sample_db$V_CALL_GENOTYPED))),
  NotInGenotype = c(mean(sample_db$V_CALL %in% not_in_genotype),
                    mean(sample_db$V_CALL_GENOTYPED %in% not_in_genotype)),
  row.names = c("Before", "After")
) %>% t() %>% round(3)