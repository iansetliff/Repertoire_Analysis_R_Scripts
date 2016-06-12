# Set example data
library(shazam)
db <- InfluenzaDb

# Create a consensus sequence for each clone to avoid over-counting mutations
# Adds the CLONAL_CONSENSUS_SEQUENCE column
db_cons <- collapseByClone(db)

# Create substitution model using silent mutations
sub_matrix <- createSubstitutionMatrix(db_cons, 
                                       sequenceColumn="CLONAL_SEQUENCE", 
                                       model="S")
# Create mutability model using silent mutations
mut_matrix <- createMutabilityMatrix(db_cons, sub_matrix, 
                                     sequenceColumn="CLONAL_SEQUENCE", 
                                     model="S")

# Extend models to include ambiguous 5-mers
sub_matrix <- extendSubstitutionMatrix(sub_matrix)
mut_matrix <- extendMutabilityMatrix(mut_matrix)

# Create targeting model matrix from substitution and mutability matrices
tar_matrix <- createTargetingMatrix(sub_matrix, mut_matrix)

# Create targeting model in one step using only silent mutations
model <- createTargetingModel(db, model="S")

# Generate hedgehog plot of mutability model
plotMutability(model, nucleotides="A", style="hedgehog")

plotMutability(model, nucleotides="C", style="hedgehog")

# Generate bar plot of mutability model
plotMutability(model, nucleotides="G", style="bar")

plotMutability(model, nucleotides="T", style="bar")

# Calculate distance matrix
dist <- calcTargetingDistance(model)