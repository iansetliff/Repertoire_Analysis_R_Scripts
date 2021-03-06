library(alakazam)
library(igraph)
library(dplyr)

# Load Change-O file
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
df <- readChangeoDb(file)

# Select clone
sub_df <- subset(df, CLONE == 164)

# This example data set does not have ragged ends
# Preprocess clone without ragged end masking (default)
clone <- makeChangeoClone(sub_df, text_fields=c("SAMPLE", "ISOTYPE"), 
                          num_fields="DUPCOUNT")

# Show combined annotations
clone@data[, c("SAMPLE", "ISOTYPE", "DUPCOUNT")]

# Run PHYLIP and parse output
dnapars_exec <- "/home/ian/BioinformaticsTools/phylip-3.696/exe/dnapars"
graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)

# The graph has shared annotations for the clone
data.frame(CLONE=graph$clone,
           JUNCTION_LENGTH=graph$junc_len,
           V_GENE=graph$v_gene,
           J_GENE=graph$j_gene)

# The vertices have sequence specific annotations
data.frame(SEQUENCE_ID=V(graph)$name, 
           ISOTYPE=V(graph)$ISOTYPE,
           DUPCOUNT=V(graph)$DUPCOUNT)

# Plot graph with defaults
plot(graph)

# Modify graph and plot attributes
V(graph)$color <- "lightgrey"
V(graph)$color[V(graph)$name == "Germline"] <- "black"
V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
V(graph)$label <- V(graph)$ISOTYPE
E(graph)$label <- ""

# Remove large default margins
par(mar=c(0, 0, 0, 0) + 0.1)

# Define a tree layout with the Germline at the top
ly <- layout_as_tree(graph, root="Germline", circular=F, flip.y=T)

# Plot graph
plot(graph, layout=ly, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=50)

# Add legend
legend("topleft", c("Germline", "Inferred", "Sample"), 
       fill=c("black", "white", "grey80"), cex=0.75)