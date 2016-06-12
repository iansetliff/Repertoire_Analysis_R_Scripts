library(alakazam)

# Load Change-O file
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
db <- readChangeoDb(file)
db <- db[db$SAMPLE == "RL01", ]

db_props <- aminoAcidProperties(db, seq="JUNCTION", nt=TRUE, trim=TRUE, 
                                label="CDR3")

# The full set of properties are calculated by default
dplyr::select(db_props[1:3, ], starts_with("CDR3"))

# Plot a subset of the properties
tmp_theme <- theme_bw() + theme(legend.position="bottom")
g1 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_LENGTH)) + tmp_theme +
  ggtitle("CDR3 length") + 
  xlab("Isotype") + ylab("Amino acids") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
g2 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_GRAVY)) + tmp_theme + 
  ggtitle("CDR3 hydrophobicity") + 
  xlab("Isotype") + ylab("GRAVY") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
g3 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_BASIC)) + tmp_theme +
  ggtitle("CDR3 basic residues") + 
  xlab("Isotype") + ylab("Basic residues") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
g4 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_ACIDIC)) + tmp_theme +
  ggtitle("CDR3 acidic residues") + 
  xlab("Isotype") + ylab("Acidic residues") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
g5 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_BULK)) + tmp_theme +
  ggtitle("CDR3 average bulkiness") + 
  xlab("Isotype") + ylab("Bulkiness") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
g6 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_ALIPHATIC)) + tmp_theme +
  ggtitle("CDR3 normalized aliphatic index") + 
  xlab("Isotype") + ylab("Aliphatic Index") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
g7 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_POLARITY)) + tmp_theme +
  ggtitle("CDR3 average polarity") + 
  xlab("Isotype") + ylab("Average Polarity") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
g8 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_CHARGE)) + tmp_theme +
  ggtitle("CDR3 charge") + 
  xlab("Isotype") + ylab("Normalized Net Charge") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
g9 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_AROMATIC)) + tmp_theme +
  ggtitle("CDR3 aromatic side chain content") + 
  xlab("Isotype") + ylab("Aromatic Side Chain Content") + scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=ISOTYPE))
multiggplot(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol=3)