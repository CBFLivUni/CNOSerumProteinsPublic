

# Plot
## Not filtering serum QC-only proteins
boxPlt <- ggplot(perSmplNACross_dropProts, aes(x=factor(PERCENTAGE, levels=c("Unfiltered","25%","20%","15%","10%","5%","0%")), y=PROP_MISSING)) + 
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=2) +
  geom_violin(fill="#99E6FF",
              width=0.75,
              alpha=0.75) + 
  geom_boxplot(width=0.25,
               outlier.size=4,
               fill="white") +
  ylim(0, 1) +
  # geom_point(size=1, alpha=0.85) +
  geom_hline(yintercept=maxPropMissing, linetype="dashed", color="darkgrey", linewidth=2) +
  ylab("Proportion of Proteins") + xlab("Threshold of Missing Values") +
  theme_classic(base_size=48) +
  theme(axis.title=element_text(size=56))
boxPltByPlate <- ggplot(perSmplNACross_dropProts, aes(x=factor(PERCENTAGE, levels=c("Unfiltered","25%","20%","15%","10%","5%","0%")), y=PROP_MISSING)) + 
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=2) +
  geom_violin(fill="#99E6FF",
              width=0.75,
              alpha=0.75) + 
  geom_boxplot(width=0.25,
               outlier.size=4,
               fill="white") +
  ylim(0, 1) +
  facet_grid(cols=vars(PLATE)) +
  # geom_point(size=1, alpha=0.85) +
  geom_hline(yintercept=maxPropMissing, linetype="dashed", color="darkgrey", linewidth=2) +
  ylab("Proportion of Proteins") + xlab("Threshold of Missing Values") +
  theme_classic(base_size=60) +
  theme(axis.title=element_text(size=56))
barPlt <- ggplot(totalNProtsCross, aes(x=factor(PERCENTAGE, levels=c("Unfiltered","25%","20%","15%","10%","5%","0%")), y=FREQ)) + 
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=2) +
  geom_bar(stat="identity", color="black", fill="#99E6FF") + 
  geom_text(aes(label = FREQ), vjust = -0.85, size=12) +
  ylim(0, 1350) +
  ylab("Number of Proteins") + xlab("Threshold of Missing Values") +
  theme_classic(base_size=48) +
  theme(axis.title=element_text(size=56),
        axis.text.x=element_text(size=32))

## Filtering by serum QC proteins
boxPlt_serumQCOnly <- ggplot(perSmplNACross_dropProts_serumQCOnly, aes(x=factor(PERCENTAGE, levels=c("Unfiltered","25%","20%","15%","10%","5%","0%")), y=PROP_MISSING)) + 
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=2) +
  geom_violin(fill="#99E6FF",
              width=0.75,
              alpha=0.75) + 
  geom_boxplot(width=0.25,
               outlier.size=4,
               fill="white") +
  ylim(0, 1) +
  # geom_point(size=1, alpha=0.85) +
  geom_hline(yintercept=maxPropMissing, linetype="dashed", color="darkgrey", linewidth=2) +
  ylab("Proportion of Proteins") + xlab("Threshold of Missing Values") +
  # ggtitle("Only Proteins Dected in all Serum QC Samples") +
  theme_classic(base_size=48) +
  theme(axis.title=element_text(size=56),
        axis.text.x=element_text(size=32))
boxPltByPlate_serumQCOnly <- ggplot(perSmplNACross_dropProts_serumQCOnly, aes(x=factor(PERCENTAGE, levels=c("Unfiltered","25%","20%","15%","10%","5%","0%")), y=PROP_MISSING)) + 
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=2) +
  geom_violin(fill="#99E6FF",
              width=0.75,
              alpha=0.75) + 
  geom_boxplot(width=0.25,
               outlier.size=4,
               fill="white") +
  ylim(0, 1) +
  facet_grid(cols=vars(PLATE)) +
  # geom_point(size=1, alpha=0.85) +
  geom_hline(yintercept=maxPropMissing, linetype="dashed", color="darkgrey", linewidth=2) +
  ylab("Proportion of Proteins") + xlab("Threshold of Missing Values") +
  # ggtitle("Only Proteins Dected in all Serum QC Samples") +
  theme_classic(base_size=48) +
  theme(axis.title=element_text(size=56),
        axis.text.x=element_text(size=32))
barPlt_serumQCOnly <- ggplot(totalNProtsCross_serumQC, aes(x=factor(PERCENTAGE, levels=c("Unfiltered","25%","20%","15%","10%","5%","0%")), y=FREQ)) + 
  geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=2) +
  geom_bar(stat="identity", color="black", fill="#99E6FF") + 
  geom_text(aes(label = FREQ), vjust = -0.85, size=12) +
  ylim(0, 1350) +
  ylab("Number of Proteins") + xlab("Threshold of Missing Values") +
  # ggtitle("Only Proteins Dected in all Serum QC Samples") +
  theme_classic(base_size=60) +
  theme(axis.title=element_text(size=56),
        axis.text.x=element_text(size=32))
