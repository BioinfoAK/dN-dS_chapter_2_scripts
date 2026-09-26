# ---- Setup ----
setwd("/Users/nastya/Downloads")
pkgs <- c("ape", "caper", "dplyr", "readr", "broom")

# Load ape
if (!require("ape")) install.packages("ape")
library(ape)

# --- 1. Load Your Data ---
# tree.nwk: your phylogenetic tree topology
# ds_data.csv: Column 1: "Species" (matches tip names), Column 2: "dS" 
# (Important: dS must be the distance from the common ancestor to that tip)

tree <- read.tree("cereal_tree_tree.nwk")
ds_table <- read.csv("ds_data_cereal.csv", header = TRUE)
plot(tree)
# --- 2. Map dS to the Tree ---
# We replace the existing branch lengths with your transcriptome-wide dS values.
# This assumes the dS in your table represents the terminal branch leading to the species.
tree$edge.length <- ds_table$dS[match(tree$tip.label[tree$edge[,2]], ds_table$Species)]
tree$edge.length[is.na(tree$edge.length)] <- 0.001

# --- 3. Calibrate using PM-TU split ---
node_pm_tu <- getMRCA(tree, tip = c("PM", "TU"))

# If Cicer is your outgroup, the MRCA of Cicer and a Medicago tip is the tree root
node_as_tu <- getMRCA(tree, tip = c("AS", "TU"))

# 3. Create a DUAL calibration
# This anchors the internal node AND puts a 'ceiling' on the root
my_calib <- makeChronosCalib(tree, 
                             node = c(node_pm_tu, node_as_tu), 
                             age.min = c(22.4, 2.8), 
                             age.max = c(31.8, 5.3))

# 4. Run chronos with both anchors
phylo_time <- chronos(tree, model = "relaxed", calibration = my_calib)
print(my_calib)
phylo_time_clock <- chronos(tree, model="clock", calibration = my_calib)

# --- 5. Extract Per-Species Rates ---
# We isolate the terminal branches (those leading directly to your species/tips)
terminal_edges <- which(tree$edge[,2] <= Ntip(tree))

# Rate = Synonymous substitutions per site (dS) / Time (Million Years)
rates_per_my <- tree$edge.length[terminal_edges] / phylo_time_clock$edge.length[terminal_edges]
species_names <- tree$tip.label[tree$edge[terminal_edges, 2]]

# --- 6. Output Final Dataframe ---
results <- data.frame(
  Species = species_names,
  dS_Value = tree$edge.length[terminal_edges],
  Time_MY = phylo_time_clock$edge.length[terminal_edges],
  Rate_per_MY = rates_per_my,
  Rate_per_Year = rates_per_my / 1e6
)
results$omega <- ds_table$dNdS[match(results$Species, ds_table$Species)]

print(results)
write.csv(results, "cereal_mutationrates_clockwise.csv")
# Optional: Plot to see the time-calibrated tree
plot(phylo_time)
axisPhylo()

library(dplyr)
library(nlme)
library(caper)

# Assuming 'results' is the dataframe from the previous script
# Manually assign groups based on your knowledge
results <- results %>%
  mutate(Use = case_when(
    Species %in% c("AS", "TMA", "HVS", "SV") ~ "Progenitor",
    Species %in% c("LR", "PA", "PM", "BT", "HM", "HMU", "ET", "ED", "SS", "TCM", "AT", "TU") ~ "NDW"
  ))
results$Use <-as.factor(results$Use)
results$dN <- ds_table$dN
results$dNdS <- ds_table$dnds

# Check your new table
print(results)
summary(results)
# Perform the Mann-Whitney U test
wilcox_result <- wilcox.test(dN ~ Use, data = results, alternative = "two.sided")

# View the test statistics and p-value
print(wilcox_result)

# --- Step B: Fix Tree Topology for PGLS ---
# 1. Explicitly root the tree using Phalaris (or your outgroup)
# resolve.root=TRUE fixes the basal polytomy error
tree_rooted <- root(phylo_time, outgroup = "LR", resolve.root = TRUE)

# 2. Fix any remaining polytomies elsewhere in the tree
tree_rooted <- multi2di(tree_rooted)

# 3. Ensure the tree and data match perfectly
# (PGLS will fail if there is even one extra tip in the tree)
tree_rooted <- keep.tip(tree_rooted, as.character(results$Species))
plot(tree_rooted)
# --- Step C: Non-Phylogenetic Test ---
# Wilcoxon Rank Sum test
wilcox_res <- wilcox.test(Rate_per_Year ~ Use, data = results)
print(wilcox_res)
tree_rooted$node.label <- NULL
# --- Step D: PGLS (The Main Event) ---
# 1. Create the comparative data object
tree_rooted$node.label <- NULL
comp_data <- comparative.data(phy = tree_rooted, 
                              data = results, 
                              names.col = "Species", 
                              vcv = TRUE)
print(comp_data)

mytree<-comp_data$phy
mydata<-comp_data$data
print(mydata)
library(phytools)
continuous_vector <- mydata$Rate_per_Year
names(continuous_vector) <- rownames(mydata)

binary_vector <- mydata$Use
names(binary_vector) <- rownames(mydata)

# Run the phylogenetic ANOVA
anova_res <- phylANOVA(mytree, binary_vector, continuous_vector, nsim = 1000)
summary(anova_res)

library(phytools)
library(OUwie)

# Step A: Create the specific 3-column data frame for OUwie
ouwie_data <- data.frame(
  Species = rownames(mydata),
  Regime  = mydata$Use,
  Trait   = mydata$Rate_per_Year
)

# Step B: Generate a Stochastic Character Map (simmap) for the tree
# This "paints" the binary history onto the branches so OUwie knows where regimes shifted
simmap_tree <- make.simmap(mytree, binary_vector, model = "ARD", nsim = 1)

# Step C: Fit the OU models 
# Let's fit a BMS model (different rates of random walk)
bms_model <- OUwie(simmap_tree, ouwie_data, model = "BMS", simmap.tree = TRUE)

# Let's fit an OUM model (different selective optima)
oum_model <- OUwie(simmap_tree, ouwie_data, model = "OUM", simmap.tree = TRUE)

oumv_model <- OUwie(simmap_tree, ouwie_data, model = "OUMV", simmap.tree = TRUE)

# Compare them using AICc to see which evolutionary model fits best
bms_model$AICc
oum_model$AICc
oumv_model$AICc

OUwie.contour(oum_model)
# Run the model
pgls_omega <- pgls(dNdS ~ Use, data = comp_data, lambda = "ML")
hist(residuals(pgls_omega))
summary(pgls_omega)
# 2. Run the PGLS model
# This tests if 'Use' predicts 'Rate' while estimating phylogenetic signal (lambda)
pgls_model <- pgls(dS_Value ~ Use, data = comp_data, lambda = "ML")

summary(pgls_model)


library(ape)
library(ggtree)
library(ggplot2)
library(patchwork)

# Example: Define your groups and new names for tree_legume
# 'label' must match the current tip labels in your tree object
tree_info_legume <- data.frame(
  label = tree_legume$tip.label,
  new_name = c("P. sativum", "P. fulvum", "L. sativus", "L. cicera", "V. peregrina",
               "L. nigricans", "L. culinaris subsp orientalis", "M. polymorpha", "C. reticulatum",
               "C. judaicum"), # Your new names
  Use = c("Progenitor", "NDW", "Progenitor","NDW", "NDW", "NDW", "Progenitor", "NDW", "Progenitor", "NDW")    # Your group assignments
)

# Repeat for the second tree object
tree_info_tree <- data.frame(
  label = tree$tip.label,
  new_name = c("L. rigidum", "P. minor", "P. arundaceae", "Tr. monoccoccum",
               "Tr. urartu", "A. tauschii", "A. speltoides", "Ta. caput-medusae", "S. vavilovii", 
               "S. striticum", "E. distans", "E. triticeum", "H. vulgare subsp spontanaeum", 
               "H. marinum subsp glaucum", "H. murinum subsp gussoneaum", "B. tectorum"),
  Use = c("NDW", "NDW", "NDW", "Progenitor", "NDW", "NDW", "Progenitor", "NDW", "Progenitor",
            "NDW", "NDW", "NDW", "Progenitor", "NDW", "NDW", "NDW")
)


# 2. Update your colors to match these exact strings
use_colors <- c("Progenitor" = "#00BFC4", "NDW" = "#F8766D" )

# 3. Fixed Plotting Function
plot_publication_tree <- function(tr, info_df, root_node) {
  
  # Root the tree
  tr_rooted <- root(tr, outgroup = root_node, resolve.root = TRUE)
  
  # Create the plot
  # We use isTip to ensure node labels only show on internal nodes (bootstraps)
  ggtree(tr_rooted) %<+% info_df +
    geom_text2(aes(label = label, subset = !isTip), 
               hjust = 1.2, vjust = -0.5, size = 3, color = "grey40") + 
    
    # Map color to 'Use' (matching the dataframe)
    geom_tiplab(aes(label = new_name, color = Use), 
                size = 4, 
                fontface = "italic",
                offset = 0.01) +
    
    theme_tree2() +
    scale_color_manual(values = use_colors) +
    hexpand(0.3, direction = 1) + 
    
    theme(
      legend.position = "right",
      text = element_text(family = "sans"),
      legend.title = element_text(face = "bold")
    )
}

#4. Generate and Combine
# (Make sure tree_info_legume also has a 'Use' column!)
p_legume <- plot_publication_tree(tree_legume, tree_info_legume, "CJ")
p_tree   <- plot_publication_tree(tree, tree_info_tree, "LR")

final_tree_figure <- (p_legume / p_tree) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

print(final_tree_figure)

# 4. Generate and Combine
# (Make sure tree_info_legume also has a 'Use' column!)
p_legume <- plot_publication_tree(tree_legume, tree_info_legume, "CJ")
p_tree   <- plot_publication_tree(tree, tree_info_tree, "LR")

final_tree_figure <- (p_legume / p_tree) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

print(final_tree_figure)






# ---- Setup ----
setwd("/Users/nastya/Downloads")
pkgs <- c("ape", "caper", "dplyr", "readr", "broom")

# Load ape
if (!require("ape")) install.packages("ape")
library(ape)

# --- 1. Load Your Data ---
# tree.nwk: your phylogenetic tree topology
# ds_data.csv: Column 1: "Species" (matches tip names), Column 2: "dS" 
# (Important: dS must be the distance from the common ancestor to that tip)

tree_legume <- read.tree("legumes_rooted.nwk")
ds_legume <- read.csv("ds_data_legumes.csv", header = TRUE)


# --- 1. Pruning the Tree ---
# Replace 'Species_to_Remove' with the exact label of the branch you want to cut
tree_legume <- drop.tip(tree_legume, "LEC")
plot(tree_legume)
# --- 2. Map dS and Placeholder Internal Lengths ---
tree_legume$edge.length <- ds_legume$dS[match(tree_legume$tip.label[tree_legume$edge[,2]], ds_legume$Species)]
tree_legume$edge.length[is.na(tree_legume$edge.length)] <- 0.001
# Correct way to show only branch (bootstrap) labels
ggtree(tree_legume) +
  geom_tiplab() +  # This handles the species names at the ends
  geom_text2(aes(label = label, subset = !isTip), 
             hjust = 1.2, 
             vjust = -0.5, 
             size = 3)
# --- 3. Calibrate using Medicago-Vicia split ---
node_med_vic <- getMRCA(tree_legume, tip = c("MP", "VP"))

# 2. Identify the Root Node (Cicer vs the rest)
# If Cicer is your outgroup, the MRCA of Cicer and a Medicago tip is the tree root
node_root <- getMRCA(tree_legume, tip = c("CJ", "MP"))

# 3. Create a DUAL calibration
# This anchors the internal node AND puts a 'ceiling' on the root
my_calib <- makeChronosCalib(tree_legume, 
                             node = c(node_med_vic, node_root), 
                             age.min = c(21.94, 24.87), 
                             age.max = c(37.77, 51))

# 4. Run chronos with both anchors
phylo_legume <- chronos(tree_legume, model = "relaxed", calibration = my_calib)

phylo_legume_clock <- chronos(tree_legume, model ="clock", calibration = my_calib)
# --- 4. Build results_legume (matching cereal format) ---
terminal_edges <- which(phylo_legume$edge[,2] <= Ntip(phylo_legume_clock))

results_legume <- data.frame(
  Species = phylo_legume$tip.label[phylo_legume$edge[terminal_edges, 2]],
  dS_Distance = tree_legume$edge.length[terminal_edges],
  Branch_Time_MY = phylo_legume$edge.length[terminal_edges],
  Rate_per_MY = tree_legume$edge.length[terminal_edges] / phylo_legume$edge.length[terminal_edges]
)

results_legume <- results_legume %>%
  mutate(Rate_per_Year = Rate_per_MY / 1e6)
results_legume$omega <- ds_table$dNdS[match(results_legume$Species, ds_legume$Species)]
write.csv(results_legume_clock, "legume_mutationrates_clock.csv")
# --- 5. Manual Grouping (matching your categories) ---
# Adjust species names to your actual legume list
results_legume <- results_legume %>%
  mutate(Use = if_else(Species %in% c("CR", "LCO", "LS", "PS"), "Progenitor", "NDW"))

# Root and clean for PGLS
tree_final <- root(phylo_legume, outgroup = "CJ", resolve.root = TRUE)
tree_final <- multi2di(tree_final)
tree_final$node.label <- NULL
print(results_legume)
comp_data_legume <- comparative.data(phy = tree_final, data = results_legume, names.col = "Species", vcv = TRUE)
pgls_legume <- pgls(dS_Distance ~ Use, data = comp_data_legume, lambda = "ML")

summary(pgls_legume)

pgls_omega <- pgls(log(omega) ~ Use, data = comp_data_legume, lambda = "ML")
hist(residuals(pgls_omega))
summary(pgls_omega)

genome_sizes_list <- c(4312, 4606, 8232, 6860, 9310, 4125, 4125, 544, 1274, 882) 

results_legume <- results_legume %>%
  mutate(genome_size = genome_sizes_list)

genome_sizes_list <- c(2692, 4410, 5586, 6076, 4802, 4060,4889, 4312, 8036, 8624, 5390, 5390, 5390, 5390, 5390, 3234) 

results <- results %>%
  mutate(genome_size = genome_sizes_list)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork) # The best library for combining plots

# 1. Update the function to be flexible with y-axis labels
create_publication_plot_final <- function(data, y_var, y_label) {
  ggplot(data, aes(x = genome_size, y = .data[[y_var]])) +
    geom_point(aes(color = Use), alpha = 0.7, size = 3) + 
    geom_text_repel(aes(label = Species), 
                    size = 3.8, 
                    box.padding = 0.5, 
                    max.overlaps = 15,
                    segment.color = 'grey50') +
    geom_smooth(method = "lm", color = "firebrick", linetype = "dashed", se = TRUE) +
    stat_cor(aes(label = ..rr.label..), label.x.npc = "left", label.y.npc = "top") +
    theme_bw() +
    labs(x = "Genome size (Mb)", y = y_label) +
    theme(
      text = element_text(size = 14, family = "sans"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14, color = "black"),
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "none" # Hide legend on individual plots to save space
    )
}

# 2. Generate the two individual plots
p1 <- create_publication_plot_final(results_legume, "dS_Distance", "dS")
p2 <- create_publication_plot_final(results_legume, "Rate_per_Year", "Rate per Year")

# 3. Combine them using patchwork
# We use 'guides = "collect"' to create a single shared legend for both plots
combined_plot <- (p1 | p2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A') & # Adds (A) and (B) labels automatically
  theme(legend.position = "bottom")

# 4. Display the result
print(combined_plot)

library(ggplot2)
library(scales)
library(ggpubr)
library(ggrepel) # Essential for non-overlapping labels

ggplot(results_legume_clock, aes(x = dS_Distance, y = Rate_per_Year)) +
  # 1. Add points
  geom_point(aes(color = Use), alpha = 0.7, size = 5) +
  
  # 2. Add Species Labels
  # box.padding helps keep the text from touching the dots
  # max.overlaps prevents the plot from becoming a mess if you have 100+ species
  geom_text_repel(aes(label = Species), 
                  size = 7, 
                  box.padding = 0.5, 
                  max.overlaps = 20,
                  segment.color = 'grey50') +
  
  # 3. Add Line of Best Fit
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  
  # 4. Log Scale and Ticks
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks(sides = "l") +
  
  # 5. R-squared Statistic
  stat_cor(aes(label = ..rr.label..), 
           method = "pearson", 
           label.y.npc = "top", 
           label.x.npc = "left") +
  
  # 4. FIXED AXIS RANGES
  # x: 0.2 to 0.35 | y: 1e-9 to 1e-6
  coord_cartesian(xlim = c(0.22, 0.35), ylim = c(1e-9, 1e-6)) +
  
  # 6. Aesthetics
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Species-Level Divergence vs. Mutation Rate",
    x = expression(italic(dS)~~Value),
    y = "Mutation Rate (substitutions/site/year)",
    color = "Species Category"
  ) +
  theme(legend.position = "bottom")

# 1. Define a range of lambda values to test (usually on a log scale)
lambdas <- 10^(-1:10) 

# 2. Run chronos with cross-validation
# 'model' can be "relaxed" or "discrete"
cv_results <- chronos(tree_rooted, model = "relaxed", 
                      lambda = lambdas, 
                      control = chronos.control(nb.rate.cat = 10, 
                                                dual.iterate = TRUE))

# 3. Check the attributes to see the CV scores
attr(cv_results, "index") # The estimated dates
attr(cv_results, "rates")           # The estimated rates per branch



results_legume <- results_legume %>% 
  rename(dS_Value = dS_Distance )

library(ggplot2)
library(ggrepel)
library(scales)
library(ggpubr)
library(patchwork) # For side-by-side plots


# 1. Custom function to ensure consistent formatting across both datasets
create_scientific_plot <- function(data, x_label = expression(italic(dS)~~Value)) {
  ggplot(data, aes(x = dS_Value, y = Rate_per_Year)) +
    # Use a neutral color if 'Use' column is missing, or map it if available
    geom_point(aes(color = Use), alpha = 0.7, size = 3) + # Reduced point size for side-by-side
    
    # Text size 3.8mm roughly equals 11pt font in Word
    geom_text_repel(aes(label = Species), 
                    size = 3.8, 
                    box.padding = 0.5, 
                    max.overlaps = 15,
                    segment.color = 'grey50') +
    
    geom_smooth(method = "lm", color = "firebrick", linetype = "dashed", se = TRUE) +
    
    # Log transformation with scientific notation
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "l") +
    
    # Statistical label
    stat_cor(aes(label = ..rr.label..), 
             method = "pearson", 
             label.y.npc = "top", 
             label.x.npc = "left",
             size = 6) + 
    
    # Fixed coordinate ranges as requested
    coord_cartesian(xlim = c(0.22, 0.35), ylim = c(1e-9, 1e-6)) +
    
    theme_bw() +
    labs(x = x_label, y = "Mutation Rate (subst/site/year)") +
    theme(
      # Standardizing font sizes for publication (10-12pt)
      text = element_text(size = 14, family = "sans"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14, color = "black"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# 2. Generate individual plot objects
p1 <- create_scientific_plot(results_legume)
p2 <- create_scientific_plot(results)

# 3. Combine using patchwork with A) and B) tagging
final_figure <- p1 + p2 + 
  plot_annotation(
    tag_levels = 'A', 
    tag_suffix = '.',
    theme = theme(plot.tag = element_text(size = 20, face = "bold"))
  )

# 4. Display and Save
print(final_figure)

# Save with dimensions suitable for a Word document (Full Page Width)
ggsave("Figure_1_Divergence_Analysis.png", final_figure, width = 10, height = 5, dpi = 300)
# Function to create the plot to avoid repeating code
create_divergence_plot <- function(data, plot_title) {
  ggplot(data, aes(x = dS_Value, y = Rate_per_Year)) +
    geom_point(aes(color = Use), alpha = 0.7, size = 3) + # Reduced point size for side-by-side
    
    # Text size ~10-12pt in Word
    geom_text_repel(aes(label = Species), 
                    size = 3.8, 
                    box.padding = 0.5, 
                    max.overlaps = 20,
                    segment.color = 'grey50') +
    
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
    
    scale_y_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "l") +
    
    stat_cor(aes(label = ..rr.label..), 
             method = "pearson", 
             label.y.npc = "top", 
             label.x.npc = "left",
             size = 4) + # Size of the R2 text
    
    coord_cartesian(xlim = c(0.22, 0.35), ylim = c(1e-9, 1e-6)) +
    
    theme_bw() +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = plot_title,
      x = expression(italic(dS)~~Value),
      y = "Mutation Rate (subst/site/year)",
      color = "Category"
    ) +
    theme(
      legend.position = "bottom",
      # Set specific font sizes for Word (11pt)
      text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      plot.title = element_text(size = 12, hjust = 0.5)
    )
}

# 1. Generate the two plots
# Replace 'dataset1' and 'dataset2' with your actual data frames
p1 <- create_divergence_plot(results, "Cereals")
p2 <- create_divergence_plot(results_legume, "Legume")

# 2. Combine them side-by-side
# guides = "collect" moves the legend to one shared spot if they are the same
combined_plot <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# 3. View the plot
print(combined_plot)

# 4. EXPORTING FOR WORD
# Using a width of 10-12 inches for a 2-panel plot ensures clarity
ggsave("mutation_plots_combined.png", combined_plot, width = 12, height = 6, dpi = 300)

wilcox_res <- wilcox.test(Rate_per_Year ~ Use, data = results_legume)
print(wilcox_res)

# 1. Extract the rates calculated by chronos
# chronos stores rates in an attribute called "rates"
all_rates <- attr(phylo_legume, "rates")

# 2. Identify which branches belong to which "Use" group
# This requires mapping your tip data to the edge matrix of the tree
library(phytools)

# Let's say you want to compare the terminal branches (the tips)
edge_indices <- match(1:length(phylo_legume$tip.label), phylo_legume$edge[,2])
terminal_rates <- all_rates[edge_indices]

# Combine into a data frame for testing
rate_comparison <- data.frame(
  Species = phylo_legume$tip.label,
  Calculated_Rate = terminal_rates
)

# Merge with your 'Use' data
rate_comparison <- merge(rate_comparison, 
                         comp_data_legume$data, 
                         by.x = "Species", 
                         by.y = "row.names")
# 3. Perform a t-test or Wilcoxon test on the rates themselves
t.test(Calculated_Rate ~ Use, data = rate_comparison)
wilcox_test <- wilcox.test(Calculated_Rate ~ Use, data = rate_comparison)
print(wilcox_test)
colnames(rate_comparison)
colnames(comp_data_legume$data)

# Plot the tree with branch rates
# We scale the colors based on the rates found in chronos
library(viridis)
rescaled_rates <- (all_rates - min(all_rates)) / (max(all_rates) - min(all_rates))
cols <- viridis(100)[as.numeric(cut(all_rates, breaks = 100))]

plot(phylo_legume, edge.color = cols, edge.width = 2, main = "Evolutionary Rates across Lineages")

