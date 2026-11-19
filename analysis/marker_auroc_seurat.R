# Marker AUROC of Seurat label transfer with the same-species reference, every dataset in one panel.
# Values: figures/cross_tissue/marker_auroc_all.csv (analysis/marker_auroc_figure.R), annotator == 'seurat'.
# One point per cell type Seurat assigns (size = % of the dataset's labelled cells given that label); black tick + value = mean over called types.
# Labelling rule: every cell type with AUROC < label_below is named. Types never assigned (NA) are dropped.
# Excluded: mouse breast cancer (its reference has no cancer cells; excluded from the subsection for the same reason).
#   Rscript --vanilla analysis/marker_auroc_seurat.R
# Output: figures/cross_tissue/marker_auroc_seurat.{png,pdf}
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(ggrepel); library(ggbeeswarm)})

# ---------------- config ----------------
tissue_labels  = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
species_labels = c(human = 'Human', mouse = 'Mouse')
mean_row       = 'Mean (all called types)'
species_cols   = c(human = '#4477AA', mouse = '#EE6677')   # as analysis/visualisations.R
exclude        = data.frame(tissue = 'breast_cancer', species = 'mouse')
label_below    = 0.7
out            = 'figures/cross_tissue/marker_auroc_seurat'
# -----------------------------------------

d = read_csv('figures/cross_tissue/marker_auroc_all.csv', show_col_types = FALSE) |> filter(annotator == 'seurat') |>
  anti_join(exclude, by = c('tissue', 'species'))
n_missing = d |> filter(celltype != mean_row) |> group_by(tissue, species) |>
  summarise(n_called = sum(!is.na(auroc)), n_types = n(), .groups = 'drop')
d = d |> left_join(n_missing, by = c('tissue', 'species')) |>
  mutate(tissue = factor(tissue_labels[tissue], levels = tissue_labels),
         row = factor(species_labels[species], levels = rev(species_labels)))
pts   = d |> filter(celltype != mean_row, !is.na(auroc)) |> mutate(celltype = gsub('_', ' ', celltype)) |>
  group_by(tissue, species) |> mutate(pct_cells = 100 * n_cells / sum(n_cells)) |> ungroup()   # share of the dataset's labelled cells
means = d |> filter(celltype == mean_row)
chk = pts |> group_by(tissue, species) |> summarise(m = mean(auroc), .groups = 'drop') |> inner_join(means, by = c('tissue', 'species'))
stopifnot(nrow(means) == 7, nrow(chk) == 7, all(abs(chk$m - chk$auroc) < 1e-9))

# labels alternate below/above the row in AUROC order, so neighbouring low types don't collide
labs_df = pts |> filter(auroc < label_below) |> group_by(tissue, species) |> arrange(auroc, .by_group = TRUE) |>
  mutate(nudge = ifelse(row_number() %% 2 == 1, -0.32, 0.32)) |> ungroup()

p = ggplot(pts, aes(auroc, row)) +
  geom_vline(xintercept = 0.5, linetype = 'dashed', colour = 'grey55', linewidth = 0.4) +
  geom_beeswarm(aes(size = pct_cells, fill = species, colour = species), shape = 21, stroke = 0.5,
                cex = 2.2, method = 'compactswarm') +   # opaque outline in species colour keeps the smallest types visible
  geom_errorbar(data = means, aes(xmin = auroc, xmax = auroc), width = 0.75, linewidth = 0.9, colour = 'black') +
  geom_label(data = means, aes(label = sprintf('mean %.2f', auroc)), nudge_y = 0.45, size = 2.6, fontface = 'bold',
             linewidth = 0, label.padding = unit(0.08, 'lines'), fill = alpha('white', 0.85)) +
  geom_text_repel(data = labs_df, aes(label = celltype), size = 2.5, colour = 'grey25', nudge_y = labs_df$nudge,
                  direction = 'both', force = 3, force_pull = 0.5, min.segment.length = 0, segment.colour = 'grey60',
                  segment.size = 0.3, box.padding = 0.2, max.overlaps = Inf, seed = 1) +
  facet_grid(tissue ~ ., switch = 'y', scales = 'free_y', space = 'free_y') +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.75))) +
  scale_fill_manual(values = alpha(species_cols, 0.7), labels = species_labels, name = 'Species',
                    guide = guide_legend(order = 1, override.aes = list(size = 4, colour = unname(species_cols)))) +
  scale_colour_manual(values = species_cols, guide = 'none') +   # outline colour shown via the fill legend's override
  guides(size = guide_legend(order = 2, override.aes = list(fill = 'grey70', colour = 'grey40'))) +
  scale_size(range = c(1.6, 8), limits = c(0, 100), breaks = c(1, 10, 25, 50, 75),
             labels = \(x) paste0(x, '%'), name = '% of dataset\'s labelled cells') +
  scale_x_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1), expand = expansion(add = 0.01)) +
  labs(x = 'Marker AUROC (Seurat, same-species reference)', y = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey93'),
        strip.placement = 'outside', strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = 'bold', size = 9, hjust = 1),
        axis.text.y = element_text(size = 8),
        legend.position = 'bottom', legend.box = 'vertical', legend.spacing.y = unit(2, 'pt'))

ggsave(paste0(out, '.png'), p, width = 7, height = 5.6, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = 7, height = 5.6, bg = 'white', device = cairo_pdf)
print(as.data.frame(means |> transmute(tissue, species, n_called, n_types, mean_auroc = round(auroc, 3))), row.names = FALSE)
cat('labelled (AUROC <', label_below, '):\n'); print(as.data.frame(pts |> filter(auroc < label_below) |> transmute(tissue, species, celltype, auroc = round(auroc, 3), n_cells, pct = round(pct_cells, 2))), row.names = FALSE)
cat('wrote', paste0(out, '.{png,pdf}\n'))
