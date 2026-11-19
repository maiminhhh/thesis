# Integration performance across tissues (dissected datasets), relative to each method's
# unintegrated source:  Rscript analysis/compare_tissues_dissected.R
#   PCA-based methods, LIGER, SAMap, scVI, scANVI -> Unintegrated (PCA of shared genes)
#   * UINMF (incl. LIGER UINMF)                   -> Unintegrated UINMF (raw H matrix)
#   * scVI-latent                                 -> Unintegrated scVI-latent
# Overall = 0.6 bio + 0.4 batch; bio uses the four metrics every tissue has (NMI, ARI,
# ASW_label, cLISI; the isolated-label metrics only exist where a cell type is species-specific).
# Output: figures/cross_tissue/delta_overall_dissected.png (+ _components.png, .csv)
source('library.R')

scores_dirs = c(bone_marrow   = 'bone_marrow/data/05_evaluated_score/dissected_mouse/',
                colon         = 'colon/data/05_evaluated_score/dissected_mouse/',
                breast_cancer = 'breast_cancer/data/dissected_mouse/',
                fetal_liver   = 'fetal_liver/data/05_evaluated_score/dissected_mouse/')
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon',
                  breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI')

read_metrics = function(f) {
  df = read_csv(f, show_col_types = FALSE)
  tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
}

scores = imap_dfr(scores_dirs, function(dir, tissue) {
  list.files(dir, pattern = '^05\\..*\\.csv$', full.names = TRUE) |>
    map_dfr(read_metrics) |>
    mutate(group = case_when(metric %in% batch_metrics ~ 'batch',
                             metric %in% bio_metrics   ~ 'bio')) |>
    filter(!is.na(group), !is.na(value)) |>
    group_by(method, group) |>
    summarise(score = mean(value), .groups = 'drop') |>
    pivot_wider(names_from = group, values_from = score) |>
    mutate(overall = 0.6 * bio + 0.4 * batch, tissue = tissue)
})

input_of = function(m) case_when(str_detect(m, 'UINMF')       ~ 'UINMF',
                                 str_detect(m, 'scVI-latent') ~ 'scVI-latent',
                                 TRUE                         ~ 'PCA / own')
source_of = c('PCA / own' = 'Unintegrated', 'UINMF' = 'Unintegrated UINMF',
              'scVI-latent' = 'Unintegrated scVI-latent')

baselines = scores |> filter(str_starts(method, 'Unintegrated')) |>
  select(tissue, source = method, bio0 = bio, batch0 = batch, overall0 = overall)
delta = scores |>
  filter(!str_starts(method, 'Unintegrated')) |>
  mutate(input = input_of(method), source = source_of[input]) |>
  left_join(baselines, by = c('tissue', 'source')) |>
  mutate(d_overall = overall - overall0, d_bio = bio - bio0, d_batch = batch - batch0,
         input = factor(input, levels = c('PCA / own', 'UINMF', 'scVI-latent')),
         tissue = factor(tissue_labels[tissue], levels = tissue_labels))
stopifnot(!anyNA(delta$overall0), all(table(delta$method) == length(scores_dirs)))

dir.create('figures/cross_tissue', showWarnings = FALSE, recursive = TRUE)
write_csv(delta, 'figures/cross_tissue/delta_overall_dissected.csv')

# categorical tissue colours (validated palette), with shapes as the secondary encoding
tissue_cols   = setNames(c('#2a78d6', '#eb6834', '#1baf7a', '#eda100'), tissue_labels)
tissue_shapes = setNames(c(16, 17, 15, 18), tissue_labels)

# methods ordered by median gain within each input panel
order_by = function(d, col) d |> group_by(method) |> mutate(med = median(.data[[col]])) |>
  ungroup() |> mutate(method = fct_reorder(method, med))

box_plot = function(d, col, xlab) {
  ggplot(d, aes(.data[[col]], method)) +
    geom_vline(xintercept = 0, colour = 'grey55', linewidth = 0.5) +
    geom_boxplot(fill = 'grey92', colour = 'grey45', width = 0.6, outlier.shape = NA, linewidth = 0.4) +
    geom_point(aes(colour = tissue, shape = tissue), size = 2.6, alpha = 0.9,
               position = position_jitter(height = 0.15, width = 0, seed = 1)) +
    facet_grid(input ~ ., scales = 'free_y', space = 'free_y') +
    scale_colour_manual(values = tissue_cols, name = 'Tissue (dissected)') +
    scale_shape_manual(values = tissue_shapes, name = 'Tissue (dissected)') +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = 'grey92'),
          strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
          strip.text.y = element_text(angle = 0, face = 'bold'),
          legend.position = 'top', plot.title.position = 'plot')
}

p = box_plot(order_by(delta, 'd_overall'), 'd_overall',
             'Change in overall score vs unintegrated source (0.6 bio + 0.4 batch)') +
  labs(title = 'Integration gain over the unintegrated input, across four dissected tissues',
       subtitle = 'Each point is one tissue. Grey line: no change from the method\'s own unintegrated input.')
ggsave('figures/cross_tissue/delta_overall_dissected.png', p, width = 10, height = 9, dpi = 300, bg = 'white')

# same, split into the bio and batch components
long = delta |> select(method, input, tissue, d_bio, d_batch) |>
  pivot_longer(c(d_bio, d_batch), names_to = 'component', values_to = 'delta') |>
  mutate(component = factor(recode(component, d_bio = 'Bio conservation', d_batch = 'Batch removal'),
                            levels = c('Batch removal', 'Bio conservation')))
long = long |> left_join(delta |> group_by(method) |> summarise(med = median(d_overall)), by = 'method') |>
  mutate(method = fct_reorder(method, med))
p2 = box_plot(long, 'delta', 'Change vs unintegrated source') +
  facet_grid(input ~ component, scales = 'free_y', space = 'free_y') +
  labs(title = 'Batch and bio components of the gain, across four dissected tissues',
       subtitle = 'Methods ordered by median overall gain (as in the main figure).')
ggsave('figures/cross_tissue/delta_components_dissected.png', p2, width = 13, height = 9, dpi = 300, bg = 'white')

print(baselines |> mutate(across(where(is.numeric), \(x) round(x, 3))), n = Inf)
delta |> group_by(input, method) |>
  summarise(median = median(d_overall), min = min(d_overall), max = max(d_overall),
            n_positive = sum(d_overall > 0), .groups = 'drop') |>
  arrange(desc(median)) |> print(n = Inf)
