# Figures on what EBCC gives beyond a majority vote. Run from the repository root.
# Reads the saved 02 annotation objects, cached EBCC fits and 05 label-sensitivity CSVs;
# no integration, EBCC fit or 05 scoring is rerun.
#   1. coverage:    % of cells each labelling labels (vote ties and SingleR's pruned calls are unlabelled)
#   2. confidence:  EBCC's per-cell max posterior vs annotator agreement and marker AUROC (by quintile)
#   3. sensitivity: Kendall tau of integration rankings under EBCC vs best-annotator labels (05 CSVs)
#   4. hard cells:  evaluation-half cells where the vote is tied, scored for the labellings that
#                   label them, incl. the vote with ties broken (to same-species Seurat, or at random);
#                   paired bootstrap difference to EBCC (chosen f), as in analysis/vote_vs_ebcc.R
# Inputs come from the 02 h5ad (analysis/ebcc_support_export.py; cached in h5ad_inputs/), not the Seurat .rds:
# Seurat no longer loads under the r45 env. Same cells, order, data layer and markers.
# Usage:  Rscript analysis/ebcc_support.R <tissue> <species>   one dataset -> figures/cross_tissue/ebcc_support/<tissue>_<species>_*.csv
#         Rscript analysis/ebcc_support.R                      figures (figure 3 needs no per-dataset step)

# helpers from analysis/vote_vs_ebcc.R (and through it analysis/visualisations.R): tissues, annotators, majority
# helpers, fit_file, fit_labels, auc_rank, score_set, B, key, tissue_order, save_plot, species_cols
# (Seurat-free: analysis/visualisations.R's config/helpers without its source()/library() calls, then analysis/vote_vs_ebcc.R's helpers)
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
vis = parse('analysis/visualisations.R'); vis_h = sapply(vis, \(e) deparse(e)[1])
for (i in seq_len(which(startsWith(vis_h, 'for (tissue in names(tissues))')) - 1))
  if (!grepl('^(source|library)\\(|run_only', vis_h[i])) eval(vis[[i]], globalenv())
ex = parse('analysis/vote_vs_ebcc.R'); ex_h = sapply(ex, \(e) deparse(e)[1])
for (i in (which(startsWith(ex_h, 'for (e in ex'))[1] + 1):(which(startsWith(ex_h, 'args = commandArgs')) - 1))
  eval(ex[[i]], globalenv())
py_scib = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'scib/bin/python')
inputs_dir = 'figures/cross_tissue/ebcc_support/h5ad_inputs'
# one dataset's metadata (rownames = cells) and marker score per cell (cells x marker groups with >= 1 panel gene)
h5ad_inputs = function(tissue, species, a) {
  dir.create(inputs_dir, showWarnings = FALSE, recursive = TRUE)
  stem = file.path(inputs_dir, sprintf('%s_%s', tissue, species))
  if (!file.exists(paste0(stem, '_score.csv.gz'))) {
    jsonlite::write_json(annotation_markers(a), paste0(stem, '_markers.json'))
    stopifnot(system2(py_scib, c('-I', 'analysis/ebcc_support_export.py', sub('\\.rds$', '.h5ad', a$obj),
                                 paste0(stem, '_markers.json'), stem)) == 0)
  }
  meta = read_csv(paste0(stem, '_obs.csv.gz'), col_types = cols(.default = 'c', ebcc_maxpost = 'd')) |>
    as.data.frame() |> column_to_rownames('cell')
  score = as.matrix(read_csv(paste0(stem, '_score.csv.gz'), col_types = cols(.default = 'd'), name_repair = 'minimal'))
  list(meta = meta, score = score)
}

out_dir = 'figures/cross_tissue/ebcc_support'
species_scale = scale_colour_manual(values = c(human = species_cols[['hsapiens']], mouse = species_cols[['mmusculus']]),
                                    name = NULL)
tissue_name = \(t) str_to_sentence(str_replace(t, '_', ' '))

compute_one = function(tissue, species) {
  t0 = Sys.time()
  tick = function(msg) message(sprintf('[%s %s %5.1f min] %s', tissue, species,
                                       as.numeric(Sys.time() - t0, units = 'mins'), msg))
  a = tissues[[tissue]]$annotation[[species]]
  inp = h5ad_inputs(tissue, species, a)
  meta = inp$meta
  cells = rownames(meta)
  tick('inputs loaded')
  set.seed(1)
  held = !seq_len(length(cells)) %in% sample(length(cells), floor(length(cells) / 2))
  score = inp$score[held, , drop = FALSE]                                      # held cells x types
  rm(inp); gc()
  ch = read_csv(file.path(dirname(a$obj), str_glue('{species}_ebcc_frac_auroc_chosen.csv')), show_col_types = FALSE)

  # the four annotators; vote = most frequent label, with the tied labels kept for tie-breaking
  ann = sapply(annotators[c('seurat', 'singleR', 'seurat_cross', 'singleR_cross')], \(col) as.character(meta[[col]]))
  top = apply(ann, 1, function(x) { tab = table(x[!is.na(x)]); names(tab)[tab == max(tab)] }, simplify = FALSE)
  tied = lengths(top) > 1
  seurat = ann[, 'seurat']
  set.seed(3)
  raw = list('Seurat'                  = seurat,
             'SingleR'                 = ann[, 'singleR'],
             'Seurat (cross)'          = ann[, 'seurat_cross'],
             'SingleR (cross)'         = ann[, 'singleR_cross'],
             'Majority vote'           = ifelse(tied, NA_character_, vapply(top, `[`, '', 1)),
             'Vote, ties to Seurat'    = ifelse(tied, seurat, vapply(top, `[`, '', 1)),
             'Vote, ties at random'    = vapply(top, \(x) if (length(x) > 1) sample(x, 1) else x, ''),
             'EBCC (f = 0)'            = fit_labels(fit_file(tissue, species, 0, ch$M), cells),
             'EBCC (chosen f)'         = as.character(meta$ebcc_label))
  tick('labels ready')
  info = \(d) mutate(d, tissue = tissue, species = species, chosen_frac = ch$frac, .before = 1)

  # 1. coverage over all cells
  coverage = tibble(labelling = names(raw), n_cells = length(cells),
                    labelled = vapply(raw, \(l) mean(!is.na(l)), numeric(1))) |> info()

  labs = lapply(raw, \(l) key(l[held]))
  # 2. EBCC confidence quintiles on the evaluation half: annotator agreement and marker AUROC
  conf = meta$ebcc_maxpost[held]
  bin = dplyr::ntile(conf, 5)
  unanimous = apply(ann[held, , drop = FALSE], 1, \(x) { x = x[!is.na(x)]; length(x) > 0 && all(x == x[1]) })
  mean_auc_in = function(l, idx) {
    types = colnames(score)[vapply(colnames(score), \(k) { p = l[idx] == key(k); any(p) && !all(p) }, logical(1))]
    mean(vapply(types, \(k) auroc(score[idx, k], l[idx] == key(k)), numeric(1)))
  }
  confidence = map_dfr(1:5, function(b) {
    idx = which(bin == b)
    tibble(bin = b, n = length(idx), maxpost_min = min(conf[idx]), maxpost_max = max(conf[idx]),
           maxpost_mean = mean(conf[idx]), unanimous = mean(unanimous[idx]), vote_tied = mean(tied[held][idx]),
           auroc_ebcc = mean_auc_in(labs[['EBCC (chosen f)']], idx), auroc_seurat = mean_auc_in(labs[['Seurat']], idx))
  }) |> info()
  tick('confidence done')

  # 4. hard (vote tied) and all evaluation cells, labellings that label them; paired bootstrap
  hard_sets = c('Seurat', 'Seurat (cross)', 'Vote, ties to Seurat', 'Vote, ties at random', 'EBCC (f = 0)', 'EBCC (chosen f)')
  hard_idx = which(tied[held])
  all_idx = which(Reduce(`&`, lapply(labs[hard_sets], \(l) !is.na(l))))
  set.seed(2)
  hard = bind_rows(score_set(score, labs[hard_sets], hard_idx, 'hard'),
                   score_set(score, labs[hard_sets], all_idx, 'all')) |>
    info() |> mutate(n_held = sum(held))
  # SingleR prunes calls, so each SingleR labelling is scored on the cells it labels, paired with
  # EBCC (chosen f) on those cells (after the above, so its resamples are unchanged)
  singler = purrr::map_dfr(c('SingleR', 'SingleR (cross)'), \(s)
    purrr::imap_dfr(list(hard = hard_idx, all = all_idx), \(idx, set_name)
      score_set(score, labs[c(s, 'EBCC (chosen f)')], idx[!is.na(labs[[s]][idx])], set_name) |> filter(labelling == s)))
  hard = bind_rows(hard, info(singler) |> mutate(n_held = sum(held)))
  tick('hard cells done')

  list(coverage = coverage, confidence = confidence, hard = hard)
}

args = commandArgs(TRUE)
if (length(args) == 2) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  r = compute_one(args[1], args[2])
  for (n in names(r)) write_csv(r[[n]], file.path(out_dir, sprintf('%s_%s_%s.csv', args[1], args[2], n)))
  message('done ', args[1], ' ', args[2])
  quit(save = 'no')
}

# ---- 3. label sensitivity (05 CSVs; colon skipped by design) ----
# alternative labels per tissue: data/05_evaluated_score/label_sensitivity_choice.json
alt_labels = c(bone_marrow   = 'human Seurat, mouse EBCC',
               fetal_liver   = 'human Seurat, mouse EBCC',
               breast_cancer = 'human Seurat, mouse Seurat (cross)')
read_sens = function(t) {
  d = read_csv(file.path(t, 'data/05_evaluated_score/label_sensitivity.csv'), show_col_types = FALSE)
  if ('kendall_tau' %in% names(d)) {          # bone_marrow: one row per dataset and score
    split3 = \(x) str_split(x, ',\\s*')
    d = d |> transmute(dataset, score, tau = kendall_tau,
                       same_top3 = map2_lgl(split3(top3_ebcc), split3(top3_alt), \(a, b) setequal(a, b)))
  } else {
    d = bind_rows(d |> transmute(dataset, score = 'overall', tau = tau_overall, same_top3 = same_top3_overall),
                  d |> transmute(dataset, score = 'bio', tau = tau_bio, same_top3 = same_top3_bio))
  }
  d |> mutate(tissue = t, dataset = str_remove(dataset, '_mouse$'))
}
sens = map_dfr(names(alt_labels), read_sens) |>
  mutate(dataset = factor(dataset, levels = c('full', 'dissected', 'genes2', 'genes4', 'genes6', 'genes8')),
         score = factor(score, levels = c('overall', 'bio'), labels = c('Overall score (0.6 bio + 0.4 batch)', 'Bio-conservation score')),
         tissue = factor(tissue, levels = tissue_order),
         top3 = if_else(same_top3, 'same top 3', 'top 3 differ'))
write_csv(sens, 'figures/cross_tissue/label_sensitivity.csv')
tissue_cols = c(colon = '#009E73', bone_marrow = '#E69F00', fetal_liver = '#56B4E9', breast_cancer = '#CC79A7')
p3 = ggplot(sens, aes(dataset, tau, colour = tissue, group = tissue)) +
  geom_line(linewidth = 0.5, alpha = 0.6) +
  geom_point(aes(shape = top3), size = 2.8, stroke = 1) +
  facet_wrap(~ score) +
  scale_colour_manual(values = tissue_cols, labels = \(x) sprintf('%s (alternative: %s)', tissue_name(x), alt_labels[x]),
                      name = NULL, drop = TRUE) +
  scale_shape_manual(values = c('same top 3' = 16, 'top 3 differ' = 1), name = NULL) +
  scale_y_continuous(limits = c(0.6, 1)) +
  labs(x = 'Integration case', y = 'Kendall tau: method ranking under EBCC vs alternative labels',
       caption = paste('26 integration methods ranked by the 05 scores under EBCC labels and under the best single annotator per species (by marker AUROC);',
                       'tau = 1: identical ranking. Filled: same top-3 methods. Colon not run: EBCC was already the best annotator in both species.',
                       sep = '\n')) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top', legend.direction = 'vertical', panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = 'bold', hjust = 0),
        axis.text.x = element_text(angle = 30, hjust = 1), plot.caption = element_text(hjust = 0, size = 9))
save_plot(p3, 'figures/cross_tissue/ebcc_label_sensitivity.png', width = 11, height = 6.5)
ggsave('figures/cross_tissue/ebcc_label_sensitivity.pdf', p3, width = 11, height = 6.5)
message('done figure 3')

# ---- figures 1, 2, 4: need every per-dataset CSV ----
read_all = function(what) map_dfr(tissue_order, \(t) map_dfr(c('human', 'mouse'), \(sp) {
  f = file.path(out_dir, sprintf('%s_%s_%s.csv', t, sp, what))
  if (file.exists(f)) read_csv(f, show_col_types = FALSE)
}))
coverage = read_all('coverage')
if (nrow(coverage) < 8 * 9) { message('per-dataset CSVs incomplete: figures 1, 2, 4 skipped'); quit(save = 'no') }
confidence = read_all('confidence'); hard = read_all('hard')
write_csv(coverage, 'figures/cross_tissue/ebcc_coverage.csv')
write_csv(confidence, 'figures/cross_tissue/ebcc_confidence.csv')
write_csv(hard, 'figures/cross_tissue/ebcc_hard_cells.csv')

# ---- 1. coverage ----
cov_sets = c('Seurat', 'SingleR', 'Seurat (cross)', 'SingleR (cross)', 'Majority vote', 'EBCC (chosen f)')
d1 = coverage |> filter(labelling %in% cov_sets) |>
  mutate(labelling = factor(labelling, levels = rev(cov_sets)), tissue = factor(tissue, levels = tissue_order))
p1 = ggplot(d1, aes(100 * labelled, labelling, fill = species)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf('%.0f%%', 100 * labelled)), position = position_dodge(width = 0.75),
            hjust = 1.1, size = 3, colour = 'white') +
  facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(tissue_name)) +
  scale_fill_manual(values = c(human = species_cols[['hsapiens']], mouse = species_cols[['mmusculus']]), name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0, 0)) +
  labs(x = '% of cells labelled', y = NULL,
       caption = 'Unlabelled: SingleR pruned calls; majority vote ties (no single most frequent label among the four annotators). EBCC labels every cell.') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top', panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = 'bold', hjust = 0), plot.caption = element_text(hjust = 0, size = 9),
        panel.spacing = unit(1, 'lines'))
save_plot(p1, 'figures/cross_tissue/ebcc_coverage.png', width = 13, height = 4.5)
ggsave('figures/cross_tissue/ebcc_coverage.pdf', p1, width = 13, height = 4.5)

# ---- 2. confidence ----
d2 = confidence |> mutate(tissue = factor(tissue, levels = tissue_order))
conf_panel = function(y, ylab, extra = NULL) {
  ggplot(d2, aes(bin, .data[[y]], colour = species)) +
    extra +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(tissue_name)) +
    species_scale +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(legend.position = 'top', panel.grid.minor = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = 'bold', hjust = 0))
}
p2a = conf_panel('unanimous', '% cells: all annotators agree') + scale_y_continuous(labels = scales::percent)
p2b = conf_panel('auroc_ebcc', 'Mean marker AUROC',
                 geom_line(aes(y = auroc_seurat), linetype = 'dashed', alpha = 0.6)) +
  labs(x = 'EBCC confidence (max posterior), quintile: 1 = lowest 20% of cells, 5 = highest') +
  guides(colour = 'none')                                                # legend once, from panel A
p2 = patchwork::wrap_plots(p2a, p2b, ncol = 1, guides = 'collect') +
  patchwork::plot_annotation(
    tag_levels = 'A',
    caption = paste('Evaluation half of cells, split into fifths by EBCC\'s max posterior per cell (quintiles within each dataset).',
                    'A: share of cells where all four annotators give the same label (EBCC\'s posterior is driven by annotator agreement, so A partly holds by construction).',
                    'B: solid = EBCC labels, dashed = same-species Seurat on the same cells. Breast cancer: EBCC\'s max posterior is nearly constant (mouse: most cells 0.539),',
                    'so some fifths hold a single EBCC label and have no AUROC (missing points).',
                    sep = '\n'),
    theme = theme(plot.caption = element_text(hjust = 0, size = 9))) &
  theme(legend.position = 'top')
save_plot(p2, 'figures/cross_tissue/ebcc_confidence.png', width = 13, height = 8)
ggsave('figures/cross_tissue/ebcc_confidence.pdf', p2, width = 13, height = 8)

# ---- 4. all, consensus and hard cells ----
consensus = read_csv('figures/cross_tissue/vote_vs_ebcc.csv', show_col_types = FALSE) |> filter(set == 'common')
# A: coverage; B: all evaluation-half cells; C: consensus cells (vote not tied, SingleR labels; = heldout_auroc.png);
# D: hard cells (vote tied). Vote ties to Seurat left out: on tied cells it is same-species Seurat.
cell_levels = c('Seurat', 'SingleR', 'Seurat (cross)', 'SingleR (cross)', 'Majority vote', 'EBCC (f = 0)', 'EBCC (chosen f)')
fig4_levels = c('Majority vote', 'EBCC (f = 0)', 'EBCC (chosen f)')   # rows of auroc_by_cell_set (the single annotators: best_annotator.png)
# thin rules between labelling rows, so each row's two species points/bars read as one group
row_dividers = geom_hline(yintercept = seq_along(fig4_levels)[-1] - 0.5, colour = 'grey75', linewidth = 0.3)
tissue_title = element_text(face = 'bold', hjust = 0, size = 14)   # tissue names, panels A-D
# one row of panels B-D: a plot per tissue with the tissue as a large title (cells used and chosen f: caption file);
# y labels on the first tissue only
cell_row = function(d) {
  d = filter(d, labelling != 'Vote, ties to Seurat') |> mutate(labelling = recode(labelling, 'Vote, ties at random' = 'Majority vote')) |>
    filter(labelling %in% fig4_levels)
  r = range(c(d$lo, d$hi, d$mean_auroc), na.rm = TRUE); w = diff(r)
  xlim = c(r[1] - 0.05 * w, r[2] + 0.18 * w)              # one x range per row, as the shared facet axis was (room for values)
  purrr::imap(tissue_order, function(t, i) {
    x = filter(d, tissue == t)
    p = plot_auroc(x, fig4_levels, xlab = NULL) + row_dividers +
      coord_cartesian(xlim = xlim, ylim = c(0.5, length(fig4_levels) + 0.5), expand = FALSE) +
      labs(title = tissue_name(t)) +
      theme(panel.grid.major.x = element_line(colour = 'grey94', linewidth = 0.25),   # faint: rows (dividers) stay the main structure
            strip.text = element_blank(), plot.title = tissue_title,
            plot.subtitle = element_blank())
    if (i > 1) p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p
  })
}
x_title = function(label) ggplot() + annotate('text', x = 0.5, y = 0.5, label = label, size = 12 / .pt) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_void()
# coverage: % of all cells each labelling labels (majority vote: ties unlabelled), rows as in the AUROC panels
cov_panel = function(title) {
  d = coverage |> filter(labelling %in% fig4_levels) |>
    mutate(labelling = factor(labelling, levels = rev(fig4_levels)), tissue = factor(tissue, levels = tissue_order))
  row = which(rev(fig4_levels) == 'EBCC (chosen f)')
  ggplot(d, aes(100 * labelled, labelling, fill = species)) +
    annotate('rect', xmin = -Inf, xmax = Inf, ymin = row - 0.5, ymax = row + 0.5, fill = '#FCEFC7') +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    row_dividers +
    geom_text(aes(label = sprintf('%.0f%%', 100 * labelled)), position = position_dodge(width = 0.75),
              hjust = 1.1, size = 2.7, colour = 'white') +
    facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(tissue_name)) +
    scale_fill_manual(values = c(human = species_cols[['hsapiens']], mouse = species_cols[['mmusculus']]), guide = 'none') +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0, 0)) +
    labs(title = title, x = '% of cells labelled', y = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          strip.text = tissue_title, panel.spacing = unit(1.8, 'lines'),
          axis.text.y = element_text(face = ifelse(rev(fig4_levels) == 'EBCC (chosen f)', 'bold', 'plain')),
          plot.title = element_text(face = 'bold'))
}
# flat grid: A spans the width; B-D are 4 tissue plots each, with one x-axis title row below
rows = list(cell_row(filter(hard, set == 'all')), cell_row(consensus), cell_row(filter(hard, set == 'hard')))
plots = c(list(cov_panel(NULL)), unlist(lapply(rows, \(r) c(r, list(x_title('Mean marker AUROC')))), recursive = FALSE))
layout = do.call(c, c(list(patchwork::area(1, 1, 1, 4)),
                      unlist(lapply(0:2, \(k) c(lapply(1:4, \(j) patchwork::area(2 + 2 * k, j)),
                                                list(patchwork::area(3 + 2 * k, 1, 3 + 2 * k, 4)))), recursive = FALSE)))
tags = c('A', unlist(lapply(c('B', 'C', 'D'), \(l) c(l, '', '', '', ''))))
p4 = patchwork::wrap_plots(plots, design = layout, guides = 'collect', heights = c(0.8, 1, 0.1, 1, 0.1, 1, 0.1)) +
  patchwork::plot_annotation(tag_levels = list(tags)) &
  theme(legend.position = 'top', plot.tag = element_text(face = 'bold', size = 14))
# caption for the thesis text, not drawn
writeLines(paste('A: coverage, the share of all cells each labelling labels; unlabelled = SingleR pruned calls, majority vote ties.',
                 'B: all cells; C: consensus cells (the four annotators give a majority and SingleR labels the cell);',
                 'D: hard cells (the four annotators tie, no majority).',
                 'B-D: evaluation half of cells (not used to choose f).',
                 'Mean AUROC over the types each labelling calls. Majority vote: in B and D ties are broken by drawing one tied label at random (seed 3), C has no ties.',
                 'Rows: majority vote, EBCC with f = 0 (no trusted prior) and EBCC with the chosen f (the single annotators: best_annotator.png).', sep = '\n'),
           'figures/cross_tissue/auroc_by_cell_set_caption.txt')
# cells used per panel (% of the evaluation half) and chosen f, formerly drawn under each tissue
cell_counts = bind_rows(filter(hard, set == 'all') |> mutate(panel = 'B'), consensus |> mutate(panel = 'C'),
                        filter(hard, set == 'hard') |> mutate(panel = 'D')) |>
  filter(labelling == 'EBCC (chosen f)') |>
  transmute(panel, tissue, species, chosen_frac, n_cells, pct_of_half = round(100 * n_cells / n_held, 1)) |>
  arrange(panel, factor(tissue, levels = tissue_order), desc(species))
write_csv(cell_counts, 'figures/cross_tissue/auroc_by_cell_set_cells.csv')
save_plot(p4, 'figures/cross_tissue/auroc_by_cell_set.png', width = 14, height = 10.5)
ggsave('figures/cross_tissue/auroc_by_cell_set.pdf', p4, width = 14, height = 10.5, device = cairo_pdf)
# ---- 5. best labelling per dataset: all cells (= auroc_by_cell_set panel B), rank within dataset ----
annotator_sets = c('Seurat', 'SingleR', 'Seurat (cross)', 'SingleR (cross)')
best = hard |> filter(set == 'all', labelling != 'Vote, ties to Seurat') |>
  mutate(labelling = recode(labelling, 'Vote, ties at random' = 'Majority vote')) |>
  group_by(tissue, species) |>
  mutate(rank = rank(-mean_auroc, ties.method = 'min'),
         best_annotator = labelling %in% annotator_sets &
           mean_auroc == max(mean_auroc[labelling %in% annotator_sets])) |>
  ungroup() |>
  mutate(tissue = factor(tissue, levels = tissue_order), labelling = factor(labelling, levels = rev(cell_levels)))
write_csv(best, 'figures/cross_tissue/best_labelling.csv')
p5 = ggplot(best, aes(species, labelling, fill = rank)) +
  geom_tile(colour = 'white', linewidth = 1.2) +
  geom_tile(data = filter(best, best_annotator), aes(colour = 'Best of the four annotators'),
            fill = NA, linewidth = 1.1, width = 0.92, height = 0.88) +
  geom_text(aes(label = sprintf('%.3f', mean_auroc), colour = rank <= 3), size = 3.3, show.legend = FALSE) +
  facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(tissue_name)) +
  scale_fill_gradient(low = '#1F4E79', high = '#EEF2F6', breaks = 1:7, name = 'Rank within dataset\n(1 = highest AUROC)',
                      guide = guide_colourbar(reverse = TRUE, barheight = unit(8, 'lines'))) +
  scale_colour_manual(values = c('Best of the four annotators' = '#E8A33D', 'TRUE' = 'white', 'FALSE' = 'grey20'),
                      breaks = 'Best of the four annotators', name = NULL) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       caption = paste('Mean marker AUROC on all cells of the evaluation half (not used to choose f); values as in auroc_by_cell_set.png panel B.',
                       'Majority vote: ties broken at random (seed 3). SingleR scored on the cells it labels (95-100%).',
                       'Not distinguishable (paired bootstrap intervals): colon human SingleR vs SingleR (cross); fetal liver mouse Seurat (cross) vs EBCC (difference 0.002, 95% CI -0.029 to 0.030).',
                       sep = '\n')) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = 'bold', hjust = 0),
        axis.text.y = element_text(face = ifelse(rev(cell_levels) == 'EBCC (chosen f)', 'bold', 'plain')),
        legend.position = 'right', plot.caption = element_text(hjust = 0, size = 9), panel.spacing = unit(0.8, 'lines'))
save_plot(p5, 'figures/cross_tissue/best_labelling.png', width = 11, height = 5)
ggsave('figures/cross_tissue/best_labelling.pdf', p5, width = 11, height = 5)
message('done figures')
