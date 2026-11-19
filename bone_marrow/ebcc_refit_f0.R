# bone_marrow EBCC at f = 0 (not cached by section 3b). With frac = 0 every trusted
# diagonal is 4 + 0 = 4, so make_beta(M) in the 02 scripts equals make_beta_prior().
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'bone_marrow'))
suppressMessages({source('../library.R'); source('../source.R')})
ann_cols = c('celltype_seurat', 'celltype_singleR', 'celltype_seurat_cross', 'celltype_singleR_cross')
dir.create('data/02_annotated_data/ebcc_frac_fits', showWarnings = FALSE)
for (sp in c('human', 'mouse')) {
  M = c(human = 5, mouse = 3)[[sp]]
  meta = readRDS(sprintf('data/02_annotated_data/02_%s_annotated.rds', sp))@meta.data
  classes = sort(unique(na.omit(unlist(lapply(ann_cols, function(a) as.character(meta[[a]]))))))
  t0 = Sys.time()
  fit = run_ebcc(meta, ann_cols, M = M, R = 10, beta = make_beta_prior(classes, ann_cols), n_cores = 5L)
  out = sprintf('data/02_annotated_data/ebcc_frac_fits/%s_frac0.rds', sp)
  saveRDS(list(label = setNames(as.character(fit$label), rownames(meta)), M = M, frac = 0,
               elbo = fit$elbo), out)
  message(sp, ' f = 0 fit saved (', format(Sys.time() - t0), '): ', out)
}
