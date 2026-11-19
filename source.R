# threads for the parallel steps below (FastMNN, Harmony, LIGER, RunUMAP); set N_THREADS to change.
# All of these give identical results to their serial versions.
default_threads = function() as.integer(Sys.getenv('N_THREADS', 16))
options(Seurat.nthreads = default_threads())   # RunUMAP neighbour search

# qc
get_or_na = function(x) if (is.null(x)) NA_character_ else as.character(x)

bg_filter = function(nCount_probe_controls, nCount_transcripts, cell_areas,
                     n_genes, n_total_probes, alpha = 0.01) {
  
  if (sum(nCount_probe_controls) == 0) {
    warning("No negative control counts: background undetectable, all cells kept")
    return(rep(TRUE, length(nCount_transcripts)))
  }
  scale  = n_genes / n_total_probes
  lambda = sum(nCount_probe_controls) / sum(cell_areas) 
  mu = lambda * cell_areas * scale  
  p = ppois(nCount_transcripts - 1, mu, lower.tail = FALSE)
  padj = p.adjust(p, method = "BH")
  
  return(padj < alpha)
}

# function for label transfer using seurat
LabelTransfer = function(data, ref,
                         celltype_column = 'cell_type_regroup',
                         dims = 1:20) {

  # pull raw counts from whatever assay/format the object is in
  get_counts = function(obj) {
    a = DefaultAssay(obj)
    if (inherits(obj[[a]], 'Assay5')) obj[[a]] = JoinLayers(obj[[a]])
    GetAssayData(obj, assay = a, layer = 'counts')
  }

  # fresh single-assay objects: only "RNA", only counts
  q = CreateSeuratObject(counts = get_counts(data)) |>
    NormalizeData()

  r = CreateSeuratObject(counts = get_counts(ref)) |>
    NormalizeData()
  r$celltype = ref[[celltype_column]]

  features = intersect(rownames(q), rownames(r))
  r = r |>
    ScaleData(features = features) |>
    RunPCA(features = features)

  anchors = FindTransferAnchors(reference = r, query = q,
                                features = features,
                                reference.reduction = 'pca',
                                dims = dims)
  pred = TransferData(anchorset = anchors, refdata = r$celltype, dims = dims)

  list(cell_type_transfer = pred$predicted.id,
       score_transfer     = pred$prediction.score.max)
}

# function for label transfer using SingleR
LabelSingleR = function(data, ref, tech = 'Xenium', celltype_column = 'cell_type_regroup', ref_tech = 'RNA',
                        BPPARAM = BiocParallel::SerialParam(),
                        num.threads = default_threads()) {
  # only v5 assays have layers to join (e.g. the Tabula Muris reference is a v3 Assay)
  if (inherits(ref[[ref_tech]], 'Assay5')) ref[[ref_tech]] = JoinLayers(ref[[ref_tech]])
  if (inherits(data[[tech]], 'Assay5')) data[[tech]] = JoinLayers(data[[tech]])
  data_layer = data |> 
    LayerData(layer = 'data', assay = tech)
  ref_layer = ref |> 
    LayerData(layer = 'data', assay = ref_tech)
  singleR = SingleR(test = data_layer,
                    ref = ref_layer,
                    de.method = 't',
                    labels = ref@meta.data[,celltype_column],
                    BPPARAM = BPPARAM,
                    num.threads = num.threads)   # identical results to 1 thread
  return(list(cell_type_singleR = singleR$pruned.labels,
              delta_singleR = singleR$delta.next))
}

# marker dot plot grouped by a label column (NA labels dropped), saved as png
save_dotplot = function(obj, features, group_by, file, dir = 'figures/02_annotation') {
  obj = obj[, !is.na(obj[[group_by]][, 1])]
  # unclipped, smaller strip labels so marker-group names fit over narrow groups
  p = DotPlot(obj, features = features, group.by = group_by) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
          strip.clip = 'off', strip.text = element_text(size = 9),
          panel.spacing.x = unit(0.8, 'lines'))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, file), p, bg = 'white',
         width = 4 + 0.2 * length(unlist(features)),
         height = 2 + 0.4 * length(unique(obj[[group_by]][, 1])))
  p
}

# function to ortholog mapping & concatenate cross-species data
concatenate = function(human_data, mouse_data, homology_tbl, orthology = 'o2o', 
                       tech_human = 'Xenium', tech_mouse = 'Xenium',
                       change_rownames = FALSE) {
  if (change_rownames) {
    
  }
  
  homology_tbl$hsapiens_homolog_associated_gene_name = homology_tbl$external_gene_name
  homology_tbl$hsapiens_homolog_ensembl_gene = homology_tbl$ensembl_gene_id
  homology_tbl$hsapiens_homolog_chromosome = homology_tbl$chromosome_name
  homology_tbl$hsapiens_homolog_chrom_start = homology_tbl$start_position
  homology_tbl$hsapiens_homolog_chrom_end = homology_tbl$end_position
  if (orthology == 'o2o') {
    one2one = homology_tbl |> 
      filter_at(vars(ends_with('homolog_orthology_type')),
                all_vars(. == 'ortholog_one2one')) |> 
      distinct(.data[['hsapiens_homolog_ensembl_gene']],
               .keep_all = TRUE)
  }
  seurats = c(human_data, mouse_data)
  names(seurats) = c('hsapiens', 'mmusculus')
  seurats_o2o = list()
  seurats_o2o[['mmusculus']] = seurats[['mmusculus']][
    (seurats[['mmusculus']] |> 
       rownames() |> 
       tolower()) %in% (one2one$mmusculus_homolog_ensembl_gene |> 
                          tolower()),]
  seurats_o2o[['hsapiens']] = seurats[['hsapiens']][
    tolower(seurats[['hsapiens']] |> rownames()) %in% tolower(one2one$hsapiens_homolog_ensembl_gene),
  ]
  one2one = one2one |> 
    filter(tolower(mmusculus_homolog_ensembl_gene) %in% 
             tolower(seurats_o2o[['mmusculus']] |> rownames())) |> 
    distinct(hsapiens_homolog_ensembl_gene, .keep_all = TRUE)
  seurats_o2o[['mmusculus']] = seurats_o2o[['mmusculus']][
    tolower(seurats_o2o[['mmusculus']] |> rownames()) %in%
      tolower(one2one$mmusculus_homolog_ensembl_gene),
  ]
  one2one = one2one[
    match(tolower(seurats_o2o[['mmusculus']] |> rownames()),
          tolower(one2one$mmusculus_homolog_ensembl_gene)),
  ]
  seurats_o2o[['mmusculus']][[tech_mouse]][['mmusculus_homolog_ensembl_gene']] = seurats_o2o[['mmusculus']] |> rownames()
  hsapiens_id = one2one$hsapiens_homolog_ensembl_gene
  names(hsapiens_id) = one2one$mmusculus_homolog_ensembl_gene
  seurats_o2o[['mmusculus']][[tech_mouse]][['hsapiens_homolog_ensembl_gene']] = hsapiens_id
  rownames(seurats_o2o[['mmusculus']]) = hsapiens_id
  seurats_o2o[['mmusculus']] = seurats_o2o[['mmusculus']][!duplicated(seurats_o2o[['mmusculus']] |> rownames()),]
  genes_intersection = intersect(
    rownames(seurats_o2o[['hsapiens']]),
    rownames(seurats_o2o[['mmusculus']])
  )
  seurats_o2o[['mmusculus']] = seurats_o2o[['mmusculus']][rownames(seurats_o2o[['mmusculus']]) %in% genes_intersection,]
  seurats_o2o[['hsapiens']] = seurats_o2o[['hsapiens']][rownames(seurats_o2o[['hsapiens']]) %in% genes_intersection,]
  seurat_o2o = merge(seurats_o2o[['mmusculus']], seurats_o2o[['hsapiens']],
                      add.cell.ids = names(seurats_o2o))
  o2o_species = lapply(
    seurat_o2o@meta.data |> rownames(),
    FUN = function(x) str_split(x, '_')[[1]][1]
  )
  names(o2o_species) = seurat_o2o@meta.data |> rownames()
  seurat_o2o@meta.data$orig.ident = o2o_species |> unlist()
  return(seurat_o2o)
}

# put a reference onto the query's gene symbols through one-to-one orthologs
# ref_species is the reference's species ('hsapiens' or 'mmusculus'); the query
# needs gene_id/gene_name in its Xenium feature metadata (01_raw_data_process.R)
cross_species_ref = function(ref, query, homology_tbl, ref_species,
                             celltype_column = 'celltype', assay = 'RNA') {
  o2o = homology_tbl |>
    filter(mmusculus_homolog_orthology_type == 'ortholog_one2one')
  if (ref_species == 'hsapiens') {
    map = data.frame(ref_sym = o2o$external_gene_name,
                     query_id = o2o$mmusculus_homolog_ensembl_gene)
  } else {
    map = data.frame(ref_sym = o2o$mmusculus_homolog_associated_gene_name,
                     query_id = o2o$ensembl_gene_id)
  }
  map = map |>
    filter(!is.na(ref_sym), ref_sym != '', !is.na(query_id), query_id != '') |>
    distinct()
  # drop symbols / ids that are not uniquely paired
  map = map |>
    filter(!ref_sym %in% ref_sym[duplicated(ref_sym)],
           !query_id %in% query_id[duplicated(query_id)])
  genes = query[['Xenium']][[]]
  map$query_name = genes$gene_name[match(map$query_id, genes$gene_id)]
  map = map |>
    filter(!is.na(query_name), ref_sym %in% rownames(ref[[assay]]))
  counts = LayerData(ref, assay = assay, layer = 'counts')[map$ref_sym, ]
  data = LayerData(ref, assay = assay, layer = 'data')[map$ref_sym, ]
  rownames(counts) = map$query_name
  rownames(data) = map$query_name
  out = CreateSeuratObject(counts = counts, assay = assay,
                           meta.data = ref@meta.data[, celltype_column, drop = FALSE])
  LayerData(out, assay = assay, layer = 'data') = data
  out
}

# one marker dot plot per annotation column, saved as png
# file_names: the part after '<prefix>_dotplot_' for each label column (default: the column names)
save_marker_dotplots = function(obj, markers, label_cols, prefix, out_dir = 'plots',
                                file_names = label_cols) {
  stopifnot(length(file_names) == length(label_cols))
  file_names = setNames(file_names, label_cols)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # keep genes on the panel, each gene once (first group wins)
  present = lapply(markers, intersect, rownames(obj))
  for (k in seq_along(present)[-1]) {
    present[[k]] = setdiff(present[[k]], unlist(present[seq_len(k - 1)]))
  }
  present = present[lengths(present) > 0]
  missing = setdiff(unlist(markers), rownames(obj))
  if (length(missing) > 0) message('markers not on panel: ', paste(missing, collapse = ', '))
  # light copy: DotPlot only needs the normalised data and the labels
  lite = CreateSeuratObject(counts = LayerData(obj, layer = 'counts'),
                            meta.data = obj@meta.data[, label_cols, drop = FALSE])
  LayerData(lite, layer = 'data') = LayerData(obj, layer = 'data')
  n_genes = length(unlist(present))
  for (col in label_cols) {
    labs = as.character(lite@meta.data[[col]])
    # rows follow the marker order, top to bottom
    lite@meta.data[[col]] = factor(labs, levels = rev(c(
      intersect(names(markers), labs),
      sort(setdiff(unique(na.omit(labs)), names(markers))))))
    p = DotPlot(lite[, !is.na(labs)], features = present, group.by = col, dot.scale = 4) +
      RotatedAxis() +
      theme(axis.text.x = element_text(size = 7),
            strip.text.x = element_text(size = 7, angle = 90),
            panel.spacing = unit(1, 'pt')) +
      labs(x = NULL, y = NULL, title = col)
    ggsave(file.path(out_dir, sprintf('%s_dotplot_%s.png', prefix, file_names[[col]])), p,
           width = max(8, 0.16 * n_genes + 4), height = 6, dpi = 300, limitsize = FALSE)
  }
}

# fastmnnn
## helper function
l2 = function(M) {
  n = sqrt(rowSums(M^2))
  n[n == 0] = 1                 
  M / n
}

add_reduction = function(data, emb, name, key, assay, loadings) {
  emb = as.matrix(emb)
  missing = setdiff(colnames(data), rownames(emb))
  if (length(missing) > 0)
    stop(length(missing), " cells in `data` not found in embedding - check barcode names")
  emb = emb[colnames(data),]        
  colnames(emb) = paste0(key, seq_len(ncol(emb)))
  data[[name]] = CreateDimReducObject(embeddings = emb, key = key, 
                                      assay = assay, loadings = loadings)
  data
}


integrate_fastmnn = function(data, batch_key = "orig.ident",
                             assay = NULL, max_features = 2000,
                             H = NULL, W = NULL, k = 20,
                             seed = 1, half_dims = TRUE,
                             BPPARAM = BiocParallel::MulticoreParam(default_threads())) {
  set.seed(seed)
  if (is.null(assay)) assay = DefaultAssay(data)
  start = Sys.time()

  if (is.null(data@misc$runtime_sec)) {
    data@misc$runtime_sec = list()
  }

    # ---- 1. FastMNN on UINMF loadings ----
  if (!is.null(H)) {
    # Keep cells that are in data AND have a usable (finite, non-zero) H row
    keep_rows = function(h) {
      h = h[intersect(rownames(h), colnames(data)), , drop = FALSE]
      h[complete.cases(h) & rowSums(abs(h)) > 0, , drop = FALSE]
    }
    H$hsapiens  = keep_rows(H$hsapiens)
    H$mmusculus = keep_rows(H$mmusculus)

    keep = c(rownames(H$hsapiens), rownames(H$mmusculus))
    message(sprintf("FastMNN UINMF: keeping %d / %d cells (%d dropped: missing or all-zero H)",
                    length(keep), ncol(data), ncol(data) - length(keep)))
    data = subset(data, cells = keep)

    # Use all UINMF factors (both species share the same factors)
    stopifnot(ncol(H$hsapiens) == ncol(H$mmusculus))
    if (ncol(H$hsapiens) < 2L) stop("UINMF reduction has too few dimensions.")

    out = batchelor::reducedMNN(
      as.matrix(H$hsapiens),
      as.matrix(H$mmusculus),
      k = k,
      BPPARAM = BPPARAM
    )

    # Put corrected rows back in the same order as the cells in data
    corrected = out$corrected
    rownames(corrected) = c(rownames(H$hsapiens), rownames(H$mmusculus))
    corrected = corrected[colnames(data), , drop = FALSE]
    stopifnot(all(is.finite(corrected)))

    data = add_reduction(
      data, corrected, "mnn_uinmf", "mnnuinmf_",
      assay, loadings = W
    )

    elapsed = as.numeric(difftime(Sys.time(), start, units = "secs"))
    data@misc$runtime_sec$fastmnn_uinmf = elapsed
    message(sprintf("FastMNN UINMF elapsed: %.1f sec", elapsed))
    return(data)
  }

  # ---- 2. PCA per batch, then FastMNN on low-dimensional embeddings ----
  object.list = SplitObject(data, split.by = batch_key)

  features = SelectIntegrationFeatures(
    object.list = object.list,
    nfeatures = min(max_features, nrow(data[[assay]]))
  )

  features = Reduce(
    intersect,
    lapply(object.list, rownames)
  ) |> intersect(features)

  if (length(features) < 2L) {
    stop("Too few shared features for FastMNN.")
  }

  # log-normalised expression per batch (genes x cells)
  mats = lapply(object.list, function(x) {
    x = NormalizeData(x, assay = assay, verbose = FALSE)
    x[[assay]] = JoinLayers(x[[assay]])
    GetAssayData(x, assay = assay, layer = "data")[features, , drop = FALSE]
  })

  n_dims = min(20L, length(features) - 1L)
  # half_dims = FALSE: use every available dimension (up to 20)
  d = min(20L, if (half_dims) floor(n_dims / 2) else n_dims)
  if (d < 1L) stop("PCA reductions have too few dimensions.")

  # MNN needs every batch in the SAME low-dimensional space:
  # multiBatchPCA computes one shared PCA (batches weighted equally),
  # instead of a separate PCA per batch
  pcs = as.list(do.call(batchelor::multiBatchPCA, c(mats, list(d = d, BPPARAM = BPPARAM))))
  for (b in seq_along(pcs)) rownames(pcs[[b]]) = colnames(mats[[b]])

  out = do.call(
    batchelor::reducedMNN,
    c(pcs, list(k = k, BPPARAM = BPPARAM))
  )
  rownames(out$corrected) = unlist(lapply(mats, colnames))

  data = add_reduction(
    data, out$corrected, "mnn", "mnn_",
    assay, loadings = matrix(numeric(0), nrow = 0, ncol = 0)
  )

  elapsed = as.numeric(difftime(Sys.time(), start, units = "secs"))
  data@misc$runtime_sec$fastmnn = elapsed
  message(sprintf(
    "FastMNN elapsed: %.1f sec (%.1f min)",
    elapsed, elapsed / 60
  ))

  data
}


# cca

integrate_cca = function(data, batch_key = "orig.ident",
                         H = NULL, W = NULL, seed = 1,
                         k.filter = NA,      # NA = skip anchor filtering
                         k.weight = 100, half_dims = TRUE) {
  set.seed(seed)
  start = Sys.time()
  assay = DefaultAssay(data)

  # Retry with a smaller k.weight if Seurat says there aren't enough anchor cells
  run_with_kweight = function(fun, k.weight) {
    tryCatch(
      fun(k.weight),
      error = function(e) {
        msg = conditionMessage(e)
        if (grepl("less than k.weight", msg)) {
          max_k = as.integer(sub(".*less than (\\d+).*", "\\1", msg))
          new_k = max(5L, max_k - 1L)
          message(sprintf("Too few anchor cells; retrying with k.weight = %d", new_k))
          fun(new_k)
        } else stop(e)
      }
    )
  }

  # Check the batch column has >= 2 groups
  check_batches = function(data, where) {
    b = droplevels(factor(data[[batch_key, drop = TRUE]]))
    message(where, ": ", paste(names(table(b)), table(b), sep = "=", collapse = ", "))
    if (nlevels(b) < 2L || anyNA(b))
      stop(where, ": need >= 2 batches and no NA in '", batch_key, "'")
    b
  }

  if (!is.null(H)) {
    H = rbind(H$hsapiens, H$mmusculus)

    # keep cells in data with a finite, non-zero H row
    ok = rownames(H)[complete.cases(H) & rowSums(abs(H)) > 0]
    common = intersect(colnames(data), ok)
    message(sprintf("CCA UINMF: keeping %d / %d cells", length(common), ncol(data)))

    data = subset(data, cells = common)
    H = as.matrix(H[colnames(data), , drop = FALSE])

    data[[batch_key]] = check_batches(data, "After subset")

    # split layers by batch (join first in case they were already split)
    data[[assay]] = JoinLayers(data[[assay]])
    data[[assay]] = split(data[[assay]], f = data[[batch_key, drop = TRUE]])
    message("Layers: ", paste(Layers(data[[assay]]), collapse = ", "))

    data = add_reduction(data, H, name = "uinmf_H", key = "uinmfraw_",
                         assay = assay, loadings = W)

    n_dims = ncol(Embeddings(data, "uinmf_H"))
    d = min(20L, if (half_dims) floor(n_dims / 2) else n_dims)
    if (d < 1L) stop("UINMF reduction has too few dimensions.")

    data = run_with_kweight(function(kw) {
      IntegrateLayers(
        object = data,
        method = CCAIntegration,
        orig.reduction = "uinmf_H",
        new.reduction = "cca_uinmf",
        assay = assay,
        dims = seq_len(d),
        k.filter = k.filter,
        k.weight = kw
      )
    }, k.weight)

    data = JoinLayers(data)

  } else {
    data[[batch_key]] = check_batches(data, "Before split")

    object.list = SplitObject(data, split.by = batch_key)
    features = SelectIntegrationFeatures(object.list = object.list,
                                         nfeatures = min(2000L, nrow(data)))
    features = Reduce(intersect, lapply(object.list, rownames)) |>
      intersect(features)
    if (length(features) < 2L) stop("Too few shared features for CCA integration.")

    d = min(20L, length(features) - 1L)

    anchors = FindIntegrationAnchors(
      object.list = object.list,
      anchor.features = features,
      reduction = "cca",
      dims = seq_len(d),
      k.filter = k.filter
    )

    # Set k.weight from the anchor cells actually available
    a = anchors@anchors
    n_anchor_cells = min(length(unique(a$cell1)), length(unique(a$cell2)))
    kw = min(k.weight, max(5L, n_anchor_cells - 1L))
    message(sprintf("Anchors: %d, anchor cells: %d, k.weight: %d",
                    nrow(a), n_anchor_cells, kw))

    data = run_with_kweight(function(kw) {
      IntegrateData(anchorset = anchors, dims = seq_len(d), k.weight = kw)
    }, kw)

    # IntegrateData returns an assay, not a reduction: build one from it
    DefaultAssay(data) = "integrated"
    int_feats = rownames(data[["integrated"]])
    VariableFeatures(data) = int_feats
    n_pcs = max(1L, min(20L, length(int_feats) - 1L))
    data = data |>
      ScaleData(features = int_feats, verbose = FALSE) |>
      RunPCA(features = int_feats, npcs = n_pcs,
             reduction.name = "cca", reduction.key = "cca_", verbose = FALSE)
  }

  elapsed = as.numeric(difftime(Sys.time(), start, units = "secs"))
  message(sprintf("CCA elapsed: %.1f sec (%.1f min)", elapsed, elapsed / 60))
  data@misc$runtime_sec = elapsed
  data
}

# rpca

integrate_rpca = function(data, batch_key = "orig.ident",
                          H = NULL, W = NULL, seed = 1, half_dims = TRUE) {
  start = Sys.time()
  set.seed(seed)

  assay = DefaultAssay(data)

  if (!is.null(H)) {
    # Combine UINMF embeddings
    H = rbind(H$hsapiens, H$mmusculus)

    # keep cells in data with a finite, non-zero H row
    ok = rownames(H)[complete.cases(H) & rowSums(abs(H)) > 0]
    common = intersect(colnames(data), ok)
    message(sprintf("RPCA UINMF: keeping %d / %d cells", length(common), ncol(data)))

    data = subset(data, cells = common)
    H = H[common, , drop = FALSE]

    data = Split_Layers(
      data,
      assay = assay,
      split.by = batch_key
    )

    data = add_reduction(
      data, H,
      name = "uinmf_H",
      key = "uinmfraw_",
      assay = assay,
      loadings = W
    )

    # Use only dimensions actually available
    n_dims = ncol(Embeddings(data, "uinmf_H"))
    d = min(20L, if (half_dims) floor(n_dims / 2) else n_dims)
    # RPCAIntegration also runs a PCA per batch on the shared genes, which has at most
    # (genes - 1) components: cap d there (only binds with very few genes, e.g. 9)
    d = min(d, nrow(data[[assay]]) - 1L)

    if (d < 1L) {
      stop("UINMF reduction has too few dimensions.")
    }

    data = IntegrateLayers(
      object = data,
      method = RPCAIntegration,
      orig.reduction = "uinmf_H",
      new.reduction = "rpca_uinmf",
      assay = assay,
      dims = seq_len(d)
    )

    data = JoinLayers(data)

  } else {
    object.list = SplitObject(data, split.by = batch_key)
    features = rownames(data)

    object.list = lapply(object.list, function(x) {
      RunPCA(
        x,
        features = features,
        npcs = min(20L, length(features) - 1L),
        verbose = FALSE,
        reduction.name = "rpca"
      )
    })

    # Use the number of dimensions available in every object
    n_dims = vapply(
      object.list,
      function(x) ncol(Embeddings(x, "rpca")),
      integer(1)
    )

    d = min(20L, if (half_dims) floor(min(n_dims) / 2) else min(n_dims))

    if (d < 1L) {
      stop("PCA reductions have too few dimensions.")
    }

  anchors = FindIntegrationAnchors(
    object.list = object.list,
    anchor.features = features,
    reduction = "rpca",
    dims = seq_len(d),
    scale = FALSE
  )

  data = IntegrateData(
    anchorset = anchors,
    dims = seq_len(d)
  )

  # IntegrateData returns an assay, not a reduction: the integrated embedding
  # is a PCA of the "integrated" assay (the stacked per-species PCAs are
  # in different spaces and are NOT an integration)
  DefaultAssay(data) = "integrated"
  int_feats = rownames(data[["integrated"]])
  VariableFeatures(data) = int_feats
  data = data |>
    ScaleData(features = int_feats, verbose = FALSE) |>
    RunPCA(features = int_feats, npcs = max(2L, min(20L, length(int_feats) - 1L)),
           reduction.name = "rpca", reduction.key = "RPCA_", verbose = FALSE)
  }

  elapsed = as.numeric(
    difftime(Sys.time(), start, units = "secs")
  )

  message(sprintf(
    "RPCA elapsed: %.1f sec (%.1f min)",
    elapsed, elapsed / 60
  ))

  data@misc$runtime_sec = elapsed
  return(data)
}

# uinmf
integrate_harmony = function(data, batch_key = 'orig.ident', H = NULL, W = NULL, seed = 1,
                             max_pcs = 50L, ncores = default_threads()) {
  set.seed(seed)
  start = Sys.time()

  if (!is.null(H)) {
    H_raw = rbind(H$mmusculus, H$hsapiens)

    # keep cells in data with a finite, non-zero H row
    ok = rownames(H_raw)[complete.cases(H_raw) & rowSums(abs(H_raw)) > 0]
    common = intersect(colnames(data), ok)
    message(sprintf("Harmony UINMF: keeping %d / %d cells", length(common), ncol(data)))

    data = subset(data, cells = common)
    H_raw = as.matrix(H_raw[colnames(data), , drop = FALSE])

    # explicit, aligned batch vector
    batch = droplevels(factor(data@meta.data[colnames(data), batch_key]))
    print(table(batch, useNA = "ifany"))
    if (nlevels(batch) < 2L || anyNA(batch)) {
      stop("Harmony needs >= 2 batches and no NA in '", batch_key, "'")
    }
    data@meta.data[[batch_key]] = batch

    data = add_reduction(data, H_raw, name = 'uinmf_H', key = 'uinmfraw_',
                         assay = DefaultAssay(data), loadings = W)

    # call Harmony on the matrix directly (cells x factors)
    emb = harmony::RunHarmony(
      data_mat  = H_raw,
      meta_data = data.frame(batch = batch, row.names = colnames(data)),
      vars_use  = "batch",
      ncores    = ncores,
      verbose   = TRUE
    )
    rownames(emb) = colnames(data)
    colnames(emb) = paste0("harmony_", seq_len(ncol(emb)))

    data[['harmony']] = CreateDimReducObject(
      embeddings = emb, key = "harmony_", assay = DefaultAssay(data)
    )

    elapsed = as.numeric(difftime(Sys.time(), start, units = "secs"))
    message(sprintf("Harmony elapsed: %.1f sec (%.1f min)", elapsed, elapsed/60))
    data@misc$runtime_sec = elapsed
    return(data)
  }

  # gene-based branch: cap PCs by number of genes (you only have 2-8)
  n_pcs = max(2L, min(max_pcs, nrow(data) - 1L))
  data = data |>
    RunPCA(features = rownames(data), npcs = n_pcs, verbose = FALSE) |>
    RunHarmony(group.by.vars = batch_key, ncores = ncores, verbose = TRUE)

  elapsed = as.numeric(difftime(Sys.time(), start, units = "secs"))
  message(sprintf("Harmony elapsed: %.1f sec (%.1f min)", elapsed, elapsed/60))
  data@misc$runtime_sec = elapsed
  data
}

# liger
integrate_liger = function(data, batch_key = 'orig.ident',
                          k = 20, n_cores = default_threads()) {
  start = Sys.time()
  
  data = data |>
    ScaleData(rownames(data), 
              split.by = batch_key, do.center = FALSE)
  
  data = data |>
    runINMF(k = k, datasetVar = batch_key, layer = 'scale.data', nCores = n_cores) |>
    quantileNorm(split.by = batch_key)
  
  elapsed = as.numeric(difftime(Sys.time(), start, units = "secs"))
  message(sprintf("LIGER elapsed: %.1f sec (%.1f min)", elapsed, elapsed/60))
  data@misc$runtime_sec = elapsed
  return(data)
}

# liger uinmf
integrate_liger_uinmf = function(human, mouse, k = 20, lambda = 5, ref = 'hsapiens',
                                 add_prefix = 'auto', n_cores = default_threads()) {
  start = Sys.time()
  obj_list = list(hsapiens = human, mmusculus = mouse)

  # add_prefix = 'auto' names cells '<dataset>_<barcode>', matching concatenate();
  # rliger's addPrefix = TRUE does NOT add the prefix (only "auto" does)
  lig = createLiger(
    sapply(obj_list, function(x) LayerData(x, layer = "counts"), simplify = FALSE),
    addPrefix = add_prefix
  )
  
  lig = lig |>
    normalize() |>
    selectGenes(useUnsharedDatasets = c("hsapiens", "mmusculus")) |>
    scaleNotCenter() |>
    runUINMF(k = k, lambda = lambda, nIteration = 30, nCores = n_cores) |>
    quantileNorm(reference = ref) |>
    runCluster() |>
    runUMAP()
  
  obj = ligerToSeurat(lig)
  if ("UMAP" %in% names(obj@reductions)) {
    obj[["umap"]] = obj[["UMAP"]]
    obj[["UMAP"]] = NULL
  }
  
  elapsed = as.numeric(difftime(Sys.time(), start, units = "secs"))
  message(sprintf("LIGER UINMF elapsed: %.1f sec (%.1f min)", elapsed, elapsed/60))
  list(seurat = obj, liger = lig)
}

h5closeAll()   # the 'open HDF5 file handle' warning: release earlier handles first

fix_nullable_strings = function(file) {
  ls = h5ls(file, recursive = TRUE)
  paths = ifelse(ls$group == "/", paste0("/", ls$name), paste0(ls$group, "/", ls$name))
  grps = paths[ls$otype == "H5I_GROUP"]
  fixed = character()

  for (g in grps) {
    a = h5readAttributes(file, g)
    if (!identical(a[["encoding-type"]], "nullable-string-array")) next

    vals = as.character(h5read(file, paste0(g, "/values")))
    mask = as.logical(as.character(h5read(file, paste0(g, "/mask"))))
    vals[mask %in% TRUE] = ""          # string-array has no NA; masked entries become ""

    h5delete(file, g)
    h5write(vals, file, g, variableLengthString = TRUE)

    fid = H5Fopen(file); did = H5Dopen(fid, g)
    h5writeAttribute("string-array", did, "encoding-type",    asScalar = TRUE)
    h5writeAttribute("0.2.0",        did, "encoding-version", asScalar = TRUE)
    H5Dclose(did); H5Fclose(fid)
    fixed = c(fixed, g)
  }
  h5closeAll()
  fixed
}

attach_raw_meta <- function(obj, raw) {
  if (colnames(obj)[1] == 1) {
    colnames(obj) = colnames(raw)
    rownames(obj) = rownames(raw)
    obj@meta.data <- raw@meta.data
    obj}
  obj = obj[,colnames(raw)]
  obj@meta.data = raw@meta.data
  obj
}

# ==============================================================================
# EBCC: Enhanced Bayesian Classifier Combination
# Li, Rubinstein & Cohn (2019), ICML. https://proceedings.mlr.press/v97/li19i/li19i.pdf
#
# Each true class k has M latent subtypes; annotators are conditionally independent
# given (class, subtype) rather than given class alone, so annotators whose accuracy
# varies together across subtypes (e.g. two methods sharing a reference) may be
# correlated within a class. Inference: mean-field variational Bayes with R random
# restarts, keeping the restart with the highest ELBO.
#
# Requires: matrixStats (rowMaxs)
# ==============================================================================

library(matrixStats)

# ---- input helpers ----

# optional relabelling onto a shared vocabulary
harmonise = function(x, map = NULL) {
  x = as.character(x)
  if (!is.null(map)) x = ifelse(x %in% names(map), map[x], x)
  x
}

# NA-safe indicator: a missing annotation is an all-zero row (marginalised out)
indicator = function(vec, classes) {
  m = matrix(0, length(vec), length(classes), dimnames = list(NULL, classes))
  ok = !is.na(vec)
  m[cbind(which(ok), match(vec[ok], classes))] = 1
  m
}

# ---- Dirichlet helpers (row-wise) ----

log_beta_rows    = function(A) rowSums(lgamma(A)) - lgamma(rowSums(A))
elog_rows        = function(A) digamma(A) - digamma(rowSums(A))   # E[log x] under Dir
dir_entropy_rows = function(A) {
  a0 = rowSums(A)
  log_beta_rows(A) + (a0 - ncol(A)) * digamma(a0) - rowSums((A - 1) * digamma(A))
}

# ---- annotator- and cell-type-specific priors ----

#' Default confusion-matrix priors: one K x K matrix per annotator
#' (rows = true class, columns = label given), beta_diag on the diagonal, beta_off elsewhere.
#' Edit the returned matrices to encode prior knowledge about specific annotators, e.g.
#'   beta$celltype_singleR["Granulocyte", "Granulocyte"] = 500   # trust SingleR on granulocytes
#'   beta$celltype_seurat["Granulocyte", "Monocyte"]     = 500   # Seurat often calls granulocytes Monocyte
#' Values are pseudo-counts per subtype, so a class's total prior weight is M times the value;
#' compare it with the number of cells in that class to judge how strong it is.
make_beta_prior = function(classes, ann_cols, beta_diag = 4, beta_off = 1) {
  b = matrix(beta_off, length(classes), length(classes), dimnames = list(classes, classes))
  diag(b) = beta_diag
  setNames(rep(list(b), length(ann_cols)), ann_cols)
}

# ---- model fit ----

#' Fit EBCC to several categorical annotations of the same cells.
#'
#' @param meta       data frame with one row per cell (e.g. seurat@meta.data)
#' @param ann_cols   names of the annotation columns in `meta`
#' @param M          subtypes per class (paper default 10; effective number is usually smaller)
#' @param a_pi       Dirichlet prior on subtype weights (paper: 0.1 / M, encourages few subtypes)
#' @param beta_diag,beta_off  Dirichlet prior on each confusion row (paper: 4 diagonal, 1 off),
#'                   i.e. "annotators are better than random"; used for any annotator not in `beta`
#' @param beta       optional named list of K x K prior matrices, one per annotator
#'                   (see make_beta_prior); overrides beta_diag/beta_off for those annotators
#' @param alpha      "empirical": class-proportion prior = majority-vote counts (paper's EBCC-Emp);
#'                   "flat": Dir(1)
#' @param mask       restrict each cell to classes at least one annotator gave
#'                   (restricts the variational family, so the ELBO stays a valid bound)
#' @param label_map  optional named vector mapping labels onto a shared vocabulary
#' @param R          random restarts; the solution with the highest ELBO is kept
#' @param maxit,tol  iteration cap and relative ELBO tolerance per restart
#' @param keep_rho   also return the cell x (class, subtype) posterior (large: cells x K*M)
#' @return list with label, posterior (cells x classes), par (tau, pi, v, conf,
#'         subtype_share), elbo, elbo_restarts, iterations, classes, obs_idx, M,
#'         and rho if keep_rho = TRUE
run_ebcc = function(meta, ann_cols, M = 10, a_pi = 0.1 / M,
                    beta_diag = 4, beta_off = 1, beta = NULL,
                    alpha = c("empirical", "flat"), mask = TRUE,
                    label_map = NULL, R = 10, maxit = 2000, tol = 1e-8,
                    seed = 1, keep_rho = FALSE, verbose = TRUE, n_cores = 1L) {
  alpha = match.arg(alpha)
  X = lapply(ann_cols, function(col) harmonise(meta[[col]], label_map))
  names(X) = ann_cols
  J = length(X); I = nrow(meta)
  classes = sort(unique(na.omit(unlist(X))))
  K = length(classes); KM = K * M
  Z = lapply(X, indicator, classes = classes)
  obs_idx = lapply(X, match, table = classes)

  # subtype columns are ordered (k1 m1, ..., k1 mM, k2 m1, ...)
  cls_of = rep(seq_len(K), each = M)
  B = matrix(0, KM, K); B[cbind(seq_len(KM), cls_of)] = 1   # sums subtypes into their class

  votes = Reduce(`+`, Z); n_votes = rowSums(votes)
  support = votes > 0; support[n_votes == 0, ] = TRUE
  support_km = support[, cls_of]

  # majority-vote initialisation of q(z_i)
  gamma0 = votes; gamma0[n_votes == 0, ] = 1
  gamma0 = gamma0 / rowSums(gamma0)

  alpha_vec = if (alpha == "empirical") colSums(gamma0) else rep(1, K)
  # confusion-matrix prior per annotator (K x K, rows = true class, cols = label)
  Beta_default = make_beta_prior(classes, "default", beta_diag, beta_off)[[1]]
  Beta_list = lapply(ann_cols, function(a) {
    if (is.null(beta[[a]])) return(Beta_default)
    b = beta[[a]]
    if (!all(classes %in% rownames(b)) || !all(classes %in% colnames(b)))
      stop("prior matrix for ", a, " must have row and column names covering all classes")
    b[classes, classes]
  })
  names(Beta_list) = ann_cols
  Beta_km_list = lapply(Beta_list, function(b) b[cls_of, ])  # KM x K per annotator

  prior_const = -log_beta_rows(matrix(alpha_vec, 1)) -
    K * log_beta_rows(matrix(a_pi, 1, M)) -
    M * sum(vapply(Beta_list, function(b) sum(log_beta_rows(b)), numeric(1)))

  # parameter step: variational Dirichlet parameters from rho = q(z_i = k, g_i = m)
  update_global = function(rho) {
    cs = colSums(rho)
    list(nu  = alpha_vec + as.vector(cs %*% B),                       # class proportions
         eta = matrix(a_pi + cs, K, M, byrow = TRUE),                  # subtype weights, K x M
         mu  = Map(function(Zj, Bkm) Bkm + crossprod(rho, Zj),         # per annotator, KM x K
                   Z, Beta_km_list))
  }

  # cell step: rho_ikm ∝ exp(E log tau_k + E log pi_km + sum_j E log v_{j,km,y_ij})
  update_local = function(g) {
    Elog_tau = elog_rows(matrix(g$nu, 1))[1, ]
    Elog_pi  = elog_rows(g$eta)
    base = rep(Elog_tau, each = M) + as.vector(t(Elog_pi))
    logR = matrix(base, I, KM, byrow = TRUE)
    for (j in seq_len(J)) logR = logR + Z[[j]] %*% t(elog_rows(g$mu[[j]]))
    if (mask) logR[!support_km] = -Inf
    rho = exp(logR - rowMaxs(logR))
    rho / rowSums(rho)
  }

  # ELBO as in Section 4.4 of the paper (valid when g has just been updated from rho)
  elbo = function(rho, g) {
    Elog_tau = elog_rows(matrix(g$nu, 1))[1, ]
    Elog_pi  = elog_rows(g$eta)
    nz = rho > 0
    sum((g$nu - 1) * Elog_tau) +
      sum((g$eta - 1) * Elog_pi) +
      sum(vapply(g$mu, function(mu) sum((mu - 1) * elog_rows(mu)), numeric(1))) +
      prior_const +
      dir_entropy_rows(matrix(g$nu, 1)) +
      sum(dir_entropy_rows(g$eta)) +
      sum(vapply(g$mu, function(mu) sum(dir_entropy_rows(mu)), numeric(1))) -
      sum(rho[nz] * log(rho[nz]))
  }

  # one restart; replaying the earlier restarts' draws gives restart r the same
  # initialisation it gets in the sequential loop, so n_cores does not change results
  fit_restart = function(r) {
    set.seed(seed)
    for (s in seq_len(r - 1)) invisible(rexp(I * KM))
    # paper's initialisation: MV class probabilities times a Dir(1_M) draw per cell and class
    w   = matrix(rexp(I * KM), I, KM)
    rho = gamma0[, cls_of] * w / (w %*% B)[, cls_of]
    rm(w)

    msgs = character(0)
    prev = -Inf
    for (it in seq_len(maxit)) {
      g = update_global(rho)
      L = elbo(rho, g)
      if (L < prev - 1e-6 * abs(prev))
        msgs = c(msgs, sprintf("restart %d: ELBO decreased at iteration %d", r, it))
      if (abs(L - prev) < tol * abs(L)) break
      rho  = update_local(g)
      prev = L
    }
    if (it == maxit) msgs = c(msgs, sprintf("restart %d did not converge within %d iterations", r, maxit))
    if (verbose) cat(sprintf("restart %2d | iterations %3d | ELBO %.2f\n", r, it, L))
    list(rho = rho, g = g, elbo = L, iterations = it, msgs = msgs)
  }

  if (n_cores > 1L) {
    fits = parallel::mclapply(seq_len(R), fit_restart, mc.cores = min(n_cores, R),
                              mc.preschedule = FALSE)
    bad = !vapply(fits, function(f) is.list(f) && !is.null(f$elbo), logical(1))
    if (any(bad)) stop("EBCC restarts failed: ", paste(which(bad), collapse = ", "), "\n",
                       paste(unlist(fits[bad]), collapse = "\n"))
  } else {
    fits = vector("list", R)
    for (r in seq_len(R)) {
      fits[[r]] = fit_restart(r)
      if (r > 1 && fits[[r]]$elbo <= max(vapply(fits[seq_len(r - 1)], `[[`, numeric(1), "elbo")))
        fits[[r]]$rho = NULL   # keep only the best restart's cells x subtypes matrix in memory
    }
  }
  for (f in fits) for (m in f$msgs) warning(m)
  elbos = vapply(fits, `[[`, numeric(1), "elbo")
  best = fits[[which.max(elbos)]]
  rm(fits)

  g = best$g
  gamma = best$rho %*% B                                   # q(z_i = k), subtypes summed out
  dimnames(gamma) = list(rownames(meta), classes)
  labels = max.col(gamma, ties.method = "first")

  sub_names = paste(classes[cls_of], rep(seq_len(M), K), sep = "#")
  tau = setNames(g$nu / sum(g$nu), classes)
  pi  = g$eta / rowSums(g$eta); dimnames(pi) = list(classes, paste0("sub", seq_len(M)))
  v   = lapply(g$mu, function(mu) {
    p = mu / rowSums(mu); dimnames(p) = list(sub_names, classes); p })
  # class-level confusion matrices (subtypes averaged by their weights), comparable to DS
  conf = lapply(v, function(vj) {
    cf = t(B * as.vector(t(pi))) %*% vj; dimnames(cf) = list(classes, classes); cf })
  share = g$eta - a_pi; share = share / rowSums(share)    # posterior subtype usage per class

  out = list(label = factor(classes[labels], levels = classes),
             posterior = gamma,
             par = list(tau = tau, pi = pi, v = v, conf = conf, subtype_share = share),
             elbo = best$elbo, elbo_restarts = elbos, iterations = best$iterations,
             classes = classes, obs_idx = obs_idx, M = M)
  if (keep_rho) {
    dimnames(best$rho) = list(rownames(meta), sub_names)
    out$rho = best$rho
  }
  out
}

# ---- posterior predictive check ----

# simulate annotator labels from fitted EBCC parameters, copying the observed missingness
simulate_ebcc = function(par, I, observed = NULL) {
  K = length(par$tau); M = ncol(par$pi)
  z = sample.int(K, I, replace = TRUE, prob = par$tau)
  g = integer(I)
  for (k in unique(z)) {
    idx = which(z == k)
    g[idx] = sample.int(M, length(idx), replace = TRUE, prob = par$pi[k, ])
  }
  row = (z - 1) * M + g
  lapply(seq_along(par$v), function(j) {
    cs = t(apply(par$v[[j]], 1, cumsum))
    l  = pmin(rowSums(runif(I) > cs[row, ]) + 1L, K)
    if (!is.null(observed)) l[is.na(observed[[j]])] = NA
    l
  })
}

# pairwise agreement between annotators (proportion of cells with the same label)
pair_agree = function(L) {
  J = length(L); m = matrix(NA_real_, J, J)
  for (a in seq_len(J)) for (b in seq_len(J)) m[a, b] = mean(L[[a]] == L[[b]], na.rm = TRUE)
  m
}

#' Compare observed pairwise agreement with agreement simulated from the fitted model.
#' Positive excess = a pair agrees more than the model explains (unmodelled dependence).
#' @return list of matrices: observed, expected, excess, z
ppc_ebcc = function(fit, nsim = 100, seed = 1) {
  set.seed(seed)
  I   = length(fit$obs_idx[[1]])
  obs = pair_agree(fit$obs_idx)
  sim = replicate(nsim, pair_agree(simulate_ebcc(fit$par, I, fit$obs_idx)), simplify = "array")
  m = apply(sim, 1:2, mean); s = apply(sim, 1:2, sd)
  nm = names(fit$obs_idx)
  lapply(list(observed = obs, expected = m, excess = obs - m, z = (obs - m) / s),
         function(x) { dimnames(x) = list(nm, nm); x })
}

fix_scale_data = function(obj) {
  for (a in Assays(obj)) {
    x = obj[[a]]
    if (inherits(x, "Assay") && !inherits(x, "Assay5")) {
      sd = slot(x, "scale.data")
      if (!is.matrix(sd)) {
        message(sprintf("Fixing invalid scale.data in assay '%s' (was %s)", a, class(sd)[1]))
        slot(x, "scale.data", check = FALSE) = matrix(numeric(0), 0, 0)
        obj[[a]] = x
      }
    }
  }
  obj
}

write_h5ad_direct = function(obj, file, assay = DefaultAssay(obj),
                             layer = "data", reductions = Reductions(obj)) {
  ad = reticulate::import("anndata", convert = FALSE)

  X = GetAssayData(obj, assay = assay, layer = layer)
  if (!inherits(X, "dgCMatrix")) X = as(as.matrix(X), "dgCMatrix")

  meta = obj@meta.data
  meta[] = lapply(meta, function(col) if (is.factor(col)) as.character(col) else col)

  adata = ad$AnnData(
    X   = reticulate::r_to_py(Matrix::t(X)),         # cells x genes
    obs = reticulate::r_to_py(meta),
    var = reticulate::r_to_py(data.frame(row.names = rownames(X)))
  )

  for (r in reductions) {
    emb = Embeddings(obj, r)[colnames(X), , drop = FALSE]
    adata$obsm$`__setitem__`(paste0("X_", r), reticulate::r_to_py(unname(emb)))
  }

  adata$write_h5ad(file)
  message("Wrote ", file, " (", ncol(X), " cells x ", nrow(X), " genes; obsm: ",
          paste0("X_", reductions, collapse = ", "), ")")
  invisible(file)
}
# scVI latent (04.18, trained without a species key) in the same shape as the UINMF H
# list, so integrate_fastmnn / _cca / _harmony / _rpca can run on it (04.19-04.24)
read_scvi_latent = function(out_path) {
  Z = as.matrix(read.csv(file.path(out_path, 'scvi_latent.csv'), row.names = 1, check.names = FALSE))
  list(hsapiens  = Z[startsWith(rownames(Z), 'hsapiens_'), , drop = FALSE],
       mmusculus = Z[startsWith(rownames(Z), 'mmusculus_'), , drop = FALSE])
}

# the integrate_* helpers name their H-based reductions after UINMF: rename one
rename_reduction = function(data, from, to) {
  emb = Embeddings(data, from)
  key = paste0(gsub('[^A-Za-z]', '', to), '_')
  colnames(emb) = paste0(key, seq_len(ncol(emb)))
  data[[to]] = CreateDimReducObject(embeddings = emb, key = key, assay = DefaultAssay(data))
  data[[from]] = NULL
  data
}
