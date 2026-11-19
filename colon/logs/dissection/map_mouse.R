# Section 5b of ../EBCC_PIPELINE.md, step 2: bin the mouse section (~250 x 250 squares) and score each bin
# for target (colon wall) EBCC labels and label-independent off-target marker sets. Saves the per-cell
# table and the bin table for choosing the cut; no pipeline outputs are changed.
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
suppressMessages({ source('../library.R') })
mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
rownames(mouse) = mouse[['Xenium']]@meta.data$gene_name
X = LayerData(mouse, assay = 'Xenium', layer = 'data')
coords = GetTissueCoordinates(mouse)
cells = data.frame(cell = coords$cell, x = coords$x, y = coords$y)
cells$label = as.character(mouse$ebcc_label[match(cells$cell, colnames(mouse))])
offtarget = list(
  squamous  = c('Krt13', 'Krt16', 'Ivl', 'Lor', 'Lce1m', 'Sbsn', 'Cnfn', 'Krt79', 'Dsc3', 'Calm4', 'Asprv1', 'Crct1', 'Col17a1', 'Trim29', 'Serpinb5'),
  skeletal  = c('Actn3', 'Myf6', 'Myoz1', 'Pygm', 'Pgam2', 'Trdn', 'Tnnc1', 'Myoz2', 'Mylk3'),
  adipose   = c('Car3', 'Aqp7', 'Retn', 'Pck1', 'Cd36'),
  cartilage_bone = c('Col2a1', 'Acan', 'Comp', 'Ucma', 'Runx2'),
  urothelium = c('Upk1a', 'Upk1b', 'Upk3a'),
  colon_epi = c('Cdx1', 'Gpa33', 'Guca2b', 'Cdhr5', 'Car2'))   # on-target check, not off-target
stopifnot(all(unlist(offtarget) %in% rownames(X)))
for (k in names(offtarget)) cells[[k]] = colMeans(as.matrix(X[offtarget[[k]], cells$cell, drop = FALSE]))
B = 250
cells$bx = floor(cells$x / B) * B; cells$by = floor(cells$y / B) * B
saveRDS(list(cells = cells, offtarget = offtarget, B = B), 'logs/dissection/mouse_cells.rds')
message('saved logs/dissection/mouse_cells.rds: ', nrow(cells), ' cells')
