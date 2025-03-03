#!/usr/bin/env Rscript
# Michael VanInsberghe 2024-07-19
# output figures for WEE1 analysis

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(data.table)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(ggplot2)
library(ggpointdensity)
library(scales)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(forcats)
library(RColorBrewer)
library(viridis)
library(uwot)
library(ggrepel)
library(dendsort)
library(DESeq2)
library(fgsea)
library(msigdbr)

theme_set(theme_cowplot())

reduceCountMatrix <- function(x,y){
  allgenes <- unique(c(rownames(x), rownames(y)))
  notInX <- allgenes[!(allgenes %in% rownames(x))]
  notInY <- allgenes[!(allgenes %in% rownames(y))]

  if(length(notInX > 0)){
    x <- rbind(x,
               Matrix(0, 
                      nrow = length(notInX),
                      ncol = ncol(x),
                      dimnames = list(notInX, colnames(x)),
                      sparse = TRUE))
  }
  if(length(notInY > 0)){
    y <- rbind(y,
               Matrix(0, 
                      nrow = length(notInY),
                      ncol = ncol(y),
                      dimnames = list(notInY, colnames(y)),
                      sparse = TRUE))
  }

  mm <- merge.Matrix(x, y, 
                     by.x = rownames(x),
                     by.y = rownames(y),
                     all.x = TRUE,
                     all.y = TRUE,
                     out.class = "CsparseMatrix")
  colnames(mm) <- c(colnames(x), colnames(y))
  return(mm)
}

readMergeCounts <- function(samples){
  require(Matrix.utils)
  counts <- mapply(readCounts,
                   folder = samples$counts,
                   sname = samples$names)
  allcounts <- Reduce(reduceCountMatrix, counts)
  #rownames(allcounts) <- gsub("_", "-", rownames(allcounts))

  return(allcounts)
}

readCounts <- function(folder, sname = ""){
  require(Matrix)
  fn_mtx <- list.files(folder, pattern = "*matrix.mtx*")
  if(length(fn_mtx) > 1){
    stop(paste0("Too many count matrix *matrix.mtx* files found in folder: ", folder))
  }
  if(length(fn_mtx) < 1){
    stop(paste0("Cannot find count matrix *matrix.mtx* in folder: ", folder))
  }

  fn_features <- list.files(folder, pattern = "*features.tsv*")
  if(length(fn_features) > 1){
    stop(paste0("Too many feature files *features.tsv* found in folder: ", folder))
  }
  if(length(fn_features) < 1){
    stop(paste0("Cannot find feature tsv *features.tsv* in folder: ", folder))
  }

  fn_barcodes <- list.files(folder, pattern = "*barcodes.tsv*")
  if(length(fn_barcodes) > 1){
    stop(paste0("Too many barcode files *barcodes.tsv* found in folder: ", folder))
  }
  if(length(fn_barcodes) < 1){
    stop(paste0("Cannot find barcode tsv *barcodes.tsv* in folder: ", folder))
  }


  counts <- readMM(file.path(folder,fn_mtx))
  rownames(counts) <- read.table(file.path(folder, fn_features), 
                                 header = FALSE,
                                 #sep = "\t",
                                 stringsAsFactors = FALSE)$V1
  if(nchar(sname) > 0){
    colnames(counts) <- paste(sname, 
                              read.table(file.path(folder, fn_barcodes),
                                         header = FALSE,
                                         #sep = "\t",
                                         stringsAsFactors = FALSE)$V1,
                              sep = "_")
  }else{
    colnames(counts) <- read.table(file.path(folder, fn_barcodes),
                                   header = FALSE,
                                   #sep = "\t",
                                   stringsAsFactors = FALSE)$V1
  }

  return(counts)
}

wiggleRegionDT <- function(reads, site = 'cut5', reference = 'cds_start', mindist = -40, maxdist = 60, by_var = 'CB'){
  group_var <- c('csd', by_var)
  if(length(by_var) > 0){
    cj <- reads[, do.call(CJ, c(c(list(csd = full_seq(mindist:maxdist,1)), .SD), unique = TRUE)), .SDcols = by_var]
  }else{
    cj <- CJ(csd=full_seq(mindist:maxdist,1), unique = TRUE)
  }

  region <- reads[, csd := .SD[[site]] - .SD[[reference]]
                  ][csd >= mindist & csd <= maxdist,
                  ][, .(n = .N), by = group_var
                  ][cj, on = group_var
                  ][, n := lapply(.SD, nafill, fill = 0), .SDcols = 'n']
  return(region)
}

scaleDFDT <- function(reads, site = 'cut5', utr5 = 100, cds = 200, utr3 = 100, by_var = 'CB'){
  group_var <- c('pos', by_var)
  if(length(by_var) > 0){
    cj <- reads[, do.call(CJ, c(c(list(pos = full_seq(c(1,utr5+cds+utr3),1)), .SD), unique = TRUE)), .SDcols = by_var]
  }else{
    cj <- CJ(pos = full_seq(c(1,utr5+cds+utr3),1), unique = TRUE)
  }

  scaled <- reads[, pos := fcase( .SD[[site]] < l_utr5 , round( (.SD[[site]]/l_utr5)*(utr5-1) )+1,
                                 ((.SD[[site]] >= l_utr5) & (.SD[[site]] < (l_utr5+l_cds) )) , round( ((.SD[[site]]-l_utr5)/l_cds)*(cds-1) ) + 1 + utr5,
                                 (.SD[[site]] >= (l_utr5 + l_cds)) , round( ((.SD[[site]] - l_utr5 - l_cds)/l_utr3)*(utr3-1) ) + utr5 + cds + 1 )
                 ][, .(n = .N), by = group_var
                 ][cj, on = group_var
                 ][, n := lapply(.SD, nafill, fill = 0), .SDcols = 'n']
  return(scaled)
}






figurePrefix <- 'wee1_figures/wee1_figures'



data_dir <- 'analysis/WEE1/data'







annot <- fread(file.path(data_dir, 'hsa_annotations.csv.gz'))
canonical_pc <- annot[
                      ][set == 'canonical' & transcript_type == 'protein_coding',
                      ][transcript_id != 'ENST00000437139.7',
                      ][!grepl("^PCDH", gene_name),
                      ]

rpfMeta <- readRDS(file.path(data_dir, 'rpf', 'compiled', 'meta.rds'))
rpfCounts <- readRDS(file.path(data_dir, 'rpf', 'compiled', 'counts.rds'))

vasaMeta <- readRDS(file.path(data_dir, 'vasa', 'compiled', 'meta.rds'))
vasaCounts <- readRDS(file.path(data_dir, 'vasa', 'compiled', 'counts.rds'))

allMeta <- setDT(bind_rows(rpfMeta %>%
                     #select(intersect(colnames(rpfMeta), colnames(vasaMeta))) %>%
                     mutate(type = "RPF"),
                     vasaMeta %>%
                       #select(intersect(colnames(rpfMeta), colnames(vasaMeta))) %>%
                       mutate(type = "VASA") ))
#rownames(allMeta) <- allMeta$CB
allMeta$treatment <- allMeta$sort_population

allCounts <- Reduce(reduceCountMatrix, list(rpfCounts, vasaCounts))
allCounts <- allCounts[, allMeta$CB]
allCounts <- allCounts[apply(allCounts>0, 1, sum)>3,]

ENSG_to_name <- data.frame(count_names = rownames(allCounts), stringsAsFactors = FALSE) %>%
  inner_join(annot, by = c("count_names" = "gene_id")) %>% 
  mutate(gene_out = paste(count_names, gene_name, sep = "-")) %>%
  select(count_names, gene_out, gene_name) %>%
  distinct() %>%
  column_to_rownames(var = "count_names")
ENSG_to_name$gene_id = rownames(ENSG_to_name)

rownames(allCounts) <- ENSG_to_name[rownames(allCounts), "gene_out"]
allCounts <- allCounts[!grepl("-PCDH", rownames(allCounts)),]

## output meta and counts for GEO
#outMeta <- allMeta
#outMeta[, sort_population := NULL]
#fwrite(outMeta, 'GEO/prep/metadata.csv.gz')

#outCounts <- as.data.table(allCounts, keep.rownames = 'gene')
#fwrite(outCounts, 'GEO/prep/counts.csv.gz')

longCounts <- as.data.frame(t(allCounts))
longCounts$CB <- rownames(longCounts)
longCounts <- longCounts %>%
  pivot_longer(cols = starts_with("ENSG"), names_to = "gene_id", values_to = "counts")

setDT(longCounts)


# tRNA reads to get size factor
soloSamples <- data.frame(names <- c("WEE11"),
                          counts = c(file.path(data_dir, "rpf", "solo", "RPFv4B-WEE11_Solo.out", "Gene", "raw")))

soloCounts <- readMergeCounts(soloSamples)[, colnames(rpfCounts)]

sumFactors <- bind_rows(data.frame(CB = colnames(soloCounts),
                                   SF = scran::calculateSumFactors(soloCounts[grepl('^tRNACluster', rownames(soloCounts)),], positive = FALSE)),
                        data.frame(CB = pull(allMeta %>% filter(type == 'VASA') %>% select(CB)),
                                   SF = scran::calculateSumFactors(allCounts[, pull(allMeta %>% filter(type == 'VASA') %>% select(CB))], positive = FALSE) ))
rownames(sumFactors) <- sumFactors$CB

setDT(sumFactors)


## QC
# correlation
sumFactors <- as.data.frame(sumFactors)
rownames(sumFactors) <- sumFactors$CB

allMeta <- as.data.frame(allMeta)
rownames(allMeta) <- allMeta$CB

normCounts <- sweep(allCounts, 2, sumFactors[colnames(allCounts), 'SF'], FUN = "/")

countsCor <- cor(as.matrix(normCounts), method = 'spearman')


cell_annotation <- HeatmapAnnotation(type       = as.character(allMeta[rownames(countsCor), 'type']),
                                     treatment  = as.character(allMeta[rownames(countsCor), 'treatment']),
                                     replicate  = as.character(allMeta[rownames(countsCor), 'replicate']),
                                     col = list(type       = setNames(c('#e41a1c', '#377eb8'), c('RPF', 'VASA')),
                                                treatment  = setNames(brewer.pal(length(unique(allMeta$treatment)), 'Set2')[1:length(unique(allMeta$treatment))], unique(allMeta$treatment)),
                                                replicate  = setNames(brewer.pal(length(unique(allMeta$replicate)), 'Pastel1')[1:length(unique(allMeta$replicate))], unique(allMeta$replicate))
                                                ),
                                     which = 'column')

cell_annotation_r <- HeatmapAnnotation(type       = as.character(allMeta[rownames(countsCor), 'type']),
                                     treatment  = as.character(allMeta[rownames(countsCor), 'treatment']),
                                     replicate  = as.character(allMeta[rownames(countsCor), 'replicate']),
                                     col = list(type       = setNames(c('#e41a1c', '#377eb8'), c('RPF', 'VASA')),
                                                treatment  = setNames(brewer.pal(length(unique(allMeta$treatment)), 'Set2')[1:length(unique(allMeta$treatment))], unique(allMeta$treatment)),
                                                replicate  = setNames(brewer.pal(length(unique(allMeta$replicate)), 'Pastel1')[1:length(unique(allMeta$replicate))], unique(allMeta$replicate))
                                                ),
                                       which = 'row',
                                       show_legend = FALSE)

row_dend <- dendsort(hclust(dist(as.matrix((countsCor)), method = 'euclidean'), method = 'ward.D2'), type = 'average', isReverse = FALSE)

pdf(paste0(figurePrefix, '_all_spearman.pdf'), height = 9, width = 10)
dc <- Heatmap(as.matrix(countsCor),
        top_annotation = cell_annotation,
        left_annotation = cell_annotation_r,
        col = colorRamp2(seq(0.8,1, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        na_col = "grey0",
        name = "spearman",
        cluster_rows = row_dend,
        show_row_dend = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        cluster_columns = row_dend,
        clustering_distance_columns = "pearson",
        clustering_method_columns = "ward.D2",
        show_column_dend = TRUE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_labels = gsub("_", "\\-", allMeta[rownames(countsCor), 'sort_population']),
        column_labels = gsub("_", "\\-", allMeta[rownames(countsCor), 'sort_population']),
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 3),
        column_title_rot = 90,
        border = FALSE,
        use_raster = FALSE,
        raster_resize_mat= TRUE,
        raster_device = "CairoPNG",
        raster_by_magick = FALSE,
        raster_quality = 15)
dc
dev.off()

# find variable genes for PCA

cv <- bind_rows(data.frame(mean = apply(normCounts[, pull(allMeta %>% filter(type == 'VASA') %>% select(CB))], 1, mean),
                           sd = apply(normCounts[, pull(allMeta %>% filter(type == 'VASA') %>% select(CB))],1,sd),
                           type = 'VASA',
                           gene_out = rownames(normCounts)),
                data.frame(mean = apply(normCounts[, pull(allMeta %>% filter(type == 'RPF') %>% select(CB))], 1, mean),
                           sd = apply(normCounts[, pull(allMeta %>% filter(type == 'RPF') %>% select(CB))],1,sd),
                           type = 'RPF',
                           gene_out = rownames(normCounts)) ) %>%
      mutate(logcv = log10(sd/mean),
             logme = log10(mean))

cv_fits <- lapply(split(cv, cv$type), function(x) {
                    robustbase::ltsReg(logcv ~ logme, data = x, nsamp = 1000, seed = set.seed(516), adjust = TRUE, alpha = 0.99)
       } )

fit_coefs <- rbindlist(lapply(cv_fits, function(x) {
                              data.frame(intercept = x$coefficients[1],
                                         slope     = x$coefficients[2])
                           } ), 
                       idcol = 'type')

plt_cv <- ggplot(cv, aes(x=log10(mean), y=log10(sd/mean), colour = type))+
  geom_point(size=0.5, alpha = 0.05)+
  geom_abline(data = fit_coefs, aes(slope = slope, intercept = intercept, colour = type), linetype = 'dashed')+
  facet_wrap(~type)
#dontsave save_plot(paste0(figurePrefix, '_plt_cv_mean.pdf'), plt_cv)

res <- rbindlist(mapply(function(dat,fit) {
                dat$res <- dat$logcv - (fit$coefficients[2]*dat$logme + fit$coefficients[1])
                return(dat)
                       },
              split(cv,cv$type), cv_fits, SIMPLIFY = FALSE) )


plt_residuals <- ggplot(res, aes(x=res))+
  geom_histogram(bins = 200)+
  facet_wrap(~type)+
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'black')+
  geom_vline(xintercept = 0.20, linetype = 'dashed', colour = 'grey70')
#dontsave save_plot(paste0(figurePrefix, '_cvfit_residual_hist.pdf'), plt_residuals, base_width = 8, base_height = 5)


res <- res %>%
  group_by(type) %>%
  mutate(res_rank = dense_rank(-res)) %>% 
  ungroup() %>%
  mutate(variable = res > .20 & logme > 1.5 )

plt_residuals <- ggplot(res, aes(x=logme, y=res, colour = variable))+
  geom_point(size=0.5, alpha = 0.1)+
  facet_wrap(~type)+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_hline(yintercept = 0.20, linetype = 'dashed', colour = 'grey70')+
  geom_vline(xintercept = 1.5, linetype = 'dashed', colour = 'grey70')
#dontsave save_plot(paste0(figurePrefix, '_cvfit_residual.pdf'), plt_residuals, base_width = 8, base_height = 5)


## all PCA
lnCounts <- log10(normCounts+1)
all_pca <- prcomp(t(lnCounts[unique(c(pull(res %>% filter(type == 'RPF',  res_rank <= 2000) %>% select(gene_out)),
                                       pull(res %>% filter(type == 'VASA', res_rank <= 2000) %>% select(gene_out)))), 
                             pull(allMeta %>% select(CB)) ]))
all_perVar <- data.frame(pct_var = all_pca$sdev^2 / sum(all_pca$sdev^2),
                         pcn     = colnames(all_pca$x),
                         PC      = 1:ncol(all_pca$x))
all_pcad <- as.data.frame(all_pca$x) %>%
  rownames_to_column('CB') %>%
  inner_join(allMeta, by = 'CB') %>%
  mutate(type_replicate = paste(treatment, replicate, sep = '_')) 

plt_all_pca_perVar <- ggplot(all_perVar %>% filter(PC < 20), aes(x=PC,y=100*pct_var))+
  geom_point()+
  ylab('% Variance')
#dontsave save_plot(paste0(figurePrefix, '_all_perVar.pdf'), plt_all_pca_perVar)

plt_all_pca_type <- ggplot(all_pcad, aes(x=PC1, y=PC2, colour = type, shape = treatment))+
  geom_point()+
  xlab(sprintf("PC1: %.1f %%", 100*pull(all_perVar %>% filter(pcn == 'PC1') %>% select(pct_var))))+
  ylab(sprintf("PC2: %.1f %%", 100*pull(all_perVar %>% filter(pcn == 'PC2') %>% select(pct_var))))+
  coord_equal()
plt_all_pca_treatment <- ggplot(all_pcad, aes(x=PC1, y=PC2, colour = treatment))+
  geom_point()+
  xlab(sprintf("PC1: %.1f %%", 100*pull(all_perVar %>% filter(pcn == 'PC1') %>% select(pct_var))))+
  ylab(sprintf("PC2: %.1f %%", 100*pull(all_perVar %>% filter(pcn == 'PC2') %>% select(pct_var))))+
  coord_equal()
plt_all_pca_replicate <- ggplot(all_pcad, aes(x=PC1, y=PC2, colour = as.factor(replicate)))+
  geom_point()+
  xlab(sprintf("PC1: %.1f %%", 100*pull(all_perVar %>% filter(pcn == 'PC1') %>% select(pct_var))))+
  ylab(sprintf("PC2: %.1f %%", 100*pull(all_perVar %>% filter(pcn == 'PC2') %>% select(pct_var))))+
  coord_equal()

#dontsave save_plot(paste0(figurePrefix, 'all_pca_12.pdf'), plot_grid(plt_all_pca_type, plt_all_pca_treatment, plt_all_pca_replicate), base_width = 10, base_height = 7)
save_plot(paste0(figurePrefix, 'all_pca_PC12_var.pdf'), plot_grid(plt_all_pca_type, plt_all_pca_perVar, rel_widths = c(1,.65) ), base_width = 8, base_height = 4)



## ribo QC

rpfReads <- fread(file.path(data_dir, 'rpf', 'compiled', 'reads.csv.gz'))#[length >= 30,]

RPFwiggle <- rbindlist(list(start = wiggleRegionDT(rpfReads,
                                                   mindist = -40, maxdist = 60,
                                                   reference = 'cds_start',
                                                   by_var = c('CB'))[, display_rank := dense_rank(csd)],
                            cds = wiggleRegionDT(rpfReads,
                                                 mindist = 100, maxdist = 200,
                                                 reference = 'cds_start',
                                                 by_var = c('CB'))[, display_rank := dense_rank(csd)+1E6],
                            stop = wiggleRegionDT(rpfReads,
                                                  mindist = -60, maxdist = 40,
                                                  reference = 'cds_end',
                                                  by_var = c('CB'))[, display_rank := dense_rank(csd)+2E6]
                            ), 
                       idcol = 'region')[rpfMeta[, c('CB', 'cell_total')], on = 'CB', nomatch = NULL
                                         ][sumFactors, on = 'CB', nomatch = NULL
                                         ][, pct := n/SF
                                         #][, pct := n/cell_total
                                         ][, display_position := dense_rank(display_rank)
                                         ]


wiggleHM <- as.data.frame(dcast(RPFwiggle, CB ~ display_position, value.var = "pct"))
rownames(wiggleHM) <- wiggleHM[,'CB']
wiggleHM <- wiggleHM[,!(colnames(wiggleHM) == 'CB')]

cell_annotation <- HeatmapAnnotation(
                                     treatment  = as.character(allMeta[rownames(wiggleHM), 'treatment']),
                                     replicate  = as.character(allMeta[rownames(wiggleHM), 'replicate']),
                                     col = list(
                                                treatment  = setNames(brewer.pal(length(unique(allMeta$treatment)), 'Set2')[1:length(unique(allMeta$treatment))], unique(allMeta$treatment)),
                                                replicate  = setNames(brewer.pal(length(unique(allMeta$replicate)), 'Pastel1')[1:length(unique(allMeta$replicate))], unique(allMeta$replicate))
                                                ),
                                     which = 'row')


row_dend <- dendsort(hclust(dist(wiggleHM, method = 'euclidean'), method = 'ward.D2'), type = 'min', isReverse = FALSE)


pdf(paste0(figurePrefix, '_CB_wigglehm_cut5.pdf'), height = 25, width = 40)
labelby = 20
length_annotation <- HeatmapAnnotation(pos = anno_mark(at = c(seq(1,101,by=labelby),
                                                              seq(102,202, by=labelby),
                                                              seq(203,303, by=labelby)),
                                                       labels = c(as.character(seq(-40,60,labelby)),
                                                                  as.character(seq(100,200,labelby)),
                                                                  as.character(seq(-60,40,labelby))),
                                                       which = "column", side = "bottom", labels_rot = 0, link_height = unit(3, "mm")))
Heatmap(as.matrix(wiggleHM),
        #col = colorRamp2(seq(0, 10, length = 11), rev(brewer.pal(11, "RdYlBu"))), # log fold change
        bottom_annotation = length_annotation,
        left_annotation = cell_annotation,
        na_col = "grey0",
        name = "psites",
        cluster_rows = row_dend,
        clustering_distance_rows = "spearman",
        clustering_method_rows = "ward.D2",
        show_row_dend = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = factor(c(rep("start_codon", length(-40:60)), rep("cds", length(100:200)), rep("stop_codon", length(-40:60))), levels = c("start_codon", "cds", "stop_codon"), ordered = TRUE),
        column_title = NULL,
        border = FALSE,
        use_raster = TRUE,
        raster_resize_mat= TRUE,
        raster_quality = 2)
dev.off()



expanded_canonical <- canonical_pc[
                                   ][, `:=` (lr_utr5 = l_utr5 - 25, lr_cds = l_cds + 25 + 25, lr_utr3 = as.numeric(l_utr3-25))
                                   ][lr_utr5 < 0, lr_utr5 := 0
                                   ][lr_cds < 0, lr_cds := 0
                                   ][lr_utr3 < 0, lr_utr3 := 0
                                   ]

regionBases <- melt(rpfReads[
                    #][length >= 27,
                        ][expanded_canonical[, c('transcript_id', 'lr_utr5', 'lr_cds', 'lr_utr3')], on = 'transcript_id', nomatch = NULL
                        ][, .(utr5 = sum(lr_utr5), cds = sum(lr_cds), utr3 = sum(lr_utr3)), by = c('CB')],
                    id.vars = c('CB'),
                    measure.vars = c('utr5', 'cds', 'utr3'),
                    variable.name = 'region', value.name = 'bp_sep')[
                                                                     #][melt(unique(rpfReads[length >= 27,][, c('CB', 'transcript_id')])[
                                                                     ][melt(unique(rpfReads[, c('CB', 'transcript_id')])[
                                                                            ][expanded_canonical[, c('transcript_id', 'lr_utr5', 'lr_cds', 'lr_utr3')], on = 'transcript_id', nomatch = NULL
                                                                            ][, .(utr5 = sum(lr_utr5), cds = sum(lr_cds), utr3 = sum(lr_utr3)), by = c('CB')],
                                                                     id.vars = c('CB'),
                                                                     measure.vars = c('utr5', 'cds', 'utr3'),
                                                                     variable.name = 'region', value.name = 'bp_one'), on = c('CB', 'region'), nomatch = NULL
                                                                     ]


regionTPM <- rpfReads[
                      #][length >= 27,
                      ][expanded_canonical[, c('transcript_id', 'lr_utr5', 'lr_cds', 'lr_utr3')], on = 'transcript_id', nomatch = NULL
                      ][, read_class := fcase(cut5 < lr_utr5, 'utr5',
                                              (cut5 >= lr_utr5) & (cut5 < (lr_utr5 + lr_cds)), 'cds',
                                              cut5 >= (lr_utr5 + lr_cds), 'utr3')
                      ][, .(reads_per_region = .N), by = c('CB', 'read_class')
                      ][regionBases, on = c('CB' = 'CB', 'read_class' = 'region')
                      ][, rph_one := reads_per_region/bp_one
                      ][, rph_sep := reads_per_region/bp_sep
                      ][, `:=` (norm_one = sum(rph_one)/100, norm_sep = sum(rph_sep)/100), by = 'CB'
                      ][, tph_one := rph_one/norm_one
                      ][, tph_sep := rph_sep/norm_sep
                      ]

plt_region <- melt(regionTPM[
                    ][allMeta[, c('CB', 'treatment', 'plate', 'replicate')], on = 'CB', nomatch = NULL
                    ][, c('CB', 'tph_one', 'tph_sep', 'read_class', 'treatment', 'replicate', 'plate')],
                   id.vars = c('CB', 'read_class', 'treatment', 'replicate', 'plate'),
                   measure.vars = c('tph_one', 'tph_sep'),
                   variable.name = 'tph_type', value.name = 'tph') %>%
  mutate(regionorder = case_when(read_class == "utr5" ~ 0,
                                 read_class == "cds"  ~ 1,
                                 read_class == "utr3" ~ 2)) %>%
  mutate(read_class = fct_reorder(read_class, regionorder)) %>%
  ggplot(aes(x=read_class, y=tph, colour = treatment))+
  geom_point(position = position_jitterdodge(jitter.width = 0.3), size = 1) +
  facet_wrap(~tph_type, nrow = 1) +
  ylab('Reads per region length per hundred')
save_plot(paste0(figurePrefix, "_plt_region_tpm.pdf"), plt_region, base_width = 8, base_height = 5)


frameSummary <- rbindlist(list(rpfReads[
                                        ][(cut5 > l_utr5 + 15) &
                                          (cut3 < (l_utr5 + l_cds - 15)),
                                        ][, frame := (psite - cds_start) %% 3
                                        ][, .(n = .N), by = c('CB', 'frame')
                                        ][, type := 'frameP'
                                        ][, cell_total := sum(n), by= 'CB'
                                        ],
                               rpfReads[
                                        ][(cut5 > l_utr5 + 15) &
                                          (cut3 < (l_utr5 + l_cds - 15)),
                                        ][, frame := (cut5 - cds_start) %% 3
                                        ][, .(n = .N), by = c('CB', 'frame')
                                        ][, type := 'frame5'
                                        ][, cell_total := sum(n), by= 'CB'
                                        ]))

plt_frame <- frameSummary[
                      ][allMeta[, c('CB', 'treatment', 'plate', 'replicate')], on = 'CB', nomatch = NULL
                      ] %>% 
                        ggplot(aes(x=frame, y=n/cell_total, colour = treatment))+
                        geom_point(position = position_jitterdodge(jitter.width = .3, seed = 516), size = 1)+
                        facet_wrap(~type, nrow = 1)+
                        ylab('Fraction of reads')
save_plot(paste0(figurePrefix, '_frame.pdf'), plt_frame, base_width = 6, base_height = 3)

save_plot(paste0(figurePrefix, '_region_frame.pdf'), plot_grid(plt_region, plt_frame, nrow = 1), base_width = 12, base_height = 5)



## translation
# global decrease
setDT(sumFactors)
setDT(allMeta)
total_comparison <- sumFactors[
                   ][allMeta[type=='RPF',][, c('CB', 'sort_population', 'cds', 'utr5', 'utr3', 'cell_total')], on = 'CB', nomatch = NULL
                   ][, norm_factor := pull(sumFactors[
                                             ][allMeta[type=='RPF',][, c('CB', 'sort_population', 'cds')], on = 'CB', nomatch = NULL
                                             ][sort_population == 'DMSO',
                                             ][, .( norm_factor = mean(cds/SF))
                                             ])
                   ]
plt_total <- ggplot(total_comparison, aes(x=sort_population, y=(cds/SF)/norm_factor, colour = sort_population))+
  geom_point(position = position_jitterdodge(jitter.width = 0.3, seed=516), size = 1)+
  xlab('')+
  ylab('relative proportion of CDS reads')
save_plot(paste0(figurePrefix, '_total_rpf.pdf'), plt_total, base_width = 4)

total_comparison[
                 ][, proportion := (cds/SF)/norm_factor
                 ][, .(mean_prop = mean(proportion), sem_prop = sd(proportion)/sqrt(.N)), by = 'sort_population'
                 ]



# metagene coverage 
lr_means <- regionBases[
                        ][, .(mean_bp_sep = median(bp_sep), 
                              mean_bp_one = median(bp_one)), by = region
                        ]

lr_cds <- 300
lr_utr3 <- ceiling(lr_cds*(pull(lr_means[region == 'utr3', 'mean_bp_one'])/pull(lr_means[region == 'cds', 'mean_bp_one'])))
lr_utr5 <- ceiling(lr_cds*(pull(lr_means[region == 'utr5', 'mean_bp_one'])/pull(lr_means[region == 'cds', 'mean_bp_one'])))

# ENSG00000062716.13 has a weird spike in the 3' UTR, remove it
covReads <- scaleDFDT(rpfReads[gene_id != 'ENSG00000062716.13',], site='cut5', utr5 = lr_utr5, cds = lr_cds, utr3 = lr_utr3, by_var = 'CB')[
                                                                                                                                            ][sumFactors, on = 'CB', nomatch = NULL
                                                                                                                                            ]

covReadsAv <- covReads[
                       ][rpfMeta[, c('CB', 'sort_population')], on = 'CB'
                       ][, .(mean_n = mean(n/SF), sem_n = sd(n/SF)/sqrt(.N)), by = c('pos', 'sort_population')
                       ]

plt_coverage_average <- ggplot(covReadsAv, aes(x=pos, y=mean_n, colour = sort_population))+
  geom_ribbon(aes(ymin=mean_n - sem_n, ymax = mean_n + sem_n, fill = sort_population), alpha = 0.5, colour = NA)+
  geom_line(linewidth = 0.5)+
  #geom_ribbon(aes(ymin=mean_norm_counts - sem_norm_counts, ymax=mean_norm_counts + sem_norm_counts, fill = ucond), alpha = 0.5)+
  scale_x_continuous(breaks=c(0,(lr_utr5/2),lr_utr5,(lr_utr5+lr_cds/2),(lr_utr5+lr_cds),(lr_utr5+lr_cds+(lr_utr3/2)),(lr_utr5+lr_cds+lr_utr3)),
                     labels=c("","5' UTR","","CDS","","3' UTR",""))+
  #facet_wrap(~sort_population)+
  xlab("position along scaled metagene")+
  ylab("normalized reads")
save_plot(paste0(figurePrefix, "_metagene_coverage_repaverage.pdf"), plt_coverage_average, base_width = 5.936, base_height = 3.71)


# dTE with tRNA SF
sumFactors <- as.data.frame(sumFactors)
rownames(sumFactors) <- sumFactors$CB

expressedGenes <- ( (rowSums(allCounts[, allMeta %>% filter(type == "VASA") %>% select(CB) %>% pull()] >= 5) >= 3) & (rowSums(allCounts[, allMeta %>% filter(type == "RPF") %>% select(CB) %>% pull()] >= 5) >= 3) )
expressedGenes <- names(which(expressedGenes))

testMeta <- as_tibble(allMeta) %>%
  select(CB, type, treatment) %>%
  mutate(typeorder = case_when(type == 'VASA' ~ 0,
                               type == 'RPF'  ~ 1)) %>%
  mutate(type = fct_reorder(type, typeorder)) %>%
  mutate(treatorder = case_when(treatment == 'DMSO' ~ 0,
                                treatment == 'WEE1i' ~ 1)) %>%
  mutate(treatment = fct_reorder(treatment, treatorder)) %>%
  mutate(group = paste(type, treatment, sep = ".")) %>%
  select(-ends_with('order')) %>%
  column_to_rownames('CB') %>%
  mutate_if(is.character, as.factor)

testCounts <- allCounts[expressedGenes, rownames(testMeta)]

rpfFullMat <- model.matrix( ~ type + treatment + type:treatment, data = testMeta)
rpfReduMat <- model.matrix( ~ type + treatment,                  data = testMeta)

tDds <- DESeqDataSetFromMatrix(countData = testCounts,
                               colData   = testMeta,
                               design    = rpfFullMat)
sizeFactors(tDds) <- sumFactors[colnames(tDds), 'SF']
tDds <- DESeq(tDds, test = 'LRT', full = rpfFullMat, reduced = rpfReduMat)

rpfRes <- results(tDds, independentFiltering = TRUE, format = 'DataFrame', alpha = 0.1) %>%
  data.frame() %>%
  rownames_to_column('gene_out')

# coefs <- coef(tDds, SE = FALSE) %>%
#   data.frame() %>%
#   rownames_to_column('gene_out')

plt_volcano_tSF <- rpfRes %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj)))+
  geom_pointdensity(size = 0.5)+
  scale_color_viridis()+
  geom_hline(yintercept = -log10(0.01), linetype = 'dashed', colour = 'grey70', alpha = .5)#+
save_plot(paste0(figurePrefix, '_tsf_dTE_volcano.pdf'), plt_volcano_tSF, base_height = 6/1.5, base_width = 9/1.5)



fgsea_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  mutate(gs_name = gsub('^REACTOME_', '', gs_name)) %>%
  split(x = .$gene_symbol, f = .$gs_name)

rpfResRank <- rpfRes %>% 
  mutate(stat = stat * sign(log2FoldChange)) %>%
  inner_join(ENSG_to_name, by = 'gene_out') %>%
  select(gene_name, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene_name) %>%
  summarize(stat = mean(stat)) %>%
  ungroup() %>%
  deframe()


rpfReact_fgsea_res <- fgsea(fgsea_reactome, stats = rpfResRank)

pdf(paste0(figurePrefix, '_tsf_dTE_gsea_table.pdf'), width = 15, height = 10)
plotGseaTable(pathways = fgsea_reactome[c(rpfReact_fgsea_res %>% filter(padj < 0.01) %>% filter(NES > 0) %>% slice_max(NES, n=20) %>% select(pathway) %>% pull(),
                                          rpfReact_fgsea_res %>% filter(padj < 0.01) %>% filter(NES < 0) %>% slice_min(NES, n=20) %>% select(pathway) %>% pull() %>% rev())],
              stats = rpfResRank,
              fgseaRes = rpfReact_fgsea_res,
              gseaParam = 1.0)
dev.off()


# dTE with RPF SF

rpfSumFactor <- bind_rows(data.frame(CB = pull(allMeta %>% filter(type == 'RPF') %>% select(CB)),
                                     SF = scran::calculateSumFactors(allCounts[, pull(allMeta %>% filter(type == 'RPF') %>% select(CB))], positive = FALSE)),
                          data.frame(CB = pull(allMeta %>% filter(type == 'VASA') %>% select(CB)),
                                     SF = scran::calculateSumFactors(allCounts[, pull(allMeta %>% filter(type == 'VASA') %>% select(CB))], positive = FALSE) ))
rownames(rpfSumFactor) <- rpfSumFactor$CB

rDds <- DESeqDataSetFromMatrix(countData = testCounts,
                               colData   = testMeta,
                               design    = rpfFullMat)
sizeFactors(rDds) <- rpfSumFactor[colnames(rDds), 'SF']
rDds <- DESeq(rDds, test = 'LRT', full = rpfFullMat, reduced = rpfReduMat)

rpfRFRes <- results(rDds, independentFiltering = TRUE, format = 'DataFrame', alpha = 0.1) %>%
  data.frame() %>%
  rownames_to_column('gene_out')


plt_volcano_rSF <- rpfRFRes %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj)))+
  geom_pointdensity(size = 0.5)+
  scale_color_viridis()+
  geom_hline(yintercept = -log10(0.01), linetype = 'dashed', colour = 'grey70')#+
save_plot(paste0(figurePrefix, '_rsf_dTE_volcano.pdf'), plt_volcano_rSF, base_height = 6/1.5, base_width = 9/1.5)


rpfRFResRank <- rpfRFRes %>% 
  mutate(stat = stat * sign(log2FoldChange)) %>%
  inner_join(ENSG_to_name, by = 'gene_out') %>%
  select(gene_name, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene_name) %>%
  summarize(stat = mean(stat)) %>%
  ungroup() %>%
  deframe()


rpfRFReact_fgsea_res <- fgsea(fgsea_reactome, stats = rpfRFResRank)

pdf(paste0(figurePrefix, '_rsf_dTE_gsea_table.pdf'), width = 15, height = 10)
plotGseaTable(pathways = fgsea_reactome[c(rpfRFReact_fgsea_res %>% filter(padj < 0.01) %>% filter(NES > 0) %>% slice_max(NES, n=20) %>% select(pathway) %>% pull(),
                                          rpfRFReact_fgsea_res %>% filter(padj < 0.01) %>% filter(NES < 0) %>% slice_min(NES, n=20) %>% select(pathway) %>% pull() %>% rev())],
              stats = rpfRFResRank,
              fgseaRes = rpfRFReact_fgsea_res,
              gseaParam = 1.0)
dev.off()


### Pausing

longNormCounts <- as.data.frame(t(normCounts))
longNormCounts$CB <- rownames(longNormCounts)
longNormCounts <- longNormCounts %>%
  pivot_longer(cols = starts_with('ENSG'), names_to = 'gene_out', values_to = 'counts')
setDT(longNormCounts)

avNC <- longNormCounts[
                       ][allMeta[, c('CB', 'treatment')], on = 'CB', nomatch = NULL
                       ][, .(mean_counts = mean(counts)), by = c('treatment', 'gene_out')
                       ][, lq := quantile(mean_counts, .10), by = 'treatment'
                       ][, uq := quantile(mean_counts, .90), by = 'treatment'
                       ][, inrange := (mean_counts >= lq) & (mean_counts <= uq)
                       ][, .(n_inrange = .N), by = 'gene_out'
                       ][n_inrange == 2
                       ][gene_out %in% expressedGenes
                       ][ENSG_to_name[, c('gene_out', 'gene_id')], on = 'gene_out', nomatch = NULL
                       ][canonical_pc[, c('gene_id', 'transcript_id')], on = 'gene_id', nomatch = NULL
                       ]


codons <- fread(file.path(data_dir, 'hsa_codons.csv.gz'))[
                                                          ][transcript_id %in% canonical_pc$transcript_id,
                                                          #][transcript_id %in% avNC$transcript_id
                                                          ][nchar(codon) == 3,
                                                          ][, start_codon_rank := dense_rank(codon_position), by = 'transcript_id'
                                                          ][, stop_codon_rank  := dense_rank(-codon_position), by = 'transcript_id'
                                                          ][start_codon_rank > 15 & stop_codon_rank > 15,
                                                          ]

codon_reads <- rpfReads[
                        ][id %in% pull(rpfReads[, .(n = .N), by = 'id'][n==1, 'id'])
                        ][, c('CB', 'transcript_id', 'psite'), with = FALSE
                        ][, readpos := .I
                        ][data.table(readpos = rep(1:.N, 7),
                                     site = rep(c('em2', 'em1', 'e', 'p', 'a', 'ap1', 'ap2'), .N),
                                     psite_diff =  rep(c(-9, -6, -3, 0, 3, 6, 9), .N) ), on = 'readpos'
                        ][, site_position := psite + psite_diff
                        ][, psite := NULL
                        ][codons, on = c('transcript_id' = 'transcript_id', 'site_position' = 'codon_position'), nomatch = NULL
                        ][three != 'Ter',
                        #][, sites_per_read := .N, by = 'readpos'
                        #][sites_per_read == 7
                        #][, sites_per_read := NULL
                        ][, readpos := NULL
                        ][, psite_diff := NULL
                        ][, start_codon_rank := NULL
                        ][, stop_codon_rank := NULL
                        ]

site_order <- data.table(site = c("em2", "em1", "e", "p", "a", "ap1", "ap2"),
                         site_order = c(1, 2, 3, 4, 5, 6, 7))




baseline_counts <- codon_reads[
                               ][!three %in% c('Ter'),
                               ][, codon_display := paste(gsub('T', 'U', codon), three, sep = '_')
                               ][, .(n = .N), by = c('CB', 'codon_display', 'site')
                               ]

total_baseline <- baseline_counts[
                                  #][CB %in% pull(allMeta %>% filter(type == 'RPF' & treatment == 'DMSO') %>% select(CB)),
                                  ][, .(codon_counts = sum(n)), by = c('CB', 'codon_display')
                                  ][, cell_total_fraction := codon_counts/sum(codon_counts), by = 'CB'
                                  ][, .(baseline_total_fraction = mean(cell_total_fraction)), by = 'codon_display'
                                  ]

site_baseline <- baseline_counts[
                                 #][CB %in% pull(allMeta %>% filter(type == 'RPF' & treatment == 'DMSO') %>% select(CB)),
                                 ][, cell_site_fraction := n/sum(n), by = c('CB', 'site')
                                 ][, .(baseline_site_fraction = mean(cell_site_fraction)), by = c('codon_display', 'site')
                                 ]

site_change <- codon_reads[
                           ][three != 'Ter',
                           ][, codon_display := paste(gsub('T', 'U', codon), three, sep = '_')
                           ][, .(n = .N), by = c('CB', 'site', 'codon_display')
                           ][, cell_site_fraction := n/sum(n), by = c('CB', 'site')
                           ][, cell_total_fraction := n/sum(n), by = 'CB'
                           ][allMeta[, c('CB', 'treatment')], on = 'CB', nomatch = NULL
                           ][total_baseline, on = 'codon_display', nomatch = NULL
                           ][site_baseline,  on = c('codon_display', 'site'), nomatch = NULL
                           ][, total_fc := cell_total_fraction/(baseline_total_fraction/nrow(site_order))
                           ][, site_fc := cell_site_fraction/baseline_site_fraction
                           ][site_order, on = 'site', nomatch = NULL
                           ][, site := fct_reorder(site, site_order)
                           ]


all_codons <- unique(site_change$codon_display)
all_codons <- all_codons[order(gsub(".*_(.*?)$","\\1",all_codons,perl=TRUE))]
all_codons <- all_codons[nchar(all_codons)==7]

# average fc across all sites to cluster codons
allHMa <- as.data.frame(dcast(site_change[
                              ][codon_display %in% all_codons,
                              ][site %in% c("em2", "em1", "e", "p", "a", "ap1", "ap2"),
                              ][, .(mean_site_fc = mean(site_fc), mean_total_fc = mean(total_fc)), by = c('CB', 'codon_display')
                              ],
                              CB ~ codon_display, value.var = 'mean_site_fc', fill = 0))
rownames(allHMa) <- allHMa$CB
allHMa <- allHMa[, colnames(allHMa) != 'CB']

allHMa_hc <- dendsort(hclust(dist(t(as.matrix(log2(allHMa))), method = 'euclidean'), method = 'ward.D2'), type = 'min')

all_codons <- all_codons[allHMa_hc$order]

allHM <- as.data.frame(dcast(site_change[
                             ][codon_display %in% all_codons,
                             ][, site_out := paste(codon_display, site, sep = '_'),
                             ],
                             CB ~ site_out, value.var = 'site_fc', fill = 0))
rownames(allHM) <- allHM$CB
hm_colnames <- c(t(outer(all_codons, pull(site_order %>% arrange(site_order) %>% select(site)), FUN = 'paste', sep = '_')))
hm_colnames <- hm_colnames[hm_colnames %in% colnames(allHM)]
allHM <- allHM[, hm_colnames]

selected_to_quant <- site_change[
                                 ][site %in% c('a'),
                                 ][, l2fc := log2(site_fc)
                                 ]
quantile(abs(selected_to_quant$l2fc), probs = 1-.01)
pausing_limit <- 0.5

allMeta <- as.data.frame(allMeta)
rownames(allMeta) <- allMeta$CB


pdf(paste0(figurePrefix, "_all_codon_stalling_sitefc_l2fc.pdf"), height = 3, width = 17)
site_labels <- HeatmapAnnotation(site = anno_mark(at = grep("^AUA_Ile", colnames(allHM)),
                                                  labels = c("-2", "-1", "E", "P", "A", "+1", "+2"),
                                                  side = 'bottom'), 
                                 which = "column" )
cell_annotation <- HeatmapAnnotation(
                                     treatment  = as.character(allMeta[rownames(allHM), 'treatment']),
                                     replicate  = as.character(allMeta[rownames(allHM), 'replicate']),
                                     col = list(
                                                treatment  = setNames(brewer.pal(length(unique(allMeta$treatment)), 'Set2')[1:length(unique(allMeta$treatment))], unique(allMeta$treatment)),
                                                replicate  = setNames(brewer.pal(length(unique(allMeta$replicate)), 'Pastel1')[1:length(unique(allMeta$replicate))], unique(allMeta$replicate))
                                                ),
                                       which = 'row')
Heatmap(log2(as.matrix((allHM))),
        bottom_annotation = site_labels,
        left_annotation = cell_annotation,
        col = colorRamp2(seq(-pausing_limit, pausing_limit, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        na_col = "grey0",
        name = "log2 fold change",
        cluster_rows = TRUE,
        clustering_distance_rows = 'euclidean',
        clustering_method_rows = 'ward.D2',
        column_title_gp = gpar(fontsize = 14),
        column_title_rot = 90,
        cluster_columns = FALSE,
        cluster_row_slices = TRUE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = factor(rep(all_codons, each = length(levels(site_change$site))), levels = all_codons, ordered = TRUE),
        row_title_rot = 0,
        border = FALSE,
        use_raster = FALSE,
        #raster_device = "CairoPNG",
        raster_quality = 10)
dev.off()




## transcriptional response

trnMeta <- testMeta[testMeta$type == 'VASA',] %>%
  mutate_if(is.factor, droplevels)

trnCount <- testCounts[, rownames(trnMeta)]

trnFull <- model.matrix( ~ treatment, data = trnMeta)

trnDds <- DESeqDataSetFromMatrix(countData = trnCount,
                                 colData   = trnMeta,
                                 design    = trnFull)
sizeFactors(trnDds) <- sumFactors[colnames(trnDds), 'SF']
trnDds <- DESeq(trnDds)

trnRes <- results(trnDds, independentFiltering = TRUE, format = 'DataFrame', alpha = 0.1) %>%
  data.frame() %>%
  rownames_to_column('gene_out')

plt_volcano_trn <- trnRes %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj)))+
  geom_pointdensity(size = 0.5)+
  scale_color_viridis()+
  geom_hline(yintercept = -log10(0.01), linetype = 'dashed', colour = 'grey70', alpha = .5)#+
save_plot(paste0(figurePrefix, '_trn_dE_volcano.pdf'), plt_volcano_trn, base_height = 6/1.5, base_width = 9/1.5)



trnResRank <- trnRes %>% 
  inner_join(ENSG_to_name, by = 'gene_out') %>%
  select(gene_name, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene_name) %>%
  summarize(stat = mean(stat)) %>%
  ungroup() %>%
  deframe()


trn_react_res <- fgsea(fgsea_reactome, stats = trnResRank)


pdf(paste0(figurePrefix, '_trn_dE_gsea_table.pdf'), width = 15, height = 10)
plotGseaTable(pathways = fgsea_reactome[c(trn_react_res %>% filter(padj < 0.01) %>% filter(NES > 0) %>% slice_max(NES, n=20) %>% select(pathway) %>% pull(),
                                          trn_react_res %>% filter(padj < 0.01) %>% filter(NES < 0) %>% slice_min(NES, n=20) %>% select(pathway) %>% pull() %>% rev())],
              stats = trnResRank,
              fgseaRes = trn_react_res,
              gseaParam = 1.0)
dev.off()


## example genes dE
exampleGenes <- c('ATF3', 'DDIT3', 'DDIT4', 'PPP1R15A')
# GADD34 is PPP1R15A
# SESN2 is on chr1 and KZ208904.1; reads likely align to both and thus multimap. If important then will re-map to major contigs only.

exampleOut <- ENSG_to_name %>%
  filter(gene_name %in% exampleGenes) %>%
  select(gene_out) %>%
  pull()


normLongCounts <- longCounts[
                             ][gene_id %in% pull(ENSG_to_name %>% filter(gene_name %in% exampleGenes) %>% select(gene_out)),
                             ][allMeta[, c('CB', 'type', 'treatment')], on = 'CB', nomatch = NULL
                             ][sumFactors, on = 'CB', nomatch = NULL
                             ][, gene_display := gsub("^.*?-(.*?)$", "\\1", gene_id, perl = TRUE)
                             ][trnRes, on = c('gene_id' = 'gene_out'), nomatch = NULL
                             ][, gene_display := paste0(gene_display, ' ', sprintf("padj: %1.2E", padj), ' ', sprintf("l2fc: %1.2f", log2FoldChange))
                             ][, dr := dense_rank(-log2FoldChange)
                             ][, gene_display := fct_reorder(gene_display, dr)
                             ][, dr := NULL
                             ][type == 'VASA', tr := 0,
                             ][type == 'RPF' , tr := 1,
                             ][, type := fct_reorder(type, tr)
                             ][, tr := NULL
                             ]

plt_check <- ggplot(normLongCounts[type == 'VASA'], aes(x=treatment, y=log2(counts/SF), colour = treatment))+
  geom_point(position = position_jitterdodge(jitter.width = 0.3, seed = 112), size = 1) +
  facet_wrap(~gene_display, scales = 'free_y')
save_plot(paste0(figurePrefix, '_plt_count_check.pdf'), plt_check, base_width = 8/1.5, base_height = 6/1.5)



avCounts <- longCounts[
                       ][gene_id %in% expressedGenes,
                       ][allMeta %>% select(CB, type, treatment), on = 'CB', nomatch = NULL
                       ][sumFactors, on = 'CB', nomatch = NULL
                       ][, .(mean_norm_counts = mean(counts/SF), sem_norm_counts = sd(counts/SF)/sqrt(.N)), by = c('gene_id', 'type', 'treatment')
                       ]

foldChange <- dcast(dcast(avCounts[,c('gene_id', 'type', 'treatment', 'mean_norm_counts')],
                    gene_id + type ~ treatment, value.var = 'mean_norm_counts')[
                                                                                ][, fc := WEE1i/DMSO
                                                                                ],
                    gene_id ~ type, value.var = 'fc')

plt_fc <- ggplot(foldChange, aes(x=log2(VASA), y=log2(RPF)))+
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'black', alpha = 1.0)+
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black', alpha = 1.0)+
  geom_pointdensity(size = 0.5)+
  scale_color_viridis()+
  xlab('log2 f.c. mean VASA WEE1i/DMSO')+
  ylab('log2 f.c. mean RPF WEE1i/DMSO')+
  coord_equal()
save_plot(paste0(figurePrefix, '_tsf_WEE1I-DMSO_fc.pdf'), plt_fc, base_width = 7, base_height = 7)

avCounts <- longCounts[
                       ][gene_id %in% expressedGenes,
                       ][allMeta %>% select(CB, type, treatment), on = 'CB', nomatch = NULL
                       ][rpfSumFactor, on = 'CB', nomatch = NULL
                       ][, .(mean_norm_counts = mean(counts/SF), sem_norm_counts = sd(counts/SF)/sqrt(.N)), by = c('gene_id', 'type', 'treatment')
                       ]

foldChange <- dcast(dcast(avCounts[,c('gene_id', 'type', 'treatment', 'mean_norm_counts')],
                    gene_id + type ~ treatment, value.var = 'mean_norm_counts')[
                                                                                ][, fc := WEE1i/DMSO
                                                                                ],
                    gene_id ~ type, value.var = 'fc')

plt_fc <- ggplot(foldChange, aes(x=log2(VASA), y=log2(RPF)))+
  geom_vline(xintercept = 0, linetype = 'dashed', colour = 'black', alpha = 1.0)+
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'black', alpha = 1.0)+
  geom_pointdensity(size = 0.5)+
  scale_color_viridis()+
  xlab('log2 f.c. mean VASA WEE1i/DMSO')+
  ylab('log2 f.c. mean RPF WEE1i/DMSO')+
  #geom_point()+
  coord_equal()
save_plot(paste0(figurePrefix, '_rsf_WEE1I-DMSO_fc.pdf'), plt_fc, base_width = 7, base_height = 7)
