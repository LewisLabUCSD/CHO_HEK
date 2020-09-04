# ggsave2 <- function(file.name, path.name = 'Figures_CHOHEK/', height, width, plot = last_plot()) {
#   ggsave(plot = plot, paste0(path.name, file.name, '.png'),height = height, width = width, dpi = 150)
#   ggsave(plot = plot, paste0(path.name, file.name, '.pdf'),height = height, width = width)
# } ## moved to ppi_assist


trend.plot <- function(dt, compare = c('CHO', '293F'), text.displacement = 0,
                       limits=c(-5, 5)) {
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  dt <- dt[cellLine%in%compare]
  cellLine1 = compare[1]
  cellLine2 = compare[2]
  # cellLine1 <- dt$cellLine %>% unique() %>% sort %>% .[1]
  # cellLine2 <- dt$cellLine %>% unique() %>% sort %>% .[2]
  
  dt[, improvement.log := log2(.SD[cellLine==cellLine2 ,1+`titer (ug/ml)`] /.SD[cellLine==cellLine1, 1+`titer (ug/ml)`]), 
     by = 'sample_ID']
  dt[, improvement.tpm.log := log2(.SD[cellLine==cellLine2 ,1+TPM] /.SD[cellLine==cellLine1, 1+TPM]), 
     by = 'sample_ID']
  
  scaffold.plot <-   ggplot(dt, aes(1+TPM, 1+`titer (ug/ml)`,
                                    group = sample_ID, label = sample_ID,
                                    color = improvement.log)) + theme_bw()
  pg <- layer_data(scaffold.plot) %>% as.data.table()
  text.dt <-  pg[, .(sample_ID = label, 
                     x = gm_mean(x),
                     y = gm_mean(y)
  ),
  by = 'label'] %>% 
    merge(dt[, .(improvement.log = unique(improvement.log),
                 improvement.tpm.log = unique(improvement.tpm.log)), by = 'sample_ID'], by = 'sample_ID')
  ## controlling for dynamic range
  yrange <- text.dt[, log(max(y)/min(y))] 
  xrange <- text.dt[, log(max(x)/min(x))]
  text.dt[, text.rad:= atan(improvement.log/yrange /(improvement.tpm.log/xrange)) ]
  text.dt[is.nan(text.rad), text.rad:=0]
  text.dt[, text.degree:= text.rad* (180/pi)]
  text.dt[, x_d:= x / exp(text.displacement * sin(text.rad))]
  text.dt[, y_d:= y * exp(text.displacement * cos(text.rad))]
  
  return(scaffold.plot + 
           geom_path(aes(color = improvement.log),alpha = .5,
                     arrow = arrow(ends = 'last', length = unit(.07, 'in'))) + 
           scale_colour_gradient2(low = 'blue',high = 'red', mid = 'grey', midpoint = 0, 
                                  limits=limits,
                                  name = "Titer foldchange")+
           geom_text(data = text.dt, 
                     aes(x_d, y_d, color = improvement.log,
                         angle = text.degree)
           ) + 
           scale_x_log10() + scale_y_log10() +
           ggtitle(sprintf('Titer improvement from %s to %s', 
                           (dt$cellLine %>% unique() %>% sort %>% .[2]),
                           (dt$cellLine %>% unique() %>% sort %>% .[1]))) 
  )
  
}


mosaic.plot.simple <- function(total.dt, tmd.dt, override.sec.alpha=F, mosaic.filter=F) {
  library(ggmosaic)
  # total.dt <- transcriptome.usage.dt[, -c('ensgid', 'geneID', 'Subsystem')]
  mosaic.dt <- merge(total.dt, tmd.dt, by.x = 'humanSymbol', by.y = 'gene.name') %>% setnames('TPM.mean', 'tpm')
  mosaic.dt[, l1:=factor(l1, levels = tmd.dt$l1 %>% unique, ordered = T) %>% droplevels()]
  mosaic.dt[, l2:=factor(l2, levels = tmd.dt$l2 %>% unique, ordered = T) %>% droplevels()]
  mosaic.dt[, l3:=factor(l3, levels = tmd.dt$l3 %>% unique, ordered = T) %>% droplevels()]
  if(mosaic.filter) mosaic.dt <- mosaic.dt[,.SD[mean(tpm)>10],by=l2]
  
  text.color.dt <- mosaic.dt[, .(color = color.dict[l2], tpm.sum = sum(tpm, na.rm = T)), by = 'l2'] %>% 
    .[, .(l2, color, frac = tpm.sum/max(tpm.sum))] %>% 
    .[order(l2), adjustcolor(color, alpha.f = min(5*frac, 1)), by = 'l2']
  if(override.sec.alpha) text.color.dt[l2=='secretory pathway', V1:='#000000FF']
  
  mosaic.dt %>% ggplot() + ggmosaic::geom_mosaic(aes(weight = tpm, x = product(cellLine, l2),
                                                     fill = l2),
                                                 offset = 0) +
    scale_fill_manual(values = color.dict,
                      guide=FALSE) + 
    coord_flip() +
    theme_pubclean() %+replace%
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_text(color = text.color.dt$V1,
                                     hjust = 1),
          axis.text.x = element_blank(),
          axis.title.x=element_blank(),
          # axis.text.x=element_text(angle = 35, vjust = .8),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(angle = 0)
    ) 
}


fig.lfc.comparison <- function(all.res, n.sig = 20, h.density = c(.05, .05),
                               geom_point.alpha = .1) {
  
  plt.dt <- all.res[baseMean.CHO>0|baseMean.HEK>0]
  plt.dt[is.na(log2FoldChange.CHO), log2FoldChange.CHO:=0]
  plt.dt[is.na(log2FoldChange.HEK), log2FoldChange.HEK:=0]
  plt.dt[, lfc.r2:= log2FoldChange.HEK**2 + log2FoldChange.CHO**2]
  plt.dt[, lfc.diff.abs:= abs(log2FoldChange.HEK-log2FoldChange.CHO)]
  plt.dt[, lfc.diff.hyperbola:=  log2FoldChange.CHO *log2FoldChange.HEK]
  # plt.dt[rank(-lfc.r2)<20, group:='extreme']
  # plt.dt[rank(-lfc.diff.abs)<=n.sig, group:='differential']
  plt.dt[rank(lfc.diff.hyperbola)<=n.sig, group:='differential']
  plt.dt[ humanSymbol %in%secM.amir$humanSymbol, secretory.pathway:=T ]
  # plt.dt[rank(--lfc.r2)<20 & rank(-lfc.diff.abs)<20, group:='extreme differential']
  
  ggplot(all.res, aes(log2FoldChange.CHO,log2FoldChange.HEK, label = humanSymbol)) +
    geom_point(alpha = geom_point.alpha) + 
    stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", alpha = .3, h = h.density) + 
    scale_fill_viridis_c(guide=FALSE) + 
    scale_color_discrete(guide = F)+
    ggrepel::geom_label_repel(data = plt.dt[!is.na(group)], 
                              aes(log2FoldChange.CHO,log2FoldChange.HEK, label = humanSymbol, color = secretory.pathway)) + 
    # coord_fixed() + 
    xlab('activation LFC in CHO') +
    ylab('activation LFC in HEK')
}

corr.pheatmap.plot <- function(., n.genes = 50,
                               corr.method = 'pearson',
                               cluster_rows =F, 
                               cluster_cols = F,
                               pheatmat.sale = 'none', anno.secM = F){
  cor.dt <- .[, .( lfc.titer_correlation = cor(log2FoldChange, lfc.titer, method = corr.method)), by = 'symbol'] 
  cor.dt[, `titer correlation`:= ifelse(lfc.titer_correlation>0, 'positive', 'negative')]
  cor.dt[, `titer correlation`:= factor(`titer correlation`, levels = c( 'negative', 'positive'))]
  if(anno.secM){
    cor.dt <- merge(cor.dt, secM.amir[, .(symbol = humanSymbol, functionalGroup)], by = 'symbol', all.x = T)
    annotation_row <- cor.dt[, .(symbol, `secretory function`=functionalGroup,
                                 `correlation w/ titer` = `titer correlation`)] %>% dt.to.df()
  }else{
    annotation_row <- cor.dt[, .(symbol, `correlation w/ titer` = `titer correlation`)] %>% dt.to.df()
  }
  genes.to.plot <- cor.dt[ order(-abs(lfc.titer_correlation))][1:n.genes, symbol]
  ## reorder:
  
  symbol.order <- 
    cor.dt[symbol%in%genes.to.plot][order(-lfc.titer_correlation), symbol]
  sample_ID.order <- 
    .[, .(lfc.titer, sample_ID)] %>% unique %>% .[order(-lfc.titer), sample_ID]
  
  mat.data.m <- .[symbol%in%genes.to.plot] %>% 
    .[, symbol:=factor(symbol, levels = symbol.order)] %>% 
    .[, sample_ID:=factor(sample_ID, levels = sample_ID.order)]
  
  mat.data <-   dcast(mat.data.m,
                      symbol~sample_ID, value.var = 'log2FoldChange')
  # mat_breaks <- quantile_breaks(as.matrix(dt.to.df(.)), n = 11)
  pheatmap( dt.to.df(mat.data),
            cluster_rows = cluster_rows,
            cluster_cols = cluster_cols,
            annotation_col = transgene.count.yield.dt.lfc[, .(sample_ID, `titer LFC` = lfc.titer)] %>% dt.to.df(),
            annotation_row =annotation_row,
            # color = RColorBrewer::brewer.pal(length(mat_breaks)-1, 'RdYlBu'),
            # breaks            = mat_breaks,
            # annotation_colors = list(lfc.titer_correlation = 'Spectral'),
            border_color      = NA, 
            scale = pheatmat.sale
            
  )
}


fig.PCA <- function(vsd, to.color, label.name = 'sample_ID') {
  color.discrete = !(vsd[[to.color]] %>% is.numeric())
  
  plotPCA(vsd, intgroup = c('sample_ID', 'producer', 'cellLine', to.color), returnData = T) %>% {
    percentVar <- round(100 * attr(., "percentVar"))
    ggplot(., aes_string('PC1', 'PC2', color=to.color, label = label.name, shape = 'cellLine')) +
      geom_point(size=3) +
      # viridis::scale_color_viridis(option = 'inferno', discrete = color.discrete, na.value="grey")+
      ggrepel::geom_text_repel()+
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()
  }
}

fig.TSNE <- function(vsd, to.color) {
  color.discrete = !(vsd[[to.color]] %>% is.numeric())
  
  
  vsd %>% assay() %>% t %>% Rtsne::Rtsne(perplexity = 5) %>% {
    tsne.dt <- cbind(vsd.cho@colData %>% as.data.frame %>% as.data.table(keep.rownames = 'sample_ID_unique'),
                     data.table(TSNE1 = .$Y[,1], TSNE2 = .$Y[, 2]))
    
    tsne.dt %>% ggplot(aes_string('TSNE1', 'TSNE2', color=to.color, label = 'sample_ID', shape = 'cellLine')) +
      geom_point(size=3) +
      viridis::scale_color_viridis(option = 'inferno', discrete = color.discrete, na.value="grey")+
      ggrepel::geom_text_repel()
  }
  
}