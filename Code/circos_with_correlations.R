library(ComplexHeatmap)
library(circlize)
library(viridis)
library(glue)

today = Sys.Date()
cell = 'B'
yl = 50
df = read.csv(glue('../Data/for_circos_{cell}.csv'), row.names = 'index')
df = df[, -grep('All', colnames(df))]
df = df[, c(6, 4, 5, 1:3, 7:ncol(df))]
row.names(df) = gsub(' Cluster ', '', row.names(df))
ids = unique(df$id)
df$id = plyr::mapvalues(df$id, from=ids, to=letters[1:length(ids)])
df[, c('Healthy', 'MPN', 'MF')] = df[, c('Healthy', 'MPN', 'MF')]*100
x_lens = as.numeric(table(df$id))
x = c()
w = c()
for (i in seq(length(ids))) {
  x = append(x, c(1:x_lens[i]))
  w = append(w, seq(0, x_lens[i], length = x_lens[i]))
}

sectors_ = df$id
col_Healthy = rgb(red = 237/255, green = 65/255, blue = 53/255, alpha = 0.3)
col_MPN = rgb(red = 37/255, green = 33/255, blue = 96/255, alpha = 0.3)
col_MF = rgb(red = 5/255, green = 207/255, blue = 194/255, alpha = 0.3)

tiff(glue("../Graphs/Circos/{today}_{cell}_circos_corr.tiff"), units="in", width=4, height=4, res=600)

circos.clear()
circos.par(start.degree = 90, gap.degree = 7, "track.height" = 0.15)
# col_fun1 = colorRamp2(c(-1, 0, 1), c("royalblue3", "white", "red3"))
col_fun1 = colorRamp2(seq(-1, 1, 0.2), rev(RColorBrewer::brewer.pal(11, 'Spectral')))
circos.heatmap.initialize(dplyr::select(df,contains("coef")), split = sectors_)
circos.heatmap(dplyr::select(df,contains("coef")), 
               col = col_fun1, rownames.side = "inside", rownames.cex = 0.6, dend.side = 'outside')

circos.par("track.height" = 0.55)
circos.track(ylim = c(0, yl), sectors = sectors_, bg.col = 'snow2', bg.border='snow1')

abundance_track = 4
for (i in seq(length(ids))){
  circos.track(ylim = c(0, yl), sectors = unique(sectors_)[i], panel.fun = function(x, y) {
  z = seq(0.5, x_lens[i], by=1)
  y = seq(5, yl, by = yl/5)
  sector_xlims = rep(x_lens[i], each=5)
  
  circos.segments(z, 0, z, yl, col = 'snow1')
  circos.segments(0, y, sector_xlims, y, col = 'snow1')
  }, track.index=abundance_track, bg.border='snow1')

}
circos.track(ylim = c(0, yl), sectors = sectors_,
             panel.fun = function(x, y) {
               y = df$MF[CELL_META$subset]
               y = y[CELL_META$row_order]
               circos.lines(seq_along(y) - 0.5, y, col = col_MF, area = TRUE, border = '#05cfc2')
             }, cell.padding = c(0.02, 0, 0.02, 0), track.index=abundance_track, bg.border='white')
circos.track(ylim = c(0, yl), sectors = sectors_,
             panel.fun = function(x, y) {
               y = df$MPN[CELL_META$subset]
               y = y[CELL_META$row_order]
               circos.lines(seq_along(y) - 0.5, y, col = col_MPN, area = TRUE, border='#252160')
             }, cell.padding = c(0.02, 0, 0.02, 0), track.index=abundance_track, bg.border='white')
circos.track(ylim = c(0, yl), sectors = sectors_,
             panel.fun = function(x, y) {
               y = df$Healthy[CELL_META$subset]
               y = y[CELL_META$row_order]
               circos.lines(seq_along(y) - 0.5, y, col = col_Healthy, area = TRUE, border='#ed4135')
             }, cell.padding = c(0.02, 0, 0.02, 0), track.index=abundance_track, bg.border='white')
circos.yaxis(sector.index = 'a', side='left', track.index=abundance_track, labels.cex = 0.5, 
             col = 'snow2', tick.length = 0.2,labels.col = 'snow4')

dev.off()
