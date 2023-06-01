library(ComplexHeatmap)
library(circlize)
library(glue)

cell = 'T'
yl = 25
df = read.csv('./Data/230120_for_circos_no_correlation.csv', row.names = 'X')
df$cell = gsub(' Cluster.*?$', '', rownames(df))
df = df[df$cell == cell, ]
df = df[c(2,1, 3:nrow(df)),]
df[, c('Healthy', 'MPN', 'MF')] = df[, c('Healthy', 'MPN', 'MF')]*100
max(df[, c('Healthy', 'MPN', 'MF')])
df$metacluster = 'a'
df$coef = 0
row.names(df) = gsub(' Cluster ', '', row.names(df))
ids = unique(df$metacluster)

x_lens = as.numeric(table(df$metacluster))
x = c(1:x_lens[1])
w = seq(0, x_lens[1], length = x_lens[1])

sectors_ = df$metacluster
col_Healthy = rgb(red = 237/255, green = 65/255, blue = 53/255, alpha = 0.3)
col_MPN = rgb(red = 37/255, green = 33/255, blue = 96/255, alpha = 0.3)
col_MF = rgb(red = 5/255, green = 207/255, blue = 194/255, alpha = 0.3)

tiff(glue("~/Documents/VIC/Hobbs flu/Graphs/Circos/{cell}_circos_no_corr.tiff"), units="in", width=4, height=4, res=500)

circos.clear()
circos.par(start.degree = 90, gap.degree = 0, "track.height" = 0.1)
col_fun1 = colorRamp2(c(-0.7, 0, 0.7), c("royalblue3", "white", "red3"))
circos.heatmap.initialize(df[, c('coef'), drop=FALSE], split = sectors_)
circos.heatmap(df[, c('coef'), drop=FALSE], col = col_fun1, rownames.side = "inside", rownames.cex = 0.6)
circos.par("track.height" = 0.7)

circos.track(ylim = c(0, yl), panel.fun = function(x, y) {
  x = seq(0.5, nrow(df)+ 0.5, by = 1)
  y = seq(0, yl, by = 5)
  
  circos.segments(x, 0, x, yl, col = 'snow1')
  circos.segments(0, y, nrow(df)+1, y, col = 'snow1')}, track.index=3, bg.border='snow1', bg.col = 'snow2')
  
circos.track(ylim = c(0, yl), sectors = sectors_,
             panel.fun = function(x, y) {
               y = df$MF[CELL_META$subset]
               y = y[CELL_META$row_order]
               circos.lines(seq_along(y) - 0.5, y, col = col_MF, area = TRUE, border = '#05cfc2')
             }, cell.padding = c(0.02, 0, 0.02, 0), track.index=3, bg.border='white')
circos.track(ylim = c(0, yl), sectors = sectors_,
             panel.fun = function(x, y) {
               y = df$MPN[CELL_META$subset]
               y = y[CELL_META$row_order]
               circos.lines(seq_along(y) - 0.5, y, col = col_MPN, area = TRUE, border='#252160')
             }, cell.padding = c(0.02, 0, 0.02, 0), track.index=3, bg.border='white')
circos.track(ylim = c(0, yl), sectors = sectors_,
             panel.fun = function(x, y) {
               y = df$Healthy[CELL_META$subset]
               y = y[CELL_META$row_order]
               circos.lines(seq_along(y) - 0.5, y, col = col_Healthy, area = TRUE, border='#ed4135')
             }, cell.padding = c(0.02, 0, 0.02, 0), track.index=3, bg.border='white')
circos.yaxis(sector.index = 'a', side='right', track.index=3, labels.cex = 0.5, col = 'snow2', tick.length = -0.1, 
             labels.col = 'snow4')

dev.off()
