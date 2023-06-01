library(viridis)
color.bar <- function(lut, min, max=-min, nticks=3, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  # dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

tiff("~/Documents/VIC/Hobbs flu/Graphs/fig4_colorbar.tiff", units="in", width=1.5, height=3, res=300)
new_pal = (RColorBrewer::brewer.pal(11, 'Spectral'))
color.bar(colorRampPalette(new_pal)(200), 1)
dev.off()
