# clean up plots
clean_plot <- function(...){
  clean_plot <- theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line.x = element_line(color = 'black'),
                      axis.line.y = element_line(color = 'black'),
                      legend.key = element_rect(fill = 'white'),
                      text = element_text(size = 15),
                      line = element_line(linewidth = 1),
                      axis.ticks = element_line(linewidth = 1),
                      ...)
  
  return(clean_plot)
}

