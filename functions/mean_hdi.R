mean_hdi <- function(x, hdi = .95){
  y = mean(x)
  ymin = quantile(x, probs = (1 - hdi)/2)
  ymax = quantile(x, probs = 1 - (1 - hdi)/ 2)
  
  return(data.frame(
    y = y,
    ymin = ymin,
    ymax = ymax
  )
  )
}
