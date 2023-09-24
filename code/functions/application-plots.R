#'
#' Useful plot functions for the application
#' 
#' 

plotMSE_pw <- function(MSE.vec, MSE.ita18 = rep(0,37), t.points, name_dir)
{
  
  library(wesanderson)
  library(latex2exp)
  pal  <- wes_palette("Cavalcanti1")
  
  df <- data.frame(t.points, MSE.vec, MSE.ita18)
  x.ticks <- seq(-3,1,by=1)
  
  t.plot <- t.points
  
  ggplot(df) +
    geom_line(aes(x = t.plot, y = MSE.vec,
                  color = "black"), size=1) +
    geom_line(aes(x = t.plot, y = MSE.ita18,
                  color = "grey70"), size=1) +
    scale_color_manual(name="Model:",values=c("black","grey70"), labels = c("F-ITA18", "ITA18")) +
    scale_x_continuous(breaks=x.ticks, labels=10^x.ticks) +
    scale_y_continuous(breaks=seq(0.08,0.18, by=0.02)) +
    theme_bw() +
    labs(x="Period [s]", y="MSE", title='(a) Comparison of point-wise MSE') +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
          text = element_text(size = 22),
          axis.text.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 19))
  
  ggsave(filename = paste0("MSE_pointwise.pdf"),
         plot = last_plot(),
         width = 9,
         height = 5,
         units = "in",
         device = NULL,
         path = name_dir,
         scale = 1,
         limitsize = TRUE,
         dpi = 300)
}

plot_sigma <- function(sigma.t, sigma.ITA18, t.points, name_dir)
{
  library(wesanderson)
  pal <- wes_palette('Cavalcanti1')
  rect.col1 <- adjustcolor(pal[3], alpha=0.2)
  rect.col2 <- adjustcolor(pal[2], alpha=0.2)

  df <- data.frame(t.points, sigma.ITA18, sigma.t)
  x.ticks <- seq(-3,1,by=1)
  
  t.plot <- t.points
  
  ggplot(df) +
    geom_line(aes(x = t.plot, y = sigma.t,
                  color = "black"), size=1) +
    geom_line(aes(x = t.plot, y = sigma.ITA18,
                  color = "grey70"), size=1) +
    scale_color_manual(name="Model:",values=c("black","grey70"), labels = c("F-ITA18", "ITA18")) +
    scale_x_continuous(breaks=x.ticks, labels=10^x.ticks) +
    scale_y_continuous(breaks=seq(0.32,0.4, by=0.04)) +
    theme_bw() +
    labs(x="Period [s]", y=TeX(r'($\hat{\sigma}$)'), title='Comparison of residuals standard deviation') +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
          text = element_text(size = 22),
          axis.text.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 19)) +
    annotate("rect", xmin = t.points[21], xmax = t.points[37], ymin = 0.31, ymax = 0.43,
             alpha = .1,fill = pal[2])
  #annotate("rect", xmin = t.points[32], xmax = t.points[37], ymin = 0.31, ymax = 0.43,
  #         alpha = .2,fill = pal[2])
  
  ggsave(filename = paste0("sigma_comparison.pdf"),
         plot = last_plot(),
         width = 10,
         height = 5,
         units = "in",
         device = NULL,
         path = name_dir,
         scale = 1,
         limitsize = TRUE,
         dpi = 300)
  
}



