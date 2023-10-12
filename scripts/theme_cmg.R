# theme_cmg()
theme_cmg <- function (base_size = 16, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.text = element_text(size = rel(0.8)),
          axis.ticks = element_line(colour = "black"),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.key.size = unit(1.2, "lines"),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "black", size = 1),
          panel.grid.major = element_line(colour = "grey90", size = 0.5),
          panel.grid.minor = element_line(colour = "grey98", size = 0.25),
          strip.background = element_rect(fill = "grey80", colour = "grey50"))
}
