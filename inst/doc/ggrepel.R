## ----setup, echo=FALSE, results='hide', warning=FALSE, error=FALSE, message=FALSE, cache=FALSE----
library(knitr)
opts_chunk$set(
  cache = TRUE,
  autodep = TRUE,
  echo = FALSE,
  warning = FALSE,
  error = FALSE,
  message = FALSE,
  out.width = 700,
  fig.width = 12,
  fig.height = 8,
  dpi = 144,
  cache.path = "cache/ggrepel/",
  fig.path = "figures/ggrepel/",
  concordance = TRUE
)

## ----geom_text_repel, echo=TRUE------------------------------------------
library(ggrepel)
set.seed(42)
ggplot(mtcars) +
  geom_point(aes(wt, mpg), color = 'red') +
  geom_text_repel(aes(wt, mpg, label = rownames(mtcars))) +
  theme_classic(base_size = 16)

## ----geom_text_repel_options, echo=TRUE----------------------------------
set.seed(42)
ggplot(mtcars) +
  geom_point(aes(wt, mpg), color = 'red') +
  geom_text_repel(
    aes(
      wt, mpg,
      color = factor(cyl),
      label = rownames(mtcars)
    ),
    size = 5,
    fontface = 'bold',
    segment.color = 'red',
    box.padding = unit(0.3, 'lines'),
    force = 2,
    max.iter = 1e4,
    expand = TRUE
  ) +
  scale_color_discrete(name = 'cyl') +
  theme_classic(base_size = 16)

## ----geom_label_repel, echo=TRUE-----------------------------------------
set.seed(42)
ggplot(mtcars) +
  geom_point(aes(wt, mpg)) +
  geom_label_repel(
    aes(wt, mpg, fill = factor(cyl), label = rownames(mtcars)),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.25, "lines")
  ) +
  theme_classic(base_size = 16)

## ----session_info, echo=TRUE---------------------------------------------
sessionInfo()

