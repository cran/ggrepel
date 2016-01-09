---
title: "ggrepel Usage Examples"
author: "Kamil Slowikowski"
date: "2016-01-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{ggrepel Usage Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# ggrepel

## Motivation

Some text labels overlap each other in plots created with [geom_text]:


```r
library(ggplot2)
ggplot(mtcars) +
  geom_point(aes(wt, mpg), color = 'red') +
  geom_text(aes(wt, mpg, label = rownames(mtcars))) +
  theme_classic(base_size = 16)
```

<img src="https://github.com/slowkow/ggrepel/blob/master/vignettes/figures/ggrepel/geom_text-1.png" title="" alt="" width="700" />

## Algorithm

`ggrepel` implements functions to repel overlapping text labels away from
each other and away from the data points that they label. The algorithm
works as follows:

- For each box:
    - Move the box into the allowed plotting area.
    - If the bounding box overlaps other boxes:
        - Repel the overlapping boxes from each other.
    - If the bounding box overlaps data points:
        - Repel the box away from the data points.
- Repeat until all overlaps are resolved, up to a preset limit
  of iterations.

## Usage Examples

### geom_text_repel

We can repel the text labels away from each other by loading `ggrepel` and
using `geom_text_repel` instead:


```r
library(ggrepel)
set.seed(42)
ggplot(mtcars) +
  geom_point(aes(wt, mpg), color = 'red') +
  geom_text_repel(aes(wt, mpg, label = rownames(mtcars))) +
  theme_classic(base_size = 16)
```

<img src="https://github.com/slowkow/ggrepel/blob/master/vignettes/figures/ggrepel/geom_text_repel-1.png" title="" alt="" width="700" />

#### Options

All options available for [geom_text] such as `size` and
`fontface` are also available for `geom_text_repel`.

However, the following parameters are not supported:

- `hjust`
- `vjust`
- `nudge_x`
- `nudge_y`
- `position`
- `check_overlap`

`ggrepel` provides additional parameters for `geom_text_repel` and `geom_label_repel`:

- `segment.color` is the line segment color
- `box.padding` is the padding surrounding the text bounding box
- `force` is the force of repulsion between overlapping text labels
- `max.iter` is the maximum number of iterations to attempt to resolve overlaps
- `expand` the text will be arranged in the expanded plot area if TRUE, or else
  the text will be arranged within the range of the data points


```r
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
```

<img src="https://github.com/slowkow/ggrepel/blob/master/vignettes/figures/ggrepel/geom_text_repel_options-1.png" title="" alt="" width="700" />

### geom_label_repel

`geom_label_repel` is based on [geom_label].


```r
set.seed(42)
ggplot(mtcars) +
  geom_point(aes(wt, mpg)) +
  geom_label_repel(
    aes(wt, mpg, fill = factor(cyl), label = rownames(mtcars)),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.25, "lines")
  ) +
  theme_classic(base_size = 16)
```

<img src="https://github.com/slowkow/ggrepel/blob/master/vignettes/figures/ggrepel/geom_label_repel-1.png" title="" alt="" width="700" />

## R Session Info


```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.10.5 (Yosemite)
## 
## locale:
## [1] C/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggrepel_0.3   ggplot2_2.0.0 knitr_1.12   
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.2      codetools_0.2-14 digest_0.6.9     plyr_1.8.3      
##  [5] grid_3.2.3       gtable_0.1.2     formatR_1.2.1    magrittr_1.5    
##  [9] evaluate_0.8     scales_0.3.0     stringi_1.0-1    rmarkdown_0.9.2 
## [13] labeling_0.3     tools_3.2.3      stringr_1.0.0    munsell_0.4.2   
## [17] yaml_2.1.13      colorspace_1.2-6 htmltools_0.3
```

[geom_text]: http://docs.ggplot2.org/current/geom_text.html
[geom_label]: http://docs.ggplot2.org/current/geom_text.html
