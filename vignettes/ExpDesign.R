---
  title: "A Gentle Guide to the Grammar of Graphics<br>with `ggplot2`"
subtitle: "Tampa R Users Meetup"
author: "Garrick Aden-Buie<br>@grrrck"
date: "2018-01-23<br><br>Follow along: **bit.ly/trug-ggplot2**"
output:
  xaringan::moon_reader:
  lib_dir: libs
css: ["default", "metropolis", "metropolis-fonts", "trug-ggplot2.css"]
nature:
  ratio: 16:9
highlightStyle: github
highlightLines: true
countIncrementalSlides: false
slideNumberFormat: "%current%"
editor_options:
  chunk_output_type: console
---

  class: fullscreen, inverse, top, center, text-white
background-image: url("images/letter-g.jpg")

.font150[**Brought to you by the letter...**]

```{r setup, include=FALSE}
knitr::opts_chunk$set(fig.width=4.25, fig.height=3.5, fig.retina=3,
                      message=FALSE, warning=FALSE, cache = TRUE,
                      autodep = TRUE, hiline=TRUE)
knitr::opts_hooks$set(fig.callout = function(options) {
  if (options$fig.callout) {
    options$echo <- FALSE
    options$out.height <- "99%"
    options$fig.width <- 16
    options$fig.height <- 8
  }
  options
})
hook_source <- knitr::knit_hooks$get('source')
knitr::knit_hooks$set(source = function(x, options) {
  if (!is.null(options$hiline) && options$hiline) {
    x <- stringr::str_replace(x, "^ ?(.+)\\s?#<<", "*\\1")
  }
  hook_source(x, options)
})
options(htmltools.dir.version = FALSE, width = 90)
as_table <- function(...) knitr::kable(..., format='html', digits = 3)
```

---
  layout: true
# Why *ggplot2*?
---
