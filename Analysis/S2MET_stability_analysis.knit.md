
<!-- rnb-text-begin -->

---
title: "S2MET Stability Analysis"
output: html_notebook
bibliography: C:/Users/Jeff/Documents/Literature/library.bib
---

This notebook will outline the stability analysis of the S2MET phenotype data


## Introduction

Load packages and set directories


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeSh0aWR5dmVyc2UpXG5gYGAifQ== -->

```r
library(tidyverse)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiXHUwMDFiWzMwbS0tIFx1MDAxYlsxbUF0dGFjaGluZyBwYWNrYWdlc1x1MDAxYlsyMm0gLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tIHRpZHl2ZXJzZSAxLjIuMCAtLVx1MDAxYlszOW1cblx1MDAxYlszMG1cdTAwMWJbMzJtdlx1MDAxYlszMG0gXHUwMDFiWzM0bWdncGxvdDJcdTAwMWJbMzBtIDIuMi4xICAgICBcdTAwMWJbMzJtdlx1MDAxYlszMG0gXHUwMDFiWzM0bXB1cnJyICBcdTAwMWJbMzBtIDAuMi4zXG5cdTAwMWJbMzJtdlx1MDAxYlszMG0gXHUwMDFiWzM0bXRpYmJsZSBcdTAwMWJbMzBtIDEuMy40ICAgICBcdTAwMWJbMzJtdlx1MDAxYlszMG0gXHUwMDFiWzM0bWRwbHlyICBcdTAwMWJbMzBtIDAuNy40XG5cdTAwMWJbMzJtdlx1MDAxYlszMG0gXHUwMDFiWzM0bXRpZHlyICBcdTAwMWJbMzBtIDAuNy4xICAgICBcdTAwMWJbMzJtdlx1MDAxYlszMG0gXHUwMDFiWzM0bXN0cmluZ3JcdTAwMWJbMzBtIDEuMi4wXG5cdTAwMWJbMzJtdlx1MDAxYlszMG0gXHUwMDFiWzM0bXJlYWRyICBcdTAwMWJbMzBtIDEuMS4xICAgICBcdTAwMWJbMzJtdlx1MDAxYlszMG0gXHUwMDFiWzM0bWZvcmNhdHNcdTAwMWJbMzBtIDAuMi4wXHUwMDFiWzM5bVxuXHUwMDFiWzMwbS0tIFx1MDAxYlsxbUNvbmZsaWN0c1x1MDAxYlsyMm0gLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tIHRpZHl2ZXJzZV9jb25mbGljdHMoKSAtLVxuXHUwMDFiWzMxbXhcdTAwMWJbMzBtIFx1MDAxYlszNG1kcGx5clx1MDAxYlszMG06Olx1MDAxYlszMm1maWx0ZXIoKVx1MDAxYlszMG0gbWFza3MgXHUwMDFiWzM0bXN0YXRzXHUwMDFiWzMwbTo6ZmlsdGVyKClcblx1MDAxYlszMW14XHUwMDFiWzMwbSBcdTAwMWJbMzRtZHBseXJcdTAwMWJbMzBtOjpcdTAwMWJbMzJtbGFnKClcdTAwMWJbMzBtICAgIG1hc2tzIFx1MDAxYlszNG1zdGF0c1x1MDAxYlszMG06OmxhZygpXHUwMDFiWzM5bVxuIn0= -->

```
[30m-- [1mAttaching packages[22m --------------------------------------- tidyverse 1.2.0 --[39m
[30m[32mv[30m [34mggplot2[30m 2.2.1     [32mv[30m [34mpurrr  [30m 0.2.3
[32mv[30m [34mtibble [30m 1.3.4     [32mv[30m [34mdplyr  [30m 0.7.4
[32mv[30m [34mtidyr  [30m 0.7.1     [32mv[30m [34mstringr[30m 1.2.0
[32mv[30m [34mreadr  [30m 1.1.1     [32mv[30m [34mforcats[30m 0.2.0[39m
[30m-- [1mConflicts[22m ------------------------------------------ tidyverse_conflicts() --
[31mx[30m [34mdplyr[30m::[32mfilter()[30m masks [34mstats[30m::filter()
[31mx[30m [34mdplyr[30m::[32mlag()[30m    masks [34mstats[30m::lag()[39m
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShyZWFkeGwpXG5saWJyYXJ5KGxtZTQpXG5gYGAifQ== -->

```r
library(readxl)
library(lme4)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiTG9hZGluZyByZXF1aXJlZCBwYWNrYWdlOiBNYXRyaXhcblxuQXR0YWNoaW5nIHBhY2thZ2U6IOOkvOO4sU1hdHJpeOOkvOO4slxuXG5UaGUgZm9sbG93aW5nIG9iamVjdCBpcyBtYXNrZWQgZnJvbSDjpLzjuLFwYWNrYWdlOnRpZHly46S847iyOlxuXG4gICAgZXhwYW5kXG4ifQ== -->

```
Loading required package: Matrix

Attaching package: 㤼㸱Matrix㤼㸲

The following object is masked from 㤼㸱package:tidyr㤼㸲:

    expand
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShicm9vbSlcbmxpYnJhcnkoc3RyaW5ncilcbmxpYnJhcnkoRlcpXG5gYGAifQ== -->

```r
library(broom)
library(stringr)
library(FW)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiVGhlIEZXIHBhY2thZ2UgcmVjZWl2ZWQgZmluYW5jaWFsIHN1cHBvcnQgZnJvbSBOSUggZ3JhbnRzIFIwMUdNMTAxMjE5IGFuZCBSMDFHTTA5OTk5MiBhbmQgZnJvbSBBcnZhbGlzXG4ifQ== -->

```
The FW package received financial support from NIH grants R01GM101219 and R01GM099992 and from Arvalis
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShnZ3JpZGdlcylcbmxpYnJhcnkoZ2dmb3JjZSlcbmxpYnJhcnkocGJyKVxuYGBgIn0= -->

```r
library(ggridges)
library(ggforce)
library(pbr)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiTG9hZGluZyBwYnI6IGFncmlkYXRcbiJ9 -->

```
Loading pbr: agridat
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBQcm9qZWN0IGFuZCBvdGhlciBkaXJlY3Rvcmllc1xuc291cmNlKFwiQzovVXNlcnMvSmVmZi9Hb29nbGUgRHJpdmUvQmFybGV5IExhYi9Qcm9qZWN0cy9TMk1FVF9NYXBwaW5nL3NvdXJjZS5SXCIpXG4jIFJlbW92ZSB0aGUgZW52aXJvbm1lbnRzIGluIHdoaWNoIHRoZSB2cCB3YXMgb25seSBvYnNlcnZlZFxuUzJfTUVUX0JMVUVzX3VzZSA8LSBTMl9NRVRfQkxVRXMgJT4lXG4gIGdyb3VwX2J5KHRyYWl0LCBlbnZpcm9ubWVudCkgJT4lXG4gIGZpbHRlcihzdW0obGluZV9uYW1lICVpbiUgdHApID4gMSkgJT4lXG4gIHVuZ3JvdXAoKVxuIyBOdW1iZXIgb2YgZW52aXJvbm1lbnRzIGFuZCBlbnRyaWVzXG5lbnYgPC0gcHVsbChkaXN0aW5jdChTMl9NRVRfQkxVRXNfdXNlLCBlbnZpcm9ubWVudCkpXG5uX2VudiA8LSBuX2Rpc3RpbmN0KGVudilcbm5fZW50cmllcyA8LSBuX2Rpc3RpbmN0KFMyX01FVF9CTFVFc191c2UkbGluZV9uYW1lKVxuIyBMb2FkIGZ3IGRhdGFcbmxvYWQoZmlsZS5wYXRoKHJlc3VsdF9kaXIsIFwiUzJNRVRfcGhlbm9fZndfcmVncmVzc2lvbl9yZXN1bHRzLlJEYXRhXCIpKVxuIyBMb2FkIHRoZSBtYXJrZXIgZWZmZWN0IEZXIHJlc3VsdHNcbmxvYWQoZmlsZS5wYXRoKHJlc3VsdF9kaXIsIFwiUzJNRVRfcGhlbm9fbWFyX2VmZl9md19yZXN1bHRzLlJEYXRhXCIpKVxuYGBgIn0= -->

```r
# Project and other directories
source("C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping/source.R")
# Remove the environments in which the vp was only observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  group_by(trait, environment) %>%
  filter(sum(line_name %in% tp) > 1) %>%
  ungroup()
# Number of environments and entries
env <- pull(distinct(S2_MET_BLUEs_use, environment))
n_env <- n_distinct(env)
n_entries <- n_distinct(S2_MET_BLUEs_use$line_name)
# Load fw data
load(file.path(result_dir, "S2MET_pheno_fw_regression_results.RData"))
# Load the marker effect FW results
load(file.path(result_dir, "S2MET_pheno_mar_eff_fw_results.RData"))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



### Finlay-Wilkinson Regression

Use the environment means to calculate the FW stability coefficient. The model for the regression coefficients looks like:

$$
y_{ij} = \mu + G_i + (1 + b_i)E_j + \epsilon_{ij},
$$

where $y_{ij}$ is the BLUE of the *i*th genotype in the *j*th environment, $\mu$ is the overall mean, $g_i$ is the effect of the *i*th genotype, $E_j$ is the mean of the *j*th environment, and $b_i$ is the regression coefficient of the response on the mean of the *j*th environment, and $\epsilon_{ij}$ is the random error.

To fit FW, we first fit a regular quantitative genetic model:

$$
y_{ij} = \mu + G_i + E_j + \epsilon_{ij}
$$

to obtain the environmental effects $E_j$. We then fit the following model for every genotype

$$
y_j = \mu + (1 + b)E_j + \epsilon_{ij}
$$

## Results

### Phenotypic Stability

#### Summary

Ranges and variances


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-frame-begin eyJtZXRhZGF0YSI6eyJjbGFzc2VzIjpbImdyb3VwZWRfZGYiLCJ0YmxfZGYiLCJ0YmwiLCJkYXRhLmZyYW1lIl0sIm5jb2wiOjMsIm5yb3ciOjZ9LCJyZGYiOiJINHNJQUFBQUFBQUFCZ3R5aVREaWl1QmlZR0JnWW1CbVlXSmdZZ1l5Z1RRREF5TURDd01uaUs1Z1lHQVdCakpBTWdKQXpBYVY0SEl2U3N6TWk4eE16VW5CS2NMdGtacVlrcG1YN3BKWWtvcEhLQ0FuTWEvRUl6VXpQYU1FbXhDcXZZeEpVQVpyU21wT1NTSkpvZ3dNZkNDVDdLZHZmaThvb05Sa3YzdnlhK2JTTmgvN1dSeDdacnJtbkxZL3hMSDk5NU82WHZ1NXFoTlhaVS9vc3Q5citWMUh0UHdMV3BDd0p1Y2tGaGREblFVVDVFcEpMRW5VU3l0S3pFMUZWNTRIRklNcFo0WUpsZ0REQ3VaYnZ1S1N4S1RNbk15U3l2aVMxS0pjcUNoellYSUttbEdjUmZubGVqRGplRUd4MWdBay92Ly8vd3RJL1VOVHpKNWZVSktabndkVXlpUU1EajlVQnpNV2dWd05Gb0FBWkVuTzFMd3l2YUxFdlBSVU5FVzhVRFliRXBzWllqSHpmN1NZWWt2TlM4L01nOFV5YTJwUlVYNFJqSk9UbUpTYUEvTSswRmRnVCtrVkZHWG13UUtGQ3loYXJGZVNYNUlJVThlVm5KOERFd0Y3bCtFZkFGZXFTUFRBQWdBQSJ9 -->

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["trait"],"name":[1],"type":["chr"],"align":["left"]},{"label":["stability_term"],"name":[2],"type":["chr"],"align":["left"]},{"label":["qcd"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"GrainYield","2":"b","3":"0.02314733"},{"1":"GrainYield","2":"delta","3":"0.10772580"},{"1":"HeadingDate","2":"b","3":"0.02542395"},{"1":"HeadingDate","2":"delta","3":"0.14089107"},{"1":"PlantHeight","2":"b","3":"0.02846363"},{"1":"PlantHeight","2":"delta","3":"0.11416573"}],"options":{"columns":{"min":{},"max":[10],"total":[3]},"rows":{"min":[10],"max":[10],"total":[6]},"pages":{}}}
  </script>
</div>

<!-- rnb-frame-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




Plot the FW regression (regular version)

This is typical FW regression (i.e. using the phenotypic means)


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG4jIFBsb3QgaGlzdG9ncmFtcyBvZiB0aGUgcmVncmVzc2lvbiBjb2VmZmljaWVudHMgYW5kIE1TRXNcbmdfcGhlbm9fZndfZGVucyA8LSBTMl9NRVRfcGhlbm9fZncgJT4lXG4gIGdncGxvdChhZXMoeCA9IGVzdGltYXRlKSkgKyBcbiAgZ2VvbV9kZW5zaXR5KGFlcyhmaWxsID0gXCJibHVlXCIpKSArIFxuICBmYWNldF93cmFwKHRyYWl0IH4gc3RhYmlsaXR5X3Rlcm0sIG5jb2wgPSAyLCBzY2FsZXMgPSBcImZyZWVcIikgK1xuICB4bGFiKFwiRXN0aW1hdGVcIikgK1xuICBsYWJzKHRpdGxlID0gXCJTdGFiaWxpdHkgTWVhc3VyZXNcIikgK1xuICBzY2FsZV9maWxsX2JyZXdlcihwYWxldHRlID0gXCJTZXQyXCIsIGd1aWRlID0gRkFMU0UpICtcbiAgdGhlbWVfYncoKSArXG4gIHRoZW1lKGF4aXMudGl0bGUueSA9IGVsZW1lbnRfYmxhbmsoKSlcblxuIyBKdXN0IHBsb3QgbGluZWFyIHN0YWJpbGl0eVxuZ19waGVub19md19kZW5zX2IgPC0gUzJfTUVUX3BoZW5vX2Z3ICU+JVxuICBmaWx0ZXIoc3RhYmlsaXR5X3Rlcm0gPT0gXCJiXCIpICU+JVxuICBnZ3Bsb3QoYWVzKHggPSBlc3RpbWF0ZSkpICsgXG4gIGdlb21fZGVuc2l0eShhZXMoZmlsbCA9IFwiYmx1ZVwiKSkgKyBcbiAgZmFjZXRfd3JhcCggfiB0cmFpdCwgbmNvbCA9IDEpICtcbiAgeGxhYihcIkVzdGltYXRlXCIpICtcbiAgbGFicyh0aXRsZSA9IFwiTGluZWFyIFN0YWJpbGl0eVwiKSArXG4gIHNjYWxlX2ZpbGxfYnJld2VyKHBhbGV0dGUgPSBcIlNldDJcIiwgZ3VpZGUgPSBGQUxTRSkgK1xuICB0aGVtZV9idygpICtcbiAgdGhlbWUoYXhpcy50aXRsZS55ID0gZWxlbWVudF9ibGFuaygpKVxuXG4jIEp1c3QgcGxvdCBub24tbGluZWFyIHN0YWJpbGl0eVxuZ19waGVub19md19kZW5zX2RlbHRhIDwtIFMyX01FVF9waGVub19mdyAlPiVcbiAgZmlsdGVyKHN0YWJpbGl0eV90ZXJtID09IFwiZGVsdGFcIikgJT4lXG4gIGdncGxvdChhZXMoeCA9IGVzdGltYXRlKSkgKyBcbiAgZ2VvbV9kZW5zaXR5KGFlcyhmaWxsID0gXCJibHVlXCIpKSArIFxuICBmYWNldF93cmFwKCB+IHRyYWl0LCBuY29sID0gMSwgc2NhbGVzID0gXCJmcmVlXCIpICtcbiAgeGxhYihcIkVzdGltYXRlXCIpICtcbiAgbGFicyh0aXRsZSA9IFwiTm9uLUxpbmVhciBTdGFiaWxpdHlcIikgK1xuICBzY2FsZV9maWxsX2JyZXdlcihwYWxldHRlID0gXCJTZXQyXCIsIGd1aWRlID0gRkFMU0UpICtcbiAgdGhlbWVfYncoKSArXG4gIHRoZW1lKGF4aXMudGl0bGUueSA9IGVsZW1lbnRfYmxhbmsoKSlcblxuIyBDb252ZXJ0IHRvIGdyb2JzXG5ncm9iX2xpc3QgPC0gbGlzdChiID0gZ19waGVub19md19kZW5zX2IsIGRlbHRhID0gZ19waGVub19md19kZW5zX2RlbHRhKSAlPiVcbiAgbWFwKGdncGxvdEdyb2IpXG4jIENvd3Bsb3RcbmdfZndfZGlzdCA8LSBjb3dwbG90OjpwbG90X2dyaWQocGxvdGxpc3QgPSBncm9iX2xpc3QpXG5cbnNhdmVfZmlsZSA8LSBmaWxlLnBhdGgoZmlnX2RpciwgXCJzdGFiaWxpdHlfZXN0aW1hdGVfZGlzdHJpdWJ0aW9ucy5qcGdcIilcbmdnc2F2ZShmaWxlbmFtZSA9IHNhdmVfZmlsZSwgcGxvdCA9IGdfZndfZGlzdCwgaGVpZ2h0ID0gNS42LCB3aWR0aCA9IDUsIGRwaSA9IDUwMClcblxuYGBgIn0= -->

```r

# Plot histograms of the regression coefficients and MSEs
g_pheno_fw_dens <- S2_MET_pheno_fw %>%
  ggplot(aes(x = estimate)) + 
  geom_density(aes(fill = "blue")) + 
  facet_wrap(trait ~ stability_term, ncol = 2, scales = "free") +
  xlab("Estimate") +
  labs(title = "Stability Measures") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_blank())

# Just plot linear stability
g_pheno_fw_dens_b <- S2_MET_pheno_fw %>%
  filter(stability_term == "b") %>%
  ggplot(aes(x = estimate)) + 
  geom_density(aes(fill = "blue")) + 
  facet_wrap( ~ trait, ncol = 1) +
  xlab("Estimate") +
  labs(title = "Linear Stability") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_blank())

# Just plot non-linear stability
g_pheno_fw_dens_delta <- S2_MET_pheno_fw %>%
  filter(stability_term == "delta") %>%
  ggplot(aes(x = estimate)) + 
  geom_density(aes(fill = "blue")) + 
  facet_wrap( ~ trait, ncol = 1, scales = "free") +
  xlab("Estimate") +
  labs(title = "Non-Linear Stability") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_blank())

# Convert to grobs
grob_list <- list(b = g_pheno_fw_dens_b, delta = g_pheno_fw_dens_delta) %>%
  map(ggplotGrob)
# Cowplot
g_fw_dist <- cowplot::plot_grid(plotlist = grob_list)

save_file <- file.path(fig_dir, "stability_estimate_distriubtions.jpg")
ggsave(filename = save_file, plot = g_fw_dist, height = 5.6, width = 5, dpi = 500)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



Plot the stability estimates as lines against g_i vs h_j


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuclxuXG5gYGAifQ== -->

```r
r

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiRXJyb3I6IG9iamVjdCAncicgbm90IGZvdW5kXG4ifQ== -->

```
Error: object 'r' not found
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



Plot the relationship between the genotypic effect and the sensitivity


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBQbG90IHRoZSBsaW5lYXIgYW5kIG5vbi1saW5lYXIgc3RhYmlsaXR5IHRvZ2V0aGVyXG5nX3N0YWJpbGl0eV9hbmRfbWVhbiA8LSBzdGFiaWxpdHlfbWVhbl9jb3JyICU+JVxuICBnZ3Bsb3QoYWVzKHggPSBnLCB5ID0gZXN0aW1hdGUsIGNvbCA9IHByb2dyYW0sIGdyb3VwID0gRkFMU0UpKSArXG4gIGZhY2V0X3dyYXAodHJhaXQgfiBzdGFiaWxpdHlfdGVybSwgc2NhbGVzID0gXCJmcmVlXCIsIG5jb2wgPSAyKSArXG4gIGxhYnModGl0bGUgPSBcIlJlbGF0aW9uc2hpcCBCZXR3ZWVuIEdlbm90eXBlIE1lYW4gYW5kIFNlbnNpdGl2aXR5XCIpICtcbiAgZ19BZGRcblxuYGBgIn0= -->

```r
# Plot the linear and non-linear stability together
g_stability_and_mean <- stability_mean_corr %>%
  ggplot(aes(x = g, y = estimate, col = program, group = FALSE)) +
  facet_wrap(trait ~ stability_term, scales = "free", ncol = 2) +
  labs(title = "Relationship Between Genotype Mean and Sensitivity") +
  g_Add

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiRXJyb3I6IG9iamVjdCAnZ19BZGQnIG5vdCBmb3VuZFxuIn0= -->

```
Error: object 'g_Add' not found
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



What is the repeatability of the estimates of linear and non-linear stability?

Use a jackknife resampling procedure to estimate the variation around these estimates


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBHcm91cCBieSB0cmFpdFxucGhlbm9zX3RvbW9kZWwgJT4lIFxuICBncm91cF9ieSh0cmFpdCkgJT4lIFxuICBkbyh7XG4gICAgIyBFeHRyYWN0IGRhdGFcbiAgICBkZiA8LSAuIFxuICAgIFxuICAgICMgVmVjdG9yIG9mIGVudmlyb25tZW50IFxuICAgIGVudnMgPC0gdW5pcXVlKGRmJGVudmlyb25tZW50KVxuICAgIFxuICAgICMgSXRlcmF0ZSBvdmVyIGVudmlyb25tZW50cywgcmVtb3ZlIGRhdGEgZnJvbSB0aGF0IGVudmlyb25tZW50LCBhbmQgZml0IHRoZSBGVyByZWdyZXNzaW9uXG4gICAgZndfZml0cyA8LSBlbnZzICU+JVxuICAgICAgbWFwKH5maWx0ZXIoZGYsIGVudmlyb25tZW50ICE9IC4pICU+JVxuICAgICAgICAgICBncm91cF9ieShsaW5lX25hbWUpICU+JSBcbiAgICAgICAgICAgZG8oZml0ID0gbG0odmFsdWUgfiBoLCAuKSkgJT4lIFxuICAgICAgICAgICB1bmdyb3VwKCkgJT4lIFxuICAgICAgICAgICBtdXRhdGUoYiA9IG1hcF9kYmwoZml0LCB+Y29lZiguKVsyXSksIFxuICAgICAgICAgICAgICAgICAgZGVsdGEgPSBtYXBfZGJsKGZpdCwgfm1lYW4ocmVzaWQoLileMikpKSApXG4gICAgXG4gICAgIyBJdGVyYXRlIG92ZXIgbGluZSBuYW1lcyBhbmQgdGVzdCBlbnZpcm9ubWVudHNcbiAgICBmd19jb2VmIDwtIGxpc3QoZndfZml0cywgZW52cykgJT4lIFxuICAgICAgcG1hcF9kZih+bXV0YXRlKC54LCB0ZXN0X2VudiA9IC55KSlcbiAgICBcbiAgICAjIEFOT1ZBIG9mIGNvZWZmaWNpZW50c1xuICAgIGZ3X2NvZWZfcmVwX2ZpdCA8LSBmd19jb2VmICU+JSBcbiAgICAgIGdhdGhlcihjb2VmLCBlc3RpbWF0ZSwgLWxpbmVfbmFtZSwgLWZpdCwgLXRlc3RfZW52KSAlPiUgXG4gICAgICBncm91cF9ieShjb2VmKSAlPiUgXG4gICAgICBkbyhmaXQgPSBsbShlc3RpbWF0ZSB+IGxpbmVfbmFtZSwgLikpXG4gICAgXG4gICAgZndfY29lZl9yZXAgJT4lIFxuICAgICAgdW5ncm91cCgpICU+JSBcbiAgICAgIG11dGF0ZShhbm92YSA9IG1hcChmaXQsIH50aWR5KGFub3ZhKC4pKSkpICU+JSBcbiAgICAgIHVubmVzdChhbm92YSkgJT4lIFxuICAgICAgZ3JvdXBfYnkoY29lZikgJT4lIFxuICAgICAgc3VtbWFyaXplKHJlcGVhdGFiaWxpdHkgPSBtZWFuc3FbMV0gLyBzdW0obWVhbnNxKSkgfSlcblxuYGBgIn0= -->

```r
# Group by trait
phenos_tomodel %>% 
  group_by(trait) %>% 
  do({
    # Extract data
    df <- . 
    
    # Vector of environment 
    envs <- unique(df$environment)
    
    # Iterate over environments, remove data from that environment, and fit the FW regression
    fw_fits <- envs %>%
      map(~filter(df, environment != .) %>%
           group_by(line_name) %>% 
           do(fit = lm(value ~ h, .)) %>% 
           ungroup() %>% 
           mutate(b = map_dbl(fit, ~coef(.)[2]), 
                  delta = map_dbl(fit, ~mean(resid(.)^2))) )
    
    # Iterate over line names and test environments
    fw_coef <- list(fw_fits, envs) %>% 
      pmap_df(~mutate(.x, test_env = .y))
    
    # ANOVA of coefficients
    fw_coef_rep_fit <- fw_coef %>% 
      gather(coef, estimate, -line_name, -fit, -test_env) %>% 
      group_by(coef) %>% 
      do(fit = lm(estimate ~ line_name, .))
    
    fw_coef_rep %>% 
      ungroup() %>% 
      mutate(anova = map(fit, ~tidy(anova(.)))) %>% 
      unnest(anova) %>% 
      group_by(coef) %>% 
      summarize(repeatability = meansq[1] / sum(meansq)) })

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiXG58PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PSAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB8IDMzJSB+MjQgcyByZW1haW5pbmcgICAgXG58PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB8IDY3JSB+MTIgcyByZW1haW5pbmcgICAgXG58PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT09PT18MTAwJSB+MCBzIHJlbWFpbmluZyAgICAgXG4ifQ== -->

```

|==================================                                                                      | 33% ~24 s remaining    
|=====================================================================                                   | 67% ~12 s remaining    
|========================================================================================================|100% ~0 s remaining     
```



<!-- rnb-output-end -->

<!-- rnb-frame-begin eyJtZXRhZGF0YSI6eyJjbGFzc2VzIjpbImdyb3VwZWRfZGYiLCJ0YmxfZGYiLCJ0YmwiLCJkYXRhLmZyYW1lIl0sIm5jb2wiOjMsIm5yb3ciOjZ9LCJyZGYiOiJINHNJQUFBQUFBQUFCcFZSVFVyRFFCUit5U1NWUmxLRTNpT2dGK2ltWUZjaXJ1cHltcnlrWThlWk1CbFFkNTdDbFJ0UDRqMDhnOUNkU0JlTkwya0cydUFQRHN5OGI3NzM4WDZ2cHZPemFCNEJnQThzOE1GbkJNa0NlQkRBc0xIM0FHeE1vUEdjMEIxMGp1amNjS0d1QmNyc1IrWjRoandUcXBoeWk3OVFsNUlyTzBOUkxPMTMxR0ZlYjlHQk1FTnArYjlZZ0ZFVGFiTCtlSDVsVHkrVDllYno5TzM5NHE5L2J5UmhLbmxWZFdVNU1zcTQ1VWx1K0MzMjVZbzRKMmVPdERRcjEyMlFhc3c3SEJzc2tVSXRoQlQyb1JkcWFQUmQ0c0xGemRZZTZhbnJla05tMnhNZjZkSUtyVWpxajl2NUhSYnNtYWJxbHRpZGZlZW9vZ3FTRzU2dVZrcmsyRlBHSFI3c1liYkx6dXJldWdhb0NxSGNxa00wUmh2M2tYeUIwbVdrMXRyT2t0SUk1U1lURVZzbFZsdnVkRkdxcFdQYW5tSDdCVlIveGh6RkFnQUEifQ== -->

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["trait"],"name":[1],"type":["chr"],"align":["left"]},{"label":["coef"],"name":[2],"type":["chr"],"align":["left"]},{"label":["repeatability"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"GrainYield","2":"b","3":"0.9988540"},{"1":"GrainYield","2":"delta","3":"0.9993859"},{"1":"HeadingDate","2":"b","3":"0.9988540"},{"1":"HeadingDate","2":"delta","3":"0.9993859"},{"1":"PlantHeight","2":"b","3":"0.9988540"},{"1":"PlantHeight","2":"delta","3":"0.9993859"}],"options":{"columns":{"min":{},"max":[10],"total":[3]},"rows":{"min":[10],"max":[10],"total":[6]},"pages":{}}}
  </script>
</div>

<!-- rnb-frame-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- Plot stability across ECs -->

<!-- ```{r plot.fw.ec} -->

<!-- # Load fw data -->
<!-- load(file.path(result_dir, "S2MET_ec_fw_regression_results.RData")) -->

<!-- # Combine the populations into one df -->
<!-- S2_MET_ec_oneyear_fw_comb <- S2_MET_ec_oneyear_fw %>% -->
<!--   list(., names(.)) %>% -->
<!--   pmap_df(~mutate(.x, population = .y)) -->

<!-- S2_MET_ec_multiyear_fw_Comb <- S2_MET_ec_multiyear_fw %>% -->
<!--   list(., names(.)) %>% -->
<!--   pmap_df(~mutate(.x, population = .y)) -->

<!-- ## Plot the value of each genotype in each environment over the environmental -->
<!-- ## gradient based on ECs -->
<!-- g_ec_oneyear_fw_tp <- S2_MET_ec_oneyear_fw_comb %>% -->
<!--   filter(stability_term == "b", population == "pop_tp") %>% -->
<!--   ggplot(aes(x = h_ec, y = value, col = program, group = line_name)) + -->
<!--   geom_point() + -->
<!--   geom_abline(aes(slope = estimate, intercept = g, col = program), alpha = 0.15) + -->
<!--   facet_grid_paginate(trait ~ variable, nrow = 3, ncol = 5, page = 1, scales = "free") + -->
<!--   theme_bw() -->

<!-- ## Iterate through the pages and save the graphs -->
<!-- for (i in seq(n_pages(g_ec_oneyear_fw_tp))) { -->
<!--   save_file <- file.path(fig_dir, str_c("oneyear_ec_fw_tp", i, ".jpg")) -->
<!--   g <- g_ec_oneyear_fw_tp + facet_grid_paginate(trait ~ variable, nrow = 3, ncol = 5, page = i, scales = "free") -->
<!--   ggsave(filename = save_file, plot = g, width = 9, height = 7) -->
<!-- } -->


<!-- # Plot using all data -->
<!-- g_ec_oneyear_fw_all <- S2_MET_ec_oneyear_fw_comb %>% -->
<!--   filter(stability_term == "b", population == "pop_all") %>% -->
<!--   ggplot(aes(x = h_ec, y = value, col = program, group = line_name)) + -->
<!--   geom_point() + -->
<!--   geom_abline(aes(slope = estimate, intercept = g, col = program), alpha = 0.15) + -->
<!--   facet_grid_paginate(trait ~ variable, nrow = 3, ncol = 5, page = 1, scales = "free") + -->
<!--   theme_bw() -->

<!-- ## Iterate through the pages and save the graphs -->
<!-- for (i in seq(n_pages(g_ec_oneyear_fw_all))) { -->
<!--   save_file <- file.path(fig_dir, str_c("oneyear_ec_fw_all", i, ".jpg")) -->
<!--   g <- g_ec_oneyear_fw_all + facet_grid_paginate(trait ~ variable, nrow = 3, ncol = 5, page = i, scales = "free") -->
<!--   ggsave(filename = save_file, plot = g, width = 9, height = 7) -->
<!-- } -->





<!-- ## Multi-year ECs -->
<!-- g_ec_multiyear_fw_tp <- S2_MET_ec_multiyear_fw_Comb %>% -->
<!--   filter(stability_term == "b", population == "pop_tp") %>% -->
<!--   ggplot(aes(x = h_ec, y = value, col = program, group = line_name)) + -->
<!--   geom_point() + -->
<!--   geom_abline(aes(slope = estimate, intercept = g, col = program), alpha = 0.15) + -->
<!--   facet_grid_paginate(trait ~ variable, nrow = 3, ncol = 5, page = 1, scales = "free") + -->
<!--   theme_bw() -->

<!-- ## Iterate through the pages and save the graphs -->
<!-- for (i in seq(n_pages(g_ec_multiyear_fw_tp))) { -->
<!--   save_file <- file.path(fig_dir, str_c("multiyear_ec_fw_tp", i, ".jpg")) -->
<!--   g <- g_ec_multiyear_fw_tp + facet_grid_paginate(trait ~ variable, nrow = 3, ncol = 5, page = i, scales = "free") -->
<!--   ggsave(filename = save_file, plot = g, width = 9, height = 7) -->
<!-- } -->


<!-- # Plot using all data -->
<!-- g_ec_multiyear_fw_all <- S2_MET_ec_multiyear_fw_Comb %>% -->
<!--   filter(stability_term == "b", population == "pop_all") %>% -->
<!--   ggplot(aes(x = h_ec, y = value, col = program, group = line_name)) + -->
<!--   geom_point() + -->
<!--   geom_abline(aes(slope = estimate, intercept = g, col = program), alpha = 0.15) + -->
<!--   facet_grid_paginate(trait ~ variable, nrow = 3, ncol = 5, page = 1, scales = "free") + -->
<!--   theme_bw() -->

<!-- ## Iterate through the pages and save the graphs -->
<!-- for (i in seq(n_pages(g_ec_multiyear_fw_all))) { -->
<!--   save_file <- file.path(fig_dir, str_c("multiyear_ec_fw_all", i, ".jpg")) -->
<!--   g <- g_ec_multiyear_fw_all + facet_grid_paginate(trait ~ variable, nrow = 3, ncol = 5, page = i, scales = "free") -->
<!--   ggsave(filename = save_file, plot = g, width = 9, height = 7) -->
<!-- } -->



<!-- ``` -->


### Marker Effect Stability

Visualize the correlation across environments of marker effects


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG4jIEZpcnN0IHNldCB1cCB0d28td2F5IHRhYmxlIG9mIG1hcmtlciBlZmZlY3RzIGFuZCBlbnZpcm9ubWVudHNcblMyX01FVF9tYXJrZXJfZWZmX3RhYmxlIDwtIFMyX01FVF9tYXJrZXJfZWZmX3BoZW5vX2Z3ICU+JSBcbiAgc2VsZWN0KC4sIHRyYWl0LCBtYXJrZXIsIGVudmlyb25tZW50LCBtYXJfZWZmZWN0KSAlPiUgXG4gIGRpc3RpbmN0KCkgJT4lIFxuICBzcHJlYWQoZW52aXJvbm1lbnQsIG1hcl9lZmZlY3QpXG5cbiMgQ2FsY3VsYXRlIHRoZSBwYWlyd2lzZSBjb3JyZWxhdGlvbiBvZiBtYXJrZXIgZWZmZWN0cyBiZXR3ZWVuIGVudmlyb25tZW50c1xuUzJfTUVUX21hcmtlcl9lZmZfY29yciA8LSBTMl9NRVRfbWFya2VyX2VmZl90YWJsZSAlPiVcbiAgZ3JvdXBfYnkoLiwgdHJhaXQpICU+JSBcbiAgZG8oe1xuICAgIGRmIDwtIC5cbiAgICBzZWxlY3QoZGYsIC10cmFpdCwgLW1hcmtlcikgJT4lIFxuICAgICAgY29yKCkgJT4lIFxuICAgICAgYXMuZGF0YS5mcmFtZSgpICU+JSBcbiAgICAgIHJvd25hbWVzX3RvX2NvbHVtbihcImVudmlyb25tZW50MVwiKSAlPiUgXG4gICAgICBnYXRoZXIoZW52aXJvbm1lbnQyLCBjb3JyLCAtZW52aXJvbm1lbnQxKSB9KVxuXG4jIEZpbmQgdGhlIG1pbmltdW0gYW5kIG1heGltdW0gY29ycmVsYXRpb25cblMyX01FVF9tYXJrZXJfZWZmX2NvcnIgJT4lIFxuICBmaWx0ZXIoLiwgZW52aXJvbm1lbnQxICE9IGVudmlyb25tZW50MikgJT4lIFxuICBzdW1tYXJpemVfYXQodmFycyhjb3JyKSwgZnVucyhtaW4sIG1heCwgbWVhbiksIG5hLnJtID0gVClcblxuIyBQbG90XG5TMl9NRVRfbWFya2VyX2VmZl9jb3JyICU+JSBcbiAgZ2dwbG90KGFlcyh4ID0gZW52aXJvbm1lbnQxLCB5ID0gZW52aXJvbm1lbnQyLCBmaWxsID0gY29ycikpICsgXG4gIGdlb21fdGlsZSgpICtcbiAgZmFjZXRfZ3JpZCguIH4gdHJhaXQpICtcbiAgc2NhbGVfZmlsbF9ncmFkaWVudDIoKSArXG4gIHRoZW1lX2J3KCkgK1xuICB0aGVtZShheGlzLnRleHQgPSBlbGVtZW50X3RleHQoc2l6ZSA9IDYpLFxuICAgICAgICBheGlzLnRleHQueCA9IGVsZW1lbnRfdGV4dChhbmdsZSA9IDQ1LCBoanVzdCA9IDEpLFxuICAgICAgICBheGlzLnRpdGxlID0gZWxlbWVudF9ibGFuaygpIClcblxuXG5gYGAifQ== -->

```r

# First set up two-way table of marker effects and environments
S2_MET_marker_eff_table <- S2_MET_marker_eff_pheno_fw %>% 
  select(., trait, marker, environment, mar_effect) %>% 
  distinct() %>% 
  spread(environment, mar_effect)

# Calculate the pairwise correlation of marker effects between environments
S2_MET_marker_eff_corr <- S2_MET_marker_eff_table %>%
  group_by(., trait) %>% 
  do({
    df <- .
    select(df, -trait, -marker) %>% 
      cor() %>% 
      as.data.frame() %>% 
      rownames_to_column("environment1") %>% 
      gather(environment2, corr, -environment1) })

# Find the minimum and maximum correlation
S2_MET_marker_eff_corr %>% 
  filter(., environment1 != environment2) %>% 
  summarize_at(vars(corr), funs(min, max, mean), na.rm = T)

# Plot
S2_MET_marker_eff_corr %>% 
  ggplot(aes(x = environment1, y = environment2, fill = corr)) + 
  geom_tile() +
  facet_grid(. ~ trait) +
  scale_fill_gradient2() +
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank() )

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



Plot the distribution of linear and non-linear stability estimates for markers


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG4jIEdldCB0aGUgZGlzdGluY3Qgb2JzZXJ2YXRpb25zIG9mIHRoZSBzdGFiaWxpdHkgbWVhc3VyZXNcblMyX01FVF9tYXJrZXJfZWZmX3BoZW5vX2Z3X3VuaXEgPC0gUzJfTUVUX21hcmtlcl9lZmZfcGhlbm9fZncgJT4lIFxuICBzZWxlY3QodHJhaXQsIG1hcmtlcjpwb3MsIGJfc3RkX2Vycm9yLCBkZiwgc3RhYmlsaXR5X3Rlcm0sIGVzdGltYXRlKSAlPiUgXG4gIGRpc3RpbmN0KClcblxuIyBQbG90IGhpc3RvZ3JhbXMgb2YgdGhlIHJlZ3Jlc3Npb24gY29lZmZpY2llbnRzIGFuZCBNU0VzXG5nX21hcl9md19kZW5zIDwtIFMyX01FVF9tYXJrZXJfZWZmX3BoZW5vX2Z3X3VuaXEgJT4lXG4gIGdncGxvdChhZXMoeCA9IGVzdGltYXRlKSkgKyBcbiAgZ2VvbV9kZW5zaXR5KGFlcyhmaWxsID0gXCJibHVlXCIpKSArIFxuICBmYWNldF93cmFwKHRyYWl0IH4gc3RhYmlsaXR5X3Rlcm0sIG5jb2wgPSAyLCBzY2FsZXMgPSBcImZyZWVcIikgK1xuICB4bGFiKFwiRXN0aW1hdGVcIikgK1xuICBsYWJzKHRpdGxlID0gXCJNYXJrZXIgRWZmZWN0IFN0YWJpbGl0eSBNZWFzdXJlc1wiKSArXG4gIHNjYWxlX2ZpbGxfYnJld2VyKHBhbGV0dGUgPSBcIlNldDJcIiwgZ3VpZGUgPSBGQUxTRSkgK1xuICB0aGVtZV9idygpICtcbiAgdGhlbWUoYXhpcy50aXRsZS55ID0gZWxlbWVudF9ibGFuaygpKVxuXG5zYXZlX2ZpbGUgPC0gZmlsZS5wYXRoKGZpZ19kaXIsIFwibWFya2VyX3N0YWJpbGl0eV9lc3RpbWF0ZV9kaXN0cml1YnRpb25zLmpwZ1wiKVxuZ2dzYXZlKGZpbGVuYW1lID0gc2F2ZV9maWxlLCBwbG90ID0gZ19tYXJfZndfZGVucywgaGVpZ2h0ID0gNiwgd2lkdGggPSA1KVxuXG5gYGAifQ== -->

```r

# Get the distinct observations of the stability measures
S2_MET_marker_eff_pheno_fw_uniq <- S2_MET_marker_eff_pheno_fw %>% 
  select(trait, marker:pos, b_std_error, df, stability_term, estimate) %>% 
  distinct()

# Plot histograms of the regression coefficients and MSEs
g_mar_fw_dens <- S2_MET_marker_eff_pheno_fw_uniq %>%
  ggplot(aes(x = estimate)) + 
  geom_density(aes(fill = "blue")) + 
  facet_wrap(trait ~ stability_term, ncol = 2, scales = "free") +
  xlab("Estimate") +
  labs(title = "Marker Effect Stability Measures") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_blank())

save_file <- file.path(fig_dir, "marker_stability_estimate_distriubtions.jpg")
ggsave(filename = save_file, plot = g_mar_fw_dens, height = 6, width = 5)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




Plot a genomic map of marker effect stability


Determine outlier markers using a 2.5% empirical threshold


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG4jIFdoYXQgc2hvdWxkIGJlIHRoZSBzaWduaWZpY2FuY2UgbGV2ZWxcbmFscGhhIDwtIDAuMDVcblxuIyBGb3IgZWFjaCB0cmFpdCwgY2FsY3VsYXRlIGVtcGlyaWNhbCB0aHJlc2hvbGRzIGZvciBzaWduaWZpY2FuY2VcblMyX01FVF9tYXJrZXJfZWZmX3BoZW5vX2Z3X3NpZyA8LSBTMl9NRVRfbWFya2VyX2VmZl9waGVub19md191bmlxICU+JVxuICBmaWx0ZXIoc3RhYmlsaXR5X3Rlcm0gPT0gXCJiXCIpICU+JSBcbiAgZ3JvdXBfYnkodHJhaXQpICU+JSBcbiAgIyBtdXRhdGUoZXN0aW1hdGUgPSBzY2FsZShlc3RpbWF0ZSkpICU+JVxuICBtdXRhdGUobG93ZXJfcGVyYyA9IHF1YW50aWxlKGVzdGltYXRlLCBhbHBoYSAvIDIpLCBcbiAgICAgICAgIHVwcGVyX3BlcmMgPSBxdWFudGlsZShlc3RpbWF0ZSwgMSAtIChhbHBoYSAvIDIpKSkgJT4lXG4gIHVuZ3JvdXAoKVxuXG4jIFBsb3QgdGhlIGxpbmVhciBzdGFiaWxpdHkgb3ZlciBnZW5vbWljIGRpc3RhbmNlXG5nX21hcl9lZmZfZndfYiA8LSBTMl9NRVRfbWFya2VyX2VmZl9waGVub19md19zaWcgJT4lIFxuICBnZ3Bsb3QoYWVzKHggPSBwb3MgLyAxMDAwMDAwLCB5ID0gZXN0aW1hdGUpKSArIFxuICBnZW9tX2FibGluZShzbG9wZSA9IDAsIGludGVyY2VwdCA9IDApICtcbiAgIyBUaHJlc2hvbGQgbGluZXNcbiAgZ2VvbV9obGluZShhZXMoeWludGVyY2VwdCA9IGxvd2VyX3BlcmMpLCBsdHkgPSAyKSArXG4gIGdlb21faGxpbmUoYWVzKHlpbnRlcmNlcHQgPSB1cHBlcl9wZXJjKSwgbHR5ID0gMikgK1xuICBnZW9tX3BvaW50KCkgKyBcbiAgZmFjZXRfZ3JpZCh0cmFpdCB+IGNocm9tLCBzcGFjZSA9IFwiZnJlZV94XCIsIHNjYWxlcyA9IFwiZnJlZVwiLCBzd2l0Y2ggPSBcInhcIikgK1xuICBzY2FsZV9jb2xvcl9ncmFkaWVudDIoKSArXG4gIHlsYWIoXCJTdGFiaWxpdHkgQ29lZmZpY2llbnRcIikgK1xuICB4bGFiKFwiUG9zaXRpb24gKE1icClcIikgK1xuICB0aGVtZV9idygpICtcbiAgdGhlbWUoYXhpcy50ZXh0LnggPSBlbGVtZW50X3RleHQoYW5nbGUgPSA0NSwgaGp1c3QgPSAxKSxcbiAgICAgICAgcGFuZWwuc3BhY2luZy54ID0gdW5pdCh4ID0gMCwgdW5pdHMgPSBcImluXCIpKSArXG4gIGxhYnModGl0bGUgPSBcIk1hcmtlciBMaW5lYXIgU3RhYmlsaXR5L1NlbnNpdGl2aXR5XCIpXG5cbiMgU2F2ZVxuc2F2ZV9maWxlIDwtIGZpbGUucGF0aChmaWdfZGlyLCBcIm1hcmtlcl9saW5lYXJfc3RhYmlsaXR5LmpwZ1wiKVxuZ2dzYXZlKGZpbGVuYW1lID0gc2F2ZV9maWxlLCBwbG90ID0gZ19tYXJfZWZmX2Z3X2IsIGhlaWdodCA9IDcsIHdpZHRoID0gOSlcblxuXG4jIyBQbG90IG5vbi1saW5lYXIgc3RhYmlsaXR5XG4jIEZvciBlYWNoIHRyYWl0LCBjYWxjdWxhdGUgZW1waXJpY2FsIHRocmVzaG9sZHMgZm9yIHNpZ25pZmljYW5jZVxuUzJfTUVUX21hcmtlcl9lZmZfcGhlbm9fZndfZGVsdGFfc2lnIDwtIFMyX01FVF9tYXJrZXJfZWZmX3BoZW5vX2Z3X3VuaXEgJT4lXG4gIGZpbHRlcihzdGFiaWxpdHlfdGVybSA9PSBcImRlbHRhXCIpICU+JSBcbiAgZ3JvdXBfYnkodHJhaXQpICU+JSBcbiAgIyBtdXRhdGUoZXN0aW1hdGUgPSBzY2FsZShlc3RpbWF0ZSkpICU+JVxuICBtdXRhdGUodXBwZXJfcGVyYyA9IHF1YW50aWxlKGVzdGltYXRlLCAxIC0gKGFscGhhIC8gMikpKSAlPiVcbiAgdW5ncm91cCgpXG5cbmdfbWFyX2VmZl9md19kZWx0YSA8LSBTMl9NRVRfbWFya2VyX2VmZl9waGVub19md19kZWx0YV9zaWcgJT4lIFxuICBnZ3Bsb3QoYWVzKHggPSBwb3MgLyAxMDAwMDAwLCB5ID0gZXN0aW1hdGUpKSArIFxuICBnZW9tX2hsaW5lKGFlcyh5aW50ZXJjZXB0ID0gdXBwZXJfcGVyYyksIGx0eSA9IDIpICtcbiAgZ2VvbV9wb2ludCgpICsgXG4gIGZhY2V0X2dyaWQodHJhaXQgfiBjaHJvbSwgc3BhY2UgPSBcImZyZWVfeFwiLCBzY2FsZXMgPSBcImZyZWVcIiwgc3dpdGNoID0gXCJ4XCIpICtcbiAgeWxhYihcIlN0YWJpbGl0eSBDb2VmZmljaWVudFwiKSArXG4gIHhsYWIoXCJQb3NpdGlvbiAoTWJwKVwiKSArXG4gIHRoZW1lX2J3KCkgK1xuICB0aGVtZShheGlzLnRleHQueCA9IGVsZW1lbnRfdGV4dChhbmdsZSA9IDQ1LCBoanVzdCA9IDEpLFxuICAgICAgICBwYW5lbC5zcGFjaW5nLnggPSB1bml0KHggPSAwLCB1bml0cyA9IFwiaW5cIikpICArXG4gIGxhYnModGl0bGUgPSBcIk1hcmtlciBOb24tTGluZWFyIFN0YWJpbGl0eS9TZW5zaXRpdml0eVwiKVxuXG5zYXZlX2ZpbGUgPC0gZmlsZS5wYXRoKGZpZ19kaXIsIFwibWFya2VyX25vbmxpbmVhcl9zdGFiaWxpdHkuanBnXCIpXG5nZ3NhdmUoZmlsZW5hbWUgPSBzYXZlX2ZpbGUsIHBsb3QgPSBnX21hcl9lZmZfZndfZGVsdGEsIGhlaWdodCA9IDcsIHdpZHRoID0gOSlcblxuXG5cbmBgYCJ9 -->

```r

# What should be the significance level
alpha <- 0.05

# For each trait, calculate empirical thresholds for significance
S2_MET_marker_eff_pheno_fw_sig <- S2_MET_marker_eff_pheno_fw_uniq %>%
  filter(stability_term == "b") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(lower_perc = quantile(estimate, alpha / 2), 
         upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup()

# Plot the linear stability over genomic distance
g_mar_eff_fw_b <- S2_MET_marker_eff_pheno_fw_sig %>% 
  ggplot(aes(x = pos / 1000000, y = estimate)) + 
  geom_abline(slope = 0, intercept = 0) +
  # Threshold lines
  geom_hline(aes(yintercept = lower_perc), lty = 2) +
  geom_hline(aes(yintercept = upper_perc), lty = 2) +
  geom_point() + 
  facet_grid(trait ~ chrom, space = "free_x", scales = "free", switch = "x") +
  scale_color_gradient2() +
  ylab("Stability Coefficient") +
  xlab("Position (Mbp)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing.x = unit(x = 0, units = "in")) +
  labs(title = "Marker Linear Stability/Sensitivity")

# Save
save_file <- file.path(fig_dir, "marker_linear_stability.jpg")
ggsave(filename = save_file, plot = g_mar_eff_fw_b, height = 7, width = 9)


## Plot non-linear stability
# For each trait, calculate empirical thresholds for significance
S2_MET_marker_eff_pheno_fw_delta_sig <- S2_MET_marker_eff_pheno_fw_uniq %>%
  filter(stability_term == "delta") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup()

g_mar_eff_fw_delta <- S2_MET_marker_eff_pheno_fw_delta_sig %>% 
  ggplot(aes(x = pos / 1000000, y = estimate)) + 
  geom_hline(aes(yintercept = upper_perc), lty = 2) +
  geom_point() + 
  facet_grid(trait ~ chrom, space = "free_x", scales = "free", switch = "x") +
  ylab("Stability Coefficient") +
  xlab("Position (Mbp)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing.x = unit(x = 0, units = "in"))  +
  labs(title = "Marker Non-Linear Stability/Sensitivity")

save_file <- file.path(fig_dir, "marker_nonlinear_stability.jpg")
ggsave(filename = save_file, plot = g_mar_eff_fw_delta, height = 7, width = 9)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




Use the slope estimates and standard errors to perform some hypothesis testing


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG4jIyBXaGF0IGhhcHBlbnMgaWYgd2UgdGVzdCB0aGUgbWFya2VyIHNsb3BlcyB2aWEgYSB0LXRlc3Q/XG5tYXJfZndfbGluZWFyX2FkaiA8LSBTMl9NRVRfbWFya2VyX2VmZl9waGVub19md191bmlxICAlPiVcbiAgZmlsdGVyKHN0YWJpbGl0eV90ZXJtID09IFwiYlwiKSAlPiVcbiAgZ3JvdXBfYnkodHJhaXQpICU+JVxuICBtdXRhdGUodF9zdGF0ID0gZXN0aW1hdGUgLyBiX3N0ZF9lcnJvcixcbiAgICAgICAgIHBfbm90X3plcm8gPSAyICogcHQoYWJzKHRfc3RhdCksIGRmLCBsb3dlci50YWlsID0gRkFMU0UpLFxuICAgICAgICAgcF9ndF96ZXJvID0gcHQodF9zdGF0LCBkZiwgbG93ZXIudGFpbCA9IEZBTFNFKSxcbiAgICAgICAgIHBfbHRfemVybyA9IDEgLSBwX2d0X3plcm8pICU+JVxuICAjIENvcnJlY3QgZm9yIG11bHRpcGxlIHRlc3RpbmdcbiAgbXV0YXRlX2F0KHZhcnMoc3RhcnRzX3dpdGgoXCJwX1wiKSksIGZ1bnMocXZhbCA9IHF2YWx1ZSguKSRxdmFsdWUpKSAlPiVcbiAgbXV0YXRlX2F0KHZhcnMoY29udGFpbnMoXCJxdmFsXCIpKSwgZnVucyhuZWdfbG9nID0gLWxvZzEwKC4pKSlcblxuIyBDb21iaW5lIHRoZSBzdGFiaWxpdHkgdnMgc2Vuc2l0aXZpdHkgcC12YWx1ZXNcbiMgZ3RsdCA9IGdyZWF0ZXIgdGhhbiAvIGxlc3MgdGhhblxubWFyX2Z3X2xpbmVhcl9hZGpfZ3RsdCA8LSBtYXJfZndfbGluZWFyX2FkaiAlPiVcbiAgbXV0YXRlKC4sIHBfbHRfemVyb19xdmFsX25lZ19sb2cgPSAtMSAqIHBfbHRfemVyb19xdmFsX25lZ19sb2cpICU+JVxuICBzZWxlY3QodHJhaXQ6cG9zLCBwX2d0X3plcm9fcXZhbF9uZWdfbG9nLCBwX2x0X3plcm9fcXZhbF9uZWdfbG9nKSAlPiVcbiAgZ2F0aGVyKHRlc3RfdHlwZSwgbmVnX2xvZ19wLCAtdHJhaXQ6LXBvcykgJT4lXG4gIG11dGF0ZShuZWdfbG9nMTBfZmRyMDUgPSBpZl9lbHNlKHRlc3RfdHlwZSA9PSBcInBfZ3RfemVyb19xdmFsX25lZ19sb2dcIiwgLWxvZzEwKDAuMDUpLCBsb2cxMCgwLjA1KSksXG4gICAgICAgICBuZWdfbG9nMTBfZmRyMTAgPSBpZl9lbHNlKHRlc3RfdHlwZSA9PSBcInBfZ3RfemVyb19xdmFsX25lZ19sb2dcIiwgLWxvZzEwKDAuMTApLCBsb2cxMCgwLjEwKSkpXG5cbiMgUGxvdHMgZm9yIGJvdGggZ3QgYW5kIGx0XG5nX21hcl9md19uZWdwb3NfcGxvdHMgPC0gbWFyX2Z3X2xpbmVhcl9hZGpfZ3RsdCAlPiVcbiAgbXV0YXRlKC4sIGNocm9tID0gYXMuZmFjdG9yKGNocm9tKSkgJT4lXG4gIGdncGxvdChhZXMoeCA9IHBvcyAvIDEwMDAwMDAsIHkgPSBuZWdfbG9nX3AsIGdyb3VwID0gY2hyb20sIGNvbCA9IGNocm9tKSkgK1xuICBnZW9tX3BvaW50KGFlcyhzaGFwZSA9IHRlc3RfdHlwZSkpICtcbiAgZ2VvbV9obGluZShhZXMoeWludGVyY2VwdCA9IG5lZ19sb2cxMF9mZHIwNSwgbHR5ID0gXCJGRFIgMDUlXCIpKSArXG4gIGdlb21faGxpbmUoYWVzKHlpbnRlcmNlcHQgPSBuZWdfbG9nMTBfZmRyMTAsIGx0eSA9IFwiRkRSIDEwJVwiKSkgK1xuICBmYWNldF9ncmlkKHRyYWl0IH4gY2hyb20sIHN3aXRjaCA9IFwieFwiLCBzY2FsZXMgPSBcImZyZWVcIiwgc3BhY2UgPSBcImZyZWVfeFwiKSArXG4gIHNjYWxlX2NvbG9yX21hbnVhbChndWlkZSA9IEZBTFNFLFxuICAgICAgICAgICAgICAgICAgICAgdmFsdWVzID0gc2V0X25hbWVzKHJlcChjKFwiYmxhY2tcIiwgXCJncmV5NzVcIiksIGxlbmd0aC5vdXQgPSA3KSwgc2VxKDEsIDcpKSkgK1xuICBzY2FsZV9zaGFwZV9tYW51YWwodmFsdWVzID0gYygxNiwxKSwgbGFiZWxzID0gYyhcIlNlbnNpdGl2aXR5XCIsIFwiU3RhYmlsaXR5XCIpLCBuYW1lID0gXCJUZXN0XCIpICtcbiAgc2NhbGVfbGluZXR5cGVfZGlzY3JldGUoZ3VpZGUgPSBGQUxTRSkgK1xuICB5bGFiKFwiLWxvZyhxKVwiKSArXG4gIHhsYWIoXCJQb3NpdGlvbiAoTWJwKVwiKSArXG4gIGxhYnModGl0bGUgPSBcIkdlbm9tZXdpZGUgVGVzdCBmb3IgTWFya2VyIExpbmVhciBTdGFiaWxpdHlcIikgK1xuICB0aGVtZV9idygpICtcbiAgdGhlbWUocGFuZWwuc3BhY2luZy54ID0gdW5pdCh4ID0gMCwgdW5pdHMgPSBcImluXCIpLFxuICAgICAgICBwYW5lbC5ib3JkZXIgPSBlbGVtZW50X2JsYW5rKCksXG4gICAgICAgIGF4aXMudGV4dC54ID0gZWxlbWVudF90ZXh0KGFuZ2xlID0gNDUsIGhqdXN0ID0gMSkpXG5cblxuc2F2ZV9maWxlIDwtIGZpbGUucGF0aChmaWdfZGlyLCBcIm1hcmtlcl9saW5lYXJfc3RhYmlsaXR5X2d0bGVfaHlwdGVzdC5qcGdcIilcbmdnc2F2ZShmaWxlbmFtZSA9IHNhdmVfZmlsZSwgcGxvdCA9IGdfbWFyX2Z3X25lZ3Bvc19wbG90cywgaGVpZ2h0ID0gNywgd2lkdGggPSA5KVxuXG5cbiMjIFBsb3QganVzdCB0aGUgdGVzdCBvZiB0aGUgc2xvcGUgYmVpbmcgZGlmZmVyZW50IGZyb20gemVyb1xubWFyX2Z3X2xpbmVhcl9hZGpfbnogPC0gbWFyX2Z3X2xpbmVhcl9hZGogJT4lXG4gIG11dGF0ZSguLCBuZWdfbG9nMTBfZmRyMDUgPSAtbG9nMTAoMC4wNSksIG5lZ19sb2cxMF9mZHIxMCA9IC1sb2cxMCgwLjEwKSlcblxuZ19tYXJfZndfbnpfcGxvdHMgPC0gbWFyX2Z3X2xpbmVhcl9hZGpfbnogJT4lXG4gIG11dGF0ZSguLCBjaHJvbSA9IGFzLmZhY3RvcihjaHJvbSkpICU+JVxuICBnZ3Bsb3QoYWVzKHggPSBwb3MgLyAxMDAwMDAwLCB5ID0gcF9ub3RfemVyb19xdmFsX25lZ19sb2cpKSArXG4gIGdlb21fcG9pbnQoKSArXG4gIGdlb21faGxpbmUoYWVzKHlpbnRlcmNlcHQgPSBuZWdfbG9nMTBfZmRyMDUsIGx0eSA9IFwiRkRSIDA1JVwiKSkgK1xuICBnZW9tX2hsaW5lKGFlcyh5aW50ZXJjZXB0ID0gbmVnX2xvZzEwX2ZkcjEwLCBsdHkgPSBcIkZEUiAxMCVcIikpICtcbiAgZmFjZXRfZ3JpZCh0cmFpdCB+IGNocm9tLCBzd2l0Y2ggPSBcInhcIiwgc2NhbGVzID0gXCJmcmVlXCIsIHNwYWNlID0gXCJmcmVlX3hcIikgK1xuICBzY2FsZV9jb2xvcl9tYW51YWwoZ3VpZGUgPSBGQUxTRSxcbiAgICAgICAgICAgICAgICAgICAgIHZhbHVlcyA9IHNldF9uYW1lcyhyZXAoYyhcImJsYWNrXCIsIFwiZ3JleTc1XCIpLCBsZW5ndGgub3V0ID0gNyksIHNlcSgxLCA3KSkpICtcbiAgc2NhbGVfbGluZXR5cGVfZGlzY3JldGUoZ3VpZGUgPSBGQUxTRSkgK1xuICB5bGFiKFwiLWxvZyhwKVwiKSArXG4gIHhsYWIoXCJQb3NpdGlvbiAoTWJwKVwiKSArXG4gICMgbGFicyh0aXRsZSA9IFwiR2Vub21ld2lkZSBBc3NvY2F0aW9uIEFuYWx5c2lzIGZvciBQaGVub3R5cGljIFN0YWJpbGl0eVwiLFxuICAjICAgICAgY2FwdGlvbiA9IFwibiA9IDE4MyBsaW5lcyB1c2VkIHRvIGNhbGN1bGF0ZSBzdGFiaWxpdHkgY29lZmZpY2llbnRzIGFuZCBuID0gMTc1IGxpbmVzIHVzZWRcXG5mb3IgYXNzb2NpYXRpb24gYW5hbHlzaXMuIFRoZSBib2xkIGxpbmUgaW5kaWNhdGVzIHRoZSBnZW5vbWV3aWRlIEZEUiB0aHJlc2hvbGQgYXQgNSU7XFxudGhlIGRhc2hlZCBsaW5lIGF0IDEwJS5cIikgK1xuICB0aGVtZV9idygpICtcbiAgdGhlbWUoXG4gICAgcGFuZWwuc3BhY2luZy54ID0gdW5pdCh4ID0gMCwgdW5pdHMgPSBcImluXCIpLFxuICAgIHBhbmVsLmJvcmRlciA9IGVsZW1lbnRfYmxhbmsoKSxcbiAgICBheGlzLnRleHQueCA9IGVsZW1lbnRfdGV4dChhbmdsZSA9IDQ1LCBoanVzdCA9IDEpXG4gIClcblxuc2F2ZV9maWxlIDwtIGZpbGUucGF0aChmaWdfZGlyLCBcIm1hcmtlcl9saW5lYXJfc3RhYmlsaXR5X256X2h5cHRlc3QuanBnXCIpXG5nZ3NhdmUoZmlsZW5hbWUgPSBzYXZlX2ZpbGUsIHBsb3QgPSBnX21hcl9md19uel9wbG90cywgaGVpZ2h0ID0gNywgd2lkdGggPSA5KVxuXG5cbmBgYCJ9 -->

```r

## What happens if we test the marker slopes via a t-test?
mar_fw_linear_adj <- S2_MET_marker_eff_pheno_fw_uniq  %>%
  filter(stability_term == "b") %>%
  group_by(trait) %>%
  mutate(t_stat = estimate / b_std_error,
         p_not_zero = 2 * pt(abs(t_stat), df, lower.tail = FALSE),
         p_gt_zero = pt(t_stat, df, lower.tail = FALSE),
         p_lt_zero = 1 - p_gt_zero) %>%
  # Correct for multiple testing
  mutate_at(vars(starts_with("p_")), funs(qval = qvalue(.)$qvalue)) %>%
  mutate_at(vars(contains("qval")), funs(neg_log = -log10(.)))

# Combine the stability vs sensitivity p-values
# gtlt = greater than / less than
mar_fw_linear_adj_gtlt <- mar_fw_linear_adj %>%
  mutate(., p_lt_zero_qval_neg_log = -1 * p_lt_zero_qval_neg_log) %>%
  select(trait:pos, p_gt_zero_qval_neg_log, p_lt_zero_qval_neg_log) %>%
  gather(test_type, neg_log_p, -trait:-pos) %>%
  mutate(neg_log10_fdr05 = if_else(test_type == "p_gt_zero_qval_neg_log", -log10(0.05), log10(0.05)),
         neg_log10_fdr10 = if_else(test_type == "p_gt_zero_qval_neg_log", -log10(0.10), log10(0.10)))

# Plots for both gt and lt
g_mar_fw_negpos_plots <- mar_fw_linear_adj_gtlt %>%
  mutate(., chrom = as.factor(chrom)) %>%
  ggplot(aes(x = pos / 1000000, y = neg_log_p, group = chrom, col = chrom)) +
  geom_point(aes(shape = test_type)) +
  geom_hline(aes(yintercept = neg_log10_fdr05, lty = "FDR 05%")) +
  geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")) +
  facet_grid(trait ~ chrom, switch = "x", scales = "free", space = "free_x") +
  scale_color_manual(guide = FALSE,
                     values = set_names(rep(c("black", "grey75"), length.out = 7), seq(1, 7))) +
  scale_shape_manual(values = c(16,1), labels = c("Sensitivity", "Stability"), name = "Test") +
  scale_linetype_discrete(guide = FALSE) +
  ylab("-log(q)") +
  xlab("Position (Mbp)") +
  labs(title = "Genomewide Test for Marker Linear Stability") +
  theme_bw() +
  theme(panel.spacing.x = unit(x = 0, units = "in"),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


save_file <- file.path(fig_dir, "marker_linear_stability_gtle_hyptest.jpg")
ggsave(filename = save_file, plot = g_mar_fw_negpos_plots, height = 7, width = 9)


## Plot just the test of the slope being different from zero
mar_fw_linear_adj_nz <- mar_fw_linear_adj %>%
  mutate(., neg_log10_fdr05 = -log10(0.05), neg_log10_fdr10 = -log10(0.10))

g_mar_fw_nz_plots <- mar_fw_linear_adj_nz %>%
  mutate(., chrom = as.factor(chrom)) %>%
  ggplot(aes(x = pos / 1000000, y = p_not_zero_qval_neg_log)) +
  geom_point() +
  geom_hline(aes(yintercept = neg_log10_fdr05, lty = "FDR 05%")) +
  geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")) +
  facet_grid(trait ~ chrom, switch = "x", scales = "free", space = "free_x") +
  scale_color_manual(guide = FALSE,
                     values = set_names(rep(c("black", "grey75"), length.out = 7), seq(1, 7))) +
  scale_linetype_discrete(guide = FALSE) +
  ylab("-log(p)") +
  xlab("Position (Mbp)") +
  # labs(title = "Genomewide Assocation Analysis for Phenotypic Stability",
  #      caption = "n = 183 lines used to calculate stability coefficients and n = 175 lines used\nfor association analysis. The bold line indicates the genomewide FDR threshold at 5%;\nthe dashed line at 10%.") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(x = 0, units = "in"),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_file <- file.path(fig_dir, "marker_linear_stability_nz_hyptest.jpg")
ggsave(filename = save_file, plot = g_mar_fw_nz_plots, height = 7, width = 9)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Filter outliers and determine windows

We will continue only with the TP data and for the slope of the regression curve. We define stable marker outliers as those below the lower threshold and sensitive markers as those above the upper threshold




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBEZXRlcm1pbmUgdGhlIG51bWJlciBvZiBzdGFibGUgb3Igc2Vuc2l0aXZlIG1hcmtlcnMgZm9yIGVhY2ggdHJhaXRcbm1hcl9lZmZfZndfYl90cCA8LSBTMl9NRVRfbWFya2VyX2VmZl9waGVub19md19zaWckcG9wX3RwICU+JSBcbiAgc2VsZWN0KHRyYWl0LCBtYXJrZXI6Y01fcG9zLCBlc3RpbWF0ZTp1cHBlcl9wZXJjKSAlPiUgXG4gIG11dGF0ZShvdXRsaWVyID0gY2FzZV93aGVuKGVzdGltYXRlID49IHVwcGVyX3BlcmMgfiBcInNlbnNpdGl2ZVwiLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgZXN0aW1hdGUgPD0gbG93ZXJfcGVyYyB+IFwic3RhYmxlXCIsIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICBUUlVFIH4gXCJub3Rfb3V0bGllclwiKSkgJT4lXG4gIGRpc3RpbmN0KClcblxuIyBGaWx0ZXIgdGhlIG91dGxpZXJzXG5tYXJfc3RhYmlsaXR5X291dGxpZXJzIDwtIG1hcl9lZmZfZndfYl90cCAlPiVcbiAgZmlsdGVyKG91dGxpZXIgIT0gXCJub3Rfb3V0bGllclwiKVxuXG4jIFN1bW1hcml6ZVxubWFyX3N0YWJpbGl0eV9vdXRsaWVycyAlPiUgXG4gIGdyb3VwX2J5KHRyYWl0LCBjaHJvbSwgb3V0bGllcikgJT4lIFxuICBzdW1tYXJpemUobl9vdXRsaWVycyA9IG4oKSlcblxuXG5cblxuIyMgRm9yIGVhY2ggdHJhaXQsIGNocm9tb3NvbWUsIGFuZCBvdXRsaWVyIHR5cGUsIGNyZWF0ZSBkaWZmZXJlbnQgd2luZG93cyBmcm9tIGVhY2hcbiMjIG91dGxpZXIgU05QIHRvIGNsdXN0ZXIgb3RoZXIgb3V0bGllciBTTlBTXG5cbiMgVmVjdG9yIG9mIHdpbmRvdyBzaXplcyAoaW4gYnApXG53aW5kX3NpemVfbGlzdCA8LSBjKDUwMDAwLCAxMDAwMDAsIDUwMDAwMCwgMTAwMDAwMCwgNTAwMDAwMDApXG5cbm1hcl9zdGFiaWxpdHlfb3V0bGllcnNfd2luZG93IDwtIG1hcl9zdGFiaWxpdHlfb3V0bGllcnMgJT4lIFxuICBzcGxpdChsaXN0KC4kdHJhaXQsIC4kY2hyb20sIC4kb3V0bGllcikpICU+JSBcbiAgLlttYXBfZGJsKC4sIG5yb3cpID4gMF0gJT4lXG4gIGxpc3QoLiwgbWFwKC4sIGZ1bmN0aW9uKGRmKSB7XG4gICAgIyBJdGVyYXRlIG92ZXIgd2luZG93IHNpemVzXG4gICAgZ3JvdXBzIDwtIHdpbmRfc2l6ZV9saXN0ICU+JSBcbiAgICAgIG1hcChmdW5jdGlvbih3aW5kX3NpemUpIHtcbiAgICAgICAgIyBDcmVhdGUgYW4gSVJhbmdlcyBvYmplY3RcbiAgICAgICAgZGZfaXJhbmdlcyA8LSBJUmFuZ2VzKHN0YXJ0ID0gZGYkcG9zIC0gd2luZF9zaXplLCBlbmQgPSBkZiRwb3MgKyB3aW5kX3NpemUpXG4gICAgICAgIFxuICAgICAgICAjIFJlZHVjZSBhbmQgZmluZCB1bmlxdWUgZ3JvdXBzXG4gICAgICAgIGRmX2lyYW5nZXMgJT4lIFxuICAgICAgICAgIHJlZHVjZSgpICU+JSBcbiAgICAgICAgICBmaW5kT3ZlcmxhcHMocXVlcnkgPSAuLCBzdWJqZWN0ID0gZGZfaXJhbmdlcykgJT4lIFxuICAgICAgICAgIHF1ZXJ5SGl0cygpICU+JVxuICAgICAgICAgIHN0cl9jKFwiY2hyXCIsIGRmJGNocm9tLCBcIl9cIiwgLikgfSlcbiAgICBcbiAgICAjIFJldHVybiBhIGRhdGEuZnJhbWVcbiAgICBncm91cHMgJT4lIFxuICAgICAgc2V0X25hbWVzKHN0cl9jKFwid2luZG93X1wiLCB3aW5kX3NpemVfbGlzdCkpICU+JSBcbiAgICAgIGFzX2RhdGFfZnJhbWUoKSB9KSkgJT4lXG4gIHBtYXBfZGYofmJpbmRfY29scygueCwgLnkpKVxuXG5cbmBgYCJ9 -->

```r
# Determine the number of stable or sensitive markers for each trait
mar_eff_fw_b_tp <- S2_MET_marker_eff_pheno_fw_sig$pop_tp %>% 
  select(trait, marker:cM_pos, estimate:upper_perc) %>% 
  mutate(outlier = case_when(estimate >= upper_perc ~ "sensitive", 
                             estimate <= lower_perc ~ "stable", 
                             TRUE ~ "not_outlier")) %>%
  distinct()

# Filter the outliers
mar_stability_outliers <- mar_eff_fw_b_tp %>%
  filter(outlier != "not_outlier")

# Summarize
mar_stability_outliers %>% 
  group_by(trait, chrom, outlier) %>% 
  summarize(n_outliers = n())




## For each trait, chromosome, and outlier type, create different windows from each
## outlier SNP to cluster other outlier SNPS

# Vector of window sizes (in bp)
wind_size_list <- c(50000, 100000, 500000, 1000000, 50000000)

mar_stability_outliers_window <- mar_stability_outliers %>% 
  split(list(.$trait, .$chrom, .$outlier)) %>% 
  .[map_dbl(., nrow) > 0] %>%
  list(., map(., function(df) {
    # Iterate over window sizes
    groups <- wind_size_list %>% 
      map(function(wind_size) {
        # Create an IRanges object
        df_iranges <- IRanges(start = df$pos - wind_size, end = df$pos + wind_size)
        
        # Reduce and find unique groups
        df_iranges %>% 
          reduce() %>% 
          findOverlaps(query = ., subject = df_iranges) %>% 
          queryHits() %>%
          str_c("chr", df$chrom, "_", .) })
    
    # Return a data.frame
    groups %>% 
      set_names(str_c("window_", wind_size_list)) %>% 
      as_data_frame() })) %>%
  pmap_df(~bind_cols(.x, .y))

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->






<!-- rnb-text-end -->

