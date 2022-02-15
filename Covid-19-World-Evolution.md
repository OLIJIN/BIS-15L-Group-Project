---
title: "Group Project"
author: "Qianyu Jin & Zayna"
date: "2/15/2022"
output: 
  html_document: 
    keep_md: yes
---


```r
library(tidyverse)
library(RColorBrewer)
library(paletteer)
library(janitor)
library(here)
```
Loading Covid-19 World Evolution Data set.

```r
covid_world_evolution <- read_csv(here("data", "covid-variants.csv"))
covid_world_evolution
```

```
## # A tibble: 100,416 × 6
##    location date       variant   num_sequences perc_sequences num_sequences_tot…
##    <chr>    <date>     <chr>             <dbl>          <dbl>              <dbl>
##  1 Angola   2020-07-06 Alpha                 0              0                  3
##  2 Angola   2020-07-06 B.1.1.277             0              0                  3
##  3 Angola   2020-07-06 B.1.1.302             0              0                  3
##  4 Angola   2020-07-06 B.1.1.519             0              0                  3
##  5 Angola   2020-07-06 B.1.160               0              0                  3
##  6 Angola   2020-07-06 B.1.177               0              0                  3
##  7 Angola   2020-07-06 B.1.221               0              0                  3
##  8 Angola   2020-07-06 B.1.258               0              0                  3
##  9 Angola   2020-07-06 B.1.367               0              0                  3
## 10 Angola   2020-07-06 B.1.620               0              0                  3
## # … with 100,406 more rows
```
