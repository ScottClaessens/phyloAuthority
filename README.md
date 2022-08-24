# The Coevolution of Religious and Political Authority in Austronesian Societies

Using novel phylogenetic methods to study the co-evolution of religious and political authority across Austronesia.

## Getting Started

### Installing

To run this code, you will need to [install R](https://www.r-project.org/) and the following R packages:

```
install.packages(c("ape", "bayestestR", "brms", "colorspace", 
                   "cowplot", "future", "future.callr", "geosphere", 
                   "ggrepel", "HDInterval", "phangorn", "phaseR", 
                   "phytools", "readxl", "rstan", "targets", 
                   "tarchetypes", "tidybayes", "tidyverse"))
```

You will also require the `rethinking` package, which you can install with the following code (see [here](https://github.com/rmcelreath/rethinking) for more details):

```
# from https://github.com/rmcelreath/rethinking
install.packages(c("coda", "mvtnorm", "devtools", "loo", "dagitty"))
devtools::install_github("rmcelreath/rethinking")
```

Finally, you will require the `ggtree` package, which you can install with the following code (see [here](https://bioconductor.org/packages/release/bioc/html/ggtree.html) for more details):

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

### Executing code

1. Set the working directory to this code repository
2. Load the `targets` package with `library(targets)`
3. To run all analyses, run `tar_make()` (or `tar_make_future(workers = N)` to run the pipeline in parallel on your local machine, with N workers)
4. To load individual targets into your environment, run `tar_load(targetName)`

## Help

Any issues, please email scott.claessens@gmail.com.

## Authors

Scott Claessens, scott.claessens@gmail.com
