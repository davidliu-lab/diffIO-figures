# diffIO-quarto

Author: Amy Huang

Quarto book to capture all the analyses for the differential IO project. 

## Set up

[Install Quarto](https://quarto.org/docs/get-started/). 

Run `quarto render`. 

## Installation notes

I was having trouble installing `xml2` which is a pretty common dependency, because, for whatever reason, where brew installed libxml2 wasn't being recognized by `install.packages`. This was the solution: 

```
install.packages("xml2", configure.vars = c("INCLUDE_DIR=/path/to/local/libxml2"))
```

For instance, my installation as at `/usr/local/Cellar/libxml2`. 

I also had an issue installing `igraph` which is dependency for Seurat. My current solution has been to install it as follows (based on advice from here: https://github.com/igraph/rigraph/issues/275): 

```
install.packages("igraph", INSTALL_opts = "--no-clean-on-error")
```

and I don't compile from source. 

Also, Matrix.utils was archived off of the CRAN repository as off 2022. 
But you can still install it through (see (this GitHub issue from monocle)[https://github.com/cole-trapnell-lab/monocle3/issues/627]): 

```
remotes::install_github("cvarrichio/Matrix.utils")
```

A few packages are also from my GitHub: 

```
remotes::install_github("amyh25/RFunctions")
remotes::install_github("amyh25/Rsc")
```
