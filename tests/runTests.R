require("dagLogo") || stop("unable to load Package:dagLogo")
require("Biostrings") || stop("unable to load Package:Biostrings")
require("biomaRt") || stop("unable to load Package:biomaRt")
data(seq.example)
data(proteome.example)
library(testthat)

test_check("dagLogo")