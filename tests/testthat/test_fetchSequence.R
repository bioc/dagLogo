test_that("test_fetchSequence", {
    fasta_file <- system.file("extdata", "sample_seq.fasta", 
        package = "dagLogo")
    seq_set <- readAAStringSet(fasta_file)
    names(seq_set) <- sapply(strsplit(names(seq_set), "|", fixed = TRUE), 
        function(chunks) chunks[[2]])
    dag.proteome <- prepareProteome(fasta = seq_set)
    dag.peptides <- fetchSequence(c("P07900", "Q15185"), type = "uniprotswissprot", 
        proteome = dag.proteome, anchorPos = c("Q10", "Q84"), 
        upstreamOffset = 5, downstreamOffset = 5)
    expect_equal(paste(dag.peptides@peptides[1, ], collapse = ""), 
        "TQTQDQPMEEE")
    if(Sys.getenv("USER")=="jianhongou"){
      mart <- useMart("ensembl")
      mart <- useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
      protein <- fetchSequence(IDs = unique(as.character(1:2)),
                             type = "entrezgene",
                             mart = mart,
                             anchorPos = c("L11", "L11"),
                             upstreamOffset = 5,
                             downstreamOffset = 5)
    }
})

