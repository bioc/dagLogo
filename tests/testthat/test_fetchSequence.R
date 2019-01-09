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
})

