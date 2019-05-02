## Create two cached environments to store those package-wide accessible constants.
#### cachedEnv stores grouping schemes
cachedEnv <- new.env(parent = emptyenv())

#### .globalEnv stores motifStack plotting configuration
.globalEnv <- new.env(parent = emptyenv())

## color and character symbol encodings and grouping for amino acids
## no grouping, colored using 20 distinct colors
no <- auto <- list(
  color =  c('A'='#CCFF00', 'C'='#FFFF00', 'D'='#FF0000', 'E'='#FF0066', 
             'F'='#00FF66', 'G'='#FF9900', 'H'='#0066FF', 'I'='#66FF00', 
             'K'='#6600FF', 'L'='#33FF00', 'M'='#00FF00', 'N'='#CC00FF', 
             'P'='#FFCC00', 'Q'='#FF00CC', 'R'='#0000FF', 'S'='#FF3300', 
             'T'='#FF6600', 'V'='#99FF00', 'W'='#00CCFF', 'Y'='#00FFCC'),
  symbol = c("Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D",
             "Cys" = "C", "Glu" = "E", "Gln" = "Q", "Gly" = "G",
             "His" = "H", "Ile" = "I", "Leu" = "L", "Lys" = "K",
             "Met" = "M", "Phe" = "F", "Pro" = "P", "Ser" = "S",
             "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V"),
  group = NULL)

#### no grouping, colored based on consensus similarities of 34 individual 
#### grouping schemes
# D consensus similarity index
# R PMID:23743923 
# A Stephenson JD and Freeland SJ
# T Unearthing the root of amino acid similarity
# J J. Mol. Evol. 77:159-169 (2013) 
consensus_similarity_SF <- list(
  color = c(
    D = "#0082c8", E = "#0082c8", N = "#0082c8", K = "#fabebe", R = "#fabebe",
    Q = "#fabebe", S = "#808000", T = "#808000", H = "#f032e6", G = "#46f0f0",
    P = "#46f0f0", A = "#e6194b", C = '#3cb44b', W = "#808080", Y = "#808080",
    F = "#808080", M = "#e6beff", L = "#e6beff", I = "#e6beff", V = "#e6beff"),
  symbol = c("Asp" = "D", "Glu" = "E", "Asn" = "N", "Lys" = "K", "Arg" = "R",
             "Gln" = "Q", "Ser" = "S", "Thr" = "T", "His" = "H", "Gly" = "G",
             "Pro" = "P", "Ala" = "A", "Cys" = "C", "Tyr" = "Y", "Trp" = "W",
             "Phe" = "F", "Met" = "M", "Leu" = "L", "Ile" = "I", "Val" = "V"),
  group = NULL)


# H KYTJ820101
# D Hydropathy index (Kyte-Doolittle, 1982)
# R PMID:7108955
# A Kyte, J. and Doolittle, R.F.
# T A simple method for displaying the hydropathic character of a protein
# J J. Mol. Biol. 157, 105-132 (1982)
#### A type of consensus scale based on a combinations of experimental 
#### observations from other studies.

hydrophobicity_KD <- list(
  color = c("R" = "#0000FF", "K" = "#0C00F2", "N" = "#1900E5",
            "D" = "#1900E5", "Q" = "#1900E5", "E" = "#1900E5",
            "H" = "#2200DC", "P" = "#5200AC", "Y" = "#5A00A4",
            "W" = "#670097", "S" = "#670097", "T" = "#6C0092",
            "G" = "#74008A", "A" = "#B50049", "M" = "#B50049",
            "C" = "#C60038", "F" = "#CF002F", "L" = "#ED0011",
            "V" = "#FA0004", "I" = "#FF0000"),
  symbol = c("Arg" = "R", "Lys" = "K", "Asn" = "N", "Asp" = "D",
             "Gln" = "Q", "Glu" = "E", "His" = "H", "Pro" = "P",
             "Tyr" = "Y", "Trp" = "W", "Ser" = "S", "Thr" = "T",
             "Gly" = "G", "Ala" = "A", "Met" = "M", "Cys" = "C",
             "Phe" = "F", "Leu" = "L", "Val" = "V", "Ile" = "I"),
  group = NULL)

# H HOPT810101
# D Hydrophilicity value (Hopp-Woods, 1981)
# R PMID:6167991
# A Hopp, T.P. and Woods, K.R.
# T Prediction of protein antigenic determinants from amino acid sequecces
# J Proc. Natl. Acad. Sci. USA 78, 3824-3828 (1981)
#### A hydrophilicity scale based on the water solubility of individual amino 
#### acids 
hydrophobicity_HW <- list(
  color = c("W" = "#0000FF", "F" = "#2200DC", "Y" = "#2B00D3",
            "I" = "#3C00C2", "L" = "#3C00C2", "V" = "#4900B5",
            "M" = "#5200AC", "C" = "#5F009F", "A" = "#74008A",
            "H" = "#74008A", "T" = "#790085", "G" = "#850079",
            "P" = "#850079", "N" = "#8E0070", "Q" = "#8E0070",
            "S" = "#92006C", "R" = "#FF0000", "D" = "#FF0000",
            "E" = "#FF0000", "K" = "#FF0000"),
  symbol = c("Trp" = "W", "Phe" = "F", "Tyr" = "Y", "Ile" = "I",
             "Leu" = "L", "Val" = "V", "Met" = "M", "Cys" = "C",
             "Ala" = "A", "His" = "H", "Thr" = "T", "Gly" = "G",
             "Pro" = "P", "Asn" = "N", "Gln" = "Q", "Ser" = "S",
             "Arg" = "R", "Asp" = "D", "Glu" = "E", "Lys" = "K"),
  group = NULL)

# H GRAR740102
# D Polarity (Grantham, 1974)
# R PMID:4843792
# A Grantham, R.
# T Amino acid difference formula to help explain protein evolution
# J Science 185, 862-864 (1974)
polarity_Grantham <- list(
  color = c("L" = "#0000FF", "I" = "#0800F6", "F" = "#0800F6",
            "W" = "#0C00F2", "C" = "#1100ED", "M" = "#1500E9",
            "V" = "#1E00E0", "Y" = "#2600D8", "P" = "#5F009F",
            "A" = "#63009B", "T" = "#74008A", "G" = "#81007D",
            "S" = "#850079", "H" = "#AC0052", "R" = "#B1004D",
            "Q" = "#B1004D", "K" = "#CB0033", "N" = "#D3002B",
            "E" = "#E90015", "D" = "#FF0000"),
  symbol = c("Leu" = "L", "Ile" = "I", "Phe" = "F", "Trp" = "W",
             "Cys" = "C", "Met" = "M", "Val" = "V", "Tyr" = "Y",
             "Pro" = "P", "Ala" = "A", "Thr" = "T", "Gly" = "G",
             "Ser" = "S", "His" = "H", "Arg" = "R", "Gln" = "Q",
             "Lys" = "K", "Asn" = "N", "Glu" = "E", "Asp" = "D"),
  group = NULL)

# H ZIMJ680104
# D Isoelectric point (Zimmerman et al., 1968)
# R PMID:5700434
# A Zimmerman, J.M., Eliezer, N. and Simha, R.
# T The characterization of amino acid sequences in proteins by statistical 
# methods
# J J. Theor. Biol. 21, 170-201 (1968)
isoelectric_point_Zimmerman <- list(
  color = c("D" = "#0000FF", "E" = "#0C00F2", "C" = "#4900B5",
            "N" = "#5200AC", "F" = "#5600A8", "Q" = "#5A00A4",
            "T" = "#5A00A4", "Y" = "#5A00A4", "S" = "#5A00A4",
            "M" = "#5F009F", "W" = "#63009B", "V" = "#63009B",
            "G" = "#670097", "L" = "#670097", "A" = "#670097",
            "I" = "#670097", "P" = "#70008E", "H" = "#9B0063",
            "K" = "#E0001E", "R" = "#FF0000"),
  symbol = c("Asp" = "D", "Glu" = "E", "Cys" = "C", "Asn" = "N",
             "Phe" = "F", "Gln" = "Q", "Thr" = "T", "Tyr" = "Y",
             "Ser" = "S", "Met" = "M", "Trp" = "W", "Val" = "V",
             "Gly" = "G", "Leu" = "L", "Ala" = "A", "Ile" = "I",
             "Pro" = "P", "His" = "H", "Lys" = "K", "Arg" = "R"),
  group = NULL)

# H ZIMJ680102
# D Bulkiness (Zimmerman et al., 1968)
# R PMID:5700434
# A Zimmerman, J.M., Eliezer, N. and Simha, R.
# T The characterization of amino acid sequences in proteins by statistical 
# methods
# J J. Theor. Biol. 21, 170-201 (1968)
bulkiness_Zimmerman <- list(
  color = c("G" = "#0000FF", "S" = "#5200AC", "A" = "#70008E",
            "D" = "#74008A", "N" = "#81007D", "C" = "#8E0070",
            "E" = "#8E0070", "H" = "#8E0070", "R" = "#970067",
            "Q" = "#9B0063", "K" = "#AC0052", "T" = "#AC0052",
            "M" = "#B50049", "P" = "#C60038", "Y" = "#CF002F",
            "F" = "#E50019", "I" = "#FF0000", "L" = "#FF0000",
            "V" = "#FF0000", "W" = "#FF0000"),
  symbol = c("Gly" = "G", "Ser" = "S", "Ala" = "A", "Asp" = "D",
             "Asn" = "N", "Cys" = "C", "Glu" = "E", "His" = "H",
             "Arg" = "R", "Gln" = "Q", "Lys" = "K", "Thr" = "T",
             "Met" = "M", "Pro" = "P", "Tyr" = "Y", "Phe" = "F",
             "Ile" = "I", "Leu" = "L", "Val" = "V", "Trp" = "W"),
  group = NULL)

# H BIGC670101
# D Residue volume (Bigelow, 1967)
# R PMID:6048539
# A Bigelow, C.C.
# T On the average hydrophobicity of proteins and the relation between it and 
# protein structure
# J J. Theor. Biol. 16, 187-211 (1967) (Asn Gln 5.0)

volume_Bigelow <- list(
  color = c("G" = "#0000FF", "A" = "#2600D8", "S" = "#2F00CF",
            "C" = "#5200AC", "D" = "#5200AC", "T" = "#5A00A4",
            "P" = "#5F009F", "N" = "#63009B", "E" = "#7D0081",
            "V" = "#7D0081", "Q" = "#8A0074", "H" = "#8E0070",
            "M" = "#9F005F", "I" = "#A80056", "L" = "#A80056",
            "K" = "#B1004D", "R" = "#BE0040", "F" = "#C60038",
            "Y" = "#CF002F", "W" = "#FF0000"),
  symbol = c("Gly" = "G", "Ala" = "A", "Ser" = "S", "Cys" = "C",
             "Asp" = "D", "Thr" = "T", "Pro" = "P", "Asn" = "N",
             "Glu" = "E", "Val" = "V", "Gln" = "Q", "His" = "H",
             "Met" = "M", "Ile" = "I", "Leu" = "L", "Lys" = "K",
             "Arg" = "R", "Phe" = "F", "Tyr" = "Y", "Trp" = "W"),
  group = NULL)

## chemistry property similarity of individual AAs (Mahler)
chemistry_property_Mahler <- list(
  color = c(
    D = "#0082c8", E = "#0082c8", K = "#fabebe", R = "#fabebe", H = "#fabebe",
    Q = "#800000", N = "#800000", S = "#808000", T = "#808000", P = "#00FF00",
    C = '#3cb44b', M = '#3cb44b', W = "#808080", Y = "#808080", F = "#808080",
    G = "#46f0f0", A = "#46f0f0", L = "#46f0f0", I = "#46f0f0", V = "#46f0f0"),
  symbol = c(
    "Asp" = "D", "Glu" = "E", "Lys" = "K", "Arg" = "R", "His" = "H",
    "Gln" = "Q", "Asn" = "N", "Ser" = "S", "Thr" = "T", "Pro" = "P", 
    "Cys" = "C", "Met" = "M", "Trp" = "W", "Tyr" = "Y", "Phe" = "F",
    "Gly" = "G", "Ala" = "A", "Leu" = "L", "Ile" = "I", "Val" = "V"),
  group = NULL)


## substitution index (Dayhoff)
sequence_alignment_Dayhoff <- list(
  color = c(
    D = "#0082c8", E = "#0082c8", N = "#0082c8", Q = "#0082c8", K = "#fabebe",
    R = "#fabebe", H = "#fabebe", S = "#808000", T = "#808000", G = "#808000",
    P = "#808000", A = "#808000", C = '#3cb44b', W = "#808080", Y = "#808080",
    F = "#808080", M = "#e6beff", L = "#e6beff", I = "#e6beff", V = "#e6beff"),
  symbol = c(
    "Asp" = "D", "Glu" = "E", "Asn" = "N", "Gln" = "Q", "Lys" = "K", 
    "Arg" = "R", "His" = "H", "Ser" = "S", "Thr" = "T", "Gly" = "G",
    "Pro" = "P", "Ala" = "A", "Cys" = "C", "Trp" = "W", "Tyr" = "Y", 
    "Phe" = "F", "Met" = "M", "Leu" = "L", "Ile" = "I", "Val" = "V"),
  group = NULL)


## structural alignment (Mirny)
structure_alignments_Mirny <- list(
  color = c(
    D = "#0082c8", E = "#0082c8", K = "#fabebe", R = "#fabebe",
    N = '#aa6e28', Q = '#aa6e28', S = '#aa6e28', T = '#aa6e28',
    G = "#46f0f0", P = "#46f0f0", H = "#f032e6", W = "#f032e6",
    Y = "#f032e6", F = "#f032e6", A = "#e6194b", C = "#e6194b",
    M = "#e6194b", L = "#e6194b", I = "#e6194b", V = "#e6194b"),
  symbol = c(
    "Asp" = "D", "Glu" = "E", "Lys" = "K", "Arg" = "R", "Asn" = "N", 
    "Gln" = "Q", "Ser" = "S", "Thr" = "T", "Gly" = "G", "Pro" = "P",
    "His" = "H", "Trp" = "W", "Tyr" = "Y", "Phe" = "F", "Ala" = "A", 
    "Cys" = "C", "Met" = "M", "Leu" = "L", "Ile" = "I", "Val" = "V"),
  group = NULL)

## spatial frequency (Maiorov)
contact_potential_Maiorov <- list(
  color = c(
    D = "#0082c8", E = "#0082c8", N = "#0082c8", Q = "#0082c8",
    K = "#fabebe", R = "#fabebe", G = '#46f0f0', P = "#00FF00",
    A = "#e6194b", V = "#e6194b", S = "#808000", T = "#808000", 
    H = "#808000", W = "#808000", Y = "#808000", C = "#3cb44b", 
    F = "#3cb44b", M = "#3cb44b", L = "#3cb44b", I = "#3cb44b"),
  symbol = c(
    "Asp" = "D", "Glu" = "E", "Asn" = "N", "Gln" = "Q", "Lys" = "K", 
    "Arg" = "R", "Gly" = "G", "Pro" = "P", "Ala" = "A", "Val" = "V", 
    "Ser" = "S", "Thr" = "T", "His" = "H", "Trp" = "W", "Tyr" = "Y",
    "Cys" = "C", "Phe" = "F", "Met" = "M", "Leu" = "L", "Ile" = "I"),
  group = NULL)

## Hydrophobicity (Kyte and Doolittle, 1982)
hydrophobicity_KD_group <- list(
  color = c(
    hydrophilic = "#000000",
    neutral = "#00811B",
    hydrophobic  = "#2000C7"),
  symbol = c(hydrophilic = 'I',
             neutral = 'U',
             hydrophobic = 'O'),
  group = list(hydrophilic = c("D", "E", "K", "R", "N", "Q", "H"),
               neutral = c("G", "P", "S", "T", "W", "Y"),
               hydrophobic = c("A", "C", "F", "I", "L", "M", "V")))

hydrophobicity_HW_group <- list(
  color = c(hydrophobic  = "#000000",
            neutral = "#00811B",
            hydrophilic = "#2000C7"),
  symbol = c(hydrophilic = 'I',
             neutral = 'U',
             hydrophobic = 'O'),
  group = list(hydrophobic = c("W", "F", "Y", "L", "I", "V", "M", "C"),
               neutral = c("H", "A", "T", "P", "G", "N", "Q", "S"),
               hydrophilic = c("R", "K", "D", "E")))

## charge 
charge_group <- list(
  color = c(
    positive = "#FFB32C",
    neutral = "#2000C7",
    negative = "#CCCCCC"),
  symbol = c(
    positive = "P",
    neutral = "U",
    negative = "N"),
  group = list(
    positive = c("H", "K", "R"),
    neutral = c("A", "C", "F", "G", "I", "L", "M", "N", "P", 
                "Q", "S", "T", "V","W","Y"),
    negative = c("D", "E"))) 

## polarity group
polarity_Grantham_group <- list(
  color = c( 
    nonpolar_aliphatic = "#000000",
    nonpolar_aromatic = "#2000C7",
    polar_uncharged = "#00811B",
    polar_positively_charged = "#800080",
    polar_negatively_charged = "#FFB32C"),
  symbol = c(
    nonpolar_aliphatic = "M",
    nonpolar_aromatic = "A",
    polar_uncharged = "U",
    polar_positively_charged = "P",
    polar_negatively_charged = "N"),
  group = list(
    nonpolar_aliphatic = c("A", "V", "L", "G", 
                           "P", "I", "M"),
    nonpolar_aromatic = c("W", "F", "Y"),
    polar_uncharged = c("S", "T", "Q", "C", "N"),
    polar_positively_charged = c("K", "R", "H"),
    polar_negatively_charged = c("D", "E")))

## Bulkiness
bulkiness_Zimmerman_group <- list(
  color = c(tiny = "#0000FF",
            small = "#74008A", 
            medium = "#970067", 
            large = "#B50049",
            gigantic = "#FF0000"),
  symbol = c(tiny = "G", 
             small = "D", 
             medium = "R", 
             large = "M", 
             gigantic = "I"),
  group = list(
    tiny = c("G", "S", "A"),
    small = c("D", "N", "C", "E", "H"),
    medium = c("R", "Q", "K", "T"),
    large = c("M", "P", "Y", "F"),
    gigantic = c("I", "L", "V", "W")))

## Residue volume
volume_Bigelow_group <- list(
  color = c(tiny = "#0000FF", small = "#5200AC", 
            medium = "#7D0081", large = "#9F005F",
            gigantic = "#FF0000"),
  symbol = c(tiny = "G", small = "C", 
             medium = "E", large = "M", 
             gigantic = "F"),
  group = list(tiny = c("G", "A", "S"),
               small = c("C", "D", "T", "P", "N"),
               medium = c("E", "V", "Q", "H"),
               large = c("M", "I", "L", "K", "R"),
               gigantic = c("F", "Y", "W")))

## side chain chemistry
side_chain_chemistry_group <- list(
  color = c( 
    hydrophobic_alipathic = "#000000",
    cyclic = "#5F009F",
    no = "#81007D",
    hydrophobic_aromatic = "#0C00F2",
    polar_neutral = "#800080",
    polar_basic = "#2000C7",
    polar_acidic = "#D00001"),
  symbol = c(
    hydrophobic_alipathic = "H",
    cyclic = "H",
    no = "G",
    hydrophobic_aromatic = "A",
    polar_neutral = "U",
    polar_basic = "B",
    polar_acidic = "A"),
  group = list(
    hydrophobic_alipathic = c("A", "V", "L", "I", "M"),
    cyclic = "P",
    no = "G",
    hydrophobic_aromatic = c("W", "F", "Y"),
    polar_neutral = c("S", "T", "Q", "C", "N"),
    polar_basic = c("K", "R", "H"),
    polar_acidic = c("D", "E")))

## Spatial frequency
contact_potential_Maiorov_group <- list(
  color = c(
    DENQ = "#0082c8",
    KR = "#fabebe",
    G = '#46f0f0',
    P = "#00FF00",
    AV = "#e6194b",
    STHWY = "#808000",
    CFMLI = "#3cb44b"),
  symbol = c(
    DENQ = "D",
    KR = "K",
    G = "G",
    P = "P",
    AV = "A",
    STHWY = "S",
    CFMLI = "C"),
  group = list(
    DENQ = c("D", "E", "N", "Q"), 
    KR = c("K", "R"),
    G = "G",
    P = "P",
    AV = c("A", "V"),
    STHWY = c("S", "T", "H", "W", "Y"),
    CFMLI = c("C", "F", "M", "L", "I"))) 

##  structure alignment
structure_alignments_Mirny_group <- list(
  color = c(
    DE = "#0082c8",
    KR = "#fabebe",
    NQST = '#aa6e28',
    GP = "#46f0f0",
    HWYF = "#f032e6",
    ACMLIV = "#e6194b"),
  symbol = c(
    DE = "D",
    KR = "K",
    NQST = "N",
    GP = "G",
    HWYF = "H",
    ACMLIV = "A"),
  group = list(
    DE = c("D","E"),
    KR = c("K","R"),
    NQST = c("N", "Q", "S", "T"),
    GP = c("G", "P"),
    HWYF = c("H", "W", "Y", "F"),
    ACMLIV = c("A", "C", "M", "L", "I", "V"))) 

## substitution index
sequence_alignment_Dayhoff_group <- list(
  color = c(
    DENQ = "#0082c8",
    KRH = "#fabebe",
    STGPA = "#808000",
    C = '#3cb44b',
    WYF = "#808080",
    MLIV = "#e6beff"),
  symbol = c(
    DENQ = "D",
    KRH = "K",
    STGPA = "S",
    C = "C",
    WYF = "W",
    MLIV = "M"),
  group = list(
    DENQ = c("D","E", "N", "Q"),
    KRH = c("K","R", "H"),
    STGPA  = c("S", "T", "G", "P", "A"),
    C = "C",
    WYF = c("W", "Y", "F"),
    MLIV = c("M", "L", "I", "V"))) 

## chemistry property similarity of individual AAs
chemistry_property_Mahler_group <- list(
  color = c(
    DE = "#0082c8",
    KRH = "#fabebe",
    QN = "#800000",
    ST = "#808000",
    P = "#00FF00",
    CM = '#3cb44b',
    WYF = "#808080",
    GALIV = "#46f0f0"),
  symbol = c(
    DE = "D",
    KRH = "K",
    QN = "Q",
    ST = "S",
    P = "P",
    CM = "C",
    WYF = "W",
    GALIV = "G"),
  group = list(
    DE= c("D","E"),
    KRH = c("K", "R", "H"),
    QN = c("Q", "N"),
    ST  = c("S", "T"),
    P = "P",
    CM = c("C", "M"),
    WYF = c("W", "Y", "F"),
    GALIV = c("G", "A", "L", "I", "V")))

## consensus AA similarity
consensus_similarity_SF_group <- list(
  color = c(
    DEN = "#0082c8",
    KRQ = "#fabebe",
    ST = "#808000",
    H = "#f032e6",
    GP = "#46f0f0",
    A = "#e6194b",
    C = '#3cb44b',
    WYF = "#808080",
    MLIV = "#e6beff"),
  symbol = c(
    DEN = "D",
    KRQ = "K",
    ST = "S",
    H = "H",
    GP = "G",
    A = "A",
    C = "C",
    WYF = "W",
    MLIV = "M"),
  group = list(
    DEN= c("D","E", "N"),
    KRQ = c("K","R", "Q"),
    ST  = c("S", "T"),
    H = "H",
    GP = c("G", "P"),
    A = "A",
    C = "C",
    WYF = c("W", "Y", "F"),
    MLIV = c("M", "L", "I", "V")))

## no grouping, coloring based on properties of individual AAs
cachedEnv$auto <- auto
cachedEnv$no <- no
cachedEnv$hydrophobicity_KD <- hydrophobicity_KD
cachedEnv$hydrophobicity_HW <- hydrophobicity_HW
cachedEnv$polarity_Grantham <- polarity_Grantham
cachedEnv$isoelectric_point_Zimmerman <- isoelectric_point_Zimmerman
cachedEnv$bulkiness_Zimmerman <- bulkiness_Zimmerman
cachedEnv$volume_Bigelow <- volume_Bigelow
cachedEnv$chemistry_property_Mahler <- chemistry_property_Mahler
cachedEnv$consensus_similarity_SF <- consensus_similarity_SF

## no grouping, coloring based on properties of AAs withing proteins
## (substitution and spatial frequency)
cachedEnv$contact_potential_Maiorov <- contact_potential_Maiorov
cachedEnv$structure_alignments_Mirny <- structure_alignments_Mirny
cachedEnv$sequence_alignment_Dayhoff <- sequence_alignment_Dayhoff

## grouping based on properties of individual AAs
cachedEnv$hydrophobicity_KD_group <- hydrophobicity_KD_group
cachedEnv$hydrophobicity_HW_group <- hydrophobicity_HW_group
cachedEnv$polarity_Grantham_group <- polarity_Grantham_group
cachedEnv$charge_group <- charge_group
cachedEnv$bulkiness_Zimmerman_group <- bulkiness_Zimmerman_group
cachedEnv$volume_Bigelow_group <- volume_Bigelow_group
cachedEnv$chemistry_property_Mahler_group <- chemistry_property_Mahler_group
cachedEnv$consensus_similarity_SF_group <- consensus_similarity_SF_group

## grouping based on properties of AAs within proteins 
## (substitution and spatial frequency)
cachedEnv$contact_potential_Maiorov_group <- contact_potential_Maiorov_group
cachedEnv$structure_alignments_Mirny_group <- structure_alignments_Mirny_group
cachedEnv$sequence_alignment_Dayhoff_group <- sequence_alignment_Dayhoff_group


