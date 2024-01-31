suppressPackageStartupMessages(library("argparse"))


# Command line arguments
parser <- ArgumentParser()

parser$add_argument("--nperm", default=500, type="integer",
                    help = "Number of permutations for each test [default %(default)s]")

parser$add_argument("--span", default=0.05, type="double",
                    help="The span parameter defining the width of local kernel [default %(default)s]")

parser$add_argument("exp_file", nargs=1, help="Input TSV file of gene expressions")

parser$add_argument("coord_file", nargs=1, help="Input TSV file of barcode coordinates")

parser$add_argument("out_path", nargs=1, help="Path to output TSV file")

args <- parser$parse_args()


# Load input data
suppressPackageStartupMessages(library("scHOT"))
suppressPackageStartupMessages(library("tidyverse"))

bc_exp <- read_tsv(args$exp_file)
bc_pos <- read_tsv(args$coord_file) %>% as.data.frame

bc_exp_mat <- t(bc_exp[,-1])
colnames(bc_exp_mat) <- bc_exp$barcode


# Build scHOT object
scHOT_spatial <- scHOT_buildFromMatrix(
    bc_exp_mat, 
    bc_pos, 
    "spatial", 
    c("imagerow","imagecol")
)

pairs <- t(combn(rownames(bc_exp_mat),2))
rownames(pairs) <- apply(pairs,1,paste0,collapse = "_")


# Run testing
set.seed(2021)
scHOT_spatial <- scHOT(
    scHOT_spatial,
    testingScaffold = pairs,
    positionType = "spatial",
    positionColData = c("imagerow","imagecol"),
    nrow.out = NULL,
    higherOrderFunction = weightedSpearman,
    higherOrderFunctionType = "weighted",
    #numberPermutations = args$nperm,
    #numberPermutations = 10,
    numberPermutations = 1000,
    higherOrderSummaryFunction = sd,
    parallel = FALSE,
    verbose = FALSE,
    span = args$span
)


# Write output
df <- as.data.frame(scHOT_spatial@scHOT_output)
drop <- c("higherOrderSequence","permutations")
df = df[,!(names(df) %in% drop)]
write_tsv(df, args$out_path)

#scHOT_spatial@scHOT_output %>% as.data.frame %>% write_tsv(args$out_path)
