args <- commandArgs(trailingOnly = TRUE)
input_csv <- args[1]
output_txt <- args[2]

if (length(args) < 2) {
  stop("Usage: Rscript extract_5ss.R <input_csv> <output_txt>")
}

input_csv <- args[1]
output_txt <- args[2]

suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

regions <- read.csv(input_csv, stringsAsFactors = FALSE)

required_cols <- c(
  "row_id",
  "Chromosome.scaffold.name",
  "Strand",
  "X5ss_start",
  "X5ss_end"
)

missing_cols <- setdiff(required_cols, colnames(regions))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

regions$Strand <- ifelse(regions$Strand == 1 | regions$Strand == "1", "+", "-")

gr <- GRanges(
  seqnames = regions$Chromosome.scaffold.name,
  ranges = IRanges(start = regions$X5ss_start, end = regions$X5ss_end),
  strand = regions$Strand
)

genome <- BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(gr) <- seqlevelsStyle(genome)

sequences <- getSeq(genome, gr)

output_df <- data.frame(
  row_id = regions$row_id,
  sequence = as.character(sequences),
  stringsAsFactors = FALSE
)

write.table(
  output_df,
  file = output_txt,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)