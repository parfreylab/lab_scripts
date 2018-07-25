# Normalization of count data
library(tools)
library(optparse)

# like argparse from Python
options.list <- list(make_option(c("-f",
                                   "--file_to_normalize"),
                                 help="Path to input count data to normalize, string."),
                     make_option(c("-c",
                                   "--category_column"),
                                 help="Name of the column header in the input file which contains protein groups, string. Case sensitive."),
                     make_option(c("-o",
                                   "--output_path"),
                                 default="./",
                                 help="Path to output folder, string. Default: current working directory."))

opt.parser <- OptionParser(option_list = options.list)
opt <- parse_args(opt.parser)

path <- opt$file_to_normalize
category <- opt$category_column
out <- opt$output_path

# check if any of the arguments are null. if so, stop the script and report that they are required
if (is.null(path)) {
  stop("count data file must be specified. see script usage (--help)")
} 
if (is.null(category)) {
  stop("protein description category from input file must be specified. see script usage (--help)")
}
if (!endsWith(out, "/")) {
  out <- paste(out, "/", sep = "")
}

# get the input filename without the extension
file.toread.noext <- file_path_sans_ext(basename(path))

# Read in the count data
to.norm <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE)

# Start a count with 1 to iterate through the columns
count <- 1

# Loop through each sample column and normalize the count data
for (i in colnames(to.norm)) {
  if (!grepl(category, i)) { # if the current column is not for protein descriptions/families/categories
    to.norm[, count] <- as.numeric(as.character(to.norm[, count])) / sum(as.numeric(as.character(to.norm[, count])))
  }
  count <- count + 1
}

write("Successfully normalized count data!", stdout())

# Write a tab-delimited text file with the normalized data
write.table(to.norm, file = paste(out, file.toread.noext, "_norm.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
