#!/usr/bin/env Rscript
# Build an ampliseq-notebook-compatible phyloseq object from Emu rel-abundance tables.
# Research use only. Not for clinical or diagnostic use.

parse_args <- function(args) {
    parsed <- list(
        abundances = character(0),
        metadata = NULL,
        output = "emu_phyloseq.rds"
    )
    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        value <- if (i < length(args)) args[[i + 1]] else NA
        if (key == "--abundances") {
            files <- character(0)
            j <- i + 1
            while (j <= length(args) && !startsWith(args[[j]], "--")) {
                files <- c(files, args[[j]])
                j <- j + 1
            }
            parsed$abundances <- files
            i <- j
            next
        } else if (key == "--metadata") {
            parsed$metadata <- value
        } else if (key == "--output") {
            parsed$output <- value
        } else {
            stop("Unknown argument: ", key)
        }
        i <- i + 2
    }
    parsed
}

sample_id_from_path <- function(path) {
    name <- basename(path)
    name <- sub("_rel-abundance\\.tsv$", "", name)
    name <- sub("-rel-abundance\\.tsv$", "", name)
    tools::file_path_sans_ext(name)
}

col_named <- function(df, name) {
    idx <- which(tolower(names(df)) == tolower(name))
    if (length(idx) == 0) {
        return(NULL)
    }
    df[[idx[[1]]]]
}

is_unassigned <- function(tax_id) {
    value <- trimws(as.character(tax_id))
    tolower(value) %in% c("", "na", "nan", "none", "null")
}

read_emu_table <- function(path) {
    table <- tryCatch(
        utils::read.delim(path, check.names = FALSE, comment.char = "#", stringsAsFactors = FALSE),
        error = function(e) NULL
    )
    if (is.null(table) || nrow(table) == 0) {
        return(NULL)
    }
    table
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package phyloseq is required")
}

abundance_files <- args$abundances
abundance_files <- abundance_files[file.exists(abundance_files)]
if (length(abundance_files) == 0) {
    stop("No Emu *_rel-abundance.tsv files were provided")
}

rank_map <- c(
    superkingdom = "Kingdom",
    phylum = "Phylum",
    class = "Class",
    order = "Order",
    family = "Family",
    genus = "Genus",
    species = "Species"
)
rank_names <- unname(rank_map)

counts <- list()
taxonomy <- list()
sample_ids <- character(0)

for (path in abundance_files) {
    table <- read_emu_table(path)
    if (is.null(table)) {
        next
    }
    tax_ids <- col_named(table, "tax_id")
    if (is.null(tax_ids)) {
        stop("Emu report is missing tax_id: ", path)
    }
    metric <- col_named(table, "estimated counts")
    if (is.null(metric)) {
        stop("Emu report is missing 'estimated counts' (run emu abundance --keep-counts): ", path)
    }
    sid <- sample_id_from_path(path)
    sample_ids <- c(sample_ids, sid)
    metric <- as.numeric(metric)
    metric[is.na(metric)] <- 0
    for (i in seq_len(nrow(table))) {
        tax_id <- as.character(tax_ids[[i]])
        if (is_unassigned(tax_id)) {
            next
        }
        if (is.null(counts[[tax_id]])) {
            counts[[tax_id]] <- list()
        }
        prev <- counts[[tax_id]][[sid]]
        counts[[tax_id]][[sid]] <- if (is.null(prev)) metric[[i]] else prev + metric[[i]]
        if (is.null(taxonomy[[tax_id]])) {
            lineage <- character(length(rank_names))
            names(lineage) <- rank_names
            for (src in names(rank_map)) {
                values <- col_named(table, src)
                if (!is.null(values)) {
                    value <- as.character(values[[i]])
                    if (!is.na(value) && !tolower(value) %in% c("", "nan", "none", "null")) {
                        lineage[[rank_map[[src]]]] <- value
                    } else {
                        lineage[[rank_map[[src]]]] <- NA_character_
                    }
                } else {
                    lineage[[rank_map[[src]]]] <- NA_character_
                }
            }
            taxonomy[[tax_id]] <- lineage
        }
    }
}

sample_ids <- unique(sample_ids)
tax_ids <- names(counts)
if (length(tax_ids) == 0 || length(sample_ids) == 0) {
    stop("No taxa with estimated counts were available to build phyloseq")
}

otu_mat <- matrix(0, nrow = length(tax_ids), ncol = length(sample_ids), dimnames = list(tax_ids, sample_ids))
for (tax_id in tax_ids) {
    for (sid in names(counts[[tax_id]])) {
        otu_mat[tax_id, sid] <- counts[[tax_id]][[sid]]
    }
}
otu_mat <- round(otu_mat)
storage.mode(otu_mat) <- "integer"
otu <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat <- matrix(NA_character_, nrow = length(tax_ids), ncol = length(rank_names), dimnames = list(tax_ids, rank_names))
for (tax_id in tax_ids) {
    tax_mat[tax_id, ] <- taxonomy[[tax_id]]
}
tax <- phyloseq::tax_table(tax_mat)

sample_data_obj <- NULL
if (!is.null(args$metadata) && nzchar(args$metadata) && file.exists(args$metadata) && file.info(args$metadata)$size > 0) {
    md <- utils::read.delim(args$metadata, check.names = FALSE, comment.char = "#", stringsAsFactors = FALSE)
    if (ncol(md) >= 1 && nrow(md) >= 1) {
        id_col <- names(md)[[1]]
        if (tolower(id_col) %in% c("id", "sample")) {
            names(md)[[1]] <- "ID"
            id_col <- "ID"
        }
        rownames(md) <- as.character(md[[id_col]])
        common <- intersect(rownames(md), colnames(otu_mat))
        if (length(common) > 0) {
            md <- md[common, , drop = FALSE]
            otu <- phyloseq::prune_samples(common, otu)
            sample_data_obj <- phyloseq::sample_data(md)
        }
    }
}

# phy_tree is omitted: ttfampliseq v2.0.0 falls back from weighted UniFrac to Bray-Curtis.
components <- list(otu, tax)
if (!is.null(sample_data_obj)) {
    components <- c(components, list(sample_data_obj))
}

ps <- do.call(phyloseq::phyloseq, components)
dir.create(dirname(args$output), showWarnings = FALSE, recursive = TRUE)
saveRDS(ps, args$output)
message("Wrote ", args$output, " with ", phyloseq::ntaxa(ps), " taxa and ", phyloseq::nsamples(ps), " samples")
