#!/usr/bin/Rscript

# Functions to read, parse and write files in the gff format

library(IRanges)
library(GenomicRanges)

rd2df <- function (rd)
{   # Convert a RangedData into a data.frame to ease later conversion into
    # gff format.
    # It really just sets the the "space" field to "seqname", and removes the
    # field "width". Leaves everything else unchanged.
    df <- as.data.frame(rd)
    names(df) <- sub("space", "seqname", names(df))
    df[["width"]] <- NULL
    df
}

df2rd <- function (df)
{   # Convert from data.frame to RangedData
    names(df) <- sub("seqname", "space", names(df))
    RangedData(df)
}

df2gff <- function (df, ...)
{   # Convert a data.frame into a form that is easily saved in disk as a gff.
    # It tries to match gff fields to the row names of the data.frame.
    # gff fields not found in the data.frame will be set to `.` and fields in
    # the data.frame not found in the gff fields will be stored accordingly in
    # the `attribute` field, separated by `;`.
    # It also accepts optional keyword arguments to specify gff fields that
    # are not present in the data.frame.
    kwargs <- list(...)
    fields <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame")
    out.df <- data.frame(matrix(nrow=nrow(df),
                                ncol=length(fields) + 1,
                                dimnames=list(c(),
                                              c(fields, "attribute"))))

    for (f in fields) {
        if (f %in% colnames(df)) {
            out.df[[f]] <- as.vector(df[[f]])
        } else if (f %in% names(kwargs)) {
            if (length(kwargs[[f]]) == 1) {
                out.df[[f]] <- rep(kwargs[[f]], nrow(out.df))
            } else {
                out.df[[f]] <- kwargs[[f]]
            }
        } else {
            out.df[[f]] <- rep(".", nrow(out.df))
        }
    }
    nonfield.columns <- colnames(df)[!colnames(df) %in% fields]
    attrVal <- function(i, df) sprintf("%s=%s", i, df[[i]])
    attrs <- do.call(paste,
                     c(lapply(nonfield.columns,
                              attrVal,
                              df),
                       sep=";"))
    if (length(attrs)) {
        out.df[["attribute"]] <- attrs
    } else {
        out.df[["attribute"]] <- rep(".", nrow(out.df))
    }
    if (nrow(out.df) > 0) {
        out.df
    } else {
        tmp <- rep(".", length(out.df))
        names(tmp) <- names(out.df)
        as.data.frame(as.list(tmp))
    }
}

writeGff <- function (df, outpath)
{
    # Use this to write the output of df2gff to disk.
    write.table(df,
                file      = outpath,
                quote     = FALSE,
                sep       = "\t",
                row.names = FALSE,
                col.names = FALSE)
}

readGff <- function (fname, load.attributes=TRUE)
{   # Read a gff file
    cols <- c("seqname",
              "source",
              "feature",
              "start",
              "end",
              "score",
              "strand",
              "frame",
              "attribute")

    df <- read.csv(fname,
                   header           = FALSE,
                   col.names        = cols,
                   stringsAsFactors = FALSE,
                   sep              = "\t")

    for (i in colnames(df)) {
        if (all(df[[i]] == ".")) {
            df[[i]] <- NULL
        }
    }

    if (load.attributes) {
        attrs <- df$attribute

        parseRowAttrs <- function(x) {
            pairs <- strsplit(x, "=")  # separate name-value pairs
            not.two <- sapply(pairs, length) < 2  # some values might be missing
            pairs[not.two] <- lapply(pairs[not.two], c, NA)  # set them to NA
            ns <- sapply(pairs, `[[`, 1)  # get the names
            vs <- lapply(pairs, `[[`, 2)  # get the values
            names(vs) <- ns
            vs
        }

        if (!is.null(attrs)) {
            row.attrs <- lapply(strsplit(attrs, ";"), parseRowAttrs)
            # make sure they all have the same attributes, even if some are NA
            params <- unique(unlist(lapply(row.attrs, names)))
            for (i in seq_along(row.attrs)) {
                xnames <- names(row.attrs[[i]])
                na.xnames <- params[!(params %in% xnames)]
                row.attrs[[i]][na.xnames] <- NA
            }

            m <- do.call(rbind, row.attrs)
            df.tmp <- as.data.frame(m)
            attrs.df <- as.data.frame(lapply(df.tmp, unlist),
                                      stringsAsFactors=FALSE)
            df <- cbind(df, attrs.df)
        }
    }

    df$attribute <- NULL

    if (ncol(df) > 0) {
        df
    } else {
        data.frame()
    }
}
