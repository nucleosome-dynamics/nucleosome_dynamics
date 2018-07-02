#!/usr/bin/Rscript

library(ggplot2)
library(IRanges)
library(plotly)
library(plyr)
library(reshape2)

dyadPos <- function(ran)
    # Find the center points of an IRanges
    round((start(ran) + end(ran)) / 2)

makeEqual <- function (xs) {
    maxlen <- max(sapply(xs, length))
    lapply(
        xs,
        function (x)
            c(x, rep(0, maxlen - length(x)))
    )
}

buildArrowDf <- function (xs, ys, middle, par.ypc, type="left") {
    if (length(xs) == 0 | length(ys) == 0) {
        data.frame(x0=vector(),
                   x1=vector(),
                   y0=vector(),
                   y1=vector(),
                   variable=factor(c(), levels=type))
    } else {
        start <- dyadPos(xs)
        end <- dyadPos(ys)
        db <- disjointBins(IRanges(start = mapply(min, start, end),
                                   end   = mapply(max, start, end)))
        if (type == "left") {
            vpos <- middle - par.ypc*db
        } else if (type == "right") {
            vpos <- middle + par.ypc*db
        }
        data.frame(x0=start,
                   x1=end,
                   y0=vpos,
                   y1=vpos,
                   variable=type)
    }
}

arrowHead <- function (x0, y0, l, par.ypc, forward=FALSE) {
    a <- sin(pi/4) * l
    if (forward) {
        list(list(x0=x0, y0=y0, x1=x0 - 10, y1=y0 + a*par.ypc),
             list(x0=x0, y0=y0, x1=x0 - 10, y1=y0 - a*par.ypc))
    } else {
        list(list(x0=x0, y0=y0, x1=x0 + 10, y1=y0 + a*par.ypc),
             list(x0=x0, y0=y0, x1=x0 + 10, y1=y0 - a*par.ypc))
    }
}

addArrows <- function (sdf, par.ypc, forward=TRUE) {
    if (nrow(sdf) > 0) {
        adply(
            sdf,
            1,
            function (row) {
                head <- arrowHead(row$x1,
                                  row$y1,
                                  0.5,
                                  par.ypc,
                                  forward=forward)
                var <- as.vector(row$variable)
                head.rows <- lapply(lapply(head,
                                           c,
                                           variable=var),
                                    as.data.frame)
                do.call(rbind, c(list(row), head.rows))
            }
        )
    } else {
        data.frame(x0=integer(),
                   y0=integer(),
                   x1=integer(),
                   y1=integer(),
                   variable=factor(levels=levels(sdf$variable)))
    }
}

colors <- c("setA"  = "dimgray",
            "setB"  = "gray",
            "ins"   = "green",
            "dels"  = "red",
            "left"  = "darkred",
            "right" = "darkblue")

buildGgplot <- function (mdf, shdf, plot.start, plot.end) {
    ggplot(mdf, aes(x=x, y=value)) +
        xlim(plot.start, plot.end) +
        xlab("") +
        ylab("") +
        geom_area(aes(fill=variable,
                      color=variable,
                      alpha=variable),
                  position="identity") +
        geom_segment(data=shdf,
                     mapping=aes(x=x0, y=y0, xend=x1, yend=y1,
                                 color=variable)) +
        scale_color_manual(name="", values=colors) +
        scale_fill_manual(name="", values=colors) +
        scale_alpha_manual(name="",
                           values=c("setA"=1,
                                    "setB"=0.7,
                                    "ins"=0.6,
                                    "dels"=0.6,
                                    "left"=1,
                                    "right"=1)) +
        theme_bw()
}

fixShifts <- function (xs, name) {
    parseShiftTxt <- function (txt)
        lapply(strsplit(txt, split="<br>"),
               function (x) {
                   pairs <- strsplit(x, split=": ")
                   names <- sapply(pairs, `[`, 1)
                   vals <- sapply(pairs, `[`, 2)
                   names(vals) <- names
                   vals
               })

    txt.ls <- parseShiftTxt(xs[["text"]])
    n <- length(xs[["text"]])
    for (i in seq(from=1, to=n, by=9)) {
        x0 <- txt.ls[[i]]["x0"]
        x1 <- txt.ls[[i]]["x1"]
        for (j in i:(i+8)) {
            if (j <= n) {
                xs[["text"]][j] <- sprintf("%s<br>from: %s<br>to: %s",
                                           name, x0, x1)
            }
        }
    }
    xs[["name"]] <- paste0(name, "s")
    xs
}

fixNonShifts <- function (xs, name) {
    xs[["name"]] <- name
    xs[["text"]] <- sprintf("%s<br>position: %d<br>coverage: %d",
                            name, xs[["x"]], xs[["y"]])
    xs
}

ggplot2widget <- function (p) {
    pdf(NULL)
    pl <- plotly_build(p)

    names <- c("coverage 1", "Coverage 2", "Deletions", "Insertions")
    for (i in seq_along(names)) {
        pl[["data"]][[i]] <- fixNonShifts(pl[["data"]][[i]], names[i])
    }

    pl[["data"]][[5]] <- fixShifts(pl[["data"]][[5]], "Upstream shift")
    pl[["data"]][[6]] <- fixShifts(pl[["data"]][[6]], "Downstream shift")

    as.widget(pl)
}

makeShifts <- function (dyn, middle, par.ypc) {
    lsdf <- buildArrowDf(dyn[[1]]$left.shifts,
                         dyn[[2]]$left.shifts,
                         middle,
                         par.ypc,
                         "left")
    rsdf <- buildArrowDf(dyn[[1]]$right.shifts,
                         dyn[[2]]$right.shifts,
                         middle,
                         par.ypc,
                         "right")
    lsdf <- addArrows(lsdf, par.ypc, forward=FALSE)
    rsdf <- addArrows(rsdf, par.ypc, forward=TRUE)
    rbind(lsdf, rsdf)
}

simplifier <- function (df, i=2)
    # remove some points from the plot for efficiency
    if (nrow(df)) {
        df[1:nrow(df) %% i == 1, ]
    } else {
        df
    }

makeCovs <- function (dyn) {
    cols <- list(setA=dyn[[1]]$originals,
                 setB=dyn[[2]]$originals,
                 dels=dyn[[1]]$indels,
                 ins=dyn[[2]]$indels)
    covs <- lapply(cols, function(x) as.vector(coverage(x)))
    valsdf <- as.data.frame(makeEqual(covs))
    n <- nrow(valsdf)
    if (n) {
        valsdf$x <- 1:n
    } else {
        valsdf$x <- integer()
    }
    melt(simplifier(valsdf), id.vars="x")
}

buildPlot <- function (dyn, start, end) {
    mdf <- makeCovs(dyn)
    if (nrow(mdf)) {
        maxy <- max(mdf$value)
    } else {
        maxy <- 0
    }
    middle <- maxy / 2
    par.ypc <- maxy / 100
    shdf <- makeShifts(subdyn, middle, par.ypc)
    gg <- buildGgplot(mdf, shdf, start, end)
    ggplot2widget(gg)
}
