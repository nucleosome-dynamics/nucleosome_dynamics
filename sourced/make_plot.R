#!/usr/bin/Rscript

library(IRanges)

makePlot <- function(dyn, plot.start, plot.end)
{   # Given an already subsetted dynamic, and its start and end, plot it

    initPlot(x0 = plot.start,
             y0 = 0,
             x1 = plot.end,
             y1 = getMaxVal(dyn))

    addCovs(dyn[[1]]$originals, dyn[[2]]$originals)
    addIndels(dyn[[1]]$indels, dyn[[2]]$indels)

    par.usr <- par("usr")
    par.ypc <- (par.usr[4] - par.usr[3]) / 100
    middle <- (par.usr[4] - par.usr[3]) * 0.5

    addLeft(dyn[[1]]$left.shifts,
            dyn[[2]]$left.shifts,
            middle,
            par.ypc)
    addRight(dyn[[1]]$right.shifts,
             dyn[[2]]$right.shifts,
             middle,
             par.ypc)

    makeLegend()
}

getMaxVal <- function(xs)
    # Get the maximum value that will be plotted
    max(unlist(lapply(xs, function(x) as.vector(coverage(x$originals)))))

initPlot <- function(x0, y0, x1, y1)
    # Initiate an empty plot with the specified size
    plot(1,
         type='n',
         xlim=c(x0, x1),
         ylim=c(y0, y1),
         xlab="",
         ylab="",
         main="")

addCovs <- function(x, y)
{   # Plot coverages for the first and second experiment.
    # First experiment will be plotted as a grey area, and second experiment
    # wil be plotted as a black dotted line.
    lines(coverage(x), col="grey", type='h', lty=1, lwd=2)
    lines(coverage(y), col="black", lty=3, lwd=2)
}

addIndels <- function(del, ins)
{   # Plot the coverages for insertions and delitons.
    # Insertions will be plotted as a green area and deletions as a red one.
    lines(coverage(del), type='h', lwd=2, col="red")
    lines(coverage(ins), type='h', lwd=2, col="green")
}

dyadPos <- function(ran)
    # Find the center points of an IRanges
    round((start(ran) + end(ran)) / 2)

addShifts <- function(xs, ys, middle, par.ypc, col, fromMid)
{   # Add shifts
    start <- dyadPos(xs)
    end <- dyadPos(ys)

    db <- disjointBins(IRanges(start = mapply(min, start, end),
                               end   = mapply(max, start, end)))

    vpos <- fromMid(middle, par.ypc*db)

    arrows(x0     = start,
           y0     = vpos,
           x1     = end,
           y1     = vpos,
           length = 0.05,
           col    = col,
           angle  = 30,
           code   = 2,
           lwd    = 1)
}

addLeft <- function(xs, ys, middle, par.ypc)
    # Add left shifts
    addShifts(xs        = xs,
              ys        = ys,
              middle    = middle,
              par.ypc   = par.ypc,
              col       = "darkred",
              fromMid   = `-`)

addRight <- function(xs, ys, middle, par.ypc)
    # Add right shifts
    addShifts(xs        = xs,
              ys        = ys,
              middle    = middle,
              par.ypc   = par.ypc,
              col       = "darkblue",
              fromMid   = `+`)

makeLegend <- function()
    # Add a legend
    legend("topleft",
           legend = c("Coverage 1",
                      "Coverage 2",
                      "Inserted reads",
                      "Removed reads",
                      "Upstream shifts",
                      "Downstream shifts"),
           col    = c("grey",
                      "black",
                      "green",
                      "red",
                      "darkred",
                      "darkblue"),
           lty    = c(1, 3, 1, 1, 1, 1),
           lwd    = 3,
           bty    = "n",
           0.8)
