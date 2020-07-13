# Circos plot tutorial

rm(list = ls())

library(ggbio)

data("CRC", package = "biovizBase")

autoplot(hg19sub, layout = "circle")

gr.crc1 <- crc.gr[values(crc.gr)$individual == "CRC-1"]

ggbio() +
  circle(mut.gr, geom = "rect", color = "steelblue") +
  circle(hg19sub, geom = "ideo") +
  circle(hg19sub, geom = "scale", size = 2) +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 2) +
  circle(gr.crc1, geom = "point", aes(x = score, y = tumreads), color = "red", grid = TRUE, radius = 30) +
  circle(gr.crc1, geom = "link", linked.to = "to.gr", aes(color = rearrangements), radius = 23) +
  scale_size(range = c(1, 2.5))

class(mut.gr)
class(hg19sub)

hg19sub
mut.gr
gr.crc1

p <- ggplot(mtcars, aes(x = mpg, y = cyl)) + geom_point()
ggbio(p)

N <- 100
library(GenomicRanges)
## ======================================================================
##  simmulated GRanges
## ======================================================================
gr <- GRanges(
  seqnames = sample(c("chr1", "chr2", "chr3"), size = N, replace = TRUE),
  ranges = IRanges(
    start = sample(1:300, size = N, replace = TRUE),
    width = sample(70:75, size = N,replace = TRUE)
  ),
  strand = sample(c("+", "-", "*"), size = N, replace = TRUE),
  value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
  sample = sample(c("Normal", "Tumor"), size = N, replace = TRUE),
  pair = sample(letters, size = N, replace = TRUE)
)

ggbio() + circle(gr)

seqlengths(gr) <- c(400, 500, 700)
values(gr)$to.gr <- gr[sample(1:length(gr), size = length(gr))]

## doesn't pass gr to the ggplot
ggplot() + layout_circle(gr, geom = "ideo", fill = "gray70", radius = 7, trackWidth = 3) +
  layout_circle(gr, geom = "bar", radius = 10, trackWidth = 4, aes(fill = score, y = score)) +
  layout_circle(gr, geom = "point", color = "red", radius = 14,
                trackWidth = 3, grid = TRUE, aes(y = score)) +
  layout_circle(gr, geom = "link", linked.to = "to.gr", radius = 6,
                trackWidth = 1)

## more formal API
ggplot(gr) + layout_circle(geom = "ideo", fill = "gray70", radius = 7, trackWidth = 3) +
  layout_circle(geom = "bar", radius = 10, trackWidth = 4, aes(fill = score, y = score)) +
  layout_circle(geom = "point", color = "red", radius = 14,
                trackWidth = 3, grid = TRUE, aes(y = score)) +
  layout_circle(geom = "link", linked.to = "to.gr", radius = 6, trackWidth = 1)
