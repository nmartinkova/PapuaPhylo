for (i in c("phytools", "ape", "geiger", "adephylo")) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}

if (!dir.exists("Results")) dir.create("Results")
if (!dir.exists("Data")) dir.create("Data")


# input data
subory <- dir("Data", full.names = T)

vzorky <- read.table("Data/TableS1.txt",
  header = T, sep = "\t",
  quote = "\"", comment.char = ""
)


pdf(
  paste0(
    "Results/FigureS2",
    # switch(ktore,    Ratt = "-Rattini",    Hydro = "-Hydromyini"  ),
    ".pdf"
  ),
  width = 2 * 4, height = 3.5 * 3
)
layout(matrix(c(1, 1, 2, 3, 3, 4), ncol = 2))
par(mar = c(4, 4, 2, 0) + .3)


for (ktore in c("Hydro", "Ratt")) {
  # load phylogeny
  strom <- read.tree(subory[grepl(ktore, subory) & grepl("treefile", subory)])

  # find tree root - Cizkova email 251126
  # plot(strom, cex = .4)
  # identify(strom)
  cisloKorena <- ifelse(ktore == "Ratt", 569, 526)

  strom2 <- phytools::reroot(strom, node.number = cisloKorena)

  # delete introduced Rattus
  if (ktore == "Ratt") {
    ktoreTips <- phytools::getDescendants(strom, cisloKorena)
    ktoreTips <- ktoreTips[ktoreTips <= Ntip(strom)]
    strom2 <- keep.tip(strom2,
      tip = strom$tip.label[ktoreTips]
    )
  }


  # save rooted tree
  write.tree(strom2, file = paste0("Data/PNGrooted", switch(ktore,
    Ratt = "Rattini",
    Hydro = "Hydromyini"
  ), ".tre"))

  # rescale to 1 for histogram
  strom3 <- rescale(strom2, model = "depth", depth = 1)

  # plot phylogeny
  plot(strom2, show.tip.label = F, show.node.label = F)
  add.scale.bar()
  mtext(switch(ktore,
    Hydro = "A",
    Ratt = "B"
  ), side = 3, adj = 0, font = 2)

  # plot root-to-tip distance histogram
  tipDist <- adephylo::distRoot(strom3, method = "patristic")
  hist(tipDist, main = "", las = 1, xlab = "Root-to-tip distance", breaks = 15)
  box()
}
dev.off()
