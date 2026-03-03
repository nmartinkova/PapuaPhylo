for (i in c("phytools", "ape", "geiger", "OUwie")) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}

if (!dir.exists("Results")) dir.create("Results")


ouwie <- TRUE

# input files
subory <- dir("Data", full.names = T)

vzorky <- read.table("Data/TableS1.txt",
  header = T, sep = "\t",
  quote = "\"", comment.char = ""
)




for (ktore in c("Ratt", "Hydro")) {
  strom <- read.tree(subory[grepl(ktore, subory) & grepl("ini", subory)])
  # calibration
  calib <- switch(ktore,
    # Rowe et al. 2019, Table 4, DOI: 10.1111/jbi.13720
    Hydro = makeChronosCalib(strom, node = "root", age.min = 6.72, age.max = 7.74, soft.bounds = F),
    Ratt = makeChronosCalib(strom, node = "root", age.min = 0.87, age.max = 1.53, soft.bounds = F)
  )
  calib$age.start <- switch(ktore,
    Hydro = mean(c(6.72, 7.74)),
    Ratt = mean(c(0.87, 1.53))
  )

  # ultrametric tree
  strom <- chronos(strom,
    calibration = calib,
    model = "discrete",
    control = chronos.control(nb.rate.cat = 3)
  )
  write.tree(strom,
    file = paste0("Results/PNGchronogram", switch(ktore,
      Ratt = "Rattini",
      Hydro = "Hydromyini"
    ), ".tre")
  )
  capture.output(print(str(strom)), file = paste0("Results/PNGchronogram", switch(ktore,
    Ratt = "Rattini",
    Hydro = "Hydromyini"
  ), ".txt"))

  strom <- read.tree(
    file = paste0("Results/PNGchronogram", switch(ktore,
      Ratt = "Rattini",
      Hydro = "Hydromyini"
    ), ".tre")
  )

  # variables to test
  elevation <- setNames(
    vzorky$Elevation[match(sub("_.+", "", strom$tip.label), vzorky$ID)],
    strom$tip.label
  )

  region <- setNames(
    vzorky$Region[match(sub("_.+", "", strom$tip.label), vzorky$ID)],
    strom$tip.label
  )
  region[grepl("Bai|Wan", region)] <- "Lowlands"
  region[grep("Finis", region)] <- "Saruwaged"

  # elevation
  ancElevation <- phytools::fastAnc(strom, elevation, vars = TRUE, CI = TRUE)
  plotElevation <- phytools::contMap(strom, elevation)
  plotElevation <- phytools::setMap(plotElevation, colors = rev(hcl.colors(7, "inferno")))


  # regions
  ancRegion <- ape::ace(region, strom, type = "discrete")


  ####### Region uncertainty


  # check which model fits best
  fitER <- fitDiscrete(strom, region, model = "ER")
  fitSYM <- fitDiscrete(strom, region, model = "SYM")
  fitARD <- fitDiscrete(strom, region, model = "ARD")
  bestModel <- switch(ktore,
    Hydro = "ARD",
    Ratt = "ER"
  )

  # stochastic map simulations
  simulacie <- phytools::make.simmap(strom, region, nsim = 100, model = bestModel, message = FALSE)
  suhrn <- phytools::describe.simmap(simulacie)
  capture.output(print(AIC(fitER, fitSYM, fitARD)), print(suhrn),
    file = paste0("Results/simmap-", switch(ktore,
      Hydro = "Hydromyini",
      Ratt = "Rattini"
    ), ".txt")
  )


  # draw trees
  pdf(paste0("Results/Figure", switch(ktore,
    Hydro = "2",
    Ratt = "3"
  ), ".pdf"), height = 12, width = 12)
  par(mfrow = c(1, 2))


  # plot elevation and scale bar
  plot(plotElevation,
    type = "phylogram", legend = FALSE,
    fsize = c(0.2, 0.9), outline = F, lwd = 1.4,
    col.text = "white"
  )

  udajeStromu <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  figLeft <- par("fig")


  axisPhylo(line = -1.8, cex = .4)
  add.color.bar(.25 * max(nodeHeights(strom)),
    cols = rev(hcl.colors(50, "inferno")),
    lims = range(elevation),
    outline = FALSE,
    title = "",
    subtitle = "Elevation",
    fsize = 1.1,
    prompt = F,
    x = 0.03 * par("usr")[2], y = ifelse(ktore == "Ratt", .95, .05) * par()$usr[4]
  )

  # taxon labels prep
  druhy <- data.frame(ID = sub("_.+", "", strom$tip.label))
  druhy$druh <- vzorky$Taxon[match(druhy$ID, vzorky$ID)]
  pocty <- rle(druhy$druh)
  endPos <- cumsum(pocty$lengths)
  startPos <- endPos - pocty$lengths + 1L
  ySegLeft <- matrix(grconvertY(c(udajeStromu$yy[startPos], udajeStromu$yy[endPos]), from = "user", to = "ndc"), ncol = 2)
  yNdcLeft <- grconvertY(rowMeans(cbind(udajeStromu$yy[startPos], udajeStromu$yy[endPos])), from = "user", to = "ndc")

  # plot region and legend
  plot(suhrn,
    colors = farby <- setNames(hcl.colors(3, "inferno"), c("Lowlands", "Mt. Wilhelm", "Saruwaged Range")),
    fsize = 0.2,
    direction = "leftwards",
    show.tip.label = F
  )

  figRight <- par("fig")
  axisPhylo(line = -1.8, cex = .4)
  udajeStromu <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  ySegRight <- grconvertY(c(udajeStromu$yy[startPos], udajeStromu$yy[endPos]), from = "user", to = "ndc")

  add.simmap.legend(
    leg = sub("_", " ", colnames(ancRegion$lik.anc)), colors = farby[sort(unique(region))], prompt = FALSE,
    x = .6 * par()$usr[2],
    y = ifelse(ktore == "Ratt", .92, .05) * par()$usr[4],
    fsize = 1.1,
    border = NA,
    title = "Region",
    adj = 1
  )

  # taxon labels
  xCenterNdc <- (figLeft[2] + figRight[1]) / 2


  par(fig = c(0, 1, 0, 1), new = TRUE)
  plot.new()
  par(usr = c(0, 1, 0, 1))

  # overlay ID labels (comment for debugging)
  rect(.424, 0.035, .5685, 1, col = "white", border = NA)

  # taxon lines
  segments(xCenterNdc - .07 + rep(c(0, 0.013), floor(nrow(ySegLeft) / 2)),
    y0 = ySegLeft[, 1], y1 = ySegLeft[, 2],
    lwd = 3
  )
  segments(xCenterNdc + .065 - rep(c(0, 0.013), floor(nrow(ySegLeft) / 2)),
    y0 = ySegLeft[, 1], y1 = ySegLeft[, 2],
    lwd = 3
  )

  # taxon text
  text(
    x = xCenterNdc, y = yNdcLeft, labels = pocty$values,
    font = 3, cex = 0.5, adj = 0.5, xpd = NA
  )
  dev.off()


  if (ouwie) {
    ######## Fenogram
    pdf(paste0("Results/Figure", switch(ktore,
      Hydro = "S3",
      Ratt = "S4"
    ), ".pdf"), height = 12, width = 10)

    # calculate and plot phenogram
    labelX <- phenogram(strom, elevation,
      ylim = c(0, 100 + max(vzorky$Elevation, na.rm = T)), xpd = NA,
      fsize = .4, xlab = "Time (Mya)", ylab = "Elevation (m asl)",
      las = 1, quiet = TRUE
    )
    udajeFenogramu <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    box()

    # draw one label per species
    # cover existing
    rect(
      xleft = nodeheight(strom, 1), ybottom = 0, xright = par("usr")[2], ytop = par("usr")[4],
      col = "white", border = NA
    )

    # calculate label heights
    fenoDruhy <- unique(data.frame(druh = druhy$druh, elevation = elevation))
    fenoStrom <- keep.tip(strom, row.names(fenoDruhy))
    fenoDruhy <- fenoDruhy[fenoStrom$tip.label, ]
    labelY <- phytools:::spreadlabels(fenoStrom,
      x = setNames(fenoDruhy$elevation, fenoDruhy$druh),
      fsize = .6, cost = c(5, 1), range = range(fenoDruhy$elevation)
    )
    text(
      x = rep(labelX[1, "x"], nrow(fenoDruhy)), y = labelY,
      labels = fenoDruhy$druh,
      cex = .6, font = 3, adj = 0
    )
    segments(nodeheight(strom, 1), fenoDruhy$elevation,
      labelX[1, "x"], labelY,
      lty = 2
    )

    # plot confidence intervals
    vyskaUzlov <- udajeFenogramu$xx[as.character((Ntip(strom) + 1):(Ntip(strom) + Nnode(strom)))]
    segments(
      x0 = vyskaUzlov,
      y0 = ancElevation$CI95[, 1],
      x1 = vyskaUzlov,
      y1 = ancElevation$CI95[, 2],
      col = hcl.colors(3, "inferno")[2]
    )
    dev.off()


    ##### OUwie

    # prepare node ancestral states and traits dataframe
    strom2 <- strom
    strom2$node.label <- unname(apply(ancRegion$lik.anc, 1, \(x) colnames(ancRegion$lik.anc)[which.max(x)]))
    traitDF <- data.frame(vzorka = names(elevation), lokalita = region, elevation = elevation / 1000)


    # region evolution along gradients
    fit <- lapply(c("BM1", "BMS", "OU1", "OUM", "OUMVA", "OUMA"), \(x) OUwie(strom2, traitDF, model = x, simmap.tree = FALSE))
    names(fit) <- c("BM1", "BMS", "OU1", "OUM", "OUMVA", "OUMA")
    summarizeOUwie <- function(x) {
      c(
        model = x$model,
        lnL = round(as.numeric(x$loglik), 1),
        AIC = round(as.numeric(x$AIC), 1),
        AICc = round(as.numeric(x$AICc), 1),
        BIC = round(as.numeric(x$BIC), 1),
        alpha = paste(round(as.numeric(x$solution[1, ]), 3), collapse = ";"),
        sigmaSq = paste(round(as.numeric(x$solution[2, ]), 3), collapse = ";"),
        theta = paste(round(as.numeric(x$theta[1:ifelse(x$model %in% c("BM1", "BMS"), 1, nlevels(x$tot.states)), 1]), 3), collapse = ";")
      )
    }
    suhrnOUwie <- do.call(rbind, lapply(fit, summarizeOUwie))
    suhrnOUwie <- suhrnOUwie[order(as.numeric(suhrnOUwie[, "AICc"])), ]

    write.table(suhrnOUwie,
      file = paste0("Results/OUwie-", switch(ktore,
        Hydro = "Hydromyini",
        Ratt = "Rattini"
      ), ".txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    capture.output(fit, file = paste0("Results/OUwie-", switch(ktore,
      Hydro = "Hydromyini",
      Ratt = "Rattini"
    ), ".txt"), append = TRUE)
  }
}
