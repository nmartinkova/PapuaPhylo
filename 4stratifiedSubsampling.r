for (i in c("phytools", "ape", "OUwie")) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}

if (!dir.exists("Results")) dir.create("Results")

# best model from 3ancestralStates
bestModel <- switch(ktore,
  Hydro = "ARD",
  Ratt = "ER"
)

# input data
subory <- dir("Data", full.names = T)

vzorky <- read.table("Data/TableS1.txt",
  header = T, sep = "\t",
  quote = "\"", comment.char = ""
)

if (TRUE) {
  for (ktore in c("Hydro", "Ratt")) {
    # chronogram
    strom <- read.tree(
      file = paste0("Results/PNGchronogram", switch(ktore,
        Ratt = "Rattini",
        Hydro = "Hydromyini"
      ), ".tre")
    )

    # trait data
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
    regiony <- sort(unique(region))
    najmensiRegion <- min(table(region))

    # table for storing results
    res <- matrix(NA, ncol = 1 + factorial(length(regiony)) + length(regiony) + 1, nrow = 100)
    resOU <- matrix(NA, ncol = 8 * 6, nrow = 100)

    # stratified subsampling
    for (i in seq_len(100)) {
      cat(i, "\n")
      # sample balanced subsamples
      ktoreTips <- NULL
      for (r in regiony) {
        ktoreTips <- c(ktoreTips, sample(which(region == r), size = najmensiRegion, replace = FALSE))
      }

      stromTips <- keep.tip(strom, names(ktoreTips))
      elevationTips <- elevation[ktoreTips]
      regionTips <- region[ktoreTips]


      ####### Region uncertainty

      simulacie <- phytools::make.simmap(stromTips, regionTips, nsim = 100, model = bestModel, message = FALSE)
      suhrn <- phytools::describe.simmap(simulacie)
      res[i, ] <- c(round(colMeans(suhrn$count), 2), round(colMeans(suhrn$times), 2))


      strom2 <- stromTips
      strom2$node.label <- unname(apply(suhrn$ace[1:Nnode(strom2), ], 1, \(x) colnames(suhrn$ace)[which.max(x)]))
      traitDF <- data.frame(vzorka = names(elevation[ktoreTips]), lokalita = region[ktoreTips], elevation = elevation[ktoreTips] / 1000)

      # region evolution along gradients
      fit <- lapply(c("BM1", "BMS", "OU1", "OUM", "OUMVA", "OUMA"), \(x) OUwie(strom2, traitDF,
        model = x, simmap.tree = FALSE,
        algorithm = "invert",
        quiet = T
      ))
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
      suhrnOUwie <- do.call(cbind, lapply(fit, summarizeOUwie))
      resOU[i, ] <- suhrnOUwie


      colnames(res) <- c(colnames(suhrn$count), colnames(suhrn$times))
      write.table(res,
        sep = "\t", row.names = F,
        file = paste0("Results/sensitivitySimmap-", switch(ktore,
          Hydro = "Hydromyini",
          Ratt = "Rattini"
        ), ".txt")
      )

      colnames(resOU) <- rep(rownames(suhrnOUwie), 6)

      write.table(resOU,
        file = paste0("Results/sensitivityOUwie-", switch(ktore,
          Hydro = "Hydromyini",
          Ratt = "Rattini"
        ), ".txt"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
}


## summarise simmap

suhrn <- NULL
hlavicka <- unname(unlist(read.table("Results/sensitivitySimmap-Hydromyini.txt", nrows = 1, sep = "\t")))
hlavicka <- sub(",", ".", hlavicka)


for (ktore in c("Hydro", "Ratt")) {
  dat <- read.table(
    paste0(
      "Results/sensitivitySimmap-",
      switch(ktore,
        Hydro = "Hydromyini",
        Ratt = "Rattini"
      ), ".txt"
    ),
    header = T, sep = "\t"
  )
  if (is.null(suhrn)) {
    # Hydromyini
    suhrn <- paste(apply(dat, 2, \(x) round(mean(x, na.rm = T), 2)),
      apply(dat, 2, \(x) round(sd(x, na.rm = T), 2)),
      sep = "\\pm"
    )
  } else {
    # Rattini
    x <- paste(apply(dat, 2, \(x) round(mean(x, na.rm = T), 2)),
      apply(dat, 2, \(x) round(sd(x, na.rm = T), 2)),
      sep = "\\pm"
    )
    suhrn <- matrix(c(suhrn, rep(NA, length(suhrn))), nrow = 2, byrow = T, dimnames = list(c("Hydromyini", "Rattini"), hlavicka))
    suhrn[2, match(colnames(dat), hlavicka)] <- x
  }
}
write.table(suhrn, file = "Results/sensitivitySimmap.txt", sep = " & ", row.names = T)


## summarise OUwie
for (ktore in c("Hydro", "Ratt")) {
  dat <- read.table(
    paste0(
      "Results/sensitivityOUwie-",
      switch(ktore,
        Hydro = "Hydromyini",
        Ratt = "Rattini"
      ), ".txt"
    ),
    header = T, sep = "\t"
  )

  subor <- paste0("Results/sensitivitySummaryOUwie-", switch(ktore,
    Hydro = "Hydromyini",
    Ratt = "Rattini"
  ), ".txt")

  tab <- table(unname(unlist(dat[1, seq(1, ncol(dat), by = 8)])[apply(dat[, seq(4, ncol(dat), by = 8)], 1, which.min)]))
  najModel <- which(dat[1, ] == "BMS")
  najDat <- dat[, c(najModel:(najModel + 7))]


  splitToMat <- function(x, meno) {
    if (grepl(";", x[1])) {
      res <- t(sapply(strsplit(as.character(x), ";", fixed = TRUE), as.numeric))
    } else {
      res <- data.frame(x)
    }
    colnames(res) <- paste0(meno, 1:ncol(res))
    return(res)
  }
  sigmaSqMat <- splitToMat(najDat[[grep("^sigmaSq", names(najDat))]], "sigmaSq")
  alphaMat <- splitToMat(najDat[[grep("^alpha", names(najDat))]], "alpha")
  thetaMat <- splitToMat(najDat[[grep("^theta", names(najDat))]], "theta")

  najDat <- cbind(najDat[, 2:5], alphaMat, sigmaSqMat, thetaMat)

  cat("Best model: ", names(tab)[which.max(tab)], "\n\nTable of all best models across 100 subsamples:\n", file = subor, append = F)
  capture.output(tab, file = subor, append = T)
  cat("\n\nBMS model summary:\n", file = subor, append = T)
  cat(colnames(najDat[, 5:ncol(najDat)]), sep = " & ", file = subor, append = T)
  cat("\n", file = subor, append = T)
  cat(paste(apply(najDat[, 5:ncol(najDat)], 2, \(x) round(mean(x, na.rm = T), 2)),
    apply(najDat[, 5:ncol(najDat)], 2, \(x) round(sd(x, na.rm = T), 2)),
    sep = "\\pm"
  ), sep = " & ", file = subor, append = T, quote = F)
}
