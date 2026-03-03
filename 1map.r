for (i in c("terra", "geodata")) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE)
    library(i, character.only = TRUE)
  }
}

if (!dir.exists("Data")) dir.create("Data")
if (!dir.exists("Results")) dir.create("Results")

# Cizkova email 260123
sites <- data.frame(
  region = c(
    rep("Mt_Wilhelm", 8),
    rep("Saruwaged", 8),
    "Lowlands", "Lowlands", "Lowlands"
  ),
  elevation = c(
    3700, 3200, 2700, 2200, 1700, 1200, 700, 200,
    2200, 200, 700, 1200, 1700, 2700, 3200, 3700,
    200, 45, 10
  ),
  site = c(
    "MW3700", "MW3200", "MW2700", "MW2200", "MW1700", "MW1200", "MW700", "MW200",
    "FR2200", "FR200", "FR700", "FR1200", "FR1700", "FR2700", "FR3200", "FR3700",
    "WG", "BAI", "NAG"
  ),
  lat = c(
    -5.788567, -5.7998, -5.814567, -5.758883, -5.7658, -5.722217, -5.724333, -5.737167,
    -6.027355, -5.924077, -5.948569, -5.965452, -6.032105, -6.110915, -6.183584, -6.188394,
    -5.255425, -5.138183, -5.156342
  ),
  lon = c(
    145.060833, 145.069267, 145.1589, 145.186183, 145.201717, 145.270083, 145.259533, 145.327967,
    146.833739, 146.856214, 146.844851, 146.834488, 146.811055, 146.895523, 146.895794, 146.90671,
    145.147331, 145.771967, 145.7949
  ),
  label = c(
    "Mt. Wilhelm MW3700", NA, NA, NA, NA, NA, NA, "Mt. Wilhelm MW200",
    NA, "Saruwaged FR200", NA, NA, NA, NA, NA, "Saruwaged FR3700",
    "Wanang", "Baitabag", "Nagada"
  ),
  stringsAsFactors = FALSE
)
# Cizkova email 260123 - delete Nagada + reorder Finisterres by elevation
sites <- sites[c(1:8, 10:13, 9, 14:18), ]
rownames(sites) <- NULL

# download elevation data for the relief map
dem <- geodata::elevation_30s(country = "PNG", path = "Data")
dem2 <- geodata::elevation_3s(lon = 145.0608, lat = -5.7885, path = "Data")
dem <- rast("Data/elevation/PNG_elv_msk.tif") # 30s
dem2 <- rast("Data/elevation/srtm_66_14.tif") # 3s
gpx <- vect(as.matrix(sites[, c("lon", "lat")]), type = "points")

# extend map border
lims <- ext(gpx) + c(1.2, 1.2, 1, 1)

# terrain shading
slp <- terrain(dem, v = "slope", unit = "radians")
asp <- terrain(dem, v = "aspect", unit = "radians")
hs <- shade(slp, asp)

# plot a small map for Figure 4
pdf("Results/map.pdf", width = 6, height = 4.6)
plot(crop(hs, lims),
  col = gray.colors(256), legend = FALSE, axes = TRUE, las = 1,
  pax = list(retro = TRUE)
)
plot(gpx, add = T)
ktore <- c(1, 8, 9, 16:18)
text(sites[ktore, "lon"], sites[ktore, "lat"], labels = sites[ktore, "label"], pos = c(1, 3, 3, 1, 3, 3))
dev.off()


extendLine <- function(w) {
  wUtm <- terra::project(w, "EPSG:32755") # meters
  xy <- terra::crds(wUtm)
  p0 <- xy[1, ]
  p1 <- xy[nrow(xy), ]

  v <- p1 - p0
  u <- v / sqrt(sum(v * v))

  extendM <- 5000 # e.g. 2 km on each end
  p0e <- p0 - extendM * u
  p1e <- p1 + extendM * u

  wExt <- terra::vect(rbind(p0e, xy, p1e), type = "lines", crs = terra::crs(wUtm))
  wExt <- terra::project(wExt, terra::crs(dem))
  return(wExt)
}


pdf("Results/FigureS1.pdf", height = 2 * 3.7, width = 2 * 5)
par(mfrow = c(2, 2), mar = c(4, 4, 0, 0) + .3)

# New Guinea
plot(hs,
  col = gray.colors(256), legend = FALSE, axes = TRUE, las = 1, mar = NA,
  pax = list(retro = TRUE)
)
lines(x = c(lims[2], par("usr")[2]), y = c(lims[3], par("usr")[3]))
lines(x = c(lims[2], par("usr")[2]), y = c(lims[4], par("usr")[4]))
plot(lims, add = T)

# sites
plot(crop(hs, lims),
  col = gray.colors(256), legend = FALSE, axes = TRUE, las = 1, mar = NA,
  pax = list(retro = TRUE)
)
plot(gpx, add = T)
text(sites[!is.na(sites$label), "lon"], sites[!is.na(sites$label), "lat"],
  labels = sites[!is.na(sites$label), "label"], pos = c(1, 3, 3, 1, 1, 3, 3)
)

# transects
for (i in c("Mt_Wilhelm", "Saruwaged")) {
  w <- extendLine(vect(as.matrix(sites[sites$region == i, c("lon", "lat")]), type = "lines", crs = crs(dem2)))

  # extract cell elevation between points
  prof <- terra::extractAlong(dem2, w, xy = TRUE, online = TRUE, bilinear = FALSE)

  # distances between cells with elevation (1km)
  prof$distance <- c(0, cumsum(terra::distance(
    prof[-nrow(prof), 2:3], prof[-1, 2:3],
    lonlat = terra::is.lonlat(dem2), pairwise = TRUE
  ))) / 1000

  ## profile plot
  plot(prof$distance, prof[, 4],
    type = "n", axes = TRUE, las = 1,
    ylab = "Elevation m asl", xlab = sub("_", ". ", i), ylim = c(0, 4100)
  )

  # shaded under curve
  polygon(c(prof$distance, rev(prof$distance)),
    c(prof[, 4], rep(0, nrow(prof))),
    border = NA, col = gray(0.8)
  )

  # draw the profile line on top
  lines(prof$distance, prof[, 4], lwd = 2)

  # add site points and names
  sitePts <- vect(as.matrix(sites[sites$region == i, c("lon", "lat")]),
    type = "points", crs = crs(dem)
  )
  profPts <- vect(as.matrix(prof[, c("x", "y")]), type = "points", crs = crs(dem))
  nearest <- apply(terra::distance(sitePts, profPts), 1L, which.min)

  points(prof$distance[nearest], prof[, 4][nearest], pch = 19)
  text(prof$distance[nearest], prof[, 4][nearest],
    labels = sites[sites$region == i, "site"],
    pos = 3, cex = 0.8
  )
}
dev.off()
