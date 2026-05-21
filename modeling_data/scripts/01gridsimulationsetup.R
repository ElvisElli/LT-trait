## Grid simulation setup for cropland in Arkansas
## Fetches weather (IEM + NASA POWER) and SSURGO soil profiles,
## then runs APSIM for each grid cell defined in locations.updated.xlsx.
rm(list = ls())

library(apsimx)
library(ggplot2)
library(stars)
library(sf)
library(dplyr)
library(rstudioapi)
library(readxl)

## ── Paths ────────────────────────────────────────────────────────────────────
apsimx_options(exe.path = "C:\\Users\\efelli\\AppData\\Local\\Programs\\APSIM2025.3.7681.0\\bin\\Models.exe")
setwd(dirname(getActiveDocumentContext()$path))

TEMPLATE_FILE <- "_corn-5-15-2026.apsimx"
DATE_RANGE    <- c("1985-01-01", "2025-12-31")

sim.grid <- read_excel("locations.updated.xlsx", sheet = "Sheet2")
print(as.data.frame(sim.grid), digits = 6)

## ── Constants (fixed across all grid cells) ──────────────────────────────────
KL <- c(0.08, 0.08, 0.08, 0.08, 0.07, 0.07, 0.07, 0.07,
        0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04,
        0.03, 0.03, 0.02, 0.02)

XF <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)

## ── Helper: IEM weather download with correct day-of-month format ─────────────
get_iem_fixed <- function(lonlat, dates, wrt.dir = ".", filename = "noname.met") {

  if (!grepl("met$", filename)) stop("filename should end in .met")

  d1 <- as.Date(dates[1]);  d2 <- as.Date(dates[2])
  yr1   <- format(d1, "%Y");  yr2   <- format(d2, "%Y")
  mnth1 <- format(d1, "%m");  mnth2 <- format(d2, "%m")
  dy1   <- as.integer(format(d1, "%d"))   # day-of-month, not %j
  dy2   <- as.integer(format(d2, "%d"))

  lon <- as.numeric(lonlat[1]);  lat <- as.numeric(lonlat[2])
  pts    <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = 1:2, crs = 4326)
  states <- sf::st_transform(spData::us_states, crs = 3857)
  pts    <- sf::st_transform(pts, crs = 3857)
  ii     <- as.integer(sf::st_intersects(pts, states))
  stt0   <- states[["NAME"]][[ii]]
  state  <- datasets::state.abb[which(stt0 == datasets::state.name)]
  stt.climate <- paste0(state, "CLIMATE")

  ftrs <- jsonlite::fromJSON(paste0(
    "http://mesonet.agron.iastate.edu/geojson/network.php?network=", stt.climate
  ))$features

  start.year <- sapply(ftrs$properties$time_domain, function(x)
    as.numeric(gsub("(", "", strsplit(x, "-")[[1]][1], fixed = TRUE)))
  end.year0  <- sapply(ftrs$properties$time_domain, function(x)
    gsub("(", "", strsplit(x, "-")[[1]][2], fixed = TRUE))
  end.year   <- ifelse(end.year0 == "Now)",
                       as.numeric(format(Sys.Date(), "%Y")),
                       suppressWarnings(as.numeric(gsub(")", "", end.year0, fixed = TRUE))))

  wch.yr <- which(start.year <= as.numeric(yr1) & end.year >= as.numeric(yr2))
  ftrs   <- ftrs[wch.yr, ]

  geo.lon <- sapply(ftrs$geometry$coordinates, `[[`, 1)
  geo.lat <- sapply(ftrs$geometry$coordinates, `[[`, 2)
  station.coords <- sf::st_as_sf(
    data.frame(lon = geo.lon, lat = geo.lat, station = ftrs$id),
    coords = c("lon", "lat"), crs = 4326
  )
  station.coords <- sf::st_transform(station.coords, crs = 3857)
  near.idx <- sf::st_nearest_feature(pts, station.coords)
  station  <- ftrs[near.idx, "id"]

  url <- paste0(
    "http://mesonet.agron.iastate.edu/cgi-bin/request/coop.py?network=", stt.climate,
    "&stations=", station,
    "&year1=", yr1,  "&month1=", mnth1, "&day1=", dy1,
    "&year2=", yr2, "&month2=", mnth2, "&day2=", dy2,
    "&vars%5B%5D=apsim&what=view&delim=comma&gis=no"
  )

  iem0 <- readLines(url)
  write(iem0, file = file.path(wrt.dir, filename))
  iem.dat <- read_apsim_met(filename, src.dir = wrt.dir, verbose = FALSE)

  if (filename == "noname.met") {
    unlink(file.path(wrt.dir, filename))
    return(iem.dat)
  } else {
    return(invisible(iem.dat))
  }
}

## ── Main loop ────────────────────────────────────────────────────────────────
results <- vector("list", nrow(sim.grid))

for (i in seq_len(nrow(sim.grid))) {

  loc     <- sim.grid$location[i]
  lonlat  <- unlist(sim.grid[i, c("x", "y")])
  cat("Processing cell", i, "/", nrow(sim.grid), "—", loc, "\n")

  ## -- Weather ----------------------------------------------------------------
  met.file <- paste0(loc, ".met")

  iem <- try(get_iem_fixed(lonlat = lonlat, dates = DATE_RANGE,
                           wrt.dir = ".", filename = met.file), silent = TRUE)
  if (inherits(iem, "try-error")) {
    warning("IEM download failed for ", loc, ": ", conditionMessage(attr(iem, "condition")))
    next
  }

  pwr <- try(get_power_apsim_met(lonlat = lonlat, dates = DATE_RANGE), silent = TRUE)
  if (inherits(pwr, "try-error")) {
    warning("NASA POWER download failed for ", loc, ": ", conditionMessage(attr(pwr, "condition")))
    next
  }

  iem$radn <- pwr$radn
  apsimx::write_apsim_met(iem, wrt.dir = ".", filename = met.file)

  ## -- Soil -------------------------------------------------------------------
  spp <- try(get_ssurgo_soil_profile(lonlat = lonlat, nlayers = 20, soil.bottom = 400),
             silent = TRUE)

  if (inherits(spp, "try-error")) {
    warning("SSURGO failed for ", loc, ": ", conditionMessage(attr(spp, "condition")))
    next
  }

  saveRDS(spp, paste0(loc, ".rds"))   # no extra list() wrapper

  soils <- spp[[1]]
  soils <- apsimx:::fix_apsimx_soil_profile(soils, verbose = FALSE)
  soils$initialwater <- initialwater_parms(
    Depth        = soils$soil$Depth,
    Thickness    = soils$soil$Thickness,
    InitialValues = soils$soil$DUL
  )
  soils$crops <- "Maize"

  ## -- Build per-cell APSIM file ----------------------------------------------
  par.sim.file <- paste0("par-sim-", loc, ".apsimx")

  # Copy directly from template — no intermediate file needed
  file.copy(TEMPLATE_FILE, par.sim.file, overwrite = TRUE)

  edit_apsimx_replace_soil_profile(
    file         = par.sim.file,
    src.dir      = ".",
    wrt.dir      = ".",
    soil.profile = soils,
    verbose      = FALSE,
    overwrite    = TRUE
  )

  edit_apsimx(file = par.sim.file, src.dir = ".", wrt.dir = ".",
              node = "Soil", soil.child = "Physical",
              parm = "KL", value = KL,
              verbose = FALSE, overwrite = TRUE)

  edit_apsimx(file = par.sim.file, src.dir = ".", wrt.dir = ".",
              node = "Soil", soil.child = "Physical",
              parm = "XF", value = XF,
              verbose = FALSE, overwrite = TRUE)

  edit_apsimx(file = par.sim.file, src.dir = ".", wrt.dir = ".",
              node    = "Weather",
              value   = met.file,
              verbose = FALSE, overwrite = TRUE)

  ## -- Run simulation ---------------------------------------------------------
  sim <- try(apsimx(file = par.sim.file, src.dir = ".", cleanup = TRUE), silent = TRUE)

  # Remove per-cell APSIM file after run
  unlink(par.sim.file)

  if (inherits(sim, "try-error")) {
    warning("APSIM run failed for ", loc, ": ", conditionMessage(attr(sim, "condition")))
    next
  }

  sim$location <- loc
  results[[i]] <- sim
}

## ── Combine and save all results ─────────────────────────────────────────────
results_all <- dplyr::bind_rows(Filter(Negate(is.null), results))
saveRDS(results_all, "grid_simulation_results.rds")
cat("Done. Results saved to grid_simulation_results.rds\n")
