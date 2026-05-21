## Standalone script: generate .met and .rds for all locations in
## locations.updated.xlsx (Sheet2).  Run from the repo root or from
## the modeling_data/scripts/ folder.
## Outputs go to the same folder as this script.

suppressPackageStartupMessages({
  library(apsimx)
  library(sf)
  library(dplyr)
  library(readxl)
  library(jsonlite)
  library(spData)
})

## ── Working directory: always the folder that contains this script ─────────────
script_dir <- normalizePath(
  if (nzchar(Sys.getenv("RSCRIPT_DIR"))) Sys.getenv("RSCRIPT_DIR")
  else dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))
  ], mustWork = FALSE))
)
setwd(script_dir)
cat("Working directory:", getwd(), "\n\n")

## ── Settings ──────────────────────────────────────────────────────────────────
DATE_RANGE <- c("1985-01-01", "2025-12-31")
NLAYERS    <- 20
SOIL_BOT   <- 400   # cm

## ── Load locations ────────────────────────────────────────────────────────────
sim.grid <- read_excel("locations.updated.xlsx", sheet = "Sheet2")
cat("Locations to process:\n")
print(as.data.frame(sim.grid), digits = 6)
cat("\n")

## ── Helper: safe error message ────────────────────────────────────────────────
try_msg <- function(e) trimws(as.character(e))

## ── Helper: IEM download (day-of-month URL fix + scalar station id) ───────────
get_iem_fixed <- function(lonlat, dates, wrt.dir = ".", filename = "noname.met") {

  if (!grepl("met$", filename)) stop("filename should end in .met")

  d1 <- as.Date(dates[1]);  d2 <- as.Date(dates[2])
  yr1 <- format(d1, "%Y");  yr2 <- format(d2, "%Y")
  m1  <- format(d1, "%m");  m2  <- format(d2, "%m")
  dy1 <- as.integer(format(d1, "%d"))
  dy2 <- as.integer(format(d2, "%d"))

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
  end.year <- ifelse(end.year0 == "Now)",
                     as.numeric(format(Sys.Date(), "%Y")),
                     suppressWarnings(as.numeric(gsub(")", "", end.year0, fixed = TRUE))))

  wch <- which(start.year <= as.numeric(yr1) & end.year >= as.numeric(yr2))
  ftrs <- ftrs[wch, ]

  geo.lon <- sapply(ftrs$geometry$coordinates, `[[`, 1)
  geo.lat <- sapply(ftrs$geometry$coordinates, `[[`, 2)
  sta.sf  <- sf::st_as_sf(data.frame(lon = geo.lon, lat = geo.lat, station = ftrs$id),
                           coords = c("lon","lat"), crs = 4326)
  sta.sf  <- sf::st_transform(sta.sf, crs = 3857)
  station <- ftrs$id[sf::st_nearest_feature(pts, sta.sf)]   # scalar

  url <- paste0(
    "http://mesonet.agron.iastate.edu/cgi-bin/request/coop.py?network=", stt.climate,
    "&stations=", station,
    "&year1=", yr1, "&month1=", m1, "&day1=", dy1,
    "&year2=", yr2, "&month2=", m2, "&day2=", dy2,
    "&vars%5B%5D=apsim&what=view&delim=comma&gis=no"
  )

  iem0    <- readLines(url)
  write(iem0, file = file.path(wrt.dir, filename))
  iem.dat <- read_apsim_met(filename, src.dir = wrt.dir, verbose = FALSE)

  if (filename == "noname.met") {
    unlink(file.path(wrt.dir, filename))
    return(iem.dat)
  }
  invisible(iem.dat)
}

## ── Main loop ─────────────────────────────────────────────────────────────────
status <- data.frame(location = sim.grid$location,
                     met = NA_character_,
                     rds = NA_character_,
                     stringsAsFactors = FALSE)

for (i in seq_len(nrow(sim.grid))) {

  loc    <- sim.grid$location[i]
  lonlat <- unlist(sim.grid[i, c("x", "y")])
  cat(sprintf("── [%d/%d] %s  (%.4f, %.4f)\n",
              i, nrow(sim.grid), loc, lonlat[1], lonlat[2]))

  met.file <- paste0(loc, ".met")
  rds.file <- paste0(loc, ".rds")

  ## ---- Weather: IEM + NASA POWER radn ----------------------------------------
  if (file.exists(met.file)) {
    cat("  .met already exists — skipping download\n")
    status$met[i] <- "skipped (exists)"
  } else {
    iem <- try(get_iem_fixed(lonlat = lonlat, dates = DATE_RANGE,
                             wrt.dir = ".", filename = met.file), silent = TRUE)
    if (inherits(iem, "try-error")) {
      warning("  IEM failed: ", try_msg(iem))
      status$met[i] <- paste0("ERROR: ", try_msg(iem))
    } else {
      pwr <- try(get_power_apsim_met(lonlat = lonlat, dates = DATE_RANGE), silent = TRUE)
      if (inherits(pwr, "try-error")) {
        warning("  NASA POWER failed: ", try_msg(pwr))
        status$met[i] <- paste0("ERROR (radn): ", try_msg(pwr))
      } else {
        ## Date-matched radn merge — safe against row-count differences
        pwr_df  <- data.frame(year = pwr$year, day = pwr$day, radn = pwr$radn)
        iem_key <- paste(iem$year, iem$day)
        iem$radn <- pwr_df$radn[match(iem_key, paste(pwr_df$year, pwr_df$day))]
        if (anyNA(iem$radn)) {
          n_na <- sum(is.na(iem$radn))
          warning("  ", n_na, " days missing POWER radn — filling with 15 MJ/m2/day")
          iem$radn[is.na(iem$radn)] <- 15
        }
        apsimx::write_apsim_met(iem, wrt.dir = ".", filename = met.file)
        cat("  .met written:", met.file, "\n")
        status$met[i] <- "OK"
      }
    }
  }

  ## ---- Soil: SSURGO ----------------------------------------------------------
  if (file.exists(rds.file)) {
    cat("  .rds already exists — skipping download\n")
    status$rds[i] <- "skipped (exists)"
  } else {
    spp <- try(get_ssurgo_soil_profile(lonlat = lonlat,
                                       nlayers  = NLAYERS,
                                       soil.bottom = SOIL_BOT), silent = TRUE)
    if (inherits(spp, "try-error")) {
      warning("  SSURGO failed: ", try_msg(spp))
      status$rds[i] <- paste0("ERROR: ", try_msg(spp))
    } else {
      saveRDS(spp, rds.file)
      cat("  .rds written:", rds.file, "\n")
      status$rds[i] <- "OK"
    }
  }

  cat("\n")
}

## ── Summary ───────────────────────────────────────────────────────────────────
cat("══ Download summary ══════════════════════════════════\n")
print(status)
write.csv(status, "download_status.csv", row.names = FALSE)
cat("\nStatus saved to download_status.csv\n")
