## Headless APSIM runner — cloud / Linux compatible.
## Reads existing .met/.rds from outputs/, builds sim files, and runs APSIM.
## Outputs: sim-<loc>.apsimx (updated), grid_simulation_results.rds in outputs/.

suppressPackageStartupMessages({
  library(apsimx)
  library(dplyr)
  library(readxl)
})

## ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR <- normalizePath(
  if (nzchar(Sys.getenv("RSCRIPT_DIR"))) Sys.getenv("RSCRIPT_DIR")
  else dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))
  ], mustWork = FALSE))
)

OUT_DIR      <- normalizePath(file.path(SCRIPT_DIR, "../outputs"))
TEMPLATE     <- normalizePath(file.path(SCRIPT_DIR, "../apsim_templates/_corn-5-15-2026.apsimx"))
LOCATIONS_XL <- file.path(SCRIPT_DIR, "locations.updated.xlsx")

apsimx_options(exe.path = "/usr/local/lib/apsim/2025.3.7681.0/bin/Models")
setwd(OUT_DIR)
cat("Working directory:", getwd(), "\n")
cat("APSIM executable: /usr/local/lib/apsim/2025.3.7681.0/bin/Models\n\n")

## ── KL/XF templates (20 layers; trimmed to actual layers inside loop) ─────────
KL_template <- c(0.08, 0.08, 0.08, 0.08, 0.07, 0.07, 0.07, 0.07,
                 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04,
                 0.03, 0.03, 0.02, 0.02)
XF_template <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)

try_msg <- function(e) trimws(as.character(e))

## ── Load locations ────────────────────────────────────────────────────────────
sim.grid <- read_excel(LOCATIONS_XL, sheet = "Sheet2")
cat("Locations:\n"); print(as.data.frame(sim.grid), digits = 6); cat("\n")

## ── Status table ──────────────────────────────────────────────────────────────
status <- data.frame(location = sim.grid$location,
                     apsimx   = NA_character_,
                     sim      = NA_character_,
                     stringsAsFactors = FALSE)

results <- vector("list", nrow(sim.grid))

for (i in seq_len(nrow(sim.grid))) {

  loc      <- sim.grid$location[i]
  met.file <- paste0(loc, ".met")
  rds.file <- paste0(loc, ".rds")
  sim.file <- paste0("sim-", loc, ".apsimx")

  cat(sprintf("── [%d/%d] %s\n", i, nrow(sim.grid), loc))

  ## ---- Guard: require weather and soil ----------------------------------------
  if (!file.exists(met.file)) {
    cat("  SKIP — missing", met.file, "\n\n"); status$apsimx[i] <- "SKIP (no .met)"; next
  }
  if (!file.exists(rds.file)) {
    cat("  SKIP — missing", rds.file, "\n\n"); status$apsimx[i] <- "SKIP (no .rds)"; next
  }

  ## ---- Build sim file (always rebuild to ensure fresh soil/weather paths) ------
  spp   <- readRDS(rds.file)
  soils <- spp[[1]]
  soils <- apsimx:::fix_apsimx_soil_profile(soils, verbose = FALSE)
  soils$initialwater <- initialwater_parms(
    Depth         = soils$soil$Depth,
    Thickness     = soils$soil$Thickness,
    InitialValues = soils$soil$DUL
  )
  soils$crops <- "Maize"

  nlayers <- nrow(soils$soil)
  if (nlayers > length(KL_template))
    stop(loc, ": soil has ", nlayers, " layers but KL template only has ", length(KL_template))
  KL <- KL_template[seq_len(nlayers)]
  XF <- XF_template[seq_len(nlayers)]

  build_ok <- tryCatch({
    file.copy(TEMPLATE, sim.file, overwrite = TRUE)

    edit_apsimx_replace_soil_profile(
      file = sim.file, src.dir = OUT_DIR, wrt.dir = OUT_DIR,
      soil.profile = soils, verbose = FALSE, overwrite = TRUE
    )
    edit_apsimx(file = sim.file, src.dir = OUT_DIR, wrt.dir = OUT_DIR,
                node = "Soil", soil.child = "Physical",
                parm = "KL", value = KL, verbose = FALSE, overwrite = TRUE)
    edit_apsimx(file = sim.file, src.dir = OUT_DIR, wrt.dir = OUT_DIR,
                node = "Soil", soil.child = "Physical",
                parm = "XF", value = XF, verbose = FALSE, overwrite = TRUE)
    edit_apsimx(file = sim.file, src.dir = OUT_DIR, wrt.dir = OUT_DIR,
                node = "Weather", value = met.file,
                verbose = FALSE, overwrite = TRUE)
    TRUE
  }, error = function(e) { warning("  Build failed: ", try_msg(e)); FALSE })

  if (!build_ok) { status$apsimx[i] <- "BUILD ERROR"; next }
  cat("  .apsimx built\n")
  status$apsimx[i] <- "OK"

  ## ---- Run simulation ----------------------------------------------------------
  sim <- try(apsimx(file = sim.file, src.dir = OUT_DIR, cleanup = FALSE), silent = TRUE)

  if (inherits(sim, "try-error")) {
    cat("  APSIM run FAILED:", try_msg(sim), "\n\n")
    status$sim[i] <- paste0("ERROR: ", try_msg(sim))
    next
  }

  sim$location <- loc
  results[[i]] <- sim
  status$sim[i] <- paste0("OK (", nrow(sim), " rows)")
  cat("  simulation OK —", nrow(sim), "rows\n\n")
}

## ── Save combined results ─────────────────────────────────────────────────────
results_all <- dplyr::bind_rows(Filter(Negate(is.null), results))
saveRDS(results_all, file.path(OUT_DIR, "grid_simulation_results.rds"))
cat("Results saved to", file.path(OUT_DIR, "grid_simulation_results.rds"), "\n")

## ── Print summary ─────────────────────────────────────────────────────────────
cat("\n══ Simulation summary ══════════════════════════════════\n")
print(status)
write.csv(status, file.path(OUT_DIR, "simulation_status.csv"), row.names = FALSE)
cat("Status saved to simulation_status.csv\n")
