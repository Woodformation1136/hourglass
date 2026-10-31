# -----
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
# -----

get_n_workers <- function(n_tasks, requested) {
  slurm_cpus <- suppressWarnings(
    as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_))
  )
  
  available <- if (is.na(slurm_cpus)) {
    max(1L, parallel::detectCores(logical = FALSE) - 1L)
  } else {
    slurm_cpus
  }
  
  max(
    1L, 
    min(as.integer(requested), as.integer(n_tasks), available)
  )
}

make_parallel_progress <- function(n) {
  pb <- utils::txtProgressBar(
    min = 0,
    max = n,
    style = 3,
    file = stderr()
  )
  
  list(
    pb = pb,
    snow_options = list(
      progress = function(completed) {
        utils::setTxtProgressBar(
          pb,
          completed
        )
      },
      
      # Required for progress to update as individual tasks finish.
      preschedule = FALSE
    )
  )
}