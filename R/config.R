#' Create LINGUINE Configuration
#'
#' @description Generates a centralized configuration object for the Ancestral Genome Reconstruction Pipeline.
#' Instead of relying on global variables, this object is passed to downstream functions.
#'
#' @param dataset Character. Name of the dataset.
#' @param min_chromosome_length_bp Numeric. Minimum chromosome/scaffold length filter (in base pairs).
#' @param base_dir Character. Root directory for inputs/outputs. Defaults to current working directory.
#' @param hmm_params List. Transition and emission probabilities.
#' @param thresholds List. Thresholds for classification, graph clustering, and plotting.
#'
#' @return A structured list of class `linguine_config` containing all pipeline parameters and paths.
#' @export
create_linguine_config <- function(
    dataset = "LINGUINE",
    min_chromosome_length_bp = 0,
    base_dir = getwd(),
    hmm_params = list(
      prob_self_syn_hmm1 = 0.995,
      prob_to_other_syn_hmm1 = 0.005,
      prob_emit_correct_syn_hmm1_value = 0.85,
      prob_emit_non_syn_obs_from_syn_hmm1_value = 0.08,
      prob_emit_no_ortholog_from_syn_hmm1_value = 0.05,
      prob_emit_multiple_b_chrs_from_syn_hmm1_value = 0.01
    ),
    thresholds = list(
      window_size = 10,
      purity_threshold = 0.80,
      min_final_block_size = 10,
      fisher_p = 0.01,
      fisher_odds = 2.0,
      min_ogs_by_block = 2,
      og_abundance = 90,
      fisher_p_lgs = 1e-20,
      fisher_odds_lgs = 4,
      paralogy_p = 1e-20,
      paralogy_odds = 10.0,
      outgroup_dominance = 0.05,
      parent_assignment = 0.60,
      ribbon_max_links = 5
    )
) {

  # 1. Generate dynamic run ID based on inputs
  run_id <- paste0(dataset, "_min_chr_size_", format(min_chromosome_length_bp, scientific = FALSE), "bp")

  # 2. Build the configuration object
  config <- list(
    dataset = dataset,
    run_id = run_id,
    min_chromosome_length_bp = min_chromosome_length_bp,
    orthology_type = "OGs",
    orthology_filename = "Orthogroups.tsv",
    paths = list(
      raw_data = file.path(base_dir, "raw_data", dataset),
      processed_data = file.path(base_dir, "runs", run_id, "processed_data"),
      intermediate_data = file.path(base_dir, "runs", run_id, "intermediate_data"),
      results = file.path(base_dir, "runs", run_id, "results"),
      plots = file.path(base_dir, "runs", run_id, "plots")
    ),
    hmm = hmm_params,
    thresholds = thresholds
  )

  # 3. Assign an S3 class so we can easily validate this object later if needed
  class(config) <- "linguine_config"

  return(config)
}
