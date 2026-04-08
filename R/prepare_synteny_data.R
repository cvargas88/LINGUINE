#' Prepare Synteny Data and HMM Emission Sequences
#'
#' @description Ingests global orthology data, filters for pairwise comparisons,
#' and maps spatial gene coordinates to discrete HMM emission symbols
#' (e.g., "ON_Chr1", "NO_ORTHOLOG", "MULTIPLE_B_CHRS").
#'
#' @param ref_species Character. Name of the reference lineage.
#' @param comp_species Character. Name of the comparison lineage.
#' @param comparison_type Character. Architecture of the comparison ("tip_vs_tip", "tip_vs_INode", "INode_vs_INode").
#' @param parent_node Character. Name of the parent node being reconstructed.
#' @param comp_daughters Character vector. Names of the comparison node's daughters (required for tip_vs_INode).
#' @param config A `linguine_config` object.
#'
#' @return Character path to the saved HMM observations RDS file, or NULL if skipped.
#' @export
prepare_synteny_data <- function(ref_species, comp_species, comparison_type, parent_node, comp_daughters = NULL, config) {

  if (comparison_type == "INode_vs_INode") {
    message("Skipping preliminary HMM emission mapping for structural Node-vs-Node analysis.")
    return(NULL)
  }

  message(sprintf("\n--- Integrating Genomic and Orthology Data (%s vs %s) ---", ref_species, comp_species))

  # ----------------------------------------------------------------------------
  # 1. Global Orthology Ingestion (Run once per dataset)
  # ----------------------------------------------------------------------------
  ortho_type <- config$orthology_type
  ortho_path <- file.path(config$paths$processed_data, paste0("ortholog_data_", ortho_type, ".rds"))

  if (file.exists(ortho_path)) {
    orthogroups_final <- readRDS(ortho_path)
  } else {
    message("Parsing global orthology data for the first time...")
    # NOTE: To keep this file clean, we rely on a helper function in utils.R
    # to handle the massive HOGs pivot_longer logic.
    orthogroups_final <- parse_global_orthology(config)
    saveRDS(orthogroups_final, ortho_path)
  }

  # ----------------------------------------------------------------------------
  # 2. Pairwise Filtration & Spatial Join
  # ----------------------------------------------------------------------------
  filtered_ortho_path <- file.path(config$paths$processed_data, paste0(ref_species, "_vs_", comp_species, "_ortholog_data_filtered.rds"))

  if (!file.exists(filtered_ortho_path)) {
    id_col_name <- if (ortho_type == "OGs") "Orthogroup" else paste0("HOG_", parent_node)
    if (!id_col_name %in% colnames(orthogroups_final)) id_col_name <- "Orthogroup" # Fallback for N0

    # Isolate shared orthogroups
    filtered_df <- orthogroups_final |>
      dplyr::select(Gene_ID, Species, OG_ID = dplyr::sym(id_col_name)) |>
      dplyr::filter(!is.na(OG_ID)) |>
      dplyr::group_by(OG_ID) |>
      dplyr::filter(dplyr::n_distinct(Species) > 1) |>
      dplyr::ungroup()

    # Map Reference Coordinates
    ref_genes <- readRDS(file.path(config$paths$processed_data, paste0(ref_species, "_genes_df.rds")))

    ref_df <- filtered_df |>
      dplyr::filter(Species == ref_species) |>
      dplyr::select(ref_gene_id = Gene_ID, OG_ID)

    if (comparison_type == "tip_vs_tip") {
      comp_genes <- readRDS(file.path(config$paths$processed_data, paste0(comp_species, "_genes_df.rds")))
      comp_df <- filtered_df |>
        dplyr::filter(Species == comp_species) |>
        dplyr::select(comp_gene_id = Gene_ID, OG_ID)

      ortholog_data_filtered <- ref_df |>
        dplyr::left_join(comp_df, by = "OG_ID", relationship = "many-to-many") |>
        dplyr::inner_join(dplyr::select(ref_genes, ref_gene_id = gene_id, ref_chromosome = chromosome, ref_start = start, ref_end = end), by = "ref_gene_id") |>
        dplyr::inner_join(dplyr::select(comp_genes, comp_gene_id = gene_id, comp_chromosome = chromosome, comp_start = start, comp_end = end), by = "comp_gene_id")

    } else { # tip_vs_INode
      ortholog_data_filtered <- ref_df |>
        dplyr::inner_join(dplyr::select(ref_genes, ref_gene_id = gene_id, ref_chromosome = chromosome, ref_start = start, ref_end = end), by = "ref_gene_id")
    }

    if (nrow(ortholog_data_filtered) == 0) stop("CRITICAL ERROR: No viable ortholog pairs remained post-spatial mapping.")
    saveRDS(ortholog_data_filtered, filtered_ortho_path)
  } else {
    ortholog_data_filtered <- readRDS(filtered_ortho_path)
  }

  # ----------------------------------------------------------------------------
  # 3. HMM Emission Sequence Generation
  # ----------------------------------------------------------------------------
  hmm_obs_path <- file.path(config$paths$intermediate_data, paste0("hmm_observations_detailed_", ref_species, "_ref_", comp_species, "_comp.rds"))

  if (file.exists(hmm_obs_path)) {
    message("Info: HMM emission sequence already exists. Skipping generation.")
    return(hmm_obs_path)
  }

  message("Computing Structural HMM Emission Sequences...")
  ref_genes_df <- readRDS(file.path(config$paths$processed_data, paste0(ref_species, "_genes_df.rds")))

  if (comparison_type == "tip_vs_tip") {
    # Resolve multiple-mapping conflicts via majority consensus
    gene_comp_chroms <- ortholog_data_filtered |>
      dplyr::group_by(ref_gene_id) |>
      dplyr::count(comp_chromosome, name = "count_chr") |>
      dplyr::filter(count_chr == max(count_chr)) |>
      dplyr::summarise(primary_comp_chr = sample(comp_chromosome, 1), .groups = 'drop')

    gene_ortho_summary <- ortholog_data_filtered |>
      dplyr::group_by(ref_gene_id) |>
      dplyr::summarise(num_ortholog_hits = dplyr::n(), all_comp_chromosomes = list(unique(comp_chromosome)), Orthogroup = unique(OG_ID)[1], .groups = 'drop')

    hmm_observation_df <- ref_genes_df |>
      dplyr::left_join(gene_comp_chroms, by = c("gene_id" = "ref_gene_id")) |>
      dplyr::left_join(gene_ortho_summary, by = c("gene_id" = "ref_gene_id")) |>
      dplyr::mutate(
        num_ortholog_hits = dplyr::coalesce(num_ortholog_hits, 0L),
        all_comp_chromosomes = purrr::map(all_comp_chromosomes, function(x) if (is.null(x)) character(0) else x)
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        temp_observation = dplyr::case_when(
          num_ortholog_hits == 0 ~ "NO_ORTHOLOG",
          length(all_comp_chromosomes) > 1 ~ "MULTIPLE_B_CHRS",
          length(all_comp_chromosomes) == 1 ~ paste0("ON_", primary_comp_chr),
          TRUE ~ "UNKNOWN_CLASSIFICATION"
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::select(gene_id, chromosome, start, end, primary_comp_chr, temp_observation, all_comp_chromosomes, Orthogroup)

  } else { # tip_vs_INode
    # Find the correct ancestral map (either A_B or B_A)
    anc_path <- file.path(config$paths$results, paste0("ancestral_genome_", comp_daughters[1], "_", comp_daughters[2], ".rds"))
    if (!file.exists(anc_path)) anc_path <- file.path(config$paths$results, paste0("ancestral_genome_", comp_daughters[2], "_", comp_daughters[1], ".rds"))
    if (!file.exists(anc_path)) stop("Ancestral genome mapping missing for structural evaluation.")

    anc_df <- readRDS(anc_path)

    # HOG mapping logic wrapped safely
    if (ortho_type == "HOGs") {
      hog_col <- paste0("HOG_", parent_node)
      anc_df <- anc_df |>
        dplyr::mutate(
          Ancestor_Full_OGs = purrr::map(Ancestor_Full_Genes, map_genes_to_hogs, mapping_df = dplyr::select(orthogroups_final, Gene_ID, Target_HOG = dplyr::sym(hog_col)))
        )
    }

    inode_og_lookup <- anc_df |>
      dplyr::select(ancestral_lg_name, Ancestor_Full_OGs) |>
      tidyr::unnest(Ancestor_Full_OGs) |>
      dplyr::rename(Orthogroup = Ancestor_Full_OGs) |>
      dplyr::distinct() |>
      dplyr::group_by(Orthogroup) |>
      dplyr::summarise(relevant_lgs = list(unique(ancestral_lg_name)), .groups = 'drop')

    hmm_observation_df <- ref_genes_df |>
      dplyr::left_join(dplyr::distinct(dplyr::select(ortholog_data_filtered, gene_id = ref_gene_id, Orthogroup = OG_ID)), by = "gene_id") |>
      dplyr::left_join(inode_og_lookup, by = "Orthogroup") |>
      dplyr::mutate(
        all_comp_chromosomes = dplyr::case_when(
          is.na(Orthogroup) | purrr::map_lgl(relevant_lgs, is.null) ~ list(character(0)),
          TRUE ~ relevant_lgs
        ),
        lg_count = purrr::map_dbl(all_comp_chromosomes, length),
        temp_observation = dplyr::case_when(
          lg_count == 0 ~ "NO_ORTHOLOG",
          lg_count > 1 ~ "MULTIPLE_B_CHRS",
          lg_count == 1 ~ purrr::map_chr(all_comp_chromosomes, ~ paste0("ON_", .x[1])),
          TRUE ~ "UNKNOWN_CLASSIFICATION"
        ),
        primary_comp_chr = ifelse(grepl("^ON_", temp_observation), sub("ON_", "", temp_observation), NA_character_)
      ) |>
      dplyr::select(gene_id, chromosome, start, end, primary_comp_chr, temp_observation, all_comp_chromosomes, Orthogroup)
  }

  saveRDS(hmm_observation_df, hmm_obs_path)
  message("HMM sequence assignment successfully completed.")
  return(hmm_obs_path)
}
