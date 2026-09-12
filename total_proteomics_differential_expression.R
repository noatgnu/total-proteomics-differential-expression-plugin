library(SummarizedExperiment)
library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(svglite)
library(scales)
library(writexl)

annotate  <- ggplot2::annotate
filter    <- dplyr::filter
select    <- dplyr::select
rename    <- dplyr::rename
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
count     <- dplyr::count

detect_delimiter <- function(filepath) {
  ext <- tolower(tools::file_ext(filepath))
  if (ext == "csv") "," else "\t"
}

basename_of <- function(x) sub(".*[/\\\\]", "", x)

replace_special_with_dot <- function(str) gsub("[^[:alnum:]]+", ".", str)

# Sanitizes a group/condition label into a valid R name for makeContrasts.
safe_name <- function(str) {
  clean <- replace_special_with_dot(str)
  if (grepl("^[0-9]", clean)) clean <- paste0("X", clean)
  clean
}

save_all_formats <- function(p, path_no_ext, width = 9, height = 8) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height, dpi = 300)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300)
  ggsave(paste0(path_no_ext, ".svg"), p, width = width, height = height)
}

read_annotation <- function(annotation_file) {
  sep <- detect_delimiter(annotation_file)
  df <- read.table(annotation_file, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("Sample" %in% colnames(df)) || !("Condition" %in% colnames(df))) {
    stop("Annotation file must have Sample and Condition columns", call. = FALSE)
  }
  df$Sample <- trimws(df$Sample)
  df$Condition <- trimws(df$Condition)
  df
}

# Matches each annotation Sample against a matrix column: exact string first, basename fallback.
match_sample_columns <- function(matrix_colnames, annotation_df) {
  matched <- rep(NA_character_, nrow(annotation_df))
  for (i in seq_len(nrow(annotation_df))) {
    s <- annotation_df$Sample[i]
    if (s %in% matrix_colnames) {
      matched[i] <- s
      next
    }
    hit <- matrix_colnames[basename_of(matrix_colnames) == basename_of(s)]
    if (length(hit) >= 1) matched[i] <- hit[1]
  }
  matched
}

# Ported verbatim from the original analysis: pools each Condition group's own samples
# into one mean/sd, then downshifts. A group-aware MNAR variant, not Perseus's literal
# per-column default (see perseus_default_impute below).
group_downshift_impute <- function(mat, groups, downshift, width) {
  imputed <- mat
  for (g in levels(groups)) {
    cols <- which(groups == g)
    sub_vals <- mat[, cols, drop = FALSE]
    observed <- sub_vals[!is.na(sub_vals)]
    if (length(observed) == 0) next
    mu <- mean(observed); sdv <- sd(observed)
    na_idx <- which(is.na(sub_vals), arr.ind = TRUE)
    if (nrow(na_idx) > 0)
      sub_vals[na_idx] <- rnorm(nrow(na_idx), mean = mu - downshift * sdv, sd = width * sdv)
    imputed[, cols] <- sub_vals
  }
  imputed
}

# Genuine Perseus-default "Replace missing values from normal distribution", mode
# "Separately for each column": mean/sd computed from each sample's own observed values.
perseus_default_impute <- function(mat, downshift, width) {
  imputed <- mat
  for (j in seq_len(ncol(mat))) {
    col_vals <- mat[, j]
    observed <- col_vals[!is.na(col_vals)]
    if (length(observed) == 0) next
    mu <- mean(observed); sdv <- sd(observed)
    na_idx <- which(is.na(col_vals))
    if (length(na_idx) > 0)
      col_vals[na_idx] <- rnorm(length(na_idx), mean = mu - downshift * sdv, sd = width * sdv)
    imputed[, j] <- col_vals
  }
  imputed
}

run_differential_expression <- function(normalized_matrix_file, annotation_file, comparison_file,
                                         output_folder, impute_method = "none",
                                         impute_downshift = 1.8, impute_width = 0.3,
                                         fdr_cutoff = 0.05, fc_threshold_log2 = log2(1.5),
                                         target_genes = NULL, top_n_genes = 15) {

  for (sub in c("DE", "curtain")) dir.create(file.path(output_folder, sub), showWarnings = FALSE, recursive = TRUE)

  up_col <- "#E8604C"; dn_col <- "#7FB3D3"; ns_col <- "grey75"

  message("Loading normalised matrix...")
  mat_sep <- detect_delimiter(normalized_matrix_file)
  mat_df <- read.table(normalized_matrix_file, sep = mat_sep, header = TRUE,
                        na.strings = c("NA", "NaN", "N/A", "#VALUE!"),
                        check.names = FALSE, stringsAsFactors = FALSE)

  meta_cols <- c("Protein.Group", "Genes", "primary_gene", "Category")
  missing_meta <- setdiff(meta_cols, colnames(mat_df))
  if (length(missing_meta) > 0) {
    stop(paste0("Normalised matrix is missing expected column(s): ", paste(missing_meta, collapse = ", ")), call. = FALSE)
  }
  gene_map <- mat_df[, meta_cols]

  message("Loading annotation file...")
  annotation_df <- read_annotation(annotation_file)
  sample_cols <- match_sample_columns(colnames(mat_df), annotation_df)
  if (any(is.na(sample_cols))) {
    stop(paste0("Could not match annotation samples to matrix columns: ",
                paste(annotation_df$Sample[is.na(sample_cols)], collapse = ", ")), call. = FALSE)
  }

  group_levels <- unique(annotation_df$Condition)
  sample_annotation <- tibble::tibble(
    sample_name = sample_cols,
    group = factor(annotation_df$Condition, levels = group_levels)
  )
  annot <- as.data.frame(sample_annotation)
  rownames(annot) <- annot$sample_name

  norm_mat <- as.matrix(mat_df[, sample_cols])
  mode(norm_mat) <- "numeric"
  rownames(norm_mat) <- mat_df$Protein.Group
  colnames(norm_mat) <- sample_cols

  message("Loading comparison file...")
  comp_sep <- detect_delimiter(comparison_file)
  comparison_df <- read.table(comparison_file, sep = comp_sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  required_comp_cols <- c("condition_A", "condition_B", "comparison_label")
  missing_comp <- setdiff(required_comp_cols, colnames(comparison_df))
  if (length(missing_comp) > 0) {
    stop(paste0("Comparison file is missing expected column(s): ", paste(missing_comp, collapse = ", ")), call. = FALSE)
  }

  # ---- imputation (always after normalisation - norm_mat is already normalised) ----
  group_vec <- annot[colnames(norm_mat), "group"]
  set.seed(42)
  imputed_mat <- switch(
    tolower(impute_method),
    "group_downshift" = group_downshift_impute(norm_mat, group_vec, impute_downshift, impute_width),
    "perseus_default" = perseus_default_impute(norm_mat, impute_downshift, impute_width),
    norm_mat
  )
  rownames(imputed_mat) <- rownames(norm_mat)
  message(paste("Missing values after imputation:", sum(is.na(imputed_mat))))

  # ---- PCA ----
  if (sum(is.na(imputed_mat)) == 0 && nrow(imputed_mat) > 1) {
    pca <- prcomp(t(imputed_mat), scale. = TRUE)
    var_e <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    pca_scores <- as.data.frame(pca$x[, 1:2]) |>
      tibble::rownames_to_column("sample_name") |>
      dplyr::left_join(sample_annotation, by = "sample_name")

    p_pca <- ggplot(pca_scores, aes(PC1, PC2, colour = group)) +
      geom_point(size = 4, alpha = 0.9) +
      labs(title = "PCA", x = paste0("PC1 (", var_e[1], "%)"), y = paste0("PC2 (", var_e[2], "%)"), colour = "Group") +
      theme_bw(base_size = 12)
    save_all_formats(p_pca, file.path(output_folder, "DE", "PCA_after_imputation"), width = 8, height = 7)
  } else {
    message("Skipping PCA: matrix still has missing values (choose an imputation method to enable it)")
  }

  # ---- design + comparisons (cell-means limma) ----
  message("Preparing limma design...")
  group_f <- factor(annot[colnames(imputed_mat), "group"], levels = group_levels)
  clean_levels <- vapply(group_levels, safe_name, character(1))
  levels(group_f) <- clean_levels
  design <- model.matrix(~ 0 + group_f)
  colnames(design) <- clean_levels

  fit <- lmFit(imputed_mat, design)

  comparisons <- list()
  for (i in seq_len(nrow(comparison_df))) {
    label <- comparison_df$comparison_label[i]
    comparisons[[label]] <- c(comparison_df$condition_A[i], comparison_df$condition_B[i])
  }

  de_list <- list(); ud_summary <- list()
  for (comp_label in names(comparisons)) {
    condition_A <- safe_name(comparisons[[comp_label]][1])
    condition_B <- safe_name(comparisons[[comp_label]][2])
    if (!(condition_A %in% clean_levels) || !(condition_B %in% clean_levels)) {
      message(paste0("[SKIP] '", comp_label, "': condition not found among annotation groups (",
                     paste(group_levels, collapse = ", "), ")"))
      next
    }

    contrast_formula <- paste0(condition_A, "-", condition_B)
    contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
    fit2 <- eBayes(contrasts.fit(fit, contrast_matrix), trend = TRUE)

    de <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P") |>
      tibble::rownames_to_column("Protein.Group") |>
      dplyr::left_join(gene_map, by = "Protein.Group") |>
      dplyr::mutate(label = dplyr::if_else(!is.na(primary_gene) & primary_gene != "", primary_gene, Protein.Group),
                    Sig = adj.P.Val < fdr_cutoff & abs(logFC) > fc_threshold_log2,
                    Direction = dplyr::case_when(
                      Sig & logFC >  fc_threshold_log2 ~ "Up",
                      Sig & logFC < -fc_threshold_log2 ~ "Down",
                      TRUE ~ "NS"))
    de_list[[comp_label]] <- de
    n_up <- sum(de$Direction == "Up"); n_dn <- sum(de$Direction == "Down")
    ud_summary[[comp_label]] <- data.frame(contrast = comp_label, up = n_up, down = n_dn)
    message(sprintf("[DE] %-30s up %d  down %d", comp_label, n_up, n_dn))

    readr::write_tsv(de, file.path(output_folder, "DE", paste0("DE_", comp_label, ".tsv")))

    # ---- volcano ----
    df_ns  <- dplyr::filter(de, Direction == "NS")
    df_sig <- dplyr::filter(de, Direction != "NS")
    target_hits <- if (!is.null(target_genes) && length(target_genes) > 0) dplyr::filter(de, label %in% target_genes) else de[0, ]

    p_vol <- ggplot(mapping = aes(logFC, -log10(adj.P.Val))) +
      geom_point(data = df_ns, colour = ns_col, alpha = 0.5, size = 1.3) +
      geom_point(data = df_sig, aes(colour = Direction), alpha = 0.75, size = 1.6) +
      geom_vline(xintercept = c(-fc_threshold_log2, fc_threshold_log2), linetype = "dashed", colour = "grey50") +
      geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", colour = "grey50") +
      scale_colour_manual(values = c("Up" = up_col, "Down" = dn_col, "NS" = "grey60"), name = NULL) +
      ggrepel::geom_text_repel(data = dplyr::slice_min(df_sig, adj.P.Val, n = 12),
                               aes(label = label), size = 2.1, max.overlaps = 20, colour = "grey30") +
      { if (nrow(target_hits) > 0) ggrepel::geom_text_repel(data = target_hits, aes(label = label),
               size = 2.6, fontface = "bold", colour = "black", max.overlaps = Inf,
               box.padding = 0.6, point.padding = 0.4, min.segment.length = 0) } +
      { if (nrow(target_hits) > 0) geom_point(data = target_hits, shape = 21, size = 2.6,
               fill = NA, colour = "black", stroke = 0.9) } +
      annotate("text", x = Inf,  y = Inf, label = paste0("Up: ", n_up),   hjust = 1.1, vjust = 1.5, size = 3.5, colour = up_col, fontface = "bold") +
      annotate("text", x = -Inf, y = Inf, label = paste0("Down: ", n_dn), hjust = -0.1, vjust = 1.5, size = 3.5, colour = dn_col, fontface = "bold") +
      labs(title = comp_label,
           subtitle = paste0("FDR < ", fdr_cutoff, "  |  |log2FC| > ", round(fc_threshold_log2, 3), "  |  limma eBayes(trend)"),
           x = "log2 fold change", y = "-log10 adjusted P") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
    save_all_formats(p_vol, file.path(output_folder, "DE", paste0("volcano_", comp_label)), width = 9, height = 8)

    # ---- top up/down genes (gene-deduped, most significant per gene) ----
    sig <- de |> dplyr::filter(Sig) |> dplyr::arrange(adj.P.Val) |> dplyr::distinct(label, .keep_all = TRUE)
    if (nrow(sig) > 0) {
      up <- sig |> dplyr::filter(logFC > 0) |> dplyr::slice_max(logFC, n = top_n_genes)
      dn <- sig |> dplyr::filter(logFC < 0) |> dplyr::slice_min(logFC, n = top_n_genes)
      top <- dplyr::bind_rows(dn, up) |> dplyr::arrange(logFC) |>
        dplyr::mutate(label = factor(label, levels = label), dir = ifelse(logFC > 0, "Up", "Down"))
      p_top <- ggplot(top, aes(logFC, label, fill = dir)) +
        geom_col(colour = "grey40", linewidth = 0.3) +
        geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
        scale_fill_manual(values = c("Up" = up_col, "Down" = dn_col), guide = "none") +
        labs(title = paste0("Top up / down genes | ", comp_label), x = "log2 fold change", y = NULL) +
        theme_bw(base_size = 11)
      save_all_formats(p_top, file.path(output_folder, "DE", paste0("top_genes_", comp_label)),
                       width = 8, height = max(5, 0.32 * nrow(top)))
    }
  }

  if (length(de_list) == 0) {
    stop("No comparison could be matched to an annotation group - check comparison_file against annotation_file", call. = FALSE)
  }

  # ---- combined up/down bar across comparisons ----
  ud <- dplyr::bind_rows(ud_summary) |>
    dplyr::mutate(contrast = factor(contrast, levels = names(de_list))) |>
    tidyr::pivot_longer(c(up, down), names_to = "dir", values_to = "n") |>
    dplyr::mutate(dir = factor(ifelse(dir == "up", "Up", "Down"), levels = c("Up", "Down")))
  p_ud <- ggplot(ud, aes(contrast, n, fill = dir)) +
    geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey40", linewidth = 0.3) +
    geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.3, size = 3, fontface = "bold") +
    scale_fill_manual(values = c("Up" = up_col, "Down" = dn_col), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Significant proteins up / down by contrast",
         subtitle = paste0("FDR < ", fdr_cutoff, ", |log2FC| > ", round(fc_threshold_log2, 3)),
         x = NULL, y = "Protein groups") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
  save_all_formats(p_ud, file.path(output_folder, "DE", "DE_up_down_by_contrast"), width = 11, height = 6)

  writexl::write_xlsx(setNames(de_list, substr(names(de_list), 1, 31)),
                      file.path(output_folder, "DE", "DE_all_contrasts.xlsx"))

  # ---- target-gene marker panels (optional) ----
  if (!is.null(target_genes) && length(target_genes) > 0) {
    gm <- gene_map
    gm$n_valid <- rowSums(!is.na(norm_mat))[match(gm$Protein.Group, rownames(norm_mat))]
    best_pg <- gm |> dplyr::filter(primary_gene %in% target_genes) |>
      dplyr::group_by(primary_gene) |> dplyr::slice_max(n_valid, n = 1, with_ties = FALSE) |> dplyr::ungroup()

    if (nrow(best_pg) > 0) {
      marker_long <- as.data.frame(imputed_mat[best_pg$Protein.Group, , drop = FALSE]) |>
        tibble::rownames_to_column("Protein.Group") |>
        tidyr::pivot_longer(-Protein.Group, names_to = "sample_name", values_to = "log2_intensity") |>
        dplyr::left_join(sample_annotation, by = "sample_name") |>
        dplyr::mutate(gene = best_pg$primary_gene[match(Protein.Group, best_pg$Protein.Group)],
                      gene = factor(gene, levels = target_genes[target_genes %in% best_pg$primary_gene]))

      p_markers <- ggplot(marker_long, aes(group, log2_intensity, fill = group)) +
        geom_boxplot(width = 0.55, alpha = 0.6, outlier.shape = NA) +
        geom_jitter(width = 0.12, size = 2, alpha = 0.85, shape = 21, colour = "grey20") +
        facet_wrap(~ gene, scales = "free_y", nrow = 1) +
        labs(title = "Target gene markers", subtitle = "points = replicates", x = NULL, y = "log2 Intensity", fill = "Group") +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
              axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
      save_all_formats(p_markers, file.path(output_folder, "DE", "target_markers_panel"), width = max(5, 3 * nrow(best_pg)), height = 5)

      marker_de <- dplyr::bind_rows(lapply(names(de_list), function(cn)
        de_list[[cn]] |> dplyr::filter(label %in% target_genes) |>
          dplyr::transmute(comparison = cn, gene = label, logFC, P.Value, adj.P.Val)))
      readr::write_tsv(marker_de, file.path(output_folder, "DE", "target_markers_DE.tsv"))
    } else {
      message("No target genes found in the normalised matrix - skipping marker panels")
    }
  }

  # ---- curtain exports ----
  curtain_intensity <- cbind(gene_map, as.data.frame(imputed_mat)) |>
    dplyr::rename(`Protein IDs` = Protein.Group, `Gene names` = Genes) |>
    dplyr::select(-primary_gene, -Category)
  readr::write_tsv(curtain_intensity, file.path(output_folder, "curtain", "curtain_intensity_matrix.tsv"))

  for (comp_label in names(de_list)) {
    de_list[[comp_label]] |>
      dplyr::rename(`Protein IDs` = Protein.Group, `Gene names` = Genes,
                    `log2 fold change` = logFC, `p value` = P.Value, `p value adjusted` = adj.P.Val) |>
      dplyr::select(`Protein IDs`, `Gene names`, primary_gene, Category,
                    `log2 fold change`, `p value`, `p value adjusted`, AveExpr, t, B, Sig, Direction) |>
      readr::write_tsv(file.path(output_folder, "curtain", paste0("curtain_DE_", comp_label, ".tsv")))
  }

  message("Differential expression complete.")
}

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (startsWith(arg, "--")) {
      key <- substring(arg, 3)
      if (i < length(args) && !startsWith(args[i + 1], "--")) {
        parsed[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        parsed[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  parsed
}

params <- parse_args(args)

normalized_matrix_file <- params$normalized_matrix_file
annotation_file <- params$annotation_file
comparison_file <- params$comparison_file
output_folder <- params$output_folder
impute_method <- ifelse(is.null(params$impute_method), "none", params$impute_method)
impute_downshift <- ifelse(is.null(params$impute_downshift), 1.8, as.numeric(params$impute_downshift))
impute_width <- ifelse(is.null(params$impute_width), 0.3, as.numeric(params$impute_width))
fdr_cutoff <- ifelse(is.null(params$fdr_cutoff), 0.05, as.numeric(params$fdr_cutoff))
fc_threshold_log2 <- ifelse(is.null(params$fc_threshold_log2), log2(1.5), as.numeric(params$fc_threshold_log2))
target_genes <- if (is.null(params$target_genes) || params$target_genes == "") {
  NULL
} else {
  trimws(unlist(strsplit(params$target_genes, ",")))
}
top_n_genes <- ifelse(is.null(params$top_n_genes), 15, as.numeric(params$top_n_genes))

if (is.null(normalized_matrix_file) || is.null(annotation_file) || is.null(comparison_file) || is.null(output_folder)) {
  stop("Missing required arguments: normalized_matrix_file, annotation_file, comparison_file, output_folder", call. = FALSE)
}

run_differential_expression(
  normalized_matrix_file = normalized_matrix_file,
  annotation_file = annotation_file,
  comparison_file = comparison_file,
  output_folder = output_folder,
  impute_method = impute_method,
  impute_downshift = impute_downshift,
  impute_width = impute_width,
  fdr_cutoff = fdr_cutoff,
  fc_threshold_log2 = fc_threshold_log2,
  target_genes = target_genes,
  top_n_genes = top_n_genes
)
