# ============================================================
# suly_visualization_standardized.R
# 与 suly_visualization.R 相同的三图可视化，但先对每个扫描做强度标准化：
#   median 版: x / median(脑内非零体素) - 1   -> 脑内中位数处 = 0
#   mean   版: x / mean(脑内非零体素)   - 1   -> 脑内均值处   = 0
# 统计量按每个扫描各自计算；脑外 (原始值 0) 设为 NA。
#
# 色标: 蓝 (负, 比参考值暗) - 白 (0) - 红 (正, 比参考值亮)，脑外浅灰。
# 上下限以 0 对称，取脑内 1%/99% 分位绝对值的较大者。
#
# 输出目录结构与 viz_png 一致:
#   SuLY_EDA/viz_png_standardized_median/<subject>/<magnet>/<param>/
#   SuLY_EDA/viz_png_standardized_mean/<subject>/<magnet>/<param>/
#     - <scan>_three_views.png     sagittal / coronal / axial 三正交中间层
#     - <scan>_axial_montage.png   axial 方向多层拼图 (montage)
#     - <scan>_mip_axis3.png       沿 axis 3 的最大强度投影 (MIP)
#
# 用法:
#   Rscript suly_visualization_standardized.R                   # 全量, median + mean
#   Rscript suly_visualization_standardized.R sub-0011          # 只跑一个被试
#   Rscript suly_visualization_standardized.R "*" median        # 只跑 median 版
# ============================================================

suppressPackageStartupMessages({
  library(RNifti)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# 1. 路径 (以脚本所在目录的上一级作为项目根目录)
# ----------------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
proj_dir <- if (length(script_arg) == 1) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  "C:/Users/yaoya/OneDrive - CUHK-Shenzhen/桌面/Huaxi"
}
data_dir <- file.path(proj_dir, "SuLY-MPRAGE-Handoff")
eda_dir  <- file.path(proj_dir, "SuLY_EDA")

args    <- commandArgs(trailingOnly = TRUE)
subpat  <- if (length(args) >= 1) args[1] else "*"
methods <- if (length(args) >= 2) args[2] else c("median", "mean")
stopifnot(all(methods %in% c("median", "mean")))

out_dirs <- setNames(file.path(eda_dir, paste0("viz_png_standardized_", methods)),
                     methods)

files <- Sys.glob(file.path(data_dir, subpat, "*", "*", "*", "*", "*.nii.gz"))
files <- files[!dir.exists(files)]

stopifnot(length(files) > 0)

cat("Found", length(files), "scans matching subject pattern:", subpat,
    "| methods:", paste(methods, collapse = ", "), "\n")


# ============================================================
# Helper functions (与 suly_visualization.R 一致)
# ============================================================

strip_ext <- function(f) {
  n <- sub("\\.nii(\\.gz)?$", "", basename(f), ignore.case = TRUE)
  sub("_native_coreg_mskd$", "", n)
}

get_volume <- function(x, volume = 1) {
  d <- dim(x)
  if (length(d) == 2) return(array(x, dim = c(d[1], d[2], 1)))
  if (length(d) == 3) return(x)
  if (length(d) == 4) {
    if (volume < 1 || volume > d[4]) stop("volume index out of range")
    return(x[, , , volume, drop = TRUE])
  }
  stop("Only 2D, 3D, or 4D NIfTI arrays are supported.")
}

slice_matrix <- function(vol, plane = c("axial", "coronal", "sagittal"),
                         slice = NULL) {
  plane <- match.arg(plane)
  d <- dim(vol)
  if (length(d) != 3) stop("Input volume must be 3D.")

  max_slice <- switch(plane, axial = d[3], coronal = d[2], sagittal = d[1])
  if (is.null(slice)) slice <- ceiling(max_slice / 2)
  slice <- as.integer(slice)
  if (slice < 1 || slice > max_slice) stop("slice index out of range")

  mat <- switch(plane,
                axial    = vol[, , slice],
                coronal  = vol[, slice, ],
                sagittal = vol[slice, , ])

  list(mat = mat, slice = slice, plane = plane)
}

matrix_to_df <- function(mat) {
  df <- expand.grid(x = seq_len(nrow(mat)), y = seq_len(ncol(mat)))
  df$intensity <- as.vector(mat)
  df
}

clip_values <- function(v, clip_q = c(0.01, 0.99)) {
  if (is.null(clip_q)) return(v)
  z <- v[is.finite(v)]
  if (length(z) < 2) return(v)
  qs <- quantile(z, probs = clip_q, na.rm = TRUE, names = FALSE)
  if (is.finite(qs[1]) && is.finite(qs[2]) && qs[2] > qs[1]) {
    v <- pmin(pmax(v, qs[1]), qs[2])
  }
  v
}

# 标准化后脑外是 NA、脑内有正有负，所以 mask 用 is.finite 而不是 > 0
brain_range <- function(vol, axis = 3) {
  prof <- apply(is.finite(vol), axis, sum)
  idx  <- which(prof > 0)
  if (length(idx) == 0) return(c(1L, dim(vol)[axis]))
  range(idx)
}

# 新增：按脑内 (非零、有限) 体素的中位数或均值做比值标准化，再减 1
#   0 = 等于参考值, 0.3 = 比参考值高 30%, -0.5 = 低 50%
# 脑外 (原始值 = 0) 设为 NA，画成灰色背景
standardize <- function(arr, method = c("median", "mean")) {
  method <- match.arg(method)
  inside <- arr > 0 & is.finite(arr)
  ref <- switch(method, median = median(arr[inside]), mean = mean(arr[inside]))
  out <- arr / ref - 1
  out[!inside] <- NA
  list(arr = out, ref = ref)
}

# 新增：以 0 为中心的对称色标上限 (脑内 1%/99% 分位绝对值的较大者)
sym_limits <- function(arr, clip_q = c(0.01, 0.99)) {
  qs <- quantile(arr, clip_q, na.rm = TRUE, names = FALSE)
  m  <- max(abs(qs))
  c(-m, m)
}

# 新增：蓝 (负) - 白 (0) - 红 (正) 发散色标，脑外 NA 为浅灰
diverging_fill <- function(limits) {
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = limits, oob = scales::squish,
                       na.value = "grey85")
}


# ============================================================
# Plot one slice
# ============================================================

plot_slice_array <- function(arr, plane = c("axial", "coronal", "sagittal"),
                             slice = NULL, volume = 1, title = NULL,
                             clip_q = c(0.01, 0.99), limits = NULL) {

  vol <- get_volume(arr, volume = volume)
  sm  <- slice_matrix(vol, plane = plane, slice = slice)

  df <- matrix_to_df(sm$mat)
  df$intensity <- clip_values(df$intensity, clip_q = clip_q)

  ttl <- if (is.null(title)) {
    sprintf("%s slice = %d", sm$plane, sm$slice)
  } else {
    sprintf("%s slice = %d", title, sm$slice)
  }

  ggplot(df, aes(x = x, y = y, fill = intensity)) +
    geom_raster(interpolate = FALSE) +
    coord_equal() +
    scale_y_reverse() +
    diverging_fill(limits) +
    labs(title = ttl, fill = "Intensity") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "right")
}


# ============================================================
# Three orthogonal middle views
# (传入已读好的数组 arr 和标签 label，避免重复读盘)
# ============================================================

plot_three_views <- function(arr, label, volume = 1, clip_q = c(0.01, 0.99)) {

  lims <- sym_limits(arr, clip_q)

  p_sag <- plot_slice_array(arr, plane = "sagittal", volume = volume,
                            title = "Sagittal", clip_q = clip_q, limits = lims)
  p_cor <- plot_slice_array(arr, plane = "coronal", volume = volume,
                            title = "Coronal", clip_q = clip_q, limits = lims)
  p_axi <- plot_slice_array(arr, plane = "axial", volume = volume,
                            title = "Axial", clip_q = clip_q, limits = lims)

  (p_sag | p_cor | p_axi) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(label, "   volume = ", volume),
      theme = theme(plot.background = element_rect(fill = "white", color = NA))
    ) &
    theme(legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))
}


# ============================================================
# Montage of multiple slices (限制在脑 mask 的 bounding box 内)
# ============================================================

plot_montage <- function(arr, label, plane = c("axial", "coronal", "sagittal"),
                         nslices = 25, volume = 1, clip_q = c(0.01, 0.99)) {

  plane <- match.arg(plane)

  vol <- get_volume(arr, volume = volume)

  axis <- switch(plane, sagittal = 1, coronal = 2, axial = 3)
  rng  <- brain_range(vol, axis = axis)

  nslices <- min(nslices, rng[2] - rng[1] + 1)
  slices  <- unique(round(seq(rng[1], rng[2], length.out = nslices)))
  slice_levels <- paste0("slice ", slices)

  dfs <- lapply(slices, function(s) {
    sm <- slice_matrix(vol, plane = plane, slice = s)
    df <- matrix_to_df(sm$mat)
    df$slice <- factor(paste0("slice ", s), levels = slice_levels)
    df
  })

  df <- do.call(rbind, dfs)
  lims <- sym_limits(arr, clip_q)

  ggplot(df, aes(x = x, y = y, fill = intensity)) +
    geom_raster(interpolate = FALSE) +
    coord_equal() +
    scale_y_reverse() +
    diverging_fill(lims) +
    facet_wrap(~ slice, ncol = ceiling(sqrt(length(slices)))) +
    labs(title = paste0(label, " | ", plane,
                        " montage | slices ", rng[1], "-", rng[2]),
         fill = "Intensity") +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          strip.text = element_text(size = 8),
          legend.position = "right")
}


# ============================================================
# Maximum intensity projection (MIP)
# ============================================================

plot_mip <- function(arr, label, axis = 3, volume = 1, clip_q = c(0.01, 0.99)) {

  axis <- as.integer(axis)
  stopifnot(axis %in% c(1, 2, 3))

  vol <- get_volume(arr, volume = volume)

  keep_dims <- setdiff(1:3, axis)
  # 整条射线都在脑外时 max 为 -Inf，改回 NA 画成灰色背景
  mat <- suppressWarnings(apply(vol, keep_dims, max, na.rm = TRUE))
  mat[!is.finite(mat)] <- NA

  df <- matrix_to_df(mat)
  lims <- sym_limits(mat, clip_q)

  ggplot(df, aes(x = x, y = y, fill = intensity)) +
    geom_raster(interpolate = FALSE) +
    coord_equal() +
    scale_y_reverse() +
    diverging_fill(lims) +
    labs(title = paste0(label, " | MIP over axis ", axis),
         fill = "Intensity") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "right")
}


# ============================================================
# Batch export
# ============================================================

for (d in out_dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()

for (i in seq_along(files)) {

  f    <- files[i]
  name <- strip_ext(f)

  message(sprintf("[%d/%d] %s", i, length(files), name))

  parts   <- strsplit(normalizePath(f, winslash = "/"), "/")[[1]]
  np      <- length(parts)
  subject <- parts[np - 5]   # sub-XXXX
  magnet  <- parts[np - 4]   # 机型
  param   <- parts[np - 2]   # ParamNN

  raw <- as.array(readNifti(f))

  for (m in methods) {

    st    <- standardize(raw, m)
    arr   <- st$arr
    label <- sprintf("%s\nx/%s - 1 (%s = %.1f)", name, m, m, st$ref)

    file_dir <- file.path(out_dirs[[m]], subject, magnet, param)
    dir.create(file_dir, showWarnings = FALSE, recursive = TRUE)

    ggsave(file.path(file_dir, paste0(name, "_three_views.png")),
           plot = plot_three_views(arr, label, volume = 1),
           width = 12, height = 4, dpi = 200)

    ggsave(file.path(file_dir, paste0(name, "_axial_montage.png")),
           plot = plot_montage(arr, label, plane = "axial", nslices = 25, volume = 1),
           width = 10, height = 10, dpi = 200)

    ggsave(file.path(file_dir, paste0(name, "_mip_axis3.png")),
           plot = plot_mip(arr, label, axis = 3, volume = 1),
           width = 7, height = 5, dpi = 200)
  }
}

message("Done in ", round(difftime(Sys.time(), t0, units = "mins"), 1),
        " min. PNG files saved to: ", paste(out_dirs, collapse = ", "))
