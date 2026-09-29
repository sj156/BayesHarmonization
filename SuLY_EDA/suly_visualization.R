# ============================================================
# suly_visualization.R
# 把 bucket_data_visualization.R 的三图可视化，应用到
# SuLY-MPRAGE-Handoff/ 里的 MPRAGE 扫描上。
#
# 对每个扫描导出三种图，按被试 / 机型分目录存放:
#   SuLY_EDA/viz_png/<subject>/<magnet>/<param>/
#     - <scan>_three_views.png     sagittal / coronal / axial 三正交中间层
#     - <scan>_axial_montage.png   axial 方向多层拼图 (montage)
#     - <scan>_mip_axis3.png       沿 axis 3 的最大强度投影 (MIP)
#
# 与 bucket 版的区别：
#   1. 输入是 5 层嵌套目录，扫描名从文件名推出 (去掉 _native_coreg_mskd)
#   2. 数据是 skull-stripped 脑 (非零体素约 9%)，montage 只在 mask 的
#      bounding box 内取层，避免拼图里一半是全黑
#   3. 色标 limits 仍按每个扫描自身的脑内 1%-99% 分位来定 —— 刻意不做
#      跨机器归一化，好让色标条上的数字直接暴露厂商间的量级差
#      (GE 中位数 ~1716, Vida/Prisma ~253, UIH ~152)
#
# 用法:
#   Rscript suly_visualization.R            # 全量 276 个扫描
#   Rscript suly_visualization.R sub-0011   # 只跑一个被试
# ============================================================

library(RNifti)
library(ggplot2)
library(patchwork)

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
out_dir  <- file.path(proj_dir, "SuLY_EDA", "viz_png")

args   <- commandArgs(trailingOnly = TRUE)
subpat <- if (length(args) >= 1) args[1] else "*"

files <- Sys.glob(file.path(data_dir, subpat, "*", "*", "*", "*", "*.nii.gz"))
files <- files[!dir.exists(files)]

stopifnot(length(files) > 0)

cat("Found", length(files), "scans matching subject pattern:", subpat, "\n")


# ============================================================
# Helper functions (与 bucket_data_visualization.R 一致)
# ============================================================

# 扫描名：去掉 .nii.gz，再去掉 _native_coreg_mskd 后缀
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

# 新增：脑 mask 沿某个轴的 bounding box (第一个/最后一个有非零信号的层)
brain_range <- function(vol, axis = 3) {
  prof <- apply(vol > 0, axis, sum)
  idx  <- which(prof > 0)
  if (length(idx) == 0) return(c(1L, dim(vol)[axis]))
  range(idx)
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
    scale_fill_viridis_c(option = "magma", na.value = "black",
                         limits = limits, oob = scales::squish) +
    labs(title = ttl, fill = "Intensity") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "right")
}


# ============================================================
# Three orthogonal middle views
# ============================================================

plot_three_views <- function(file, volume = 1, clip_q = c(0.01, 0.99)) {

  arr <- as.array(readNifti(file))

  # 色标 limits 按该扫描自身的脑内强度分布 —— 不做跨机器归一化
  b <- arr[arr > 0 & is.finite(arr)]
  lims <- as.numeric(quantile(b, c(0.01, 0.99)))

  p_sag <- plot_slice_array(arr, plane = "sagittal", volume = volume,
                            title = "Sagittal", clip_q = clip_q, limits = lims)
  p_cor <- plot_slice_array(arr, plane = "coronal", volume = volume,
                            title = "Coronal", clip_q = clip_q, limits = lims)
  p_axi <- plot_slice_array(arr, plane = "axial", volume = volume,
                            title = "Axial", clip_q = clip_q, limits = lims)

  (p_sag | p_cor | p_axi) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(strip_ext(file), "   volume = ", volume),
      theme = theme(plot.background = element_rect(fill = "white", color = NA))
    ) &
    theme(legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))
}


# ============================================================
# Montage of multiple slices (限制在脑 mask 的 bounding box 内)
# ============================================================

plot_montage <- function(file, plane = c("axial", "coronal", "sagittal"),
                         nslices = 25, volume = 1, clip_q = c(0.01, 0.99)) {

  plane <- match.arg(plane)

  arr <- as.array(readNifti(file))
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
  df$intensity <- clip_values(df$intensity, clip_q = clip_q)

  ggplot(df, aes(x = x, y = y, fill = intensity)) +
    geom_raster(interpolate = FALSE) +
    coord_equal() +
    scale_y_reverse() +
    scale_fill_viridis_c(option = "magma", na.value = "black") +
    facet_wrap(~ slice, ncol = ceiling(sqrt(length(slices)))) +
    labs(title = paste0(strip_ext(file), " | ", plane,
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

plot_mip <- function(file, axis = 3, volume = 1, clip_q = c(0.01, 0.99)) {

  axis <- as.integer(axis)
  stopifnot(axis %in% c(1, 2, 3))

  arr <- as.array(readNifti(file))
  vol <- get_volume(arr, volume = volume)

  keep_dims <- setdiff(1:3, axis)
  mat <- apply(vol, keep_dims, max, na.rm = TRUE)

  df <- matrix_to_df(mat)
  df$intensity <- clip_values(df$intensity, clip_q = clip_q)

  ggplot(df, aes(x = x, y = y, fill = intensity)) +
    geom_raster(interpolate = FALSE) +
    coord_equal() +
    scale_y_reverse() +
    scale_fill_viridis_c(option = "magma", na.value = "black") +
    labs(title = paste0(strip_ext(file), " | MIP over axis ", axis),
         fill = "Intensity") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "right")
}


# ============================================================
# Batch export
# ============================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()

for (i in seq_along(files)) {

  f    <- files[i]
  name <- strip_ext(f)

  message(sprintf("[%d/%d] %s", i, length(files), name))

  # 输出目录按被试 / 机型 / 参数分开：
  #   viz_png/<subject>/<magnet>/<param>/<scan>_*.png
  # 从输入路径的 5 层结构里取 subject / magnet / param
  parts   <- strsplit(normalizePath(f, winslash = "/"), "/")[[1]]
  np      <- length(parts)
  subject <- parts[np - 5]   # sub-XXXX
  magnet  <- parts[np - 4]   # 机型
  param   <- parts[np - 2]   # ParamNN

  file_dir <- file.path(out_dir, subject, magnet, param)
  dir.create(file_dir, showWarnings = FALSE, recursive = TRUE)

  ggsave(file.path(file_dir, paste0(name, "_three_views.png")),
         plot = plot_three_views(f, volume = 1),
         width = 12, height = 4, dpi = 200)

  ggsave(file.path(file_dir, paste0(name, "_axial_montage.png")),
         plot = plot_montage(f, plane = "axial", nslices = 25, volume = 1),
         width = 10, height = 10, dpi = 200)

  ggsave(file.path(file_dir, paste0(name, "_mip_axis3.png")),
         plot = plot_mip(f, axis = 3, volume = 1),
         width = 7, height = 5, dpi = 200)
}

message("Done in ", round(difftime(Sys.time(), t0, units = "mins"), 1),
        " min. PNG files saved to: ", out_dir)
