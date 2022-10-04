# From https://github.com/satijalab/seurat/blob/master/R/preprocessing.R and https://github.com/satijalab/seurat/blob/master/R/convenience.R
require(httr)
require(Matrix)

ReadMtx <- function(
  mtx,
  cells,
  features,
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
) {
  all.files <- list(
    "expression matrix" = mtx,
    "barcode list" = cells,
    "feature list" = features
  )
  for (i in seq_along(along.with = all.files)) {
    uri <- tryCatch(
      expr = {
        con <- url(description = all.files[[i]])
        close(con = con)
        all.files[[i]]
      },
      error = function(...) {
        return(normalizePath(path = all.files[[i]], winslash = '/'))
      }
    )
    err <- paste("Cannot find", names(x = all.files)[i], "at", uri)
    uri <- build_url(url = parse_url(url = uri))
    if (grepl(pattern = '^[A-Z]?:///', x = uri)) {
      uri <- gsub(pattern = '^://', replacement = '', x = uri)
      if (!file.exists(uri)) {
        stop(err, call. = FALSE)
      }
    } else {
      if (!Online(url = uri, seconds = 2L)) {
        stop(err, call. = FALSE)
      }
      if (file_ext(uri) == 'gz') {
        con <- url(description = uri)
        uri <- gzcon(con = con, text = TRUE)
      }
    }
    all.files[[i]] <- uri
  }
  cell.barcodes <- read.table(
    file = all.files[['barcode list']],
    header = FALSE,
    sep = cell.sep,
    row.names = NULL,
    skip = skip.cell
  )
  feature.names <- read.table(
    file = all.files[['feature list']],
    header = FALSE,
    sep = feature.sep,
    row.names = NULL,
    skip = skip.feature
  )
  # read barcodes
  bcols <- ncol(x = cell.barcodes)
  if (bcols < cell.column) {
    stop(
      "cell.column was set to ",
      cell.column,
      " but ",
      cells,
      " only has ",
      bcols,
      " columns.",
      " Try setting the cell.column argument to a value <= to ",
      bcols,
      "."
    )
  }
  cell.names <- cell.barcodes[, cell.column]
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  # read features
  fcols <- ncol(x = feature.names)
  if (fcols < feature.column) {
    stop(
      "feature.column was set to ",
      feature.column,
      " but ",
      features,
      " only has ",
      fcols, " column(s).",
      " Try setting the feature.column argument to a value <= to ",
      fcols,
      "."
    )
  }
  if (any(is.na(x = feature.names[, feature.column]))) {
    na.features <- which(x = is.na(x = feature.names[, feature.column]))
    replacement.column <- ifelse(test = feature.column == 2, yes = 1, no = 2)
    if (replacement.column > fcols) {
      stop(
        "Some features names are NA in column ",
        feature.column,
        ". Try specifiying a different column.",
        call. = FALSE
        )
    } else {
      warning(
        "Some features names are NA in column ",
        feature.column,
        ". Replacing NA names with ID from column ",
        replacement.column,
        ".",
        call. = FALSE
        )
    }
    feature.names[na.features, feature.column] <- feature.names[na.features, replacement.column]
  }
  feature.names <- feature.names[, feature.column]
  if (unique.features) {
    feature.names <- make.unique(names = feature.names)
  }
  data <- readMM(file = all.files[['expression matrix']])
  if (mtx.transpose) {
    data <- t(x = data)
  }
  if (length(x = cell.names) != ncol(x = data)) {
    stop(
      "Matrix has ",
      ncol(data),
      " columns but found ", length(cell.names),
      " barcodes. ",
      ifelse(
        test = length(x = cell.names) > ncol(x = data),
        yes = "Try increasing `skip.cell`. ",
        no = ""
      ),
      call. = FALSE
      )
  }
  if (length(x = feature.names) != nrow(x = data)) {
    stop(
      "Matrix has ",
      nrow(data),
      " rows but found ", length(feature.names),
      " features. ",
      ifelse(
        test = length(x = feature.names) > nrow(x = data),
        yes = "Try increasing `skip.feature`. ",
        no = ""
      ),
      call. = FALSE
      )
  }

  colnames(x = data) <- cell.names
  rownames(x = data) <- feature.names
  data <- as(data, Class = "dgCMatrix")
  return(data)
}

ReadSTARsolo <- function(data.dir, ... ) {
  mtx <- file.path(data.dir, "matrix.mtx")
  cells <- file.path(data.dir, "barcodes.tsv")
  features <- file.path(data.dir, "features.tsv")
  return(ReadMtx(mtx = mtx, cells = cells, features = features, ...))
}