dirnameScript <- function(){
    # get full directory path of current script located
    cmd = commandArgs(trailingOnly = FALSE)
    scriptName = sub('--file=', "", cmd[grep('^--file=', cmd)])
    if (length(scriptName) > 0) {
        path = normalizePath(scriptName)
        dirname = dirname(path)
    } else {
        print('Warning: not a runtime environment, using current directory instead.')
        dirname = getwd()
    }
    return(dirname)
}

print(dirnameScript())


get_color <- function (x, n = 100,palette = "Paired", palcolor = NULL, type = "auto", 
          matched = FALSE, reverse = FALSE, NA_keep = FALSE, NA_color = "grey80") 
{
  palette_list <- readRDS(str_glue("{dirnameScript()}/COLORS/colorDB.rds"))
  if (missing(x)) {
    x <- 1:n
    type <- "continuous"
  }
  if (!palette %in% names(palette_list)) {
    stop("The palette name (", palette, ") is invalid! You can check the available palette names with 'show_palettes()'. Or pass palette colors via the 'palcolor' parameter.")
  }
  if (is.list(palcolor)) {
    palcolor <- unlist(palcolor)
  }
  if (all(palcolor == "")) {
    palcolor <- palette_list[[palette]]
  }
  if (is.null(palcolor) || length(palcolor) == 0) {
    palcolor <- palette_list[[palette]]
  }
  if (!is.null(names(palcolor))) {
    if (all(x %in% names(palcolor))) {
      palcolor <- palcolor[intersect(names(palcolor), x)]
    }
  }
  pal_n <- length(palcolor)
  if (!type %in% c("auto", "discrete", "continuous")) {
    stop("'type' must be one of 'auto','discrete' and 'continuous'.")
  }
  if (type == "auto") {
    if (is.numeric(x)) {
      type <- "continuous"
    }
    else {
      type <- "discrete"
    }
  }
  if (type == "discrete") {
    if (!is.factor(x)) {
      x <- factor(x, levels = unique(x))
    }
    n_x <- nlevels(x)
    if (isTRUE(attr(palcolor, "type") == "continuous")) {
      color <- colorRampPalette(palcolor)(n_x)
    }
    else {
      color <- ifelse(rep(n_x, n_x) <= pal_n, palcolor[1:n_x], 
                      colorRampPalette(palcolor)(n_x))
    }
    names(color) <- levels(x)
    if (any(is.na(x))) {
      color <- c(color, setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      color <- color[x]
      color[is.na(color)] <- NA_color
    }
  }
  else if (type == "continuous") {
    if (!is.numeric(x) && all(!is.na(x))) {
      stop("'x' must be type of numeric when use continuous color palettes.")
    }
    if (all(is.na(x))) {
      values <- as.factor(rep(0, n))
    }
    else if (length(unique(na.omit(as.numeric(x)))) == 1) {
      values <- as.factor(rep(unique(na.omit(as.numeric(x))), 
                              n))
    }
    else {
      if (isTRUE(matched)) {
        values <- cut(x, breaks = seq(min(x, na.rm = TRUE), 
                                      max(x, na.rm = TRUE), length.out = n + 1), 
                      include.lowest = TRUE)
      }
      else {
        values <- cut(1:100, breaks = seq(min(x, na.rm = TRUE), 
                                          max(x, na.rm = TRUE), length.out = n + 1), 
                      include.lowest = TRUE)
      }
    }
    n_x <- nlevels(values)
    color <- ifelse(rep(n_x, n_x) <= pal_n, palcolor[1:n_x], 
                    colorRampPalette(palcolor)(n_x))
    names(color) <- levels(values)
    if (any(is.na(x))) {
      color <- c(color, setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      if (all(is.na(x))) {
        color <- NA_color
      }
      else if (length(unique(na.omit(x))) == 1) {
        color <- color[as.character(unique(na.omit(x)))]
        color[is.na(color)] <- NA_color
      }
      else {
        color <- color[as.character(values)]
        color[is.na(color)] <- NA_color
      }
    }
  }
  if (isTRUE(reverse)) {
    color <- rev(color)
  }
  if (!isTRUE(NA_keep)) {
    color <- color[names(color) != "NA"]
  }
  return(color)
}