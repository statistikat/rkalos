seats_from_tramo <- function(series, ssq, effective, estimated, d = 0L, bd = 0L,
    phi = numeric(), theta = numeric(), bphi = numeric(), btheta = numeric(),
    log = FALSE, mean_correction = TRUE, outliers = NULL, forecast = -2L, backcast = 0L) {
  spec <- as.list(environment())[setdiff(names(formals(seats_from_tramo)), c("series", "outliers"))]
  domain <- seats_model_domain(series, spec)
  spec$effective <- native_horizon(effective, "effective")
  spec$estimated <- native_horizon(estimated, "estimated")
  if (!is.numeric(ssq) || is.complex(ssq) || !is.null(dim(ssq)) ||
      length(ssq) != 1L || !is.finite(ssq) || ssq <= 0 ||
      spec$effective > length(series) || spec$effective <= spec$estimated ||
      !is.null(dim(effective)) || !is.null(dim(estimated)))
    stop("The fitted likelihood requires positive ssq and 0 <= estimated < effective <= length(series).",
         call. = FALSE)
  if (is.null(outliers)) outliers <- data.frame(type = character(), position = integer(), coefficient = numeric())
  if (!is.data.frame(outliers) || !setequal(names(outliers), c("type", "position", "coefficient")) ||
      length(names(outliers)) != 3L || !is.character(outliers$type) || !is.null(dim(outliers$type)) ||
      any(!outliers$type %in% c("AO", "LS")) ||
      !is.numeric(outliers$position) || is.complex(outliers$position) || !is.null(dim(outliers$position)) ||
      any(!is.finite(outliers$position)) || any(outliers$position != floor(outliers$position)) ||
      any(outliers$position < 1L | outliers$position > length(series)) ||
      !is.numeric(outliers$coefficient) || is.complex(outliers$coefficient) ||
      !is.null(dim(outliers$coefficient)) || any(!is.finite(outliers$coefficient)))
    stop("outliers must have AO or LS type, one-based observed position and finite coefficient.", call. = FALSE)
  outliers$position <- as.integer(outliers$position)
  result <- .native_seats_from_tramo(spec, domain, as.numeric(series), outliers)
  seats_date_components(result, domain$frequency, c("initial", "final", "merged"))
}
