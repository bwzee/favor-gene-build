
library(parallel)


#' Title
#'
#' @param df
#' @param fl
#' @param chunks
#'
#' @return
#' @export
#'
#' @examples
write_fast_json <- function(df, fl, chunks=5000) {
    df1 = Map(function(x, nm) setNames(x, rep(nm, length(x))), df, names(df))
    df1 <- lapply(df1, function(x) {
        if (is(x, "factor")) as.character(x) else x
    })
    lol = unname(do.call(Map, c(list, df1)))
    idx <- splitIndices(length(lol), chunks)
    json <- sapply(idx, function(i) jsonlite::toJSON(lol[i],pretty=T,auto_unbox = T))
    substring(json[-length(json)], nchar(json)[-length(json)]) <- ","
    substring(json[-1], 1, 1) <- ""
    writeLines(json, fl)
}

