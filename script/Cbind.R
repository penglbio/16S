padNA <- function (mydata, rowsneeded, first = TRUE) 
{
  temp1 = colnames(mydata)
  rowsneeded = rowsneeded - nrow(mydata)
  temp2 = setNames(
    data.frame(matrix(rep(NA, length(temp1) * rowsneeded), 
                      ncol = length(temp1))), temp1)
  if (isTRUE(first)) rbind(mydata, temp2)
  else rbind(temp2, mydata)
}

dotnames <- function(...) {
  vnames <- as.list(substitute(list(...)))[-1L]
  vnames <- unlist(lapply(vnames,deparse), FALSE, FALSE)
  vnames
}

Cbind <- function(..., first = TRUE) {
  Names <- dotnames(...)
  datalist <- setNames(list(...), Names)
  nrows <- max(sapply(datalist, function(x) 
    ifelse(is.null(dim(x)), length(x), nrow(x))))
  datalist <- lapply(seq_along(datalist), function(x) {
    z <- datalist[[x]]
    if (is.null(dim(z))) {
      z <- setNames(data.frame(z), Names[x])
    } else {
      if (is.null(colnames(z))) {
        colnames(z) <- paste(Names[x], sequence(ncol(z)), sep = "_")
      } else {
        colnames(z) <- paste(Names[x], colnames(z), sep = "_")
      }
    }
    padNA(z, rowsneeded = nrows, first = first)
  })
  do.call(cbind, datalist)
}