#library(data.table)
#library(foreign)
#library(formula.tools)
#library(gtools)
#library(lmtest)
#library(sandwich)
#library(stargazer)


## System functions
object.sizes <- function(){
  return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name)
    object.size(get(object.name))))))
}

plot.object.sizes <- function() {
  barplot(object.sizes(),
          main="Memory usage by object", ylab="Bytes", xlab="Variable name",
          col=heat.colors(length(object.sizes())))
}


##LISP-like functions
all <- function (x) {
  Reduce(`&`, x)
}

some <- function (x) {
  Reduce(`|`, x)
}

#' 
#' @export
map <- function(x, fn, combine=c, n.cores=1) {
  if (n.cores > 1) {
    registerDoParallel(n.cores)
    foreach(j=x, .combine=combine) %dopar% {
      fn(j)
    }
  } else {
    foreach(j=x, .combine=combine) %do% {
      fn(j)
    }
  }
}



parse.formula <- function (formula, ...)
    # copied from package mosaic
{
  op <- formula[[1]]
  condition <- NULL
  if (length(formula) == 2) {
    rhs <- formula[[2]]
    lhs <- NULL
  }
  else if (length(formula) == 3) {
    rhs <- formula[[3]]
    lhs <- formula[[2]]
  }
  else {
    stop("Invalid formula type.")
  }
  if (inherits(rhs, "call") && rhs[[1]] == "|") {
    condition <- rhs[[3]]
    rhs <- rhs[[2]]
  }
  return(structure(list(op = op, lhs = lhs, rhs = rhs, condition = condition), 
                   class = "parsedFormula"))
}

get.terms <- function(fmla) {
  if(is.null(fmla)) {
    return(c())
  } else if(!(class(fmla) %in% c('call', 'formula'))) {
    return(paste(deparse(fmla), collapse=""))
  } else {
    parsed.formula <- parse.formula(fmla)
    if(!(is.null(parsed.formula$op) | paste(parsed.formula$op) %in% c('+', '~'))) {
      return(paste(deparse(fmla), collapse=""))
    } else {
      lhs.terms <- get.terms(parsed.formula$lhs)
      rhs.terms <- get.terms(parsed.formula$rhs)
      condition.terms <- get.terms(parsed.formula$condition)
      return(c(lhs.terms, rhs.terms, condition.terms))
    }
  }
}


#' @export
#' @aliases as.lm
#' @rdname as.lm
as.lm <- function(object, ...) UseMethod("as.lm")

#' @name as.lm
#' @aliases as.lm.nls
#' @export
#' @author Walmes Zeviani, \email{walmes@@ufr.br}
#' @title Converts nls class objects to lm class
#' @description This is a modified version of Walmes Zeviani, which is
#'     a modified version of the code
#'     available at the \code{nls2} package on github
#'     \url{https://github.com/ggrothendieck/nls2}. This function is no
#'     longer available in the current version of \code{nls2} package.
#' @param object An object of class \code{nls}.
#' @param ... Currently not used.
#' @return This function returns an object of class \code{lm} created
#'     from the \code{nls} object. The linear model design matrix used
#'     is the the partial derivatives of the \code{nls} model function
#'     with relation to the model parameters, i.e. the gradient matrix.
#' @details This function is useful to get the residuals plot for the
#'     fitted nls model. The dependence on the \code{as.proto.list}
#'     function in the \code{proto} package was removed (also this
#'     function no longer exists). Thanks for the original author of
#'     this function G. Grothendieck.
#' @examples
#'
#' # An simple nls fit.
#' n0 <- nls(dist ~ A + B * speed^C,
#'           data = cars,
#'           start = c(A = 0, B = 1, C = 1))
#' summary(n0)
#'
#' # The results.
#' plot(dist ~ speed, data = cars)
#' with(as.list(coef(n0)), {
#'     curve(A + B * speed^C, xname = "speed", add = TRUE)
#' })
#'
#' # Residual analysis.
#' par(mfrow = c(2, 2))
#' plot(as.lm(n0))
#' layout(1)
#'
# if object is an "nls" object then its often used like this:
# predict(as.lm(object), ...) where ... are any predict.lm args
#
# as.lm.nls effectively just does this:
# lm(lhs ~ gradient - 1, offset = fitted(object),
#   list(gradient = object$m$gradient(), lhs = object$m$lhs()))
# so most of the code is just to get the names right.
#
as.lm.nls <- function(object, ...) {
    if (!inherits(object, "nls")) {
        w <- paste(
            "expected object of class nls but got object of class:",
            paste(class(object), collapse = " "))
        warning(w)
    }

    gradient <- object$m$gradient()
    if (is.null(colnames(gradient))) {
        colnames(gradient) <- names(object$m$getPars())
    }

    response.name <- if (length(stats::formula(object)) == 2) {
                         "0"
                     } else {
                         as.character(stats::formula(object)[[2]])
                     }
    lhs <- object$m$lhs()

    L <- data.frame(lhs, gradient)
    names(L)[1] <- response.name

    fo <- sprintf("%s ~ %s - 1",
                  response.name,
                  paste(colnames(gradient), collapse = " + "))
    fo <- stats::as.formula(fo)

    m <- lm(fo, #RB: I changes this 
           #    offset = substitute(fitted(object)),
            data = L)
    
    #RB: Make sure the coefficients are exactly the same. For some reason they are not always the same
    m$coefficients <- object$m$getPars()
    m$call <- object$call
    m$call[[1]] <- quote(lm)
    return(m)
}


## data.table functions


# Anonymize some columns of a data.frame.
# Replaces the original values with strings starting with the column name
# followed by a number.
# df can be a data frame or a data table.
# Returns : A data.frame with the selected columns anonymized.
anonymize.cols <- function(df, cols) {
  df <- data.frame(df)
  for(col in cols) {
    anon.col <- match(df[, col], sample(unique(df[, col])), incomparables =c(NA))
    df[, col] <- ifelse(is.na(anon.col), NA, paste0(col, "_", anon.col))
  }
  df
}


##Statistics functions

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' 
#' @export
cl <- function(fm, cluster, vcov=FALSE, fast=FALSE){
    cluster <- cluster[!(1:length(cluster) %in% fm$na.action)]
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- fm$rank
    dfc <- (M/(M-1))*((N-1)/(N-K))
    if (fast==TRUE) {
      # Alternative (runs faster in some computers, but seems to be wrong sometimes): 
      dt.uj <- data.table(cbind(cluster, estfun(fm)))
      dt.uj <- dt.uj[, lapply(.SD, sum), by=cluster]
      dt.uj[, cluster := NULL]
      uj <- as.matrix(dt.uj)
    } else {
      uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum))
    }
    vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
    if(vcov) {
      vcovCL
    } else{
      coeftest(fm, vcovCL)
      }
}

#' 
#' @export
cl.se <- function(fm, cluster, fast=FALSE){
    cl(fm, cluster, fast=fast)[, "Std. Error"]
}

#' 
#' @export
robust.se <- function(fm) {
    cov.fm <- vcovHC(fm, type = "HC")
    rob.std.err.fm <- sqrt(diag(cov.fm))
    rob.std.err.fm
}

#' 
#' @export
naive.se <- function(fm) {
    sqrt(diag(vcov(fm)))
}

#' 
#' @export
summ.stats <- function(vec.list, na.rm=TRUE) {
    list(
         Var = names(vec.list),
         N = sapply(vec.list, FUN=function(x) sum(!is.na(x))),
         Mean = sapply(vec.list, mean, na.rm=na.rm),
         Median = as.double(sapply(vec.list, median, na.rm=na.rm)),
         S.D. = sapply(vec.list, sd, na.rm=na.rm),
         #S.E. = sapply(vec.list, sd, na.rm=na.rm)/sqrt(sapply(vec.list, FUN=function(x) sum(!is.na(x)))),
         Min = sapply(vec.list, min, na.rm=na.rm),
         Max = sapply(vec.list, max, na.rm=na.rm))
}

#' 
#' @export
summ <- gtools::defmacro(df.name, DOTS, by=list(), out=NULL, expr={
    .varlist = quote(c(...))
    if(length(.varlist) == 1) {
        .varlist=names(df.name)
    } else {
        .varlist=paste(.varlist)[2:length(.varlist)]
    }
    stargazer(data.frame(map(.varlist,
                             function(var) df.name[, summ.stats(.SD), .SDcols=var, by=by],
                             combine=rbind,
                             n.cores=4)),
              summary=FALSE, type="text", out=out, float = FALSE)
})

compress.dt <- function(dt) {
   data.table(as.data.frame(lapply(dt ,function(x) type.convert(as.character(x)))))
}

# Concentration indices

herf.coeff <- function(vec) {
  vec <- vec[!is.na(vec)]
  tot <- sum(vec)
  N <- length(vec)
  H <- sum((vec/tot)^2)
  H.star <- (H-1/N)/(1-1/N)
  H.star
}

gini.coeff <- function(vec) {
  vec <- vec[!is.na(vec)]
  tot <- sum(vec)
  cum.vec <- cumsum(sort(vec))
  unif.dist <- sapply(1:length(vec), function(x) {x/length(vec) * tot})
  sum(unif.dist - cum.vec) / tot / length(vec) * 2
}

theil.coeff <- function(vec) {
  avg <- mean(vec, na.rm=TRUE)
  mean((vec/avg) * log(vec/avg), na.rm=TRUE)
}


# Bootstrap

bootstrap <- function(x, iterations, fn, n.cores=1) {
    set.seed(1979)
    map(1:iterations, function(i) fn(sample(x, replace=TRUE)), n.cores=n.cores)
}

