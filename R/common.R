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

cl.se <- function(fm, cluster, fast=FALSE){
    cl(fm, cluster, fast=fast)[, "Std. Error"]
}

robust.se <- function(fm) {
    cov.fm <- vcovHC(fm, type = "HC")
    rob.std.err.fm <- sqrt(diag(cov.fm))
    rob.std.err.fm
}

naive.se <- function(fm) {
    sqrt(diag(vcov(fm)))
}

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
