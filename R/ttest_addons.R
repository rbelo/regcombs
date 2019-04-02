#library(data.table)
#library(formula.tools)
#library(Hmisc)

#install.packages("devtools")
#library("devtools")
# #devtools::install_github("klutometis/roxygen")
#library(roxygen2)
#setwd("~/research")
#create("regcombs")

#source("common.R")


#' Creates an ANOVA table
#' @export
#' @examples
#' my.table <- anova.table(education + age ~ treated | month_id, dt.data, out="tables/anova_table_1.tex")
anova.table <- function(fmla, data, ...) {
  dep.vars <- get.terms(parse.formula(fmla)$lhs)
  indep.vars <- get.terms(parse.formula(fmla)$rhs)
  control.vars <- get.terms(parse.formula(fmla)$condition)
  if(is.null(control.vars)) {
    control.vars <- c()
  }
  indep.var <- indep.vars[1] # use only one variable for testing
  dt.result <-
    map(1:length(dep.vars),
        function(i) {
          dep.var <- dep.vars[i]
          dt.tmp <- data[, list(var= dep.var,
                                " avg"= mean(get(dep.var), na.rm=TRUE),
                                " sd" = sd(get(dep.var) , na.rm=TRUE),
                                n  = as.double(sum(!is.na(get(dep.var))))),
                         by=c(control.vars, indep.var)]
          dt.tmp <- melt(dt.tmp,
                         measure.vars=c(" avg", " sd", "n"),
                         variable.name = "stat")
          dt.tmp <- dcast.data.table(dt.tmp, as.formula(paste(paste(c("var","stat", control.vars), collapse=" + ")," ~ ", indep.var)),
                                     value.var=c("value"))
          dt.tmp <- merge(dt.tmp,
                          data[, c(var=dep.var, stat=c(" avg"),
                                   summary(aov(eval(as.formula(paste(dep.var, " ~ ", indep.var)))))[[1]][1,c(4,5)]),
                               by=control.vars],
                          by=c("var", control.vars, "stat"), all.x=TRUE)
          dt.tmp
        },
        combine=rbind)

  dt.result[, n_row := 1:.N, by=c("var", control.vars)]
  dt.result[n_row != 1, c("var", control.vars) := NA]
  dt.result[, n_row := 1:.N, by=c("var")]
  dt.result[n_row != 1, c("var") := NA]
  dt.result[, n_row := NULL]

  stargazer(dt.result, summary=FALSE,
            type="text",
            float=FALSE,
            digits=2,
            rownames=FALSE, ...)

  return(dt.result)
}

#' Create a ttest table
#' @export
#' @examples
#' my.table <- ttest.table(education + age ~ treated | month_id, dt.data, out="tables/ttest_table_1.tex")
ttest.table <- function(fmla, data, add.stats=c(), ...) {

  dep.vars <- get.terms(parse.formula(fmla)$lhs)
  indep.vars <- get.terms(parse.formula(fmla)$rhs)
  control.vars <- get.terms(parse.formula(fmla)$condition)
  if(is.null(control.vars)) {
    control.vars <- c()
  }
  indep.var <- indep.vars[1] # use only one variable for testing

  indep.var.values <- sort(data[, unique(get(indep.var))])
  indep.var.comparisons <- combn(indep.var.values, 2)
  indep.var.comparisons.list <- lapply(1:ncol(indep.var.comparisons), function(x) indep.var.comparisons[, x])

  dt.result <-
    map(1:length(dep.vars),
        function(i) {
          dep.var <- dep.vars[i]
          map(indep.var.comparisons.list,
              function(comparison) {
                data[, list(var_order = i,
                            var=dep.var,
                            test = paste(comparison[1], "vs.", comparison[2]),
                            mean_1=mean(get(dep.var)[get(indep.var)== comparison[1]], na.rm=TRUE),
                            sd_1=sd(get(dep.var)[get(indep.var)== comparison[1]], na.rm=TRUE),
                            n_1 = as.double(sum(!is.na(get(dep.var)[get(indep.var)==comparison[1]]))),
                            min_1 = min(get(dep.var)[get(indep.var)== comparison[1]], na.rm=TRUE),
                            max_1 = max(get(dep.var)[get(indep.var)== comparison[1]], na.rm=TRUE),
                            mean_2=mean(get(dep.var)[get(indep.var)==comparison[2]], na.rm=TRUE),
                            sd_2=sd(get(dep.var)[get(indep.var)== comparison[2]], na.rm=TRUE),
                            n_2 = as.double(sum(!is.na(get(dep.var)[get(indep.var)==comparison[2]]))),
                            min_2 = min(get(dep.var)[get(indep.var)== comparison[2]], na.rm=TRUE),
                            max_2 = max(get(dep.var)[get(indep.var)== comparison[2]], na.rm=TRUE),
                            diff = mean(get(dep.var)[get(indep.var)==comparison[1]], na.rm=TRUE) -
                              mean(get(dep.var)[get(indep.var)== comparison[2]], na.rm=TRUE),
                            t_stat=t.test(get(dep.var)[get(indep.var)==comparison[1]],
                                          get(dep.var)[get(indep.var)==comparison[2]])["statistic"][[1]],
                            pval= t.test(get(dep.var)[get(indep.var)==comparison[1]],
                                         get(dep.var)[get(indep.var)==comparison[2]])["p.value"][[1]]),
                     by=control.vars]
              }, combine=rbind)
        },
        combine=rbind)

  if (length(add.stats) == 0) {
    stats <- c("avg", "SD", "n")
  } else {
    stats <- c("avg", "SD", add.stats, "n")
  }

  if (length(indep.var.values) == 2) {
    test.var <- c()
  } else {
    test.var <- "test"
  }

  dt.tmp.1 <- dt.result[, c("var", control.vars, test.var, paste0(c("mean", "sd", add.stats, "n"), "_1")), with=FALSE]
  setnames(dt.tmp.1, paste0(c("mean", "sd", add.stats, "n"), "_1"), capitalize(stats))
  dt.tmp.2 <- dt.result[, c("var", control.vars, test.var, paste0(c("mean", "sd", add.stats, "n"), "_2")), with=FALSE]
  setnames(dt.tmp.2, paste0(c("mean", "sd", add.stats, "n"), "_2"), capitalize(stats))

  dt.result <-
    merge(
      merge(
        melt(dt.tmp.1,
             measure.vars = capitalize(stats),
             variable.name="stat",
             value.name = "group_1")[order(get(c("var", control.vars, test.var)))],
        melt(dt.tmp.2,
             measure.vars = capitalize(stats),
             variable.name="stat",
             value.name = "group_2")[order(get(c("var", control.vars, test.var)))],
        by=c(control.vars, "var", test.var, "stat")),
      dt.result[, c("var_order", "var", control.vars, test.var, "diff", "t_stat", "pval"), with=FALSE],
      by=c("var", control.vars, test.var), all.x=TRUE)[order(var_order, get(c("var", control.vars, test.var)))]

  dt.result[, var_order := NULL]
  dt.result[stat != "Avg",
            `:=`(test = NA,
                 diff = NA,
                 t_stat = NA,
                 pval = NA)]
  dt.result[, n_row := 1:.N, by=c("var", control.vars)]
  dt.result[n_row != 1, c("var", control.vars) := NA]
  dt.result[, n_row := 1:.N, by=c("var")]
  dt.result[n_row != 1, c("var") := NA]
  dt.result[, n_row := NULL]
  setnames(dt.result,
           c("var", "stat", "diff", "t_stat", "pval"),
           c("Variable", "Statistic", "Difference", "t-Value", "p-Value"))
    if (length(indep.var.values) == 2) {
    setnames(dt.result, c("group_1", "group_2"), paste(indep.var.values))
    dt.result[, test := NULL]
  }

  stargazer(dt.result, summary=FALSE,
            type="text",
            float=FALSE,
            digits=2,
            rownames = FALSE,
            ...)

  return(dt.result)
}

