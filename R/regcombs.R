#library(data.table)
#library(foreign)
#library(formula.tools)
#library(gtools)
#library(lmtest)
#library(sandwich)
#library(stargazer)

add.to.models <- function(dt.models, var.name, fmla) {
  dt.models <- copy(dt.models)
  dt.models[, k := 1]
  if ("formula" %in% class(fmla)) {
    to.repeat <- get.terms(parse.formula(fmla)$lhs)
    to.crossjoin <- get.terms(parse.formula(fmla)$rhs)
    to.recycle <- get.terms(parse.formula(fmla)$condition)
    if (length(to.recycle) > 0) {
      dt.models <- cbind(dt.models,
                         data.table(tmp_var = to.recycle))
    } else if(length(to.repeat) > 0) {
      dt.models <- merge(dt.models,
                         data.table(dep_var = dt.models[, unique(dep_var)],
                                    tmp_var = to.repeat),
                         by="dep_var", all.x = TRUE)
    } else { # length(to.crossjoin) is assumed greater than 0
      dt.models <- merge(dt.models,
                         data.table(k = 1,
                                    tmp_var = to.crossjoin),
                         by="k", allow.cartesian = TRUE)
    }
  } else {
    dt.models <-
      merge(dt.models,
            data.table(k = 1,
                       tmp_var = paste(deparse(substitute(fmla)), collapse="")),
            by="k", all.x = TRUE)
  }
  dt.models[, substitute(var.name) := tmp_var]
  dt.models[, k := NULL]
  dt.models[, tmp_var := NULL]
  copy(dt.models)
}


reg.combs.models <- function(fmla, data, reg.fn = ~ felm, reg.params = ~ 0,
                             fe = ~ 0, iv = ~ 0, cl = ~ 0, w = ~ 0) {
  dep.vars <- get.terms(parse.formula(fmla)$lhs)
  indep.vars <- get.terms(parse.formula(fmla)$rhs)
  controls <- get.terms(parse.formula(fmla)$condition)
  if (length(controls) == 0) {
    controls <- "1"
  }
  dt.models <- data.table(expand.grid(indep_vars = indep.vars, # indep_vars first because expand.grid groups the second variable instead of the first.
                                      dep_var = dep.vars))[, list(dep_var, indep_vars)]
  dt.models <- add.to.models(dt.models, "controls",
                             as.formula(paste( "~", paste0(controls, collapse=" + "))))
  dt.models[controls == "1", controls := ""]
  dt.models <- add.to.models(dt.models, "reg_fn", reg.fn)
  dt.models <- add.to.models(dt.models, "reg_params", reg.params)
  dt.models <- add.to.models(dt.models, "data", data)
  dt.models <- add.to.models(dt.models, "fe", fe)
  dt.models <- add.to.models(dt.models, "iv", iv)
  dt.models <- add.to.models(dt.models, "cl", cl)
  dt.models <- add.to.models(dt.models, "w", w)
  copy(dt.models)
  }

#' Reg Combs
#'
#' This function allows you to perform regression combinations.
#' @import data.table
#' @import doParallel
#' @import foreach
#' @import foreign
#' @import formula.tools
#' @import gtools
#' @import Hmisc
#' @import iterators
#' @import lmtest
#' @import parallel
#' @import sandwich
#' @import stargazer
#' @param fmla a formula
#' @param data data
#' @keywords regression combinations
#' @export
#' @examples
#' reg.combs(a1 + a2 ~ b1 + b2 + (b1 + b2) | c1 + (c1 + c2),
#'          fe = ~ (f1 + f2),
#'          iv = ~
#'            (e1 | e2 ~ i1 + i2) +
#'            (e3 | e4 ~ i3 + i4),
#'          cl = ~ cl1 + cl2,
#'          w = ~ w1 + w2,
#'          data = dt[1] + dt[2] ~ .,
#'          reg.fn = ~ felm + logit,
#'          test = TRUE)
reg.combs <- function(fmla, data, reg.fn = ~ felm, reg.params = ~ 0,
                      fe = ~ 0, iv = ~ 0, cl = ~ 0, w = ~ 0,
                      robust = FALSE,
                      omit.stat=c("f", "ser"),
                      n.cores=1,
                      model.summaries = TRUE,
                      test = FALSE,
                      ...) {
  dt.models <- reg.combs.models(fmla = fmla,
                                data = data,
                                reg.fn = reg.fn,
                                reg.params = reg.params,
                                fe = fe, iv = iv, cl = cl, w = w)
  commands <- map(1:nrow(dt.models), function(iter) {
    weights.string <- ifelse(dt.models[iter, w] == "0", "",
                             paste0(", weights = ",
                                    dt.models[iter, data], "[, ",
                                    dt.models[iter, w], "]"))
    if (dt.models[iter, reg_fn]  == "felm") {
      reg.fn.string <- "felm"
      reg.fn.params.string <- paste0(" | ", dt.models[iter, fe],
                                     " | ", dt.models[iter, iv],
                                     " | ", dt.models[iter, cl],
                                     weights.string)
    } else if (dt.models[iter, reg_fn] == "logit") {
      reg.fn.string <- "glm"
      reg.fn.params.string <-
        paste0(ifelse(dt.models[iter, fe] == "0", "", paste0(" + ", dt.models[iter, fe])),
               weights.string, ", family=binomial(link='logit')")
    } else if (dt.models[iter, reg_fn] == "probit") {
      reg.fn.string <- "glm"
      reg.fn.params.string <-
        paste0(ifelse(dt.models[iter, fe] == "0", "", paste0(" + ", dt.models[iter, fe])),
               weights.string, ", family=binomial(link='probit')")
    } else {
      reg.fn.string <- dt.models[iter, reg_fn]
      reg.fn.params.string <- ifelse(dt.models[iter, reg_params] == "0", "",
                                     paste(",", dt.models[iter, reg_params]))
    }
    formula.string <- paste0(dt.models[iter, dep_var], " ~ ",
                             dt.models[iter, indep_vars],
                             ifelse(dt.models[iter, controls] == "", "", " + "),
                             dt.models[iter, controls])
    command <- paste0(reg.fn.string, "(formula = ", formula.string,
                      reg.fn.params.string, ", data = ", dt.models[iter, data],
                      ")")
    command
  }, combine=c)
  if (test) {
    return(commands)
  }
  models <- map(1:length(commands),
                function(i) {
                  message(paste("Time:", Sys.time()))
                  message(paste("Model:", commands[i]))
                  list(eval(parse(text=commands[i])))
                },
                combine = c,
                n.cores = n.cores)
  ## Make sure that nls models are transformed to lm objects so that stargazer can display them
  models <- lapply(models, function(model) {
    if (class(model) == "nls") {
      as.lm(model)
    } else {
      model
    }
  })
  if (dt.models[reg_fn != "felm" & cl != "0", .N] > 0) {
    ## Calculate SE using cl.se
    models.se <-
      map(1:length(commands), function(iter) {
        message(paste("Time:", Sys.time()))
        if (dt.models[iter, reg_fn] == "felm" | dt.models[iter, cl] == "0") {
          cl.cmd <- paste0("naive.se(models[[", iter, "]])")
        } else {
        cl.cmd <- paste0("cl.se(models[[", iter, "]],",
                         dt.models[iter, data], "[, ", dt.models[iter, cl], "])")
        }
        message(paste("Model SE:", cl.cmd))
        list(eval(parse(text=cl.cmd)))
      }, combine = c, n.cores = n.cores)
    stargazer(models, se=models.se, type="text", omit.stat=omit.stat, ...)
    if(model.summaries == TRUE) {
      models <- lapply(models, summary)
    }
    list(models=models,
         models.se=models.se)
  } else if (robust == TRUE ) {
    ## Calculate SE using robust.se
    models.se <-
      map(1:length(commands), function(iter) {
        message(paste("Time:", Sys.time()))
        robust.cmd <- paste0("robust.se(models[[", iter, "]])")
        message(paste("Model SE:", robust.cmd))
        list(eval(parse(text=robust.cmd)))
      }, combine = c, n.cores = n.cores)
    stargazer(models, se=models.se, type="text", omit.stat=omit.stat, ...)
    if(model.summaries == TRUE) {
      models <- lapply(models, summary)
    }
    list(models=models,
         models.se=models.se)
  } else {
    stargazer(models, type="text", omit.stat=omit.stat, ...)
    if(model.summaries == TRUE) {
      models <- lapply(models, summary)
    }
    models
  }
}


