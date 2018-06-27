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
  if (class(fmla) == "formula") {
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
    } else {
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


reg.combs.models <- function(fmla, data, reg.fn,
                           fe, iv, cl, w,
                           controls.alone,
                           include.all) {
  dep.vars <- get.terms(parse.formula(fmla)$lhs)
  indep.vars <- get.terms(parse.formula(fmla)$rhs)
  if (include.all == TRUE) {
    indep.vars.comb <- unique(c(indep.vars, paste(indep.vars, collapse=" + ")))
  } else {
    indep.vars.comb <- indep.vars
  }
  if (controls.alone == TRUE) {
    indep.vars.comb <- c(1, indep.vars.comb)
  }
  controls <- get.terms(parse.formula(fmla)$condition)
  dt.models <- CJ(dep_var = dep.vars,
                  indep_vars = indep.vars.comb,
                  controls = paste0(controls, collapse=" + "))
  dt.models <- add.to.models(dt.models, "reg_fn", reg.fn)
  dt.models <- add.to.models(dt.models, "data", data)
  dt.models <- add.to.models(dt.models, "fe", fe)
  dt.models <- add.to.models(dt.models, "iv", iv)
  dt.models <- add.to.models(dt.models, "cl", cl)
  dt.models <- add.to.models(dt.models, "w", w)
  copy(dt.models)
  }

#' Reg Combs
#'
#' This function allows you to perform regression combinations
#' @param fmla a formula
#' @param data data
#' @keywords cats
#' @export
#' @examples
#' reg.combs(a1 + a2 ~ b1 + b2 | c1 + c2 ,
#'          fe = ~ (f1 + f2),
#'          iv = ~
#'            (e1 | e2 ~ i1 + i2) +
#'            (e3 | e4 ~ i3 + i4),
#'          cl = ~ (cl1 + cl2),
#'          w = ~ w1 + w2,
#'          data = dt[1] + dt[2] ~ .,
#'          reg.fn = ~ felm + logit,
#'          include.all = TRUE)
reg.combs <- function(fmla, data, reg.fn="felm",
                      fe = ~ 0, iv = ~ 0, cl = ~ 0, w = ~ 0,
                      controls.alone=FALSE,
                      include.all=FALSE,
                      omit.stat=c("f", "ser"),
                      n.cores=1,
                      model.summaries = TRUE,
                      test = FALSE,
                      ...) {
  dt.models <- reg.combs.models(fmla = fmla,
                                data = data,
                                reg.fn = reg.fn,
                                fe = fe, iv = iv, cl = cl, w = w,
                                controls.alone = controls.alone,
                                include.all = include.all)
  commands <- map(1:nrow(dt.models), function(iter) {
    weights.string <- ifelse(w == "0", "",
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
      reg.fn.string <- paste(deparse(substitute(reg.fn)), collapse="")
      reg.fn.params.string <- ""
    }
    formula.string <- paste0(dt.models[iter, dep_var], " ~ ",
                             dt.models[iter, indep_vars], " + ",
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
                  print(paste("Time:", Sys.time()))
                  print(paste("Model:", commands[i]))
                  list(eval(parse(text=commands[i])))
                },
                combine = c,
                n.cores = n.cores)
  if (dt.models[reg_fn != "felm" & cl != "0", .N] > 0) {
    ## Calculate SE using cl.se
    models.se <-
      map(1:length(commands), function(iter) {
        print(paste("Time:", Sys.time()))
        if (dt.models[iter, reg_fn] == "felm" | dt.models[iter, cl] == "0") {
          cl.cmd <- paste0("naive.se(models[[", iter, "]])")
        } else {
        cl.cmd <- paste0("cl.se(models[[", iter, "]],",
                         dt.models[iter, data], "[, ", dt.models[iter, cl], "])")
        }
        print(paste("Model SE:", cl.cmd))
        list(eval(parse(text=cl.cmd)))
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

reg.combs.old <- function(fmla, data, reg.fn="felm",
                           fe = 0, iv = 0, cl = 0, w = 0,
                           controls.alone=FALSE,
                           include.all=FALSE,
                           omit.stat=c("f", "ser"),
                           n.cores=1,
                           model.summaries = TRUE,
                           ...) {

  data.string <- paste(deparse(substitute(data)), collapse="")
  dep.vars <- get.terms(parse.formula(fmla)$lhs)
  indep.vars <- get.terms(parse.formula(fmla)$rhs)
  fe.vars        <- ifelse(fe == 0, 0,
                           paste(get.terms(parse.formula(fe)$rhs), collapse=" + "))
  iv.fmla.string <- ifelse(iv == 0, 0, paste0("(", iv, ")"))
  cl.vars        <- ifelse(cl == 0, 0,
                           paste(get.terms(parse.formula(cl)$rhs), collapse=" + "))
  weights.string <- ifelse(w == 0, "",
                           paste0(", weights = ", data.string, "[, ", paste(get.terms(parse.formula(w)$rhs), collapse=""), "]"))

  ## Prepare additional parameter depending on the regression
  if (class(reg.fn) == "character") {
    if (reg.fn  == "felm") {
      reg.fn.string <- "felm"
      reg.fn.params.string <- paste0(" | ", paste0(fe.vars, collapse=" + "),
                                     " | ", iv.fmla.string,
                                     " | ", paste0(cl.vars, collapse=" + "),
                                     weights.string)
    } else if (reg.fn == "logit") {
      reg.fn.string <- "glm"
      reg.fn.params.string <- paste0(ifelse(fe.vars == 0, "", paste0(" + ", fe.vars)),
                                     weights.string, ", family=binomial(link='logit')")
    } else if (reg.fn == "probit") {
      reg.fn.string <- "glm"
      reg.fn.params.string <- paste0(ifelse(fe.vars == 0, "", paste0(" + ", fe.vars)),
                                     weights.string, ", family=binomial(link='probit')")
    }
  } else {
    reg.fn.string <- paste(deparse(substitute(reg.fn)), collapse="")
    reg.fn.params.string <- ""
  }

  if (include.all == TRUE) {
    indep.vars.comb <- unique(c(indep.vars, paste(indep.vars, collapse=" + ")))
  } else {
    indep.vars.comb <- indep.vars
  }
  if (controls.alone == TRUE) {
    indep.vars.comb <- c(1, indep.vars.comb)
  }
  controls <- get.terms(parse.formula(fmla)$condition)
  formulas <- map(dep.vars, function(dep.var) {
    map(indep.vars.comb, function(indep.var.comb) {
      iter.fmla <- as.formula(paste0(dep.var, " ~ ", indep.var.comb, ifelse(is.null(controls), " ", " + "), paste0(controls, collapse=" + ")))
      list(iter.fmla)}, combine=c)}, combine=c)
  models <-
    map(1:length(formulas), function(i) {
      print(paste("Time:", Sys.time()))
      reg.cmd <- paste0(reg.fn.string, "(formula = ", formulas[i],
                          reg.fn.params.string, ", data = ", data.string,
                            ")")
      print(paste("Model:", reg.cmd))
      list(eval(parse(text=reg.cmd)))
    }, combine = c, n.cores = n.cores)
  if (!(reg.fn.string == "felm") & cl.vars != 0 ) {
    ## Calculate SE using cl.se
    models.se <-
      map(1:length(formulas), function(i) {
        print(paste("Time:", Sys.time()))
        cl.cmd <- paste0("cl.se(models[[", i, "]],",
                         data.string, "[, ", cl.vars[1], "])")
        print(paste("Model SE:", cl.cmd))
        list(eval(parse(text=cl.cmd)))
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

reg.combs.felm <- reg.combs.old

