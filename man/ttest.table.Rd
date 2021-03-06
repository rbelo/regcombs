% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ttest_addons.R
\name{ttest.table}
\alias{ttest.table}
\title{Create a ttest table}
\usage{
ttest.table(fmla, data, add.stats = c(), ...)
}
\description{
Create a ttest table
}
\examples{
my.table <- ttest.table(education + age ~ treated | month_id, dt.data, out="tables/ttest_table_1.tex")
}
