% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/robusta.R
\name{cluster}
\alias{cluster}
\title{Setting cluster robust standard error to be reported in regression summary}
\usage{
cluster(reg, clustervar)
}
\arguments{
\item{reg}{regression objects, for now only works for lm, glm}

\item{clustervar}{the clustering variable in the data frame that is used for regression. The name of data frame need not to be mentioned since it will be read from the object's call.}
}
\description{
Summary command will use the cluster variable attached by this command to the regression object to calculate cluster robust standard error
}
\examples{
lmo = lm(mpg ~ carb + cyl, mtcars)
texreg::screenreg(list(lmo, cluster(lmo, am)))
}

