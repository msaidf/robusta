#' Specifying the type of robust standard error to be reported in regression summary
#' 
#' Summary command will identify which type of standard error to use by reading robust attribute attached by this command
#' @param reg regression objects, for now only works for lm, glm, and plm objects 
#' @param se type of robust standard error. The default is HC1 to make it easy to compare it with stata robust standard error. The default heteroskedasticity consistent (HC) covariance matrix in sandwich package is instead HC3, that is less affected by outlier in small sample. All types of robust standard error can be referred by the name of author who propose it. 'hc1' is equivalent to 'Hinkley', 'hc0' is equivalent to 'White', 'hc2' is equivalent to 'Horn', 'hc3' is equivalent to 'Efron', 'hc4' is equivalent to 'Cribari'. For data containing serial correlation, the available heteroscedasticity and autocorrelation consistent (HAC) standard error is 'NeweyWest', 'Andrew' that is equivalent to 'hac', and 'Lumley'. All HC types and 'arrellano' method are available for plm object in panel data analysis.  
#' @examples 
#' lmo = lm(mpg ~ carb + cyl, mtcars) 
#' texreg::screenreg(list(lmo, robust(lmo), robust(lmo, 'hc3'), robust(lmo, 'White')))
#' @export 
robust = function(reg, se = c('Hinkley', 'hc1', 'stata', 'White', 'hc0','HornDuncan', 'hc2', 'Efron', 'hc3', 'CribariNeto', 'hc4', 'NeweyWest', 'Andrews', 'hac', 'Lumley', 'Arellano'))
{
  if(length(se) > 1) se = 'hc1'
  switch(se 
	, 'Hinkley' = {attributes(reg)$robust <- 'HC1'}
	, 'hc1'  = {attributes(reg)$robust <- 'HC1'}
	, 'White'  = {attributes(reg)$robust <- 'HC0'}
	, 'hc0'  = {attributes(reg)$robust <- 'HC0'}
	, 'Horn'  = {attributes(reg)$robust <- 'HC2'}
	, 'hc2'  = {attributes(reg)$robust <- 'HC1'}
	, 'Efron'  = {attributes(reg)$robust <- 'HC3'}
	, 'hc3'  = {attributes(reg)$robust <- 'HC3'}
	, 'Cribari'  = {attributes(reg)$robust <- 'HC4'}
	, 'hc4'  = {attributes(reg)$robust <- 'HC4'}
	, 'NeweyWest'  = {attributes(reg)$robust <- 'NeweyWest'}
	, 'Andrews'  = {attributes(reg)$robust <- 'kernHAC'}
	, 'hac'  = {attributes(reg)$robust <- 'kernHAC'}
	, 'Lumley'  = {attributes(reg)$robust <- 'weave'}
	, 'Arellano'  = {attributes(reg)$robust <- 'Arellano'}
  )
	reg
}

#' Setting cluster robust standard error to be reported in regression summary
#' 
#' Summary command will use the cluster variable attached by this command to the regression object to calculate cluster robust standard error
#' @param reg regression objects, for now only works for lm, glm
#' @param clustervar the clustering variable in the data frame that is used for regression. The name of data frame need not to be mentioned since it will be read from the object's call.
#' @examples 
#' lmo = lm(mpg ~ carb + cyl, mtcars) 
#' texreg::screenreg(list(lmo, cluster(lmo, am)))
#' @export 
cluster = function(reg, clustervar)
{
  data = as.character(reg$call)[3]
  data = eval(parse(text = data), envir = .GlobalEnv)
  reg$cluster <- data[[deparse(substitute(clustervar))]]
  reg
}

coeftest_cluster <- function(model) {
	cluster <- model$cluster
	M <- length(unique(cluster))
	N <- length(cluster)
	K <- model$rank
	dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
	uj <- apply(sandwich::estfun(model), 2, function(x) tapply(x, cluster, sum))
	rcse.cov <- sandwich::sandwich(model, meat = crossprod(uj)/N) * dfc
	lmtest::coeftest.default(model, rcse.cov)
}

#' @export 
summary.lm <- function(reg) {
	sumreg = stats::summary.lm(reg)
	if(!is.null(attributes(reg)$robust)) {
		switch(attributes(reg)$robust  
		, 'HC1' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC1')))}
		, 'HC0' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC0')))}
		, 'HC2' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC2')))}
		, 'HC3' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC3')))}
		, 'HC4' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC4')))}
		, 'NeweyWest' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::NeweyWest))}
		, 'kernHAC' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::kernHAC))}
		, 'weave' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::weave))}
		)
	}
	if(!is.null(reg$cluster)) sumreg$coefficients = coeftest_cluster(reg)
	sumreg
}

#' @export 
summary.glm <- function(reg) {
	sumreg = stats::summary.lm(reg)
	if(!is.null(attributes(reg)$robust)) {
		switch(attributes(reg)$robust  
		, 'HC1' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC1')))}
		, 'HC0' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC0')))}
		, 'HC2' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC2')))}
		, 'HC3' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC3')))}
		, 'HC4' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::vcovHC(reg, type = 'HC4')))}
		, 'NeweyWest' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::NeweyWest))}
		, 'kernHAC' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::kernHAC))}
		, 'weave' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = sandwich::weave))}
		)
	}
	if(!is.null(reg$cluster)) sumreg$coefficients = coeftest_cluster(reg)
	sumreg
}

summary.plm <- function(reg) {
	sumreg = summary(reg)
	if(!is.null(attributes(reg)$robust)) {
		switch(attributes(reg)$robust  
		, 'HC1' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = plm::vcovBK(reg, type = 'HC1')))}
		, 'HC0' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = plm::vcovBK(reg, type = 'HC0')))}
		, 'HC2' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = plm::vcovBK(reg, type = 'HC2')))}
		, 'HC3' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = plm::vcovBK(reg, type = 'HC3')))}
		, 'HC4' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = plm::vcovBK(reg, type = 'HC4')))}
		, 'Arellano' = {sumreg$coefficients = as.matrix(lmtest::coeftest.default(reg, vcov = plm::vcovHC(reg, method = 'arellano')))}
		)
	sumreg
	}
}

