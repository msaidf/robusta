#' @export 
robust = function(reg, se = c('Hinkley', 'hc1', 'stata', 'White', 'hc0','HornDuncan', 'hc2', 'Efron', 'hc3', 'CribariNeto', 'hc4', 'NeweyWest', 'Andrews', 'hac', 'LumleyHeagerty', 'Arellano', 'bk'))
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
	prcse.cov <- sandwich::sandwich(model, meat = crossprod(uj)/N) * dfc
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

