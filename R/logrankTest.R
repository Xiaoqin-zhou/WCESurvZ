#' logrankTest
#'
#' Perform Log-Rank Test for Survival Data
#'
#' This function performs a log-rank test to compare survival curves between different groups.
#'
#' @param data
#'
#' @return NA
#' @export
#'
#' @examples NA
logrankTest <- function(data){

  data_temp = data %>% select(id,event,time,group,status,weight,n.risk_1,n.risk_2)

  nsize = nrow(data_temp)
  ngroup = length(unique(data_temp$group))
  ord_time = drop(data_temp$time)
  ord_status = drop(data_temp$status)
  ord_group = drop(data_temp$group)
  ord_weight = drop(data_temp$weight)
  ord_n.risk_1 = drop(data_temp$n.risk_1)
  ord_n.risk_2 = drop(data_temp$n.risk_2)
  ord_n.risk_all = ord_n.risk_1+ord_n.risk_2



      xx <- .C(Cwcelogrank, as.integer(nsize),
		   as.integer(ngroup),
		   as.double(ord_time),  #time
		   as.integer(ord_status),  #status
		   as.integer(ord_group),  #group
		   as.double(ord_weight),  #weight
		   as.double(ord_n.risk_1),  #n.risk_1
		   as.double(ord_n.risk_2),  #n.risk_2
		   as.double(ord_n.risk_all),  #n.risk_all
		   observed = double(ngroup),
		   expected = double(ngroup),
		   tempted = double(ngroup),
		   var.e    = double(1),
		   double(ngroup))

        otmp= xx$observed
        etmp= xx$expected
      	df   <- (etmp >0)

  	    temp2 <- ((otmp - etmp)[df])[-1]
  	    vv <- xx$var.e
  	    chi <- sum(solve(vv, temp2) * temp2)
        df <- (sum(1*(etmp>0))) -1
        pvalue= pchisq(chi, df, lower.tail=FALSE)

    list(expected = etmp,
         observed = otmp,
  			 chi = chi,
  			 pvalue  =pvalue)


}
