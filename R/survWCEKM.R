survWCEKM <- function(data) {

  data_kmtemp <- data %>%
      group_by(group)  %>%
      mutate(
        n.censor.all = ave(n.censor, time, FUN = sum)
      ) %>%
      distinct(time, .keep_all = TRUE) %>%
      arrange(group, time,id) %>% ungroup()


  ntime = nrow(data_kmtemp)
  group_lv = levels(data_kmtemp$group)
  ngroup <- length(group_lv)
  sort_time = drop(data_kmtemp$time)
  sort_group = drop(data_kmtemp$group)
  sort_n.event = data_kmtemp$n.event
  sort_n.risk = data_kmtemp$n.risk
  n.censor.all = data_kmtemp$n.censor.all


  # Call the external C function to perform the log-rank test
  xxx <- .C(CsurvWCEKM,
           ntime = as.integer(ntime),
           ngroup = as.integer(ngroup),
           time = as.double(sort_time),       # time
           group = as.integer(sort_group),     # group
           n.event = as.double(sort_n.event),     # n.event
           n.risk = as.double(sort_n.risk),   # n.risk
           survival = double(ntime),
           std.err = double(ntime),
           cum.haz = double(ntime),
           std.chaz = double(ntime)
           )


  # Return the results as a dataframe
  data.frame(
    time = xxx$time,
    group = xxx$group,
    n.event = xxx$n.event,
    n.censor = n.censor.all,
    n.risk = xxx$n.risk,
    survival = xxx$survival,
    std.err = xxx$std.err,
    cumhaz = xxx$cum.haz,
    std.chaz = xxx$std.chaz
    )
}
