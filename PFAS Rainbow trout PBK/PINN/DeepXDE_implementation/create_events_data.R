events <- create.events(params)
events = events$data[events$data$var=='M_lumen',]

events = events[, c('time', 'value')]

build_up_events <- events[1,]
for (i in 2:dim(events)[1]) {
  
  time_grid <- seq(events[i-1,1] + 5/60, events[i,1] - 5/60, length.out = 50)
  intermediate_grid = cbind(time_grid, rep(0, length(time_grid)))
  colnames(intermediate_grid) = colnames(build_up_events)
  
  build_up_events <- rbind(build_up_events, intermediate_grid, events[i,])
}

extra_time_grid <- seq(build_up_events[dim(build_up_events)[1],1] + 5/60, 650, length.out=20)
extra_grid <- cbind(extra_time_grid, rep(0, length(extra_time_grid)))

colnames(extra_grid) <- colnames(build_up_events)
build_up_events <- rbind(build_up_events, extra_grid)

write.csv(build_up_events, file = '/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/events.csv',
          row.names = F)

#===============================================================================
# Just write the events to a csv file
events = events$data[events$data$var=='M_lumen',]
events = events[, c('time', 'value')]

write.csv(events, file = '/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/events_simplified.csv',
          row.names = F)
