if (!isNamespace(topenv())) {
    source('R/4domain_reactor.R')
    source('R/ANALOG_GATE_LIB.R')
    source('R/analysis.R')
    source('R/crn_reactor.R')
    source('R/dsd.R')
    source('R/ELECTRO_LIB.R')
    source('R/GATE_LIB.R')
    source('R/io.R')
    source('R/neuron_hjelmfelt.R')
    source('R/parser.R')
    source('R/util_functions.R')
}

#library(binaryLogic)

# Auxiliary function to join strings
jn <- function(...) { paste(..., sep = '') }

make_rgbyclock <- function(name) {
  gate <- list(
    name      = name,
    species   = c(),
    reactions = list(),
    ki        = list(),
    ci        = list()
  )
  
  gate$name      = name
  gate$species   = list( r = jn(name,'_r'), value0 = jn(name,'_value0'),
                         g = jn(name,'_g'), G = jn(name,'_G'),
                         b = jn(name,'_b'), value1 = jn(name,'_value1'),
                         y = jn(name,'_y'), Y = jn(name,'_Y')
  )
  
  gate$ci        = list(0, 1,
                        0, 0,
                        0, 0,
                        0, 0)
  
  gate$reactions = list(jn('0 -> ', gate$species$r),
                        jn(gate$species$value0, ' + ', gate$species$r, ' -> ', gate$species$value0),
                        jn(gate$species$value0, ' + ', gate$species$y, ' -> ', gate$species$G, ' + ', gate$species$y),
                        jn(gate$species$value0, ' + 2', gate$species$G, ' -> 3', gate$species$G),
                        
                        jn('0 -> ', gate$species$g),
                        jn(gate$species$G, ' + ', gate$species$g, ' -> ', gate$species$G),
                        jn(gate$species$G, ' + ', gate$species$r, ' -> ', gate$species$value1, ' + ', gate$species$r),
                        jn(gate$species$G, ' + 2', gate$species$value1, ' -> 3', gate$species$value1),
                        
                        jn('0 -> ', gate$species$b),
                        jn(gate$species$value1, ' + ', gate$species$b, ' -> ', gate$species$value1),
                        jn(gate$species$value1, ' + ', gate$species$g, ' -> ', gate$species$Y, ' + ', gate$species$g),
                        jn(gate$species$value1, ' + 2', gate$species$Y, ' -> 3', gate$species$Y),
                        
                        jn('0 -> ', gate$species$y),
                        jn(gate$species$Y, ' + ', gate$species$y, ' -> ', gate$species$Y),
                        jn(gate$species$Y, ' + ', gate$species$b, ' -> ', gate$species$value0, ' + ', gate$species$b),
                        jn(gate$species$Y, ' + 2', gate$species$value0, ' -> 3', gate$species$value0)
                        
  )
  
  fast <- 0.1
  slow <- fast / 1000
  
  gate$ki        = c(slow, fast, slow, fast,
                     slow, fast, slow, fast,
                     slow, fast, slow, fast,
                     slow, fast, slow, fast)#rep(1, times = 12)
  
  return(gate)
}

make_jkff <- function(name, jValue, kValue, eValue) {
  entity <- make_circuit(0)
  
  # JK Flip-Flop interface
  interface <- make_generic3to1_gate(name, 'JKFF', jValue, kValue, eValue)
  
  entity$inputs <- list(interface$species$input1,
                        interface$species$input2,
                        interface$species$input3)
  
  entity$output <- c(interface$species$output)
  
  # Instantiate the sub-entities
  jand <- make_and_gate(jn(name, '_jand'), 0, 0)
  kand <- make_and_gate(jn(name, '_kand'), 0, 0)
  jkor <- make_or_gate(jn(name, '_jkor'), 0, 0)
  ffd  <- make_flipflopd(jn(name, '_bit'), 0, 0)
  
  # Insert sub-entities
  entity <- circuit_insert_gate(entity, interface)
  entity <- circuit_insert_gate(entity, jand)
  entity <- circuit_insert_gate(entity, kand)
  entity <- circuit_insert_gate(entity, jkor)
  entity <- circuit_insert_gate(entity, ffd)
  
  # Sub-entities interconnections
  entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input1, jand$species$input2))
  entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input2, kand$species$input1, negated = TRUE))
  entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input3, ffd$inputs[[2]]))
  
  entity <- circuit_insert_gate(entity, make_link_gate(jand$species$output, jkor$species$input1))
  entity <- circuit_insert_gate(entity, make_link_gate(kand$species$output, jkor$species$input2))
  entity <- circuit_insert_gate(entity, make_link_gate(jkor$species$output, ffd$inputs[[1]]))
  
  entity <- circuit_insert_gate(entity, make_link_gate(ffd$output, interface$species$output))
  entity <- circuit_insert_gate(entity, make_link_gate(ffd$output, jand$species$input1, negated = TRUE))
  entity <- circuit_insert_gate(entity, make_link_gate(ffd$output, kand$species$input2))
  
  # Combine all the sub-entities in the final entity
  entity <- circuit_compile(entity)
  
  return(entity)
}

make_sync_counter <- function(name, tValue, eValue, n = 4) {
  entity <- make_circuit(0)
  
  # Instantiate the sub-entities
  interface <- make_generic2to1_gate(name, 'counter', tValue, eValue)
  entity$inputs <- list(interface$species$input1,
                        interface$species$input2)
  
  entity <- circuit_insert_gate(entity, interface)
  
  t <- list()
  entity$outputs <- list()
  
  if (n >= 1)
  {
    t <- rlist::list.append(t, make_jkff(jn(name, '_fft1'), 0, 0, 0))
    entity$outputs <- rlist::list.append(entity$outputs, t[[1]]$output)
    
    entity <- circuit_insert_gate(entity, t[[1]])
    entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input1, t[[1]]$inputs[[1]]))
    entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input1, t[[1]]$inputs[[2]]))
    entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input2, t[[1]]$inputs[[3]]))
  }
  
  if (n >= 2)
  {
    t <- rlist::list.append(t, make_jkff(jn(name, '_fft2'), 0, 0, 0))
    entity$outputs <- rlist::list.append(entity$outputs, t[[2]]$output)
    
    entity <- circuit_insert_gate(entity, t[[2]])
    entity <- circuit_insert_gate(entity, make_link_gate(t[[1]]$output, t[[2]]$inputs[[1]]))
    entity <- circuit_insert_gate(entity, make_link_gate(t[[1]]$output, t[[2]]$inputs[[2]]))
    entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input2, t[[2]]$inputs[[3]]))
  }
  
  if (n >= 3)
  {
    t <- rlist::list.append(t, make_jkff(jn(name, '_fft3'), 0, 0, 0))
    entity$outputs <- rlist::list.append(entity$outputs, t[[3]]$output)
    
    previous <- 2
    thirdToLast <- 1
    
    entity <- circuit_insert_gate(entity, t[[3]])
    and <- make_and_gate(jn(name, '_and', thirdToLast, previous), 0, 0)
    entity <- circuit_insert_gate(entity, and)
    entity <- circuit_insert_gate(entity, make_link_gate(t[[thirdToLast]]$output, and$species$input1))
    entity <- circuit_insert_gate(entity, make_link_gate(t[[previous]]$output, and$species$input2))
    entity <- circuit_insert_gate(entity, make_link_gate(and$species$output, t[[3]]$inputs[[1]]))
    entity <- circuit_insert_gate(entity, make_link_gate(and$species$output, t[[3]]$inputs[[2]]))
    entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input2, t[[3]]$inputs[[3]]))
    andPrev <- and
  }
  
  if (n >= 4)
  {
    for (i in seq(4, n))
    {
      t <- rlist::list.append(t, make_jkff(jn(name, '_fft', i), 0, 0, 0))
      entity$outputs <- rlist::list.append(entity$outputs, t[[i]]$output)
      
      previous <- (i-1)
      thirdToLast <- (i-2)
      
      entity <- circuit_insert_gate(entity, t[[i]])
      and <- make_and_gate(jn(name, '_and', thirdToLast, previous), 0, 0)
      entity <- circuit_insert_gate(entity, and)
      entity <- circuit_insert_gate(entity, make_link_gate(andPrev$species$output, and$species$input1))
      entity <- circuit_insert_gate(entity, make_link_gate(t[[previous]]$output,   and$species$input2))
      entity <- circuit_insert_gate(entity, make_link_gate(and$species$output, t[[i]]$inputs[[1]]))
      entity <- circuit_insert_gate(entity, make_link_gate(and$species$output, t[[i]]$inputs[[2]]))
      entity <- circuit_insert_gate(entity, make_link_gate(interface$species$input2, t[[i]]$inputs[[3]]))
      andPrev <- and
    }
  }
  
  entity <- circuit_compile(entity)
  
  return(entity)
}

plot_circuit_facet <- function (results, signals, newColNames, legend = TRUE) {
  filteredResults <- results[, c('time', signals)]
  colnames(filteredResults) <- newColNames
  
  groupedResults <- reshape2::melt(filteredResults,
                                   id.vars = c('time'),
                                   measure.vars = dimnames(filteredResults)[[2]][-1])
  
  g <- ggplot2::ggplot(groupedResults, ggplot2::aes(time, value, color = variable)) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::labs(x = 'Time (a.u.)', y = 'Concentration (a.u.)', color = 'Species') +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::facet_wrap(variable ~ ., ncol=1, strip.position = 'right') +
    ggplot2::theme( strip.background = ggplot2::element_blank(),
                    strip.text.x = ggplot2::element_blank() ) +
    ggplot2::geom_line(size = 1.3) +
    ggplot2::scale_y_continuous(breaks=seq(0,1,1))
  
  if (legend == FALSE)
  {
    g <- g + ggplot2::theme(legend.position = "none")
  }
  
  return(g)
}

## Example for manual execution with source().
if (!isNamespace(topenv())) {
    timing  <- seq(0, 36000, 200)
    circuit <- DNArLogic::make_circuit(timing)

    counter <- make_sync_counter('ctr', 1, 0, n=2)
    jkff <- make_jkff('jkff', 0, 0, 0)

    circuit <- circuit_insert_gate(circuit, counter)
    circuit <- circuit_insert_gate(circuit, jkff)

    clock <- make_rgbyclock('clock')
    circuit <- circuit_insert_gate(circuit, clock)
    circuit <- circuit_insert_gate(circuit, make_link_gate(clock$species, counter$inputs[[2]]))
    circuit <- circuit_insert_gate(circuit, make_link_gate(clock$species, jkff$inputs[[3]], negated = TRUE))
    circuit <- circuit_insert_gate(circuit, make_link_gate(counter$outputs[[1]], jkff$inputs[[1]])) # J
    circuit <- circuit_insert_gate(circuit, make_link_gate(counter$outputs[[2]], jkff$inputs[[2]])) # K

    compiled <- DNArLogic::circuit_compile(circuit)
    nreactions_4bcounter <- length(compiled$reactions)

    results <- react(
      species   = compiled$species,
      ci        = compiled$ci,
      reactions = compiled$reactions,
      ki        = compiled$ki,
      t         = compiled$t
    )

    speciesToPlot <- c(counter$outputs[[1]]$value1,
                       counter$outputs[[2]]$value1,
                       clock$species$value0,
                       jkff$output$value1)

    newColNames <- c('time', 'Value J', 'Value K', 'Clock', 'Value Q')
    ctr_g <- plot_circuit_facet(results, speciesToPlot, newColNames, legend = FALSE)
    ctr_g <- ctr_g + ggplot2::geom_vline(xintercept=8000)
    ctr_g <- ctr_g + ggplot2::geom_vline(xintercept=16000)
    ctr_g <- ctr_g + ggplot2::geom_vline(xintercept=24100)
    ctr_g <- ctr_g + ggplot2::geom_vline(xintercept=32100)
    print(ctr_g)
}
