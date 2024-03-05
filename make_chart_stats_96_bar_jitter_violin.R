make_chart_stats_96_bar_jitter_violin <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, logplot = F, label_adjust = 10000, baradjust = 0, baradjust_bottom = 0, loglabeladjust = 0, aceto_bottom = F, yaxis_perc = NULL, perc_errorbars = NULL, add_jitter = NULL, add_violin = NULL) {
  
  bac$testcol <- unlist(unname(bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
#    print(levels(factor(bac$testcol)))
  }
  
  if(logplot == F) {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  } else {
    if(aceto_bottom==F) {
      bac3 <- 1
      bac3 <- bac %>% 
        group_by(testcol) %>%
        dplyr::summarize(cfuLmean = mean(log10(cfuL+1), na.rm=T),  totalmean = mean(log10(cfuA+cfuL+1), na.rm=T)) %>%
        ungroup() %>%
        mutate(cfuAmean = totalmean-cfuLmean) %>%
        dplyr::select(testcol, cfuL = cfuLmean, total = totalmean, cfuA = cfuAmean) %>%
        reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
        mutate(tempvar = paste0(testcol,"_",variable)) %>%
        inner_join(bac %>% 
                     group_by(testcol) %>%
                     dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                     ungroup() %>%
                     reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                     mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                     mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
    } else {
  #    print("hi")
      bac3 <- bac %>% 
        group_by(testcol) %>%
        dplyr::summarize(cfuAmean = mean(log10(cfuA+1), na.rm=T),  totalmean = mean(log10(cfuA+cfuL+1), na.rm=T)) %>%
        ungroup() %>%
        mutate(cfuLmean = totalmean-cfuAmean) %>%
        dplyr::select(testcol, cfuL = cfuLmean, total = totalmean, cfuA = cfuAmean) %>%
        reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
        mutate(tempvar = paste0(testcol,"_",variable)) %>%
        inner_join(bac %>% 
                     group_by(testcol) %>%
                     dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                     ungroup() %>%
                     reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                     mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                     mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
    }
  }
  
  bac3$value.y = ifelse(is.na(bac3$value.y),0,bac3$value.y) ## added 2023-08-12 - maybe remove if things get broken
  
  print(paste0("N:",dim(bac)[1]))
  bac3$addcol <- bac3$value.x
  numvals <- length(table(list(bac$testcol %>% droplevels)))
  if (aceto_bottom == T) {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
    for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
  } else {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
    for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
  }
  
  
  try(rm(ghi, mid_test),T)
  cat("Stats on AAB abundances")
 # print(kruskal.test(bac$cfuA ~ bac$testcol))
  mid_test <- kruskal.test(bac$cfuA ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = T)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
 #     print(ghi)
    } else {
      cfuAclds <- c("a","b")
    }
 #   cat("AAB Plot_values:")
  } else {
    cfuAclds <- c(rep("a", numvals))
  }
  
  try(rm(ghi, mid_test),T)
  cat("\n\nStats on LAB abundances")
#  print(kruskal.test(bac$cfuL ~ bac$testcol))
  mid_test <- kruskal.test(bac$cfuL ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
  #    print(ghi)
    } else {
      cfuLclds <- c("a","b")
    }
 #   cat("LAB Plot_values:")
  } else {
    cfuLclds <- c(rep("a", numvals))
  }
  
#  print(bac3)
  if(aceto_bottom != T) {
    bac3plot <-1 
    bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
      geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
      scale_fill_manual(values = c("red","blue"))+
      theme(legend.position = "bottom") + 
      geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuAclds), vjust = 1)) + #, inherit.aes = T)) +
      geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuLclds), vjust = 1)) + #, inherit.aes = T)+
      theme_cowplot()
  } else {
    bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
      geom_bar(stat = "identity") + 
      geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
      scale_fill_manual(values = c("blue","red"))+
      theme(legend.position = "bottom") +  
      geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuLclds), vjust = 1)) +
      geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuAclds), vjust = 1)) + 
      theme_cowplot()
  }
  
  if(!is.null(yaxis)) {    
    bac3plot <- bac3plot + ylim(yaxis)
  }  
  
  
  ## make the mean and raw data data frames
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    #mutate(percA = cfuA/total, percL = cfuL/total) %>%
    group_by(testcol) %>%
    dplyr::summarize(meanA = mean(cfuA, na.rm=T), meanL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    #    dplyr::summarize(meanA = mean(cfuA, na.rm=T), meanL = mean(percL, na.rm=T), semA = sd(percA)/sqrt(dplyr::n()), semL = sd(percL)/sqrt(dplyr::n())) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = meanA + meanL) %>%
    mutate(aperc = meanA/totalcol, lperc = meanL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc")) 
  bac4$addcol <- bac4$value
  if (aceto_bottom == T) {
    bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
  } else {
    bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
  }
  
  bac4A <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    mutate(percA = cfuA/total, percL = cfuL/total) %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  bac4A$addcol <- bac4A$value
  if (aceto_bottom == T) {
    bac4A$variable = factor(bac4A$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4A$addcol[i] = bac4A$value[i-numvals]}
  } else {
    bac4A$variable = factor(bac4A$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4A$addcol[i] = bac4A$value[i+numvals]}
  }
  
  ## run the stats
  try(rm(ghi, mid_test),T)
  cat("\n\nStats on relative abundances")
 # print(kruskal.test(bac$aperc ~ bac$testcol))
  mid_test <- kruskal.test(bac$aperc ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
  #    print(kruskal.test(bac$aperc ~ bac$testcol))
      cat("AAB Plot_values:")
   #   print(ghi)
    } else {
      apercclds <- c("a","b")
    }
  } else {
    apercclds <- c(rep("a",numvals))
  }
  
  
  
#  print(bac4)
  if(!is.null(yaxis_perc)) {    
    if(aceto_bottom == F) {
      bac4g <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc") %>% mutate(value = value - (1-yaxis_perc[2])))
 #     print(bac4g)
 #     print(1)
      bac4plot <- ggplot(data = NULL) + 
        geom_bar(data = bac4g, aes(x = testcol, y = value, fill = variable), stat="identity") + 
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        ylim(yaxis_perc*1.001) +
        theme_cowplot() 
    } else {
  #    print(2)
      bac4g <- bac4 %>% filter(variable == "aperc") %>% rbind(bac4 %>% filter(variable == "lperc") %>% mutate(value = value - (1-yaxis_perc[2])))
      bac4g$variablerev = factor(bac4g$variable, levels = rev(levels(bac4g$variable)))
      bac4plot <- ggplot(data = NULL) + 
        geom_bar(data = bac4g, aes(x = testcol, y = value, fill = variablerev), stat="identity") + 
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        ylim(yaxis_perc*1.001) +
        theme_cowplot() 
    }
  } else {
    bac4g <- bac4
    if(aceto_bottom == F) {
 #     print(3)
      bac4plot <- ggplot(data = NULL) + 
        geom_bar(data = bac4g, aes(x = testcol, y = value, fill = variable), stat="identity") + 
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        theme_cowplot() 
    } else {
  #    print(4)
      bac4c <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc"))
      bac4c$variablerev = factor(bac4c$variable, levels = rev(levels(bac4c$variable)))
      bac4plot <- ggplot(data = NULL) + 
        geom_bar(data = bac4g, aes(x = testcol, y = value, fill = variablerev), stat="identity") + 
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        theme_cowplot() 
    }
  }
 # print("gothere5")
  
  if(!is.null(add_jitter)) {    
    set.seed(43)
    bac4plot <- bac4plot + geom_jitter(data = bac4A %>% filter(variable == "lperc"), aes(x = testcol, y=value, group = testcol), width = 0.1, alpha = .4, size = .85, color = "white")
  }
  
  if(!is.null(add_violin)) {    
    set.seed(43)
    bac4plot <- bac4plot + geom_violin(data = bac4A, aes(x = testcol, y=value), alpha = 0, width = 0.6)
  }
  
 # print(bac4plot)
  ifelse(plotout == 1, 
         return(bac3plot),
         ifelse(plotout == 2, 
                return(bac4plot),
                ifelse(plotout == 3, 
                       return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))
  
  
  
}
