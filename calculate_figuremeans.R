calculate_figuremeans <- function(datain, humidityvar, sexvar) {
  
  ##
  # perc of means
  ##
  
  figure4asam_means = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    group_by(Geography) %>% 
    dplyr::summarize(meanA = mean(cfuA), 
                     semA = sd(cfuA)/sqrt(dplyr::n()), 
                     meanL = mean(cfuL), 
                     semL = sd(cfuL)/sqrt(dplyr::n()), 
                     count = dplyr::n()
    ) %>% 
    ungroup() %>% 
    mutate(percA = meanA/(meanA+meanL), percL = meanL/(meanA+meanL), total = meanA+meanL) %>% 
    droplevels()
  
  figure4asam_means <- figure4asam_means %>%
    mutate(latitude = c(1:length(table(figure4asam_means$Geography)))) 
  
  cat("\nperc of means\n")
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanA, method = "spearman")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$total, method = "spearman")
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanL, method = "spearman")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$percL, method = "spearman")
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
  plot(figure4asam_means$latitude, figure4asam_means$percL)  
  
}
