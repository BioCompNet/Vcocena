plot_correlation_distributions <- function(integrated_output, data) {
  test <- integrated_output$integrated_edgelist
  test1 <- dplyr::filter(test, V1 %in% rownames(data$set1_counts) & V2 %in% rownames(data$set1_counts))
  test2 <- dplyr::filter(test, V1 %in% rownames(data$set2_counts) & V2 %in% rownames(data$set2_counts))
  test3 <- dplyr::filter(
    test,
    (V1 %in% rownames(data$set1_counts) & V2 %in% rownames(data$set2_counts)) |
      (V1 %in% rownames(data$set2_counts) & V2 %in% rownames(data$set1_counts))
  )
  par(mfrow = c(1, 3))
  hist(test1$weight, probability = TRUE)
  hist(test2$weight, probability = TRUE)
  hist(test3$weight, probability = TRUE)

  plot(density(test1$weight), col = "red", ylim = c(0, 5))
  lines(density(test2$weight), col = "blue")
  lines(density(test3$weight), col = "green")
  lines(density(test$weight), col = "purple")
}
