library(tidyverse)
library(readxl)
library(pROC)
library(randomForest)

animals <- read_excel("speciesTable.xlsx")
# head(animals)

print(unique(animals$`NCBI class`))

mammals <- animals %>%
    filter(`NCBI class` == "Mammalia") %>%
    pull(`NCBI scientific name`) %>%
    unique()

# head(mammals)
# length(mammals)

methyl <- read_tsv("liver_methyl_matrix.txt") %>%
    rename(gene = `...1`) %>%
    filter(!(is.na(gene))) %>%
    mutate(gene = make.unique(gene)) %>%
    column_to_rownames(var = "gene")

# head(methyl)
# dim(methyl)

# separating imprinted gene list
lines <- readLines("ImpGenesKnown.txt")
imprinted_list <- lapply(lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  parts <- parts[parts != ""]
  list(species = parts[1], genes = parts[-1])
})
imprinted_df <- do.call(rbind, lapply(imprinted_list, function(x) {
  if (length(x$genes) == 0) return(NULL)
  data.frame(species = x$species, gene = x$genes, stringsAsFactors = FALSE)
}))
imprinted_genes_all <- unique(imprinted_df$gene)

# add another column in gene matrix to determine whether imprinted
methyl$is_imprinted <- rownames(methyl) %in% imprinted_genes_all
# head(methyl)
# table(methyl$is_imprinted)

######## FIGURE 1 #########
# rotate the table in order to make histogram
# each instance of gene in each differen species is a row
methyl_mammal <- methyl %>%
    rownames_to_column("gene") %>%
    select(gene, is_imprinted, any_of(mammals))
# head(methyl_mammal)
dim(methyl_mammal)

methyl_mammal_long <- methyl_mammal %>%
    pivot_longer(-c(gene, is_imprinted), names_to = "species", values_to = "methylation") %>%
    filter(!is.na(methylation))
# head(methyl_mammal_long)

# summarize methylation by gene
gene_summary_mammal <- methyl_mammal_long %>%
    group_by(gene, is_imprinted) %>%
    summarise(
        mean_methyl = mean(methylation),
        median_methyl = median(methylation),
        var_methyl = var(methylation),
        n_species = n(), # how many species w/ methyl data
        .groups = "drop" # ungroup the data
    ) %>%

    # difference from prediction
    mutate(dist_from_50 = (mean_methyl - 50)) 

# print(gene_summary_mammal)

# make histogram!
fig1 <- ggplot(gene_summary_mammal, aes(x = mean_methyl, fill = is_imprinted)) +
    geom_histogram(binwidth = 5, position = "identity", alpha = 0.6) +
    scale_fill_manual(values = c("black", "purple"), 
    labels = c("Non-imprinted", "Imprinted")) +
    facet_wrap(~is_imprinted, labeller = label_both, scales = "free_y") +
    labs(
        title = "Promoter Methylation Distribution in Mammalian Liver",
        subtitle = "Known imprinted vs non-imprinted genes",
        x = "Mean Promoter Methylation Across Mammalian Species (%)",
        y = "Frequency of Genes"
    ) +
    theme_minimal()

# print(fig1)

fig1_5 <- ggplot(gene_summary_mammal, aes(x = dist_from_50, fill = is_imprinted)) +
    geom_histogram(binwidth = 5, position = "identity", alpha = 0.6) +
    scale_fill_manual(values = c("yellow", "green"), 
    labels = c("Non-imprinted", "Imprinted")) +
    facet_wrap(~is_imprinted, labeller = label_both, scales = "free_y") +
    labs(
        title = "Promoter Methylation Distribution in Mammalian Liver",
        subtitle = "Known imprinted vs non-imprinted genes",
        x = "Distance of Mean Promoter from 50% (%)",
        y = "Frequency of Genes"
    ) +
    theme_minimal()
# print(fig1_5)

# t-test to ensure not imprinted mean is different enough from 50%
# dependent ~ independent
print("\n--- Figure 1: Distance from 50% t-test --- \n")
print(t.test(dist_from_50 ~ is_imprinted, data = gene_summary_mammal))

######## FIGURE 2 #########
#violin plot + boxplot
fig2 <- ggplot(gene_summary_mammal, aes(x = is_imprinted, y = var_methyl, fill = is_imprinted)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    labs(
        title = "Methylation Variance Across Mammalian Species",
        x = NULL,
        y = "Variance in Methylation"
    ) +
    scale_x_discrete(label = c("Non-Imprinted", "Imprinted")) +
    theme_minimal() + 
    theme(legend.position = "none")
# print(fig2)

# print("\n--- Figure 2: Methylation variance t-test ---\n")
# print(t.test(var_methyl ~ is_imprinted, data = gene_summary_mammal))

# linear model using cooks distance (better than pairwise t-test)
varModel <- lm(gene_summary_mammal$var_methyl ~ gene_summary_mammal$is_imprinted)
cooks <- cooks.distance(varModel)

rowNum <- c(1:17370)
cookdf <- data.frame(rowNum, cooks)
cooks_plot <- ggplot(cookdf, aes(x = rowNum, y = cooks)) +
    geom_point() +
    labs(
        title = "Cooks Distance of Variance Linear Model",
        x = NULL,
        y = "Cook's Distance"
    ) +
    theme_minimal()
# print(cooks_plot)

######## FIGURE 3 #########
# logistic regression
model_data <- gene_summary_mammal %>%
    filter(!is.na(var_methyl)) %>%
    mutate(is_imprinted = as.factor(is_imprinted))

# compute inverse-frequency weights
# balance data, bc there are so many more not imprinted genes
n_total <- nrow(model_data)
n_imp <- sum(model_data$is_imprinted == TRUE)
n_not <- sum(model_data$is_imprinted == FALSE)
model_data <- model_data %>%
    mutate(weight = ifelse(is_imprinted == TRUE, n_total/(2 * n_imp), n_total/(2 * n_not)))

model <- glm(is_imprinted ~ var_methyl + mean_methyl, data = model_data, family = "binomial", weights = weight)
# summary(model)

model_data$predicted_response <- predict(model, type = "response")
roc_curve <- roc(model_data$is_imprinted, model_data$predicted_response)
cat("\nAUC:", round(auc(roc_curve), 3), "\n")

fig3 <- ggroc(roc_curve, color = "steelblue", size = 1) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.3, y = 0.2,
            label = paste0("AUC = ", round(auc(roc_curve), 3)),
            size = 5, color = "steelblue") +
    labs(
    title = "Logistic Regression: Predicting Imprinting Status from Methylation",
    subtitle = "Trained on mammalian promoter methylation data",
    x = "Specificity", y = "Sensitivity"
    ) +
  theme_minimal()
# print(fig3)


######## FIGURE 4 #########
# random forest for better prediction?
rf_model <- randomForest(
    is_imprinted ~ var_methyl + mean_methyl,
    data = model_data,
    ntree = 1000,
    mtry = 1, # tune this: typically sqrt(# features)
    class_weight = weight, # adjust for imbalance
    importance = TRUE
)
# summary(rf_model)

model_data$predicted_prob <- predict(rf_model, type = "prob")[, "TRUE"]
roc_curve <- roc(model_data$is_imprinted, model_data$predicted_prob)
cat("RF AUC:", round(auc(roc_curve), 3), "\n")

fig4 <- ggroc(roc_curve, color = "steelblue", size = 1) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.3, y = 0.2,
            label = paste0("AUC = ", round(auc(roc_curve), 3)),
            size = 5, color = "steelblue") +
    labs(
    title = "Random Forest: Predicting Imprinting Status from Methylation",
    subtitle = "Trained on mammalian promoter methylation data",
    x = "Specificity", y = "Sensitivity"
    ) +
    theme_minimal()
print(fig4)