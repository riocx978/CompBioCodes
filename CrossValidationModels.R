# Load libraries
library(caret)
library(mgcv)
library(ggplot2)
library(gam)
setwd("C:/Users/riocx/Documents/Masters new/Spring 2024/Data Science in PH/HW5/")
df <- na.omit(read.csv("Week4_HW_data.csv"))
df$viral_status <- factor(ifelse(df$viral_status == "SC2", "SC2", "other virus"))
df$gender <- factor(df$gender)
# Define models and cross-validation function
crossValidateLogistic <- function(train, test) {
  # Fit logistic regression model
  fit <- glm(viral_status ~ gender + age + IL1B + IFI6 + IL1R2, data = train, family = "binomial")
  pred_prob <- predict(fit, newdata = test, type = "response")
  pred_class <- ifelse(pred_prob > 0.5, "SC2", "other virus")
  return(data.frame(Prediction = pred_class, Reference = test$viral_status))
}
crossValidateGAM <- function(train, test) {
  # Fit GAM model
  fit <- gam(viral_status ~ gender + s(age,df=4) + s(IL1B,df=4) + s(IFI6,df=4) + s(IL1R2,df=4), data = train, family = "binomial")
  pred_prob <- predict(fit, newdata = test, type = "response")
  pred_class <- ifelse(pred_prob > 0.5, "SC2", "other virus")
  return(data.frame(Prediction = pred_class, Reference = test$viral_status))
}
# Perform cross-validation for logistic regression
all_logistic_predictions <- data.frame(Prediction = character(), Reference = character(), stringsAsFactors = FALSE)
for (i in 1:10) {
  set.seed(i)
  folds <- createFolds(df$viral_status, k = 10)
  train <- df[-folds[[i]], ]
  test <- df[folds[[i]], ]
  predictions <- crossValidateLogistic(train, test)
  predictions$Prediction <- factor(predictions$Prediction, levels = levels(test$viral_status))
  predictions$Reference <- factor(predictions$Reference, levels = levels(test$viral_status))
  if (i == 1) {
    all_logistic_predictions <- predictions
  } else {
    all_logistic_predictions <- rbind(all_logistic_predictions, predictions)
  }
}
# Calculate confusion matrix for logistic regression
logistic_conf_matrix <- confusionMatrix(all_logistic_predictions$Prediction, all_logistic_predictions$Reference)
# Perform cross-validation for GAM
all_gam_predictions <- data.frame(Prediction = character(), Reference = character(), stringsAsFactors = FALSE)
for (i in 1:10) {
  set.seed(i) 
  folds <- createFolds(df$viral_status, k = 10)
  train <- df[-folds[[i]], ]
  test <- df[folds[[i]], ]
  predictions <- crossValidateGAM(train, test)
  predictions$Prediction <- factor(predictions$Prediction, levels = levels(test$viral_status))
  predictions$Reference <- factor(predictions$Reference, levels = levels(test$viral_status))
  if (i == 1) {
    all_gam_predictions <- predictions
  } else {
    all_gam_predictions <- rbind(all_gam_predictions, predictions)
}
}
# Calculate confusion matrix for GAM
gam_conf_matrix <- confusionMatrix(all_gam_predictions$Prediction, all_gam_predictions$Reference)
print(logistic_conf_matrix)
print(gam_conf_matrix)
# Performance metrics for each model
logistic_accuracy<- mean(all_logistic_predictions$Prediction == all_logistic_predictions$Reference)
gam_accuracy <- mean(all_gam_predictions$Prediction == all_gam_predictions$Reference)
# Compare average accuracies
if (logistic_accuracy > gam_accuracy) {
  cat("Logistic regression model performs better with average accuracy:", logistic_accuracy, "\n")
} else if (logistic_accuracy < gam_accuracy) {
  cat("GAM model performs better with average accuracy:", gam_accuracy, "\n")
} else {
  cat("Both models perform equally with average accuracy:", logistic_accuracy, "\n")
}
#Graphical representation
model_accuracies <- data.frame(Model = c("Logistic Regression", "GAM"),
                               Accuracy = c(logistic_accuracy, gam_accuracy))
ggplot(model_accuracies, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Comparative Performance of Models",
       x = "Model",
       y = "Accuracy") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))