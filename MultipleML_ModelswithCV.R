library(MLDataR)
library(dplyr)
library(tidyr)
library(tidymodels)
library(data.table)
library(ConfusionTableR)
library(OddsPlotty)
library(ROCR)
library(ggplot2)
library(glmnet)
library(caret)
library(randomForest)
library(caTools)
require(gbm)
library(gbm)
library(pROC)
data(package = "MLDataR")
df <- heartdisease
df$Sex <- ifelse(df$Sex == "F", 1, 0)
df$RestingECG <- ifelse(df$RestingECG == "ST", 1, 0)
df$Angina <- ifelse(df$Angina == "Y", 1, 0)
#Explore data
# Boxplot for MaxHR with heart disease status
ggplot(df, aes(x = factor(HeartDisease), y = MaxHR, fill = factor(HeartDisease))) +
  geom_boxplot() +
  labs(x = "Heart Disease", y = "MaxHR", fill = "Heart Disease") +
  theme_minimal()

# Boxplot for HeartPeakReading with heart disease status
ggplot(df, aes(x = factor(HeartDisease), y = HeartPeakReading, fill = factor(HeartDisease))) +
  geom_boxplot() +
  labs(x = "Heart Disease", y = "HeartPeakReading", fill = "Heart Disease") +
  theme_minimal()
# Correlation matrix
# Calculate correlation matrix
correlation_matrix <- cor(df)

# Extract correlations with HeartDisease
heart_disease_correlations <- correlation_matrix["HeartDisease", -c(1, ncol(correlation_matrix))]

# Print the correlations
print(heart_disease_correlations)

# Assuming df_selected contains the relevant variables
df_selected <- select(df, MaxHR, Angina, HeartPeakReading, HeartDisease)
correlation_matrix <- cor(df_selected)

# Create a heatmap
ggplot(data = as.data.frame(correlation_matrix), aes(x = rownames(correlation_matrix), y = rownames(correlation_matrix), fill = HeartDisease)) +
  geom_tile() +
  geom_text(aes(label = round(HeartDisease, 2)), color = "black") +
  scale_fill_gradient(low = "orange", high = "red") +  # Customize the orange color scale
  theme_minimal() +
  labs(title = "Correlation Matrix with Heart Disease") +
  scale_x_discrete(name = "Disease Variables") +
  scale_y_discrete(name = "Disease Variables")

#Split into test data
df <- heartdisease
n <- nrow(df)
set.seed(200)
ntest <- trunc(n / 3)
testid <- sample(1:n, ntest)
#mean absolute prediction error = 0.3040682
# Fit logistic regression model
model <- glm(HeartDisease ~ .,data = df[-testid, ], family = binomial)
# Calculate predicted probabilities for the test data
lpred <- predict(model, newdata = df[testid, ], type = "response")

# Convert probabilities to binary predictions based on a threshold (0.5)
pred <- ifelse(lpred > 0.5, 1, 0)

# Compute accuracy metrics
accuracy <- mean(pred == df$HeartDisease[testid])
sensitivity <- sum(pred == 1 & df$HeartDisease[testid] == 1) / sum(df$HeartDisease[testid] == 1)
specificity <- sum(pred == 0 & df$HeartDisease[testid] == 0) / sum(df$HeartDisease[testid] == 0)
precision <- sum(pred == 1 & df$HeartDisease[testid] == 1) / sum(pred == 1)
recall <- sensitivity
f1_score <- 2 * (precision * recall) / (precision + recall)
# Summarize the model
summary(model)
# Evaluate confusion matrix
pred <- ifelse(lpred > 0.5, 1, 0)  # Convert probabilities to binary predictions
conf_matrix <- table(Predicted = pred, Actual = df$HeartDisease[testid])
print(conf_matrix)
# Print accuracy metrics
cat("Accuracy:", accuracy, "\n")
cat("Sensitivity:", sensitivity, "\n")
cat("Specificity:", specificity, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1-score:", f1_score, "\n")

# Make predictions with Naive Bayes model
df <- na.omit(heartdisease)
df$HeartDisease <- as.factor(df$HeartDisease)
#Naïve Bayes model to predict SC2 vs. no_virus vs. other_virus on the training set and calculate correct prediction rate (accuracy) on the validation set
set.seed(100)
dt = sort(sample(nrow(df), nrow(df)*.7))
train <- df[dt, ]
test <- df[-dt, ]
trctrl <- trainControl(method = "cv", number = 10, savePredictions=TRUE)
nb_fit <- train(factor(HeartDisease) ~ ., data = df, method = "naive_bayes", trControl=trctrl, tuneLength = 0)
nb_fit
# Make predictions with Naive Bayes model
nb.pred <- predict(nb_fit, test)
# Calculate accuracy for Naive Bayes model
confusionMatrix(test$HeartDisease, nb.pred)
nb.acc <- mean(nb.pred == test$HeartDisease)
nb.acc

#LASSO
#Use 10-fold CV and LASSO to find best set of features
# Load data
df <- heartdisease
x <- model.matrix(HeartDisease ~ ., data = df)[, -1]
y <- df$HeartDisease

# Fit Lasso model
fit.lasso <- cv.glmnet(x, y, alpha = 1)

# Plot Lasso coefficient paths and CV plot
plot(fit.lasso, xvar = "lambda", label = TRUE)
coef(fit.lasso)

# Determine best lambda and retrieve coefficients
best_lambda <- fit.lasso$lambda.min
best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coef(best_model)

# Cross-validation
set.seed(123)
folds <- sample(rep(1:10, length = nrow(df)))
cv.errors <- matrix(NA, 10, 19)
test_errors <- numeric(10)

for (k in 1:10) {
  # Split data into training and test sets based on the fold
  x_train <- x[folds != k, ]
  y_train <- y[folds != k]
  x_test <- x[folds == k, ]
  y_test <- y[folds == k]
  
  # Perform LASSO regression on the training set
  cv.lasso <- cv.glmnet(x_train, y_train, alpha = 1)
  
  # Evaluate the model on the test set
  pred <- predict(cv.lasso, newx = x_test, s = "lambda.min")
  test_errors[k] <- mean((y_test - pred)^2)
  
  # Perform cross-validation on the training set
  for (i in 1:19) {
    pred_cv <- predict(cv.lasso, s = cv.lasso$lambda[i], newx = x_train)
    cv.errors[k, i] <- mean((y_train - pred_cv)^2)
  }
}

# Compute average test error
mean_test_error <- mean(test_errors)

# Plot ROC curve and calculate AUC
pred_roc <- prediction(pred, y_test)
perf_roc <- performance(pred_roc, "tpr", "fpr")
plot(perf_roc, colorize = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray")
auc <- performance(pred_roc, measure = "auc")
auc_value <- auc@y.values[[1]]
cat("Area under the ROC curve (AUC):", auc_value, "\n")

# Compute predictions on the test set
pred <- predict(best_model, newx = x_test, s = best_lambda, type = "response")

# Convert probabilities to binary predictions based on a threshold (0.5 by default)
binary_pred <- ifelse(pred > 0.5, 1, 0)

# Compute confusion matrix
conf_matrix <- table(Actual = y_test, Predicted = binary_pred)

# Compute accuracy metrics
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
specificity <- conf_matrix[1, 1] / sum(conf_matrix[1, ])
precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
recall <- sensitivity
f1_score <- 2 * precision * recall / (precision + recall)

# Print accuracy metrics
cat("Accuracy:", accuracy, "\n")
cat("Sensitivity:", sensitivity, "\n")
cat("Specificity:", specificity, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1-score:", f1_score, "\n")

#Random Forest
df <- heartdisease
df$HeartDisease <- as.factor(df$HeartDisease)
df$Sex <- as.factor(df$Sex)
df$RestingECG <- as.factor(df$RestingECG)
df$Angina <- as.factor(df$Angina)
df <- na.omit(df)
set.seed(222)
dt = sort(sample(nrow(df), nrow(df)*.7))
train <- df[dt, ]
test <- df[-dt, ]
rf <- randomForest(HeartDisease ~ ., data=train, proximity=TRUE)
p1 <- predict(rf, train)
confusionMatrix(p1, train$HeartDisease)
p2 <- predict(rf, test)
confusionMatrix(p2, test$HeartDisease)
plot(rf)

#Boosting
df <- heartdisease
# Convert HeartDisease to numeric (assuming it's binary)
df$HeartDisease <- as.numeric(df$HeartDisease)

# Convert categorical variables to factors
df$Sex <- as.factor(df$Sex)
df$RestingECG <- as.factor(df$RestingECG)
df$Angina <- as.factor(df$Angina)

# Split data into train and test sets
set.seed(123)  # For reproducibility
dt <- sample(nrow(df), 0.7 * nrow(df))
train <- df[dt, ]
test <- df[-dt, ]

# Train the boosting model
boost <- gbm(HeartDisease ~ ., data = train, distribution = "bernoulli",
             n.trees = 10000, shrinkage = 0.01, interaction.depth = 4)
summary(boost)
# Predict class probabilities on the test set
predmat <- predict(boost, newdata = test, n.trees = 10000, type = "response")

# Calculate accuracy
pred_labels <- ifelse(predmat > 0.5, 1, 0)
accuracy <- mean(pred_labels == test$HeartDisease)
cat("Accuracy:", accuracy, "\n")

# Calculate AUC-ROC
rocobj <- roc(test$HeartDisease, predmat)
auc <- round(auc(rocobj), 4)
cat("AUC-ROC:", auc, "\n")
# Calculate precision, recall, and F1-score
conf_matrix <- confusionMatrix(factor(pred_labels), factor(test$HeartDisease))
precision <- conf_matrix$byClass["Pos Pred Value"]
recall <- conf_matrix$byClass["Sensitivity"]
f1_score <- conf_matrix$byClass["F1"]
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1-score:", f1_score, "\n")

# Plot ROC curve
plot(rocobj, main = paste0("ROC Curve (AUC = ", auc, ")"))
