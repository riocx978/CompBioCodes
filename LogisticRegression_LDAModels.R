setwd("C:/Users/riocx/Documents/Masters new/Spring 2024/Data Science in PH/HW2/")
df <- read.csv("Week4_HW_data.csv")
library(ISLR)

#Build a logistic regression model to predict SC2 vs. non-SC2 (no_virus+other_virus); SC2=0, nonSC2=1
df$viral_status <- ifelse(df$viral_status == "no_virus" | df$viral_status == "other_virus", 1, 0)
df$viral_status <- as.factor(df$viral_status)
glm.fit=glm(viral_status~gender+age+IL1B+IFI6+IL1R2, data=df,family=binomial)
summary(glm.fit)
glm.probs=predict(glm.fit,type="response") 
glm.probs[1:5]
glm.pred=ifelse(glm.probs>0.5, 1, 0)
table(glm.pred, df$viral_status)
mean(glm.pred== df$viral_status)
# Split the data into a training set and a testing set
set.seed(200)
sample <- sample(c(TRUE, FALSE), nrow(df), replace = TRUE, prob = c(0.7, 0.3))
train <- df[sample, ]
test <- df[!sample, ]
model <- glm(viral_status ~ gender + age + IL1B + IFI6 + IL1R2, family = "binomial", data = train)
summary(model)
probs <- predict(model, test, type = "response")
preds <- ifelse(probs > 0.5, 1, 0)
table(preds, test$viral_status)
mean(preds == test$viral_status)

#Build LDA model to predict SC2 vs. no_virus vs. other_virus
library(MASS)
df <- read.csv("Week4_HW_data.csv")
df$viral_status <- as.factor(df$viral_status)
lda_model <- lda(viral_status ~ gender + age + IL1B + IFI6 + IL1R2, data = df)
print(lda_model)
predictions <- predict(lda_model, df)
table(predictions$class, df$viral_status)

#Build KNN model to predict SC2 vs. no_virus vs. other_virus
library(class)
df <- read.csv("Week4_HW_data.csv")
df$gender <- ifelse(df$gender == "male", 1, 0)
df$gender <- as.factor(df$gender)
df$viral_status <- as.factor(df$viral_status)
Xlag <- df[, c("gender", "age", "IL1B", "IFI6", "IL1R2")]
target <- df$viral_status
train <- sample(nrow(df), nrow(df)*0.7)
k <- 1
knn.pred <- knn(Xlag[train,], Xlag[-train,], target[train], k = k)
table(knn.pred, target[-train])
mean(knn.pred == target[-train])
