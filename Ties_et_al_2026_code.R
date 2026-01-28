#Swift Fox maxent

#Code to get times new roman for plot
install.packages("extrafont")
font_import(prompt = FALSE)
loadfonts(device = "win")

#load the required packages
library(dismo)
library(extrafont)
library(raster)
library(sp)
library(ggplot2)
library(sf)
library(ggpubr)
library(rJava)
library(Metrics)

#setwd
presence_background_points <- read.csv("swiftfox_presence_absence_occurences_with_enviromental_variables.csv")

#Run MaxEnt
maxent_full <- maxent(x = subset(presence_background_points, select = -presence), p = presence_background_points$presence)
maxent_full

#### Cross Validation ####
sdm_csv_raw <- read.csv("swiftfox_presence_absence_occurences_with_enviromental_variables.csv")
summary(sdm_csv_raw)
sdm_csv <- na.omit(sdm_csv_raw)

BackgroundIndices <- (1:10000)
BackgroundFold1 <- sample(BackgroundIndices, 2500)
BackgroundFold2 <- sample(BackgroundIndices[-BackgroundFold1], 2500)
any(BackgroundFold1%in%BackgroundFold2)
BackgroundFold3 <- sample(BackgroundIndices[-c(BackgroundFold1, BackgroundFold2)], 2500)
BackgroundFold4 <- BackgroundIndices[-c(BackgroundFold1, BackgroundFold2, BackgroundFold3)]

NPresence <- nrow(sdm_csv[sdm_csv$presence == 1, ])
PresenceIndices <- 1:NPresence
NPresenceFolds <- floor(NPresence / 4)
PresenceFold1 <- sample(PresenceIndices, NPresenceFolds)
PresenceFold2 <- sample(PresenceIndices[-PresenceFold1], NPresenceFolds)
PresenceFold3 <- sample(PresenceIndices[-c(PresenceFold1, PresenceFold2)], NPresenceFolds)
PresenceFold4 <- PresenceIndices[-c(PresenceFold1, PresenceFold2, PresenceFold3)]

BackgroundFoldList <- list(BackgroundFold1, BackgroundFold2, BackgroundFold3, BackgroundFold4)
PresenceFoldList <- list(PresenceFold1, PresenceFold2, PresenceFold3, PresenceFold4)

PresenceDF <- sdm_csv[sdm_csv$presence == 1, ]
BackgroundDF <- sdm_csv[sdm_csv$presence == 0, ]
AUC_test <- c()
AUC_train <- c()

for (i in 1:4) {
  PTI_i <- unlist(PresenceFoldList[1:4][-i])
  BTI_i <- unlist(BackgroundFoldList[1:4][-i])
  PTestI_i <- unlist(PresenceFoldList[i])
  BTestI_i <- unlist(BackgroundFoldList[i])
  TestingPresencePoints_i <- PresenceDF[PTestI_i, ]
  TestingBackgroundPoints_i <- BackgroundDF[BTestI_i, ]
  AllTestPoints_i <- na.omit(rbind(TestingPresencePoints_i, TestingBackgroundPoints_i))
  
  TrainingPresencePoints_i <- PresenceDF[PTI_i, ]
  TrainingBackgroundPoints_i <- BackgroundDF[BTI_i, ]
  AllTrainingPoints_i <- rbind(TrainingPresencePoints_i, TrainingBackgroundPoints_i)
  
  MaxEnt_i <- maxent(x = subset(AllTrainingPoints_i, select = -presence), p = AllTrainingPoints_i$presence)
  
  # Get predictions for the test set
  Predictions_i <- predict(MaxEnt_i, AllTestPoints_i)
  AUC_test_i <- auc(AllTestPoints_i$presence, Predictions_i)
  AUC_test[i] <- AUC_test_i
  
  # Get predictions for the training set
  Predictions_i_train <- predict(MaxEnt_i, AllTrainingPoints_i)
  if (length(Predictions_i_train) != length(AllTrainingPoints_i$presence)) {
    Predictions_i_train <- Predictions_i_train[1:length(AllTrainingPoints_i$presence)]
  }
  AUC_train_i <- auc(AllTrainingPoints_i$presence, Predictions_i_train)
  AUC_train[i] <- AUC_train_i
}

AUC_test
AUC_train
mean(AUC_test)
mean(AUC_train)
summary(AllTestPoints_i)

#### Making a Grid Map of Maxent Predictions ####
#create point map
#"Enivormental_variables_gridpoints_TexasPandhandle.csv" to large to add to github please contact author directly if desired
grid_csv <- read.csv("Enivormental_variables_gridpoints_TexasPanhandle.csv")

prediction_1 <- predict(maxent_full, grid_csv)
grid_csv$prediction <- prediction_1
xyz <- data.frame(x=grid_csv$x_coo, y=grid_csv$y_coo, z=grid_csv$prediction)

map_points <- st_as_sf(xyz, coords = c("x", "y"))
plot(map_points)

# Convert to an sf object
map_points_2 <- st_as_sf(xyz, coords = c("x", "y"), crs = 4326)  # Adjust CRS as needed

# Save the sf object as a shapefile
st_write(map_points_2, "swiftfox_HabitatSuitability.shp", delete_layer = TRUE)  # Use delete_layer = TRUE if overwriting

#### Making on combined plot ####

original_data <- read.csv("swiftfox_presence_absence_occurences_with_enviromental_variables_unstandarized.csv")

#####
# Sequence in RAW 0–1 space
woody_seq_raw <- seq(
  min(original_data$Woody, na.rm = TRUE),
  max(original_data$Woody, na.rm = TRUE),
  length.out = 100
)

# Standardize because the model was trained on z-scores
woody_mean <- mean(original_data$Woody, na.rm = TRUE)
woody_sd   <- sd(original_data$Woody, na.rm = TRUE)

woody_seq_z <- (woody_seq_raw - woody_mean) / woody_sd

# Convert 0–1 → %
woody_seq_percent <- woody_seq_raw * 100

woody_pred <- data.frame(
  Woody = woody_seq_z,
  High_Developed = median(presence_background_points$High_Developed, na.rm = TRUE),
  PFG = median(presence_background_points$PFG, na.rm = TRUE),
  AFG = median(presence_background_points$AFG, na.rm = TRUE),
  PPT = median(presence_background_points$PPT, na.rm = TRUE),
  claycontent = median(presence_background_points$claycontent, na.rm = TRUE)
)

woody_pred$Prediction <- predict(maxent_full, woody_pred)

woody_pred$Woody_percent <- woody_seq_percent

#Woody Cover

# --- Woody Cover Prediction ---
# Sequence on z-scale (model input)
woody_seq_z <- seq(
  min(presence_background_points$Woody, na.rm = TRUE),
  max(presence_background_points$Woody, na.rm = TRUE),
  length.out = 100
)

# Back-transform to percent cover (assuming raw data in 0-1 scale)
woody_mean <- mean(original_data$Woody, na.rm = TRUE)
woody_sd   <- sd(original_data$Woody, na.rm = TRUE)
woody_seq_raw <- woody_seq_z * woody_sd + woody_mean
woody_seq_percent <- woody_seq_raw * 100  # convert 0-1 → 0-100

# Build prediction dataframe
woody_pred <- data.frame(
  Woody = woody_seq_z,
  High_Developed = median(presence_background_points$High_Developed, na.rm = TRUE),
  PFG = median(presence_background_points$PFG, na.rm = TRUE),
  AFG = median(presence_background_points$AFG, na.rm = TRUE),
  PPT = median(presence_background_points$PPT, na.rm = TRUE),
  claycontent = median(presence_background_points$claycontent, na.rm = TRUE)
)

woody_pred$Prediction <- predict(maxent_full, woody_pred, se.fit = TRUE)
woody_pred$Woody_percent <- woody_seq_percent

plot_woody <- ggplot(woody_pred, aes(x = Woody_percent, y = Prediction)) +
  geom_line(linewidth = 1.2, color = "black") +
  theme_minimal(base_size = 15, base_family = "Times New Roman") +
  labs(
    x = "Woody Cover (%)",
    y = NULL,
    title = "B."
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

#####
#PFG

# --- PFG Cover Prediction ---
# Sequence on z-scale (model input)
pfg_seq_z <- seq(
  min(presence_background_points$PFG, na.rm = TRUE),
  max(presence_background_points$PFG, na.rm = TRUE),
  length.out = 100
)

# Back-transform to percent cover (assuming raw data in 0-1 scale)
pfg_mean <- mean(original_data$PFG, na.rm = TRUE)
pfg_sd   <- sd(original_data$PFG, na.rm = TRUE)
pfg_seq_raw <- pfg_seq_z * pfg_sd + pfg_mean
pfg_seq_percent <- pfg_seq_raw * 100  # convert 0-1 → 0-100

# Build prediction dataframe
pfg_pred <- data.frame(
  PFG = pfg_seq_z,
  High_Developed = median(presence_background_points$High_Developed, na.rm = TRUE),
  Woody = median(presence_background_points$Woody, na.rm = TRUE),
  AFG = median(presence_background_points$AFG, na.rm = TRUE),
  PPT = median(presence_background_points$PPT, na.rm = TRUE),
  claycontent = median(presence_background_points$claycontent, na.rm = TRUE)
)

pfg_pred$Prediction <- predict(maxent_full, pfg_pred, se.fit = TRUE)
pfg_pred$PFG_percent <- pfg_seq_percent

# Plot
plot_pfg <- ggplot(pfg_pred, aes(x = PFG_percent, y = Prediction)) +
  geom_line(linewidth = 1.2, color = "black") +
  theme_minimal(base_size = 15, base_family = "Times New Roman") +
  labs(
    x = "PFG Cover (%)",
    y = NULL,
    title = "C."
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
#####
#AFG
# --- AFG Cover Prediction ---
# Sequence on z-scale (model input)
afg_seq_z <- seq(
  min(presence_background_points$AFG, na.rm = TRUE),
  max(presence_background_points$AFG, na.rm = TRUE),
  length.out = 100
)

# Back-transform to percent cover (assuming raw data in 0-1 scale)
afg_mean <- mean(original_data$AFG, na.rm = TRUE)
afg_sd   <- sd(original_data$AFG, na.rm = TRUE)
afg_seq_raw <- afg_seq_z * afg_sd + afg_mean
afg_seq_percent <- afg_seq_raw * 100  # convert 0-1 → 0-100

# Build prediction dataframe
afg_pred <- data.frame(
  AFG = afg_seq_z,
  High_Developed = median(presence_background_points$High_Developed, na.rm = TRUE),
  PFG = median(presence_background_points$PFG, na.rm = TRUE),
  Woody = median(presence_background_points$Woody, na.rm = TRUE),
  PPT = median(presence_background_points$PPT, na.rm = TRUE),
  claycontent = median(presence_background_points$claycontent, na.rm = TRUE)
)

afg_pred$Prediction <- predict(maxent_full, afg_pred, se.fit = TRUE)
afg_pred$AFG_percent <- afg_seq_percent

# Plot
plot_afg <- ggplot(afg_pred, aes(x = AFG_percent, y = Prediction)) +
  geom_line(linewidth = 1.2, color = "black") +
  theme_minimal(base_size = 15, base_family = "Times New Roman") +
  labs(
    x = "AFG Cover (%)",
    y = NULL,
    title = "D."
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
#####
#PPT
# --- PPT Prediction ---
# Sequence on z-scale (model input)
ppt_seq_z <- seq(
  min(presence_background_points$PPT, na.rm = TRUE),
  max(presence_background_points$PPT, na.rm = TRUE),
  length.out = 100
)

# Back-transform to original scale if needed
ppt_mean <- mean(original_data$PPT, na.rm = TRUE)
ppt_sd   <- sd(original_data$PPT, na.rm = TRUE)
ppt_seq_raw <- ppt_seq_z * ppt_sd + ppt_mean

# Build prediction dataframe
ppt_pred <- data.frame(
  PPT = ppt_seq_z,
  High_Developed = median(presence_background_points$High_Developed, na.rm = TRUE),
  PFG = median(presence_background_points$PFG, na.rm = TRUE),
  Woody = median(presence_background_points$Woody, na.rm = TRUE),
  AFG = median(presence_background_points$AFG, na.rm = TRUE),
  claycontent = median(presence_background_points$claycontent, na.rm = TRUE)
)

ppt_pred$Prediction <- predict(maxent_full, ppt_pred, se.fit = TRUE)
ppt_pred$PPT_value <- ppt_seq_raw  # keep original units

# Plot
plot_ppt <- ggplot(ppt_pred, aes(x = PPT_value, y = Prediction)) +
  geom_line(linewidth = 1.2, color = "black") +
  theme_minimal(base_size = 15, base_family = "Times New Roman") +
  labs(
    x = "Average Precipitation (cm)",
    y = NULL,
    title = "A."
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
#####
#High_Developed
# --- High_Developed Prediction ---
# Sequence on z-scale (model input)
hd_seq_z <- seq(
  min(presence_background_points$High_Developed, na.rm = TRUE),
  max(presence_background_points$High_Developed, na.rm = TRUE),
  length.out = 100
)

# Back-transform to percent cover (assuming raw data in 0-1 scale)
hd_mean <- mean(original_data$High_Developed, na.rm = TRUE)
hd_sd   <- sd(original_data$High_Developed, na.rm = TRUE)
hd_seq_raw <- hd_seq_z * hd_sd + hd_mean
hd_seq_percent <- hd_seq_raw * 100  # convert 0-1 → 0-100

# Build prediction dataframe
hd_pred <- data.frame(
  High_Developed = hd_seq_z,
  PFG = median(presence_background_points$PFG, na.rm = TRUE),
  Woody = median(presence_background_points$Woody, na.rm = TRUE),
  AFG = median(presence_background_points$AFG, na.rm = TRUE),
  PPT = median(presence_background_points$PPT, na.rm = TRUE),
  claycontent = median(presence_background_points$claycontent, na.rm = TRUE)
)

hd_pred$Prediction <- predict(maxent_full, hd_pred, se.fit = TRUE)
hd_pred$High_Developed_percent <- hd_seq_percent

# Plot
plot_hd <- ggplot(hd_pred, aes(x = High_Developed_percent, y = Prediction)) +
  geom_line(linewidth = 1.2, color = "black") +
  theme_minimal(base_size = 15, base_family = "Times New Roman") +
  labs(
    x = "Highly Developed Cover (%)",
    y = NULL,
    title = "F."
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
#####
#Clay
# --- Clay Content Prediction ---
# Sequence on z-scale (model input)
clay_seq_z <- seq(
  min(presence_background_points$claycontent, na.rm = TRUE),
  max(presence_background_points$claycontent, na.rm = TRUE),
  length.out = 100
)

# Back-transform to original units (0-1 scale) and convert to percent
clay_mean <- mean(original_data$claycontent, na.rm = TRUE)
clay_sd   <- sd(original_data$claycontent, na.rm = TRUE)
clay_seq_raw <- clay_seq_z * clay_sd + clay_mean
clay_seq_percent <- clay_seq_raw * 100  # convert 0-1 → 0-100

# Build prediction dataframe
clay_pred <- data.frame(
  claycontent = clay_seq_z,
  PFG = median(presence_background_points$PFG, na.rm = TRUE),
  Woody = median(presence_background_points$Woody, na.rm = TRUE),
  AFG = median(presence_background_points$AFG, na.rm = TRUE),
  PPT = median(presence_background_points$PPT, na.rm = TRUE),
  High_Developed = median(presence_background_points$High_Developed, na.rm = TRUE)
)

clay_pred$Prediction <- predict(maxent_full, clay_pred, se.fit = TRUE)
clay_pred$clay_percent <- clay_seq_percent  # now in 0-100

# Plot
plot_clay <- ggplot(clay_pred, aes(x = clay_percent, y = Prediction)) +
  geom_line(linewidth = 1.2, color = "black") +
  theme_minimal(base_size = 15, base_family = "Times New Roman") +
  labs(
    x = "Clay Content (%)",
    y = NULL,
    title = "E."
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
#####

plot_clay <- plot_clay + 
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 35))  # l = left margin in points

combine_plot <- ggarrange(plot_ppt, plot_woody, 
                          plot_pfg, plot_afg, plot_clay, plot_hd,
                          ncol = 2, nrow = 3, align = "v")

combined_plot_with_title <- annotate_figure(
  combine_plot,
  top = text_grob(
    "Effect of Environmental Variables on Swift Fox Habitat Suitability",
    color = "black", face = "bold", size = 18,
    family = "Times New Roman"  # Set font
  ),
  left = text_grob(
    "Relative Likelihood of Suitability",
    color = "black", rot = 90, face = "bold", size = 16,
    family = "Times New Roman"  # Set font
  )
)
combined_plot_with_title

file_location <- "*"

ggsave(filename = paste0(file_location, "combined_plots_with_title.jpg"), plot = combined_plot_with_title, width = 10, height = 12, dpi = 300)
