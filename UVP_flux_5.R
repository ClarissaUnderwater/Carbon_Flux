setwd("/Users/clarissakarthauser/Documents/WHOI/EN688")

library(dplyr)
library(stringr)
library(ggplot2)
library("plotrix")
library(zoo)

################################################################
######## reading datasets
################################################################

# read example UVP data
UVP_ex <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_examples.csv",header = TRUE, sep = ";")

CTD_EN8 <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/CTD/EN688_008.csv",header = TRUE, sep = ";")
CTD_EN9 <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/CTD/EN688_009.csv",header = TRUE, sep = ";")
CTD_EN11 <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/CTD/EN688_011.csv",header = TRUE, sep = ";")
# Add a column indicating the data table of origin
CTD_EN8 <- CTD_EN8 %>% mutate(Station = "EN4")
CTD_EN9 <- CTD_EN9 %>% mutate(Station = "EN5")
CTD_EN11 <- CTD_EN11 %>% mutate(Station = "EN6")

CTD_all <- rbind(CTD_EN8, CTD_EN9, CTD_EN11)

# read abundance data
UVP_abundance <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_EN688/abundance_456.csv",header = TRUE, sep = ";")

UVP_abundance$number_from_256 <- (UVP_abundance$LPM..256.512.µm...no.l.1. 
                                  + UVP_abundance$LPM..0.512.1.02.mm...no.l.1. 
                                  + UVP_abundance$LPM..1.02.2.05.mm...no.l.1.
                                  + UVP_abundance$LPM..2.05.4.1.mm...no.l.1.)

UVP_abundance$number_from_512 <- ( UVP_abundance$LPM..0.512.1.02.mm...no.l.1. 
                                   + UVP_abundance$LPM..1.02.2.05.mm...no.l.1.
                                   + UVP_abundance$LPM..2.05.4.1.mm...no.l.1.)

UVP_abundance$number_from_128 <- ( UVP_abundance$LPM..128.256.µm...no.l.1.
                                   + UVP_abundance$LPM..256.512.µm...no.l.1. 
                                   + UVP_abundance$LPM..0.512.1.02.mm...no.l.1. 
                                   + UVP_abundance$LPM..1.02.2.05.mm...no.l.1.
                                   + UVP_abundance$LPM..2.05.4.1.mm...no.l.1.)

# read POC data
POC <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_EN688/POC.csv",header = TRUE, sep = ";")

POC$ID <- paste(POC$Station, "_",POC$Depth)

POC_1 <- subset(POC, poresize_um == 1)
POC_51 <- subset(POC, poresize_um == 51)
POC <- merge(POC_1, POC_51, by.x = "ID", by.y = "ID")
POC$totalC <- POC$umol_PIC_L.1.y + POC$umol_PIC_L.1.x + POC$umol_POC._L.1.x + POC$umol_POC._L.1.y
POC$Station <- POC$Station.x

# read fluxes
Th_flux <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_EN688/Th_flux.csv",header = TRUE, sep = ";")

 Th_flux <- Th_flux %>%
  mutate(Station = case_when(
    Event == 8 ~ "EN4",
    Event == 10 ~ "EN5",
    Event == 11 ~ "EN6",
    TRUE ~ NA_character_))

 Th_flux_all <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_EN688/Th_flux.csv",header = TRUE, sep = ";")
  
Trap_flux <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_EN688/Trap_flux.csv",header = TRUE, sep = ";")

# read traditional flux
UVP_Flux_Kiko <- read.table("/Users/ckarthauser/Documents/WHOI/Particle_properties/UVP_flux.csv",header = TRUE, sep = ";")

### read ecotaxa data
UVP5 <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_EN688/EN_5_ecotaxa_clean_nodupl.tsv",header = TRUE, sep = "\t")
UVP4 <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_EN688/EN_4_ecotaxa_clean_nodupl.tsv",header = TRUE, sep = "\t")
UVP6 <- read.table("/Users/ckarthauser/Documents/WHOI/EN688/UVP_EN688/EN_6_ecotaxa.tsv",header = TRUE, sep = "\t")

# Add a column indicating the data table of origin
UVP4 <- UVP4 %>% mutate(Station = "EN4")
UVP5 <- UVP5 %>% mutate(Station = "EN5")
UVP6 <- UVP6 %>% mutate(Station = "EN6")

UVP_all <- rbind(UVP4, UVP5, UVP6)

# calculate sizes
UVP_all$ESD <- ((UVP_all$object_area/3.14/2)^(1/2)) *2 *0.073 # 73 µm per px or # 14 px / 1 mm
UVP_all$volume <- 4/3*pi*(UVP_all$ESD/2)^3 
UVP_all$area <- UVP_all$object_area *0.073^2

# subset to only include particles
UVP_p <- UVP_all %>%
  filter(startsWith(object_annotation_hierarchy, "not-living") &
           !object_annotation_category %in% c("artefact", "reflection"))

UVP_z <- UVP_all %>%
  filter(!startsWith(object_annotation_hierarchy, "not-living"))

UVP_other <- UVP_all %>%
  filter(object_annotation_category %in% c("artefact", "reflection"))

################################################################
###### Step 1: UVP_p - Convert UVP skew to microscope skew
################################################################
ggplot() +
geom_point(size = 3, data = UVP_ex, aes(y = skew.1, x = skew)) + # skew.1 is manual one
  geom_smooth(data = UVP_ex, aes(y = skew.1, x = skew), method = "lm", se = F, color = "black" , size = 0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("skew_methods.pdf", width=6, height=3, dpi=300)

#subset_data <- na.omit(data.frame(X = UVP_ex$skew, Y = UVP_ex$skew.1))
#source("/Users/clarissakarthauser/Documents/WHOI/Particle_properties/correlation_summary.R", echo=TRUE)

# equation from correlation above:
UVP_p$skew <- -1.131318181 * UVP_p$object_skew -0.053566155

ggplot() + 
  geom_violin(data = (UVP_p), aes(y = skew, x =  Station, fill = Station)) +
  geom_boxplot(data = (UVP_p), aes(y = skew, x =  Station)) +
  # geom_point(size = 0.01, position = position_jitterdodge(),data = UVP_p, aes(y = skew, x =  Station, fill = Station)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~object_annotation_category)
ggsave("skew_UVP.pdf", width=6, height=6, dpi=300)

# skew over depth
ggplot() + 
  geom_point(data = UVP_p, aes(y = skew, x = object_depth_min, color = object_annotation_category), size = 0.1) +
  geom_smooth(data = UVP_p, aes(y = skew, x = object_depth_min), size = 0.5) +
  coord_flip() +
  scale_x_reverse() +
  facet_wrap(~Station)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("skew_depth.pdf", width=8, height=5, dpi=300)


################################################################
###### Step 2 a: UVP_p - How much carbon is stored in each particle?
################################################################

# calculate carbon from correlation with skew (eq. 3)
UVP_p$carbon <- (-2.02 * UVP_p$skew + 3.1)*UVP_p$area # carbon per particle in µg, area in mm2

# plot carbon by type 
ggplot() +
  geom_boxplot(data = UVP_p, aes(y = carbon/area, x =  object_annotation_category)) 

# carbon over area
ggplot() +
  geom_point(data = subset(UVP_p, UVP_p$ESD < 1.024), aes(y = carbon, x =  area, color = object_annotation_category)) 

# carbon over depth
ggplot() + 
  geom_point(data = UVP_p, aes(y = carbon/area, x = object_depth_min, color = object_annotation_category)) +
  coord_flip() +
  scale_x_reverse() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("C_depth.pdf", width=4, height=3, dpi=300)


################################################################
###### Step 2 b: UVP_p - How fast does each particle sink?
################################################################
# using eq. 2 for all particles
UVP_p$sinking_speed_45_all <- (-85.23 * UVP_p$skew + 131.72) 

ggplot() +
 # geom_violin(data = UVP_p, aes(y = sinking_speed_45, x =  object_annotation_category)) +
  #geom_boxplot(data = UVP_p, aes(y = sinking_speed_45, x =  object_annotation_category)) +
 geom_point(size = 0.01, position = position_jitterdodge(),data = subset(EN_harbor, Station != "harbor"), aes(y = speed_m_d, x =  type, fill = type)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Station)+
  ylim(0,400)
ggsave("sink_measured.pdf", width=6, height=6, dpi=300)
ggsave("sink_estimated.pdf", width=12, height=6, dpi=300)

ggplot() +
  geom_point(data = UVP_p, aes(y = sinking_speed_45_all, x = object_depth_min, color = object_annotation_category)) +
  coord_flip() +
  scale_x_reverse() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("sink_depth.pdf", width=4, height=3, dpi=300)

################################################################
###### Step 2 c: UVP_p - Calculating individual particle flux
################################################################

# flux of each individual vignette in umol C * m/d
UVP_p$individual_flux <- UVP_p$carbon/12 * UVP_p$sinking_speed_45_all 

ggplot() +
  geom_boxplot(data = subset(UVP_p, ESD<1.024), aes(y = individual_flux, x =  object_annotation_category)) 

ggplot() +
  geom_point(data = UVP_p, aes(y = individual_flux, x = object_depth_min, color = object_annotation_category)) +
  coord_flip() +
  scale_x_reverse() +
  facet_wrap(~Station)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("flux_depth.pdf", width=4, height=3, dpi=300)


################################################################
###### Step 3 a: UVP_p_binned - summarize data into depth bins
################################################################
# Create depth bins
create_depth_bins <- function(data) {
  data %>% mutate(depth_bin = cut(object_depth_min, breaks = seq(0, max(object_depth_min + 10, na.rm = TRUE), by = 5), 
                                  labels = seq(2.5, max(object_depth_min + 10, na.rm = TRUE) - 2.5, by = 5), include.lowest = TRUE))
}

UVP_p <- create_depth_bins(UVP_p)
UVP_all <- create_depth_bins(UVP_all)
UVP_z <- create_depth_bins(UVP_z)

# Summarize data by depth bins
UVP_z_binned <- UVP_z %>% group_by(depth_bin, Station) %>% summarize(count_z = sum(!is.na(Station)))
UVP_p_binned <- UVP_p %>% group_by(depth_bin, Station) %>% summarize(total_carbon_bin = sum(carbon, na.rm = TRUE), 
                                                                     median_carbon_bin = median(carbon, na.rm = TRUE),  
                                                                     median_flux_bin = median(individual_flux, na.rm = TRUE),  
                                                                     total_flux_bin = sum(individual_flux, na.rm = TRUE), 
                                                                     count_p = sum(!is.na(Station)))
UVP_all_binned <- UVP_all %>% group_by(depth_bin, Station) %>% summarize(count_all = sum(!is.na(Station)))

# Merge data
UVP_p_binned$Sample <- paste(UVP_p_binned$Station, "_", UVP_p_binned$depth_bin)
UVP_z_binned$Sample <- paste(UVP_z_binned$Station, "_", UVP_z_binned$depth_bin)
UVP_all_binned$Sample <- paste(UVP_all_binned$Station, "_", UVP_all_binned$depth_bin)

UVP_p_binned <- UVP_p_binned %>%
  ungroup() %>%
  select(-Station, -depth_bin) %>%  # Remove duplicate column
  left_join(UVP_z_binned, by = "Sample")

UVP_p_binned <- UVP_p_binned %>%
  ungroup() %>%
  select(-Station, -depth_bin) %>%  # Remove duplicate column
  left_join(UVP_all_binned, by = "Sample")

# Calculate the ratio of particles to all objects
UVP_p_binned$ratio_particles <- UVP_p_binned$count_p / UVP_p_binned$count_all

# Calculate the ratio of animals (UVP_z) to all objects
UVP_p_binned$ratio_z <- UVP_p_binned$count_z / UVP_p_binned$count_all

# Calculate the ratio of animals (UVP_z) to all objects
UVP_p_binned$ratio_other <- 1 - (UVP_p_binned$ratio_z + UVP_p_binned$ratio_particles)

################################################################
###### Step 3 b: UVP_p_binned_grouped - calculate percentage of particle types
################################################################

# by type
UVP_p_binned_grouped <- UVP_p %>%
  group_by(depth_bin, Station, object_annotation_category) %>%
  summarize(
    total_carbon = sum(carbon, na.rm = TRUE),  
    total_flux = sum(individual_flux, na.rm = TRUE),  
    median_carbon = median(carbon, na.rm = TRUE),  
    median_flux = median(individual_flux, na.rm = TRUE),  
    count = sum(!is.na(Station)))

UVP_p_binned_grouped$Sample <- paste(UVP_p_binned_grouped$Station, "_", UVP_p_binned_grouped$depth_bin)

# Convert depth_bin to numeric
UVP_p$depth_bin <- as.numeric(as.character(UVP_p$depth_bin))
UVP_p_binned$depth_bin <- as.numeric(as.character(UVP_p_binned$depth_bin))
UVP_p_binned_grouped$depth_bin <- as.numeric(as.character(UVP_p_binned_grouped$depth_bin))

# Merge particle types with UVP_p_binned
UVP_p_binned_grouped <- UVP_p_binned_grouped %>%
  ungroup() %>%
  select(-Station, -depth_bin) %>%  # Remove duplicate column
  left_join(UVP_p_binned, by = "Sample")

# Calculate the percentage of each particle type within each Sample
UVP_p_binned_grouped <- UVP_p_binned_grouped %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(
    ratio_type = replace_na(count / count_all, 0),
    ratio_z = replace_na(count_z / count_all, 0),  
    ratio_other = 1 - (ratio_z + ratio_particles),
    total_ratio = sum(ratio_type) + ratio_z + ratio_other
  )

UVP_p_binned_grouped_b$depth_bin <- as.numeric(as.character(UVP_p_binned_grouped_b$depth_bin))

# select order of levels
UVP_p_binned_grouped <- UVP_p_binned_grouped %>%
  mutate(object_annotation_category = factor(object_annotation_category, 
                                             levels = c("filament<detritus", "fiber<detritus", "light<detritus", "fluffy<detritus","feces", "compact-dark")))  # Replace with your order

# Figure: Percentage of particles (split by type), animals, and other objects
ggplot() +
  geom_area(position = "fill", data = UVP_p_binned_grouped, aes(x = depth_bin, y = ratio_type, fill = object_annotation_category)) +
 facet_wrap(~Station) +
  coord_flip() +
  scale_x_reverse() +
  theme_minimal() +
  labs(y = "Percentage", x = "Depth Bin", fill = "Category") 
ggsave("ratios_depth_z_types.pdf", width = 9, height = 3, dpi = 300)

median(subset(UVP_p_binned_grouped$ratio_type,UVP_p_binned_grouped$object_annotation_category=="filament<detritus"))/0.73
median(subset(UVP_p_binned_grouped$ratio_type,UVP_p_binned_grouped$object_annotation_category=="light<detritus"))/0.73
median(subset(UVP_p_binned_grouped$ratio_type,UVP_p_binned_grouped$object_annotation_category=="fiber<detritus"))/0.73
median(subset(UVP_p_binned_grouped$ratio_type,UVP_p_binned_grouped$object_annotation_category=="fluffy<detritus"))/0.73
median(subset(UVP_p_binned_grouped$ratio_type,UVP_p_binned_grouped$object_annotation_category=="feces"))/0.73
median(subset(UVP_p_binned_grouped$ratio_type,UVP_p_binned_grouped$object_annotation_category=="compact-dark"))/0.73

median(subset(UVP_p_binned_grouped$total_flux,UVP_p_binned_grouped$object_annotation_category=="filament<detritus"))/564
median(subset(UVP_p_binned_grouped$total_flux,UVP_p_binned_grouped$object_annotation_category=="light<detritus"))/564
median(subset(UVP_p_binned_grouped$total_flux,UVP_p_binned_grouped$object_annotation_category=="fiber<detritus"))/564
median(subset(UVP_p_binned_grouped$total_flux,UVP_p_binned_grouped$object_annotation_category=="fluffy<detritus"))/564
median(subset(UVP_p_binned_grouped$total_flux,UVP_p_binned_grouped$object_annotation_category=="feces"))/564
median(subset(UVP_p_binned_grouped$total_flux,UVP_p_binned_grouped$object_annotation_category=="compact-dark"))/564

################################################################
###### Step 3 c: UVP_C - merging UVP_abundance with UVP_p_binned - 
###### How much carbon from particles > 512 µm is in each liter of seawater?
################################################################

UVP_abundance$Sample <- paste(UVP_abundance$Station, "_",UVP_abundance$depth)

# Merge particle abundances with UVP_p_binned
UVP_C <- UVP_p_binned %>%
  ungroup() %>%
  select(-Station, -depth_bin) %>%  # Remove duplicate column
  left_join(UVP_abundance, by = "Sample")

UVP_C_grouped <- UVP_p_binned_grouped %>%
  ungroup() %>%
  select(-Station, -depth_bin) %>%  # Remove duplicate column
  left_join(UVP_abundance, by = "Sample")

# calculate carbon per liter based on number from abundance * µg C per particle 

# C all particles
UVP_C <- UVP_C %>% mutate(C_uM_512_mean = number_from_512 /12 * (total_carbon_bin / count_p) * ratio_particles)  # mean particles
UVP_C <- UVP_C %>% mutate(C_uM_512_median = number_from_512 /12 * (median_carbon_bin) * ratio_particles) # median particles
UVP_C <- UVP_C %>% mutate(C_uM_512_all_objects = number_from_512 /12 * (median_carbon_bin) ) # median all objects

UVP_C_grouped <- UVP_C_grouped %>% mutate(C_uM_512_median = number_from_512 / 12 * (median_carbon) * ratio_type) # median

# select order of levels
UVP_C_grouped <- UVP_C_grouped %>%
  mutate(object_annotation_category = factor(object_annotation_category, 
                                             levels = c("filament<detritus", "fiber<detritus", "light<detritus", "fluffy<detritus","feces", "compact-dark")))  # Replace with your order

# Figure 6a: Percentage of count
ggplot() +
  geom_area(position = "stack", data = UVP_C_grouped, aes(x = depth, y = count / count_p * number_from_512, fill = object_annotation_category)) +
  facet_wrap(~Station) +
  coord_flip() +
  #ylim(0,1)+
  scale_x_reverse() +
  theme_minimal() +
  labs(y = "Number per liter", x = "Depth Bin", fill = "Category") 
ggsave("ratios_depth_count.pdf", width = 9, height = 3, dpi = 300)


# Figure 6 b: Percentage of carbon
ggplot() +
  geom_area(position = "stack", data = UVP_C_grouped, aes(x = depth, y = C_uM_512_median, fill = object_annotation_category)) +
  facet_wrap(~Station) +
  coord_flip() +
 # ylim(0,1)+
  scale_x_reverse() +
  theme_minimal() +
  labs(y = "Carbon µmol per L", x = "Depth Bin", fill = "Category") 
ggsave("ratios_depth_C.pdf", width = 9, height = 3, dpi = 300)

################################################################
###### Step 3 d: UVP_C - How much carbon flux from particles > 512 µm?
################################################################

# Assuming UVP_p_binned_grouped has columns `depth_bin`, `percentage_flux`, `object_annotation_category`
# and UVP_C_grouped has column `C_from_512_vignettes`
UVP_C_grouped <- UVP_C_grouped %>%
  mutate(flux_from_512_vignettes = total_flux / UVP_C_grouped$Sampled.volume..L. ) # per m3, mmol, corrected for non-particle objects

UVP_C_grouped <- UVP_C_grouped %>%
  mutate(flux_512_median = number_from_512 * ratio_particles * 1000 / 1000 * median_flux * count / count_p) # per L, mmol, corrected for non-particle objects

# median flux in umol C * m/d
# number per L
# total mmol C m-2 d-1

# Figure 6 c: Percentage of flux
ggplot() +
  geom_area(position = "fill", data = UVP_C_grouped, aes(x = depth, y = flux_512_median, fill = object_annotation_category)) +
  facet_wrap(~Station) +
  coord_flip() +
  #ylim(0,1)+
  scale_x_reverse() +
  theme_minimal() +
  labs(y = "Flux mmol C m-2 day-1", x = "Depth Bin", fill = "Category") 
ggsave("ratios_depth_flux_rough.pdf", width = 9, height = 3, dpi = 300)

median(subset(UVP_C_grouped$flux_512_median,UVP_C_grouped$object_annotation_category=="filament<detritus"))/13.07
median(subset(UVP_C_grouped$flux_512_median,UVP_C_grouped$object_annotation_category=="light<detritus"))/13.07
median(subset(UVP_C_grouped$flux_512_median,UVP_C_grouped$object_annotation_category=="fiber<detritus"))/13.07
median(subset(UVP_C_grouped$flux_512_median,UVP_C_grouped$object_annotation_category=="fluffy<detritus"))/13.07
median(subset(UVP_C_grouped$flux_512_median,UVP_C_grouped$object_annotation_category=="feces"))/13.07
median(subset(UVP_C_grouped$flux_512_median,UVP_C_grouped$object_annotation_category=="compact-dark"))/13.07


################################################################
###### Step 3 e: UVP_C - How much carbon from small particles is in each liter of seawater?
################################################################
# we have carbon per particle in size M, it correlates with area - extrapolate to S 
##### carbon over area
ggplot(UVP_p, aes(y = carbon, x = area, color = object_annotation_category)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +  # Adds linear regression line
    stat_poly_eq(aes(label = after_stat(eq.label)), 
                 formula = y ~ 0 + x, 
                 parse = TRUE, 
                 label.x.npc = "left", label.y.npc = "top") + # Adds equation
   # facet_wrap(~object_annotation_category) +
    theme_minimal()
 ggsave("C_area_one.pdf", width=5, height=4, dpi=300)

# calculate correlation, force through zero
model <- summary(lm(lm(UVP_p$carbon ~ 0 + UVP_p$area)))
model$coefficients[, "Estimate"]
model$coefficients[, "Std. Error"]
model$r.squared
summary_model$coefficients[, 4] 

UVP_C$Slope_0 <- 3.599552
UVP_C$Intercept_0 <- 0

# biolvolume / number = average volume per particle in bin
# calculate average ESD per particle in bin, calculate C from that, multiply with number
# 73 µm per px, size classes < 128 can not be used accurately here

# C small 
UVP_C$ESR_128 <-  ((3 / (4 * pi)) * (UVP_C$LPM.biovolume..128.256.µm...mm3.l.1. / UVP_C$LPM..128.256.µm...no.l.1.))^(1/3)
UVP_C$cross_sectional_area_128 <- UVP_C$ESR_128^2 * pi

UVP_C$ESR_256 <- ((3 / (4 * pi)) * (UVP_C$LPM.biovolume..256.512.µm...mm3.l.1. / UVP_C$LPM..256.512.µm...no.l.1.))^(1/3)
UVP_C$cross_sectional_area_256 <- UVP_C$ESR_256^2 * pi

UVP_C$C_p_p_128 <- (UVP_C$cross_sectional_area_128 * UVP_C$Slope_0) 
UVP_C$C_p_p_256 <- (UVP_C$cross_sectional_area_256 * UVP_C$Slope_0)

UVP_C$C_128 <- UVP_C$LPM..128.256.µm...no.l.1. * UVP_C$C_p_p_128 * UVP_C$ratio_particles / 12
UVP_C$C_256 <- UVP_C$LPM..256.512.µm...no.l.1. * UVP_C$C_p_p_256 * UVP_C$ratio_particles / 12

# standard: 70–120 mg m−3 = 70-120 µg per L = 0.6-1 µmol L-1

### plot C per liter Wokil and UVP
ggplot() + 
  geom_line(data = UVP_C, aes(y = C_uM_512_median, x =  depth, color="512")) +
geom_line(data = subset(UVP_C), aes(y = C_uM_512_median + C_256  + C_128, x =  depth, color="all")) +
  geom_point(data = subset(POC, !(Station.x %in% c("EN1", "EN2", "EN3"))), aes(y = umol_POC._L.1.x+umol_PIC_L.1.x, x =  Depth.x, color="filter > 1 µm"), shape = 3) +
  geom_point(data = subset(POC, !(Station.x %in% c("EN1", "EN2", "EN3"))), aes(y = umol_POC._L.1.y +umol_PIC_L.1.y, x =  Depth.x, color="filter > 51 µm")) +
  coord_flip() +
  scale_x_reverse() +
  ylim(0,2)+
  facet_wrap(~Station)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("POC_UVP.pdf", width=9, height=4, dpi=300)

# compare shapes of profiles
POC$depth <- POC$Depth.x - 2.51/
POC$ID2 <- paste(POC$Station, "_",POC$depth)
POC_compare <- merge(UVP_C, POC, by.x = "Sample", by.y = "ID2")
POC_compare <- subset(POC_compare, Depth.x>30)

subset_data <- na.omit(data.frame(X = POC_compare$C_uM_512_median, Y = POC_compare$umol_POC._L.1.y + POC_compare$umol_PIC_L.1.y))
subset_data <- na.omit(data.frame(X = POC_compare$C_uM_512_median + POC_compare$C_256, Y = POC_compare$umol_POC._L.1.y + POC_compare$umol_PIC_L.1.y))
subset_data <- na.omit(data.frame(X = POC_compare$C_uM_512_median + POC_compare$C_256 + POC_compare$C_128, Y = POC_compare$umol_POC._L.1.y + POC_compare$umol_PIC_L.1.y))

mean(POC_compare$C_uM_512_median)/mean(POC_compare$umol_POC._L.1.y + POC_compare$umol_PIC_L.1.y)

source("/Users/ckarthauser/Documents/WHOI/Particle_properties/correlation_summary.R", echo=TRUE)

################################################################
###### Step 4 UVP_C - Carbon fluxes
################################################################
# for > 512 µm
UVP_C <- UVP_C %>%
  mutate(Flux_512 = number_from_512*1000 * ratio_particles * (median_flux_bin/1000) ) # per m3, mmol, corrected for non-particle objects

# particle flux: umol C * m /d 
# UVP flux: particle / m3 * umol C * m / d / particle /1000
# results: mmol C /m2 /d

################################################################
###### Step 4b: How much carbon flux from small particles?
################################################################

##### flux over area
ggplot(UVP_p, aes(y = individual_flux, x = area, color = object_annotation_category)) +
  geom_point() +
  xlim(0,1)+
  ylim(0,500)+
  geom_smooth(method = "lm", se = FALSE) +  # Adds linear regression line
  stat_poly_eq(aes(label = after_stat(eq.label)), 
               formula = y ~ 0+ x, 
               parse = TRUE, 
               label.x.npc = "left", label.y.npc = "top") + # Adds equation
  #facet_wrap(~object_annotation_category) +
  theme_minimal()
ggsave("Flux_area_one.pdf", width=5, height=4, dpi=300)

# calculate correlation 
subset_data <- na.omit(data.frame(X = UVP_p$area, Y = UVP_p$individual_flux))
source("/Users/ckarthauser/Documents/WHOI/Particle_properties/correlation_summary.R", echo=TRUE)

# force through zero not possible for flux
model <- summary(lm(UVP_p$individual_flux ~ 0 + UVP_p$area))
model$coefficients[, "Estimate"]
model$coefficients[, "Std. Error"]
model$r.squared

UVP_C$Slope_flux_0 <- 48.67169

# flux small

UVP_C$Flux_p_p_128 <- UVP_C$cross_sectional_area_128 * UVP_C$Slope_flux_0 
UVP_C$Flux_p_p_256 <- UVP_C$cross_sectional_area_256 * UVP_C$Slope_flux_0 

UVP_C$Flux_128 <- UVP_C$LPM..128.256.µm...no.l.1. * UVP_C$Flux_p_p_128 * UVP_C$ratio_particles
UVP_C$Flux_256 <- UVP_C$LPM..256.512.µm...no.l.1. * UVP_C$Flux_p_p_256 * UVP_C$ratio_particles

################################################################
###### Step 4c: Big flux summary figure - comparing with other flux data
################################################################

UVP_Flux_Kiko$Sample <- paste(UVP_Flux_Kiko$Station, "_",UVP_Flux_Kiko$z)

UVP_C <- merge(UVP_C, UVP_Flux_Kiko, by.x = "Sample", by.y = "Sample")
UVP_C$Station <- UVP_C$Station.x

### the big flux summary
ggplot() + 
  geom_vline(xintercept = 0, color = "darkgreen")+
 geom_line(data = subset(Th_flux, Station != "NA" & Depth > 70), aes(y = POC_Flux_mmol_m2_d1, x =  Depth, color="Thorium")) +
  geom_point(data = subset(Th_flux, Station != "NA" & Depth > 70), aes(y = POC_Flux_mmol_m2_d1, x =  Depth, color="Thorium"), size = 0.7) +
  geom_point(data = subset(UVP_C, depth > 70), aes(y = Flux_512, x = depth, color = "our UVP"), size = 0.7, alpha = 0.7) +
  geom_smooth(data = subset(UVP_C, depth > 70), aes(y = Flux_512, x = depth, color = "our UVP"), size = 0.5, alpha = 0.4) +
  geom_smooth(data = subset(UVP_C, depth > 70), aes(y = Flux_mgC_m2, x = depth, color = "UVP general"), size = 0.5, alpha = 0.4) +
  geom_point(data = subset(UVP_C, depth > 70), aes(y = Flux_mgC_m2, x = depth, color = "UVP general"), size = 0.7, alpha = 0.7) +
  geom_point(data = subset(Trap_flux, Station != "EN1"), aes(y = POC_Flux_mmol_C_m2_d1, x =  Depth, color="Trap")) +
  coord_flip() +
  scale_x_reverse() +
  facet_wrap(~Station)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("flux_UVP_line_all.pdf", width=9, height=4, dpi=300)


data_to_export <- data.frame(
  Flux_512 = UVP_C$Flux_512,
  Flux_256 = UVP_C$Flux_256,
  depth = UVP_C$depth, 
  station = UVP_C$Station)

write.csv(data_to_export, "UVP_C_flux_depth.csv", row.names = FALSE)

################################################################
###### Step 5: Supplement extra data Flux
################################################################

Trap_flux$Event <- Trap_flux$Station
Trap_flux$Event <- gsub("EN", "", Trap_flux$Event)

Th_flux_all <- Th_flux_all %>%
  mutate(Event = as.numeric(Event)) %>%
  arrange(Event)

# Filter out rows with missing values in relevant columns
Th_flux_all_clean <- Th_flux_all %>%
  filter(!is.na(POC_Flux_mmol_m2_d1), !is.na(Error), !is.na(Depth))

Trap_flux <- Trap_flux %>%
  mutate(Event = as.numeric(Event)) %>%
  arrange(Event)

Trap_flux$Event <- Trap_flux$Event +4

Th_flux_all$Error_n <- as.numeric(Th_flux_all$Error)

#Th_flux_all$POC_Flux_mmol_m2_d1 contains flux
#Th_flux_all$Event contains stations
#Th_flux_all$Depth contains depths
#Th_flux_all$Teff contains  
# Th_flux_all$POC_Flux_mmol_m2_d1 divided by Th_flux_all$POC_Flux_mmol_m2_d1 of the same station at 75 m

# For each station (Event), get the flux at 75 meters

# Create a copy of the table where Depth is 75 meters
flux_at_75m <- Th_flux_all[Th_flux_all$Depth == 75, c("Event", "POC_Flux_mmol_m2_d1")]
flux_at_90m <- Th_flux_all[Th_flux_all$Depth == 90, c("Event", "POC_Flux_mmol_m2_d1")]

# Rename the column for clarity
colnames(flux_at_75m)[2] <- "POC_Flux_at_75m"
colnames(flux_at_90m)[2] <- "POC_Flux_at_90m"

# Merge the original data with the 75m flux data based on the Event (station)
Th_flux_all <- merge(Th_flux_all, flux_at_75m, by = "Event", all.x = TRUE)
Th_flux_all <- merge(Th_flux_all, flux_at_90m, by = "Event", all.x = TRUE)

Th_flux_all$flux_at_Ez <- (Th_flux_all$POC_Flux_at_90m+2*Th_flux_all$POC_Flux_at_75m)/3

# Calculate Teff by dividing the flux by the flux at 75m for the same station
Th_flux_all$Teff <- Th_flux_all$POC_Flux_mmol_m2_d1 / Th_flux_all$flux_at_Ez * 100

# Replace negative values in Th_flux_all$Teff with 0
Th_flux_all$Teff_pos <- ifelse(Th_flux_all$Teff < 0, 0, Th_flux_all$Teff)

data_to_export <- data.frame(
  Teff = Th_flux_all$Teff_pos,
  depth = Th_flux_all$Depth, station = Th_flux_all$Event)

write.csv(data_to_export, "Th_Teff.csv", row.names = FALSE)

mean(subset(Th_flux_all$Teff_pos, Th_flux_all$Depth == 180 ))
sd(subset(Th_flux_all$Teff_pos, Th_flux_all$Depth == 180))

mean(subset(Th_flux_all$Teff_pos, Th_flux_all$Depth == 350 | Th_flux_all$Depth == 400))
sd(subset(Th_flux_all$Teff_pos, Th_flux_all$Depth == 350 | Th_flux_all$Depth == 400))

Th_flux_75 <- subset(Th_flux_all, Depth > 74)

ggplot()+
  geom_errorbar(data = Th_flux_75, 
                aes(y = POC_Flux_mmol_m2_d1, x = Depth,
                    ymin = POC_Flux_mmol_m2_d1 - Error_n, ymax = POC_Flux_mmol_m2_d1 + Error_n, color="Thorium"), 
               size = 0.2) +
  geom_line(data = subset(Th_flux_75), aes(y = POC_Flux_mmol_m2_d1, x =  Depth, color="Thorium")) +
  geom_point(data = subset(Th_flux_75), aes(y = POC_Flux_mmol_m2_d1, x =  Depth, color="Thorium"), size = 0.5) +
# geom_point(data = subset(Trap_flux), aes(y = POC_Flux_mmol_C_m2_d1, x =  Depth, color="Trap")) +
  coord_flip() +
  scale_x_reverse() +
  facet_wrap(~Event)+
  xlim(500,0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
ggsave("flux_Th.png", width=9, height=9, dpi=300)

################################################################
###### Step 6: Supplement extra data CTDs
################################################################

max(CTD_EN8$Cpar) # 0.01 % = 0.0955  (total max 95.5)
max(CTD_EN9$Cpar) # 0.01 % = 0.0378
max(CTD_EN11$Cpar) # 0.01 % = 0.0294

max(CTD_EN8$FlECO.AFL) # 10 % = 0.189
max(CTD_EN9$FlECO.AFL) # 10 % = 0.31
max(CTD_EN11$FlECO.AFL) # 10 % = 0.15

# < 0.1 % 
UVP_C[UVP_C$Station.x == "EN4", "EZ"] <- 18
UVP_C[UVP_C$Station.x == "EN5", "EZ"] <- 64 # 81 to 0.0378
UVP_C[UVP_C$Station.x == "EN6", "EZ"] <- 70 # 83 to 0.0294

UVP_C[UVP_C$Station.x == "EN4", "EZ_pc"] <- 18
UVP_C[UVP_C$Station.x == "EN5", "EZ_pc"] <- 81 # 81 to 0.0378
UVP_C[UVP_C$Station.x == "EN6", "EZ_pc"] <- 83 # 83 to 0.0294

UVP_C[UVP_C$Station.x == "EN4", "Chl_10"] <- 81   
UVP_C[UVP_C$Station.x == "EN5", "Chl_10"] <- 59 # 81 to 0.0378
UVP_C[UVP_C$Station.x == "EN6", "Chl_10"] <- 82 # 83 to 0.0294

UVP_C[UVP_C$Station.x == "EN4", "Sal_low"] <- 81
UVP_C[UVP_C$Station.x == "EN5", "Sal_low"] <- 59 
UVP_C[UVP_C$Station.x == "EN6", "Sal_low"] <- 82 

ggplot() +
  #geom_line(data = subset(CTD_all, Station != "NA"), aes(y = CTD_all$Potemp090C, x =  DepSM, color="T")) +
  geom_line(data = subset(CTD_all, Station != "NA"), aes(y = CTD_all$Sal00, x =  DepSM)) +
  coord_flip() +
  scale_x_reverse() +
  #xlim(240,180)+
  #ylim(35.3,35.4)+
  facet_wrap(~Station)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("CTD_T.pdf", width=9, height=4, dpi=300)


#### CTD plot
ggplot() +
  geom_line(data = subset(CTD_all, Station != "NA"), aes(y = CTD_all$FlECO.AFL, x =  DepSM, color="Chl a")) +
  #geom_line(data = subset(CTD_all, Station != "NA"), aes(y = Cpar, x =  DepSM, color="PAR")) +
  coord_flip() +
  scale_x_reverse() +
  #xlim(200,0)+
  # ylim(0,0.3)+
  facet_wrap(~Station)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("CTD_Chl.pdf", width=9, height=4, dpi=300)



################################################################
###### Step 7: Correlations and old plots
################################################################

subset_data <- na.omit(data.frame(X = subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN4" & UVP_C$depth>50),
                                             Y= subset(UVP_C$Flux_mgC_m2, UVP_C$Station =="EN4" & UVP_C$depth>50)))


subset_data <- na.omit(data.frame(X = subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN5" & UVP_C$depth>50),
                                  Y= subset(UVP_C$Flux_mgC_m2, UVP_C$Station =="EN5")))
subset_data <- na.omit(data.frame(X = subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN5"),
                                  Y= subset(UVP_C$Kiko_corr, UVP_C$Station =="EN5" )))


subset_data <- na.omit(data.frame(X = subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN6" & UVP_C$depth>50),
                                  Y= subset(UVP_C$Flux_mgC_m2, UVP_C$Station =="EN6" & UVP_C$depth>50)))


subset_data <- na.omit(data.frame(X = UVP_C$flux_from_512_median_vignettes,
                                  Y= UVP_C$Kiko_corr))
subset_data <- na.omit(data.frame(X = UVP_C$flux_from_512_median_vignettes,
                                  Y= UVP_C$Flux_mgC_m2))


source("/Users/clarissakarthauser/Documents/WHOI/Particle_properties/correlation_summary.R", echo=TRUE)


#Trap: 75, 125, 175, 375, 500
#UVP: depth bin 72.5 + 77.5 

# Martin curve: y = ln(flux) against x = ln(depth) regression
# ctds: night - day - day

ggplot() +
  geom_vline( xintercept = log(70), color = "darkgreen", alpha=0.5)+
  geom_vline( xintercept = log(200), color = "grey")+
  # geom_point(data = na.omit(subset(UVP_C, depth > EZ)), aes(y = log(flux_from_512_median_vignettes), log(depth)), size = 1, color = "red", alpha = 0.4) +
  geom_point(data = na.omit(subset(UVP_C, depth > 0)), aes(y = log(flux_from_512_median_vignettes), log(depth)), size = 1, color = "red", alpha = 0.4) +
  geom_point(data = na.omit(subset(UVP_C, depth > 0)), aes(y = log(Flux_mgC_m2), log(depth)), size = 1, color = "gray", alpha = 0.4) +
 # geom_point(data = na.omit(subset(UVP_C, depth > EZ)), aes(y = log(Flux_mgC_m2), log(depth)), size = 1, color = "red", alpha = 0.4) +
  #geom_point(data = na.omit(subset(UVP_C, depth > Chl_10)), aes(y = log(Flux_mgC_m2), log(depth)), size = 1, color = "green", alpha = 0.4) +
  #geom_point(data = subset(UVP_C, depth > 0), aes(y = log(Flux_mgC_m2), log(depth)), size = 1, color = "blue", alpha = 0.4) +
  geom_point(data = na.omit(subset(Th_flux, Depth > 0)), aes(y = log(POC_Flux_mmol_m2_d1), log(Depth)), size = 1, color = "blue", alpha = 0.4) +
geom_point(data = na.omit(subset(Trap_flux, Depth > 0)), aes(y = log(POC_Flux_mmol_C_m2_d1), log(Depth)), size = 1, color = "black") +
  facet_wrap(~Station)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
      axis.text.x = element_text(angle = 90, hjust = 1))


ggplot() +
  geom_point(data = na.omit(subset(Th_flux_all, Depth > 0)), aes(y = log(POC_Flux_mmol_m2_d1), log(Depth)), size = 1, color = "red", alpha = 0.4) +
  facet_wrap(~Event)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1))

subset_data <- na.omit(data.frame(Y = log(subset(Th_flux_all$POC_Flux_mmol_m2_d1, Th_flux_all$POC_Flux_mmol_m2_d1 > 0)),
                                  X= log(subset(Th_flux_all$Depth, Th_flux_all$POC_Flux_mmol_m2_d1 > 0))))



# UVP estimate large, cutoff EZ
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN4" & UVP_C$depth>UVP_C$EZ)),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN4" & UVP_C$depth>UVP_C$EZ))))
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN5" & UVP_C$depth>UVP_C$EZ)),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN5" & UVP_C$depth>UVP_C$EZ))))
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN6" & UVP_C$depth>UVP_C$EZ)),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN6" & UVP_C$depth>UVP_C$EZ))))
source("/Users/clarissakarthauser/Documents/WHOI/Particle_properties/correlation_summary.R", echo=TRUE)

# UVP estimate Elena
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$Flux_mgC_m2, UVP_C$Station =="EN4" & UVP_C$depth>0)),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN4" & UVP_C$depth>0))))
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$Flux_mgC_m2, UVP_C$Station =="EN5" & UVP_C$depth>0)),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN5" & UVP_C$depth>0))))
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$Flux_mgC_m2, UVP_C$Station =="EN6" & UVP_C$depth>0)),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN6" & UVP_C$depth>0))))

# UVP estimate large, cutoff 0
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN4" )),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN4" ))))
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN5" )),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN5" ))))
subset_data <- na.omit(data.frame(Y = log(subset(UVP_C$flux_from_512_median_vignettes, UVP_C$Station =="EN6" )),
                                  X= log(subset(UVP_C$depth, UVP_C$Station =="EN6" ))))
source("/Users/clarissakarthauser/Documents/WHOI/Particle_properties/correlation_summary.R", echo=TRUE)


# thorium all
subset_data <- na.omit(data.frame(Y = log(subset(Th_flux$POC_Flux_mmol_m2_d1, Th_flux$Station =="EN4"  & Th_flux$POC_Flux_mmol_m2_d1 > 0)),
                                  X= log(subset(Th_flux$Depth, Th_flux$Station =="EN4" & Th_flux$POC_Flux_mmol_m2_d1 > 0))))
subset_data <- na.omit(data.frame(Y = log(subset(Th_flux$POC_Flux_mmol_m2_d1, Th_flux$Station =="EN5" & Th_flux$POC_Flux_mmol_m2_d1 > 0)),
                                  X= log(subset(Th_flux$Depth, Th_flux$Station =="EN5" & Th_flux$POC_Flux_mmol_m2_d1 > 0))))
subset_data <- na.omit(data.frame(Y = log(subset(Th_flux$POC_Flux_mmol_m2_d1, Th_flux$Station =="EN6"  & Th_flux$POC_Flux_mmol_m2_d1 > 0)),
                                  X= log(subset(Th_flux$Depth, Th_flux$Station =="EN6"  & Th_flux$POC_Flux_mmol_m2_d1 > 0))))

# thorium below 70
subset_data <- na.omit(data.frame(Y = log(subset(Th_flux$POC_Flux_mmol_m2_d1, Th_flux$Station =="EN4" & Th_flux$Depth>70 & Th_flux$POC_Flux_mmol_m2_d1 > 0)),
                                  X= log(subset(Th_flux$Depth, Th_flux$Station =="EN4" & Th_flux$Depth>70 & Th_flux$POC_Flux_mmol_m2_d1 > 0))))
subset_data <- na.omit(data.frame(Y = log(subset(Th_flux$POC_Flux_mmol_m2_d1, Th_flux$Station =="EN5" & Th_flux$Depth>70 & Th_flux$POC_Flux_mmol_m2_d1 > 0)),
                                  X= log(subset(Th_flux$Depth, Th_flux$Station =="EN5" & Th_flux$Depth>70 & Th_flux$POC_Flux_mmol_m2_d1 > 0))))
subset_data <- na.omit(data.frame(Y = log(subset(Th_flux$POC_Flux_mmol_m2_d1, Th_flux$Station =="EN6" & Th_flux$Depth>70 & Th_flux$POC_Flux_mmol_m2_d1 > 0)),
                                  X= log(subset(Th_flux$Depth, Th_flux$Station =="EN6" & Th_flux$Depth>70 & Th_flux$POC_Flux_mmol_m2_d1 > 0))))

# trap
subset_data <- na.omit(data.frame(Y = log(subset(Trap_flux$POC_Flux_mmol_C_m2_d1, Trap_flux$Station =="EN4"& Trap_flux$Depth>0 )),
                                  X= log(subset(Trap_flux$Depth, Trap_flux$Station =="EN4" & Trap_flux$Depth>0 ))))
subset_data <- na.omit(data.frame(Y = log(subset(Trap_flux$POC_Flux_mmol_C_m2_d1, Trap_flux$Station =="EN5"& Trap_flux$Depth>0 )),
                                  X= log(subset(Trap_flux$Depth, Trap_flux$Station =="EN5" & Trap_flux$Depth>0 ))))

source("/Users/clarissakarthauser/Documents/WHOI/Particle_properties/correlation_summary.R", echo=TRUE)

############################
