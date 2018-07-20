## S2MET Mapping
## 
## Phenotypic data summary and analysis
## This notebook will provide some phenotyping data summaries for the S2MET Mapping project. It will include:

# 1. Basic model for g + e + gxe
# 2. Heritability estimates
# 3. Correlations among environments

# Repository directory
repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))

# Other packages
library(lme4)
library(ggridges)
library(ggrepel)
library(cowplot)
library(modelr)
library(pbr)


### Basic Summaries
### 
### Look at the number of lines per environment

# Find the total number of possible line x environment combinations and find
# the proportion that are observed for each trait
observed <- S2_MET_BLUEs_tp %>% 
  distinct(trait, environment, line_name) %>%
  mutate(observed = TRUE) %>%
  group_by(trait) %>%
  complete(environment, line_name, fill = list(observed = FALSE)) 

(prob_observed <- observed %>%
  # group_by(trait) %>%
  summarize(prop_obs = mean(observed)))

## For each line, calculate the balance - sort lowest to highest
line_observed <- observed %>% 
  group_by(trait, line_name) %>% 
  summarize(p_env = mean(observed)) %>% 
  arrange(p_env)
# How many lines were observed everywhere?
line_observed %>% filter(p_env == 1) %>% summarize(n_line = n_distinct(line_name))

## For each environment, calculate the balance
env_observed <- observed %>% 
  group_by(trait, environment) %>% 
  summarize(p_line = mean(observed)) %>% 
  arrange(p_line)


## Number of total environments and number of environments per trait
n_distinct(S2_MET_BLUEs_tp$environment)

S2_MET_BLUEs_tp %>%
  group_by(trait) %>%
  summarize(n_env = n_distinct(environment))

# trait       n_env
# <chr>       <int>
# 1 GrainYield     35
# 2 HeadingDate    33
# 3 PlantHeight    35



## TP and VP environments
S2_MET_BLUEs_tpvp %>%
  filter(line_name %in% vp) %>%
  group_by(trait) %>%
  summarize(n_env = n_distinct(environment))



## Interpretation: of the possible line x environment combinations in which at least one line was observed, the above is the proportion of the lines that were observed.




## Map of locations
### Plot the trial locations and summarize trial sites by the number of years
### a location was included

# First create a df to plot
trial_info_toplot <- S2_MET_BLUEs_tp %>% 
  distinct(environment) %>% 
  left_join(trial_info) %>% 
  distinct(environment, location, latitude, longitude)

## Order environments on latitude
# Sort on latitute
loc_order <- S2_MET_BLUEs_tp %>% 
  left_join(trial_info) %>% 
  distinct(location, latitude) %>%
  arrange(latitude, location) %>%
  # distinct(environment, longitude) %>% 
  # arrange(longitude, environment) %>% 
  pull(location)


# Create a gradient of colors
f_colors <- colorRampPalette(colors = umn_palette(2)[3:5])
colors_use <- set_names(f_colors(length(loc_order)), loc_order)


# Get the map data for canada
canada <- map_data("world", "Canada")

# Download map data for US by county
usa_county <- map_data(map = "county")
# Download state data
usa_state <- map_data(map = "state")

# Adjust the groups in the states
usa_state <- usa_state %>%
  mutate(group = group + max(canada$group))

# Adjust the groups in the counties
usa_county <- usa_county %>%
  mutate(group = group + max(usa_state$group))

# Tidy and combine
north_america <- bind_rows(usa_state, usa_county, canada)

# Map
g_map <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, fill = NA, color = "grey50", lwd = 0.2) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "grey50", lwd = 0.2) +
  geom_point(data = trial_info_toplot, aes(x = longitude, y = latitude, group = location, color = location), size = 0.5) +
  geom_text_repel(data = trial_info_toplot, aes(x = longitude, y = latitude, label = environment), inherit.aes = FALSE, direction = "y", 
                  hjust = 1.2, vjust = 1, segment.size = 0.1, size = 1.5, box.padding = 0.02) + 
  coord_fixed(ratio = 1.5, xlim = c(-125, -60), ylim = c(35, 50)) +
  scale_color_manual(guide = FALSE, values = colors_use) + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_pnas() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank())


# Save the figure
ggsave(filename = "trial_location_map.jpg", plot = g_map, path = fig_dir,
       width = 8.7, height = 5, units = "cm", dpi = 1000)





## Make an additional map of the locations in which the TP and VP were evaluated
# First create a df to plot
trial_info_toplot <- S2_MET_BLUEs_tpvp %>% 
  filter(line_name %in% vp) %>% 
  distinct(environment) %>% 
  left_join(trial_info) %>% 
  distinct(environment, location, latitude, longitude)

## Order environments on latitude
# Sort on latitute
loc_order <- trial_info_toplot %>% 
  distinct(location, latitude) %>%
  arrange(latitude, location) %>%
  # distinct(environment, longitude) %>% 
  # arrange(longitude, environment) %>% 
  pull(location)

# Create a gradient of colors
f_colors <- colorRampPalette(colors = umn_palette(2)[3:5])
colors_use <- set_names(f_colors(length(loc_order)), loc_order)



# Map
g_map_vp <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, fill = NA, color = "grey50", lwd = 0.2) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "grey50", lwd = 0.2) +
  geom_point(data = trial_info_toplot, aes(x = longitude, y = latitude, group = location, color = location), size = 0.5) +
  geom_text_repel(data = trial_info_toplot, aes(x = longitude, y = latitude, label = environment), inherit.aes = FALSE, direction = "y", 
                  hjust = 1.2, vjust = 1, segment.size = 0.1, size = 2, box.padding = 0.02) + 
  coord_fixed(ratio = 1.5, xlim = c(-125, -60), ylim = c(35, 50)) +
  scale_color_manual(guide = FALSE, values = colors_use) + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_pnas() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank())

# Save the figure
ggsave(filename = "trial_location_vp_map.jpg", plot = g_map_vp, path = fig_dir,
       width = 8.7, height = 5, units = "cm", dpi = 1000)








## Visualization of distributions

# Sort on grain yield environmental mean
env_order <- S2_MET_BLUEs_tp %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(environment, trait) %>% 
  mutate(env_mean = mean(value, na.rm = TRUE)) %>% 
  filter(trait == "GrainYield") %>% 
  complete(environment) %>%
  arrange(env_mean) %>%
  pull(environment) %>% 
  unique()

# Plot
(g_met_dist <- S2_MET_BLUEs_tp %>%
  mutate(environment = parse_factor(environment, levels = env_order)) %>%
  ggplot(aes(x = value, y = environment, fill = environment)) +
  geom_density_ridges() +
  facet_grid(. ~ trait, scales = "free_x") +
  scale_fill_discrete(guide = FALSE) +
  ylab("Environment") +
  xlab("") +
  labs(title = "Trait Distributions in All Environments") +
  theme_bw() )





## Plot, sorted on latitude
# Combine location colors with environment
env_colors <- distinct(S2_MET_BLUEs_tp, environment, location) %>% 
  left_join(data.frame(location = names(colors_use), color = colors_use, stringsAsFactors = FALSE), .) %>% 
  {set_names(.$color, .$environment)}


(g_met_dist <- S2_MET_BLUEs_tp %>%
    mutate(environment = parse_factor(environment, levels = names(env_colors))) %>%
    ggplot(aes(x = value, y = environment, fill = environment)) +
    geom_density_ridges() +
    facet_grid(. ~ trait, scales = "free_x") +
    scale_fill_manual(guide = FALSE, values = env_colors) +
    ylab("Environment") +
    xlab("Phenotypic Value") + 
    theme_pnas() +
    theme(panel.grid = element_blank(), axis.line.y = element_blank()) )


# Save it
ggsave(filename = "met_trait_dist.jpg", plot = g_met_dist, path = fig_dir, 
       width = 8.7, height = 10, units = "cm", dpi = 1000)


## Combine with the map
g_map_and_dist <- plot_grid(
  g_map + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank()),
  g_met_dist, 
  ncol = 1, labels = c("A", "B"), rel_heights = c(0.5, 1), label_size = 10)

# Save it
ggsave(filename = "met_trait_dist_and_map.jpg", plot = g_map_and_dist, path = fig_dir, 
       width = 8.7, height = 12, units = "cm", dpi = 1000)



### Heritability

# Use a different optimizer
library(optimx)

## Calculate heritability on a per-environment basis - just for the TP
# Use the stage-one data
per_env_herit_tp <- stage_one_data %>%
  filter(!str_detect(trial, "S2C1F4"),
         environment %in% unique(S2_MET_BLUEs_tp$environment)) %>%
  ungroup() %>%
  distinct(environment, location, year, trait, heritability)

## Add trial information and create a printable table
per_env_herit_tp_toprint <- per_env_herit_tp %>% 
  left_join(distinct(trial_info, environment, location, state, year, latitude, longitude)) %>% 
  select(environment:year, latitude, longitude, trait, heritability) %>% 
  spread(trait, heritability) %>%
  arrange(year, location) %>%
  # Round
  mutate_at(vars(latitude:PlantHeight), funs(round(., 3)))

write_csv(x = per_env_herit_tp_toprint, path = file.path(fig_dir, "per_env_tp_herit.csv"))


## Create another table with the heritability information for all environments
per_env_herit_all <- stage_one_data %>%
  filter(!str_detect(trial, "S2C1F4")) %>%   
  ungroup() %>%
  distinct(environment, location, year, trait, heritability)

## Add trial information and create a printable table
per_env_herit_all_toprint <- per_env_herit_all %>% 
  left_join(distinct(trial_info, environment, location, state, year, latitude, longitude)) %>% 
  select(environment:year, latitude, longitude, trait, heritability) %>% 
  spread(trait, heritability) %>%
  arrange(year, location)

write_csv(x = per_env_herit_all_toprint, path = file.path(fig_dir, "per_env_all_herit.csv"))



# Group by trait and fit the multi-environment model
stage_two_model_fits <- S2_MET_BLUEs_tp %>% 
  group_by(trait) %>%
  # mutate(value = scale(value)) %>% # Scale and center the observations
  rename(env = environment) %>% # Rename environment
  do({
    # Extract the df
    df <- .
    
    # List of models
    models <- formulas(~ value,
                      base = ~ 1 + (1|env) + (1|line_name),
                      full = add_predictors(base, ~ (1 | line_name:env)))
    
    # Lmer control
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                optimizer ='optimx', optCtrl=list(method='nlminb'), calc.derivs = FALSE)
    # Get the weights
    wts <- df$std_error^2
    
    # Fit the models
    fits <- models %>% 
      map(~lmer(formula = ., data = df, control = lmer_control, weights = wts))
    
    # Table of lines by environments (i.e. plots)
    plot_table <- xtabs(formula = ~ line_name + env, data = df)
    
    # Find the number of environments
    n_e <- plot_table %>%
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Now replicates
    n_r <- plot_table %>% 
      harm_mean()
    
    # Return data_frame
    data_frame(fits = list(fits), n_e = n_e, n_r = n_r) })
    


stage_two_lrt <- stage_two_model_fits %>% 
  mutate(llik = map(fits, ~map_df(., logLik))) %>% 
  unnest(llik) %>%
  select(-fits:-n_r) %>%
  mutate(lr = -2 * (base - full),
         p_value = pchisq(q = lr, df = 1, lower.tail = FALSE) / 2)



# The interaction term is significant for all traits
# The broad-sense heritability is calculated on a entry-mean basis according to the formula

# Now calculate heritability across all environments

# Extract the full model and calculate heritability
herit_out <- stage_two_model_fits %>%
  ungroup() %>%
  mutate(full_fit = map(fits, "full")) %>%
  mutate(herit = select(., full_fit, n_e, n_r) %>% as.list() %>% 
           pmap(~herit(object = ..1, n_e = ..2, n_r = ..3, exp = "line_name / (line_name + (line_name:env / n_e) + (Residual / n_r))")),
         herit = map_dbl(herit, "heritability"))

# trait       herit
# 1 GrainYield  0.289
# 2 HeadingDate 0.985
# 3 PlantHeight 0.942

# # Bootstrapping
# herit_out <- stage_two_model_fits %>%
#   ungroup() %>%
#   mutate(full_fit = map(fits, "full")) %>%
#   select(n_e:full_fit) %>%
#   pmap(~{
#     herit_boot(object = ..3, n_e = ..1, n_r = ..2, boot.reps = 100,
#                exp = "line_name / (line_name + (line_name:env / n_e) + (Residual / n_r))")
#   })

stage_two_herit <- stage_two_model_fits %>%
  ungroup() %>%
  mutate(herit = herit_out)
    

# Plot
g_herit <- stage_two_herit %>% 
  mutate(herit_correct = map_dbl(herit, "heritability")) %>%
  # mutate(herit_correct = heritability - bias) %>%
  ggplot(aes(x = trait, y = herit_correct)) +
  geom_col(aes(fill = "fill"), width = 0.80) +
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.5) +
  scale_fill_discrete(guide = FALSE) +
  ylab("Heritability") +
  xlab("Trait") +
  labs(title = "Broad-Sense Heritability") +
  theme_bw()

ggsave(filename = "heritability.jpg", plot = g_herit, path = fig_dir,
       height = 5, width = 4, dpi = 1000)


# What is the proportion of each variance component to the total phenotypic variance?
prop_varcomp <- stage_two_model_fits %>% 
  mutate(full_fit = map(fits, "full")) %>% 
  do(var_comp = as.data.frame(VarCorr(.$full_fit[[1]]))) %>% 
  unnest() %>%
  select(trait:grp, vcov) %>% 
  group_by(trait) %>%
  mutate(grp = str_replace_all(grp, c("line_name:env" = "GxE", "line_name" = "Genotype", "env" = "Environment")),
         grp = factor(grp, levels = c("Genotype", "GxE", "Environment", "Residual")),
         prop_vcov = vcov / sum(vcov))

# Display nicer
prop_varcomp %>% 
  select(-vcov) %>% 
  spread(grp, prop_vcov)

# trait       Genotype    GxE Environment    Residual
# 1 GrainYield    0.0162 0.0709       0.875 0.0379     
# 2 HeadingDate   0.223  0.111        0.666 0.000000715
# 3 PlantHeight   0.0454 0.0961       0.858 0.00000721


# Plot
g_prop_varcomp <- prop_varcomp %>%
  ggplot(aes(x = trait, y = prop_vcov, fill = grp)) +
  geom_col() +
  ylab("Proportion") +
  xlab("Trait") + 
  scale_fill_discrete(guide = guide_legend(title = "Variance\nComponent")) +
  labs(title = "Components of Phenotypic Variance") +
  theme_bw()


ggsave(filename = "var_comp_prop.jpg", plot = g_prop_varcomp, path = fig_dir,
       height = 5, width = 5, dpi = 1000)

## Output a table
prop_varcomp_table <- prop_varcomp %>% 
  mutate(vcov_prop = str_c(round(vcov, 2), " (", round(prop_vcov * 100, 2), "%)")) %>% 
  select(-vcov:-prop_vcov) %>% 
  spread(trait, vcov_prop)



### Calculate the proportion of GxE that is due to environmental genetic variance
### heterogeneity versus lack of environmental correlation

# For each environment, calculate the genetic variance via reml
env_varG <- S2_MET_BLUEs_tp %>% 
  group_by(trait, environment) %>%
  do(varG = {
    df <- .
    
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                optimizer ='optimx', optCtrl=list(method='nlminb'), calc.derivs = FALSE)
    wts <- df$std_error^2
  
    fit <- lmer(formula = value ~ 1 + (1|line_name), data = df, control = lmer_control,
                weights = wts)  
    
    as.data.frame(VarCorr(fit))[1,"vcov"]
    
  })

# Calculate the heterogeneity
env_varG_V <- env_varG %>% 
  group_by(trait) %>% 
  unnest() %>% 
  summarize(V = var(sqrt(varG)))

# Use the estimate of varGE across all environments to calculate L
env_L <- left_join(env_varG_V, subset(prop_varcomp, grp == "GxE",c(trait, vcov))) %>% 
  mutate(L = vcov - V)

# Use the estimate of genetic variance across all environments to calculate the 
# genetic correlation
env_r <- left_join(env_L, subset(prop_varcomp, grp == "Genotype", c(trait, vcov)), by = "trait") %>% 
  mutate(r_G = vcov.y / (vcov.y + L))

## Genetic correlation across environments:
## GY: 0.203
## HD: 0.726
## PH: 0.334


## What proportion do V and L make up of varGE?
## This is from Li et al 2018
## Add to the variance component table
varGE_components <- env_r %>% 
  select(trait, varGE = vcov.x, V, L) %>% 
  mutate_at(vars(V, L), funs(prop = . / varGE)) %>%
  mutate(heterogeneity = str_c(round(V, 3), " (", round(V_prop, 2) * 100, "%)"), 
         lackCorrelation = str_c(round(L, 3), " (", round(L_prop, 2) * 100, "%)")) %>% 
  select(trait, heterogeneity, lackCorrelation) %>% 
  gather(grp, value, -trait) %>% 
  spread(trait, value)

# Combine with the variance component table
prop_varcomp_table1 <- bind_rows(prop_varcomp_table, varGE_components)


write_csv(x = prop_varcomp_table1, path = file.path(fig_dir, "trait_varcomp.csv"))



