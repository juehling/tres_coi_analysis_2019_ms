# Data analyses for "Predictors and consequences of diet variation in a 
# declining generalist aerial insectivore"

# Written by: Jenny Uehling and Conor Taff
# Last updated: 9/26/2022
# Run under R Studio 4.1.0 on a Mac OS

# This is the main analysis script for analyzing data after processing them
# in the data_filtering_script.R, relative_abundance_occurrence_calculations.R,
# and diversity_calculations.R (all located in this repo).

# Install packages and load libraries ------------------------------------------

  ## Install packages for McElreath's package "rethinking" to work -------------

  # RStan:
  # remove.packages(c("StanHeaders", "rstan")) # if already installed previously
  # install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  
  # cmdstanr:
  # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  # cmdstanr::install_cmdstan()
  
  # Other necessary packages and McElreath's rethinking package
  # install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
  # devtools::install_github("rmcelreath/rethinking", branch = experimental)
  
  # devtools::install_github("rmcelreath/rethinking")

  ## Load packages and set ggplot theme ----------------------------------------
  pacman::p_load("tidyverse", "plyr", "dplyr", "here", "ggpubr", "data.table",
                 "gtools", "lme4", "lmerTest", "emmeans", "MuMIn", "DHARMa", "glmmTMB",
                 "sjPlot", "RColorBrewer", "MASS", "coda","mvtnorm","devtools","loo", 
                 "shape", "rethinking")
  
  # tidyverse, plyr, and dplyr for data wrangling
  # here for file reference by relative paths
  # ggpubr for plotting
  # data.table for creating data tables
  # gtools for using logit transformation
  # lme4 for modeling
  # lmerTest for getting p-values from summary tables
  # emmeans for calculating means, confidence intervals, etc.
  # MuMIn for calculating r-squared from GLMM
  # DHARMa for testing glmer model assumptions
  # glmmTMB for a different package for GLMM models
  # sjPlot for created html tables
  # RColorBrewer for examining color plotting options
  # MASS, coda, mvtnorm, devtools, loo, shape, rethinking for plotting confidence intervals

  # Set ggplot theme
  theme_set(theme_classic())

# Description of data frames ---------------------------------------------------

# These are the major data frames we will use:

# "aquatic" contains a row for every sample 
# (prefixes: L19N for day 6 nestlings fecal samples and P19N for all other fecal samples)
aquatic <- read.csv("2_modified_data/aquatic.csv")
  
# Classify control treatments as "Control"
aquatic$exp_treat[aquatic$exp_treat == "Color_Stress_Cont_Con"] <- "Control"
aquatic$exp_treat[aquatic$exp_treat == "Light_Light_Control"] <- "Control"

# Summarize samples in this dataset
table(aquatic$Sex[aquatic$Adult_or_Nestling == "Adult"])
table(aquatic$Age[aquatic$Adult_or_Nestling == "Nestling"])

# "aquatic_nestlingsadults" contains a row for every nestling fecal sample, and information
# about each nestling's mother in the same row.
aquatic_nestlingsadults <- read.csv("2_modified_data/aquatic_nestlingsadults.csv")

# Classify control treatments as "Control" -- nestlings
aquatic_nestlingsadults$n_exp_treat[aquatic_nestlingsadults$n_exp_treat == "Color_Stress_Cont_Con"] <- "Control"
aquatic_nestlingsadults$n_exp_treat[aquatic_nestlingsadults$n_exp_treat == "Light_Light_Control"] <- "Control"

# Classify control treatments as "Control" -- adults
aquatic_nestlingsadults$ad_exp_treat[aquatic_nestlingsadults$ad_exp_treat == "Color_Stress_Cont_Con"] <- "Control"
aquatic_nestlingsadults$ad_exp_treat[aquatic_nestlingsadults$ad_exp_treat == "Light_Light_Control"] <- "Control"

# Can we use mass and wing in the same model? Check for relationship between them.
plot(Mass~Flat_Wing, data = aquatic[aquatic$Adult_or_Nestling == "Adult" ,])
cor.test(aquatic$Flat_Wing[aquatic$Adult_or_Nestling == "Adult"], aquatic$Mass[aquatic$Adult_or_Nestling == "Adult"])
# There appears to be no relationship, so we can use them in the same models.

# Combine alpha diversity information with data frames -------------------------

# "alpha_div" contains a row for every sample ID and the Simpson's diversity index
# for each sample
alpha_div <- read.csv("2_modified_data/alpha_div.csv")
alpha_div <- subset(alpha_div, select= c(sampleID, simpson_psmelt)) # only keep sampleID and Simpson's diversity index

# Add alpha_div to aquatic
aquatic <- merge(aquatic, alpha_div, by = "sampleID")

# Rename alpha_div column so it will merge with aquatic_nestlingsadults to import
# diversity information for nestlings
alpha_div <- dplyr::rename(alpha_div, n_sampleID = sampleID, n_simpson_psmelt = simpson_psmelt)

# Merge alpha_div for nestlings to aquatic_nestlingsadults
aquatic_nestlingsadults <- merge(aquatic_nestlingsadults, alpha_div,
                                 by = "n_sampleID", all.x = TRUE, all.y = FALSE)

# What are the consequences of diet variation? ---------------------------------
  
  ## Do diet characteristics affect nestling mass? (day 12) --------------------

    ### Prepare data frame for modeling and plotting ---------------------------
    # Filter to just nestlings
    aquatic_n <- aquatic[aquatic$Adult_or_Nestling == "Nestling" ,]
    # Take out just day 12 and day 15 nestlings
    aquatic_n <- aquatic_n[aquatic_n$Age == "12" | aquatic_n$Age == "15" ,]
    # Take out just day 12 nestlings
    aquatic_n_day12 <- aquatic_n[aquatic_n$Age == "12" ,]
    
    # Check structure of each variable in the model to make sure it is correct
    str(aquatic_n_day12$Mass)
    str(aquatic_n_day12$percent_aquatic_ra)
    str(aquatic_n_day12$percent_aquatic_occ)
    str(aquatic_n_day12$percent_aquatic_occ_fam)
    str(aquatic_n_day12$simpson_psmelt)
    str(aquatic_n_day12$exp_treat)
    str(aquatic_n_day12$site_box_year)

    ### Plot mass ~ treatment for visualization --------------------------------
    p_exp_treat <- ggplot(aquatic_n_day12) +
      geom_boxplot((aes(x=exp_treat, y=Mass))) + 
      xlab("Experimental Treatment") +
      ylab("Mass") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14))
    p_exp_treat

    ### Model for mass ~ alpha diversity (day 12 nestlings) --------------------
    mod_div <- lmer(Mass ~ scale(simpson_psmelt) + (1|exp_treat) + (1|site_box_year), data = aquatic_n_day12)
    summary(mod_div, ddf = "Kenward-Roger")
    # Checking conditions
    hist(resid(mod_div))
    plot(mod_div)
    
    # Plot for mass ~ diversity
    p <- ggplot(aquatic_n_day12) +
      geom_point(aes(x=simpson_psmelt, y=Mass)) + 
      xlab("Simpson's diversity index") +
      ylab("Mass") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14))
    p
    
    ### Model for mass ~ percent aquatic, relative abundance (day 12 nestlings) -----
    mod_ra <- lmer(Mass ~ scale(percent_aquatic_ra) + (1|exp_treat) + (1|site_box_year), data = aquatic_n_day12)
    summary(mod_ra, ddf = "Kenward-Roger")
    # Checking conditions
    hist(resid(mod_ra))
    plot(mod_ra)
  
    # Plot for mass ~ percent aquatic, relative abundance
    p <- ggplot(aquatic_n_day12) +
      geom_point((aes(x=percent_aquatic_ra, y=Mass))) + 
      xlab("Proportion of nestling diet composed of aquatic insects") +
      ylab("Mass") +
      theme_classic() +
      theme(axis.title = element_text(size = 20)) + 
      theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
    p
    
    ### Model for mass ~ percent aquatic, occurrence distinct taxonomic group (day 12 nestlings) --------
    mod_occ <- lmer(Mass ~ scale(percent_aquatic_occ) + (1|exp_treat) + (1|site_box_year), data = aquatic_n_day12)
    summary(mod_occ, ddf = "Kenward-Roger")
    # Checking conditions
    hist(resid(mod_occ))
    plot(mod_occ)
    
    # Plot mass ~ percent aquatic, occurrence distinct taxonomic group (day 12 nestlings)
    p <- ggplot(aquatic_n_day12) +
      geom_point((aes(x=percent_aquatic_occ, y=Mass))) +
      geom_smooth((aes(x=percent_aquatic_occ, y=Mass)), method=lm , se=TRUE) +
      xlab("Proportion of nestling diet composed of aquatic insects") +
      ylab("Mass") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14))
    p
    
    ### Model for mass ~ percent aquatic, occurrence family (day 12 nestlings) --------
    mod_occ_fam <- lmer(Mass ~ scale(percent_aquatic_occ_fam) + (1|exp_treat) + (1|site_box_year), data = aquatic_n_day12)
    summary(mod_occ_fam, ddf = "Kenward-Roger")
    # Checking conditions
    hist(resid(mod_occ))
    plot(mod_occ)
    
    # Plot for mass ~ percent aquatic, occurrence family (day 12 nestlings)
    p <- ggplot(aquatic_n_day12) +
      geom_point((aes(x=percent_aquatic_occ_fam, y=Mass))) +
      geom_smooth((aes(x=percent_aquatic_occ_fam, y=Mass)), method=lm , se=TRUE) +
      xlab("Proportion of nestling diet composed of aquatic insects") +
      ylab("Mass") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14))
    p

  ## Do nestling diet characteristics predict fledge success? (day 12) ---------

    ### Prepare data frame for modeling and plotting ---------------------------
    # Label the unknown fate bird as died because it probably was abandoned
    aquatic_n_day12$Nestling_Fate[aquatic_n_day12$Nestling_Fate == "Unknown"] <- "Died"
    # Take out nestlings that were predated (4) because there are so few of them
    # that we can't really do anything with them in modeling efforts.
    aquatic_n_day12_binary <- aquatic_n_day12[aquatic_n_day12$Nestling_Fate == "Died" | aquatic_n_day12$Nestling_Fate == "Fledged" ,]

    # Changed "Died" to 0 and "Fledged" to 1
    aquatic_n_day12_binary$Nestling_Fate[aquatic_n_day12_binary$Nestling_Fate == "Died"] <- "0"
    aquatic_n_day12_binary$Nestling_Fate[aquatic_n_day12_binary$Nestling_Fate == "Fledged"] <- "1"
    aquatic_n_day12_binary$Nestling_Fate <- as.integer(aquatic_n_day12_binary$Nestling_Fate)
    
    # Check structure of each variable in the model to make sure it is correct
    str(aquatic_n_day12_binary$Nestling_Fate)
    str(aquatic_n_day12_binary$simpson_psmelt)
    str(aquatic_n_day12_binary$percent_aquatic_ra)
    str(aquatic_n_day12_binary$percent_aquatic_occ)
    str(aquatic_n_day12_binary$percent_aquatic_occ_fam)
    str(aquatic_n_day12_binary$exp_treat)
    str(aquatic_n_day12_binary$site_box_year)

    ### Model for nestling fate ~ alpha diversity ------------------------------
    mod_div <- glmer(Nestling_Fate ~ scale(simpson_psmelt) + (1|exp_treat) + (1|site_box_year), 
                     data = aquatic_n_day12_binary, family = binomial(link = "logit"))
    summary(mod_div)
    r.squaredGLMM(mod_div)
    # Check conditions
    simulationOutput <- simulateResiduals(fittedModel = mod_div)
    plot(simulationOutput)
    testDispersion(simulationOutput)
    
    # ### Plot nestling fate ~ alpha diversity -- FIGURE 2 -----------------------
    # commented out because this takes a long time to run
    # # Create the same model without scaling simpson for ease of plotting
    # mod_div <- glmer(Nestling_Fate ~ simpson_psmelt + (1|exp_treat) + (1|site_box_year), 
    #                  data = aquatic_n_day12_binary, family = binomial(link = "logit"),
    #                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000000)))
    # 
    # # sample from the fit model to get a 'posterior' distribution
    # # Increasing 'n' will result in more samples = smoother, but takes longer
    # post <- mvrnorm(n = 1e5, mu = fixef(mod_div), Sigma = vcov(mod_div))
    # 
    # # create a vector of the values to evaluate the model at
    # # decreasing interval will result in more samples = smoother, but takes longer
    # r <- seq(0, 1, 0.0001)
    # 
    # # evaluate the predicted likelihood of fledging at each value of r
    # # this is getting the values for plotting best fit line
    # lh_fl <- sapply(r, function(z)logistic(mean(post[, 1] + post[, 2]*z)))
    # 
    # # evaluate the confidence interval (same as above but gets high and low)
    # lh_ci <- sapply(r, function(z)logistic(HPDI(post[, 1] + post[, 2]*z)))
    # 
    # # put these together in data frame
    # fit <- data.frame(r = r,
    #                   mu = lh_fl,
    #                   lo = lh_ci[1, ],
    #                   hi = lh_ci[2, ])
    # 
    # # Plot
    # p <- ggplot(data = aquatic_n_day12_binary, mapping = aes(x = simpson_psmelt, y = Nestling_Fate)) +
    #   geom_jitter(height = 0.02, alpha = 0.7, color = "dodgerblue") +
    #   theme_bw() +
    #   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
    #         axis.title = element_text(size = 16), axis.text = element_text(size = 14)) +
    #   xlab("Simpson's diversity index") +
    #   ylab("Likelihood of fledging") +
    #   geom_ribbon(data = fit, inherit.aes = FALSE,
    #               mapping = aes(x = r, ymin = lo, ymax = hi), fill = "dodgerblue", alpha = 0.5) +
    #   geom_line(data = fit, inherit.aes = FALSE,
    #             mapping = aes(x = r, y = mu), color = "dodgerblue", size = 2)
    # p
    # ggsave(here("3_r_scripts/figs/figs_diversity/Fig_2_nd12_aq_fate_div.pdf"), p, width = 5, height = 4, device = "pdf")

    ### Model for nestling fate ~ percent aquatic, relative abundance ----------
    mod_ra <- glmer(Nestling_Fate ~ scale(percent_aquatic_ra) + (1|exp_treat) + (1|site_box_year), 
                    data = aquatic_n_day12_binary, family = binomial(link = "logit"))
    summary(mod_ra)
    r.squaredGLMM(mod_ra)
    # Check conditions
    simulationOutput <- simulateResiduals(fittedModel = mod_ra)
    plot(simulationOutput)
    testDispersion(simulationOutput)

    # Plot nestling fate ~ percent aquatic, relative abundance
    p <- ggplot(aquatic_n_day12_binary) +
      geom_boxplot(aes(x=as.character(Nestling_Fate), y=percent_aquatic_ra)) +
      theme_classic() +
      theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 16))
    p

    ### Model for nestling fate ~ percent aquatic, occurrence distinct taxonomic group ------------------
    mod_occ <- glmer(Nestling_Fate ~ scale(percent_aquatic_occ) + (1|exp_treat) + (1|site_box_year),
                    data = aquatic_n_day12_binary, family = binomial(link = "logit"))
    summary(mod_occ)
    r.squaredGLMM(mod_occ)
    # Check conditions
    simulationOutput <- simulateResiduals(fittedModel = mod_occ)
    plot(simulationOutput)
    testDispersion(simulationOutput)
    
    # Plot nestling fate ~ percent aquatic, occurrence distinct taxonomic group
    p <- ggplot(aquatic_n_day12_binary) +
      geom_boxplot(aes(x=as.character(Nestling_Fate), y=percent_aquatic_occ)) +
      theme_classic() +
      theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 16))
    p
    
    ### Model for nestling fate ~ percent aquatic, occurrence family ------------------
    mod_occ_fam <- glmer(Nestling_Fate ~ scale(percent_aquatic_occ_fam) + (1|exp_treat) + (1|site_box_year),
                     data = aquatic_n_day12_binary, family = binomial(link = "logit"))
    summary(mod_occ_fam)
    r.squaredGLMM(mod_occ_fam)
    # Check conditions
    simulationOutput <- simulateResiduals(fittedModel = mod_occ_fam)
    plot(simulationOutput)
    testDispersion(simulationOutput)
    
    # Plot nestling fate ~ percent aquatic, occurrence family
    p <- ggplot(aquatic_n_day12_binary) +
      geom_boxplot(aes(x=as.character(Nestling_Fate), y=percent_aquatic_occ_fam)) +
      theme_classic() +
      theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 16))
    p
    
  ## Do diet characteristics affect nestling mass? (day 15) --------------------
    
    ### Prepare data frame for modeling and plotting ---------------------------
    # Extract just day 15 nestlings
    aquatic_n_day15 <- aquatic_n[aquatic_n$Age == "15" ,]
    
    # Check structure of each variable in the model to make sure it is correct
    str(aquatic_n_day15$Mass)
    str(aquatic_n_day15$percent_aquatic_ra)
    str(aquatic_n_day15$percent_aquatic_occ)
    str(aquatic_n_day15$percent_aquatic_occ_fam)
    str(aquatic_n_day15$simpson_psmelt)
    str(aquatic_n_day15$exp_treat)
    str(aquatic_n_day15$site_box_year)
    
    ### Plot mass ~ treatment for visualization --------------------------------
    p_exp_treat <- ggplot(aquatic_n_day15) +
      geom_boxplot((aes(x=exp_treat, y=Mass))) + 
      xlab("Experimental Treatment") +
      ylab("Mass") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14))
    p_exp_treat
    
    ### Model for mass ~ alpha diversity (day 15 nestlings) --------------------
    mod_div <- lmer(Mass ~ scale(simpson_psmelt) + (1|exp_treat) + (1|site_box_year), data = aquatic_n_day15)
    summary(mod_div, ddf = "Kenward-Roger")
    # Checking conditions
    hist(resid(mod_div))
    plot(mod_div)
    
    # Plot for mass ~ diversity
    p <- ggplot(aquatic_n_day15) +
      geom_point(aes(x=simpson_psmelt, y=Mass)) + 
      xlab("Simpson's diversity index") +
      ylab("Mass") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14))
    p
    
    ### Model for mass ~ percent aquatic, relative abundance (day 15 nestlings) -----
    mod_ra <- lmer(Mass ~ scale(percent_aquatic_ra) + (1|exp_treat) + (1|site_box_year), data = aquatic_n_day15)
    summary(mod_ra, ddf = "Kenward-Roger")
    # Checking conditions
    hist(resid(mod_ra))
    plot(mod_ra)
    
    # Plot for mass ~ percent aquatic, relative abundance
    p <- ggplot(aquatic_n_day15) +
      geom_point((aes(x=percent_aquatic_ra, y=Mass))) + 
      xlab("Proportion of nestling diet composed of aquatic insects") +
      ylab("Mass") +
      theme_classic() +
      theme(axis.title = element_text(size = 20)) + 
      theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
    p
    
    ### Model for mass ~ percent aquatic, occurrence distinct taxonomic group (day 15 nestlings) --------
    mod_occ <- lmer(Mass ~ scale(percent_aquatic_occ) + (1|exp_treat) + (1|site_box_year), data = aquatic_n_day15)
    summary(mod_occ, ddf = "Kenward-Roger")
    # Checking conditions
    hist(resid(mod_occ))
    plot(mod_occ)
    
    # Plot for mass ~ percent aquatic, occurrence distinct taxonomic group (day 15 nestlings)
    p <- ggplot(aquatic_n_day15) +
      geom_point((aes(x=percent_aquatic_occ, y=Mass))) +
      geom_smooth((aes(x=percent_aquatic_occ, y=Mass)), method=lm , se=TRUE) +
      xlab("Nestling Percent Aquatic") +
      ylab("Mass") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14))
    p

    ### Model for mass ~ percent aquatic, occurrence family (day 15 nestlings) --------
    mod_occ_fam <- lmer(Mass ~ scale(percent_aquatic_occ_fam) + (1|exp_treat) + (1|site_box_year), data = aquatic_n_day15)
    summary(mod_occ_fam, ddf = "Kenward-Roger")
    # Checking conditions
    hist(resid(mod_occ))
    plot(mod_occ)
    
    # Plot for mass ~ percent aquatic, occurrence family (day 15 nestlings)
    p <- ggplot(aquatic_n_day15) +
      geom_point((aes(x=percent_aquatic_occ_fam, y=Mass))) +
      geom_smooth((aes(x=percent_aquatic_occ_fam, y=Mass)), method=lm , se=TRUE) +
      xlab("Nestling Percent Aquatic") +
      ylab("Mass") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14))
    p
    
  ## Do nestling diet characteristics predict fledge success? (day 15) ---------

    ### Prepare data frame for modeling and plotting ------------------------------    

    # Label the unknown fate bird as died because it probably was abandoned
    aquatic_n_day15$Nestling_Fate[aquatic_n_day15$Nestling_Fate == "Unknown"] <- "Died"
    aquatic_n_day15_binary <- aquatic_n_day15
  
    # Changed "Died" to 0 and "Fledged" to 1
    aquatic_n_day15_binary$Nestling_Fate[aquatic_n_day15_binary$Nestling_Fate == "Died"] <- "0"
    aquatic_n_day15_binary$Nestling_Fate[aquatic_n_day15_binary$Nestling_Fate == "Fledged"] <- "1"
    aquatic_n_day15_binary$Nestling_Fate <- as.integer(aquatic_n_day15_binary$Nestling_Fate)
    
    # Check structure of each variable in the model to make sure it is correct
    str(aquatic_n_day15_binary$Nestling_Fate)
    str(aquatic_n_day15_binary$simpson_psmelt)
    str(aquatic_n_day15_binary$percent_aquatic_ra)
    str(aquatic_n_day15_binary$percent_aquatic_occ)
    str(aquatic_n_day15_binary$percent_aquatic_occ_fam)
    str(aquatic_n_day15_binary$Site)
    str(aquatic_n_day15_binary$site_box_year)

    ### Model for nestling fate ~ alpha diversity --------------------------------
    mod_div <- glmer(Nestling_Fate ~ scale(simpson_psmelt) + (1|exp_treat), 
                     data = aquatic_n_day15_binary, family = binomial(link = "logit"))
    summary(mod_div)
    simulationOutput <- simulateResiduals(fittedModel = mod_div)
    plot(simulationOutput)
    testDispersion(simulationOutput)
    
    # Plot for nestling fate ~ alpha diversity
    p <- ggplot(aquatic_n_day15_binary) +
      geom_boxplot(aes(x=as.character(Nestling_Fate), y=simpson_psmelt)) + 
      xlab("Fate") +
      ylab("Simpson's diversity index (day 15 nestlings)") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14)) +
      scale_x_discrete(labels = c("0" = "Died", "1" = "Fledged"))
    p
  
    ### Model for nestling fate ~ percent aquatic, relative abundance ------------
    mod_ra <- glmer(Nestling_Fate ~ scale(percent_aquatic_ra) + (1|exp_treat), 
                    data = aquatic_n_day15_binary, family = binomial(link = "logit"))
    summary(mod_ra)
    r.squaredGLMM(mod_ra)
    simulationOutput <- simulateResiduals(fittedModel = mod_div)
    plot(simulationOutput)
    testDispersion(simulationOutput)

    # Plot for nestling fate ~ percent aquatic, relative abundance
    p <- ggplot(aquatic_n_day15_binary) +
      geom_boxplot((aes(x=as.character(Nestling_Fate), y=percent_aquatic_ra))) + 
      xlab("Fate") +
      ylab("Percent aquatic (relative abundance)") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14)) +
      scale_x_discrete(labels = c("0" = "Died", "1" = "Fledged"))
    p

    ### Model for nestling fate ~ percent aquatic, occurrence distinct taxonomic group ------------------
    mod_occ <- glmer(Nestling_Fate ~ scale(percent_aquatic_occ) + (1|exp_treat), 
                     data = aquatic_n_day15_binary, family = binomial(link = "logit"))
    summary(mod_occ)
    r.squaredGLMM(mod_occ)
    simulationOutput <- simulateResiduals(fittedModel = mod_occ)
    plot(simulationOutput)
    testDispersion(simulationOutput)

    # Plot for nestling fate ~ percent aquatic, occurrence distinct taxonomic group
    p <- ggplot(aquatic_n_day15_binary) +
      geom_boxplot((aes(x=as.character(Nestling_Fate), y=percent_aquatic_occ))) + 
      xlab("Fate") +
      ylab("Percent aquatic (occurrence)") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14)) +
      scale_x_discrete(labels = c("0" = "Died", "1" = "Fledged"))
    p
    
    ### Model for nestling fate ~ percent aquatic, occurrence family ------------------
    mod_occ_fam <- glmer(Nestling_Fate ~ scale(percent_aquatic_occ_fam) + (1|exp_treat), 
                     data = aquatic_n_day15_binary, family = binomial(link = "logit"))
    summary(mod_occ_fam)
    r.squaredGLMM(mod_occ_fam)
    simulationOutput <- simulateResiduals(fittedModel = mod_occ_fam)
    plot(simulationOutput)
    testDispersion(simulationOutput)
    
    # Plot for nestling fate ~ percent aquatic, occurrence family
    p <- ggplot(aquatic_n_day15_binary) +
      geom_boxplot((aes(x=as.character(Nestling_Fate), y=percent_aquatic_occ_fam))) + 
      xlab("Fate") +
      ylab("Percent aquatic (occurrence, family)") +
      theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
      theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
      theme(strip.text = element_text(size = 14)) +
      scale_x_discrete(labels = c("0" = "Died", "1" = "Fledged"))
    p

# Does adult female phenotype predict offspring diet variation? ----------------

  # We will use the aquatic_nestlingsadults data frame because we're looking at
  # the relationship between mom's phenotypic traits and her nestlings' diets.

  # Use gtools package to perform a logit transformation on percentage data (relative abundance).
  aquatic_nestlingsadults$n_percent_aquatic_ra_logit <- gtools::logit(aquatic_nestlingsadults$n_percent_aquatic_ra, min = -0.0001, max = 1.0001)
  
  # Use gtools package to perform a logit transformation on percentage data (occurrence).
  aquatic_nestlingsadults$n_percent_aquatic_occ_logit <- gtools::logit(aquatic_nestlingsadults$n_percent_aquatic_occ, min = -0.0001, max = 1.0001)

  # Use gtools package to perform a logit transformation on percentage data (occurrence family).
  aquatic_nestlingsadults$n_percent_aquatic_occ_fam_logit <- gtools::logit(aquatic_nestlingsadults$n_percent_aquatic_occ_fam, min = -0.0001, max = 1.0001)
  
  ## Prepare data frame for modeling and plotting ------------------------------
    
  # Check structure of each variable in the model to make sure it is correct
  str(aquatic_nestlingsadults$n_percent_aquatic_ra_logit)
  str(aquatic_nestlingsadults$n_percent_aquatic_occ_logit)
  str(aquatic_nestlingsadults$n_percent_aquatic_occ_fam_logit)
  str(aquatic_nestlingsadults$n_simpson_psmelt)
  str(aquatic_nestlingsadults$Site)
  str(aquatic_nestlingsadults$n_Age)
  aquatic_nestlingsadults$n_Age <- as.factor(aquatic_nestlingsadults$n_Age) # change to factor
  str(aquatic_nestlingsadults$n_Age)
  str(aquatic_nestlingsadults$n_exp_treat)
  str(aquatic_nestlingsadults$ad_Mass)
  str(aquatic_nestlingsadults$ad_Flat_Wing)
  str(aquatic_nestlingsadults$n_brood_size_sampling_time)

  ## Plot the relationship between various metrics and percent aquatic, relative abundance. ----
    # This could also be done subbing in percent aquatic occurrence, percent aquatic occurrence family, and diversity.
    # We're looking to make sure there aren't any big abnormalities in predictor variables.

  # Plot mom's mass
  p_mass <- ggplot(aquatic_nestlingsadults) +
    geom_point((aes(x=ad_Mass, y=n_percent_aquatic_ra))) + 
    xlab("Female mass") +
    ylab("Proportion of nestling diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_mass
  # No major outliers or biologically unrealistic values.
  # One female does not have a third capture mass: 272107803. She was captured,
  # but her mass was not measured during the third capture.

  # Plot mom's wing
  p_wing <- ggplot(aquatic_nestlingsadults) +
    geom_point((aes(x=ad_Flat_Wing, y=n_percent_aquatic_ra))) + 
    xlab("Female flatwing") +
    ylab("Proportion of nestling diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_wing
  # No major outliers or biologically unrealistic values

  p_exp_treat <- ggplot (aquatic_nestlingsadults) +
    geom_boxplot((aes(x=n_exp_treat, y=n_percent_aquatic_ra))) + 
    xlab("Experimental Treatment") +
    ylab("Proportion of nestling diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_exp_treat
  
  # Plot brood size + percent aquatic
  p_brood <- ggplot(aquatic_nestlingsadults) +
    geom_boxplot((aes(x=factor(n_brood_size_sampling_time), y=n_percent_aquatic_ra))) + 
    xlab("Brood Size") +
    ylab("Proportion of nestling diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_brood
  
  ## Model for nestling alpha diversity ~ adult phenotype -------------------------------
  mod_div <- lmer(n_simpson_psmelt ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) + 
                  n_Age + (1|n_exp_treat) + (1|site_box_year),
                  data = aquatic_nestlingsadults)
  # Checking conditions
  hist(resid(mod_div))
  plot(mod_div)
  summary(mod_div, ddf = "Kenward-Roger")
  anova(mod_div, ddf = "Kenward-Roger")
  tab_model(mod_div, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_diversity/ad_phen_nestling_diet_tab_div.html"))
  # Examine pairwise differences between variables
  emmeans(mod_div, list(pairwise ~ Site), adjust = "tukey", ddf = "Kenward-Roger")
  
  # Plot alpha diversity ~ nestling site
  mod_div <- lmer(n_simpson_psmelt ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) + 
                  n_Age + (1|n_exp_treat) + (1|site_box_year),
                  data = aquatic_nestlingsadults)
  m_em <- as.data.frame(emmeans(mod_div, "Site", lmer.df = "Kenward-Roger"))
  m_eml <- pivot_longer(m_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
  
  p_site <- ggplot(aquatic_nestlingsadults, (aes(x=Site, y=n_simpson_psmelt))) +
    geom_boxplot(width = 0.25, position = position_nudge(x = -0.4), fill = "dodgerblue") + 
    geom_jitter(width = 0.15, alpha = 0.5, shape = 16, color = "dodgerblue", size = 2) +
    xlab("Site") +
    ylab("Simpson's diversity index") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14)) + 
    theme(axis.text.y = element_text(size = 14)) +
    scale_x_discrete(labels = c("Turkey_Hill" = "Turkey Hill", "Unit_1" = "Unit 1", 
                                "Unit_2" = "Unit 2", "Unit_4" = "Unit 4"))
  p_site2 <- p_site + geom_line(data = m_eml, mapping = aes(x = Site, y = y), color = "dodgerblue", size = 1) +
    geom_point(data = m_em, mapping = aes(x = Site, y = emmean), color = "dodgerblue", size = 4, shape = 23, fill = "dodgerblue")
  p_site2
  
  ## Model for nestling percent aquatic, relative abundance ~ adult phenotype -----
  mod_ra <- lmer(n_percent_aquatic_ra_logit ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) + 
               n_Age + (1|n_exp_treat) + (1|site_box_year),
               data = aquatic_nestlingsadults)
  # Checking conditions
  hist(resid(mod_ra))
  plot(mod_ra)
  summary(mod_ra, ddf = "Kenward-Roger")
  anova(mod_ra, ddf = "Kenward-Roger")
  tab_model(mod_ra, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_relativeabundance/ad_phen_nestling_diet_mod_tab_ra.html"))
  # Examine pairwise differences between variables
  emmeans(mod_ra, list(pairwise ~ n_Age), adjust = "tukey", ddf = "Kenward-Roger")
  
  ## Plot percent aquatic, relative abundance ~ nestling age -- FIGURE 3 -------------------
  mod_ra <- lmer(n_percent_aquatic_ra_logit ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) +
                 n_Age + (1|n_exp_treat) + (1|site_box_year),
                 data = aquatic_nestlingsadults)
  m_em <- as.data.frame(emmeans(mod_ra, "n_Age", lmer.df = "Kenward-Roger"))
  m_eml <- pivot_longer(m_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
  
  p_age <- ggplot(aquatic_nestlingsadults, (aes(x=as.factor(n_Age), y=n_percent_aquatic_ra))) +
    geom_boxplot(width = 0.25, position = position_nudge(x = -0.3), fill = "dodgerblue", outlier.shape = NA) + 
    geom_jitter(width = 0.12, alpha = 0.5, shape = 16, color = "dodgerblue", size = 2) +
    xlab("Nestling age (days)") +
    ylab(paste0("Proportion of nestling diet", "\n", "composed of aquatic insects")) +
    theme(axis.title.x = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
    theme(axis.title.y = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_age2 <- p_age + geom_line(data = m_eml, mapping = aes(x = n_Age, y = inv.logit(y, min = -0.0001, max = 1.0001)), position = position_nudge(x = 0.3), col = "dodgerblue", size = 1) +
    geom_point(data = m_em, mapping = aes(x = n_Age, y = inv.logit(emmean, min = -0.0001, max = 1.0001)), size = 4, shape = 23, position = position_nudge(x = 0.3), fill = "dodgerblue")
  p_age2
  ggsave(here("3_r_scripts/figs/figs_relativeabundance/Fig_3_nestling_age_ra.pdf"), p_age2, width = 5, height = 5, device = "pdf")

  ## Model for nestling percent aquatic, occurrence distinct taxonomic group ~ adult phenotype ----------
  mod_occ <- lmer(n_percent_aquatic_occ_logit ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) +
                n_Age + (1|n_exp_treat) + (1|site_box_year),
                data = aquatic_nestlingsadults)
  # Checking conditions
  hist(resid(mod_occ))
  plot(mod_occ)
  summary(mod_occ, ddf = "Kenward-Roger")
  anova(mod_occ, ddf = "Kenward-Roger")
  tab_model(mod_occ, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_occurrence/ad_phen_nestling_diet_mod_tab_occ.html"))
  # Examine pairwise differences between variables
  emmeans(mod_occ, list(pairwise ~ Site), adjust = "tukey", ddf = "Kenward-Roger")
  emmeans(mod_occ, list(pairwise ~ n_Age), adjust = "tukey", ddf = "Kenward-Roger")
  
  # Plot percent aquatic, occurrence distinct taxonomic group ~ site
  p_occ_site <- ggplot(aquatic_nestlingsadults) +
    geom_boxplot((aes(x=Site, y=n_percent_aquatic_occ)))
  p_occ_site
  
  # Plot percent aquatic, occurrence distinct taxonomic group ~ nestling age
  m_em <- as.data.frame(emmeans(mod_occ, "n_Age", lmer.df = "Kenward-Roger"))
  m_eml <- pivot_longer(m_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
  
  p_age <- ggplot(aquatic_nestlingsadults, (aes(x=as.factor(n_Age), y=n_percent_aquatic_occ))) +
    geom_boxplot(width = 0.25, position = position_nudge(x = -0.4), fill = "dodgerblue") + 
    geom_jitter(width = 0.15, alpha = 0.5, shape = 16, color = "dodgerblue", size = 2) +
    xlab("Nestling Age") +
    ylab("Proportion of nestling diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_age2 <- p_age + geom_line(data = m_eml, mapping = aes(x = n_Age, y = inv.logit(y, min = -0.0001, max = 1.0001)), col = "dodgerblue", size = 1) +
    geom_point(data = m_em, mapping = aes(x = n_Age, y = inv.logit(emmean, min = -0.0001, max = 1.0001)), color = "dodgerblue", size = 4, shape = 23, fill = "dodgerblue")
  p_age2
  
  ## Model for nestling percent aquatic, occurrence family ~ adult phenotype ----------
  mod_occ_fam <- lmer(n_percent_aquatic_occ_fam_logit ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) +
                    n_Age + (1|n_exp_treat) + (1|site_box_year),
                    data = aquatic_nestlingsadults)
  # Checking conditions
  hist(resid(mod_occ_fam))
  plot(mod_occ_fam) # weird
  summary(mod_occ_fam, ddf = "Kenward-Roger")
  anova(mod_occ_fam, ddf = "Kenward-Roger")
  tab_model(mod_occ_fam, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_occurrence_family/ad_phen_nestling_diet_mod_tab_occ_fam.html"))
  # Examine pairwise differences between variables
  emmeans(mod_occ_fam, list(pairwise ~ n_Age), adjust = "tukey", ddf = "Kenward-Roger")
  
  # Plot percent aquatic, occurrence family ~ nestling age
  m_em <- as.data.frame(emmeans(mod_occ_fam, "n_Age", lmer.df = "Kenward-Roger"))
  m_eml <- pivot_longer(m_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
  
  p_age <- ggplot(aquatic_nestlingsadults, (aes(x=as.factor(n_Age), y=n_percent_aquatic_occ_fam))) +
    geom_boxplot(width = 0.25, position = position_nudge(x = -0.4), fill = "dodgerblue") + 
    geom_jitter(width = 0.15, alpha = 0.5, shape = 16, color = "dodgerblue", size = 2) +
    xlab("Nestling Age") +
    ylab("Proportion of nestling diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_age2 <- p_age + geom_line(data = m_eml, mapping = aes(x = n_Age, y = inv.logit(y, min = -0.0001, max = 1.0001)), col = "dodgerblue", size = 1) +
    geom_point(data = m_em, mapping = aes(x = n_Age, y = inv.logit(emmean, min = -0.0001, max = 1.0001)), color = "dodgerblue", size = 4, shape = 23, fill = "dodgerblue")
  p_age2
  
  ## Examination of day 12 and day 15 patterns with individual ID --------------
  
  # We have individual bird ID for both day 12 and day 15 nestlings. Do we see
  # the same patterns if we use that information in the model?
  day12_15 <- aquatic_nestlingsadults[aquatic_nestlingsadults$n_Age == "12" |
                                        aquatic_nestlingsadults$n_Age == "15" ,]
  
  # Check structure of band
  str(day12_15$n_Individual_Band)
  day12_15$n_Individual_Band <- as.factor(day12_15$n_Individual_Band) # change band to factor
  str(day12_15$n_Individual_Band)
  
  mod_div <- lmer(n_simpson_psmelt ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) + 
                    n_Age + (1|n_Individual_Band) + (1|n_exp_treat) + (1|site_box_year),
                  data = day12_15)
  # Checking conditions
  hist(resid(mod_div))
  plot(mod_div)
  summary(mod_div, ddf = "Kenward-Roger")
  anova(mod_div)
  
  mod_ra <- lmer(n_percent_aquatic_ra_logit ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) + 
                   n_Age + (1|n_Individual_Band) + (1|n_exp_treat) + (1|site_box_year),
                 data = day12_15)
  # Checking conditions
  hist(resid(mod_ra))
  plot(mod_ra)
  summary(mod_ra, ddf = "Kenward-Roger")
  anova(mod_ra)
  
  mod_occ <- lmer(n_percent_aquatic_occ_logit ~ Site + scale(ad_Mass) + scale(ad_Flat_Wing) + scale(n_brood_size_sampling_time) +
                    n_Age + (1|n_Individual_Band) + (1|n_exp_treat) + (1|site_box_year),
                  data = day12_15)
  # Checking conditions
  hist(resid(mod_occ))
  plot(mod_occ)
  summary(mod_occ, ddf = "Kenward-Roger")
  anova(mod_occ)
  
  # Broadly see the same patterns.

# Does an adult female’s phenotype predict her own diet variation? -------------

  ## Prepare data frame for modeling and plotting ------------------------------
  aquatic_ad <- aquatic[aquatic$Adult_or_Nestling == "Adult" ,]
  aquatic_adF <- aquatic_ad[aquatic_ad$Sex == "F" ,]

  # Use gtools package to perform a logit transformation on percentage data (relative abundance)
  aquatic_adF$percent_aquatic_ra_logit <- gtools::logit(aquatic_adF$percent_aquatic_ra, min = -0.0001, max = 1.0001)
  
  # Use gtools package to perform a logit transformation on percentage data (occurrence)
  aquatic_adF$percent_aquatic_occ_logit <- gtools::logit(aquatic_adF$percent_aquatic_occ, min = -0.0001, max = 1.0001)
  
  # Use gtools package to perform a logit transformation on percentage data (occurrence family)
  aquatic_adF$percent_aquatic_occ_fam_logit <- gtools::logit(aquatic_adF$percent_aquatic_occ_fam, min = -0.0001, max = 1.0001)

  # Some third capture females do not have wing measurements, so import the wing measurements from their first capture
  # However, some of these first wing measurements are recorded during captures when there was not a
  # fecal sample taken. We need to go back into the raw data to extract what we need.
  captures <- read.csv("~/Dropbox/tres_database_data/Captures_Hormone_Bleeding_Blood_DNA_12.09.21.csv")
  captures <- captures[!is.na(captures$Exp_Year), ]
  captures <- captures[captures$Exp_Year == "2019" ,]
  captures <- captures[captures$Adult_or_Nestling == "Adult" ,]
  captures <- captures[captures$Sex == "F" ,]
  
  captures_wing <- captures[captures$Capture_Number == "1" ,]
  
  # There were some birds that were re-nests, and so they had multiple "first captures."
  # In the following code, I remove those birds' "first" captures that actually occurred at their re-nest nestbox.
  # 278172145, cap 1 renest = 6/13/19
  captures_wing <- captures_wing[!(captures_wing$Individual_Band == "278172145" & captures_wing$Capture_Date == "13-Jun-19") ,]
  # 278172146, cap 1 renest = 6/9/19
  captures_wing <- captures_wing[!(captures_wing$Individual_Band == "278172146" & captures_wing$Capture_Date == "09-Jun-19") ,]
  # 278172231, cap 1 renest = 6/21/19
  captures_wing <- captures_wing[!(captures_wing$Individual_Band == "278172231" & captures_wing$Capture_Date == "21-Jun-19") ,]
  # 281128404, cap 1 renest = 6/17/19
  captures_wing <- captures_wing[!(captures_wing$Individual_Band == "281128404" & captures_wing$Capture_Date == "17-Jun-19") ,]
  
  cols_needed <- as.vector(c("Individual_Band", "Flat_Wing"))
  captures_wing <- captures_wing[,cols_needed]
  
  # Delete column in aquatic_adF that has wing
  aquatic_adF <- subset(aquatic_adF, select = -c(Flat_Wing))
  
  # Add in new column with wing measurements from first captures
  aquatic_adF <- left_join(aquatic_adF, captures_wing)
  
  ## Add new variable for the treatment experienced at the time of capture  ------
  
  # For Unit 1 and Unit 2 adults, at first capture, they had not experienced any treatment.
  # Indicate this in new "treat_time_of_capture" variable.
  
  aquatic_adF$treat_time_of_capture <- NA
  
  for (i in 1:length(aquatic_adF$sampleID)){
    if (aquatic_adF$Capture_Number[i] == "1"){
      if (aquatic_adF$Site[i] == "Unit_1"){
      aquatic_adF$treat_time_of_capture[i] <- "Control"
      }
      if (aquatic_adF$Site[i] == "Unit_2"){
        aquatic_adF$treat_time_of_capture[i] <- "Control"
      }
      if (aquatic_adF$Site[i] == "Turkey_Hill"){
        aquatic_adF$treat_time_of_capture[i] <- aquatic_adF$exp_treat[i]
      }
      if (aquatic_adF$Site[i] == "Unit_4"){
        aquatic_adF$treat_time_of_capture[i] <- aquatic_adF$exp_treat[i]
      }
    }
  }

  # For Unit 1 and Unit 2 adults, at second capture, some of the birds had experienced one treatment.
  # "Pred" birds had the predation treatment, and "Cont" birds received no treatment and so should be treated
  # as "Control" birds at the point of second capture.
  # Indicate this in "treat_time_of_capture" variable.
  
  for (i in 1:length(aquatic_adF$exp_treat)){
    if (aquatic_adF$Capture_Number[i] == "2"){
      if (aquatic_adF$exp_treat[i] == "Control"){
        aquatic_adF$treat_time_of_capture[i] <- aquatic_adF$exp_treat[i]
      }
      else{
        if (aquatic_adF$Nest_Experiment[i] == "Color_Stress"){
          ttt <- strsplit(aquatic_adF$exp_treat[i], "_")
          aquatic_adF$treat_time_of_capture[i] <- ttt[[1]][3] # Extract first experimental treatment from string
        }
        if (aquatic_adF$Nest_Experiment[i] != "Color_Stress"){
          aquatic_adF$treat_time_of_capture[i] <- aquatic_adF$exp_treat[i]
        }
      }
    }
  }
  
  # Re-label "Cont" treatment birds as "Control" from second capture
  for (i in 1:length(aquatic_adF$exp_treat)){
    if (aquatic_adF$Capture_Number[i] == "2"){
      if (aquatic_adF$treat_time_of_capture[i] == "Cont"){
        aquatic_adF$treat_time_of_capture[i] <- "Control"
      }
    }
  }
      
  # For final capture (3), all treatments should be listed.
  for (i in 1:length(aquatic_adF$exp_treat)){
    if (aquatic_adF$Capture_Number[i] == "3"){
      aquatic_adF$treat_time_of_capture[i] <- aquatic_adF$exp_treat[i]
    }
  }
  
  # Check structure of each variable in the model to make sure it is correct
  str(aquatic_adF$Site)
  str(aquatic_adF$Flat_Wing)
  str(aquatic_adF$Mass)
  str(aquatic_adF$Capture_Number)
  aquatic_adF$Capture_Number <- as.factor(aquatic_adF$Capture_Number)
  str(aquatic_adF$Capture_Number)
  str(aquatic_adF$treat_time_of_capture)
  table(aquatic_adF$treat_time_of_capture)
  str(aquatic_adF$site_box_year)
  str(aquatic_adF$percent_aquatic_ra_logit)
  str(aquatic_adF$percent_aquatic_occ_logit)
  str(aquatic_adF$percent_aquatic_occ_fam_logit)
  str(aquatic_adF$simpson_psmelt)

  ## Plot the relationship between various metrics and percent aquatic, relative abundance. ----
  # This could also be done subbing in percent aquatic occurrence, occurrence family, and/or diversity.
  # We're looking to make sure there aren't any big abnormalities in predictor variables here.

  # Plot mass + percent aquatic
  p_mass <- ggplot(aquatic_adF, (aes(x=Mass, y=percent_aquatic_ra))) +
    geom_point() + 
    geom_smooth(method=lm , se=TRUE, alpha = 0.5) +
    xlab("Mass") +
    ylab("Proportion of adult diet composed of aquatic insect") +
    theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 16))
  p_mass
  # No major outliers or biologically unrealistic values
  
  # Plot wing + percent aquatic
  p_wing <- ggplot(aquatic_adF, (aes(x=Flat_Wing, y=percent_aquatic_ra))) +
    geom_point() + 
    geom_smooth(method=lm , se=TRUE) +
    xlab("Female flatwing") +
    ylab("Proportion of adult diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_wing
  # No major outliers or biologically unrealistic values

  # Plot experimental treatment + percent aquatic
  p_exp <- ggplot(aquatic_adF) +
    geom_boxplot((aes(x=treat_time_of_capture, y=percent_aquatic_ra))) + 
    xlab("Treatment") +
    ylab("Proportion of adult diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_exp
  
  # Plot mass by stage
  p_mass <- ggplot(aquatic_adF) +
    geom_boxplot((aes(x=Capture_Number, y=Mass))) + 
    xlab("capture number") +
    ylab("mass") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14))
  p_mass
  # This shows that mass declines from capture 1 to 2 to 3, so we should include all
  # capture levels in the models (rather than delineating between incubation and
  # provisioning only).
  
  ## Model for adult alpha diversity ~ adult phenotype -------------------------
  mod_div <- lmer(simpson_psmelt ~ Site + scale(Flat_Wing) + scale(Mass)*Capture_Number + (1|treat_time_of_capture) + (1|site_box_year),
                  data = aquatic_adF)
  
  # Checking conditions
  hist(resid(mod_div))
  plot(mod_div)
  summary(mod_div, ddf = "Kenward-Roger")
  anova(mod_div, ddf = "Kenward-Roger")
  tab_model(mod_div, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_diversity/ad_phen_ad_diet_tab_mod_div.html"))
  # Examine pairwise differences between variables
  emtrends(mod_div, pairwise ~ Capture_Number, var = "Mass", ddf = "Kenward-Roger")
  emmip(mod_div, Capture_Number ~ Mass, cov.reduce = range, ddf = "Kenward-Roger")
  
 
  
  # ## Plot adult diversity ~ mass + stage -- FIGURE 4 ---------------------------------------
  # commented out because this takes a long time to run
  # # Create the same model without scaling simpson for ease of plotting
  # mod_div <- lmer(simpson_psmelt ~ Site + Flat_Wing + Mass*Capture_Number + (1|exp_treat) + (1|site_box_year),
  #                 data = aquatic_adF, control = lmerControl(optimizer = "bobyqa"))
  # 
  # # sample from the fit model to get a 'posterior' distribution
  # # Increasing 'n' will result in more samples = smoother, but takes longer
  # post <- mvrnorm(n = 1e5, mu = fixef(mod_div), Sigma = vcov(mod_div))
  # post <- as.data.frame(post)
  # 
  # # create a vector of the values to evaluate the model at
  # # decreasing interval will result in more samples = smoother, but takes longer
  # r.cap1 <- seq(min(subset(aquatic_adF$Mass,aquatic_adF$Capture_Number=="1")),
  #               max(subset(aquatic_adF$Mass,aquatic_adF$Capture_Number=="1")), 0.0001)
  # r.cap2 <- seq(min(subset(aquatic_adF$Mass,aquatic_adF$Capture_Number=="2")),
  #               max(subset(aquatic_adF$Mass,aquatic_adF$Capture_Number=="2")), 0.0001)
  # r.cap3 <- seq(min(subset(aquatic_adF$Mass,aquatic_adF$Capture_Number=="3")),
  #               max(subset(aquatic_adF$Mass,aquatic_adF$Capture_Number=="3")), 0.0001)
  # 
  # # Add other variables
  # # evaluate the predicted Simpson's diversity index at each value of r
  # # this is getting the values for plotting best fit line
  # lh_fl_cap1 <- sapply(r.cap1, function(z)mean(post$"(Intercept)"+post$SiteUnit_2+post$Flat_Wing*mean(aquatic_adF$Flat_Wing)+post$Mass*z))
  # lh_fl_cap2 <- sapply(r.cap2, function(z)mean(post$"(Intercept)"+post$SiteUnit_2+post$Flat_Wing*mean(aquatic_adF$Flat_Wing)+post$"Capture_Number2"+post$Mass*z+post$"Mass:Capture_Number2"*z))
  # lh_fl_cap3 <- sapply(r.cap3, function(z)mean(post$"(Intercept)"+post$SiteUnit_2+post$Flat_Wing*mean(aquatic_adF$Flat_Wing)+post$"Capture_Number3"+post$Mass*z+post$"Mass:Capture_Number3"*z))
  # 
  # # evaluate the confidence interval (same as above but gets high and low)
  # lh_ci_cap1 <- sapply(r.cap1, function(z)HPDI(post$"(Intercept)"+post$SiteUnit_2+post$Flat_Wing*mean(aquatic_adF$Flat_Wing)+post$Mass*z))
  # lh_ci_cap2 <- sapply(r.cap2, function(z)HPDI(post$"(Intercept)"+post$SiteUnit_2+post$Flat_Wing*mean(aquatic_adF$Flat_Wing)+post$"Capture_Number2"+post$Mass*z+post$"Mass:Capture_Number2"*z))
  # lh_ci_cap3 <- sapply(r.cap3, function(z)HPDI(post$"(Intercept)"+post$SiteUnit_2+post$Flat_Wing*mean(aquatic_adF$Flat_Wing)+post$"Capture_Number3"+post$Mass*z+post$"Mass:Capture_Number3"*z))
  # 
  # # put these together in data frame
  # fit_cap1 <- data.frame(r = r.cap1,
  #                   mu = lh_fl_cap1,
  #                   lo = lh_ci_cap1[1, ],
  #                   hi = lh_ci_cap1[2, ])
  # 
  # fit_cap2 <- data.frame(r = r.cap2,
  #                        mu = lh_fl_cap2,
  #                        lo = lh_ci_cap2[1, ],
  #                        hi = lh_ci_cap2[2, ])
  # 
  # fit_cap3 <- data.frame(r = r.cap3,
  #                        mu = lh_fl_cap3,
  #                        lo = lh_ci_cap3[1, ],
  #                        hi = lh_ci_cap3[2, ])
  # 
  # cap_names <- c(
  #   `1` = "Early incubation",
  #   `2` = "Late incubation",
  #   `3` = "Provisioning"
  # )
  # aquatic_adF$Capture_Number <- as.factor(aquatic_adF$Capture_Number)
  # 
  # p_mass <- ggplot(aquatic_adF, (aes(x=Mass, y=simpson_psmelt, color = Capture_Number, fill = Capture_Number))) +
  #   geom_point() +
  #   scale_fill_manual(values = c("brown", "lightpink3", "lightpink1")) +
  #   scale_color_manual(name = "Breeding stage",
  #                      labels = c("Mid incubation", "Late incubation", "Provisioning"),
  #                      values = c("brown", "lightpink3", "lightpink1")) +
  #   guides(fill = "none") +
  #   xlab("Mass (g)") +
  #   ylab("Simpson's diversity index") +
  #   theme(axis.title = element_text(size = 16)) + theme(axis.text = element_text(size = 14)) +
  #   theme(legend.title = element_text(size = 16), legend.text = element_text (size = 14)) +
  #   geom_ribbon(data = fit_cap1, inherit.aes = FALSE,
  #               mapping = aes(x = r, ymin = lo, ymax = hi), fill = "brown", alpha = 0.5) +
  #   geom_line(data = fit_cap1, inherit.aes = FALSE,
  #             mapping = aes(x = r, y = mu), color = "brown", size = 2) +
  #   geom_ribbon(data = fit_cap2, inherit.aes = FALSE,
  #             mapping = aes(x = r, ymin = lo, ymax = hi), fill = "lightpink3", alpha = 0.5) +
  #   geom_line(data = fit_cap2, inherit.aes = FALSE,
  #             mapping = aes(x = r, y = mu), color = "lightpink3", size = 2) +
  #   geom_ribbon(data = fit_cap3, inherit.aes = FALSE,
  #               mapping = aes(x = r, ymin = lo, ymax = hi), fill = "lightpink1", alpha = 0.5) +
  #   geom_line(data = fit_cap3, inherit.aes = FALSE,
  #             mapping = aes(x = r, y = mu), color = "lightpink1", size = 2)
  # 
  # p_mass
  # ggsave(here("3_r_scripts/figs/figs_diversity/Fig_4_mom_mass_mom_div.pdf"), p_mass, width = 7, height = 5, device = "pdf")

  ## Model for adult percent aquatic, relative abundance ~ adult phenotype -----
  mod_ra <- lmer(percent_aquatic_ra_logit ~ Site + scale(Flat_Wing) + scale(Mass)*Capture_Number + (1|treat_time_of_capture) + (1|site_box_year),
               data = aquatic_adF)
  # Checking conditions
  hist(resid(mod_ra))
  plot(mod_ra)
  summary(mod_ra, ddf = "Kenward-Roger")
  anova(mod_ra, ddf = "Kenward-Roger")
  tab_model(mod_ra, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_relativeabundance/ad_phen_ad_diet_tab_mod_ra.html"))
  # Examine pairwise differences between variables
  emmeans(mod_ra, list(pairwise ~ Site), adjust = "tukey", ddf = "Kenward-Roger")
  
  # Plot for adult percent aquatic, relative abundance ~ site
  p_ra_site <- ggplot(aquatic_adF) +
    geom_boxplot((aes(x=Site, y=percent_aquatic_ra_logit)))
  p_ra_site

  ## Model for adult percent aquatic, occurrence distinct taxonomic group ~ adult phenotype -------------
  mod_occ <- lmer(percent_aquatic_occ_logit ~ Site + scale(Flat_Wing) + scale(Mass)*Capture_Number + (1|treat_time_of_capture) + (1|site_box_year),
                data = aquatic_adF)
  # Checking conditions
  hist(resid(mod_occ))
  plot(mod_occ)
  summary(mod_occ, ddf = "Kenward-Roger")
  anova(mod_occ, ddf = "Kenward-Roger")
  tab_model(mod_occ, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_occurrence/ad_phen_ad_diet_tab_mod_occ.html"))

  ## Model for adult percent aquatic, occurrence family ~ adult phenotype -------------
  mod_occ_fam <- lmer(percent_aquatic_occ_fam_logit ~ Site + scale(Flat_Wing) + scale(Mass)*Capture_Number + (1|treat_time_of_capture) + (1|site_box_year),
                  data = aquatic_adF)
  # Checking conditions
  hist(resid(mod_occ_fam))
  plot(mod_occ_fam)
  summary(mod_occ_fam, ddf = "Kenward-Roger")
  anova(mod_occ_fam, ddf = "Kenward-Roger")
  tab_model(mod_occ_fam, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_occurrence_family/ad_phen_ad_diet_tab_mod_occ_fam.html"))
  
# How do male and female diets differ? -----------------------------------------

  ## Prepare data frame for modeling and plotting ------------------------------
  aquatic_prov_adF <- aquatic[aquatic$Adult_or_Nestling == "Adult" & aquatic$Sex == "F" ,] # Select females
  aquatic_prov_adF <- aquatic_prov_adF[aquatic_prov_adF$Capture_Number == "3" ,] # Select females from third capture only (during provisioning)
  aquatic_prov_adM <- aquatic[aquatic$Adult_or_Nestling == "Adult" & aquatic$Sex == "M" ,] # Select males
  
  adults_compare <- rbind(aquatic_prov_adF, aquatic_prov_adM)
  adults_compare <- adults_compare[adults_compare$Site == "Unit_1" | adults_compare$Site == "Unit_2" ,]
  
  # Use gtools package to perform a logit transformation on percentage data
  adults_compare$percent_aquatic_ra_logit <- gtools::logit(adults_compare$percent_aquatic_ra, min = -0.0001, max = 1.0001)
  adults_compare$percent_aquatic_occ_logit <- gtools::logit(adults_compare$percent_aquatic_occ, min = -0.0001, max = 1.0001)
  adults_compare$percent_aquatic_occ_fam_logit <- gtools::logit(adults_compare$percent_aquatic_occ_fam, min = -0.0001, max = 1.0001)
  
  str(adults_compare$simpson_psmelt)
  str(adults_compare$percent_aquatic_ra_logit)
  str(adults_compare$percent_aquatic_occ_logit)
  str(adults_compare$percent_aquatic_occ_fam_logit)
  str(adults_compare$Sex)
  str(adults_compare$site_box_year)
  
  ## Model for alpha diversity ~ sex -------------------------------------------
  mod_div <- lmer(simpson_psmelt ~ Sex + (1|exp_treat) + (1|site_box_year), data = adults_compare)
  summary(mod_div, ddf = "Kenward-Roger")
  anova(mod_div, ddf = "Kenward-Roger")
  # Checking conditions
  hist(resid(mod_div))
  plot(mod_div)
  
  ## Model for percent aquatic, relative abundance ~ Sex -----------------------
  mod_ra <- lmer(percent_aquatic_ra_logit ~ Sex + (1|exp_treat) + (1|site_box_year), data = adults_compare) 
  summary(mod_ra, ddf = "Kenward-Roger")
  anova(mod_ra, ddf = "Kenward-Roger")
  # Checking conditions
  hist(resid(mod_ra))
  plot(mod_ra)
  
  # Plot percent aquatic, relative abundance ~ Sex
  p_ra <- ggplot(adults_compare) +
    geom_boxplot(width = 0.5, (aes(x=Sex, y=percent_aquatic_ra))) + 
    facet_wrap(~ Site) +
    xlab("Sex") +
    ylab("Proportion of adult diet composed of aquatic insects") +
    theme_classic() +
    theme(axis.title = element_text(size = 26)) + theme(axis.text.x = element_text(size = 26), 
                                                        axis.text.y = element_text(size = 18)) +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 24)) +
    scale_fill_manual(values=c("orange", "gray"))
  p_ra

  ## Model for percent aquatic, occurrence distinct taxonomic group ~ Sex -------------------------------
  mod_occ <- lmer(percent_aquatic_occ_logit ~ Sex + (1|exp_treat) + (1|site_box_year), data = adults_compare)
  summary(mod_occ, ddf = "Kenward-Roger")
  anova(mod_occ, ddf = "Kenward-Roger")
  # Checking conditions
  hist(resid(mod_occ))
  plot(mod_occ)
  
  ## Model for percent aquatic, occurrence family ~ Sex -------------------------------
  mod_occ_fam <- lmer(percent_aquatic_occ_fam_logit ~ Sex + (1|exp_treat) + (1|site_box_year), data = adults_compare)
  summary(mod_occ_fam, ddf = "Kenward-Roger")
  anova(mod_occ_fam, ddf = "Kenward-Roger")
  # Checking conditions
  hist(resid(mod_occ_fam))
  plot(mod_occ_fam)

# Do adults and nestlings differ in their diet variation? ----------------------
  
  ## Prepare data frame for modeling and plotting ------------------------------
  
  # We will use the aquatic data frame, but use only adult samples from birds 
  # captured during provisioning.
  aquatic_prov_adF <- aquatic[aquatic$Adult_or_Nestling == "Adult" & aquatic$Sex == "F" ,] # Select females
  aquatic_prov_adF <- aquatic_prov_adF[aquatic_prov_adF$Capture_Number == "3" ,] # Select females from third capture only (during provisioning)
  aquatic_prov_adF$bird_category <- paste(aquatic_prov_adF$Adult_or_Nestling, aquatic_prov_adF$Sex, sep = "_") # Classify as adult females
  aquatic_prov_adM <- aquatic[aquatic$Adult_or_Nestling == "Adult" & aquatic$Sex == "M" ,] # Select males
  aquatic_prov_adM$bird_category <- paste(aquatic_prov_adM$Adult_or_Nestling, aquatic_prov_adM$Sex, sep = "_") # Classify as adult males
  aquatic_prov_nestlings <- aquatic[aquatic$Adult_or_Nestling == "Nestling" ,] # Select nestlings
  aquatic_prov_nestlings$bird_category <- "Nestling" # Classify as nestlings
  
  # Bind together all of these samples from provisioning
  aquatic_prov <- rbind(aquatic_prov_adF, aquatic_prov_adM, aquatic_prov_nestlings)
  
  # Use gtools package to perform a logit transformation on percentage data
  aquatic_prov$percent_aquatic_ra_logit <- gtools::logit(aquatic_prov$percent_aquatic_ra, min = -0.0001, max = 1.0001)
  aquatic_prov$percent_aquatic_occ_logit <- gtools::logit(aquatic_prov$percent_aquatic_occ, min = -0.0001, max = 1.0001)
  aquatic_prov$percent_aquatic_occ_fam_logit <- gtools::logit(aquatic_prov$percent_aquatic_occ_fam, min = -0.0001, max = 1.0001)
  
  # Check the structure of each variable before modeling
  str(aquatic_prov$bird_category)
  str(aquatic_prov$site_box_year)
  str(aquatic_prov$percent_aquatic_ra_logit)
  str(aquatic_prov$percent_aquatic_occ_logit)
  str(aquatic_prov$percent_aquatic_occ_fam_logit)
  str(aquatic_prov$simpson_psmelt)
  str(aquatic_prov$exp_treat)

  ## Model for alpha diversity ~ bird category ---------------------------------
  mod_div <- lmer(simpson_psmelt ~ bird_category + (1|exp_treat) + (1|site_box_year), data = aquatic_prov)
  summary(mod_div, ddf = "Kenward-Roger")
  anova(mod_div, ddf = "Kenward-Roger")
  tab_model(mod_div, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_diversity/bird_category_div.html"))
  # Checking conditions
  hist(resid(mod_div))
  plot(mod_div)

  ## Plot alpha diversity ~ bird category
  theme_set(theme_classic())
  
  m_em <- as.data.frame(emmeans(mod_div, c("bird_category"), lmer.df = "Kenward-Roger"))
  m_eml <- pivot_longer(m_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
  
  colors <- c()
  p_div <- ggplot(aquatic_prov, aes(x=bird_category, y=simpson_psmelt, fill = bird_category)) +
    geom_boxplot(aes(fill = factor(bird_category)), width = 0.25, position = position_nudge(x = -0.4)) +
    geom_jitter(aes(color = factor(bird_category)), width = 0.15, alpha = 0.5, shape = 16, size = 2) +
    xlab("Age") +
    ylab("Simpson's diversity index") +
    theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14)) +
    theme(legend.position="none") +
    scale_color_manual(values = c("Adult" = "lightpink", "Nestling" = "dodgerblue")) +
    scale_fill_manual(values = c("Adult" = "lightpink", "Nestling" = "dodgerblue"))
  p_div2 <- p_div + geom_line(data = m_eml, mapping = aes(x = bird_category, y = y, color = bird_category), size = 1) +
    geom_point(data = m_em, mapping = aes(x = bird_category, y = emmean), size = 4, shape = 23)
  p_div2

  ## Model for percent aquatic, relative abundance ~ bird category -------------
  mod_ra <- lmer(percent_aquatic_ra_logit ~ bird_category + (1|exp_treat) + (1|site_box_year), data = aquatic_prov)
  summary(mod_ra, ddf = "Kenward-Roger")
  anova(mod_ra, ddf = "Kenward-Roger")
  tab_model(mod_ra, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_relativeabundance/bird_category_ra.html"))
  # Checking conditions
  hist(resid(mod_ra))
  plot(mod_ra)
  # Post hoc pairwise comparisons
  emmeans(mod_ra, list(pairwise ~  bird_category), adjust = "tukey", ddf = "Kenward-Roger")
  
  # Plot percent aquatic, relative abundance ~ bird category
  theme_set(theme_classic())
  
  m_em <- as.data.frame(emmeans(mod_ra, c("bird_category"), lmer.df = "Kenward-Roger"))
  m_eml <- pivot_longer(m_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
  
  colors <- c()
  p_ra <- ggplot(aquatic_prov, aes(x=bird_category, y=inv.logit(percent_aquatic_ra_logit, min = -0.0001, max = 1.0001), fill = bird_category)) +
    geom_boxplot(aes(fill = factor(bird_category)), width = 0.25, position = position_nudge(x = -0.4)) +
    geom_jitter(aes(color = factor(bird_category)), width = 0.15, alpha = 0.5, shape = 16, size = 2) +
    xlab("Age") +
    ylab("Proportion of diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 14)) + theme(axis.text.x = element_text(size = 14)) +
    theme(legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 14)) +
    theme(legend.position="none") +
    scale_color_manual(values = c("Adult_F" = "lightpink", "Adult_M" = "orange", "Nestling" = "dodgerblue")) +
    scale_fill_manual(values = c("Adult_F" = "lightpink", "Adult_M" = "orange", "Nestling" = "dodgerblue"))
  p_ra2 <- p_ra + geom_line(data = m_eml, mapping = aes(x = bird_category, y = inv.logit(y, min = -0.0001, max = 1.0001), color = bird_category), size = 1) +
    geom_point(data = m_em, mapping = aes(x = bird_category, y = inv.logit(emmean, min = -0.0001, max = 1.0001)), size = 4, shape = 23)
  p_ra2

  ## Model for percent aquatic, occurrence distinct taxonomic group ~ bird category  -----
  mod_occ <- lmer(percent_aquatic_occ_logit ~ bird_category + (1|exp_treat) + (1|site_box_year), data = aquatic_prov)
  summary(mod_occ, ddf = "Kenward-Roger")
  anova(mod_occ, ddf = "Kenward-Roger")
  # Checking conditions
  hist(resid(mod_occ))
  plot(mod_occ)
  tab_model(mod_occ, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_occurrence/bird_category_occ.html"))
  # Post hoc pairwise comparisons
  emmeans(mod_occ, list(pairwise ~  bird_category), adjust = "tukey", ddf = "Kenward-Roger")
  
  # Plot percent aquatic, occurrence distinct taxonomic group ~ bird category
  
  theme_set(theme_classic())
  
  m_em <- as.data.frame(emmeans(mod_occ, c("bird_category"), lmer.df = "Kenward-Roger"))
  m_eml <- pivot_longer(m_em, cols = c("lower.CL", "upper.CL"), values_to = "y")
  
  colors <- c()
  p_occ <- ggplot(aquatic_prov, aes(x=bird_category, y=inv.logit(percent_aquatic_occ_logit, min = -0.0001, max = 1.0001), fill = bird_category)) +
    geom_boxplot(aes(fill = factor(bird_category)), width = 0.25, position = position_nudge(x = -0.4)) +
    geom_jitter(aes(color = factor(bird_category)), width = 0.15, alpha = 0.5, shape = 16, size = 2) +
    xlab("Age") +
    ylab("Proportion of diet composed of aquatic insects") +
    theme(axis.title = element_text(size = 14)) + theme(axis.text.x = element_text(size = 14)) +
    theme(legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
    theme(strip.text = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 14)) +
    theme(legend.position="none") +
    scale_color_manual(values = c("Adult_F" = "lightpink", "Adult_M" = "orange", "Nestling" = "dodgerblue")) +
    scale_fill_manual(values = c("Adult_F" = "lightpink", "Adult_M" = "orange", "Nestling" = "dodgerblue"))
  p_occ2 <- p_occ + geom_line(data = m_eml, mapping = aes(x = bird_category, y = inv.logit(y, min = -0.0001, max = 1.0001), color = bird_category), size = 1) +
    geom_point(data = m_em, mapping = aes(x = bird_category, y = inv.logit(emmean, min = -0.0001, max = 1.0001)), size = 4, shape = 23)
  p_occ2

  ## Model for percent aquatic, occurrence family ~ bird category  ---------------------
  mod_occ_fam <- lmer(percent_aquatic_occ_fam_logit ~ bird_category + (1|exp_treat) + (1|site_box_year), data = aquatic_prov)
  summary(mod_occ_fam, ddf = "Kenward-Roger")
  anova(mod_occ_fam, ddf = "Kenward-Roger")
  # Checking conditions
  hist(resid(mod_occ_fam))
  plot(mod_occ_fam)
  tab_model(mod_occ_fam, p.val = "kr", show.df = TRUE, file = here("3_r_scripts/model_outputs/model_outputs_occurrence_family/bird_category_occ_fam.html"))
  # Post hoc pairwise comparisons
  emmeans(mod_occ_fam, list(pairwise ~  bird_category), adjust = "tukey", ddf = "Kenward-Roger")
  
  # Plot for percent aquatic, occurrence family ~ bird category
  p_occ <- ggplot(aquatic_prov) +
    geom_boxplot((aes(x=bird_category, y=percent_aquatic_occ_fam)))
  p_occ