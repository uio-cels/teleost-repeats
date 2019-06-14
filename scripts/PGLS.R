#install.packages('caper')
library('caper')

#### Data frame 
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.tre')
data <- read.csv('DataFrame/Repeats.tsv', sep = '\t', row.names = 'ALIAS', na.strings = 'nan')
data['ALIAS'] <- rownames(data)
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)

### genome size, dropping outliers
fishTree <- drop.tip(fishTree, 'fish_107') # outlier 
fishTree <- drop.tip(fishTree, 'fish_109') # outlier 
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)

# best model
mod.pgls <- pgls(bp / 1e6 ~ TEs + log10(strs_perc) + brood1_cat, data = compdata, lambda = 'ML')
summary(mod.pgls)


#### strs vs. aquatic habitat
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.tre')
data <- read.csv('DataFrame/Repeats.tsv', sep = '\t', row.names = 'ALIAS', na.strings = 'nan')
data['ALIAS'] <- rownames(data)
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)

mod.pgls <- pgls(log10(strs_perc) ~ water, data = compdata, lambda = 'ML')
summary(mod.pgls)

# robust to removing Ovalentarians?
ova = c('fish_137', 'fish_130', 'fish_122', 'fish_129', 'fish_136', 'fish_114', 'fish_135', 'fish_141', 'fish_112', 'fish_132', 'fish_140', 'fish_125', 'fish_126', 'fish_143', 'fish_106', 'fish_87', 'fish_88', 'fish_90', 'fish_101')
fishTree <- drop.tip(fishTree, ova)
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)
summary(pgls(log10(strs_perc) ~ water, data = compdata, lambda = 'ML'))
# yes

#### TEs vs. STRs
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.tre')
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)
mod <- pgls(log10(strs_perc) ~ TEs, data = compdata, lambda = 'ML')
summary(mod)
# No correlation

#### assembly size, individual variables
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.tre')
data <- read.csv('DataFrame/Repeats.tsv', sep = '\t', row.names = 'ALIAS', na.strings = 'nan')
data['ALIAS'] <- rownames(data)
fishTree <- drop.tip(fishTree, 'fish_107') # outlier 
fishTree <- drop.tip(fishTree, 'fish_109') # outlier 
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)

for (var in c('nTEs', 'TEs', 'Total', 'num_strs_unit10', 'strs_perc', 'log_strs')) {
  mod <- pgls(as.formula(paste("bp ~", var)), data = compdata, lambda = 'ML')
  print(summary(mod))
  print(mod$model$coef)
}

#### Diversification rates
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.families_stem.tre')
groupdata <- read.csv('DataFrame/RepeatsGrouped.tsv', sep = '\t', row.names = 'FAMILY')
groupdata['FAMILY'] <- rownames(groupdata)
compdata <- comparative.data(fishTree, groupdata, 'FAMILY', na.omit=F)

#### Full model
for (var in c('TEs', 'Total', 'strs_perc', 'log_strs')) {
  mod <- pgls(as.formula(paste("div_rate_e5_family ~", var)), data = compdata, lambda = 'ML')
  print(summary(mod))
  print(mod$model$coef)
}

