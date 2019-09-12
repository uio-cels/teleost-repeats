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

# full model
mod.pgls <- pgls(bp / 1e6 ~ TEs + 
                 log10(strs_perc) + 
                 brood1_cat + 
                 water +
                 busco_complete_single +
                 n50_contig,
                 data = compdata, lambda = 'ML')
summary(mod.pgls)

#### TEs vs all interspersed
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.tre')
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)
mod <- pgls(TEs ~ Total, data = compdata, lambda = 'ML' )
summary(mod)

#### TEs vs. STRs
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.tre')
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)
mod <- pgls(log10(strs_perc) ~ TEs, data = compdata, lambda = 'ML')
summary(mod)
# No correlation

#### genome assembly size, individual variables
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.tre')
data <- read.csv('DataFrame/Repeats.tsv', sep = '\t', row.names = 'ALIAS', na.strings = 'nan')
data['ALIAS'] <- rownames(data)
fishTree <- drop.tip(fishTree, 'fish_107') # outlier 
fishTree <- drop.tip(fishTree, 'fish_109') # outlier 
compdata <- comparative.data(fishTree, data, 'ALIAS', na.omit=F)

for (var in c('nTEs', 'num_strs_unit10', 'TEs', 'log_strs')) {
  mod <- pgls(as.formula(paste("bp / 1e6 ~", var)), data = compdata, lambda = 'ML')
  print(summary(mod))
  print(mod$model$coef)
}

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

#### Diversification rates
fishTree <- read.tree('Trees/ALIASES_NEWICK_101g_nucl_conc_28bfossils.combined.families_stem.tre')
groupdata <- read.csv('DataFrame/RepeatsGrouped.tsv', sep = '\t', row.names = 'FAMILY')
groupdata['FAMILY'] <- rownames(groupdata)
compdata <- comparative.data(fishTree, groupdata, 'FAMILY', na.omit=F)

#### Full model
summary(pgls(div_rate_e5_family ~ TEs + log_strs + I(bp / 1e6), data = compdata, lambda = 'ML'))
# Testing div rates with different extinction rates
summary(pgls(div_rate_e0_family ~ TEs + log_strs + I(bp / 1e6), data = compdata, lambda = 'ML'))
summary(pgls(div_rate_e9_family ~ TEs + log_strs + I(bp / 1e6), data = compdata, lambda = 'ML'))
# very similar

#### Full model, total interspersed
summary(pgls(div_rate_e5_family ~ Total + log_strs + I(bp / 1e6), data = compdata, lambda = 'ML'))
summary(pgls(div_rate_e0_family ~ Total + log_strs + I(bp / 1e6), data = compdata, lambda = 'ML'))
summary(pgls(div_rate_e9_family ~ Total + log_strs + I(bp / 1e6), data = compdata, lambda = 'ML'))
# very similar
