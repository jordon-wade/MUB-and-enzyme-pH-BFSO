library(plyr)
library(dplyr)
library(car)
library(tidyr)
library(ggplot2)
library(nationalparkcolors)
library(lme4)

#### READ IN DATA AND MANIPULATE ####
setwd("~/Dropbox")
MUB.pH.data <- read.csv("Wade Margenot/Water vs buffer MS/MUB and pH/MUB vs pH data.csv")
MUB.pH.data <- MUB.pH.data[-100,]
MUB.pH.data$pH.change <- MUB.pH.data$pH.t60 - MUB.pH.data$pH.t0
MUB.pH.data$pH.diff0 <- MUB.pH.data$pH.t0 - MUB.pH.data$Target.pH
MUB.pH.data$pH.diff60 <- MUB.pH.data$pH.t60 - MUB.pH.data$Target.pH
summary.data <- as.data.frame(MUB.pH.data %>% dplyr::group_by(Target.pH, Soil.ID, MUB.conc, soil.pH, Substrate.conc) %>% dplyr::summarize_at(vars(pH.t0:pH.diff60), funs(mean, sd, se=sd(.)/sqrt(n()))))
summary.data$Target.pH <- factor(summary.data$Target.pH)
MUB.pH.data$MUB.conc <- ordered(MUB.pH.data$MUB.conc, levels=c("1x", "2x", "4x"))
MUB.pH.data$Target.pH <- factor(MUB.pH.data$Target.pH)

## Define color palettes ##
pal5 <- park_palette("Everglades")
pal3 <- park_palette("Badlands", 3)

#### Figure 1: pH changes with increasing [MUB] ####
pal.5 <- c("#0E6B9A", "#35BDDB", "#ACDEDE", "#D59D4F", "#FF7F00")
pal.3 <- c("#FF0000", "#FBAE00", "#00A58B")
summary.data$Soil.ID <- mapvalues(summary.data$Soil.ID, from=c("Soil 1", "Soil 2", "Soil 3", "Soil 4"), to = c("Soil 1 (pH = 8.6)", "Soil 2 (pH = 7.8)", "Soil 3 (pH = 5.1)", "Soil 4 (pH = 4.4)"))

Fig.1a <- ggplot(data=summary.data, aes(x=MUB.conc, y=pH.diff0_mean, group=Target.pH, color=Target.pH)) + geom_abline(slope=0, intercept=0, color="black", linetype="dashed") + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1)) + facet_grid(Substrate.conc ~ Soil.ID) + theme_bw() + ylim(-1.75,1.5) + labs(y="Difference from Target pH\n(Assay pH - Target pH at t = 0 min)", x="MUB Concentration") + theme(legend.position="top", legend.box = "horizontal") + scale_color_manual(values=pal.5, name="Target pH")

Fig.1b <- ggplot(data=summary.data, aes(x=pH.change_mean, fill=MUB.conc)) + geom_density(color=NA, alpha=0.5) + theme_bw() + scale_fill_manual(values=pal.3, name="MUB Concentration") + labs(y="Density", x="pH change over time\n(60 min - 0 min)") + theme(legend.position="top", legend.box = "horizontal") + geom_vline(xintercept=0, color="black", linetype="dashed") +xlim(-0.5,0.3)

Figure1 <- ggarrange(Fig.1a, Fig.1b, labels=c("a", "b"), heights=c(1.2, 0.75), ncol=1, nrow=2)
Figure1
ggsave("Figure 1.png", width=7, height=7)

#### STATS (just in case) ####
library(lmerTest)
pH.diff.lmer <- lmer(pH.diff60 ~ Substrate.conc*MUB.conc*Target.pH + (1|Soil.ID), REML=FALSE, data=MUB.pH.data)
plot(resid(pH.diff.lmer))
shapiro.test(resid(pH.diff.lmer))
anova(pH.diff.lmer) # use only with lmerTest attached
detach("package:lmerTest", unload=TRUE)


### Figure S1: pH drift separate by substrate level
FigureS1 <- ggplot(data=summary.data, aes(x=pH.change_mean, fill=MUB.conc)) + geom_density(color=NA, alpha=0.5) + theme_bw() + scale_fill_manual(values=pal.3, name="MUB Concentration") + labs(y="Density", x="pH change over time\n(60 min - 0 min)") + theme(legend.position="top", legend.box = "horizontal") + geom_vline(xintercept=0, color="black", linetype="dashed") +xlim(-0.5,0.3) + facet_grid(Substrate.conc ~ .)
FigureS1
ggsave("Figure S1.png", width=3.5, height=5)

### Figure: response to reviewers on activity and pH (from pH optima paper) ###
### Load and wrangle the modeled data ###
modeled.data.all <- read.csv("/Users/TheDankness/Dropbox/Wade Margenot/pH optimization/Data and Quantitative Analyses/Standardized data for R.csv")
modeled.data.all$Rel.Act_est <- modeled.data.all$Activity.std_est*100
modeled.data.all$Rel.Act_min <- modeled.data.all$Activity.std_min*100
modeled.data.all$Rel.Act_max <- modeled.data.all$Activity.std_max*100

modeled.data.all$Soil.id <- factor(modeled.data.all$Soil.id, levels=c("Soil 1", "Soil 2", "Soil 3", "Soil 4", "Soil 5", "Soil 6", "Soil 7", "Soil 8", "Soil 9", "Soil 10", "Soil 11", "Soil 12", "Soil 13", "Soil 14", "Soil 15", "Soil 16", "Soil 17", "Soil 18", "Soil 19", "Soil 20", "Soil 21", "Soil 22", "Soil 23", "Soil 24", "Soil 25", "Soil 26", "Prunus dulcis", "Aspergillus niger", "Escherichia coli", "Wheat germ"))

# Unique to this dataset #
df1 <- subset(modeled.data.all, Soil.id=="Soil 1"| Soil.id=="Soil 2"|Soil.id=="Soil 14"|Soil.id=="Soil 17")
df2 <- subset(df1, Enzyme=="Phos")
df2$Soil.id <- mapvalues(df2$Soil.id, from=c("Soil 17", "Soil 1", "Soil 14", "Soil 2"), to = c("Soil 1 (pH = 8.6)", "Soil 2 (pH = 7.8)", "Soil 3 (pH = 5.1)", "Soil 4 (pH = 4.4)"))
df2$Soil.id <- factor(df2$Soil.id, levels=c("Soil 1 (pH = 8.6)", "Soil 2 (pH = 7.8)", "Soil 3 (pH = 5.1)", "Soil 4 (pH = 4.4)"))

pH.vec <- as.data.frame(c(8.6, 7.8, 5.1, 4.4))
colnames(pH.vec) <- "soil.pH"
Response.fig <- ggplot(data=df2, aes(x=pH_est, y=Rel.Act_est, color="Enzyme")) + geom_ribbon(aes(ymin=Rel.Act_min, ymax=Rel.Act_max, color=NULL), alpha=0.3) + geom_point(size=0.25) + facet_wrap(.~Soil.id, nrow=2, ncol=2) + scale_color_manual(values="#DC3220") + theme_bw() + theme(legend.position="none") + labs(y="Relative Activity (%)", x="Buffer pH") + scale_y_continuous(breaks=seq(0,100, by=25)) + scale_x_continuous(breaks=seq(3,13,2)) + coord_cartesian(ylim=c(0,100), xlim=c(3,12))
Response.fig
