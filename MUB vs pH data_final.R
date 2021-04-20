library(plyr)
library(dplyr)
library(car)
library(tidyr)
library(ggplot2)
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

#### Figure 1: pH changes with increasing [MUB] ####
pal.5 <- c("#0E6B9A", "#35BDDB", "#ACDEDE", "#D59D4F", "#FF7F00")
pal.3 <- c("#FF0000", "#FBAE00", "#00A58B")
summary.data$Soil.ID <- mapvalues(summary.data$Soil.ID, from=c("Soil 1", "Soil 2", "Soil 3", "Soil 4"), to = c("Soil 1 (pH = 8.6)", "Soil 2 (pH = 7.8)", "Soil 3 (pH = 5.1)", "Soil 4 (pH = 4.4)"))

Fig.1a <- ggplot(data=summary.data, aes(x=MUB.conc, y=pH.diff0_mean, group=Target.pH, color=Target.pH)) + geom_abline(slope=0, intercept=0, color="black", linetype="dashed") + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1)) + facet_grid(Substrate.conc ~ Soil.ID) + theme_bw() + ylim(-1.75,1.5) + labs(y="Difference from Target pH\n(Assay pH - Target pH at t = 0 min)", x="MUB Concentration") + theme(legend.position="top", legend.box = "horizontal") + scale_color_manual(values=pal.5, name="Target pH")

Fig.1b <- ggplot(data=summary.data, aes(x=pH.change_mean, fill=MUB.conc)) + geom_density(color=NA, alpha=0.5) + theme_bw() + scale_fill_manual(values=pal.3, name="MUB Concentration") + labs(y="Density", x="pH change over time\n(60 min - 0 min)") + theme(legend.position="top", legend.box = "horizontal") + geom_vline(xintercept=0, color="black", linetype="dashed") +xlim(-0.5,0.3)

Figure1 <- ggarrange(Fig.1a, Fig.1b, labels=c("a", "b"), heights=c(1.2, 0.75), ncol=1, nrow=2)
Figure1
ggsave("Figure 1.png", width=7, height=7)

#### STATS ####
library(lmerTest)
pH.diff.lmer <- lmer(pH.diff60 ~ Substrate.conc*MUB.conc*Target.pH + (1|Soil.ID), REML=FALSE, data=MUB.pH.data)
plot(resid(pH.diff.lmer))
shapiro.test(resid(pH.diff.lmer))
anova(pH.diff.lmer) # use only with lmerTest attached
detach("package:lmerTest", unload=TRUE)