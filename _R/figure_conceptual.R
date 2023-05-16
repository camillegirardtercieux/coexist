
#################
### Libraries ###
#################

library(plyr)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggplot2)
library(gridExtra)
library(scales)

#################
### Read data ###
#################

mortality <- "prop"
fecundity <- "abund"

seed_example <- 9
seed_r_example <- 1

species <- 1

load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed_example}_{seed_r_example}_sites.RData")))
load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed_example}_{seed_r_example}_perf_E_Sp.RData")))

n_axes <- ncol(sites)
Conceptual_fig <- cbind(perf_E_Sp[,1], sites)
colnames(Conceptual_fig) <- c("Perf", c(sprintf("Env_%d", 1:n_axes)))

load(here::here("outputs", "Comparison", "Species_all.RData"))

change_label <- paste0("X", 1:n_axes)

colnames_long <- paste0("Env_", 1:n_axes)

##################
## Build figure ##
##################

### Panel A: Species response to environmental dimensions (perfect vs partial knowledge)

pa <- Conceptual_fig %>%
  pivot_longer(colnames_long, names_to = "environment", values_to = "response") %>%
  mutate(environment = factor(environment, levels = colnames_long)) %>%
  ggplot(aes(x = response, y = Perf))+
  geom_point(size = 0.5, alpha = 0.5, aes(color = environment), show.legend=F)+
  scale_colour_viridis_d(option= "magma", begin=0, end=0.9090909)+
  facet_wrap(vars(environment), ncol=5, labeller = as_labeller(c(Env_=change_label))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), data = ~ subset(., environment %in% c("Env_6")), colour = "grey60", size = 1, fill = NA, inherit.aes = F) +
  geom_smooth(method ="lm", formula = y~poly(x,2), se=F, aes(color = environment), show.legend=F)+
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0", sub("^0+", "", x)))+
  xlab("Environmental dimensions")+
  ylab("Performance")+
  #ggtitle(paste("(",letters[1],")",sep="")) +
  theme_classic() +
  theme(legend.position="none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))


### Panel B: uIV and sIV in species response

## Chose environmental dimension (here: Env_6)

## Envelope showing individual variability

# Linear model with polynomial
mod <- lm(Perf ~ 1 + Env_6 +I(Env_6^2), data=Conceptual_fig)
V <- var(mod$residuals)
# Predictions
nseq <- 100
Env_6_seq <- seq(min(Conceptual_fig$Env_6), max(Conceptual_fig$Env_6), length.out=nseq)
nsamp <- 1000
posterior <- matrix(NA, nrow=nsamp, ncol=nseq)
set.seed(1234)
variability <- rnorm(nsamp, 0, sqrt(V))
beta <- matrix(mod$coefficients)
X <- cbind(1, Env_6_seq, Env_6_seq^2)
for (i in 1:nseq) {
  posterior[, i] <- c(X[i,] %*% beta) + variability
}
# Credible interval
ci <- apply(posterior, 2, quantile, c(0.025, 0.975))
ci_df <- data.frame(t(ci))
names(ci_df) <- c("q_min", "q_max")
ci_df$Env_6 <- Env_6_seq
ci_df$Perf <- rep(0, nrow(ci_df)) # Can be any value

# Plot: min and max for Perf
y_min <- floor(min(Conceptual_fig$Perf))
y_max <- ceiling(max(Conceptual_fig$Perf))

pb <- Conceptual_fig %>%
  ggplot(aes(x=Env_6, y=Perf)) +
  geom_point(size=0.5, alpha=0.5) +
  geom_smooth(method="lm", formula=y~poly(x,2), se=F, col="#00BFC4") +
  coord_cartesian(xlim=c(0, 1), ylim=c(y_min, y_max)) +
  geom_ribbon(data=ci_df, aes(ymin=q_min, ymax=q_max), fill = "grey70", alpha=0.35) +
  scale_x_continuous(breaks=seq(0, 1, length.out=5)) +
  scale_y_continuous(breaks=seq(y_min, y_max, length.out=4),
                     labels = scales::number_format(accuracy = 0.01)) +
  xlab("Environment X6") +
  ylab("Performance") +
  theme_classic() +
  theme(legend.position="none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        plot.margin = margin(r = 3,
                             l = 0))


## Distribution of performance for Env_6=1 (last value)

P_seq <- seq(y_min, y_max, length.out=100)
dens <- dnorm(P_seq, mean=c(X[nseq,] %*% beta), sd=sqrt(V))
dist_df <- data.frame(P_seq, dens)
q_min <- ci_df$q_min[nseq]
q_max <- ci_df$q_max[nseq]

p_dist <- dist_df %>%
  ggplot(aes(x=P_seq, y=dens)) +
  geom_line() +
  geom_ribbon(data=subset(dist_df, q_min < P_seq & P_seq < q_max), aes(ymin = 0, ymax = dens), fill = "grey70", alpha=0.35)+
  geom_vline(xintercept=c(q_min, q_max), linetype = "dashed") +
  scale_y_continuous(breaks=seq(0, round(max(dens), digits=1),length.out=2)) +
  scale_x_continuous(breaks=seq(y_min, y_max, length.out=4),
                     labels = scales::number_format(accuracy = 0.01)) +
  coord_flip(xlim=c(y_min, y_max)) +
  xlab("") + ylab("Density") +
  theme_classic() +
  theme(legend.position="none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(r = 5,
                             l = -10))


### Panel C: Approach and objectives

SR_pc <- Species_all[which(Species_all$Mortality==mortality&Species_all$Fecundity==fecundity),]
#SR_pc <- SR_pc[SR_pc$Nb_obs!=n_axes,]
SR_pc <- SR_pc[SR_pc$Mod!="Perf_know",]

x1 <- -6
y1 <- y2 <- 102
y3 <- y4 <- -2
x2 <- x1 + 8

b2 <- y2 + x2
x3 <- b2 - y3
x4 <- x3 - 8


SR_pc <- SR_pc%>%
  dplyr::group_by(Mod, Nb_obs)%>%
  dplyr::mutate(mean_SR=mean(N_sp))%>%
  dplyr::slice(1)%>%
  dplyr::select(Mod, Nb_obs, mean_SR)%>%
  dplyr::arrange(factor(Mod, levels = c("Part_know_IV", "Part_know", "Perf_know")), Nb_obs)%>%
  dplyr::ungroup()
SR_pc$sIV <- c(((SR_pc[1:(nrow(SR_pc)-1),]$Nb_obs)/n_axes)*100, 100)
SR_pc$uIV <- c(0)
SR_pc[SR_pc$Mod=="Part_know",]$uIV <- rep(0, nrow(SR_pc[SR_pc$Mod=="Part_know",]))
SR_pc[SR_pc$Mod=="Part_know_IV",]$uIV <- ((n_axes - SR_pc[SR_pc$Mod=="Part_know_IV",]$Nb_obs)/n_axes)*100
#SR_pc[SR_pc$Mod=="Perf_know",]$uIV <- 0

pc<-data.frame(x = SR_pc$sIV,
               y = SR_pc$uIV,
               z = SR_pc$mean_SR)%>%
  ggplot(aes(x, y, z))+
  geom_point(size=2, aes(color = z), show.legend=T)+
  ggplot2::scale_colour_viridis_c('Species richness', guide=guide_colourbar(title.position = "top"))+
  xlab("sIV (%)")+
  ylab("uIV (%)")+
  scale_x_continuous(limits=c(x1, x3), breaks=seq(0,100,20))+
  scale_y_continuous(limits=c(y3, y1), breaks=seq(0,100,20))+
  geom_polygon(data = data.frame(x=c(10.5,16,16,10.5),y=c(89,89,-2,-2)), color="grey70", fill = NA) + # objective 2
  geom_polygon(data = data.frame(x=c(64,69.5,69.5,64),y=c(35,35,-2,-2)), color="grey70", fill = NA) +
  geom_polygon(data = data.frame(x=c(x1,x2,x3,x4),y=c(y1,y2,y3,y4)), color="dodgerblue4", fill = NA) + # objective 1
  theme_classic()+
  theme(
        aspect.ratio=1,
        legend.position = c(.7, .8),
        legend.direction="horizontal",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))

pc_second_x_axis <- ggplot2::ggplot(SR_pc, ggplot2::aes(Nb_obs))+
  ggplot2::geom_blank()+
  ggplot2::theme_classic()+
  ggplot2::scale_x_continuous(name="Number of observed dimensions",
                              expand = c(0, 0),
                              breaks=c(0, 3, 6, 9, 12, 15),
                              limits=c(min(SR_pc$Nb_obs), max(SR_pc$Nb_obs)))+
  ggplot2::theme(aspect.ratio=0.1,
                 panel.background = ggplot2::element_rect(fill='transparent'),
                 plot.background = ggplot2::element_rect(fill='transparent', color=NA),
                 text=ggplot2::element_text(size = 13),
                 axis.text=ggplot2::element_text(size=10),
                 axis.line.x=element_line(colour = "deeppink3"),
                 axis.ticks.x=ggplot2::element_line(colour="deeppink3"),
                 axis.text.x=ggplot2::element_text(colour="deeppink3"),
                 axis.title.x=ggplot2::element_text(colour="deeppink3"),
                 axis.line.y=ggplot2::element_blank(),
                 axis.text.y=ggplot2::element_blank(),
                 axis.ticks.y=ggplot2::element_blank(),
                 axis.title.y=ggplot2::element_blank(),
                 panel.grid.minor.y=ggplot2::element_blank(),
                 panel.grid.major.y=ggplot2::element_blank(),
                 plot.margin = ggplot2::margin(l = 36,
                                               r = 25,
                                               t = 0,
                                               b = 0))

pc_second_y_axis <- ggplot2::ggplot(SR_pc, ggplot2::aes(y=Nb_obs))+
  ggplot2::geom_blank()+
  ggplot2::theme_classic()+
  ggplot2::scale_y_continuous(name="Number of unobserved dimensions",
                              expand = c(0, 0),
                              breaks=c(0, 3, 6, 9, 12, 15),
                              limits=c(min(SR_pc$Nb_obs), max(SR_pc$Nb_obs)))+
  ggplot2::theme(aspect.ratio=1,
                 panel.background = ggplot2::element_rect(fill='transparent'),
                 plot.background = ggplot2::element_rect(fill='transparent', color=NA),
                 text=ggplot2::element_text(size = 13),
                 axis.text=ggplot2::element_text(size=10),
                 axis.line.y=element_line(colour = "deeppink3"),
                 axis.ticks.y=ggplot2::element_line(colour="deeppink3"),
                 axis.text.y=ggplot2::element_text(colour="deeppink3"),
                 axis.title.y=ggplot2::element_text(colour="deeppink3"),
                 axis.line.x=ggplot2::element_blank(),
                 axis.text.x=ggplot2::element_blank(),
                 axis.ticks.x=ggplot2::element_blank(),
                 axis.title.x=ggplot2::element_blank(),
                 panel.grid.minor.x=ggplot2::element_blank(),
                 panel.grid.major.x=ggplot2::element_blank(),
                 plot.margin = ggplot2::margin(l = 0,
                                               r = 0,
                                               t = 0,
                                               b = 0))

pb_dist <- ggpubr::ggarrange(pb, p_dist, nrow=1, ncol=2, widths=c(4,1))

pc_axes <-cowplot::ggdraw()+
  cowplot::draw_plot(pc, y=0.1, scale=0.8)+
  cowplot::draw_plot(pc_second_x_axis, x = 0.025, y = -0.35, scale=0.37)+
  cowplot::draw_plot(pc_second_y_axis, x = -0.11, y = 0.155, scale=0.57)

#fig <- ggpubr::ggarrange(ggpubr::ggarrange(pa, pb_dist, nrow=1, ncol=2, labels=c("A", "B")), pc, nrow=2, ncol=1, labels=c("", "C"))
fig <- ggpubr::ggarrange(ggpubr::ggarrange(pa, pb_dist, nrow=1, ncol=2, labels=c("A", "B")), pc_axes, nrow=2, ncol=1, labels=c("", "C"))

ggplot2::ggsave(fig, filename=here::here("outputs", "Fig_conceptual.png"),
                width=20, height=20, units="cm", dpi=300, bg="white")
