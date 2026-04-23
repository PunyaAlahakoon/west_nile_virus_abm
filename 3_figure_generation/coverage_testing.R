
#load the default theme for figures:
theme<-source("functions/default_theme.R")
default_theme<-default_theme()

optim_prev_esti<-readRDS("optim_prev_estims_2.rds")

#write a script to check when the pool sizes are long  prevelance of the two are high:

#generate a known prevlanace:
#take pool size to be 5 
pool_sizes=50 
n_pools<-c(10,30,50,100)

#get the proprtion sample 


p_dynamic <- c(0.001, 0.01, 0.1)
prop      <- c(0.2, 0.4, 0.6, 0.8, 0.9) #proportion for dynamic

grids <- expand.grid(p_dynamic,prop, pool_sizes, n_pools)
grids$static<-(grids$Var1/grids$Var2)-grids$Var1
#grid <- expand.grid(p_static,px_grid[,1], pool_sizes, n_pools)


#grid <- expand.grid(p_static,px_grid[,1], pool_sizes, n_pools)

colnames(grids)<-c("p_dynamic","prop","pool_sizes","n_pools","p_static")
p_combi <- crossing(rep = 1:100, grids)


 n_prev<-nrow(p_combi)
#
 n_prev
#  
#  grid$true_dynamic_prop=round(grid$p_dynamic_true/(grid$p_dynamic_true+grid$p_static_true),2)
#  
#  optim_prev_esti$p_dynamic_true=p_combi$p_dynamic
#  optim_prev_esti$p_static_true=p_combi$p_static

 optim_prev_esti$true_dynamic_prop=optim_prev_esti$p_dynamic_true/(optim_prev_esti$p_dynamic_true+optim_prev_esti$p_static_true)

levs<-sort(round(unique( optim_prev_esti$true_dynamic_prop),3))
optim_prev_esti$true_dynamic_prop<-as.factor(optim_prev_esti$true_dynamic_prop)
levels(optim_prev_esti$true_dynamic_prop)<-levs



optim_prev_esti$propr_error_dynamic=optim_prev_esti$p_dynamic_true-optim_prev_esti$optim_p_dynamic

optim_prev_esti$pool_size=p_combi$pool_sizes
optim_prev_esti$n_pools<-p_combi$n_pools

lvs=paste0("Pool size = ",pool_sizes)
optim_prev_esti$pool_size=as.factor(optim_prev_esti$pool_size)
levels(optim_prev_esti$pool_size)=lvs

lv2=paste0("Productive infection \n prevalence = ",p_dynamic)
optim_prev_esti$true_Prevalence=as.factor(optim_prev_esti$p_dynamic_true)
levels(optim_prev_esti$true_Prevalence)=lv2



p1 <- ggplot(data = optim_prev_esti) +
  geom_boxplot(
    aes(
      x = as.factor(true_dynamic_prop),
      y = propr_error_dynamic,
      fill = as.factor(n_pools)
    ),position = "dodge"
  ) +
  geom_abline(
    slope = 0,
    intercept = 0,
     colour = "black",
    linetype = "dashed" 
  ) +
   default_theme+
  facet_wrap(~ as.factor(true_Prevalence),scales="free_y") +
  labs(
    x = "True productive infection proportion",
    y = " Actual- estimated prevalence",
    fill = "Number of pools"     
  )+
  ggtitle("Productive infection estimation error")+
  theme(legend.position = "bottom")

p1

#prevalence
#check if the cis include the true prevlance 
optim_prev_esti$within_ct_dynamic=between(optim_prev_esti$p_dynamic_true,optim_prev_esti$dynamic_q1,optim_prev_esti$dynamic_q2)*1

#get the counts 
sumrry_cts<-optim_prev_esti %>% 
  group_by(n_pools,p_dynamic_true,true_Prevalence,true_dynamic_prop) %>%
  summarise(perce_pooled=(sum(within_ct_dynamic)/n())*100,.groups = 'drop')

m_summary_cts=melt(sumrry_cts,id=c("n_pools","p_dynamic_true","true_Prevalence","true_dynamic_prop"))


#
p2<-ggplot(data=m_summary_cts,aes(x=as.factor(true_dynamic_prop),
                                  y=value,fill = as.factor(n_pools)))+
  geom_point(size=2,shape=24,alpha=1,stroke=1)+
  facet_wrap(~true_Prevalence,scales="free_y")+
  #scale_color_viridis(discrete=T,option="D")+
  
  # scale_colour_manual(
  #  values = c("#007160", "#008e63", "#2fab63", "#78c664", "#b8e067","#f9f871")
  # )+
  ylim(70,100)+
  labs(
    x = "True productive infection proportion",
    y = "Percentage (%)",
    fill = "Number of pools"     
  )+
  default_theme+
  theme(legend.position = "bottom")+
  ggtitle("Confidence interval coverage")

p2


ggarrange(p1,p2,ncol=1,common.legend = T, legend = "bottom")


ggsave("figures/coverage.png",last_plot(),height = 10,width = 14)




optim_prev_esti$covered<-(optim_prev_esti$dynamic_q1 <= optim_prev_esti$p_dynamic_true) & (optim_prev_esti$dynamic_q2 >= optim_prev_esti$p_dynamic_true)

optim_prev_esti <- optim_prev_esti %>%
  group_by(true_Prevalence, p_dynamic_true, true_dynamic_prop, n_pools) %>%
  arrange(p_dynamic_true, .by_group = TRUE) %>%
  mutate(
    rep_in_group = row_number()
  ) %>%
  ungroup()


coverage_summary <- optim_prev_esti %>%
  group_by(true_Prevalence, p_dynamic_true, true_dynamic_prop, n_pools) %>%
  summarise(
    coverage_rate = mean(covered, na.rm = TRUE),
    n_rep = n(),
    # median_ci_width = median(ci_width, na.rm = TRUE),
    .groups = "drop"
  )




# numeric pools
pools_to_consider <- c(10, 50,100)          # numeric
true_prop_to_consider <- c(0.2, 0.6)  # numeric

n_dd <- subset(optim_prev_esti,
               # true_Prevalence == "Prevalence = 0.001" &
                 n_pools %in% pools_to_consider &
                 true_dynamic_prop %in% true_prop_to_consider)

# make columns factors (operate on the data.frame itself)
n_dd$n_pools <- factor(n_dd$n_pools)
n_dd$true_dynamic_prop <- factor(n_dd$true_dynamic_prop)
n_dd<-droplevels(n_dd)


# Define expected raw levels (example)
pool_lvls <- sort(unique(n_dd$n_pools))
prop_lvls <- sort(unique(n_dd$true_dynamic_prop))

# Create labelled factors (keep order explicit if needed)
n_dd$n_pools <- factor(
  n_dd$n_pools,
  levels = pool_lvls,
  labels = paste0("Number of pools = ", pool_lvls)
)

n_dd$true_dynamic_prop <- factor(
  n_dd$true_dynamic_prop,
  levels = prop_lvls,
  labels = paste0("Proportion = ", prop_lvls)
)

cov_dat <- subset(coverage_summary,
                    n_pools %in% pools_to_consider &
                    true_dynamic_prop %in% true_prop_to_consider)

# factor and drop unused levels
cov_dat$n_pools <- factor(cov_dat$n_pools)
cov_dat$true_dynamic_prop <- factor(cov_dat$true_dynamic_prop)
cov_dat <- droplevels(cov_dat)   # drops any unused factor levels across the data.frame


# Define expected raw levels (example)
pool_lvls <- sort(unique(cov_dat$n_pools))
prop_lvls <- sort(unique(cov_dat$true_dynamic_prop))

# Create labelled factors (keep order explicit if needed)
cov_dat$n_pools <- factor(
  cov_dat$n_pools,
  levels = pool_lvls,
  labels = paste0("Number of pools = ", pool_lvls)
)

cov_dat$true_dynamic_prop <- factor(
  cov_dat$true_dynamic_prop,
  levels = prop_lvls,
  labels = paste0("Proportion = ", prop_lvls)
)

n_dd_1<-subset(n_dd, true_Prevalence == "Productive infection \n prevalence = 0.001")
n_dd_1<-droplevels(n_dd_1)

cov_dat_1<-subset(cov_dat,true_Prevalence == "Productive infection \n prevalence = 0.001")
cov_dat_1<-droplevels(cov_dat_1)

p_intervals_p1<- ggplot(n_dd_1, 
                              aes(x = rep_in_group, y = p_dynamic_true)) +
  geom_linerange(aes(ymin = dynamic_q1, ymax = dynamic_q2, color = covered),
                 size = 0.6, alpha = 0.9) +
  geom_point(aes(fill = covered), shape = 21, color = "black", size = 1.6) +
  # dashed line at the numeric truth for each facet (use data = coverage_summary to annotate)
  geom_hline(aes(yintercept = p_dynamic_true), data = cov_dat_1, color = "purple",
             linetype = "dashed", size = 0.5) +
  scale_color_manual(name = "Coverage", values = c("TRUE" = "#2c7bb6", "FALSE" = "#d73027")) +
  scale_fill_manual(name = "Coverage", values = c("TRUE" = "#2c7bb6", "FALSE" = "#d73027")) +
  labs(
    x = "Replication index",
    y = "Confidence interval "
  ) +
  coord_cartesian(ylim = c(0, 0.01)) +
  default_theme+
  facet_wrap(as.factor(n_pools)~true_dynamic_prop,ncol=2 )+
  # facet_grid(rows = vars(n_pools),
  #            cols = vars(true_Prevalence),
  #            scales = "free_y")
# annotate each facet with coverage rate and number of reps Confidence 
geom_label(
  data = cov_dat_1,
  aes(x = Inf, y = Inf,
      label = paste0("Coverage = ", percent(coverage_rate, accuracy = 0.1))),
  hjust = 1.05, vjust = 8, size = 3.0, inherit.aes = FALSE
) +
  ggtitle("Productive infection prevalence = 0.001")
# coord_cartesian(ylim = c(0, 0.7)) +

p_intervals_p1

n_dd_1<-subset(n_dd, true_Prevalence == "Productive infection \n prevalence = 0.01")
n_dd_1<-droplevels(n_dd_1)

cov_dat_1<-subset(cov_dat,true_Prevalence == "Productive infection \n prevalence = 0.01")
cov_dat_1<-droplevels(cov_dat_1)

p_intervals_p2<- ggplot(n_dd_1, 
                        aes(x = rep_in_group, y = p_dynamic_true)) +
  geom_linerange(aes(ymin = dynamic_q1, ymax = dynamic_q2, color = covered),
                 size = 0.6, alpha = 0.9) +
  geom_point(aes(fill = covered), shape = 21, color = "black", size = 1.6) +
  # dashed line at the numeric truth for each facet (use data = coverage_summary to annotate)
  geom_hline(aes(yintercept = p_dynamic_true), data = cov_dat_1, color = "purple",
             linetype = "dashed", size = 0.5) +
  scale_color_manual(name = "Coverage", values = c("TRUE" = "#2c7bb6", "FALSE" = "#d73027")) +
  scale_fill_manual(name = "Coverage", values = c("TRUE" = "#2c7bb6", "FALSE" = "#d73027")) +
  coord_cartesian(ylim = c(0, 0.05)) +
  labs(
    x = "Replication index",
    y = "Confidence interval"
  ) +
  default_theme+
  facet_wrap(as.factor(n_pools)~true_dynamic_prop,ncol=2 )+
  geom_label(
    data = cov_dat_1,
    aes(x = Inf, y = Inf,
        label = paste0("Coverage = ", percent(coverage_rate, accuracy = 0.1))),
    hjust = 1.05, vjust = 8, size = 3.0, inherit.aes = FALSE
  ) +
  ggtitle("Productive infection prevalence = 0.01")
p_intervals_p2



n_dd_1<-subset(n_dd, true_Prevalence == "Productive infection \n prevalence = 0.1")
n_dd_1<-droplevels(n_dd_1)

cov_dat_1<-subset(cov_dat,true_Prevalence == "Productive infection \n prevalence = 0.1")
cov_dat_1<-droplevels(cov_dat_1)

p_intervals_p3<- ggplot(n_dd_1, 
                        aes(x = rep_in_group, y = p_dynamic_true)) +
  geom_linerange(aes(ymin = dynamic_q1, ymax = dynamic_q2, color = covered),
                 size = 0.6, alpha = 0.9) +
  geom_point(aes(fill = covered), shape = 21, color = "black", size = 1.6) +
  # dashed line at the numeric truth for each facet (use data = coverage_summary to annotate)
  geom_hline(aes(yintercept = p_dynamic_true), data = cov_dat_1, color = "purple",
             linetype = "dashed", size = 0.5) +
  scale_color_manual(name = "Coverage", values = c("TRUE" = "#2c7bb6", "FALSE" = "#d73027")) +
  scale_fill_manual(name = "Coverage", values = c("TRUE" = "#2c7bb6", "FALSE" = "#d73027")) +
  coord_cartesian(ylim = c(0, 0.45)) +
  labs(
    x = "Replication index",
    y = "Confidence interval"
  ) +
  default_theme+
  facet_wrap(as.factor(n_pools)~true_dynamic_prop,ncol=2  )+
  geom_label(
    data = cov_dat_1,
    aes(x = Inf, y = Inf,
        label = paste0("Coverage = ", percent(coverage_rate, accuracy = 0.1))),
    hjust = 1.05, vjust = 8, size = 3.0, inherit.aes = FALSE
  ) +
  ggtitle("Productive infection prevalence = 0.1")
p_intervals_p3

ggarrange(p_intervals_p1,p_intervals_p2,p_intervals_p3,ncol=3,common.legend = T,legend = "bottom")

ggsave("figures/coverage_perc.png",last_plot(),height = 12,width = 18)
