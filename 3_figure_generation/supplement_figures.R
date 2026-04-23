##plot the kdes 

pool_sizes_test=c(1,10, 20, 50,100,200)
f_kde_all=lapply(1:length(pool_sizes_test),function(x) readRDS(file = paste0("pre_calculations/kde_calcs/f_kde_",pool_sizes_test[x],".rds")))
ct_vec=1:40

couts<-expand.grid(pool_sizes_test,ct_vec)
fk_sets<-NULL


for(i in 1:nrow(couts)){
  inx=which(pool_sizes_test %in% couts[i,1])
  fkde=f_kde_all[[inx]]
  xx=unlist(lapply(fkde[1:length(fkde)], function(f) f(couts[i,2])))
  fk_sets=rbind(fk_sets,data.frame("ct_value"=couts[i,2],"n_postives"=0:(length(fkde)-1),"pool_size"=couts[i,1],"fk"=xx))
  
}

fk_sets$pool_size=as.factor(fk_sets$pool_size)
levels(fk_sets$pool_size)=paste0("Pool size = ",pool_sizes_test)  

ggplot(data=fk_sets,aes(x=ct_value,y=fk,color=(n_postives),group=as.factor(n_postives)))+
  geom_line()+
  facet_wrap(~as.factor(pool_size))+
  scale_color_viridis(discrete = F,option = "turbo") +
  default_theme+
  xlab("Ct value")+
  ylab("Density")+
  theme(legend.position = "bottom")+
  labs(color="Number of positives")

ggsave("figures/kde_all.png",last_plot(),height = 6,width = 8)  

yy=subset(fk_sets,fk_sets$pool_size=="Pool size = 50")
yy$n_postives=as.factor(yy$n_postives)
levels(yy$n_postives)=paste0(0:50, "\n  positives")  


ggplot(data=yy,aes(x=ct_value,y=fk),color="black")+
  geom_line()+
  facet_wrap(~as.factor(n_postives),scales="free_y",ncol=10)+
  # scale_color_viridis(discrete = TRUE,option = "turbo") +
  default_theme+
  xlab("Ct value")+
  ylab("Density")+
  theme(legend.position = "none")

ggsave("figures/kde_50_pool_size.png",last_plot(),height = 12,width = 16)  





 ###testing the accuracy and coverage of the prevalence estimates from ct based method vs binary only method 

estimated_prev_1=readRDS("pre_calculations/estimated_prevalence/estimated_prev_all_pools_all_1.rds")
estimated_prev_2=readRDS("pre_calculations/estimated_prevalence/estimated_prev_all_pools_all_2.rds")
estimated_prev_3=readRDS("pre_calculations/estimated_prevalence/estimated_prev_all_pools_all_3.rds")
estimated_prev_4=readRDS("pre_calculations/estimated_prevalence/estimated_prev_all_pools_all_4.rds")

estimated_prev=rbind(estimated_prev_1,estimated_prev_2,estimated_prev_3,estimated_prev_4)

prev_sub<-data.frame("n_pools_test"=estimated_prev$n_pools_test,"pool_sizes_test"=estimated_prev$pool_sizes_test,
                     "prev"=estimated_prev$prev,"prev_pooled"=estimated_prev$prev_pooled,"prev_cdc"=estimated_prev$prev_cdc )

m_prev_sub<-melt(prev_sub,id=c("n_pools_test","pool_sizes_test","prev"))




sot_estimated_prev=subset(prev_sub,prev_sub$n_pools_test!=1 )
sot_estimated_prev=subset(sot_estimated_prev,sot_estimated_prev$prev %in% c(0.001,0.01 ,0.1))

pre_dat=sot_estimated_prev %>%
  group_by(n_pools_test,pool_sizes_test,prev) %>%
  reframe(diff_pooled=prev-prev_pooled,diff_cdc=prev-prev_cdc)

dat<-pre_dat %>%
  group_by(n_pools_test,pool_sizes_test,prev) %>% 
  reframe(mn_pooled=mean(diff_pooled),pooled_q1=quantile(diff_pooled,probs=0.025,na.rm = T),
          pooled_q2=quantile(diff_pooled,probs=0.975,na.rm = T),
          mn_cdc=mean(diff_cdc,na.rm=T),cdc_q1=quantile(diff_cdc,probs=0.025,na.rm = T),
          cdc_q2=quantile(diff_cdc,probs=0.975,na.rm = T) )


lvs=paste0("Pool size = ",c(1,10,20,50,100,200))
dat$pool_sizes_test=as.factor(dat$pool_sizes_test)
levels(dat$pool_sizes_test)=lvs

lv2=paste0("Prevalence = ",c(0.001,0.01,0.1))
dat$prev=as.factor(dat$prev)
levels(dat$prev)=lv2

p1<-ggplot()+
  geom_line(data=subset(dat,prev=="Prevalence = 0.001"),aes(x=as.factor(n_pools_test),y=mn_pooled,group=as.factor(pool_sizes_test)),color="black")+
  geom_ribbon(data=subset(dat,prev=="Prevalence = 0.001"),aes(x=as.factor(n_pools_test),ymin=pooled_q1, 
                                                              ymax=pooled_q2,group=as.factor(pool_sizes_test)),fill="#ee8494", alpha = 0.4)+
  facet_wrap(~as.factor(pool_sizes_test),ncol=1)+default_theme+
  geom_hline(yintercept = 0,linetype="dotted",color="red")+
  ylab("Actual prevalence - estimated prevalence")+
  xlab("Number of pools")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p1


p2<-ggplot()+
  geom_line(data=subset(dat,prev=="Prevalence = 0.01"),aes(x=as.factor(n_pools_test),y=mn_pooled,group=as.factor(pool_sizes_test)),color="black")+
  geom_ribbon(data=subset(dat,prev=="Prevalence = 0.01"),aes(x=as.factor(n_pools_test),ymin=pooled_q1, 
                                                             ymax=pooled_q2,group=as.factor(pool_sizes_test)),fill="#008e7c", alpha = 0.4)+
  facet_wrap(~as.factor(pool_sizes_test),ncol=1)+default_theme+
  geom_hline(yintercept = 0,linetype="dotted",color="red")+
  ylab("Actual prevalence - estimated prevalence")+
  xlab("Number of pools")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p2

p3<-ggplot()+
  geom_line(data=subset(dat,prev=="Prevalence = 0.1"),aes(x=as.factor(n_pools_test),y=mn_pooled,group=as.factor(pool_sizes_test)),color="black")+
  geom_ribbon(data=subset(dat,prev=="Prevalence = 0.1"),aes(x=as.factor(n_pools_test),ymin=pooled_q1, 
                                                            ymax=pooled_q2,group=as.factor(pool_sizes_test)),fill="#816997", alpha = 0.4)+
  facet_wrap(~as.factor(pool_sizes_test),ncol=1)+default_theme+
  geom_hline(yintercept = 0,linetype="dotted",color="red")+
  ylab("Actual prevalence - estimated prevalence")+
  xlab("Number of pools")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
p3


p1|p2|p3

ggsave("figures/supp_fig_1.png",last_plot(),height = 10,width = 10)

#subset pooled ct values 
pooled_esti<-data.frame("n_pools_test"=estimated_prev$n_pools_test,"pool_sizes_test"=estimated_prev$pool_sizes_test,
                        "prev"=estimated_prev$prev,
                        "pooled_q1"=estimated_prev$pooled_q1,"pooled_q2"=estimated_prev$pooled_q2)

#check if the cis include the true prevlance 
estimated_prev$within_ct=between(estimated_prev$prev,estimated_prev$pooled_q1,estimated_prev$pooled_q2)*1
estimated_prev$within_cdc=between(estimated_prev$prev,estimated_prev$cdc_q1,estimated_prev$cdc_q2)*1

#get the counts 
sumrry_cts<-estimated_prev %>% 
  group_by(n_pools_test,pool_sizes_test,prev) %>%
  summarise(perce_pooled=(sum(within_ct)/n())*100,perce_cdc=(sum(within_cdc)/n())*100,.groups = 'drop')

m_summary_cts=melt(sumrry_cts,id=c("n_pools_test","pool_sizes_test","prev"))

m_summary_cts=subset(m_summary_cts,m_summary_cts$n_pools_test!=1)


p2<-ggplot(data=m_summary_cts,aes(x=as.factor(n_pools_test),
                                  y=value,color=as.factor(pool_sizes_test),group = as.factor(pool_sizes_test)))+
  geom_point(size=2,shape=17,alpha=1,stroke=2)+
  facet_wrap(variable~prev)+
  scale_color_viridis(discrete=T,option="D")+
  
  # scale_colour_manual(
  #  values = c("#007160", "#008e63", "#2fab63", "#78c664", "#b8e067","#f9f871")
  # )+
  default_theme+
  xlab("Number of pools")+
  ylab("Percentage (%)")+
  theme(legend.position = "bottom")+
  labs(color = "Pool size")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

p2

ggsave("figures/pre_coverage.png",last_plot())

p2<-ggplot()+
  geom_line(data=dat,aes(x=as.factor(n_pools_test),y=mn_cdc,group=as.factor(pool_sizes_test),color=as.factor(prev)))+
  geom_ribbon(data=dat,aes(x=as.factor(n_pools_test),ymin=cdc_q1, 
                           ymax=cdc_q2,group=as.factor(pool_sizes_test),fill=as.factor(pool_sizes_test)), alpha = 0.5)+
  facet_wrap(as.factor(pool_sizes_test)~as.factor(prev),ncol=3,scales="free_y")+
  default_theme+
  geom_hline(yintercept = 0,linetype="dotted")+
  theme(legend.position = "none")+
  ylab("Actual prevalence-estimated prevalence")+
  xlab("Number of pools")+
  theme(legend.position = "none")+
  ggtitle(" Estimation with binarised data")

p2


ggsave("figures/pre_test.png",last_plot(),height = 8,width = 14)


###visualisation of the data:
#read the clean real datasets 
source("functions/clean_nebraska_data.R")
dat_neb<-clean_nebraska_data()
source("functions/clean_colorado_data.R")
dat_col<-clean_colorado_data()
all_dat=rbind(dat_col,dat_neb)

#estimated prevalances for the them 
combined <- readRDS("combined.RDS")
combined <- combined %>%  
  rename(  surv_year = yer,    disease_week = week  )

combined$surv_year=as.factor(combined$surv_year)

nebbb_summary <- all_dat %>%
  group_by(surv_year, disease_week,state) %>%
  summarise(
    n_pools = n(),                                 # total number of pools tested
    mean_Ct = mean(ctval, na.rm = TRUE),              # mean Ct across all pools
    min_Ct = min(ctval, na.rm = TRUE),                # minimum Ct
    max_Ct = max(ctval, na.rm = TRUE),                # maximum Ct
    n_positive = sum(WNV == "1", na.rm = TRUE),  # count of positive pools
    mean_Ct_positive = mean(ctval[WNV == "1"], na.rm = TRUE), # mean Ct only positives
    avg_pool_size = mean(pool_size, na.rm = TRUE)  # average pool size
  ) %>%
  ungroup()

nebbb_summary$surv_year=as.factor(nebbb_summary$surv_year)

# combined_neb=subset(combined,combined$state=="Nebraska")
# Keep only rows that match in BOTH dataframe
joined_inner <- nebbb_summary %>%  inner_join(combined, by = c("surv_year", "disease_week","state"))

joined_inner <- joined_inner %>%
  mutate(
    cdc_combined = paste0(
      round(cdc_prev, 3), " (",
      round(q1_cdc, 3), ", ",
      round(q2_cdc, 3), ")"
    ), 
    ct_overall_combined = paste0(
      round(prev_esti, 3), " (",
      round(q1_esti, 3), ", ",
      round(q2_esti, 3), ")"
    ), 
   optim_p_combined = paste0(
      round(optim_p_dynamic, 3), " (",
      round(dynamic_q1, 3), ", ",
      round(dynamic_q2, 3), ")"
    )
  )

joined_inner_f=cbind(joined_inner[,c(1:10)],joined_inner[,c("cdc_combined","ct_overall_combined","optim_p_combined")])

neb_set<-subset(joined_inner_f,state=="Nebraska")
col_set<-subset(joined_inner_f,state=="Colorado")
library(xtable)

col_set$avg_pool_size=round(col_set$avg_pool_size)
xtable(neb_set[,-3])
xtable(col_set[,-3])



ci_summary <- joined_inner %>%
  mutate(
    binary_prev = round(cdc_prev, 5),
    binary_width = round(q2_cdc, 5) - round(q1_cdc, 5),
    
    ct_prev = round(prev_esti, 5),
    ct_width = round(q2_esti, 5) - round(q1_esti, 5),
    
    binary_rel_width = if_else(binary_prev > 0, binary_width / binary_prev),
    ct_rel_width = if_else(ct_prev > 0, ct_width / ct_prev),
    
    binary_width_for_zeros = if_else(binary_prev == 0, binary_width)
    ct_width_for_zeros = if_else(ct_prev == 0, ct_width)
  ) %>%
  group_by(state) %>%
  summarise(
    mean_binary_width_when_binary_zero = mean(binary_width_for_zeros, na.rm = TRUE),
    mean_ct_width_when_ct_zero = mean(ct_width_for_zeros, na.rm = TRUE),
   
    
    mean_binary_rel_width_when_binary_positive =
      mean(binary_rel_width, na.rm = TRUE),
    mean_ct_rel_width_when_ct_positive =
      mean(ct_rel_width, na.rm = TRUE),
    
    .groups = "drop"
  )



p_v2<-ggplot(all_dat,aes(x=disease_week,y=ctval,colour=state))+
  geom_point()+
  facet_wrap(~(as.factor(surv_year)),ncol=3)+
  default_theme+  
  xlab("Disease week")+
  ylab("Ct value")+
  scale_color_manual(values  =c("#903920","#006a55"))+
  theme(legend.position = "bottom")

p_v2


yers=2022:2024
weks<-c(1:23,39:52)
yr_wek<-expand.grid(yers,weks)

new_dat<-data.frame("surv_year"=yr_wek$Var1,"site_name"=rep(NA),"calculated_county"= rep(NA),"calculated_state"=rep(NA),
                    "collection_date"=rep(NA),"disease_week"=yr_wek$Var2,"trap_type"=rep(NA),"lure"=rep(NA),"species"=rep(NA),
                    "num_count"=rep(NA),"WNV"=rep(NA), "ctval"=rep(NA))
all_dat=rbind(dat,new_dat)



all_dat$yr_wek=paste(all_dat$surv_year,all_dat$disease_week)



p3<-ggplot(data = all_dat,aes(x=disease_week,y=ctval,colour = surv_year))+
  geom_point()+
  default_theme+
  xlab("Disease week")+
  scale_color_manual(values  = wes_palette("AsteroidCity1", 3, type = "continuous"))+
  ylab("Ct value")+
  facet_wrap(as.factor(species)~state)+
  ggtitle("Ct values by species")+
  theme(legend.position = "bottom")

p3

ggsave("figures/all_dat_species.png",p3,height = 10,width = 10)



p4<-ggplot(data = dat_neb,aes(x=disease_week,y=ctval,colour = surv_year))+
  geom_point()+
  default_theme+
  xlab("Disease week")+
  ylab("Ct value")+
  facet_wrap(~as.factor(calculated_county))+
  ggtitle("Ct values by the calculated county (Nebraska)")
p4


lv=c("Berthud","Larimer")
col$county=as.factor(col$county)
levels(col$county)=lv
p4<-ggplot(data = col,aes(x=week,y=cq,colour = as.factor(year)))+
  geom_point()+
  default_theme+
  xlab("Disease week")+
  ylab("Ct value")+
  facet_wrap(~as.factor(county))+
  ggtitle("Ct values by the calculated county (Colorado)")
p4

