#this script will calculaet the kdes for the given pool sizes 
#read the data saved from the abm, they are saved in teh synthetic data folder 
source("functions/kde_approx.R")
all_data<-readRDS("synthetic_data/all_data.rds")


#get the ct values that are static and and dynamic 




ct_threshold=40

#KDE estimation:
#positive_ctss<-all_data %>% filter(ct_value < ct_threshold) %>% pull(ct_value)

static_moz<-subset(all_data,ct_type_method==2)
dynamic_moz<-subset(all_data,ct_type_method==3)

#KDE estimation:
positive_ctss_static<-static_moz %>% filter(ct_value < ct_threshold) %>% pull(ct_value)
positive_ctss_dynamic<-dynamic_moz %>% filter(ct_value < ct_threshold) %>% pull(ct_value)

#transform to viral loads:
ct_to_vl<-function(ct,intercept,slope){ #ct=m*log10(vl)+c0
  vl<-10^((ct-intercept)/slope)
  vl   
}

vl_to_ct<-function(slope,intercept,vl){ #vl in 10^
  ct1<-slope*log10(vl) +intercept
  #ct<-min(ct1,40)
  ct1
}

intercept=36.9
slope= -2.7 

viral_loads_vec_static=ct_to_vl(positive_ctss_static,intercept,slope)
viral_loads_vec_dynamic=ct_to_vl(positive_ctss_dynamic,intercept,slope)

#pool_sizes=c(1:50,100, 200,500)

pool_sizes=c(1:50)

for (i in 1:length(pool_sizes)) {
  f_kde_i_static<-kde_approx(B=10000,x=pool_sizes[i],viral_loads=viral_loads_vec_static,intercept=36.9, slope= -2.7)
  f_kde_i_dynamic<-kde_approx(B=10000,x=pool_sizes[i],viral_loads=viral_loads_vec_dynamic,intercept=36.9, slope= -2.7)
  fil_s=paste0("kde_calcs/f_kde_static_",pool_sizes[i],".rds")
  fil_d=paste0("kde_calcs/f_kde_dynamic_",pool_sizes[i],".rds")
  saveRDS(f_kde_i_static, file = fil_s)
  saveRDS(f_kde_i_dynamic, file = fil_d)
}



#plot them:

pool_sizes_test=c(1,10, 20, 50,100,200)
f_kde_all=lapply(1:length(pool_sizes_test),function(x) readRDS(file = paste0("kde_calcs/f_kde_",pool_sizes_test[x],".rds")))
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



