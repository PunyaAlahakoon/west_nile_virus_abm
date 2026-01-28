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


# Define range of pool sizes
pool_sizes <- 1:50


# Maximum pool size
max_pool <- 50

# Initialize list to store combinations for each pool size
all_combos <- vector("list", max_pool)


joint_kde<-function(B,obs_k_mat,viral_loads_vec_static,
                    viral_loads_vec_dynamic ){
  ct_threshold<-40 
  
  intercept=36.9
  slope= -2.7
  #  den_estimates <- NULL
  nn<-nrow(obs_k_mat)
  all_kdes <- list(length=nn)
  
  
  
  # den_estimates<-data.frame(matrix(NA,ncol=3,nrow=0))
  # colnames(den_estimates)<-c("k","x","y")
  for(i in 1:nn){
    si<-obs_k_mat[i,1]
    di<-obs_k_mat[i,2]
    x<-sum(obs_k_mat[i,]) #total pools size 
    
    #this will always be the first value
    if((si+di)==0){
      all_kdes[[i]] <- function(x,r=0.00001){
        ct_threshold=40
        y <- numeric(length(x))
        y[x < ct_threshold] <- r
        y[x >= ct_threshold] <- 1-r
        y
      }
    } else {
      
      z<-replicate(B,sum(c(sample(viral_loads_vec_static,si),
                           sample(viral_loads_vec_dynamic,di)) )) #replicate the function B times 
      
      z=z/x
      ct_values<-slope*log10(z)+intercept
      ct_values[ct_values>ct_threshold]=ct_threshold
      #kernsmpooth package:
      
      kernel_i<-adaptiveKernel(ct_values,to = ct_threshold,from=0)
      all_kdes[[i]] <- approxfun(kernel_i$x, kernel_i$y,rule=2)
      #den_estimates<-rbind(den_estimates,data.frame("k"= rep(i),"x"=kernel_i$x,"y"=kernel_i$y))
    }
    
  }
  return(all_kdes)
}



for (n in 1:max_pool) {
  # For each pool size n, generate all combinations of (static, dynamic, negatives)
  combos_n <- expand.grid(
    static = 0:n,
    dynamic = 0:n,
    negatives = 0:n
  )
  
  # Keep only rows where sum == n
  combos_n <- subset(combos_n, static + dynamic + negatives == n)
  
  # Store in list
  all_combos[[n]] <- combos_n
}

# Example: combinations for pool size 5
head(all_combos[[5]])




for (i in 1:length(pool_sizes)) {
  
  obs_k_mat<-all_combos[[i]]
  joint_fkde<-joint_kde(B=10000,obs_k_mat,viral_loads_vec_static,
  viral_loads_vec_dynamic)
  fil_s=paste0("joint_kde_calc/joint_fkde_",pool_sizes[i],".rds")

  saveRDS(joint_fkde, file = fil_s)

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



