

clean_nebraska_data<-function() {
#load the data:
data=read.csv("data/Nebraska_Mosquito_Pool_Data_22_23_24.csv")

#extract relevant col;umns 
dat<-data.frame("surv_year"=as.factor(data$surv_year),"site_name"=data$site_name,"calculated_county"= data$calculated_county,"calculated_state"=data$calculated_state,
                "collection_date"=data$collection_date,"disease_week"=data$disease_week,"trap_type"=data$trap_type,"lure"=data$lure,"species"=data$species,
                "pool_size"=data$num_count,"WNV"=data$WNV, "ctval"=data$ctval)

#if there are ct values ==0, then make the WNV zero and ct value =40 
dat$ctval[dat$ctval==0]=40
dat$WNV[dat$ctval==40]=0
dat$WNV[dat$ctval<40]=1

dat<-dat  %>% filter(species!="Culex")
dat<-dat  %>% filter(species!="Culex erraticus")
dat <-dat %>% filter(species!="Coquillettidia perturbans")
dat <-dat %>% filter(species!="Culex territans")
dat <-dat %>% filter(species!="Culiseta inornata")

dat_neb<-dat[,c(1,6,9,10,11,12)]

dat_neb$state<-rep("Nebraska")

#if ct values are 0; make them ct=40 
#dat_neb$ctval[dat_neb$WNV==0]=40

dat_neb$species[dat_neb$species=="Culex pipiens/restuans/salinarius"]="Pipiens"
dat_neb$species[dat_neb$species=="Culex tarsalis"]="Tarsalis"

return(dat_neb)

}
