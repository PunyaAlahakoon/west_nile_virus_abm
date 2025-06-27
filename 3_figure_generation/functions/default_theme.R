

default_theme<-function(){

theme=theme_minimal(15)+
  theme(plot.background=element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        #axis.title.y=element_blank(),
        #strip.text.x = element_text(size = 14),
        # strip.text.y = element_text(size = 14)
  )+
  theme(panel.grid.major.y = element_line(color = "grey",
                                          linewidth = 0.5,
                                          linetype = "dotted"))+
  theme(legend.title=element_blank(),legend.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

return(theme)

}