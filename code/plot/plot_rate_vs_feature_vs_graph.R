source("plot_functions.R")

library(ggplot2)
library(dplyr)
library(HDInterval)
library(cowplot)
library(pdftools)
library(latex2exp)
library(colorspace)
library(ggrepel)
library(Cairo)

# files
fp = "/Users/mlandis/projects/fig_model/"
fn = "out"
out_fp    = paste0(fp, "output/anolis/")
code_fp   = paste0(fp, "code/")
dat_fp    = paste0(fp, "data/anolis/")

geo_fn    = paste0(out_fp, fn, ".geo_features.log")
mdl_fn    = paste0(out_fp, fn, ".model.log")
dist_fn   = paste0(dat_fp, "/anolis_nr9.distance_km.txt" )
size_fn   = paste0(dat_fp, "/anolis_nr9.size_km2.txt" )
tip_fn    = paste0(dat_fp, "anolis_nr9_ns383.range.tsv")
lbl_fn    = paste0(dat_fp, "state_labels.csv")
plot_fn   = paste0(code_fp, "plot/fig_rate_vs_feature.pdf")

# read in datasets
dat_geo    = read.table(geo_fn, header=T)
dat_mdl    = read.table(mdl_fn, header=T)
dat_dist   = read.table(dist_fn, header=T, sep=",")
dat_size   = read.table(size_fn, header=T, sep=",")
dat_states = read.table(lbl_fn, header=T, sep=",", colClasses = c("character","character"))

# burn-in/thinning (speed-up)
f_burn   = 0.5 # 
thin_by = 5
keep_idx = max(1,floor(nrow(dat_geo)*f_burn)):nrow(dat_geo)
keep_idx = keep_idx[ seq(1,length(keep_idx),by=thin_by) ]
dat_geo  = dat_geo[ keep_idx, ]
dat_mdl  = dat_mdl[ keep_idx, ]

# # unit distance from 1 to 2 (centroid)
# unit_dist = 13.904
# # km distance from 1 to 2 (centroid)
# km_dist = 1475.423
# # unit size region 8
# unit_size = 6.953
# # km^s size of region 8
# km_size = 76192


# convert to absolute units
dat_dist = dat_dist #* (km_dist/unit_dist) #106.115 #/ dat_dist[1,2] * 2300
dat_size = t(dat_size)[,1]
dat_size = dat_size #* (km_size/unit_size) #/ dat_size[8] * 76192  # Hispaniola (region H/8) is 76192 km^2

# thresholds
sig_mult       = 5
log10_sig_mult = log(sig_mult,base=10)

# regions and colors
num_regions  = 9
reg_lbl      = LETTERS[1:num_regions]
reg_long_lbl = c('NA & Upper CA',
                 'Lower CA',
                 'Chocó & Caribe',
                 'Andes',
                 'Amazonia',
                 'LA & PR',
                 'Cuba & Jamaica',
                 'Hispaniola',
                 'Bahamas')

colors       = c(Continental="#9E248F",Insular="#00AEEF")
colors       = c("#9E248F","#00AEEF")
reg_types    = c("Continental","Insular")[ c(1,1,1,1,1,2,2,2,2) ] # 1: continental, 2: insular

# get matrix of event region codes
lbl = make_region_code_matrix( reg_lbl )


#--------------#
# Process data #
#--------------#

# simple plotting data
dat12a = make_plot_dat( dat_mdl, dat_geo, num_regions )
#dat12b = make_plot_dat_old( dat_mdl, dat_geo, num_regions )
dat12 = make_plot_dat2(dat12a)
dat12_bd = dat12[ dat12$from_region != dat12$to_region, ]
dat12_we = dat12[ dat12$from_region == dat12$to_region, ]


# linear model plotting data
lm_thin_by = 10
df_lm = make_lm_plot_dat( dat_dist, dat_size, dat12a )
summ_lm = summarize_lm_dat(df_lm)
thresh_bd_cont = summ_lm$log10_upper[ summ_lm$type=="continental" & summ_lm$ratio=="bd" ]
thresh_bd_ins = summ_lm$log10_upper[ summ_lm$type=="insular" & summ_lm$ratio=="bd" ]
thresh_we_cont = summ_lm$log10_upper[ summ_lm$type=="continental" & summ_lm$ratio=="we" ]
thresh_we_ins = summ_lm$log10_upper[ summ_lm$type=="insular" & summ_lm$ratio=="we" ]

    
# range data table for tip taxa
dat_tip = read.table(tip_fn, header=F)
colnames(dat_tip) = c("name","state")
dat_tip$range = dat_states[match(dat_tip$state, dat_states$state),2]
dat_count = make_dat_count(dat_tip, dat_size, region_names=reg_lbl)
dat_count$type = c("Continental","Insular")[ c(1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1) ]
dat_count_region = make_dat_count_region( dat_count )
dat_count$n[ dat_count$lbl == "AB" ] = dat_count$n[ dat_count$lbl == "AB" ] - 1 # rm Bplumifrons
dat_count$n[ dat_count$lbl == "E" ] = dat_count$n[ dat_count$lbl == "E" ] - 1   # rm Pmarmoratus
dat_count$n[ dat_count$lbl == "CDE" ] = dat_count$n[ dat_count$lbl == "CDE" ] = - 1 # rm Pscapulatus
dat_count$n[ dat_count$lbl == "CDE" ] = dat_count$n[ dat_count$lbl == "CDE" ] = - 1 # rm Ugallardoi
dat_count = dat_count[ -dat_count$n<=0, ]

#---------------------#
# Plotting aesthetics #
#---------------------#

# theme
p_theme = theme(legend.position='none',
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.title       = element_text(size=10),
              axis.line        = element_line(colour = "black"),
              axis.text        = element_text(size=8),
              plot.title       = element_text(hjust = 0.5),
              strip.background = element_rect(color="black",fill="white",size=1.5,linetype="solid"))

p_theme2 = theme(legend.position='none',
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.title = element_text(size=10),
              axis.line = element_line(colour = "black"),
              axis.text.y= element_text(size=8),
              axis.text.x = element_text(size=8, angle=90, hjust=1, vjust=0.5),
              plot.title = element_text(hjust = 0.5),
              strip.background=element_rect(color="black",fill="white",size=1.5,linetype="solid"))



# rate plots
point_size_color       = 1.25
point_size_white       = point_size_color*0.25
text_size              = 2.25

# histogram plot
lwd_histogram          = 0.5
text_size_histogram    = 2.8
text_y_histogram       = 73
line_dx_histogram      = 0.3


# B vs. D
my_aes1_bd = aes(x=mean_dij, y=mean_bij, color=insular, label=lbl )
my_aes1_b  = aes(x=mean_dij, xend=mean_dij, y=hdil_bij, yend=hdiu_bij, color=insular)
my_aes1_d  = aes(y=mean_bij, yend=mean_bij, x=hdil_dij, xend=hdiu_dij, color=insular)

# W vs. E
my_aes2_we = aes(x=mean_ei, y=mean_wi, color=insular, label=lbl )
my_aes2_w  = aes(x=mean_ei, xend=mean_ei, y=hdil_wi, yend=hdiu_wi, color=insular)
my_aes2_e  = aes(y=mean_wi, yend=mean_wi, x=hdil_ei, xend=hdiu_ei, color=insular)

# B/D vs. dist
my_aes3_bd  = aes(x=log10( dist_ij), y=log10(mean_bdij), color=insular,  label=lbl)
my_aes3_hpd = aes(y=log10(hdil_bdij), yend=log10(hdiu_bdij), x=log10(dist_ij), xend=log10(dist_ij), color=insular)

# W/E vs. size
my_aes4_we  = aes(x=log10(size_i), y=log10(mean_wei), color=insular,  label=lbl)
my_aes4_hpd = aes(y=log10(hdil_wei), yend=log10(hdiu_wei), x=log10(size_i), xend=log10(size_i), color=insular)

my_aes_from    = aes(x=mean_dij, y=mean_bij, color=from_insular, label=from_lbl )
my_aes_to      = aes(x=mean_dij, y=mean_bij, color=to_insular, label=to_lbl )

# DB 2
my_aes_bd_seg  = aes(y=hdil_bdij, yend=hdiu_bdij, x=dist_ij, xend=dist_ij, color=insular)

my_aes_from_bd = aes(x=dist_ij, y=mean_bdij, color=from_insular, label=from_lbl )
my_aes_to_bd   = aes(x=dist_ij, y=mean_bdij, color=to_insular, label=to_lbl )

# EW 3,4

my_aes_we2 =  aes(x=size_i, y=mean_wei, color=insular, label=lbl )
my_aes_we_seg = aes(x=size_i, xend=size_i, y=hdil_wei, yend=hdiu_wei, color=insular)

# dist vs size
my_aes_from_sd = aes(x=size_i, y=dist_ij, color=from_insular, label=from_lbl )
my_aes_to_sd   = aes(x=size_i, y=dist_ij, color=to_insular, label=to_lbl )


##############
# MAKE PLOTS #
##############



#-----------------------------------#
# p1: plot r_d vs. r_b
#-----------------------------------#


p = ggplot(dat12_bd, my_aes1_bd)
p = p + xlab( TeX("$r_d$") )
p = p + ylab( TeX("$r_b$" ))
p = p + scale_x_log10(limits=c(5e-5,5e-2))
p = p + scale_y_log10(limits=c(1e-5,2e1))
p = p + geom_abline(intercept=0, slope=1, color="black", lty=3)
p = p + geom_abline(intercept=c(-log10_sig_mult,log10_sig_mult), slope=1, color="gray", lty=3)
p = p + geom_segment( data=dat12_bd, my_aes1_b, alpha=0.5)
p = p + geom_segment( data=dat12_bd, my_aes1_d, alpha=0.5)
p = p + scale_color_manual(values=colors)
p = p + geom_point(size=point_size_color)
p = p + geom_point(size=point_size_white, color="white")
p = p + p_theme
p1 = p

#-----------------------------------#
# p2: plot r_e vs. r_w
#-----------------------------------#

#min_hdil_ei = 1e-5
#dat12_we$hdil_ei[ dat12_we$hdil_ei<min_hdil_ei ] = min_hdil_ei    

p = ggplot(dat12_we, my_aes2_we)
p = p + xlab(TeX("$r_e$"))
p = p + ylab(TeX("$r_w$"))
p = p + scale_x_log10(limits=c(5E-7,1E-1))
#p = p + scale_x_log10(limits=c(min_hdil_ei,1E-1))
p = p + scale_y_log10(limits=c(1E-2,1E-1))
p = p + geom_abline(intercept=0, slope=1, color="black", lty=3)
p = p + geom_abline(intercept=c(-log10_sig_mult,log10_sig_mult), slope=1, color="gray", lty=3)
p = p + geom_segment( data=dat12_we, my_aes2_w, alpha=0.5)
p = p + geom_segment( data=dat12_we, my_aes2_e, alpha=0.5)
p = p + scale_color_manual(values=colors)
p = p + geom_point(size=point_size_color)
p = p + geom_point(size=point_size_white, color="white")
p = p + p_theme
p2 = p
    
#--------------------------------------#
# plm: plot r_b/r_d-ratio vs. distance #
#--------------------------------------#

df_lm_thin_bd = df_lm[ df_lm$ratio=="bd", ]
df_lm_thin_bd = df_lm_thin_bd[ seq(1, nrow(df_lm_thin_bd), by=lm_thin_by), ]

p = ggplot(dat12_bd, my_aes3_bd)
p = p + scale_y_continuous( limits=c(-5,5), breaks=seq(-4,4,by=2), labels=c("0.0001", "0.01", "1", "100", "10000"))
p = p + scale_x_continuous( limits=c(1.4,4), breaks=log10(c(100,300,1000,3000)), labels=c("100","300","1000","3000"))
p = p + geom_hline( yintercept=0, lty=3 )
p = p + geom_hline( yintercept=log10(c(5,1/5)), lty=3, col="gray" )
p = p + geom_vline( xintercept=thresh_bd_cont, color=colors[1], lty=2)
p = p + geom_vline( xintercept=thresh_bd_ins, color=colors[2], lty=2)
for (i in 1:nrow(df_lm_thin_bd)) {
    p = p + geom_abline(intercept=df_lm_thin_bd$intercept_p[i], slope=df_lm_thin_bd$slope_p[i], alpha=0.05, color=colors[df_lm_thin_bd$insular[i]+1] )
}

below_mrcd = dat12_bd[ dat12_bd$dist_ij < 10^thresh_bd_cont & !dat12_bd$insular, ]
#p = p + annotate("text", x=log10(below_mrcd$dist_ij)+0.05, y=log10(below_mrcd$mean_bdij)-0.3, label=below_mrcd$lbl, color=colors[1])


p = p + scale_color_manual( values=colors )
p = p + geom_segment( data=dat12_bd, mapping=my_aes3_hpd, alpha=0.5)
p = p + geom_point(size=point_size_color)
p = p + geom_point(size=point_size_white, color="white")

p = p + geom_text_repel(data=below_mrcd,
                        mapping=aes(x=log10(dist_ij),
                        y=log10(mean_bdij)),
                        color="black",
                        size=2,
                        min.segment.length=0, box.padding=0.3, point.padding=0.3,
                        segment.shape=-1, segment.size=0.15, nudge_x=+0.37, nudge_y=-2.6)

p = p + xlab( "Distance (km)" )
p = p + ylab( TeX("$r_b / r_d$") )
p = p + p_theme
p3 = p

#--------------------------------------#
# plm: plot r_w/r_e-ratio vs. size     #
#--------------------------------------#

df_lm_thin_we = df_lm[ df_lm$ratio=="we", ]
df_lm_thin_we = df_lm_thin_we[ seq(1, nrow(df_lm_thin_we), by=lm_thin_by), ]

p = ggplot(dat12_we, my_aes4_we)
p = p + scale_y_continuous( limits=c(-2,5), breaks=seq(-4,4,by=2), labels=c("0.0001", "0.01", "1", "100", "10000"))
p = p + scale_x_continuous( limits=c(4,7), breaks=log10(c(1e4,1e5,1e6,1e7)), labels=10^seq(4,7))
for (i in 1:nrow(df_lm_thin_we)) {
    p = p + geom_abline(intercept=df_lm_thin_we$intercept_p[i], slope=df_lm_thin_we$slope_p[i], alpha=0.05, color=colors[ df_lm_thin_we$insular[i]+1] )
}
p = p + scale_color_manual( values=colors )
p = p + geom_hline( yintercept=0, lty=3 )
p = p + geom_hline( yintercept=0, lty=3 )
p = p + geom_hline( yintercept=log10(c(5,1/5)), lty=3, col="gray" )
p = p + geom_segment( data=dat12_we, mapping=my_aes4_hpd, alpha=0.5)
p = p + geom_point(size=point_size_color)
p = p + geom_point(size=point_size_white, color="white")
p = p + xlab( TeX("Size (km$^2$)") )
p = p + ylab( TeX("$r_w / r_e$") )
p = p + p_theme
p4 = p


####################
# Plot frequencies #
####################

#-------------------------#
# plot range histogram    #
#-------------------------#

p = ggplot(dat_count, aes(x=lbl, y=n, fill=type))
p = p + geom_col()
p = p + scale_fill_manual(values=colors)
p = p + ylim(c(0,82))
p = p + xlab("Range")
p = p + ylab("Species per range")
p = p + annotate(geom="text", label="186 spp.\n(49%)", x=3,   y=text_y_histogram+7, color=colors[1], size=text_size_histogram)
p = p + annotate(geom="text", label="164 spp.\n(43%)", x=7.5, y=text_y_histogram+7, color=colors[2], size=text_size_histogram)
p = p + annotate(geom="text", label="29 spp.\n(8%)",   x=13,  y=text_y_histogram+7, color=colors[1], size=text_size_histogram)

p = p + geom_segment(x= 1-line_dx_histogram, xend= 5+line_dx_histogram, y=text_y_histogram, yend=text_y_histogram,   lwd=lwd_histogram, color=colors[1])
p = p + geom_segment(x= 1-line_dx_histogram, xend= 1-line_dx_histogram, y=text_y_histogram, yend=text_y_histogram-2, lwd=lwd_histogram, color=colors[1])
p = p + geom_segment(x= 5+line_dx_histogram, xend= 5+line_dx_histogram, y=text_y_histogram, yend=text_y_histogram-2, lwd=lwd_histogram, color=colors[1])

p = p + geom_segment(x= 6-line_dx_histogram, xend= 9+line_dx_histogram, y=text_y_histogram, yend=text_y_histogram,   lwd=lwd_histogram, color=colors[2])
p = p + geom_segment(x= 6-line_dx_histogram, xend= 6-line_dx_histogram, y=text_y_histogram, yend=text_y_histogram-2, lwd=lwd_histogram, color=colors[2])
p = p + geom_segment(x= 9+line_dx_histogram, xend= 9+line_dx_histogram, y=text_y_histogram, yend=text_y_histogram-2, lwd=lwd_histogram, color=colors[2])

p = p + geom_segment(x=10-line_dx_histogram, xend=16+line_dx_histogram, y=text_y_histogram, yend=text_y_histogram,   lwd=lwd_histogram, color=colors[1])
p = p + geom_segment(x=10-line_dx_histogram, xend=10-line_dx_histogram, y=text_y_histogram, yend=text_y_histogram-2, lwd=lwd_histogram, color=colors[1])
p = p + geom_segment(x=16+line_dx_histogram, xend=16+line_dx_histogram, y=text_y_histogram, yend=text_y_histogram-2, lwd=lwd_histogram, color=colors[1])

p = p + p_theme2
pc_region = p


#----------------------------------------#
# plot number of spp per region vs. size #
#----------------------------------------#

p = ggplot(dat_count_region, aes(x=mean_size, y=n, color=type, label=lbl, alpha=single))
p = p + geom_text(hjust=0)
p = p + scale_x_log10( limits=c( min(dat_size)*(1/1.3), max(dat_size)*1.3) )
p = p + ylim(c(0, 82))
p = p + scale_color_manual(values=colors)
p = p + scale_alpha(range=c(0.4,1.0))
p = p + xlab( TeX("Size (km$^2$)") )
p = p + ylab("Species per region")
p = p + p_theme
pc_size = p


#-------------------------------------------#
# pc_sd plot distance vs. widespread ranges #
#-------------------------------------------#

dat12ws = dat12[ dat12$from_region!=dat12$to_region, ]
dat12ws$is_widespread = 0
dat12ws$is_widespread[ dat12ws$lbl %in% c("BA","BC","CD","CE","DC","DE") ] = 1

p = ggplot(dat12ws, aes(x=size_i, y=dist_ij, color=insular, label=lbl, alpha=is_widespread))
p = p + scale_x_log10( limits=c( min(dat_size)*(1/1.3), max(dat_size)*1.3) )
p = p + scale_y_log10()
p = p + scale_color_manual(values=colors)
p = p + scale_alpha(range=c(0.3,1.0))
p = p + xlab(TeX( "Size (km$^2$)") )
p = p + ylab("Distance (km)")
p = p + geom_hline( yintercept=10^thresh_bd_cont, lty=2, color=colors[1])
p = p + geom_hline( yintercept=10^thresh_bd_ins,  lty=2, color=colors[2])
p = p + annotate( x=1e7, y=10^thresh_bd_cont, hjust=1, vjust=-0.2, geom="text", label=TeX("$r_d \\geq r_b$"), color=colors[1], size=3.5)
p = p + annotate( x=1e7, y=10^thresh_bd_ins,  hjust=1, vjust=-0.2, geom="text", label=TeX("$r_d \\geq r_b$"), color=colors[2], size=3.5)

for (i in 1:nrow(dat12_bd)) {
    p = p + geom_text(data=dat12ws[i,], my_aes_from_sd, size=text_size*1.5, nudge_x=-0.06)
    p = p + geom_text(data=dat12ws[i,], my_aes_to_sd,   size=text_size*1.5, nudge_x=+0.06)
}
p = p + p_theme
pc_sd = p



#-----------------------------------#
# pl: plot legend
#-----------------------------------#
p = ggplot()
p = p + annotate( "text",
                  label=reg_long_lbl,
                  x=rep(0,9), y=(9:1)*0.3,
                  hjust=0, size=3,
                  color=colors[c(rep(1,5),rep(2,4))])
p = p + xlim(0,2)
p = p + theme(legend.position='none',
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.title       = element_blank(),
              axis.line        = element_blank(),
              axis.text        = element_blank(),
              axis.ticks       = element_blank(),
              plot.title       = element_blank(),
              strip.background = element_rect(color="black",fill="white",size=1.5,linetype="solid"))

pl = p


#-------------------------------#
# pm: plot map
#-------------------------------#
map_font_size = 4
map_fn = paste0( code_fp, "plot/anolis_nr9_geo.pdf" )
my_map = magick::image_read_pdf(map_fn, density = 600)
p = ggdraw() + draw_image(my_map, clip="on", scale=1.45, height=1, width=1)
coords = matrix( c( 0.22, 0.68,
                    0.36, 0.55,
                    0.47, 0.69,
                    0.42, 0.25,
                    0.67, 0.59,
                    0.62, 0.78,
                    0.38, 0.90,
                    0.54, 0.83,
                    0.48, 0.93), byrow=T, ncol=2)
p = p + geom_point(x=coords[,1], y=coords[,2], col=colors[c(1,1,1,1,1,2,2,2,2)], size=10)
p = p + annotate("text", label=LETTERS[1:9], x=coords[,1], y=coords[,2], size=map_font_size, col=colors[c(1,1,1,1,1,2,2,2,2)])
p = p + annotate( "text",
                  label=paste0(LETTERS[1:9], ":"),
                  x=rep(0.0,9), y=(9:1)*0.05,
                  hjust=1, size=3,
                  color=colors[c(rep(1,5),rep(2,4))])
p = p + annotate( "text",
                  label=reg_long_lbl,
                  x=rep(0.03,9), y=(9:1)*0.05,
                  hjust=0, size=3,
                  color=colors[c(rep(1,5),rep(2,4))])
p = p + theme(plot.margin = unit(c(-2, -2, -2, -2), "cm"))
pm = p


#---------------------------#
# pgg plot graphs 
#---------------------------#



xy = c(0.1, 0.3)
dxy = c(0.02, 0.05)
font_size=3
ana_fn = paste0( code_fp, "plot/fig_geo_DE.pdf" )
ana_graph = magick::image_read_pdf(ana_fn, density = 600)
p = ggdraw() + draw_image(ana_graph, clip="on", scale=1.05, height=1, width=1, x=0.0, y=-0.0)

p = p + annotate("text", label=c("Rates", TeX("Extinction, $r_e$"), TeX("Dispersal, $r_d$")),
                  x = c( xy[1], xy[1]+dxy[1], xy[1]+dxy[1]),
                  y=  c( xy[2], xy[2]-dxy[2], xy[2]-(2*dxy[2])), hjust=c(0.5,0,0), size=font_size)
p = p + annotate("point",   x=xy[1], y=xy[2]-dxy[2], colour="black", size=2)
p = p + annotate("point",   x=xy[1], y=xy[2]-dxy[2], colour="white", size=1)
p = p + annotate("segment", x=xy[1]-0.01, xend=xy[1]+0.01, y=xy[2]-(2*dxy[2]), yend=xy[2]-(2*dxy[2]), size=1.5)
p = p + theme(plot.margin = unit(rep(0,4), "cm"))
p_ana = p



cla_fn = paste0( code_fp, "plot/fig_geo_BW.pdf" )
cla_graph = magick::image_read_pdf(cla_fn, density = 600)
p = ggdraw() + draw_image(cla_graph, clip="on", scale=1.05, height=1, width=1, x=0.0, y=-0.0)
p = p + annotate("text", label=c("Rates", "Within-region", TeX("speciation, $r_w$"), "Between-region", TeX("speciation, $r_b$")),
                  x = c( xy[1], xy[1]+dxy[1], xy[1]+dxy[1], xy[1]+dxy[1], xy[1]+dxy[1]),
                  y=  c( xy[2], xy[2]-dxy[2], xy[2]-2*dxy[2], xy[2]-3*dxy[2], xy[2]-4*dxy[2]), hjust=c(0.5,0,0,0,0), size=font_size)
p = p + annotate("point",   x=xy[1], y=xy[2]-dxy[2], colour="black", size=2)
p = p + annotate("point",   x=xy[1], y=xy[2]-dxy[2], colour="white", size=1)
p = p + annotate("segment", x=xy[1]-0.01, xend=xy[1]+0.01, y=xy[2]-(3*dxy[2]), yend=xy[2]-(3*dxy[2]), size=1.5)
p = p + theme(plot.margin = unit(rep(0,4), "cm"))
p_cla = p


#-----------------------------------#
# pg: plot  them all!
#-----------------------------------#


pg3u = plot_grid(  pm, pc_region, align="hv", nrow=1, labels=c(LETTERS[1:2]), rel_widths=c(2,1.2))
pmm = plot_grid( p_ana, p_cla, align="hv", nrow=1, labels=LETTERS[3:4])
pg3l = plot_grid( p1, p3, p2, p4, align="hv", nrow=2, labels=c(LETTERS[5:8]))
pg3 = plot_grid( pg3u, pmm, pg3l, nrow=3, align="hv", rel_heights=c(1,1,2))

# pdf( plot_fn, height=11, width=7.5 )
# print(pg3)
# dev.off()

CairoPDF( plot_fn, height=11, width=7.5 )
print(pg3)
dev.off()
