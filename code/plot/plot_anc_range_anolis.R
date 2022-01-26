# libraries
library(cowplot)
library(ggtree)
library(treeio)
library(ggimage)

# helper functions
source("plot_functions.R")
source("plot_ancestral_states.R")

# Files
fp = "../../"
out_str = "out"
dat_fp  = paste0(fp, "data/anolis/")
out_fp  = paste(fp, "output/anolis/", sep="")
plot_fp = paste(fp, "code/plot/", sep="")
tree_fn = paste(out_fp, out_str, ".ase_marginal.tre", sep="")
plot_fn = paste(plot_fp, "fig_anc_range.pdf",sep="")
lbl_fn  = paste0(fp, "code/plot/state_labels_n", num_reg, ".txt")
dat_fn  = paste0(dat_fp, "anolis_nr9_ns383.range_table.tsv")

# other info
num_reg = 9; num_states=255
reg_lbl = LETTERS[1:9]
long_names =
"A: NA+Upper CA
B: Lower CA
C: Chocó+Caribe
D: Andes
E: Amazonia
F: LA+PR        
G: Cuba+Cay.+Jam.
H: Hispaniola
I: Bahamas"

# read and format tip data
tip_dat = read.csv(dat_fn, sep=",", header=T, row.names=1, colClasses=c("character"))
for (i in 1:nrow(tip_dat)) {
    for (j in 1:ncol(tip_dat)) {
        if (tip_dat[i,j] == 0) {
            tip_dat[i,j] = NA
        }
    }
}
for (i in 1:ncol(tip_dat)) {
    tip_dat[ !is.na(tip_dat[,i]), i ] = reg_lbl[i]
    tip_dat[,i] = factor(tip_dat[,i])
} 
colnames(tip_dat)=reg_lbl
tip_dat = data.frame(tip_dat)

# read tree
phy = read.nexus(tree_fn)
phy_state = read.beast(tree_fn)

# read and format region label information
dat_lbl_raw = read.csv( lbl_fn, sep="\t",  colClasses=c('character'), header=F )[ 1:num_states, 2]
dat_lbl = sapply( dat_lbl_raw, function(x) { bitstr_to_regset(x, reg_lbl)  })
names(dat_lbl) = NULL

# what states are used in the ancestral state tree?
used_states = get_used_states(phy_state, dat_lbl_raw)
#num_used_states = length(used_states) 
st_lbl = list()
for (i in 1:num_states) {
    st_lbl[[ as.character(i-1) ]] = dat_lbl[i]
}

# set up colors
bg_color_all = rep("black", 255)
names(bg_color_all) =  unlist(st_lbl )
bg_color_all[ unlist(used_states[[1]]) ] = c("firebrick3","gold","deepskyblue","darkorchid","forestgreen","blue","cyan","magenta","orange")
bg_color_all[ unlist(used_states[[2]]) ] = darken( rainbow( n=length(used_states[[2]]), start=0, end=0.95 ), amount=-0.3 )
bg_color_all[ unlist(used_states[[3]]) ] = darken( rainbow( n=length(used_states[[3]]), start=0, end=0.95 ), amount=0.3 )
bg_color_all[ unlist(used_states[[4]]) ] = darken( rainbow( n=length(used_states[[4]]), start=0, end=0.75 ), amount=0.5 )
st_colors = c( bg_color_all )
names(st_colors) = unlist( st_lbl )

# Build figure
cat("Processing...\n")
zz=plot_ancestral_states(tree_file=tree_fn,
                      include_start_states=T,
                      shoulder_label_size=0,
                      summary_statistic="PieRange",
                      state_labels=st_lbl,
                      state_colors=st_colors,
                      node_label_size=0,
                      node_size_range=c(0.5,2.0),
                      node_label_nudge_x=0.5,
                      tip_node_size=2,
                      tip_label_size=2.5,
                      tip_label_offset=2,
                      xlim_visible=c(0,70),
                      shoulder_label_nudge_x=-0.1,
                      show_posterior_legend=T,
                      node_pie_diameter=4.5,
                      tip_pie_diameter=3,
                      pie_nudge_x=0.2,
                      pie_nudge_y=0.2,
                      alpha=1)

# create copy of plotting object so in case you want to modify ggplot style
p2 = zz


# get dimensions
x_height = max(p2$data$x)
x_new_root = 90
x_offset = x_new_root - x_height
x_label = 10
dx = x_new_root %% 10
n_node = length(phy$tip.label)

# set coordinates
p2 = p2 + geom_rootedge(rootedge=5)
p2 = p2 + theme_tree2()
p2 = p2 + coord_cartesian(xlim = c(-(x_offset+3),x_height+x_label),
                          ylim=c(-7,n_node+20), expand=F)
p2 = p2 + scale_x_continuous(breaks = seq(0,x_new_root,10)-x_offset+dx,
                             labels = rev(seq(0,x_new_root,10)))
p2 = p2 + labs(x="Age (Ma)")

# correct tip labels
tip_lbl = p2$data$label
tip_idx = which(!is.na(tip_lbl))
tip_lbl = tip_lbl[tip_idx]
p2$data$label[tip_idx] = tip_lbl

# adjust tip node sizes
p2 = p2 + geom_tippoint(aes(color=end_state_1), size=1.6)

# correct legend fill
p2 = p2 + theme(plot.title = element_text(size=6, face="bold"),
                legend.position=c(0.02, 0.98),
                legend.justification=c("left", "top"),
                legend.key = element_blank(),
                legend.text=element_text(size=6),
                axis.line = element_line(colour = "black"))

p2 = p2 + guides(colour=guide_legend(title="Range", override.aes=list(size=2)))

# set geological ages
p2 = add_epoch_times(p2, x_new_root,  dy_bars=-7, dy_text=-3)

pq = p2
p2 = p2 + annotate(geom="text", label=long_names,  x=-20, y=315, hjust=0)


# Plot grid of 0/1 range character at tips
p3 = ggtree( phy_state, color="white")
p3 = gheatmap(p3, tip_dat, offset=1, #0.018,
               width=1.5, font.size=3,
               colnames_angle=0,
               colnames_level=reg_lbl,
               colnames_position="top",
               colnames_offset_y=1)
p3 = p3 + scale_fill_manual(breaks=reg_lbl, values=st_colors[1:9], na.value="transparent")
p3 = p3 + coord_cartesian(xlim = c(x_height,x_height*4),
                          ylim=c(-7,n_node+20), expand=F)
p3 = p3 + guides(fill = FALSE)
p3 = p3 + theme_tree2(axis.line.x=element_line(color=NA), axis.ticks.x=element_line(color=NA), axis.text.x=element_text(color=NA))
p3 = p3 + xlab(" ")



pg = plot_grid(p2, p3, ncol=2, rel_widths=c(1,0.2))


pdf(file=plot_fn, height=30, width=12)
print(pg)
dev.off()
