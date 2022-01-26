source("plot_functions.R")

library(rgdal)
library(maptools)


# filesystem input/output
fp          = "../../"
dat_fp      = paste0( fp, "output/anolis/" )
dat_fn      = paste0(dat_fp, "out.geo_features.log")
regions_fn  = paste0(fp, "code/plot/anolis_regions.csv")
plot_fn     = paste0(fp, "code/plot/fig_anolis_n9_graph")

# read in data
dat         = read.table(dat_fn, header=T, sep="\t")
region_info = read.csv(regions_fn, sep=",", header=T)
burn        = 0.2*nrow(dat)
dat         = dat[burn:nrow(dat), ]
cn          = colnames(dat)

# various labels and coordinates for plotting
nn   = region_info$code
nr   = length(nn)
lat  = region_info$latitude
lon  = region_info$longitude
size  = 2

# where to write labels (offset from ~centroids)
coords  = matrix( c(lon, lat), ncol=2 )
#coords_text = coords + 1.2*matrix(c(
#    -2.5,+2.5,
#    -2.5,-1.0,
#    -2.5,+2.5,
#    -3.0,-1.0,
#    +3.0,-1.0,
#    +3.0,+2.0,
#    -2.5,+2.5,
#    +3.0,+2.0,
#    +2.0,+3.0), ncol=2, byrow=T)

# colors for plotting
color_base = c("#9E248F","#00AEEF")
color_idx = matrix( 1, nrow=nr, ncol=nr )
color_idx[6:9,] = 2
color_idx[,6:9] = 2
colors = matrix( color_base[color_idx], nrow=nr, ncol=nr )

# extract posteriors for m_p functions
dat[ dat==Inf ] = 0
dat_d = dat[, grepl("m_d", cn)]
dat_e = dat[, grepl("m_e", cn)]
dat_w = dat[, grepl("m_w", cn)]
dat_b = dat[, grepl("m_b", cn)]

# get means for within/between m-functions
ld = make_m2(dat_d, nn)
le = make_m1(dat_e, nn)
lw = make_m1(dat_w, nn)
lb = make_m2(dat_b, nn)

# generate cladogenetic graph (w.r. speciation, b.r. speciation)
clado = lb$median; diag(clado) = lw$median

# generate anagenetic graph (extinction, dispersal)
ana = ld$median; diag(ana) = le$median

# read in shape files
shp_raw     = readOGR( dsn="shapefiles", layer="anolis_regions" )
shp_all_tmp = as(unionSpatialPolygons( shp_raw, c(1,2,3,4,5,6,7,8,9,10,11,12,13) ), "SpatialPolygonsDataFrame")
shp_all     = shp_all_tmp[ c(6,1,10,11,13,12,2,3,4,5,7,8,9), ]
shp         = as(unionSpatialPolygons( shp_all, c(1,1,2,3,3,4,5,6,6,8,7,7,9) ), "SpatialPolygonsDataFrame")
op          = par(mar = rep(0, 4), oma=rep(0,4))
col         = c("firebrick3","gold","deepskyblue","darkorchid","forestgreen","blue","cyan","magenta","orange")
col         = sapply( col, function(x) { t_col(x, 85) })
xrange=range(coords[,1]); xrange[1]=xrange[1]-5;xrange[2]=xrange[2]+5;
yrange=range(coords[,2]); yrange[1]=yrange[1]-5;yrange[2]=yrange[2]+5;

# generate cladogenetic plot (within-region and between-region speciation)
pdf("fig_geo_BW.pdf", height=7, width=7)
plot(shp, border="#999999", lwd=0.5, col=col, xaxs="i", yaxs="i", xlim=xrange, ylim=yrange)
par(op)
for (i in 1:9) {
    for (j in 1:9) {
        if (i > j) {
            col_ij = rgb( t(col2rgb(colors[i,j])/255), alpha=0.6 )
            segments( x0=coords[i,1], y0=coords[i,2],  x1=coords[j,1], y1=coords[j,2], col=col_ij , lwd=lb$median[i,j] ) 
        }
    }
}

points( coords[,1], coords[,2], cex=size*3, pch=21, lwd=lw$median*2, col="white", bg="white" )
points( coords[,1], coords[,2], cex=size*3, pch=21, lwd=lw$median*2, col=diag(colors) )
#text( coords_text[,1], coords_text[,2], cex=1.5, label=LETTERS[1:9], col="black" )
text( coords[,1], coords[,2], cex=2.2, label=LETTERS[1:9], col="black" )
dev.off()

# generate anagenetic plot (dispersal and extinction)
op <- par(mar = rep(0, 4), oma=rep(0,4))
pdf("fig_geo_DE.pdf", height=7, width=7)
plot(shp, border="#999999", lwd=0.5, col=col, xaxs="i", yaxs="i", xlim=xrange, ylim=yrange)
for (i in 1:9) {
    for (j in 1:9) {
        if (i > j) {
            col_ij = rgb( t(col2rgb(colors[i,j])/255), alpha=0.6 )
            segments( x0=coords[i,1], y0=coords[i,2],  x1=coords[j,1], y1=coords[j,2], col=col_ij , lwd=ld$median[i,j] ) 
        }
    }
}
points( coords[,1], coords[,2], cex=size*3, pch=21, lwd=le$median*2, col="white", bg="white" )
points( coords[,1], coords[,2], cex=size*3, pch=21, lwd=le$median*2, col=diag(colors) )
#text( coords_text[,1], coords_text[,2], cex=1.4, label=LETTERS[1:9], col="black" )
text( coords[,1], coords[,2], cex=2.2, label=LETTERS[1:9], col="black" )
dev.off()
