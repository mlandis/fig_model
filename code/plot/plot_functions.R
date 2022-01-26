library(HDInterval)
library(data.table)

gm = function(x) {
    n = length(x[ !is.na(x)])
    z = prod(x, na.rm=T)
    return(z^(1/n))
}

make_w = function(a_c, a_d, beta=0, sigma=1) {
    idx = which(a_d==0)
    a_c[idx] = a_c[idx] * sigma
    a_c = ( a_c / gm(a_c) )^beta
    return(a_c)
}

make_b = function(b_c, b_d, beta=0, sigma=1) {
    idx = which(b_d==0)
    b_c[idx] = b_c[idx] * sigma
    b_c = ( b_c / gm(b_c) )^beta
    return(b_c)
}

mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.lty  <- params("vertex", "lty")
    if (length(vertex.frame.lty) != 1 && !is.null(v)) {
      vertex.frame.lty <- vertex.frame.lty[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width, vertex.frame.lty,
         FUN=function(x, y, bg, fg, size, lwd, lty) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd, lty=lty,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

make_m1 = function(d, nn) {
    nr = length(nn)
    ret = list()
    colnames(d) = nn
    xm = sapply(d, median)
    h = apply(d, 2, HDInterval::hdi) #  function(z) {  hdi(z, allowSplit=F) })
    xl = h[1,]
    xu = h[2,]
    ret = list( median=xm, lower=xl, upper=xu)
    return(ret)
}

make_m2 = function(d, nn) {
    nr = length(nn)
    ret = list()
    xm = matrix( sapply(d, median), byrow=T, nrow=nr)
    xm = ( xm + t(xm) ) / 2 # mean
    h = apply(d, 2, HDInterval::hdi)
    xl = matrix( h[1,], byrow=T, nrow=nr ) 
    xu = matrix( h[2,], byrow=T, nrow=nr )
    diag(xm)=diag(xl)=diag(xu)=NA
    rownames(xm)=colnames(xm)=nn
    rownames(xl)=colnames(xl)=nn
    rownames(xu)=colnames(xu)=nn
    ret = list( median=xm, lower=xl, upper=xu)
    return(ret)
}

make_FIG_rate_graph = function(a, fn,
    colors=NULL,
    coords=NULL,
    node_order = NULL,
    node_size=20,
    node_frame_size = 1,
    edge_size=1,
    prop_size=T,
    add=F) {
    
    # make graph
    g = graph_from_adjacency_matrix(a, weighted=T, mode="undirected", diag=FALSE)
    
    # set colors
    node_colors="black"
    edge_colors="black"
    if (!is.null(colors)) {
        node_colors = diag(colors)
        edge_colors = c() #colors[ upper.tri(colors) ]
        for (i in 1:nrow(colors)) {
            for (j in i:ncol(colors)) {
                if (j > i) {
                    edge_colors = c( edge_colors, colors[i,j] )
                }
            }
        }
    }
    
    # set layout
    if (is.null(node_order)) {
        node_order = 1:nrow(a)
    }
    if (is.null(coords)) {
        coords <- layout.circle(g, node_order)
    }
    g$layout = coords
    
    # set graph object sizes
    V(g)$size               = node_size #* diag(a) / median(diag(a))
    V(g)$vertex.frame.width = node_frame_size * 2 * diag(a) / mean(diag(a))
    E(g)$weight             = edge_size * 1 * E(g)$weight / mean(E(g)$weight)
 
    # set graph object colors
    V(g)$color = node_colors
    E(g)$color = edge_colors
    
    # make plot
    if (!add) {
        pdf(height=3.25, width=3.25, file=fn)
        par(mar = rep(0.5, 4))
    }
    
    plot.igraph(g,
        vertex.shape="fcircle",
        vertex.color="white",
        vertex.size=V(g)$size,
        vertex.frame.width=V(g)$vertex.frame.width,
        vertex.frame.color=V(g)$color,
        vertex.lty=1,
        vertex.label=V(g)$name,
        vertex.label.color=V(g)$color,
        vertex.label.cex=0.55,
        vertex.label.family="Helvetica",
        edge.width=E(g)$weight,
        edge.color=adjustcolor( E(g)$color, alpha=0.6 ),
        edge.label.family="Helvetica",
        edge.label.color="black",
        add=add)
    
    if (!add) {
        dev.off()
    }
}
    
make_region_code_matrix = function( x ) {
    n = length(x)
    m = matrix(NA, nrow=n, ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) {
                m[i,j] = x[i]
            } else {
                m[i,j] = paste0(x[i], x[j])
            }
        }
    }
    return(m)
}

make_plot_dat = function( dat_mdl, dat_geo, num_regions ) {
    # process input data
    dat_plot = data.frame()
    
    list_idx = 1
    dat_plot_list = list()
    
    for (i in 1:nrow(dat_mdl)) {
        
        df = data.frame()
        
        dat_mdl_i = dat_mdl[i,]
        idx = dat_mdl_i$Iteration
        dat_geo_i = dat_geo[dat_geo$Iteration==idx, ]
        rho_d = dat_mdl_i$rho_d
        rho_b = dat_mdl_i$rho_b
        rho_w = dat_mdl_i$rho_w
        rho_e = dat_mdl_i$rho_e
        for (j in 1:num_regions) {
            s_w_j = paste0( "m_w.", j, "." )
            m_w_j = dat_geo_i[[ s_w_j ]]
            s_e_j = paste0( "m_e.", j, "." )
            m_e_j = dat_geo_i[[ s_e_j ]]
            for (k in 1:num_regions) {
                #if (j != k) {
                    s_b_jk = paste0( "m_b.", j, "..", k, "." )
                    m_b_jk = dat_geo_i[[ s_b_jk ]]
                    s_d_jk = paste0( "m_d.", j, "..", k, "." )
                    m_d_jk = dat_geo_i[[ s_d_jk ]]
                    tmp = c(idx, j, k, rho_w, rho_e, rho_b, rho_d, m_w_j, m_e_j, m_b_jk, m_d_jk )
                    df = rbind(df, tmp)
                    #names(tmp) =  c("iteration", "from_region", "to_region", "r_w", "r_e", "r_b", "r_d", "m_w_i", "m_e_i", "m_b_ij", "m_d_ij")
                    
                    #dat_plot_list[[list_idx]] = data.table( t(c(idx, j, k, r_w, r_e, r_b, r_d, m_w_j, m_e_j, m_b_jk, m_d_jk )) )
                    
                    #dat_plot = rbind( dat_plot, c(idx, j, k, r_w, r_e, r_b, r_d, m_w_j, m_e_j, m_b_jk, m_d_jk ))
                #}
            }
        }
        colnames(df) = c("iteration", "from_region", "to_region", "rho_w", "rho_e", "rho_b", "rho_d", "m_w_i", "m_e_i", "m_b_ij", "m_d_ij")
        #print(list_idx)
        dat_plot_list[[list_idx]] = df
        list_idx = list_idx + 1
    }
    #return(dat_plot_list)
    dat_plot = rbindlist(dat_plot_list)
    #colnames(dat_plot) = c("iteration", "from_region", "to_region", "r_w", "r_e", "r_b", "r_d", "m_w_i", "m_e_i", "m_b_ij", "m_d_ij")
    
    print("A")
    dat_plot$m_d_ij[ dat_plot$m_d_ij == Inf ] = 0
    dat_plot$m_b_ij[ dat_plot$m_b_ij == Inf ] = 0
    dat_plot$m_w_i[ dat_plot$m_w_i == Inf ] = 0
    dat_plot$m_e_i[ dat_plot$m_e_i == Inf ] = 0
    dat_plot$r_d_ij = dat_plot$rho_d*dat_plot$m_d_ij 
    dat_plot$r_b_ij = dat_plot$rho_b*dat_plot$m_b_ij
    dat_plot$r_w_i = dat_plot$rho_w*dat_plot$m_w_i 
    dat_plot$r_e_i = dat_plot$rho_e*dat_plot$m_e_i
    print("B")
    dat_plot$insular = (dat_plot$from_region >= 6) | (dat_plot$to_region >= 6)
    dat_plot$we_ratio = dat_plot$r_w_i / dat_plot$r_e_i
    dat_plot$bd_ratio = dat_plot$r_b_ij / dat_plot$r_d_ij
    dat_plot$dist_ij = 0
    dat_plot$size_i = 0
    print("C")
    #return(dat_plot)
    
    # note, flattened matrix access is column first, row second!
    from_region = dat_plot$from_region
    to_region = dat_plot$to_region
    from_to_region = (to_region-1)*nrow(dat_dist) + from_region
    
    for (i in 1:nrow(dat_dist)) {
        for (j in 1:ncol(dat_dist)) {
            k = (j-1)*nrow(dat_dist) + i
            cat(i, j, k, dat_dist[i,j], as.matrix(dat_dist)[k], "\n" )
        }
    }
    
    dat_plot$dist_ij = as.matrix(dat_dist)[ from_to_region ]
    #dist_tmp = dat_dist[ from_region, ]
    dat_plot$size_i = dat_size[from_region]
    # for (i in 1:nrow(dat_plot)) {
    #     src = dat_plot$from_region[i]
    #     dst = dat_plot$to_region[i]
    #     dat_plot$dist_ij[i] = dat_dist[src,dst]
    #     dat_plot$size_i[i] = dat_size[src]
    # }
    print("D")
    #return(from_to_region)
    return(dat_plot)   
}


make_plot_dat_old = function( dat_mdl, dat_geo, num_regions ) {
    # process input data
    dat_plot = data.frame()

    for (i in 1:nrow(dat_mdl)) {
        
        df = data.frame()
        
        dat_mdl_i = dat_mdl[i,]
        idx = dat_mdl_i$Iteration
        dat_geo_i = dat_geo[dat_geo$Iteration==idx, ]
        r_d = dat_mdl_i$r_d
        r_b = dat_mdl_i$r_b
        r_w = dat_mdl_i$r_w
        r_e = dat_mdl_i$r_e
        for (j in 1:num_regions) {
            s_w_j = paste0( "m_w.", j, "." )
            m_w_j = dat_geo_i[[ s_w_j ]]
            s_e_j = paste0( "m_e.", j, "." )
            m_e_j = dat_geo_i[[ s_e_j ]]
            for (k in 1:num_regions) {
                #if (j != k) {
                    s_b_jk = paste0( "m_b.", j, "..", k, "." )
                    m_b_jk = dat_geo_i[[ s_b_jk ]]
                    s_d_jk = paste0( "m_d.", j, "..", k, "." )
                    m_d_jk = dat_geo_i[[ s_d_jk ]]
                    # tmp = c(idx, j, k, r_w, r_e, r_b, r_d, m_w_j, m_e_j, m_b_jk, m_d_jk )
                    # df = rbind(df, tmp)
                    # #names(tmp) =  c("iteration", "from_region", "to_region", "r_w", "r_e", "r_b", "r_d", "m_w_i", "m_e_i", "m_b_ij", "m_d_ij")
                    # 
                    # #dat_plot_list[[list_idx]] = data.table( t(c(idx, j, k, r_w, r_e, r_b, r_d, m_w_j, m_e_j, m_b_jk, m_d_jk )) )
                    # 
                    dat_plot = rbind( dat_plot, c(idx, j, k, r_w, r_e, r_b, r_d, m_w_j, m_e_j, m_b_jk, m_d_jk ))
                #}
            }
        }
        #colnames(df) = c("iteration", "from_region", "to_region", "r_w", "r_e", "r_b", "r_d", "m_w_i", "m_e_i", "m_b_ij", "m_d_ij")
        #print(list_idx)
        #dat_plot_list[[list_idx]] = df
        #list_idx = list_idx + 1
    }
    #return(dat_plot_list)
    #dat_plot = rbindlist(dat_plot_list)
    colnames(dat_plot) = c("iteration", "from_region", "to_region", "r_w", "r_e", "r_b", "r_d", "m_w_i", "m_e_i", "m_b_ij", "m_d_ij")
    
    print("A")
    
    dat_plot$r_d_ij = dat_plot$r_d*dat_plot$m_d_ij 
    dat_plot$r_b_ij = dat_plot$r_b*dat_plot$m_b_ij
    dat_plot$r_w_i = dat_plot$r_w*dat_plot$m_w_i 
    dat_plot$r_e_i = dat_plot$r_e*dat_plot$m_e_i
    print("B")
    dat_plot$insular = (dat_plot$from_region >= 6) | (dat_plot$to_region >= 6)
    dat_plot$we_ratio = dat_plot$r_w_i / dat_plot$r_e_i
    dat_plot$bd_ratio = dat_plot$r_b_ij / dat_plot$r_d_ij
    dat_plot$dist_ij = 0
    dat_plot$size_i = 0
    print("C")
    #return(dat_plot)
    # 
    # from_region = dat_plot$from_region
    # from_to_region = (nrow(dat_dist)-1)*from_region + dat_plot$to_region
    # 
    # dat_plot$dist_ij = as.matrix(dat_dist)[ from_to_region ]
    # #dist_tmp = dat_dist[ from_region, ]
    # dat_plot$size_i = dat_size[from_region]
    for (i in 1:nrow(dat_plot)) {
        src = dat_plot$from_region[i]
        dst = dat_plot$to_region[i]
        dat_plot$dist_ij[i] = dat_dist[src,dst]
        dat_plot$size_i[i] = dat_size[src]
    }
    print("D")
    #return(from_to_region)
    return(dat_plot)   
}

make_plot_dat2 = function(dat_plot) {
    
    dat_plot_2 = dat_plot %>% group_by(from_region, to_region)  %>%
        summarise(
            mean_wi = mean(r_w_i, na.rm=T),
            hdil_wi = HDInterval::hdi(r_w_i)[1],
            hdiu_wi = HDInterval::hdi(r_w_i)[2],
            mean_ei = mean(r_e_i, na.rm=T),
            hdil_ei = HDInterval::hdi(r_e_i)[1],
            hdiu_ei = HDInterval::hdi(r_e_i)[2],
            mean_bij = mean(r_b_ij, na.rm=T),
            hdil_bij = HDInterval::hdi(r_b_ij)[1],
            hdiu_bij = HDInterval::hdi(r_b_ij)[2],
            mean_dij = mean(r_d_ij, na.rm=T),
            hdil_dij = HDInterval::hdi(r_d_ij)[1],
            hdiu_dij = HDInterval::hdi(r_d_ij)[2],
            mean_wei = mean(r_w_i/r_e_i, na.rm=T),
            hdil_wei = HDInterval::hdi(r_w_i/r_e_i)[1],
            hdiu_wei = HDInterval::hdi(r_w_i/r_e_i)[2],
            mean_bdij = mean( r_b_ij/r_d_ij, na.rm=T),
            hdil_bdij = HDInterval::hdi(r_b_ij/r_d_ij)[1],
            hdiu_bdij = HDInterval::hdi(r_b_ij/r_d_ij)[2]
            )
    
    
    dat_plot_2$insular      = (dat_plot_2$from_region >= 6) | (dat_plot_2$to_region >= 6)
    dat_plot_2$from_insular = (dat_plot_2$from_region >= 6)
    dat_plot_2$to_insular   = (dat_plot_2$to_region >= 6)
    dat_plot_2$from_lbl     = LETTERS[dat_plot_2$from_region] 
    dat_plot_2$to_lbl       = LETTERS[dat_plot_2$to_region]
    
    
    lbl2 = c()
    for (i in 1:nrow(dat_plot_2)) {
        lbl2[i] = lbl[ dat_plot_2$from_region[i], dat_plot_2$to_region[i] ]
    }
    dat_plot_2$lbl = lbl2
    
    
    
    dat1 = dat_plot_2[ dat_plot_2$from_region < dat_plot_2$to_region, ]
    dat2 = dat_plot_2[ dat_plot_2$from_region >= dat_plot_2$to_region, ]
    dat12 = rbind(dat1, dat2)
    
    
    dat12$dist_ij = rep(0, nrow(dat12))
    dat12$size_i  = rep(0, nrow(dat12))
    for (k in 1:nrow(dat12)) {
        i = dat12$from_region[k]
        j = dat12$to_region[k]
        #print(dat_dist[i,j])
        dat12[k,]$dist_ij = dat_dist[i,j]
        dat12[k,]$size_i = dat_size[i]
    }
    
    return(dat12)
}

make_lm_plot_dat = function( dat_dist, dat_size, dat_plot ) {
    
    # get lm_matrix
    dd = c()
    dd_cont = c()
    dd_ins = c()
    for (i in 1:nrow(dat_dist)) {
        for (j in 1:ncol(dat_dist)) {
            dd = c(dd, dat_dist[i,j]) 
            if (i <= 5 && j <= 5) {
                dd_cont = c(dd_cont, dat_dist[i,j])
            } else {
                dd_ins = c(dd_ins, dat_dist[i,j])
            }
        }
    }
    aa_cont = dat_size[1:5]
    aa_ins = dat_size[6:9]
    dd_cont = dd_cont[ dd_cont!= 0 ]
    dd_ins = dd_ins[ dd_ins!= 0 ]
    
    dat_plot_it = unique(dat_plot$iteration)
    mtx_lm = matrix(0, nrow=length(dat_plot_it), ncol=8)
    colnames(mtx_lm) = c("intercept_bd_cont", "slope_bd_cont","intercept_bd_ins", "slope_bd_ins","intercept_we_cont", "slope_we_cont","intercept_we_ins", "slope_we_ins")
    for (i in 1:length(dat_plot_it)) {
        # get iteration
        it = dat_plot_it[i]
        dat_plot_i = dat_plot[ dat_plot$iteration==it, ]
        df_dist_cont = dat_plot_i[ dat_plot_i$from_region!=dat_plot_i$to_region & !dat_plot_i$insular, ]
        df_size_cont = dat_plot_i[ dat_plot_i$from_region==dat_plot_i$to_region & !dat_plot_i$insular, ]
        df_dist_ins  = dat_plot_i[ dat_plot_i$from_region!=dat_plot_i$to_region &  dat_plot_i$insular, ]
        df_size_ins  = dat_plot_i[ dat_plot_i$from_region==dat_plot_i$to_region &  dat_plot_i$insular, ]
        
        # collect lms
        xy_bd_cont = matrix( c( log10(dd_cont), log10(df_dist_cont$bd_ratio)), nrow=nrow(df_dist_cont))
        xy_we_cont = matrix( c( log10(aa_cont), log10(df_size_cont$we_ratio)), nrow=nrow(df_size_cont))
        xy_bd_ins  = matrix( c( log10(dd_ins),  log10(df_dist_ins$bd_ratio)),  nrow=nrow(df_dist_ins))
        xy_we_ins  = matrix( c( log10(aa_ins),  log10(df_size_ins$we_ratio)),  nrow=nrow(df_size_ins))
        
        # generate linear models across posterior samples
        idx1 = 1; idx2 = 2
        #idx1 = 2; idx2 = 1
        lm_bd_cont = lm( xy_bd_cont[,idx1] ~ xy_bd_cont[,idx2] )
        mtx_lm[i, c("intercept_bd_cont","slope_bd_cont")] = coef(lm_bd_cont)
        lm_bd_ins = lm( xy_bd_ins[,idx1] ~ xy_bd_ins[,idx2] )
        mtx_lm[i, c("intercept_bd_ins","slope_bd_ins")] = coef(lm_bd_ins)
        lm_we_cont = lm( xy_we_cont[,idx1] ~ xy_we_cont[,idx2] )
        mtx_lm[i, c("intercept_we_cont","slope_we_cont")] = coef(lm_we_cont)
        lm_we_ins = lm( xy_we_ins[,idx1] ~ xy_we_ins[,idx2] )
        mtx_lm[i, c("intercept_we_ins","slope_we_ins")] = coef(lm_we_ins)
    }
    
    
    mtx_lm2 = rbind( cbind( mtx_lm[,c("intercept_bd_cont","slope_bd_cont")], "bd", "continental" ),
                     cbind( mtx_lm[,c("intercept_we_cont","slope_we_cont")], "we", "continental" ),
                     cbind( mtx_lm[,c("intercept_bd_ins","slope_bd_ins")], "bd", "insular" ),
                     cbind( mtx_lm[,c("intercept_we_ins","slope_we_ins")], "we", "insular" ) )
        
    df_lm = data.frame( mtx_lm2 )
    colnames(df_lm) = c("intercept","slope","ratio","type")
    df_lm$intercept = as.numeric(df_lm$intercept)
    df_lm$slope = as.numeric(df_lm$slope)
    df_lm$insular = df_lm$type=="insular"
    
    df_lm$x0 = df_lm$intercept #(0 - df_lm$intercept) / df_lm$slope
    df_lm$y0 = 0
    df_lm$slope_p = 1/df_lm$slope #(0-df_lm$intercept) / df_lm$slope
    df_lm$intercept_p = -df_lm$intercept/df_lm$slope
    
    return(df_lm)

}

summarize_lm_dat = function(df_lm) {
    
    med_bd_cont = mean( df_lm$x0[ df_lm$type=="continental" & df_lm$ratio=="bd" ] )
    med_bd_ins  = mean( df_lm$x0[ df_lm$type=="insular" & df_lm$ratio=="bd" ] )
    med_we_cont = mean( df_lm$x0[ df_lm$type=="continental" & df_lm$ratio=="we" ] )
    med_we_ins  = mean( df_lm$x0[ df_lm$type=="insular" & df_lm$ratio=="we" ] )
    hpd_bd_cont = hdi( df_lm$x0[ df_lm$type=="continental" & df_lm$ratio=="bd" ] )
    hpd_bd_ins  = hdi( df_lm$x0[ df_lm$type=="insular"     & df_lm$ratio=="bd" ] )
    hpd_we_cont = hdi( df_lm$x0[ df_lm$type=="continental" & df_lm$ratio=="we" ] )
    hpd_we_ins  = hdi( df_lm$x0[ df_lm$type=="insular"     & df_lm$ratio=="we" ] )
    
    cat( "MRCD [BD, cont.] = ", 10^med_bd_cont, " [ ", 10^hpd_bd_cont[1], "  ", 10^hpd_bd_cont[2], " ]\n")
    cat( "MRCD [BD, ins. ] = ", 10^med_bd_ins,  " [ ", 10^hpd_bd_ins[1],  "  ", 10^hpd_bd_ins[2],  " ]\n")
    cat( "MRCD [WE, cont.] = ", 10^med_we_cont, " [ ", 10^hpd_we_cont[1], "  ", 10^hpd_we_cont[2], " ]\n")
    cat( "MRCD [WE, ins.]  = ", 10^med_we_ins,  " [ ", 10^hpd_we_ins[1],  "  ", 10^hpd_we_ins[2],  " ]\n")
    
    ret = data.frame()
    ret = rbind( ret, c("continental", "bd", med_bd_cont, hpd_bd_cont[1], hpd_bd_cont[2]) )
    ret = rbind( ret, c("continental", "we", med_we_cont, hpd_we_cont[1], hpd_we_cont[2]) )
    ret = rbind( ret, c("insular", "bd", med_bd_ins, hpd_bd_ins[1], hpd_bd_ins[2]) )
    ret = rbind( ret, c("insular", "we", med_we_ins, hpd_we_ins[1], hpd_we_ins[2]) )
    
    df = data.frame( type=ret[,1],
                     ratio=ret[,2],
                     log10_median=as.numeric(ret[,3]),
                     log10_lower=as.numeric(ret[,4]),
                     log10_upper=as.numeric(ret[,5]) )
    #colnames(ret) = c("type","ratio","log10_median","log10_lower","log10_upper")
    #ret[,3:5] = as.numeric( ret[,3:5] )
    return(df)   
}

make_dat_count = function( dat_tip, dat_size, region_names ) {

    df = dat_tip %>% count(state)
    df$range = dat_states[match(df$state, dat_states$state),2]
    lbl = c()
    for (i in 1:nrow(df)) {
        x = df$range[i]
        lbl[i] = paste( region_names[ which( unlist(strsplit(x, split=""))==1 ) ], collapse="" )
    }
    df$lbl = lbl
    df$lbl = factor(df$lbl, ordered=T, levels=lbl[order(nchar(lbl))])
    
    df$mean_size = rep(0, nrow(df))
    reg_idx2 = lapply(df$range, function(x) { which(strsplit(x, "")[[1]]=="1") } )
    for (i in 1:length(reg_idx2)) {
        j = reg_idx2[[i]]
        df$mean_size[i] = mean( dat_size[ c(j) ] )
    }
    df$widespread = sapply( as.vector(df$lbl), nchar ) > 1
    
    return(df)
}    

make_dat_count_region = function(dat_count) {
    
    dat_count_single = dat_count[ !dat_count$widespread, ]
    dat_count_widespread = dat_count[ dat_count$widespread, ]
    for (i in 1:nrow(dat_count_single)) {
        for (j in 1:nrow(dat_count_widespread)) {
            contains_region_i = dat_count_single$lbl[i] %in% strsplit( as.vector(dat_count_widespread$lbl)[j], "" )[[1]]
            if (contains_region_i) {
                dat_count_single$n[i] = dat_count_single$n[i] + dat_count_widespread$n[j]
            }
        }
    }
    
    ret = dat_count_single
    ret$lbl = as.vector(ret$lbl)
    for (i in 1:nrow(dat_count_single)) {
        di = dat_count[ dat_count$lbl == dat_count_single$lbl[i], ]
        if ( di$n != dat_count_single$n[i] ) {
            di$lbl = paste0(di$lbl, "'")
            ret = rbind(ret, di)
        }
        ret$lbl[i] = paste0(ret$lbl[i], " ")
    }
    ret$single = unlist( lapply( strsplit(ret$lbl, ""), function(x) { if (x[2]=="'") { 0 } else { 1 } } ) )
    
    return(ret)
}


# converts bit string to region-set
# for example, converts "1100" to "A+B"
bitstr_to_regset = function(s,n) {
    rs = n[ which(as.integer(strsplit(s, split="")[[1]])==1) ]
    rs = paste(rs, collapse="+")
    return(rs)
}

# process states and ranges
get_used_states = function(phy, bitset) {
    dat = phy@data[, c("end_state_1", "end_state_2", "end_state_3", "start_state_1", "start_state_2", "start_state_3")]
    vec = as.numeric( unlist(dat) )
    vec = sort(unique(vec[ !is.na(vec) ]))
    #print(vec)
    
    # container for ranges sorted by # regions
    n_reg = length( strsplit( bitset[ vec[1]+1 ], split="" )[[1]] )
    ret = list()
    for (i in 1:n_reg) {
        ret[[i]] = list()
    }

    for (i in 1:length(vec)) {
        range_n = vec[i] + 1 # convert to base-1 range-states
        range_01 = bitset[ range_n ]
        idx = sum( as.numeric( strsplit( range_01, split="" )[[1]] ) )
        #print(idx)
        ret[[idx]] = c( ret[[idx]], range_n )
    }
    
    return(ret)
    #return(sort(unique(vec)))
}



# adds epoch boxes to plots
add_epoch_times <- function( p, max_age, dy_bars, dy_text ) {
    
    max_x = max(p$data$x)
    max_y = max(p$data$y)
    epoch_names = c("Late\nCretaceous","Paleocene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")
  

    x_pos = max_x-c(max_age, 65, 56, 48, 33.9, 23, 16, 5.3, 0)
    y_pos = rep(max_y, length(x_pos))
    x_pos_mid = ( x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)] ) / 2 

    for (k in 2:(length(x_pos))) {
        box_col = "gray92"
        if (k %% 2 == 0) box_col = "white"
        box = geom_rect( xmin=x_pos[k-1], xmax=x_pos[k], ymin=dy_bars, ymax=y_pos[k], fill=box_col )
        p = append_layers(p, box, position = "bottom")
    }
    for (k in 1:length(epoch_names)) {
        p = p + annotate( geom="text", label=epoch_names[k], x=x_pos_mid[k], y=dy_text, hjust=0.5, size=3.25)
    }
    return(p)

}

# add transparency to color
t_col <- function(color, percent = 50, name = NULL) {
    # color = color name
    # percent = % transparency
    # name = an optional name for the color

    ## Get RGB values for named color
    rgb.val <- col2rgb(color)

    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1],
                 rgb.val[2],
                 rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)

    ## Save the color
    invisible(t.col)
}
