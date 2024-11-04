# Figure 1
# a. map sample
# b. pca
# c. phylogeny
# d. ngsadm

# 0. env all ----
library(tidyverse)
library(ggrepel)
library(glue)

library(sf)
library(rnaturalearth)

source('./lib_meta.R')
meta_sample = load_meta_gdoc()
list_dirs = get_dir_list()
pop_legend = get_legend(meta_sample)


source('./lib_utility.R')

# a. map v1 ----
## 1. data ----
world_country = ne_countries(scale = "medium", returnclass = "sf")
world_coast = ne_coastline(scale = "medium", returnclass = "sf")

f_id_se   = glue('{list_dirs$prj_dir}3-ngsadm/add_MCZ/list.id.qc.rel')
df_id_se   = read_tsv(f_id_se, col_names = c("id"))

## 2. zanzibar ----
p_zanzibar =
    meta_sample %>% 
    filter(Plot_ID %in% df_id_se$id) %>% 
    filter(Plot_spp %in% c("African leopard")) %>% 
    
    ggplot() +  
    # plot map
    geom_sf(data=world_country,
            fill="#fef0d9", color = "#fef0d9", linewidth = 0.1) +
    geom_sf(data=world_coast, linewidth = 0.1) +
    
    # plot sample point
    geom_point(
        aes(x=Lon, y=Lat, shape=Plot_pop, fill=Plot_pop),
        size=1, position=position_jitter(h=0,w=0),
        stroke=0.2) +
    
    ##  coord
    coord_sf(xlim=c(38,41), ylim=c(-7,-4)) +            # zanzibar
    
    scale_shape_manual(breaks = pop_legend$Plot_pop,
                       values = pop_legend$Plot_region_shape) +
    scale_fill_manual(breaks = pop_legend$Plot_pop,
                      values = pop_legend$Plot_region_color) +
    
    scale_x_continuous(breaks = c(38,  40 )) + 
    scale_y_continuous(breaks = c(-7, -4)) +
    
    guides(shape = guide_legend(
        ncol=1, override.aes = list(size = 2),
    )) +
    
    # scale_fill_manual(values = fill9) +
    theme(
        panel.background = element_rect(
        fill = '#abd8ea', colour = "black",
        linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        axis.ticks.length = unit(1, "pt"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 3))
p_zanzibar

## 3. africa ----
p_a = meta_sample %>% 
    filter(Plot_ID %in% df_id_se$id) %>% 
    filter(Plot_spp %in% c("African leopard")) %>% 
    
    ggplot() +  
    # plot map
    geom_sf(data=world_country,
            fill="#fef0d9", color = "#fef0d9", linewidth = 0.25) +
    geom_sf(data=world_coast, linewidth = 0.25) +
    
    # plot sample point
    geom_point(
        aes(x=Lon, y=Lat, shape=Plot_pop, fill=Plot_pop),
               size=1, position=position_jitter(h=0,w=0),
        stroke=0.2) +
    
    ##  coord
    # coord_sf(xlim=c(-20,140), ylim=c(-40,50)) +       # full
    coord_sf(xlim=c(-20,50), ylim=c(-40,40)) +          # africa

    scale_shape_manual(breaks = pop_legend$Plot_pop,
                       values = pop_legend$Plot_region_shape) +
    scale_fill_manual(breaks = pop_legend$Plot_pop,
                      values = pop_legend$Plot_region_color) +
    
    guides(shape = guide_legend(
        ncol=1, override.aes = list(size = 1),
        )) +
    
    annotation_custom(
        grob = ggplotGrob(p_zanzibar),
        xmin = -20, xmax = 5,
        ymin = -30, ymax = -10) +
    
    # scale_fill_manual(values = fill9) +
    theme(
        panel.background = element_rect(
            fill = '#abd8ea', colour = "black",
            linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.key.height = unit(6, "pt"),
        legend.key.width = unit(4, "pt"),
        legend.key = element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 5))

p_a
# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 400
fname   = 'Fig1/Panel_A.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=400, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))

# a. map v2 ----
## 1. data ----
world_country = ne_countries(scale = "medium", returnclass = "sf")
world_coast = ne_coastline(scale = "medium", returnclass = "sf")

f_id_se   = glue('{list_dirs$prj_dir}3-ngsadm/add_MCZ/list.id.qc.rel')
df_id_se   = read_tsv(f_id_se, col_names = c("id"))

## 2. zanzibar ----
p_zanzibar =
    meta_sample %>% 
    filter(Plot_ID %in% df_id_se$id) %>% 
    filter(Plot_spp %in% c("African leopard")) %>% 
    
    ggplot() +  
    # plot map
    geom_sf(data=world_country,
            fill="#fef0d9", color = "#fef0d9", linewidth = 0.1) +
    geom_sf(data=world_coast, linewidth = 0.1) +
    
    # plot sample point
    geom_point(
        aes(x=Lon, y=Lat, shape=Plot_pop, fill=Plot_pop),
        size=1, position=position_jitter(h=0,w=0),
        stroke=0.2) +
    
    ##  coord
    coord_sf(xlim=c(35,41), ylim=c(-10,-1)) +            # zanzibar
    
    scale_shape_manual(breaks = pop_legend$Plot_pop,
                       values = pop_legend$Plot_region_shape) +
    scale_fill_manual(breaks = pop_legend$Plot_pop,
                      values = pop_legend$Plot_region_color) +
    
    scale_x_continuous(breaks = c(38,  40 )) + 
    scale_y_continuous(breaks = c(-7, -4)) +
    
    guides(shape = guide_legend(
        ncol=1, override.aes = list(size = 2),
    )) +
    
    # scale_fill_manual(values = fill9) +
    theme(
        panel.background = element_rect(
            fill = '#abd8ea', colour = "black",
            linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        axis.ticks.length = unit(1, "pt"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 3))
p_zanzibar

## 3. africa ----
p_a = meta_sample %>% 
    filter(Plot_ID %in% df_id_se$id) %>% 
    filter(Plot_spp %in% c("African leopard")) %>% 
    
    ggplot() +  
    # plot map
    geom_sf(data=world_country,
            fill="#fef0d9", color = "#fef0d9", linewidth = 0.25) +
    geom_sf(data=world_coast, linewidth = 0.25) +
    
    # plot sample point
    geom_point(
        aes(x=Lon, y=Lat, shape=Plot_pop, fill=Plot_pop),
        size=1, position=position_jitter(h=0,w=0),
        stroke=0.2) +
    
    ##  coord
    # coord_sf(xlim=c(-20,140), ylim=c(-40,50)) +       # full
    coord_sf(xlim=c(-20,50), ylim=c(-40,40)) +          # africa
    
    scale_shape_manual(breaks = pop_legend$Plot_pop,
                       values = pop_legend$Plot_region_shape) +
    scale_fill_manual(breaks = pop_legend$Plot_pop,
                      values = pop_legend$Plot_region_color) +
    
    guides(shape = guide_legend(
        ncol=1, override.aes = list(size = 1),
    )) +
    
    annotation_custom(
        grob = ggplotGrob(p_zanzibar),
        xmin = -20, xmax = 5,
        ymin = -30, ymax = -10) +
    
    # scale_fill_manual(values = fill9) +
    theme(
        panel.background = element_rect(
            fill = '#abd8ea', colour = "black",
            linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.key.height = unit(6, "pt"),
        legend.key.width = unit(4, "pt"),
        legend.key = element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 5))

p_a
# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 400
fname   = 'Fig1/Panel_A.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=400, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# b. pca ----
## 1. env & data ----

f_id_full  = glue('{list_dirs$prj_dir}2-pca/add_MCZ/list.pro.afr_mcz')
f_cov_full = glue('{list_dirs$prj_dir}2-pca/add_MCZ/pca.pro.proafr_mcz.cov')

df_id_full = read_tsv(f_id_full, col_names = c("id"))
v_id_full  = df_id_full$id

v_id_pro = c()

df_cov_full = read_table(f_cov_full, col_names = FALSE)

res_pro = pro_pcangsd(
    v_id_full,
    v_id_full,
    # v_id_se,
    v_id_pro,
    df_cov_full,
    v_pc = c(1:4)
)

## 2. plot ----
res_pro$vec_var_explained / sum(res_pro$vec_var_explained)

p_b = res_pro$df_evec %>%
# res %>%
    mutate(spp=meta_sample[match(id, meta_sample$Plot_ID),]$Plot_pop) %>%
    
    ggplot() + 
    geom_point(aes(x=PC1, y=PC3, fill=spp, shape=spp),
               size= 1, stroke=0.2,
               ) + 
    
    scale_shape_manual(breaks = pop_legend$Plot_pop,
                       values = pop_legend$Plot_region_shape) +
    scale_fill_manual(breaks = pop_legend$Plot_pop,
                      values = pop_legend$Plot_region_color) +
    
    labs(fill="Populations",
         shape="Populations",
         x="PC1 (variation 6.27%)",
         y="PC3 (variation 2.98%)",
    ) +
    
    # guides(shape = guide_legend(ncol=1, override.aes = list(size = 5/ .pt))) +
    
    theme(panel.background = element_rect(
        fill = 'white', colour = "black",
        linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = "none",
        # legend.text = element_text(size = 7),
        # legend.title = element_text(size = 7),
        # legend.key.height = unit(10, "pt"),
        # legend.key.width = unit(7, "pt"),
        # legend.key = element_blank(),
        # legend.background = element_blank(),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p_b

# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig1/Panel_B.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# c. phylogeny v1 ----
library(ggtree)
library(treeio)
# tree_file = './4-phylo/wg_nj/se.maf001.withbp.raxml.support.root.tre'
tree_file = './4-phylo/wg_nj/se.maf001.withbp.raxml.support.root.copy.tre'
# tree_raw = treeio::read.nexus(glue('{list_dirs$prj_dir}{tree_file}'))
tree_raw = ggtree::read.tree(glue('{list_dirs$prj_dir}{tree_file}'))

tree_root = root(
    tree_raw, outgroup = "MFN_MAM_056095")

# change tips
tree_add = treeio::as_tibble(tree_root) %>% 
    mutate(
        bp=c(rep(NA,60) ,tree_root$node.label),
        bp=as.numeric(bp),
        bootstrap=case_when(
            bp > 90 ~ "> 0.90",
            bp < 50 ~ "< 0.50",
            TRUE ~ "0.5~0.9")
    ) %>% 
    select(c("node", "bootstrap"))
# view()
tree_merge = full_join(tree_root, tree_add, by="node")


p_barsize     = 0.2
p_bar_offset1 = 0.011
p_bar_offset2 = 0.02
p_text_offset = 0.002
p_font_size   = 4 /.pt

p_c = tree_merge %>% 
    ggplot(aes(x, y)) +
    geom_tree(size = .1,)
p_c = p_c %>% 
    flip(68,69)

p_c = p_c + 
    geom_tiplab(
        mapping = aes(
            label = paste0(label),
            color = label,
            # label = paste0(node, "--",label),
        ),
        size= 2/ .pt) +
    # geom_nodelab(mapping = aes(label = node)) +

    scale_color_manual(
        breaks = meta_sample$Plot_ID,
        values = meta_sample$Plot_region_color) +
    
    coord_cartesian(clip = 'off') + 
    theme_tree(
        plot.margin=margin(t=3, r=10, l=3, b=3),
        legend.position = c(0.1, 0.65),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill=NA),
        legend.text = element_text(size= 7 / .pt),
        legend.title = element_text(size= 7/ .pt),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(2, "pt")
        ) + 
    
    geom_nodepoint(aes(fill=bootstrap),
                   shape=21,
                   size= 1 / .pt,
                   stroke=0.1) +
    scale_fill_manual(breaks = c("> 0.90", "0.5~0.9", "< 0.50"),
                      values = c("#636363", "#bdbdbd", "#f0f0f0")) +
    
    guides(
        fill = guide_legend(
            # override.aes = list(size = 1),
            title="Bootstrap support",
        ),
        color = "none",
    ) +
    
    geom_strip(
        53, 58, label="Namibia",
        angle = 270,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        hjust = 0.5,
        align=TRUE,
        offset = p_bar_offset1) + 
    
    geom_strip(
        12, 6, label="Tanzania",
        angle = 270,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        hjust = 0.5,
        align=TRUE,
        offset = p_bar_offset1) + 
    
    geom_strip(
        39, 60, label="Zambia",
        angle = 270,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        hjust = 0.5,
        offset = p_bar_offset1) + 
    
    geom_strip(
        38, 19, label="NE",
        angle = 270,
        hjust = 0.5,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset1) + 
    
    geom_strip(
        37, 34, label="NW",
        angle = 270,
        hjust = 0.5,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset1) +
    
    geom_strip(
        32, 33, label="Zanzibar",
        angle = 270,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset1) +

    geom_strip(
        26, 30, label="Asia",
        angle = 270,
        hjust = 0.5,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset2) +
    
    geom_strip(
        53, 31, label="Africa",
        angle = 270,
        hjust = 0.5,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset2) 

p_c

# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig1/Panel_C.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# c. phylogeny v2----
library(ggtree)
library(treeio)
# tree_file = './4-phylo/wg_nj/se.maf001.withbp.raxml.support.root.tre'
tree_file = './4-phylo/wg_nj/se.maf001.withbp.raxml.support.root.copy.tre'
# tree_raw = treeio::read.nexus(glue('{list_dirs$prj_dir}{tree_file}'))
tree_raw = ggtree::read.tree(glue('{list_dirs$prj_dir}{tree_file}'))

tree_root = root(
    tree_raw, outgroup = "MFN_MAM_056095")

# change tips
tree_add = treeio::as_tibble(tree_root) %>% 
    mutate(
        bp=c(rep(NA,60) ,tree_root$node.label),
        bp=as.numeric(bp),
        bootstrap=case_when(
            bp > 90 ~ "> 0.90",
            bp < 50 ~ "< 0.50",
            TRUE ~ "0.5~0.9")
    ) %>% 
    select(c("node", "bootstrap"))
# view()
tree_merge = full_join(tree_root, tree_add, by="node")


p_barsize     = 0.2
p_bar_offset_adjust = 0.025
p_text_offset_adjust = 0.0025
p_bar_offset1 = 0.011 + p_bar_offset_adjust
p_bar_offset2 = 0.02 + p_bar_offset_adjust
p_text_offset = 0.002 + p_text_offset_adjust
p_font_size   = 6 /.pt

p_c = tree_merge %>% 
    ggplot(aes(x, y)) +
    geom_tree(size = .1,)
p_c = p_c %>% 
    flip(68,69)

p_c = p_c + 
    geom_tiplab(
        mapping = aes(
            label = paste0(label),
            color = label,
            # label = paste0(node, "--",label),
        ),
        size= 5/ .pt) +
    # geom_nodelab(mapping = aes(label = node)) +
    
    scale_color_manual(
        breaks = meta_sample$Plot_ID,
        values = meta_sample$Plot_region_color) +
    
    coord_cartesian(clip = 'off') + 
    theme_tree(
        plot.margin=margin(t=3, r=10, l=3, b=3),
        legend.position = c(0.1, 0.65),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill=NA),
        legend.text = element_text(size= 14 / .pt),
        legend.title = element_text(size= 14/ .pt),
        legend.key.height = unit(7, "pt"),
        legend.key.width = unit(7, "pt")
    ) + 
    
    geom_nodepoint(aes(fill=bootstrap),
                   shape=21,
                   size= 1 / .pt,
                   stroke=0.1) +
    scale_fill_manual(breaks = c("> 0.90", "0.5~0.9", "< 0.50"),
                      values = c("#636363", "#bdbdbd", "#f0f0f0")) +
    
    guides(
        fill = guide_legend(
            # override.aes = list(size = 1),
            title="Bootstrap support",
        ),
        color = "none",
    ) +
    
    geom_strip(
        53, 58, label="Namibia",
        angle = 270,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        hjust = 0.5,
        align=TRUE,
        offset = p_bar_offset1) + 
    
    geom_strip(
        12, 6, label="Tanzania",
        angle = 270,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        hjust = 0.5,
        align=TRUE,
        offset = p_bar_offset1) + 
    
    geom_strip(
        39, 60, label="Zambia",
        angle = 270,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        hjust = 0.5,
        offset = p_bar_offset1) + 
    
    geom_strip(
        38, 19, label="NE",
        angle = 270,
        hjust = 0.5,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset1) + 
    
    geom_strip(
        37, 34, label="NW",
        angle = 270,
        hjust = 0.5,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset1) +
    
    geom_strip(
        32, 33, label="Zanzibar",
        angle = 270,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset1) +
    
    geom_strip(
        26, 30, label="Asia",
        angle = 270,
        hjust = 0.5,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset2) +
    
    geom_strip(
        53, 31, label="Africa",
        angle = 270,
        hjust = 0.5,
        offset.text = p_text_offset,
        barsize = p_barsize,
        fontsize= p_font_size,
        align=TRUE,
        offset = p_bar_offset2) 

p_c

# save plots
pwidth  = 80 # mm
pheight = 100 # mm
pdpi    = 600
fname   = 'Fig1/Panel_C.r.v2'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))



# d. ngsadm ----
## 1.data ---- 

f_id_full = glue('{list_dirs$prj_dir}3-ngsadm/add_MCZ/list.id')
f_id_se   = glue('{list_dirs$prj_dir}3-ngsadm/add_MCZ/list.id.qc.rel')

df_id_full = read_tsv(f_id_full, col_names = c("id"))
df_id_se   = read_tsv(f_id_se, col_names = c("id"))

v_id_se = df_id_full$id[df_id_full$id %in% df_id_se$id]
v_pop_se = meta_sample[match(v_id_se, meta_sample$Plot_ID),]$Plot_pop


# # quick
# K = 3
# f_Q = glue('{list_dirs$prj_dir}3-ngsadm/add_MCZ/se.thin.{K}.Q')
# color_panel = c(
#     '#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
#     '#ff7f00', '#a65628', '#f781bf', '#999999')
# 
# p_k = qplot_admixture(
#     f_Q       = f_Q,
#     K         = K,
#     v_id      = v_id_se,
#     v_pop     = v_pop_se,
#     v_color   = color_panel,
#     start_col = 1
# )
# 
# p_k

order_pop = c(
    "Namibia", 
    "Equatorial Guinea", "South Africa", "Zambia", "Malawi",
    "TanzaniaN", "TanzaniaS", "Tanzania Pemba", 
    "Zanzibar",
    "NE.Kenya", "NE.Uganda", "NE.Ethiopia", "NE.Eritrea", 
    "NW.Cameroon",  "NW.Gabon", "NW.Ghana","NW.Morocco",
    "Persian leopard", "Indochinese leopard",
    "Indian leopard", "Sri Lankan leopard",
    "Javan leopard", "Amur leopard")

color_panel = c(
    '#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
    '#ff7f00', '#a65628', '#f781bf', '#999999')


# batch
df_fQ_K = data.frame(
    K=c(2:5)) %>% 
    mutate(
        file=glue('{list_dirs$prj_dir}/3-ngsadm/add_MCZ/1706691322/aligned.files/se.thin.{K}.Q.converted'))

l_res_p_adm = plot_admixture(
    df_fQ_K      = df_fQ_K,
    v_id         = v_id_se,
    v_pop        = v_pop_se,
    K_order      = 4,
    v_pop_order  = order_pop,
    v_color      = color_panel,
    start_col    = 6
)

## 2.plot v1----
p_d = l_res_p_adm$df %>% 
    mutate(
        id=factor(id, levels = l_res_p_adm$id_order),
        pop=factor(pop, levels = order_pop),
        K=glue('K={K}')
    ) %>%
    # view()
    
    ggplot() + 
    geom_bar(aes(
        x=id,
        y=percent,
        fill=factor(group, levels=1:20)),  # change number of color
        width = 1,  # no gap in plot
        stat = 'identity') +
    facet_grid(
        rows = vars(K), 
        cols = vars(pop),
        labeller = labeller(pop = label_wrap_gen(11)),
        switch = "x", scales = "free", space = "free") +
    
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expand_scale(add = 1)) +
    scale_fill_manual(values=color_panel) + 
    
    theme(
        # axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        strip.text.y = element_text(angle = 0, size = 9 / .pt, hjust = 0),
        strip.text.x = element_text(angle = 90,
                                    size = 7 / .pt,
                                    hjust = 1,
                                    vjust = 0.5,
                                    ),
        strip.clip = "off",
        strip.background = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.spacing.y = unit(3, "pt"),
        panel.spacing.x = unit(0, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
p_d

# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig1/Panel_D.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))

## 2.plot v2----
p_d = l_res_p_adm$df %>% 
    mutate(
        id=factor(id, levels = l_res_p_adm$id_order),
        pop=factor(pop, levels = order_pop),
        K=glue('K={K}')
    ) %>%
    # view()
    
    ggplot() + 
    geom_bar(aes(
        y=id,
        x=percent,
        fill=factor(group, levels=1:20)),  # change number of color
        width = 1,  # no gap in plot
        stat = 'identity') +
    facet_grid(
        cols = vars(K), 
        rows = vars(pop),
        # labeller = labeller(pop = label_wrap_gen(11)),
        switch = "x", scales = "free", space = "free") +
    
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expand_scale(add = 1)) +
    scale_fill_manual(values=color_panel) + 

    theme(
        # axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        strip.text.x = element_text(size = 20 / .pt,
                                    vjust = 0.5),
        strip.text.y = element_text(size = 18 / .pt,
                                    angle = 0,
                                    hjust = 0,
                                    vjust = 0.5,
        ),
        strip.clip = "off",
        strip.background = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.spacing.x = unit(3, "pt"),
        panel.spacing.y = unit(0, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
p_d

# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig1/Panel_D.r.v2'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# Final plot v1----
library(cowplot)

p_abcd = plot_grid(
    p_a, 
    p_c,
    p_b,
    p_d,
    # align = "vh",
    ncol=2,
    label_size = 12 /.pt,
    labels = c('A', 'C', 'B', 'D')
)

p_abcd

# save plots
pwidth  = 160 # mm
pheight = 120 # mm
pdpi    = 600
fname   = 'Fig1/Fig1.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))



# Final plot v2----
library(cowplot)

p_abcd = plot_grid(
    p_a, 
    p_b,
    p_c,
    p_d,
    # align = "vh",
    rel_heights = c(1,1.5),
    ncol=2,
    label_size = 12 /.pt,
    labels = c('A', 'C', 'B', 'D')
)

p_abcd

# save plots
pwidth  = 160 # mm
pheight = 180 # mm
pdpi    = 600
fname   = 'Fig1/Fig1.r.v2'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


