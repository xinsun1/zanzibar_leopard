# supplementary figures


# 0. env_all ----
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


# pca: all ----
## 1. env & data ----
f_id_full  = glue('{list_dirs$prj_dir}2-pca/add_MCZ/list.id')

# f_id_se    = glue('{list_dirs$prj_dir}2-pca/add_MCZ/list.pro.afr')
# f_cov_full = glue('{list_dirs$prj_dir}2-pca/add_MCZ/pca.pro.proafr.pro_sites.cov')

f_id_se    = glue('{list_dirs$prj_dir}2-pca/add_MCZ/list.pro.all')
f_cov_full = glue('{list_dirs$prj_dir}2-pca/add_MCZ/pca.pro.proall.pro_sites.cov')

# f_id_se    = glue('{list_dirs$prj_dir}2-pca/add_MCZ/list.pro.afr_mcz')
# f_cov_full = glue('{list_dirs$prj_dir}2-pca/add_MCZ/pca.pro.proafr_mcz.pro_sites.cov')



df_id_full = read_tsv(f_id_full, col_names = c("id"))
v_id_full  = df_id_full$id
df_id_se   = read_tsv(f_id_se, col_names = c("id"))
v_id_se    = df_id_se$id

# v_id_pro = c()
# v_id_pro = c("MCZ40953")

v_id_se = v_id_se[! v_id_se %in% c("MCZ36709")]
v_id_pro = c("MCZ40953", "MCZ36709")


df_cov_full = read_table(f_cov_full, col_names = FALSE)

res_pro = pro_pcangsd(
    v_id_full,
    # v_id_full,
    v_id_se,
    v_id_pro,
    df_cov_full,
    v_pc = c(1:4)
)


## 2. plot ----
res_pro$vec_var_explained / sum(res_pro$vec_var_explained)

res_pro$df_evec %>%
    # res %>%
    mutate(spp=meta_sample[match(id, meta_sample$Plot_ID),]$Plot_pop) %>%
    
    ggplot() + 
    geom_point(aes(x=PC1, y=PC2, fill=spp, shape=spp),
               size= 2, stroke=0.2,
    ) + 
    
    scale_shape_manual(breaks = pop_legend$Plot_pop,
                       values = pop_legend$Plot_region_shape) +
    scale_fill_manual(breaks = pop_legend$Plot_pop,
                      values = pop_legend$Plot_region_color) +
    
    labs(fill="Populations",
         shape="Populations",
         x="PC1 (variation 12.2%)",
         y="PC2 (variation 3.79%)",
    ) +
    
    guides(shape = guide_legend(ncol=1, override.aes = list(size = 5/ .pt))) +
    
    theme(panel.background = element_rect(
        fill = 'white', colour = "black", linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.height = unit(7, "pt"),
        legend.key.width = unit(7, "pt"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.title = element_text(size = 9),
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plots
pwidth  = 120 # mm
pheight = 80 # mm
pdpi    = 600
fname   = 'Spp/PCA_All.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# pca: afr 1-2 ----
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

res_pro$df_evec %>%
    # res %>%
    mutate(spp=meta_sample[match(id, meta_sample$Plot_ID),]$Plot_pop) %>%
    
    ggplot() + 
    geom_point(aes(x=PC1, y=PC2, fill=spp, shape=spp),
               size= 2, stroke=0.2,
    ) + 
    
    scale_shape_manual(breaks = pop_legend$Plot_pop,
                       values = pop_legend$Plot_region_shape) +
    scale_fill_manual(breaks = pop_legend$Plot_pop,
                      values = pop_legend$Plot_region_color) +
    
    labs(fill="Populations",
         shape="Populations",
         x="PC1 (variation 6.27%)",
         y="PC2 (variation 3.66%)",
    ) +
    
    guides(shape = guide_legend(ncol=1, override.aes = list(size = 5/ .pt))) +
    
    theme(panel.background = element_rect(
        fill = 'white', colour = "black",
        linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(7, "pt"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plots
pwidth  = 120 # mm
pheight = 80 # mm
pdpi    = 600
fname   = 'Spp/PCA_afr_PC12.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# IBD ----
## 1. data ----
# localNGSrealate

f_ibd_pair = glue('{list_dirs$prj_dir}10-ibd/list_to_compare')
df_ibd_pair = read_table(f_ibd_pair, col_names = c('p1', 'p2'))

df_ibd_pair = df_ibd_pair %>%
    # add region information
    mutate(
        sp1=meta_sample[match(p1, meta_sample$Plot_ID),]$Plot_region,
        sp2=meta_sample[match(p2, meta_sample$Plot_ID),]$Plot_region,
        ) %>%
    
    # generate name
    # MCZ36709--ZMUC3490.IBDtractinference.gz
    mutate(
        in_file=glue('{list_dirs$prj_dir}10-ibd/{p1}--{p2}.IBDtractinference.gz'),
        pop=glue("{sp1}--{sp2}"),
        )

df_res_all = c()
for(i in 1:nrow(df_ibd_pair)){
    df_tmp_res = read_table(df_ibd_pair[i,]$in_file)
    df_tmp_res = df_tmp_res %>% 
        mutate(
            g = df_ibd_pair[i,]$pop,
        )
    df_res_all = rbind(df_res_all, df_tmp_res)
}

## 2. plot ----

df_res_all %>%
    
    filter(
        ! g %in% c(
            "Zanzibar--Malawi"
        )) %>%
    
    mutate(
        Chr=substr(Chr, 4,5)
    ) %>% 
    # view()
    
    ggplot() +
    geom_point(
        aes(x=Pos, y=Viterbi, color=factor(Viterbi)),
        shape=20, size=0.5,
        ) +
    facet_grid(
        rows = vars(g),
        cols = vars(Chr),
        switch = "x",
        space = "free_x",
        scales = "free_x") + 
    scale_y_continuous(breaks = c(0,1,2)) +

    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,
             linewidth=0.25)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,
             linewidth=0.25)+
    
    theme(
        panel.background = element_blank(),
        line = element_line(linewidth = 0.25),
        # axis.line=element_line(),
        # aspect.ratio = 16/4,
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        strip.clip = "off",
        strip.background = element_blank(),
        strip.text.x = element_text(size=5),
        strip.text.y = element_text(size=6, angle=0, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=7),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5)
        ) +
    labs(color="IBD sharing",
         y="IBD-Viterbi",
    )

# save plots
pwidth  = 120 # mm
pheight = 90 # mm
pdpi    = 400
fname   = 'Spp/IBD.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=pdpi, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))



# admK ----
## 1. data ----
f_ld_full = glue('{list_dirs$prj_dir}3-ngsadm/add_MCZ/best_likelihood.all.se.thin')
df_ld_full = read_delim(f_ld_full, delim=" ",
                      col_names = c("K", "rep", "ld"))

## 2. plot ----
df_ld_full %>% 
    
    # get highest ld per K
    group_by(K) %>% 
    summarise(ld=max(ld)) %>% 
    
    # view()
    
    ggplot() +
    
    geom_line(aes(x=K, y=ld)) + 
    
    labs(x="K for NGSadmix",
         y="Maximum log likelihood for each K",
    ) +
    
    theme(panel.background = element_rect(
        fill = 'white', colour = "black", linetype = "solid"),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.height = unit(7, "pt"),
        legend.key.width = unit(7, "pt"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.title = element_text(size = 9),
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plots
pwidth  = 120 # mm
pheight = 120 # mm
pdpi    = 600
fname   = 'Spp/K_adm.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))




## 2.2 plot range ----
df_ld_full %>% 
    
    # get highest ld per K
    group_by(K) %>% 
    
    ggplot() +
    
    geom_boxplot(aes(group=K, y=ld))

# f3 Pemba ----
## 1. data ----
dir_mount = '~/Downloads/test/'
f_rdata = '14-admixtools2/res.f3.mac3.nopoly.rdata' # MCZ1
f_rdata = '14-admixtools2/res.f3.mac3.nopoly.mcz2.rdata' # MCZ2
load(glue('{dir_mount}{f_rdata}'))
res_f3_old = res_f3

f_rdata = glue('{list_dirs$prj_dir}8-admixtools/res.f3.mac3.nopoly.pemba.rdata')
load(f_rdata)


## 2. plot ----
res_f3 %>%
    
    # change pop3 name
    mutate(
        id  = pop3,
        spp = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_spp,
        dp  = meta_sample[match(id, meta_sample$Plot_ID),]$DP_mem_noS,
        err = meta_sample[match(id, meta_sample$Plot_ID),]$Filter_error_rate,
        region = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_pop,
    ) %>%
    
    filter(! pop3 %in% c("MFN_MAM_056356", "MCZ40953")) %>%
    
    arrange(est) %>%
    mutate(id=factor(id, levels=id)) %>%
    
    ggplot() +
    geom_pointrange(
        aes(x=id,
            y=est, ymin=est-se, ymax=est+se,
            shape=region,
            fill=region),
        linewidth = 0.25,
        size= 0.5,
        stroke=0.2) +
    
    scale_shape_manual(breaks = pop_legend$Plot_pop,
                       values = pop_legend$Plot_region_shape) +
    scale_fill_manual(breaks = pop_legend$Plot_pop,
                      values = pop_legend$Plot_region_color) +
    
    # guides(shape = guide_legend(
    #     ncol=1, override.aes = list(size = 0.25),
    # )) +
    
    theme(
        panel.background = element_rect(
            fill = NA, colour = "black",
            linetype = 1, linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key = element_blank(),
        legend.box.spacing = unit(1, "pt"),
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 5)) +
    
    coord_flip() 

# save plots
pwidth  = 80 # mm
pheight = 120 # mm
pdpi    = 600
fname   = 'Spp/f3_pemba.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=pdpi, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


