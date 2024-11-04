# Figure 2
# a. f3
# b. hPSMC
# c. qpAdm

# 0. env all ----
library(tidyverse)
library(ggrepel)
library(glue)

source('./lib_meta.R')
meta_sample = load_meta_gdoc()
list_dirs = get_dir_list()
pop_legend = get_legend(meta_sample)

source('./lib_utility.R')


# a. f3 ----
## 1. data ----
f_rdata = '8-admixtools/res.f3.mac3.nopoly.rdata' # MCZ1
f_rdata = '8-admixtools/res.f3.mac3.nopoly.mcz2.rdata' # MCZ2

load(glue('{list_dirs$prj_dir}{f_rdata}'))
res_f3_old = res_f3


## 2. plot ----
p_a = res_f3 %>%
    
    # change pop3 name
    mutate(
        id  = pop3,
        spp = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_spp,
        dp  = meta_sample[match(id, meta_sample$Plot_ID),]$DP_mem_noS,
        err = meta_sample[match(id, meta_sample$Plot_ID),]$Filter_error_rate,
        region = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_pop,
    ) %>%
    
    filter(! region %in% c("Zanzibar")) %>%
    
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
fname   = 'Fig2/Panel_A.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=pdpi, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# b. hPSMC ----
## 1. data ----
f_res = glue('{list_dirs$prj_dir}6-hpsmc/plot_hpsmc.X.all_g7.5u1e8')
# f_res = glue('{list_dirs$prj_dir}6-hpsmc/plot_hpsmc.X.all_g5u1.1e9')
# f_res = glue('{list_dirs$prj_dir}6-hpsmc/plot_hpsmc.X.N.tv.all_g7.5u1e8')
# f_res = glue('{list_dirs$prj_dir}6-hpsmc/plot_hpsmc.X.N.tv.all_g5u1.1e9')

df_res = read_tsv(
    f_res,
    col_names = c('time', 'Ne', 'type', 'rep', 'indi'),
)

## 2. plot ----
color_list = c(
    '#4daf4a', '#984ea3', '#e41a1c', '#377eb8',
    '#ff7f00', '#f781bf','#999999')

p_b = df_res %>% 
    separate(indi, c("g", "g2"), sep = "-vs-", remove = F) %>%
    mutate(g=case_when(g=="4343_Tanzania" ~ "Africa-Asia",
                       g=="MFN_MAM_056356" ~ "Africa-Africa",
                       g=="MCZ36709" &
                           g2 %in% c("4343_Tanzania", "L033-L0665a", "MFN_MAM_056356") ~ "Zanzibar-Africa" ,
                       g=="MCZ36709" &
                           g2 %in% c("Amurleopard_PPO1", "Bhagya", "MFN_MAM_050746") ~ "Zanzibar-Asia")) %>%
    
    ggplot() + 
    geom_step(
        aes(x=time, y=Ne, group=factor(indi), color=factor(g)),
        size= 0.25) +
    geom_vline(xintercept = 2e5, linetype=2, linewidth=0.25) +
    geom_hline(yintercept = 3, linetype=2, linewidth=0.25) +
    
    scale_x_log10(
        # breaks = scales::trans_breaks(
        #     "log10", function(x) 10^x, n=2)(c(1e4,1e6)),
        # breaks = c(1e4,1e5,1e6),
        breaks = c(1e5,1e6),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + 
    
    scale_color_manual(values = color_list) +
    annotation_logticks(sides = "b", size=0.25, linewidth=0.25) + 
    theme(
        panel.background = element_rect(
            fill = NA, colour = "black",
            linetype = 1, linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.major = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        axis.text = element_text(size=5),
        axis.title=element_text(size=7),
        legend.text=element_text(size=7),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.key = element_blank(),
        legend.position = c(0.7,0.7),
        legend.background = element_rect(fill=NA),
        
        legend.title = element_blank()
    ) +
    
    ylab(expression("Effective population size (x10"^4*")")) +
    xlab(expression("g=7.5, mutation rate=1x10"^-8)) +
    # xlab(expression("g=5, mutation rate=1.1x10"^-9)) +
    coord_cartesian(ylim=c(0, 40), xlim=c(9e4,4e6))

p_b

# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig2/Panel_B.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=pdpi, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# c. qpAdm ----
## 1. data ----
f_res_qpadm = glue('{list_dirs$prj_dir}/8-admixtools/res.qpadm.rdata')
load(f_res_qpadm)

pop_legend_qpadm = pop_legend %>% 
    mutate(
        pop_adm = case_when(
            Plot_pop == "TanzaniaS" ~ "Tanzania",
            TRUE ~ Plot_pop
        )
    )
                   
## 2. plot ----
p_c = best_df_all %>% 
    filter(id %in% c("MCZ36709")) %>% 
    
    mutate(
        source = case_when(
            source == "SMNH582373" ~ "NW.Morocco",
            source == "SMNH595313" ~ "NE.Eritrea",
            source == "X4343_Tanzania" ~ "Tanzania",
            source == "ZMUC4446" ~ "Zambia",
            source == "X7942_Namibia" ~ "Namibia",
            source == "Bhagya" ~ "Indian leopard",
            TRUE ~ source
        )
    ) %>% 
    
    filter(pat == "00") %>% 
    
    # sort ids
    arrange(id, source, desc(value)) %>% 
    mutate(id = factor(id, levels=unique(id))) %>%  
    
    ggplot() + 
    geom_bar(aes(x=midx,
                 y=value,
                 fill=factor(source)),  # change number of color
             width = 1,  # no gap in plot
             stat = 'identity') +
    facet_grid(
        cols = vars(id),
        # rows = vars(pop),
        switch = "x", scales = "free", space = "free") +
    
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expand_scale(add = 0)) +
    scale_fill_manual(
        breaks = pop_legend_qpadm$pop_adm,
        values = pop_legend$Plot_region_color) +
    
    theme(
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 5),
        strip.text.x = element_text(angle = 0, size = 7),
        strip.clip = "off",
        strip.background = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.height = unit(7, "pt"),
        # panel.spacing.y = unit(13, "pt"),
        panel.spacing.x = unit(3, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    guides(
        fill = guide_legend(
            title=str_wrap("qpAdm source populations", width = 15),
        ),
    )

p_c

# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig2/Panel_C.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=pdpi, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))



# Final plot ----
library(cowplot)

p_bc = plot_grid(
    p_b,
    p_c,
    ncol=1,
    label_size = 12 /.pt,
    labels = c('B', 'C')
)
p_bc

p_abc = plot_grid(
    p_a, 
    p_bc,
    ncol=2,
    label_size = 12 /.pt,
    labels = c('A', '')
)

p_abc

# save plots
pwidth  = 160 # mm
pheight = 120 # mm
pdpi    = 600
fname   = 'Fig2/Fig2.r.v2'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))





