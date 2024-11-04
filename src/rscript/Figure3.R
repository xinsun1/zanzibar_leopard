# Figure 3
# a. point plot, het, roh, gload
# b. het comparision
# c. roh distribution

# 0. env all ----
library(tidyverse)
library(ggrepel)
library(glue)

source('./lib_meta.R')
meta_sample = load_meta_gdoc()
list_dirs = get_dir_list()
pop_legend = get_legend(meta_sample)

source('./lib_utility.R')


# a. point plot ----
## 1. data ----
f_roh  = glue('{list_dirs$prj_dir}5-het_roh/dp5.het3.maf05.mis50.tv.mis20.hom.indiv')
f_roh_all = glue('{list_dirs$prj_dir}5-het_roh/dp5.het3.maf05.mis50.tv.mis20.hom')

f_het  = glue('{list_dirs$prj_dir}5-het_roh/sample_het_tv_noS')

f_load_pha = glue('{list_dirs$prj_dir}9-load/gload.phastcon100.t0.75')
f_load_phy = glue('{list_dirs$prj_dir}9-load/gload.phylop100.t1.5')
f_load_snf = glue('{list_dirs$prj_dir}9-load/ann.sum.load')
f_load_snp = glue('{list_dirs$prj_dir}9-load/gload.pha.phy.totalsnp')

# roh
df_roh = read_table(f_roh)
len_wg_noXMT=2281341621

m_roh = df_roh %>%
    group_by(IID) %>%
    filter(KB>1000) %>%
    summarise(
        nroh_1mb = n(),
        sroh = sum(KB)/1e3,
        froh = sum(KB)*1e3/len_wg_noXMT) %>% 
    
    mutate(
        spp = meta_sample[match(IID, meta_sample$Plot_ID),]$Plot_spp,
        dp  = meta_sample[match(IID, meta_sample$Plot_ID),]$DP_mem_noS,
        err = meta_sample[match(IID, meta_sample$Plot_ID),]$Filter_error_rate,
        region = meta_sample[match(IID, meta_sample$Plot_ID),]$Plot_region) %>%
    
    filter(dp > 10) %>%
    
    mutate(plot_spp=case_when(
        region == "Zanzibar" ~ region,
        TRUE ~ spp)) %>% 
    
    mutate(
        type  = "Inbreeding coefficent",
        value = froh,
        pop   = plot_spp,
        id    = IID,
    ) %>% 
    select(type, value, pop, id)


m_roh %>% 
    ggplot() +
    geom_boxplot(aes(x=pop, y=value, color=pop),
                 outlier.shape = NA) +
    geom_point(aes(x=pop, y=value))


# het
df_het = read_table(f_het,
                    col_names = c('id', 'het'))


m_het = df_het %>%
    mutate(
        spp = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_spp,
        dp  = meta_sample[match(id, meta_sample$Plot_ID),]$DP_mem_noS,
        err = meta_sample[match(id, meta_sample$Plot_ID),]$Filter_error_rate,
        region = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_region) %>%
    # view()
    filter(dp > 10) %>%
    filter(err %in% c("No")) %>%
        
    mutate(plot_spp=case_when(
        region == "Zanzibar" ~ region,
        TRUE ~ spp)) %>% 
    
    mutate(
        type  = "Heterozygosity",
        value = het,
        pop   = plot_spp,
        id    = id,
    ) %>% 
    select(type, value, pop, id)

m_het %>% 
    ggplot() +
    geom_boxplot(aes(x=pop, y=value, color=pop),
                 outlier.shape = NA) +
    geom_point(aes(x=pop, y=value))

# load phy
df_load_snp = read_table(
    f_load_snp,
    col_names = c("id", "tSNP")
)
df_load_pha = read_table(
    f_load_pha,
    col_names = c("id", "NSNP", "load", "r_load"))

df_load_phy = read_table(
    f_load_phy,
    col_names = c("id", "NSNP", "load", "r_load"))

df_load_snf = read_table(
    f_load_snf,
    col_names = c("id",  "lof/syn", "non/syn"))


m_load_pha = df_load_pha %>% 
    mutate(
        tSNP=df_load_snp[match(id, df_load_snp$id), ]$tSNP,
        load=load/tSNP,
    ) %>% 
    
    mutate(
        spp = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_spp,
        dp  = meta_sample[match(id, meta_sample$Plot_ID),]$DP_mem_noS,
        err = meta_sample[match(id, meta_sample$Plot_ID),]$Filter_error_rate,
        region = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_region) %>%
    # view()
    filter(dp > 10) %>%

    mutate(plot_spp=case_when(
        region == "Zanzibar" ~ region,
        TRUE ~ spp)) %>% 
    
    mutate(
        type  = "load_pha",
        value = load,
        pop   = plot_spp,
        id    = id,
    ) %>% 
    select(type, value, pop, id)

m_load_phy = df_load_phy %>% 
    mutate(
        tSNP=df_load_snp[match(id, df_load_snp$id), ]$tSNP,
        load=load/tSNP,
    ) %>% 
    
    mutate(
        spp = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_spp,
        dp  = meta_sample[match(id, meta_sample$Plot_ID),]$DP_mem_noS,
        err = meta_sample[match(id, meta_sample$Plot_ID),]$Filter_error_rate,
        region = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_region) %>%
    # view()
    filter(dp > 10) %>%
    
    mutate(plot_spp=case_when(
        region == "Zanzibar" ~ region,
        TRUE ~ spp)) %>% 
    
    mutate(
        type  = "load_phy",
        value = load,
        pop   = plot_spp,
        id    = id,
    ) %>% 
    select(type, value, pop, id)

m_load_snf = df_load_snf %>% 
    pivot_longer(
        cols = c("lof/syn", "non/syn"),
        names_to = "type",
        values_to = "value") %>%
    
    mutate(
        spp = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_spp,
        dp  = meta_sample[match(id, meta_sample$Plot_ID),]$DP_mem_noS,
        err = meta_sample[match(id, meta_sample$Plot_ID),]$Filter_error_rate,
        region = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_region) %>%
    filter(dp > 10) %>%
    
    mutate(plot_spp=case_when(
        region == "Zanzibar" ~ region,
        TRUE ~ spp)) %>% 
    
    mutate(
        pop   = plot_spp,
    ) %>% 
    select(type, value, pop, id)

m_load = rbind(m_load_pha, m_load_phy, m_load_snf)

# merge
m_all = rbind(m_load, m_het, m_roh)


## 2. plot ----
order_pop = c(
    "Zanzibar",
    "African leopard",
    "Persian leopard", "Indochinese leopard",
    "Indian leopard", "Sri Lankan leopard",
    "Javan leopard", "Amur leopard")

p_a = m_all %>%
    filter(
        type %in% c(
            "Heterozygosity",
            "Inbreeding coefficent",
            "load_phy", 
            "lof/syn"),
    ) %>% 
    
    mutate(
        type = case_when(
            type == "load_phy" ~ 'Genetic load (PhyloP score)',
            type == "lof/syn" ~ "Genetic load (Lof/synonoymous)",
            type == "Heterozygosity" ~ "Heterozygosity (transversions)",
            type == "Inbreeding coefficent" ~  "Inbreeding coefficent",
            TRUE ~ type
        )
    ) %>% 
    
    mutate(
        pop = str_wrap(pop, width=10),
        pop = factor(pop, levels=str_wrap(order_pop, width=10)),
    ) %>% 

    ggplot() +
    geom_boxplot(
        aes(x=pop,y=value),
        outlier.shape = NA,
        fill='grey', linewidth=0.2) +
    
    geom_jitter(
        aes(x=pop,y=value),
        width = 0.1,
        size= 0.5, stroke=0.1,) +
    
    facet_wrap(vars(type),
               ncol=2, scales = "free_y") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,
             linewidth=0.1)+
    
    theme_classic() +
    theme(
        line = element_line(linewidth = 0.1),
        aspect.ratio = 9/16,
        axis.line=element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=16/ .pt),
        axis.text.x = element_text(
            angle = 90, hjust=1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
p_a


# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig3/Panel_A.r.v2'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))

# b. roh: cumulative ----
## 1. load data ----
f_roh_all = glue('{list_dirs$prj_dir}5-het_roh/dp5.het3.maf05.mis50.tv.mis20.hom')
# f_roh_all = glue('{list_dirs$prj_dir}5-het_roh/test.mis10.hom')
df_roh_all = read_table(f_roh_all)

## 2. plot ----
raw_meta = df_roh_all %>%
    
    mutate(id=IID) %>%
    mutate(
        spp = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_spp,
        dp  = meta_sample[match(id, meta_sample$Plot_ID),]$DP_mem_noS,
        err = meta_sample[match(id, meta_sample$Plot_ID),]$Filter_error_rate,
        region = meta_sample[match(id, meta_sample$Plot_ID),]$Plot_pop) %>% 
    
    mutate(pop = case_when(
        spp == "African leopard" ~ region,
        TRUE ~ spp,
    ))


len_wg_noXMT=2281341621
# !!! changed for rev1
# x: max_roh_len, y: csum_roh
list_len = seq(0.5e6, 6e6, length.out=100) # 500kb, 30 mb
froh_cumlen = c()
for(l in list_len){
    tmp_froh = raw_meta %>% 
        
        filter(KB <= (l/1e3)) %>% 
        
        group_by(pop) %>%
        summarise(
            n = length(unique(id)),
            s_roh = sum(KB),
            f_roh = s_roh*1e3/len_wg_noXMT/n,
            pop = unique(pop)
        ) %>%
        mutate(max_len = l)
    froh_cumlen = rbind(froh_cumlen, tmp_froh)
}

# reverse length order for rev1
list_len = seq(0.5e6, 6e6, length.out=100) # 500kb, 30 mb
list_len = rev(list_len)
froh_cumlen = c()
for(l in list_len){
    tmp_froh = raw_meta %>% 
        
        filter(KB > (l/1e3)) %>% 
        
        group_by(pop) %>%
        summarise(
            n = length(unique(id)),
            s_roh = sum(KB),
            f_roh = s_roh*1e3/len_wg_noXMT/n,
            pop = unique(pop)
        ) %>%
        mutate(max_len = l)
    froh_cumlen = rbind(froh_cumlen, tmp_froh)
}



froh_cumlen_label = froh_cumlen %>% 
    filter(max_len == max(list_len)) %>% 
    mutate(name=glue('{pop} ({n})'))

p_b = froh_cumlen %>%
    ggplot() +
    
    geom_line(
        aes(x=max_len,
            y=f_roh,
            group = pop,
            color = pop
        ),
        linewidth=0.5
    ) +
    scale_color_manual(
        breaks = pop_legend$Plot_pop,
        values = pop_legend$Plot_region_color
    ) +
    xlim(0.5e6, 2e7) +

    scale_x_continuous(trans = c("log10", "reverse")) +
    # scale_x_log10(
    #     breaks = c(1e6,1e7),
    #     labels = scales::trans_format("log10", scales::math_format(10^.x))
    # ) +

    annotation_logticks(sides = "b",
                        size=0.25,
                        short = unit(.5,"mm"),
                        mid = unit(1,"mm"),
                        long = unit(1.5,"mm")) +
    geom_text_repel(
        data=froh_cumlen_label,
        aes(x=max_len,
            y=f_roh, label=name),
        size = 3 / .pt,
        force  = 2,
        nudge_x = 0.3,
        direction = "y",
        hjust = 0,
        max.overlaps = 3, 
        max.iter = 1e6,
        # seed = 1234,
        # segment.curvature = 0.1,
        min.segment.length = 0,
        segment.size = 0.1
    ) +
    theme(panel.background = element_rect(
        fill = 'white', colour = "black",
        linetype = "solid", linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.position = "none",
        aspect.ratio = 1,
        axis.text = element_text(size=5),
        axis.title = element_text(size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    
    xlab("Maximum ROH length (bp)") +
    ylab("Proportion of genome in ROH")
p_b


# save plots
pwidth  = 80 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig3/Panel_B.r.v2'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))


# Final plot v1 ----
library(cowplot)

p_ab = plot_grid(
    p_a,
    p_b,
    nrow=1,
    rel_widths = c(2, 2),
    label_size = 12 /.pt,
    labels = c('A', 'B')
)
p_ab

# save plots
pwidth  = 150 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig3/Fig3.r.v1'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))

    
# Final plot v2 ----
library(cowplot)

p_a
# save plots
pwidth  = 120 # mm
pheight = 60 # mm
pdpi    = 600
fname   = 'Fig3/Fig3.r.v2'
save_plots(pwidth=pwidth,
           pheight=pheight, # mm
           pdpi=600, 
           fname_prefix=glue('{list_dirs$plot_dir}{fname}'))    
    
    
    
