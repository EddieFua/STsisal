rm(list = ls())
library(ggsignif)
library(ggpubr)
# load('/Users/fuyinghao/Documents/STsisal/exp_res/sim_res_8_22.RData')
# load('/Users/fuyinghao/Documents/STsisal/Revision/res/sim_res.RData')
load('/Users/fuyinghao/Documents/STsisal/exp_res/sim_res_updated.RData')
library(ggplot2)
color = RColorBrewer::brewer.pal(5,'Set3')[1:3]
##rmse
rmse =  rbind(STsisal = cbind(unlist(i_STsisal_rmse),rep('STsisal',60)),
              STdeconvolve = cbind(unlist(i_STdeconvolve_rmse),rep('STdeconvolve',60)),
              CARD = cbind(unlist(i_CARD_rmse),rep('CARD',60)),
              RCTD = cbind(unlist(i_RCTD_rmse),rep('RCTD',60)))
rmse = cbind(rmse,noise = rep(c(rep(0,20),rep(0.3,20),rep(0.6,20)),4))
colnames(rmse) = c('RMSE','Method','Heterozygosity rate')
rmse = as.data.frame(rmse)
rmse$Method = as.factor(rmse$Method)
rmse$`Heterozygosity rate` = as.factor(rmse$`Heterozygosity rate`)
rmse$RMSE = as.numeric(rmse$RMSE)
p_rmse = ggplot(rmse, aes(x=Method, y=RMSE, fill=`Heterozygosity rate`)) + 
  geom_boxplot()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_brewer(palette="BuPu")+labs(title = "RMSE")+
  theme(legend.position = "none")

##corr
corr =  rbind(STsisal = cbind(unlist(i_STsisal_corr),rep('STsisal',60)),
              STdeconvolve = cbind(unlist(i_STdeconvolve_corr),rep('STdeconvolve',60)),
              CARD = cbind(unlist(i_CARD_corr),rep('CARD',60)),
              RCTD = cbind(unlist(i_RCTD_corr),rep('RCTD',60)))
corr = cbind(corr,noise = rep(c(rep(0,20),rep(0.3,20),rep(0.6,20)),4))
colnames(corr) = c('Corr','Method','Heterozygosity rate')
corr = as.data.frame(corr)
corr$Method = as.factor(corr$Method)
corr$`Heterozygosity rate` = as.factor(corr$`Heterozygosity rate`)
corr$Corr = as.numeric(corr$Corr)
p_corr = ggplot(corr, aes(x=Method, y=Corr, fill=`Heterozygosity rate`)) + 
  geom_boxplot()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_brewer(palette="BuPu")+labs(title = "Corr")+
  theme(legend.position = "none")

##jsd
jsd =  rbind(STsisal = cbind(unlist(i_STsisal_jsd),rep('STsisal',60)),
             STdeconvolve = cbind(unlist(i_STdeconvolve_jsd),rep('STdeconvolve',60)),
             CARD = cbind(unlist(i_CARD_jsd),rep('CARD',60)),
             RCTD = cbind(unlist(i_RCTD_jsd),rep('RCTD',60)))
jsd = cbind(jsd,noise = rep(c(rep(0,20),rep(0.3,20),rep(0.6,20)),4))
colnames(jsd) = c('Jsd','Method','Heterozygosity rate')
jsd = as.data.frame(jsd)
jsd$Method = as.factor(jsd$Method)
jsd$`Heterozygosity rate` = as.factor(jsd$`Heterozygosity rate`)
jsd$Jsd = as.numeric(jsd$Jsd)
p_jsd = ggplot(jsd, aes(x=Method, y=Jsd, fill=`Heterozygosity rate`)) + 
  geom_boxplot()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_brewer(palette="BuPu")+labs(title = "JSD")+
  theme(legend.position = "none")

##MAE
MAE =  rbind(STsisal = cbind(unlist(i_STsisal_MAE),rep('STsisal',60)),
             STdeconvolve = cbind(unlist(i_STdeconvolve_MAE),rep('STdeconvolve',60)),
             CARD = cbind(unlist(i_CARD_MAE),rep('CARD',60)),
             RCTD = cbind(unlist(i_RCTD_MAE),rep('RCTD',60)))
MAE = cbind(MAE,noise = rep(c(rep(0,20),rep(0.3,20),rep(0.6,20)),4))
colnames(MAE) = c('MAE','Method','Heterozygosity rate')
MAE = as.data.frame(MAE)
MAE$Method = as.factor(MAE$Method)
MAE$`Heterozygosity rate` = as.factor(MAE$`Heterozygosity rate`)
MAE$MAE = as.numeric(MAE$MAE)
p_MAE = ggplot(MAE, aes(x=Method, y=MAE, fill=`Heterozygosity rate`)) + 
  geom_boxplot()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_brewer(palette="BuPu")+labs(title = "MAE")+
  theme(legend.position = "none")

t.test(corr$Corr[corr$Method=='STsisal'&corr$`Heterozygosity rate`==0.3],corr$Corr[corr$Method=='STdeconvolve'&corr$`Heterozygosity rate`==0.3])
library(cowplot)

# 设置图像大小和图例位置
plot_size <- 5  # 图像大小为 5 英寸
legend_position <- "top"  # 图例位置设为上方

# 调整图像大小和间距
p_rmse <- p_rmse + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
p_corr <- p_corr + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
p_jsd <- p_jsd + theme(plot.margin = margin(1, 1, 1, 1, "cm"))
p_MAE <- p_MAE + theme(plot.margin = margin(1, 1, 1, 1, "cm"))

# 组合图像

uniform_theme <- theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA)  # Add borders
  )

# Adjust individual plots
p_rmse <- p_rmse + uniform_theme
p_corr <- p_corr + uniform_theme
p_jsd <- p_jsd + uniform_theme
p_MAE <- p_MAE + uniform_theme

# Combine plots into a single row with shared legend
combined_plot <- plot_grid(p_rmse, p_corr, p_jsd, p_MAE, align = "hv", nrow = 1)
plot_grid(p_rmse, p_corr, p_jsd, p_MAE, nrow = 2, ncol = 2, align = "hv") +
  theme(legend.position = legend_position) +
  cowplot::get_legend(p_rmse)



# 保存图像
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/sim_res.pdf', width = plot_size * 2, height = plot_size*2)





library(tidyr)
library(dplyr)
colnames(rmse)[1] = 'Value'
colnames(corr)[1] = 'Value'
colnames(jsd)[1] = 'Value'
colnames(MAE)[1] = 'Value'
long_rmse <- rmse %>% mutate(Metric = "RMSE")
long_corr <- corr %>% mutate(Metric = "Corr")
long_jsd <- jsd %>% mutate(Metric = "JSD")
long_MAE <- MAE %>% mutate(Metric = "MAE")

# 合并数据
combined_data_long <- bind_rows(long_rmse, long_corr, long_jsd, long_MAE)
combined_data_long$Col_Group <- factor(combined_data_long$`Heterozygosity rate`, levels = c("0", "0.3", "0.6"))
combined_data_long$Row_Group <- factor(combined_data_long$Metric, levels = c("RMSE", "Corr", "JSD", "MAE"))
custom_labeller <- labeller(
  Col_Group = c(
    `0` = "Heterogeneity rate: 0",
    `0.3` = "Heterogeneity rate: 0.3",
    `0.6` = "Heterogeneity rate: 0.6"
  )
)
combined_plot <- ggplot(combined_data_long, 
                        aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = "dodge") +  # 隐藏离群值，调整箱线图宽度
  facet_grid(Row_Group ~ Col_Group, scales = "free_y", labeller = custom_labeller) +  # 每个指标根据自己的范围调整y轴
  scale_fill_brewer(palette = "BuPu") +  # 使用柔和且高对比度的配色
  theme_minimal(base_size = 15) +  # 设置基础字体大小
  theme(
    legend.position = "bottom",  # 图例放在底部
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    axis.title.x = element_blank(),  # 删除横坐标标题
    axis.text.x = element_blank(),   # 删除横坐标标签
    axis.ticks.x = element_blank(),  # 删除横坐标刻度线
    axis.text.y = element_text(size = 16, face = "bold"),  # 调整Y轴文字
    axis.title = element_text(size = 20, face = "bold"),  # 调整轴标题
    strip.text = element_text(size = 20, face = "bold"),  # 调整分面标题
    strip.background = element_rect(fill = "grey80", color = NA),  # 设置灰色阴影背景
    panel.spacing.y = unit(0.5, "lines"),  # 增加面板之间的垂直间距
    panel.spacing.x = unit(2, "lines"),  # 保持水平间距
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold")  # 调整图形标题
  ) +
  labs(fill = "Method") +  # 自定义图例标题
  ylab("Metric Value")

# 添加显著性标注，并手动调整位置
combined_plot <- combined_plot +  
  geom_signif(
    comparisons = list(
      c("STsisal", "STdeconvolve"), 
      c("STsisal", "CARD"), 
      c("STsisal", "RCTD")
    ),
    map_signif_level = TRUE,
    margin_top = 0.15,  # 增加顶部边距
    tip_length = 0.05,
    textsize = 5,
    test = "t.test",
    step_increase = 0.15,  # 增加显著性标注的步长，确保标注不重叠
    vjust = 0.4  # 垂直调整显著性标注的位置
  )


# 打印并保存图像
print(combined_plot)
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/sim_res_with_significance_rate.pdf', combined_plot, width = 16, height = 14, dpi = 300)
