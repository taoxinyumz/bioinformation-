## 代码流程
~~~
data.dir <- 'data.dir/'
list.files(data.dir)
expression <- Read10X(data.dir = data.dir)
yzjp_1 <- CreateSeuratObject(counts = expression,project = 'yzjp',min.cells = 3,min.features = 200,names.field = 1,names.delim = '_') #用于样本的分组，将barcode拆分成不同信息，names.field = 1表示要提取分隔符分割后的第1个元素，names.delim = "_"表示元信息的分隔符为"_"。names.field参数指定了将barcode拆分后取哪个元素作为meta.data中orig.ident列的属性。
yzjp_1
head(yzjp_1@meta.data) #meta.data是一个数据框，存储有关每个单元（cell）的元信息，例如样本来源、批次、是否为质控单元等。
yzjp_1[['percent.mt']] <- PercentageFeatureSet(yzjp_1,pattern = '^MT[-]')  #计算以mt开头的基因，即线粒体基因，并计算它们在细胞中的比例
yzjp_1[['percent.virus']] <- PercentageFeatureSet(yzjp_1,pattern = 'M1|M2|HA|NA|NP|NS1|NEP|PB1|PB1-F2|PA-X|PA|PB2') #正则表达式 "M1|M2|HA|NA|NP|NS1|NEP|PB1|PB1-F2|PA-X|PA|PB2" 匹配了这些病毒基因的名称，因此该函数会计算每个细胞中这些基因的表达量占所有基因表达量的百分比。
head(yzjp_1@meta.data$percent.mt)  #查看名为percent.mt的元数据
head(yzjp_1@meta.data$percent.virus)  #查看名为percent》virus的元件数
VlnPlot(yzjp_1, features = c('nFeature_RNA',"nCount_RNA",'percent.mt'),ncol = 3)
yzjp_1$log10GenesPerUMI <- log10(yzjp_1$nFeature_RNA)/log10(yzjp_1$nCount_RNA)    #计算每一个UMI所覆盖的基因数量，值越高表示测序质量越高
yzjp_1$mitoRatio <- PercentageFeatureSet(object = yzjp_1,pattern = "^MT")
yzjp_1$mitoRatio <- yzjp_1$mitoRatio/100  #计算MT基因所占的比例
metadata <- yzjp_1@meta.data   #创建metadata数据框
metadata$cells <- rownames(metadata)   #将细胞的ID添加到metadata中
metadata$sample <- NA  #创建一个名为sample的column
metadata$sample[which(str_detect(metadata$cells, "^d5_1"))] <- "d5_1"  #这段代码的含义是将metadata数据框中cells列以"d5_1"开头的元素的对应sample列的值设为"d5_1",其中which(str_detect(metadata$cells, "^d5_1"))这段代码会返回一个布尔向量，表示哪些元素（细胞）的cell名称以d5_1开头，然后使用which函数将这些元素的位置（下标）提取出来返回一个整数向量。
~~~

# 可视化每个样本中的细胞计数
~~~
pdf("d5_cell_counts.pdf",width = 6,height = 4) #函数中的参数"d5_cell_counts.pdf"指定了输出文件的名称，而width和height参数则分别指定了PDF文件的宽度和高度。
metadata %>%      #“将metadata传递给下一个函数
  ggplot(aes(x=sample, fill=sample)) +    #其中x=sample表示以数据集中的sample变量为x轴，fill=sample表示以sample变量的不同取值来填充（fill）颜色
  geom_bar(width =0.3) +      #geom_bar() 函数创建一个条形图，其中width = 0.3表示每个条形的宽度为0.3
  theme_classic() +      #theme_classic() 是 ggplot2 中的一个函数，它提供了一种基础主题样式，可以使绘图更加美观。
  geom_text(stat='count',aes(label=..count..), vjust=1.6, color="white", size=3.5)+   #通过 geom_text 函数来实现。其中，stat='count' 表示使用 ggplot2 统计每个 x 变量出现的次数，aes(label=..count..) 表示使用计数值作为文本标签，vjust=1.6 表示文本在垂直方向上偏移 1.6 个单位（相对于原来位置），color="white" 表示文本颜色为白色，size=3.5 表示文本大小为 3.5。
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +   #axis.text.x是指定X轴的文本标签样式。element_text是指定文本标签的样式类型。angle设置文本标签与X轴之间的旋转角度。vjust设置文本标签的垂直对齐方式。hjust设置文本标签的水平对齐方式。
  theme(plot.title = element_text(hjust=0.5, face="bold")) +  #theme(plot.title = ) 表示设置图形标题的样式；element_text(hjust=0.5, face="bold") 表示将文本标签的水平对齐方式设为0.5，即居中对齐，同时将字体样式设置为粗体。
  ggtitle("NCells")      #ggtitle("New Title")会将主标题改为"New Title"
dev.off()   #dev.off() 是 R 语言中的一个函数，用于关闭绘图设备，将绘制好的图形保存到文件或显示到屏幕上
~~~

# 可视化每个细胞的UMI/转录本数量
~~~
pdf("d5_UMIs_counts.pdf",width = 6,height = 4)         #在R语言中用于生成一个宽度为6英寸、高度为4英寸的PDF文件，命名为"d5_UMIs_counts.pdf"，后续图形绘制的结果会输出到这个文件中。
metadata %>%                                           #metadata %>%将metadata数据框作为管道的输入，将其传递给后续的函数。这种方法可以使代码更加简洁和易读，尤其是在多个函数调用和数据转换的情况下。
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +    #ggplot()创建了一个空白的坐标系,aes()指定了数据中每个变量对应到图形中的哪个视觉元素。在这里，color=sample将sample变量映射到点的颜色，x=nUMI将nUMI变量映射到点的位置上。
  geom_density(alpha = 0.4) +                          #geom_density()用于绘制核密度图.alpha参数控制了点的透明度。值为0表示完全透明，1表示完全不透明
  scale_x_log10() +                                    #
  theme_classic() +                                    #
  ylab("Cell density") +                               #可以使用 ylab("Cell density") 来将 y 轴的标签设置为 "Cell density"
  geom_vline(xintercept = 500)                         #在 ggplot2 绘图中添加一条竖直位置在 X 轴上 500 位置处的参考线
dev.off()                                              #与PDF指令相对应
**横座标是nUMI，纵座标是cell density的图片，该怎么看**:
图中纵坐标的cell density表示的是在某个nUMI区间内，有多少个细胞。如果某个nUMI区间的cell density较高，则表示这个区间内的细胞数量比较多，这个基因表达水平的细胞数量相对较多。
~~~

# 可视化每个细胞检测到的基因分布
~~~
pdf("yzjp_per_cell.pdf",width = 6,height = 4)                 #
metadata %>%                                                  #
  ggplot(aes(color=sample, x=nGene, fill= sample)) +          #表示用nGene作为横坐标，以样本为区分颜色并填充。
  geom_density(alpha = 0.2) +                                 #geom_density函数，表示将密度曲线添加到图表上，alpha=0.2表示透明度为0.2。
  theme_classic() +                                           #
  scale_x_log10() +                                           #
  geom_vline(xintercept = 300)                                #
dev.off()                                                     #
**x轴是nGene y轴是geom_density的图，其中曲线较高的地方有什么含义**:
曲线较高的地方表示在该nGene值下有更多的单个细胞，即单个细胞在这个基因集上的表达水平更高。这可能是由于这些基因在细胞类型或状态中的特异性表达或增强的功能，或者与细胞功能密切相关。
~~~

## 通过箱线图可视化每个细胞检测到的基因分布
~~~
pdf("yzjp_detected_per_cell_boxplot.pdf",width = 5,height = 4)
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +  #x轴表示不同的样本，y轴表示每个单元的log10基因数，fill参数指定根据不同的样本填充颜色
  geom_boxplot(notch=TRUE,notchwidth = 0.3,outlier.size = 1) +  #添加一个箱线图，其中notch参数用于绘制缺口，notchwidth参数指定缺口的宽度，outlier.size参数指定异常值的大小
  stat_summary(fun = mean, geom = "point", shape = 23, size=4)+  #添加一个点图，其中stat_summary函数用于计算每个样本的均值，geom参数指定为点图，shape参数指定点的形状，size参数指定点的大小；
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +  #使用theme_classic函数设定绘图的主题为经典风格，axis.text.x参数用于设置x轴标签的旋转角度和位置；
  theme(plot.title = element_text(hjust=0.5, face="bold")) +  #
  ggtitle("NCells vs NGenes")  #
dev.off()
~~~

## color = sample和fill = sample 的区别是什么？
~~~
color = sample是指定每个数据点的颜色，它所对应的绘图元素是数据点的边框。fill = sample是指定每个数据点的填充颜色，它所对应的绘图元素是数据点的内部。
~~~

## Vjust和Hjust的意义是什么？
~~~
Vjust表示垂直对齐方式。vjust=1（文本沿底部与图形对齐） vjust=0（文本沿顶部与图形对齐）vjust=0.5（文本垂直居中对齐）
Hjust表示水平对齐方式。Hjust=0（表示文本左对齐） Hjust=1(表示文本右对齐) Hjust=0.5(文本居中对齐)
~~~

## 可视化检测到的基因与 UMI 数量之间的相关性#并确定是否大量存在基因/UMI 数量少的细胞
~~~
  pdf("d5_genes_detected_and_number_of_UMIs.pdf",width = 6,height = 4)
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +        # mitoRatio表示一个细胞中线粒体基因的比例，通常用于评估RNA测序数据的质量。高的mitoRatio值可能表示细胞死亡或处理过程中的损伤，因为线粒体的RNA会更加稳定，相对于核基因而言，在细胞死亡或受到损伤后更加持久。
  geom_point() +      # 添加一个散点图，横座标是x变量，纵座标是Y变量
  scale_colour_gradient(low = "gray90", high = "black") +   # 
  stat_smooth(method=lm) +       # ggplot2 中，这个函数通常与 aes 函数一起使用，用于设定图形中的 X 轴和 Y 轴的变量映射，从而对数据进行可视化。
  scale_x_log10() +              
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +   #  画线x
  geom_hline(yintercept = 250) +   #  画线y
  facet_wrap(~sample)  # 按照sample进行分组展示
dev.off()
~~~

## 可视化每个细胞检测到的线粒体基因表达分布
~~~
pdf("d5_mitochondrial_gene_expression.pdf",width = 6,height = 4)
metadata %>%        
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +   # color指定样本分组的颜色映射，x指定x轴上显示的变量，fill指定样本分组的填充颜色映射。
  geom_density(alpha = 0.2) +    # 在核密度估计图中，横坐标表示变量的取值范围，纵坐标表示概率密度，曲线下的面积为1。
  scale_x_log10() +   # 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()
~~~

# 将基因表达的总体复杂性可视化
## 可视化每个UMI检测到的基因数
~~~
pdf("d5_genes_detected_per_UMI.pdf",width = 6,height = 4)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +  # aes映射，其中 x 轴为 log10GenesPerUMI，颜色和填充颜色根据 sample 变量进行映射。
  geom_density(alpha = 0.2) +                                       
  theme_classic() +
  geom_vline(xintercept = 0.8)       # xintercept 参数指定垂直线的位置，此处设置为0.8，即在x轴上的0.8处添加一条垂直线。
dev.off()
~~~

# 查看病毒基因情况
## 检查四分位数值
~~~
summary(d5_merged_seurat@meta.data$percent.virus)      # 这段代码是对 d5_merged_seurat 数据集的元数据中的 percent.virus 变量进行汇总统计，并返回其基本的描述性统计信息，如最大值、最小值、中位数、平均值和四分位数等。
View(d5_merged_seurat@meta.data)                       # 查看元数据
~~~
# 基于四分位数值将 mitoRatio 转换为分类因子向量
~~~
d5_merged_seurat@meta.data$percent_virus <- cut(d5_merged_seurat@meta.data$percent.virus,      # 
                                                breaks=c(-Inf, 0,0.05,0.5,5, Inf), 
                                                labels=c("Per=0%,n=12757","Per<0.05%,n=8482",
                                                         "Per0.05-0.5%,n=10349","Per0.5-5,n=1180",
                                                         "Per>5,n=417"))
Freq <- as.data.frame(table(d5_merged_seurat@meta.data$percent_virus))
pdf("d5_FeatureScatter_percent_virus.pdf",width = 5,height = 4)
FeatureScatter(object = d5_merged_seurat, feature1 = 'nGene', feature2 = 'percent.virus',
               group.by = "sample")
dev.off()
~~~

## 一些错误及解决办法
~~~
1.显示失去了wt权限

~~~




