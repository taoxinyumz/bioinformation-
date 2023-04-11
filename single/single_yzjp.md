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

# 可视化每个样本中的细胞计数
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

# 


~~~

## 一些错误及解决办法
~~~
1.显示失去了wt权限

~~~




