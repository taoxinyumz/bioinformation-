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

## 一些错误及解决办法
~~~
1.显示失去了wt权限

~~~




