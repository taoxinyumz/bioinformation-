## 获取SRA数据库中感兴趣样本的Run ID的方法
+ 1.搜索SRA库(样本多)。NCBI的SRA栏-搜物种名-Public-DNA-genome-PacBio SMRT-fastq-"Send results to Run selector"(或进入某个筛选结果并点击"All runs")
-下载Metadata (SraRunTable.txt)。**Metadata是SRA数据的元数据，包括关于数据来源、样本、实验、测序平台等详细信息的描述**。
+ 2.获取"SraRunTable.txt"后 (含Run ID、样本的各种属性)，将文件后缀**改为.csv** ，用Excel打开、自定义排序、根据研究兴趣筛选特定样本(增加或删减)。按照以上方法，可
按照这些列汇集不同的表格(注意顺序)：Sample、Platform、Run、Collection_Date、Country、Location、Assay_Type... 手动(或编写自动化程序)汇集成大的表格，另存为"文
本文件(制表符分隔)(*.txt)文件"，文件名：metadata.txt，上传至`~/wgs/fastq/bak/`目录下
 根据Run的编号下载fastq至`~/wgs/fastq/bak/`目录，代码如下：
~~~
 # 提前用Excel按照“Sample、Platform、Run ...”将Run放在第3列
sed 's/\r//g' ~/wgs/fastq/bak/metadata.txt | tail -n +2 | cut -f 3 | \
   ~/wgs/public/bin/rush -j 25 "~/wgs/soft/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump \
   -v --split-3 --gzip {}
~~~

## 包含SRA的accession号的csv文件，怎么根据该文件下载序列
~~~
while read sra; do
    fastq-dump --split-files $sra
done < sra_file.csv
~~~
















