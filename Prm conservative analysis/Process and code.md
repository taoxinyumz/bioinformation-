# 进行蛋白突变分析以确定蛋白序列的保守性
## 下载 entrez-direct
~~~
conda install -c bioconda entrez-direct # 在下载之前需要先激活conda环境
~~~

## 收集同源蛋白序列
+ NCBI中搜索prm蛋白质的氨基酸序列：打开NCBI网站，进入“Protein”数据库，在搜索框中输入Zika virus prM protein。得到完整蛋白的accession Q32ZE1.1 
+ 在Enter Query Sequence栏中输入蛋白质的序列号，并在“Choose Search Set" 中选择您要搜索的数据库。通常建议选择 "nr" 数据库，该数据库包含了NCBI收集的全部非冗余的蛋白质序列。随后点击BLAST。
+ 下载所有的accession号，并用entrez-direct在服务器上下载。
~~~
grep "^>" seqdump.txt | sed 's/>//g' > output.txt # 这个命令会从名为seqdump.txt的文本文件中提取所有的序列号，并将结果保存到名为output.txt的文本文件中
efetch -db protein -id "$(cat accession.txt)" -format fasta > zk.fasta # 根据所获取的序列号，下载相关的蛋白质的序列，并保存到zk.fasta文件中




























