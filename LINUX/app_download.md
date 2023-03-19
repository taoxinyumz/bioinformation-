## 下载rz/sz安装及使用方法
~~~
 首先通过sftp工具把安 装文件上传到tmp目录下
 cd app
 tar zxvf lrzsz-0.12.20.tar.gz && cd lrzsz-0.12.20
  ./configure --prefix=$HOME/app/lrzsz && make && make install   #指定安装路径并安装该软件
 cd /data/home/taoxy/app/lrzsz/bin
 ln -s /data/home/taoxy/app/lrzsz/bin/lrz rz
 ln -s /data/home/taoxy/app/lrzsz/bin/lsz sz 
 PATH=/data/home/taoxy/app/lrzsz/bin:$PATH   #添加到环境变量中去
 source ~/.bashrc  #保存修改过后的bashrc文件
 ~~~
  










