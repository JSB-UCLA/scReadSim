# Install required softwares of scReadSim

This tutorial introduces the bash commands used to install the following required softwares for scReadSim:
- [MACS3 >= 3.0.0](https://github.com/macs3-project/MACS)
- [samtools >= 1.12](http://www.htslib.org/)
- [bedtools >= 2.29.1](https://bedtools.readthedocs.io/en/latest/)
- [seqtk >= 1.3](https://github.com/lh3/seqtk)
- [fgbio >= 2.0.1](https://github.com/fulcrumgenomics/fgbio)

Specify the absolute path where the softwares will be installed. **Note that user needs to use their own absolute path.**  

```{code-block} console
$ Tools_path=/path/to/your/folder 
```

## MACS3
For installation details, please check [MACS3 website](https://github.com/macs3-project/MACS). The final command returns the path to macs3 software, which is a required input of scReadSim.

```{code-block} console
$ pip install macs3
$ macs3_executefile=$(which macs3)
$ echo ${macs3_executefile%/macs3}
```

## samtools
For installation details, please check [samtools website](http://www.htslib.org/). The final command returns the path to samtools software, which is a required input of scReadSim.

```{code-block} console
$ cd ${Tools_path}
$ sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
$ wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 -O samtools.tar.bz2
$ tar -xjvf samtools.tar.bz2
$ cd samtools-1.12
$ make
$ sudo make install
$ export PATH=$PATH:${Tools_path}/samtools-1.12
$ echo ${Tools_path}/samtools-1.12
```


## bedtools
For installation details, please check [bedtools website](https://bedtools.readthedocs.io/en/latest/). The final command returns the path to bedtools software, which is a required input of scReadSim.

```{code-block} console
$ cd ${Tools_path}
$ wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
$ tar -zxvf bedtools-2.29.1.tar.gz
$ cd bedtools2
$ make
$ export PATH=$PATH:${Tools_path}/bedtools2/bin
$ echo ${Tools_path}/bedtools2/bin
```


## seqtk
For installation details, please check [seqtk website](https://github.com/lh3/seqtk). The final command returns the path to seqtk software, which is a required input of scReadSim.

```{code-block} console
$ cd ${Tools_path}
$ git clone https://github.com/lh3/seqtk.git
$ cd seqtk
$ make
$ export PATH=$PATH:${Tools_path}/seqtk
$ echo ${Tools_path}/seqtk
```



## fgbio
For installation details, please check [fgbio website](https://github.com/fulcrumgenomics/fgbio). 

The installation of fgbio needs build tool [sbt](https://www.scala-sbt.org/download.html) to compile. First install sbt.


```{code-block} console
$ echo "deb https://repo.scala-sbt.org/scalasbt/debian all main" | sudo tee /etc/apt/sources.list.d/sbt.list
$ echo "deb https://repo.scala-sbt.org/scalasbt/debian /" | sudo tee /etc/apt/sources.list.d/sbt_old.list
$ curl -sL "https://keyserver.ubuntu.com/pks/lookup?op=get&search=0x2EE0EA64E40A89B84B2DF73499E82A75642AC823" | sudo apt-key add
$ sudo apt-get update
$ sudo apt-get install sbt 
$ sbt -h
```

Now install fgbio. The final command returns the absolute path of fgbio java script, which is a required input of scReadSim. **Note**: user needs to specify the corresponding version for the final command.


```{code-block} console
$ cd ${Tools_path}
$ git clone https://github.com/fulcrumgenomics/fgbio.git
$ cd fgbio
$ sbt assembly
$ echo ${Tools_path}/fgbio/target/scala-{version}/fgbio-{version}-SNAPSHOT.jar
```