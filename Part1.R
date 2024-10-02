#Diana Escalante A01253468
#Libia Flores A01253425
library(ggplot2)
library(seqinr)
library(gridExtra)

#=============Funciones==============
#Función Tamaño
calculartam <- function(secuencia){
  return(length(secuencia[[1]]))
}
#Contar nucleótidos
contarnuc2 <- function(secuencia){
  s <- calculartam(secuencia)
  cA <- 0
  cC <- 0
  cG <- 0
  cT <- 0
  for (i in secuencia[[1]]){
    if(i=="a"){
      cA=cA+1
    }
    if(i=="c"){
      cC=cC+1
    }
    if(i=="g"){
      cG=cG+1
    }
    if(i=="t"){
      cT=cT+1
    }
  }
  
  return(data.frame(a =cA, c=cC, g=cG, t=cT))
}
#Conteo de GC
contarGC <- function(secuencia){
  s <- calculartam(secuencia)
  cC <- 0
  cG <- 0
  for (i in secuencia[[1]]){
    if(i=="c"){
      cC=cC+1
    }
    if(i=="g"){
      cG=cG+1
    }
  }
  return (((cG+cC)*100)/s)
}
#=======Leer Bases de Datos==========
setwd("C:/Users/libia/OneDrive/Escritorio/2DOSEMESTRE/AnalisisdeBiologiaComputacional/Evidencia1_genomas")

#=======Hemisferio Norte=============
usa1<- read.fasta(file = "usa1.fasta")
usa2<- read.fasta(file = "usa2.fasta")
usa3<- read.fasta(file = "usa3.fasta")
canada1<- read.fasta(file = "canada1.fasta")
canada2<- read.fasta(file = "canada2.fasta")
canada3<- read.fasta(file = "canada3.fasta")
#========Hemisferio Sur===============
brasil1<- read.fasta(file = "brasil1.fasta")
brasil2<- read.fasta(file = "brasil2.fasta")
brasil3<- read.fasta(file = "brasil3.fasta")
uruguay1<- read.fasta(file = "uruguay1.fasta")
uruguay2<- read.fasta(file = "uruguay2.fasta")
uruguay3<- read.fasta(file = "uruguay3.fasta")

#================IDs===================
#USA
idusa1 <- names(usa1)[1]
idusa2 <- names(usa2)[1]
idusa3 <- names(usa3)[1]
#Canada
idcanada1 <- names(canada1)[1]
idcanada2 <- names(canada2)[1]
idcanada3 <- names(canada3)[1]
#Brasil
idbrasil1 <- names(brasil1)[1]
idbrasil2 <- names(brasil2)[1]
idbrasil3 <- names(brasil3)[1]
#Uruguay
iduruguay1 <- names(uruguay1)[1]
iduruguay2 <- names(uruguay2)[1]
iduruguay3 <- names(uruguay3)[1]

#=======Tamaño de cada secuencia=======
#USA
tamusa1 <- calculartam(usa1)
tamusa2 <- calculartam(usa2)
tamusa3 <- calculartam(usa3)
tamusa1
tamusa2
tamusa3

#Canada
tamcanada1 <- calculartam(canada1)
tamcanada2 <- calculartam(canada2)
tamcanada3 <- calculartam(canada3)
tamcanada1
tamcanada2
tamcanada3

#Brasil
tambrasil1 <- calculartam(brasil1)
tambrasil2 <- calculartam(brasil2)
tambrasil3 <- calculartam(brasil3)
tambrasil1
tambrasil2
tambrasil3

#Uruguay
tamuruguay1 <- calculartam(uruguay1)
tamuruguay2 <-calculartam(uruguay2)
tamuruguay3 <-calculartam(uruguay3)
tamuruguay1
tamuruguay2
tamuruguay3

#======Gráfica Longitud de Genomas=====
dataframelong <- data.frame(ID=c(idusa1, idusa2, idusa3, idcanada1, idcanada2, idcanada3, idbrasil1, idbrasil2, idbrasil3, iduruguay1, iduruguay2, iduruguay3),
                            Longitud=c(tamusa1, tamusa1, tamusa3, tamcanada1, tamcanada2, tamcanada3, tambrasil1, tambrasil2, tambrasil3, tamuruguay1, tamuruguay2, tamuruguay3))
graficalong <- ggplot(dataframelong, aes(x=ID, y=Longitud, fill=ID))+geom_bar(stat="identity")
graficalong

#=====Composición de Nucleótidos=======
#USA
cnusa1 <- contarnuc2(usa1)
cnusa2 <- contarnuc2(usa2)
cnusa3 <- contarnuc2(usa3)
cnusa1
cnusa2
cnusa3

#Canada
cncanada1 <- contarnuc2(canada1)
cncanada2 <- contarnuc2(canada2)
cncanada3 <- contarnuc2(canada3)
cncanada1
cncanada2
cncanada3

#Brasil
cnbrasil1 <- contarnuc2(brasil1)
cnbrasil2 <- contarnuc2(brasil2)
cnbrasil3 <- contarnuc2(brasil3)
cnbrasil1
cnbrasil2
cnbrasil3

#Uruguay
cnuruguay1 <- contarnuc2(uruguay1)
cnuruguay2 <- contarnuc2(uruguay2)
cnuruguay3 <- contarnuc2(uruguay3)
cnuruguay1
cnuruguay2
cnuruguay3


#===Gráfica Composición de Nucleótidos===
#USA
dataframeusa1 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnusa1$a, cnusa1$c, cnusa1$g, cnusa1$t))
dataframeusa2 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnusa2$a, cnusa2$c, cnusa2$g, cnusa2$t))
dataframeusa3 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnusa3$a, cnusa3$c, cnusa3$g, cnusa3$t))
grafusa1 <- ggplot(dataframeusa1, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idusa1)
grafusa2 <- ggplot(dataframeusa2, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idusa2)
grafusa3 <- ggplot(dataframeusa3, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idusa3)
#CANADA
dataframecanada1 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cncanada1$a, cncanada1$c, cncanada1$g, cncanada1$t))
dataframecanada2 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cncanada2$a, cncanada2$c, cncanada2$g, cncanada2$t))
dataframecanada3 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cncanada3$a, cncanada3$c, cncanada3$g, cncanada3$t))
grafcanada1 <- ggplot(dataframecanada1, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idcanada1)
grafcanada2 <- ggplot(dataframecanada2, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idcanada2)
grafcanada3 <- ggplot(dataframecanada3, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idcanada3)
#Brasil
dataframebrasil1 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnbrasil1$a, cnbrasil1$c, cnbrasil1$g, cnbrasil1$t))
dataframebrasil2 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnbrasil2$a, cnbrasil2$c, cnbrasil2$g, cnbrasil2$t))
dataframebrasil3 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnbrasil3$a, cnbrasil3$c, cnbrasil3$g, cnbrasil3$t))
grafbrasil1 <- ggplot(dataframebrasil1, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idbrasil1)
grafbrasil2 <- ggplot(dataframebrasil2, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idbrasil2)
grafbrasil3 <- ggplot(dataframebrasil3, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=idbrasil3)

#Uruguay
dataframeuruguay1 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnuruguay1$a, cnuruguay1$c, cnuruguay1$g, cnuruguay1$t))
dataframeuruguay2 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnuruguay2$a, cnuruguay2$c, cnuruguay2$g, cnuruguay3$t))
dataframeuruguay3 <- data.frame(Nucleótido=c("a", "c", "g", "t"),Freq=c(cnuruguay3$a, cnuruguay3$c, cnuruguay3$g, cnuruguay3$t))
grafuruguay1 <- ggplot(dataframeuruguay1, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=iduruguay1)
grafuruguay2 <- ggplot(dataframeuruguay2, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=iduruguay2)
grafuruguay3 <- ggplot(dataframeuruguay3, aes(x=Nucleótido, y=Freq, fill=Nucleótido))+geom_bar(stat="identity")+labs(title=iduruguay3)

graficacn <- grid.arrange(grafusa1, grafusa2, grafusa3, grafcanada1, grafcanada2, grafcanada3, grafbrasil1, grafbrasil2, grafbrasil3, grafuruguay1, grafuruguay2, grafuruguay3, ncol=4)
graficacn

#========Porcentaje de GC=======

#USA
cgusa1 <- contarGC(usa1)
cgusa2 <- contarGC(usa2)
cgusa3 <- contarGC(usa3)
cgusa1
cgusa2
cgusa3

#Canada
cgcanada1 <- contarGC(canada1)
cgcanada2 <- contarGC(canada2)
cgcanada3 <- contarGC(canada3)
cgcanada1
cgcanada2
cgcanada3

#Brasil
cgbrasil1 <- contarGC(brasil1)
cgbrasil2 <- contarGC(brasil2)
cgbrasil3 <- contarGC(brasil3)
cgbrasil1
cgbrasil2
cgbrasil3

#Uruguay
cguruguay1 <- contarGC(uruguay1)
cguruguay2 <- contarGC(uruguay2)
cguruguay3 <- contarGC(uruguay3)
cguruguay1
cguruguay2
cguruguay3

#===Gráfica Porcentaje de GC===
#USA
dataframegcusa1 <- data.frame(ID=idusa1, ContenidoGC=cgusa1)
dataframegcusa2 <- data.frame(ID=idusa2, ContenidoGC=cgusa2)
dataframegcusa3 <- data.frame(ID=idusa3, ContenidoGC=cgusa3)
grafgcusa1 <- ggplot(dataframegcusa1, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idusa1)
grafgcusa2 <- ggplot(dataframegcusa2, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idusa2)
grafgcusa3 <- ggplot(dataframegcusa3, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idusa3)

#Canada
dataframegccanada1 <- data.frame(ID=idcanada1, ContenidoGC=cgcanada1)
dataframegccanada2 <- data.frame(ID=idcanada2, ContenidoGC=cgcanada2)
dataframegccanada3 <- data.frame(ID=idcanada3, ContenidoGC=cgcanada3)
grafgccanada1 <- ggplot(dataframegccanada1, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idcanada1)
grafgccanada2 <- ggplot(dataframegccanada2, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idcanada2)
grafgccanada3 <- ggplot(dataframegccanada3, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idcanada3)

#Brasil
dataframegcbrasil1 <- data.frame(ID=idbrasil1, ContenidoGC=cgbrasil1)
dataframegcbrasil2 <- data.frame(ID=idbrasil2, ContenidoGC=cgbrasil2)
dataframegcbrasil3 <- data.frame(ID=idbrasil3, ContenidoGC=cgbrasil3)
grafgcbrasil1 <- ggplot(dataframegcbrasil1, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idbrasil1)
grafgcbrasil2 <- ggplot(dataframegcbrasil2, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idbrasil2)
grafgcbrasil3 <- ggplot(dataframegcbrasil3, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=idbrasil3)

#Uruguay
dataframegcuruguay1 <- data.frame(ID=iduruguay1, ContenidoGC=cguruguay1)
dataframegcuruguay2 <- data.frame(ID=iduruguay2, ContenidoGC=cguruguay2)
dataframegcuruguay3 <- data.frame(ID=iduruguay3, ContenidoGC=cguruguay3)
grafgcuruguay1 <- ggplot(dataframegcuruguay1, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=iduruguay1)
grafgcuruguay2 <- ggplot(dataframegcuruguay2, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=iduruguay2)
grafgcuruguay3 <- ggplot(dataframegcuruguay3, aes(x=ID, y=ContenidoGC, fill=ContenidoGC))+geom_bar(stat="identity")+labs(title=iduruguay3)

graficagc <- grid.arrange(grafgcusa1, grafgcusa2, grafgcusa3, grafgccanada1, grafgccanada2, grafgccanada3, grafgcbrasil1, grafgcbrasil2, grafgcbrasil3, grafgcuruguay1, grafgcuruguay2, grafgcuruguay3, ncol=4)
graficagc

#============Marco de Datos==============
marcodedatos <- data.frame(IDs=c(idusa1, idusa2, idusa3, idcanada1, idcanada2, idcanada3, idbrasil1, idbrasil2, idbrasil3, iduruguay1, iduruguay2, iduruguay3),
                           Longitud=c(tamusa1, tamusa1, tamusa3, tamcanada1, tamcanada2, tamcanada3, tambrasil1, tambrasil2, tambrasil3, tamuruguay1, tamuruguay2, tamuruguay3),
                           ContenidoA=c(cnusa1$a, cnusa2$a, cnusa3$a, cncanada1$a, cncanada2$a, cncanada3$a, cnbrasil1$a, cnbrasil2$a, cnbrasil3$a, cnuruguay1$a, cnuruguay2$a, cnuruguay3$a),
                           ContenidoC=c(cnusa1$c, cnusa2$c, cnusa3$c, cncanada1$c, cncanada2$c, cncanada3$c, cnbrasil1$c, cnbrasil2$c, cnbrasil3$c, cnuruguay1$c, cnuruguay2$c, cnuruguay3$c),
                           ContenidoG=c(cnusa1$g, cnusa2$g, cnusa3$g, cncanada1$g, cncanada2$g, cncanada3$g, cnbrasil1$g, cnbrasil2$g, cnbrasil3$g, cnuruguay1$g, cnuruguay2$g, cnuruguay3$g),
                           ContenidoT=c(cnusa1$t, cnusa2$t, cnusa3$t, cncanada1$t, cncanada2$t, cncanada3$t, cnbrasil1$t, cnbrasil2$t, cnbrasil3$t, cnuruguay1$t, cnuruguay2$t, cnuruguay3$t)
                           )
marcodedatos
