#require("mvtnorm")
gear.names = c("Tuna HL","Troll","Longline","Bottom/inshore HL","Aku Boat",
                "Klingon")
sgn = c("THL","Troll","LL","BHL","Aku","K")
have.xssams.R = FALSE
have.issams.R = TRUE
start.year = 1952


plot.diagnostics=function(dat=NULL,file="diagnostics.dat",dt,ngear,
                 sdlogPop, sdlogYield, sdlogF, sdlogQ, K, r, B1, indexed,
                 plot.Fmort,plot.prod, plot.impact,
                 devices,block)
{
   print("plot.diagnostics, issams.R")
   print(paste("Block: ",block))
   if (is.null(dat))
     dat = read.table(file=file,header=TRUE)

   if (dat$t[1] == 1)
   {
      dat$t = (start.year-0.4*dt +  dat$t*dt)
   #  print(names(dat))
   }

   wy5 = which((floor(dat$t%%5))==0)
   ncol = ncol(dat)
#  print("Names of all variables:")
#  print(names(dat))
#  print(head(dat))
#  print(tail(dat))
   gear.col = 4
   dd = c(2,(gear.col+1):ncol)
   print("Names of log transformed variables:")
   print(names(dat)[dd])
   dat[,dd] = exp(dat[,dd])
#  print(head(dat))
#  print(tail(dat))

   #gear.names = c("TunaHL","Troll","Longline","Bottom/inshore HL","AkuBoat")
   title.line = -1
   lwd = 3
   #sd.lwd = 3
   #sd.lty = "dotted"
   #old.par = par(no.readonly = TRUE) 

   F.ndx=grep("F",names(dat)) 
   predC.ndx=grep("predC",names(dat)) 
   obsC.ndx=grep("obsC",names(dat)) 

   # biomass plots
   d = 1
   xpos = 0
   ypos = 0
   if (devices[d] > 0)
   {
      s = dev.set(devices[d])
   }
   else
   {
      width = 9.0
      height = 9.0
      xpos = 100
      ypos = 0
      x11(width=width,height=height,xpos=xpos,ypos=ypos,title="Biomass")
      devices[d] = dev.cur()
   }

   par(mar=c(4,5,0,4)+0.1)
   ntime = length(dat$t)
   x = dat$t
   y = matrix(nrow=ntime,ncol=1)
   y[,1]=dat$pop
#  y[,2]=dat$forcing
   print(paste("sdlogPop",sdlogPop))
   plot.biomass(dat$t,y,sd=sdlogPop,block=block,K=dat$K,B1=B1,
                forcing=dat$forcing,indexed=indexed)

#  options(scipen=6)
#  xrange=nice.ts.plot(x,y,ylab="Biomass (mt)")

#  lines(x,dat$K,lwd=2,lty="dotdash",col="blue",xlim=xrange)
#  text(x[ntime],dat$K[ntime],"  K",adj=c(0,0.5),col="blue")

#  plot.error(dat$t,dat$pop,sdlogPop,bcol="blue",fcol="lightblue")
#  double.lines(dat$t,dat$pop,bcol="blue",fcol="lightblue",lwd=5)
#  text(x[ntime],dat$pop[ntime],adj=c(0,1),col="blue",labels="  N")

#  if (max(dat$forcing,na.rm=TRUE) > 1.0)
#  {
#     plot.error(x,dat$forcing,sdlogQ,bcol="purple4",fcol="purple1")
#     double.lines(x,dat$forcing,bcol="purple4",fcol="purple1",lwd=5) 
#     text(x[ntime],dat$forcing[ntime],adj=c(0,0),col="purple4",labels="  I")
#  }
#  show.block.number(block,dat$t[1])

   # catch plots
   d = 2
   if (devices[d] > 0)
   {
      s = dev.set(devices[d])
   }
   else
   {
     width = 9.0
     height =11.0
     xpos = xpos + 50
     ypos = ypos + 50
     x11(width=width,height=height,xpos=xpos,ypos=ypos,title="Catch")
     devices[d] = dev.cur()
   }

   plot.catches(dat$t,dat[,c(gear.col+2*ngear+1):(gear.col+3*ngear)],
                      dat[,c(gear.col+ngear+1)  :(gear.col+2*ngear)],
                      sdlogYield,block) 

#  plot.Fmort = TRUE
   if (plot.Fmort)
   {
       d = d + 1
       if (devices[d] > 0)
       {
          s = dev.set(devices[d])
       }
       else
       {
          width = 9.0
          height =11.0
          xpos = xpos + 50
          ypos = ypos + 50
          x11(width=width,height=height,xpos=xpos,ypos=ypos,title="Fishing Mortality")
          devices[d] = dev.cur()
       }
       par(mar=c(3,4.5,0,0)+0.1)
       np = ngear
       lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
       layout.show(lm)
       for (g in 1:ngear)
       {
       #  labelsY=parse(text=paste(seq(50,100,10), "^o ", "*N", sep="")) 
          Flab = parse(text=paste("F~(","y^-1",")",sep=""))
          nice.ts.plot(dat$t,dat[,(gear.col+g)], ylab=Flab,
                bcol="orange4",fcol="orange", lwd=lwd)

          plot.error(dat$t,dat[,(gear.col+g)],sdlogF,
                    bcol="orange4",fcol="orange")
          lines(dat$t,dat[,(gear.col+g)],col="orange4",lwd=lwd+2)
          title(main=gear.names[g],line=title.line)
       }
       if (!is.null(block))
          show.block.number(block,dat$t[1],line=2)
   } #if (plot.Fmort)

   if (plot.prod)
   {
       d = d + 1
       if (devices[d] > 0)
       {
          s = dev.set(devices[d])
       }
       else
       {
          width = 9.0
          height =9.0
          xpos = xpos + 50
          ypos = ypos + 50
          x11(width=width,height=height,xpos=xpos,ypos=ypos,
               title="Production")
          devices[d] = dev.cur()
       }
       Fmort = rowSums(dat[,F.ndx])
       obsC = rowSums(dat[,obsC.ndx])
       predC = rowSums(dat[,predC.ndx])
       plot.production(Fmort,obsC,predC,dat$t, r,K,block)


   } #if (plot.prod)

   if (plot.impact)
   {
       d = d + 1
       if (devices[d] > 0)
       {
          s = dev.set(devices[d])
       }
       else
       {
          width = 9.0
          height =9.0
          xpos = xpos + 50
          ypos = ypos + 50
          x11(width=width,height=height,xpos=xpos,ypos=ypos,
               title="Impact")
          devices[d] = dev.cur()
       }
       par(mar=c(4,5,0,4)+0.1)
       ntime = length(dat$t)
       x = dat$t
       Fmort = rowSums(dat[,F.ndx])
       y = matrix(nrow=ntime,ncol=3)
       y[,3]=dat$pop
       y[1,1] = y[1,3]
       y[1,2] = y[1,3]
       for (t in 2:ntime)
       {
          for (p in 1:2)
          {
             FF = Fmort[t]
             if (p == 2)
                FF = 0.0

             y[t,p]  = (K*(r-FF))/((((K*(r-FF))/y[t-1,p])*exp(-(r-FF))) - r*exp(-(r-FF))  + r) # p5
          } 
      }

      legend = c(" F (est)"," F= 0", " B (est)")
      xrange=nice.ts.plot(x,y[,1:2],ylab="Biomass (mt)",legend=legend[1:2])

   #  impact = 100.0*(1.0 - y[,1]/y[,2]) # PNAS
   #  impact = (1.0 - y[,1]/y[,2])
      impact = y[,1]/y[,2]  # WCPFC
   #  print(impact)
      par("new"=TRUE)
      plot(x,impact,type='l',xlim=xrange,lwd=3,lty="dashed",col="red",
           ann=FALSE,axes=FALSE)
      axis(4,col="red",ylab="p",col.axis="red")
      #mtext("p",side=4,col="red",line=0.1)

   } #if (plot.impact)
   new.devices = devices
   return(new.devices)
}

log.diagnostics=function(file="issams_program.log",ntime=61,dt=1,ngear=5,
                         plot.Fmort=FALSE,plot.prod=FALSE,plot.impact=FALSE)
{
      
   print(paste("Scanning file",file))
   log = scan(file,what="character")
   res = grep("Residuals:",log)

   max.counter = length(res)
   counter = max.counter

   c = 'n'
   dev.list = vector(mode="numeric",length=5)
#  dev.file.names=c("tmp/est_pop","tmp/est_catch","tmp/est_F","tmp/prod")
   dev.file.names=c("est_pop","est_catch","est_F","prod")
#  while (c != 'q')
   while ( (c != 'q') && (c != 'x') )
   {
      tmp=get.diagnostics(log,ntime=ntime,dt=dt,ngear=ngear,block=counter,mtype="i")

      new.devices = plot.diagnostics(tmp$resid,dt=dt,ngear=ngear,
                    sdlogPop=tmp$sdlogPop,
                    sdlogYield=tmp$sdlogYield,
                    sdlogF=tmp$sdlogF,
                    sdlogQ=tmp$sdlogQ,
                    K=tmp$K, r=tmp$r, B1=tmp$B1, indexed=tmp$indexed,
                    plot.Fmort=plot.Fmort,
                    plot.prod=plot.prod,
                    plot.impact=plot.impact,
                    devices=dev.list,block=counter)


      dev.list=new.devices
     
      c = readline("next, back, save, quit or exit? (enter n,b,s,q,x):")
      print(paste(c," entered"))
      if (c == 'n')
         counter = counter + 1
      else if (c == 'b')
         counter = counter - 1
      else if (c == 's')
      {
         print(dev.list)
         for (d in 1:length(dev.list))
         {
            if (dev.list[d] > 0)
            {
               tname = paste(dev.file.names[d],"B",counter,sep="")
               print(paste("Saving device",d,"to file",tname))
               dev.set(dev.list[d])
               din = par("din")
               print(din)
               save.png.plot(tname,width=din[1],height=din[2])
            }
            else
               print(paste("Skipping device",d))
          }
      }
      if (counter < 1)
         counter = max.counter
      else if (counter > max.counter)
         counter = 1
   } #while
   graphics.off()
   if (c == 'x')
      q("no")
#  return(counter)
}

show.block.number=function(block.number,x,line=3)
{
   mtext(text=paste("(",block.number,")",sep=""),side=1,line=line,
          at=c(x,0),cex=0.8)
}

extract.value=function(text,text.list)
{
   name = paste("^",text,sep="")
#  print(paste(text,name))
   ndx = grep(name,text.list)[1]
#  print(paste("ndx",ndx))
#  print(paste("nchar",nchar(name)))
   colon = grep(":",name)
#  if (text=="K")
#  {
#     print(paste("ndx",ndx))
#     print(paste("colon",colon,length(colon),is.numeric(colon),as.numeric(colon),sep="|"))
#  }
   offset = 2
   if(length(colon)>0)
      offset = 1
#  print(paste("offset",offset))
#  print(text.list[ndx])
#  print(text.list[ndx+offset])
   x = as.numeric(text.list[ndx+offset])
#  print(paste(text,x))
   return(x)
}


read.rep=function(file="issams.rep",ntime=61,dt=1,ngear=5)

{
#  print(paste("Scanning file",file))
   rep = scan(file,what="character",quiet=TRUE)
   ret = list()
   ret$nll = extract.value("nll",rep)
   ret$nvar = extract.value("nvar",rep)
   ret$logr = extract.value("logr",rep)
   ret$r = extract.value("r",rep)
   ret$logK = extract.value("logK",rep)
   ret$K = extract.value("K",rep)
   ret$logsdlogF = extract.value("logsdlogF:",rep)
   ret$sdlogF = extract.value("sdlogF:",rep)
   ret$logsdlogPop = extract.value("logsdlogPop:",rep)
   ret$sdlogPop = extract.value("sdlogPop:",rep)
   ret$logsdlogYield = extract.value("logsdlogYield:",rep)
   ret$sdlogYield = extract.value("sdlogYield:",rep)
   ret$LQ = extract.value("LQ:",rep)
   ret$Q = extract.value("Q:",rep)
   ret$logsdlogQ = extract.value("logsdlogQ:",rep)
   ret$sdlogQ = extract.value("sdlogQ:",rep)
   ret$Lpcon = extract.value("Lpcon",rep)
   ret$pcon = extract.value("pcon",rep)
   ret$KF = extract.value("klingon_multiplier",rep)
   return(ret)
}

read.rep.files=function(path.list)# read.rep.files(dir(path=".",pattern="CV*"))
{
   print(path.list)
   have.names=FALSE
   csv = "rep_summary.csv"
   nvar = 0
   rep.list=list() 
   r = 0
   for (p in path.list)
   {
      path = paste(".",p,"issams.rep",sep="/")
   #  print(path)
      rep = read.rep(file=path)
      if (nvar <= 0)
         nvar = length(rep)
   #  print(rep)
      r = r+1
      rep.list = c(rep.list,rep)
      if (!have.names)
      {
         cat(paste(names(rep)[1],", ",sep=""),file=csv,append=FALSE)
         for (i in 2:nvar)
         {
            print(paste(i,names(rep)[i]))
         #  cat(names(rep)[i],sep=",",file=csv,append=TRUE)
            cat(names(rep)[i],file=csv,append=TRUE)
            if (i < nvar)
               cat(", ",file=csv,append=TRUE)
         }
         cat("\n",file=csv,append=TRUE)
         have.names=TRUE
      }

      for (i in 1:nvar)
      {
      #  print(paste(i,rep[i]))
      #  cat(paste(rep[i],", ",sep=""),file=csv,append=TRUE)
         cat(rep[[i]],file=csv,append=TRUE)
         if (i < nvar)
            cat(", ",file=csv,append=TRUE)
      }
      cat("\n",file=csv,append=TRUE)
   }
   
   rep.matrix = matrix(unlist(rep.list),nrow=length(path.list),byrow=TRUE)
   colnames(rep.matrix)=names(rep)
   return(as.data.frame(rep.matrix))

}

nice.xrange=function(x,nchar=3)
{
   xr = c(min(x,na.rm=TRUE),max(x,na.rm=TRUE))
   if (nchar > 0)
      xr[2] = xr[2]*1.05 #*nchar

   return(xr)
}

plot.nll.K=function(path.list)
{
   rep = read.rep.files(path.list)
   width = 9.0
   height = 6.0
   x11(width=width,height=height)
   par(mar=c(4,4,0,4)+0.1)
   xrange = nice.xrange(rep$sdlogF,5)
   nx = length(rep$sdlogF)
#  plot(rep$sdlogF,rep$nll,type='n',xlab="Process Error SD",ylab="Negative Log Likelihood")
   plot(xrange,range(rep$nll),type='n',xlab="Process Error SD",ylab="Negative Log Likelihood")
   double.lines(rep$sdlogF,rep$nll,bcol="blue",fcol="lightblue",lwd=5)
   points(rep$sdlogF,rep$nll,col="blue",pch=16)
   text(rep$sdlogF[nx],rep$nll[nx]," NLL",pos=4,offset=0.25,cex=1.0)

   par("new"=TRUE)
#  plot(rep$sdlogF,rep$K,type='n',axes=FALSE,ann=FALSE)
   plot(xrange,range(rep$K),type='n',axes=FALSE,ann=FALSE)
   double.lines(rep$sdlogF,rep$K,bcol="orange4",fcol="orange",lwd=5)
   points(rep$sdlogF,rep$K,col="orange4",pch=16)
   text(rep$sdlogF[nx],rep$K[nx]," K",pos=4,offset=0.25,cex=1.0)
   axis(4) # ,col="orange4",col.axis="orange4")
   mtext("K",side=4,line=3) #,col="orange")

   save.png.plot("NLL-K",width=width,height=height)
   

}

comp.nll.K.r=function()
{ 
   path.0.4 = "CVruns-r0.4"
   dirs.0.4 = dir(path=path.0.4,pattern="CV*")
#  print(dirs.0.4)
   rep.0.4 = read.rep.files(paste(path.0.4,dirs.0.4,sep="/"))
#  print(rep.0.4)

   path.est = "CVruns-rest"
   dirs.est = dir(path=path.est,pattern="CV*")
#  print(path.est)
   rep.est = read.rep.files(paste(path.est,dirs.est,sep="/"))
#  print(rep.est)

   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   par(mar=c(4,4,0,0)+0.1)
   lm = layout(matrix(c(1:3),ncol=1,byrow=TRUE))
   layout.show(lm)
   xrange = nice.xrange(rep.0.4$sdlogF,5)
   nx = length(rep.0.4$sdlogF)

   yrange = range(c(rep.0.4$nll,rep.est$nll))
   plot(xrange,yrange,type='n',xlab="",ylab="Negative Log Likelihood")
   double.lines(rep.0.4$sdlogF,rep.0.4$nll,bcol="blue",fcol="lightblue",lwd=5)
   points(rep.0.4$sdlogF,rep.0.4$nll,col="blue",pch=16)
   text(rep.0.4$sdlogF[1],rep.0.4$nll[1]," r=0.4",pos=4,offset=0.25,cex=1.0)
   double.lines(rep.est$sdlogF,rep.est$nll,bcol="darkgreen",fcol="lightgreen",lwd=5)
   points(rep.est$sdlogF,rep.est$nll,col="darkgreen",pch=16)
   text(rep.est$sdlogF[1],rep.est$nll[1]," est r",pos=4,offset=0.25,cex=1.0)

   yrange = range(c(rep.0.4$r,rep.est$r))
   plot(xrange,yrange,type='n',xlab="",ylab="Growth rate, r")
   double.lines(rep.0.4$sdlogF,rep.0.4$r,bcol="blue",fcol="lightblue",lwd=5)
   points(rep.0.4$sdlogF,rep.0.4$r,col="blue",pch=16)
   text(rep.0.4$sdlogF[nx],rep.0.4$r[nx]," r=0.4",pos=4,offset=0.25,cex=1.0)
   double.lines(rep.est$sdlogF,rep.est$r,bcol="darkgreen",fcol="lightgreen",lwd=5)
   points(rep.est$sdlogF,rep.est$r,col="darkgreen",pch=16)
   text(rep.est$sdlogF[nx],rep.est$r[nx]," est r",pos=4,offset=0.25,cex=1.0)
   
   yrange = range(c(rep.0.4$K,rep.est$K))
   plot(xrange,yrange,type='n',xlab="Process Error SD",ylab="Asymptotic Size, K")
   double.lines(rep.0.4$sdlogF,rep.0.4$K,bcol="blue",fcol="lightblue",lwd=5)
   points(rep.0.4$sdlogF,rep.0.4$K,col="blue",pch=16)
   text(rep.0.4$sdlogF[nx],rep.0.4$K[nx]," r=0.4",pos=4,offset=0.25,cex=1.0)
   double.lines(rep.est$sdlogF,rep.est$K,bcol="darkgreen",fcol="lightgreen",lwd=5)
   points(rep.est$sdlogF,rep.est$K,col="darkgreen",pch=16)
   text(rep.est$sdlogF[nx],rep.est$K[nx]," est r",pos=4,offset=0.25,cex=1.0)
   

   save.png.plot("nll-K-r-comp",width=width,height=height)
}


read.fit.files=function(path.list)
{
   nfit = length(path.list)
   print(nfit)
   print(path.list)
   csv = "std_summary.csv"
   nvar = 0
   fit.list=list() 
   all.pnames=c("ar", "aK", "asdlogF", "asdlogPop", "asdlogYield", "aQ", "asdlogQ", "apcon")
   pnames=c("ar", "aK", "asdlogYield", "aQ")

   var.names = c("mult", "nlogl", "r","s.r","cv.r", "K","s.K","cv.K",
                         "sY","s.sY","cv.sY", "Q","s.Q","cv.Q")
   nvar=length(var.names)
   print(var.names)
   print(nvar)
   est.mat = matrix(nrow=nfit,ncol=nvar)
   colnames(est.mat) = var.names 
#  for (p in path.list)
   for (j in 1:nfit)
   {
      p = path.list[j]
      path = paste(".",p,"issams",sep="/")
      fit = read.fit(file=path)
      pndx = which(match(fit$names,pnames,nomatch=0)>0)
      cv =fit$std[pndx]/fit$est[pndx]

      est.mat[j,1] = as.numeric(strsplit(p,"KF")[[1]][2])
      est.mat[j,2] = fit$nlogl
      c = 2
      for (v in 1:4)
      {
         c = c + 1
         est.mat[j,c] = fit$est[pndx[v]]
         c = c + 1
         est.mat[j,c] = fit$std[pndx[v]]
         c = c + 1
         est.mat[j,c] = cv[v]
      }
   }
   est.mat = as.data.frame(est.mat)
   return(est.mat[order(est.mat$mult),])

}

get.cor.table=function(file)
{
   fit = read.fit(file)
   fixed = which(fit$names!="U")
   cor.mat = fit$cor[fixed,fixed]
   colnames(cor.mat) = fit$names[fixed]
   rownames(cor.mat) = fit$names[fixed]
#  print(cor.mat)
   csv = paste(file,"-cor.csv",sep="")
   cat(paste("#",file,"\n",sep=" "),file=csv,append=FALSE)
   cat(paste("#","nll =",fit$nlogl,"\n",sep=" "),file=csv,append=TRUE)
#  cat(", ",file=csv,append=TRUE)
   cat("Name, Est., S.D., ",file=csv,append=TRUE)
   n = length(fixed)
   for (i in 1:n)
   {
      cat(colnames(cor.mat)[i],file=csv,append=TRUE)
   #  cat(paste(colnames(cor.mat)[i],fit$est[fixed[i]],fit$std[fixed[i]],
   #            sep=", "),file=csv,append=TRUE)
      if (i < n)
         cat(", ",file=csv,append=TRUE)
      else
         cat("\n",file=csv,append=TRUE)
   }       
   
   for (i in 1:n)
   {
   #  cat(paste(rownames(cor.mat)[i],", ",sep=""),file=csv,append=TRUE)
      cat(paste(rownames(cor.mat)[i],", ",fit$est[fixed[i]],", ",fit$std[fixed[i]],", ",
                sep=""),file=csv,append=TRUE)
      for (j in 1:n)
      {
      #  print(paste(i,j))
         if (j <= i)
            cat(cor.mat[i,j],file=csv,append=TRUE)
         else
            cor.mat[i,j] = NA

         if (j <= i)
            cat(", ",file=csv,append=TRUE)

      }
      cat("\n",file=csv,append=TRUE)
   }
#  print(cor.mat)
   return(fit)

}

plot.U=function(file,ntime=61,dt=1,ngear=5)
{
   fit = read.fit(file)
   uu = which(fit$names=="U")
   utPop = ngear*ntime
   est = fit$est[uu]
   std = fit$std[uu]
   tt = seq(1,ntime,dt)
   #  dat$t = (start.year-0.4*dt +  dat$t*dt)
   t = (start.year-0.4*dt +  tt*dt)
   P = exp(est[tt+utPop])
   sdP = std[tt+utPop]

   options(scipen=6)
   nice.ts.plot(t,P)
   plot.error(t,P,sdP,bcol="blue",fcol="lightblue",mult=1)
#  return(std[tt+utPop])
}


#con <- pipe(paste0("cut -f1,2,3,4,5 -d, ","/home/test/Test.csv"))
#x <- matrix(scan(con,skip=1,sep=","),ncol=5)
#close(con)
awkread = function(file,name) 
{
   tname = paste("#a",name,sep="")
   awk = paste("awk '(x)&&(NF==0){exit}",
         "(!x)&&(/^",tname,"/){x=1;next}",
         "(x){print}' ",file,sep="")
#  dat = as.data.frame(matrix(scan(pipe(awk),quiet=TRUE),byrow=TRUE,ncol=2))
   con = pipe(awk)
   dat = as.data.frame(matrix(scan(con,quiet=TRUE),byrow=TRUE,ncol=2))
   close(con)
   names(dat) <- c("x","p")
   return(dat)
}  


#  fit = read.fit(file)
hst.plot=function(model="issams")
{
#  hst.names=c("logMSY","MSY","logFmsy","Fmsy","r","K","sdlogProc",
#              "sdlogYield","logQ","Q")
#  print(length(hst.names))
   hst.names=c("logFmsy","logMSY","sdlogProc","sdlogYield","logQ",
               "Fmsy","MSY","K","r","Q")
#  print(length(hst.names))
   dat = scan(paste(model,".dat",sep=""),what="character")
#  print(dat)
   wp = which(dat == "r_prior")
   r.prior.use = as.numeric(dat[wp+4])
   r.prior = as.numeric(dat[wp+5])
   sd.r.prior = as.numeric(dat[wp+6])
   print(paste(wp,r.prior.use,r.prior,sd.r.prior))

   nvar = length(hst.names)
   nrc = ceiling(sqrt(nvar))
   print(paste(nvar,nrc))
#  lm = layout(matrix(c(1:(nrc*nrc)),ncol=nrc,nrow=nrc,byrow=TRUE))
#  print(lm)
#  layout.show(lm)
   xpos = 0
   ypos = 50
   incr = 50
   lwd = 2
   hst.file = paste(model,".hst",sep="")
   fit = read.fit(model)
#  print(fit)
   for (v in 1:nvar)
   {
      w = which(fit$names==paste("a",hst.names[v],sep=""))
   #  x1 = fit$est[w]-3.0*fit$std[w]
   #  x2 = fit$est[w]+3.0*fit$std[w]
      dist = awkread(hst.file,hst.names[v])
   #  dist$p = dist$p/sum(dist$p,na.rm=TRUE)
      x1 = min(dist$x,na.rm=TRUE)
      x2 = max(dist$x,na.rm=TRUE)
      x1 = fit$est[w]-5*fit$std[w]
      x2 = fit$est[w]+5*fit$std[w]
   #  print(" ")
      print(paste(v,fit$names[w],fit$est[w],fit$std[w],x1,x2))
      xx = seq(x1,x2,(x2-x1)/100.0)
      yy = dnorm(xx,mean=fit$est[w],sd=fit$std[w])
      pp = dist$p
   #  print(paste(sum(pp),sum(yy)))
   #  yy = yy/sum(yy,na.rm=TRUE)
   #  pp = pp/(sum(pp,na.rm=TRUE))
   #  print(paste(sum(pp),sum(yy)))
   #  print(cbind(xx,yy))
      if (nrow(dist) > 4)
      {
         x11(xpos=xpos,ypos=ypos,title=hst.names[v])
         xrange = c(max(min(xx),min(dist$x)),min(max(xx),max(dist$x))) #range(xx,dist$x)
         xrange = c(x1,x2) #range(dist$x)
   #     print(paste(range(pp),range(yy)))
         yrange = range(pp,yy) #range(yy,dist$p)
         ylab=paste("p(",hst.names[v],")",sep="") 
         plot(xrange,yrange,xlab=hst.names[v],ylab=ylab,type='n')
         lines(dist$x,pp,type='b',lwd=lwd)
         abline(v=fit$est[w],lty="dotdash",col="blue",lwd=lwd)
         abline(v=fit$est[w]-2.0*fit$std[w],lty="dotted",col="blue",lwd=lwd)
         abline(v=fit$est[w]+2.0*fit$std[w],lty="dotted",col="blue",lwd=lwd)
         lines(xx,yy,col="blue",lwd=lwd)
         if ((hst.names[v] == "r" ) && r.prior.use > 0)
         {
            pyy = dnorm(xx,mean=r.prior,sd=sd.r.prior)
            lines(xx,pyy,col="red",lwd=lwd)
            abline(v=r.prior,lty="dotdash",col="red",lwd=lwd)
         }
         #   ar = 2.0*mfexp(logFmsy);
         if ((hst.names[v] == "logFmsy" ) && r.prior.use > 0)
         {
            logFmsy.prior = log(r.prior)-log(2)
            pyy = dnorm(xx,mean=logFmsy.prior,sd=sd.r.prior)
            lines(xx,pyy,col="red",lwd=lwd)
            abline(v=logFmsy.prior,lty="dotdash",col="red")
         }
         if ((hst.names[v] == "Fmsy" ) && r.prior.use > 0)
         {
            Fmsy.prior = r.prior/2.0
            sd.Fmsy.prior = sd.r.prior/2.0
            pyy = dnorm(xx,mean=Fmsy.prior,sd=sd.Fmsy.prior)
            lines(xx,pyy,col="red",lwd=lwd)
            abline(v=Fmsy.prior,lty="dotdash",col="red")
         }
         #   aK = 4.0*mfexp(logMSY)/(1.0e-20+ar);
         if ((hst.names[v] == "K" ) && r.prior.use > 0)
         {
            wMSY = which(fit$names=="aMSY")
            MSY = fit$est[wMSY]
            K.prior = 4.0*MSY/r.prior
            sd.K.prior = 4.0*MSY/r.prior
            pyy = dnorm(xx,mean=K.prior,sd=sd.K.prior)
            lines(xx,pyy,col="red",lwd=lwd)
            abline(v=K.prior,lty="dotdash",col="red")
         }

         xpos = xpos + incr
         ypos = ypos + incr
      }
   }

}

# logit(N1/(N1+N2)) = log(N1)-log(N2)
logit<-function(p)
{
   return(log(p/(1-p)))
}

alogit=function(alpha)
{
   return(1/(1+exp(-alpha)))
}


plot.Qprop.prior=function(propL=0.01)
{
   p = seq(0.01,0.99,0.01)
   p = c(0.0001,p,0.9999)
   wp = which(p == propL)
   
   plot(c(0.0,1.0),c(0.0,1.0),type='n',xlab="x",ylab="p(L(x))")
   abline(v=propL,col="red",lty="dotdash")
   for (sd in c(0.1,0.2,0.3,0.4,5,0.6,0.7,0.8,0.9,0.95,0.99,0.999,0.9999))
#  sd = 0.4
   {
      pp = dnorm(logit(p),mean=logit(propL),sd=logit(sd))
   #  print(paste("sd",sd,logit(sd)))
   #  print(pp)
      double.lines(p,pp,lwd=5,fcol="lightblue",bcol="blue")
      text(propL,pp[wp],sd,col="blue")
   }

   save.png.plot("Qprop_prior",width=6.5,height=6.5)
}

TMB.diagnostics=function(data,obj,opt,
                         plot.Fmort=TRUE,plot.prod=TRUE,plot.impact=FALSE)
{
   dt = data$dt
   ntime = data$ntime
   ngear = data$ngear
   start.year = 1952
   yy = seq(1,ntime,1)
   yy = (start.year-0.4*dt +  yy*dt)
#  print(yy)
   resid = as.data.frame(cbind(yy,obj$report()$resid))
   resid.names =c("t","pop","K","forcing","F1","F2","F3","F4","F5","predC1","predC2","predC3","predC4","predC5","obsC1","obsC2","obsC3","obsC4","obsC5")
   colnames(resid)=resid.names
   print(head(resid))


   ests = opt$par
#  exp(junk[which(names(junk)=="logMSY")])
   MSY = exp(ests[which(names(ests)=="logMSY")])
   Fmsy = exp(ests[which(names(ests)=="logFmsy")])
   print(paste("TMB.diagnostics",MSY,Fmsy))
   r = 2.0*Fmsy
   K = 4.0*MSY/(1.0e-20+r)
   print(paste("TMB.diagnostics",r,K))
   B1 = NULL
   d = NULL
   sdlogProc = exp(ests[which(names(ests)=="logsdlogProc")])
   sdlogYield = exp(ests[which(names(ests)=="logsdlogYield")])

   max.counter = 1
   counter = max.counter
   c = 'n'
   dev.list = vector(mode="numeric",length=5)
#  dev.file.names=c("tmp/est_pop","tmp/est_catch","tmp/est_F","tmp/prod")
   dev.file.names=c("est_pop","est_catch","est_F","prod")
#  while (c != 'q')
   while ( (c != 'q') && (c != 'x') )
   {
   #  tmp=get.diagnostics(log,ntime=ntime,dt=dt,ngear=ngear,block=counter,mtype="i")
G
      new.devices = plot.diagnostics(resid,dt=dt,ngear=ngear,
                    sdlogPop=sdlogProc,
                    sdlogYield=sdlogYield,
                    sdlogF=sdlogProc,
                    sdlogQ=sdlogProc,
                    K=K, r=r, B1=B1, indexed=data$use_Q,
                    plot.Fmort=plot.Fmort,
                    plot.prod=plot.prod,
                    plot.impact=plot.impact,
                    devices=dev.list,block=NULL)


      dev.list=new.devices
     
      c = readline("next, back, save, quit or exit? (enter n,b,s,q,x):")
      print(paste(c," entered"))
      if (c == 'n')
         counter = counter + 1
      else if (c == 'b')
         counter = counter - 1
      else if (c == 's')
      {
         print(dev.list)
         for (d in 1:length(dev.list))
         {
            if (dev.list[d] > 0)
            {
               tname = paste(dev.file.names[d],"B",counter,sep="")
               print(paste("Saving device",d,"to file",tname))
               dev.set(dev.list[d])
               din = par("din")
               print(din)
               save.png.plot(tname,width=din[1],height=din[2])
            }
            else
               print(paste("Skipping device",d))
          }
      }
      if (counter < 1)
         counter = max.counter
      else if (counter > max.counter)
         counter = 1
   } #while
   graphics.off()
   if (c == 'x')
      q("no")
#  return(counter)
}


TMB.sd.comp=function(TMB.path="/home/jsibert/Projects/xssa/TMB/issams-fit.Rdata",
                    ADMB.path="/home/jsibert/Projects/xssa/run-issams/use_r_prior/r2/Q1/issams")
{
   varnames = c("logFmsy","logMSY","logsdlogProc","logsdlogYield",
                "logQ","sdlogProc","sdlogYield","Q",
                "MSY","Fmsy","r","K")
#  print(varnames)
   nnames = length(varnames)

   TMB.fit = load(TMB.path)
   TMB.std = std
#  print(TMB.fit)

   ADMB.std.file = paste(ADMB.path,".std",sep="")
   ADMB.std = matrix(ncol=2,nrow=nnames)

   rownames(ADMB.std) = varnames
   colnames(ADMB.std) = c("Estimate","Std. Error")
   ttt=try(scan(file=ADMB.std.file,what="character"),TRUE)
#  print(paste("class(std))",class(std)))
   if (class(std) == "try-error")
   {
      print(paste(ADMB.std.file,"missing"))
   }
   else
   {
      for (i in 1:nnames)
      {
         w = which(ttt == varnames[i])

         if (length(w) < 1)
         {
            tname = paste("a",varnames[i],sep="")
            w = which(ttt == tname)
         #  print(paste(i,tname,w,length(w),ttt[w+1],ttt[w+2]))
       
         }
         ADMB.std[i,1] = as.numeric(ttt[w+1])
         ADMB.std[i,2] = as.numeric(ttt[w+2])
      }  
   #  print(ADMB.std)
   }

   fit.list = list(ADMB=ADMB.std,TMB=TMB.std)
#  return(fit.list)

   get.sd.range=function(var,fits)
   {
      for (v in var)
      {
         smax = (-Inf)
         smin =   Inf 
      #  print(paste(smin,smax))
         for (f in 1:length(fits))
         {
            ff = as.data.frame(fits[[f]])
            smin = min(smin,ff$Estimate[v]-ff$"Std. Error"[v])
            smax = max(smax,ff$Estimate[v]+ff$"Std. Error"[v])
         }
      }
      return(range(smin,smax))
   } 
#  print(get.sd.range(1,fit.list))
#  print(get.sd.range(c(3,4),fit.list))
#  print(get.sd.range(12,fit.list))

   plot.sd.range=function(var,fits)
   { 
      yinc = c(1.5/5,3.5/5)
      vcount = 0
      cols=c("blue","red")
      for (v in var)
      {
         xrange = get.sd.range(v,fits)
         print(varnames[var])
         print(xrange)
         plot(xrange,c(0,1),type='n',ylab="",axes=FALSE,ann=FALSE,col="black")
         for (f in 1:length(fits))
         {
            vcount = vcount + 1
            vndx   = 1 + vcount%%2
         #  print(paste(vcount,vndx))
            ff = as.data.frame(fits[[f]])
            x = c(ff$Estimate[v]-ff$"Std. Error"[v],
                  ff$Estimate[v]+ff$"Std. Error"[v])
            y = c(yinc[vndx],yinc[vndx])
            print(x)
         #  print(y)
            lines(x,y,lwd=3,col=cols[f])
            points(ff$Estimate[v],y[1],pch=16,col=cols[f],cex=2)         
            title(main=varnames[var],line=-0.5,cex.main=0.9)
            axis(1)
         }
      }
   }   


   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   par(mar=c(2,2,0,0)+0.1)
   lm = layout(matrix(c(1:nnames),ncol=1,byrow=TRUE))
   layout.show(lm)

   plot.sd.range(1,fit.list)
   plot.sd.range(2,fit.list)
   plot.sd.range(3,fit.list)
   plot.sd.range(4,fit.list)
   plot.sd.range(5,fit.list)
   plot.sd.range(6,fit.list)
   plot.sd.range(7,fit.list)
   plot.sd.range(8,fit.list)
   plot.sd.range(9,fit.list)
   plot.sd.range(10,fit.list)
   plot.sd.range(11,fit.list)
   plot.sd.range(12,fit.list)

   save.png.plot("ADMB-TMB-sd-comp",width=width,height=height)
}
