

#####constant fitter

constantobjective = function(N,A){
  value = N*(log(N/A))-N
  value[is.nan(value)]<-0
  return(value)
}

#######linear fitter
linfit <-function(N,t){
  
  if (length(t)==0){ # if no data is recieved, set cost to -inf
    return(list("a"=0,"b"=0,"cost"=-Inf))
  }  
  if (length(t)==1){ #if only one point is supplied, treat as constant intensity
    return(list("a"=0,"b"=N/t[1],"cost"=log(N/t[1])-N))
  }
  if (length(t)==2){ #if only two points are supplied, fit it perfectly
    int1=1/t[1]
    int2=1/(t[2]-t[1])
    return(list("a"=(int2-int1)/t[2],"b"=int1,"cost"=N[1]*log(N[1]/t[1])+N[2]*log(N[2]/(t[2]-t[1]))-sum(N)))
  }
  
  t=rep(t,N)
  
  epsilon  <- 1
  n        <- length(t)
  i        <- 1
  M        <- t[n]
  S        <- M^2/2
  
  #start with some initial values for a and b
  coef     <- matrix(c(0,length(t)/M),nrow=2)
  
  new.cost <- -Inf
  
  while (abs(epsilon) >10e-4 && i<100){
    old.cost <- new.cost
    
    #first derivatives
    f=matrix(c(sum(t/(coef[1]*t+coef[2]))-S,sum(1/(coef[1]*t+coef[2]))-M),nrow=2)
    
    #hessian values
    fa <- -sum((t/(coef[1]*t+coef[2]))^2)
    fb <- -sum(t/(coef[1]*t+coef[2])^2)
    ga <- fb
    gb <- -sum(1/(coef[1]*t+coef[2])^2)
    
    #create the inverse of the hessian
    invhess <- 1/(fa*gb-fb*ga)*matrix(c(gb,-fb,-ga,fa),nrow=2,ncol=2,byrow=TRUE) 
    
    #run newtons method
    coef=coef-invhess%*%f 
    
    if(coef[2]<0 || (M*coef[1]+coef[2])<0){#if best line has a negative intensity, dont use linear
      
      return(list("a"=0,"b"=length(t)/t[length(t)],"cost"=constantobjective(length(t),t[length(t)])))
      
    }
    #calculate the new cost
    logsum=coef[1]*t+coef[2]
    logsum=logsum[logsum>0]
    new.cost <- sum(log(logsum))-coef[1]*S-coef[2]*M
    epsilon  <- new.cost-old.cost
    
    i <- i+1
  }
  list("a"=coef[1],"b"=coef[2],"cost"=new.cost)
}

"plot.BB" <- function(x,show=c("hist","blocks"),binwidth=NULL,bins=NULL,ylim=NULL,xlim=NULL,xact=NULL,main="Bayesian Blocks",xlab="Time",ylab="Intensity") {
  

  data       <- x$data
  n          <- length(data)
  intensity  <- x$N/x$A
  legend     <- vector()
  legend.col <- vector()

  #check if xlim is supplied
  if (!is.null(xlim)){
    lowerbound=xlim[1]
    upperbound=xlim[2]
  }
  else{
    lowerbound <- data[1]-(data[2]-data[1])/2
    upperbound <- data[n]+(data[n]-data[n-1])/2
  }
  
  #check if ylim is supplied
  if(is.null(ylim)){
    ylim=c(0,max(intensity)*0.75)
  }
  
  #check if the hist bins are supplied
  if (is.null(binwidth) & is.null(bins)){
    binwidth=data[n]/100
    bins=100
  }
  else if (is.null(bins)){
    bins   <- round(data[n]/binwidth,0)
  }
  else{
    binwidth=data[n]/bins
  }
  
  xx=c(lowerbound,(data[1:(length(data)-1)] + data[2:(length(data))])/2,upperbound)
  intensity=1/diff(xx)
  
  #initialize plot
  plot(NA, ylim = ylim, xlim = c(lowerbound,upperbound),xaxt=xact,bty="n",
       main=main,xlab=xlab,ylab=ylab)
  grid(NA,NULL)#add grid for y-axis
  
  
  
  #plot histogram if requested
  if ("hist" %in% show){
    histdata <- rep(data,x$N)
    end      <- histdata[length(histdata)]
    i        <- 0:(end/binwidth)
    
    height   <- hist(histdata, breaks=c(i,bins+1)*binwidth,plot=FALSE)$counts/binwidth
    rect(xleft=i*binwidth,ybottom=0,xright= (i+1)*binwidth, ytop=height,col="grey70",border="grey40")
    legend=c(legend,"Binned Data")
    legend.col=c(legend.col,"grey70")
  }
  
  
  #plot individual intensities if requested
  if ("points" %in% show){
    segments(x0=c(xx[1: (length(xx)-1)]), y0=intensity,
             x1=xx[2:length(xx)], y1=intensity,
             col = "cornflowerblue", lwd = 3)
    
    legend=c(legend,"Data Points")
    legend.col=c(legend.col,"cornflowerblue")
  }
  
  
  #plot blocks if requested
  if ("blocks" %in% show){
    
    cpoints=c(x$left,length(x$data)+1)
      
    
    x0=NULL
    y0=NULL
    for(i in 1:length(x$type)){
      x0=c(x0,xx[cpoints[i]])
      y0=c(y0,x$params[[i]]$b)
      x0=c(x0,xx[cpoints[i+1]])
      y0=c(y0,x$params[[i]]$b+ifelse(is.null(x$params[[i]]$a),0,x$params[[i]]$a*(xx[cpoints[i+1]]-xx[cpoints[i]])))
    }
    lines(x0,y0,lwd=2,col="red")
    
    legend=c(legend,"Bayesian Blocks")
    legend.col=c(legend.col,"red")
    
  }
  legend=c(legend,"True Intensity")
  legend.col=c(legend.col,"black")

  
  #plot legend
  legend(x=x$data[1],y=ylim[2]*0.98,legend= legend, lty=c(rep(1,length(legend)-1),2),
         lwd=rep(2,length(legend)), col=legend.col)
}

#S3 generic plot function 


# OPTINTERVAL FUNCTION 
# An O(N^2) algorithm for finding the optimal partition of N data points on an interval
# INPUT:
# data is the vector ofcells
# N is the number of "observations" per point. This is the same "N" we've been using all along
# c is the block penalty term
# OUTPUT:
# a list of objects pertaining to the "BB" object

optinterval = function(data,N,c,type=c("constant"),alpha=0.05,pen="Likelihood Ratio Test",verbose=FALSE){
  start.time <- Sys.time()
  
  n          <- length(N)
  percent    <- floor(n/100)
  lowerbound <- data[1]-(data[2]-data[1])/2
  upperbound <- data[n]+(data[n]-data[n-1])/2
  xx         <- c(lowerbound,(data[1:(length(data)-1)] + data[2:(length(data))])/2,upperbound) # voronoi cell vertices
  A          <- diff(xx) # length of each voronoi cell
  
  if(length(type)==1 && type=="linear"){
    alpha=1
  }
  
  if (pen=="Likelihood Ratio Test"){
    chi.lin    <- qchisq(1-alpha,df=1)/2
  
  }
  else if(pen=="AIC"){
    c=c+1
    chi.lin    <- 1
  }
  else if(pen=="BIC"){
    c=c+0.5*log(n)
    chi.lin    <- log(n)
  }
  

  opt           <- rep(0,n+1)
  lastchange    <- rep(1,n)
  changeA       <- rep(0,n) 
  changeN       <- rep(0,n)
  optint        <- matrix(-Inf,nrow=n,ncol=length(type))
  last          <- rep(0,length(type))
  lasttype      <- rep("None",n)
  lastparams    <- list()
  unpruned <- NULL
  endobj   <- rep(0,n)
  
  
  
  #begin looping through each point
  for (i in 1:n){
    
    unpruned          <- c(unpruned,i)
    
    ##### constant blocks ##########
    if ("constant" %in% type){
      changeA[unpruned]  <- changeA[unpruned] + A[i]
      changeN[unpruned]  <- changeN[unpruned] + N[i]
      optint[unpruned,which(type=="constant")]  <- opt[unpruned] + constantobjective(changeN[unpruned],changeA[unpruned])-c
      last[which(type=="constant")]             <- which.max(optint[unpruned,which(type=="constant")])
    }
    ################################
    
    ##### linear blocks ############
    if ("linear" %in% type){
      linblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        linblocks[[j]]   <- linfit(N[j:i],x)
        optint[j,which(type=="linear")]       <- opt[j] + linblocks[[j]][["cost"]]-c-chi.lin
      }
      last[which(type=="linear")]              <- which.max(optint[unpruned,which(type=="linear")])
    }
    
    ################################
    
    bestshape<-which.max(apply(optint,2,max))


    lastchange[i]     <- which.max(optint[,bestshape])-1

    
    if(("constant" %in% type) && (bestshape==which(type=="constant"))){ #constant block is best
      lasttype[i]       <- "constant"
      lastparams[[i]]   <- list("b"    = sum(N[(lastchange[i]+1):(i)])/sum(A[(lastchange[i]+1):(i)]),
                                "cost" = constantobjective(sum(N[(lastchange[i]+1):(i+1)]),sum(A[(lastchange[i]+1):(i+1)])))
    }
    if (("linear" %in% type) && (bestshape==which(type=="linear"))){#linear block is best
      lasttype[i]       <- "linear"
      lastparams[[i]]   <- linblocks[[unpruned[which(type=="linear")]]]
    }
    
    
    opt[i+1]          <- max(optint[unpruned,bestshape])
    unpruned          <- unpruned[((optint[unpruned,bestshape]+c-opt[i+1])>0)]
    
    if((verbose==TRUE) && (i %% percent==0)){#print out the progress of the algorithm
      cat("\n",round(100*i/n,0),"% of points completed",round(100*(i-length(unpruned))/i,0),"% pruned ")
    }
    
  }
  BBdata   <- list("data"         = data, 
                   "N"            = N,
                   "A"            = A,
                   "opt"          = opt,
                   "lastchange"   = lastchange,
                   "lastparams"   = lastparams,
                   "lasttype"     = lasttype,
                   "pruned"       = setdiff(1:i,unpruned))
  
  
  
  print(Sys.time()-start.time)
  model <- getStep(BBdata,n)
  summary(model)
  return(model)
}

#S3 Generic summary function
summary.BB <- function(x) {
  cat(length(x$data), 'data points\n\n',
      length(x$left),'total blocks:\n\t',
      sum(x$type=='constant'),'constant blocks\n\t',
      sum(x$type=='linear'),'linear blocks\n\t',
      floor((length(x$pruned)*100)/length(x$data)),'% of the points were pruned\n\n')
}


#this function returns a BB object for the optimal blocks at any given cell location "n" where 1 < n < N
getStep <- function(BBobj,n){
  
  lasttype      <- BBobj$lasttype[1:n]
  lastchange    <- BBobj$lastchange[1:n]
  lastparams    <- BBobj$lastparams[1:n]
  
  
  left          <- vector()
  params        <- list()
  type          <- vector()
  type[1]       <- lasttype[n]
  left[1]       <- lastchange[n]
  params[[1]]   <- lastparams[[n]]
  i=1
  while (left[i] > 1) {
    
    left[i+1]      <- lastchange[left[i]]
    type[i+1]      <- lasttype[left[i]]
    params[[i+1]]  <- lastparams[[ left[i] ]]
    
    i <- i+1
    
  }
  
  left      <- rev(left)+1
  type      <- rev(type)
  params    <- rev(params)
  
  if (length(left) == 1){
    right <- n
  }
  else {
    right <- c(left[2:length(left)]-1,n)
  }
  
  BBdata   <- list("data"         = BBobj$data[1:n],
                   "N"            = BBobj$N[1:n],
                   "A"            = BBobj$A[1:n],
                   "left"         = left,
                   "right"        = right,
                   "type"         = type,
                   "params"       = params,
                   "opt"          = BBobj$opt[1:n],
                   "lastchange"   = lastchange,
                   "lastparams"   = lastparams,
                   "lasttype"     = lasttype,
                   "pruned"       = BBobj$pruned)
  
  BBobject <- structure(BBdata, class = "BB")
  
  return(BBobject)
}


#Non-homogeneous Poisson Process Generator!!
NHPois <- function(time_length, equation){
  
  expr<-function(x){
    x<-x
    eval(equation)
  }
  current <- 0
  i<-1
  data <- vector()
  maxrate = optim(par=time_length/2,fn=expr,method='L-BFGS-B',lower=0,upper=time_length,control = list(fnscale = -1))$value
  while(current < time_length){
    current<- current+ rexp(n=1,rate=maxrate)
    u <- runif(1)
    
    if (u < expr(current) / ( maxrate)){
      data[i] <- current
      i<-i+1
    }
    
  }
  return(data[c(-length(data))])
}



#### create some simulated data
set.seed(10)
equation1 <- expression(75)
equation2 <- expression(300-120*x)
equation3 <- expression(60)
x1 <- NHPois(time_length =2,equation = equation1)
x2 <- NHPois(time_length =2,equation = equation2)+x1[length(x1)]
x3 <- NHPois(time_length =1,equation = equation3)+x2[length(x2)]
x=c(x1,x2,x3)
N <- rep(1,length(x))#let all cells be of size 1
N<-as.data.frame(table(x))$Freq





##################INTERACTIVITY CODE 

#server function for the Shiny app
function(input, output,session) {

 runmodel <- reactive({model <- optinterval(x, N,input$c,type=rev(tolower(input$type)),alpha=input$alpha,input$c_choice,verbose=FALSE)
})
  # Fill in the spot we created for a plot
  output$plot <- renderPlot({
    #create model
        plot.BB(runmodel(),show=input$show,binwidth=input$binwidth/1000,xlab="Seconds",ylim=c(0,600))
      if(input$true) { 
        lines(x=c(0,2),y=c(75,75),lwd=2,lty='dashed')
        lines(x=c(2,2),y=c(75,300),lwd=2,lty='dashed')
        lines(x=c(2,4),y=c(300,60),lwd=2,lty='dashed')
        lines(x=c(4,5),y=c(60,60),lwd=2,lty='dashed')
      }
    # Render a barplot
  },height = 700)
  
  output$constants <- renderText({ 
    paste0("Constant blocks: ",sum(runmodel()$type=='constant'))
  })
  output$linears <- renderText({
    paste0("Linear blocks: ",sum(runmodel()$type=='linear'))
  })
  output$pruned <- renderText({
    paste0( floor((length(runmodel()$pruned)*100)/length(runmodel()$data)),"% pruned")
  })
  
}

