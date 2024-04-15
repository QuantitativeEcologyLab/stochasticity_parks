## attempt at a beta gamlss family, partiallyy using auto-generated derivative pipeline, for
## michael.noonan@ubc.ca, stefano.mezzini@ubc.ca, rekhamarcus@gmail.com
## (c) Simon N. Wood 2023

betals <- function(link=list("identity","identity"),eps=.Machine$double.eps*1000) {
## General family for beta location scale model...
## First linear predictor gives logit of mean, mu.
## Second linear predictor gives logiy of phi, where var(y) = phi(1-mu)*mu
## For beta, phi is bounded between 0 and 1.
## The usual beta distribution parameters are
## a = mu(1/phi-1) and b = (1-mu)(1/phi-1).
## NOTE: null deviance STILL MISSING
## 1. get derivatives wrt mu, rho and xi.
## 2. get required link derivatives and tri indices.
## 3. transform derivs to derivs wrt eta (gamlss.etamu).
## 4. get the grad and Hessian etc for the model
##    via a call to gamlss.gH  
  
  ## first deal with links and their derivatives...
  if (length(link)!=2) stop("betals requires 2 links specified as character strings")
  okLinks <- list("identity", "identity")
  stats <- list()
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- make.link(link[[i]]) else 
    stop(link[[i]]," link not available for betals")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
           mu.eta=stats[[i]]$mu.eta),
           class="family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }

  eval(parse(text=paste("preinitialize <- function(G) { eps <-",
  eps,";y <- G$y;y[y >= 1-eps] <- 1 - eps;y[y<= eps] <- eps; return(list(y=y));}")))
      
  

  eval(parse(text=paste("residuals <- function(object,type=c(\"deviance\",\"pearson\",\"response\")) {\n",
      "mu <- object$fitted[,1]; phi <- object$fitted[,2];\n",
      "y <- object$y; type <- match.arg(type)\n",
      "if (type==\"deviance\") { a <- mu*(1/phi-1); b <- (1-mu)*(1/phi-1);\n",
      "mode <- rep(0.5,length(a)); ii <- a>1 & b>1;\n mode[ii] <- (a[ii]-1)/(a[ii]+b[ii]-2);\n",
      "mode[a<=1 & b>1] <-", eps,"; mode[b<=1 & a>1] <-",1-eps,";\n",
      "rsd <- 2*(dbeta(mode,a,b,log=TRUE) - dbeta(y,a,b,log=TRUE))*sign(y-mu)\n",
      "} else if (type==\"pearson\") {\n",
      "sd <- sqrt(phi*mu*(1-mu));rsd <- (y-mu)/sd\n",
      "} else { rsd <- y-mu };\n",
      "rsd } ## betals residuals")))

   predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                beta=NULL,off=NULL,Vb=NULL) {
  ## optional function to give predicted values - idea is that 
  ## predict.gam(...,type="response") will use this, and that
  ## either eta will be provided, or {X, beta, off, Vb}. family$data
  ## contains any family specific extra information. 
  ## if se = FALSE returns one item list containing matrix otherwise 
  ## list of two matrices "fit" and "se.fit"... 

    if (is.null(eta)) {
      discrete <- is.list(X) 
      lpi <- attr(X,"lpi") 
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      }
      nobs <- if (discrete) nrow(X$kd) else nrow(X)
      eta <- matrix(0,nobs,2)
      ve <- matrix(0,nobs,2) ## variance of eta 
      for (i in 1:2) {
        if (discrete) {
	  eta[,i] <- Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]])
        } else {
          Xi <- X[,lpi[[i]],drop=FALSE]
          eta[,i] <- Xi%*%beta[lpi[[i]]] ## ith linear predictor
        } 
        if (!is.null(off[[i]])) eta[,i] <- eta[,i] + off[[i]]
        if (se) ve[,i] <- if (discrete) diagXVXd(X$Xd,Vb,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,nthreads=1,lt=X$lpid[[i]]) else
	                  drop(pmax(0,rowSums((Xi%*%Vb[lpi[[i]],lpi[[i]]])*Xi)))
      }
    } else { 
      se <- FALSE
    }
    bino <- binomial()
    bino$linkinv -> logis
    gamma <- cbind(logis(eta[,1]),logis(eta[,2]))
   
    if (se) { ## need to loop to find se of probabilities...
      mu.eta <- bino$mu.eta
      vp <- gamma
      vp[,1] <- sqrt(ve[,1])*abs(mu.eta(eta[,1]))
      vp[,2] <- sqrt(ve[,2])*abs(mu.eta(eta[,2]))
      return(list(fit=gamma,se.fit=vp))
    } ## if se
    list(fit=gamma)
  } ## betals predict


  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    object$fitted <- binomial()$linkinv(object$fitted) ## put on mean and phi scale
    object$null.deviance <- NA
    
  })

  ncv <- function(X,y,wt,nei,beta,family,llf,H=NULL,Hi=NULL,R=NULL,offset=NULL,dH=NULL,db=NULL,deriv=FALSE,nt=1) {
    gamlss.ncv(X,y,wt,nei,beta,family,llf,H=H,Hi=Hi,R=R,offset=offset,dH=dH,db=db,deriv=deriv,nt=nt)
  } ## ncv 

  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
  ## function defining the gamlss beta model log lik. 
  ## deriv: 0 - eval
  ##        1 - grad and Hess
  ##        2 - diagonal of first deriv of Hess
  ##        3 - first deriv of Hess
  ##        4 - everything.
    jj <- attr(X,"lpi") ## extract linear predictor index
    eta <- X[,jj[[1]],drop=FALSE]%*%coef[jj[[1]]]
    m <- family$linfo[[1]]$linkinv(eta) ## mean
    etav <- X[,jj[[2]],drop=FALSE]%*%coef[jj[[2]]] ## log sigma
    v <- family$linfo[[2]]$linkinv(etav) ## log sigma
   
    n <- length(y)
    l1 <- matrix(0,n,2)
    
    ## note that the derivative code is largely auto-generated, and
    ## auto-simplified. Optimized Maxima derivs exported as Maxima
    ## code, translated to R code in R, then processed in R to
    ## remove redundant auxiliary variables and their definitions.
    gamma3 <- function(x) psigamma(x,deriv=2);
    gamma4 <- function(x) psigamma(x,deriv=3);
    em1 <- exp(m)+1
    emv <- exp(m-v)
    aa1 <- 1/em1;
    aa2 <- aa1*emv;
    aa3 <- exp(-v);
    aa4 <- aa1*aa3;
    logy <- log(y);log1y <- log(1-y)
    l0 <-  (aa2-1)*logy+(aa4-1)*log1y - lgamma(aa4) + lgamma(aa3) - lgamma(aa2);
    l <- sum(l0)
    if (deriv>0) {
      ## first derivatives m, v...
      bb1 <- em1;
      bb2 <- 1/bb1;
      bb3 <- -v;
      bb4 <- exp(bb3+m);
      bb5 <- bb2*bb4;
      bb6 <- 1/bb1^2;
      bb7 <- bb5-bb6*exp(bb3+2*m);
      dgabb5 <- -digamma(bb5)
      l1[,1]  <-  bb7*logy-bb6*bb4*log1y + bb6*bb4*digamma(bb2/exp(v))+bb7*dgabb5;
      cc2 <- emv;
      dgaaa3 <- -digamma(aa3)
      l1[,2]  <-  (-aa1*cc2*logy)-aa1*aa3*log1y + aa1*aa3*digamma(aa1*aa3)+aa3*dgaaa3 +
                   aa1*cc2*digamma(aa1*cc2);

      ## the second derivatives mm mv vv
    
      l2 <- matrix(0,n,3)
      dd07 <- 2*m;
      dd08 <- exp(bb3+dd07);
      dd09 <- 1/bb1^3;
      dd10 <- 2*dd09*exp(bb3+3*m)-3*bb6*dd08+bb5;
      dd11 <- bb2/exp(v);
      dd12 <- -digamma(dd11);
      dd13 <- log1y;
      tgadd11 <- -trigamma(dd11)
      tgabb5 <- -trigamma(bb5)
      l2[,1]  <-  dd10*logy+2*dd09*dd08*dd13-bb6*bb4*dd13+(exp(dd07-2*v)*tgadd11)/bb1^4+
                  2*dd09*dd08*dd12-bb6*bb4*dd12+(bb5-bb6*dd08)^2*tgabb5+dd10*dgabb5;
      ee6 <- exp(bb3+2*m);
      ee7 <- bb6*ee6-bb2*bb4;
      l2[,2]  <-  ee7*logy+bb6*bb4*log1y+(exp(m-2*v)*tgadd11)/bb1^3+bb6*bb4*dd12-
                  bb2*bb4*(bb5-bb6*ee6)*tgabb5+ee7*dgabb5;
      ff4 <- bb2*cc2;
      ff7 <- exp(-2*v);
      ff8 <- bb2*aa3;
      tgaff8 <- -trigamma(ff8)
      dgaff8 <- -digamma(ff8)
      tgaaa3 <- -trigamma(aa3)
      tgaff4 <- -trigamma(ff4)
      dgaff4 <- -digamma(ff4)
      l2[,3]  <-  bb2*cc2*logy+bb2*aa3*log1y+bb6*ff7*tgaff8+bb2*aa3*dgaff8-
                  ff7*tgaaa3-aa3*dgaaa3+bb6*exp(2*m-2*v)*tgaff4+bb2*cc2*dgaff4;
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),family$linfo[[2]]$mu.eta(etav))
      g2 <- cbind(family$linfo[[1]]$d2link(m),family$linfo[[2]]$d2link(v))
    } ## if (deriv>0)
 
    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      ## the third derivatives
      ## order mmm mmv mvv vvv 
      l3 <- matrix(0,n,4)
      gg10 <- 3*m;
      gg11 <- exp(bb3+gg10);
      gg12 <- 1/bb1^4;
      gg13 <- (-6*gg12*exp(bb3+4*m))+12*dd09*gg11-7*bb6*dd08+bb5;
      gg14 <- bb5-bb6*dd08;
      gg17 <- -2*v;
      gg18 <- tgadd11;
      g3add11 <- -gamma3(dd11)
      g3abb5 <- -gamma3(bb5)
      l3[,1]  <-  gg13*logy-6*gg12*gg11*dd13+6*dd09*dd08*dd13-bb6*bb4*dd13-
              (exp(gg10-3*v)*g3add11)/bb1^6-(6*exp(gg17+gg10)*gg18)/bb1^5+
	      3*gg12*exp(gg17+dd07)*gg18-6*gg12*gg11*dd12+6*dd09*dd08*dd12-
	      bb6*bb4*dd12+gg14^3*g3abb5+3*gg14*(2*dd09*gg11-3*bb6*dd08+bb5)*tgabb5+gg13*dgabb5;
      hh05 <- -bb2*bb4;
      hh10 <- exp(bb3+3*m);
      hh11 <- (-2*dd09*hh10)+3*bb6*dd08+hh05;
      hh14 <- tgabb5;
      l3[,2]  <-  hh11*logy-2*dd09*dd08*dd13+bb6*bb4*dd13-(exp(dd07-3*v)*g3add11)/bb1^5-
         (4*exp(gg17+dd07)*gg18)/bb1^4+dd09*exp(gg17+m)*gg18-2*dd09*dd08*dd12+
	 bb6*bb4*dd12-bb2*bb4*gg14^2*g3abb5-bb2*bb4*(2*dd09*hh10-3*bb6*dd08+bb5)*hh14+
	 2*gg14*(bb6*dd08+hh05)*hh14+hh11*dgabb5;
      l3[,3]  <-  gg14*log(y)-bb6*bb4*log(1-y)-(exp(m-3*v)*g3add11)/bb1^4-
         (3*exp(gg17+m)*tgadd11)/bb1^3-bb6*bb4*dd12+bb6*exp(gg17+dd07)*gg14*g3abb5-
	 2*bb2*bb4*(bb6*dd08-bb2*bb4)*hh14+bb2*bb4*gg14*hh14+gg14*dgabb5;
      jj09 <- exp(-3*v);
      g3aff8 <- -gamma3(ff8)
      g3aaa3 <- -gamma3(aa3)
      g3aff4 <- -gamma3(ff4)
      l3[,4]  <-  (-bb2*cc2*logy)-bb2*aa3*log1y-dd09*jj09*g3aff8-3*bb6*ff7*tgaff8-
         bb2*aa3*dgaff8+jj09*g3aaa3+3*ff7*tgaaa3+aa3*dgaaa3-dd09*exp(3*m-3*v)*g3aff4-
	 3*bb6*exp(2*m-2*v)*tgaff4-bb2*cc2*dgaff4;

      g3 <- cbind(family$linfo[[1]]$d3link(m),family$linfo[[2]]$d3link(v))
    }

    if (deriv>3) {
      ## the fourth derivatives
      ## mmmm mmmv mmvv mvvv vvvv
      l4 <- matrix(0,n,5)
      kk13 <- 4*m;
      kk14 <- exp(bb3+kk13);
      kk15 <- 1/bb1^5;
      kk16 <- 24*kk15*exp(bb3+5*m)-60*gg12*kk14+50*dd09*gg11-15*bb6*dd08+bb5;
      kk17 <- 2*dd09*gg11-3*bb6*dd08+bb5;
      kk24 <- 1/bb1^6;
      kk25 <- -3*v;
      kk26 <- g3add11;
      g4add11 <- -gamma4(dd11)
      g4abb5 <- -gamma4(bb5)
      l4[,1]  <-  kk16*logy+24*kk15*kk14*dd13-36*gg12*gg11*dd13+14*dd09*dd08*dd13-
      bb6*bb4*dd13+(exp(kk13-4*v)*g4add11)/bb1^8+(12*exp(kk25+kk13)*kk26)/bb1^7-
      6*kk24*exp(kk25+gg10)*kk26+36*kk24*exp(gg17+kk13)*gg18-36*kk15*exp(gg17+gg10)*gg18+
      7*gg12*exp(gg17+dd07)*gg18+24*kk15*kk14*dd12-36*gg12*gg11*dd12+14*dd09*dd08*dd12-
      bb6*bb4*dd12+gg14^4*g4abb5+6*gg14^2*kk17*g3abb5+4*gg14*((-6*gg12*kk14)+12*dd09*gg11-
      7*bb6*dd08+bb5)*hh14+3*kk17^2*hh14+kk16*dgabb5;
      ll13 <- exp(bb3+4*m);
      ll14 <- 6*gg12*ll13-12*dd09*gg11+7*bb6*dd08+hh05;
      ll18 <- bb6*dd08+hh05;
      ll20 <- g3abb5;
      l4[,2]  <-  ll14*logy+6*gg12*gg11*dd13-6*dd09*dd08*dd13+bb6*bb4*dd13+
      (exp(gg10-4*v)*g4add11)/bb1^7+(9*exp(kk25+gg10)*kk26)/bb1^6-3*kk15*exp(kk25+dd07)*kk26+
      18*kk15*exp(gg17+gg10)*gg18-12*gg12*exp(gg17+dd07)*gg18+dd09*exp(gg17+m)*gg18+
      6*gg12*gg11*dd12-6*dd09*dd08*dd12+bb6*bb4*dd12-bb2*bb4*gg14^3*g4abb5-
      3*bb2*bb4*gg14*kk17*ll20+3*gg14^2*ll18*ll20-bb2*bb4*((-6*gg12*ll13)+12*dd09*gg11-
      7*bb6*dd08+bb5)*hh14+3*ll18*kk17*hh14+3*gg14*((-2*dd09*gg11)+3*bb6*dd08+hh05)*hh14+ll14*dgabb5;
      mm11 <- 2*dd09*hh10-3*bb6*dd08+bb5; 
      mm13 <- gg14^2;
      mm19 <- exp(gg17+dd07);
      l4[,3] <-  mm11*logy+2*dd09*dd08*dd13-bb6*bb4*dd13+(exp(dd07-4*v)*g4add11)/bb1^6+
      (7*exp(kk25+dd07)*kk26)/bb1^5-gg12*exp(kk25+m)*kk26+10*gg12*mm19*gg18-
      3*dd09*exp(gg17+m)*gg18+2*dd09*dd08*dd12-bb6*bb4*dd12+bb6*mm19*mm13*g4abb5+
      bb6*mm19*mm11*ll20+bb2*bb4*mm13*ll20-4*bb2*bb4*gg14*ll18*ll20+bb2*bb4*mm11*hh14-
      2*bb2*bb4*((-2*dd09*hh10)+3*bb6*dd08+hh05)*hh14+2*ll18^2*hh14+2*mm13*hh14+mm11*dgabb5;
      nn08 <- bb6*dd08-bb2*bb4;
      l4[,4]  <-  nn08*logy+bb6*bb4*log1y+(exp(m-4*v)*g4add11)/bb1^5+
      (6*exp(kk25+m)*g3add11)/bb1^4+7*dd09*exp(gg17+m)*tgadd11+bb6*bb4*dd12-
      dd09*exp(kk25+3*m)*gg14*g4abb5+3*bb6*mm19*nn08*ll20-3*bb6*mm19*gg14*ll20+
      3*bb2*bb4*nn08*hh14-4*bb2*bb4*gg14*hh14+nn08*dgabb5;
      oo11 <- exp(-4*v);
      l4[,5]  <-  bb2*cc2*logy+bb2*aa3*log1y - gg12*oo11*gamma4(ff8)+6*dd09*jj09*g3aff8+
      7*bb6*ff7*tgaff8+bb2*aa3*dgaff8 + oo11*gamma4(aa3)-6*jj09*g3aaa3-7*ff7*tgaaa3-
      aa3*dgaaa3 - gg12*exp(4*m-4*v)*gamma4(ff4)+6*dd09*exp(3*m-3*v)*g3aff4+
      7*bb6*exp(2*m-2*v)*tgaff4+bb2*cc2*dgaff4;

      g4 <- cbind(family$linfo[[1]]$d4link(m),family$linfo[[2]]$d4link(v))
    }
    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
   
      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- mgcv:::gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- mgcv:::gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                      d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D) 
    } else ret <- list()
    ret$l0 <- l0;ret$l <- l; ret
  } ## end ll betals

  sandwich <- function(y,X,coef,wt,family,offset=NULL) {
  ## compute filling for sandwich estimate of cov matrix
    ll(y,X,coef,wt,family,offset=NULL,deriv=1,sandwich=TRUE)$lbb
  }

 
  initialize <- expression({
  ## start out with mean given by least squres fit to logit(y)
  ## where y pre-initialized to never be quite 0 or 1.
  ## phi initially 0.5.
      n <- rep(1, nobs)
      ## should E be used unscaled or not?..
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      if (is.null(start)) {
        jj <- attr(x,"lpi")
        start <- rep(0,ncol(x))
        yt1 <- y #binomial()$linkfun(y) ## y preinitialized to have no actual 0s or 1s
        x1 <- x[,jj[[1]],drop=FALSE]
        e1 <- E[,jj[[1]],drop=FALSE] ## square root of total penalty
        if (use.unscaled) {
          qrx <- qr(rbind(x1,e1))
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(yt1,rep(0,nrow(E))))
          startji[!is.finite(startji)] <- 0       
        } else startji <- pen.reg(x1,e1,yt1)
        start[jj[[1]]] <- startji
        ## initialization sets phi=.5 initially
      }
  }) ## initialize Betals

  structure(list(family="betals",ll=ll,link=paste(link),nlp=2,ncv=ncv,sandwich=sandwich,predict=predict,
    tri = mgcv:::trind.generator(2), ## symmetric indices for accessing derivative arrays
    initialize=initialize,postproc=postproc,residuals=residuals,preinitialize=preinitialize,
    linfo = stats, ## link information list
    d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done    
    ls=1, ## signals that ls not needed here
    available.derivs = 2 ## can use full Newton here
    ),class = c("general.family","extended.family","family"))
} ## end betals


library(mgcv)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2 * x^11 * (10 * (1 - x))^6 + 10 * 
            (10 * x)^3 * (1 - x)^10
n <- 400
x0 <- runif(n, 0, 1)
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)

mu <- binomial()$linkinv(f0(x0)+f2(x2)/2-3)
phi <- binomial()$linkinv(f1(x1)-2)
a <- mu*(1/phi-1); b <- (1-mu)*(1/phi-1)
#a <- .5; b <- .5
y <- rbeta(n,a,b)
hist(y)
dat <- data.frame(y=y,x0=x0,x1=x1,x2=x2)

#m0 <- gam(list(y~1,~1),family=betals(),data=dat)

m1 <- gam(list(y~s(x0)+s(x2),~s(x1)),family=betals(),data=dat)
plot(m1,pages=1,scale=0,scheme=1)

gam.check(m1,type="pearson")

