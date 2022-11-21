HWE_B <- function(){
  Version <- "HWE_Bayes_Forensic v.1.0"
  
  if(require(tcltk)){
    cat('"tcltk" is loaded correctly \n')
  }else{
    cat('Trying to install "tcltk"... \n')
    install.packages("tcltk")
    if(require("tcltk")){
      cat('"tcltk" installed and loaded \n')
    }else{
      stop(tk_messageBox("The package \"tcltk\" could not be installed.\n Please execute the following command in the console and check the error message.\n >install.packages(\"tcltk\").", icon = "info", type = "ok"))
    }
  } #check tkltk package
  
  if(require(tcltk2)){
    cat('"tcltk2" is loaded correctly \n')
  }else{
    cat('trying to install "tcltk2"... \n')
    install.packages("tcltk2")
    if(require("tcltk2")){
      cat('"tcltk2" installed and loaded \n')
    }else{
      stop(tk_messageBox("The package \"tcltk2\" could not be installed.\n Please execute the following command in the console and check the error message.\n >install.packages(\"tcltk2\")", icon = "info", type = "ok"))
    }
  } #check tcltk2 package
  
  if(require(bridgesampling)){
    cat('"bridgesampling" is loaded correctly \n')
  }else{
    cat('trying to install "bridgesampling"... \n')
    install.packages("bridgesampling")
    if(require("bridgesampling")){
      cat('"bridgesampling" installed and loaded \n')
    }else{
      stop(tk_messageBox("The package \"bridgesampling\" could not be installed.\n Please execute the following command in the console and check the error message.\n >install.packages(\"bridgesampling\")", icon = "info", type = "ok"))
    }
  } #check bridgesampling package
  
  if(require(rstan)){
    cat('"rstan" is loaded correctly \n')
  }else{
    stop(tk_messageBox("Please install \"RStan\" by referring to \"RStan Getting Started\" on the official website, and then run this program again after confirming the operation.", icon = "info", type = "ok"))
  } #check RStan environment
  
  if(require(ggplot2)){
    cat('"ggplot2" is loaded correctly \n')
  }else{
    cat('Trying to install "ggplot2"... \n')
    install.packages("ggplot2")
    if(require("ggplot2")){
      cat('"ggplot2" installed and loaded \n')
    }else{
      stop(tk_messageBox("The package \"ggplot2\" could not be installed.\n Please execute the following command in the console and check the error message.\n >install.packages(\"ggplot2\").", icon = "info", type = "ok"))
    }
  } #check ggplot2 package
  
  Sm_f <- "
        data {
          real N;
          int<lower=0> Ka;
          int<lower=0> Kg;
          int<lower=0> Dg[Kg];
          vector<lower=0>[Ka] Da; 
          real m;
        }
    
        parameters {
          simplex[Ka] pa;
          real<lower=m,upper=1> f;
        }
    
        transformed parameters {
          vector[Kg] pg;
          {
            real pmin;
            int c;	
            c = 0;
            for(i in 1:Ka){
        			for(j in 1:i){
      				c += 1;
      				if(i == j)
      					if(pa[i] < 0.5)
			      			pmin = pa[i];
			      		else
			      			pmin = 1-pa[i];
	      				if(f < (-pmin)/(1-pmin))
		      				pg[c] = 1/(2*N)^2 ;
			      		else
				      		pg[c] = pa[i]^2 + f*pa[i]*(1-pa[i]) + 1/(2*N)^2;
			      	if(i != j)
				      		pg[c] = (1-f)*2*pa[i]*pa[j] - 2/(2*N)^2;
		        	}
	        	}
          pg = pg/sum(pg);
          }
        }
        
        model {
          target += dirichlet_lpdf(pa|Da+1);
          target += normal_lpdf(f|(m+1)/2,(1-m)/2);
          target += multinomial_lpmf(Dg|pg);
        }
      "
  Sm <- "
        data {
          real N;
          int<lower=0> Ka;
          int<lower=0> Kg;
          int<lower=0> Dg[Kg];
          vector<lower=0>[Ka] Da; 
          real m;
        }
    
        parameters {
          simplex[Ka] pa;
        }
    
        transformed parameters {
          simplex[Kg] pg;
          {
            int c;	
            c = 0;
        
            for(i in 1:Ka){
              for(j in 1:i){
                c += 1;
                pg[c] = (i == j) 
                ? pa[i]^2 
                : 2*pa[i]*pa[j];
              }
            }
          }
        }
      
        model {
          target += dirichlet_lpdf(pa|Da+1);
          target += multinomial_lpmf(Dg|pg);
        }
      "
  if(exists("F_stanmodel")){cat('Loading F_stanmodel\n')} else {
    cat('Compiling F_stanmodel. Please wait for a few minutes.\n')
    F_stanmodel <<- stan_model(model_code = Sm_f)}
  
  if(exists("nF_stanmodel")){cat('Loading nF_stanmodel\n')} else {
    cat('Compiling nF_stanmodel. Please wait for a few minutes.\n')
    nF_stanmodel <<- stan_model(model_code = Sm)}
  
  Opendirectory<- function(dp, var, Cursor){
    DirectoryName <- tkchooseDirectory(initialdir = tclvalue(dp))
    tclvalue(dp) <- DirectoryName
    tclvalue(var) <- DirectoryName
    tclvalue(Cursor) <- "hand2"
  } #choose directory
  
  Opencsv <- function(fp, var, top, filestate, Cursor){
    FileName <- tclvalue(tkgetOpenFile(parent = top, initialdir = tclvalue(fp), multiple = "true", filetypes = "{{CSV Files} {.csv}}"))
    if(!nchar(FileName)){
      tkmessageBox(message = "No file was selected!", icon = "error", type = "ok")
    }else{
      tmp <- sub("\\}", FileName, replacement = "")
      tmp2 <- sub("\\{", tmp, replacement = "")
      tclvalue(fp) <- tmp2
      foo3 <- strsplit(tmp2, "/")[[1]]
      tclvalue(var) <- strsplit(foo3[length(foo3)], "\\.csv")[[1]][1]
      tclvalue(filestate) <- "normal"
      tclvalue(Cursor) <- "hand2"
    }
  } #choose data.csv
  
  Popdata.make <- function(Popdata.filepath){
    Popdata.File <- read.csv(Popdata.filepath, header=T)
  } #adjuast_data.csv
  
  Alleledata.make <- function(pop){
    if(ncol(pop)==1){
      Alleledata.File <- list(table(pop[,1]))
    }else{
      Alleledata.File <- lapply(pop,table)
    }
    names(Alleledata.File) <- names(pop)
    return(Alleledata.File)
  } #aggregate observed allele number
  
  Genotypedata.make <- function(pop, an){
    gt <- list()
      for (i in 1:ncol(pop)) {
        gx <- matrix(0, length(as.vector(an[[i]])), length(as.vector(an[[i]])))
        rownames(gx) <- names(an[[i]])
        colnames(gx) <- names(an[[i]])
        pop_rn <- pop[[i]][!is.na(pop[[i]])]
        for (j in 1:(sum(as.vector(an[[i]]))/2)) {
          if(pop_rn[2*j-1] <= pop_rn[2*j]){
            gx[which(rownames(gx)==pop_rn[2*j-1]),which(colnames(gx)==pop_rn[2*j])]<-gx[which(rownames(gx)==pop_rn[2*j-1]),which(colnames(gx)==pop_rn[2*j])]+1
          }else{
            gx[which(rownames(gx)==pop_rn[2*j]),which(colnames(gx)==pop_rn[2*j-1])]<-gx[which(rownames(gx)==pop_rn[2*j]),which(colnames(gx)==pop_rn[2*j-1])]+1
          } 
        }
        gt <- c(gt,list(gx))
      }
    names(gt) <- names(an)
    return(gt)
  } #aggregate obseved genotype number
  
  Fis.calc <- function(pop, an, gn){
    Fis <- matrix(0,ncol(pop),3)
      for(i in 1:ncol(pop)){
        he <- 1-sum((an[[i]]/sum(an[[i]]))^2)#Expected Hz
        ho <- 1-(sum(diag(gn[[i]]))/sum(gn[[i]]))#Observed Hz
        hf <- round((he-ho)/he,3)
        Fis[i,1] <- hf
        Fis[i,2] <- ho
    }
    return(Fis)
  } #calculate Fis by heterozygosity
    
  Standata.make <- function(pop, an, gn){
    N <- sum(as.vector(an))
    Ka <- length(as.vector(an))
    Kg <- Ka*(Ka+1)/2
    Da <- as.vector(an)
    Dg <- as.vector(gn[upper.tri(gn, diag = TRUE)])
    m <- round(1/(1-Ka),3)
    Sdata <- list(N = N, Ka = Ka, Kg = Kg, Dg = Dg, Da = Da, m = m)
    return(Sdata)
  } #make input data for stan
  
  Ratio.calc <- function(f,hz){
      no <- 1/(1-hz)
      ro <- round(2*(1-f)/(f*(no-1)+1),3)
    return(ro)
  } #calculate homo:hetero ratio for estimated f
  
  
  RunningStan <- function(model, data, ch, it, wa, th){
    RS <- sampling(model, data, chains = ch, iter = it, warmup = wa, thin = th, control=list(adapt_delta=0.99))
  } #run stan code
  
  Saveresult <- function(pop, locus, index, Min, First, Med, Mode, Third, Max, Rhat, Ratio, Fis, BF01, BF10){
    Resultcsv <- matrix("", ncol(pop), 11)
    rownames(Resultcsv) <- locus
    colnames(Resultcsv) <- index
    Resultcsv[,1] <- Min
    Resultcsv[,2] <- First
    Resultcsv[,3] <- Med
    Resultcsv[,4] <- Mode
    Resultcsv[,5] <- Third
    Resultcsv[,6] <- Max
    Resultcsv[,7] <- Rhat
    Resultcsv[,8] <- Ratio
    Resultcsv[,9] <- Fis[,1]
    Resultcsv[,10] <- BF01
    Resultcsv[,11] <- BF10
    
    Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
    if(tclvalue(Save.as) != ""){
      if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
        ReportName <- tclvalue(Save.as)
      }else{
        ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
      }
      write.csv(Resultcsv, file = ReportName, row.names = TRUE)
    }
  } #save results
  
  Fv <- function(){
    Tab1 <- function(){
      tkdestroy(frame1)
      frame1 <<- tkframe(tab1)
      top.label <- tklabel(frame1, text = "This program provide Bayesian approach of test for Hardy-Weinberg equilibrium of Forensic population data.")
      tkgrid(top.label, padx = 20, pady = 20, sticky = "w")
      tkgrid(frame1)
    } #top   
    
    Tab2 <- function(){
      tkdestroy(frame2)
      frame2 <<- tkframe(tab2)
      frame2.label <- tklabel(frame2, text = "Select the directory where the data is located. The results will also be saved here.")
      tkgrid(frame2.label, padx = 20, pady = 20 ,columnspan=3)
      
      #choose dir
      Directory.label <- tklabel(frame2, text = "Directory")
      dataname0.label <- tklabel(frame2, textvariable = Directory.var, width = 40, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      Directory.butt <- tkbutton(frame2, text = "Browse", cursor = "hand2", command = function() Opendirectory(Directorypath, Directory.var, Directory.Cursor))
      tkgrid(Directory.label, dataname0.label, Directory.butt, padx = 20, pady = 20, sticky = "w")
      
      #choose data.csv
      Popdata.label <- tklabel(frame2, text = "Population data")
      dataname.label <- tklabel(frame2, textvariable = Popdata.var, width = 40, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      Popdata.butt <- tkbutton(frame2, text = "Load", cursor = "hand2", command = function() Opencsv(Popdata.filepath, Popdata.var, tf, Popdata.butt.state, Popdata.butt.cursor))
      tkgrid(Popdata.label, dataname.label, Popdata.butt, padx = 20, pady = 20, sticky = "w")
      
      tkgrid(frame2)
    } #load file
    
    Tab3 <- function(){
      tkdestroy(frame3)
      frame3 <<- tkframe(tab3)
      frame3.label <- tklabel(frame3, text = "Parameters settings")
      tkgrid(frame3.label, padx = 20, pady = 20, sticky = "w")
      
      Chains.label <- tklabel(frame3, text = "Chain (Default = 4)")
      Chains.entry <- tkentry(frame3, textvariable = Chains, width = 25, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      tkgrid(Chains.label, Chains.entry, padx = 20, pady = 20, sticky = "w")
      
      Iter.label <- tklabel(frame3, text = "Iteration (Default = 10000)")
      Iter.entry <- tkentry(frame3, textvariable = Iter, width = 25, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      tkgrid(Iter.label, Iter.entry, padx = 20, pady = 20, sticky = "w")
      
      Warmup.label <- tklabel(frame3, text = "Warmup (Default = 500)")
      Warmup.entry <- tkentry(frame3, textvariable = Warmup, width = 25, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      tkgrid(Warmup.label, Warmup.entry, padx = 20, pady = 20, sticky = "w")
      
      Thin.label <- tklabel(frame3, text = "Thinning (Default = 1)")
      Thin.entry <- tkentry(frame3, textvariable = Thin, width = 25, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      tkgrid(Thin.label, Thin.entry, padx = 20, pady = 20, sticky = "w")
      
      tkgrid(frame3)
    } #parameters for MCMC sampling
    
    Tab4 <- function(){
      tkdestroy(frame4)
      frame4 <<- tkframe(tab4)
      
      option.label <- tklabel(frame4, text = "Check the box for additional output files to save.")
      tkgrid(option.label, padx = 20, pady = 20, sticky = "w")
      
      fitoption <- tk2checkbutton(frame4, text="fit_f_[Locus name] .Rda")
      tkconfigure(fitoption, variable = cbValuef)
      tkgrid(fitoption, padx = 20, pady = 20, sticky = "w")
      
      fitoption2 <- tk2checkbutton(frame4, text="fit_0_[Locus name] .Rda")
      tkconfigure(fitoption2, variable = cbValue0)
      tkgrid(fitoption2, padx = 20, pady = 20, sticky = "w")
      
      pngoption <- tk2checkbutton(frame4, text="EstmatedF_[Locus name] .png")
      tkconfigure(pngoption, variable = cbValuep)
      tkgrid(pngoption, padx = 20, pady = 20, sticky = "w")
      
      estimate.butt <- tkbutton(frame4, text = "Run", cursor = "hand2", command = function() Tab5())
      tkgrid(estimate.butt, padx = 20, pady = 20, sticky = "w")
      
      tkgrid(frame4)
    } #select output files
    
    Tab5 <- function(){
      tkdestroy(frame5)
      frame5 <<- tkframe(tab5)
      if(exists('Finish.label')){tkdestroy(Finish.label)}
      
      if(tclvalue(Popdata.filepath) != ""){
        Popdata <- Popdata.make(tclvalue(Popdata.filepath))
        setwd(tclvalue(Directorypath))
        
        pb.label <- tk2label(frame4, text = "The locus count for which calculation has been completed: 0")
        pb <- tk2progress(frame4, length = 200)
        tkconfigure(pb, value = 0, maximum = ncol(Popdata))
        tkgrid(pb.label, padx = 20, pady = 20, sticky = "w")
        tkgrid(pb, padx = 20, pady = 20, sticky = "w")
        tcl("update")
        
        Allele <- Alleledata.make(Popdata)
        Genotype <- Genotypedata.make(Popdata, Allele)
        Fis <- Fis.calc(Popdata, Allele, Genotype)
        Locus <- c(colnames(Popdata))
        Index <- c("Min", "2.5%", "Median", "Mode", "97.5%", "Max", "Rhat", "Ratio", "F_IS","BF01","BF10")
        Min <- c()
        First <- c()
        Med <- c()
        Mode <- c()
        Third <- c()
        Max <- c()
        Rhat <- c()
        Ratio <- c()
        BF01 <- c()
        BF10 <- c()
        
        options(mc.cores = parallel::detectCores())
        rstan_options(auto_write = TRUE)
        
        for(n in 1:ncol(Popdata)){
          Standata <- Standata.make(Popdata, Allele[[n]], Genotype[[n]])

          fitH0 <- RunningStan(nF_stanmodel, Standata, as.numeric(tclvalue(Chains)), as.numeric(tclvalue(Iter)), as.numeric(tclvalue(Warmup)), as.numeric(tclvalue(Thin)))
          fitH1 <- RunningStan(F_stanmodel, Standata, as.numeric(tclvalue(Chains)), as.numeric(tclvalue(Iter)), as.numeric(tclvalue(Warmup)), as.numeric(tclvalue(Thin)))
          
          H0 <- bridge_sampler(fitH0, method = "warp3", silent = TRUE)
          H1 <- bridge_sampler(fitH1, method = "warp3", silent = TRUE)
          
          BF0 <- bf(H0,H1)
          BF1 <- bf(H1,H0)
          BF01 <- c(BF01, round(BF0[[1]],3)) 
          BF10 <- c(BF10, round(BF1[[1]],3))
          
          Ms_f <- rstan::extract(fitH1)
          ms95 <- quantile(Ms_f$f, probs=c(0, 0.025, 0.5, 0.975, 1))
          rh <- all(summary(fitH1)$summary[,"Rhat"] <= 1.10, na.rm = T)
          mo <- density(Ms_f$f)$x[which.max(density(Ms_f$f)$y)]
          rt <- Ratio.calc(mo,Fis[n,2])
          
          Min <- c(Min, round(ms95[1], 3))
          First <- c(First, round(ms95[2], 3))
          Med <- c(Med, round(ms95[3], 3))
          Mode <- c(Mode, round(mo,3))
          Third <- c(Third, round(ms95[4], 3))
          Max <- c(Max, round(ms95[5], 3))
          Rhat <- c(Rhat, rh)
          Ratio <- c(Ratio,rt)
          
          if(tclvalue(cbValuef) == 1){
            ffile.name <- sprintf("fit_f_%s.Rda", Locus[[n]])
            saveRDS(fitH1, file = ffile.name)
          }
          
          if(tclvalue(cbValue0) == 1){
            ffile.name <- sprintf("fit_0_%s.Rda", Locus[[n]])
            saveRDS(fitH0, file = ffile.name)
          }
          
          if(tclvalue(cbValuep) == 1){
            pfile.name <- sprintf("EstimatedF_%s.png", Locus[[n]])
            dens <- density(Ms_f$f)
            d <- data.frame(X=dens$x, Y=dens$y)
            v <- data.frame(V=c(mo,ms95[3]),Values=c("Mode","Median"))
            g1 <-  ggplot()+
              geom_ribbon(data=subset(d, X>=ms95[2] & X<=ms95[4]),
                          aes(x=X, ymin=0, ymax=Y),fill="LightGray")+
              geom_line(data=d, aes(x=X, y=Y),size=1)+
              geom_vline(data=v,aes(xintercept=V,colour=Values),
                         show.legend = TRUE,size=2)+
              geom_vline(xintercept=0,linetype="dashed",size=1)+
              theme_minimal(base_size = 50)+
              xlab("Estemated f value")+
              ylab("Density")+
              theme(axis.title.y=element_blank(),axis.text.y=element_blank())+
              scale_y_continuous(breaks = NULL)+
              ggtitle(sprintf("Density of estimated f at %s", Locus[[n]]))
            g2 <- stan_dens(fitH1,pars="f",separate_chains = T)+
              theme(axis.title.x=element_text(size=50),
                    axis.text.x = element_text(size=40),
                    legend.title = element_text(size=50),
                    legend.text = element_text(size=40))
            G1 <- ggplot_gtable(ggplot_build(g1)) 
            G2 <- ggplot_gtable(ggplot_build(g2))
            png(pfile.name, width = 1200, height = 1200)
              grid::grid.draw(rbind(G1,G2))
            dev.off()

          }
          
          rm(fitH0)
          rm(fitH1)
          rm(Ms_f)
          
          tkconfigure(pb.label, text = paste("The locus count for which calculation has been completed: ",n))
          tkconfigure(pb, value = n)
          tcl("update")
        }
        
        tkdestroy(pb)
        tkdestroy(pb.label)
        
        Finish.label <<- tklabel(frame4, text = "Done. Check the Results tab.")
        tkgrid(Finish.label, padx = 20, pady = 20, sticky = "w")
        Success.label <- tklabel(frame5, text = "Save the result CSV file.")
        tkgrid(Success.label, padx = 20, pady = 20, sticky = "w")
        Save.butt <- tkbutton(frame5, text = " Save ", cursor = "hand2", command = function() Saveresult(Popdata, Locus, Index, Min, First, Med, Mode, Third, Max, Rhat, Ratio, Fis, BF01, BF10))
        tkgrid(Save.butt, padx = 20, pady = 20, sticky = "w")
      }
      tkgrid(frame5)
    } #run & save
    
    #definition
    Directorypath <- tclVar("")
    Directory.var <- tclVar("")
    Directory.Cursor <- tclVar("")
    Popdata.filepath <- tclVar("")
    Popdata.var <- tclVar("")
    Popdata.butt.state <- tclVar("disabled")
    Popdata.butt.cursor <- tclVar("arrow")
    Chains <- tclVar("4")
    Iter <- tclVar("10000")
    Warmup <- tclVar("500")
    Thin <- tclVar("1")
    cbValuef <- tclVar("0")
    cbValue0 <- tclVar("0")
    cbValuep <- tclVar("0")
    Locus <- tclVar("")
    Index <- tclVar("")
    co <- tclVar("")
    Min <- tclVar("")
    First <- tclVar("")
    Med <- tclVar("")
    Mode <- tclVar("")
    Third <- tclVar("")
    Max <- tclVar("")
    Rhat <- tclVar("")
    ffile.name <- tclVar("")
    mfile.name <- tclVar("")
    
    #Tabs
    tf <- tktoplevel()
    tkwm.title(tf, Version)
    tabs <- tk2notebook(tf, tabs = c("Top","Step 1", "Step 2", "Step 3", "Results"))
    tkpack(tabs, fill = "both", expand = 1)
    tab1 <- tk2notetab(tabs, "Top")
    tab2 <- tk2notetab(tabs, "Step 1")
    tab3 <- tk2notetab(tabs, "Step 2")
    tab4 <- tk2notetab(tabs, "Step 3")
    tab5 <- tk2notetab(tabs, "Results")
    frame1 <- tkframe(tab1)
    Tab1()
    frame2 <- tkframe(tab2)
    Tab2()
    frame3 <- tkframe(tab3)
    Tab3()
    frame4 <- tkframe(tab4)
    Tab4()
    frame5 <- tkframe(tab5)
    tcl("wm", "attributes", tf, topmost=TRUE)
  }
  Fv()
}
HWE_B()

