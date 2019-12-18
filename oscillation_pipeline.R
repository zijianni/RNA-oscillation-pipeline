############################
#Step 0: configure parameters
#Change the following parameters based on your onw machine
#Also change the raw .xlsx file name because this name will be used among
#the whole pipeline and for generating final results
############################
rm(list=ls())
#directory of raw .xlsx files
raw_dir <- "~/Box/HES7oscillation/Human_HES7_data/H1-HES7mutationTest/01_C73T_testRcode"

#directory of output results
output_dir <- "~/Box/HES7oscillation/Human_HES7_data/H1-HES7mutationTest/01_C73T_testRcode"

#file position of Perl script "bandpass.pl"
bp_script <- "~/Box/HES7oscillation/Human_HES7_data/H1-HES7mutationTest/TBX6mutation/bandpass.pl"

#first bandwidth parameter
bw1 <- 45

#second bandwidth parameter
bw2 <- 3


###########################
#Step 1: load packages and functions
###########################
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("chron", quietly = TRUE)) install.packages("chron")

library(ggplot2)
library(readxl)   
library(chron)

dat_modify <- function(dat){    #modify data structure to fit the input format
    dat[,1] <- seq(from=0,to=5*(nrow(dat)-1),length.out = nrow(dat))
    colnames(dat)[1] <- "mins"
    return(dat)
}


find_peak_valley <- function(exp,time,degree=20,interval=c(100,1500)){   #function to call peak/valley positions and values
    obj <- lm(exp~poly(time, degree))
    x <- min(time):max(time)
    F.V <- predict(obj,data.frame(time=x))
    dif.F.V <- sign(diff(F.V))
    #dif.mat <- cbind((min(time)+1):max(time),dif.F.V)
    peak <- c()
    valley <- c()
    for(i in 1:(length(dif.F.V)-1)){
        is.valley <- all(dif.F.V[i:(i+1)]==c(-1,1))||all(dif.F.V[i:(i+1)]==c(-1,0))
        is.peak <- all(dif.F.V[i:(i+1)]==c(1,-1))||all(dif.F.V[i:(i+1)]==c(1,0))
        if(is.valley)
            valley <- c(valley,x[i+1])
        if(is.peak)
            peak <- c(peak,x[i+1])
    }
    p_res <- peak[(peak>=interval[1])&(peak<=interval[2])]
    v_res <- valley[(valley>=interval[1])&(valley<=interval[2])]
    return(list(peak=p_res,valley=v_res, peak_diff=diff(p_res),valley_diff=diff(v_res),fitted_value=F.V))
}

range_11 <- function(x){x/(max(x)-min(x))}  #for data scaling

plot_osc <- function(exp,time,plot_raw=F,raw=NULL,degree=20,interval=c(100,1500),ask=T,...){  #for plotting
    pv <- find_peak_valley(exp=exp,time=time,degree=degree,interval=interval)
    
    x <- min(time):max(time)
    obj <- lm(exp~poly(time,degree))
    F.V <- predict(obj,data.frame(time=x))
    
    plot(exp~time,type="l",ylab="Detrended expression",xlab="mins", 
         main="Detrended data and polynomial-fitted curve",...)
    lines(F.V~x,col="red")
    for(i in pv$peak){
        abline(v=i,col="blue",lty=2)
    }
    
    for(i in pv$valley){
        abline(v=i,col="green",lty=2)
    }
    legend("topright",lty=c(1,1,2,2),col=c("black","red","blue","green"),
           legend=c("detrended data","fitted value","peak","valley"))
    
    par(ask=ask)
    
    plot(range_11(F.V)~x,type="l",xlab="mins",ylab="scaled expression",col="red",
         main="Polynomial-fitted curve and lagged difference",...)
    lines(range_11(diff(F.V))~x[-1],col="pink")
    abline(0,0,lty="dashed")
    legend("topright",lty=c(1,1,2,2),col=c("red","pink","blue","green"),
           legend=c("fitted value","lagged diff","peak","valley"))
    for(i in pv$peak){
        abline(v=i,col="blue",lty=2)
    }
    
    for(i in pv$valley){
        abline(v=i,col="green",lty=2)
    }
    
    if(plot_raw){
        plot(raw~time,type="l",xlab="mins",ylab="raw expression",col="black",
             main="Raw data without detrending",...)
        #lines(range_11(diff(F.V))~x[-1],col="pink")
        #abline(0,0,lty="dashed")
        legend("topright",lty=c(1,2,2),col=c("black","blue","green"),
               legend=c("raw","peak","valley"))
        for(i in pv$peak){
            abline(v=i,col="blue",lty=2)
        }
        
        for(i in pv$valley){
            abline(v=i,col="green",lty=2)
        }
        
    }
    
    
    par(ask=F)
}

###########################
#Step 2: pre-processing raw .xlsx files
###########################
raw_dir <- gsub("/$", "", raw_dir)
raw_files <- list.files(raw_dir,pattern = ".xlsx")


for(i in seq_along(raw_files)){
    raw_temp <- read_xlsx(file.path(raw_dir,raw_files[i]),sheet = "Results Table",skip = 6)
    raw_temp <- raw_temp[,5:(which(colnames(raw_temp)=="Part")-1)]
    colnames(raw_temp)[1] <- "#Steps"
    tm <- round(raw_temp$`#Steps`,0)
    raw_temp[,1] <- paste0(tm%/%60,":",ifelse(tm%%60<10,paste0(0,tm%%60),tm%%60),":00")
    write.table(raw_temp,file = paste0(raw_dir,"/",gsub(".xlsx\\>","_raw.txt",raw_files[i])),
                quote = F, row.names = F,col.names = T,sep = "\t")
}

##########################
#Step 3: Run bandpass.pl
##########################
bp_file <- list.files(raw_dir,pattern = "_raw.txt")
bp_out <- gsub("_raw.txt","_detrended.txt",bp_file)

for(i in seq_along(bp_file)){
    cmd <- paste(bp_script, "<", file.path(raw_dir,bp_file[i]), bw1, bw2, ">", file.path(raw_dir,bp_out[i]), sep = " ")
    system(cmd)
}

#########################
#Step 4: Polynomial fitting and final results
#########################
Exp <- list()
Ind <- gsub(pattern = "(.*)_raw.txt" ,replacement ="\\1", x= bp_file)


for(i in seq_along(bp_file)){
    raw_dat <- read.delim(file.path(raw_dir,bp_file[i]),sep = "\t",header = T)
    detrend_dat <- read.delim(file.path(raw_dir,bp_out[i]),sep = "\t",header = T)
    Exp[[i]] <- list(raw=raw_dat,detrended=detrend_dat)
}

names(Exp) <- Ind

for(j in 1:length(Exp)){
    
    ###modify data to fit input format
    isRD <- is.data.frame(Exp[[j]][[1]])
    if(isRD){
        raw_temp <- Exp[[j]][[1]]
        raw_temp <- dat_modify(raw_temp)
        dat_temp <- Exp[[j]][[2]]
        dat_temp <- dat_modify(dat_temp)
    }else{
        dat_temp <- Exp[[j]]
        if(!is.factor(dat_temp[,1])){
            dat_temp[[1]] <- NULL
        }
        
        dat_temp <- dat_modify(dat_temp)
    }
    
    ###get peak/valley information (saved in `pvpv`)
    if(dim(dat_temp)[2]==2){
        pvpv <- list()
        pvpv[[1]] <- find_peak_valley(dat_temp[,2],time=dat_temp[,1],
                                      degree=25,interval = c(100,dat_temp[0.9*nrow(dat_temp),1]))
        names(pvpv) <- colnames(dat_temp)[2]
    }else{
        pvpv <- apply(dat_temp[,-1],2,find_peak_valley,time=dat_temp[,1],
                      degree=25,interval = c(100,dat_temp[0.9*nrow(dat_temp),1]))
    }
    
    ###pre-define empty variables to save peak/value information
    Exp_Peak <- matrix(NA,10,length(pvpv))
    colnames(Exp_Peak) <- names(pvpv)
    Exp_Peakval <- Exp_Peak
    Exp_Valley <- matrix(NA,10,length(pvpv))
    colnames(Exp_Valley) <- names(pvpv)
    Exp_Valleyval <- Exp_Valley
    Ave_Period <- numeric(length(pvpv))
    names(Ave_Period) <- names(pvpv)
    Exp_FV <- matrix(NA,max(dat_temp$mins)-min(dat_temp$mins)+1,length(pvpv)+1)
    Exp_FV[,1] <- min(dat_temp$mins):max(dat_temp$mins)
    colnames(Exp_FV) <- c("time",names(pvpv))
    
    
    ###extract information from `pvpv`
    for(i in 1:length(pvpv)){
        p_temp <- pvpv[[i]]$peak
        v_temp <- pvpv[[i]]$valley
        period_temp <- c(diff(p_temp),diff(v_temp))
        ###peak position
        Exp_Peak[1:length(p_temp),i] <- p_temp
        ###valley position
        Exp_Valley[1:length(v_temp),i] <- v_temp
        Ave_Period[i] <- mean(period_temp)
        ###fitted value of polynomial regression
        Exp_FV[,(i+1)] <- pvpv[[i]]$fitted_value
    }
    
    
    for(d1 in 1:10){   #assume no more than 10 peaks(valleys)
        ###peak values
        for(d2 in 1:ncol(Exp_Peakval)){
            if(!is.na(Exp_Peak[d1,d2])){
                Exp_Peakval[d1,d2] <- Exp_FV[which(Exp_FV[,1]==Exp_Peak[d1,d2]),d2+1]
            }
        }
        
        ###valley values
        for(d3 in 1:ncol(Exp_Valleyval)){
            if(!is.na(Exp_Valley[d1,d3])){
                Exp_Valleyval[d1,d3] <- Exp_FV[which(Exp_FV[,1]==Exp_Valley[d1,d3]),d3+1]
            }
        }
    }
    
    ###peak oscillation period
    Exp_Peak_Period <- apply(Exp_Peak,2,diff)
    rownames(Exp_Peak_Period) <- paste0("Period ",1:nrow(Exp_Peak_Period))
    ###valley oscillation period
    Exp_Valley_Period <- apply(Exp_Valley,2,diff)
    rownames(Exp_Valley_Period) <- paste0("Period ",1:nrow(Exp_Valley_Period))
    
    
    ###save results to output directory
    dir.create(file.path(output_dir,Ind[j]))
    
    write.table(Exp_Valley,file=paste0(output_dir,"/",Ind[j],"/",Ind[j],"_Valley.csv"),quote = F, na="",col.names = T,row.names = F,sep = ",")
    write.table(Exp_Peak,file=paste0(output_dir,"/",Ind[j],"/",Ind[j],"_Peak.csv"),quote = F, na="",col.names = T,row.names = F,sep = ",")
    
    write.table(Exp_Valleyval,file=paste0(output_dir,"/",Ind[j],"/",Ind[j],"_Valleyval.csv"),quote = F, na="",col.names = T,row.names = F,sep = ",")
    write.table(Exp_Peakval,file=paste0(output_dir,"/",Ind[j],"/",Ind[j],"_Peakval.csv"),quote = F, na="",col.names = T,row.names = F,sep = ",")
    
    write.table(Exp_Peak_Period,file=paste0(output_dir,"/",Ind[j],"/",Ind[j],"_Peak_Period.csv"),quote = F, na="",col.names = T,row.names = F,sep = ",")
    write.table(Exp_Valley_Period,file=paste0(output_dir,"/",Ind[j],"/",Ind[j],"_Valley_Period.csv"),quote = F, na="",col.names = T,row.names = F,sep = ",")

    write.table(Exp_FV,file=paste0(output_dir,"/",Ind[j],"/",Ind[j],"_fitted_values.csv"),quote = F, na="",col.names = T,row.names = F,sep = ",")
    
    ###generate plot to output directory
    fname <- paste0(output_dir,"/",Ind[j],"/",Ind[j],".pdf")
    pdf(fname,onefile = T,width = 10,height = 7)
    
    for(i in 1:(ncol(dat_temp)-1)){
        
        if(isRD){
            plot_osc(dat_temp[,i+1],dat_temp[,1],plot_raw = T,
                     raw = raw_temp[,i+1],degree = 25,interval = c(100,dat_temp[0.9*nrow(dat_temp),1]),ask=F,sub=colnames(dat_temp)[i+1])
            
        }else{
            plot_osc(dat_temp[,i+1],dat_temp[,1],plot_raw = F,
                     degree = 25,interval = c(100,dat_temp[0.9*nrow(dat_temp),1]),ask=F,sub=colnames(dat_temp)[i+1])
            
        }
    }
    dev.off()
}
