#!/usr/bin/Rscript

cor_ci <- function(x = NULL, y = NULL, r = NULL, n = NULL, alpha = 0.05){

   # if (length(args) != 2){
   #	stop("put in R^2 value and then n value")}

    trans <- function(r){
        return(log((1+r)/(1-r))/2)}

    if (is.null(n)){
    n <- length(x)}  #check y is same lenght
    if (is.null(r)){
    r <- cor(x,y)}

    z <- trans(r)
    l <- z - (trans(1-alpha/2))/(sqrt(n-3))
    u <- z + (trans(1-alpha/2))/(sqrt(n-3))
    l <- (exp(2*l)-1)/(exp(2*l)+1)
    u <- (exp(2*u)-1)/(exp(2*u)+1)
    ans <- c(round(l,4), round(u,4))
    return(ans)}


pro_dna <- function(tab, title){
    if (nrow(tab) > 1){
    	plot(tab$proEnt, tab$dnaEnt,
    	    main = paste(title, "Pro_DNA"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$proEnt),]$proEnt),
    	    ylim = c(0, tab[which.max(tab$dnaEnt),]$dnaEnt),
    	    xlab="Protein Entropy (bits)", ylab="DNA Entropy (bits)")

    	lmDNA <- lm(tab$dnaEnt ~ tab$proEnt, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 <- format(DNAsum$r.squared, digits = 4)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot coult not be created")}
}

#this makes R plots for DNA_Protein + puts them all in a pdf
dna_pro <- function(tab, title){
    if (nrow(tab) > 1){
    	plot(tab$dnaEnt, tab$proEnt,
    	    main = paste(title, "DNA_Pro"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$dnaEnt),]$dnaEnt),
    	    ylim = c(0, tab[which.max(tab$proEnt),]$proEnt),
    	    xlab="DNA Entropy (bits)", ylab="Protein Entropy (bits)")

    	lmDNA <- lm(tab$proEnt ~ tab$dnaEnt, data=tab)
    	DNAsum <- summary(lmDNA)
    	abline(lmDNA, col = "black")
    	R2 = format(DNAsum$r.squared, digits = 4) 
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

pro_cod <- function(tab, title){
    if (nrow(tab) > 1){
    	plot(tab$proEnt, tab$codEnt,
    	    main = paste(title, "Pro_Cod"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$proEnt),]$proEnt),
    	    ylim = c(0, tab[which.max(tab$codEnt),]$codEnt),
    	    xlab="Protein Entropy (bits)", ylab="Codon Entropy (bits)")

    	lmDNA <- lm(tab$codEnt ~ tab$proEnt, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 = format(DNAsum$r.squared, digits = 4)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

cod_pro <- function(tab, title){
    if (nrow(tab) > 1){
    	plot(tab$codEnt, tab$proEnt,
    	    main = paste(title, "Cod_Pro"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$codEnt),]$codEnt),
    	    ylim = c(0, tab[which.max(tab$proEnt),]$proEnt),
    	    xlab="Codon Entropy (bits)", ylab="Protein Entropy (bits)")

    	lmDNA <- lm(tab$proEnt ~ tab$codEnt, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 = format(DNAsum$r.squared, digits = 4)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

dna_cod <- function(tab, title){
    if (nrow(tab) > 1){
    	plot(tab$dnaEnt, tab$codEnt,
    	    main = paste(title, "DNA_Cod"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$dnaEnt),]$dnaEnt),
    	    ylim = c(0, tab[which.max(tab$codEnt),]$codEnt),
    	    xlab="DNA Entropy (bits)", ylab="Codon Entropy (bits)")

    	lmDNA <- lm(tab$codEnt ~ tab$dnaEnt, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 = format(DNAsum$r.squared, digits = 4)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

cod_dna <- function(tab, title){
    if (nrow(tab) > 1){
    	plot(tab$codEnt, tab$dnaEnt,
    	    main = paste(title, "Cod_DNA"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$codEnt),]$codEnt),
    	    ylim = c(0, tab[which.max(tab$dnaEnt),]$dnaEnt),
    	    xlab="Codon Entropy (bits)", ylab="DNA Entropy (bits)")

    	lmDNA <- lm(tab$dnaEnt ~ tab$codEnt, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 = format(DNAsum$r.squared, digits = 4)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}



ExecuteScript <- function(fileR){

    if (length(args) != 1){
	stop("must indicated and filename (which contains other filenames)")}

    print("started running")
    f <- readLines(fileR)[1]
    print(f)
    direction <- substr(f, nchar(f)-15, nchar(f)-9) 
    print(direction)

    if (direction == "pro_dna"){
        filename <- paste0("period_PDplot.pdf")}
    else if (direction == "dna_pro"){
        filename <- paste0("period_DPplot.pdf")}
    else if (direction == "pro_cod"){
        filename <- paste0("period_PCplot.pdf")}
    else if (direction == "cod_pro"){
        filename <- paste0("period_CPplot.pdf")}
    else if (direction == "dna_cod"){
        filename <- paste0("period_DCplot.pdf")}    
    else if (direction == "cod_dna"){
        filename <- paste0("period_CDplot.pdf")}

    df <- data.frame("parameters", "rep", "period", "seq", "numdatapoints", "R2", "low", "high")
    pdf(filename)
    for (file in readLines(fileR)){
        tab <- read.table(file, header=TRUE)
        tab <- na.omit(tab)
        num.obs <- as.character(nrow(tab))
        print(num.obs)
        parameters <- substr(file, nchar(file)-24, nchar(file)-17)  #change for each organism
	reppat <- regexpr("_[0-9][0-9][0-9]_", file)  #for if patterns are id at the protein level
	repnum <- substr(file, reppat+1, reppat+3)
	organism <- substr(file, reppat+5, nchar(file)-25)
	
	seqtype <- grep("dna_repeats", file)   #if patterns are identified at the DNA level
	if (length(seqtype) == 1){
	    reppat <- regexpr("_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_", file)
	    repnum <- substr(file, reppat+1, reppat+10)
	    organism <- substr(file, reppat+12, nchar(file)-25)
		}
	print(repnum)
	print(parameters)
	print(organism)

	#this is to feed the plot the title
	title <- paste(organism, parameters, repnum) 
	w <- grep("wholeseq", file)
	noper <- grep("noperiodicity", file)
	if (length(noper) == 1){
	    title <- paste(title, "noperiod")
	    per <- "FALSE"
	    	if (length(w) == 1){
	    	    title <- paste(title, "wholeseq")
	            w <- "wholeseq"}
		
		else {
	    	    title <- paste(title, "partseq")
	            w <- "partseq"}
		}
	else {
	    title <- paste(title, "period")
	    per <- "TRUE"
	    w <- "  "}
	 	

        if (direction == "pro_dna"){
            R2 <- pro_dna(tab, title)}      #plots will be made and R2 returned
        else if (direction == "dna_pro"){
            R2 <- dna_pro(tab, title)}
        else if (direction == "pro_cod"){
            R2 <- pro_cod(tab, title)}
        else if (direction == "cod_pro"){
            R2 <- cod_pro(tab, title)}
        else if (direction == "dna_cod"){
            R2 <- dna_cod(tab, title)}
        else if (direction == "cod_dna"){
            R2 <- cod_dna(tab, title)}
	CI <- cor_ci(r=sqrt(as.numeric(R2)), n=as.numeric(num.obs))
	low <- CI[1]
	high <- CI[2]
        df <- rbind(df, c(parameters, repnum, per, w, num.obs, R2, low, high))
    }
    print("writing table")
    write.table(df, file = paste0("period_numdatapoints_", organism, "_", as.character(parameters),"_", repnum,"_", direction), sep="\t")
    print("num_datapoints written")
    dev.off()
}

#ExecuteScript("test")
args <- commandArgs(trailingOnly = TRUE)
ExecuteScript(args[1]) 
#second argument will be rep number (exp. 422)
