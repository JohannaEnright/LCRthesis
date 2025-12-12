#!/usr/bin/Rscript


#confidence interval calc
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
    ans <- c(round(l,3), round(u,3))
    return(ans)}


#this script is meant to make a bunch of R plots from the entropy data
#takes a file of filenames (one for each line) (use ls | more > filename)
#this makes the R plots from Pro_DNA and puts them all in a pdf
#also will write the number of points for each plot of to a file

pro_dna <- function(tab, title, parameters){
    print(parameters)
    if (nrow(tab) > 1){
    	plot(tab$EntropyProtein, tab$EntropyDNA,
    	    main = paste(title, parameters, "Entropy Pro_DNA"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$EntropyProtein),]$EntropyProtein),
    	    ylim = c(0, tab[which.max(tab$EntropyDNA),]$EntropyDNA),
    	    xlab="Protein Entropy (bits)", ylab="DNA Entropy (bits)")

    	lmDNA <- lm(tab$EntropyDNA ~ tab$EntropyProtein, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 <- format(DNAsum$r.squared, digits = 3)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot coult not be created")}
}

#this makes R plots for DNA_Protein + puts them all in a pdf
dna_pro <- function(tab, title, parameters){
    print(parameters)
    if (nrow(tab) > 1){
    	plot(tab$EntropyDNA, tab$EntropyProtein,
    	    main = paste(title, parameters, "Entropy DNA_Pro"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$EntropyDNA),]$EntropyDNA),
    	    ylim = c(0, tab[which.max(tab$EntropyProtein),]$EntropyProtein),
    	    xlab="DNA Entropy (bits)", ylab="Protein Entropy (bits)")

    	lmDNA <- lm(tab$EntropyProtein ~ tab$EntropyDNA, data=tab)
    	DNAsum <- summary(lmDNA)
    	abline(lmDNA, col = "black")
    	R2 = format(DNAsum$r.squared, digits = 3)
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

pro_cod <- function(tab, title, parameters){
    print(parameters)
    if (nrow(tab) > 1){
    	plot(tab$EntropyProtein, tab$EntropyCodon,
    	    main = paste(title, parameters, "Entropy Pro_Cod"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$EntropyProtein),]$EntropyProtein),
    	    ylim = c(0, tab[which.max(tab$EntropyCodon),]$EntropyCodon),
    	    xlab="Protein Entropy (bits)", ylab="Codon Entropy (bits)")

    	lmDNA <- lm(tab$EntropyCodon ~ tab$EntropyProtein, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 = format(DNAsum$r.squared, digits = 3)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

cod_pro <- function(tab, title, parameters){
    print(parameters)
    if (nrow(tab) > 1){
    	plot(tab$EntropyCodon, tab$EntropyProtein,
    	    main = paste(title, parameters, "Entropy Cod_Pro"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$EntropyCodon),]$EntropyCodon),
    	    ylim = c(0, tab[which.max(tab$EntropyProtein),]$EntropyProtein),
    	    xlab="Codon Entropy (bits)", ylab="Protein Entropy (bits)")

    	lmDNA <- lm(tab$EntropyProtein ~ tab$EntropyCodon, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 = format(DNAsum$r.squared, digits = 3)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

dna_cod <- function(tab, title, parameters){
    print(parameters)
    if (nrow(tab) > 1){
    	plot(tab$EntropyDNA, tab$EntropyCodon,
    	    main = paste(title, parameters, "Entropy DNA_Cod"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$EntropyDNA),]$EntropyDNA),
    	    ylim = c(0, tab[which.max(tab$EntropyCodon),]$EntropyCodon),
    	    xlab="DNA Entropy (bits)", ylab="Codon Entropy (bits)")

    	lmDNA <- lm(tab$EntropyCodon ~ tab$EntropyDNA, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 = format(DNAsum$r.squared, digits = 3)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

cod_dna <- function(tab, title, parameters){
    print(parameters)
    if (nrow(tab) > 1){
    	plot(tab$EntropyCodon, tab$EntropyDNA,
    	    main = paste(title, parameters, "Entropy Cod_DNA"),
    	    pch = 21, cex = 0.5, col = "green", bg="black", lwd=1,
    	    xlim = c(0, tab[which.max(tab$EntropyCodon),]$EntropyCodon),
    	    ylim = c(0, tab[which.max(tab$EntropyDNA),]$EntropyDNA),
    	    xlab="Codon Entropy (bits)", ylab="DNA Entropy (bits)")

    	lmDNA <- lm(tab$EntropyDNA ~ tab$EntropyCodon, data=tab)
    	DNAsum <- summary(lmDNA)
    	R2 = format(DNAsum$r.squared, digits = 3)
    	abline(lmDNA, col = "black")
    	legend("bottomright", legend=paste("R2=", R2))
    	print("plot created")
    	return(R2)}
    else{
	return("plot could not be created")}
}

#this is the function you want to give parameters to every time
#direction is Pro_DNA or DNA_Pro, fileR is a list of the filenames
#containing the Entropy Data
#the file contining the file names does not like spaces at the end
#t is file title, title is plot title
ExecuteScript <- function(direction, fileR, t=""){

    if (length(args) != 3){
	stop("must indicated 'direction' and filename (which contains other filenames) and title")}

    if (direction != "pro_dna" && direction != "dna_pro" &&
	direction != "pro_cod" && direction != "cod_pro" &&
	direction != "dna_cod" && direction != "cod_dna"){
	stop("args[1] must be 'pro_dna' or 'dna_pro' or 
	     'cod_pro' or 'pro_cod' or 'cod_dna' or 'dna_cod'")}

    print("started running")
    #print(title)
    print(readLines(fileR)[1])  #print the file

    if (direction == "pro_dna"){
        filename <- paste0("PDplot_", t, ".pdf")}
    else if (direction == "dna_pro"){
        filename <- paste0("DPplot_", t, ".pdf")}
    else if (direction == "pro_cod"){
        filename <- paste0("PCplot_", t, ".pdf")}
    else if (direction == "cod_pro"){
        filename <- paste0("CPplot_", t, ".pdf")}
    else if (direction == "dna_cod"){
        filename <- paste0("DCplot_", t, ".pdf")}    
    else if (direction == "cod_dna"){
        filename <- paste0("CDplot_", t, ".pdf")}
    else {
        stop("options for directions are 'pro_dna' or 'dna_pro'\
            or 'pro_cod' or 'cod_pro' or 'cod_dna' or 'dna_cod'")}

    df <- data.frame("parameters", "R2", "r", "numdatapoints", "lowB", "highB")
    pdf(filename)
    for (file in readLines(fileR)){
        tab <- read.table(file, header=TRUE)
        tab <- na.omit(tab)
        num.obs <- as.character(nrow(tab))
        print(num.obs)
        parameters <- substr(file, nchar(file)-22, nchar(file)-15)  #change for each title
	lastslash <- gregexpr("/", file)
	lastslash <- lastslash[[1]][length(lastslash[[1]])]
    	title <- substr(file, lastslash+1, nchar(file)-33)
	print(title)

        if (direction == "pro_dna"){
            R2 <- pro_dna(tab, title, parameters)} #plots will be made and R2 returned
        else if (direction == "dna_pro"){
            R2 <- dna_pro(tab, title, parameters)}
        else if (direction == "pro_cod"){
            R2 <- pro_cod(tab, title, parameters)}
        else if (direction == "cod_pro"){
            R2 <- cod_pro(tab, title, parameters)}
        else if (direction == "dna_cod"){
            R2 <- dna_cod(tab, title, parameters)}
        else if (direction == "cod_dna"){
            R2 <- cod_dna(tab, title, parameters)}
	CI <- cor_ci(r=sqrt(as.numeric(R2)), n=as.numeric(num.obs))
	low <- CI[1]
	high <- CI[2]
	r <- as.character(round(sqrt(as.numeric(R2)),3))
        df <- rbind(df, c(parameters, R2, r, num.obs, low, high))
    }
    print("writing table")
    write.table(df, file= paste0("numdatapoints_", t, "_", direction), sep="\t")
    print("num_datapoints written")
    dev.off()
}
args <- commandArgs(trailingOnly = TRUE)
ExecuteScript(args[1], args[2], args[3])  #"dna_pro", "proj1/Human_Output/DNA_Pro/ran0729/allfilenames"

