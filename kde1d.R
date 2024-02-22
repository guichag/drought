### INSTALL R PACKAGES WITH RPY2  ###
# from rpy2.robjects.packages import importr
# utils = importr('utils')
# utils.install_packages('nomdupackage')


# Define function to fit KDE and get corresponding CDF

require('kde1d')


make_kde1d <- function(values){
    out = kde1d(values)
    return(out)
}


make_pkde1d <- function(values, tmp){
    out = pkde1d(values, tmp)
    return(out)
}
