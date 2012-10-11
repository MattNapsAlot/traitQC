HDF5ReadAllData <- function(in.file,
                            log=FALSE)
{
	require(bat)
    ret.object <- NULL

    file.sum <- HDF5Summary(in.file)

    #
    # Do some checks to make sure the file looks good.
    #
    if(!file.exists(in.file))
        stop(">> The files does not exist\n")
    if(length(file.sum$groupsummary) > 1)
        stop(">> There is too many data groups in the hdf5 file\n")
    if(length(file.sum$groupsummary) == 0)
        stop(">> There is no data groups in the hdf5 file\n")

    fields.to.create <- strsplit(file.sum$datasetsummary,
                                 paste(file.sum$groupsummary,"/",sep=""))

    for(x in 1:length(fields.to.create))
    {
        if(log)
            cat(">> Reading the field: ",file.sum$datasetsummary[x],"\n",sep="")
        ret.object[[ fields.to.create[[x]][2] ]] <- HDF5ReadData(in.file, 
                                                                 file.sum$datasetsummary[x])
    }
    return(ret.object)
}
