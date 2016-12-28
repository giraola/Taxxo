#' getOS
#'
#' Helper function that identifies the Operating System.
#' @keywords OS
#' @export
#' @examples
#' .getOS()

.getOS<-function()

	{
  
  	sysinf<-Sys.info()
  
  	if (!is.null(sysinf)) {
    	
    	os<-sysinf['sysname']
  	}
    
    if (os=='Darwin') {
    
		os<-"darwin"

	} else {

		os<-.Platform$OS.type
	}
    
    if (grepl("^darwin",R.version$os)) {

		os<-"darwin"
    }
	
	if (grepl("linux-gnu",R.version$os)) {
	
		os<-"linux"
  	}

	tolower(os)
}
