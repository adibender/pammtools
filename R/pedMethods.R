# @Author: andreas.bender@stat.uni-muenchen.de
# @Date:   2017-03-20 17:57:02
# @Last Modified by:   andreas.bender@stat.uni-muenchen.de
# @Last Modified time: 2017-03-20 19:16:16


#' @importFrom dplyr filter
#' @export
filter_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- filter_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr slice
#' @export
slice_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- slice_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr arrange
#' @export
arrange_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- arrange_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr select
#' @export
select_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- select_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr rename
#' @export
rename_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- rename_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}


#' @importFrom dplyr distinct
#' @export
distinct_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- distinct_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}


#' @importFrom dplyr mutate
#' @export
mutate_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- mutate_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

} 

#' @importFrom dplyr transmute
#' @export
transmute_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- transmute_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}


#' @importFrom dplyr summarise
#' @export
summarise_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- summarise_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}


#' @importFrom dplyr sample_n
#' @export
sample_n_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- sample_n_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

} 

#' @importFrom dplyr sample_frac
#' @export
sample_frac_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- sample_frac_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}