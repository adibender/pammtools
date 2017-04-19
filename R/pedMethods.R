# @Author: andreas.bender@stat.uni-muenchen.de
# @Date:   2017-03-20 17:57:02
# @Last Modified by:   andreas.bender@stat.uni-muenchen.de
# @Last Modified time: 2017-04-19 09:54:30


#' @importFrom dplyr filter
#' @export
filter.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- filter(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr filter_
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
slice.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- slice(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr slice_
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
arrange.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- arrange(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr arrange_
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
mutate.ped <- function(ped, keep.attributes=TRUE, ...) {

	if(keep.attributes) {
		ped.attr   <- attributes(ped)
		attr.names <- setdiff(names(ped.attr), c("class", "row.names", "names"))
	}
	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- mutate(ped, ...)
	class(ped) <- c(ped.class, class(ped))
	if(keep.attributes) {
		attributes(ped) <- c(attributes(ped), ped.attr[attr.names])
	}

	return(ped)

} 

#' @importFrom dplyr mutate_
#' @export
mutate_.ped <- function(ped, keep.attributes=TRUE, ...) {

	if(keep.attributes) {
		ped.attr   <- attributes(ped)
		attr.names <- setdiff(names(ped.attr), c("class", "row.names", "names"))
	}
	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- mutate_(ped, ...)
	class(ped) <- c(ped.class, class(ped))
	if(keep.attributes) {
		attributes(ped) <- c(attributes(ped), ped.attr[attr.names])
	}

	return(ped)

} 

#' @importFrom dplyr transmute
#' @export
transmute.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- transmute(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr transmute_
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


#' @import dplyr 
#' @export
sample_n.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- sample_n(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

} 

#' @import dplyr
#' @export
sample_frac.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- sample_frac(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr left_join
#' @export
left_join.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- left_join(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}


#' @importFrom dplyr group_by
#' @export
group_by_.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- group_by_(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}

#' @importFrom dplyr ungroup
#' @export
ungroup.ped <- function(ped, ...) {

	ped.class  <- class(ped)[1]
	class(ped) <- class(ped)[-1]
	ped        <- ungroup(ped, ...)
	class(ped) <- c(ped.class, class(ped))

	return(ped)

}