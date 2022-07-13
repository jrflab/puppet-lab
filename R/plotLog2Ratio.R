'plotLog2Ratio' <- function(facets, num_levels = NA, purity = NA, ploidy = NA)
{
	log2_ratio = facets$jointseg %>%
		     dplyr::as_tibble() %>%
		     dplyr::select(Chromosome = chrom,
				   Position = maploc,
				   Log2Ratio = cnlr)
	segmented = facets$cncf$cncf %>%
		    dplyr::as_tibble() %>%
		    dplyr::select(Chromosome = chrom,
				  Start = start,
				  End = end,
				  Log2Ratio = cnlr.median,
				  N = num.mark)
	seg = puppetlabprune_(x = seg %>% data.frame(), n)
		seg = seg %>%
		      dplyr::as_tibble() %>%
		      dplyr::mutate(n = cumsum(n))
		purity = ifelse(is.na(purity), 1, purity)
		ploidy = ifelse(is.na(ploidy), 2, ploidy)
		data(CytoBand)
   		end = dplyr::as_tibble(CytoBand) %>%
		      dplyr::filter(Chromosome %in% 1:23) %>%
		      dplyr::group_by(Chromosome) %>%
		      dplyr::summarize(end = max(End)) %>%
		      .[["end"]]
   		end = cumsum(end)
   		start = c(1, end[1:22]+1)
   		CytoBand = dplyr::tibble(Chromosome = 1:23,
					 Start = start,
					 End = end) %>%
			   dplyr::mutate(Centromere = .5*(Start + End))
   		index = NULL
   		for (i in 1:23) {
   			index = c(index, seq(from = CytoBand[i, "Start"][[1]], to=CytoBand[i, "End"][[1]], length=sum(cna$chrom==i)))
   		}
		plot_ = dplyr::tibble(x = index,
				      y = cna$log2,
				      title = gsub(pattern = "_", replacement = " ", x = opt$sample_name)) %>%
			ggplot(mapping = aes(x = x, y = y)) +
			geom_point(stat = "identity", shape = 1, size = 1, color = "black") +
			xlab("\n\nChromosome\n") +
                	ylab(expression(Log[2]~"Ratio")) +
                	scale_x_continuous(breaks = CytoBand[["Centromere"]],
					   labels = c(1:22, "X"),
					   minor_breaks = c(CytoBand[["Start"]], CytoBand[["End"]][23])) +
                	scale_y_continuous(limits = c(-5, 5)) +
                	theme(plot.margin = unit(c(1, 1, 3, 1), "lines")) +
			facet_wrap(~title)
		
		plot(index, cna$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=FALSE, xlab="", ylab="", main="", ylim=c(-4.5,4.5))
 		for (j in 1:nrow(seg)) {
 			if (j == 1) {
 				lines(x=c(1, index[seg[j,"n"]]), y=rep(seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
 			} else {
 				lines(x=c(index[seg[j-1,"n"]], index[seg[j,"n"]]), y=rep(seg[j,"log2"],2), lty=1, lwd=2.75, col="red")
 			}
  		}
  		axis(side=1, at=c(CytoBand[,"start"],CytoBand[nrow(CytoBand),"end"]), labels=rep("", nrow(CytoBand)+1), tcl=.5)
		axis(side=1, at=apply(CytoBand[,c("start", "end"),drop=FALSE], 1, mean), labels=c(1:22, "X"), lwd=0, lwd.ticks=1, tcl=-.35)
  		axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1, lwd = 1.15, lwd.ticks = 0.95, line = -.5)
		mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
		points(c(0-.05*max(index),max(index)+.01*max(index)), c(0,0), type="l", col="black", lwd=1)
		title(main = paste0(title, " | alpha = ", signif(purity, 3), " | psi = ", signif(ploidy, 3)), cex.main=.75, font.main=1)
	}