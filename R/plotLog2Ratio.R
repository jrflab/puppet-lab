'plotlog2ratio' <- function(facets, num_levels = 10, purity = 1, ploidy = 2)
{
	log2ratio = facets$jointseg %>%
		    dplyr::as_tibble() %>%
		    dplyr::mutate(maploc = 1:nrow(.)) %>%
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
	pruned = puppetlab::prune(x = segmented, n = num_levels)
	pruned = pruned %>%
		 dplyr::mutate(N = cumsum(N)) %>%
		 dplyr::mutate(Start = c(1, 1+N[1:(nrow(.)-1)])) %>%
		 dplyr::mutate(End = N)
	cytoband = log2ratio %>%
		   dplyr::group_by(Chromosome) %>%
		   dplyr::summarize(Start = min(Position),
				    End = max(Position)) %>%
		   dplyr::mutate(Centromere = .5*(Start + End))
	
	plot_ = log2ratio %>%
		ggplot(mapping = aes(x = Position, y = Log2Ratio)) +
		geom_point(stat = "identity", shape = 1, size = 1, color = "black") +
		geom_segment(mapping = aes(x = Start, y = Log2Ratio, xend = End, yend = Log2Ratio), data = pruned, color = "red", size = 1.5, inherit.aes = FALSE) +
		xlab("\n\nChromosome\n") +
                ylab(expression(Log[2]~"Ratio")) +
                scale_x_continuous(breaks = cytoband[["Centromere"]],
				   labels = c(1:22, "X"),
				   minor_breaks = c(cytoband[["Start"]], cytoband[["End"]][23])) +
		scale_y_continuous(limits = c(-5, 5),
				   sec.axis = sec_axis(trans = ~(((2^(.))*((purity*ploidy) + (1-purity)*2)) - 2*(1-purity))/purity,
				  		       breaks = c(1, 2, 3, 4, 5, 10, 15, 20, 25, 40, 50)),
				   		       name = expression(Log[2]~"Ratio")) +
                theme(plot.margin = unit(c(1, 1, 3, 2), "lines"))
		
	return(invisible(plot_))
}
