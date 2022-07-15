'plot_cncf' <- function(facets, purity = 1, ploidy = 2)
{
	log2ratio = facets$jointseg %>%
		    dplyr::as_tibble() %>%
		    dplyr::mutate(maploc = 1:nrow(.)) %>%
		    dplyr::select(Chromosome = chrom,
				  Position = maploc,
				  Log2Ratio = cnlr,
				  LogOR = valor)
	segmented = facets$cncf$cncf %>%
		    dplyr::as_tibble() %>%
		    dplyr::select(Chromosome = chrom,
				  Start = start,
				  End = end,
				  Log2Ratio = cnlr.median,
				  LogOR = mafR,
				  TCN = tcn.em,
				  LCN = lcn.em,
				  CF = cf.em,
				  N = num.mark)
	pruned = segmented %>%
		 dplyr::mutate(N = cumsum(N)) %>%
		 dplyr::mutate(Start = c(1, 1+N[1:(nrow(.)-1)])) %>%
		 dplyr::mutate(End = N) %>%
		 dplyr::mutate(LogOR = sqrt(abs(LogOR)))
	cytoband = log2ratio %>%
		   dplyr::group_by(Chromosome) %>%
		   dplyr::summarize(Start = min(Position),
				    End = max(Position)) %>%
		   dplyr::mutate(Centromere = .5*(Start + End))
	
	pa = log2ratio %>%
	     ggplot(mapping = aes(x = Position, y = Log2Ratio)) +
	     geom_point(stat = "identity", shape = 1, size = .75, color = "grey10") +
	     geom_segment(mapping = aes(x = Start, y = Log2Ratio, xend = End, yend = Log2Ratio), data = pruned, color = "red", size = 1.5, inherit.aes = FALSE) +
	     xlab("") +
	     ylab(expression(Log[2]~"Ratio")) +
	     scale_x_continuous(breaks = cytoband[["Centromere"]],
				labels = rep("", 23),
				minor_breaks = c(cytoband[["Start"]], cytoband[["End"]][23])) +
	     scale_y_continuous(limits = c(-5, 5)) +
	     theme(plot.margin = unit(c(2, 1, 0, 2), "lines"),
		   axis.text.x = element_blank(),
		   axis.ticks.x = element_blank(),
		   axis.text.y = element_text(size = 8))
	
	pb = log2ratio %>%
	     ggplot(mapping = aes(x = Position, y = abs(LogOR))) +
	     geom_point(stat = "identity", shape = 1, size = .75, color = "grey10") +
	     geom_point(stat = "identity", mapping = aes(x = Position, y = -abs(LogOR)), shape = 1, size = .75, color = "grey10") +
	     geom_segment(mapping = aes(x = Start, y = LogOR, xend = End, yend = LogOR), data = pruned, color = "red", size = 1.5, inherit.aes = FALSE) +
	     geom_segment(mapping = aes(x = Start, y = -LogOR, xend = End, yend = -LogOR), data = pruned, color = "red", size = 1.5, inherit.aes = FALSE) +
	     xlab("") +
	     ylab(expression(Log~"OR")) +
	     scale_x_continuous(breaks = cytoband[["Centromere"]],
				labels = rep("", 23),
				minor_breaks = c(cytoband[["Start"]], cytoband[["End"]][23])) +
	     scale_y_continuous(limits = c(-5, 5)) +
	     theme(plot.margin = unit(c(0, 1, 0, 2), "lines"),
		   axis.text.x = element_blank(),
		   axis.ticks.x = element_blank(),
		   axis.text.y = element_text(size = 8))
	
	pc = log2ratio %>%
	     ggplot(mapping = aes(x = Position, y = abs(LogOR))) +
	     geom_point(stat = "identity", shape = 1, size = .75, color = NA) +
	     geom_segment(mapping = aes(x = Start, y = LCN, xend = End, yend = LCN), data = pruned, color = "red", size = 1.5, inherit.aes = FALSE) +
     	     geom_segment(mapping = aes(x = Start, y = TCN, xend = End, yend = TCN), data = pruned, color = "grey10", size = 1.5, inherit.aes = FALSE) +
	     xlab("") +
	     ylab(expression("Copy number")) +
	     scale_x_continuous(breaks = cytoband[["Centromere"]],
				labels = rep("", 23),
				minor_breaks = c(cytoband[["Start"]], cytoband[["End"]][23])) +
	     scale_y_continuous(limits = c(0, 5),
			        labels = paste0("     ", 0:5)) +
	     theme(plot.margin = unit(c(0, 1, 1, 2), "lines"),
		   axis.text.x = element_blank(),
		   axis.ticks.x = element_blank(),
		   axis.text.y = element_text(size = 8))
	
	pd = log2ratio %>%
	     ggplot(mapping = aes(x = Position, y = abs(LogOR))) +
	     geom_point(stat = "identity", shape = 1, size = .75, color = NA) +
	     geom_segment(mapping = aes(x = Start, y = 1, xend = End, yend = 1, color = CF), data = pruned, size = 3.5, inherit.aes = FALSE) +
	     xlab("") +
	     ylab(expression(" ")) +
	     scale_colour_viridis_c() +
	     scale_x_continuous(breaks = cytoband[["Centromere"]],
				labels = c(1:22, "X"),
				minor_breaks = c(cytoband[["Start"]], cytoband[["End"]][23])) +
	     scale_y_continuous(limits = c(0.5, 1.5),
				breaks = c(0, 1.5),
			        labels = paste0("    ", c(" ", " "))) +
	     theme(plot.margin = unit(c(-1.5, 1, 2, 3.95), "lines"),
		   axis.text.y = element_blank(),
		   axis.ticks.y = element_blank()) +
	     guides(color = FALSE)
	
	
	
	return(invisible(list(pa, pb, pc, pd)))
	
}
