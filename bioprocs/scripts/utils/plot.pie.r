require ('ggplot2')
if (!exists('plotPie')) {
	plotPie = function (m, filename, ggs = list(), devpars = list(res=300, width=2000, height=2000)) {
		percent = function(x) paste0(format(round(x*100, 1), nsmall = 1), "%")
		do.call(png, c(list(filename=filename), devpars))
		N         = nrow(m)
		Group     = colnames(m)
		oValue    = if (N == 1) m[1,] else colSums(m)
		oValue    = as.vector(as.matrix(oValue))
		GroupName = paste0(Group, " (", oValue, ")")
		Value     = 100*oValue/sum(oValue)
		Labels    = percent(Value / 100)
		df        = data.frame (Group = Group, Value = Value)
		df$Group  = factor(df$Group, levels = df$Group)
		#if (length(Group) == 2) { 
		#	ytext = 100 - cumsum(rev(Value)) + rev(Value) / 2
		#} else {
			ytext = cumsum(rev(Value)) - rev(Value) / 2
		#}
		pie =   ggplot(df, aes(x="", y=Value, fill=Group)) +
				geom_bar(width = 1, stat = "identity") +
				geom_text(aes(y = ytext), label = rev(Labels)) +
				coord_polar("y") +
				scale_fill_discrete(labels = GroupName, breaks = Group)

		if (length(ggs) > 0) {
			for (i in 1:length(ggs)) {
				pie = pie + ggs[[i]]
			}
		}
		print(pie)
		dev.off()
	}
}