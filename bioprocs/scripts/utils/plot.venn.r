if (!exists('plotVenn')) {
	plotVenn = function(mat, filename, params) {
		library(VennDiagram)
		nc   = ncol(mat)
		cats = colnames(mat)
		csum = colSums(mat)
		rsum = rowSums(mat)
		png (filename, res=300, width=2000, height=2000)
		if (nc == 1) {
			v = do.call(draw.single.venn, c(list(area = csum[1,1], category = cats[1]), params))
		} else if (nc == 2) {
			cross.area = length(which(rsum == nc))
			v = do.call(draw.pairwise.venn, c(list(area1 = csum[1], area2 = csum[2], cross.area = cross.area, category = cats), params))
		} else if (nc == 3) {
			n12  = length(which(rowSums(mat[, c(1, 2)]) == 2))
			n23  = length(which(rowSums(mat[, c(2, 3)]) == 2))
			n13  = length(which(rowSums(mat[, c(1, 3)]) == 2))
			n123 = length(which(rsum == nc))
			v = do.call(draw.triple.venn, c(list(area1 = csum[1], area2 = csum[2], area3 = csum[3], n12 = n12, n23 = n23, n13 = n13, n123 = n123, category = cats), params))
		} else if (nc == 4) {
			n12  = length(which(rowSums(mat[, c(1, 2)]) == 2))
			n13  = length(which(rowSums(mat[, c(1, 3)]) == 2))
			n14  = length(which(rowSums(mat[, c(1, 4)]) == 2))
			n23  = length(which(rowSums(mat[, c(2, 3)]) == 2))
			n24  = length(which(rowSums(mat[, c(2, 4)]) == 2))
			n34  = length(which(rowSums(mat[, c(3, 4)]) == 2))
			n123 = length(which(rowSums(mat[, c(1,2,3)]) == 3))
			n124 = length(which(rowSums(mat[, c(1,2,4)]) == 3))
			n134 = length(which(rowSums(mat[, c(1,3,4)]) == 3))
			n234 = length(which(rowSums(mat[, c(2,3,4)]) == 3))
			n1234 = length(which(rsum == nc))
			v = do.call(draw.quad.venn, c(list(
				area1 = csum[1], 
				area2 = csum[2], 
				area3 = csum[3], 
				area4 = csum[4], 
				n12   = n12, 
				n13   = n13, 
				n14   = n14, 
				n23   = n23, 
				n24   = n24, 
				n34   = n34, 
				n123  = n123, 
				n124  = n124, 
				n134  = n134, 
				n234  = n234, 
				n1234 = n1234, 
				category = cats), params))
		} else if (nc == 5) {
			n12    = length(which(rowSums(mat[, c(1, 2)]) == 2))
			n13    = length(which(rowSums(mat[, c(1, 3)]) == 2))
			n14    = length(which(rowSums(mat[, c(1, 4)]) == 2))
			n15    = length(which(rowSums(mat[, c(1, 5)]) == 2))
			n23    = length(which(rowSums(mat[, c(2, 3)]) == 2))
			n24    = length(which(rowSums(mat[, c(2, 4)]) == 2))
			n25    = length(which(rowSums(mat[, c(2, 5)]) == 2))
			n34    = length(which(rowSums(mat[, c(3, 4)]) == 2))
			n35    = length(which(rowSums(mat[, c(3, 5)]) == 2))
			n45    = length(which(rowSums(mat[, c(4, 5)]) == 2))
			n123   = length(which(rowSums(mat[, c(1, 2, 3)]) == 3))
			n124   = length(which(rowSums(mat[, c(1, 2, 4)]) == 3))
			n125   = length(which(rowSums(mat[, c(1, 2, 5)]) == 3))
			n134   = length(which(rowSums(mat[, c(1, 3, 4)]) == 3))
			n135   = length(which(rowSums(mat[, c(1, 3, 5)]) == 3))
			n145   = length(which(rowSums(mat[, c(1, 4, 5)]) == 3))
			n234   = length(which(rowSums(mat[, c(2, 3, 4)]) == 3))
			n235   = length(which(rowSums(mat[, c(2, 3, 5)]) == 3))
			n245   = length(which(rowSums(mat[, c(2, 4, 5)]) == 3))
			n345   = length(which(rowSums(mat[, c(3, 4, 5)]) == 3))
			n1234  = length(which(rowSums(mat[, c(1, 2, 3, 4)]) == 4))
			n1235  = length(which(rowSums(mat[, c(1, 2, 3, 5)]) == 4))
			n1245  = length(which(rowSums(mat[, c(1, 2, 4, 5)]) == 4))
			n1345  = length(which(rowSums(mat[, c(1, 3, 4, 5)]) == 4))
			n2345  = length(which(rowSums(mat[, c(2, 3, 4, 5)]) == 4))
			n12345 = length(which(rsum == nc))
			v = do.call(draw.quintuple.venn, c(list(
				area1  = csum[1], 
				area2  = csum[2], 
				area3  = csum[3], 
				area4  = csum[4], 
				area5  = csum[5], 
				n12    = n12, 
				n13    = n13, 
				n14    = n14, 
				n15    = n15, 
				n23    = n23, 
				n24    = n24, 
				n25    = n25, 
				n34    = n34, 
				n35    = n35, 
				n45    = n45, 
				n123   = n123, 
				n124   = n124, 
				n125   = n125, 
				n134   = n134, 
				n135   = n135, 
				n145   = n145, 
				n234   = n234, 
				n235   = n235, 
				n245   = n245, 
				n345   = n345, 
				n1234  = n1234, 
				n1235  = n1235, 
				n1245  = n1245, 
				n1345  = n1345, 
				n2345  = n2345, 
				n12345 = n12345, 
				category = cats), params))
		} else {
			stop('Too many columns to draw a venn diagram, please use UpSetR.')
		}
		print (v)
		dev.off()
	}
}