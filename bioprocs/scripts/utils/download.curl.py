if 'downloadCurl' not in vars() or not callable(downloadCurl):
	def downloadCurl(url, outfile, username = '', password = '', curl = 'curl'):
		"curl -Lfv -o filename.zip -u EMAIL:PASSWORD https://www.drugbank.ca/releases/5-0-7/downloads/target-approved-uniprot-links"
		cmd  = curl + ' -L '
		if username and password:
			cmd += '-u {}:{}'.format(username, password)
		elif username:
			cmd += '-u {}'.format(username)
		cmd += ' -o "{}" "{}"'.format(outfile, url)
		runcmd (cmd)