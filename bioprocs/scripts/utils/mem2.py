if 'mem2' not in vars() or not callable (mem2):
	
	def _autoUnit(num):
		if num % (1024 * 1024) == 0:
			return num / (1024*1024), 'G'
		elif num % 1024 == 0:
			return num / 1024, 'M'
		else:
			return num, 'K'

	# unit = 'auto'/'G'/'M'/'K'/'java'
	def mem2 (mem, unit = 'auto'):
		mem   = str(mem)
		ounit = mem[-1].upper()
		if ounit == 'G':
			num   = int(mem[:-1]) * 1024 * 1024
		elif ounit == 'M':
			num   = int(mem[:-1]) * 1024
		elif ounit == 'K':
			num   = int(mem)
		elif not ounit.isdigit():
			raise ValueError('Unknown memory unit: ' + ounit)
		else:
			num   = int(mem)
		
		unit = unit.upper()
		retn, retu = _autoUnit(num)
		if unit == 'G':
			if retu == 'M':   retn /= 1024.
			elif retu == 'K': retn /= (1024. * 1024.)
			retu = 'G'
		elif unit == 'M':
			if retu == 'G':   retn *= 1024
			elif retu == 'K': ret /= 1024.
			retu = 'M'
		elif unit == 'K':
			if retu == 'G':   retn *= 1024 * 1024
			elif retu == 'M': retn *= 1024 
			retu = 'K'
		
		if unit == 'JAVA':
			xmx = "-Xmx" + str(retn) + retu
			n, u = _autoUnit(num / 8)
			return '-Xms' + str(n) + u + ' ' + xmx
		else:
			return str(retn) + retu
