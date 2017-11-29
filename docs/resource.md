{% raw %}

## pTxt

### description
	Download CSV format files.

### input
#### `in`:
 The name of the resource  

### output
#### `outfile:file`:
 The output file  

### args
#### `cols`:
      Select the columns to keep. Default: '' (all cols)  
#### `rowfilter`:
 Filter rows. For example, to filter out rows not start with 'Chr':  
		- `"lambda x: not x[0].startswith('Chr')"`
		- Note that rowfilter applied before cols filter.
#### `urls`:
      Available resources and their urls.  
#### `gz`:
        Whether to gzip the output file.  
{% endraw %}
