<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pTxt](#ptxt)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


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
