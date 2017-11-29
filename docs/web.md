<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pDownloadForm](#pdownloadform)
- [pDownloadGet](#pdownloadget)
- [pDownloadPost](#pdownloadpost)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pDownloadForm

### description
	Download results by submitting a form, supporting pagination.

### input
#### `url`   :
 the URL contains the form  
#### `data`  :
 the data used to fill the form (JSON string or transformed from dict by json.dumps).  
#### `submit`:
 the submit button to submit the form (use Xpath).  
#### `next`  :
 the button for next page (use Xpath)  

### output
#### `outdir:file`:
 The directory saves the results  

### args
#### `interval`:
 seconds to wait between fetching each page. Default: 1  

## pDownloadGet

### description
	Download results by urls.

### input
#### `url`:
 the URLs to download  

## pDownloadPost

### description
	Download results by POST.

### input
#### `url` :
 the URLs to download  
#### `data`:
 the POST data.  
