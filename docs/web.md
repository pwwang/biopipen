
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
