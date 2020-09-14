"""Script for web.pDownloadGet"""
# pylint: disable=invalid-name,undefined-variable

import requests
from bioprocs.utils import logger

logger.info('Beginning file download with requests')

url = {{i.url | quote}}
resp = requests.get(url)

with open({{o.outfile | quote}}, 'wb') as f:
    f.write(resp.content)

# Retrieve HTTP meta-data
logger.info(str(resp.status_code))
logger.info(str(resp.headers['content-type']))
logger.info(resp.encoding)
