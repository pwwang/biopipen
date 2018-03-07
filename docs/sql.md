# sql
<!-- toc -->
{% raw %}

## pCreateTable

### description
Create tables in the database

### input
#### `dsn`:: The dsn to connect to the database  
	- currently support `sqlite:file=...`
#### `schema:file`:: The schema file  
	- could be a pure schema file:
	```
	Field	Type	Statement
	ID	INT	PRIMARY KEY
	...
	```
	- or a data file with header

### output
#### `dsn`:: The dsn  

## pImportData

### description
Create tables and import the data

### input
#### `dsn`:: The dsn to connect to the database  
	- currently support `sqlite:file=...`
#### `datafile:file`:: The schema file  
	- must have header

### output
#### `dsn`:: The dsn  

## pUpdateTable

### description
Update table using sql.

### input
#### `dsn`:: The dsn to connect to the database  
	- currently support `sqlite:file=...`

### output
#### `dsn`:: The dsn  

## pSelectTable

### description
Select data from table and dump it.

### input
#### `dsn`:: The dsn to connect to the database  
	- currently support `sqlite:file=...`

### output
#### `outfile:file`:: The dumped file  
{% endraw %}
