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

### args
#### `intype`:: The input file schema file or a data file. Default: `schema`  
#### `drop`::  Force creating the table (drop the pre-existing table)  
#### `delimit`::The delimit of input file. Default: `\\t`  

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

### args
#### `delimit`::The delimit of input file. Default: `\\t`  

## pUpdateTable

### description
Update table using sql.

### input
#### `dsn`:: The dsn to connect to the database  
	- currently support `sqlite:file=...`

### output
#### `dsn`:: The dsn  

### args
#### `sql`:: The sql to update the table (list)  

## pSelectTable

### description
Select data from table and dump it.

### input
#### `dsn`:: The dsn to connect to the database  
	- currently support `sqlite:file=...`

### output
#### `outfile:file`:: The dumped file  

### args
#### `sql`:: The sql to select data from the table (list)  
{% endraw %}
