
# This **readme** documents the use of `error_annotation.xlsx` 

This can be used to manually annotate plate specific errors that occur or are noticed during the bench protocols. 

The file `error_annotation.xlsx` has seven columns: 

--- 

> column name: **assay_name**   
> description: This is the assay identifier and must exactly match the filename that it references (exclude the file type: `.xlsx`).  
> allowed values: String; any ASCII characters  

--- 

> column name: **plate_number**   
> description: The plate number identifier marking the error. Refer to plate map for specifics.    
> allowed values: Integer; 1-6

--- 

This format allows for contiguous regions on the plate to be annotated easily, for non-contiguous annotations, please use separate annotations. **To indicate only a single well**, input the same value for XXX_from and XXX_to (eg. `col_from` = 3, `col_to` = 3).

> column name: **col_from**   
> description: Start column (inclusive) 
> allowed values: Integer; 1-24   

> column name: **col_to**   
> description: End column (inclusive)  
> allowed values: Integer; 1-24  

> column name: **row_from**  
> description: Start row (inclusive)  
> allowed values: String; A-P   

> column name: **row_to**   
> description: End row (inclusive)    
> allowed values: String; A-P     

--- 

> column name: **annotations**   
> description: comments and description of the error that occured; reason for removal/flag     
> allowed values: String; any ASCII characters    

--- 

> column name: **remove**   
> description: True if data should be removed/filtered; False will merely falg the data and add the annotation to the `comments` column of the output file.   
> allowed values: Boolean; True/False

--- 

## Example Annotation

An example annotation observation in the `error_annotation.xlsx` file: 

| assay_name                                                             | plate_number | col_from | col_to | row_from | row_to | annotations                                            | remove |
|------------------------------------------------------------------------|--------------|----------|--------|----------|--------|--------------------------------------------------------|--------|
| lab_id=10004-norm=Blank490-plate_version_id=OHSU_HNSCC_derm002-note=NA | 1            | 9        | 14     | A        | A      | EXAMPLE: raw data looks entirely dead, maybe an issue. | FALSE  |

This indicates an error on plate 1, from column 9 to 14 (inclusive!) in the first row (A), the annotation describes the issue and we see that the data is not supposed to be removed. **Please Note**: It is important that the assay_name match the file exactly it's filename (excluding `.xlsx`). 


## To add annotation to `HNSCC_all_functional_data.csv` file... 

To add these annotations to the relevant dose-response assay, run the file: 

```bash
$ python add_manual_annotations.py
```

Optional arguments: 
```bash
$ python add_manual_annotations.py --help
usage: add_manual_annotations.py [-h] [--no-verbose] [--annotation-path <pth>]
                                 [--data-path <pth>] [--assayID-path <pth>]
                                 [--output-path <pth>]

This script adds annotations from `../annotations/error_annotation.xlsx` to
the summary output file `HNSCC_all_functional_data.csv`. For more information,
please refer to `../annotations/readme.md`

optional arguments:
  -h, --help            show this help message and exit
  --no-verbose          whether to print output messages to console
  --annotation-path <pth>
                        path to annotations
  --data-path <pth>     path to all data output file
  --assayID-path <pth>  path to .assayid file
  --output-path <pth>   path to output file
```

You should see the output similar to... 

(no errors)
```bash
annotation file path: ../annotation/error_annotation.xlsx
output data file path: ../output/HNSCC_all_functional_data.csv
adding annotations... 18 annotations added.
saving data...data saved.
script completed.
```

(with error)
```bash
$ python add_manual_annotations.py
annotation file path: ../annotation/error_annotation.xlsx
output data file path: ../output/HNSCC_all_functional_data.csv
adding annotations...
Assay ID not found: lab_id=99999;A-MISTAKE;-norm=Blank490-plate_version_id=OHSU_HNSCC_dermMISTAKE-note=NA
Please double check that the assay name is correct and has been processed.

-----------------------------------------------------------
failed to add annotation:
assay_name      lab_id=99999;A-MISTAKE;-norm=Blank490-plate_ve...
plate_number                                                    2
col_from                                                        9
col_to                                                         14
row_from                                                        A
row_to                                                          C
annotations     EXAMPLE: raw data looks entirely dead, maybe a...
remove                                                      False
Name: 1, dtype: object

Exception traceback:
No panel ID found.
-----------------------------------------------------------
18 annotations added.
saving data...data saved.
script completed.

```

**Please Note**: The script will attempt to add all annotations, errors will be displayed in output but not halt file execution. For this reason, carefully check that the script output does not indicate any failures. If there are failures, follow exception description to debug... most likely a user error. 

An annotated version of `HNSCC_all_functional_data.csv` will be saved to ./outputs/ as `HNSCC_all_functional_data-annotated.csv` unless otherwise specified by the `--output-path` argument.

The annotated output file will have two new columns, `manual_flag` which will be True for all columns that have annotations and `manual_remove` which will be True if the annotation specifies `remove` = True. 