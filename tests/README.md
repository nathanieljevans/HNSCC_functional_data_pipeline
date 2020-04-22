# HNSCC pipeline Test Plan 

- [HNSCC pipeline Test Plan](#hnscc-pipeline-test-plan)
  - [overview](#overview)
- [Test 001](#test-001)
- [Test 002](#test-002)
- [Test 003](#test-003)


## overview 

The [requirements document](../docs/design_documents/requirementsDocument.md) specifies the system requirements that must be met for sufficiency in this project. 

# Test 001 

This requirement is meant to test the accurate mapping of data from photospectrometer output data and linking to inhibitor, lab_id and optical_density using the provided plate maps. 

We'll do this by setting up a test plate map and dataset, where mapped data fidelity is easily tested.

Execute this test by navigating to `./test001/` and running: 

```bash
$ python run_test001.py
```

output should look like: 

```
path lab_id=tst01-norm=blank490-plate_version_id=test001_platemap-note=test001.xlsx
---------------------------------------------------------
please double check file name parsing is accurate
file name: lab_id=tst01-norm=blank490-plate_version_id=test001_platemap-note=test001.xlsx
lab_id: tst01
norm: blank490
version_id: test001_platemap
notes: test001
---------------------------------------------------------

--------------------------------------------------------------------------
This is the message log for:
        lab_id: tst01
        version_id: test001_platemap
        notes: test001
--------------------------------------------------------------------------
Data has already been normalized by positive controls: True
This assay has 1 plates
mapping data... [test001_platemap]
mapping complete.
         mapped data shape: (384, 11)
         column names: [['lab_id', 'norm_type', 'plate_num', 'plate_row', 'assay_version_id', 'note', 'plate_col', 'optical_density', 'conc', 'inhibitor', 'map_version_id']]

-------------------------------------------------------------------------
processing: lab_id=tst01-norm=blank490-plate_version_id=test001_platemap-note=test001.xlsx
total number of rows in data after mapping: 384
number of unique inhibitors: 384
platemap being used: ./HNSCC_plate_map-version_id=test001_platemap.xlsx
-------------------------------------------------------------------------

    inhibitor   conc optical_density
0         I-1    C-1            OD-1
1        I-25   C-25           OD-25
2        I-49   C-49           OD-49
3        I-73   C-73           OD-73
4        I-97   C-97           OD-97
..        ...    ...             ...
379     I-288  C-288          OD-288
380     I-312  C-312          OD-312
381     I-336  C-336          OD-336
382     I-360  C-360          OD-360
383     I-384  C-384          OD-384

[384 rows x 3 columns]
testing mapping...
mapping is verified correct by test001

```

# Test 002 

# Test 003 