# External data

This folder contains large or licensed datasets that are not distributed with this repository.

## Instructions

Download the datasets listed below and place them in this folder using the exact filenames specified.

---

## Datasets

### 1. Environment Agency biological survey data
- **File(s):** OPEN_DATA_TAXON_INFO.rds; INV_OPEN_DATA_TAXA.rds; INV_OPEN_DATA_SITE.rds; INV_OPEN_DATA_METRICS.rds; FW_Fish_Counts.rds
- **Source:** https://environment.data.gov.uk/ecology/explorer/downloads/
- **Access:** Open Government Licence v3.0
- **Notes:**

Users should:
1. Download compressed bulk download files for Freshwater fish counts (NFPD), Freshwater river macroinvertebrate surveys (Biosys), and Biosys Taxon Info
2. Extract comma-separated values files, import them into R and save them in RDS format in data/external using the original filenames

---

### 2. Environment Agency water quality data
- **File(s):** Comma-separated values files with the naming convention <year>.csv
- **Source:** https://environment.data.gov.uk/dataset/8034d47a-ba8a-4978-aca0-bd5d6e870536
- **Access:** Open Government Licence v3.0
- **Notes:** 

Bulk downloads of these data are no longer available. Users can request the wrangled data from Dr Martin Wilkes ('m.wilkes@essex.ac.uk') or follow the steps below to reconstruct the bulk downloads.

To reconstruct the bulk downloads, users should:
1. Visit the Water Quality Explorer map homepage available via the source above
2. Click Open Filters and under Area Type select "River Basin District"
3. The user must now select each River Basin District in turn
4. For each River Basin District, click View Observations > Open Filters
5. Set dates from 01/01/2002 to 31/12/2002
6. Click Closer Filters then Download
7. For each year, bind the rows of comma-separate values files for each Rievr Basin District and name the file <year>.csv
8. Place the processed comma-searate values files in data/external

---

### 3. HadUK-Grid rainfall data
- **File(s):** rainfall_hadukgrid_uk_5km_*.nc
- **Source:** https://www.metoffice.gov.uk/
- **Access:** Open Government Licence v3.0
- **Notes:** Monthly rainfall data (5 km resolution)

Users should:
1. Download monthly (5 km resolution) rainfall files covering 2002 to 2023 from https://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.1.ceda/5km/rainfall/mon/v20250415
2. Download long-term monthly mean (5 kmn resolution) rainfall files for 1961 to 1990 inclusive and 1991 to 2020 inclusive from https://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.1.ceda/5km/rainfall/mon-30y/v20250415
3. Place all downloaded files in data/external

---

### 4. OS Open Rivers river network
- **File(s):** WatercourseLink.*
- **Source:** Ordnance Survey
- **Access:** Open Government Licence v3.0

Users should:
1. Download the data in ESRI&copy; Shapefile format
2. Extract WatercoureLink.shp and companion files to data/external

---

### 5. Land Cover Map (UKCEH)
- **File(s):** gb2020lcm1km_percentage_target.tif &copy;UK Centre for Ecology & Hydrology (UKCEH)
- **Source:** UK Centre for Ecology & Hydrology
- **Access:** Subject to UKCEH data licence (free for research use; registration may be required)  

Users should:
1. Download the compressed file from https://catalogue.ceh.ac.uk/documents/304e5acd-bea3-46f2-abc0-20bcc174e15b
2. Extract gb2020lcm1km_percentage_target.tif into data/external

---

### 6. Agricultural Sediment Risk (ASR) and Channel Resectioning Index (CRI)
- **File(s):** ASR_CRI.xlsx
- **Source:** Dr Marc Naura
- **Access:** Not licensed for redistribution

Users should request the data from Dr Marc Naura ('Marc.J.Naura@cranfield.ac.uk') and place the file in data/external

---

### 7. AMBER Barrier Atlas
- **File(s):** ASR_CRI.xlsx
- **Source:** https://amber.international/european-barrier-atlas/
- **Access:** Creative Commons Attribution 4.0 International
- **Attribution statement:** Jones, Joshua (2019), “Barrier database (AMBER-GB) for: A Global Assessment of Stream Fragmentation in Great Britain”, Mendeley Data, v1 http://dx.doi.org/10.17632/trsrvw4swg.1

Users should:
1. Download the comma-separated values file from the source
2. The data can be cropped to the UK, if desired
3. Rename the file 'atlas-country-United-Kingdom.csv' and place in data/external

---

### 8. CHESS-SCAPE monthly mean precipitation projections
- **File(s):** rainfall_hadukgrid_uk_5km_*.nc
- **Source:** https://catalogue.ceda.ac.uk/uuid/8194b416cbee482b89e0dfbe17c5786c/
- **Access:** Open Government Licence v3.0
- **Notes:** Future projections of precipitation at 1 km resolution for the United Kingdom 1980-2080 derived from UK Climate Projections 2018

Users should download the following files and place them in data/external
1. https://dap.ceda.ac.uk/badc/deposited2021/chess-scape/data/rcp26_bias-corrected/01/monthly/chess-scape_rcp26_bias-corrected_01_pr_uk_1km_monthly_19801201-20801130.nc?download=1
2. https://dap.ceda.ac.uk/badc/deposited2021/chess-scape/data/rcp60_bias-corrected/01/monthly/chess-scape_rcp60_bias-corrected_01_pr_uk_1km_monthly_19801201-20801130.nc?download=1
3. https://dap.ceda.ac.uk/badc/deposited2021/chess-scape/data/rcp85_bias-corrected/01/monthly/chess-scape_rcp85_bias-corrected_01_pr_uk_1km_monthly_19801201-20801130.nc?download=1

---

## Notes

These datasets are not included in the repository due to file size constraints or licensing restrictions.

All analyses can be reproduced once these datasets are downloaded and placed in this folder.