# DIA-SpecLib_Result_Summarizer
Extract the identified peptides and proteins from spectral library search of DIA proteomic data in Skyline

## How to use DIA-SpecLib_Result_Summarizer in Windows/Linux
```
>python extract_skyline_result.py -h
usage: extract_skyline_result.py [-h] -ip [-ip ...]

Extract and summarize Skyline results from DIA data analysis against spectral
library

positional arguments:
  -ip         Skyline output from DIA-Spectral Library search

optional arguments:
  -h, --help  show this help message and exit

```

# iq input generator from Skyline search results

## How to use iq_input_gen.py in Windows/Linux
```
>python iq_input_gen.py -h
usage: iq_input_gen.py [-h] -i [-i ...] -c [-c ...]

Generate a input matrix for the calculation of protein abundance using iq R
scrip

positional arguments:
  -i          Skyline output from DIA-Spectral Library search
  -c          A .txt file consisting the details of condition specific raw
              files used for the search

optional arguments:
  -h, --help  show this help message and exit
 ```
