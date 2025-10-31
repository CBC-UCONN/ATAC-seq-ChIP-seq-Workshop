
[Day 1 Slides](https://CBC-UCONN.github.io/ATAC-seq-ChIP-seq-Workshop/slides/day1.html)

[Day 2 Slides](https://CBC-UCONN.github.io/ATAC-seq-ChIP-seq-Workshop/slides/day2.html)

[Mamba Installation](https://CBC-UCONN.github.io/ATAC-seq-ChIP-seq-Workshop/slides/mamba.html)

[Day 3 Slides](https://CBC-UCONN.github.io/ATAC-seq-ChIP-seq-Workshop/slides/day3.html)


## Getting Workshop Code
Code for the workshop can be downloaded from https://github.com/CBC-UCONN/Single-Cell-Transcriptomics/releases.

Copy the url for the version you wish to download and use `wget` or `curl` to download it.

To download the most recent version:

```bash
wget https://github.com/CBC-UCONN/Single-Cell-Transcriptomics/archive/refs/tags/2025.10.zip
```
Once downloaded, unzip the file:

```bash
unzip 2025.10.zip
``` 

On the Xanadu cluster, the results directory can be populated with pre-computed results using the `mock_run.sh` command as follows:
```bash 
./mock_run.sh <script name>.sh
```
