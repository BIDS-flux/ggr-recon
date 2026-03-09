<div align="center"><img src="ggr-recon.png" width="50%" height="50%"></div>

<div align="center"><a href="https://scholar.harvard.edu/files/suiyao/files/sui_miccai_2019.pdf">Paper 1</a> | <a href="https://scholar.harvard.edu/files/suiyao/files/sui_deepgg_miccai_2020.pdf">Paper 2</a></div>

# GGR-recon
A deconvolution-based MRI super-resolution reconstruction method with a gradient guidance regularization (GGR).

The reconstruction comprises two steps: 1) ***preprocessing*** and
2) ***deconvolution***, corresponding to the two python scripts,
respectively, ***preprocess.py*** and ***recon.py***. In the
preprocessing step, the algorithm deals with image alignment for motion
compensation, computes automatically the geometric properties of the
high-res reconstruction, creates the filters (for slice profiles and
downsamplings) used in the deconvolution step, and offers a gradient
guidance reference for the regularization of the deconvolution. The
preporcessing also provides a resampling mode, which is usually used
to determine the geometric properties of the high-res reconstruction
manually. In the deconvolution step, the regularization is created and
then the regularized deconvolution is performed in the Fourier domain.

## Dependencies
### Docker mode
- **Docker**: https://www.docker.com/

### Python mode
- **NumPy**: https://numpy.org/
- **Scipy**: https://www.scipy.org/
- **SimpleITK**: https://simpleitk.org/
- **Rich**: https://rich.readthedocs.io/en/stable/introduction.html
- **pybids**: https://bids-standard.github.io/pybids/
- **CRKIT**: http://crl.med.harvard.edu/software/

## Getting started
GGR-recon can be run in either ***docker*** or ***python*** mode. We *strongly* recommend using GGR-recon in the ***docker*** mode, as the issues about the environment configuration and version conflicts can be maximally mitigated.

### Docker mode
If you are using a proxy in your network, you need to configure your docker enviroment with the proxy. The configuration is presented in the Appendix section below. If not, you could directly build the docker image and use it with docker.

#### Build docker image
```console
cd /path/to/your/code/folder
docker build -t crl/ggr-recon:latest .
```
In general, the tag is set in the form of foo/bar:x.x.x_\*.\*.\*, where x.x.x denotes the docker image version and \*.\*.\* denotes the software version, e.g., crl/ggr:0.0.1_0.9.0

### Python mode
#### Configure python environment
Run the following command to install the python libraries
```console
python -m pip install -r ./requirements.txt
```
Install the CRKIT software suit and configure it in your *~/.bashrc*
```console
wget http://crl.med.harvard.edu/CRKIT/CRKIT-1.6.0-RHEL6.tar.gz
tar -xf CRKIT-1.6.0-RHEL6.tar.gz

export BUNDLE=/path/to/crkit/crkit-1.6.0
export PATH=$PATH:$BUNDLE/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BUNDLE/Frameworks/InsightToolkit:$BUNDLE/Frameworks/vtk-6.1:$BUNDLE/Frameworks/qt-5.3.2/lib:$BUNDLE/lib:$BUNDLE/bin
export QT_PLUGIN_PATH=$BUNDLE/Frameworks/qt-5.3.2/plugins
export DYLD_LIBRARY_PATH=""
```

## Usage
### View help
#### Docker mode
```console
docker run -it --rm --name ggr crl/ggr-recon preprocess.py -h
```
#### Docker mode
```console
docker run -it --rm --name ggr crl/ggr-recon recon.py -h
```

#### Docker mode
```console
docker run -it --rm --name ggr crl/ggr-recon pipeline.py -h
```

### Repository example data
This repository ships two example layouts with the same phantom images:

```text
data/
├── acr-axial/
│   ├── ax_t2_phantom.nii.gz
│   └── ax_t2_phantom.json
├── acr-coronal/
│   ├── cor_t2_phantom.nii.gz
│   └── cor_t2_phantom.json
└── acr-sagittal/
    ├── sag_t2_phantom.nii.gz
    └── sag_t2_phantom.json

data-bids-example/
├── dataset_description.json
└── sub-phantom/
    └── ses-01/
        └── anat/
            ├── sub-phantom_ses-01_acq-ax_T2w.nii.gz
            ├── sub-phantom_ses-01_acq-ax_T2w.json
            ├── sub-phantom_ses-01_acq-cor_T2w.nii.gz
            ├── sub-phantom_ses-01_acq-cor_T2w.json
            ├── sub-phantom_ses-01_acq-sag_T2w.nii.gz
            └── sub-phantom_ses-01_acq-sag_T2w.json
```

### Input and output
`preprocess.py` supports:

1. BIDS auto-discovery using `pybids` (`--path` + optional `--bids-filter KEY=VALUE`)
2. Explicit input files with `-f/--filenames` (legacy non-BIDS mode)

For BIDS mode, pass the dataset root with `--path`. A reconstructable group requires a complete `acq-{sag,cor,ax}` set with suffix `T2w` (`.nii` or `.nii.gz`) within a single BIDS entity group. Additional entities (`run`, `rec`, `desc`, and others) are supported and used for grouping.

Intermediate files are written to `--temp_path` (default: `/opt/GGR-recon/temp/`).
Final outputs are written under `--out_path` in the matching anatomical folder:

`sub-xxx[/ses-xx]/anat/`

Output naming:

- `acq-*` and existing `rec-*` are removed from the input stem
- `rec-superesolution` is inserted before `_T2w`
- a sidecar JSON is created next to the output NIfTI

Example:

- input: `sub-001_ses-01_run-1_acq-sag_T2w.nii.gz`
- output: `sub-001_ses-01_run-1_rec-superesolution_T2w.nii.gz`

If source inputs are split by `rec-*` (for example `rec-filtered` and `rec-orig`), the source `rec` value is copied into `acq-*` to avoid output collisions:

- `..._acq-sag_rec-filtered_T2w.nii.gz` -> `..._acq-filtered_rec-superesolution_T2w.nii.gz`
- `..._acq-sag_rec-orig_T2w.nii.gz` -> `..._acq-orig_rec-superesolution_T2w.nii.gz`

Run preprocessing only:

```console
docker run -it --rm --name ggr-recon \
  -v /your/bids:/bids \
  -v /your/temp:/temp \
  crl/ggr-recon preprocess.py \
  --path /bids --temp_path /temp --out_path /bids \
  --bids-filter subject=001 --bids-filter session=01
```

Run reconstruction only:

```console
docker run --rm -it \
  -v /your/temp:/temp \
  -v /your/bids:/bids \
  crl/ggr-recon recon.py \
  --temp_path /temp --out_path /bids --ggr -w 0.03
```

Run preprocessing + reconstruction with `pipeline.py`:

```console
docker run --rm -it \
  -v /your/bids:/bids \
  -v /your/temp:/temp \
  crl/ggr-recon pipeline.py \
  --path /bids --temp_path /temp --out_path /bids \
  --bids-filter subject=2983 --bids-filter session=1a \
  -- --ggr -w 0.03
```

`pipeline.py` also supports direct BIDS-App style selectors:

- `--participant-label` (maps to subject filter)
- `--session-label` (maps to session filter)

These can be used instead of `--bids-filter subject=...` and `--bids-filter session=...`:

```console
docker run --rm -it \
  -v /your/bids:/bids \
  -v /your/temp:/temp \
  crl/ggr-recon pipeline.py \
  --path /bids --temp_path /temp --out_path /bids \
  --participant-label 2983 --session-label 1a \
  -- --ggr -w 0.03
```

Both label options accept multiple values (space-separated and/or repeated). If label options and equivalent `--bids-filter` entries are both provided, labels take precedence for subject/session. In explicit filename mode (`-f/--filenames`), labels are ignored and `pipeline.py` prints a warning.

When no explicit `-f/--filenames` is provided, `pipeline.py` reconstructs all complete BIDS groups that match your filters. If you do not pass a `rec` filter, all matching `rec-*` groups are processed.

### BIDS phantom example
Run the bundled BIDS example directly from this repository:

```console
docker run --rm -it \
  -v `pwd`/data-bids-example:/bids \
  -v `pwd`/temp:/temp \
  crl/ggr-recon pipeline.py \
  --path /bids --temp_path /temp --out_path /bids \
  --bids-filter subject=phantom --bids-filter session=01 \
  -- --ggr -w 0.03
```

This writes:

- `data-bids-example/sub-phantom/ses-01/anat/sub-phantom_ses-01_rec-superesolution_T2w.nii.gz`
- `data-bids-example/sub-phantom/ses-01/anat/sub-phantom_ses-01_rec-superesolution_T2w.json`

### Non-BIDS phantom example
The original `data/` folder is still usable with explicit `-f` inputs:

```console
docker run --rm -it --volume `pwd`:/data crl/ggr-recon \
  preprocess.py --temp_path /data/temp --out_path /data/recons \
  -f /data/data/acr-axial/ax_t2_phantom.nii.gz /data/data/acr-coronal/cor_t2_phantom.nii.gz /data/data/acr-sagittal/sag_t2_phantom.nii.gz

docker run --rm -it \
  --volume `pwd`/temp:/opt/GGR-recon/temp \
  --volume `pwd`/recons:/opt/GGR-recon/recons \
  crl/ggr-recon recon.py --temp_path /opt/GGR-recon/temp --out_path /opt/GGR-recon/recons --ggr -w 0.03
```

### Baseline implementation
In the **deconvolution** step, a total variation (TV) regularization is also implemented with the Tikhonov criterion, for the comparison to our gradient guidance regularization (GGR). To enable the TV regularization instead of GGR, use the option *--tik* when running *recon.py*

## Appendix
### Data acquisition protocol
We recommend acquiring three low-res images in the three complementary planes respectively. Each low-res image comprises high in-plane resolution and thick slices. For example, we acquire T2 TSE images with an in-plane resolution of 0.5mm x 0.5mm and thickness of 2mm, and reconstruct the high-res image at the isotropic resolution of 0.5mm. It takes two minutes to acquire such an image on our scanner. With this protocol, GGR-recon never enhances in-plane resoltuion but reduces slice thickness.

### Configure docker environment
The configuration for docker is mainly about the proxy setting. The proxy needs to be set by two steps: 1) declare in the Dockerfile and 2) set in the system level.

In the *Dockerfile*, add the environment variables by
```Docker
ENV http_proxy "http://your.proxy.edu:3128"
ENV no_proxy "127.0.0.1,localhost"
```
For the system level configuration, edit the file */etc/systemd/system/docker.service.d/http-proxy.conf* and insert the following lines
```
[Service]
Environment="HTTP_PROXY=http://your.proxy.edu:3128"
Environment="NO_PROXY=localhost,127.0.0.1"
```
Flush changes and restart Docker
```console
 sudo systemctl daemon-reload
 sudo systemctl restart docker
```
Verify that the configuration has been loaded and matches the changes you made, for example:
```console
sudo systemctl show --property=Environment docker
```
```
Environment=HTTP_PROXY=http://your.proxy.edu:3128 NO_PROXY=localhost,127.0.0.1
```

The proxy setting can also be set in the file *~/.docker/config.json* with the following lines
```js
{
  "proxies":
  {
    "default":
    {
      "httpProxy": "http://your.proxy.edu:3128",
      "noProxy": "127.0.0.1,localhost"
    }
  }
}
```

Then, build your docker image with the proxy arguemnt by
```
docker build --build-arg http_proxy=http://your.proxy.edu:3128 -t your-ggr-tag .
```

or build your docker image with the build computer network:
```
docker build --network=host  -t crl/ggr-recon .
```


Build command without proxy:
```
docker build  -t crl/ggr-recon .
```


## References
  1. Yao Sui, Onur Afacan, Ali Gholipour, and Simon K. Warfield. 2019. “**Isotropic MRI Super-Resolution Reconstruction with Multi-Scale Gradient Field Prior**.” *International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI)*. Shen Zhen, China. <a href="https://scholar.harvard.edu/files/suiyao/files/sui_miccai_2019.pdf">PDF</a>
  2. Yao Sui, Onur Afacan, Ali Gholipour, and Simon K. Warfield. 2020. “**Learning a Gradient Guidance for Spatially Isotropic MRI Super-Resolution Reconstruction**.” *International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI)*. Lima, Peru. <a href="https://scholar.harvard.edu/files/suiyao/files/sui_deepgg_miccai_2020.pdf">PDF</a>


## Applications
1. "**Super-resolution reconstruction of T2-weighted thick-slice neonatal brain MRI scans.**" Incebak et al. 2022. PMID: 34506677 PMCID: PMC8752487 DOI: 10.1111/jon.12929  https://pubmed.ncbi.nlm.nih.gov/34506677/

   This paper examined the use of this type of super-resolution for imaging of neonates. Qualitative and quantitative assessments showed that 3D SRR of several LR images produces images that are of comparable quality to standard 2D HR image acquisition for healthy neonatal imaging without loss of anatomical details with similar edge definition allowing the detection of fine anatomical structures and permitting comparable morphometric measurement.
