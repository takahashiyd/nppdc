# Data

| Folder | File | Source | Used by |
|---|---|---|---|
| `allen/` | `subset_movie.mat` (37 MB): 8 V1 + 8 LM LFP channels, 1250 Hz, 600 s of `natural_movie_three`, Allen Brain Observatory Neuropixels session 715093703 | extracted with `demos/extract_allen_subset.m` from the public session and probe NWB files | `demos/allen_lfp_analysis.m`, `demos/hodge_demo.m` |
| `allen/` | `sessions.csv`, `probes.csv`, `channels.csv` | Allen metadata tables | reference for the session, probe and channel identifiers |

The raw Allen data are public; only this small subset is redistributed here so
that every figure can be regenerated without a download.

## Source, attribution and terms of use

The files in `allen/` are derived from the **Allen Brain Observatory – Visual
Coding Neuropixels** dataset, © Allen Institute for Brain Science. They are
redistributed for reproducibility under the Allen Institute
[Terms of Use](https://alleninstitute.org/legal/terms-of-use) and are **not**
covered by the MIT license of the code. If you use them, follow the Allen
[Citation Policy](https://alleninstitute.org/legal/citation-policy) and cite:

> Siegle, J. H., Jia, X., Durand, S., et al. (2021). Survey of spiking in the
> mouse visual system reveals functional hierarchy. *Nature*, 592, 86–92.
> https://doi.org/10.1038/s41586-020-03171-x

## Rebuilding `subset_movie.mat` from the raw NWB files

`demos/extract_allen_subset.m` reads three NWB files (about 2–3 GB each) from the
Allen public S3 bucket. Download them into `data/allen/` **under the local names
the script expects**:

| Local name expected by the script | Allen file | Probe |
|---|---|---|
| `session_715093703.nwb` | `session_715093703.nwb` | — |
| `probe_C_V1_lfp.nwb` | `probe_810755801_lfp.nwb` | probeC (V1, VISp) |
| `probe_D_LM_lfp.nwb` | `probe_810755803_lfp.nwb` | probeD (LM, VISl) |

```sh
base=https://allen-brain-observatory.s3.us-west-2.amazonaws.com/visual-coding-neuropixels/ecephys-cache/session_715093703
cd data/allen
curl -o session_715093703.nwb $base/session_715093703.nwb
curl -o probe_C_V1_lfp.nwb   $base/probe_810755801_lfp.nwb
curl -o probe_D_LM_lfp.nwb   $base/probe_810755803_lfp.nwb
```

Then run `demos/extract_allen_subset.m` from MATLAB. It uses only the built-in
`h5read` (no MatNWB needed). The probe IDs come from `probes.csv`
(`ecephys_session_id = 715093703`).
