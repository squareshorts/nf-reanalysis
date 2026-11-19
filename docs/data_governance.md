# Data governance

The public repository intentionally excludes:

- raw EEG recordings (`.edf`, `.mat`, and related source files);
- the participant-level derived EEG table used by the scripts;
- BrainAvatar score workbooks;
- cognitive-test workbooks and operational records;
- participant crosswalks, local ID inventories, and source-document excerpts;
- individual-level PCA coordinates and the PCA score plot.

The derived EEG table contains participant-level rows and must remain in authorized local storage. The analysis scripts accept an explicit local path and write aggregate tables. The local output directory and PCA score plot should not be committed.

Only aggregate results are included in `results/aggregate/`. The corresponding author may provide the derived table upon reasonable request, subject to institutional approval and applicable data-governance requirements. The raw recordings remain restricted by institutional and occupational confidentiality constraints.
